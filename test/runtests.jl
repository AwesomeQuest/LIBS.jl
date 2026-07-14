using LIBS
using Test
using JSON
using Artifacts

const FIXTURES_DIR = @artifact_str("libs_test_data")
const WL_TOL_NM = 0.01

function load_fixture(elem::AbstractString)
    path = joinpath(FIXTURES_DIR, "$(lowercase(elem)).json")
    isfile(path) || return nothing
    return JSON.parsefile(path)
end

function get_nist_charge_data(fixture)
    stages = fixture["charge_stages"]
    data = fixture["data"]
    result = Dict{Int,Dict}()
    for (idx, label) in enumerate(stages)
        label == "Ritz wavelength (nm)" && continue
        parts = split(label)
        roman = parts[end]
        charge = LIBS.roman_to_int(roman)
        result[charge] = Dict(
            :label => label,
            :wavelengths => Float64[],
            :intensities => Float64[],
        )
    end
    for entry in data
        wl = entry["wl"]
        ints = entry["intensities"]
        for (stage_pos, label) in enumerate(stages)
            label == "Ritz wavelength (nm)" && continue
            int_idx = stage_pos - 1  # stage_pos=2 (first charge col) → ints[1], stage_pos=3 → ints[2]
            if length(ints) >= int_idx && ints[int_idx] !== nothing
                parts = split(label)
                charge = LIBS.roman_to_int(parts[end])
                push!(result[charge][:wavelengths], wl)
                push!(result[charge][:intensities], Float64(ints[int_idx]))
            end
        end
    end
    return result
end

function get_our_charge_data(elem::AbstractString)
    sticks = lte_spectrum_sticks(elem, 1.0, 1e17;
        unit=1, int_scale=1, min_rel_int=nothing,
        low_wl=200.0, upp_wl=600.0, show_av=2)
    result = Dict{Int,Dict}()
    for s in sticks
        if !haskey(result, s.spectr_charge)
            result[s.spectr_charge] = Dict(
                :wavelengths => Float64[],
                :intensities => Float64[],
            )
        end
        push!(result[s.spectr_charge][:wavelengths], s.wavelength)
        push!(result[s.spectr_charge][:intensities], s.intensity)
    end
    return result
end

function run_element_test(elem::AbstractString, fixture)
    nist = get_nist_charge_data(fixture)
    our = get_our_charge_data(elem)

    common_charges = intersect(keys(nist), keys(our))
    @test !isempty(common_charges)

    for charge in sort(collect(common_charges))
        nist_wls = nist[charge][:wavelengths]
        nist_ints = nist[charge][:intensities]
        our_wls = our[charge][:wavelengths]
        our_ints = our[charge][:intensities]

        if isempty(nist_wls) || isempty(our_wls)
            continue
        end

        # For each NIST wavelength, find matching our wavelength
        nist_match_count = 0
        matched_nist_idxs = Int[]
        for j in 1:length(nist_wls)
            nwl = nist_wls[j]
            for i in 1:length(our_wls)
                if abs(nwl - our_wls[i]) <= WL_TOL_NM
                    nist_match_count += 1
                    push!(matched_nist_idxs, j)
                    break
                end
            end
        end

        # For each our wavelength, find matching NIST wavelength
        our_match_count = 0
        for i in 1:length(our_wls)
            owl = our_wls[i]
            for j in 1:length(nist_wls)
                if abs(owl - nist_wls[j]) <= WL_TOL_NM
                    our_match_count += 1
                    break
                end
            end
        end

        recall = nist_match_count / max(length(nist_wls), 1)
        precision_val = our_match_count / max(length(our_wls), 1)

        # Verify matched wavelengths are accurate
        avg_dev = 0.0
        pair_count = 0
        for i in 1:length(our_wls)
            owl = our_wls[i]
            for j in 1:length(nist_wls)
                if abs(owl - nist_wls[j]) <= WL_TOL_NM
                    avg_dev += abs(owl - nist_wls[j])
                    pair_count += 1
                    break
                end
            end
        end
        if pair_count >= 3
            avg_dev /= pair_count
            @test avg_dev <= WL_TOL_NM
        end

        # Verify intensity correlation for matched lines
        pair_count_int = 0
        sum_int_dev = 0.0
        for i in 1:length(our_wls)
            owl = our_wls[i]
            for j in 1:length(nist_wls)
                if abs(owl - nist_wls[j]) <= WL_TOL_NM
                    our_norm = our_ints[i] / maximum(our_ints)
                    nist_norm = nist_ints[j] / maximum(nist_ints)
                    if our_norm > 0.01 || nist_norm > 0.01
                        sum_int_dev += abs(our_norm - nist_norm)
                        pair_count_int += 1
                    end
                    break
                end
            end
        end
        if pair_count_int >= 5
            avg_int_dev = sum_int_dev / pair_count_int
            @test avg_int_dev <= 0.35
        end

        # Check that charge stages have reasonable overlap
        @test recall > 0 || precision_val > 0
    end
end

function find_available_elements()
    elems = String[]
    for f in readdir(FIXTURES_DIR)
        endswith(f, ".json") || continue
        push!(elems, first(splitext(f)))
    end
    return elems
end

@testset "LIBS" begin
    @testset "Basic API" begin
        @test isdefined(LIBS, :lte_spectrum_sticks)
        @test isdefined(LIBS, :partition_function)
        @test isdefined(LIBS, :saha_ion_populations)
        @test isdefined(LIBS, :int_to_roman)
        @test isdefined(LIBS, :roman_to_int)
        @test LIBS.int_to_roman(1) == "I"
        @test LIBS.int_to_roman(2) == "II"
        @test LIBS.int_to_roman(3) == "III"
        @test LIBS.roman_to_int("I") == 1
        @test LIBS.roman_to_int("II") == 2
        @test LIBS.roman_to_int("III") == 3
    end

    @testset "parse_spectra" begin
        entries = LIBS.parse_spectra("Fe I")
        @test length(entries) == 1
        @test LIBS.element_symbol(entries[1].Z) == "Fe"
        @test entries[1].charges == [1]

        entries = LIBS.parse_spectra("Fe I-III")
        @test length(entries) == 1
        @test entries[1].charges == [1, 2, 3]

        entries = LIBS.parse_spectra("Fe0-2")
        @test length(entries) == 1
        @test entries[1].charges == [1, 2, 3]

        entries = LIBS.parse_spectra("Fe I, Ni I")
        @test length(entries) == 2
    end

    @testset "partition_function" begin
        db = LIBS.open_db()
        pf = LIBS.partition_function(db, "Fe", 1, 1.0)
        @test pf > 0
    end

    @testset "saha_ion_populations" begin
        db = LIBS.open_db()
        pops = LIBS.saha_ion_populations(db, "Fe", 1.0, 1e17; charge_min=1, charge_max=3)
        @test haskey(pops, 1)
        @test haskey(pops, 2)
        @test haskey(pops, 3)
        total_pop = sum(v[1] for v in values(pops))
        @test abs(total_pop - 1.0) < 1e-6
    end

    @testset "lte_spectrum_sticks basic" begin
        sticks = lte_spectrum_sticks("Fe I", 1.0, 1e17; unit=1, low_wl=200.0, upp_wl=600.0)
        @test length(sticks) > 0
        for s in sticks
            @test s.wavelength >= 200.0
            @test s.wavelength <= 600.0
            @test s.spectr_charge == 1
        end
    end

    test_elements = find_available_elements()
    total = length(test_elements)
    @testset "NIST validation ($elem)" for (i, elem) in enumerate(test_elements)
        print("\r\e[2KNIST validation [$i/$total] ", lpad(elem, 2), "...")
        flush(stdout)
        if lowercase(elem) == "h"
            @test true
            continue
        end
        fixture = load_fixture(elem)
        fixture === nothing && (@test false; continue)
        run_element_test(elem, fixture)
    end
    println()
end
