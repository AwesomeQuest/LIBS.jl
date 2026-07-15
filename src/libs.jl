struct SpectrumOverlay
    sticks::Vector{LIBSStickLine}
    spectrum::Vector{DopplerGridPoint}
end

const WL_UNIT_TO_ANGSTROM = [1.0, 0.1, 0.0001]

function to_vacuum_angstrom(wl::Real, unit::Int, show_av::Int)
    unit = clamp(unit, 0, 2)
    wl_A = wl / WL_UNIT_TO_ANGSTROM[unit + 1]
    in_air = if show_av == 0
        wl_A >= 2000.0 && wl_A <= 10000.0
    elseif show_av == 2
        wl_A >= 2000.0 && wl_A <= 20000.0
    elseif show_av == 4
        wl_A >= 1850.0
    else
        false
    end
    return in_air ? air_to_vac(wl_A) : wl_A
end

function lte_spectrum_sticks(spectra::AbstractString, temp_eV::Real, eden::Real;
    composition::Union{Dict{String,Float64},Nothing}=nothing,
    int_scale::Int=0,
    min_rel_int::Union{Real,Nothing}=nothing,
    unit::Int=0,
    show_av::Int=2,
    low_wl::Union{Real,Nothing}=nothing,
    upp_wl::Union{Real,Nothing}=nothing,
    db::Union{LIBSDB,Nothing}=nothing)

    db = _v(db, open_db())
    entries = parse_spectra(spectra)
    isempty(entries) && return LIBSStickLine[]

    # Expand bare elements to all available charges for that element
    all_entries = SpectrumEntry[]
    for e in entries
        if isempty(e.charges)
            elem = element_symbol(e.Z)
            tbl = get(db.lines, elem, nothing)
            if tbl === nothing
                push!(all_entries, SpectrumEntry(e.Z, e.isotope, [0]))
                continue
            end
            charges = Int[]
            for row in Tables.rows(tbl)
                if !(row.spectr_charge in charges)
                    push!(charges, row.spectr_charge)
                end
            end
            isempty(charges) && push!(charges, 0)
            push!(all_entries, SpectrumEntry(e.Z, e.isotope, sort(charges)))
        else
            push!(all_entries, e)
        end
    end

    # Compute Saha populations for each unique element
    saha_cache = Dict{String,Dict{Int,Tuple{Float64,Float64}}}()
    for e in all_entries
        elem = element_symbol(e.Z)
        min_c = minimum(e.charges)
        max_c = maximum(e.charges)
        saha_cache[elem] = saha_ion_populations(db, elem, temp_eV, eden; charge_min=min_c, charge_max=max_c)
    end

    # Composition fractions
    elem_names = unique([element_symbol(e.Z) for e in all_entries])
    if composition === nothing
        composition = Dict(name => 1.0 / length(elem_names) for name in elem_names)
    end
    total_frac = sum(v for v in values(composition); init=0.0)
    total_frac <= 0 && throw(ArgumentError("composition fractions must be positive"))
    composition = Dict(k => v / total_frac for (k, v) in composition)

    low_wl_A = low_wl === nothing ? nothing : to_vacuum_angstrom(low_wl, unit, show_av)
    upp_wl_A = upp_wl === nothing ? nothing : to_vacuum_angstrom(upp_wl, unit, show_av)
    unit_clamped = clamp(unit, 0, 2)
    wl_unit_factor = WL_UNIT_TO_ANGSTROM[unit_clamped + 1]

    # Collect lines per element from per-element Arrow files
    sticks = LIBSStickLine[]
    for e in all_entries
        elem = element_symbol(e.Z)
        abundance = get(composition, elem, 0.0)
        abundance <= 0 && continue
        saha = get(saha_cache, elem, Dict{Int,Tuple{Float64,Float64}}())
        isempty(saha) && continue

        charge_set = Set(e.charges)
        tbl = get(db.lines, elem, nothing)
        tbl === nothing && continue

        for row in Tables.rows(tbl)
            row.spectr_charge in charge_set || continue
            vwl = _v(row.vac_wl_num, 0.0)
            if low_wl_A !== nothing && vwl < low_wl_A
                continue
            end
            if upp_wl_A !== nothing && vwl > upp_wl_A
                continue
            end

            intensity = lte_line_intensity(row, saha, temp_eV; int_scale=int_scale, abundance=abundance)
            intensity <= 0 && continue

            calc_wl = _v(row.calc_wl_num, 0.0)
            vac_wl = _v(row.vac_wl_num, 0.0)
            wl_in_air = _v(row.wl_in_air, false)

            if calc_wl > 0
                wl_angstrom = calc_wl
            elseif wl_in_air == true
                wl_angstrom = vac_to_air(vac_wl)
            else
                wl_angstrom = vac_wl
            end
            wl_angstrom > 0 || continue

            wl = wl_angstrom * wl_unit_factor
            push!(sticks, LIBSStickLine(
                elem, row.spectr_charge, wl, intensity,
                _v(row.low_conf, ""),
                _v(row.low_term, ""),
                _v(row.upp_conf, ""),
                _v(row.upp_term, "")
            ))
        end
    end

    if min_rel_int !== nothing && min_rel_int > 0
        max_per_charge = Dict{Int,Float64}()
        for s in sticks
            cur = get(max_per_charge, s.spectr_charge, 0.0)
            if s.intensity > cur
                max_per_charge[s.spectr_charge] = s.intensity
            end
        end
        filter!(s -> s.intensity >= min_rel_int * get(max_per_charge, s.spectr_charge, 0.0), sticks)
    end
    sort!(sticks, by=s -> s.wavelength)
    sticks
end

function lte_spectrum_data(spectra::AbstractString, temp_eV::Real, eden::Real, resolution::Real;
    composition::Union{Dict{String,Float64},Nothing}=nothing,
    int_scale::Int=0,
    min_rel_int::Union{Real,Nothing}=nothing,
    unit::Int=0,
    show_av::Int=2,
    low_wl::Union{Real,Nothing}=nothing,
    upp_wl::Union{Real,Nothing}=nothing,
    db::Union{LIBSDB,Nothing}=nothing,
    kwargs...)

    sticks = lte_spectrum_sticks(spectra, temp_eV, eden;
        composition=composition, int_scale=int_scale, min_rel_int=min_rel_int,
        unit=unit, show_av=show_av, low_wl=low_wl, upp_wl=upp_wl, db=db, kwargs...)
    spectrum = doppler_spectrum(sticks, resolution)
    return SpectrumOverlay(sticks, spectrum)
end

function lte_spectrum_sticks(db::LIBSDB, spectra::AbstractString, temp_eV::Real, eden::Real; kwargs...)
    lte_spectrum_sticks(spectra, temp_eV, eden; db=db, kwargs...)
end

function lte_spectrum_data(db::LIBSDB, spectra::AbstractString, temp_eV::Real, eden::Real, resolution::Real; kwargs...)
    lte_spectrum_data(spectra, temp_eV, eden, resolution; db=db, kwargs...)
end
