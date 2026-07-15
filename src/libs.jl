struct SpectrumOverlay
    sticks::Vector{LIBSStickLine}
    spectrum::Vector{DopplerGridPoint}
end

function lte_spectrum_sticks(spectra::AbstractString, temp, density;
    composition=nothing, int_scale::Integer=0, min_rel_int=nothing,
    low_wl=nothing, upp_wl=nothing, db=nothing)

    db = unwrap_or(db, open_db())
    entries = parse_spectra(spectra)
    isempty(entries) && return LIBSStickLine[]

    temp_eV = to_internal_temp(temp)
    eden = to_internal_density(density)

    # Expand bare elements to all available charges
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

    saha_cache = Dict{String,Dict{Int,Tuple{Float64,Float64}}}()
    for e in all_entries
        elem = element_symbol(e.Z)
        min_c = minimum(e.charges)
        max_c = maximum(e.charges)
        saha_cache[elem] = saha_ion_populations(db, elem, temp_eV, eden; charge_min=min_c, charge_max=max_c)
    end

    elem_names = unique(element_symbol(e.Z) for e in all_entries)
    if composition === nothing
        composition = Dict(name => 1.0 / length(elem_names) for name in elem_names)
    end
    total_frac = sum(values(composition); init=0.0)
    total_frac <= 0 && throw(ArgumentError("composition fractions must be positive"))
    if total_frac != 1.0
        composition = Dict(k => v / total_frac for (k, v) in composition)
    end

    # Wavelength bounds → internal Å for DB filtering
    low_wl_A = low_wl === nothing ? nothing : to_internal_length(low_wl)
    upp_wl_A = upp_wl === nothing ? nothing : to_internal_length(upp_wl)

    # Output unit inferred from input bounds
    out_unit = output_unit(low_wl, upp_wl)

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
            vwl = unwrap_or(row.vac_wl_num, 0.0)
            if low_wl_A !== nothing && vwl < low_wl_A
                continue
            end
            if upp_wl_A !== nothing && vwl > upp_wl_A
                continue
            end

            intensity = lte_line_intensity(row, saha, temp_eV; int_scale=int_scale, abundance=abundance)
            intensity <= 0 && continue

            calc_wl = unwrap_or(row.calc_wl_num, 0.0)
            vac_wl = unwrap_or(row.vac_wl_num, 0.0)

            if calc_wl > 0
                wl_angstrom = calc_wl
            elseif unwrap_or(row.wl_in_air, false)
                wl_angstrom = vac_to_air(vac_wl)
            else
                wl_angstrom = vac_wl
            end
            wl_angstrom > 0 || continue

            push!(sticks, LIBSStickLine(
                elem, row.spectr_charge, wrap_output(wl_angstrom, out_unit), intensity,
                unwrap_or(row.low_conf, ""),
                unwrap_or(row.low_term, ""),
                unwrap_or(row.upp_conf, ""),
                unwrap_or(row.upp_term, "")
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
    sort!(sticks)
end

function lte_spectrum_data(spectra::AbstractString, temp, density, resolution::Real; kwargs...)
    sticks = lte_spectrum_sticks(spectra, temp, density; kwargs...)
    spectrum = doppler_spectrum(sticks, resolution)
    SpectrumOverlay(sticks, spectrum)
end

lte_spectrum_sticks(db::LIBSDB, args...; kwargs...) =
    lte_spectrum_sticks(args...; db=db, kwargs...)

lte_spectrum_data(db::LIBSDB, args...; kwargs...) =
    lte_spectrum_data(args...; db=db, kwargs...)
