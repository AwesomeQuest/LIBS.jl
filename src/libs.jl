# ---------------------------------------------------------------------------
# Top-level LTE spectrum synthesis.
#
# `lte_spectrum_sticks` is the primary public API: given a spectrum
# specification string, a temperature, and an electron density, it returns
# a vector of `LIBSStickLine` sorted by wavelength.
#
# `lte_spectrum_data` wraps sticks with a Doppler-broadened continuum
# into a `SpectrumOverlay` for easy plotting.
# ---------------------------------------------------------------------------

"""
    SpectrumOverlay

Combines raw stick lines and a Doppler-broadened spectrum for plotting.
Fields:
- `sticks`   — `Vector{LIBSStickLine}`  (discrete transitions)
- `spectrum` — `Vector{DopplerGridPoint}`  (broadened continuum)
"""
struct SpectrumOverlay
    sticks::Vector{LIBSStickLine}
    spectrum::Vector{DopplerGridPoint}
end

"""
    lte_spectrum_sticks(spectra, temp, density; composition=nothing,
                        int_scale=0, min_rel_int=nothing,
                        low_wl=nothing, upp_wl=nothing, db=nothing)

Compute LTE line stick spectrum for one or more elements/charge stages.

# Arguments
- `spectra`   — spectrum specification string (e.g. "Fe I", "Fe I-III, Ni I",
                "Fe0-2", "Fe I, Ni I").  See `parse_spectra` for full syntax.
- `temp`      — plasma temperature: `Quantity` in eV or K, or bare Float64 (eV)
- `density`   — electron number density: `Quantity` in cm⁻³, or bare Float64
- `composition` — optional `Dict{String, Float64}` of elemental abundances
                (e.g. `Dict("Fe" => 0.7, "Ni" => 0.3)`).  Defaults to equal
                abundance for all elements in `spectra`.
- `int_scale` — if 1, weight intensities by photon energy (E_u − E_l)
- `min_rel_int` — discard lines with intensity < min_rel_int × max intensity
                per charge stage (e.g. 0.01 to keep only strong lines).
                `nothing` keeps all lines.
- `low_wl`    — lower wavelength bound (`Quantity` or nothing)
- `upp_wl`    — upper wavelength bound (`Quantity` or nothing)
- `db`        — optional `LIBSDB` instance (uses the singleton if not provided)

# Returns
`Vector{LIBSStickLine}` sorted by wavelength.  Each stick has:
- `wavelength` as a `Quantity` in the inferred output unit (default nm)
- `intensity` in arbitrary relative units

# Algorithm
1. Parse the spectrum specification into element/charge entries.
2. For each element, compute Saha–Boltzmann ionisation populations.
3. For each transition in the database matching the requested elements,
   charges, and wavelength bounds, compute the LTE line intensity.
4. Optionally filter by minimum relative intensity and sort by wavelength.
"""
function lte_spectrum_sticks(spectra::AbstractString, temp, density;
    composition=nothing, int_scale::Integer=0, min_rel_int=nothing,
    low_wl=nothing, upp_wl=nothing, db=nothing)

    db = unwrap_or(db, open_db())
    entries = parse_spectra(spectra)
    isempty(entries) && return LIBSStickLine[]

    temp_eV = to_internal_temp(temp)
    eden = to_internal_density(density)

    # Expand bare element entries (empty charges) to all available
    # charge stages in the database.
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

    # Cache Saha populations per element across shared charge stages
    saha_cache = Dict{String,Dict{Int,Tuple{Float64,Float64}}}()
    for e in all_entries
        elem = element_symbol(e.Z)
        min_c = minimum(e.charges)
        max_c = maximum(e.charges)
        saha_cache[elem] = saha_ion_populations(db, elem, temp_eV, eden; charge_min=min_c, charge_max=max_c)
    end

    # Normalise elemental abundances
    elem_names = unique(element_symbol(e.Z) for e in all_entries)
    if composition === nothing
        composition = Dict(name => 1.0 / length(elem_names) for name in elem_names)
    end
    total_frac = sum(values(composition); init=0.0)
    total_frac <= 0 && throw(ArgumentError("composition fractions must be positive"))
    if total_frac != 1.0
        # Renormalise to sum = 1
        composition = Dict(k => v / total_frac for (k, v) in composition)
    end

    # Convert wavelength bounds to internal Å for DB filtering
    low_wl_A = low_wl === nothing ? nothing : to_internal_length(low_wl)
    upp_wl_A = upp_wl === nothing ? nothing : to_internal_length(upp_wl)

    # Determine the output unit from the input bounds (default: nm)
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

            # Filter by wavelength bounds (vacuum wavelengths in Å)
            vwl = unwrap_or(row.vac_wl_num, 0.0)
            if low_wl_A !== nothing && vwl < low_wl_A
                continue
            end
            if upp_wl_A !== nothing && vwl > upp_wl_A
                continue
            end

            # Compute LTE line intensity
            intensity = lte_line_intensity(row, saha, temp_eV; int_scale=int_scale, abundance=abundance)
            intensity <= 0 && continue

            # Determine the output wavelength:
            #   calc_wl → theoretical (HFR/MCDF) wavelength (most accurate)
            #   vac_wl  → tabulated vacuum wavelength (fallback)
            #   if the database marks wl_in_air, convert vac_wl → air wavelength
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

    # Optional filtering: keep only lines within a factor of min_rel_int
    # of the maximum-intensity line in each charge stage.
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

"""
    lte_spectrum_data(spectra, temp, density, resolution; kwargs...)

Compute a full LTE spectrum (sticks + Doppler broadened) in one call.
Shorthand for:
    sticks = lte_spectrum_sticks(spectra, temp, density; kwargs...)
    spectrum = doppler_spectrum(sticks, resolution)
    SpectrumOverlay(sticks, spectrum)

Keyword arguments are forwarded to `lte_spectrum_sticks`.
"""
function lte_spectrum_data(spectra::AbstractString, temp, density, resolution::Real; kwargs...)
    sticks = lte_spectrum_sticks(spectra, temp, density; kwargs...)
    spectrum = doppler_spectrum(sticks, resolution)
    SpectrumOverlay(sticks, spectrum)
end

# Convenience methods accepting a pre-opened database
lte_spectrum_sticks(db::LIBSDB, args...; kwargs...) =
    lte_spectrum_sticks(args...; db=db, kwargs...)

lte_spectrum_data(db::LIBSDB, args...; kwargs...) =
    lte_spectrum_data(args...; db=db, kwargs...)
