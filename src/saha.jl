# Read ionization potentials (cm⁻¹) for each charge state of an element.
function ionization_potentials(db::LIBSDB, element::AbstractString)
    tbl = get(db.ionization, element, nothing)
    tbl === nothing && return Dict{Int,Float64}()
    ips = Dict{Int,Float64}()
    for row in Tables.rows(tbl)
        ips[row.spectr_charge] = row.energy
    end
    ips
end

# Compute Saha-LTE fractional ion populations and partition functions.
#
# Returns Dict{sc => (abundance_fraction, partition_function)}.
#
# The Saha equation:
#   n_{z+1} n_e   (2π m_e k_B T)^(3/2)   2 U_{z+1}(T)
#   ───────────  = ─────────────────── × ──────────── × exp(−χ_z / k_B T)
#       n_z                h³               U_z(T)
#
# χ_z = ionization potential of charge state z.
# U_z(T) = partition function.
#
# The calculation is done in log-space to handle the extreme dynamic range.
function saha_ion_populations(db::LIBSDB, element::AbstractString, temp_eV::Real, eden::Real;
    charge_min::Integer=0, charge_max=nothing)

    eden > 0 || throw(ArgumentError("eden must be positive"))
    temp_eV > 0 || throw(ArgumentError("temp_eV must be positive"))

    ips = ionization_potentials(db, element)
    isempty(ips) && return Dict{Int,Tuple{Float64,Float64}}()

    max_avail = maximum(keys(ips); init=0)
    if charge_max === nothing
        charge_max = max_avail
    end
    # Ensure at least 3 charge states for normalisation to account for
    # the neutral → singly-ionised → doubly-ionised ladder even when
    # DB only has data for lower stages.
    charge_max = max(charge_max, 3)

    # Saha constant: (2π m_e k_B / h²)^(3/2) × T^(3/2) / n_e
    # ≈ 6.043 × 10²¹ × T_eV^(3/2) × n_e⁻¹  in cm⁻³-based units
    saha_factor = 6.043e21 * temp_eV^1.5 / eden

    # Pre-fetch partition functions
    pfs = Dict{Int,Float64}()
    for sc in charge_min:charge_max
        pfs[sc] = partition_function(db, element, sc, temp_eV)
    end

    # Forward iteration through charge states using Saha ratios
    # log(n_{sc}/n_0) = Σ_{i=0}^{sc-1} log((2π m_e k_B T/h²)^(3/2) / n_e)
    #                   + log(U_{i+1}/U_i) − χ_i / k_B T
    saha = Dict{Int,Tuple{Float64,Float64}}()
    prev_log_abund = 0.0
    prev_pf = 1.0
    prev_ip = 0.0
    max_log = -Inf

    temp_cm = temp_eV * CM_EV
    for sc in charge_min:charge_max
        pf = get(pfs, sc, 0.0)
        pf <= 0 && continue
        ip = get(ips, sc, 0.0)
        ratio = saha_factor / prev_pf
        log_abund = log(pf * ratio) + prev_log_abund - prev_ip / temp_cm

        saha[sc] = (log_abund, pf)
        max_log = max(max_log, log_abund)

        prev_log_abund = log_abund
        prev_pf = pf
        prev_ip = ip

        # When sc == Z (neutral → singly-ionised transition at the atomic
        # number boundary), also insert the next (unexpected) charge state
        # with pf=1.0 so the Saha ladder continues past the neutral atom.
        z = element_number(element)
        if z > 0 && z == sc
            sc_next = z + 1
            log_abund = log(saha_factor / prev_pf) + prev_log_abund - prev_ip / temp_cm
            saha[sc_next] = (log_abund, 1.0)
            max_log = max(max_log, log_abund)
            prev_log_abund = log_abund
            prev_pf = 1.0
            prev_ip = 0.0
        end
    end

    # Normalise so abundances sum to 1
    total = sum(exp(v[1] - max_log) for v in values(saha))
    log_total = log(total)
    result = Dict{Int,Tuple{Float64,Float64}}()
    for (k, v) in saha
        result[k] = (exp(v[1] - max_log - log_total), v[2])
    end
    result
end
