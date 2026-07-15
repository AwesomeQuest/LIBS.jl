# ---------------------------------------------------------------------------
# Saha ionisation and Saha-Boltzmann population fractions.
#
# The Saha equation relates the relative number densities of successive
# ionisation stages in a plasma in local thermodynamic equilibrium (LTE):
#
#     n_{z+1} · n_e   2 · Z_{z+1}(T) · (2π m_e k_B T)^{3/2}
#     ─────────────  =  ────────────────────────────────────── · exp(−χ_z / k_B T)
#         n_z                       Z_z(T) · h³
#
# where
#   n_z   = number density of ionisation stage z (z = 1 for neutral)
#   n_e   = electron number density (cm⁻³)
#   Z_z   = internal partition function of stage z
#   χ_z   = ionisation potential from stage z to z+1 (eV or cm⁻¹)
#   T     = temperature (eV internally)
#   m_e   = electron mass
#   k_B   = Boltzmann constant
#   h     = Planck constant
#
# In logarithmic form, the Saha decrement from stage z to z+1 is:
#
#     log(n_{z+1} / n_z) = log(S_T · Z_{z+1} / Z_z) − χ_z / k_B T
#
# with the Saha prefactor:
#     S_T = 2 · (2π m_e / h²)^{3/2} · (k_B T)^{3/2} / n_e
#         = 6.043 × 10²¹ · T_{eV}^{3/2} / n_e    [in CGS units]
#
# After computing the relative Boltzmann factors, populations are normalized
# so that Σ_z n_z = 1 (fractional abundances).
# ---------------------------------------------------------------------------

"""
    ionization_potentials(db, element)

Return a `Dict{Int, Float64}` mapping ionisation stage → ionisation potential
(in cm⁻¹) for a given element, as read from the NIST ionization tables.
"""
function ionization_potentials(db::LIBSDB, element::AbstractString)
    tbl = get(db.ionization, element, nothing)
    tbl === nothing && return Dict{Int,Float64}()
    ips = Dict{Int,Float64}()
    for row in Tables.rows(tbl)
        ips[row.spectr_charge] = row.energy
    end
    ips
end

"""
    saha_ion_populations(db, element, temp, eden; charge_min=0, charge_max=nothing)

Compute fractional ionisation stage populations for a plasma in LTE using the
Saha–Boltzmann equation.

# Arguments
- `db`        — `LIBSDB` instance
- `element`   — element name (e.g. "Fe")
- `temp`      — temperature (eV, K, or bare Float64 in eV)
- `eden`      — electron number density (cm⁻³ or compatible Quantity)
- `charge_min`– minimum ionisation stage to consider (default 0)
- `charge_max`– maximum stage (default: auto-detect from available IPs, min 3)

# Returns
A `Dict{Int, Tuple{Float64, Float64}}` mapping ionisation stage →
    (fractional_abundance, partition_function).
The abundances sum to 1.0.

# Algorithm
1. Retrieve partition functions Z_z(T) for each charge stage.
2. Compute the Saha factor S_T = 6.043e21 · T^{3/2} / n_e.
3. Recursively apply the Saha decrement from the lowest to the highest stage:
       log_abund(z+1) = log_abund(z) + log(S_T · Z_{z+1} / Z_z) − IP_z / k_B T
4. For z = Z (nuclear charge), insert a "bare nucleus" stage where Z = 1
   (fully stripped ion).
5. Normalize so that Σ exp(log_abund(z) − log_max) = 1.
"""
function saha_ion_populations(db::LIBSDB, element::AbstractString, temp, eden;
    charge_min::Integer=0, charge_max=nothing)

    temp_eV = to_internal_temp(temp)
    eden_cm3 = to_internal_density(eden)

    eden_cm3 > 0 || throw(ArgumentError("eden must be positive"))
    temp_eV > 0 || throw(ArgumentError("temp must be positive"))

    # Retrieve ionisation potentials (cm⁻¹) for each charge stage
    ips = ionization_potentials(db, element)
    isempty(ips) && return Dict{Int,Tuple{Float64,Float64}}()

    # Determine the maximum charge stage to consider
    max_avail = maximum(keys(ips); init=0)
    if charge_max === nothing
        charge_max = max_avail
    end
    # Ensure we have at least 3 stages for multi-stage LTE to work
    charge_max = max(charge_max, 3)

    # Saha prefactor: S_T = 2 · (2π m_e k_B / h²)^{3/2} · (k_B T)^{3/2} / n_e
    #   = 6.043 × 10²¹ · T_{eV}^{3/2} / n_e   (in CGS)
    saha_factor = 6.043e21 * temp_eV^1.5 / eden_cm3

    # Pre-compute partition functions for all requested stages
    pfs = Dict{Int,Float64}()
    for sc in charge_min:charge_max
        pfs[sc] = partition_function(db, element, sc, temp_eV)
    end

    saha = Dict{Int,Tuple{Float64,Float64}}()
    prev_log_abund = 0.0   # log-abundance of the previous (lower) stage
    prev_pf = 1.0          # Z_{z-1}  (Z of neutral atom = 1 by convention)
    prev_ip = 0.0           # IP_{z-1}
    max_log = -Inf

    # Temperature in cm⁻¹: T_cm = T_eV × 8065.54393734921
    temp_cm = temp_eV * 8065.54393734921

    for sc in charge_min:charge_max
        pf = get(pfs, sc, 0.0)
        pf <= 0 && continue
        ip = get(ips, sc, 0.0)

        # Saha decrement: log(n_{z+1}/n_z) = log(S_T · Z_{z+1}/Z_z) − IP_z / k_B T
        ratio = saha_factor / prev_pf
        log_abund = log(pf * ratio) + prev_log_abund - prev_ip / temp_cm

        saha[sc] = (log_abund, pf)
        max_log = max(max_log, log_abund)

        prev_log_abund = log_abund
        prev_pf = pf
        prev_ip = ip

        # Special case: when we reach the stage equal to the nuclear charge Z,
        # insert a "bare nucleus" stage (fully stripped ion with Z = 1).
        z = element_number(element)
        if z > 0 && z == sc
            sc_next = z + 1
            # For the bare nucleus, the Saha decrement uses Z_next = 1
            log_abund = log(saha_factor / prev_pf) + prev_log_abund - prev_ip / temp_cm
            saha[sc_next] = (log_abund, 1.0)
            max_log = max(max_log, log_abund)
            prev_log_abund = log_abund
            prev_pf = 1.0
            prev_ip = 0.0
        end
    end

    # Normalise: compute Σ exp(log_abund − max_log), then divide each by total
    total = sum(exp(v[1] - max_log) for v in values(saha))
    log_total = log(total)
    result = Dict{Int,Tuple{Float64,Float64}}()
    for (k, v) in saha
        # exp(v[1] - max_log - log_total) = exp(v[1]) / total  [numerically stable]
        result[k] = (exp(v[1] - max_log - log_total), v[2])
    end
    result
end
