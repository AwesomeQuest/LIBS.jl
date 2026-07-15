# ---------------------------------------------------------------------------
# Line intensity computation and the LIBSStickLine data structure.
#
# In LTE, the spectral line intensity for a transition from upper level u to
# lower level l is proportional to:
#
#     I_{ul} ∝ (E_u − E_l) · g_u · A_{ul} · exp(−E_u / k_B T) · n_z(T) / Z_z(T)
#
# where
#   A_{ul}  = Einstein A coefficient (s⁻¹)
#   g_u     = statistical weight of the upper level
#   E_u     = energy of the upper level (cm⁻¹)
#   n_z(T)  = relative abundance of ionisation stage z (from Saha)
#   Z_z(T)  = internal partition function of stage z
# ---------------------------------------------------------------------------

"""
    LIBSStickLine

A single spectral line produced by `lte_spectrum_sticks`.
Fields:
- `element`       — element symbol (e.g. "Fe")
- `spectr_charge` — ionisation stage (1 = neutral, 2 = singly ionised, …)
- `wavelength`    — transition wavelength as a `Quantity` (default output: nm)
- `intensity`     — relative line intensity (arbitrary units)
- `low_conf`      — lower level configuration string
- `low_term`      — lower level term string
- `upp_conf`      — upper level configuration string
- `upp_term`      — upper level term string

Sticks are sorted by wavelength via `Base.isless`.
"""
struct LIBSStickLine
    element::String
    spectr_charge::Int
    wavelength::Quantity  # output unit (nm by default); internally stored via wrap_output
    intensity::Float64
    low_conf::String
    low_term::String
    upp_conf::String
    upp_term::String
end

Base.isless(a::LIBSStickLine, b::LIBSStickLine) = a.wavelength < b.wavelength

"""
    lte_line_intensity(row, saha, temp_eV; int_scale=0, abundance=1.0)

Compute the LTE line intensity for a single transition.

# Arguments
- `row`       — a database row (Arrow table entry) with fields:
                spectr_charge, A, upp_g, upp_energy, low_energy,
                line_str, vac_wl_num, calc_wl_num
- `saha`      — output of `saha_ion_populations`: a Dict mapping
                stage → (fractional_abundance, partition_function)
- `temp_eV`   — temperature in eV (bare Float64)
- `int_scale` — if 1, multiply by (E_u − E_l) to get photon-energy-weighted
                intensity (more realistic for detector response)
- `abundance` — elemental abundance fraction (0…1)

# Returns
A Float64 intensity value (arbitrary relative scale).  Returns 0.0 if
the transition has no valid A-value or line strength.

# Intensity formula
    I = abundance · n_z / Z_z · g_u · A_{ul} · exp(−E_u / k_B T) · f(int_scale)

where f(int_scale) = (E_u − E_l) if int_scale == 1, else 1.
If A · g_u is not available, an attempt is made to compute g_u · A from
the line strength S (in atomic units) via:
    g_u · A = S · 2.0261269 × 10¹⁸ / λ³     [λ in Å, A in s⁻¹]
"""
function lte_line_intensity(row, saha::Dict{Int,Tuple{Float64,Float64}}, temp_eV::Real; int_scale::Integer=0, abundance::Real=1.0)
    sc = row.spectr_charge
    entry = get(saha, sc, nothing)
    entry === nothing && return 0.0
    ion_abund, pf = entry
    pf <= 0 && return 0.0

    # Upper level energy (cm⁻¹)
    upp_e = unwrap_or(row.upp_energy, 0.0)
    low_e = unwrap_or(row.low_energy, 0.0)
    upp_e <= 0 && return 0.0

    # Try to get gA = g_u · A_{ul} directly from the database
    A = unwrap_or(row.A, 0.0)
    upp_g = unwrap_or(row.upp_g, 0)
    gA = A > 0 && upp_g > 0 ? A * upp_g : 0.0

    if gA <= 0
        # Fallback: compute gA from line strength S (in a.u.)
        #   g_u · A_{ul} = S · 2π² e² / (ε₀ m_e c λ³)
        #                 = S · 2.0261269 × 10¹⁸ / λ³    [λ in Å]
        ls = unwrap_or(row.line_str, 0.0)
        if ls > 0
            wl = unwrap_or(row.vac_wl_num, 0.0)
            wl = wl > 0 ? wl : unwrap_or(row.calc_wl_num, 0.0)
            gA = ls * 2.0261269e18 / wl^3 / upp_g
        end
    end
    gA <= 0 && return 0.0

    # Convert eV → cm⁻¹ for Boltzmann factor: 1 eV = 8065.54393734921 cm⁻¹
    temp_cm = temp_eV * 8065.54393734921

    # Boltzmann factor: exp(−E_u / k_B T)
    boltz = exp(-upp_e / temp_cm)

    # Optional photon-energy weighting (E_u − E_l) for response scaling
    factor = int_scale == 1 ? (upp_e - low_e) : 1.0

    # I = abundance · n_z / Z_z · g_u · A_{ul} · exp(−E_u / k_B T) · factor
    factor * gA * boltz * ion_abund / pf * abundance
end
