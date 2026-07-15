struct LIBSStickLine
    element::String
    spectr_charge::Int
    wavelength::Float64
    intensity::Float64
    low_conf::String
    low_term::String
    upp_conf::String
    upp_term::String
end

# Order by wavelength (used by sort! in lte_spectrum_sticks)
Base.isless(a::LIBSStickLine, b::LIBSStickLine) = a.wavelength < b.wavelength

# Compute LTE line intensity for a single transition.
#
# Formula:
#   I ∝ abundance × (N_sc / U_sc) × gA × exp(−E_up / k_B T)
#   where gA = A_ki × g_k (if A available) or estimated from line strength.
#
# int_scale=0: raw intensity
# int_scale=1: multiplied by transition energy ΔE (closer to photon flux)
function lte_line_intensity(row, saha::Dict{Int,Tuple{Float64,Float64}}, temp_eV::Real; int_scale::Integer=0, abundance::Real=1.0)
    sc = row.spectr_charge
    entry = get(saha, sc, nothing)
    entry === nothing && return 0.0
    ion_abund, pf = entry
    pf <= 0 && return 0.0

    upp_e = _v(row.upp_energy, 0.0)
    low_e = _v(row.low_energy, 0.0)
    upp_e <= 0 && return 0.0

    # gA = A × g_upp (Einstein coefficient × upper statistical weight)
    A = _v(row.A, 0.0)
    upp_g = _v(row.upp_g, 0)
    gA = A > 0 && upp_g > 0 ? A * upp_g : 0.0
    if gA <= 0
        # Fallback: estimate gA from line strength
        ls = _v(row.line_str, 0.0)
        if ls > 0
            wl = _v(row.vac_wl_num, 0.0)
            wl = wl > 0 ? wl : _v(row.calc_wl_num, 0.0)
            # gA = S × (2π m_e c / e² × 10⁸) / λ³ / g_upp ≈ S × 2.026e18 / λ³ / g_upp
            gA = ls * 2.0261269e18 / wl^3 / upp_g
        end
    end
    gA <= 0 && return 0.0

    temp_cm = temp_eV * CM_EV
    boltz = exp(-upp_e / temp_cm)
    factor = int_scale == 1 ? (upp_e - low_e) : 1.0
    factor * gA * boltz * ion_abund / pf * abundance
end
