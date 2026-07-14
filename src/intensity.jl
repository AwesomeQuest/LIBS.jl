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

Base.isless(a::LIBSStickLine, b::LIBSStickLine) = a.wavelength < b.wavelength

function lte_line_intensity(row, saha::Dict{Int,Tuple{Float64,Float64}}, temp_eV::Real; int_scale::Int=0, abundance::Real=1.0)
    sc = row.spectr_charge
    entry = get(saha, sc, nothing)
    entry === nothing && return 0.0
    ion_abund, pf = entry
    pf <= 0 && return 0.0

    upp_e = _v(row.upp_energy, 0.0)
    low_e = _v(row.low_energy, 0.0)
    upp_e <= 0 && return 0.0

    A = _v(row.A, 0.0)
    upp_g = _v(row.upp_g, 0)
    gA = A > 0 && upp_g > 0 ? A * upp_g : 0.0
    if gA <= 0
        ls = _v(row.line_str, 0.0)
        if ls > 0
            wl = _v(row.vac_wl_num, 0.0)
            wl = wl > 0 ? wl : _v(row.calc_wl_num, 0.0)
            gA = ls * 2.0261269e18 / wl^3 / upp_g
        end
    end
    gA <= 0 && return 0.0

    temp_cm = temp_eV * CM_EV
    boltz = exp(-upp_e / temp_cm)
    factor = int_scale == 1 ? (upp_e - low_e) : 1.0
    return factor * gA * boltz * ion_abund / pf * abundance
end
