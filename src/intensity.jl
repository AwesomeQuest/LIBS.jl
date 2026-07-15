struct LIBSStickLine
    element::String
    spectr_charge::Int
    wavelength::Quantity  # internal base: Å
    intensity::Float64
    low_conf::String
    low_term::String
    upp_conf::String
    upp_term::String
end

Base.isless(a::LIBSStickLine, b::LIBSStickLine) = a.wavelength < b.wavelength

function lte_line_intensity(row, saha::Dict{Int,Tuple{Float64,Float64}}, temp_eV::Real; int_scale::Integer=0, abundance::Real=1.0)
    sc = row.spectr_charge
    entry = get(saha, sc, nothing)
    entry === nothing && return 0.0
    ion_abund, pf = entry
    pf <= 0 && return 0.0

    upp_e = unwrap_or(row.upp_energy, 0.0)
    low_e = unwrap_or(row.low_energy, 0.0)
    upp_e <= 0 && return 0.0

    A = unwrap_or(row.A, 0.0)
    upp_g = unwrap_or(row.upp_g, 0)
    gA = A > 0 && upp_g > 0 ? A * upp_g : 0.0
    if gA <= 0
        ls = unwrap_or(row.line_str, 0.0)
        if ls > 0
            wl = unwrap_or(row.vac_wl_num, 0.0)
            wl = wl > 0 ? wl : unwrap_or(row.calc_wl_num, 0.0)
            gA = ls * 2.0261269e18 / wl^3 / upp_g
        end
    end
    gA <= 0 && return 0.0

    temp_cm = temp_eV * 8065.54393734921
    boltz = exp(-upp_e / temp_cm)
    factor = int_scale == 1 ? (upp_e - low_e) : 1.0
    factor * gA * boltz * ion_abund / pf * abundance
end
