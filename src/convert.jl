# k_B in eV → cm⁻¹: 1 eV = 8065.54393734921 cm⁻¹
const CM_EV = 8065.54393734921

# Coalesce — return x if not nothing/missing, else default
_v(x, default) = x === nothing || ismissing(x) ? default : x

# Dispersion formula for refractive index of air (Peck & Reeder, 1972).
# σ = 1/λ in μm⁻¹; nair1 = (n − 1) × 10⁸
function nair1(sigma::Real)
    σ² = float(sigma) * float(sigma)
    8060.51e-8 + 2480990 / (132.274e8 - σ²) + 17455.7 / (39.32957e8 - σ²)
end

# Convert vacuum wavelength (Å) to air wavelength (Å)
function vac_to_air(λ_vac::Real)
    L = float(λ_vac) / (1.0 + nair1(1e8 / float(λ_vac)))
    return L
end

# Convert air wavelength (Å) to vacuum wavelength (Å).
# Uses fixed-point iteration for λ > 1500 Å; below that the correction is negligible.
function air_to_vac(λ_air::Real)
    λ = float(λ_air)
    λ <= 1500.0 && return λ
    tol = 1e-12
    Lv = λ
    for _ in 1:10
        Lv_prev = Lv
        # Refractive index depends on vacuum wavelength — iterate to converge
        Lv = λ * (1.0 + nair1(1e8 / Lv_prev))
        abs(Lv - Lv_prev) <= tol && break
    end
    Lv
end
