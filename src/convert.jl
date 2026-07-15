"""
    unwrap_or(x, default)

Return `x` if it is not `nothing` and not `missing`; otherwise return `default`.
Used throughout the codebase to safely extract nullable database columns.
"""
unwrap_or(x, default) = x === nothing || ismissing(x) ? default : x

"""
    nair1(sigma)

Refractive index of air at standard conditions (15 °C, 760 mmHg) from
the Edlén (1966) formula, where `sigma` = 1e8 / λ_vac with λ_vac in Å.

Returns the quantity (n - 1) × 10⁸, the refractivity in units of 10⁻⁸.

Reference: Edlén, B. (1966). The refractive index of air. *Metrologia*, 2(2), 71.
The formula is:
    (n - 1) × 10⁸ = 8060.51 + 2480990 / (132.274 × 10⁸ - σ²)
                        + 17455.7 / (39.32957 × 10⁸ - σ²)
where σ = 10⁸ / λ_vac, with λ_vac in Ångströms.
"""
function nair1(sigma::Real)
    # σ² = (10⁸ / λ_vac)² — wavenumber squared in (cm⁻¹)²
    σ² = float(sigma) * float(sigma)
    8060.51e-8 + 2480990 / (132.274e8 - σ²) + 17455.7 / (39.32957e8 - σ²)
end

"""
    vac_to_air(λ_vac)

Convert a vacuum wavelength `λ_vac` (in Å) to the corresponding
air wavelength (in Å) using the Edlén dispersion formula.

The conversion is:
    λ_air = λ_vac / (1 + nair1(10⁸ / λ_vac))

where nair1(sigma) is the refractivity of air at wavenumber sigma = 10⁸/λ.
"""
function vac_to_air(λ_vac::Real)
    float(λ_vac) / (1.0 + nair1(1e8 / float(λ_vac)))
end

"""
    air_to_vac(λ_air)

Convert an air wavelength `λ_air` (in Å) back to the corresponding
vacuum wavelength (in Å) by iteratively inverting the Edlén formula.

Uses a fixed-point iteration that typically converges within 2–3 iterations.
Returns `λ_air` unchanged if it is ≤ 1500 Å (below which air dispersion is negligible).
"""
function air_to_vac(λ_air::Real)
    λ = float(λ_air)
    # For λ ≤ 1500 Å, the air/vacuum shift is < 10⁻⁴ Å; skip conversion.
    λ <= 1500.0 && return λ
    tol = 1e-12
    # Fixed-point iteration: Lv = λ_air * (1 + nair1(10⁸ / Lv))
    Lv = λ
    for _ in 1:10
        Lv_prev = Lv
        Lv = λ * (1.0 + nair1(1e8 / Lv_prev))
        abs(Lv - Lv_prev) <= tol && break
    end
    Lv
end
