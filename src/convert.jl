unwrap_or(x, default) = x === nothing || ismissing(x) ? default : x

function nair1(sigma::Real)
    σ² = float(sigma) * float(sigma)
    8060.51e-8 + 2480990 / (132.274e8 - σ²) + 17455.7 / (39.32957e8 - σ²)
end

function vac_to_air(λ_vac::Real)
    float(λ_vac) / (1.0 + nair1(1e8 / float(λ_vac)))
end

function air_to_vac(λ_air::Real)
    λ = float(λ_air)
    λ <= 1500.0 && return λ
    tol = 1e-12
    Lv = λ
    for _ in 1:10
        Lv_prev = Lv
        Lv = λ * (1.0 + nair1(1e8 / Lv_prev))
        abs(Lv - Lv_prev) <= tol && break
    end
    Lv
end
