const CM_EV = 8065.54393734921

_v(x, default) = x === nothing || ismissing(x) ? default : x

function nair1(sigma::Float64)
    σ² = sigma * sigma
    n = 8060.51e-8 + 2480990 / (132.274e8 - σ²) + 17455.7 / (39.32957e8 - σ²)
    return n
end

function vac_to_air(λ_vac::Float64)
    L = λ_vac / (1.0 + nair1(1e8 / λ_vac))
    return L
end

function air_to_vac(λ_air::Float64)
    tol = 1e-18
    if λ_air > 1500.0
        Lv1 = λ_air
        for _ in 1:10
            Lv0 = 1.000281383 * λ_air
            Lv1 = λ_air * (1.0 + nair1(1e8 / Lv0))
            if abs(Lv1 - Lv0) <= tol
                break
            end
        end
        return Lv1
    else
        return λ_air
    end
end
