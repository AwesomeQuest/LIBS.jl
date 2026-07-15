struct DopplerGridPoint
    wavelength::Float64
    intensity::Float64
end

function doppler_spectrum(sticks::Vector{LIBSStickLine}, resolution::Real; grid::Union{Vector{Float64},Nothing}=nothing, pad::Real=6)
    isempty(sticks) && return DopplerGridPoint[]
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    wl_min = minimum(s.wavelength for s in sticks)
    wl_max = maximum(s.wavelength for s in sticks)
    wl_min <= wl_max || return DopplerGridPoint[]

    if grid === nothing
        sigma_min = wl_min / resolution
        sigma_max = wl_max / resolution
        sigma_edge = max(sigma_min, sigma_max)
        wl_lo = wl_min - pad * sigma_edge
        wl_hi = wl_max + pad * sigma_edge
        dl = sigma_min / 4
        npoints = max(2, Int(ceil((wl_hi - wl_lo) / dl)))
        grid = range(wl_lo, wl_hi; length=npoints)
    end

    sqpi = sqrt(π)
    intensities = zeros(Float64, length(grid))
    for s in sticks
        width = s.wavelength / resolution
        width2 = width * width
        for (i, wl) in enumerate(grid)
            if abs(wl - s.wavelength) <= 6 * width
                intensities[i] += s.intensity / width / sqpi * exp(-(s.wavelength - wl)^2 / width2)
            end
        end
    end

    [DopplerGridPoint(wl, intensities[i]) for (i, wl) in enumerate(grid)]
end
