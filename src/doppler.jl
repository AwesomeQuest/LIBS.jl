struct DopplerGridPoint
    wavelength::Float64
    intensity::Float64
end

function doppler_spectrum(sticks::Vector{LIBSStickLine}, resolution::Real; grid::Union{Vector{Float64},Nothing}=nothing)
    isempty(sticks) && return DopplerGridPoint[]
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    wl_min = minimum(s.wavelength for s in sticks)
    wl_max = maximum(s.wavelength for s in sticks)
    wl_min <= wl_max || return DopplerGridPoint[]

    if grid === nothing
        fwhm = wl_min / resolution
        dl = fwhm / 6
        npoints = max(2, Int(ceil((wl_max - wl_min) / dl)))
        grid = range(wl_min, wl_max; length=npoints)
    end

    sqpi = sqrt(π)
    fwhm_to_sigma = 1 / (2 * sqrt(log(2)))
    intensities = zeros(Float64, length(grid))
    for s in sticks
        sigma = s.wavelength / resolution * fwhm_to_sigma
        sigma2 = sigma * sigma
        for (i, wl) in enumerate(grid)
            if abs(wl - s.wavelength) <= 6 * sigma
                intensities[i] += s.intensity / sigma / sqpi * exp(-(s.wavelength - wl)^2 / sigma2)
            end
        end
    end

    [DopplerGridPoint(wl, intensities[i]) for (i, wl) in enumerate(grid)]
end
