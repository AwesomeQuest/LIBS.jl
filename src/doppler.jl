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
        min_width = wl_min / resolution
        dl = min_width / 4
        npoints = max(1, Int(ceil((wl_max - wl_min) / dl)))
        grid = range(wl_min, wl_max; length=npoints)
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
