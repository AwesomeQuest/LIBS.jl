struct DopplerGridPoint
    wavelength::Quantity  # internal base: Å
    intensity::Float64
end

function doppler_spectrum(sticks::AbstractVector{LIBSStickLine}, resolution::Real; grid=nothing, pad::Real=6)
    isempty(sticks) && return DopplerGridPoint[]
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    # Work in raw Å for performance
    wls = [ustrip(s.wavelength) for s in sticks]
    ints = [s.intensity for s in sticks]
    out_unit = unit(first(sticks).wavelength)

    wl_min = minimum(wls)
    wl_max = maximum(wls)
    wl_min <= wl_max || return DopplerGridPoint[]

    if grid === nothing
        sigma_min = wl_min / resolution
        sigma_max = wl_max / resolution
        sigma_edge = max(sigma_min, sigma_max)
        wl_lo = wl_min - pad * sigma_edge
        wl_hi = wl_max + pad * sigma_edge
        dl = sigma_min / 4
        npoints = max(2, Int(ceil((wl_hi - wl_lo) / dl)))
        grid_A = range(wl_lo, wl_hi; length=npoints)
    else
        grid_A = ustrip.(uconvert.(INTERNAL_LENGTH, grid))
    end

    sqpi = sqrt(π)
    intensities = zeros(Float64, length(grid_A))
    for idx in eachindex(wls)
        wl = wls[idx]
        intensity = ints[idx]
        width = wl / resolution
        width2 = width * width
        for (i, gw) in enumerate(grid_A)
            if abs(gw - wl) <= pad * width
                intensities[i] += intensity / width / sqpi * exp(-(wl - gw)^2 / width2)
            end
        end
    end

    [DopplerGridPoint(wrap_output(gw, out_unit), intensities[i]) for (i, gw) in enumerate(grid_A)]
end
