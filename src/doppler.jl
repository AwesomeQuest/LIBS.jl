# ---------------------------------------------------------------------------
# Doppler broadening.
#
# In a plasma in LTE, each atomic transition line is broadened by the
# thermal motion of the emitting atoms.  The line profile is Gaussian:
#
#     ϕ(λ) = 1 / (σ √{2π}) · exp(−(λ − λ₀)² / (2σ²))    [in wavelength]
#
# with the Doppler width σ (1σ Gaussian width):
#
#     σ = λ₀ · √(k_B T / m c²)   =   λ₀ / R
#
# where R = c / v_th is the resolving power corresponding to the
# thermal velocity v_th = √(k_B T / m).  In this code the user specifies
# a resolution R = λ₀ / σ directly, which captures the combined effect
# of temperature, atomic mass, and instrument broadening.
#
# The broadened spectrum is the convolution:
#
#     I(λ) = Σ_i I_i · ϕ(λ − λ_{0,i})
#
# where each stick line at λ_{0,i} with intensity I_i is replaced by
# a unit-area Gaussian of width σ = λ_{0,i} / R.
# ---------------------------------------------------------------------------

"""
    DopplerGridPoint

A single point in a Doppler-broadened spectrum.
Fields:
- `wavelength` — wavelength as a `Quantity` (same unit as the input sticks)
- `intensity`  — convolved spectral intensity at this wavelength
"""
struct DopplerGridPoint
    wavelength::Quantity  # output unit (default: nm, inherited from input sticks)
    intensity::Float64
end

"""
    doppler_spectrum(sticks, resolution; grid=nothing, pad=6)

Convolve a set of stick lines with a Doppler (Gaussian) profile to produce
a smooth spectrum.

# Arguments
- `sticks`     — vector of `LIBSStickLine` from `lte_spectrum_sticks`
- `resolution` — resolving power R = λ / Δλ_FWHM, where Δλ_FWHM is the
                 full-width at half-maximum of the Gaussian kernel.
                 The relationship is: σ = λ / R, FWHM = 2σ√(2 ln 2) ≈ 2.355 σ.
- `grid`       — optional explicit wavelength grid (Quantity vector).
                 If `nothing`, an adaptive grid is constructed.
- `pad`        — number of σ to extend the grid beyond the extreme stick
                 wavelengths (default: 6, which captures >99.999% of the area).

# Returns
`Vector{DopplerGridPoint}` representing the broadened spectrum on a
uniform wavelength grid (in the same unit as the input sticks).

# Algorithm
1. Convert all stick wavelengths to internal Å units.
2. Build an adaptive uniform grid spanning [λ_min − pad·σ_max, λ_max + pad·σ_max]
   with spacing σ_min / 4 (where σ_min = λ_min / R, σ_max = λ_max / R).
3. For each stick, add Gaussian contributions to all grid points within
   pad·σ of the line centre:
       I_g += I_i / (σ√π) · exp(−(λ_i − λ_g)² / σ²)
   Each kernel has unit area: ∫ I_i / (σ√π) · exp(−(Δλ/σ)²) dλ = I_i.
4. Wrap output wavelengths back to the user's requested unit.
"""
function doppler_spectrum(sticks::AbstractVector{LIBSStickLine}, resolution::Real; grid=nothing, pad::Real=6)
    isempty(sticks) && return DopplerGridPoint[]
    resolution > 0 || throw(ArgumentError("resolution must be positive"))

    # Work in internal Å for consistent numerics regardless of output unit.
    # Convert stick wavelengths to Å explicitly.
    wls = [ustrip(INTERNAL_LENGTH, s.wavelength) for s in sticks]
    ints = [s.intensity for s in sticks]
    out_unit = unit(first(sticks).wavelength)

    wl_min = minimum(wls)
    wl_max = maximum(wls)
    wl_min <= wl_max || return DopplerGridPoint[]

    if grid === nothing
        # Gaussian 1σ widths at the extremes (in Å)
        sigma_min = wl_min / resolution
        sigma_max = wl_max / resolution
        # Use the larger σ for padding to ensure full coverage
        sigma_edge = max(sigma_min, sigma_max)
        # Extend the grid by pad·σ on each side
        wl_lo = wl_min - pad * sigma_edge
        wl_hi = wl_max + pad * sigma_edge
        # Grid spacing: resolve the narrowest line with ~4 points per σ
        dl = sigma_min / 4
        npoints = max(2, Int(ceil((wl_hi - wl_lo) / dl)))
        grid_A = range(wl_lo, wl_hi; length=npoints)
    else
        grid_A = ustrip.(uconvert.(INTERNAL_LENGTH, grid))
    end

    # Pre-factor: the Gaussian kernel is exp(−(Δλ/σ)²) / (σ√π)
    # ∫ (1 / (σ√π)) · exp(−(Δλ/σ)²) dλ = 1
    sqpi = sqrt(π)
    intensities = zeros(Float64, length(grid_A))

    # Convolve: sum Gaussian contributions from each stick
    for idx in eachindex(wls)
        wl = wls[idx]
        intensity = ints[idx]
        width = wl / resolution  # σ for this line
        width2 = width * width
        # Only compute contributions within pad·σ of the line centre
        for (i, gw) in enumerate(grid_A)
            if abs(gw - wl) <= pad * width
                intensities[i] += intensity / width / sqpi * exp(-(wl - gw)^2 / width2)
            end
        end
    end

    # Convert the grid wavelengths back to the output unit
    [DopplerGridPoint(wrap_output(gw, out_unit), intensities[i]) for (i, gw) in enumerate(grid_A)]
end
