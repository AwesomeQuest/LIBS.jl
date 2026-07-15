# ---------------------------------------------------------------------------
# Unitful integration: internal base units, type aliases, and conversion
# helpers.  All computation is done in CGS-based atomic units internally;
# user-facing quantities can be in any Unitful-compatible unit.
# ---------------------------------------------------------------------------

using Unitful: @u_str, ustrip, uconvert, Quantity

# -- Internal base units ----------------------------------------------------
# All wavelength/frequency data is stored in the NIST database in
# Ångströms.  Temperature is converted to electron-volts for Boltzmann
# factors.  Number densities are in cm⁻³ (standard LIBS/CGS convention).
# Energy levels in the database are in cm⁻¹.
const INTERNAL_LENGTH  = u"Å"      # 1 Å = 10⁻¹⁰ m
const INTERNAL_TEMP    = u"eV"     # 1 eV ≈ 11604.5 K
const INTERNAL_DENSITY = u"cm^-3"  # particles per cm³
const INTERNAL_ENERGY  = u"cm^-1"  # wavenumber energy units

# -- Temperature: kelvin ↔ eV conversion factor -----------------------------
# k_B = 8.617333262 × 10⁻⁵ eV/K  (CODATA 2018)
const _KB_EV_PER_K = 8.617333262e-5

# Dimensional anchor for detecting temperature-in-kelvin inputs.
const THERMAL_DIM = dimension(1.0u"K")

"""
    to_internal_temp(t)

Convert a temperature to internal units (eV).  Accepts:
- A `Quantity` in eV → returned as bare Float64.
- A `Quantity` in K → converted via k_B × T (Boltzmann constant).
- A bare `Real` → assumed already in eV.

Returns a bare `Float64` (unitless) for use in Boltzmann factors.
"""
function to_internal_temp(t::Quantity)
    # If the quantity has dimensions of temperature (K), convert K → eV
    if dimension(t) == THERMAL_DIM
        return Float64(ustrip(t) * _KB_EV_PER_K)
    end
    # Otherwise assume it's already in energy units (eV) and convert
    uconvert(INTERNAL_TEMP, t) |> ustrip
end
to_internal_temp(t::Real) = Float64(t)

"""
    to_internal_density(d)

Convert a number density to internal units (cm⁻³).
Accepts a `Quantity` (any compatible unit) or a bare `Real` (assumed cm⁻³).
"""
to_internal_density(d::Quantity) = uconvert(INTERNAL_DENSITY, d) |> ustrip
to_internal_density(d::Real) = Float64(d)

"""
    to_internal_length(wl)

Convert a length/wavelength to internal units (Å).
Accepts a `Quantity` (any compatible unit) or a bare `Real` (assumed Å).
Used internally for database queries and Doppler grid computation.
"""
to_internal_length(wl::Quantity) = uconvert(INTERNAL_LENGTH, wl) |> ustrip
to_internal_length(wl::Real) = Float64(wl)

"""
    output_unit(low_wl, upp_wl)

Determine the output length unit for stick/spectrum wavelengths.
- If both bounds are `nothing` → default to nm.
- If one or both bounds are `Quantity` → use the unit of the first one.
"""
output_unit(::Nothing, ::Nothing) = u"nm"
output_unit(low::Quantity, ::Nothing) = unit(low)
output_unit(::Nothing, upp::Quantity) = unit(upp)
output_unit(low::Quantity, ::Quantity) = unit(low)

"""
    wrap_output(wl, u)

Convert an internal Å value `wl` into a `Quantity` in the requested unit `u`.
If `u` is `nothing`, returns a bare Float64 (unitless, for backward compatibility).
"""
wrap_output(wl::Real, u::Unitful.Units) = uconvert(u, Quantity(Float64(wl), INTERNAL_LENGTH))
wrap_output(wl::Real, ::Nothing) = Float64(wl)
