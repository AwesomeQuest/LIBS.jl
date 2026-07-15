using Unitful: @u_str, ustrip, uconvert, Quantity

# Internal base units for computation
const INTERNAL_LENGTH = u"Å"
const INTERNAL_TEMP = u"eV"
const INTERNAL_DENSITY = u"cm^-3"
const INTERNAL_ENERGY = u"cm^-1"

# Convert user-facing temperature to internal eV.
# Accepts Quantity in eV or K, or bare Float64 (assumed eV).
# K → eV: k_B × T with k_B = 8.617333262e-5 eV/K.
const THERMAL_DIM = dimension(1.0u"K")

function to_internal_temp(t::Quantity)
    if dimension(t) == THERMAL_DIM
        return Float64(ustrip(t) * 8.617333262e-5)
    end
    uconvert(INTERNAL_TEMP, t) |> ustrip
end
to_internal_temp(t::Real) = Float64(t)

# Convert user-facing density to internal cm⁻³.
to_internal_density(d::Quantity) = uconvert(INTERNAL_DENSITY, d) |> ustrip
to_internal_density(d::Real) = Float64(d)

# Convert user-facing length to internal Å.
to_internal_length(wl::Quantity) = uconvert(INTERNAL_LENGTH, wl) |> ustrip
to_internal_length(wl::Real) = Float64(wl)

# Infer output length unit from input bounds.
# If neither bound given, default to nm.
output_unit(::Nothing, ::Nothing) = u"nm"
output_unit(low::Quantity, ::Nothing) = unit(low)
output_unit(::Nothing, upp::Quantity) = unit(upp)
output_unit(low::Quantity, ::Quantity) = unit(low)

# Wrap an internal Å value into the given unit.
wrap_output(wl::Real, u::Unitful.Units) = uconvert(u, Quantity(Float64(wl), INTERNAL_LENGTH))
wrap_output(wl::Real, ::Nothing) = Float64(wl)
