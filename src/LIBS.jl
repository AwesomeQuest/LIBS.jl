"""Saha-LTE LIBS spectrum calculations with Unitful.jl quantities and Arrow.jl data tables."""
module LIBS

using Arrow
using Tables
using Unitful

export open_db, partition_function, ionization_potentials, saha_ion_populations
export lte_line_intensity, lte_spectrum_sticks, lte_spectrum_data, doppler_spectrum
export LIBSStickLine, DopplerGridPoint, SpectrumOverlay

include("elements.jl")
include("unitfull.jl")
include("db.jl")
include("convert.jl")
include("partition.jl")
include("saha.jl")
include("intensity.jl")
include("doppler.jl")
include("libs.jl")

end
