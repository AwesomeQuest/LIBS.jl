"""
Saha-LTE LIBS spectrum calculations with Unitful.jl quantities and Arrow.jl data tables.

# Overview
This package computes LTE (Local Thermodynamic Equilibrium) emission spectra
for Laser-Induced Breakdown Spectroscopy (LIBS) using atomic line data from
the NIST Atomic Spectra Database.

The primary workflow is:

    db = LIBS.open_db()
    sticks = lte_spectrum_sticks("Fe I", 1.0u"eV", 1e17u"cm^-3";
                                 low_wl=200.0u"nm", upp_wl=600.0u"nm")
    so = lte_spectrum_data("Fe I", 1.0u"eV", 1e17u"cm^-3", 2000;
                           low_wl=200.0u"nm", upp_wl=600.0u"nm")

All temperature and density inputs accept `Unitful.jl` quantities
(e.g. `1.0u"eV"`, `12000.0u"K"`, `1e17u"cm^-3"`).  Bare Float64 values
are assumed to be in internal units (eV, cm⁻³).

# Physics
1. **Saha–Boltzmann populations** (`saha_ion_populations`):
   Compute the fractional abundance of each ionization stage in LTE
   given temperature and electron density using the Saha equation.

2. **Partition functions** (`partition_function`):
   Z(T) = Σ g_i exp(−E_i / k_B T), summed over all NIST energy levels.

3. **Line intensities** (`lte_line_intensity`, `lte_spectrum_sticks`):
   I_{ul} ∝ n_z/Z_z · g_u · A_{ul} · exp(−E_u / k_B T) · (E_u − E_l)ⁱ

4. **Doppler broadening** (`doppler_spectrum`, `lte_spectrum_data`):
   Convolve sticks with a Gaussian kernel of width σ = λ / R.

# Exported types
- `LIBSStickLine` — a single spectral line
- `DopplerGridPoint` — a point in a Doppler-broadened spectrum
- `SpectrumOverlay` — combined sticks + broadened spectrum for plotting

# Extensions
When `Plots.jl` is loaded, recipe methods are activated for plotting
`LIBSStickLine` vectors (scatter), `DopplerGridPoint` vectors (line),
and `SpectrumOverlay` (overlay).
"""
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
