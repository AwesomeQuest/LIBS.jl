# LIBS.jl

LTE emission spectrum synthesis for laser-induced breakdown spectroscopy,
using atomic line data from the NIST Atomic Spectra Database.

All physical inputs accept [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)
quantities (eV, K, cm⁻³, nm, Å) with automatic conversion to internal CGS-based
atomic units.

## Installation

```julia
julia> import Pkg; Pkg.add(url="https://github.com/AwesomeQuest/LIBS.jl")
```

The NIST line data is packaged as Julia artifacts and downloads automatically
on first use.

## Quick start

```julia
using LIBS
using Unitful

db = open_db()
```

**Stick spectrum** — discrete emission lines:

```julia
sticks = lte_spectrum_sticks("Fe I", 1.0u"eV", 1e17u"cm^-3";
    low_wl=200.0u"nm", upp_wl=600.0u"nm")
```

**Doppler-broadened spectrum** — sticks convolved with a Gaussian kernel:

```julia
so = lte_spectrum_data("Fe I", 1.0u"eV", 1e17u"cm^-3", 2000;
    low_wl=200.0u"nm", upp_wl=600.0u"nm")
```

**Multi-element mixtures:**

```julia
sticks = lte_spectrum_sticks("Fe I, Ni I, Cr I", 1.2u"eV", 5e16u"cm^-3";
    composition=Dict("Fe" => 0.7, "Ni" => 0.2, "Cr" => 0.1))
```

Temperature can also be given in Kelvin:

```julia
sticks = lte_spectrum_sticks("Fe I", 12000.0u"K", 1e17u"cm^-3")
```

## API

### Database

| Function | Description |
|---|---|
| `open_db()` | Return the module-level LIBSDB singleton (Arrow tables) |

### Plasma statics

| Function | Description |
|---|---|
| `partition_function(db, elem, charge, temp)` | Internal partition function Z(T) = Σ g·exp(−E/k_B T) |
| `ionization_potentials(db, elem)` | Dict mapping charge stage → IP (cm⁻¹) |
| `saha_ion_populations(db, elem, temp, eden; ...)` | Saha–Boltzmann fractional abundances for each ionisation stage |

### Spectrum synthesis

| Function | Description |
|---|---|
| `lte_spectrum_sticks(spectra, temp, density; ...)` | Discrete line stick spectrum |
| `doppler_spectrum(sticks, resolution; ...)` | Convolve sticks with Gaussian Doppler profile |
| `lte_spectrum_data(spectra, temp, density, resolution; ...)` | Sticks + Doppler in one call |

### Spectrum specification syntax

Examples passed to `lte_spectrum_sticks` / `lte_spectrum_data`:

| Input | Meaning |
|---|---|
| `"Fe I"` | Fe neutral (charge 1) |
| `"Fe I-III"` | Fe stages I, II, III (charge 1–3) |
| `"Fe0-2"` | Same, 0-indexed numeric form |
| `"Fe"` | All available charge states for Fe |
| `"Fe I, Ni I"` | Multi-element |
| `"56Fe I"` | Isotope-specific (⁵⁶Fe) |
| `"all spectra"` | No-op (return empty) |

## Physics

### Partition function

Z(T) = Σ_i g_i · exp(−E_i / k_B T)

summed over all NIST energy levels with valid energy and degeneracy.
Temperatures are converted to eV internally; Boltzmann factors use
the conversion 1 eV = 8065.54 cm⁻¹.

### Saha–Boltzmann ionisation

The Saha equation in LTE:

    n_{z+1} · n_e     Z_{z+1}(T) · (2π m_e k_B T)^{3/2}
    ─────────────  =  ────────────────────────────────── · exp(−χ_z / k_B T)
        n_z                    Z_z(T) · h³

In log form with the CGS prefactor S_T = 6.043 × 10²¹ · T_{eV}^{3/2} / n_e:

    log(n_{z+1}/n_z) = log(S_T · Z_{z+1}/Z_z) − IP_z / k_B T

Populations are normalised so Σ_z n_z = 1.

### Line intensity

I_{ul} ∝ (E_u − E_l)^i · g_u · A_{ul} · exp(−E_u / k_B T) · n_z / Z_z

where i = int_scale (0 or 1).  If the Einstein A-value is unavailable,
g_u · A is computed from the line strength S:

    g_u · A = S · 2.026 × 10¹⁸ / λ³    (λ in Å, A in s⁻¹)

### Doppler broadening

Each stick line at λ_0 is replaced by a unit-area Gaussian:

    ϕ(λ) = exp(−(λ − λ_0)² / σ²) / (σ√π)

with σ = λ_0 / R (R = user-specified resolving power).
The adaptive grid spans λ_min/max ± 6σ_max with spacing σ_min / 4.

## Plotting

When `Plots.jl` is loaded, recipe methods are activated:

```julia
using Plots

# Stick lines (scatter)
plot(sticks)

# Broadened spectrum (line)
plot(so.spectrum)

# Overlay (sticks + spectrum)
plot(so)
```

## Data

Atomic line data is sourced from the NIST Atomic Spectra Database (ASD)
and packaged as Arrow tables per element.  The build pipeline is in
`build/` and includes scripts to download and convert the NIST ASD
export files.

## Testing

```bash
julia --project -e 'import Pkg; Pkg.test()'
```

Tests validate stick wavelengths against NIST reference data across
14 elements, and verify Doppler spectrum area conservation, monotonicity,
and resolution scaling.

## References

- NIST Atomic Spectra Database: https://physics.nist.gov/asd
- Edlén, B. (1966). The refractive index of air. *Metrologia*, 2(2), 71.
- Saha, M. N. (1920). Ionisation in the solar chromosphere. *Phil. Mag.*, 40, 472.
- Griem, H. R. (1997). *Principles of Plasma Spectroscopy*. Cambridge.
