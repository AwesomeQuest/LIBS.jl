using LIBS
using Plots
using Unitful

sticks = lte_spectrum_sticks("Fe I", 1.0u"eV", 1e17u"cm^-3";
    low_wl=370.0u"nm", upp_wl=380.0u"nm", min_rel_int=0.05)

so = lte_spectrum_data("Fe I", 1.0u"eV", 1e17u"cm^-3", 2000;
    low_wl=370.0u"nm", upp_wl=380.0u"nm", min_rel_int=0.05)

# Sticks alone — discrete delta-like intensities
p1 = plot(sticks; title="Fe I sticks (discrete transitions)", size=(600,300))

# Broadened spectrum alone — continuous line
p2 = plot(so.spectrum; title="Fe I Doppler spectrum (R=2000)", size=(600,300))

# Overlay — both on the same axes
p3 = plot(so; title="Fe I overlay (sticks + spectrum)", size=(600,300))

# Combined figure
plot(p1, p2, p3; layout=(3,1), size=(700,800))
savefig("plotdemo.png")
println("Saved plotdemo.png")
