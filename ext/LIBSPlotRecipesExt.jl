module LIBSPlotRecipesExt

using LIBS
using Plots

@recipe function f(sticks::AbstractVector{LIBSStickLine})
    seriestype := :scatter
    xguide --> "Wavelength (nm)"
    yguide --> "Intensity"
    label --> "sticks"
    linecolor := :steelblue
    markercolor := :steelblue
    markerstrokecolor := :steelblue
    markershape := :circle
    markersize := 3
    x := [s.wavelength for s in sticks]
    y := [s.intensity for s in sticks]
end

@recipe function f(spectrum::AbstractVector{DopplerGridPoint})
    seriestype := :path
    xguide --> "Wavelength (nm)"
    yguide --> "Intensity"
    label --> "spectrum"
    linecolor := :crimson
    linewidth := 2
    x := [p.wavelength for p in spectrum]
    y := [p.intensity for p in spectrum]
end

@recipe function f(so::SpectrumOverlay)
    xguide --> "Wavelength (nm)"
    yguide --> "Intensity"

    @series begin
        seriestype := :scatter
        label --> "sticks"
        linecolor := :steelblue
        markercolor := :steelblue
        markerstrokecolor := :steelblue
        markershape := :circle
        markersize := 3
        x := [s.wavelength for s in so.sticks]
        y := [s.intensity for s in so.sticks]
    end

    @series begin
        seriestype := :path
        label --> "spectrum"
        linecolor := :crimson
        linewidth := 2
        x := [p.wavelength for p in so.spectrum]
        y := [p.intensity for p in so.spectrum]
    end
end

end
