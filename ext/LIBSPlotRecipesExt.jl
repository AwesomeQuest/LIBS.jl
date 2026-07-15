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
    getproperty.(sticks, :wavelength), getproperty.(sticks, :intensity)
end

@recipe function f(spectrum::AbstractVector{DopplerGridPoint})
    seriestype := :path
    xguide --> "Wavelength (nm)"
    yguide --> "Intensity"
    label --> "spectrum"
    linecolor := :crimson
    linewidth := 2
    getproperty.(spectrum, :wavelength), getproperty.(spectrum, :intensity)
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
        getproperty.(so.sticks, :wavelength), getproperty.(so.sticks, :intensity)
    end

    @series begin
        seriestype := :path
        label --> "spectrum"
        linecolor := :crimson
        linewidth := 2
        getproperty.(so.spectrum, :wavelength), getproperty.(so.spectrum, :intensity)
    end
end

end
