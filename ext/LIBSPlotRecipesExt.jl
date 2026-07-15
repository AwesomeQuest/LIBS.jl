module LIBSPlotRecipesExt

using LIBS
using Plots
using Unitful: unit, ustrip

@recipe function f(sticks::AbstractVector{LIBSStickLine})
    seriestype := :scatter
    label --> "sticks"
    linecolor := :steelblue
    markercolor := :steelblue
    markerstrokecolor := :steelblue
    markershape := :circle
    markersize := 3
    ws = ustrip.(unit(first(sticks).wavelength), getproperty.(sticks, :wavelength))
    xguide --> "Wavelength ($(unit(first(sticks).wavelength)))"
    yguide --> "Intensity"
    ws, getproperty.(sticks, :intensity)
end

@recipe function f(spectrum::AbstractVector{DopplerGridPoint})
    seriestype := :path
    label --> "spectrum"
    linecolor := :crimson
    linewidth := 2
    ws = ustrip.(unit(first(spectrum).wavelength), getproperty.(spectrum, :wavelength))
    xguide --> "Wavelength ($(unit(first(spectrum).wavelength)))"
    yguide --> "Intensity"
    ws, getproperty.(spectrum, :intensity)
end

@recipe function f(so::SpectrumOverlay)
    u = unit(first(so.sticks).wavelength)
    xguide --> "Wavelength ($u)"
    yguide --> "Intensity"
    stick_wls = ustrip.(u, getproperty.(so.sticks, :wavelength))
    xlims --> (minimum(stick_wls), maximum(stick_wls))

    @series begin
        seriestype := :scatter
        label --> "sticks"
        linecolor := :steelblue
        markercolor := :steelblue
        markerstrokecolor := :steelblue
        markershape := :circle
        markersize := 3
        stick_wls, getproperty.(so.sticks, :intensity)
    end

    @series begin
        seriestype := :path
        label --> "spectrum"
        linecolor := :crimson
        linewidth := 2
        ustrip.(u, getproperty.(so.spectrum, :wavelength)), getproperty.(so.spectrum, :intensity)
    end
end

end
