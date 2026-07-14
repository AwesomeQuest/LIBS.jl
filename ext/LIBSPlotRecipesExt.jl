module LIBSPlotRecipesExt

using LIBS
using Plots

@recipe function f(::Type{<:Union{LIBSStickLine, AbstractVector{LIBSStickLine}}}, sticks::AbstractVector{LIBSStickLine})
    seriestype --> :stem
    xguide --> "Wavelength"
    yguide --> "Intensity"
    label --> "sticks"
    linecolor --> :steelblue
    markercolor --> :steelblue
    markerstrokecolor --> :steelblue
    markershape --> :circle
    markersize --> 3

    [s.wavelength for s in sticks], [s.intensity for s in sticks]
end

@recipe function f(::Type{<:Union{DopplerGridPoint, AbstractVector{DopplerGridPoint}}}, spectrum::AbstractVector{DopplerGridPoint})
    seriestype --> :path
    xguide --> "Wavelength"
    yguide --> "Intensity"
    label --> "spectrum"
    linecolor --> :crimson
    linewidth --> 2

    [p.wavelength for p in spectrum], [p.intensity for p in spectrum]
end

@recipe function f(::Type{<:Union{SpectrumOverlay, AbstractVector{SpectrumOverlay}}}, so::SpectrumOverlay)
    xguide --> "Wavelength"
    yguide --> "Intensity"

    @series begin
        seriestype --> :stem
        label --> "sticks"
        linecolor --> :steelblue
        markercolor --> :steelblue
        markerstrokecolor --> :steelblue
        markershape --> :circle
        markersize --> 3
        [s.wavelength for s in so.sticks], [s.intensity for s in so.sticks]
    end

    @series begin
        seriestype --> :path
        label --> "spectrum"
        linecolor --> :crimson
        linewidth --> 2
        [p.wavelength for p in so.spectrum], [p.intensity for p in so.spectrum]
    end
end

end
