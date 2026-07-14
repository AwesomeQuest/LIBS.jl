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

@userplot SpectrumPlot

@recipe function f(sp::SpectrumPlot)
    args = sp.args
    if length(args) == 1 && args[1] isa Tuple
        sticks, spectrum = args[1]
    elseif length(args) == 2
        sticks, spectrum = args
        sticks isa AbstractVector{LIBSStickLine} || error("First argument must be Vector{LIBSStickLine}")
        spectrum isa AbstractVector{DopplerGridPoint} || error("Second argument must be Vector{DopplerGridPoint}")
    else
        error("SpectrumPlot requires (sticks, spectrum) or a tuple thereof")
    end

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
        [s.wavelength for s in sticks], [s.intensity for s in sticks]
    end

    @series begin
        seriestype --> :path
        label --> "spectrum"
        linecolor --> :crimson
        linewidth --> 2
        [p.wavelength for p in spectrum], [p.intensity for p in spectrum]
    end
end

end
