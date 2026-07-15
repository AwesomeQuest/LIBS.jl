const ELEMENTS = [
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
]

const ELEMENT_TO_Z = Dict{String,Int}()
for i in 1:length(ELEMENTS)
    ELEMENT_TO_Z[lowercase(ELEMENTS[i])] = i
end

element_symbol(z::Integer) = (1 <= z <= length(ELEMENTS)) ? ELEMENTS[z] : ""

function element_number(elem::AbstractString)
    s = replace(elem, r"^\d+" => "")
    get(ELEMENT_TO_Z, lowercase(s), 0)
end

function roman_to_int(roman::AbstractString)::Int
    r = uppercase(roman)
    d = Dict('I'=>1,'V'=>5,'X'=>10,'L'=>50,'C'=>100,'D'=>500,'M'=>1000)
    result = 0
    prev = 0
    for i in length(r):-1:1
        val = get(d, r[i], 0)
        val == 0 && return 0
        result += val < prev ? -val : val
        prev = val
    end
    result
end

struct SpectrumEntry
    Z::Int
    isotope::Union{Int,Nothing}
    charges::Vector{Int}
end

function parse_spectra(input::AbstractString)::Vector{SpectrumEntry}
    s = strip(input)
    isempty(s) && return SpectrumEntry[]
    lowercase(s) == "all spectra" && return SpectrumEntry[]

    entries = SpectrumEntry[]

    for group in split(s, r"\s*,\s*")
        group = strip(group)
        isempty(group) && continue

        m = match(r"^(\d+|\*)?([A-Za-z]+)\s*(.*)$", group)
        m === nothing && continue
        iso_str = m.captures[1]
        elem_part = m.captures[2]
        rest = strip(m.captures[3])

        has_isotope = iso_str !== nothing
        isotope = has_isotope ? (iso_str == "*" ? nothing : parse(Int, iso_str)) : nothing
        z = get(ELEMENT_TO_Z, lowercase(elem_part), 0)
        z == 0 && continue

        if isempty(rest)
            push!(entries, SpectrumEntry(z, isotope, Int[]))
            continue
        end

        rest = strip(rest)
        if occursin(r"^\d+(-\d*)?$", rest)
            m2 = match(r"^(\d+)(?:-(\d*))?$", rest)
            a = parse(Int, m2[1]) + 1
            b = (m2[2] !== nothing && !isempty(m2[2])) ? parse(Int, m2[2]) + 1 : a
            push!(entries, SpectrumEntry(z, isotope, collect(a:b)))
        elseif all(c -> c in "IVLCivlc", rest)
            stage = roman_to_int(rest)
            stage > 0 && push!(entries, SpectrumEntry(z, isotope, [stage]))
        elseif all(c -> c == '-', rest)
            push!(entries, SpectrumEntry(z, isotope, [-length(rest)]))
        elseif occursin(r"^[IVLCivlc]+(-[IVLCivlc]+)?$", rest)
            parts = split(rest, r"\s*-\s*")
            if length(parts) == 2
                a = roman_to_int(parts[1])
                b = roman_to_int(parts[2])
                if a > 0 && b >= a
                    push!(entries, SpectrumEntry(z, isotope, collect(a:b)))
                end
            end
        else
            push!(entries, SpectrumEntry(z, isotope, Int[]))
        end
    end

    return entries
end


