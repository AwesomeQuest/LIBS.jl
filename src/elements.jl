# ---------------------------------------------------------------------------
# Periodic table data, element lookups, Roman numeral conversion, and
# spectrum specification parsing.
# ---------------------------------------------------------------------------

# Elements Z = 1 … 118 (IUPAC names, sorted by atomic number).
const ELEMENTS = [
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
]

# Reverse lookup: lowercase symbol → atomic number Z.
const ELEMENT_TO_Z = Dict{String,Int}()
for (i, sym) in enumerate(ELEMENTS)
    ELEMENT_TO_Z[lowercase(sym)] = i
end

"""
    element_symbol(z)

Return the standard IUPAC symbol for atomic number `z` (1 … 118).
"""
element_symbol(z::Integer) = ELEMENTS[z]

"""
    element_number(elem)

Return the atomic number for a given element string (case-insensitive,
ignoring leading digits used for isotope notation, e.g. "56Fe" → 26).
Returns 0 if the element is not found.
"""
function element_number(elem::AbstractString)
    # Strip leading isotope number, e.g. "56Fe" → "Fe"
    s = replace(elem, r"^\d+" => "")
    get(ELEMENT_TO_Z, lowercase(s), 0)
end

# Roman numeral → integer mapping (I=1, V=5, X=10, L=50, C=100, D=500, M=1000).
const ROMAN_VALUES = Dict('I'=>1, 'V'=>5, 'X'=>10, 'L'=>50, 'C'=>100, 'D'=>500, 'M'=>1000)

"""
    roman_to_int(roman)

Convert a Roman numeral (e.g. "I", "II", "III", "IV") to an integer.
Used to parse NIST ionization stages (I = neutral, II = singly ionized, ...).
Standard subtractive notation is supported (IV → 4, IX → 9, etc.).
"""
function roman_to_int(roman::AbstractString)
    r = uppercase(roman)
    result = 0
    prev = 0
    # Traverse right-to-left: subtract smaller-before-larger, add otherwise.
    for c in reverse(r)
        val = get(ROMAN_VALUES, c, 0)
        val == 0 && return 0
        result += val < prev ? -val : val
        prev = val
    end
    result
end

"""
    SpectrumEntry

A parsed element/isotope/charge-range specification, e.g. from "Fe I-III"
or "56Fe 0-2".  Fields:
- `Z`        — atomic number
- `isotope`  — isotope mass number, or `nothing` for natural abundance
- `charges`  — list of ionization stages (1 = neutral, 2 = singly ionized, …),
              or empty vector meaning "all available charge states"
"""
struct SpectrumEntry
    Z::Int
    isotope::Union{Int,Nothing}
    charges::Vector{Int}
end

"""
    parse_spectra(input)

Parse a spectrum specification string into a vector of `SpectrumEntry`.

Accepted formats (comma-separated groups):
  - `"Fe I"`         → Fe, charge 1 (neutral)
  - `"Fe I-III"`     → Fe, charges 1, 2, 3
  - `"Fe0-2"`        → Fe, charges 1, 2, 3  (numeric 0-indexed form)
  - `"Fe"`           → Fe, all available charges (empty charges vector)
  - `"Fe I, Ni I"`   → two elements
  - `"56Fe I"`       → isotope ⁵⁶Fe, charge 1
  - `"all spectra"`  → empty (catch-all)
"""
function parse_spectra(input::AbstractString)
    s = strip(input)
    isempty(s) && return SpectrumEntry[]
    lowercase(s) == "all spectra" && return SpectrumEntry[]

    entries = SpectrumEntry[]
    for group in split(s, r"\s*,\s*")
        group = strip(group)
        isempty(group) && continue

        # Group pattern: optional isotope (digits or *), then element symbol, then charge spec.
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
            # Bare element name → load all available charge states (determined later).
            push!(entries, SpectrumEntry(z, isotope, Int[]))
            continue
        end

        rest = strip(rest)

        # Numeric form: "n" or "n-m" (0-indexed: 0=I, 1=II, …)
        if occursin(r"^\d+(-\d*)?$", rest)
            m2 = match(r"^(\d+)(?:-(\d*))?$", rest)
            a = parse(Int, m2[1]) + 1
            b = (m2[2] !== nothing && !isempty(m2[2])) ? parse(Int, m2[2]) + 1 : a
            push!(entries, SpectrumEntry(z, isotope, collect(a:b)))
        # Roman numeral form: "I", "III", "I-III", etc.
        elseif occursin(r"^[IVLCivlc]+(-[IVLCivlc]+)?$", rest)
            parts = split(rest, r"\s*-\s*")
            if length(parts) == 1
                stage = roman_to_int(rest)
                stage > 0 && push!(entries, SpectrumEntry(z, isotope, [stage]))
            elseif length(parts) == 2
                a = roman_to_int(parts[1])
                b = roman_to_int(parts[2])
                if a > 0 && b >= a
                    push!(entries, SpectrumEntry(z, isotope, collect(a:b)))
                end
            end
        # Dashes-only → negative charge offset (e.g. "---" = stage −3).
        # This is used internally for special cases like electron energies.
        elseif all(c -> c == '-', rest)
            push!(entries, SpectrumEntry(z, isotope, [-length(rest)]))
        end
    end

    entries
end
