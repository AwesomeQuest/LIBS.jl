#!/usr/bin/env julia
"""Fetch NIST ASD LIBS reference data for a set of elements."""

using Downloads
using JSON

const BASE_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl" *
    "?composition={elem}%3A100" *
    "&mytext%5B%5D={elem}" *
    "&myperc%5B%5D=100" *
    "&spectra={elem}{lo}-{hi}" *
    "&low_w=200" *
    "&limits_type=0" *
    "&upp_w=600" *
    "&show_av=2" *
    "&unit=1" *
    "&resolution=1000" *
    "&temp=1" *
    "&eden=1e17" *
    "&maxcharge={maxcharge}" *
    "&min_rel_int=0.001" *
    "&int_scale=1" *
    "&libs=1"

const ELEMENTS = [
    ("H", 1, 0, 0),
    ("He", 2, 0, 1),
    ("Li", 3, 0, 1),
    ("C", 6, 0, 2),
    ("N", 7, 0, 2),
    ("O", 8, 0, 1),
    ("Al", 13, 0, 2),
    ("Si", 14, 0, 2),
    ("Ar", 18, 0, 1),
    ("Fe", 26, 0, 2),
    ("Ni", 28, 0, 2),
    ("Cu", 29, 0, 1),
    ("U", 92, 0, 2),
    ("Kr", 36, 0, 1),
    ("Zn", 30, 0, 1),
]

function parse_data_sticks_array(html::AbstractString)
    m = match(r"var dataSticksArray\s*=\s*(\[.*?\])\s*;", html)
    m === nothing && return nothing
    arr_str = m.captures[1]

    rows = String[]
    current = IOBuffer()
    depth = 0
    for ch in arr_str
        if ch == '['
            if depth > 0
                write(current, ch)
            end
            depth += 1
        elseif ch == ']'
            depth -= 1
            if depth == 0
                push!(rows, String(take!(current)))
            else
                write(current, ch)
            end
        else
            if depth >= 1
                write(current, ch)
            end
        end
    end

    isempty(rows) && return nothing

    header = rows[1]
    charge_labels = [m.match for m in eachmatch(r"label:'([^']+)'", header)]

    data = []
    for row_str in rows[2:end]
        parts = split(row_str, ",")
        isempty(parts) && continue
        wl = parse(Float64, strip(parts[1]))
        ints = [p == "null" ? nothing : parse(Float64, p) for p in parts[2:end]]
        push!(data, Dict("wl" => wl, "intensities" => ints))
    end

    stages = String[]
    pops = Float64[]
    for label in charge_labels
        m_stage = match(r"(.+?)\s*\(([\d.eE+-]+)\)", label)
        if m_stage !== nothing && label != "Ritz wavelength (nm)"
            push!(stages, strip(m_stage[1]))
            push!(pops, parse(Float64, m_stage[2]))
        else
            push!(stages, label)
            push!(pops, NaN)
        end
    end

    return Dict(
        "charge_stages" => stages,
        "saha_populations" => pops,
        "data" => data,
    )
end

function main()
    script_dir = @__DIR__
    fixtures_dir = joinpath(script_dir, "fixtures")
    mkpath(fixtures_dir)

    success = 0
    fail = 0

    for (symbol, z, lo, hi) in ELEMENTS
        url = replace(BASE_URL,
            "{elem}" => symbol,
            "{lo}" => string(lo),
            "{hi}" => string(hi),
            "{maxcharge}" => string(hi))
        fixture_path = joinpath(fixtures_dir, "$(lowercase(symbol)).json")
        html_path = joinpath(fixtures_dir, "$(lowercase(symbol)).html")

        print("Fetching $symbol (Z=$z)... ")
        flush(stdout)

        try
            html = String(Downloads.download(url; headers=Dict(
                "User-Agent" => "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
            )))

            write(html_path, html)

            parsed = parse_data_sticks_array(html)
            if parsed === nothing
                println("FAIL (no dataSticksArray)")
                fail += 1
                continue
            end

            if isempty(parsed["data"])
                println("WARN (empty data)")
            else
                nstages = count(s -> !occursin("nm", s), parsed["charge_stages"])
                println("OK ($(length(parsed["data"])) lines, $nstages stages)")
            end

            parsed["element"] = symbol
            parsed["url_params"] = Dict(
                "temp" => 1,
                "eden" => 1e17,
                "resolution" => 1000,
                "low_w" => 200,
                "upp_w" => 600,
                "min_rel_int" => 0.001,
            )

            write(fixture_path, JSON.json(parsed, indent=2))

            success += 1
            sleep(1.5)

        catch e
            println("FAIL ($e)")
            fail += 1
        end
    end

    println("\nDone: $success OK, $fail FAIL")
    return fail == 0 ? 0 : 1
end

main()
