using SQLite
using Arrow

const SOURCE_DB = "/home/torfi/nistlibs/newnist/NISTLIBS/asd.sqlite"
const DATA_DIR = joinpath(@__DIR__, "..", "data")

_val(x, default) = x === nothing || ismissing(x) ? default : x

function get_elements(db)
    rows = DBInterface.execute(db, """
        SELECT DISTINCT element FROM ASD_Lines
        WHERE element != '' AND int_rel_calc > 0
        ORDER BY element
    """)
    [r.element for r in rows]
end

function build_element_lines(db, elem)
    sql = """
    SELECT
        t1.spectr_charge,
        t1.vac_wl_num,
        t1.calc_wl_num,
        t1.wl_in_air,
        CAST(NULLIF(t1.A, '') AS REAL) AS A,
        CAST(NULLIF(t1.line_str, '') AS REAL) AS line_str,
        t1.int_rel_calc,
        COALESCE(CAST(t2.energy AS REAL) + 0, 0.0) AS low_energy,
        COALESCE(CAST(t3.energy AS REAL) + 0, 0.0) AS upp_energy,
        COALESCE(t2.g, 0) AS low_g,
        COALESCE(t3.g, 0) AS upp_g,
        COALESCE(t2.conf, '') AS low_conf,
        COALESCE(t2.term, '') AS low_term,
        COALESCE(t3.conf, '') AS upp_conf,
        COALESCE(t3.term, '') AS upp_term
    FROM ASD_Lines t1
    LEFT JOIN ASD_Levels t2 ON t1.low_level_id = t2.level_id
    LEFT JOIN ASD_Levels t3 ON t1.upp_level_id = t3.level_id
    WHERE t1.element = ? AND t1.int_rel_calc > 0
    ORDER BY t1.spectr_charge, t1.vac_wl_num
    """

    rows = DBInterface.execute(db, sql, [elem])
    data = [(;
        spectr_charge = _val(r.spectr_charge, 0),
        vac_wl_num = _val(r.vac_wl_num, 0.0),
        calc_wl_num = _val(r.calc_wl_num, 0.0),
        wl_in_air = _val(r.wl_in_air, 0) == 1,
        A = _val(r.A, 0.0),
        line_str = _val(r.line_str, 0.0),
        int_rel_calc = _val(r.int_rel_calc, 0.0),
        low_energy = _val(r.low_energy, 0.0),
        upp_energy = _val(r.upp_energy, 0.0),
        low_g = _val(r.low_g, 0),
        upp_g = _val(r.upp_g, 0),
        low_conf = _val(r.low_conf, ""),
        low_term = _val(r.low_term, ""),
        upp_conf = _val(r.upp_conf, ""),
        upp_term = _val(r.upp_term, ""),
    ) for r in rows]

    isempty(data) && return
    out = joinpath(DATA_DIR, "lines_$(elem).arrow")
    Arrow.write(out, Arrow.toarrow(data))
    print(" $(length(data))")
end

function build_element_levels(db, elem)
    sql = """
    SELECT
        spectr_charge,
        COALESCE(CAST(energy AS REAL) + 0, 0.0) AS energy,
        g
    FROM ASD_Levels
    WHERE element = ? AND energy != '' AND term != 'Limit' AND g > 0
    ORDER BY spectr_charge, energy
    """

    rows = DBInterface.execute(db, sql, [elem])
    data = [(; spectr_charge = r.spectr_charge, energy = r.energy, g = r.g) for r in rows]
    isempty(data) && return
    out = joinpath(DATA_DIR, "levels_$(elem).arrow")
    Arrow.write(out, Arrow.toarrow(data))
    print(" $(length(data))")
end

function build_element_ionization(db, elem)
    sql = """
    SELECT
        spectr_charge,
        COALESCE(CAST(energy AS REAL) + 0, 0.0) AS energy
    FROM ASD_Levels
    WHERE element = ? AND level_id LIKE '%Lim001%' AND energy != ''
    ORDER BY spectr_charge
    """

    rows = DBInterface.execute(db, sql, [elem])
    data = [(; spectr_charge = r.spectr_charge, energy = r.energy) for r in rows]
    isempty(data) && return
    out = joinpath(DATA_DIR, "ionization_$(elem).arrow")
    Arrow.write(out, Arrow.toarrow(data))
    print(" $(length(data))")
end

mkpath(DATA_DIR)
db = SQLite.DB(SOURCE_DB)
elems = get_elements(db)
println("Building per-element Arrow files for $(length(elems)) elements...")
for (i, elem) in enumerate(elems)
    print("$i/$(length(elems)) $(elem): lines")
    build_element_lines(db, elem)
    print(" levels")
    build_element_levels(db, elem)
    print(" ionization")
    build_element_ionization(db, elem)
    println()
end
println("Build complete.")
