function partition_function(db::LIBSDB, element::AbstractString, spectr_charge::Int, temp)
    temp_eV = to_internal_temp(temp)
    tbl = get(db.levels, element, nothing)
    tbl === nothing && return 0.0
    temp_cm = temp_eV * 8065.54393734921
    z = 0.0
    for row in Tables.rows(tbl)
        if row.spectr_charge == spectr_charge && row.energy >= 0 && row.g > 0
            z += row.g * exp(-row.energy / temp_cm)
        end
    end
    z
end
