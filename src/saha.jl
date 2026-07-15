function ionization_potentials(db::LIBSDB, element::AbstractString)
    tbl = get(db.ionization, element, nothing)
    tbl === nothing && return Dict{Int,Float64}()
    ips = Dict{Int,Float64}()
    for row in Tables.rows(tbl)
        ips[row.spectr_charge] = row.energy
    end
    ips
end

function saha_ion_populations(db::LIBSDB, element::AbstractString, temp, eden;
    charge_min::Integer=0, charge_max=nothing)

    temp_eV = to_internal_temp(temp)
    eden_cm3 = to_internal_density(eden)

    eden_cm3 > 0 || throw(ArgumentError("eden must be positive"))
    temp_eV > 0 || throw(ArgumentError("temp must be positive"))

    ips = ionization_potentials(db, element)
    isempty(ips) && return Dict{Int,Tuple{Float64,Float64}}()

    max_avail = maximum(keys(ips); init=0)
    if charge_max === nothing
        charge_max = max_avail
    end
    charge_max = max(charge_max, 3)

    saha_factor = 6.043e21 * temp_eV^1.5 / eden_cm3

    pfs = Dict{Int,Float64}()
    for sc in charge_min:charge_max
        pfs[sc] = partition_function(db, element, sc, temp_eV)
    end

    saha = Dict{Int,Tuple{Float64,Float64}}()
    prev_log_abund = 0.0
    prev_pf = 1.0
    prev_ip = 0.0
    max_log = -Inf

    temp_cm = temp_eV * 8065.54393734921
    for sc in charge_min:charge_max
        pf = get(pfs, sc, 0.0)
        pf <= 0 && continue
        ip = get(ips, sc, 0.0)
        ratio = saha_factor / prev_pf
        log_abund = log(pf * ratio) + prev_log_abund - prev_ip / temp_cm

        saha[sc] = (log_abund, pf)
        max_log = max(max_log, log_abund)

        prev_log_abund = log_abund
        prev_pf = pf
        prev_ip = ip

        z = element_number(element)
        if z > 0 && z == sc
            sc_next = z + 1
            log_abund = log(saha_factor / prev_pf) + prev_log_abund - prev_ip / temp_cm
            saha[sc_next] = (log_abund, 1.0)
            max_log = max(max_log, log_abund)
            prev_log_abund = log_abund
            prev_pf = 1.0
            prev_ip = 0.0
        end
    end

    total = sum(exp(v[1] - max_log) for v in values(saha))
    log_total = log(total)
    result = Dict{Int,Tuple{Float64,Float64}}()
    for (k, v) in saha
        result[k] = (exp(v[1] - max_log - log_total), v[2])
    end
    result
end
