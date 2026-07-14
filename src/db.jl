using Artifacts

struct LIBSDB
    lines::Dict{String,Arrow.Table}
    levels::Dict{String,Arrow.Table}
    ionization::Dict{String,Arrow.Table}
end

function LIBSDB(data_dir::AbstractString)
    lines = Dict{String,Arrow.Table}()
    levels = Dict{String,Arrow.Table}()
    ionization = Dict{String,Arrow.Table}()

    for f in readdir(data_dir)
        path = joinpath(data_dir, f)
        if startswith(f, "lines_") && endswith(f, ".arrow")
            elem = f[7:end-6]
            lines[elem] = Arrow.Table(path)
        elseif startswith(f, "levels_") && endswith(f, ".arrow")
            elem = f[8:end-6]
            levels[elem] = Arrow.Table(path)
        elseif startswith(f, "ionization_") && endswith(f, ".arrow")
            elem = f[12:end-6]
            ionization[elem] = Arrow.Table(path)
        end
    end

    LIBSDB(lines, levels, ionization)
end

function _find_data_dir()
    local_dir = joinpath(@__DIR__, "..", "data")
    isdir(local_dir) && return local_dir
    return artifact"libs_data"
end

const _DB = Ref{LIBSDB}()

function __init__()
    _DB[] = LIBSDB(_find_data_dir())
end

function open_db()
    _DB[]
end

function _test_fixtures_dir()
    return @artifact_str("libs_test_data")
end
