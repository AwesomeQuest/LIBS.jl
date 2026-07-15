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
        endswith(f, ".arrow") || continue
        path = joinpath(data_dir, f)
        if startswith(f, "lines_")
            elem = f[length("lines_")+1:end-length(".arrow")]
            lines[elem] = Arrow.Table(path)
        elseif startswith(f, "levels_")
            elem = f[length("levels_")+1:end-length(".arrow")]
            levels[elem] = Arrow.Table(path)
        elseif startswith(f, "ionization_")
            elem = f[length("ionization_")+1:end-length(".arrow")]
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
