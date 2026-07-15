# ---------------------------------------------------------------------------
# Database access layer.
#
# NIST atomic line data is pre-compiled into Arrow tables
# (one file per element: lines_<elem>.arrow, levels_<elem>.arrow,
# ionization_<elem>.arrow).  The tables are memory-mapped on first access
# via a module-level singleton.
# ---------------------------------------------------------------------------

using Artifacts

"""
    LIBSDB

Container for the three Arrow table collections:
- `lines`      — radiative transition data (A-values, wavelengths, energies, …)
- `levels`     — energy level data (term energies, degeneracies, configurations)
- `ionization` — ionization potentials per charge state

Each collection is a `Dict{String, Arrow.Table}` keyed by lowercase element symbol.
"""
struct LIBSDB
    lines::Dict{String,Arrow.Table}
    levels::Dict{String,Arrow.Table}
    ionization::Dict{String,Arrow.Table}
end

"""
    LIBSDB(data_dir)

Construct an `LIBSDB` by scanning `data_dir` for files matching
`lines_*.arrow`, `levels_*.arrow`, and `ionization_*.arrow`.
The element name is extracted from the filename.
"""
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

"""
    _find_data_dir()

Locate the data directory.  First checks a local `data/` directory next to
the package source (for development); otherwise falls back to the
`libs_data` artifact (for installed packages).
"""
function _find_data_dir()
    local_dir = joinpath(@__DIR__, "..", "data")
    isdir(local_dir) && return local_dir
    return artifact"libs_data"
end

# Module-level singleton database reference (lazily initialized in __init__).
const _DB = Ref{LIBSDB}()

function __init__()
    _DB[] = LIBSDB(_find_data_dir())
end

"""
    open_db()

Return the module-level `LIBSDB` singleton, loading it on first access.
This is the primary entry point for all database operations.
"""
function open_db()
    _DB[]
end

"""
    _test_fixtures_dir()

Return the path to the `libs_test_data` artifact directory containing
JSON test fixtures for NIST validation.
"""
function _test_fixtures_dir()
    return @artifact_str("libs_test_data")
end
