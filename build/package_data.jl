#!/usr/bin/env julia
"""Package the Arrow data files into a tarball and compute hashes for Artifacts.toml."""

using Pkg.Artifacts
using SHA
using Tar

const PKG_DIR = joinpath(@__DIR__, "..")
const DATA_DIR = joinpath(PKG_DIR, "data")
const ARTIFACT_NAME = "libs_data"
const TARBALL_NAME = "libs_data.tar.gz"

function main()
    isdir(DATA_DIR) || error("data directory not found: $DATA_DIR")

    # Create a temporary directory for the artifact
    hash = create_artifact() do artifact_dir
        for f in readdir(DATA_DIR)
            cp(joinpath(DATA_DIR, f), joinpath(artifact_dir, f))
        end
    end
    git_tree_sha1 = string(hash)
    println("git-tree-sha1: $git_tree_sha1")

    # Create the tarball
    archive_path = joinpath(PKG_DIR, TARBALL_NAME)
    archive_artifact(hash, archive_path)
    println("Tarball created: $archive_path ($(round(filesize(archive_path)/1e6, digits=1))) MB")

    # Compute SHA256 of tarball
    tarball_sha256 = open(archive_path) do io
        bytes2hex(sha256(io))
    end
    println("tarball sha256: $tarball_sha256")

    # Generate Artifacts.toml content
    println()
    println("="^60)
    println("Artifacts.toml entry:")
    println("="^60)
    println("[$ARTIFACT_NAME]")
    println("git-tree-sha1 = \"$git_tree_sha1\"")
    println()
    println("    [[$ARTIFACT_NAME.download]]")
    println("    sha256 = \"$tarball_sha256\"")
    println("    url = \"https://github.com/YOUR_USER/LIBS.jl/releases/download/v0.1.0/$TARBALL_NAME\"")

    println()
    println("="^60)
    println("To write directly to Artifacts.toml, run with --write flag")

    if "--write" in ARGS
        toml_path = joinpath(PKG_DIR, "Artifacts.toml")
        if isfile(toml_path)
            # Update existing
            bind_artifact!(toml_path, ARTIFACT_NAME, hash;
                download_info=[("https://github.com/YOUR_USER/LIBS.jl/releases/download/v0.1.0/$TARBALL_NAME", tarball_sha256)],
                force=true)
        else
            # Create new
            bind_artifact!(toml_path, ARTIFACT_NAME, hash;
                download_info=[("https://github.com/YOUR_USER/LIBS.jl/releases/download/v0.1.0/$TARBALL_NAME", tarball_sha256)])
        end
        println("Artifacts.toml written to $toml_path")
    end

    return hash
end

main()
