# Generate the reference fixtures that pin PRESTO's C `coherentsearch` to
# CoherentSearch.jl (v1.0.0 production code; `src/` is unchanged since then).
#
# Run from this directory, with CoherentSearch.jl's project:
#
#     julia --project=$HOME/git/CoherentSearch.jl -t 1 gen_fixtures.jl
#
# It reads `synth_cs.fft` (see make_synth_fft.py) and rewrites `fixtures/`.  See
# README.md for what is pinned and the file layout.  Everything here calls the
# Julia package's own internals, so the numbers are the production path's, not a
# re-derivation.
#
# Two flavours of every per-trial number are written:
#   * "exact": Float64 interpolation weights and a Float64 profile stage.  This is
#     the Julia test suite's own pin precision (~1e-10 from Eqn. 30 evaluated point
#     by point), so it is the truth both the numpy oracle and the C port are
#     compared against.
#   * "prod": Float32 weights and a Float32 profile stage, i.e. what a default
#     `coherent_search` run actually computes.  Its distance from "exact" is the
#     scale of rounding a single-precision C port should also show.

using CoherentSearch
using FFTW
using LinearAlgebra
using Printf
using TOML

const CS = CoherentSearch

const HERE = @__DIR__
const FFTPATH = joinpath(HERE, "synth_cs.fft")
const OUTDIR = joinpath(HERE, "fixtures")

# The `coherent_search` CLI defaults.  (`search()`'s own keyword default for
# `hifreq` is 100; the CLI passes 125.)
const NHARMS = 60
const M = 16
const HIDR = 0.5
const MAXDECIM = 6
const LOFREQ = 0.1
const HIFREQ = 125.0
const THRESHOLD = 8.0
const BLOCKSIZE = 2048
const SIGMA_CAP = CS._BOXCAR_SIGMA_SAMPLES   # production σ̂ buffers are always this long

# (name, first global trial t0, number of trials n).  Every t0 is deliberately
# NOT a multiple of DIRECT_GROUP_V = 32, so both partial end groups are masked.
const CASES = [
    ("low", 0, 160),            # bottom of the band: harmonics 1-2 do not fit
    ("slowpsr", 23574, 160),    # r = 200.37 pulsar, every harmonic filled
    ("fastpsr", 348861, 160),   # r = 2911.13 pulsar, harmonics 12-60 past Nyquist
    ("nyqedge", 654896, 160),   # k = 6 decimated trials cross Nyquist mid-chunk
]

function writebin(path, a)
    open(path, "w") do io
        write(io, a)
    end
end

function git_describe(dir)
    try
        return strip(read(`git -C $dir describe --tags --always --dirty`, String))
    catch
        return "unknown"
    end
end

# All per-k quantities for one chunk, from a filled workspace.
function per_k(ws, params, ft, k, rstart, lodr, n, ::Type{P}) where {P}
    if k == 1
        prof = ws.profs
        Hk = params.nharms
        widths = ws.bcwidths
    else
        db = ws.decims[findfirst(d -> d.k == k, ws.decims)]
        mul!(db.dprofs, db.brfftplan, db.src)
        prof = db.dprofs
        Hk = db.Hk
        widths = db.bcwidths
    end
    nbins = 2Hk
    # Decimated trials whose fundamental k*r reaches Nyquist are skipped, and they
    # form a suffix of the chunk (see `decim_pass!`).  The base pass scores all n.
    nyq = ft.N / 2
    nvalid = k == 1 ? n : count(j -> k * (rstart + (j - 1) * lodr) < nyq, 1:n)
    sig_a = CS._analytic_sigma(ws.filled, k, Hk)
    buf = Vector{P}(undef, SIGMA_CAP)
    sig_m = Float64(CS._block_sigma(prof, nbins, nvalid, buf))
    shape = zeros(Float64, length(widths))
    CS._boxcar_shape!(shape, nbins, widths, ws.filled, k, Hk)
    psum = Vector{P}(undef, nbins + widths[end] + 1)
    inv(s) = s > 0 ? one(P) / P(s) : zero(P)
    ma = [Float64(CS._profile_boxcar(prof, j, psum, widths, nbins, inv(sig_a), shape))
          for j in 1:nvalid]
    mm = [Float64(CS._profile_boxcar(prof, j, psum, widths, nbins, inv(sig_m), shape))
          for j in 1:nvalid]
    return (; Hk, nbins, widths, nvalid, sig_a, sig_m, shape, ma, mm, prof)
end

function chunk_case(ft, name, t0, n)
    lodr = HIDR / NHARMS
    r_lo = LOFREQ * ft.T
    rstart = r_lo + t0 * lodr
    decims = decimation_set(NHARMS, MAXDECIM)
    meta = Dict{String,Any}("name" => name, "t0" => t0, "n" => n, "r_lo" => r_lo,
                            "lodr" => lodr, "rstart" => rstart)
    dir = joinpath(OUTDIR, name)
    mkpath(dir)
    for (flavour, P) in (("exact", Float64), ("prod", Float32))
        params = SearchParams(nharms=NHARMS, m=M, hidr=HIDR, decimations=decims,
                              precision=(P === Float64 ? :f64 : :f32))
        dplans = CS.build_direct_plans(P, params, r_lo)
        ws = CS.with_plan_rigor(FFTW.ESTIMATE) do
            CS.Workspace(P, params, n)
        end
        CS.fill_chunk_profiles!(ws, dplans, ft, params, rstart, lodr, n; t0=t0)
        if flavour == "exact"
            meta["filled"] = [h for h in 1:NHARMS if ws.filled[h]]
            # (nharms+1, n) column-major: each trial's harmonic stack (DC row
            # first, always zero) is contiguous.  Stored as complex64.
            writebin(joinpath(dir, "ftprofs.c64"), ComplexF32.(ws.ftprofs[:, 1:n]))
        end
        for k in decims
            r = per_k(ws, params, ft, k, rstart, lodr, n, P)
            key = "k$k"
            if flavour == "exact"
                meta[key] = Dict{String,Any}(
                    "Hk" => r.Hk, "nbins" => r.nbins, "widths" => r.widths,
                    "nvalid" => r.nvalid, "sigma_analytic" => r.sig_a,
                    "sigma_measured" => r.sig_m, "shape" => r.shape)
                # (nbins, n) column-major: each trial's profile contiguous.  Only
                # k = 1: a decimated stack is the stride-k rows of `ftprofs`, so its
                # profiles add bytes but no information.
                k == 1 &&
                    writebin(joinpath(dir, "profs_$key.f32"), Float32.(r.prof[:, 1:n]))
                writebin(joinpath(dir, "metric_$(key)_analytic.f64"), r.ma)
                writebin(joinpath(dir, "metric_$(key)_measured.f64"), r.mm)
            else
                meta[key]["sigma_measured_prod"] = r.sig_m
                writebin(joinpath(dir, "metric_$(key)_analytic_prod.f64"), r.ma)
            end
        end
    end
    open(joinpath(dir, "meta.toml"), "w") do io
        TOML.print(io, meta; sorted=true)
    end
    return nothing
end

function search_case(ft)
    dir = joinpath(OUTDIR, "search")
    mkpath(dir)
    decims = decimation_set(NHARMS, MAXDECIM)
    params = SearchParams(nharms=NHARMS, m=M, hidr=HIDR, decimations=decims)
    lodr = HIDR / NHARMS
    r_lo = LOFREQ * ft.T
    r_hi = HIFREQ * ft.T
    total = max(0, floor(Int, (r_hi - r_lo) / lodr) + 1)
    cstarts = CS.chunk_starts(CS.build_direct_plans(params, r_lo), ft, params, total,
                              BLOCKSIZE)
    # Raw (uncollapsed) candidates at a lower threshold, to pin the trials the
    # gate lets through, then the CLI's own default output.
    raw = search(ft, params; lofreq=LOFREQ, hifreq=HIFREQ, threshold=6.0,
                 blocksize=BLOCKSIZE, remove=false, harm_remove=false, wisdom=false)
    sort!(raw; by=c -> (c.r, c.nharm))
    open(joinpath(dir, "raw_cands_thresh6.txt"), "w") do io
        println(io, "# r(bins)  freq(Hz)  metric  nharm  -- uncollapsed, threshold 6")
        for c in raw
            @printf(io, "%.17g %.17g %.17g %d\n", c.r, c.freq, c.metric, c.nharm)
        end
    end
    meta = Dict{String,Any}("total" => total, "r_lo" => r_lo, "r_hi" => r_hi,
                            "lodr" => lodr, "blocksize" => BLOCKSIZE,
                            "chunk_starts" => cstarts, "n_raw_thresh6" => length(raw))
    open(joinpath(dir, "meta.toml"), "w") do io
        TOML.print(io, meta; sorted=true)
    end
    # The candidate file exactly as the CLI writes it with every default.
    CS.main(["--nowisdom", "--noprogress", "-o", joinpath(dir, "default.cohout"),
             FFTPATH])
    return nothing
end

function main()
    ft = FFTFile(FFTPATH)
    mkpath(OUTDIR)
    for (name, t0, n) in CASES
        chunk_case(ft, name, t0, n)
    end
    search_case(ft)
    prov = Dict{String,Any}(
        "coherentsearch_jl" => git_describe(pkgdir(CoherentSearch)),
        "julia" => string(VERSION),
        "fftw" => string(FFTW.version),
        "input" => "synth_cs.fft",
        "nharms" => NHARMS, "m" => M, "hidr" => HIDR, "maxdecim" => MAXDECIM,
        "lofreq" => LOFREQ, "hifreq" => HIFREQ, "threshold" => THRESHOLD,
        "sigma_cap" => SIGMA_CAP,
        "cases" => [c[1] for c in CASES])
    open(joinpath(OUTDIR, "provenance.toml"), "w") do io
        TOML.print(io, prov; sorted=true)
    end
end

main()
