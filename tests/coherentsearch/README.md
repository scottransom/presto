# Reference numbers for the `coherentsearch` port

PRESTO's `coherentsearch` is a C port of the production search in
[CoherentSearch.jl](https://github.com/scottransom/CoherentSearch.jl) v1.0.0. This directory
pins it to the Julia code. Everything here is small and committed, so the checks run in the
ordinary `pixi run test` on every push.

| file | what it is |
|---|---|
| `synth_cs.fft`, `synth_cs.inf` | the input: 2^16 points (dt = 0.5 ms, T = 32.768 s) of white noise plus two pulsars, normalised to unit mean power. One pulsar is at r = 200.37 bins (6.115 Hz, all 60 harmonics below Nyquist). The other is at r = 2911.13 bins (88.84 Hz, only harmonics 1–11 below Nyquist). |
| `make_synth_fft.py` | how that input was made (NumPy, fixed seed). The bytes are committed because NumPy does not promise identical random streams across versions. |
| `gen_fixtures.jl` | writes `fixtures/` by calling CoherentSearch.jl's own functions. |
| `fixtures/` | the reference numbers (below). |
| `cs_oracle.py` | a plain NumPy statement of the algorithm, stage by stage. It is the specification the C code is written against, and its docstring lists the conventions that are easy to get wrong. |
| `cs_dump.c` | a test helper (built by meson as `build-pixi/src/cs_dump`, not installed) that runs the C code on one chunk and writes its stages in the fixture layout. |
| `../test_coherentsearch.py` | checks the oracle, and the C code through `cs_dump`, against every fixture. |

## The fixtures

All of them use the `coherent_search` CLI defaults: `nharms = 60`, `m = 16`, `hidr = 0.5`,
`maxdecim = 6` (decimations k = 1–6), `lofreq = 0.1`, `hifreq = 125`, analytic σ,
`blocksize = 2048` and threshold 8. The global trial grid starts at `r_lo = 0.1 T` and steps by
`lodr = hidr/nharms = 1/120` bin.

**Four fixed chunks** of 160 trials, each in its own directory. Each one starts off a 32-trial
group boundary, so the interpolator's masked end groups are exercised:

| chunk | t0 | what it exercises |
|---|---|---|
| `low` | 0 | bottom of the band: harmonics 1–2 do not fit, so the filled set is not a prefix |
| `slowpsr` | 23574 | the strong pulsar, every harmonic filled |
| `fastpsr` | 348861 | the fast pulsar, only harmonics 1–11 filled |
| `nyqedge` | 654896 | only harmonics 1–5 filled. For k = 6 the decimated trials cross Nyquist after 71 of the 160, and none of its stride-6 harmonics is filled, so σ = 0 and every metric is 0 |

Each chunk directory holds:
- `meta.toml`: `t0`, `n`, `r_lo`, `lodr`, `rstart`, the filled harmonics, and per decimation
  (`[k1]` … `[k6]`) `Hk`, `nbins`, the pruned boxcar `widths`, `nvalid`, `sigma_analytic`,
  `sigma_measured`, `sigma_measured_prod` and the per-width `shape` factors.
- `ftprofs.c64`: the interpolated harmonic stacks, complex64, `n × (nharms+1)`. Each trial's
  stack is contiguous and row 0 (DC) is zero.
- `profs_k1.f32`: the k = 1 profiles (unnormalised c2r inverse FFT), float32, `n × 120`.
- `metric_k{k}_analytic.f64`, `metric_k{k}_measured.f64`: the peak boxcar metric of every
  valid trial (`nvalid` of them) with the analytic and the measured σ, from the **exact** path
  (Float64 weights and profiles; the Julia test suite's own pin precision).
- `metric_k{k}_analytic_prod.f64`: the same from the **production** path (Float32 weights and
  profiles), i.e. what a default Julia run computes. It sits within 2.4e-6 of the exact values,
  which is the scale a single-precision C port should match.

All binary files are raw little-endian with no header.

**One full search**, in `fixtures/search/`:
- `meta.toml`: the trial count and the 295 `chunk_starts`. These are the 2048-trial grid plus a
  forced split wherever a harmonic starts or stops fitting.
- `raw_cands_thresh6.txt`: all 1095 trials above threshold 6 across every decimation, before
  any collapsing (`r`, frequency, metric, `nharm`).
- `default.cohout`: the CLI's output with every default. It holds the two pulsars, S/N 17.95
  (k = 2) and 11.60 (k = 3).

## Regenerating

Do this only when the Julia code or the input changes:

    cd tests/coherentsearch
    julia --project=$HOME/git/CoherentSearch.jl -t 1 gen_fixtures.jl   # ~10 s
    pixi run pytest tests/test_coherentsearch.py                       # from the top level

`fixtures/provenance.toml` records the CoherentSearch.jl commit, Julia and FFTW versions used.
