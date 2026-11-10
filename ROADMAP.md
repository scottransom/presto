# PRESTO Roadmap

Notes on significant planned changes. This is a living, aspirational document — items may
change, be reordered, or be dropped. For stable "how to work in this repo" guidance, see
[CLAUDE.md](CLAUDE.md); for what has already shipped, see [CHANGELOG.md](CHANGELOG.md) and the
git history; for tracked work, see GitHub issues.

## Where things stand (late 2026)

The 2026 goal — **PRESTO on conda-forge** — is done. The package is
[`presto-pulsar`](https://github.com/conda-forge/presto-pulsar-feedstock) (the plain `presto`
name is taken on bioconda by an unrelated tool). It shipped with v6.0.1 (2026-08-18) and is
built for linux-64, linux-aarch64, osx-64 and osx-arm64 (native arm runners since v6.1.1).
`RELEASE.md` is the release-to-conda-forge checklist.

The prerequisite work that made clean packaging possible is also complete:

- **No required environment variables.** Runtime data installs to `{prefix}/share/presto`;
  `PRESTO` is now only an optional override (`src/datadir.c`, `python/presto/_datadir.py`).
  `makewisdom` is installed and writes wisdom to the same place.
- **Fortran-free compiled code.** FFTFIT is pure Python (`classic` and `aarchiba` methods),
  and the last Fortran (`dgels`) was replaced with GSL. PGPLOT is still an external Fortran
  dependency.
- **No TEMPO.** Barycentering uses ERFA; polycos come from `tempo2 -tempo1`.
- **ERFA and GSL in place of vendored/NR code:** SLALIB date routines, `amoeba.c`, `median.c`,
  randlib, and several `misc_utils.c` helpers.
- **Man pages** for every clig-based tool, installed to `{prefix}/share/man/man1`.
- **Dead code removed:** `dftfold`, `a2x`, `cand_output.c`, `solvopt.c`/`apprgrdn.c`, the
  unbuildable legacy C tests.

The next big goal, a C port of CoherentSearch.jl, **shipped in v6.2.0** as `coherentsearch`
(plus `coherent_sift.py`), merged from the `coherent-search` branch on 2026-11-10. What is
left of it is optional: Highway kernels (phase 4's open items and phase 5 below) and
backporting the C version's additions to Julia (see "Other open items").

## The coherent harmonic search in PRESTO (shipped in v6.2.0)

### What and why

[CoherentSearch.jl](https://github.com/scottransom/CoherentSearch.jl) (v1.0.0, 2026-10-08; the
version the submitted paper describes) is a periodicity search over PRESTO `.fft` files that
keeps the *phase* of every harmonic. At each trial frequency it Fourier-interpolates up to 60
harmonics (Eqn. 30 of Ransom, Eikenberry & Middleditch 2002), inverse-FFTs them into a pulse
profile, and matched-filters that profile against a boxcar bank. A decimation ladder reuses the
same amplitudes for shallower folds. In its injection Monte Carlo it detects 76% of white-noise
injections against `accelsearch -zmax 0`'s 42% and riptide's 50–71%, with a calculable
false-alarm rate. Single-threaded it is faster than riptide's `rseek`.

It belongs in PRESTO because its input *is* PRESTO's pipeline (`realfft` → `rednoise` →
normalized `_red.fft`; the analytic noise scale assumes exactly that), and its candidates go
to `prepfold`. Most PRESTO users will never install Julia. Porting also removes Julia's ~1 s
per-invocation start-up and the sysimage/precompile machinery.

**This is a reach and packaging play, not a speed play.** The Julia code is already near
LLVM-quality, and ~45% of its runtime is FFTW, which a C port calls identically. Expect speed
parity (roughly ±20%), not a win.

Julia stays the reference implementation, the research platform, and the home of the GPU
(CUDA.jl) path. The C port follows it; it does not replace it.

### Language: C, with the SIMD kernels in C++ (Highway)

Plain C for everything except the two or three hand-shaped SIMD kernels. Templates were the
main reason to consider C++, and the port only needs the production single-precision path
(the `Float64` path exists for Julia's equivalence pins, which stay in Julia). A fully
independent C++ code was considered and rejected: PRESTO already provides the file formats,
the `.inf` reader, FFTW/wisdom handling, clig/pyclig, PGPLOT, meson, and a conda-forge
package with the users.

### SIMD on conda-forge (benefits all of PRESTO, not just this)

conda-forge's x86-64 builds target a baseline CPU (no AVX2/FMA), whereas Julia JIT-compiles for
the host. FFTW already dispatches at runtime, so only our own loops are affected. Options, in
the order to pursue them:

0. **Measure the ceiling first (free).** Build PRESTO at the default flags and again with
   `-march=x86-64-v3`, then time `accelsearch`, `prepsubband`, `rfifind` and `prepfold` on real
   data. That bounds what any SIMD work could buy the existing tools.
   - **accelsearch, measured 2026-11-09.** On the Ter5 test file (DM 242, downsampled ×2,
     5.28 M points, `_red.fft`), on an i7-10510U (Comet Lake, AVX2, no AVX-512), with
     conda-forge's flags (`-march=nocona -mtune=haswell -O2`), with and without
     `-march=x86-64-v3`, 1 thread, best of 2–3 runs (`hyperfine`):

     | search | baseline (6.1.1) | baseline + `fmaf` fix | x86-64-v3 + fix |
     |---|---|---|---|
     | `-zmax 200` | 73.2 s | 29.9 s (2.45×) | 26.4 s (2.8×) |
     | `-zmax 50 -wmax 200` | 297.1 s | 126.3 s (2.35×) | ≈124 s |

     **Most of the first-measured v3 gain was a bug, not the CPU level.** A gperftools profile
     (conda-forge `gperftools` via `LD_PRELOAD`, `go tool pprof`; no root needed) put ~70% of
     the baseline's time in `subharm_fderivs_vol` and libm: its inner loops called `fmaf()`
     explicitly, and on a target without hardware FMA (x86-64 baseline) each call is an
     exact-rounding libm routine that also blocks vectorization. Plain multiply-adds (which
     the compiler still fuses to FMA where available) fixed it in 6.1.2. ARM always has FMA,
     so it was never affected. All ACCEL files are byte-identical across the four builds.
   - **What v3 still buys accelsearch is ~10%.** With the fix, the profile at 1 thread
     (~27 s) is roughly: FFTW 8.5 s, `search_ffdotpows` 6.9 s, `add_ffdotpows` 5.6 s,
     `subharm_fderivs_vol` 5.9 s. `gcc -fopt-info-vec` shows that neither `add_ffdotpows`
     (the `rinds[kk]` gather in the harmonic sum) nor `search_ffdotpows` (the threshold scan,
     an `omp parallel for collapse(3)` that recomputes indices per element) vectorizes at
     any level. Those two (~45% of the remaining time) are the next targets; restructuring
     the scan's loop nest may be enough without explicit SIMD.
   - **The other tools gain little from v3** (same host and file, 1 thread): `rfifind -time 2`
     ~20.5 s both ways; `prepsubband` (64 DMs, masked) ~17 s both ways (I/O-bound);
     `prepfold` on the raw data with a DM/P search 96 → 85 s (1.13×). Outputs agree to float
     rounding (≤ 2e-7 relative; the rfifind mask is byte-identical).
   - Watch out when benchmarking: the tools use `DT_RPATH` with the conda env's `lib/` first,
     which overrides `LD_LIBRARY_PATH`, so install each build into its own prefix or both
     will silently load the same `libpresto`.
1. **conda-forge microarchitecture variants — tried and dropped.** Adding the
   `x86_64-microarch-level` metapackage to the feedstock's build requirements would produce
   v1 and v3 builds, with the solver choosing via the `__archspec` virtual package
   (conda-forge/presto-pulsar-feedstock#7 did this). After the `fmaf` fix, v3 gains ~10% on
   accelsearch and less elsewhere, which does not justify doubling the x86 builds; conda-forge
   also prefers runtime dispatch to variants. Closed unmerged. Revisit only if a profile shows
   a large auto-vectorization gap that runtime-dispatched kernels (below) do not cover.
2. **Google Highway for hand-written kernels.** [Highway](https://github.com/google/highway) is
   portable explicit SIMD with built-in runtime dispatch (`HWY_DYNAMIC_DISPATCH`) across SSE4,
   AVX2, AVX3 (Skylake-X), AVX3_DL, AVX3_ZEN4, AVX3_SPR, NEON, SVE/SVE2 and RVV. On conda-forge
   as `libhwy` (1.4.0, all six platforms); also in Debian and Homebrew.
   - *Why it fits:* the CoherentSearch.jl optimization log is largely a record of the
     auto-vectorizer surprising us. The worst case was one `vscatterqps` in the boxcar tile
     transpose, which dropped Skylake-SP to AVX-512 licence level 2 for the *whole* search
     (2.84 → 2.06 GHz). With explicit `Load`/`MulAdd`/`Store`, a scatter appears only if
     written. Target selection can exclude plain AVX3 (Skylake-SP) while keeping the Zen 4 and
     Sapphire Rapids targets. SVE reaches Graviton/Grace cluster nodes, which neither the
     conda-forge baseline nor NEON-only code serves well. `contrib/` adds VQSort and vector
     math.
   - *Costs:* C++17. Confine it to a `src/simd/` directory of `.cc` files with `extern "C"`
     entry points so all data structures and libpresto's ABI stay C (libpresto would then
     link libstdc++, and meson needs `project(..., 'c', 'cpp')`). Highway code is verbose,
     each kernel compiles once per target, and scalable vectors (SVE/RVV) complicate the
     fixed 32-trial accumulator layout (use `FixedTag`/`CappedTag`, or unroll on `Lanes(d)`).
3. **Fallback if C++ proves unwelcome: plain-C dispatch.** Compile the kernel files two or
   three times in meson with different flags (baseline, `-mavx2 -mfma`) and pick function
   pointers once at startup via `__builtin_cpu_supports`. Works everywhere (including macOS,
   unlike `target_clones`, which needs ifunc), but each variant still relies on the
   auto-vectorizer.

**Plan:** (0) is done and (1) dropped; do (2) for the new
search's hot kernels: the grouped interpolation matvec and the transposed boxcar
gate/prefix-sum scan. Then adopt Highway kernels in `accelsearch` *only where profiling says
the non-FFTW loops matter* (see `fourier_interpolation_tradeoffs.md` §5 in CoherentSearch.jl).
The profile now names two `accelsearch` targets: the harmonic sum in `add_ffdotpows`, whose
`inpows[rinds[kk]]` gather is a near-arithmetic index sequence that could become contiguous
loads, and the threshold scan in `search_ffdotpows`.

### Scope of the port

**In:** the production search from CoherentSearch.jl v1.0.0, single precision throughout:
- direct interpolation plans (`DirectPlan`, `build_group_table`) and the V=32 group kernel;
- base and stride-`k` decimated inverse FFTs (`fftwf_plan_many_dft_c2r` with `istride = k`);
- analytic σ (default), `--sigma measured` (quickselect MAD — the red-noise fallback), and the
  band-limited boxcar normalization (`_boxcar_shape!`) plus the σ sanity check;
- the `Float32` batched boxcar gate with exact rescore above `exactcut`, the ladder-pruned
  width bank, and chunk starts / Nyquist give-up;
- duplicate and harmonic collapse, a candidate output file, and multi-file invocations;
- OpenMP over whole chunks, round-robin (`schedule(static,1)`), one private workspace per
  thread, all FFTW planning before the parallel region;
- **PGPLOT candidate plots** (`/CPS`, `/PNG`), drawn per file as each search finishes, plus a
  small replot tool that reads a saved candidate file. Profile reconstruction
  (`candidate_profile`) is microseconds per candidate.
- `bin/sift_candidates.py` (pure stdlib) moves over unchanged.

**Out:** the metric-statistics/normalize histograms, the reference/oracle paths, the
`Float64` profile path, the toy search, and the GPU.

### Validation

The Julia correctness chain (Python oracle → Julia reference → Julia production) gets one more
link: Julia v1.0.0 production → C.

- **Not byte-identical to Julia, by design:** FMA contraction and summation order differ, and
  FFTW `MEASURE` plans differ between processes.
- **Pin:** per-trial metrics from Julia's `chunk_metrics` on fixed chunks, committed as small
  fixtures, against a C dump mode at about 1e-6 relative. Candidate sets must be identical
  above threshold plus the gate margin on a few real files (e.g. `PM0063`, NGC6624).
- **Within C,** keep the Julia standard: deterministic and chunk-invariant output, checked by
  `diff` on the candidate file.
- **Per kernel:** re-sweep the constants tuned under Julia/LLVM (`V = 32`, `_BC_BATCH = 128`,
  `BJ = 8`) on fitzroy (Skylake-SP), foops (Comet Lake) and an aarch64 host, and audit the
  disassembly of every hot kernel for scatters and 512-bit code.
- **Watch item:** FFTW's c2r overwrites its input by default, and every decimated pass reads
  `ftprofs` after the base pass. Use `FFTW_PRESERVE_INPUT` (supported for rank-1 c2r) or a
  separate buffer, and check what FFTW.jl does rather than assume.

### Prerequisite: CI and a real test suite (do this first, on `master`)

PRESTO has had no CI since the stale Travis setup was removed in 2026-07. The port needs a
harness anyway, and that harness is most of the CI. This work is independent of the port and
benefits everything, so it goes straight to `master`:

1. **[DONE]** GitHub Actions + pixi (`.github/workflows/ci.yml`, `prefix-dev/setup-pixi`)
   runs `pixi run test` on linux-64, linux-aarch64 and osx-arm64. The test list is explicit in
   `pixi.toml`: `test_barycenter.py`, `test_binary_calcs.py` (rewritten from a plotting demo
   into two real checks), `test_maximize_rz.py`, the fftfit suite, and
   `test_presto_python.py`.
   - Still to triage, item by item: the stale scripts in `tests/` (`dedisp.py`,
     `test_max_rzw_harmonics.py`, `tests/python/`), the old f2py-era demo
     `python/fftfit_src/test_fftfit.py` (superseded by the real fftfit tests beside it), and the
     legacy C `.mak` test programs.
     Convert each to a real test, keep it as a demo in `examplescripts/`, or delete it.
2. **pytest as the single runner**, driving executables through `subprocess`.
3. **[DONE]** End-to-end tests on synthetic data generated in-test
   (`tests/test_pipeline_synthetic.py`, generator `tests/synthpulsar.py`): a 23.4567 Hz,
   10%-duty, S/N 12, DM 50 pulsar in a 64-channel 8-bit filterbank, through prepdata and
   prepsubband (4 DM trials), realfft, rednoise, `accelsearch -zmax 0` and prepfold (on the
   `.dat` and on the raw filterbank). Results are checked against the injected values with
   tolerances, not bytes. 1-, 2- and 4-bit files are checked channel by channel against both
   readers, and the pulsar must still be found at those depths. Runs in ~6 s. PSRFITS,
   rfifind and `accelsearch` with `-zmax > 0` are not covered here, but the real-data tests
   (item 5) cover all three.
4. **[DONE]** The aarch64 signed-8-bit check is a permanent test (the round-trip of every
   sample format through `prepdata -dm 0`), and CI also builds on linux-64 with
   `-funsigned-char`. With the pre-6.1.1 `char` code put back, that job fails exactly the
   signed test, so the bug class is caught on x86 too.
5. **[DONE] Scheduled real-data tests** (`.github/workflows/realdata.yml`, `pixi run
   test-realdata`, ~7 min): weekly and on demand, on linux-64, linux-aarch64 and osx-arm64.
   The 173 MB Terzan 5 PSRFITS file lives on the repo's `test-data-v1` release, not in git
   (`tests/prepfold/get_test_data.sh` fetches and checksums it; new data go under a new tag).
   `tests/test_realdata_search.py` runs rfifind → prepsubband → realfft → rednoise →
   accelsearch → sifting → prepfold and `single_pulse_search.py`, checked against the
   Ter5A/Ter5N timing solutions (Doppler-shifted f and z, DM, and Ter5A's single pulses
   phase-locked to its ephemeris). Then come the `prepfold_tests.sh` and
   `prepfold_multi_tests.sh` regression suites. Possible additions: a second
   dataset in another format (e.g. the 4-bit `GBT_Lband_PSR.fil` from the tutorial, on the
   same release), `-zmax > 0` folding of a `.dat` from an accelerated candidate, and the
   `-ncpus` determinism check from `tests/prepfold/determinism_check.sh`.

The feedstock's own tests remain the release-time check; CI catches problems before they get
there.

### Workflow: a long-lived feature branch

This departs from the usual direct-to-`master` workflow (see CLAUDE.md), because it is exactly
the "very large-scale change" that rule reserves branches for: a new executable, a new
dependency (Highway), and C++ in the build for the first time.

- CI (above) and the `-march=x86-64-v3` measurement went to `master` as usual (done, and
  the feedstock microarch change was dropped) — neither needed the port.
- The port lives on a branch, e.g. `coherent-search`, rebased on `master` periodically so it
  runs under CI from day one.
- Merge into `master` via a PR once phase 3 below passes its pins on all CI platforms; a
  release with the new tool follows. **Done:** fast-forwarded into `master` on 2026-11-10
  (0050861) once CI and the real-data tests passed on every platform, and released as
  v6.2.0.

### Phases

0. **[DONE] Reference numbers from Julia** (`tests/coherentsearch/`, `tests/test_coherentsearch.py`,
   in `pixi run test`). A committed 256 KB synthetic `.fft` with two pulsars;
   `gen_fixtures.jl` dumps CoherentSearch.jl's own numbers stage by stage for four 160-trial
   chunks at the CLI defaults: filled harmonics, interpolated stacks, profiles, and for each k
   = 1–6 the width banks, analytic and measured σ, shape factors and every per-trial metric
   (exact Float64 path, plus the Float32 production path, which is within 2.4e-6). It also
   dumps a full search's 295 chunk starts, its 1095 raw candidates above threshold 6 and the
   default `.cohout`. `cs_oracle.py`, a plain NumPy statement of the algorithm, reproduces
   all of it (metrics to 3e-14). It is the specification the C code follows, and the C
   dump mode will be checked against the same fixtures.
1. **[DONE] Single-threaded scalar C for one chunk** (`include/coherentsearch.h`,
   `src/cs_interp.c`, `src/cs_metric.c`, `src/cs_chunk.c`; for now linked into executables,
   not libpresto). It covers the exact-rational interpolation plans (and the non-rational
   fallback), the batched FFTW c2r with `FFTW_PRESERVE_INPUT`, and every decimation k with
   its pruned width bank. It computes the analytic and measured σ, the shape factors and the
   per-trial boxcar metric, all in float32 like Julia's production path. More than planned:
   the decimated passes, pruned widths and measured σ came along because they are only a
   strided view of the same stacks. The `cs_dump` test helper writes one chunk's stages;
   `tests/test_coherentsearch.py` checks them against the Julia fixtures. Integers, analytic
   σ and shapes match exactly; stacks match to 2.8e-7, measured σ to 2e-7 and metrics to
   2.9e-6 (Julia's own float32 path is 2.4e-6 from exact). Clean under valgrind.
   **Not done, deliberately:** the batched float32 *gate*. In scalar C it would be the same
   arithmetic as the exact rescore, so it only pays as a SIMD kernel and moves to phase 4.
   The plain per-trial interpolation loop also replaces Julia's trials-axis group kernel
   (`build_group_table`); that is a SIMD layout, likewise phase 4.
2. **[DONE] The full search** (`src/cs_search.c`, `src/coherentsearch.c`,
   `clig/coherentsearch_cmd.cli`, `docs/coherentsearch.1`, `src/Makefile`). It covers the
   chunk loop with `chunk_starts` (forced splits where a harmonic starts or stops fitting),
   the analytic-σ sanity check, and candidate collapse (duplicates, then harmonics, both
   with stable sorts as in Julia). Duty cycles come from an m = 64 refold (`measure_ducy`),
   and the output is the `.cohout` file. Candidates are a plain growable array, not
   accelsearch's `GSList`: the search collects them and then sorts and collapses them in
   passes, and phase 3's per-thread lists merge more simply. On the synthetic file the
   chunk starts are identical, so are the 1095 raw candidates above threshold 6
   (bit-identical r, metrics within 3.8e-6), and the default `.cohout` is
   **byte-identical** to the Julia CLI's. 1.3 s single-threaded (foops) for 491k trials ×
   6 folds. Clean under valgrind.
   Differences from the Julia CLI: it always writes `<input>.cohout` (beside the input, or
   in `-outdir`; `-o` for one file) instead of stdout for a single file, PRESTO-style
   single-dash options, `-measuredsigma` for `--sigma measured`, and `-lobin` takes effect
   whenever it is given rather than whenever it is not 100. Left out on purpose:
   `--precision f64`, `--metricstats`, `--normalize` and the GPU (see "Out" above).
   **Known-pulsar labels (added after phase 4's first round):** `cs_identify_known()`
   runs `comp_psr_to_cand()` on each reported candidate, with accelsearch's matching (a
   catalogue pulsar's harmonics 1–40 within 5 bins, or 0.1% for binaries, plus ±20-bin
   "SL" sidelobe matches, within 10 beam diameters of the pointing). A "Known PSR"
   column appears in the `.cohout` only when something matches (`sift_candidates.py`
   reads just the first six columns), and the plots show the name in bold red.
   On PM0063 the candidate is labelled PSR J1020-6026. On the Terzan 5 test
   observation (four DMs, 236–242) it labels Ter5A and its 2nd–6th harmonics, Ter5C
   and Ter5N. With short observations the sidelobe matches are loose (±20 bins is
   ±0.6 Hz at T = 33 s). That Ter5 data also trips the analytic-σ warning (a factor
   of 2.35 at 0.1 Hz), with spurious 0.11–0.12 Hz candidates. Julia reports the
   identical numbers, so it is the rednoised spectrum's lowest bins, not the port.
3. **[DONE] Throughput and plots.**
   - **OpenMP** (`-ncpus`) over whole chunks, round-robin, with one workspace per thread
     built up front. Each chunk keeps its own candidates and they are joined in chunk
     order, so output is byte-identical for any thread count (checked by a test).
   - **Multi-file runs** reuse the workspaces and FFTW plans.
   - **Wisdom:** the batched c2r plans come from FFTW wisdom when it exists, and
     `makewisdom` now stores them for the default `-nharms`/`-blocksize`, sharing
     `include/cs_fftplan.h`. Otherwise they use FFTW_ESTIMATE, never run-time MEASURE,
     whose timing-dependent choices would change results in the last bits from run to
     run. Wisdom makes the transforms 1.8x faster.
   - **Plots:** `-plot` draws PGPLOT pages of the candidates' profiles, refolded at full
     depth with an m = 64 kernel, 15 to a portrait page; this is `src/cs_plot.c`.
     When Ghostscript is available it also renders the PostScript into PNGs, as prepfold does,
     because PGPLOT's PNG driver draws poorly (and valgrind finds uninitialised reads
     inside it). `-replot` redraws from an existing `.cohout`, which stands in for a
     separate replot tool.
   - **Timings (foops):** on PM0063 (8.4 M points, 0.1–125 Hz), C and Julia give the
     same single candidate (S/N 11.65, 7.1185 Hz) from the same 108 trials. C takes
     **92.9 s on 1 thread and 40.6 s on 4**, against Julia's 35.2 s and 17.9 s. The C
     profile is 45% boxcar scan, 30% interpolation and 24% FFTW. Those are exactly the
     three places Julia uses SIMD layouts (the batched float32 gate, the trials-axis group
     kernel) and PATIENT wisdom, so phase 4 is where the C catches up.
   *Merge point* once phase 4 brings the speed to at least Julia's.
4. **SIMD.** *In progress, plain C first (Highway deferred to see how far that gets).*
   - **Done:** the two batched layouts from Julia, written as plain C loops that gcc
     auto-vectorizes at conda-forge's baseline flags (SSE). Both give results
     **bit-identical** to the scalar code, since each lane does the same float32 operations
     in the same order, so no gate or exact rescore is needed.
     - The boxcar scan transposes `CS_BC_BATCH = 64` profiles into a tile and runs the
       recurrence across them (2.5x on that kernel).
     - The interpolation uses trials-axis groups of `CS_GROUP_V = 16` (Julia's
       `build_group_table`): one contiguous matrix-vector product per group.
   - `makewisdom` now always plans these transforms with PATIENT. EXHAUSTIVE takes 20x
     longer to plan and is no faster.
   - **PM0063, full band (foops):** 92.9 → **41.2 s** on 1 thread (Julia 35.2 s), and
     40.6 → **19.6 s** on 4 (Julia 17.9 s). The `.cohout` is unchanged.
   - **Measured and rejected:**
     - V = 32 groups (slower with both SSE and AVX2);
     - a single V-float GCC vector type, or an array of 4-float vector types, for the
       group accumulators (slower; gcc lowers them badly);
     - phase-outer/width-inner boxcar loops;
     - `-O3` (1.65x *slower*).
     An `x86-64-v3` build is only ~10% faster, so these kernels are not mainly limited by
     vector width.
   - **Profile now:** boxcar ~37%, interpolation ~31%, FFTW ~29% (at PATIENT, so near
     its floor).
   - **eiger** (Xeon w5-3433, Sapphire Rapids, AVX-512; `-march=native` system gcc 14.2),
     PM0063 full band. Identical `.cohout` from both codes at every thread count:

     | threads | C | Julia |
     |---|---|---|
     | 1 | 20.2 s | 19.4 s |
     | 4 | 5.80 s | 5.74 s |
     | 8 | 3.25 s | 3.54 s |
     | 16 | 1.87 s | 2.54 s |

     Interleaved runs (0.1–30 Hz, 1 thread, ±1%) on that machine: baseline flags 6.43 s
     and native 6.44 s, with batch sizes 32 / 64 / 128 at 6.59 / 6.44 / 7.56 s, and
     `CS_GROUP_V = 32` at 5.72 s. So wider vectors alone buy nothing, but AVX-512's 32
     registers fit 32-trial groups. `CS_GROUP_V` is now 32 when `__AVX512F__` is
     defined and 16 otherwise; conda-forge's baseline builds keep 16.
   - **fitzroy** (Xeon Silver 4114, Skylake-SP, AVX-512; `-march=native` gcc 14.2; a
     RAID resync was running, so ±5%), 0.1–30 Hz, 1 thread: baseline flags 12.6 s;
     native 11.2 s with V = 16 and 10.9 s with V = 32; batch sizes 32 and 128 at 11.9
     and 12.2 s. So V = 32 under AVX-512 and a batch size of 64 hold here too. Full band
     against Julia, identical output:

     | threads | C | Julia |
     |---|---|---|
     | 1 | 38.6 s | 45.3 s |
     | 10 | 4.94 s | 6.38 s |
     | 20 | 2.81 s | 5.08 s |

   - Still open: an aarch64 host (NEON has 32 registers, so V = 32 may win there too),
     disassembly audits, and Highway only if a kernel still needs it.
5. **Afterwards, optional:** profile `accelsearch` and move any loops worth it onto Highway.

### Open questions

- The candidate-file format: `.cohout` for now (what CoherentSearch.jl and
  `sift_candidates.py` use); something `ACCEL_sift`/`prepfold`-friendly could follow. The
  tool's name is decided: **`coherentsearch`**.
- Whether to accept C++ in libpresto at all, or keep the Highway kernels in a separate small
  library linked only by the executables that use them.
- Whether `--sigma measured` should warn, or switch automatically, when the σ sanity check
  fires (Julia warns and does not switch).
- A cheaper alternative to rule out first: Julia 1.12's `juliac --trim` for a self-contained
  executable. Likely fragile with FFTW.jl and ArgParse, and it still does not slot into
  PRESTO's meson build or the feedstock, so it is not the plan.

## Other open items

### Backport coherentsearch additions to CoherentSearch.jl

Four things the C port gained (November 2026) belong in the Julia reference
(`~/git/CoherentSearch.jl`) too:
- **Periodic-signal (comb) removal** (`cs_remove_periodic()` in `src/cs_search.c`): a strong
  sinusoid lands on one high harmonic of many slow trials, giving combs of identical peaks.
  Drop a candidate whose strongest harmonic is the 7th or higher and more than 8x the largest
  of harmonics 1-4 (`CS_PERIODIC_*` in `include/coherentsearch.h`; `-noperiodic` disables).
- **The FFT normalisation check** (`cs_check_normalisation()`): √(median power / ln 2) over
  bins 10000-19999 (the top half of a shorter FFT), straight from the file before the search,
  warning if off by more than 10%. It replaces the chunk-based `_sigma_sanity_check` that
  Julia still uses. That check samples the lowest chunk (about 17 fundamental bins at 0.1 Hz)
  and takes the worst of ~18 comparisons, so it warned falsely on well-normalised Ter5 data,
  and on white noise where only a few harmonics fit.
- **The `-deep` pass:** a search for slow, very narrow pulses. It defaults to `-nharms 120
  -maxdecim 12 -hifreq 5` (`CS_DEEP_*` in `include/cs_fftplan.h`; the paper's
  `coherent_tier`), but explicitly given options win. It writes `<input>_deep.cohout`, so it
  never overwrites the default search's output.
- **Sifting a deep pass against the default `.cohout`** (`cs_remove_harmonics_of()`): run
  after a default search, `-deep` reads `<input>.cohout` and drops its candidates that are
  harmonically related to stronger ones there, such as subharmonics of a bright pulsar above
  the deep band (Ter5A/2 on the Ter5 data). `-noharmremove` skips it, and a note is printed
  if no default `.cohout` exists.
- **Sifter improvements** made porting `sift_candidates.py` to `bin/coherent_sift.py`:
  reading the Known PSR column (a `known` class and flag, and `--no-known`), reading
  `*_deep.cohout` files and flagging deep-pass signals (`--no-deep`), and classing
  single-observation runs on their DM behaviour instead of calling everything `weak`.

### libpresto headers and soversion

`libpresto` installs as an unversioned `.so` with **no headers** (only `presto.pc`), so nothing
can link against it from the conda package, and the recipe declares no `run_exports`. To open
that up: install the headers under `include/presto/` (names like `mask.h`, `median.h`,
`vectors.h` are too generic for a shared `include/`), give the library a `soversion:`, and only
then add `run_exports` in the feedstock. Relevant to the C++/Highway question above.

### Make the pulsar-database code more robust

`src/database.c` (`read_database`, `comp_psr_to_cand`, `comp_bin_to_cand`,
`comp_rawbin_to_cand`) is very handy for searches of fields with many known pulsars, such as
globular clusters, but it is old code. It keeps the whole catalogue in a static array of
fixed size `NP`, matches by pointer-keyed statics (one such cache was removed for
coherentsearch), writes into caller buffers of implied size, and exits if `pulsars.cat`
is missing. Worth a cleanup: a loaded-catalogue struct instead of statics, explicit
buffer lengths, an error return instead of `exit`, and perhaps an option to drop the
loose sidelobe matches.

### OpenBLAS threads spin in every PRESTO tool

libpresto links GSL, whose conda-forge build needs `libcblas.so.3`, which is OpenBLAS
(pthreads) by default. OpenBLAS starts one thread per core at load time, even though PRESTO
never calls BLAS, and each busy-waits about 0.1 s before sleeping. Measured 2026-11-10 on an
8-core machine: `coherentsearch` used ~0.7 s of extra CPU on both a 0.7 s and an 18 s run,
so it is a fixed cost per process, not continuous. It matters for many short runs, for one
job per core over many DMs (the reason `-ncpus` defaults to 1), and on many-core nodes,
where each process starts that many threads. `OPENBLAS_NUM_THREADS=1` removes it. The pixi
environment now sets it on activation, and README.md recommends it.

Still open, for the conda-forge package:
- An `etc/conda/activate.d` (and `deactivate.d`) script in the recipe could set it, but it
  would apply to the whole environment, making NumPy/SciPy linear algebra single-threaded
  for users who rely on it, and it does nothing for binaries run without activation.
- Better: stop loading OpenBLAS at all. libgsl itself needs a CBLAS, but conda-forge's gsl
  also ships its own `libgslcblas`; linking libpresto against that ahead of `libcblas` may
  not be enough, since libgsl's own NEEDED entry still pulls in `libcblas.so.3`. Check
  whether anything in libpresto needs a fast BLAS, and whether the GSL routines it uses
  touch BLAS at all.

### FFTFIT uncertainty calibration

Both FFTFIT methods under-cover their 1-sigma errors by ~5–15%, worst for broad profiles (e.g.
kappa = 1: aarchiba 0.656, classic 0.625 against 0.683). It is inherited from the Fortran, not a
regression. `aarchiba` is the better starting point: resolve Anne Archibald's
`std`-normalization `FIXME`, look at linearization bias for broad profiles, or validate an
empirical inflation factor against the coverage tests.

### Small issues found while writing the synthetic tests

- `sifting.candlist_from_candfile` requires `_DM<float>` in the filename (`_DM50.00`, not
  `_DM50`) and raises otherwise.
- Raw-data readers work in fixed blocks (2400 spectra for SIGPROC), and prepdata drops the
  final partial block, padding it with the mean and recording that in the `.inf`. For a file
  of only one block, the padding value is 0 instead of the mean. Harmless, but odd.

### `prepfold_multi` generated command-line files

`src/prepfold_multi_cmd.[ch]` need **manual** edits after `pyclig` generates them from
`clig/prepfold_multi_cmd.cli`. Fix the generator or the `.cli` so it regenerates cleanly like
every other tool. Until then, do not blindly overwrite the hand edits.

### Python docstrings, type hints and ruff

Continue adding numpy-style docstrings and type hints across `python/presto`, especially
`psr_utils.py`, `prepfold.py`, `psrfits.py`, `spectra.py`, `filterbank.py` and `sifting.py`, and
ruff-format the files touched.

### Remaining replaceable code

- `python/presto/parfile.py` optionally imports the old `slalib` package (`sla_ecleq`/
  `sla_eqecl`/`sla_eqgal`). Move these to astropy (already required), minding TEMPO's
  ecliptic-obliquity conventions.
- Parts of `misc_utils.c` still originate in Numerical Recipes; replace them with GSL/ERFA
  where possible.
- Further dead-code removal stays item-by-item, with each removal approved first.

## Postponed, potentially indefinitely

### Drop support for retired pulsar backends

Remove the code for old backends (Spigot, BCPM, WAPP, and others), believed to affect only the
`readfile` code path. Postponed: the code compiles without warnings and is unlikely to be
touched again, and old files from those machines still turn up, where `readfile` is very useful.
