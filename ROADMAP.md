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

## Next major goal: a coherent harmonic search in PRESTO

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
   data. That bounds what any SIMD work could buy the existing tools. Expect modest gains,
   since FFTW dominates most of them.
1. **conda-forge microarchitecture variants (no code changes).** Add the
   `x86_64-microarch-level` metapackage to the feedstock's build requirements to produce v1 and
   v3 builds; the solver chooses one via the `__archspec` virtual package. That gives every
   auto-vectorized loop in PRESTO AVX2/FMA. Costs: double the x86 build configs; CI cannot
   guarantee level 4, so no AVX-512 — arguably a feature, given the AVX-512 licence-downclock
   history below.
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

**Plan:** do (0) now, (1) on the feedstock independently of the port, and (2) for the new
search's hot kernels: the grouped interpolation matvec and the transposed boxcar
gate/prefix-sum scan. Then adopt Highway kernels in `accelsearch` *only where profiling says
the non-FFTW loops matter* (see `fourier_interpolation_tradeoffs.md` §5 in CoherentSearch.jl).
The one plausible `accelsearch` target is the harmonic sum in `add_ffdotpows`, whose
`inpows[rinds[kk]]` gather is a near-arithmetic index sequence that could become contiguous
loads.

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

- CI (above), the `-march=x86-64-v3` measurement, and the feedstock microarch change go to
  `master` (or the feedstock) as usual — none of them needs the port.
- The port lives on a branch, e.g. `coherent-search`, rebased on `master` periodically so it
  runs under CI from day one.
- Merge into `master` via a PR once phase 3 below passes its pins on all CI platforms; a
  release with the new tool follows.

### Phases

1. **k = 1 only, single-threaded, scalar C:** interpolation, base transform, analytic σ,
   gate and exact rescore. Pinned against Julia per-trial metrics.
2. **The full search:** decimation ladder, pruned widths, measured σ, candidate collapse,
   output file, clig `.cli` and man page.
3. **Throughput:** OpenMP, multi-file invocations, wisdom; then PGPLOT plots and the replot
   tool. *Merge point.*
4. **SIMD:** Highway kernels with dispatch, constant re-tuning, disassembly audits; meson and
   feedstock changes for `libhwy` and C++.
5. **Afterwards, optional:** profile `accelsearch` and move any loops worth it onto Highway.

### Open questions

- The tool's name (`coherentsearch`, `cohsearch`, …) and its candidate-file format (keep
  `.cohout`, or something `ACCEL_sift`/`prepfold`-friendly?).
- Whether to accept C++ in libpresto at all, or keep the Highway kernels in a separate small
  library linked only by the executables that use them.
- Whether `--sigma measured` should warn, or switch automatically, when the σ sanity check
  fires (Julia warns and does not switch).
- A cheaper alternative to rule out first: Julia 1.12's `juliac --trim` for a self-contained
  executable. Likely fragile with FFTW.jl and ArgParse, and it still does not slot into
  PRESTO's meson build or the feedstock, so it is not the plan.

## Other open items

### libpresto headers and soversion

`libpresto` installs as an unversioned `.so` with **no headers** (only `presto.pc`), so nothing
can link against it from the conda package, and the recipe declares no `run_exports`. To open
that up: install the headers under `include/presto/` (names like `mask.h`, `median.h`,
`vectors.h` are too generic for a shared `include/`), give the library a `soversion:`, and only
then add `run_exports` in the feedstock. Relevant to the C++/Highway question above.

### Raise the minimum Python version to 3.12 (decided, next release)

`python/pyproject.toml` declares `requires-python = '>=3.9'`, and `INSTALL.md` still says
">=3.8" in one place. NumPy, SciPy, astropy and matplotlib follow
[SPEC 0](https://scientific-python.org/specs/spec-0000/), which retires a Python version three
years after release, so 3.11 (October 2022) is already past due and 3.12 puts PRESTO in step
with its dependencies. conda-forge will drop 3.10 at the next rerender regardless. For the next
release:
  * set `requires-python = '>=3.12'` in `python/pyproject.toml`,
  * fix the stale ">=3.8" line in `INSTALL.md` (and check for others),
  * add `python >=3.12` to `requirements/host` and `requirements/run` in the feedstock recipe —
    a bare `python` builds for whatever conda-forge's global pinning offers, so the recipe, not
    `requires-python`, is what narrows the matrix,
  * confirm nothing in the package still guards on older versions.

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
