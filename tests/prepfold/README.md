# prepfold / prepfold_multi test suite

This directory holds the regression and validation tests for `prepfold` and
`prepfold_multi`. All scripts are run from inside `tests/prepfold/` and most
expect the test data file `Ter5_080912_short2bits.fits` to be present (it is
downloaded automatically if missing). Prepend a fresh build's `src/` dir to
`PATH` before running anything so the scripts pick up the binaries you just
built, e.g.:

```sh
export PATH=/home/rlynch/sandboxes/presto/build/src:$PATH
cd tests/prepfold
```

## Test scripts

- **`prepfold_tests.sh`** -- the single-candidate `prepfold` regression gate.
  Runs 8 folds covering topocentric/barycentric time series, PRESTO subbands,
  raw PSRFITS (with and without an `rfifind` mask), and barycentric events,
  comparing each `.bestprof` against the reference copies in `goodfolds/` via
  `compare_bestprof.py`. It also runs the `reader_equiv_check` raw-block
  reader-equivalence harness twice (with and without a mask), and finishes
  with `check_cmd_layout.sh`.

- **`prepfold_multi_tests.sh`** -- the `prepfold_multi` equivalence gate.
  Proves the defining correctness property, "multi == many singles": folding
  N candidates (each with its own DM) in a single pass over the raw data must
  reproduce, field for field, what N independent `prepfold` runs produce.
  Checked for both supported inputs: raw PSRFITS (`multi_cands.txt`, 3
  candidates) and PRESTO subbands (`multi_cands_sub.txt`, 3 candidates). Also
  reports the wall-clock of the one multi run vs. the sum of the single runs,
  to make the "the raw data is read once" claim measured rather than just
  structural.

- **`check_cmd_layout.sh`** -- a source-level drift guard, not a data test.
  `prepfold_multi` reuses `prepfold`'s shared fold pipeline, so its parsed
  options must populate the *same* `struct s_Cmdline` memory layout that
  `prepfold` uses (`include/prepfold_multi_cmd.h` is `prepfold_cmd.h` with one
  `candfile` field appended after `full_cmd_line`). This script extracts the
  `s_Cmdline` body from both generated headers and asserts that
  `prepfold_cmd.h`'s body is a line-for-line prefix of
  `prepfold_multi_cmd.h`'s, catching any reorder, rename, type change, or
  removal of a shared field before it silently corrupts a shared-struct read.
  Pure header check; no build or test data required.

- **`compare_bestprof.py`** -- the shared tolerance-based `.bestprof`
  comparator used by every script above. Checks folded sample count exactly,
  averages/stddevs/chi^2/spin parameters within relative or error-based
  tolerances, best DM against the DM step (read from the `.pfd` if given via
  `--pfd`), and the profile shape to `1e-3`. See the docstring in the script
  for the full tolerance table.

- **`reader_equiv_check.c`** -- a small standalone C harness (not built by
  default; build with `meson compile -C build reader_equiv_check`) that
  proves the factored read-once/dedisperse-many raw-block reader produces
  bit-identical subbands to the original `read_subbands()`, with and without
  an `rfifind` mask applied.

### prepfold_multi validation scripts

These extend the equivalence gate above with determinism, memory-safety, and
performance evidence for `prepfold_multi`'s OpenMP candidate-loop
parallelism (`-ncpus`):

- **`determinism_check.sh`** -- folds the RAWDATA and insubs candidate sets at
  `-ncpus` 1, 4, and 8 into distinct output roots, then md5-compares each
  candidate's `.bestprof` across thread counts. Bit-identical across thread
  counts means the output is both deterministic and unaffected by `-ncpus`.
  Also tolerance-compares the `-ncpus 4` output directly against the
  single-`prepfold` references left behind by `prepfold_multi_tests.sh`.

- **`valgrind_check.sh`** -- memory and race checks on a short data slice
  (`-start`/`-end`, small `-npart`) so the instrumented runs finish quickly
  while still exercising setup, the parallel fold pass, the parallel
  optimize/output stage, and cleanup. Runs `memcheck` at `-ncpus 1` (serial)
  and `-ncpus 2` (parallel; memcheck serializes threads, so this still
  catches out-of-bounds/uninitialized access deterministically), then
  `helgrind` at `-ncpus 2` for data races.

- **`perf_check.sh`** -- two measurements: (A) read-once I/O win, one 3-candidate
  multi pass vs. three separate single `prepfold` runs over the same raw
  file; (B) OpenMP candidate-loop scaling, one 8-candidate multi run
  (`multi_cands_perf.txt`) at `-ncpus` 1, 2, and 4, plus an md5 check that the
  `-ncpus 1` and `-ncpus 4` results stay bit-identical.

## Candidate files

`multi_cands.txt`, `multi_cands_sub.txt`, and `multi_cands_perf.txt` are
`-candfile` inputs for `prepfold_multi`: plain-text tables, one candidate per
row, each with its own DM and period (or frequency) plus optional
derivatives. See the `prepfold_multi(1)` "CANDIDATE FILE FORMAT" section for
the full format spec; `multi_cands.txt` is reproduced there as the worked
example. `multi_cands_sub.txt` covers the PRESTO-subband path (residual DMs
around the subbands' DM), and `multi_cands_perf.txt` is an 8-candidate file
sized for the OpenMP scaling measurement in `perf_check.sh`.

## What to expect

Based on the most recent validation run:

- `prepfold_tests.sh`: all pass -- 8 folds, `reader_equiv_check` x2 (no mask /
  with mask), and the command-layout check.
- `prepfold_multi_tests.sh`: RAWDATA 3/3 and insubs 3/3 pass -- every
  `prepfold_multi` candidate matches its single-`prepfold` counterpart field
  for field.
- `determinism_check.sh`: `.bestprof` files are bit-identical across `-ncpus`
  1, 4, and 8 for every candidate, in both RAWDATA and insubs.
- `valgrind_check.sh`: memcheck reports 0 errors / 0 leaks for both the
  serial and parallel runs. helgrind reports a number of "possible data
  race" lines, but they are all libgomp-internal false positives rooted in
  PRESTO's pre-existing OpenMP region in `get_PSRFITS_subint`
  (`psrfits.c:924`), triggered during the serial raw read -- none of them
  reference `prepfold_multi`'s own parallel regions. The authoritative
  thread-safety evidence is the bit-identical determinism above combined
  with the clean memcheck run, not the helgrind log.
- `perf_check.sh` (illustrative numbers from one machine -- expect different
  absolute numbers on other hardware, not exact reproductions): read-once
  speedup around 1.12x for 3 candidates; the 8-candidate OpenMP run went from
  about 189s at `-ncpus 1` to about 93.6s at `-ncpus 4` (about 2.0x at 4
  cores), Amdahl-bounded by the serial per-block raw read that every
  candidate still shares.
