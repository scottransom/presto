# Making a PRESTO release

This is the repeatable checklist for cutting a tagged PRESTO release and keeping the
conda-forge package in sync. Replace `X.Y.Z` with the new version throughout.

## 1. Prepare the release commit (on `master`)

1. Make sure `master` is clean and up to date and the build/tests pass
   (`pixi run build && pixi run test`).
2. Stamp the version into the three files that must agree
   (`meson.build`, `python/pyproject.toml`, `python/presto/__init__.py`):

   ```sh
   python determine_version.py --set X.Y.Z --write
   ```

   `determine_version.py` with no `--set` derives a development version from
   `git describe` (e.g. `X.Y.Z.dev12`); `--set` is only for real releases.
3. Update `CHANGELOG.md`: rename the top `## Development (unreleased, since ...)` heading
   to `## Version X.Y.Z:` (keeping its accumulated bullets, with a short lead summary), and
   open a fresh empty `## Development (unreleased, since vX.Y.Z):` section above it.
4. Add a `## Version X.Y.Z:` highlights block to the top of `README.md`.
5. Commit: `git commit -am "Release PRESTO X.Y.Z"`.

## 2. Tag, push, and publish on GitHub

```sh
git tag -a vX.Y.Z -m "PRESTO vX.Y.Z"
git push origin master
git push origin vX.Y.Z
gh release create vX.Y.Z --title "PRESTO vX.Y.Z" --notes-file <changelog excerpt>
```

After tagging, `python determine_version.py` (no args) should print exactly `X.Y.Z`.

## 3. Get the source hash for packaging

conda-forge hashes GitHub's auto-generated tag tarball:

```sh
curl -sL https://github.com/scottransom/presto/archive/refs/tags/vX.Y.Z.tar.gz | sha256sum
```

## 4. Update conda-forge

- **First release (bootstrapping the feedstock):** submit a recipe to
  [`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes) under
  `recipes/presto-pulsar/` (`recipe.yaml` + `build.sh`, in the v1/rattler-build recipe
  format). A copy of the recipe is kept in this repo under `conda-recipe/` for reference.
  Notes:
  - ERFA has no conda-forge feedstock, so the recipe vendors the ERFA release tarball as a
    second `source:` unpacked into `subprojects/erfa-2.0.1/` (meson builds it with no
    network); its BSD-3 `LICENSE` is bundled alongside PRESTO's `COPYING`.
  - First release targets `linux-64` + `osx-64`; `osx-arm64` is skipped until `tempo2` is
    packaged for it (tempo2 is only needed at runtime for `prepfold -timing` polycos).
  - Test the recipe offline with staged-recipes' `python build-locally.py` before opening
    the PR.
- **Subsequent releases:** once `presto-pulsar-feedstock` exists, conda-forge's
  **regro-cf-autotick-bot** usually opens the version-bump PR automatically within a day of
  the GitHub release. Just review/merge it (bump `build: number: 0`, re-pin deps if needed).
  To do it manually, edit `recipe/recipe.yaml`'s `version` + `sha256` in the feedstock and
  open a PR.

## Notes

- The ERFA tarball version tracked by the recipe is whatever `subprojects/erfa.wrap`
  pins; keep the two in sync when bumping ERFA.
- Version numbers live in exactly the three files `determine_version.py` rewrites &mdash;
  don't edit them by hand.
