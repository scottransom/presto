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

PRESTO is on conda-forge as **`presto-pulsar`**, built by the feedstock
[`conda-forge/presto-pulsar-feedstock`](https://github.com/conda-forge/presto-pulsar-feedstock).
The feedstock's `recipe/recipe.yaml` is the authoritative recipe; `conda-recipe/` in this repo
is only a reference copy (see [`conda-recipe/README.md`](conda-recipe/README.md)).

### The easy way: let the bot do it

Within a day or so of the GitHub release, conda-forge's **regro-cf-autotick-bot** notices the
new tag and opens a version-bump PR on the feedstock ("presto-pulsar v X.Y.Z"). All you have
to do is:

1. Check the diff: `version:` and `sha256:` match the new tag (compare against the hash from
   step 3), and `build: number:` is back to `0`.
2. Wait for the Azure CI to go green for linux-64, osx-64, and osx-arm64.
3. Merge it. The new packages hit the `conda-forge` channel roughly half an hour later.

If the recipe needs anything beyond the version/hash (a new dependency, a changed build step,
a new tool in `package_contents`), just push those extra commits onto the bot's PR branch
before merging. If you change dependencies or anything in `conda-forge.yml`, ask for a
re-render by commenting on the PR:

    @conda-forge-admin, please rerender

### The manual way: bump it yourself

Worth doing if the bot is slow or gets confused, or if the release needs recipe changes
anyway. Working from the fork at `~/git/presto-pulsar-feedstock` (add the upstream remote once
with `git remote add upstream https://github.com/conda-forge/presto-pulsar-feedstock.git`):

```sh
cd ~/git/presto-pulsar-feedstock
git fetch upstream
git checkout -b vX.Y.Z upstream/main       # feedstock's default branch is 'main'
# In recipe/recipe.yaml, edit:
#   context: version:  ->  "X.Y.Z"
#   source: sha256:    ->  the hash from step 3
#   build: number:     ->  0   (only matters if the previous release bumped it)
git commit -am "presto-pulsar vX.Y.Z"
git push -u origin vX.Y.Z
gh pr create --repo conda-forge/presto-pulsar-feedstock \
    --base main --head scottransom:vX.Y.Z \
    --title "presto-pulsar vX.Y.Z" --body "Update to PRESTO vX.Y.Z"
```

Then the same as above: watch CI, `@conda-forge-admin, please rerender` if the dependencies or
`conda-forge.yml` changed, and merge when green. To try a build before pushing, run
`python build-locally.py` in the feedstock (it uses docker for the linux targets).

### Keep the in-repo copy in sync

If the *packaging* changed (not just version/hash), copy the changes back into
`conda-recipe/recipe.yaml` and `conda-recipe/build.sh` in this repo as part of the release
commit, and check for drift with the `diff` commands in `conda-recipe/README.md`.

### First-time bootstrapping (historical)

The initial submission went through
[`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes) under
`recipes/presto-pulsar/`, in the v1/rattler-build recipe format, and was merged as PR #34070
in August 2026. That path is no longer used now that the feedstock exists.

## Notes

- The recipe takes ERFA from conda-forge's `liberfa`, so it no longer has to track the
  version pinned by `subprojects/erfa.wrap` (that wrap still serves source builds without a
  system ERFA).
- Version numbers live in exactly the three files `determine_version.py` rewrites &mdash;
  don't edit them by hand.
