# conda-forge recipe (reference copy)

This directory holds a reference copy of PRESTO's conda-forge recipe. It is **not** used by
the normal meson/pixi build; it documents how PRESTO is packaged for conda-forge.

- `recipe.yaml` / `build.sh` — the recipe. It builds `libpresto` + the C tools with meson,
  then the Python package with pip, taking ERFA from conda-forge's `liberfa` package.

The recipe uses the **v1 (rattler-build) recipe format** defined by
[CEP-13](https://rattler.build/latest/reference/recipe_file/) — `recipe.yaml`, not the older
conda-build `meta.yaml`. Differences to remember when editing: `${{ ... }}` jinja, a `context:`
block instead of `{% set %}`, `if: <selector>` / `then:` lists instead of `# [selector]`
comments, `target_directory:` instead of `folder:`, and a `tests:` list instead of `test:`.
`conda-recipe-manager convert meta.yaml` does most of a v0 → v1 conversion automatically.

## First submission

Copy these two files into a fork of
[`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes) at
`recipes/presto-pulsar/`, fill in the `sha256` of the v6.0.0 GitHub tarball
(`curl -sL https://github.com/scottransom/presto/archive/refs/tags/v6.0.0.tar.gz | sha256sum`),
test locally with `python build-locally.py`, and open a PR. See `../RELEASE.md` for the full
checklist.

## After the feedstock exists

Once `presto-pulsar-feedstock` is created, its `recipe/recipe.yaml` becomes authoritative and
conda-forge's autotick bot proposes version bumps automatically. Keep this copy roughly in
sync when the packaging itself changes (deps, build steps).
