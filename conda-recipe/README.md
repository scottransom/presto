# conda-forge recipe (reference copy)

PRESTO is on conda-forge as **[`presto-pulsar`](https://anaconda.org/conda-forge/presto-pulsar)**.
This directory holds a *reference copy* of the recipe. It is **not** used by the normal
meson/pixi build, and it is **not** what conda-forge builds from &mdash; it just documents how
PRESTO is packaged, and gives somewhere to prototype packaging changes.

- `recipe.yaml` / `build.sh` — the recipe. It builds `libpresto` + the C tools with meson,
  then the Python package with pip, taking ERFA from conda-forge's `liberfa` package.

The authoritative copy lives in the feedstock,
[`conda-forge/presto-pulsar-feedstock`](https://github.com/conda-forge/presto-pulsar-feedstock),
under `recipe/`. **Keep the two in sync** whenever the packaging itself changes (dependencies,
build steps, tests). The `version:`/`sha256:` here will lag the feedstock between releases;
that is fine and expected. To check for drift:

    diff <(tail -n +14 conda-recipe/recipe.yaml) \
         <(tail -n +5 /path/to/presto-pulsar-feedstock/recipe/recipe.yaml)
    diff conda-recipe/build.sh /path/to/presto-pulsar-feedstock/recipe/build.sh

(the `tail` offsets just skip the differing header comments).

## Recipe format

The recipe uses the **v1 (rattler-build) recipe format** defined by
[CEP-13](https://rattler.build/latest/reference/recipe_file/) — `recipe.yaml`, not the older
conda-build `meta.yaml`. Differences to remember when editing: `${{ ... }}` jinja, a `context:`
block instead of `{% set %}`, `if: <selector>` / `then:` lists instead of `# [selector]`
comments, `target_directory:` instead of `folder:`, and a `tests:` list instead of `test:`.
`conda-recipe-manager convert meta.yaml` does most of a v0 → v1 conversion automatically.

## Releasing a new version

See [`../RELEASE.md`](../RELEASE.md). In short: tag the release on GitHub, then wait a day for
conda-forge's **regro-cf-autotick-bot** to open the version-bump PR on the feedstock and merge
it once CI is green (or do the bump by hand from a fork of the feedstock if you are impatient
or the bot cannot figure it out).
