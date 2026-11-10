#!/bin/sh

# Fetch the fixtures for the real-data tests into this directory (if they are
# not already here) and verify their SHA-256 checksums:
#   Ter5_080912_short2bits.fits         GBT PSRFITS of Terzan 5 (173 MB), for
#                                       the search pipeline and prepfold suites
#   PM0063_034C1_DM445.0_red.{fft,inf}  a rednoised Parkes Multibeam spectrum
#                                       (34 MB) holding the slow pulsar
#                                       J1020-6026, for coherentsearch
#
# The data are too big for the git repo, so they live as assets on a dedicated
# GitHub release of scottransom/presto.  The tag is versioned: to change or add
# a fixture, upload the whole set under a new tag (test-data-v3, ...) and
# update TAG and the checksums below.  Never replace an asset in place, so old
# commits keep fetching the data their references were made from.  The NRAO
# copies are a fallback.
#
# Exits non-zero if any file cannot be fetched or its checksum is wrong.

TAG=test-data-v2
# file  sha256
FILES="Ter5_080912_short2bits.fits 4a20519afbd1dcbce5c906ef11ffc92b2692ab513efa17d1b43df07b121d2581
PM0063_034C1_DM445.0_red.fft dacfd7d998dc8f19b21673ea76ef7199c033f68c16223cff42fbae7789c06782
PM0063_034C1_DM445.0_red.inf e31459907039d5bcfe76e40d945f4bb8f5b1f8b7dd3d4610a1cdde29b1b49c5b"

cd "$(dirname "$0")" || exit 2

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

# fetch FILE SHA256: download FILE unless it is already here and intact.
fetch() {
    if [ -f "$1" ] && [ "$(sha256_of "$1")" = "$2" ]; then
        return 0
    fi
    for url in "https://github.com/scottransom/presto/releases/download/${TAG}/$1" \
               "https://www.cv.nrao.edu/~sransom/$1"; do
        echo "Downloading $url"
        if curl -fsSL --retry 3 -o "$1.part" "$url" && \
           [ "$(sha256_of "$1.part")" = "$2" ]; then
            mv "$1.part" "$1"
            return 0
        fi
        echo "  failed (download error or checksum mismatch)"
    done
    rm -f "$1.part"
    echo "ERROR: could not fetch $1 with SHA-256 $2"
    return 1
}

status=0
for f in $(echo "$FILES" | cut -d' ' -f1); do
    fetch "$f" "$(echo "$FILES" | grep "^$f " | cut -d' ' -f2)" || status=1
done
exit $status
