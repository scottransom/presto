#!/bin/sh

# Fetch the real-data fixture for the prepfold suites into this directory
# (if it is not already here) and verify its SHA-256 checksum.
#
# The data are too big for the git repo (173 MB), so they live as assets on a
# dedicated GitHub release of scottransom/presto.  The tag is versioned: if the
# fixture ever changes, upload it under a new tag (test-data-v2, ...) and
# update TAG and the checksum below, so old commits keep fetching the data
# their reference folds were made from.  The original NRAO copy is a fallback.
#
# Exits non-zero if the file cannot be fetched or its checksum is wrong.

TAG=test-data-v1
FILE=Ter5_080912_short2bits.fits
SHA256=4a20519afbd1dcbce5c906ef11ffc92b2692ab513efa17d1b43df07b121d2581
URLS="https://github.com/scottransom/presto/releases/download/${TAG}/${FILE}
https://www.cv.nrao.edu/~sransom/${FILE}"

cd "$(dirname "$0")" || exit 2

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

if [ -f "$FILE" ] && [ "$(sha256_of "$FILE")" = "$SHA256" ]; then
    exit 0
fi

for url in $URLS; do
    echo "Downloading $url"
    if curl -fsSL --retry 3 -o "$FILE.part" "$url" && \
       [ "$(sha256_of "$FILE.part")" = "$SHA256" ]; then
        mv "$FILE.part" "$FILE"
        exit 0
    fi
    echo "  failed (download error or checksum mismatch)"
done
rm -f "$FILE.part"
echo "ERROR: could not fetch $FILE with SHA-256 $SHA256"
exit 1
