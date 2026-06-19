#!/usr/bin/env python3
"""Tolerance-based .bestprof comparator for prepfold validation.

Usage:
    compare_bestprof.py <new>.bestprof <ref>.bestprof [--pfd FILE]

Exit 0 (PASS) if all checked fields agree within tolerance.
Exit 1 (FAIL) if any field exceeds its tolerance.

Field tolerances:
    Data Folded              exact integer
    Data Avg / Data StdDev   relative 1e-5
    Profile Avg / StdDev     relative 1e-5
    Reduced chi-sqr          relative 1e-2
    P / P' / P''             within max(err_ref, err_new); zero derivs skipped
    Best DM                  abs <= one DM step (from --pfd) or 0.01 fallback
    Profile shape            max |delta_normalized| <= 1e-3

Ignored: Input file, Candidate, Telescope, Epoch_*, T_sample, Prob(Noise).
"""

import sys
import os
import argparse
try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required. Install it in your PRESTO Python environment: pip install numpy")


# Add parent python/ dir so presto can be imported when --pfd is used.
_here = os.path.dirname(os.path.abspath(__file__))
_presto_python = os.path.join(_here, '..', 'python')
if os.path.isdir(_presto_python) and _presto_python not in sys.path:
    sys.path.insert(0, _presto_python)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_bestprof(filename):
    """Return (header_dict, profile_array).

    header_dict maps normalized key (no unit suffix) -> raw value string.
    profile_array is a float64 array of the profile bins.
    """
    header = {}
    profile = []
    in_profile = False

    with open(filename) as f:
        for line in f:
            line = line.rstrip('\n')
            if in_profile:
                parts = line.split()
                if parts:
                    profile.append(float(parts[-1]))
            elif line.startswith('#'):
                if line.startswith('######'):
                    in_profile = True
                    continue
                body = line[1:].strip()
                if '=' in body:
                    key_part, val_part = body.split('=', 1)
                elif '<' in body:
                    # e.g. "Prob(Noise)      <  2.2e-18 (~8.7 sigma)"
                    key_part, val_part = body.split('<', 1)
                    val_part = '<' + val_part
                else:
                    continue
                key = key_part.strip()
                # Strip unit suffix "(ms)", "(s/s)", etc.
                if '(' in key:
                    key = key[:key.index('(')].strip()
                header[key] = val_part.strip()

    return header, np.array(profile, dtype=np.float64)


def parse_val_err(raw):
    """Parse 'value +/- err' string.  Returns (float, float) or (None, None) for N/A."""
    raw = raw.strip()
    if raw in ('N/A', ''):
        return None, None
    if '+/-' in raw:
        v, e = raw.split('+/-', 1)
        return float(v.strip()), float(e.strip())
    return float(raw), None


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def rel_dev(ref, new):
    if ref == 0.0:
        return abs(new - ref)          # absolute when ref is zero
    return abs(new - ref) / abs(ref)


def fmt_row(field, ref_v, new_v, deviation, status):
    return f"  {field:<35s}  ref={str(ref_v):<22s}  new={str(new_v):<22s}  dev={str(deviation):<20s}  {status}"


# Required header keys for comparison
_REQUIRED_KEYS = ('Data Folded', 'Data Avg', 'Data StdDev',
                  'Profile Avg', 'Profile StdDev', 'Reduced chi-sqr')


def _load_file(path, label, rows, failures):
    """Open and parse a .bestprof file; record FAIL rows on error.

    Returns (hdr, prof) on success, or (None, None) on failure.
    """
    if not os.path.exists(path):
        rows.append(f'  ERROR: {label} file not found: {path}')
        failures.append(f'{label}_missing')
        return None, None
    try:
        hdr, prof = parse_bestprof(path)
    except Exception as exc:
        rows.append(f'  ERROR: could not parse {label} file ({path}): {exc}')
        failures.append(f'{label}_parse_error')
        return None, None
    missing = [k for k in _REQUIRED_KEYS if k not in hdr]
    if missing:
        rows.append(f'  ERROR: {label} file is missing required header keys: {missing}')
        rows.append(f'  Keys found: {sorted(hdr)}')
        rows.append(f'  File: {path}')
        failures.append(f'{label}_missing_keys')
        return None, None
    return hdr, prof


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def compare(new_file, ref_file, pfd_file=None):
    """Compare new vs ref.  Returns (passed: bool, report_lines: list[str])."""
    rows = []
    failures = []

    ref_hdr, ref_prof = _load_file(ref_file, 'ref',  rows, failures)
    new_hdr, new_prof = _load_file(new_file, 'new', rows, failures)

    if ref_hdr is None or new_hdr is None:
        return False, rows

    def check(name, ok, ref_v, new_v, dev_str):
        status = 'OK' if ok else 'FAIL'
        rows.append(fmt_row(name, ref_v, new_v, dev_str, status))
        if not ok:
            failures.append(name)

    # --- Data Folded: exact integer ---
    ref_n = int(float(ref_hdr['Data Folded']))
    new_n = int(float(new_hdr['Data Folded']))
    check('Data Folded', ref_n == new_n, ref_n, new_n,
          '0' if ref_n == new_n else str(abs(new_n - ref_n)))

    # --- Stats: relative 1e-5 ---
    for field in ('Data Avg', 'Data StdDev', 'Profile Avg', 'Profile StdDev'):
        rv = float(ref_hdr[field])
        nv = float(new_hdr[field])
        dev = rel_dev(rv, nv)
        check(field, dev <= 1e-5, rv, nv, f'{dev:.3e}')

    # --- Reduced chi-sqr: relative 1e-2 ---
    rv = float(ref_hdr['Reduced chi-sqr'])
    nv = float(new_hdr['Reduced chi-sqr'])
    dev = rel_dev(rv, nv)
    check('Reduced chi-sqr', dev <= 1e-2, rv, nv, f'{dev:.3e}')

    # --- Period / derivative comparisons ---
    period_compared = False
    deriv_keys = [
        (0, 'P'),
        (1, "P'"),
        (2, "P''"),
    ]
    for frame in ('topo', 'bary'):
        for order, prefix in deriv_keys:
            key = f'{prefix}_{frame}'      # e.g. "P_bary", "P'_topo"
            ref_raw = ref_hdr.get(key, 'N/A')
            new_raw = new_hdr.get(key, 'N/A')

            ref_v, ref_e = parse_val_err(ref_raw)
            new_v, new_e = parse_val_err(new_raw)

            if ref_v is None:
                continue                   # N/A in reference — skip

            if order > 0 and ref_v == 0.0:
                continue                   # zero derivative — skip per spec

            if new_v is None:
                check(key, False, ref_raw, 'N/A', 'new missing value')
                continue

            err_tol = max(ref_e or 0.0, new_e or 0.0)
            if err_tol > 0.0:
                diff = abs(new_v - ref_v)
                ok = diff <= err_tol
                dev_str = f'{diff:.3e} (tol {err_tol:.3e})'
            else:
                dev = rel_dev(ref_v, new_v)
                ok = dev <= 1e-5
                dev_str = f'{dev:.3e} (rel)'

            check(key, ok, ref_raw, new_raw, dev_str)
            if order == 0:
                period_compared = True

    if not period_compared:
        rows.append('  ERROR: no period field present in reference — treating as FAIL')
        failures.append('period_present')

    # --- Best DM (only when reference has it) ---
    ref_dm_raw = ref_hdr.get('Best DM')
    if ref_dm_raw is not None:
        ref_dm = float(ref_dm_raw)
        new_dm_raw = new_hdr.get('Best DM')
        if new_dm_raw is None:
            check('Best DM', False, ref_dm, 'MISSING', 'absent in new file')
        else:
            new_dm = float(new_dm_raw)
            dm_tol = _dm_tolerance(pfd_file, rows)
            diff = abs(new_dm - ref_dm)
            check('Best DM', diff <= dm_tol, ref_dm, new_dm,
                  f'{diff:.4f} (tol {dm_tol:.4f})')

    # --- Profile shape ---
    if len(ref_prof) > 0 and len(new_prof) > 0:
        ref_avg = float(ref_hdr['Profile Avg'])
        ref_std = float(ref_hdr['Profile StdDev'])
        new_avg = float(new_hdr['Profile Avg'])
        new_std = float(new_hdr['Profile StdDev'])
        if ref_std > 0.0 and new_std > 0.0:
            ref_norm = (ref_prof - ref_avg) / ref_std
            new_norm = (new_prof - new_avg) / new_std
            max_diff = float(np.max(np.abs(ref_norm - new_norm)))
            tol = 1e-3
            check('Profile shape (max|Δnorm|)', max_diff <= tol,
                  '0', f'{max_diff:.3e}', f'{max_diff:.3e} (tol {tol:.0e})')

    passed = len(failures) == 0
    return passed, rows


def _dm_tolerance(pfd_file, rows):
    """Return abs DM tolerance: one DM step from pfd, or 0.01 fallback."""
    default = 0.01
    if pfd_file is None:
        return default
    try:
        from presto.prepfold import pfd
        p = pfd(pfd_file)
        if p.numdms > 1:
            dms = np.atleast_1d(p.dms)
            if len(dms) > 1:
                return float(abs(np.diff(dms).mean()))
    except Exception as exc:
        rows.append(f'  [DM] could not load pfd ({exc}); using tol={default}')
    return default


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description='Compare .bestprof files with numerical tolerances.')
    ap.add_argument('new', help='New .bestprof file (produced by current build)')
    ap.add_argument('ref', help='Reference .bestprof file')
    ap.add_argument('--pfd', metavar='FILE',
                    help='Path to new .pfd binary (enables precise DM-step tolerance)')
    args = ap.parse_args()

    print(f'Comparing : {args.new}')
    print(f'Reference : {args.ref}')
    print()

    passed, rows = compare(args.new, args.ref, pfd_file=args.pfd)

    for r in rows:
        print(r)

    print()
    if passed:
        print('PASS')
        sys.exit(0)
    else:
        print('FAIL')
        sys.exit(1)


if __name__ == '__main__':
    main()
