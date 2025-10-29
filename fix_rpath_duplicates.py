#!/usr/bin/env python3
"""
Post-install script to fix duplicate RPATH entries in installed binaries.
This works around a bug where gcc15 and gfortran both add /opt/local/lib/libgcc.
"""

import subprocess
import sys
import os
from pathlib import Path

def get_rpaths(binary):
    """Get list of RPATHs from a Mach-O binary."""
    try:
        result = subprocess.run(
            ['otool', '-l', str(binary)],
            capture_output=True,
            text=True,
            check=True
        )

        rpaths = []
        lines = result.stdout.split('\n')
        for i, line in enumerate(lines):
            if 'cmd LC_RPATH' in line:
                # The path is 2 lines down
                if i + 2 < len(lines):
                    path_line = lines[i + 2]
                    if 'path' in path_line:
                        # Extract path (format: "         path /some/path (offset 12)")
                        path = path_line.split('path')[1].split('(offset')[0].strip()
                        rpaths.append(path)
        return rpaths
    except subprocess.CalledProcessError:
        return []

def remove_duplicate_rpaths(binary):
    """Remove duplicate RPATH entries from a binary."""
    rpaths = get_rpaths(binary)
    if not rpaths:
        return False

    # Find duplicates
    seen = set()
    duplicates = []
    for rpath in rpaths:
        if rpath in seen:
            duplicates.append(rpath)
        else:
            seen.add(rpath)

    if not duplicates:
        return False

    print(f"Fixing {binary.name}: removing {len(duplicates)} duplicate RPATH(s)")

    # Remove all duplicates (install_name_tool removes first occurrence)
    for dup_rpath in duplicates:
        try:
            subprocess.run(
                ['install_name_tool', '-delete_rpath', dup_rpath, str(binary)],
                check=True,
                capture_output=True
            )
        except subprocess.CalledProcessError as e:
            # May fail if already removed, that's ok
            pass

    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: fix_rpath_duplicates.py <install_dir>")
        sys.exit(1)

    install_dir = Path(sys.argv[1])
    bin_dir = install_dir / 'bin'

    if not bin_dir.exists():
        print(f"Warning: {bin_dir} does not exist")
        return

    fixed_count = 0
    for binary in bin_dir.iterdir():
        if binary.is_file() and os.access(binary, os.X_OK):
            if remove_duplicate_rpaths(binary):
                fixed_count += 1

    if fixed_count > 0:
        print(f"Fixed {fixed_count} binaries")
    else:
        print("No duplicate RPATHs found")

if __name__ == '__main__':
    main()
