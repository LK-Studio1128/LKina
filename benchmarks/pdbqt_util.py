#!/usr/bin/env python3
"""Shared PDBQT I/O helpers for LKina benchmarks.

CRITICAL (v3 fix): PDBQT is column-strict.  The atom-type token must sit at
0-indexed columns 77-78 and coordinates at 30-37/38-45/46-53, exactly as
AutoDockTools / LKina / Vina 1.2.7 expect.  Earlier hand-rolled formatters
shifted x one column right (trailing space) which Vina 1.2.7 rejects with
"Coordinate ... is not valid", and put 2-char types at 78-79, truncating
NA/OA/Zn -> N/O/Z and silently removing metals from the grid.

This module builds lines by ABSOLUTE column positions:
  rec 0-5  serial 6-10  ' ' 11  name 12-15  resn 16-18  ' ' 19
  chain 20  resi 21-24  '     ' 25-29  x 30-37  y 38-45  z 46-53
  occ 54-59  b 60-65  '    ' chg(6) ' ' -> 66-76  type 77-78  ' ' 79
Verified to parse under BOTH LKina and Vina 1.2.7.
"""
import math

def atom_line(record, serial, name, resn, chain, resi, x, y, z, q, typ):
    """Standard PDBQT atom line (80 chars).

    Absolute 0-indexed columns: rec 0-5, serial 6-10, ' ' 11, name 12-15,
    altloc 16, resn 17-19, ' ' 20, chain 21, resi 22-25, '    ' 26-29,
    x 30-37, y 38-45, z 46-53, occ 54-59, b 60-65, '    '+chg+' ' 66-76,
    type 77-78, ' ' 79.
    Verified under LKina (grid + reactive atom resolution) and Vina 1.2.7.
    """
    return (f"{record:<6s}{serial:5d} {name:>4s} {resn:>3s} {chain:1s}{resi:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"    {q:>6s} {typ:<2s} ")

def norm(v):
    l = math.sqrt(sum(c*c for c in v))
    return tuple(c/l for c in v)

def write_receptor(path, tok, donors):
    """Metal at origin + donor shell + one C backbone atom."""
    lines = []
    s = 1
    lines.append(atom_line("ATOM", s, "CA", "ALA", "A", 301, -3.0, 0.0, 0.0, "+0.100", "C")); s += 1
    for i, dd in enumerate(donors):
        typ = "NA" if i % 2 == 0 else "OA"
        name = "ND1" if typ == "NA" else "OD1"
        resn = "HIS" if typ == "NA" else "ASP"
        lines.append(atom_line("ATOM", s, name, resn, "A", 301, *dd, "+0.000", typ)); s += 1
    lines.append(atom_line("HETATM", s, "M", tok, "A", 302, 0.0, 0.0, 0.0, "+0.000", tok))
    lines.append("END")
    open(path, "w").write("\n".join(lines) + "\n")

def write_ligand(path, enter, dist):
    nx, ny, nz = enter[0]*dist, enter[1]*dist, enter[2]*dist
    open(path, "w").write("\n".join([
      "REMARK  probe ligand (NA donor)",
      "ROOT",
      atom_line("ATOM", 1, "N", "LIG", "A", 1, nx, ny, nz, "-0.350", "NA"),
      atom_line("ATOM", 2, "C", "LIG", "A", 1, nx, ny, nz - 1.45, "+0.100", "C"),
      atom_line("ATOM", 3, "H", "LIG", "A", 1, nx + 0.9, ny, nz - 1.45, "+0.100", "HD"),
      "ENDROOT", "TORSDOF 0"]) + "\n")
