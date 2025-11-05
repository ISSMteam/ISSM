#!/usr/bin/env python3
"""
patch_bin.py — edit any ISSM .bin record.

Examples
--------
# change a scalar double in place:
patch_bin.py -f run.bin -n md.timestepping.final_time -v 3.1536e7

# replace a 10×10 matrix with a different-sized one (auto-rebuilds file):
patch_bin.py -f run.bin -n md.geometry.surface \
             -v '[ [0,1],[2,3] ]' --rebuild
"""
from __future__ import annotations
import argparse, json, struct, pathlib, tempfile, shutil, sys
import numpy as np
FMT = {"Boolean":1,"Integer":2,"Double":3,"String":4,
       "BooleanMat":5,"IntMat":6,"DoubleMat":7,"MatArray":8,"StringArray":9}
CODE2FMT = {v:k for k,v in FMT.items()}

# ─────────────────────────── helper: iterate over records ────────────────────
def iterate_records(fp: pathlib.Path):
    """Yield (name,pos,reclen,code,data_pos) for every record in *fp*."""
    with fp.open("rb") as f:
        while True:
            raw = f.read(4)
            if not raw: break
            namelen = struct.unpack("i", raw)[0]
            name = f.read(namelen).decode()
            reclen = struct.unpack("q", f.read(8))[0]
            code   = struct.unpack("i", f.read(4))[0]
            data_pos = f.tell()
            yield name, f.tell() - 4, reclen, code, data_pos
            f.seek(reclen-4, 1)

# ─────────────────────────── decode / encode payloads ────────────────────────
def read_payload(f, code, reclen):
    if code in (FMT["Boolean"], FMT["Integer"]):
        return struct.unpack("i", f.read(4))[0]
    if code == FMT["Double"]:
        return struct.unpack("d", f.read(8))[0]
    if code == FMT["String"]:
        strlen = struct.unpack("i", f.read(4))[0]
        return f.read(strlen).decode()
    if code in (FMT["BooleanMat"], FMT["IntMat"], FMT["DoubleMat"]):
        mattype = struct.unpack("i", f.read(4))[0]
        n,m = struct.unpack("ii", f.read(8))
        data = np.frombuffer(f.read(n*m*8), dtype=np.float64).reshape(n,m)
        return mattype, data
    if code == FMT["MatArray"]:
        # skip detailed parsing, return raw bytes
        return f.read(reclen-4)
    if code == FMT["StringArray"]:
        return f.read(reclen-4)
    raise TypeError("unknown code")

def encode_payload(val, code):
    """Return bytes for *val* plus the 'code' (4B) and new reclen."""
    if code in (FMT["Boolean"], FMT["Integer"]):
        payload = struct.pack("i", int(val))
    elif code == FMT["Double"]:
        payload = struct.pack("d", float(val))
    elif code == FMT["String"]:
        b = str(val).encode()
        payload = struct.pack("i", len(b)) + b
    elif code in (FMT["BooleanMat"], FMT["IntMat"], FMT["DoubleMat"]):
        mattype, arr = val                         # expect (mattype, np.ndarray)
        n,m = arr.shape
        payload = struct.pack("i", int(mattype))
        payload += struct.pack("ii", n, m)
        payload += arr.astype(np.float64).tobytes(order="C")
    elif code == FMT["MatArray"]:
        payload = bytes(val)                       # raw bytes
    elif code == FMT["StringArray"]:
        payload = bytes(val)
    else:
        raise TypeError("unsupported code")
    reclen = len(payload)+4                        # +4 for the code itself
    return struct.pack("i", code) + payload, reclen

# ─────────────────────────── patch logic ─────────────────────────────────────
def patch(file: pathlib.Path, field: str, new_val, rebuild=False):
    records = list(iterate_records(file))
    match = next((r for r in records if r[0]==field), None)
    if not match:
        raise KeyError(f"{field} not found")
    name, name_pos, old_reclen, code, data_pos = match
    new_payload, new_reclen = encode_payload(new_val, code)

    if not rebuild and new_reclen==old_reclen:
        # simple in-place overwrite
        with file.open("r+b") as f:
            f.seek(data_pos-4)           # back up to code
            f.write(new_payload)
        return

    # otherwise rebuild the entire file
    tmp = pathlib.Path(tempfile.mktemp(suffix=".bin", dir=file.parent))
    with file.open("rb") as src, tmp.open("wb") as dst:
        for n, pos, reclen, c, dpos in records:
            src.seek(pos)
            chunk = src.read(4 + len(n.encode()) + 8 + reclen)  # full record
            if n==field:
                # slice header up to code (everything before code)
                header = chunk[:4+len(n.encode())+8]
                dst.write(header + new_payload)
            else:
                dst.write(chunk)
    shutil.move(tmp, file)

# ─────────────────────────── CLI ─────────────────────────────────────────────
def cli():
    p = argparse.ArgumentParser()
    p.add_argument("-f","--file", required=True)
    p.add_argument("-n","--name", required=True)
    p.add_argument("-v","--value", required=True,
                   help="Scalar, JSON array, or JSON object depending on type.")
    p.add_argument("--rebuild", action="store_true",
                   help="Force full-file rewrite even if size matches.")
    args = p.parse_args()

    # value parsing: try JSON first, fall back to raw string / float
    try:
        val = json.loads(args.value)
        if isinstance(val, list):  # assume matrix if nested list
            val = (0, np.array(val, dtype=np.float64))
    except Exception:
        try: val = float(args.value)
        except ValueError: val = args.value

    patch(pathlib.Path(args.file), args.name, val, rebuild=args.rebuild)
    print(f"✔ Patched {args.name} in {args.file}")

if __name__ == "__main__": cli()
