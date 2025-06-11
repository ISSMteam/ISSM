#!/usr/bin/env python3
import argparse, struct, pathlib
from typing import Union

FMT2CODE = {"bool":1,"int":2,"double":3}

def locate(fp, field):
    with fp.open('rb') as f:
        while (raw:=f.read(4)):
            namelen = struct.unpack('i', raw)[0]
            name    = f.read(namelen).decode()
            reclen  = struct.unpack('q', f.read(8))[0]
            codepos = f.tell(); code = struct.unpack('i', f.read(4))[0]
            datapos = f.tell()
            if name == field: return datapos, code
            f.seek(reclen-4,1)
    raise KeyError(field)

def patch(fp: pathlib.Path, field: str, value: Union[str,int,float]):
    datapos, code = locate(fp, field)
    with fp.open('r+b') as f:
        f.seek(datapos)
        if code==FMT2CODE["double"]:
            f.write(struct.pack('d', float(value)))
        elif code==FMT2CODE["int"]:
            f.write(struct.pack('i', int(value)))
        elif code==FMT2CODE["bool"]:
            f.write(struct.pack('i', 1 if value.lower()=='true' else 0))
        else:
            raise TypeError("Only scalar bool/int/double supported in-place")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-f","--file", required=True)
    p.add_argument("-n","--name", required=True)
    p.add_argument("-v","--value", required=True)
    args = p.parse_args()
    patch(pathlib.Path(args.file), args.name, args.value)
