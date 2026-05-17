"""
write_netCDF  –  Save an ISSM model to a NetCDF4 file.

Usage
-----
    from write_netCDF import write_netCDF
    write_netCDF(md, 'model.nc')
    write_netCDF(md, 'model.nc', verbose=True)

The file can be read back with read_netCDF.py (Python) or read_netCDF.m (MATLAB).

NetCDF layout
-------------
Each attribute of md (mesh, mask, geometry, …) becomes a NetCDF group.
Attributes within those objects are stored as variables inside the group.
Nested sub-classes get nested sub-groups.

Every group carries a 'classtype' global attribute with the Python class name
so that read_netCDF can reconstruct the correct object on load.

Results storage
---------------
md.results.<SolutionName> is a solution() instance whose .steps list holds
solutionstep() objects.  Each step is stored as a sub-group named 'step_<k>'
inside the solution sub-group.
"""

import os
import time
import numpy as np
from netCDF4 import Dataset


def write_netCDF(md, filename: str, verbose: bool = False) -> None:
    """Save model *md* to *filename* (NetCDF4).  Overwrites silently."""
    if verbose:
        print('write_netCDF v2.0  (Python)')

    # Overwrite silently
    if os.path.exists(filename):
        os.remove(filename)
        if verbose:
            print(f'Existing file {filename} deleted.')

    nc = Dataset(filename, 'w', format='NETCDF4')
    try:
        nc.history     = 'Created by write_netCDF.py on ' + time.ctime()
        nc.Conventions = 'ISSM-NetCDF-2.0'

        # Global dimensions reused throughout the file
        nc.createDimension('scalar', 1)
        nc.createDimension('Unlim',  None)   # unlimited

        # Walk every top-level attribute of the model
        for gname, val in md.__dict__.items():
            if val is None:
                if verbose:
                    print(f'  md.{gname} is None – skipped')
                continue

            grp = nc.createGroup(gname)
            grp.classtype = type(val).__name__

            _write_obj(val, grp, nc, f'md.{gname}', verbose)

    except Exception:
        nc.close()
        raise

    nc.close()
    if verbose:
        print('Model successfully saved as NetCDF4.')


# =====================================================================
#  _write_obj  –  write all attributes of a Python object into a group
# =====================================================================
def _write_obj(obj, grp, root_nc, path: str, verbose: bool) -> None:
    """Iterate over obj.__dict__ and write each attribute into *grp*."""
    try:
        items = obj.__dict__.items()
    except AttributeError:
        return   # not an object with attributes

    for fname, val in items:
        if fname in ('errlog', 'outlog'):
            continue
        _write_value(val, fname, grp, root_nc, f'{path}.{fname}', verbose)


# =====================================================================
#  _write_value  –  dispatch a single value to the right writer
# =====================================================================
def _write_value(val, vname: str, grp, root_nc, path: str, verbose: bool) -> None:

    # ---- None / empty -----------------------------------------------
    if val is None:
        _write_scalar_nan(vname, grp, root_nc)
        if verbose:
            print(f'    [None]    {path}')
        return

    # ---- numpy ndarray ----------------------------------------------
    if isinstance(val, np.ndarray):
        _write_numpy(val, vname, grp, root_nc, path, verbose)
        return

    # ---- bool (must come before int because bool is a subclass of int)
    if isinstance(val, (bool, np.bool_)):
        _write_bool(np.array([int(val)], dtype=np.int16), vname, grp, root_nc)
        if verbose:
            print(f'    [bool]    {path}')
        return

    # ---- int / float scalars ----------------------------------------
    if isinstance(val, (int, np.integer)):
        _write_scalar(int(val), vname, 'i8', grp, root_nc)
        if verbose:
            print(f'    [int]     {path}')
        return

    if isinstance(val, (float, np.floating)):
        _write_scalar(float(val), vname, 'f8', grp, root_nc)
        if verbose:
            print(f'    [float]   {path}')
        return

    # ---- str --------------------------------------------------------
    if isinstance(val, str):
        _write_string(val, vname, grp)
        if verbose:
            print(f'    [string]  {path}')
        return

    # ---- list -------------------------------------------------------
    if isinstance(val, list):
        _write_list(val, vname, grp, root_nc, path, verbose)
        return

    # ---- dict / OrderedDict -----------------------------------------
    if isinstance(val, dict):
        _write_dict(val, vname, grp, root_nc, path, verbose)
        return

    # ---- ISSM sub-class (has __dict__) ------------------------------
    if hasattr(val, '__dict__'):
        sub_grp = grp.createGroup(vname)
        sub_grp.classtype = type(val).__name__
        _write_obj(val, sub_grp, root_nc, path, verbose)
        return

    # ---- fallback ---------------------------------------------------
    if verbose:
        print(f'    [SKIP]    {path} – unsupported type {type(val).__name__}')


# =====================================================================
#  _write_list  –  dispatch a Python list
# =====================================================================
def _write_list(val: list, vname: str, grp, root_nc, path: str, verbose: bool) -> None:
    if len(val) == 0:
        # Empty list: store as a scalar NaN sentinel
        _write_scalar_nan(vname, grp, root_nc)
        return

    first = val[0]

    # list of ISSM objects → cell_of_objects group
    if hasattr(first, '__dict__') and not isinstance(first, dict):
        _write_cell_of_objects(val, vname, grp, root_nc, path, verbose)
        return

    # list of str
    if isinstance(first, str):
        _write_cell_strings(val, vname, grp)
        if verbose:
            print(f'    [cellstr] {path}')
        return

    # list of solution-step objects (results.solution.steps)
    # This is handled by _write_solution_steps when the parent is a solution();
    # if we land here generically, treat as cell_of_objects
    if hasattr(first, '__dict__'):
        _write_cell_of_objects(val, vname, grp, root_nc, path, verbose)
        return

    # Plain numeric list → convert to numpy array
    try:
        arr = np.array(val)
        _write_numpy(arr, vname, grp, root_nc, path, verbose)
    except Exception:
        if verbose:
            print(f'    [SKIP]    {path} – list of mixed/unsupported types')


# =====================================================================
#  Solution / solutionstep special handling
# =====================================================================
def _write_obj_with_results(obj, grp, root_nc, path: str, verbose: bool) -> None:
    """Like _write_obj but handles solution().steps lists as step_k subgroups."""
    try:
        items = obj.__dict__.items()
    except AttributeError:
        return

    for fname, val in items:
        if fname in ('errlog', 'outlog'):
            continue

        # Detect a list of solutionstep-like objects
        if (isinstance(val, list) and len(val) > 0
                and hasattr(val[0], '__dict__')
                and not isinstance(val[0], dict)):
            _write_solution_steps(val, fname, grp, root_nc, path, verbose)
        else:
            _write_value(val, fname, grp, root_nc, f'{path}.{fname}', verbose)


def _write_solution_steps(steps: list, vname: str, grp, root_nc, path: str, verbose: bool) -> None:
    """Write a list of solutionstep objects as step_1, step_2, … sub-groups."""
    sub_grp = grp.createGroup(vname)
    sub_grp.classtype = 'struct'          # MATLAB struct array convention
    sub_grp.nsteps    = len(steps)

    for k, step in enumerate(steps, start=1):
        sg = sub_grp.createGroup(f'step_{k}')
        for sfield, sval in step.__dict__.items():
            if sfield in ('errlog', 'outlog'):
                continue
            _write_value(sval, sfield, sg, root_nc, f'{path}.step_{k}.{sfield}', verbose)

    if verbose:
        print(f'    [struct]  {path}.{vname}  ({len(steps)} steps)')


# =====================================================================
#  Leaf writers
# =====================================================================

def _write_numpy(arr: np.ndarray, vname: str, grp, root_nc, path: str, verbose: bool) -> None:
    if arr.size == 0:
        _write_scalar_nan(vname, grp, root_nc)
        return

    # Booleans → int16
    if arr.dtype == bool:
        arr = arr.astype(np.int16)
        is_bool = True
    else:
        is_bool = False

    # Squeeze trailing size-1 dimensions; keep 0-d as scalar
    arr = np.squeeze(arr)
    if arr.ndim == 0:
        arr = arr.reshape(1)

    # Determine NetCDF dtype
    if np.issubdtype(arr.dtype, np.integer):
        nc_dtype = 'i8'
    else:
        nc_dtype = 'f8'
        arr = arr.astype(np.float64)

    dims = _make_dims(arr.shape, grp, root_nc)

    var = grp.createVariable(vname, nc_dtype, dims, zlib=(arr.size > 1))
    if is_bool:
        var.type_is = 'bool'
    var[:] = arr

    if verbose:
        print(f'    [numeric] {path}  shape={arr.shape}  dtype={arr.dtype}')


def _write_scalar(val, vname: str, nc_dtype: str, grp, root_nc) -> None:
    var = grp.createVariable(vname, nc_dtype, ('scalar',))
    var[:] = val


def _write_scalar_nan(vname: str, grp, root_nc) -> None:
    var = grp.createVariable(vname, 'f8', ('scalar',))
    var[:] = np.nan
    var.type_is = 'empty'


def _write_string(s: str, vname: str, grp) -> None:
    if len(s) == 0:
        dim_name = 'char0'
        try:
            grp.createDimension(dim_name, 0)
        except Exception:
            pass
        v = grp.createVariable(vname, 'S1', (dim_name,))
        v.type_is = 'string'
        return
    n        = len(s)
    dim_name = f'char{n}'
    try:
        grp.createDimension(dim_name, n)
    except Exception:
        pass
    v = grp.createVariable(vname, 'S1', (dim_name,))
    from netCDF4 import stringtochar
    v[:] = stringtochar(np.array([s], dtype=f'S{n}'))
    v.type_is = 'string'


def _write_cell_strings(strings: list, vname: str, grp) -> None:
    if not strings:
        dim_name = 'char0'
        try:
            grp.createDimension(dim_name, 0)
        except Exception:
            pass
        v = grp.createVariable(vname, 'S1', (dim_name,))
        v.type_is = 'cell_of_strings'
        return

    nrows   = len(strings)
    max_len = max(len(s) for s in strings)
    if max_len == 0:
        max_len = 1

    rows_name = f'strrows_{vname}'
    cols_name = f'strcols_{vname}'
    try:
        grp.createDimension(rows_name, nrows)
    except Exception:
        pass
    try:
        grp.createDimension(cols_name, max_len)
    except Exception:
        pass

    v = grp.createVariable(vname, 'S1', (rows_name, cols_name))
    v.type_is = 'cell_of_strings'

    arr = np.chararray((nrows, max_len))
    for i, s in enumerate(strings):
        padded = s.ljust(max_len)
        arr[i] = list(padded.encode('utf-8'))
    v[:] = arr


def _write_bool(arr: np.ndarray, vname: str, grp, root_nc) -> None:
    if arr.size == 1:
        v = grp.createVariable(vname, 'i2', ('scalar',))
        v[:] = arr[0]
    else:
        n = arr.size
        dim_name = f'booldim{n}'
        try:
            grp.createDimension(dim_name, n)
        except Exception:
            pass
        v = grp.createVariable(vname, 'i2', (dim_name,))
        v[:] = arr.flatten()
    v.type_is = 'bool'


def _write_cell_of_objects(val: list, vname: str, grp, root_nc, path: str, verbose: bool) -> None:
    sub_grp = grp.createGroup(vname)
    sub_grp.classtype = 'cell_of_objects'
    sub_grp.nrows     = 1
    sub_grp.ncols     = len(val)

    for c, elem in enumerate(val, start=1):
        eg = sub_grp.createGroup(f'item_1_{c}')
        eg.classtype = type(elem).__name__
        if hasattr(elem, '__dict__'):
            _write_obj(elem, eg, root_nc, f'{path}[{c}]', verbose)

    if verbose:
        print(f'    [cell]    {path}  ({len(val)} objects)')


def _write_dict(val: dict, vname: str, grp, root_nc, path: str, verbose: bool) -> None:
    """Store a dict as a sub-group with one string variable per key."""
    sub_grp = grp.createGroup(vname)
    sub_grp.classtype = 'dict'
    for k, v in val.items():
        _write_value(v, str(k), sub_grp, root_nc, f'{path}.{k}', verbose)


# =====================================================================
#  Dimension helper
# =====================================================================

def _make_dims(shape: tuple, grp, root_nc) -> tuple:
    """Return a tuple of dimension-name strings, creating dims as needed."""
    dim_names = []
    for n in shape:
        if n == 1:
            dim_names.append('scalar')
        else:
            dname = f'dim{n}'
            # Try to create in local group; if it exists anywhere, reuse it
            for target in (grp, root_nc):
                try:
                    target.createDimension(dname, n)
                    break
                except Exception:
                    pass   # already exists – that's fine
            dim_names.append(dname)
    return tuple(dim_names)
