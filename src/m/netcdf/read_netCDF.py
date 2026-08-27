"""
read_netCDF  –  Load an ISSM model from a NetCDF4 file written by
                write_netCDF.py (Python) or write_netCDF.m (MATLAB).

Usage
-----
    from read_netCDF import read_netCDF
    md = read_netCDF('model.nc')
    md = read_netCDF('model.nc', verbose=True)

Every call returns an independent model() instance (fully re-entrant,
no module-level global state).
"""

import sys
import warnings
from collections import OrderedDict
import numpy as np
from netCDF4 import Dataset
from model import model
from results import results, solution, solutionstep
from toolkits import toolkits


def read_netCDF(filename: str, verbose: bool = False):
    """Read *filename* and return a populated ISSM model object."""
    if verbose:
        print('read_netCDF v2.0  (Python)')

    from os.path import exists
    if not exists(filename):
        raise FileNotFoundError(f'read_netCDF: file not found: {filename}')

    md = model()

    nc = Dataset(filename, 'r')
    nc.set_auto_mask(False)
    try:
        md = _read_all_groups(nc, md, verbose)
    finally:
        nc.close()

    if verbose:
        print('Model successfully loaded from NetCDF4.')
    return md


def _read_all_groups(nc, md, verbose: bool):  # {{{
    # =====================================================================
    #  _read_all_groups  –  iterate over every top-level group
    # =====================================================================
    for gname, grp in nc.groups.items():
        if verbose:
            print(f'  Reading group: {gname}')

        ct = _get_classtype(grp)

        # results group: walk solution sub-groups and build results()/solution() objects
        if gname == 'results':
            md.results = _read_results_group(grp, verbose)
            continue

        # toolkits group: each analysis sub-group holds an OrderedDict of solver options
        if gname == 'toolkits':
            md.toolkits = _read_toolkits_group(grp, verbose)
            continue

        # Instantiate the correct sub-class for polymorphic model fields
        md = _instantiate_subclass(md, gname, ct, verbose)

        # Populate all fields
        try:
            target = getattr(md, gname)
        except AttributeError:
            if verbose:
                print(f'    [SKIP] md.{gname} does not exist in model()')
            continue

        target = _read_group_into_obj(grp, target, verbose)
        setattr(md, gname, target)

    return md
# }}}
def _read_results_group(grp, verbose: bool):  # {{{
    # =====================================================================
    #  _read_results_group  –  builds results() / solution() / solutionstep()
    # =====================================================================
    # Reconstruct md.results as a proper results()/solution()/solutionstep()
    # hierarchy, matching what loadresultsfromdisk produces.
    res = results()

    for sol_name, sol_grp in grp.groups.items():
        ct = _get_classtype(sol_grp)

        if ct == 'struct':
            # Each step_k sub-group → one solutionstep
            steps = []
            for sg_name in sorted(sol_grp.groups.keys(),
                                  key=lambda n: int(n.split('_')[-1])
                                  if n.startswith('step_') else 0):
                sg   = sol_grp.groups[sg_name]
                step = solutionstep()
                for vname, var in sg.variables.items():
                    setattr(step, vname, _read_variable(sg, vname, var, verbose))
                # Recurse into any nested groups within a step (rare)
                for sub_name, sub_sg in sg.groups.items():
                    setattr(step, sub_name, _read_group_into_obj(sub_sg, {}, verbose))
                steps.append(step)

            if not steps:
                # Empty struct – store as empty solution
                setattr(res, sol_name, solution([]))
            elif len(steps) == 1 and sol_name != 'TransientSolution':
                # Single-step non-transient: store the solutionstep directly,
                # wrapped in a solution so attribute access still works
                sol_obj = solution(steps)
                setattr(res, sol_name, sol_obj)
            else:
                setattr(res, sol_name, solution(steps))

            if verbose:
                print(f'    [results] {sol_name}: {len(steps)}-step solution')
        else:
            # Unexpected sub-group inside results – fall back to dict
            if verbose:
                print(f'    [results] {sol_name}: unrecognised classtype "{ct}", stored as dict')
            setattr(res, sol_name, _read_group_into_obj(sol_grp, {}, verbose))

    return res
# }}}
def _read_toolkits_group(grp, verbose: bool):  # {{{
    # =====================================================================
    #  _read_toolkits_group  –  rebuilds toolkits() with OrderedDict per analysis
    # =====================================================================
    # Reconstruct md.toolkits as a toolkits() instance whose analysis
    # attributes are OrderedDicts of solver options (matching mumpsoptions() etc.)
    tk = toolkits.__new__(toolkits)   # bypass __init__ / setdefaultparameters
    tk.DefaultAnalysis  = None
    tk.RecoveryAnalysis = None

    for aname, agrp in grp.groups.items():
        ct = _get_classtype(agrp)
        opts = OrderedDict()

        if ct == 'struct':
            # Variables are stored one level deeper in step_1
            step_grps = agrp.groups
            src = step_grps.get('step_1', agrp) if step_grps else agrp
        else:
            src = agrp

        for vname, var in src.variables.items():
            opts[vname] = _read_variable(src, vname, var, verbose)

        setattr(tk, aname, opts)
        if verbose:
            print(f'    [toolkits] {aname}: {list(opts.keys())}')

    return tk
# }}}
def _read_group_into_obj(grp, obj, verbose: bool):  # {{{
    # =====================================================================
    #  _read_group_into_obj  –  populate obj from a NetCDF group, recursively
    # =====================================================================
    # Read variables at this level
    for vname, var in grp.variables.items():
        data = _read_variable(grp, vname, var, verbose)
        obj  = _set_attr(obj, vname, data, verbose)

    # Recurse into sub-groups
    for sg_name, sg in grp.groups.items():
        ct = _get_classtype(sg)

        if ct == 'struct':
            # Struct array (results.TransientSolution etc.)
            nsteps = int(getattr(sg, 'nsteps', len(sg.groups)))
            data   = _read_struct_array(sg, nsteps, verbose)

        elif ct == 'cell_of_objects':
            data = _read_cell_of_objects(sg, verbose)

        elif ct == 'dict':
            data = _read_dict(sg, verbose)

        else:
            # Regular sub-class or unknown: read into existing attr or new dict
            try:
                sub_obj = getattr(obj, sg_name)
            except AttributeError:
                sub_obj = {}

            # If the classtype tells us which class to instantiate
            if ct and ct not in ('', 'struct', 'cell_of_objects', 'dict'):
                sub_obj = _instantiate_class(ct, sub_obj, verbose)

            data = _read_group_into_obj(sg, sub_obj, verbose)

        obj = _set_attr(obj, sg_name, data, verbose)

    return obj
# }}}
def _read_struct_array(grp, nsteps: int, verbose: bool):  # {{{
    # =====================================================================
    #  _read_struct_array  –  reconstruct a list of dicts from step_k groups
    # =====================================================================
    steps = []
    for sg_name in sorted(grp.groups.keys(),
                          key=lambda n: int(n.split('_')[-1]) if n.startswith('step_') else 0):
        sg   = grp.groups[sg_name]
        step = {}
        for vname, var in sg.variables.items():
            step[vname] = _read_variable(sg, vname, var, verbose)
        # Recurse into any nested groups within a step (rare)
        for sub_name, sub_grp in sg.groups.items():
            step[sub_name] = _read_group_into_obj(sub_grp, {}, verbose)
        steps.append(step)

    if verbose:
        print(f'    [struct] loaded {len(steps)}-step struct array')
    return steps
# }}}
def _read_cell_of_objects(grp, verbose: bool):  # {{{
    # =====================================================================
    #  _read_cell_of_objects  –  reconstruct a nested list of ISSM objects
    # =====================================================================
    nrows = int(getattr(grp, 'nrows', 1))
    ncols = int(getattr(grp, 'ncols', len(grp.groups)))
    result = [[None] * ncols for _ in range(nrows)]

    for ig_name, ig in grp.groups.items():
        import re
        m = re.match(r'^item_(\d+)_(\d+)$', ig_name)
        if not m:
            continue
        r   = int(m.group(1)) - 1
        c   = int(m.group(2)) - 1
        ct  = _get_classtype(ig)
        obj = _instantiate_class(ct, {}, verbose)
        obj = _read_group_into_obj(ig, obj, verbose)
        result[r][c] = obj

    if nrows == 1:
        result = result[0]  # return flat list for 1-D case
    if verbose:
        print(f'    [cell]  loaded {nrows}x{ncols} cell of objects')
    return result
# }}}
def _read_dict(grp, verbose: bool) -> dict:  # {{{
    d = {}
    for vname, var in grp.variables.items():
        d[vname] = _read_variable(grp, vname, var, verbose)
    for sg_name, sg in grp.groups.items():
        d[sg_name] = _read_group_into_obj(sg, {}, verbose)
    return d
# }}}
def _read_variable(grp, vname: str, var, verbose: bool):  # {{{
    # =====================================================================
    #  _read_variable  –  read one NetCDF variable and convert to a Python type
    # =====================================================================
    type_is = getattr(var, 'type_is', '')

    # Empty sentinel (numeric empty or empty cell/list – both come back as [])
    if type_is in ('empty', 'cell_empty'):
        return []

    raw = var[:]

    # Bool
    if type_is == 'bool':
        return np.array(raw, dtype=bool)

    # String / cell-of-strings
    if type_is == 'string' or (hasattr(var, 'dimensions') and
                                any('char' in d for d in var.dimensions)):
        return _decode_string(raw)

    if type_is == 'cell_of_strings':
        return _decode_cell_strings(raw)

    # Numpy array / scalar
    data = np.array(raw)

    # Squeeze trailing size-1 dims
    data = np.squeeze(data)

    # Single-element array → Python scalar
    if data.ndim == 0:
        v = data.item()
        # Return int if value is a whole number (guard against nan/inf first)
        if isinstance(v, float) and np.isfinite(v) and v == int(v):
            return int(v)
        return v

    return data
# }}}
def _decode_string(raw) -> str:  # {{{
    """Convert raw NetCDF char data to a Python str."""
    try:
        return raw.tobytes().decode('utf-8').strip('\x00').strip()
    except Exception:
        try:
            return ''.join(c.decode('utf-8') if isinstance(c, bytes) else c
                           for c in np.ndarray.flatten(raw))
        except Exception:
            return str(raw)
# }}}
def _decode_cell_strings(raw) -> list:  # {{{
    """Convert a 2-D char array to a list of strings."""
    if raw.ndim == 1:
        return [_decode_string(raw)]
    result = []
    for row in raw:
        s = _decode_string(row)
        result.append(s)
    return result
# }}}
def _instantiate_subclass(md, gname: str, ct: str, verbose: bool):  # {{{
    # =====================================================================
    #  _instantiate_subclass  –  create the right class for polymorphic fields
    # =====================================================================
    if not ct:
        return md

    # Only act for fields whose type can vary
    if gname == 'inversion':
        if ct == 'm1qn3inversion':
            from m1qn3inversion import m1qn3inversion
            md.inversion = m1qn3inversion()
        elif ct == 'taoinversion':
            from taoinversion import taoinversion
            md.inversion = taoinversion()
        elif ct == 'inversion':
            from inversion import inversion
            md.inversion = inversion()
        else:
            warnings.warn(f"inversion class '{ct}' not supported yet")

    elif gname == 'smb':
        if ct == 'SMBforcing':
            from SMBforcing import SMBforcing
            md.smb = SMBforcing()
        elif ct == 'SMBpdd':
            from SMBpdd import SMBpdd
            md.smb = SMBpdd()
        elif ct == 'SMBd18opdd':
            from SMBd18opdd import SMBd18opdd
            md.smb = SMBd18opdd()
        elif ct == 'SMBgradients':
            from SMBgradients import SMBgradients
            md.smb = SMBgradients()
        elif ct == 'SMBcomponents':
            from SMBcomponents import SMBcomponents
            md.smb = SMBcomponents()
        elif ct == 'SMBmeltcomponents':
            from SMBmeltcomponents import SMBmeltcomponents
            md.smb = SMBmeltcomponents()
        elif ct == 'SMBgradientscomponents':
            from SMBgradientscomponents import SMBgradientscomponents
            md.smb = SMBgradientscomponents()
        elif ct == 'SMBgradientsela':
            from SMBgradientsela import SMBgradientsela
            md.smb = SMBgradientsela()
        elif ct == 'SMBhenning':
            from SMBhenning import SMBhenning
            md.smb = SMBhenning()
        elif ct == 'SMBgemb':
            from SMBgemb import SMBgemb
            md.smb = SMBgemb()
        elif ct == 'SMBpddSicopolis':
            from SMBpddSicopolis import SMBpddSicopolis
            md.smb = SMBpddSicopolis()
        elif ct == 'SMBsemic':
            from SMBsemic import SMBsemic
            md.smb = SMBsemic()
        elif ct == 'SMBdebrisEvatt':
            from SMBdebrisEvatt import SMBdebrisEvatt
            md.smb = SMBdebrisEvatt()
        else:
            warnings.warn(f"SMB class '{ct}' not supported yet")

    elif gname == 'friction':
        known = [
            'friction', 'frictioncoulomb', 'frictioncoulomb2',
            'frictionhydro', 'frictionjosh', 'frictionpism',
            'frictionregcoulomb', 'frictionregcoulomb2', 'frictionschoof',
            'frictionshakti', 'frictiontsai', 'frictionwaterlayer',
            'frictionweertman', 'frictionweertmantemp',
        ]
        if ct in known:
            mod = __import__(ct, fromlist=[ct])
            md.friction = getattr(mod, ct)()
        else:
            warnings.warn(f"friction class '{ct}' not supported yet")

    elif gname == 'hydrology':
        known = [
            'hydrologyshreve', 'hydrologydc', 'hydrologyglads',
            'hydrologypism', 'hydrologyshakti', 'hydrologytws', 'hydrologyarmapw',
        ]
        if ct in known:
            mod = __import__(ct, fromlist=[ct])
            md.hydrology = getattr(mod, ct)()
        else:
            warnings.warn(f"hydrology class '{ct}' not supported yet")

    elif gname == 'mesh':
        if ct == 'mesh3dprisms':
            from mesh3dprisms import mesh3dprisms
            md.mesh = mesh3dprisms()
        elif ct == 'mesh2d':
            from mesh2d import mesh2d
            md.mesh = mesh2d()
        else:
            warnings.warn(f"mesh class '{ct}' not supported yet")

    # All other groups keep their default instance; we just overwrite fields
    return md
# }}}
def _instantiate_class(ct: str, fallback, verbose: bool):  # {{{
    """Try to instantiate *ct* by importing it; return *fallback* on failure."""
    if not ct or ct in ('struct', 'cell_of_objects', 'dict', ''):
        return fallback
    try:
        mod = __import__(ct, fromlist=[ct])
        cls = getattr(mod, ct)
        return cls()
    except Exception:
        pass
    # Try to find ct in already-imported modules
    for mod_name, mod in list(sys.modules.items()):
        if mod is not None and hasattr(mod, ct):
            cls = getattr(mod, ct)
            if isinstance(cls, type):
                try:
                    return cls()
                except Exception:
                    pass
    if verbose:
        print(f'    [WARN] could not instantiate class {ct}, using fallback')
    return fallback
# }}}
def _set_attr(obj, name: str, value, verbose: bool):  # {{{
    # =====================================================================
    #  _set_attr  –  safely set an attribute on either a class or dict
    # =====================================================================
    try:
        if isinstance(obj, dict):
            obj[name] = value
        else:
            setattr(obj, name, value)
    except Exception as e:
        if verbose:
            print(f'    [WARN] could not set {name}: {e}')
    return obj
# }}}
def _get_classtype(grp) -> str:  # {{{
    return getattr(grp, 'classtype', '')
# }}}
