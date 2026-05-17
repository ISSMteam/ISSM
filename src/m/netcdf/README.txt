write_netCDF / read_netCDF  –  Cross-language ISSM model serialisation
======================================================================

These four files let you save an ISSM model (md) to NetCDF4 and reload it
in either MATLAB or Python, interchangeably.

    write_netCDF.m   –  save  from MATLAB
    read_netCDF.m    –  load  into MATLAB
    write_netCDF.py  –  save  from Python
    read_netCDF.py   –  load  into Python

Version: 2.0


USAGE
-----

Python:
    from write_netCDF import write_netCDF
    from read_netCDF  import read_netCDF

    write_netCDF(md, 'model.nc')
    write_netCDF(md, 'model.nc', verbose=True)

    md2 = read_netCDF('model.nc')
    md2 = read_netCDF('model.nc', verbose=True)

MATLAB:
    write_netCDF(md, 'model.nc')
    write_netCDF(md, 'model.nc', 'verbose', true)

    md2 = read_netCDF('model.nc')
    md2 = read_netCDF('model.nc', 'verbose', true)

Cross-language round-trip:
    %% MATLAB → Python
    write_netCDF(md, 'model.nc');        % MATLAB
    md2 = read_netCDF('model.nc')        # Python

    ## Python → MATLAB
    write_netCDF(md, 'model.nc')         # Python
    md2 = read_netCDF('model.nc');       % MATLAB


DEPENDENCIES
------------
Python:  numpy, netCDF4, and the ISSM model() class (plus any sub-class
         modules that your model uses, e.g. m1qn3inversion, SMBpdd, etc.)
MATLAB:  built-in netcdf package (no extra toolboxes required), plus
         the ISSM model() class and any sub-class .m files.


FILE FORMAT  (ISSM-NetCDF-2.0)
------------------------------
Global attributes:
    Conventions = 'ISSM-NetCDF-2.0'
    history     = creation timestamp

Layout:
    /mesh/          ← group per top-level md field
        classtype   ← NC_GLOBAL attribute: Python/MATLAB class name
        x           ← variable
        y           ← variable
        …
    /geometry/
        …
    /results/
        /TransientSolution/   ← sub-group (classtype='struct')
            nsteps            ← NC_GLOBAL attribute: number of steps
            /step_1/
                Vel           ← variable
                …
            /step_2/
                …
    /inversion/
        classtype = 'm1qn3inversion'   ← enables correct class reconstruction
        …

Key design decisions:
  - Every group carries a 'classtype' global attribute so readers can
    reconstruct the exact Python/MATLAB class without heuristics.
  - Struct arrays (results.TransientSolution) are stored as step_1 … step_n
    sub-groups, one per time step.
  - The file is overwritten silently if it already exists (no interactive
    prompts – safe for use in scripts and Jupyter notebooks).
  - No module-level global state in Python (fully re-entrant; calling
    read_netCDF twice returns two independent model objects).
  - No MATLAB 'persistent' variable bugs (the old code skipped writing
    inversion/smb/friction/hydrology class names on the second call).


KNOWN LIMITATIONS
-----------------
- md.qmu  is not yet supported (complex OrderedDict-based structure).
- Cell arrays whose elements are themselves struct arrays are not yet
  supported.
- MATLAB classes that have no Python equivalent (e.g. SMBgemb) can be
  saved from MATLAB and round-tripped back to MATLAB, but the Python
  reader will fall back to a plain dict for those fields.
- Very large arrays (> a few GB) work fine thanks to NetCDF4 chunking,
  but you may want to increase the zlib compression level in write_netCDF.py
  for highly compressible fields.


OTHER FILES IN THIS DIRECTORY
------------------------------
export_netCDF.m / export_netCDF.py
    An older, alternative serialisation that stacks results variables by
    time dimension rather than using per-step sub-groups.  It writes a
    different file format and does not have a matching read_netCDF.
    Kept for backward compatibility.

restable.m
    Helper used by export_netCDF.m.
