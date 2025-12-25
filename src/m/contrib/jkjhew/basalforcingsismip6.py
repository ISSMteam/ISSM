import numpy as np

# ISSM helper imports (same APIs as MATLAB equivalents)
from project3d import project3d
from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class basalforcingsismip6:
    """
    Python port of ISSM's basalforcingsismip6 class.

    Usage:
        bf = basalforcingsismip6()
        # or
        bf = basalforcingsismip6(existing_obj_dict)  # structtoobj equivalent
    """

    # --------------------- Constructor ---------------------
    def __init__(self, *args):
        # Defaults (match MATLAB)
        self.num_basins               = 0
        self.basin_id                 = np.nan
        self.gamma_0                  = 0.0
        self.tf                       = np.nan
        self.tf_depths                = np.nan
        self.delta_t                  = np.nan
        self.islocal                  = False
        self.geothermalflux           = np.nan
        self.groundedice_melting_rate = np.nan
        self.melt_anomaly             = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        elif len(args) == 1:
            self.setdefaultparameters()
            # structtoobj: copy matching keys
            src = args[0]
            if isinstance(src, dict):
                for k, v in src.items():
                    if hasattr(self, k):
                        setattr(self, k, v)
            else:
                # Best-effort: copy attributes
                for k in dir(src):
                    if not k.startswith("_") and hasattr(self, k):
                        setattr(self, k, getattr(src, k))
        else:
            raise RuntimeError("constructor not supported")

    # --------------------- Methods -------------------------
    def extrude(self, md):
        """
        Match MATLAB:
          basin_id: element, layer=1
          tf: node (time series OK, handled as 2D [nNodes+1, nTimes])
          geothermalflux: element, layer=1
          groundedice_melting_rate: node, layer=1
          melt_anomaly: element, layer=1
        """
        # basin_id is element-wise scalar
        self.basin_id = project3d(md, "vector", self.basin_id, "type", "element", "layer", 1)

        # tf: list/array over depths; each item is (nNodes+1, nTimes).
        # Project node-wise vectors; keep the time axis intact.
        if isinstance(self.tf, (list, tuple)):
            new_tf = []
            for arr in self.tf:
                if isinstance(arr, np.ndarray) and arr.ndim == 2:
                    # Last row is time; separate before projection
                    nodal = arr[:-1, :]
                    years = arr[-1:, :]
                    nodal_proj = project3d(md, "vector", nodal, "type", "node")
                    new_tf.append(np.vstack([nodal_proj, years]))
                else:
                    new_tf.append(arr)
            self.tf = new_tf
        elif isinstance(self.tf, np.ndarray) and self.tf.ndim == 3:
            # Optional: depth-major array; treat last row as years per slice
            tf_list = []
            for iz in range(self.tf.shape[0]):
                arr = self.tf[iz, :, :]
                nodal = arr[:-1, :]
                years = arr[-1:, :]
                nodal_proj = project3d(md, "vector", nodal, "type", "node")
                tf_list.append(np.vstack([nodal_proj, years]))
            self.tf = tf_list

        # geothermalflux: element, bedrock only (layer=1)
        if isinstance(self.geothermalflux, np.ndarray):
            self.geothermalflux = project3d(md, "vector", self.geothermalflux, "type", "element", "layer", 1)

        # groundedice_melting_rate: node, layer=1
        if isinstance(self.groundedice_melting_rate, np.ndarray):
            self.groundedice_melting_rate = project3d(
                md, "vector", self.groundedice_melting_rate, "type", "node", "layer", 1
            )

        # melt_anomaly: element, layer=1
        if isinstance(self.melt_anomaly, np.ndarray):
            self.melt_anomaly = project3d(md, "vector", self.melt_anomaly, "type", "element", "layer", 1)

        return self

    def initialize(self, md):
        # gamma_0 default if 0
        if self.gamma_0 == 0:
            self.gamma_0 = 14477
            print("      no basalforcings.gamma_0 specified: value set to 14477 m/yr")
        # groundedice_melting_rate default to zeros if NaN
        if isinstance(self.groundedice_melting_rate, float) and np.isnan(self.groundedice_melting_rate):
            self.groundedice_melting_rate = np.zeros(md.mesh.numberofvertices)
            print("      no basalforcings.groundedice_melting_rate specified: values set as zero")
        return self

    def setdefaultparameters(self):
        self.gamma_0 = 14477.0  # m/yr
        self.islocal = False
        return self

    def checkconsistency(self, md, solution=None, analyses=None):
        # num_basins
        md = checkfield(md, fieldname="basalforcings.num_basins", numel=1, NaN=1, Inf=1, _gt=0)

        # basin_id ∈ [0, num_basins], element-sized
        md = checkfield(
            md,
            fieldname="basalforcings.basin_id",
            Inf=1,
            _ge=0,
            _le=md.basalforcings.num_basins,
            size=[md.mesh.numberofelements, 1],
        )

        # gamma_0 > 0
        md = checkfield(md, fieldname="basalforcings.gamma_0", numel=1, NaN=1, Inf=1, _gt=0)

        # tf_depths: row vector, <= 0 (depths are non-positive)
        md = checkfield(md, fieldname="basalforcings.tf_depths", NaN=1, Inf=1, size=[1, None], _le=0)

        # delta_t: length = num_basins, row vector
        md = checkfield(
            md,
            fieldname="basalforcings.delta_t",
            NaN=1,
            Inf=1,
            numel=md.basalforcings.num_basins,
            size=[1, md.basalforcings.num_basins],
        )

        # islocal: boolean {0,1}
        md = checkfield(md, fieldname="basalforcings.islocal", values=[0, 1])

        # geothermalflux: timeseries on elements
        md = checkfield(md, fieldname="basalforcings.geothermalflux", NaN=1, Inf=1, _ge=0, timeseries=1)

        # groundedice_melting_rate: timeseries on nodes
        md = checkfield(md, fieldname="basalforcings.groundedice_melting_rate", NaN=1, Inf=1, timeseries=1)

        # melt_anomaly: only check if array is used (length>1)

        if np.size(self.melt_anomaly) > 1:
            md = checkfield(md, fieldname="basalforcings.melt_anomaly", NaN=1, Inf=1, timeseries=1)


        # tf cell structure: [1,1,num_depths], each cell (nNodes+1 x nTime)
        depths = np.size(self.tf_depths)
        md = checkfield(md, fieldname="basalforcings.tf", size=[1, 1, depths])
        # Validate each depth's matrix
        def _matrix_ok(mat):
            return (
                isinstance(mat, np.ndarray)
                and mat.ndim == 2
                and mat.shape[0] == md.mesh.numberofvertices + 1
            )

        if isinstance(self.tf, (list, tuple)) and len(self.tf) == depths:
            for i in range(depths):
                mat = self.tf[i]
                if not _matrix_ok(mat):
                    raise ValueError(
                        f"basalforcings.tf[{i}] must be (nVertices+1) x nTime with last row = years."
                    )
                # MATLAB also checks NaN/Inf/>=0 and timeseries; rely on WriteData handling plus basic numeric sanity:
                if not np.all(np.isfinite(mat)):
                    raise ValueError(f"basalforcings.tf[{i}] contains NaN/Inf.")
                if np.any(mat[:-1, :] < 0):
                    raise ValueError(f"basalforcings.tf[{i}] has negative TF values (excluding last time row).")
        else:
            raise ValueError("basalforcings.tf must be a list of length numel(tf_depths).")

        return md

    def disp(self):
        print("   ISMIP6 basal melt rate parameterization:")
        fielddisplay(self, "num_basins", "number of basins the model domain is partitioned into [unitless]")
        fielddisplay(self, "basin_id", "basin number assigned to each node (unitless)")
        fielddisplay(self, "gamma_0", "melt rate coefficient (m/yr)")
        fielddisplay(self, "tf_depths", "elevation of vertical layers in ocean thermal forcing dataset")
        fielddisplay(self, "tf", "thermal forcing (ocean temperature minus freezing point) (degrees C)")
        fielddisplay(self, "delta_t", "Ocean temperature correction per basin (degrees C)")
        fielddisplay(self, "islocal", "boolean to use local ISMIP6 melt-rate parameterization (default false)")
        fielddisplay(self, "geothermalflux", "geothermal heat flux (W/m^2)")
        fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) (m/yr)")
        fielddisplay(self, "melt_anomaly", "floating ice basal melt anomaly (m/yr)")

    def marshall(self, prefix, md, fid):
        """
        Binary marshalling identical to MATLAB version.
        Notes:
          - basin_id is written 0-indexed (subtract 1).
          - gamma_0 is scaled by 1/yts to convert from m/yr to m/s in the file.
          - Time-series arrays use timeserieslength consistent with MATLAB:
              tf, delta_t, groundedice_melting_rate: nVertices+1
              geothermalflux, melt_anomaly: nElements+1
        """
        yts = md.constants.yts

        # Model id = 7 (ISMIP6 basal melt parameterization)
        WriteData(fid, prefix, name="md.basalforcings.model", data=7, format="Integer")

        WriteData(fid, prefix, obj=self, fieldname="num_basins", format="Integer")
        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="basin_id",
            data=(self.basin_id - 1),  # 0-indexed
            name="md.basalforcings.basin_id",
            format="IntMat",
            mattype=2,
        )
        WriteData(fid, prefix, obj=self, fieldname="gamma_0", format="Double", scale=1.0 / yts)
        WriteData(fid, prefix, obj=self, fieldname="tf_depths", format="DoubleMat", name="md.basalforcings.tf_depths")

        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="tf",
            format="MatArray",
            name="md.basalforcings.tf",
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=yts,
        )
        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="delta_t",
            format="DoubleMat",
            name="md.basalforcings.delta_t",
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=yts,
        )
        WriteData(fid, prefix, obj=self, fieldname="islocal", format="Boolean")

        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="geothermalflux",
            format="DoubleMat",
            name="md.basalforcings.geothermalflux",
            mattype=1,
            timeserieslength=md.mesh.numberofelements + 1,
            yts=yts,
        )
        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="groundedice_melting_rate",
            format="DoubleMat",
            mattype=1,
            scale=1.0 / yts,
            timeserieslength=md.mesh.numberofvertices + 1,
            yts=yts,
        )
        WriteData(
            fid,
            prefix,
            obj=self,
            fieldname="melt_anomaly",
            format="DoubleMat",
            mattype=1,
            scale=1.0 / yts,
            timeserieslength=md.mesh.numberofelements + 1,
            yts=yts,
        )
