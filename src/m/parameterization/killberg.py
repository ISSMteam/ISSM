### October 2025: A fast implementation of iceberg removal using pure NumPy
### Justin Hew 

import numpy as np

def killberg_fast(md):
    """
    Remove isolated floating ice patches (icebergs) from md.mask.ice_levelset,
    using an edge-connected flood fill on the ice-only subgraph.

    Requirements on `md` (ISSM-style object):
      md.mesh.elements            : (nE, nVperElem) int array (0- or 1-based)
      md.mesh.numberofelements    : int
      md.mesh.numberofvertices    : int
      md.mask.ice_levelset        : (nV,) float array  (ice if < 0)
      md.mask.ocean_levelset      : (nV,) float array  (grounded if > 0)
    """
    elems = md.mesh.elements
    if elems.min() == 1:    # convert 1-based (MATLAB) to 0-based (Python)
        elems = elems - 1

    nE = md.mesh.numberofelements
    nV = md.mesh.numberofvertices

    ice_ls   = md.mask.ice_levelset
    ocean_ls = md.mask.ocean_levelset

    print("Looking for isolated patches of floating ice (icebergs) [fast]")

    # ---- Element classification ------------------------------------------------
    # Ice elements: at least one vertex with ice (<0)
    is_ice_elem = (ice_ls[elems].min(axis=1) < 0)

    # Grounded elements: >2 vertices with ocean_levelset > 0
    is_grounded_elem = (ocean_ls[elems] > 0).sum(axis=1) > 2

    # Seed vertices = vertices of grounded elements, but only those that are ice vertices
    seed_vertices = np.unique(elems[is_grounded_elem].ravel())
    seed_vertices = seed_vertices[ice_ls[seed_vertices] < 0]

    # If no seeds (no grounded-ice contact), nothing to propagate; everything ice stays as-is
    if seed_vertices.size == 0:
        print("No grounded-ice seeds found; No iceberg removed.")
        return ice_ls.copy()

    # ---- Build edge list on the ICE submesh -----------------------------------
    # Keep only edges from elements that are ice elements
    ice_elems = elems[is_ice_elem]
    if ice_elems.size == 0:
        print("No ice elements found on the mesh; nothing to remove.")
        return ice_ls.copy()

    # Construct undirected edges from polygon/triangle elements:
    # connect consecutive vertices (v0-v1, v1-v2, ..., v_{k-1}-v0)
    k = ice_elems.shape[1]
    v_roll = np.roll(ice_elems, -1, axis=1)
    edges = np.stack([ice_elems.ravel(), v_roll.ravel()], axis=1)

    # Normalize edge endpoints so (a,b) with a<b (to deduplicate)
    a = edges.min(axis=1)
    b = edges.max(axis=1)
    edges = np.stack([a, b], axis=1)

    # Drop self-loops and duplicates
    edges = edges[a != b]
    if edges.size == 0:
        print("Ice submesh has no edges; nothing to remove.")
        return ice_ls.copy()
    edges = np.unique(edges, axis=0)

    # We only want to traverse through ICE vertices
    ice_vertices_mask = (ice_ls < 0)
    # Filter edges to those whose BOTH endpoints are ice vertices
    keep = ice_vertices_mask[edges[:, 0]] & ice_vertices_mask[edges[:, 1]]
    edges = edges[keep]
    if edges.size == 0:
        print("No traversable edges within ice; nothing to remove.")
        return ice_ls.copy()

    # ---- Flood fill (BFS) on the ice-only graph from grounded seeds -----------
    # Induced subgraph on ice vertices: we’ll do a queue-based BFS using numpy sets.
    # Start frontier with the seed ice vertices that exist in the ice-only graph:
    # (i.e., that are incident to at least one kept edge)
    incident = np.zeros(nV, dtype=bool)
    incident[edges[:, 0]] = True
    incident[edges[:, 1]] = True

    frontier = seed_vertices[incident[seed_vertices]]
    if frontier.size == 0:
        # Seeds exist but none touch the ice-only graph — all floating → will be removed
        iceberg_vertices = np.where(ice_ls < 0)[0]
        ice_levelset = ice_ls.copy()
        if iceberg_vertices.size:
            print(f"REMOVING {iceberg_vertices.size} vertices on icebergs")
            ice_levelset[iceberg_vertices] = +1
        else:
            print("No iceberg found!")
        return ice_levelset

    visited = np.zeros(nV, dtype=bool)
    visited[frontier] = True

    # Build adjacency lists from edge list for quick expansion
    # (array-of-lists approach in pure NumPy)
    # For performance, build CSR-like slices:
    deg = np.bincount(edges.ravel(), minlength=nV)
    offsets = np.zeros(nV + 1, dtype=int)
    np.cumsum(deg, out=offsets[1:])
    nbrs = np.empty(edges.shape[0] * 2, dtype=edges.dtype)
    # Fill both directions
    # position cursors
    cursor = offsets.copy()
    for u, v in edges:
        nbrs[cursor[u]] = v; cursor[u] += 1
        nbrs[cursor[v]] = u; cursor[v] += 1

    # BFS
    q = frontier.tolist()
    while q:
        u = q.pop()
        # neighbors slice
        start, end = offsets[u], offsets[u+1]
        neigh = nbrs[start:end]
        # only move to ICE vertices (already ensured by edges filter) and not visited
        next_nodes = neigh[~visited[neigh]]
        if next_nodes.size:
            visited[next_nodes] = True
            q.extend(next_nodes.tolist())

    # visited ∧ ice_vertices_mask are the ice vertices connected (via edges) to grounded ice
    connected_to_grounded = visited & ice_vertices_mask

    # Any remaining ice vertices are disconnected icebergs → set to +1
    iceberg_vertices = np.where(ice_vertices_mask & (~connected_to_grounded))[0]
    ice_levelset = ice_ls.copy()
    if iceberg_vertices.size:
        n_remove = iceberg_vertices.size
        print(f"REMOVING {n_remove} vertex{'es' if n_remove > 1 else ''} on icebergs")
        ice_levelset[iceberg_vertices] = +1
    else:
        print("No iceberg found!")

    return ice_levelset

