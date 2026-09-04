import numpy as np
import sys
import yaml

MAX_GRADIENT = 0.03
R_EARTH = 6371.229
_N_SAMPLES = 360


def load_config(config_file="mesh.yaml"):
    with open(config_file, "r") as f:
        return yaml.safe_load(f)["mesh"]


def finest_resolution(config_file="mesh.yaml"):
    """Return the finest (smallest) resolution across all refinement regions."""
    config = load_config(config_file)
    regions, coarse_res = _build_regions(config)
    if not regions:
        return coarse_res
    return min(r["resolution"] for r in regions)


def _profile_max_slope(profile):
    """Return the maximum normalized slope of a transition profile."""
    if profile == "linear":
        return 1.0
    if profile == "smoothstep":
        return 30.0 / 16.0
    sys.exit(f"Error: Unknown transition profile '{profile}'")


def transition_width(h_fine, h_coarse, scale=1.0, profile="linear"):
    """Compute transition zone width for the given profile."""
    return scale * (h_coarse - h_fine) * _profile_max_slope(profile) / MAX_GRADIENT


# ── Tangent plane projection ──────────────────────────────────────────────────


def _small_circle(clon_rad, clat_rad, radius_km, angles):
    """Return (lons_rad, lats_rad) for a small circle on the sphere."""
    alpha = radius_km / R_EARTH
    sin_c, cos_c = np.sin(clat_rad), np.cos(clat_rad)
    sin_a, cos_a = np.sin(alpha), np.cos(alpha)
    lats = np.arcsin(np.clip(sin_c * cos_a + cos_c * sin_a * np.cos(angles), -1, 1))
    lons = clon_rad + np.arctan2(
        np.sin(angles) * sin_a * cos_c,
        cos_a - sin_c * np.sin(lats)
    )
    return lons, lats


def _bearing(clon, clat, lons, lats):
    """Compute initial bearing (azimuth) from center to each point, in radians."""
    dlon = lons - clon
    x = np.sin(dlon) * np.cos(lats)
    y = np.cos(clat) * np.sin(lats) - np.sin(clat) * np.cos(lats) * np.cos(dlon)
    return np.arctan2(x, y)


def _ellipse_small_circle(clon, clat, a, b, orient, angles):
    """Return (lons_rad, lats_rad) for a spherical ellipse at given azimuths."""
    r_bnd = a * b / np.sqrt(
        (b * np.cos(angles - orient)) ** 2 + (a * np.sin(angles - orient)) ** 2
    )
    return _small_circle(clon, clat, r_bnd, angles)


# ── Signed distance functions ─────────────────────────────────────────────────

def _signed_distance(shape, lons_rad, lats_rad, max_dist=None):
    """Signed distance (km) from shape boundary: negative inside, positive outside."""
    stype = shape["type"]
    if stype == "circle":
        return _sd_circle(shape, lons_rad, lats_rad)
    elif stype == "ellipse":
        return _sd_ellipse(shape, lons_rad, lats_rad, max_dist)
    elif stype == "polygon":
        return _sd_polygon(shape, lons_rad, lats_rad, max_dist)
    else:
        sys.exit(f"Error: Unknown shape type '{stype}'")


def _sd_circle(shape, lons_rad, lats_rad):
    clat = np.radians(shape["center_latitude"])
    clon = np.radians(shape["center_longitude"])
    radius = shape["radius"]

    x_c, y_c, z_c = geo_to_cartesian(clon, clat)
    p_center = np.array([x_c, y_c, z_c])

    x, y, z = geo_to_cartesian(lons_rad.flatten(), lats_rad.flatten())
    p = np.column_stack((x, y, z))
    r = R_EARTH * unit_sphere_distance(p_center, p)

    return (r - radius).reshape(lons_rad.shape)


def _sd_ellipse(shape, lons_rad, lats_rad, max_dist=None):
    """Signed distance to spherical ellipse via Newton nearest-point."""
    clon = np.radians(shape["center_longitude"])
    clat = np.radians(shape["center_latitude"])
    a = shape["semi_major"]
    b = shape["semi_minor"]
    orient = np.radians(shape.get("orientation", 0.0))

    flat_lons = lons_rad.flatten()
    flat_lats = lats_rad.flatten()

    # Angular distance and bearing from center to each point
    x_c, y_c, z_c = geo_to_cartesian(clon, clat)
    px, py, pz = geo_to_cartesian(flat_lons, flat_lats)
    P = np.column_stack((px, py, pz))
    n_pts = len(P)
    p_center = np.array([x_c, y_c, z_c])
    delta = np.arccos(np.clip(P @ p_center, -1.0, 1.0))
    beta = _bearing(clon, clat, flat_lons, flat_lats)

    # Pre-filter: skip points too far from ellipse to matter
    if max_dist is not None:
        cutoff_ang = a / R_EARTH + max_dist / R_EARTH + 0.01
        candidate = delta < cutoff_ang
        P_full, n_pts_full = P, n_pts
        delta = delta[candidate]
        beta = beta[candidate]
        flat_lons = flat_lons[candidate]
        flat_lats = flat_lats[candidate]
        P = P[candidate]
        n_pts = len(P)

    # Inside test: radial distance vs boundary radius at point's bearing
    alpha_p = beta - orient
    r_at_beta = a * b / np.sqrt((b * np.cos(alpha_p))**2 + (a * np.sin(alpha_p))**2)
    inside = (R_EARTH * delta) <= r_at_beta

    # Coarse sweep to find best Newton start (avoids wrong-extremum convergence)
    a2, b2, ab = a**2, b**2, a * b
    diff2 = a2 - b2
    n_sweep = 36
    sweep_angles = np.linspace(0, 2 * np.pi, n_sweep, endpoint=False)
    sweep_D = b2 * np.cos(sweep_angles - orient)**2 + a2 * np.sin(sweep_angles - orient)**2
    sweep_lons, sweep_lats = _small_circle(clon, clat, ab / np.sqrt(sweep_D), sweep_angles)
    sx, sy, sz = geo_to_cartesian(sweep_lons, sweep_lats)
    S = np.column_stack((sx, sy, sz))
    best = np.argmax(P @ S.T, axis=1)
    theta = sweep_angles[best]

    cos_d, sin_d = np.cos(delta), np.sin(delta)

    for _ in range(8):
        al = theta - orient
        cos_al, sin_al = np.cos(al), np.sin(al)
        D = b2 * cos_al**2 + a2 * sin_al**2
        sqrt_D = np.sqrt(D)
        rho = ab / (sqrt_D * R_EARTH)

        Dp = 2 * diff2 * sin_al * cos_al
        rhop = -ab * Dp / (2 * D * sqrt_D * R_EARTH)

        Dpp = 2 * diff2 * (cos_al**2 - sin_al**2)
        rhopp = -ab * (Dpp * D - 1.5 * Dp**2) / (2 * D**2 * sqrt_D * R_EARTH)

        cos_r, sin_r = np.cos(rho), np.sin(rho)
        gamma = beta - theta
        cos_g, sin_g = np.cos(gamma), np.sin(gamma)

        F = cos_d * cos_r + sin_d * sin_r * cos_g
        A_val = sin_d * cos_r * cos_g - cos_d * sin_r
        g = rhop * A_val + sin_d * sin_r * sin_g
        gp = (rhopp * A_val - rhop**2 * F
              + 2 * rhop * sin_d * cos_r * sin_g
              - sin_d * sin_r * cos_g)
        gp = np.where(np.abs(gp) < 1e-12, np.copysign(1e-12, gp), gp)
        theta -= g / gp

    # Final GC distance to nearest boundary point
    al = theta - orient
    D = b2 * np.cos(al)**2 + a2 * np.sin(al)**2
    rho = ab / (np.sqrt(D) * R_EARTH)
    gamma = beta - theta
    cos_pb = cos_d * np.cos(rho) + sin_d * np.sin(rho) * np.cos(gamma)
    gc_dist = R_EARTH * np.arccos(np.clip(cos_pb, -1.0, 1.0))

    sd = np.where(inside, -gc_dist, gc_dist)

    if max_dist is not None:
        result = np.full(n_pts_full, max_dist * 1.1)
        result[candidate] = sd
        return result.reshape(lons_rad.shape)
    return sd.reshape(lons_rad.shape)

def _sd_polygon(shape, lons_rad, lats_rad, max_dist=None):
    """Signed distance to convex polygon using great-circle arcs on the sphere."""
    vertices = np.array(shape["vertices"])
    vlats = np.radians(vertices[:, 0])
    vlons = np.radians(vertices[:, 1])

    vx, vy, vz = geo_to_cartesian(vlons, vlats)
    V = np.column_stack((vx, vy, vz))
    nv = len(V)

    # Precompute edge normals (depend only on polygon vertices)
    V_next = np.roll(V, -1, axis=0)
    edge_normals = np.cross(V, V_next)
    edge_lens = np.linalg.norm(edge_normals, axis=1, keepdims=True)
    edge_normals /= np.maximum(edge_lens, 1e-15)

    # Centroid for inside test and pre-filter
    centroid = V.mean(axis=0)
    centroid /= np.linalg.norm(centroid)

    px, py, pz = geo_to_cartesian(lons_rad.flatten(), lats_rad.flatten())
    P = np.column_stack((px, py, pz))
    n_pts = len(P)

    # Pre-filter: skip points too far from polygon to matter
    if max_dist is not None:
        # Check edge-arc extrema (arcs can bulge beyond vertices)
        poly_radius = np.max(np.arccos(np.clip(V @ centroid, -1.0, 1.0)))
        for i in range(nv):
            j = (i + 1) % nv
            omega = np.arccos(np.clip(V[i] @ V[j], -1.0, 1.0))
            if omega < 1e-12:
                continue
            a = centroid @ V[i]
            b = centroid @ V[j]
            phi = np.arctan2((b - a * np.cos(omega)) / np.sin(omega), a)
            t_ext = (phi + np.pi) / omega
            if 0 < t_ext < 1:
                so = np.sin(omega)
                p = (np.sin((1 - t_ext) * omega) * V[i]
                     + np.sin(t_ext * omega) * V[j]) / so
                poly_radius = max(poly_radius, np.arccos(np.clip(centroid @ p, -1, 1)))
        cutoff_ang = poly_radius + max_dist / R_EARTH + 0.01
        cos_cutoff = np.cos(cutoff_ang)
        cos_cent = P @ centroid
        candidate = cos_cent > cos_cutoff
        P_full, n_pts_full = P, n_pts
        P = P[candidate]
        n_pts = len(P)

    # Pass 1: unsigned distance to nearest edge (cosine domain, single arccos at end)
    max_cos = np.full(n_pts, -1.0)
    best_edge = np.zeros(n_pts, dtype=np.intp)

    for i in range(nv):
        j = (i + 1) % nv
        A, B = V[i], V[j]

        if edge_lens[i, 0] < 1e-15:
            continue
        edge_n = edge_normals[i]

        cos_dA = P @ A
        cos_dB = P @ B
        cos_dist_edge = np.maximum(cos_dA, cos_dB)

        side = P @ edge_n
        proj = P - side[:, np.newaxis] * edge_n
        proj_norms = np.linalg.norm(proj, axis=1, keepdims=True)
        safe = (proj_norms > 1e-15).flatten()
        proj_n = np.divide(proj, proj_norms, where=safe[:, np.newaxis],
                           out=np.zeros_like(proj))

        cross_AP = np.cross(A, proj_n)
        cross_PB = np.cross(proj_n, B)
        on_arc = safe & (cross_AP @ edge_n >= -1e-10) & (cross_PB @ edge_n >= -1e-10)

        cos_dgc = np.sqrt(np.clip(1.0 - side**2, 0.0, 1.0))
        cos_dist_edge[on_arc] = cos_dgc[on_arc]

        closer = cos_dist_edge > max_cos
        best_edge[closer] = i
        max_cos = np.maximum(max_cos, cos_dist_edge)

    min_dist = np.arccos(np.clip(max_cos, -1.0, 1.0))

    # Pass 2: nearest-edge sign test (no projection or hemisphere limit)
    max_ang = (max_dist / R_EARTH + 0.01) if max_dist else np.pi
    idx = np.where(min_dist < max_ang)[0]

    inside = np.zeros(n_pts, dtype=bool)

    if len(idx) > 0:
        # Orient edge normals inward (toward centroid)
        centroid_signs = np.sign(centroid @ edge_normals.T)
        inward_n = edge_normals * centroid_signs[:, np.newaxis]
        # Point is inside if it's on the inward side of its nearest edge
        be = best_edge[idx]
        side_nearest = np.sum(P[idx] * inward_n[be], axis=1)
        inside[idx] = side_nearest > -1e-10

    min_dist_km = R_EARTH * min_dist
    sd = np.where(inside, -min_dist_km, min_dist_km)

    if max_dist is not None:
        result = np.full(n_pts_full, max_dist * 1.1)
        result[candidate] = sd
        return result.reshape(lons_rad.shape)
    return sd.reshape(lons_rad.shape)



# ── Shape parsing ─────────────────────────────────────────────────────────────

def _parse_shape(ref):
    stype = ref.get("shape", "circle")
    shape = {"type": stype}

    if stype == "circle":
        shape["center_latitude"] = ref["center_latitude"]
        shape["center_longitude"] = ref["center_longitude"]
        shape["radius"] = ref["radius"]
    elif stype == "ellipse":
        shape["center_latitude"] = ref["center_latitude"]
        shape["center_longitude"] = ref["center_longitude"]
        shape["semi_major"] = ref["semi_major"]
        shape["semi_minor"] = ref["semi_minor"]
        shape["orientation"] = ref.get("orientation", 0.0)
    elif stype == "polygon":
        shape["vertices"] = ref["vertices"]
        _check_convex(ref["vertices"])
    else:
        sys.exit(f"Error: Unknown shape type '{stype}'")

    return shape


def _check_convex(vertices_deg):
    """Exit with error if polygon vertices are not convex (spherical geometry)."""
    verts = np.array(vertices_deg)
    n = len(verts)
    if n < 3:
        sys.exit("Error: Polygon must have at least 3 vertices.")

    vlats = np.radians(verts[:, 0])
    vlons = np.radians(verts[:, 1])
    vx, vy, vz = geo_to_cartesian(vlons, vlats)
    V = np.column_stack((vx, vy, vz))

    sign = None
    centroid = V.mean(axis=0)
    c_len = np.linalg.norm(centroid)
    if c_len < 1e-15:
        sys.exit("Error: Polygon vertices cancel out (antipodal or degenerate).")
    centroid /= c_len
    for i in range(n):
        j = (i + 1) % n
        edge_n = np.cross(V[i], V[j])
        if np.linalg.norm(edge_n) < 1e-15:
            continue
        s = centroid @ edge_n
        if abs(s) < 1e-15:
            continue
        if sign is None:
            sign = s > 0
        elif (s > 0) != sign:
            sys.exit("Error: Polygon vertices do not form a convex polygon.")
    if sign is None:
        sys.exit("Error: Polygon has no valid edges (coincident or degenerate vertices).")


def _shape_center(shape):
    """Return (latitude, longitude) in degrees for a shape's center."""
    if shape["type"] in ("circle", "ellipse"):
        return shape["center_latitude"], shape["center_longitude"]
    vertices = np.array(shape["vertices"])
    return float(np.mean(vertices[:, 0])), float(np.mean(vertices[:, 1]))


# ── Region building ───────────────────────────────────────────────────────────

def _build_regions(config):
    coarse_res = config["resolution"]

    refinements = config.get("refinements", [])
    if not refinements:
        return [], coarse_res

    scale = config.get("transition_scale", 1.0)
    profile = config.get("transition_profile", "linear")

    regions = []
    _flatten_tree(refinements, coarse_res, coarse_res, regions, None, scale,
                  profile)

    return regions, coarse_res


def _flatten_tree(refinements, coarse_res, parent_res, regions, parent_idx,
                  scale, profile):
    """Recursively flatten a hierarchical refinement tree into a flat list."""
    for ref in refinements:
        shape = _parse_shape(ref)
        is_nested = parent_idx is not None
        target = parent_res
        tw = transition_width(ref["resolution"], target, scale, profile)

        idx = len(regions)
        regions.append({
            "resolution": ref["resolution"],
            "shape": shape,
            "nested": is_nested,
            "parent_index": parent_idx,
            "transition_target": target,
            "transition_width": tw,
            "profile": profile,
        })

        children = ref.get("refinements", [])
        if children:
            _flatten_tree(children, coarse_res, ref["resolution"], regions, idx,
                          scale, profile)


# ── Validation ────────────────────────────────────────────────────────────────

def _polygon_boundary(vlons_rad, vlats_rad, n_samples=_N_SAMPLES):
    """Sample a closed polygon boundary with great-circle arc interpolation."""
    vx, vy, vz = geo_to_cartesian(vlons_rad, vlats_rad)
    V = np.column_stack((vx, vy, vz))
    nv = len(V)
    spe = max(10, n_samples // nv)
    pts = []
    for i in range(nv):
        j = (i + 1) % nv
        A, B = V[i], V[j]
        omega = np.arccos(np.clip(A @ B, -1, 1))
        if omega < 1e-10:
            pts.append(A)
            continue
        sin_omega = np.sin(omega)
        for t in np.linspace(0, 1, spe, endpoint=False):
            P = (np.sin((1-t)*omega) * A + np.sin(t*omega) * B) / sin_omega
            pts.append(P)
    pts.append(pts[0].copy())  # close the loop
    pts = np.array(pts)
    lons = np.arctan2(pts[:, 1], pts[:, 0])
    lats = np.arcsin(np.clip(pts[:, 2], -1, 1))
    return lons, lats


def _offset_boundary(bnd_lons, bnd_lats, center_lon, center_lat, offset_km):
    """Offset a closed boundary curve by offset_km in the outward normal direction."""
    bx, by, bz = geo_to_cartesian(bnd_lons, bnd_lats)
    B = np.column_stack((bx, by, bz))
    n = len(B)

    # Central finite difference for tangent (periodic)
    closed = np.linalg.norm(B[0] - B[-1]) < 1e-10
    T = np.empty_like(B)
    if closed:
        T[0] = B[1] - B[-2]
        T[-1] = T[0].copy()
    else:
        T[0] = B[1] - B[0]
        T[-1] = B[-1] - B[-2]
    if n > 2:
        T[1:-1] = B[2:] - B[:-2]

    # Normal in tangent plane: N = B x T, normalized
    N = np.cross(B, T)
    norms = np.linalg.norm(N, axis=1, keepdims=True)
    N /= np.where(norms < 1e-15, 1, norms)

    # Ensure outward (away from center)
    cx, cy, cz = geo_to_cartesian(center_lon, center_lat)
    C = np.array([np.asarray(cx).item(), np.asarray(cy).item(), np.asarray(cz).item()])
    dots_C = B @ C
    away = B * dots_C[:, np.newaxis] - C[np.newaxis, :]
    flip = np.sum(N * away, axis=1) < 0
    N[flip] *= -1

    alpha = offset_km / R_EARTH
    P = np.cos(alpha) * B + np.sin(alpha) * N
    P /= np.linalg.norm(P, axis=1, keepdims=True)

    return np.arctan2(P[:, 1], P[:, 0]), np.arcsin(np.clip(P[:, 2], -1, 1))


def _sample_outer_boundary(shape, tw, n_samples=_N_SAMPLES):
    """Sample points on the outer edge of a shape's transition zone."""
    stype = shape["type"]
    angles = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)

    if stype == "circle":
        clat = np.radians(shape["center_latitude"])
        clon = np.radians(shape["center_longitude"])
        return _small_circle(clon, clat, shape["radius"] + tw, angles)

    elif stype == "ellipse":
        clat = np.radians(shape["center_latitude"])
        clon = np.radians(shape["center_longitude"])
        a = shape["semi_major"]
        b = shape["semi_minor"]
        orient = np.radians(shape.get("orientation", 0.0))
        r_bnd = a * b / np.sqrt(
            (b * np.cos(angles - orient)) ** 2 + (a * np.sin(angles - orient)) ** 2
        )
        bnd_lons, bnd_lats = _small_circle(clon, clat, r_bnd, angles)
        return _offset_boundary(bnd_lons, bnd_lats, clon, clat, tw)

    elif stype == "polygon":
        vertices = np.array(shape["vertices"])
        vlats = np.radians(vertices[:, 0])
        vlons = np.radians(vertices[:, 1])
        vx, vy, vz = geo_to_cartesian(vlons, vlats)
        V = np.column_stack((vx, vy, vz))
        nv = len(V)

        centroid = V.mean(axis=0)
        centroid /= np.linalg.norm(centroid)

        alpha = tw / R_EARTH
        sin_a, cos_a = np.sin(alpha), np.cos(alpha)

        outward_normals = []
        for i in range(nv):
            j = (i + 1) % nv
            n = np.cross(V[i], V[j])
            n_len = np.linalg.norm(n)
            if n_len > 1e-15:
                n /= n_len
            outward_normals.append(-np.sign(centroid @ n) * n)

        pts = []
        spe = max(5, n_samples // (2 * nv))
        spc = max(3, n_samples // (2 * nv))

        for i in range(nv):
            j = (i + 1) % nv
            A, B = V[i], V[j]
            outward = outward_normals[i]

            omega = np.arccos(np.clip(A @ B, -1, 1))
            sin_omega = np.sin(omega) if omega > 1e-10 else 1.0
            for t in np.linspace(0, 1, spe, endpoint=False):
                if omega < 1e-10:
                    Q = A.copy()
                else:
                    Q = (np.sin((1-t)*omega) * A + np.sin(t*omega) * B) / sin_omega
                Q /= np.linalg.norm(Q)
                P = cos_a * Q + sin_a * outward
                P /= np.linalg.norm(P)
                pts.append(P)

            # Corner arc at vertex j
            outward_next = outward_normals[j]
            u1 = outward - (outward @ V[j]) * V[j]
            u1_len = np.linalg.norm(u1)
            u2 = outward_next - (outward_next @ V[j]) * V[j]
            u2_len = np.linalg.norm(u2)
            if u1_len < 1e-15 or u2_len < 1e-15:
                continue
            u1 /= u1_len
            u2 /= u2_len

            ang = np.arccos(np.clip(u1 @ u2, -1, 1))
            if np.cross(u1, u2) @ V[j] < 0:
                ang = -ang

            if ang >= 0:
                # Concave vertex on sphere: no corner arc needed.
                # Just bridge to the next edge's offset direction.
                P = cos_a * V[j] + sin_a * u2
                P /= np.linalg.norm(P)
                pts.append(P)
            else:
                # Convex vertex: sweep a rounded corner arc
                vj_cross_u1 = np.cross(V[j], u1)
                for s in np.linspace(0, ang, spc, endpoint=True):
                    u_rot = u1 * np.cos(s) + vj_cross_u1 * np.sin(s)
                    P = cos_a * V[j] + sin_a * u_rot
                    P /= np.linalg.norm(P)
                    pts.append(P)

        pts = np.array(pts)
        lons = np.arctan2(pts[:, 1], pts[:, 0])
        lats = np.arcsin(np.clip(pts[:, 2], -1, 1))
        return lons, lats


def _is_ancestor(regions, ancestor_idx, descendant_idx):
    """Check if ancestor_idx is an ancestor of descendant_idx in the tree."""
    current = descendant_idx
    while regions[current]["parent_index"] is not None:
        current = regions[current]["parent_index"]
        if current == ancestor_idx:
            return True
    return False


def _minimum_width(shape):
    """Return the minimum width (km) across the narrowest point of a shape."""
    stype = shape["type"]
    if stype == "circle":
        return 2.0 * shape["radius"]
    elif stype == "ellipse":
        return 2.0 * min(shape["semi_major"], shape["semi_minor"])
    elif stype == "polygon":
        vertices = np.array(shape["vertices"])
        vlats = np.radians(vertices[:, 0])
        vlons = np.radians(vertices[:, 1])
        vx, vy, vz = geo_to_cartesian(vlons, vlats)
        V = np.column_stack((vx, vy, vz))
        nv = len(V)
        min_w = np.inf
        for i in range(nv):
            j = (i + 1) % nv
            n = np.cross(V[i], V[j])
            n_len = np.linalg.norm(n)
            if n_len < 1e-15:
                continue
            n /= n_len
            dists = np.abs(V @ n)
            dists[i] = 0.0
            dists[j] = 0.0
            w = R_EARTH * np.arcsin(np.clip(dists.max(), 0, 1))
            min_w = min(min_w, w)
        return min_w


def _validate_regions(regions):
    """Validate containment for nested regions and separation for non-nested."""
    n = len(regions)

    min_cells = 100
    for reg in regions:
        width = _minimum_width(reg["shape"])
        n_cells = width / reg["resolution"]
        if n_cells < min_cells:
            clat, clon = _shape_center(reg["shape"])
            sys.exit(
                f"Error: Refinement region is too narrow.\n"
                f"  Region ({reg['resolution']} km) at ({clat}, {clon})\n"
                f"  is only {n_cells:.0f} cells wide at its narrowest point.\n"
                f"  Minimum required: {min_cells} cells"
                f" ({min_cells * reg['resolution']:.0f} km)."
            )

    for i in range(n):
        for j in range(i + 1, n):
            ri, rj = regions[i], regions[j]

            # Direct parent-child: validate containment
            if ri["parent_index"] == j or rj["parent_index"] == i:
                if ri["parent_index"] == j:
                    child, parent = ri, rj
                else:
                    child, parent = rj, ri

                outer_lons, outer_lats = _sample_outer_boundary(
                    child["shape"], child["transition_width"]
                )
                parent_sd = _signed_distance(
                    parent["shape"], outer_lons, outer_lats
                ).flatten()

                if np.any(parent_sd > 0):
                    clat, clon = _shape_center(child["shape"])
                    plat, plon = _shape_center(parent["shape"])
                    sys.exit(
                        f"Error: Nested refinement region does not fit"
                        f" within its parent.\n"
                        f"  Inner region ({child['resolution']} km)"
                        f" at ({clat}, {clon})"
                        f" with transition"
                        f" {child['transition_width']:.1f} km\n"
                        f"  does not fit within outer region"
                        f" ({parent['resolution']} km)"
                        f" at ({plat}, {plon})."
                    )
                continue

            # Skip indirect ancestor-descendant pairs
            if _is_ancestor(regions, i, j) or _is_ancestor(regions, j, i):
                continue

            # Non-related regions: check for overlap
            outer_i_lons, outer_i_lats = _sample_outer_boundary(
                ri["shape"], ri["transition_width"]
            )
            outer_j_lons, outer_j_lats = _sample_outer_boundary(
                rj["shape"], rj["transition_width"]
            )

            sd_i_at_j = _signed_distance(
                rj["shape"], outer_i_lons, outer_i_lats
            ).flatten()
            sd_j_at_i = _signed_distance(
                ri["shape"], outer_j_lons, outer_j_lats
            ).flatten()

            overlap = (np.any(sd_i_at_j < rj["transition_width"]) or
                       np.any(sd_j_at_i < ri["transition_width"]))
            if overlap:
                clat_i, clon_i = _shape_center(ri["shape"])
                clat_j, clon_j = _shape_center(rj["shape"])
                sys.exit(
                    f"Error: Refinement regions are too close together.\n"
                    f"  Region ({ri['resolution']} km)"
                    f" at ({clat_i}, {clon_i})"
                    f" and region ({rj['resolution']} km)"
                    f" at ({clat_j}, {clon_j})"
                    f" have overlapping transition zones."
                )


# ── h function ────────────────────────────────────────────────────────────────

def _h_single(d, resolution, transition_target, tw, nested, profile="linear"):
    """Compute h from signed distance d to shape boundary."""
    if nested:
        ret = np.full_like(d, np.inf)
    else:
        ret = np.full_like(d, transition_target)

    ret[d <= 0] = resolution

    if tw > 0:
        mask = (d > 0) & (d < tw)
        if profile == "linear":
            ret[mask] = resolution + (transition_target - resolution) * d[mask] / tw
        elif profile == "smoothstep":
            t = d[mask] / tw
            ret[mask] = resolution + (transition_target - resolution) * (
                6.0 * t ** 5 - 15.0 * t ** 4 + 10.0 * t ** 3)

    return ret


def get_hfun(longitude, latitude, config_file="mesh.yaml"):
    config = load_config(config_file)
    regions, coarse_res = _build_regions(config)
    _validate_regions(regions)

    shape = longitude.shape
    h_values = np.full(shape, coarse_res, dtype=float)

    for reg in regions:
        d = _signed_distance(reg["shape"], longitude, latitude,
                             max_dist=reg["transition_width"]).flatten()

        h_region = _h_single(
            d,
            reg["resolution"],
            reg["transition_target"],
            reg["transition_width"],
            reg["nested"],
            reg["profile"],
        )

        h_values = np.minimum(h_values.flatten(), h_region).reshape(shape)

    return h_values


# ── Utilities ─────────────────────────────────────────────────────────────────

def geo_to_cartesian(lam, phi):
    x = np.cos(lam) * np.cos(phi)
    y = np.sin(lam) * np.cos(phi)
    z = np.sin(phi)
    return (x, y, z)


def unit_sphere_distance(p, q_arr):
    return np.arccos(np.clip(q_arr @ p, -1.0, 1.0))


def get_region_boundaries(config_file="mesh.yaml", **kwargs):
    """Extract boundary contours from hfun values on a regular grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config = load_config(config_file)
    regions, base_res = _build_regions(config)
    if not regions:
        return []

    # Compute hfun on a regular grid
    lons_deg = np.linspace(-180, 180, 720)
    lats_deg = np.linspace(-90, 90, 360)
    longrid, latgrid = np.meshgrid(np.radians(lons_deg), np.radians(lats_deg))
    hvals = get_hfun(longrid, latgrid, config_file)

    # Build (contour_level, label, kind) tuples; offset to separate shared levels
    specs = []
    seen = set()
    for reg in regions:
        res, tgt = reg["resolution"], reg["transition_target"]
        if ("inner", res) not in seen:
            seen.add(("inner", res))
            specs.append((res + 0.01, res, "boundary"))
        if ("outer", tgt) not in seen:
            seen.add(("outer", tgt))
            clevel = tgt - 0.5 if tgt >= base_res else tgt - 0.01
            specs.append((clevel, tgt, "transition"))

    boundaries = []
    fig, ax = plt.subplots()
    for clevel, label, kind in specs:
        cs = ax.contour(lons_deg, lats_deg, hvals, levels=[clevel])
        for seg in cs.allsegs[0]:
            if len(seg) >= 10:
                boundaries.append({
                    "lons_deg": seg[:, 0],
                    "lats_deg": seg[:, 1],
                    "resolution": label,
                    "kind": kind,
                })
        ax.cla()
    plt.close(fig)
    return boundaries
