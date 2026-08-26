import numpy as np
import sys
import yaml

MAX_GRADIENT = 0.03
R_EARTH = 6371.229


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


def transition_width(h_fine, h_coarse):
    """Compute transition zone width so the gradient does not exceed MAX_GRADIENT."""
    return (h_coarse - h_fine) / MAX_GRADIENT


def _gc_distance_deg(lat1, lon1, lat2, lon2):
    """Great circle distance (km) between two points specified in degrees."""
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)
    x1 = np.cos(lon1) * np.cos(lat1)
    y1 = np.sin(lon1) * np.cos(lat1)
    z1 = np.sin(lat1)
    x2 = np.cos(lon2) * np.cos(lat2)
    y2 = np.sin(lon2) * np.cos(lat2)
    z2 = np.sin(lat2)
    return R_EARTH * np.arccos(np.clip(x1*x2 + y1*y2 + z1*z2, -1.0, 1.0))


def _build_regions(config):
    """Parse refinement regions and compute transition parameters.

    For each region, determines whether it is nested inside a coarser region
    by checking if its center falls within the coarser region's radius.
    Nested regions transition from their resolution to their parent's resolution.
    Non-nested regions transition from their resolution to the coarse resolution.
    """
    coarse_res = config["coarse_resolution"]

    refinements = config.get("refinements", [])
    if not refinements:
        ref = config.get("refinement")
        if ref:
            refinements = [ref]

    if not refinements:
        return [], coarse_res

    # Sort coarsest resolution first
    refinements = sorted(refinements, key=lambda r: -r["resolution"])

    regions = []
    for ref in refinements:
        regions.append({
            "resolution": ref["resolution"],
            "center_latitude": ref["center_latitude"],
            "center_longitude": ref["center_longitude"],
            "radius": ref["radius"],
            "nested": False,
            "transition_target": coarse_res,
            "transition_width": transition_width(ref["resolution"], coarse_res),
        })

    # A region is nested if its center is inside a coarser region's radius.
    # Among all such candidates, pick the finest (tightest parent).
    for i in range(len(regions)):
        ri = regions[i]
        best_parent_idx = None
        best_parent_res = None

        for j in range(len(regions)):
            if i == j:
                continue
            rj = regions[j]
            if rj["resolution"] <= ri["resolution"]:
                continue

            dist = _gc_distance_deg(
                ri["center_latitude"], ri["center_longitude"],
                rj["center_latitude"], rj["center_longitude"],
            )

            if dist < rj["radius"]:
                if best_parent_idx is None or rj["resolution"] < best_parent_res:
                    best_parent_idx = j
                    best_parent_res = rj["resolution"]

        if best_parent_idx is not None:
            ri["nested"] = True
            ri["transition_target"] = regions[best_parent_idx]["resolution"]
            ri["transition_width"] = transition_width(
                ri["resolution"], ri["transition_target"]
            )

    return regions, coarse_res


def _validate_regions(regions):
    """Validate that refinement region transition zones do not overlap.

    For nested regions, the inner region plus its transition zone must be fully
    contained within the outer region's uniform area.  For non-nested regions,
    their total footprints (radius + transition) must not overlap.
    """
    n = len(regions)
    for i in range(n):
        for j in range(i + 1, n):
            ri, rj = regions[i], regions[j]

            dist = _gc_distance_deg(
                ri["center_latitude"], ri["center_longitude"],
                rj["center_latitude"], rj["center_longitude"],
            )

            ri_outer = ri["radius"] + ri["transition_width"]
            rj_outer = rj["radius"] + rj["transition_width"]

            # Determine if one is nested in the other
            ri_nested_in_rj = ri["nested"] and dist < rj["radius"]
            rj_nested_in_ri = rj["nested"] and dist < ri["radius"]

            if ri_nested_in_rj:
                if dist + ri_outer > rj["radius"]:
                    sys.exit(
                        f"Error: Nested refinement region does not fit within"
                        f" its parent.\n"
                        f"  Inner region at"
                        f" ({ri['center_latitude']}, {ri['center_longitude']})"
                        f" with radius {ri['radius']} km"
                        f" + transition {ri['transition_width']:.1f} km\n"
                        f"  does not fit within outer region at"
                        f" ({rj['center_latitude']}, {rj['center_longitude']})"
                        f" with radius {rj['radius']} km.\n"
                        f"  Inner footprint extends to {dist + ri_outer:.1f} km"
                        f" from outer center,"
                        f" but must be within {rj['radius']:.1f} km."
                    )
                continue

            if rj_nested_in_ri:
                if dist + rj_outer > ri["radius"]:
                    sys.exit(
                        f"Error: Nested refinement region does not fit within"
                        f" its parent.\n"
                        f"  Inner region at"
                        f" ({rj['center_latitude']}, {rj['center_longitude']})"
                        f" with radius {rj['radius']} km"
                        f" + transition {rj['transition_width']:.1f} km\n"
                        f"  does not fit within outer region at"
                        f" ({ri['center_latitude']}, {ri['center_longitude']})"
                        f" with radius {ri['radius']} km.\n"
                        f"  Inner footprint extends to {dist + rj_outer:.1f} km"
                        f" from outer center,"
                        f" but must be within {ri['radius']:.1f} km."
                    )
                continue

            # Not nested — footprints must not overlap
            if dist < ri_outer + rj_outer:
                sys.exit(
                    f"Error: Refinement regions are too close together.\n"
                    f"  Region at ({ri['center_latitude']}, {ri['center_longitude']})"
                    f" with radius {ri['radius']} km"
                    f" + transition {ri['transition_width']:.1f} km\n"
                    f"  and region at ({rj['center_latitude']}, {rj['center_longitude']})"
                    f" with radius {rj['radius']} km"
                    f" + transition {rj['transition_width']:.1f} km\n"
                    f"  are separated by {dist:.1f} km"
                    f" but require at least {ri_outer + rj_outer:.1f} km."
                )


def _h_single(r, resolution, transition_target, radius, tw, nested):
    """Compute h for one region.

    Non-nested regions return transition_target outside their footprint.
    Nested regions return inf outside their footprint.
    """
    if nested:
        ret = np.full_like(r, np.inf)
    else:
        ret = np.full_like(r, transition_target)

    ret[r <= radius] = resolution

    if tw > 0:
        mask = (r > radius) & (r < radius + tw)
        ret[mask] = (
            resolution + (r[mask] - radius) * (transition_target - resolution) / tw
        )

    return ret


def get_hfun(longitude, latitude, config_file="mesh.yaml"):
    config = load_config(config_file)
    regions, coarse_res = _build_regions(config)
    _validate_regions(regions)

    shape = longitude.shape
    h_values = np.full(shape, coarse_res, dtype=float)

    for reg in regions:
        x_c, y_c, z_c = geo_to_cartesian(
            np.radians(reg["center_longitude"]),
            np.radians(reg["center_latitude"]),
        )
        p_center = np.array([x_c, y_c, z_c])

        x, y, z = geo_to_cartesian(longitude, latitude)
        p = np.column_stack((x.flatten(), y.flatten(), z.flatten()))

        r = R_EARTH * unit_sphere_distance(p_center, p)

        h_region = _h_single(
            r,
            reg["resolution"],
            reg["transition_target"],
            reg["radius"],
            reg["transition_width"],
            reg["nested"],
        )

        h_values = np.minimum(h_values.flatten(), h_region).reshape(shape)

    return h_values


def geo_to_cartesian(lam, phi):
    x = np.cos(lam) * np.cos(phi)
    y = np.sin(lam) * np.cos(phi)
    z = np.sin(phi)
    return (x, y, z)


def unit_sphere_distance(p, q_arr):
    return np.arccos(np.clip(q_arr @ p, -1.0, 1.0))
