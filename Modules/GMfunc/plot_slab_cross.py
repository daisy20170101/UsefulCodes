import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------ core geometry helpers ------------------
def _project_aeqd(lon, lat, lon0, lat0):
    """
    Project lon/lat to a local Azimuthal Equidistant system centered at (lon0, lat0).
    Uses pyproj if available; otherwise falls back to a spherical approximation.
    Returns x, y in meters.
    """
    try:
        from pyproj import CRS, Transformer

        aeqd = CRS.from_proj4(
            f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +x_0=0 +y_0=0 +units=m +no_defs"
        )
        tfm = Transformer.from_crs("EPSG:4326", aeqd, always_xy=True)
        x, y = tfm.transform(lon, lat)
        return np.asarray(x, float), np.asarray(y, float)
    except Exception:
        # Fallback: equirectangular about (lon0, lat0), Earth radius 6371 km
        R = 6371000.0
        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        lon0 = float(lon0)
        lat0 = float(lat0)
        lam = np.deg2rad(lon)
        phi = np.deg2rad(lat)
        lam0 = np.deg2rad(lon0)
        phi0 = np.deg2rad(lat0)
        x = R * (lam - lam0) * np.cos(phi0)
        y = R * (phi - phi0)
        return x, y


def _profile_axes(lon1, lat1, lon2, lat2, lon_pts, lat_pts):
    """
    Build along-track and cross-track coordinates (meters) for points lon_pts/lat_pts
    relative to the segment P1->P2 defined by (lon1,lat1) and (lon2,lat2).
    """
    lon1, lat1, lon2, lat2 = map(float, (lon1, lat1, lon2, lat2))
    # center for projection
    lon0 = 0.5 * (lon1 + lon2)
    lat0 = 0.5 * (lat1 + lat2)

    # project to local metric coordinates
    x1, y1 = _project_aeqd(lon1, lat1, lon0, lat0)
    x2, y2 = _project_aeqd(lon2, lat2, lon0, lat0)
    xp, yp = _project_aeqd(np.asarray(lon_pts), np.asarray(lat_pts), lon0, lat0)

    # section unit vector and length
    v = np.array([x2 - x1, y2 - y1], dtype=float)
    L = np.hypot(v[0], v[1])
    if L <= 0:
        raise ValueError("Section endpoints are identical.")
    e = v / L

    # vector from P1 to each point
    wx = xp - x1
    wy = yp - y1

    # along-track (signed) and cross-track (signed)
    s = e[0] * wx + e[1] * wy  # meters
    n = e[0] * wy - e[1] * wx  # right-hand signed distance (meters)

    return s, n, L  # along, cross, section length


# -----------------------------------------------------------


def cross_section_from_catalog(
    catalog_df: pd.DataFrame,
    line_start: tuple,  # (lat1, lon1)
    line_end: tuple,  # (lat2, lon2)
    *,
    lat_col="latitude",
    lon_col="longitude",
    depth_col="depth_km",
    mag_col="Mw",
    width_km=50.0,
) -> pd.DataFrame:
    """
    Select earthquakes within a corridor around the profile and compute along-section distances.

    Returns a new DataFrame with columns:
        's_km' (along-section, origin at line_start),
        'n_km' (signed cross-track),
        'depth_km' (positive-down),
        'Mw' (if present), plus original metadata.
    """
    lat1, lon1 = float(line_start[0]), float(line_start[1])
    lat2, lon2 = float(line_end[0]), float(line_end[1])

    if not {lat_col, lon_col, depth_col}.issubset(catalog_df.columns):
        raise KeyError(f"catalog_df must have '{lat_col}', '{lon_col}', '{depth_col}'.")

    lats = catalog_df[lat_col].values
    lons = catalog_df[lon_col].values

    s_m, n_m, L_m = _profile_axes(lon1, lat1, lon2, lat2, lons, lats)
    width_m = float(width_km) * 1e3

    mask = np.isfinite(s_m) & np.isfinite(n_m) & (np.abs(n_m) <= width_m)
    sub = catalog_df.loc[mask].copy()
    sub["s_km"] = s_m[mask] / 1e3
    sub["n_km"] = n_m[mask] / 1e3
    # ensure depth positive-down
    sub["depth_km"] = sub[depth_col].astype(float)
    if mag_col in sub.columns:
        sub["Mw"] = sub[mag_col].astype(float)
    sub.attrs["section_length_km"] = L_m / 1e3
    return sub


def plot_hikurangi_cross_section(
    section_df: pd.DataFrame,
    *,
    slab_profile=None,  # None, or DataFrame with ['distance_km','depth_km'], or (dist_km, depth_km)
    title="Hikurangi cross-section",
    size_by_mag=True,
    mag_col="Mw",
    cmap="viridis",
    station_distances_km=None,  # optional ticks/markers along section
    xlim=None,
    ylim=None,
    ax=None,
):
    """
    Plot a cross-section: earthquakes within the corridor and optional slab interface.

    Parameters
    ----------
    section_df : DataFrame from cross_section_from_catalog (needs 's_km','depth_km').
    slab_profile : None | DataFrame | tuple
        If DataFrame, expects columns 'distance_km' and 'depth_km'.
        If tuple, pass (distance_km_array, depth_km_array).
    size_by_mag : bool
        Scale marker size by Mw if available; otherwise fixed.
    """
    if not {"s_km", "depth_km"}.issubset(section_df.columns):
        raise KeyError("section_df must contain 's_km' and 'depth_km'.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4.6))
    else:
        fig = ax.figure

    # Earthquakes
    x = section_df["s_km"].to_numpy()
    y = section_df["depth_km"].to_numpy()

    if size_by_mag and (mag_col in section_df.columns):
        Mw = section_df[mag_col].to_numpy()
        # perceptual size mapping (tunable)
        sizes = 10.0 + 8.0 * (Mw - np.nanmin(Mw)) ** 2
    else:
        Mw = None
        sizes = 18.0

    sc = ax.scatter(
        x,
        y,
        s=sizes,
        c=(Mw if Mw is not None else y),
        cmap=cmap,
        alpha=0.85,
        edgecolor="k",
        linewidths=0.2,
    )

    # Slab interface (optional)
    if slab_profile is not None:
        if isinstance(slab_profile, pd.DataFrame):
            if not {"distance_km", "depth_km"}.issubset(slab_profile.columns):
                raise KeyError(
                    "slab_profile DataFrame must have 'distance_km' and 'depth_km'."
                )
            xs, ys = (
                slab_profile["distance_km"].to_numpy(),
                slab_profile["depth_km"].to_numpy(),
            )
        else:
            xs, ys = slab_profile
            xs = np.asarray(xs, float)
            ys = np.asarray(ys, float)
        order = np.argsort(xs)
        ax.plot(xs[order], ys[order], color="k", lw=2, label="Slab interface")

    # Stations/markers along section (optional)
    if station_distances_km is not None:
        for d in np.atleast_1d(station_distances_km):
            ax.axvline(float(d), color="grey", lw=0.8, ls="--", alpha=0.6)

    # Axes aesthetics
    ax.set_xlabel("Distance along section (km)")
    ax.set_ylabel("Depth (km)")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.invert_yaxis()  # depth positive-down
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    cbar = fig.colorbar(sc, ax=ax, pad=0.01)
    cbar.set_label("Magnitude Mw" if Mw is not None else "Depth (km)")

    ax.legend(loc="best")
    fig.tight_layout()
    return fig, ax


# ------------------ convenience wrapper ------------------
def plot_hikurangi_section_from_catalog(
    catalog_df: pd.DataFrame,
    line_start: tuple,
    line_end: tuple,
    *,
    lat_col="latitude",
    lon_col="longitude",
    depth_col="depth_km",
    mag_col="Mw",
    width_km=50.0,
    slab_profile=None,
    title=None,
):
    """
    One-call convenience: select events within corridor and plot the cross-section.
    """
    sub = cross_section_from_catalog(
        catalog_df,
        line_start,
        line_end,
        lat_col=lat_col,
        lon_col=lon_col,
        depth_col=depth_col,
        mag_col=mag_col,
        width_km=width_km,
    )
    if title is None:
        title = "Hikurangi cross-section"
    return plot_hikurangi_cross_section(
        sub, slab_profile=slab_profile, title=title, size_by_mag=True, mag_col=mag_col
    )
