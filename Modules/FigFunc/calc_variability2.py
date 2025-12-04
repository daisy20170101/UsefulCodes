import seissolxdmf
import numpy as np
import matplotlib.tri as tri


# Rjb calculation and GME loading code


def closest_point_to_line(point, p1, p2, p3, p4):
    """
    Compute the closest point and distance from a 3D point to a polyline
    defined by four points p1, p2, p3, p4 (connected in sequence).
    The polyline consists of three segments: p1->p2, p2->p3, p3->p4.

    Parameters
    ----------
    point : (3,) array_like
        Query point in 3D.
    p1, p2, p3, p4 : (3,) array_like
        Four points defining the polyline in order.

    Returns
    -------
    dist : float
        Euclidean distance from point to the closest point on the polyline.
    closest : (3,) ndarray
        Coordinates of the closest point on the polyline.
    segment_id : int
        Index of the segment containing the closest point (0, 1, or 2).
    t : float
        Parameter along the segment (0 to 1), where closest = p_start + t * (p_end - p_start).

    Examples
    --------
    >>> point = np.array([1.0, 1.0, 1.0])
    >>> p1 = np.array([0.0, 0.0, 0.0])
    >>> p2 = np.array([1.0, 0.0, 0.0])
    >>> p3 = np.array([1.0, 1.0, 0.0])
    >>> p4 = np.array([0.0, 1.0, 0.0])
    >>> dist, closest, seg_id, t = closest_point_to_line(point, p1, p2, p3, p4)
    """
    P = np.asarray(point, dtype=float)
    points = [np.asarray(p, dtype=float) for p in [p1, p2, p3, p4]]

    # Define the three segments
    segments = [(points[i], points[i+1]) for i in range(3)]

    min_dist = np.inf
    closest_point = None
    closest_segment = None
    closest_t = None

    for seg_id, (p_start, p_end) in enumerate(segments):
        # Vector from segment start to end
        seg_vec = p_end - p_start
        seg_length_sq = np.dot(seg_vec, seg_vec)

        if seg_length_sq < 1e-12:
            # Degenerate segment (point)
            dist = np.linalg.norm(P - p_start)
            closest = p_start
            t = 0.0
        else:
            # Project point onto the line defined by the segment
            # t = (P - p_start) · seg_vec / ||seg_vec||^2
            t = np.dot(P - p_start, seg_vec) / seg_length_sq

            # Clamp t to [0, 1] to stay on the segment
            t = np.clip(t, 0.0, 1.0)

            # Closest point on the segment
            closest = p_start + t * seg_vec
            dist = np.linalg.norm(P - closest)

        # Update if this is the closest so far
        if dist < min_dist:
            min_dist = dist
            closest_point = closest
            closest_segment = seg_id
            closest_t = t

    return min_dist


def closest_point_to_rectangle(point, p1, p2, p3, p4, eps=1e-12):
    """
    Compute the closest point and distance from a 3D point to a rectangle
    defined by four vertices p1, p2, p3, p4 (ordered around the perimeter).
    Assumes a planar rectangle (p1->p2 and p1->p4 are adjacent edges).

    Parameters
    ----------
    point : (3,) array_like
        Query point in 3D.
    p1, p2, p3, p4 : (3,) array_like
        Vertices of the rectangle in consistent order (e.g., counterclockwise).
        Edges are assumed to be (p1->p2) and (p1->p4).
    eps : float
        Tolerance for degeneracy checks.

    Returns
    -------
    dist : float
        Euclidean distance from point to the rectangle.
    closest : (3,) ndarray
        Coordinates of the closest point on (or within) the rectangle.
    uv : (2,) ndarray
        Local coordinates (u, v) in edge-length units, clamped to [0, Lu] x [0, Lv].
        u measures along (p2 - p1), v along (p4 - p1).
    inside : bool
        True if the orthogonal projection of the point falls inside the rectangle.
    """
    P = np.asarray(point, dtype=float)
    P1 = np.asarray(p1, dtype=float)
    P2 = np.asarray(p2, dtype=float)
    P4 = np.asarray(p4, dtype=float)

    # Edge vectors from p1
    u = P2 - P1
    v = P4 - P1
    Lu = np.linalg.norm(u)
    Lv = np.linalg.norm(v)

    if Lu < eps or Lv < eps:
        # Degenerate: treat as a segment or point
        if Lu < eps and Lv < eps:
            closest = P1.copy()
            return np.linalg.norm(P - closest), closest, np.array([0.0, 0.0]), False
        if Lu < eps:  # collapse to segment P1--P4
            t = np.clip(np.dot(P - P1, v) / (Lv*Lv), 0.0, 1.0)
            closest = P1 + t * v
            return np.linalg.norm(P - closest), closest, np.array([0.0, t*Lv]), (0.0 <= t <= 1.0)
        else:         # collapse to segment P1--P2
            s = np.clip(np.dot(P - P1, u) / (Lu*Lu), 0.0, 1.0)
            closest = P1 + s * u
            return np.linalg.norm(P - closest), closest, np.array([s*Lu, 0.0]), (0.0 <= s <= 1.0)

    # Build an orthonormal frame from the edges (robust to slight non-orthogonality)
    u_hat = u / Lu
    v_ortho = v - np.dot(v, u_hat) * u_hat
    Lv_ortho = np.linalg.norm(v_ortho)
    if Lv_ortho < eps:
        # Edges are (nearly) colinear; treat as segment along u
        s = np.clip(np.dot(P - P1, u) / (Lu*Lu), 0.0, 1.0)
        closest = P1 + s * u
        return np.linalg.norm(P - closest), closest, np.array([s*Lu, 0.0]), (0.0 <= s <= 1.0)
    v_hat = v_ortho / Lv_ortho

    # Express point in this local (u_hat, v_hat) basis
    w = P - P1
    su = np.dot(w, u_hat)           # signed distance along u_hat
    sv = np.dot(w, v_hat)           # signed distance along v_hat

    # Rectangle extents in this orthonormalized frame
    # (use original edge lengths along their respective directions)
    u_min, u_max = 0.0, Lu
    v_min, v_max = 0.0, Lv_ortho

    # Clamp to the rectangle
    su_c = np.clip(su, u_min, u_max)
    sv_c = np.clip(sv, v_min, v_max)

    # Closest point on the rectangle in world coordinates
    closest = P1 + su_c * u_hat + sv_c * v_hat

    # Distance and flags
    dist = np.linalg.norm(P - closest)
    inside = (u_min <= su <= u_max) and (v_min <= sv <= v_max)

    # uv returned in edge-length units (0..Lu, 0..Lv)
    uv = np.array([su_c, sv_c], dtype=float)
    return dist


def closest_point_to_rectangles(point, rectangles):
    """
    Compute the closest distance from a 3D point to multiple rectangles.

    Parameters
    ----------
    point : (3,) array_like
        Query point in 3D.
    rectangles : list of tuples
        List of rectangles, where each rectangle is a tuple of 4 vertices (p1, p2, p3, p4).
        Each vertex should be a (3,) array_like.

    Returns
    -------
    min_dist : float
        Minimum Euclidean distance from point to any of the rectangles.
    closest_rect_idx : int
        Index of the rectangle with the closest point.

    Examples
    --------
    >>> point = np.array([0.0, 0.0, 5.0])
    >>> rect1 = (np.array([0, 0, 0]), np.array([1, 0, 0]),
    ...          np.array([1, 1, 0]), np.array([0, 1, 0]))
    >>> rect2 = (np.array([2, 0, 0]), np.array([3, 0, 0]),
    ...          np.array([3, 1, 0]), np.array([2, 1, 0]))
    >>> rect3 = (np.array([4, 0, 0]), np.array([5, 0, 0]),
    ...          np.array([5, 1, 0]), np.array([4, 1, 0]))
    >>> min_dist, idx = closest_point_to_rectangles(point, [rect1, rect2, rect3])
    """
    min_dist = np.inf
    closest_rect_idx = None

    for idx, rect in enumerate(rectangles):
        p1, p2, p3, p4 = rect
        dist = closest_point_to_rectangle(point, p1, p2, p3, p4)

        if dist < min_dist:
            min_dist = dist
            closest_rect_idx = idx

    return min_dist, closest_rect_idx


def calc_rjb(x, y, p1=(1719988.0, 5407437.69), p2=(1783672.647, 5452039.0)):

    s = (x, y)
    t = ((s[0] - p1[0]) * (p2[0] - p1[0]) + (s[1] - p1[1]) * (p2[1] - p1[1])) / (
        (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2
    )

    if t < 0:
        d = np.sqrt((s[0] - p1[0]) ** 2 + (s[1] - p1[1]) ** 2)
    else:
        if t > 1:
            d = np.sqrt((s[0] - p2[0]) ** 2 + (s[1] - p2[1]) ** 2)
        else:
            c = (p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]))
            d = np.sqrt((s[0] - c[0]) ** 2 + (s[1] - c[1]) ** 2)

    return d

# hik slab surface
p1 = (1800650.75, 5504513.0)
p2 = (1888556.30, 5415316.55)
p3 = (1767281.22, 5312777.18)
p4 = (1693832.38, 5386412.84)


def calc_rjb3D(x, y, p1=(1, 1), p2=(1, 2), p3=(2, 2), p4=(2, 1)):
    """compute the shortest distance to the surfface projection of a fault plane"""

    s = (x, y)
    rjb = np.zeros(4)

    pxy = np.array([p1, p2, p3, p4])  # pxy is (4,2) in shape

    # check if s is inside the rectangle:

    i = 0

    t = (
        (s[0] - pxy[i, 0]) * (pxy[i + 1, 0] - pxy[i, 0])
        + (s[1] - pxy[i, 1]) * (pxy[i + 1, 1] - pxy[i, 1])
    ) / ((pxy[i, 0] - pxy[i + 1, 0]) ** 2 + (pxy[i + 1, 1] - pxy[i, 1]) ** 2)

    idd = 3

    q = (
        (s[0] - pxy[idd, 0]) * (pxy[0, 0] - pxy[idd, 0])
        + (s[1] - pxy[idd, 1]) * (pxy[0, 1] - pxy[idd, 1])
    ) / ((pxy[idd, 0] - pxy[0, 0]) ** 2 + (pxy[idd, 1] - pxy[0, 1]) ** 2)

    if t > 1 or t < 0 or q < 0 or q > 1:

        # if not:

        for i in range(3):
            # print('line:',i)

            t = (
                (s[0] - pxy[i, 0]) * (pxy[i + 1, 0] - pxy[i, 0])
                + (s[1] - pxy[i, 1]) * (pxy[i + 1, 1] - pxy[i, 1])
            ) / ((pxy[i, 0] - pxy[i + 1, 0]) ** 2 + (pxy[i + 1, 1] - pxy[i, 1]) ** 2)

            if t < 0:
                d = np.sqrt((s[0] - pxy[i, 0]) ** 2 + (s[1] - pxy[i, 1]) ** 2)
            else:
                if t > 1:
                    d = np.sqrt(
                        (s[0] - pxy[i + 1, 0]) ** 2 + (s[1] - pxy[i + 1, 1]) ** 2
                    )
                else:
                    c = (
                        pxy[i, 0] + t * (pxy[i + 1, 0] - pxy[i, 0]),
                        pxy[i, 1] + t * (pxy[i + 1, 1] - pxy[i, 1]),
                    )
                    d = np.sqrt((s[0] - c[0]) ** 2 + (s[1] - c[1]) ** 2)

            rjb[i] = d

        # print('line:',i+1)
        idd = i + 1

        t = (
            (s[0] - pxy[idd, 0]) * (pxy[0, 0] - pxy[idd, 0])
            + (s[1] - pxy[idd, 1]) * (pxy[0, 1] - pxy[idd, 1])
        ) / ((pxy[idd, 0] - pxy[0, 0]) ** 2 + (pxy[idd, 1] - pxy[0, 1]) ** 2)

        if t < 0:
            d = np.sqrt((s[0] - pxy[idd, 0]) ** 2 + (s[1] - pxy[idd, 1]) ** 2)
        else:
            if t > 1:
                d = np.sqrt((s[0] - pxy[0, 0]) ** 2 + (s[1] - pxy[0, 1]) ** 2)
            else:
                c = (
                    pxy[idd, 0] + t * (pxy[0, 0] - pxy[idd, 0]),
                    pxy[idd, 1] + t * (pxy[0, 1] - pxy[idd, 1]),
                )
                d = np.sqrt((s[0] - c[0]) ** 2 + (s[1] - c[1]) ** 2)
        rjb[i + 1] = d

        rjb_min = np.min(rjb)

    else:
        rjb_min = 0

    return rjb_min


def load_gm_data_slab(xdmfFilename, origin=(0, 0)):
    """
    the Rjb is calculated assuming the fault plane Ax+By+D=0

    """
    p1 = (1800650.75, 5504513.0)
    p2 = (1888556.30, 5415316.55)
    p3 = (1767281.22, 5312777.18)
    p4 = (1693832.38, 5386412.84)

    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    ndt = sx.ReadNdt()
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()

    print(surfxyz.shape, connect.shape)

    centers = (
        surfxyz[connect[:, 0]] + surfxyz[connect[:, 1]] + surfxyz[connect[:, 2]]
    ) / 3.0
    depi = np.sqrt((centers[:, 0] - origin[0]) ** 2 + (centers[:, 1] - origin[1]) ** 2)

    Rjb = []

    for i in range(len(centers[:, 1])):
        Rjb.append(calc_rjb3D(centers[i, 0], centers[i, 1], p1, p2, p3, p4))

    ##%%
    pga = sx.ReadData("PGA")
    pgv = sx.ReadData("PGV")
    sa1 = sx.ReadData("SA01.000s")
    sa3 = sx.ReadData("SA03.000s")

    return pga, pgv, sa1, sa3, depi, Rjb



def load_gm_data_crust(xdmfFilename, origin=(0, 0)):
    """
    the Rjb is calculated assuming the fault line: Ax+B=0

    """

    sx = seissolxdmf.seissolxdmf(xdmfFilename)

    ndt = sx.ReadNdt()
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()

    print(surfxyz.shape, connect.shape)

    centers = (
        surfxyz[connect[:, 0]] + surfxyz[connect[:, 1]] + surfxyz[connect[:, 2]]
    ) / 3.0
    depi = np.sqrt((centers[:, 0] - origin[0]) ** 2 + (centers[:, 1] - origin[1]) ** 2)

    Rjb = []

    for i in range(len(centers[:, 1])):
        Rjb.append(calc_rjb(centers[i, 0], centers[i, 1]))

    ##%%
    pga = sx.ReadData("PGA")
    pgv = sx.ReadData("PGV")
    sa1 = sx.ReadData("SA01.000s")
    sa3 = sx.ReadData("SA03.000s")

    return pga, pgv, sa1, sa3, depi, Rjb


def calc_gm_variability(rjb, sa1, nbins=40):
    """give the model name, number of bins in distance (km) and ground intensity index (either, pgv, sa1, or sa3)"""

    # group into bins of Rjb and compute mean and std

    rjb_bins = np.logspace(0, 2, nbins)
    mean_bins = np.zeros(nbins - 1)
    std_bins = np.zeros(nbins - 1)
   
    for ik, rjb_bin in enumerate(rjb_bins[0:-1]):
              
        mean_bins[ik] = np.mean(
            np.log10(                   
                sa1[np.where((rjb / 1e3 > rjb_bin) & (rjb / 1e3 < rjb_bins[ik + 1]))]
            )
        )
        std_bins[ik] = np.std(
            np.log10(
                sa1[np.where((rjb / 1e3 > rjb_bin) & (rjb / 1e3 < rjb_bins[ik + 1]))]
            )
        )

    return rjb_bins, mean_bins,                     std_bins
        
import pandas as pd      

def prepare_gm_table_from_model(
        site_table: pd.DataFrame | np.ndarray,
        xdmfFilename: str = 'c9s30r27A-GME-surface.xmdf',
        modelname: str='c9s30r27A',
        mag: float=7.2,
        period_list: list[str]=("pga","sa1","sa3") ,
        hypodepth: float = 20.0,
        interface_flag: bool=True,
        wairarapa_flag: bool=True,
        ) -> pd.DataFrame:

    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    surfxyz = sx.ReadGeometry()
    connect = sx.ReadConnect()

    print(surfxyz.shape, connect.shape)

    centers = (
        surfxyz[connect[:, 0]] + surfxyz[connect[:, 1]] + surfxyz[connect[:, 2]]
    ) / 3.0

    Rjb = []
    rrup = []

    # hik slab
    p1 = (1800650.75, 5504513.0,-30000.0)
    p2 = (1888556.30, 5415316.55,-10000.0)
    p3 = (1767281.22, 5312777.18,-10000.0)
    p4 = (1693832.38, 5386412.84,-30000.0)


    if interface_flag:
        for i in range(len(centers[:, 1])):
            Rjb.append(calc_rjb3D(centers[i, 0], centers[i, 1], p1, p2, p3, p4))
            rrup.append(closest_point_to_rectangle(np.array([centers[i, 0], centers[i, 1],centers[i, 2]]),np.array(p1), np.array(p2), np.array(p3),np.array(p4)))

    else:
        if wairarapa_flag:
            p1 = (1.71306e+06, 5.38688e+06 )
            p4 = (1.80433e+06, 5.4603e+06)
            p2 = (1.7438E+06 , 5.3975E+06)
            p3 = (1.7751E+06,	5.4295E+06)

            for i in range(len(centers[:, 1])):
                Rjb.append(closest_point_to_line((centers[i, 0], centers[i, 1]),p1,p2,p3,p4))
        else:
            for i in range(len(centers[:, 1])):
                Rjb.append(calc_rjb(centers[i, 0], centers[i, 1]))
    
        rrup = Rjb.copy()

    ##%%
    pga = sx.ReadData("PGA")
    sa1 = sx.ReadData("SA01.000s")
    
    try:
        sa3 = sx.ReadData("SA03.000s")
    except:
        try:
            sa3 = sx.ReadData("SA02.000s")
            print('replaced by SA 2.0')
        except:
            print('skip:' +modelname)
            pass
    # sa0_3=sx.ReadData("SA00.300s")

    # Find closest points in centers for each (x,y) in site_table
    if isinstance(site_table, pd.DataFrame):
        site_x = site_table['x'].values
        site_y = site_table['y'].values
        site_sta = site_table['sta'].values if 'sta' in site_table.columns else None
        site_vs30 = site_table['Vs30'].values if 'Vs30' in site_table.columns else None
    else:
        # Assume numpy array with columns [x, y, sta, Vs30]
        site_x = site_table[:, 0]
        site_y = site_table[:, 1]
        site_sta = site_table[:, 2] if site_table.shape[1] > 2 else None
        site_vs30 = site_table[:, 3] if site_table.shape[1] > 3 else None

    # Find nearest center for each site in site_table
    selected_indices = []
    sta_list = []
    vs30_list = []

    for site_idx in range(len(site_x)):
        # Calculate distances from this site to all centers
        distances = np.sqrt((centers[:, 0] - site_x[site_idx])**2 + (centers[:, 1] - site_y[site_idx])**2)
        closest_center_idx = np.argmin(distances)

        selected_indices.append(closest_center_idx)

        # Get sta and Vs30 from the site_table
        if site_sta is not None:
            sta_list.append(site_sta[site_idx])
        else:
            sta_list.append(f"site_{site_idx}")

        if site_vs30 is not None:
            vs30_list.append(site_vs30[site_idx])
        else:
            vs30_list.append(np.nan)

    # Filter data to only include selected centers
    selected_indices = np.array(selected_indices)
    Rjb_array = np.array(Rjb)
    rrup_array=np.array(rrup)

    table=pd.DataFrame({
        "pSA_1.0": sa1[selected_indices]/9.8,
        "pSA_3.0": sa3[selected_indices]/9.8,
        # "pSA_0.3": sa0_3[selected_indices]/9.8,
        "PGA": pga[selected_indices]/9.8,
        "r_jb": Rjb_array[selected_indices]/1e3,
        "r_rup": rrup_array[selected_indices]/1e3,
        "evid": np.full(len(selected_indices),modelname),
        "mag": np.full(len(selected_indices),mag),
        "sta": sta_list,
        "Vs30": vs30_list,
    })

    print(table.keys())

    return table

def prepare_resid_table_at22crust(
        df: pd.DataFrame,
        ) -> pd.DataFrame :

    from openquake.hazardlib import gsim

    from openquake.hazardlib.gsim.nz22 import atkinson_2022
    gsim = atkinson_2022.Atkinson2022Crust()

    print("Sites:", gsim.REQUIRES_SITES_PARAMETERS)
    print("Rupture:", gsim.REQUIRES_RUPTURE_PARAMETERS)
    print("Dists:", gsim.REQUIRES_DISTANCES)

    from openquake.hazardlib.imt import PGA,SA
    from openquake.hazardlib.const import StdDev
    from openquake.hazardlib.contexts import SitesContext, DistancesContext, RuptureContext

    all_residuals = []

    for i in range(len(df)):
        mag, rrup, vs30, sta = df.iloc[i][["mag", "r_rup", "Vs30","sta"]]

        resid_list = []
        
        # PGA
        imt = PGA()
        vs30  = np.full(1, vs30)            # m/s
        mag = np.full(1,mag)
        rrup = np.full(1,rrup)
        # hypo_depth = np.full(1,hypo_depth)

        sites = SitesContext()
        sites.sids = np.arange(1, dtype=int)            # REQUIRED int ids
        sites.vs30 = vs30.astype(float)

        rup = RuptureContext()
        rup.mag = mag.astype(float)

        dists = DistancesContext()
        dists.rrup = rrup.astype(float)
        
        mean_ln, [sigma_tot] = gsim.get_mean_and_stddevs(
                sites, rup, dists, imt, [StdDev.TOTAL]
            )
        
        pga_median = mean_ln[0]
        resid_list.append(np.log(df.iloc[i]["PGA"]) - pga_median)
        
        # Pre-filter once
        psa_columns_filtered = [col for col in df.columns 
                            if col.startswith('pSA_') and float(col.split("_")[1]) <= 10.0]
        periods = [float(col.split("_")[1]) for col in psa_columns_filtered]
        residual_cols = ['resid_PGA'] + [f'resid_{col}' for col in psa_columns_filtered]

        # PsA
        for col, period in zip(psa_columns_filtered, periods):

            imt = SA(period)
            mean_ln, [sigma_tot] = gsim.get_mean_and_stddevs(
                sites, rup, dists, imt, [StdDev.TOTAL]
            )
            pred_median = mean_ln[0]
            resid_list.append(np.log(df.iloc[i][col]) - pred_median)
        
        all_residuals.append(resid_list)
        
        if (i + 1) % 500 == 0:  # Progress indicator
            print(f"Processed {i + 1}/{len(df)} rows")


    # Create residuals DataFrame
    df_residuals = pd.DataFrame(all_residuals, columns=residual_cols, index=df.index)

    # Add ALL residual columns at once - NO fragmentation!
    df_resd_all = pd.concat([df, df_residuals], axis=1)

    print("Done! Residuals added to df_sel")

    return df_resd_all

def prepare_resid_table_at22interface(
        df: pd.DataFrame,
        ) -> pd.DataFrame :

    from openquake.hazardlib import gsim

    from openquake.hazardlib.gsim.nz22 import atkinson_2022
    gsim = atkinson_2022.Atkinson2022SInter()

    print("Sites:", gsim.REQUIRES_SITES_PARAMETERS)
    print("Rupture:", gsim.REQUIRES_RUPTURE_PARAMETERS)
    print("Dists:", gsim.REQUIRES_DISTANCES)

    from openquake.hazardlib.imt import PGA,SA
    from openquake.hazardlib.const import StdDev
    from openquake.hazardlib.contexts import SitesContext, DistancesContext, RuptureContext

    all_residuals = []

    backarc_flag=False

    for i in range(len(df)):
        mag, rrup, vs30, sta = df.iloc[i][["mag", "r_rup", "Vs30","sta"]]

        resid_list = []
        
        # PGA
        imt = PGA()
        vs30  = np.full(1, vs30)            # m/s
        mag = np.full(1,mag)
        rrup = np.full(1,rrup)
        # hypo_depth = np.full(1,hypo_depth)

        sites = SitesContext()
        sites.sids = np.arange(1, dtype=int)            # REQUIRED int ids
        sites.vs30 = vs30.astype(float)
        sites.backarc = np.full(1, backarc_flag, dtype=bool)  # REQUIRED for BC Hydro


        rup = RuptureContext()
        rup.mag = mag.astype(float)

        dists = DistancesContext()
        dists.rrup = rrup.astype(float)
        
        mean_ln, [sigma_tot] = gsim.get_mean_and_stddevs(
                sites, rup, dists, imt, [StdDev.TOTAL]
            )
        
        pga_median = mean_ln[0]
        resid_list.append(np.log(df.iloc[i]["PGA"]) - pga_median)
        
        # Pre-filter once
        psa_columns_filtered = [col for col in df.columns 
                            if col.startswith('pSA_') and float(col.split("_")[1]) <= 10.0]
        periods = [float(col.split("_")[1]) for col in psa_columns_filtered]
        residual_cols = ['resid_PGA'] + [f'resid_{col}' for col in psa_columns_filtered]

        # PsA
        for col, period in zip(psa_columns_filtered, periods):

            imt = SA(period)
            mean_ln, [sigma_tot] = gsim.get_mean_and_stddevs(
                sites, rup, dists, imt, [StdDev.TOTAL]
            )
            pred_median = mean_ln[0]
            resid_list.append(np.log(df.iloc[i][col]) - pred_median)
        
        all_residuals.append(resid_list)
        
        if (i + 1) % 500 == 0:  # Progress indicator
            print(f"Processed {i + 1}/{len(df)} rows")


    # Create residuals DataFrame
    df_residuals = pd.DataFrame(all_residuals, columns=residual_cols, index=df.index)

    # Add ALL residual columns at once - NO fragmentation!
    df_resd_all = pd.concat([df, df_residuals], axis=1)

    print("Done! Residuals added to df_sel")

    return df_resd_all