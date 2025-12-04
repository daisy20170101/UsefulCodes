import numpy as np


def calc_dist_rect(x, y, p=[(1, 1), (0, 1), (1, 0), (0, 0)]):
    """compute the shortest distance to the surfface projection of a fault plane"""

    s = (x, y)

    rjb = np.zeros(4)

    pxy = np.array(p)  # pxy is (4,2) in shape

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


def point_to_polygon_distance(point, polygon):
    """
    Calculate the shortest distance from a point to an arbitrary polygon.

    Parameters:
    -----------
    point : tuple or array-like
        (x, y) coordinates of the point
    polygon : array-like
        List or array of (x, y) coordinates defining the polygon vertices
        Shape: (n_vertices, 2)

    Returns:
    --------
    float
        Shortest distance from point to polygon (0 if point is inside)
    """
    point = np.array(point)
    polygon = np.array(polygon)

    # Check if point is inside polygon using ray casting algorithm
    def point_in_polygon(px, py, poly):
        n = len(poly)
        inside = False
        p1x, p1y = poly[0]
        for i in range(1, n + 1):
            p2x, p2y = poly[i % n]
            if py > min(p1y, p2y):
                if py <= max(p1y, p2y):
                    if px <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (py - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or px <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    # If point is inside polygon, distance is 0
    if point_in_polygon(point[0], point[1], polygon):
        return 0.0

    # Calculate minimum distance to all edges
    min_distance = float("inf")
    n_vertices = len(polygon)

    for i in range(n_vertices):
        # Get current edge (from vertex i to vertex (i+1) % n_vertices)
        v1 = polygon[i]
        v2 = polygon[(i + 1) % n_vertices]

        # Calculate distance from point to line segment
        edge_vector = v2 - v1
        point_vector = point - v1

        # Project point onto line segment
        edge_length_sq = np.dot(edge_vector, edge_vector)

        if edge_length_sq == 0:
            # Degenerate edge (point)
            distance = np.linalg.norm(point - v1)
        else:
            # Parameter t for projection onto line
            t = np.dot(point_vector, edge_vector) / edge_length_sq

            if t < 0:
                # Closest point is v1
                closest_point = v1
            elif t > 1:
                # Closest point is v2
                closest_point = v2
            else:
                # Closest point is on the line segment
                closest_point = v1 + t * edge_vector

            distance = np.linalg.norm(point - closest_point)

        min_distance = min(min_distance, distance)

    return min_distance


def point_to_polygon_distance_vectorized(points, polygon):
    """
    Calculate shortest distances from multiple points to a polygon (vectorized version).

    Parameters:
    -----------
    points : array-like
        Array of (x, y) coordinates, shape: (n_points, 2)
    polygon : array-like
        List or array of (x, y) coordinates defining the polygon vertices
        Shape: (n_vertices, 2)

    Returns:
    --------
    numpy.ndarray
        Array of shortest distances for each point
    """
    points = np.array(points)
    if points.ndim == 1:
        points = points.reshape(1, -1)

    distances = np.array(
        [point_to_polygon_distance(point, polygon) for point in points]
    )
    return distances
