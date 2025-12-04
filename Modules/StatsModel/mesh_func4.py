#!/usr/bin/env python3
"""
Mesh Topology and Discrete Curvature Computation
Complete implementation of neighbor finding and curvature operators for triangle meshes
"""

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#!/usr/bin/env python3
"""
Updated Robust Curvature Computation Functions
Handles edge cases and validates vertex indexing properly
"""

import numpy as np
from collections import defaultdict
import warnings


def find_neighbors(vertex_idx, faces):
    """
    Find all neighboring vertices of a given vertex in a triangle mesh.

    Parameters:
    -----------
    vertex_idx : int
        Index of the vertex whose neighbors we want to find
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices

    Returns:
    --------
    neighbors : list
        List of neighboring vertex indices

    Raises:
    -------
    ValueError
        If vertex_idx is not found in any face
    """
    if not isinstance(vertex_idx, (int, np.integer)):
        raise TypeError(f"vertex_idx must be integer, got {type(vertex_idx)}")

    neighbors = set()
    vertex_found = False

    # Go through all faces and find ones containing our vertex
    for face in faces:
        if vertex_idx in face:
            vertex_found = True
            # Add the other two vertices of this face as neighbors
            for v in face:
                if v != vertex_idx:
                    neighbors.add(v)

    if not vertex_found:
        raise ValueError(f"Vertex {vertex_idx} not found in any face")

    return list(neighbors)


def get_ordered_neighbors(vertex_idx, faces):
    """
    Get neighbors in order around the vertex (for angle calculations).

    Parameters:
    -----------
    vertex_idx : int
        Index of the vertex
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices

    Returns:
    --------
    ordered_neighbors : list
        Ordered list of neighbor vertices around the central vertex
    """
    # Find all faces containing this vertex
    incident_faces = []
    for i, face in enumerate(faces):
        if vertex_idx in face:
            incident_faces.append((i, face))

    if len(incident_faces) < 2:
        # Not enough faces to determine order
        return find_neighbors(vertex_idx, faces)

    # Build adjacency between faces
    face_edges = {}
    for face_idx, face in incident_faces:
        # Find the edge opposite to our vertex
        other_vertices = [v for v in face if v != vertex_idx]
        if len(other_vertices) == 2:
            edge = tuple(sorted(other_vertices))
            face_edges[edge] = face_idx

    # Try to build ordered sequence
    if not face_edges:
        return find_neighbors(vertex_idx, faces)

    # Start with first edge
    edges = list(face_edges.keys())
    ordered_neighbors = list(edges[0])
    used_edges = {edges[0]}

    # Build the ring by following connectivity
    while len(used_edges) < len(edges):
        last_vertex = ordered_neighbors[-1]
        found_next = False

        for edge in edges:
            if edge not in used_edges:
                if last_vertex in edge:
                    # Add the other vertex in this edge
                    next_vertex = edge[0] if edge[1] == last_vertex else edge[1]
                    ordered_neighbors.append(next_vertex)
                    used_edges.add(edge)
                    found_next = True
                    break

        if not found_next:
            break

    # Remove potential duplicate at the end
    if len(ordered_neighbors) > 2 and ordered_neighbors[0] == ordered_neighbors[-1]:
        ordered_neighbors = ordered_neighbors[:-1]

    return ordered_neighbors


def compute_cotangent_weights(vertex_idx, vertices, neighbors, faces):
    """
    Compute cotangent weights for discrete Laplacian operator.

    Parameters:
    -----------
    vertex_idx : int
        Index of central vertex
    vertices : array_like, shape (n_vertices, 3)
        3D coordinates of all vertices
    neighbors : list
        List of neighbor vertex indices
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices

    Returns:
    --------
    weights : dict
        Dictionary mapping neighbor index to cotangent weight
    """
    weights = {}
    center = vertices[vertex_idx]

    for neighbor_idx in neighbors:
        neighbor = vertices[neighbor_idx]
        total_weight = 0.0

        # Find all faces containing both vertex_idx and neighbor_idx
        shared_faces = []
        for face in faces:
            if vertex_idx in face and neighbor_idx in face:
                shared_faces.append(face)

        # For each shared face, compute cotangent of angle at third vertex
        for face in shared_faces:
            third_vertex_idx = None
            for v in face:
                if v != vertex_idx and v != neighbor_idx:
                    third_vertex_idx = v
                    break

            if third_vertex_idx is not None:
                third_vertex = vertices[third_vertex_idx]

                # Vectors from third vertex to the other two
                v1 = center - third_vertex
                v2 = neighbor - third_vertex

                # Compute cotangent of angle
                dot_product = np.dot(v1, v2)
                cross_product = np.cross(v1, v2)
                cross_magnitude = np.linalg.norm(cross_product)

                if cross_magnitude > 1e-10:  # Avoid division by zero
                    cot_angle = dot_product / cross_magnitude
                    total_weight += cot_angle

        weights[neighbor_idx] = max(total_weight, 1e-10)  # Avoid negative weights

    return weights


def discrete_mean_curvature(vertex_idx, vertices, faces, method="uniform"):
    """
    Compute discrete mean curvature at a vertex.

    Parameters:
    -----------
    vertex_idx : int
        Index of the vertex
    vertices : array_like, shape (n_vertices, 3)
        3D coordinates of all vertices
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices
    method : str, optional
        Method for computing weights ('uniform', 'cotangent')

    Returns:
    --------
    H : float
        Discrete mean curvature

    Raises:
    -------
    ValueError
        If vertex_idx is invalid or has insufficient neighbors
    IndexError
        If vertex_idx is out of range for vertices array
    """
    # Validate inputs
    if vertex_idx < 0 or vertex_idx >= len(vertices):
        raise IndexError(f"vertex_idx {vertex_idx} out of range [0, {len(vertices)-1}]")

    try:
        neighbors = find_neighbors(vertex_idx, faces)
    except ValueError as e:
        raise ValueError(f"Cannot compute curvature: {e}")

    if len(neighbors) < 2:
        warnings.warn(
            f"Vertex {vertex_idx} has only {len(neighbors)} neighbors, "
            "curvature may be unreliable"
        )
        return 0.0

    center = vertices[vertex_idx]

    # Compute Laplacian using different weighting schemes
    if method == "uniform":
        # Simple uniform weights
        laplacian = np.zeros(3)
        for neighbor_idx in neighbors:
            neighbor = vertices[neighbor_idx]
            laplacian += neighbor - center
        laplacian /= len(neighbors)

    elif method == "cotangent":
        # Cotangent weights (more accurate but more complex)
        try:
            weights = compute_cotangent_weights(vertex_idx, vertices, neighbors, faces)
            laplacian = np.zeros(3)
            total_weight = 0.0

            for neighbor_idx, weight in weights.items():
                neighbor = vertices[neighbor_idx]
                laplacian += weight * (neighbor - center)
                total_weight += weight

            if total_weight > 1e-10:
                laplacian /= total_weight
            else:
                # Fallback to uniform weights
                laplacian = np.zeros(3)
                for neighbor_idx in neighbors:
                    neighbor = vertices[neighbor_idx]
                    laplacian += neighbor - center
                laplacian /= len(neighbors)

        except Exception as e:
            warnings.warn(
                f"Cotangent weights failed for vertex {vertex_idx}: {e}. "
                "Using uniform weights."
            )
            # Fallback to uniform weights
            laplacian = np.zeros(3)
            for neighbor_idx in neighbors:
                neighbor = vertices[neighbor_idx]
                laplacian += neighbor - center
            laplacian /= len(neighbors)

    else:
        raise ValueError(f"Unknown method: {method}. Use 'uniform' or 'cotangent'")

    # Mean curvature is half the magnitude of the Laplacian
    H = 0.5 * np.linalg.norm(laplacian)

    return H


def discrete_gaussian_curvature(vertex_idx, vertices, faces):
    """
    Compute discrete Gaussian curvature using angle defect method.

    Parameters:
    -----------
    vertex_idx : int
        Index of the vertex
    vertices : array_like, shape (n_vertices, 3)
        3D coordinates of all vertices
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices

    Returns:
    --------
    K : float
        Discrete Gaussian curvature
    """
    # Validate inputs
    if vertex_idx < 0 or vertex_idx >= len(vertices):
        raise IndexError(f"vertex_idx {vertex_idx} out of range [0, {len(vertices)-1}]")

    try:
        neighbors = get_ordered_neighbors(vertex_idx, faces)
    except ValueError as e:
        raise ValueError(f"Cannot compute Gaussian curvature: {e}")

    if len(neighbors) < 3:
        warnings.warn(
            f"Vertex {vertex_idx} has only {len(neighbors)} neighbors, "
            "Gaussian curvature may be unreliable"
        )
        return 0.0

    center = vertices[vertex_idx]
    total_angle = 0.0

    # Compute sum of angles around the vertex
    for i in range(len(neighbors)):
        v1_idx = neighbors[i]
        v2_idx = neighbors[(i + 1) % len(neighbors)]

        v1 = vertices[v1_idx] - center
        v2 = vertices[v2_idx] - center

        # Compute angle between v1 and v2
        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)

        if v1_norm > 1e-10 and v2_norm > 1e-10:
            cos_angle = np.dot(v1, v2) / (v1_norm * v2_norm)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Numerical stability
            angle = np.arccos(cos_angle)
            total_angle += angle

    # Gaussian curvature is angle defect (2π - sum of angles)
    K = 2 * np.pi - total_angle

    return K


def validate_mesh(vertices, faces):
    """
    Validate mesh data and return information about vertex usage.

    Parameters:
    -----------
    vertices : array_like, shape (n_vertices, 3)
        3D coordinates of vertices
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices

    Returns:
    --------
    validation_info : dict
        Dictionary containing validation results and mesh statistics
    """
    vertices = np.asarray(vertices)
    faces = np.asarray(faces)

    n_vertices = len(vertices)

    # Find which vertices are used in faces
    used_vertex_indices = set()
    max_index = -1
    min_index = float("inf")

    for face in faces:
        for vertex_idx in face:
            used_vertex_indices.add(vertex_idx)
            max_index = max(max_index, vertex_idx)
            min_index = min(min_index, vertex_idx)

    # Check for invalid indices
    invalid_indices = {
        idx for idx in used_vertex_indices if idx < 0 or idx >= n_vertices
    }

    # Find unused vertices
    expected_indices = set(range(n_vertices))
    unused_vertices = expected_indices - used_vertex_indices

    # Check if indexing is consecutive
    is_consecutive = used_vertex_indices == expected_indices

    validation_info = {
        "n_vertices": n_vertices,
        "n_faces": len(faces),
        "used_vertex_indices": sorted(used_vertex_indices),
        "unused_vertices": sorted(unused_vertices),
        "invalid_indices": sorted(invalid_indices),
        "min_vertex_index": min_index if min_index != float("inf") else None,
        "max_vertex_index": max_index,
        "is_consecutive_indexing": is_consecutive,
        "is_valid": len(invalid_indices) == 0,
    }

    return validation_info


def compute_mesh_curvatures(vertices, faces, method="uniform", validate=True):
    """
    Compute curvatures for entire mesh with robust error handling.

    Parameters:
    -----------
    vertices : array_like, shape (n_vertices, 3)
        3D coordinates of vertices
    faces : array_like, shape (n_faces, 3)
        Triangle faces defined by vertex indices
    method : str, optional
        Method for computing mean curvature weights ('uniform', 'cotangent')
    validate : bool, optional
        Whether to validate mesh before computation

    Returns:
    --------
    curvatures : dict
        Dictionary with vertex indices as keys and curvature info as values
        Each value is a dict with keys: 'H' (mean), 'K' (Gaussian), 'valid'

    Raises:
    -------
    ValueError
        If mesh validation fails with invalid vertex indices
    """
    vertices = np.asarray(vertices)
    faces = np.asarray(faces)

    if validate:
        validation_info = validate_mesh(vertices, faces)

        if not validation_info["is_valid"]:
            raise ValueError(
                f"Invalid mesh: faces reference non-existent vertices "
                f"{validation_info['invalid_indices']}"
            )

        if validation_info["unused_vertices"]:
            warnings.warn(
                f"Mesh has unused vertices: {validation_info['unused_vertices']}"
            )

        print(
            f"Mesh validation: {validation_info['n_vertices']} vertices, "
            f"{validation_info['n_faces']} faces"
        )

    # Get vertices that are actually used in faces
    used_vertex_indices = set()
    for face in faces:
        used_vertex_indices.update(face)

    curvatures = {}

    # Compute curvatures for each used vertex
    for vertex_idx in used_vertex_indices:
        curvature_info = {"H": 0.0, "K": 0.0, "valid": False, "error": None}

        try:
            # Compute mean curvature
            H = discrete_mean_curvature(vertex_idx, vertices, faces, method=method)
            curvature_info["H"] = H

            # Compute Gaussian curvature
            K = discrete_gaussian_curvature(vertex_idx, vertices, faces)
            curvature_info["K"] = K

            curvature_info["valid"] = True

        except Exception as e:
            curvature_info["error"] = str(e)
            warnings.warn(f"Failed to compute curvature for vertex {vertex_idx}: {e}")

        curvatures[vertex_idx] = curvature_info

    return curvatures


def get_curvature_arrays(curvatures, n_vertices):
    """
    Convert curvature dictionary to arrays for easier plotting/analysis.

    Parameters:
    -----------
    curvatures : dict
        Curvature dictionary from compute_mesh_curvatures
    n_vertices : int
        Total number of vertices in mesh

    Returns:
    --------
    H_array : array, shape (n_vertices,)
        Mean curvature array (NaN for vertices without curvature)
    K_array : array, shape (n_vertices,)
        Gaussian curvature array (NaN for vertices without curvature)
    valid_array : array, shape (n_vertices,)
        Boolean array indicating which vertices have valid curvature
    """
    H_array = np.full(n_vertices, np.nan)
    K_array = np.full(n_vertices, np.nan)
    valid_array = np.zeros(n_vertices, dtype=bool)

    for vertex_idx, curvature_info in curvatures.items():
        if vertex_idx < n_vertices:  # Safety check
            H_array[vertex_idx] = curvature_info["H"]
            K_array[vertex_idx] = curvature_info["K"]
            valid_array[vertex_idx] = curvature_info["valid"]

    return H_array, K_array, valid_array


# Example usage and testing
if __name__ == "__main__":
    print("Testing Updated Curvature Functions")
    print("=" * 40)

    # Create simple test mesh (tetrahedron)
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.5, 0.5, 1.0]]
    )

    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])

    print("Test mesh: 4 vertices, 4 faces (tetrahedron)")

    # Validate mesh
    validation = validate_mesh(vertices, faces)
    print(f"Validation: {validation}")

    # Test individual functions
    print("\nTesting individual functions:")

    vertex_idx = 0
    print(f"\nVertex {vertex_idx}:")

    try:
        neighbors = find_neighbors(vertex_idx, faces)
        print(f"  Neighbors: {neighbors}")

        H = discrete_mean_curvature(vertex_idx, vertices, faces, method="uniform")
        print(f"  Mean curvature (uniform): {H:.6f}")

        H_cot = discrete_mean_curvature(vertex_idx, vertices, faces, method="cotangent")
        print(f"  Mean curvature (cotangent): {H_cot:.6f}")

        K = discrete_gaussian_curvature(vertex_idx, vertices, faces)
        print(f"  Gaussian curvature: {K:.6f}")

    except Exception as e:
        print(f"  Error: {e}")

    # Test full mesh computation
    print("\nComputing curvatures for entire mesh:")
    curvatures = compute_mesh_curvatures(vertices, faces, method="uniform")

    for vertex_idx, info in curvatures.items():
        print(
            f"Vertex {vertex_idx}: H={info['H']:.6f}, K={info['K']:.6f}, "
            f"valid={info['valid']}"
        )

    # Convert to arrays
    H_array, K_array, valid_array = get_curvature_arrays(curvatures, len(vertices))
    print(f"\nCurvature arrays:")
    print(f"H: {H_array}")
    print(f"K: {K_array}")
    print(f"Valid: {valid_array}")
