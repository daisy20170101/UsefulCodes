"""
Fault curvature analysis from triangular mesh data.

This module computes discrete differential geometry measures on triangular meshes,
particularly designed for fault surface analysis in geoscience applications.

Key curvature indices:
- Gaussian curvature (K): Identifies synforms, saddles, and planar regions
- Mean curvature (H): Measures overall surface bending and roughness
- Shape Index (SI): Classifies surface morphology independent of scale
"""

import numpy as np
from typing import Tuple, Dict
import warnings


def compute_vertex_normals(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """
    Compute vertex normals using area-weighted face normals.

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face

    Returns
    -------
    normals : np.ndarray, shape (n_vertices, 3)
        Unit normal vectors at each vertex
    """
    n_verts = len(vertices)
    normals = np.zeros((n_verts, 3), dtype=float)

    # Compute face normals weighted by face area
    for face in faces:
        v0, v1, v2 = vertices[face]
        # Two edge vectors
        e1 = v1 - v0
        e2 = v2 - v0
        # Cross product gives area-weighted normal
        face_normal = np.cross(e1, e2)
        # Accumulate to vertices
        for idx in face:
            normals[idx] += face_normal

    # Normalize
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-12)  # Avoid division by zero
    normals = normals / norms

    return normals


def compute_vertex_areas(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """
    Compute the Voronoi area associated with each vertex.
    Uses mixed area (1/3 of incident triangle areas) as approximation.

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face

    Returns
    -------
    areas : np.ndarray, shape (n_vertices,)
        Area associated with each vertex
    """
    n_verts = len(vertices)
    areas = np.zeros(n_verts, dtype=float)

    for face in faces:
        v0, v1, v2 = vertices[face]
        # Triangle area
        e1 = v1 - v0
        e2 = v2 - v0
        area = 0.5 * np.linalg.norm(np.cross(e1, e2))
        # Distribute 1/3 to each vertex
        for idx in face:
            areas[idx] += area / 3.0

    return areas


def compute_principal_curvatures(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute principal curvatures κ₁ and κ₂ at each vertex using discrete operators.

    Uses the Meyer et al. (2003) discrete differential geometry operators.
    Reference: "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds"

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")
        - "meyer": Meyer et al. cotangent formula
        - "simple": Simple averaging (faster but less accurate)

    Returns
    -------
    kappa_1 : np.ndarray, shape (n_vertices,)
        Maximum principal curvature at each vertex
    kappa_2 : np.ndarray, shape (n_vertices,)
        Minimum principal curvature at each vertex
    dir_1 : np.ndarray, shape (n_vertices, 3)
        Direction of maximum curvature (principal direction 1)
    dir_2 : np.ndarray, shape (n_vertices, 3)
        Direction of minimum curvature (principal direction 2)
    """
    n_verts = len(vertices)

    # Compute vertex normals
    normals = compute_vertex_normals(vertices, faces)

    # Compute vertex areas
    areas = compute_vertex_areas(vertices, faces)

    # Build adjacency information
    # For each vertex, store neighboring vertices
    neighbors = [set() for _ in range(n_verts)]
    for face in faces:
        neighbors[face[0]].update([face[1], face[2]])
        neighbors[face[1]].update([face[0], face[2]])
        neighbors[face[2]].update([face[0], face[1]])

    # Compute mean curvature normal using cotangent formula
    mean_curvature_normal = np.zeros((n_verts, 3), dtype=float)

    if method == "meyer":
        # Meyer et al. cotangent formula
        for i in range(n_verts):
            vi = vertices[i]
            for j in neighbors[i]:
                vj = vertices[j]
                edge = vj - vi

                # Find triangles containing edge (i,j)
                cotangent_sum = 0.0
                for face in faces:
                    if i in face and j in face:
                        # Find the opposite vertex
                        k = [v for v in face if v != i and v != j][0]
                        vk = vertices[k]

                        # Compute cotangent of angle at k
                        e1 = vi - vk
                        e2 = vj - vk
                        cos_angle = np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2) + 1e-12)
                        # Clamp to avoid numerical issues
                        cos_angle = np.clip(cos_angle, -1.0, 1.0)
                        sin_angle = np.sqrt(1.0 - cos_angle**2)
                        cot = cos_angle / (sin_angle + 1e-12)
                        cotangent_sum += cot

                mean_curvature_normal[i] += cotangent_sum * edge

        # Normalize by 2*area
        mean_curvature_normal /= (2.0 * areas[:, np.newaxis] + 1e-12)

    else:  # simple method
        for i in range(n_verts):
            vi = vertices[i]
            for j in neighbors[i]:
                vj = vertices[j]
                mean_curvature_normal[i] += (vj - vi)
            mean_curvature_normal[i] /= (len(neighbors[i]) + 1e-12)

    # Mean curvature H (scalar)
    H = np.einsum('ij,ij->i', mean_curvature_normal, normals) / 2.0

    # Gaussian curvature K using angle defect formula
    K = np.zeros(n_verts, dtype=float)
    for i in range(n_verts):
        angle_sum = 0.0
        # Sum angles at vertex i in all incident triangles
        for face in faces:
            if i in face:
                # Get vertices of the triangle
                v_indices = list(face)
                i_idx = v_indices.index(i)
                j_idx = (i_idx + 1) % 3
                k_idx = (i_idx + 2) % 3

                vi = vertices[v_indices[i_idx]]
                vj = vertices[v_indices[j_idx]]
                vk = vertices[v_indices[k_idx]]

                # Compute angle at vi
                e1 = vj - vi
                e2 = vk - vi
                cos_angle = np.dot(e1, e2) / (np.linalg.norm(e1) * np.linalg.norm(e2) + 1e-12)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.arccos(cos_angle)
                angle_sum += angle

        # Angle defect formula: K = (2π - sum of angles) / area
        K[i] = (2.0 * np.pi - angle_sum) / (areas[i] + 1e-12)

    # Compute principal curvatures from H and K
    # κ₁ = H + sqrt(H² - K)
    # κ₂ = H - sqrt(H² - K)
    discriminant = H**2 - K
    discriminant = np.maximum(discriminant, 0.0)  # Handle numerical errors
    sqrt_disc = np.sqrt(discriminant)

    kappa_1 = H + sqrt_disc
    kappa_2 = H - sqrt_disc

    # Principal directions (simplified - along tangent plane)
    # For a complete implementation, would need to compute shape operator eigenvectors
    # Here we provide approximate directions orthogonal to normal
    dir_1 = np.zeros((n_verts, 3), dtype=float)
    dir_2 = np.zeros((n_verts, 3), dtype=float)

    for i in range(n_verts):
        n = normals[i]
        # Choose arbitrary tangent direction
        if abs(n[0]) < 0.9:
            t1 = np.array([1.0, 0.0, 0.0])
        else:
            t1 = np.array([0.0, 1.0, 0.0])
        # Project to tangent plane
        t1 = t1 - np.dot(t1, n) * n
        t1 = t1 / (np.linalg.norm(t1) + 1e-12)
        # Second direction orthogonal to both
        t2 = np.cross(n, t1)
        t2 = t2 / (np.linalg.norm(t2) + 1e-12)

        dir_1[i] = t1
        dir_2[i] = t2

    return kappa_1, kappa_2, dir_1, dir_2


def gaussian_curvature(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> np.ndarray:
    """
    Compute Gaussian curvature K = κ₁ × κ₂ at each vertex.

    Gaussian curvature characterizes the intrinsic geometry:
    - K > 0: Elliptic points (bowl-like, synform)
    - K < 0: Hyperbolic points (saddle-like) - critical for fault mechanics
    - K ≈ 0: Parabolic points (cylindrical or planar)

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")

    Returns
    -------
    K : np.ndarray, shape (n_vertices,)
        Gaussian curvature at each vertex
    """
    kappa_1, kappa_2, _, _ = compute_principal_curvatures(vertices, faces, method)
    K = kappa_1 * kappa_2
    return K


def mean_curvature(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> np.ndarray:
    """
    Compute mean curvature H = (κ₁ + κ₂)/2 at each vertex.

    Mean curvature measures overall surface bending:
    - H > 0: Convex regions (bulges/asperities)
    - H < 0: Concave regions (valleys)
    - |H| magnitude: Intensity of bending

    Useful for identifying fault surface roughness and asperities.

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")

    Returns
    -------
    H : np.ndarray, shape (n_vertices,)
        Mean curvature at each vertex
    """
    kappa_1, kappa_2, _, _ = compute_principal_curvatures(vertices, faces, method)
    H = (kappa_1 + kappa_2) / 2.0
    return H


def shape_index(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> np.ndarray:
    """
    Compute Shape Index SI = (2/π) × arctan((κ₁ + κ₂)/(κ₁ - κ₂)).

    Shape index provides scale-independent classification of surface morphology:
    - SI = -1.0: Spherical valley (cup)
    - SI = -0.5: Trough (valley)
    - SI =  0.0: Saddle (critical for fault analysis!)
    - SI = +0.5: Ridge
    - SI = +1.0: Spherical dome (cap)

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")

    Returns
    -------
    SI : np.ndarray, shape (n_vertices,)
        Shape index at each vertex, range [-1, 1]
    """
    kappa_1, kappa_2, _, _ = compute_principal_curvatures(vertices, faces, method)

    # Handle numerical issues when κ₁ ≈ κ₂ (planar regions)
    epsilon = 1e-10
    denominator = kappa_1 - kappa_2
    denominator = np.where(np.abs(denominator) < epsilon, epsilon, denominator)

    numerator = kappa_1 + kappa_2

    # Compute shape index
    SI = (2.0 / np.pi) * np.arctan(numerator / denominator)

    return SI


def curvedness(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> np.ndarray:
    """
    Compute curvedness C = sqrt((κ₁² + κ₂²)/2).

    Curvedness is a scale-independent measure of "how much" the surface bends,
    independent of the shape type (complementary to shape index).

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")

    Returns
    -------
    C : np.ndarray, shape (n_vertices,)
        Curvedness at each vertex (always ≥ 0)
    """
    kappa_1, kappa_2, _, _ = compute_principal_curvatures(vertices, faces, method)
    C = np.sqrt((kappa_1**2 + kappa_2**2) / 2.0)
    return C


def compute_all_curvatures(
    vertices: np.ndarray,
    faces: np.ndarray,
    method: str = "meyer"
) -> Dict[str, np.ndarray]:
    """
    Compute all curvature measures at once (more efficient than separate calls).

    Parameters
    ----------
    vertices : np.ndarray, shape (n_vertices, 3)
        3D coordinates of mesh vertices
    faces : np.ndarray, shape (n_faces, 3)
        Indices of vertices forming each triangular face
    method : str
        Method for curvature computation (default: "meyer")

    Returns
    -------
    curvatures : dict
        Dictionary containing:
        - 'kappa_1': Maximum principal curvature
        - 'kappa_2': Minimum principal curvature
        - 'gaussian': Gaussian curvature K
        - 'mean': Mean curvature H
        - 'shape_index': Shape index SI
        - 'curvedness': Curvedness C
        - 'dir_1': Principal direction 1 (shape (n_vertices, 3))
        - 'dir_2': Principal direction 2 (shape (n_vertices, 3))
    """
    kappa_1, kappa_2, dir_1, dir_2 = compute_principal_curvatures(vertices, faces, method)

    # Gaussian curvature
    K = kappa_1 * kappa_2

    # Mean curvature
    H = (kappa_1 + kappa_2) / 2.0

    # Shape index
    epsilon = 1e-10
    denominator = kappa_1 - kappa_2
    denominator = np.where(np.abs(denominator) < epsilon, epsilon, denominator)
    SI = (2.0 / np.pi) * np.arctan((kappa_1 + kappa_2) / denominator)

    # Curvedness
    C = np.sqrt((kappa_1**2 + kappa_2**2) / 2.0)

    return {
        'kappa_1': kappa_1,
        'kappa_2': kappa_2,
        'gaussian': K,
        'mean': H,
        'shape_index': SI,
        'curvedness': C,
        'dir_1': dir_1,
        'dir_2': dir_2,
    }


# ============================================================================
# Visualization helpers
# ============================================================================

def classify_surface_type(SI: np.ndarray) -> np.ndarray:
    """
    Classify surface type based on shape index value.

    Parameters
    ----------
    SI : np.ndarray
        Shape index values

    Returns
    -------
    types : np.ndarray of str
        Surface type classification
    """
    types = np.empty(len(SI), dtype=object)

    types[SI < -0.75] = 'spherical_valley'
    types[(SI >= -0.75) & (SI < -0.25)] = 'trough'
    types[(SI >= -0.25) & (SI < 0.25)] = 'saddle'
    types[(SI >= 0.25) & (SI < 0.75)] = 'ridge'
    types[SI >= 0.75] = 'spherical_dome'

    return types


def curvature_statistics(curvatures: np.ndarray, name: str = "Curvature") -> Dict:
    """
    Compute summary statistics for curvature values.

    Parameters
    ----------
    curvatures : np.ndarray
        Array of curvature values
    name : str
        Name of the curvature measure

    Returns
    -------
    stats : dict
        Dictionary with statistical measures
    """
    return {
        'name': name,
        'mean': np.mean(curvatures),
        'median': np.median(curvatures),
        'std': np.std(curvatures),
        'min': np.min(curvatures),
        'max': np.max(curvatures),
        'q25': np.percentile(curvatures, 25),
        'q75': np.percentile(curvatures, 75),
    }


# ============================================================================
# Example usage
# ============================================================================

if __name__ == "__main__":
    # Example: Simple mesh (4 vertices forming 2 triangles - a bent plane)
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.1],
        [0.0, 1.0, 0.1],
        [1.0, 1.0, 0.2],
    ])

    faces = np.array([
        [0, 1, 2],
        [1, 3, 2],
    ])

    print("Computing curvatures for example mesh...")
    print(f"Vertices shape: {vertices.shape}")
    print(f"Faces shape: {faces.shape}")

    # Compute all curvatures
    curv = compute_all_curvatures(vertices, faces)

    print("\nResults:")
    print(f"Gaussian curvature (K): {curv['gaussian']}")
    print(f"Mean curvature (H): {curv['mean']}")
    print(f"Shape index (SI): {curv['shape_index']}")
    print(f"Curvedness (C): {curv['curvedness']}")

    # Statistics
    print("\nStatistics:")
    for key in ['gaussian', 'mean', 'shape_index']:
        stats = curvature_statistics(curv[key], key)
        print(f"{key}: mean={stats['mean']:.6f}, std={stats['std']:.6f}")
