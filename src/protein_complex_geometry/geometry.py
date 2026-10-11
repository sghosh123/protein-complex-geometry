"""Rigid-body geometry operations using row-vector coordinates and Angstroms."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import math
import numpy as np

Array = np.ndarray


def _coords(values: Iterable[Iterable[float]], name: str = "coordinates") -> Array:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_atoms, 3); got {arr.shape}")
    if arr.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one coordinate")
    if not np.isfinite(arr).all():
        raise ValueError(f"{name} contain non-finite values")
    return arr


def rotation_matrix(axis: str, angle_degrees: float) -> Array:
    """Return a right-handed active rotation matrix for row-vector coordinates.

    Coordinates are transformed as ``xyz @ R``. Positive angles follow the
    right-hand rule. This API deliberately accepts only Cartesian x/y/z axes.
    """
    if axis.lower() not in {"x", "y", "z"}:
        raise ValueError("axis must be one of: x, y, z")
    if not math.isfinite(float(angle_degrees)):
        raise ValueError("angle_degrees must be finite")
    theta = math.radians(float(angle_degrees))
    c, s = math.cos(theta), math.sin(theta)
    if axis.lower() == "x":
        return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])
    if axis.lower() == "y":
        return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def euler_rotation_matrix(x_degrees: float = 0.0, y_degrees: float = 0.0,
                          z_degrees: float = 0.0) -> Array:
    """Compose x, then y, then z rotations for row-vector coordinates."""
    return (rotation_matrix("x", x_degrees)
            @ rotation_matrix("y", y_degrees)
            @ rotation_matrix("z", z_degrees))


@dataclass(frozen=True)
class RigidTransform:
    """Rigid transform for row-vector coordinates: ``xyz @ rotation + translation``."""
    rotation: Array
    translation: Array

    def __post_init__(self) -> None:
        r = np.asarray(self.rotation, dtype=float)
        t = np.asarray(self.translation, dtype=float)
        if r.shape != (3, 3):
            raise ValueError("rotation must have shape (3, 3)")
        if t.shape != (3,):
            raise ValueError("translation must have shape (3,)")
        if not np.isfinite(r).all() or not np.isfinite(t).all():
            raise ValueError("transform contains non-finite values")
        if not np.allclose(r.T @ r, np.eye(3), atol=1e-7):
            raise ValueError("rotation matrix must be orthonormal")
        if not np.isclose(np.linalg.det(r), 1.0, atol=1e-7):
            raise ValueError("rotation matrix must be a proper rotation (determinant +1)")
        object.__setattr__(self, "rotation", r)
        object.__setattr__(self, "translation", t)

    def apply(self, coordinates: Iterable[Iterable[float]]) -> Array:
        xyz = _coords(coordinates)
        return xyz @ self.rotation + self.translation

    def as_dict(self) -> dict:
        return {"convention": "row_vector: xyz @ rotation + translation",
                "rotation": self.rotation.tolist(), "translation_angstrom": self.translation.tolist()}


def transform_from_euler(x_degrees: float = 0.0, y_degrees: float = 0.0,
                          z_degrees: float = 0.0,
                          origin: Iterable[float] = (0.0, 0.0, 0.0),
                          translation: Iterable[float] = (0.0, 0.0, 0.0)) -> RigidTransform:
    """Create an Euler rotation about ``origin``, followed by a translation.

    ``origin`` and ``translation`` are in Angstroms. The rotation order is x,
    then y, then z. For row-vector coordinates, the affine offset is
    ``origin - origin @ rotation + translation``.
    """
    pivot = np.asarray(origin, dtype=float)
    shift = np.asarray(translation, dtype=float)
    if pivot.shape != (3,) or shift.shape != (3,):
        raise ValueError("origin and translation must each have three components")
    if not np.isfinite(pivot).all() or not np.isfinite(shift).all():
        raise ValueError("origin and translation must be finite")
    r = euler_rotation_matrix(x_degrees, y_degrees, z_degrees)
    return RigidTransform(r, pivot - pivot @ r + shift)


def kabsch_transform(mobile: Iterable[Iterable[float]],
                     reference: Iterable[Iterable[float]]) -> tuple[RigidTransform, float]:
    """Fit ``mobile`` onto ``reference`` using the Kabsch algorithm.

    The arrays must contain corresponding points in the same order and units.
    Returns the rigid transform and post-fit RMSD in the input coordinate unit.
    Reflection is explicitly prohibited. At least three non-collinear points
    are required for a well-defined 3D fit.
    """
    mob = _coords(mobile, "mobile coordinates")
    ref = _coords(reference, "reference coordinates")
    if mob.shape != ref.shape:
        raise ValueError(f"mobile and reference shapes differ: {mob.shape} vs {ref.shape}")
    if mob.shape[0] < 3:
        raise ValueError("at least three corresponding points are required")
    mob_ctr = mob.mean(axis=0)
    ref_ctr = ref.mean(axis=0)
    p = mob - mob_ctr
    q = ref - ref_ctr
    if np.linalg.matrix_rank(p, tol=1e-8) < 2 or np.linalg.matrix_rank(q, tol=1e-8) < 2:
        raise ValueError("alignment points are collinear or degenerate")
    covariance = p.T @ q
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rot = u @ correction @ vt
    translation = ref_ctr - mob_ctr @ rot
    transform = RigidTransform(rot, translation)
    fitted = transform.apply(mob)
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - ref) ** 2, axis=1))))
    return transform, rmsd


def rmsd(a: Iterable[Iterable[float]], b: Iterable[Iterable[float]]) -> float:
    """Calculate unweighted coordinate RMSD for corresponding points."""
    aa, bb = _coords(a, "a"), _coords(b, "b")
    if aa.shape != bb.shape:
        raise ValueError(f"coordinate shapes differ: {aa.shape} vs {bb.shape}")
    return float(np.sqrt(np.mean(np.sum((aa - bb) ** 2, axis=1))))
