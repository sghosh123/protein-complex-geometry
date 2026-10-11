import numpy as np
import pytest
from protein_complex_geometry.geometry import (
    RigidTransform, euler_rotation_matrix, kabsch_transform, rmsd,
    rotation_matrix, transform_from_euler,
)

def test_rotation_about_z_follows_right_hand_rule():
    out = np.array([[1., 0., 0.]]) @ rotation_matrix("z", 90)
    np.testing.assert_allclose(out, [[0., 1., 0.]], atol=1e-7)

def test_rotation_preserves_pairwise_distance(noncollinear_points):
    before = np.linalg.norm(noncollinear_points[0] - noncollinear_points[1])
    after_xyz = noncollinear_points @ euler_rotation_matrix(19, -32, 73)
    after = np.linalg.norm(after_xyz[0] - after_xyz[1])
    assert after == pytest.approx(before)

def test_transform_rotates_about_pivot_then_translates():
    transform = transform_from_euler(z_degrees=90, origin=(1, 0, 0), translation=(2, 3, 4))
    np.testing.assert_allclose(transform.apply([[1, 0, 0]]), [[3, 3, 4]], atol=1e-7)
    np.testing.assert_allclose(transform.apply([[2, 0, 0]]), [[3, 4, 4]], atol=1e-7)

def test_kabsch_recovers_known_transform(noncollinear_points):
    rotation = euler_rotation_matrix(17, 42, -91)
    translation = np.array([12.5, -3.0, 7.25])
    mobile = noncollinear_points
    reference = mobile @ rotation + translation
    fit, final_rmsd = kabsch_transform(mobile, reference)
    np.testing.assert_allclose(fit.apply(mobile), reference, atol=1e-7)
    assert final_rmsd < 1e-7
    assert np.linalg.det(fit.rotation) == pytest.approx(1.0)

def test_kabsch_rejects_reflection_fit_as_impossible_to_match_exactly(noncollinear_points):
    reflected = noncollinear_points.copy(); reflected[:, 0] *= -1
    transform, after = kabsch_transform(noncollinear_points, reflected)
    assert np.linalg.det(transform.rotation) == pytest.approx(1.0)
    assert after > 0.01

def test_kabsch_rejects_shape_mismatch(noncollinear_points):
    with pytest.raises(ValueError, match="shapes differ"):
        kabsch_transform(noncollinear_points, noncollinear_points[:3])

def test_kabsch_requires_three_non_collinear_points():
    with pytest.raises(ValueError, match="at least three"):
        kabsch_transform([[0,0,0],[1,0,0]], [[0,0,0],[1,0,0]])
    with pytest.raises(ValueError, match="collinear"):
        kabsch_transform([[0,0,0],[1,0,0],[2,0,0]], [[0,0,0],[1,0,0],[2,0,0]])

def test_rmsd_correct():
    assert rmsd([[0,0,0],[1,0,0]], [[0,0,0],[2,0,0]]) == pytest.approx(np.sqrt(0.5))

def test_transform_rejects_improper_rotation():
    with pytest.raises(ValueError, match="proper rotation"):
        RigidTransform(np.diag([-1., 1., 1.]), np.zeros(3))
