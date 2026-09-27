import numpy as np
import pytest

from gvpy.ice_cauldron import IceCauldron
from gvpy.utilities.geometry import ellipse_perimeter, get_elliptical_cylinder_sa, s_spheroid


class TestValidators:
    """attrs.field(validator=...) checks on IceCauldron, translated from
    IceCauldron.m's declarative {mustBeX} property validators.
    """

    def test_default_construction(self):
        c = IceCauldron()
        assert c.params is not None
        assert c.E_prime is not None
        assert c.fastest_cauldron_index is not None

    def test_geometry_must_be_member(self):
        with pytest.raises(ValueError):
            IceCauldron(geometry="not-a-shape")

    def test_f_i_must_be_in_range_inclusive(self):
        with pytest.raises(ValueError):
            IceCauldron(f_i=1.5)

    def test_A_must_be_positive(self):
        with pytest.raises(ValueError):
            IceCauldron(A=-1.0)

    def test_G_n_must_be_nonnegative(self):
        with pytest.raises(ValueError):
            IceCauldron(G_n=-5.0)

    def test_on_setattr_revalidates(self):
        """attrs.define validates on every subsequent assignment by default,
        not just at construction."""
        c = IceCauldron()
        with pytest.raises(ValueError):
            c.f_i = 2.0


class TestGeometry:
    """ellipse_perimeter, get_elliptical_cylinder_sa and s_spheroid (from
    utilities/geometry/S_spheroid.m, with its edge cases corrected)."""

    def test_ellipse_perimeter_circle(self):
        assert ellipse_perimeter(3.0, 3.0) == pytest.approx(2 * np.pi * 3.0)

    def test_elliptical_cylinder_sa_uses_perimeter(self):
        assert get_elliptical_cylinder_sa(4.0, 2.0, l=5) == pytest.approx(ellipse_perimeter(4.0, 2.0) * 5 / 2)
        assert get_elliptical_cylinder_sa(4.0, 2.0, halved=False) == pytest.approx(ellipse_perimeter(4.0, 2.0))

    def test_sphere_and_disc(self):
        a = 100.0
        assert s_spheroid(a, a, halved=False)[0] == pytest.approx(4 * np.pi * a**2)
        assert s_spheroid(a, a)[0] == pytest.approx(2 * np.pi * a**2)
        assert s_spheroid(a, 0.0)[0] == pytest.approx(np.pi * a**2)  # MATLAB's c == 0 value

    def test_limits_are_continuous(self):
        a = 100.0
        assert s_spheroid(a, a * (1 - 1e-6))[0] == pytest.approx(s_spheroid(a, a)[0], rel=1e-5)  # oblate -> sphere
        assert s_spheroid(a, a * (1 + 1e-6))[0] == pytest.approx(s_spheroid(a, a)[0], rel=1e-5)  # prolate -> sphere
        assert s_spheroid(a, 1e-6)[0] == pytest.approx(s_spheroid(a, 0.0)[0], rel=1e-5)  # oblate -> disc

    def test_mixed_vector_matches_elementwise(self):
        a = np.array([100.0, 100.0, 100.0, 100.0])
        c = np.array([0.0, 50.0, 100.0, 150.0])
        S, P = s_spheroid(a, c)
        for i in range(len(a)):
            S_i, P_i = s_spheroid(a[i], c[i])
            assert S[i] == pytest.approx(S_i)
            assert P[i] == pytest.approx(P_i)

    def test_roof_curve(self):
        # Scalar a, c (errored in MATLAB) and vector a with scalar c
        for a, c in [(100.0, 50.0), (np.array([100.0, 200.0]), 50.0)]:
            _, _, x, z = s_spheroid(a, c, n=51, return_xz=True)
            a_vec = np.broadcast_to(a, (x.shape[1],))
            assert x.shape == z.shape == (51, np.size(a))
            np.testing.assert_allclose(x**2 / a_vec**2 + z**2 / c**2, 1.0)
