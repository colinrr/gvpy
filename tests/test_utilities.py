import numpy as np
import pytest

from gvpy.ice_cauldron import IceCauldron
from gvpy.utilities.geometry import ellipse_perimeter, get_elliptical_cylinder_sa, s_spheroid
from gvpy.utilities.math_helpers import SmoothTopHat


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


class TestSmoothTopHat:
    """SmoothTopHat in utilities/math_helpers.py (new, no MATLAB source)."""

    def test_default_shape(self):
        th = SmoothTopHat()
        assert th.x_end() == 5.0
        x, y = th.get_xy()
        assert len(x) == 101 and x[0] == 0.0 and x[-1] == 5.0
        assert y[0] == 0.0  # cosine rise starts exactly at 0
        assert np.isclose(y[-1], th.tol)  # logistic tail ends at tol
        assert np.all(th(np.linspace(1, 3, 11)) == 1.0)
        assert np.all(th(np.array([-1.0, 5.1])) == 0.0)

    @pytest.mark.parametrize("kernel", ["cosine", "logistic", "smoothstep"])
    def test_kernels_monotonic_and_bounded(self, kernel):
        th = SmoothTopHat(rise_kernel=kernel, tail_kernel=kernel)
        x, y = th.get_xy(n=1001)
        i_peak = np.searchsorted(x, th.rise_length)
        assert np.all(np.diff(y[:i_peak]) >= 0) and np.all(np.diff(y[i_peak:]) <= 0)
        assert y.min() >= 0 and y.max() == 1.0
        assert y[0] <= th.tol and y[-1] <= th.tol

    @pytest.mark.parametrize("kernel", ["cosine", "logistic", "smoothstep"])
    def test_area_matches_numerical_integral(self, kernel):
        th = SmoothTopHat(rise_length=0.5, flat_length=1.5, tail_length=3.0, rise_kernel=kernel, tail_kernel=kernel)
        x, y = th.get_xy(n=200001)
        # tolerance covers trapezoid error across the logistic tol jumps
        assert np.isclose(np.trapezoid(y, x), th.get_area(), rtol=1e-4)

    def test_fixed_step(self):
        x, _ = SmoothTopHat().get_xy(dx=0.1)
        assert len(x) == 51 and np.allclose(np.diff(x), 0.1) and np.isclose(x[-1], 5.0)

    def test_unknown_kernel_raises(self):
        with pytest.raises(ValueError):
            SmoothTopHat(rise_kernel="tanh")

    def test_dimensionalize(self):
        th = SmoothTopHat(y_max=5.4e7, area=4e11)
        x, y = th.get_xy(n=200001)
        assert y.max() == 5.4e7
        assert np.isclose(x[-1], 5.0 * th.x_scale)
        assert th.get_area() == pytest.approx(4e11)
        assert np.isclose(np.trapezoid(y, x), 4e11, rtol=1e-4)
        # Non-dimensional outputs are unchanged by the dimensional parameters
        assert th.get_area(dimensionalize=False) == SmoothTopHat().get_area()
        assert np.array_equal(th.get_xy(dimensionalize=False)[1], SmoothTopHat().get_xy()[1])

    def test_dimensionalize_without_area_only_scales_y(self):
        th = SmoothTopHat(y_max=2.0)
        assert th.x_scale == 1.0 and th.get_area() == 2.0 * 3.5
        assert np.array_equal(th(np.array([0.5, 2.0])), 2.0 * SmoothTopHat()(np.array([0.5, 2.0])))

    def test_x_0_shifts_dimensional_start(self):
        th = SmoothTopHat(y_max=2.0, area=14.0, x_0=-3.0)
        base = SmoothTopHat(y_max=2.0, area=14.0)
        x, y = th.get_xy()
        x_base, y_base = base.get_xy()
        assert x[0] == -3.0 and np.allclose(x, x_base - 3.0) and np.allclose(y, y_base)
        assert th(np.array([-3.5]))[0] == 0.0
        assert th.get_area() == base.get_area()
        assert th.get_xy(dimensionalize=False)[0][0] == 0.0

    def test_x_end_dimensionalize(self):
        th = SmoothTopHat(y_max=2.0, area=14.0, x_0=-3.0)
        assert th.x_end(dimensionalize=False) == 5.0
        assert th.x_end() == -3.0 + 5.0 * th.x_scale
        assert th.get_xy()[0][-1] == th.x_end()
