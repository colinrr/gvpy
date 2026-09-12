import pytest

from gvpy.ice_cauldron import IceCauldron


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
