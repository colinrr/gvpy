import attrs
import pytest

from gvpy.gv_main import gv_main
from gvpy.ice_cauldron import IceCauldron, MissingVarInfoWarning

REQUIRED_KEYS = {"units", "long_name", "symbol"}


def test_every_field_has_complete_var_info():
    for f in attrs.fields(IceCauldron):
        assert set(IceCauldron.VAR_INFO.get(f.name, {})) == REQUIRED_KEYS, f.name


def test_solver_var_tuples_have_var_info():
    for name in IceCauldron.vector_solution_vars + IceCauldron.vector_derived_vars:
        assert name in IceCauldron.VAR_INFO, name


def test_glen_A_units_follow_glen_n():
    assert IceCauldron().get_var_info("A")["units"] == "Pa^(-3) s^-1"
    assert IceCauldron(glen_n=4).get_var_info("A")["units"] == "Pa^(-4) s^-1"


def test_missing_entry_warns():
    c = IceCauldron()
    ds = c.params.assign(not_a_registered_var=c.params["G_n"] * 2)
    with pytest.warns(MissingVarInfoWarning, match="not_a_registered_var"):
        c.add_var_attrs(ds)


@pytest.mark.parametrize("n_cauldrons, G_n", [(1, 400), (2, [300, 400])])
def test_outputs_carry_units_and_long_name(n_cauldrons, G_n):
    # A MissingVarInfoWarning anywhere in the run fails this test (pytest filterwarnings config)
    result = gv_main(n_cauldrons=n_cauldrons, G_n=G_n, disable_plume_flux=True)
    for ds in (result.gv.params, result.data, result.events):
        for name in list(ds.data_vars) + list(ds.coords):
            info = IceCauldron.VAR_INFO[name]
            assert ds[name].attrs.get("long_name") == info["long_name"], name
            if info["units"] is not None:
                assert "units" in ds[name].attrs, name
