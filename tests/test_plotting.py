import warnings

import attrs
from bokeh.models import Plot
import holoviews as hv
import panel as pn
import pytest

from gvpy.gv_main import gv_main
from gvpy.plotting import plot_cauldron_geometry, plot_results, plot_time_series


@pytest.fixture(scope="module")
def r1():
    return gv_main(n_cauldrons=1, G_n=400, disable_plume_flux=True, supraglacial_drainage_mode="off")


@pytest.fixture(scope="module")
def r2():
    return gv_main(n_cauldrons=2, G_n=[300, 400], disable_plume_flux=True)


def _plots(layout):
    return list(layout.get_root().select({"type": Plot}))


def _markdown_text(layout):
    return " ".join(obj.object for obj in layout.select(pn.pane.Markdown))


@pytest.mark.parametrize("which", ["r1", "r2", "both"])
def test_full_layout_renders_linked_without_warnings(which, r1, r2):
    results = {"r1": r1, "r2": r2, "both": {"single": r1, "double": r2}}[which]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        layout = plot_results(results)
        plots = _plots(layout)
    assert len(plots) == 12  # 14 panels minus 2 placeholders
    x_groups = sorted([sum(p.x_range is q.x_range for q in plots) for p in plots], reverse=True)
    y_groups = sorted([sum(p.y_range is q.y_range for q in plots) for p in plots], reverse=True)
    assert x_groups[:11] == [11] * 11 and x_groups[11] == 1  # all time panels linked; geometry separate
    assert y_groups[:3] == [3] * 3  # the three plotted volume panels share y
    assert "Sheet Volume" in _markdown_text(layout) and "Plume Flux" in _markdown_text(layout)


def test_panel_subset(r1):
    layout = plot_results(r1, panels=["Cauldron Dimensions", "Supraglacial Flux"])
    assert len(_plots(layout)) == 2


def test_event_markers(r2):
    overlay = plot_time_series(r2, "a_n")
    markers = [el for el in overlay if isinstance(el, hv.Scatter)][0]
    assert len(markers) == r2.events.sizes["event"]  # every event here is per-cauldron


def test_missing_long_name_raises(r1):
    data = r1.data.copy()
    data["a_n"].attrs = {}
    with pytest.raises(ValueError, match="long_name"):
        plot_time_series(attrs.evolve(r1, data=data), "a_n")


def test_unlisted_units_warn(r1):
    data = r1.data.copy()
    data["a_n"].attrs = {**data["a_n"].attrs, "units": "furlong"}
    with pytest.warns(UserWarning, match="furlong"):
        plot_time_series(attrs.evolve(r1, data=data), "a_n")


def test_missing_panel_variable_placeholder(r1):
    result = attrs.evolve(r1, data=r1.data.drop_vars("f_i_n"))
    with pytest.warns(UserWarning, match="Heating Efficiency"):
        layout = plot_results(result)
    assert "missing from results: f_i_n" in _markdown_text(layout)


def test_geometry_single_closed_snapshot(r1):
    # One closed-cauldron snapshot: scalar a, c into s_spheroid (errored in MATLAB)
    overlay = plot_cauldron_geometry(r1, steps=[5])
    hv.render(overlay, backend="bokeh")
