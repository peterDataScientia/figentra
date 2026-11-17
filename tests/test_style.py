import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from figentra import new_figure, style_axes, legend_bold, save_figure, style_context

def test_canvas_and_axes():
    fig, ax = new_figure()
    assert tuple(round(x * 25.4, 3) for x in fig.get_size_inches()) == (130, 100)
    style_axes(ax, xlabel="Time [s]", ylabel="Power [µW]")
    assert ax.xaxis.label.get_fontsize() == 15
    assert ax.xaxis.label.get_fontweight() == "bold"
    assert ax.spines["left"].get_linewidth() == 2
    assert not ax.spines["top"].get_visible()
    plt.close(fig)

def test_export(tmp_path):
    fig, ax = new_figure()
    ax.plot([1, 2], [2, 3], label="Example")
    legend_bold(ax)
    paths = save_figure(fig, "Figure_01", output_dir=tmp_path, close=True)
    assert {p.suffix for p in paths} == {".png", ".pdf", ".svg"}
    assert all(p.stat().st_size > 100 for p in paths)

def test_invalid_filename(tmp_path):
    fig, _ = new_figure()
    with pytest.raises(ValueError):
        save_figure(fig, "../unsafe", output_dir=tmp_path)
    plt.close(fig)

def test_style_context_restores_settings():
    previous = matplotlib.rcParams["axes.linewidth"]
    with style_context():
        assert matplotlib.rcParams["axes.linewidth"] == 2
    assert matplotlib.rcParams["axes.linewidth"] == previous

def test_approved_unit_brackets():
    from figentra import approved_unit_label
    assert approved_unit_label("Concentration (uM)") == "Concentration [µM]"
    assert approved_unit_label("Response (%)") == "Response [%]"
    assert approved_unit_label("Time (s)") == "Time [s]"
    assert approved_unit_label("Signal (background)") == "Signal (background)"
    assert approved_unit_label("Heat [µcal]") == "Heat [µcal]"

def test_existing_axis_labels_normalized():
    fig, ax = new_figure()
    ax.set_xlabel("Concentration (uM)")
    ax.set_ylabel("Response (%)")
    style_axes(ax)
    assert ax.get_xlabel() == "Concentration [µM]"
    assert ax.get_ylabel() == "Response [%]"
    plt.close(fig)
