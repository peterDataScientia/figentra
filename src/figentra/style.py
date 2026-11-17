"""Approved figure settings extracted from SNT100_MASTER_STYLE_1_TO_11.ipynb.

Only presentation is changed; scientific observations and fits are untouched.
"""
from contextlib import contextmanager
from pathlib import Path
import re
import matplotlib as mpl
import matplotlib.pyplot as plt

APPROVED_STYLE = {
    "size_mm": (130, 100), "font_family": "DejaVu Sans",
    "base_font_pt": 12, "axis_label_pt": 15, "tick_label_pt": 12,
    "legend_pt": 11, "axes_width_pt": 2.0, "tick_width_pt": 2.0,
    "tick_length_pt": 6, "label_pad": 7, "dpi": 600,
    "formats": ("png", "pdf", "svg"),
}
_RC = {
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.labelsize": 15, "axes.labelweight": "bold",
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "axes.linewidth": 2.0, "pdf.fonttype": 42,
    "svg.fonttype": "none", "savefig.facecolor": "white",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.grid": False,
}

def apply_style():
    """Apply approved Matplotlib defaults globally."""
    mpl.rcParams.update(_RC)

@contextmanager
def style_context():
    """Apply the approved style temporarily, without global leakage."""
    with mpl.rc_context(rc=_RC):
        yield

def new_figure(*, size_mm=None, **kwargs):
    """Create figure and axes using the approved 130 × 100 mm canvas."""
    apply_style()
    w, h = size_mm or APPROVED_STYLE["size_mm"]
    fig, ax = plt.subplots(figsize=(w / 25.4, h / 25.4), **kwargs)
    style_axes(ax)
    return fig, ax

_UNIT_ALIASES = {
    "%": "%", "uM": "µM", "μM": "µM", "µM": "µM",
    "nM": "nM", "mM": "mM", "M": "M", "pM": "pM",
    "s": "s", "ms": "ms", "min": "min", "h": "h",
    "µs": "µs", "us": "µs", "μs": "µs",
    "RU": "RU", "a.u.": "a.u.", "AU": "AU",
    "°C": "°C", "K": "K", "µW": "µW", "uW": "µW",
    "kcal/mol": "kcal/mol", "Å": "Å", "nm": "nm",
}
_TRAILING_UNIT = re.compile(r"^(.*?)\s*\(([^()]*)\)\s*$")

def approved_unit_label(value):
    """Convert recognized trailing unit parentheses to approved square brackets.

    E.g. 'Concentration (uM)' -> 'Concentration [µM]'.
    Unrecognized parenthetical phrases are left untouched.
    """
    if value is None:
        return None
    match = _TRAILING_UNIT.fullmatch(str(value))
    if not match:
        return value
    label, unit = match.groups()
    normalized = _UNIT_ALIASES.get(unit.strip())
    if normalized is None:
        return value
    return f"{label.rstrip()} [{normalized}]"

def style_axes(ax, *, xlabel=None, ylabel=None):
    """Enforce approved spines, ticks, labels, and white background."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(2)
    ax.tick_params(axis="both", which="major", width=2,
                   length=6, direction="out", labelsize=12)
    ax.minorticks_off()
    ax.grid(False)
    ax.set_facecolor("white")
    # Normalize requested labels and labels previously set on the axes.
    ax.set_xlabel(approved_unit_label(
        xlabel if xlabel is not None else ax.get_xlabel()))
    ax.set_ylabel(approved_unit_label(
        ylabel if ylabel is not None else ax.get_ylabel()))
    for label in (ax.xaxis.label, ax.yaxis.label):
        label.set_fontsize(15)
        label.set_fontweight("bold")
    ax.xaxis.labelpad = 7
    ax.yaxis.labelpad = 7
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(12)
        label.set_fontweight("bold")
    return ax

def legend_bold(ax, **kwargs):
    """Apply the approved 11-point, bold, frameless legend."""
    kwargs.setdefault("frameon", False)
    kwargs.setdefault("prop", {"weight": "bold", "size": 11})
    return ax.legend(**kwargs)

def save_figure(fig, name, *, output_dir=".", formats=None, dpi=None, close=False):
    """Save PNG/PDF/SVG; returns list of Path objects."""
    name = str(name)
    if not name or Path(name).name != name or Path(name).suffix:
        raise ValueError("name must be a filename stem without directories or extension")
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    formats = tuple(formats or APPROVED_STYLE["formats"])
    if not formats or any(f not in ("png", "pdf", "svg") for f in formats):
        raise ValueError("formats must contain png, pdf or svg")
    saved = []
    for fmt in formats:
        path = folder / (name + "." + fmt)
        fig.savefig(path, dpi=dpi or APPROVED_STYLE["dpi"], facecolor="white")
        saved.append(path)
    if close:
        plt.close(fig)
    return saved
