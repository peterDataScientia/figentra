# Figentra

**Scientific figures. Consistently refined.**

Figentra is a lightweight Matplotlib styling package that reproduces the approved scientific figure settings from the `SNT100_MASTER_STYLE_1_TO_11.ipynb` research notebook.

**Figentra only changes presentation**. It does not modify data, fit models, smooth curves, normalize measurements or fabricate experimental results.

## Installation

Install from the official GitHub repository today:

```bash
pip install git+https://github.com/peterDataScientia/figentra.git
```

For local development:

```bash
git clone https://github.com/peterDataScientia/figentra.git
cd figentra
pip install -e ".[dev]"
pytest -q
```

The `pip install figentra` command will become available **after** the first successful PyPI release.

## Usage

```python
import numpy as np
from figentra import new_figure, style_axes, legend_bold, save_figure

x = np.linspace(0, 120, 200)
y = np.exp(-x / 40)  # illustration, not experimental data

fig, ax = new_figure()
ax.plot(x, y, lw=2, label="Example")
style_axes(ax, xlabel="Time [s]", ylabel="Response [a.u.]")
legend_bold(ax)
fig.subplots_adjust(left=.22, right=.92, bottom=.15, top=.81)
save_figure(fig, "Figure_01", output_dir="figures", close=True)
```

To style an existing Matplotlib figure:

```python
from figentra import apply_style, style_axes
apply_style()
# draw with Matplotlib
style_axes(ax)
```

For temporary style changes without altering global configuration, use `with style_context():`.

## Approved visual standards

| Setting | Value |
|---|---|
| Canvas | 130 × 100 mm |
| Font | DejaVu Sans |
| Base font | 12 pt |
| Axis labels | 15 pt, bold |
| Tick labels | 12 pt, bold |
| Legend | 11 pt, bold, no frame |
| Axes | Left/bottom 2 pt; top/right hidden |
| Major ticks | 2 pt wide, 6 pt long, outward |
| Minor ticks and grid | Disabled |
| Background | White |
| Export formats | PNG, PDF and SVG |
| Raster export | 600 dpi |
| Vector font settings | PDF fonttype 42; SVG fonttype none |

Approved figure-specific colors, data scales, legend positioning and `subplots_adjust` are deliberately controlled by each plot rather than forced globally.

## Python API

- `apply_style()`: enable approved global Matplotlib styling
- `style_context()`: temporary styling context
- `new_figure(size_mm=None)`: create a 130 × 100 mm figure and axes
- `style_axes(ax, xlabel=None, ylabel=None)`: standardize labels, ticks, spines
- `legend_bold(ax, **kwargs)`: frame-free bold legend
- `save_figure(fig, stem, output_dir=".", formats=None, dpi=None, close=False)`: export all formats

## Release and PyPI

GitHub Actions tests Python 3.9, 3.11 and 3.13 on push/PR. A separate release workflow builds and publishes using PyPI Trusted Publishing. To activate public `pip install figentra`, configure a trusted publisher on PyPI using GitHub owner `peterDataScientia`, repository `figentra`, workflow `publish.yml`, environment `pypi`; create the GitHub `pypi` environment and publish a GitHub Release after CI passes. Publishing requires your own PyPI authorization; simply pushing this repository does not publish the package.

## License

MIT — see [LICENSE](LICENSE).
