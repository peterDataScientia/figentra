"""Example only: these are generated demonstration values, not laboratory data."""
import numpy as np
from figentra import new_figure, style_axes, legend_bold, save_figure

x = np.linspace(0, 120, 200)
y = np.exp(-x / 40)
fig, ax = new_figure()
ax.plot(x, y, linewidth=2, label="Demonstration")
style_axes(ax, xlabel="Time [s]", ylabel="Response [a.u.]")
legend_bold(ax)
fig.subplots_adjust(left=.22, right=.92, bottom=.15, top=.81)
save_figure(fig, "Figure_01_demo", output_dir="figures", close=True)
