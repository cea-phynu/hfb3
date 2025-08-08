#!/usr/bin/env python3

# Plot the local 1-body density of a state.

# ==============================================================================
# ==============================================================================
# ==============================================================================

import sys
import hfb3

# ==============================================================================
# ==============================================================================
# ==============================================================================


def plot_matplotlib(_state, zmin, zmax, xmin, xmax):
    """plot using Matplotlib"""

    discrete = hfb3.Discrete(_state.basis, hfb3.Mesh.regular(xmin, 0, zmin, xmax, 0, zmax, 101, 1, 201))
    denst = discrete.getLocalXZ(_state.rho(hfb3.NEUTRON) + _state.rho(hfb3.PROTON), True)

    import numpy as np
    from matplotlib.pyplot import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, ax = plt.subplots()
    im = ax.imshow(np.flip(denst, axis=0), cmap=cm.inferno, extent=[zmin, zmax, xmin, xmax])
    plt.xlabel(r'$z$ [fm]')
    plt.ylabel(r'$r$ [fm]')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax)

    plt.ylabel(r'$\rho_{tot}$ [fm$^{-3}$]')
    plt.show()

# ==============================================================================
# ==============================================================================
# ==============================================================================


def plot_contour_matplotlib(_state, zmin, zmax, xmin, xmax, xstp, outputFile=None):
    """plot using Matplotlib"""
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib

    discrete = hfb3.Discrete(_state.basis, hfb3.Mesh.regular(xmin, 0, zmin, xmax, 0, zmax, 101, 1, 201))
    denst = np.flip(discrete.getLocalXZ(_state.rho(hfb3.NEUTRON) + _state.rho(hfb3.PROTON), True), axis=1)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 4))

    axes[0].plot(discrete.mesh.az.p, denst[0, :], color="black")
    axes[0].grid(True, linestyle='--')
    axes[0].set_xlim([zmin, zmax])
    axes[0].set_ylim(bottom=0)
    axes[0].set_xlabel(r'$z$ [fm]')
    axes[0].set_ylabel(r'$\rho_{tot}$ [fm$^{-3}$]')

    axes[1].contour(denst, np.arange(0, 0.2, 0.01), colors=["black",], extent=[zmin, zmax, xmin, xmax])
    plt.setp(axes[0], xticklabels=[])
    axes[1].grid(True, linestyle='--')
    axes[1].set_xlabel(r'$z$ [fm]')
    axes[1].set_ylabel(r'$r$ [fm]')

    # invert y axis
    ax = axes[1]
    ax.set_ylim(ax.get_ylim()[::-1])

    axes[1].yaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(xstp, xmax + 0.5, xstp)))

    plt.subplots_adjust(hspace=0, left=0.12, right=0.97, top=0.97, bottom=0.12)

    if outputFile:
        plt.savefig(outputFile)
        print(f"File {outputFile} generated.")
    else:
        plt.show()

    plt.close(fig)

# ==============================================================================
# ==============================================================================
# ==============================================================================


def plot_bokeh(_state, zmin, zmax, xmin, xmax):
    """plot using Bokeh"""

    from bokeh.layouts import gridplot
    from bokeh.plotting import figure, show
    from bokeh.models import ColorBar, LinearColorMapper

    # from bokeh.io import output_notebook
    # output_notebook()

    discrete = hfb3.Discrete(_state.basis, hfb3.Mesh.regular(xmin, 0, zmin, xmax, 0, zmax, 101, 1, 201))
    denst = discrete.getLocalXZ(_state.rho(hfb3.NEUTRON) + _state.rho(hfb3.PROTON), True)

    p = figure(tools="pan, reset, save, wheel_zoom",  # height=300, width=800,
               active_drag="pan", active_scroll="wheel_zoom", match_aspect=True,
               x_axis_label='z [fm]', y_axis_label='r [fm]')
    color = LinearColorMapper(palette="Inferno256")
    p.image(image=[denst], x=zmin, y=xmin, dw=zmax - zmin, dh=xmax - xmin, color_mapper=color)
    color_bar = ColorBar(color_mapper=color, label_standoff=10, location=(0, 0), width=10)
    p.add_layout(color_bar, 'right')

    # show(p)  # inside a notebook
    show(gridplot([p,], ncols=1, sizing_mode='stretch_both'))

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":

    if len(sys.argv) == 1:
        print("Trying to reproduce the figures from the manuscript.")
        try:
            plot_contour_matplotlib(hfb3.State("240Pu_deformed.msg.gz"), -15, 15, 0, 10, 2.0, "240Pu_deformed_density.eps")
        except Exception:
            print("Please generate '240Pu_deformed.msg.gz' first with 'bin/hfb3 examples/240Pu_deformed.hfb3'.")
        try:
            plot_contour_matplotlib(hfb3.State("16O_deformed.msg.gz"), -6, 6, 0, 4, 1.0, "16O_deformed_density.eps")
        except Exception:
            print("Please generate '16O_deformed.msg.gz' first with 'bin/hfb3 examples/16O_deformed.hfb3'.")

    else:
        if len(sys.argv) > 2:
            print(f"Usage: {sys.argv[0]} filename.msg[.gz]")
            sys.exit(0)

        fileName = sys.argv[1]
        state = hfb3.State(fileName)
        plot_contour_matplotlib(state, -15, 15, 0, 10, 1.0)

        # possible other visualization
        # plot_matplotlib(state, -20, 20, -10, 10)

        # possible other visualization
        # plot_bokeh(state, -20, 20, -10, 10)
