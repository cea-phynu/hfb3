#!/usr/bin/env python3

from bokeh.plotting import figure, show
from bokeh.models import HoverTool
from bokeh.palettes import Category10_10
import csv
import sys

# ==============================================================================
# ==============================================================================
# ==============================================================================


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} file.csv")
        sys.exit(-1)

    files = sys.argv[1:]

    curves = []

    for f in files:
        with open(f, newline='') as csvfile:
            data = csv.reader(csvfile, delimiter=',')
            x = []
            y = []
            for d in data:
                x.append(float(d[0]))
                y.append(float(d[1]))
            curves.append([f, x, y])

    # Adding the plot
    p = figure(title='Main Title',
               x_axis_label='x',
               y_axis_label='y'
               )

    # Rendering the graph
    color = iter(Category10_10)
    for c in curves:
        line = p.line(c[1], c[2],
                      legend_label=c[0],
                      line_width=1.5,
                      color=next(color))
        # p.circle(c[1], c[2])

    p.tools.append(HoverTool(tooltips=[('y', '@y'), ('x', '@x')],
                             renderers=[line], mode='vline'))

    p.legend.click_policy = "hide"

    # Display the results
    # show(row(p, sizing_mode='stretch_both'))
    show(p)

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == '__main__':
    main()
