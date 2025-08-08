#!/usr/bin/env python3


import pandas as pd
from bokeh.plotting import figure, show
from bokeh.models import CustomJS, HoverTool, ColumnDataSource
from bokeh.models.widgets import Select
from bokeh.layouts import row, column
# from bokeh.layouts import gridplot
from bokeh.palettes import Category10
import sys

# ==============================================================================
# ==============================================================================
# ==============================================================================


def create_plot(dataFrames, titles):

    sources = [ColumnDataSource(df) for df in dataFrames]

    render_sources = [ColumnDataSource({"x": df[df.columns[0]],
                                        "y": df[df.columns[1]],
                                        "f": [titles[i]] * len(df),
                                        }) for i, df in enumerate(dataFrames)]

    tools = "crosshair,pan,wheel_zoom,box_zoom,reset,hover"

    colors = Category10[len(dataFrames)]

    f = figure(tools=tools,
               title='',
               x_axis_label=dataFrames[0].columns[0],
               y_axis_label=dataFrames[0].columns[1],
               sizing_mode="stretch_both",
               )

    f.select_one(HoverTool).tooltips = [('f', '@f'), ('x', '$x'), ('y', '$y')]

    for i, rs in enumerate(render_sources):
        f.scatter(source=rs, x="x", y="y", color=colors[i], line_width=1, legend_label=titles[i])

    f.legend.location = "top_right"
    f.legend.click_policy="hide"

    columnNames = list(dataFrames[0].columns)

    select0 = Select(title="X axis:", options=columnNames, value=columnNames[0])
    select1 = Select(title="Y axis:", options=columnNames, value=columnNames[1])

    jsArgs = dict(
        render_sources=render_sources,
        sources=sources,
        x_selector=select0,
        y_selector=select1,
        xaxis=f.xaxis,
        yaxis=f.yaxis,)

    jsCode = """
    // Extract what we want to color by from selector
    // let colorby = colorby_selector.value;

    // View of the colors for convenience
    // let colors = render_sources[0].data['color'];

    // Convenient to have the number of data points
    // let n = colors.length;

    // New data
    for (var i = 0; i < render_sources.length; i++)
    {
      render_sources[i].data['x'] = sources[i].data[x_selector.value];
      render_sources[i].data['y'] = sources[i].data[y_selector.value];
    }

    console.log(x_selector.value, y_selector.value)

    // Update axis labels to reflect what was selected
    xaxis[0].axis_label = x_selector.value;
    yaxis[0].axis_label = y_selector.value;

    for (var i = 0; i < render_sources.length; i++)
    {
      render_sources[i].change.emit();
    }
    """

    select0.js_on_change("value", CustomJS(code=jsCode, args=jsArgs))
    select1.js_on_change("value", CustomJS(code=jsCode, args=jsArgs))

    # layout = row(
    #     f,
    #     Spacer(width=15),
    #     column( select0, Spacer(height=15), select1, Spacer(height=15),))

    # show(row(column(select0, select1, height=100, width=100, sizing_mode="fixed"), f), sizing_mode='stretch_width')
    # show(gridplot([row([select0, select1], height=100, width=100, sizing_mode="fixed"), f], ncols=1, sizing_mode='stretch_both'))
    show(row(column(select0, select1, width=100), f, sizing_mode="stretch_both"))

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":
    fileNames = sys.argv[1:]
    datas = [pd.read_csv(f) for f in fileNames]
    create_plot(datas, fileNames)
