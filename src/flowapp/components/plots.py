import pandas as pd
import plotly.graph_objects as go


def column_selectable_scatter(df: pd.DataFrame) -> go.Figure:
    default_x = df.columns[0]
    default_y = df.columns[1]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df[default_x], y=df[default_y], mode="markers", marker=dict(size=8)
        )
    )

    x_buttons = []
    for col in df.columns:
        x_buttons.append(dict(args=[{"x": [df[col]]}], label=col, method="restyle"))

    y_buttons = []
    for col in df.columns:
        y_buttons.append(dict(args=[{"y": [df[col]]}], label=col, method="restyle"))

    updatemenus = [
        dict(
            buttons=x_buttons,
            direction="down",
            showactive=True,
            y=1.15,
            pad={"r": 10, "t": 10},
            name="X-axis",
            xanchor="left",
            yanchor="top",
        ),
        dict(
            buttons=y_buttons,
            direction="down",
            showactive=True,
            y=1.15,
            pad={"r": 10, "t": 10},
            name="Y-axis",
            xanchor="right",
            yanchor="top",
        ),
    ]

    fig.update_layout(
        updatemenus=updatemenus,
        xaxis_title=default_x,
        yaxis_title=default_y,
        height=600,
    )
    return fig
