import altair as alt
import pandas as pd
import streamlit as st


def truncate_text(text):
    if isinstance(text, str) and len(text) > 10:
        return text[:10] + "..."
    return text


def column_selectable_scatter(df: pd.DataFrame) -> alt.Chart:
    df = df.map(truncate_text)

    x_dropdown = alt.param(
        name="x_axis",
        bind=alt.binding_select(
            options=df.columns, name="X-axis:", element="#column-selectable-scatter-X"
        ),
        value=df.columns[0],
    )
    y_dropdown = alt.param(
        name="y_axis",
        bind=alt.binding_select(
            options=df.columns, name="Y-axis:", element="#column-selectable-scatter-Y"
        ),
        value=df.columns[1],
    )

    scatter = (
        alt.Chart(df)
        .mark_point()
        .encode(
            x=alt.X(
                "x:Q",
                title="X-axis",
                axis=alt.Axis(
                    titleY=-10,
                    labelPadding=-25,
                    offset=35,
                ),
            ),
            y=alt.Y(
                "y:Q",
                title="Y-axis",
                axis=alt.Axis(
                    titleX=-10,
                    labelPadding=-25,
                    offset=35,
                ),
            ),
        )
        .transform_calculate(x="datum[x_axis]", y="datum[y_axis]")
        .add_params(x_dropdown, y_dropdown)
        .properties(
            height=400,
            width=500,
            padding={"left": 30, "top": 30, "right": 30, "bottom": 20},
        )
        .interactive()
    )
    st.markdown(
        """
        <style>
            #column-selectable-scatter-X, #column-selectable-scatter-Y {
                margin-bottom: 10px;
                margin-right: 10px;
                text-align: left;
            }
            select {
                font-size: 14px;
                padding: 5px;
                margin-left:10px;
                margin-top: 5px;
                margin-bottom: 5px;
                border: 1px solid #ccc;
                border-radius: 5px;
                text-align: left;
                min-width: 200px; 
            }
        </style>
        """,
        unsafe_allow_html=True,
    )
    return scatter


def create_parity_plot(
    y_train_true: pd.Series,
    y_train_pred: pd.Series,
    y_test_true: pd.Series,
    y_test_pred: pd.Series,
) -> alt.Chart:
    df_train = pd.DataFrame(
        {"actual": y_train_true, "predicted": y_train_pred, "type": "train"}
    )
    df_test = pd.DataFrame(
        {"actual": y_test_true, "predicted": y_test_pred, "type": "test"}
    )
    df = pd.concat([df_train, df_test], ignore_index=True)

    min_val = float(min(df["actual"].min(), df["predicted"].min()))
    max_val = float(max(df["actual"].max(), df["predicted"].max()))

    df_line = pd.DataFrame(
        {"actual": [min_val, max_val], "predicted": [min_val, max_val]}
    )

    points = (
        alt.Chart(df)
        .mark_point(opacity=0.5, size=60)
        .encode(
            x=alt.X(
                "predicted",
                title="Predicted",
                scale=alt.Scale(domain=[min_val, max_val]),
            ),
            y=alt.Y(
                "actual", title="Actual", scale=alt.Scale(domain=[min_val, max_val])
            ),
            color=alt.Color("type:N", title="Data"),
            tooltip=["type", "actual", "predicted"],
        )
    )

    line = (
        alt.Chart(df_line)
        .mark_line(color="gray", strokeDash=[5, 5])
        .encode(x="predicted", y="actual")
    )

    chart = (points + line).properties(width=500, height=400).interactive()

    return chart
