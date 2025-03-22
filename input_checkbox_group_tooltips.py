from shiny.express import ui, input, render


icon_title = "About tooltips"


def fa_info_circle(title: str):
    # Enhanced from https://rstudio.github.io/fontawesome/ via `fontawesome::fa(&quot;info-circle&quot;, a11y = &quot;sem&quot;, title = icon_title)`
    return ui.HTML(
        f'<svg aria-hidden="true" role="img" viewBox="0 0 512 512" style="height:1em;width:1em;vertical-align:-0.125em;margin-left:auto;margin-right:auto;font-size:inherit;fill:currentColor;overflow:visible;position:relative;"><title>{title}</title><path d="M256 512A256 256 0 1 0 256 0a256 256 0 1 0 0 512zM216 336h24V272H216c-13.3 0-24-10.7-24-24s10.7-24 24-24h48c13.3 0 24 10.7 24 24v88h8c13.3 0 24 10.7 24 24s-10.7 24-24 24H216c-13.3 0-24-10.7-24-24s10.7-24 24-24zm40-208a32 32 0 1 1 0 64 32 32 0 1 1 0-64z"/></svg>'
    )


def input_checkbox_group_tooltips(
    id: str, label: str, choices: dict, tooltips: dict | None = None
):
    # Create a checkbox group with tooltips

    # Create the input object

    selections = []

    # If tooltips are not provided, use the choice names as tooltips

    if tooltips is None:
        tooltips = choices

    for choice in choices:
        with ui.tooltip(id="btn_tooltip", placement="right"):
            ui.input_action_button("btn", "A button with a tooltip")
            "The tooltip message"

        ui.input_checkbox(id, label, choices[choice])
    # Return the input object
    return selections


@render.text
def btn_tooltip_state():
    # Return the state of the tooltip

    return f"Tooltip state: {input.btn_tooltip()}"
