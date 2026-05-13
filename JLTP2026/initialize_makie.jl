using CairoMakie, ColorSchemes, LaTeXStrings

temperature_colors = Makie.wong_colors()
mass_colors = ColorSchemes.tab10
asymptotic_colors = getindex.(Ref(ColorSchemes.tab10), 4:6) # ColorSchemes.tab20[2:2:end]
light_mass_colors = asymptotic_colors
alpha = 1.0

pt = 4/3 # px
inch = 96 # px
mm = inch/25.4


my_theme = Theme(
    fontsize=8pt,
    labelsize=9pt,
    figure_padding = 0,
    Axis = (
        xticksize = 3,
        yticksize = 3,
        xlabelsize = 9pt,
        ylabelsize = 9pt,
    ),
    Label = (
        fontsize = 9pt,
        ),
        )
set_theme!(merge(my_theme, theme_latexfonts()))
