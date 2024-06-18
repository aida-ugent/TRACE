import matplotlib
import numpy as np


def toHex(colormap):
    return [matplotlib.colors.rgb2hex(colormap[i]) for i in range(len(colormap))]


PurpleGreen = matplotlib.colors.LinearSegmentedColormap.from_list(
    "diverging_purple_green",
    [
        "#762a83",
        "#884c91",
        "#996b9f",
        "#aa8aad",
        "#baa9bb",
        "#c9c9c9",
        "#88dd95",
        "#5ab26c",
        "#318947",
        "#176022",
        "#003a00",
    ],
    N=21,
)
matplotlib.colormaps.register(PurpleGreen)

palettes = {
    "viridis": toHex(matplotlib.colormaps["viridis"](np.linspace(0, 1, 21))),
    "continuous_BrGn": toHex(matplotlib.colormaps["BrBG"](np.linspace(0, 1, 21))),
    "continuous_PRGn": toHex(
        matplotlib.colormaps["diverging_purple_green"](np.linspace(0, 1, 21))
    ),
    "continuous_PRGn_old": [
        "#762a83",
        "#884c91",
        "#996b9f",
        "#aa8aad",
        "#baa9bb",
        "#c9c9c9",
        "#88dd95",
        "#5ab26c",
        "#318947",
        "#176022",
        "#003a00",
    ],
}
