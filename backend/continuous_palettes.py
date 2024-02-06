import matplotlib.pyplot as plt
import matplotlib

def toHex(colormap):
    return [matplotlib.colors.rgb2hex(colormap(i)) for i in range(colormap.N)]

purple_green = toHex(matplotlib.colors.LinearSegmentedColormap.from_list(
    "continuous_purple_green",
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
))

# toHex(plt.cm.get_cmap("PRGn", 20)), this colormap has white in the center

palettes = {
    "viridis": toHex(plt.cm.get_cmap("viridis", 21)),
    "continuous_BrGn": toHex(plt.cm.get_cmap("BrBG", 21)),
    "continuous_PRGn": purple_green,
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
