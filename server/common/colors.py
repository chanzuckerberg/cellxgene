import re

from server.common.errors import ColorFormatException

HEX_COLOR_FORMAT = re.compile("^#[a-fA-F0-9]{6,6}$")

# https://www.w3.org/TR/css-color-4/#named-colors
CSS4_NAMED_COLORS = dict(
    aliceblue="#f0f8ff",
    antiquewhite="#faebd7",
    aqua="#00ffff",
    aquamarine="#7fffd4",
    azure="#f0ffff",
    beige="#f5f5dc",
    bisque="#ffe4c4",
    black="#000000",
    blanchedalmond="#ffebcd",
    blue="#0000ff",
    blueviolet="#8a2be2",
    brown="#a52a2a",
    burlywood="#deb887",
    cadetblue="#5f9ea0",
    chartreuse="#7fff00",
    chocolate="#d2691e",
    coral="#ff7f50",
    cornflowerblue="#6495ed",
    cornsilk="#fff8dc",
    crimson="#dc143c",
    cyan="#00ffff",
    darkblue="#00008b",
    darkcyan="#008b8b",
    darkgoldenrod="#b8860b",
    darkgray="#a9a9a9",
    darkgreen="#006400",
    darkgrey="#a9a9a9",
    darkkhaki="#bdb76b",
    darkmagenta="#8b008b",
    darkolivegreen="#556b2f",
    darkorange="#ff8c00",
    darkorchid="#9932cc",
    darkred="#8b0000",
    darksalmon="#e9967a",
    darkseagreen="#8fbc8f",
    darkslateblue="#483d8b",
    darkslategray="#2f4f4f",
    darkslategrey="#2f4f4f",
    darkturquoise="#00ced1",
    darkviolet="#9400d3",
    deeppink="#ff1493",
    deepskyblue="#00bfff",
    dimgray="#696969",
    dimgrey="#696969",
    dodgerblue="#1e90ff",
    firebrick="#b22222",
    floralwhite="#fffaf0",
    forestgreen="#228b22",
    fuchsia="#ff00ff",
    gainsboro="#dcdcdc",
    ghostwhite="#f8f8ff",
    gold="#ffd700",
    goldenrod="#daa520",
    gray="#808080",
    green="#008000",
    greenyellow="#adff2f",
    grey="#808080",
    honeydew="#f0fff0",
    hotpink="#ff69b4",
    indianred="#cd5c5c",
    indigo="#4b0082",
    ivory="#fffff0",
    khaki="#f0e68c",
    lavender="#e6e6fa",
    lavenderblush="#fff0f5",
    lawngreen="#7cfc00",
    lemonchiffon="#fffacd",
    lightblue="#add8e6",
    lightcoral="#f08080",
    lightcyan="#e0ffff",
    lightgoldenrodyellow="#fafad2",
    lightgray="#d3d3d3",
    lightgreen="#90ee90",
    lightgrey="#d3d3d3",
    lightpink="#ffb6c1",
    lightsalmon="#ffa07a",
    lightseagreen="#20b2aa",
    lightskyblue="#87cefa",
    lightslategray="#778899",
    lightslategrey="#778899",
    lightsteelblue="#b0c4de",
    lightyellow="#ffffe0",
    lime="#00ff00",
    limegreen="#32cd32",
    linen="#faf0e6",
    magenta="#ff00ff",
    maroon="#800000",
    mediumaquamarine="#66cdaa",
    mediumblue="#0000cd",
    mediumorchid="#ba55d3",
    mediumpurple="#9370db",
    mediumseagreen="#3cb371",
    mediumslateblue="#7b68ee",
    mediumspringgreen="#00fa9a",
    mediumturquoise="#48d1cc",
    mediumvioletred="#c71585",
    midnightblue="#191970",
    mintcream="#f5fffa",
    mistyrose="#ffe4e1",
    moccasin="#ffe4b5",
    navajowhite="#ffdead",
    navy="#000080",
    oldlace="#fdf5e6",
    olive="#808000",
    olivedrab="#6b8e23",
    orange="#ffa500",
    orangered="#ff4500",
    orchid="#da70d6",
    palegoldenrod="#eee8aa",
    palegreen="#98fb98",
    paleturquoise="#afeeee",
    palevioletred="#db7093",
    papayawhip="#ffefd5",
    peachpuff="#ffdab9",
    peru="#cd853f",
    pink="#ffc0cb",
    plum="#dda0dd",
    powderblue="#b0e0e6",
    purple="#800080",
    rebeccapurple="#663399",
    red="#ff0000",
    rosybrown="#bc8f8f",
    royalblue="#4169e1",
    saddlebrown="#8b4513",
    salmon="#fa8072",
    sandybrown="#f4a460",
    seagreen="#2e8b57",
    seashell="#fff5ee",
    sienna="#a0522d",
    silver="#c0c0c0",
    skyblue="#87ceeb",
    slateblue="#6a5acd",
    slategray="#708090",
    slategrey="#708090",
    snow="#fffafa",
    springgreen="#00ff7f",
    steelblue="#4682b4",
    tan="#d2b48c",
    teal="#008080",
    thistle="#d8bfd8",
    tomato="#ff6347",
    turquoise="#40e0d0",
    violet="#ee82ee",
    wheat="#f5deb3",
    white="#ffffff",
    whitesmoke="#f5f5f5",
    yellow="#ffff00",
    yellowgreen="#9acd32",
)


def convert_color_to_hex_format(unknown):
    """
    Try to convert color info to a hex triplet string https://en.wikipedia.org/wiki/Web_colors#Hex_triplet.

    The function accepts for the following formats:
    - A CSS4 color name, as supported by matplotlib https://matplotlib.org/3.1.0/gallery/color/named_colors.html
    - RGB tuple/list with values ranging from 0.0 to 1.0, as in [0.5, 0.75, 1.0]
    - RFB tuple/list with values ranging from 0 to 255, as in [128, 192, 255]
    - Hex triplet string, as in "#08c0ff"

    :param unknown: color info of unknown format
    :return: a hex triplet representing that color
    """
    try:
        if type(unknown) in (list, tuple) and len(unknown) == 3:
            if all(0.0 <= ele <= 1.0 for ele in unknown):
                tup = tuple(int(ele * 255) for ele in unknown)
            elif all(0 <= ele <= 255 and isinstance(ele, int) for ele in unknown):
                tup = tuple(unknown)
            else:
                raise ColorFormatException("Unknown color iterable format!")
            return "#%02x%02x%02x" % tup
        elif isinstance(unknown, str) and unknown.lower() in CSS4_NAMED_COLORS:
            return CSS4_NAMED_COLORS[unknown.lower()]
        elif isinstance(unknown, str) and HEX_COLOR_FORMAT.match(unknown):
            return unknown.lower()
        else:
            raise ColorFormatException("Unknown color format type!")
    except Exception as e:
        raise ColorFormatException(e)


def convert_anndata_category_colors_to_cxg_category_colors(data):
    """
    Convert color information from anndata files to the cellxgene color data format as described below:
        {
            "<category_name>": {
                "<label_name>": "<color_hex_code>",
                ...
            },
            ...
        }

    For more on the cxg color data structure, see https://github.com/chanzuckerberg/cellxgene/issues/1307.

    For more on the anndata color data structure, see
    https://github.com/chanzuckerberg/cellxgene/issues/1152#issuecomment-587276178.

    Handling of malformed data:
    - For any color info in a adata.uns[f"{category}_colors"] color array that convert_color_to_hex_format cannot
      convert to a hex triplet string, a ColorFormatException is raised
    - No category_name key group is returned for adata.uns[f"{category}_colors"] keys for which there is no
      adata.obs[f"{category}"] key

    :param data: the anndata file
    :return: cellxgene color data structure as described above
    """
    cxg_colors = dict()
    color_key_suffix = "_colors"
    for uns_key in data.uns.keys():
        # find uns array that describes colors for a category
        if not uns_key.endswith(color_key_suffix):
            continue

        # check to see if we actually have observations for that category
        category_name = uns_key[: -len(color_key_suffix)]
        if category_name not in data.obs.keys():
            continue

        # create the cellxgene color entry for this category
        cxg_colors[category_name] = dict(
            zip(data.obs[category_name].cat.categories, [convert_color_to_hex_format(c) for c in data.uns[uns_key]])
        )
    return cxg_colors
