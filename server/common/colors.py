import re

HEX_COLOR_FORMAT = re.compile("^#[a-fA-F0-9]{6,6}$")


def convert_color_to_hex_format(unknown):
    if type(unknown) in (list, tuple):
        if all(0.0 <= ele <= 1.0 for ele in unknown):
            tup = tuple(int(ele * 255) for ele in unknown)
        elif all(0 <= ele <= 255 and isinstance(ele, int) for ele in unknown):
            tup = tuple(unknown)
        else:
            raise Exception("Unknown color format!")
        return "#%02x%02x%02x" % tup
    elif isinstance(unknown, str) and HEX_COLOR_FORMAT.match(unknown):
        return unknown.lower()
    else:
        raise Exception("Unknown color format!")


def anndata_colors_to_cxg_colors(data):
    return dict(
        (k[:-7], dict(zip(data.obs[k[:-7]].cat.categories, [convert_color_to_hex_format(c) for c in data.uns[k]])))
        for k in data.uns.keys()
        if k.endswith("_colors") and k[:-7] in data.obs.keys()
    )
