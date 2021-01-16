import numpy as np


def pack_selector_from_mask(boolarray):
    """
    pack all contiguous selectors into slices.  Remember that
    tiledb multi_index requires INCLUSIVE indices.
    """

    if boolarray is None:
        return slice(None)

    assert type(boolarray) == np.ndarray
    assert boolarray.dtype == bool

    selector = np.nonzero(boolarray)[0]
    return pack_selector_from_indices(selector)


def pack_selector_from_indices(selector):

    if len(selector) == 0:
        return None

    result = []
    current = slice(selector[0], selector[0])
    for sel in selector[1:]:
        if sel == current.stop + 1:
            current = slice(current.start, sel)
        else:
            result.append(current if current.start != current.stop else current.start)
            current = slice(sel, sel)

    if len(result) == 0 or result[-1] != current:
        result.append(current if current.start != current.stop else current.start)

    return result
