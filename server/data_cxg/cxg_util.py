import numpy as np


def pack_selector(boolarray):
    """
    pack all contiguous selectors into slices.  Remember that
    tiledb multi_index requires INCLUSIVE indices.
    """

    if boolarray is None:
        return slice(None)

    assert type(boolarray) == np.ndarray

    selector = np.nonzero(boolarray)[0]

    if len(selector) == 0:
        return slice(None)

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
