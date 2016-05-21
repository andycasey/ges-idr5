
""" General utility functions. """


def safe_int(x, fill_value=-1):
    try:
        return int(x)
    except:
        return fill_value


def wg_as_int(wg):
    return int(str(wg).strip().lower().lstrip("wg"))