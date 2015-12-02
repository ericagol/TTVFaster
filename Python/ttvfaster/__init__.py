__all__ = ["run_ttvfaster"]

try:
    from ._ttvfaster import run_ttvfaster
except ImportError:
    raise ImportError(
        "It looks like you need to build the Cython extension or navigate "
        "to a different directory if you already built"
    )
