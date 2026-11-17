import numpy as np

from typing import override

from fractions import Fraction

from matplotlib.ticker import Formatter

from dexter.types import Array1
from dexter._core import _py_get_max_threads, _py_set_num_threads

# ================================================================================================


def get_max_threads() -> int:
    """Returns the number of the device's available threads."""
    return _py_get_max_threads()


def set_num_threads(num: int):
    """Sets the global number of threads.

    Parameters
    ----------
    num
        The number of threads to use. If `num` is either $0$ or greater than the device’s number
        of threads, it defaults to the device’s number of threads.

    Example
    ------
    ```python title="Configure number of threads"
    >>> dex.utils.set_num_threads(4)

    ```
    """
    _py_set_num_threads(num)


# ================================================================================================


def prepare_contour_levels(levels: Array1 | list[float]) -> Array1:
    r"""Prepares a 1D array to be passed to `plt.contour`.

    Matplotlib's contour method require the `levels` argument to be a 1D array of strictly
    increasing values. This method ensures the returned array meets this condition.

    Parameters
    ----------
    levels
        The values to be sorted.
    """
    return np.sort(list(set(levels)))


class FractionFormatter(Formatter):
    r"""Formats values as integer ratios.

    An instance of this class can be passed as the `fmt` argument to `plt.clabel` to format contour
    labels to the closest rational fraction.

    Parameters
    ----------
    denominator_limit
        The maximum value the fraction's denominator is allowed to have.
    """

    denominator_limit: int

    def __init__(self, denominator_limit: int = 50) -> None:
        self.denominator_limit = denominator_limit
        super().__init__()

    @override
    def format_data(self, value: float) -> str:
        return str(Fraction(value).limit_denominator(self.denominator_limit))

    @override
    def format_ticks(self, values: list[float]) -> list[str]:
        return [self.format_data(value) for value in values]

    @override
    def __call__(self, x: float, _: int | None = None) -> str:
        return self.format_data(x)
