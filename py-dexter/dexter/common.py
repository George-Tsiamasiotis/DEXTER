"""Final wrappers of common exported free functions."""

from dexter._core import _py_get_max_threads, _py_set_num_threads


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
    >>> dex.set_num_threads(4)

    ```
    """
    _py_set_num_threads(num)
