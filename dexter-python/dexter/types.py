import numpy as np
from typing import Literal, TypeAlias, Union


Interp1DType: TypeAlias = Literal[
    "linear", "cubic", "cubic periodic", "akima", "akima periodic", "steffen"
]
Interp2DType: TypeAlias = Literal["bilinear", "bicubic"]

NDArrayShape: TypeAlias = tuple[int, ...]
NDArray1D: TypeAlias = np.ndarray[tuple[int], np.dtype[np.float64]]
NDArray2D: TypeAlias = np.ndarray[tuple[int, int], np.dtype[np.float64]]

PoincareSection: TypeAlias = Literal["ConstTheta", "ConstZeta"]

CalculatedFrequency: TypeAlias = Union[float, None]
