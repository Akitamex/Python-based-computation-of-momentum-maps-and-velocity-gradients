from typing import TypeVar

# Function Bar Alias
Vec = TypeVar('Vec', list, tuple)  # Abstract 1D Array
Mat = TypeVar('Mat', list, tuple)  # Abstract 2D Array
Cube = TypeVar('Cube', list, tuple)  # Abstract 3D Array
Vecx = TypeVar('Vecx', list, tuple)  # Linear Range

