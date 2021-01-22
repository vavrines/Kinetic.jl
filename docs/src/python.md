# Calling from Python 

For maximum convenience, a wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate all the methods from Python.

## How to use

Let's start by cloning the repository and changing into the directory.
```bash
git clone https://github.com/vavrines/kineticpy.git
cd kineticpy
```

Next, we start `python`.
The Julia main module can be installed and initialized by
```python
>>> import kineticpy
>>> kineticpy.install()
```

The basic structs and methods are stored in the base module, and can be imported via
```python
>>> from kineticpy import base
```

## Example

We provide some quick tutorial here for kineticpy.

```python
from kineticpy import base
import numpy as np

u = np.linspace(-5, 5, 28)
prim_var = np.array([1.0, 0.0, 1.0])
M = base.maxwellian(u, prim_var) # compute Maxwellian distribution
M.view()
```

```
array([7.83543327e-12, 2.77323769e-10, 7.46041809e-09, 1.52542631e-07,
       2.37067103e-06, 2.80029217e-05, 2.51412806e-04, 1.71562923e-03,
       8.89839075e-03, 3.50793472e-02, 1.05109877e-01, 2.39379825e-01,
       4.14365469e-01, 5.45169515e-01, 5.45169515e-01, 4.14365469e-01,
       2.39379825e-01, 1.05109877e-01, 3.50793472e-02, 8.89839075e-03,
       1.71562923e-03, 2.51412806e-04, 2.80029217e-05, 2.37067103e-06,
       1.52542631e-07, 7.46041809e-09, 2.77323769e-10, 7.83543327e-12])
```