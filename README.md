# polytrope
A simple polytrope class written in Python. The `Polytrope` object is initialised with the polytropic index `n` and calculates the relevant dimensionless parameters.

# getting started.
All that is required to use this class is `numpy`, `scipy` and `matplotlib`. These can be installed using `conda`:
```
conda install numpy scipy matplotlib
```

# basic usage.
Here's how to use this class to get the dimensionless radius, mass, moment of inertia and potential:
```
from polytrope import Polytrope

P = Polytrope(1.5)

print(P.xi_1, P.m, P.i, P.Omega)
```
