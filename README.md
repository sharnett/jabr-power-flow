A Python implementation of a conic programming formulation to the power flow
problem on tree networks, as described in this paper:

```
R. A. Jabr, "Radial distribution load flow using conic programming,"
IEEE Transactions on Power Systems, vol. 21, no. 3, pp. 1458-1459, 2006.
```

It relies on MOSEK, a commercial (but free academic license) optimization
solver. Unfortunately the nicer ('fusion') version of their Python API only
works with Python 2.x.

The work I'm doing on the [power grid attack problem](https://github.com/sharnett/tree-power-flow)
is going to involve testing the power flow on grids that are far from their
normal operating conditions. I don't trust the usual Newton's method approach
(as in MATPOWER) to be sufficient for all cases; I want a certificate of
infeasibility, which this method provides.

This will also be my first project using Test Driven Development.

