A Python implementation of a conic programming formulation to the power flow
problem on tree networks, as described in the paper

```
R. A. Jabr, "Radial distribution load flow using conic programming,"
IEEE Transactions on Power Systems, vol. 21, no. 3, pp. 1458-1459, 2006.
```

It relies on Gurobi, a commercial (but free academic license) optimization
solver, which unfortunately requires python 2.

The work I'm doing on the [power grid attack problem](https://github.com/sharnett/tree-power-flow)
involves testing the power flow on grids that are far from their
normal operating conditions. The usual Newton's method approach can't be trusted
(as in MATPOWER) for all cases; I want a certificate of
infeasibility, which this method provides.
