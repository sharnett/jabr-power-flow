A Python implementation of a conic programming solution to power flow on tree networks, as described in this paper:

```
R. A. Jabr, "Radial distribution load flow using conic programming,"
IEEE Transactions on Power Systems, vol. 21, no. 3, pp. 1458-1459, 2006.
```

It relies on MOSEK, a commercial (but free academic licenses) optimization solver. Unfortunately the nicer version
of their Python API only works with Python 2.x.