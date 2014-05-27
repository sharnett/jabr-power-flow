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

Once you have gurobi installed, try running

```
python jabr.py
```

which should produce

```
  1   0.853 -27.815
  2   0.985  -5.829
  3   0.873 -20.100
  4   0.853 -27.815
  5   1.000   0.000
```

Try it on some other examples in the cases directory:

```python
answer = solve('cases/case15_v2.m')
for bus in sorted(answer.keys()):
    v, t = answer[bus]
    print('%3d %7.3f %7.3f' % (bus, v, 180/pi*t))
```

```

  1   1.000   0.000
  2   0.971   0.032
  3   0.957   0.049
  4   0.951   0.056
  5   0.950   0.068
  6   0.958   0.189
  7   0.956   0.216
  8   0.957   0.205
  9   0.968   0.072
 10   0.967   0.085
 11   0.950   0.131
 12   0.946   0.182
 13   0.945   0.198
 14   0.949   0.085
 15   0.948   0.087
 ```
