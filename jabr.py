from __future__ import division
from mosek.fusion import Matrix, DenseMatrix, Model, Expr, Domain, ObjectiveSense, SolutionError
from numpy import array, sqrt, real, imag
from math import sin, cos, asin, pi
from scipy.sparse import dok_matrix, hstack

S2 = sqrt(2)


def z2y(r, x):
    """ converts impedance Z=R+jX to admittance Y=G+jB """
    return r/(r**2+x**2), -x/(r**2+x**2)


def jabr_baby():
    r, x = .02, .2  # impedance and resistance, from matpower
    g, b = z2y(r, x)  # convert to susceptance and conductance
    A = DenseMatrix([[0., g*S2, -g, -b], [0, -b*S2, b, -g]])
    rhs = [-1., 0]
    with Model('jabr') as M:
        u = M.variable("u", 2, Domain.greaterThan(0.0))
        R = M.variable("R", 1, Domain.greaterThan(0.0))
        I = M.variable("I", 1, Domain.unbounded())
        x = Expr.vstack(u, R, I)
        M.objective("maxRsum", ObjectiveSense.Maximize, R)
        M.constraint("u0",  u.index(0).asExpr(), Domain.equalsTo(1/S2))
        M.constraint("flow",  Expr.mul(A, x), Domain.equalsTo(rhs))
        M.constraint("cone", x, Domain.inRotatedQCone())
        M.solve()
        print(u.level())
        print(R.level())
        print(I.level())


def jabr5():
    r, x = .01, .1
    g, b = z2y(r, x)
    A1 = DenseMatrix([[0., g*S2, -g, -b], [0, -b*S2, b, -g]])
    A2 = DenseMatrix([[0., 3*g*S2, -g, -b, 0, 0, -g, b, 0, 0, -g, b],
                      [0., -3*b*S2, b, -g, 0, 0, b, g, 0, 0, b, g]])
    A3 = DenseMatrix([[0., g*S2, -g, -b], [0, -b*S2, b, -g]])
    A4 = DenseMatrix([[0., g*S2, -g, -b], [0, -b*S2, b, -g]])
    rhs = [-1., 0]

    def jabr():
        with Model('jabr') as M:
            u = M.variable("u", 5, Domain.greaterThan(0.0))
            R = M.variable("R", 4, Domain.greaterThan(0.0))
            I = M.variable("I", 4, Domain.unbounded())

            x01 = Expr.vstack(Expr.vstack(u.index(0), u.index(1)), R.index(0), I.index(0))
            x02 = Expr.vstack(Expr.vstack(u.index(0), u.index(2)), R.index(1), I.index(1))
            x23 = Expr.vstack(Expr.vstack(u.index(2), u.index(3)), R.index(2), I.index(2))
            x24 = Expr.vstack(Expr.vstack(u.index(2), u.index(4)), R.index(3), I.index(3))

            M.constraint("u0",  u.index(0).asExpr(), Domain.equalsTo(1/S2))

            M.constraint("cone01", x01, Domain.inRotatedQCone())
            M.constraint("cone02", x02, Domain.inRotatedQCone())
            M.constraint("cone23", x23, Domain.inRotatedQCone())
            M.constraint("cone24", x24, Domain.inRotatedQCone())

            M.constraint("flow1", Expr.mul(A1, x01), Domain.equalsTo(rhs))
            M.constraint("flow2", Expr.mul(A2, Expr.vstack(x02, x23, x24)), Domain.equalsTo(rhs))
            M.constraint("flow3", Expr.mul(A3, x23), Domain.equalsTo(rhs))
            M.constraint("flow4", Expr.mul(A4, x24), Domain.equalsTo(rhs))

            M.objective("maxRsum", ObjectiveSense.Maximize, Expr.sum(R))
            M.solve()
            return {'u': u.level(), 'R': R.level(), 'I': I.level()}
    return jabr()


def process_answer(ans):
    r, x = .01, .1
    g, b = z2y(r, x)
    U = array(ans['u'])
    V = sqrt(sqrt(2)*U)
    R, I = ans['R'], ans['I']
    t01 = asin(I[0]/V[0]/V[1])
    t02 = asin(I[1]/V[0]/V[2])
    t23 = asin(I[2]/V[2]/V[3])
    t24 = asin(I[3]/V[2]/V[4])
    print(V)
    print([round(180/pi*x, 3) for x in [t01, t02, t23, t24]])
    print(2*sqrt(2)*g*U[0]-g*(R[0]+R[1])+b*(I[0]+I[1]))
    print(g*(V[0]**2-V[0]*V[1]*cos(t01)) - b*V[0]*V[1]*sin(-t01),
          g*(V[0]**2-V[0]*V[2]*cos(t02)) - b*V[0]*V[2]*sin(-t02))


def build_U_matrices(G, B):
    n = G.shape[0]
    Ureal = dok_matrix((n-1, n))
    Ureac = dok_matrix((n-1, n))
    S2 = 2**.5
    for i in range(1, n):
        Ureal[i-1, i] = S2*G[i, :].sum()
        Ureac[i-1, i] = -S2*B[i, :].sum()
    return Ureal, Ureac


def build_R_matrices(G, B, branch_map):
    """ rows are buses 2 to n; cols are branches  """
    n = G.shape[0]
    Rreal = dok_matrix((n-1, n-1))
    Rreac = dok_matrix((n-1, n-1))
    for fbus in range(1, n):
        for tbus in B[fbus, :].nonzero()[1]:
            branch = branch_map[(fbus, tbus)]
            Rreal[fbus-1, branch] = -G[fbus, tbus]
            Rreac[fbus-1, branch] = B[fbus, tbus]
    return Rreal, Rreac


def build_I_matrices(G, B, branch_map):
    """ rows are buses 2 to n; cols are branches  """
    n = B.shape[0]
    Ireal = dok_matrix((n-1, n-1))
    Ireac = dok_matrix((n-1, n-1))
    for fbus in range(1, n):
        for tbus in B[fbus, :].nonzero()[1]:
            branch = branch_map[(fbus, tbus)]
            s = -1 if fbus < tbus else 1
            Ireal[fbus-1, branch] = s*B[fbus, tbus]
            Ireac[fbus-1, branch] = s*G[fbus, tbus]
    return Ireal, Ireac


def build_constraint_matrix(G, B, branch_map):
    Ureal, Ureac = build_U_matrices(G, B)
    Rreal, Rreac = build_R_matrices(G, B, branch_map)
    Ireal, Ireac = build_I_matrices(G, B, branch_map)
    Areal = hstack([Ureal, Rreal, Ireal])
    Areac = hstack([Ureac, Rreac, Ireac])
    return Areal, Areac


def convert_matrix_to_mosek(M):
    """ converts coo_matrix to mosek SparseMatrix type  """
    m, n = M.shape
    rows = M.row.astype(int).tolist()
    cols = M.col.astype(int).tolist()
    vals = M.data.tolist()
    return Matrix.sparse(m, n, rows, cols, vals)


def build_mosek_model(case):
    """ the mosek Model type seems to be a bit brittle and doesn't like being
    passed around, inspecting its attributes, etc. Simplest to do this all in
    one shot rather than splitting it into smaller pieces
    """
    Areal, Areac = build_constraint_matrix(case.G, case.B, case.branch_map)
    Areal, Areac = map(convert_matrix_to_mosek, (Areal, Areac))
    P = -1 * real(case.demands)[1:]
    Q = -1 * imag(case.demands)[1:]
    branches = case.branch_list
    n = len(case.demands)
    v = case.vhat
    with Model('jabr%d' % n) as M:
        u = M.variable("u", n, Domain.greaterThan(0.0))
        R = M.variable("R", n-1, Domain.greaterThan(0.0))
        I = M.variable("I", n-1, Domain.unbounded())

        expr = lambda x: x.asExpr()

        def cone_constraint(fbus, tbus, branch):
            x = Expr.vstack(map(expr, [u.index(fbus), u.index(tbus), R.index(branch), I.index(branch)]))
            M.constraint("cone%d_%d" % (fbus, tbus), x, Domain.inRotatedQCone())

        M.constraint("u0",  u.index(0).asExpr(), Domain.equalsTo(v*v/sqrt(2)))
        for branch, (fbus, tbus) in enumerate(branches):
            cone_constraint(fbus, tbus, branch)
        M.constraint("flow_real", Expr.mul(Areal, Expr.vstack(u, R, I)), Domain.equalsTo(P))
        M.constraint("flow_reac", Expr.mul(Areac, Expr.vstack(u, R, I)), Domain.equalsTo(Q))

        M.objective("maxRsum", ObjectiveSense.Maximize, Expr.sum(R))
        out_file = open('/Users/srharnett/Dropbox/power/jabr-power-flow/mosek.log', 'a')
        M.setLogHandler(out_file)
        #M.setSolverParam("intpntCoTolPfeas", 1.0e-2)
        #M.setSolverParam("intpntCoTolDfeas", 1.0e-2)
        #M.setSolverParam("intpntCoTolRelGap", 1.0e-2)
        M.solve()
        try:
            return u.level(), R.level(), I.level()
        except SolutionError:
            status = M.getPrimalSolutionStatus()
            raise SolutionError("mosek failed to converge: %s (check log)" % status)


def recover_original_variables(u, R, I):
    """ incomplete. right now just computes the original bus voltages  """
    V = sqrt(sqrt(2) * array(u))
    return V


if __name__ == '__main__':
    process_answer(jabr5())