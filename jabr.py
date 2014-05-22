from __future__ import division, print_function
import sh
from gurobipy import Model as GurobiModel, GRB, quicksum
from mosek.fusion import Matrix, DenseMatrix, Model, Expr, Domain, ObjectiveSense, SolutionError
from numpy import array, sqrt, real, imag
from math import sin, cos, asin, pi
from scipy.sparse import dok_matrix, hstack
from datetime import datetime
from collections import defaultdict, deque

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
        out_file.write('\n' + str(datetime.now()) + '\n')
        M.setLogHandler(out_file)
        M.writeTask('mosek.opf')
        M.solve()
        try:
            raise SolutionError
            #return u.level(), R.level(), I.level()
        except SolutionError:
            out_file.write("\nPython thinks it didn't converge, trying again from command line\n")
            # this is crazy. the Python solver reports unknown solution status,
            # even though the log clearly shows it converged. run the command
            # line solver and parse the output.
            #with NamedTemporaryFile() as t:
                #sh.mosek('mosek.opf', itro=t.name)
            filename = '/Users/srharnett/Dropbox/power/jabr-power-flow/wtf.sol'
            try:
                sh.mosek('mosek.opf', '-itro', filename, '-p', 'mosek.par',
                         '-q', 'wtf.log', _out='wtf.stdout', _err='wtf.stderr')
            except sh.ErrorReturnCode_22:
                out_file.write("\nreturn code 22, whatever that means\n")
            problem_status, solution_status, u, R, I = parse_mosek_output_file(filename)
            if solution_status not in ('OPTIMAL', 'NEAR_OPTIMAL'):
                raise SolutionError("mosek failed to converge: %s (check log)" % solution_status)
            return u, R, I


def parse_mosek_output_file(filename):
    u, R, I = [], [], []
    with open(filename) as f:
        line = f.readline()  # NAME
        problem_status = f.readline().split()[3].strip()
        solution_status = f.readline().split()[3].strip()
        while line != 'VARIABLES':
            line = f.readline().strip()
        f.readline()
        line = f.readline().split()
        var, val = line[1], float(line[3])
        while var[0] != 'R':
            u.append(val)
            line = f.readline().split()
            var, val = line[1], float(line[3])
        while var[0] != 'I':
            R.append(val)
            line = f.readline().split()
            var, val = line[1], float(line[3])
        while var[0] != 'c':
            I.append(val)
            line = f.readline().split()
            var, val = line[1], float(line[3])
    return problem_status, solution_status, u, R, I


def build_gurobi_model(case):
    G, B = case.G, case.B
    P = real(case.demands)
    Q = imag(case.demands)
    branches = case.branch_list
    n = len(case.demands)
    vhat = case.vhat
    s2 = 2**.5
    gens = {bus: gen.v for bus, gen in case.gens.items()}
    del gens[0]

    m = GurobiModel("jabr")
    u = [m.addVar(name='u_%d'%i) for i in range(n)]
    R = {(i, j): m.addVar(name='R_%d_%d' % (i, j)) for i, j in branches}
    I = {(i, j): m.addVar(lb=-GRB.INFINITY, name='I_%d_%d' % (i, j)) for i, j in branches}
    for i, j in branches:
        R[j, i] = R[i, j]
        I[j, i] = I[i, j]
    m.update()
    m.addConstr(u[0] == vhat*vhat/s2, 'u0')
    for gen, v in gens.iteritems():
        m.addConstr(u[gen] == v*v/s2, 'u%d' % gen)
    for i, j in branches:
        m.addQConstr(2*u[i]*u[j] >= R[i,j]*R[i,j] + I[i,j]*I[i,j], 'cone_%d_%d' % (i, j))
    k = lambda i: (j for j in B[i, :].nonzero()[1])
    s = lambda i, j: 1 if i < j else -1
    for i in range(1, n):
        m.addConstr(-s2*u[i]*G[i, :].sum() + quicksum(G[i,j]*R[i,j] + B[i,j]*s(i,j)*I[i,j] for j in k(i)) == P[i],
                    'real_flow_%d_%d' % (i, j))
        if i in gens:
            continue
        m.addConstr(s2*u[i]*B[i, :].sum() + quicksum(-B[i,j]*R[i,j] + G[i,j]*s(i,j)*I[i,j] for j in k(i)) == Q[i],
                    'reac_flow_%d_%d' % (i, j))
    m.setObjective(quicksum(R[i,j] for i, j in branches), sense=GRB.MAXIMIZE)
    m.params.outputFlag = 0
    m.optimize()
    if m.status != 2:
        raise SolutionError("gurobi failed to converge: %s (check log)" % m.status)
    u_opt = [x.getAttr('x') for x in u]
    R_opt = {(i, j): x.getAttr('x') for (i, j), x in R.items()}
    I_opt = {(i, j): x.getAttr('x') for (i, j), x in I.items()}
    #for (i, j), t in I_opt.items():
        #if j < i:
            #I_opt[(i, j)] = -t
    return u_opt, R_opt, I_opt


def recover_original_variables(u, I):
    """ given Jabr variables u and I, return bus voltages and angles """
    V = sqrt(sqrt(2) * array(u))
    theta_branch = {(i, j): asin(I[(i, j)]/(V[i]*V[j])) for (i, j) in I}
    theta_bus = recover_bus_angles(theta_branch)
    return V, theta_bus


def recover_bus_angles(theta_branch):
    """ assumes 0 is the root. traverses the tree and converts branch voltage
        angles to bus voltage angles """
    A = defaultdict(list) # adjacency matrix
    theta_bus = {0: 0}
    for i, j in theta_branch.keys():
        A[i] += [j]
    q = deque([0])
    while q:
        parent = q.popleft()
        for child in A[parent]:
            if child not in theta_bus:
                theta_ij = theta_branch[(parent, child)]
                if parent > child:
                    theta_ij *= -1
                theta_bus[child] = theta_bus[parent] - theta_ij
                q.append(child)
    return theta_bus


if __name__ == '__main__':
    process_answer(jabr5())