from __future__ import division, print_function
from gurobipy import Model as GurobiModel, GRB, quicksum
from numpy import array, sqrt, real, imag, pi
from math import asin
from scipy.sparse import dok_matrix, hstack
from collections import defaultdict, deque
from loadcase import load_case


def build_U_matrices(G, B):
    S2 = sqrt(2)
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
    #m.params.barQCPConvTol = 5e-10
    m.optimize()
    if m.status != 2:
        raise ValueError("gurobi failed to converge: %s (check log)" % m.status)
    u_opt = [x.getAttr('x') for x in u]
    R_opt = {(i, j): x.getAttr('x') for (i, j), x in R.items()}
    I_opt = {(i, j): x.getAttr('x') for (i, j), x in I.items()}
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
    # TODO clarify what you're doing here
    A = defaultdict(list) # adjacency matrix
    theta_bus_dict = {0: 0}
    for i, j in theta_branch.keys():
        A[i] += [j]
    q = deque([0])
    while q:
        parent = q.popleft()
        for child in A[parent]:
            if child not in theta_bus_dict:
                theta_ij = theta_branch[(parent, child)]
                if parent > child:
                    theta_ij *= -1
                theta_bus_dict[child] = theta_bus_dict[parent] - theta_ij
                q.append(child)
    theta_bus = [theta_bus_dict[i] for i in range(max(theta_bus_dict.keys())+1)]
    return theta_bus


def solve(casefile):
    """ given a matpower casefile, solves the power flow using the Jabr method
        and returns a dictionary mapping bus number to
        (voltage magnitude, voltage angle) tuples. angles are in degrees
    """
    case = load_case('cases/case5_renumber_tree.m')
    u, R, I = build_gurobi_model(case)
    V, theta = recover_original_variables(u, I)
    i2e = case.i2e
    answer = {i2e[i]: (v, t) for (i, (v, t)) in enumerate(zip(V, theta))}
    return answer


if __name__ == '__main__':
    answer = solve('cases/case5_renumber_tree.m')
    for bus in sorted(answer.keys()):
        v, t = answer[bus]
        print('%3d %7.3f %7.3f' % (bus, v, 180/pi*t))
