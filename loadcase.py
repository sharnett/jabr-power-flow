from jabr import z2y
from scipy.sparse import dok_matrix


def load_buses(casefileobj):
    """ returns a dictionary of bus demands, the root bus, and its voltage
    uses external bus numbering
    """
    line = ''
    demand_dict = {}
    root = 1
    vhat = 1
    while line.find("mpc.bus = [") == -1:
        line = casefileobj.readline()
    line = casefileobj.readline()
    while line.find("];") == -1:
        words = line.split()
        bus, bus_type, p, q = (int(words[0]), int(words[1]), float(words[2]),
                               float(words[3]))
        if bus_type == 3:
            root = bus
            vhat = float(words[7])
        demand_dict[bus] = p + q*1j
        line = casefileobj.readline()
    return demand_dict, root, vhat


def renumber_buses(demand_dict, root):
    """ creates a map of external to internal bus numbering, and a list with the
    reverse map. also returns the demands as a list using internal numbering
    """
    e2i = {root: 0}
    i2e = [0]*len(demand_dict)
    demands = [0]*len(demand_dict)
    i2e[0] = root
    i = 1
    for bus in sorted(demand_dict.keys()):
        if bus != root:
            e2i[bus] = i
            i2e[i] = bus
            i += 1
    for bus, d in demand_dict.iteritems():
        i = e2i[bus]
        demands[i] = d
    return e2i, i2e, demands


def load_gens(casefileobj, e2i):
    """ returns a dictionary of generator buses and their power output. internal
    numbering """
    line = ''
    gens = {}
    while line.find("mpc.gen = [") == -1:
        line = casefileobj.readline()
    line = casefileobj.readline()
    while line.find("];") == -1:
        words = line.split()
        bus, p, q = e2i[int(words[0])], float(words[1]), float(words[2])
        gens[bus] = p+q*1j
        line = casefileobj.readline()
    return gens


def adjust_demands(demands, gens):
    """ subtracts off any power generated at non-slack buses. internal numbering
    """
    for (bus, p) in gens.iteritems():
        if bus != 0:
            demands[bus] -= p


def load_branches(casefileobj, e2i):
    """ returns 'dictionary of keys' sparse matrices for branch conductance and
    susceptance. also returns a list of (from_bus, to_bus) connections, and a
    dictionary that inverts the list. uses internal numbering

    assumes there are exactly n-1 branches
    """
    line = ''
    n = len(e2i)
    G = dok_matrix((n, n))
    B = dok_matrix((n, n))
    branch_list = [None]*(n-1)
    branch_map = {}
    i = 0
    while line.find("mpc.branch = [") == -1:
        line = casefileobj.readline()
    line = casefileobj.readline()
    while line.find("];") == -1:
        words = line.split()
        fbus, tbus = e2i[int(words[0])], e2i[int(words[1])]
        r, x = float(words[2]), float(words[3])
        g, b = z2y(r, x)
        #if g < 1e-8:
            #g = 1e-8
        G[fbus, tbus] = g
        G[tbus, fbus] = g
        B[fbus, tbus] = b
        B[tbus, fbus] = b
        branch_list[i] = (fbus, tbus) if fbus < tbus else (tbus, fbus)
        line = casefileobj.readline()
        i += 1
    assert i == n-1, "it doesn't look like there are exactly n-1 branches"
    branch_list.sort()
    for i, (fbus, tbus) in enumerate(branch_list):
        branch_map[(fbus, tbus)] = i
        branch_map[(tbus, fbus)] = i
    return G, B, branch_list, branch_map


class Case(object):
    pass


def load_case(casefile):
    """ returns list of demands, conductance and susceptance matrices, list of
    branches, map of branches to list index, root voltage, and internal to
    external numbering map. uses internal numbering
    """
    casefileobj = open(casefile)
    demand_dict, root, vhat = load_buses(casefileobj)
    e2i, i2e, demands = renumber_buses(demand_dict, root)
    gens = load_gens(casefileobj, e2i)
    adjust_demands(demands, gens)
    G, B, branch_list, branch_map = load_branches(casefileobj, e2i)
    c = Case()
    c.demands = demands
    c.G = G
    c.B = B
    c.branch_list = branch_list
    c.branch_map = branch_map
    c.vhat = vhat
    c.i2e = i2e
    return c

