from jabr import z2y


def load_buses(casefileobj):
    line = ''
    demand_dict = {}
    root = 1
    vhat = 1
    while line.find("mpc.bus = [") == -1:
        line = casefileobj.readline()
    line = casefileobj.readline()
    while line.find("];") == -1:
        words = line.split()
        bus, bus_type, d = int(words[0]), int(words[1]), float(words[2])
        if bus_type == 3:
            root = bus
            vhat = float(words[7])
        demand_dict[bus] = d
        line = casefileobj.readline()
    return demand_dict, root, vhat


def renumber_buses(demand_dict, root):
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


def load_branches(casefileobj, e2i):
    line = ''
    susceptances = {}
    while line.find("mpc.branch = [") == -1:
        line = casefileobj.readline()
    line = casefileobj.readline()
    while line.find("];") == -1:
        words = line.split()
        fbus, tbus = e2i[int(words[0])], e2i[int(words[1])]
        r, x = float(words[2]), float(words[3])
        if tbus < fbus:
            fbus, tbus = tbus, fbus
        g, b = z2y(r, x)
        susceptances[(fbus, tbus)] = g+b*1j
        line = casefileobj.readline()
    return susceptances
