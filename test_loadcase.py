import pytest
from loadcase import *
from jabr import z2y
from scipy.sparse import dok_matrix
from numpy import real, imag
from numpy.testing import assert_almost_equal

CASE_DIRECTORY = '/Users/srharnett/Dropbox/power/jabr-power-flow/cases/'

@pytest.fixture
def case5():
    return open(CASE_DIRECTORY + 'case5_renumber_tree.m')


@pytest.fixture
def case14():
    return open(CASE_DIRECTORY + 'case14_tree.m')


def test_load_buses_small(case5):
    demand_dict = {1: 1, 2: 1, 3: 1, 4: 1, 5: 0}
    root = 5
    vhat = 1
    assert load_buses(case5) == (demand_dict, root, vhat)


def test_load_buses_medium(case14):
    demand_dict = {1: 0, 2: 21.7+12.7j, 3: 94.2+19j, 4: 47.8-3.9j, 5: 7.6+1.6j,
                   6: 11.2+7.5j, 7: 0, 8: 0, 9: 29.5+16.6j, 10: 9+5.8j,
                   11: 3.5+1.8j, 12: 6.1+1.6j, 13: 13.5+5.8j, 14: 14.9+5j}
    root = 1
    vhat = 1.06
    assert load_buses(case14) == (demand_dict, root, vhat)


def test_renumber_buses_small(case5):
    demand_dict, root, _ = load_buses(case5)
    e2i = {5: 0, 1: 1, 2: 2, 3: 3, 4: 4}
    i2e = [5, 1, 2, 3, 4]
    demands = [0, 1, 1, 1, 1]
    assert renumber_buses(demand_dict, root) == (e2i, i2e, demands)


def test_renumber_buses_medium(case14):
    demand_dict, root, _ = load_buses(case14)
    e2i = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10,
           12: 11, 13: 12, 14: 13}
    i2e = range(1, 15)
    demands = [0, 21.7+12.7j, 94.2+19j, 47.8-3.9j, 7.6+1.6j, 11.2+7.5j, 0, 0,
               29.5+16.6j, 9+5.8j, 3.5+1.8j, 6.1+1.6j, 13.5+5.8j, 14.9+5j]
    assert renumber_buses(demand_dict, root) == (e2i, i2e, demands)


def test_load_gens_small(case5):
    demand_dict, root, _ = load_buses(case5)
    e2i, _, _ = renumber_buses(demand_dict, root)
    gens = {0: (4.1590465, 1.0)}
    gens_hat = load_gens(case5, e2i)
    for bus, gen in gens.items():
        assert gens_hat[bus].p == gen[0]
        assert gens_hat[bus].v == gen[1]


def test_load_gens_med(case14):
    demand_dict, root, _ = load_buses(case14)
    e2i, _, _ = renumber_buses(demand_dict, root)
    gens = {0: (238.679924, 1.06),
            1: (40, 1.045),
            2: (0, 1.01),
            5: (0, 1.07),
            7: (0, 1.09)}
    gens_hat = load_gens(case14, e2i)
    for bus, gen in gens.items():
        assert gens_hat[bus].p == gen[0]
        assert gens_hat[bus].v == gen[1]


def test_adjust_demands_small(case5):
    adjusted_demands = [0, 1, 1, 1, 1]
    demand_dict, root, _ = load_buses(case5)
    e2i, _, demands = renumber_buses(demand_dict, root)
    gens = load_gens(case5, e2i)
    adjust_demands(demands, gens)
    assert demands == adjusted_demands


def test_adjust_demands_medium(case14):
    adjusted_demands = [0, -18.3+12.7j, 94.2+19j, 47.8-3.9j,
                        7.6+1.6j, 11.2+7.5j, 0, 0, 29.5+16.6j,
                        9+5.8j, 3.5+1.8j, 6.1+1.6j, 13.5+5.8j, 14.9+5j]
    demand_dict, root, _ = load_buses(case14)
    e2i, _, demands = renumber_buses(demand_dict, root)
    gens = load_gens(case14, e2i)
    adjust_demands(demands, gens)
    assert_almost_equal(demands, adjusted_demands)


def test_load_branches_small(case5):
    demand_dict, root, _ = load_buses(case5)
    e2i, _, _ = renumber_buses(demand_dict, root)
    r, x = .01, .1
    g, b = z2y(r, x)
    n = len(e2i)
    Ghat = dok_matrix((n, n))
    Bhat = dok_matrix((n, n))
    keys = [(0, 2), (0, 3), (1, 3), (3, 4)]
    for i, j in keys:
        Ghat[i, j] = g
        Ghat[j, i] = g
        Bhat[i, j] = b
        Bhat[j, i] = b
    branch_list_hat = list(sorted(keys))
    branch_map_hat = {}
    for i, (fbus, tbus) in enumerate(branch_list_hat):
        branch_map_hat[(fbus, tbus)] = i
        branch_map_hat[(tbus, fbus)] = i
    G, B, branch_list, branch_map = load_branches(case5, e2i)
    assert_almost_equal(G.todense(), Ghat.todense())
    assert_almost_equal(B.todense(), Bhat.todense())
    assert branch_list == branch_list_hat
    assert branch_map == branch_map_hat


def test_load_branches_medium(case14):
    demand_dict, root, _ = load_buses(case14)
    e2i, _, _ = renumber_buses(demand_dict, root)
    n = len(e2i)
    Ghat = dok_matrix((n, n))
    Bhat = dok_matrix((n, n))
    s_dict = {(0, 1): (499.9131600798035-1526.3086523179554j),
              (0, 4): (102.58974549701888-423.4983682334831j),
              (1, 2): (113.50191923073959-478.1863151757718j),
              (3, 4): (684.0980661495671-2157.855398169159j),
              (3, 6): -478.1943381790359j,
              (4, 5): -396.79390524561546j,
              (5, 10): (195.50285631772607-409.4074344240442j),
              (5, 11): (152.59674404509738-317.5963965029401j),
              (5, 12): (309.89274038379875-610.2755448193116j),
              (6, 7): -567.6979846721543j,
              (6, 8): -909.0082719752751j,
              (8, 9): (390.2049552447428-1036.5394127060915j),
              (8, 13): (142.4005487019931-302.90504569306034j)}
    for (i, j), y in s_dict.iteritems():
        Ghat[i, j] = y.real
        Ghat[j, i] = y.real
        Bhat[i, j] = y.imag
        Bhat[j, i] = y.imag
    branch_list_hat = list(sorted(s_dict.keys()))
    branch_map_hat = {}
    for i, (fbus, tbus) in enumerate(branch_list_hat):
        branch_map_hat[(fbus, tbus)] = i
        branch_map_hat[(tbus, fbus)] = i
    G, B, branch_list, branch_map = load_branches(case14, e2i)
    assert_almost_equal(G.todense(), Ghat.todense())
    assert_almost_equal(B.todense(), Bhat.todense())
    assert branch_list == branch_list_hat
    assert branch_map == branch_map_hat


def test_load_case():
    casefile = CASE_DIRECTORY + 'case5_renumber_tree.m'
    vhat = 1
    i2e = [5, 1, 2, 3, 4]
    demands = [0, 1, 1, 1, 1]
    r, x = .01, .1
    g, b = z2y(r, x)
    n = len(demands)
    G = dok_matrix((n, n))
    B = dok_matrix((n, n))
    keys = [(0, 2), (0, 3), (1, 3), (3, 4)]
    for i, j in keys:
        G[i, j] = g
        G[j, i] = g
        B[i, j] = b
        B[j, i] = b
    c = load_case(casefile)
    assert c.demands == demands
    assert c.G == G
    assert c.B == B
    assert c.vhat == vhat
    assert c.i2e == i2e
