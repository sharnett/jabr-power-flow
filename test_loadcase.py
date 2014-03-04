import pytest
from loadcase import load_buses, renumber_buses, load_branches
from jabr import z2y

@pytest.fixture
def case5():
    return open('/Users/srharnett/Dropbox/power/jabr-power-flow/case5_renumber_tree.m')


@pytest.fixture
def case14():
    return open('/Users/srharnett/Dropbox/power/jabr-power-flow/case14_tree.m')


def test_load_buses_small(case5):
    demand_dict = {1: 1, 2: 1, 3: 1, 4: 1, 5: 0}
    root = 5
    vhat = 1
    assert load_buses(case5) == (demand_dict, root, vhat)


def test_load_buses_medium(case14):
    demand_dict = {1: 0, 2: 21.7, 3: 94.2, 4: 47.8, 5: 7.6, 6: 11.2, 7: 0, 8: 0,
                   9: 29.5, 10: 9, 11: 3.5, 12: 6.1, 13: 13.5, 14: 14.9}
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
    demands = [0, 21.7, 94.2, 47.8, 7.6, 11.2, 0, 0, 29.5, 9, 3.5, 6.1, 13.5, 14.9]
    assert renumber_buses(demand_dict, root) == (e2i, i2e, demands)


def test_load_branches_small(case5):
    demand_dict, root, _ = load_buses(case5)
    e2i, _, _ = renumber_buses(demand_dict, root)
    r, x = .01, .1
    g, b = z2y(r, x)
    y = g + 1j*b
    branches = {(1, 3): y, (0, 2): y, (3, 4): y, (0, 3): y}
    assert load_branches(case5, e2i) == branches


def test_load_branches_medium(case14):
    demand_dict, root, _ = load_buses(case14)
    e2i, _, _ = renumber_buses(demand_dict, root)
    r, x = .01, .1
    g, b = z2y(r, x)
    y = g + 1j*b
    branches = {(1, 3): y, (0, 2): y, (3, 4): y, (0, 3): y}
    assert load_branches(case5, e2i) == branches
