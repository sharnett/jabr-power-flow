import pytest
from jabr import *
from loadcase import load_case
from numpy.testing import assert_almost_equal
from scipy.sparse import coo_matrix
from numpy import array


@pytest.fixture
def case5():
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case5_renumber_tree.m'
    return load_case(c)


@pytest.fixture
def case5q():
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case5_with_q.m'
    return load_case(c)


@pytest.fixture
def case9():
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case9_tree.m'
    return load_case(c)


@pytest.fixture
def case14():
    #c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case14_tree.m'
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case14_matpower.m'
    return load_case(c)


def test_z2y_simple():
    r, x = 0, 1
    assert z2y(r, x) == (0, -1)


def test_z2y_normal():
    r, x = .02, .2
    g, b = z2y(r, x)
    ghat, bhat = 0.495049504950495, -4.9504950495049505
    assert_almost_equal(g, ghat)
    assert_almost_equal(b, bhat)


def test_build_U_matrices(case5):
    G, B = case5.G, case5.B
    r, x = .01, .1
    g, b = z2y(r, x)
    S2 = 2**.5
    Ureal = dok_matrix(
        [[0., g*S2,    0,    0,       0],
         [0.,    0, g*S2,    0,       0],
         [0.,    0,    0, 3*g*S2,     0],
         [0.,    0,    0,    0,    g*S2]])
    Ureac = dok_matrix(
        [[0., -b*S2,     0,       0,     0],
         [0.,     0, -b*S2,       0,     0],
         [0.,     0,     0, -3*b*S2,     0],
         [0.,     0,     0,       0, -b*S2]])
    Ureal_hat, Ureac_hat = build_U_matrices(G, B)
    assert_almost_equal(Ureal_hat.todense(), Ureal.todense())
    assert_almost_equal(Ureac_hat.todense(), Ureac.todense())


def test_build_R_matrices(case5):
    G, B, branch_map = case5.G, case5.B, case5.branch_map
    r, x = .01, .1
    g, b = z2y(r, x)
    Rreal = dok_matrix(
        [[0,  0, -g,  0],
         [-g, 0,  0,  0],
         [0, -g, -g, -g],
         [0,  0,  0, -g]])
    Rreac = dok_matrix(
        [[0, 0, b, 0],
         [b, 0, 0, 0],
         [0, b, b, b],
         [0, 0, 0, b]])
    Rreal_hat, Rreac_hat = build_R_matrices(G, B, branch_map)
    assert_almost_equal(Rreal_hat.todense(), Rreal.todense())
    assert_almost_equal(Rreac_hat.todense(), Rreac.todense())


def test_build_I_matrices(case5):
    G, B, branch_map = case5.G, case5.B, case5.branch_map
    r, x = .01, .1
    g, b = z2y(r, x)
    Ireal = dok_matrix(
        [[0,  0, -b,  0],
         [b, 0,  0,  0],
         [0,  b,  b, -b],
         [0,  0,  0,  b]])
    Ireac = dok_matrix(
        [[0, 0, -g,  0],
         [g, 0,  0,  0],
         [0, g,  g, -g],
         [0, 0,  0,  g]])
    Ireal_hat, Ireac_hat = build_I_matrices(G, B, branch_map)
    assert_almost_equal(Ireal_hat.todense(), Ireal.todense())
    assert_almost_equal(Ireac_hat.todense(), Ireac.todense())


def test_build_contraint_matrix(case5):
    G, B, branch_map = case5.G, case5.B, case5.branch_map
    r, x = .01, .1
    g, b = z2y(r, x)
    S2 = 2**.5
    Areal = coo_matrix(
        [[0.,  g*S2,     0,      0,      0,  0,  0, -g,  0,  0,  0, -b,  0],
         [0.,     0,  g*S2,      0,      0, -g,  0,  0,  0,  b,  0,  0,  0],
         [0.,     0,     0, 3*g*S2,      0,  0, -g, -g, -g,  0,  b,  b, -b],
         [0.,     0,     0,      0,   g*S2,  0,  0,  0, -g,  0,  0,  0,  b]])
    Areac = coo_matrix(
        [[0., -b*S2,     0,       0,     0,  0,  0,  b,  0,  0,  0, -g,  0],
         [0.,     0, -b*S2,       0,     0,  b,  0,  0,  0,  g,  0,  0,  0],
         [0.,     0,     0, -3*b*S2,     0,  0,  b,  b,  b,  0,  g,  g, -g],
         [0.,     0,     0,       0, -b*S2,  0,  0,  0,  b,  0,  0,  0,  g]])
    Areal_hat, Areac_hat = build_constraint_matrix(G, B, branch_map)
    assert_almost_equal(Areal_hat.todense(), Areal.todense())
    assert_almost_equal(Areac_hat.todense(), Areac.todense())


def test_build_mosek_model(case5):
    U = [0.70710678, 0.51489306, 0.68559881, 0.53884305, 0.51489306]
    R = [0.97958314, 0.81977994, 0.73816875, 0.73816875]
    I = [0.1, 0.3, -0.1, 0.1]
    Uhat, Rhat, Ihat = build_mosek_model(case5)
    assert_almost_equal(U, Uhat)
    assert_almost_equal(R, Rhat)
    assert_almost_equal(I, Ihat)


def test_full_case5(case5):
    V = {5: 1, 1: 0.85332805, 2: 0.98467413, 3: 0.87294854, 4: 0.85332805}
    u, R, I = build_mosek_model(case5)
    Vhat = recover_original_variables(u, R, I)
    i2e = case5.i2e
    for bus, v in enumerate(Vhat):
        assert_almost_equal(V[i2e[bus]], v)


def test_full_case5q(case5q):
    V = {5: 1, 1: 0.79806, 2: 0.97430, 3: 0.82022, 4: 0.78499}
    u, R, I = build_mosek_model(case5q)
    Vhat = recover_original_variables(u, R, I)
    i2e = case5q.i2e
    for bus, v in enumerate(Vhat):
        assert_almost_equal(V[i2e[bus]], v, decimal=4)


def test_full_case9(case9):
    V = [0, 1, 1, 1, 0.95627, 0.93106, 0.98732, 0.91403, 0.95389, 0.90988]
    u, R, I = build_mosek_model(case9)
    Vhat = recover_original_variables(u, R, I)
    i2e = case9.i2e
    for bus, v in enumerate(Vhat):
        assert_almost_equal(V[i2e[bus]], v, decimal=4)


def test_full_case14(case14):
    V = [0, 1.06, 1.045, 1.01, 0.88095, 0.8999, 0.92693, 0.87004, 0.93897,
         0.82911, 0.8196, 0.9194, 0.91413, 0.90867, 0.78673]
    u, R, I = build_mosek_model(case14)
    Vhat = recover_original_variables(u, R, I)
    i2e = case14.i2e
    for bus, v in enumerate(Vhat):
        assert_almost_equal(V[i2e[bus]], v, decimal=4)
