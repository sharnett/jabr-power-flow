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
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case14_tree.m'
    return load_case(c)


@pytest.fixture
def case118():
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case118_v2.m'
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


def test_parse_mosek_output_file(filename='/Users/srharnett/Dropbox/power/jabr-power-flow/case14.out'):
    problem_status, solution_status, u_hat, _, _ = parse_mosek_output_file(filename)
    v = [1.06, 1.045, 1.01, 0.88095, 0.8999, 0.92693, 0.87004, 0.93897,
         0.82911, 0.8196, 0.9194, 0.91413, 0.90867, 0.78673]
    u = array(v)**2/sqrt(2)
    assert problem_status == 'PRIMAL_AND_DUAL_FEASIBLE'
    assert solution_status == 'OPTIMAL'
    assert_almost_equal(u, u_hat, decimal=4)


def test_full_case5(case5):
    V = {5: 1, 1: 0.85332805, 2: 0.98467413, 3: 0.87294854, 4: 0.85332805}
    u, R, I = build_mosek_model(case5)
    Vhat = recover_original_variables(u, R, I)
    i2e = case5.i2e
    for bus, v in enumerate(Vhat):
        assert_almost_equal(V[i2e[bus]], v, decimal=4)


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


def test_full_case118(case118):
    V = [0, 0.955, 0.98047, 0.96243, 0.998, 1.0032, 0.99, 0.98869, 1.015, 1.023,
         1.05, 0.98107, 0.99, 0.96077, 0.98618, 0.97, 0.95675, 0.98835, 0.973,
         0.962, 0.91971, 0.9262, 0.94582, 0.99291, 0.992, 1.05, 1.015, 0.968,
         0.94859, 0.93822, 0.96681, 0.967, 0.963, 0.96438, 0.984, 0.97959, 0.98,
         0.98799, 0.94192, 0.96662, 0.97, 0.95377, 0.985, 0.88745, 0.92024,
         0.94488, 1.005, 1.0021, 1.016, 1.025, 1.0118, 0.97702, 0.97017,
         0.93377,  0.955, 0.952, 0.954, 1.002, 0.97162, 0.985, 0.99246, 0.995,
         0.998, 0.96269, 0.97828, 1.005, 1.05, 1.0368, 0.99932, 1.035, 0.984,
         0.98626, 0.98, 0.991, 0.958, 0.95473, 0.943, 1.006, 1.0001, 1.0115,
         1.04, 0.98881, 0.96363, 0.96244, 0.97699, 0.985, 0.98392, 1.015,
         0.99053, 1.005, 0.985, 0.98, 0.99, 0.98075, 0.96795, 0.9478, 0.99708,
         1.0291, 1.0231, 1.01, 1.017, 0.95209, 0.97567, 1.01, 0.971, 0.965,
         0.94646, 0.952, 0.95931, 0.95753, 0.973, 0.98, 0.975, 0.993, 0.95995,
         0.9587, 1.005, 0.97135, 0.94188]
    u, R, I = build_mosek_model(case118)
    Vhat = recover_original_variables(u, R, I)
    Vhat_i = [0]*(len(Vhat)+1)
    i2e = case118.i2e
    for bus, v in enumerate(Vhat):
        Vhat_i[i2e[bus]] = v
    assert_almost_equal(V, Vhat_i)
