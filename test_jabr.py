import pytest
from jabr import *
from loadcase import load_case
from numpy.testing import assert_almost_equal
from mosek.fusion import DenseMatrix


@pytest.fixture
def case5():
    c = '/Users/srharnett/Dropbox/power/jabr-power-flow/case5_renumber_tree.m'
    demands, susceptances, vhat, i2e = load_case(c)
    return demands, susceptances, vhat, i2e


def test_z2y_simple():
    r, x = 0, 1
    assert z2y(r, x) == (0, -1)


def test_z2y_normal():
    r, x = .02, .2
    g, b = z2y(r, x)
    ghat, bhat = 0.495049504950495, -4.9504950495049505
    assert_almost_equal(g, ghat)
    assert_almost_equal(b, bhat)


def test_build_U_matrix(case5):
    susceptances = case5[1]
    r, x = .01, .1
    g, b = z2y(r, x)
    S2 = 2**.5
    Ureal = DenseMatrix(
        [[0.,  g*S2,       0,     0,     0],
         [0.,     0,  3*g*S2,     0,     0],
         [0.,     0,       0,  g*S2,     0],
         [0.,     0,       0,     0,  g*S2]])
    assert build_U_matrix(susceptances) == Ureal


def test_build_contraint_matrix(case5):
    susceptances = case5[1]
    r, x = .01, .1
    g, b = z2y(r, x)
    S2 = 2**.5
    Areal = DenseMatrix(
        [[0.,  g*S2,       0,     0,     0, -g,  0,  0,  0, -b,  0,  0,  0],
         [0.,     0,  3*g*S2,     0,     0,  0, -g, -g, -g,  0, -b,  b,  b],
         [0.,     0,       0,  g*S2,     0,  0,  0, -g,  0,  0,  0, -b,  0],
         [0.,     0,       0,     0,  g*S2,  0,  0,  0, -g,  0,  0,  0, -b]])
    assert build_constraint_matrix(susceptances) == Areal
