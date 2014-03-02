from jabr import z2y
from numpy.testing import assert_almost_equal

def test_true():
    assert True


def test_false():
    assert False


def test_z2y_simple():
    r, x = 0, 1
    assert z2y(r, x) == (0, -1)


def test_z2y_normal():
    r, x = .02, .2
    g, b = z2y(r, x)
    ghat, bhat = 0.495049504950495, -4.9504950495049505
    assert_almost_equal(g, ghat)
    assert_almost_equal(b, bhat)

