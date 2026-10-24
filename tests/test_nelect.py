"""Tests for common._next_nelect, the pure Newton-step helper extracted
from set_pot's optimization loop. See __init__.py for the branch logic
this is meant to mirror exactly (equal-nelect escape, tiny-gradient
escape, two-point vs. three-point gradient estimate, and the +/-0.75
clamp on oversized steps).
"""
import pytest

import common


def test_equal_last_two_nelect_newton_step_fallback_slope():
    # post-relax duplicate, no distinct earlier point: fallback 3.2 e/V,
    # U above target -> more electrons (U falls as N rises)
    nel_data = {'nelect': [100.0, 100.0], 'potential': [-0.3, -0.3]}
    assert common._next_nelect(nel_data, -0.5) == pytest.approx(100.0 + 3.2*0.2)
    # U below target -> fewer electrons
    assert common._next_nelect(nel_data, -0.1) == pytest.approx(100.0 - 3.2*0.2)


def test_equal_last_two_nelect_uses_distinct_slope():
    # last distinct pair: (101, 0.06) vs (100, 0.30) -> dN/dU = 1/(-0.24)
    nel_data = {'nelect': [100.0, 101.0, 101.0], 'potential': [0.30, 0.05, 0.06]}
    dNdU = 1.0/(0.06-0.30)
    expected = 101.0 - (0.06-0.0)*dNdU
    assert common._next_nelect(nel_data, 0.0) == pytest.approx(expected)


def test_equal_last_two_nelect_slope_out_of_bounds_uses_fallback():
    # dN/dU = -0.5 e/V (too small) -> fallback 3.2; custom fallback honoured
    nel_data = {'nelect': [100.0, 101.0, 101.0], 'potential': [0.30, -1.70, -1.70]}
    assert common._next_nelect(nel_data, -1.9) == pytest.approx(101.0 + 3.2*0.2)
    assert common._next_nelect(nel_data, -1.9, fallback_slope=2.0) == pytest.approx(101.0 + 2.0*0.2)


def test_tiny_gradient_denominator_uses_fallback_step():
    # two points, but nelect barely changed -> denominator ~ 0
    nel_data = {'nelect': [100.0, 100.00001], 'potential': [-0.30, -0.31]}
    assert common._next_nelect(nel_data, -0.5) == pytest.approx(100.00001 + 3.2*0.19)


def test_two_point_newton_step():
    nel_data = {
        'nelect': [100.0, 101.0],
        'potential': [-0.30, -0.40],  # dV/dN = -0.10 per electron
    }
    desired_U = -0.50
    grad = (-0.30 - -0.40) / (100.0 - 101.0)  # matches inline formula: grad_numer/grad_denom
    y = nel_data['potential'][-1] - desired_U
    expected_diff = abs(y) ** 2 / (y * grad)
    expected = nel_data['nelect'][-1] - expected_diff

    new_nel = common._next_nelect(nel_data, desired_U)
    assert new_nel == pytest.approx(expected)


def test_three_point_newton_step_uses_linear_fit():
    # exact linear V(N) relationship over 3 points: slope -0.1, intercept 9.7
    nelect = [95.0, 100.0, 101.0]
    potential = [0.0, -0.5, -0.6]
    nel_data = {'nelect': nelect, 'potential': potential}
    desired_U = -0.55

    xax, yax, grad, intercept = common.get_line(nelect, potential)
    y = potential[-1] - desired_U
    diff = abs(y) ** 2 / (y * grad)
    expected = nelect[-1] - diff

    new_nel = common._next_nelect(nel_data, desired_U)
    assert new_nel == pytest.approx(expected)


def test_clamp_positive_diff():
    # slope grad = 0.1; with desired_U far below the last potential, the
    # raw Newton step diff = y/grad = 11.0, well above the >5.0 clamp
    # threshold, so it should be replaced with the fixed 0.75 step.
    nel_data = {
        'nelect': [100.0, 101.0],
        'potential': [0.0, 0.1],
    }
    desired_U = -1.0
    new_nel = common._next_nelect(nel_data, desired_U)
    assert new_nel == pytest.approx(101.0 - 0.75)


def test_clamp_negative_diff():
    # same slope, but desired_U far above the last potential now makes
    # diff = y/grad = -19.0, below the <-5.0 clamp threshold, so it
    # should be replaced with the fixed -0.75 step.
    nel_data = {
        'nelect': [100.0, 101.0],
        'potential': [0.0, 0.1],
    }
    desired_U = 2.0
    new_nel = common._next_nelect(nel_data, desired_U)
    assert new_nel == pytest.approx(101.0 + 0.75)


def test_min_step_guard():
    # a huge fitted gradient would give a ~1e-5 e step; floor to 1e-3 e
    nel_data = {'nelect': [100.0, 101.0], 'potential': [0.0, 1000.0]}
    new_nel = common._next_nelect(nel_data, 999.97)
    assert new_nel - 101.0 == pytest.approx(1e-3)  # y > 0 -> more electrons
    nel_data = {'nelect': [100.0, 101.0], 'potential': [0.0, 1000.0]}
    assert 101.0 - common._next_nelect(nel_data, 1000.03) == pytest.approx(1e-3)


def test_nan_and_zero_gradient_edge_cases():
    import math
    # flat potential -> zero gradient -> fallback Newton step, finite
    nel_data = {'nelect': [100.0, 101.0, 102.0], 'potential': [0.3, 0.3, 0.3]}
    assert common._next_nelect(nel_data, 0.0) == pytest.approx(102.0 + 3.2*0.3)
    # NaN potential propagates (set_pot exits on NaN)
    nel_data = {'nelect': [100.0, 101.0], 'potential': [0.3, float('nan')]}
    assert math.isnan(common._next_nelect(nel_data, 0.0))


def test_identical_history_never_returns_last():
    nel_data = {'nelect': [100.0, 100.0, 100.0], 'potential': [0.3, 0.3, 0.3]}
    for U in (0.0, 0.29, 0.31, 0.6):
        assert abs(common._next_nelect(nel_data, U)-100.0) >= 1e-3


def _simulate(slope, N0, U0, target, hist):
    """Linear U(N) = U0 + slope*(N-N0) (V per e, slope < 0); run the Newton
    loop from a post-relax duplicate entry until within 0.02 V."""
    nd = {'nelect': list(hist[0]), 'potential': list(hist[1])}
    for it in range(10):
        if abs(nd['potential'][-1]-target) < 0.02:
            return it, nd
        n = common._next_nelect(nd, target)
        nd['nelect'].append(n)
        nd['potential'].append(U0 + slope*(n-N0))
    return 10, nd


@pytest.mark.parametrize('slope', [-0.314, -0.48, -0.2])
@pytest.mark.parametrize('side', [+0.15, -0.15])
def test_linear_model_converges_from_post_relax_duplicate(slope, side):
    N0, U0, target = 1000.0, 1.0, 1.2
    n1 = N0 + (target+side-U0)/slope
    hist = ([N0, n1, n1], [U0, target+side, target+side])
    it, nd = _simulate(slope, N0, U0, target, hist)
    assert it <= 2
    assert len(set(nd['nelect'])) == len(nd['nelect']) - 1  # only the relax duplicate
