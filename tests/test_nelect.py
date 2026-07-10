"""Tests for common._next_nelect, the pure Newton-step helper extracted
from set_pot's optimization loop. See __init__.py for the branch logic
this is meant to mirror exactly (equal-nelect escape, tiny-gradient
escape, two-point vs. three-point gradient estimate, and the +/-0.75
clamp on oversized steps).
"""
import pytest

import common


def test_equal_last_two_nelect_perturbs_by_point_one():
    nel_data = {
        'nelect': [100.0, 100.0],
        'potential': [-0.3, -0.3],
    }
    new_nel = common._next_nelect(nel_data, desired_U=-0.5)
    assert new_nel == pytest.approx(100.1)


def test_tiny_gradient_denominator_perturbs_by_point_one():
    # two points, but nelect barely changed -> denominator ~ 0
    nel_data = {
        'nelect': [100.0, 100.00001],
        'potential': [-0.30, -0.31],
    }
    new_nel = common._next_nelect(nel_data, desired_U=-0.5)
    assert new_nel == pytest.approx(100.00001 + 0.1)


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
