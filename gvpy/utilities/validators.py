"""Generic value-validation helpers, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: @IceCauldron/IceCauldron.m's local validator
functions (mustBePosOrNan, mustBeInRangeInclusive), plus MATLAB's built-in
property validators used throughout that same properties block
(mustBeNonnegative, mustBePositive).
"""

import numpy as np


def must_be_pos_or_nan(x):
    """Value must be positive or NaN."""
    x = np.asarray(x)
    if not np.all((x > 0) | np.isnan(x)):
        raise ValueError("Value must be positive or NaN.")


def must_be_nonnegative(x):
    """Value must be non-negative."""
    x = np.asarray(x)
    if not np.all(x >= 0):
        raise ValueError("Value must be non-negative.")


def must_be_positive(x):
    """Value must be positive."""
    x = np.asarray(x)
    if not np.all(x > 0):
        raise ValueError("Value must be positive.")


def must_be_in_range_inclusive(x, a, b):
    """Value must be in range [a, b]."""
    x = np.asarray(x)
    if not np.all((x >= a) & (x <= b)):
        raise ValueError(f"Value must be in range [{a:.2f}, {b:.2f}].")


# TRANSLATION NOTE: the four functions above are the direct equivalents of MATLAB's declarative
# per-field validators ({mustBePosOrNan}, {mustBeNonnegative}, {mustBePositive},
# {mustBeInRangeInclusive(...)}) used throughout IceCauldron.m's properties block (~28 fields).
# The adapters below bridge them into attrs' validator signature (instance, attribute, value)
# so they can be attached per-field via attrs.field(validator=...) on IceCauldron, preserving
# the same declarative, per-field structure (and line-for-line traceability to the MATLAB
# source) that hand-written checks in __post_init__ would not. See gvpy/ice_cauldron.py.
def attrs_validate_pos_or_nan(instance, attribute, value):
    must_be_pos_or_nan(value)


def attrs_validate_nonnegative(instance, attribute, value):
    must_be_nonnegative(value)


def attrs_validate_positive(instance, attribute, value):
    must_be_positive(value)


def attrs_range_inclusive(a, b):
    """Build an attrs validator enforcing value in range [a, b]."""

    def validator(instance, attribute, value):
        must_be_in_range_inclusive(value, a, b)

    return validator


def allow_none(validate):
    """Wrap an attrs validator so it skips fields that are still unset (None) -
    for computed fields (e.g. E_prime) that start as None and are only assigned
    a real value later, in __attrs_post_init__.
    """

    def validator(instance, attribute, value):
        if value is None:
            return
        validate(instance, attribute, value)

    return validator
