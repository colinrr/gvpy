"""Generic value-validation helpers, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: @IceCauldron/IceCauldron.m's local validator
functions (mustBePosOrNan, mustBeInRangeInclusive), plus MATLAB's built-in
property validators used throughout that same properties block
(mustBeNonnegative, mustBePositive).
"""

from collections.abc import Callable
from typing import Any

import attrs
import numpy as np
from numpy.typing import ArrayLike

# attrs validator signature: (instance, attribute, value) -> None
AttrsValidator = Callable[[Any, attrs.Attribute, Any], None]


def must_be_pos_or_nan(x: ArrayLike) -> None:
    """Value must be positive or NaN."""
    x = np.asarray(x)
    if not np.all((x > 0) | np.isnan(x)):
        raise ValueError("Value must be positive or NaN.")


def must_be_nonnegative(x: ArrayLike) -> None:
    """Value must be non-negative."""
    x = np.asarray(x)
    if not np.all(x >= 0):
        raise ValueError("Value must be non-negative.")


def must_be_positive(x: ArrayLike) -> None:
    """Value must be positive."""
    x = np.asarray(x)
    if not np.all(x > 0):
        raise ValueError("Value must be positive.")


def must_be_in_range_inclusive(x: ArrayLike, a: float, b: float) -> None:
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
def attrs_validate_pos_or_nan(instance: object, attribute: attrs.Attribute, value: ArrayLike) -> None:
    must_be_pos_or_nan(value)


def attrs_validate_nonnegative(instance: object, attribute: attrs.Attribute, value: ArrayLike) -> None:
    must_be_nonnegative(value)


def attrs_validate_positive(instance: object, attribute: attrs.Attribute, value: ArrayLike) -> None:
    must_be_positive(value)


def attrs_range_inclusive(a: float, b: float) -> AttrsValidator:
    """Build an attrs validator enforcing value in range [a, b]."""

    def validator(instance: object, attribute: attrs.Attribute, value: ArrayLike) -> None:
        must_be_in_range_inclusive(value, a, b)

    return validator
