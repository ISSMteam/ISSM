# compare_models.py
from __future__ import annotations
import numbers
import numpy as np
from typing import Any, Iterable, Mapping

def compare_models(md1: Any, md2: Any, *, quiet: bool = False) -> list[str]:
    """
    Recursively compare two ISSM model objects (or nested objects/arrays).

    - Enforces same class for nested objects (skips subtree if different).
    - For numeric arrays/scalars: requires identical shape and values.
      NaNs at identical positions are considered equal (MATLAB-like logic).
    - Handles lists/tuples and dicts recursively.
    - Returns a list of human-readable difference messages.
      If quiet=False (default), also prints them.

    Example:
        diffs = compare_models(md1, md2)
        if diffs:
            print("\\n".join(diffs))
    """
    diffs: list[str] = []

    def _is_numeric_like(x: Any) -> bool:
        return isinstance(x, (numbers.Number, np.ndarray)) or (
            hasattr(x, "dtype") and hasattr(x, "shape")
        )

    def _as_ndarray(x: Any) -> np.ndarray:
        if isinstance(x, np.ndarray):
            return x
        return np.asarray(x)

    def _get_fields(obj: Any) -> dict[str, Any]:
        """
        Safely list data attributes of an ISSM Python object or plain object.
        - Skips dunder/private names and callables.
        - Uses __dict__ when available (faster, fewer false positives).
        - Falls back to dir() + getattr() for SWIG/wrapped objects.
        """
        fields: dict[str, Any] = {}
        # Prefer __dict__ if present
        src: Iterable[str]
        if hasattr(obj, "__dict__") and isinstance(obj.__dict__, dict):
            src = obj.__dict__.keys()
        else:
            src = [n for n in dir(obj) if not (n.startswith("_"))]

        for name in src:
            if name.startswith("_"):
                continue
            # Avoid bound methods / functions / properties that raise
            try:
                val = getattr(obj, name)
            except Exception:
                continue
            if callable(val):
                continue
            # Some ISSM wrappers expose numpy arrays via properties; keep those
            fields[name] = val
        return fields

    def _compare(path: str, a: Any, b: Any):
        # Case 1: numeric scalars/arrays — MATLAB-like behavior
        if _is_numeric_like(a) and _is_numeric_like(b):
            A = _as_ndarray(a)
            B = _as_ndarray(b)
            if A.shape != B.shape:
                diffs.append(f"{path} do not have the same size: {A.shape} vs {B.shape}")
                return

            # Equality with NaN-at-same-indices treated as equal
            # Build a mask where both are NaN; set those to zero for comparison
            a_nan = np.isnan(A)
            b_nan = np.isnan(B)
            # If NaN positions differ, it's a diff
            if A.dtype.kind in "fc" or B.dtype.kind in "fc":
                if not np.array_equal(a_nan, b_nan):
                    diffs.append(f"{path} differs (NaN pattern not identical)")
                    return
                # Replace NaNs by zero and compare
                A_cmp = np.where(a_nan, 0, A)
                B_cmp = np.where(b_nan, 0, B)
            else:
                # Non-float dtypes: just compare directly
                A_cmp, B_cmp = A, B

            if not np.array_equal(A_cmp, B_cmp):
                diffs.append(f"{path} differs")
            return

        # Case 2: strings (exact compare)
        if isinstance(a, str) and isinstance(b, str):
            if a != b:
                diffs.append(f"{path} differs: {a!r} != {b!r}")
            return

        # Case 3: mappings (dict-like)
        if isinstance(a, Mapping) and isinstance(b, Mapping):
            a_keys = set(a.keys()); b_keys = set(b.keys())
            only_a = sorted(a_keys - b_keys)
            only_b = sorted(b_keys - a_keys)
            for k in only_a:
                diffs.append(f"{path}.{k} present only in first")
            for k in only_b:
                diffs.append(f"{path}.{k} present only in second")
            for k in sorted(a_keys & b_keys):
                _compare(f"{path}.{k}", a[k], b[k])
            return

        # Case 4: lists/tuples — compare length then elementwise
        if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
            if len(a) != len(b):
                diffs.append(f"{path} do not have the same size: {len(a)} vs {len(b)}")
                return
            for i, (ai, bi) in enumerate(zip(a, b)):
                _compare(f"{path}[{i}]", ai, bi)
            return

        # Case 5: objects — recurse on attributes if same class, else skip (MATLAB behavior)
        if hasattr(a, "__class__") and hasattr(b, "__class__"):
            if a.__class__ is not b.__class__:
                diffs.append(
                    f"Skipping '{path}' because classes are not consistent "
                    f"({a.__class__.__name__} vs {b.__class__.__name__})"
                )
                return
            # Recurse into fields
            a_fields = _get_fields(a)
            b_fields = _get_fields(b)
            keys = sorted(set(a_fields) | set(b_fields))
            for k in keys:
                in_a = k in a_fields
                in_b = k in b_fields
                if not in_a or not in_b:
                    which = "first" if in_a else "second"
                    diffs.append(f"{path}.{k} present only in {which}")
                    continue
                _compare(f"{path}.{k}", a_fields[k], b_fields[k])
            return

        # Fallback: strict equality for other types
        if a != b:
            diffs.append(f"{path} differs: {a!r} != {b!r}")

    _compare("md", md1, md2)

    if not quiet:
        for d in diffs:
            print(d)
    return diffs
