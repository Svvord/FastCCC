"""Shared helpers for pathway-classification report panels."""

from __future__ import annotations

import pandas as pd


UNANNOTATED_CLASSIFICATIONS = {
    "",
    "unknown",
    "unannotated",
    "nan",
    "none",
    "na",
    "n/a",
}


def keep_annotated_classifications(
    df: pd.DataFrame,
    column: str = "classification",
) -> pd.DataFrame:
    """Return rows with a real pathway classification annotation."""
    if column not in df.columns:
        return df.iloc[0:0].copy()

    out = df.copy()
    labels = out[column].astype("string").str.strip()
    mask = labels.notna() & ~labels.str.lower().isin(UNANNOTATED_CLASSIFICATIONS)
    out = out.loc[mask].copy()
    out[column] = labels.loc[mask].astype(str)
    return out
