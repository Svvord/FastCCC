import pandas as pd

def filter_empty_genes(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Remove gene with all counts values to zero
    """
    if counts.empty:
        return counts
    
    select_columns = counts.apply(lambda col: col.sum() > 0)
    select_columns = counts.columns[select_columns]
    filtered_counts = counts[select_columns]

    return filtered_counts