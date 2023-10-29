import glob
import os
import pickle
import numpy as np
from collections import Counter
import pandas as pd


def get_sheet_name(result_type: str) -> str:
    """Determine the sheet name based on the result type."""
    return {
        "sum": "sum",
        "tetrad": "tetrad",
        "triad": "triad"
    }.get(result_type, "triad")

def update_sheet_with_results(sheet: pd.DataFrame, prefix: str, result_dict: Dict[str, Any], detected_tetrad: Dict[str, Any], output_extent: str) -> pd.DataFrame:
    """
    Update the provided sheet with results.
    
    Args:
    - sheet: DataFrame representing the sheet to be updated.
    - prefix: Prefix string to identify rows.
    - result_dict: Dictionary containing result values.
    - detected_tetrad: Dictionary containing tetrad values.
    - output_extent: Extent specification for output.
    
    Returns:
    - Updated DataFrame.
    """
    
    # Ensure 'Prefix' column exists in the sheet
    if 'Prefix' not in sheet.columns:
        sheet['Prefix'] = ''
    # Map of limited columns and their corresponding verbose names
    limited_columns = ['cM_rb', 'cM_ry', 'cM_yb','significance_map_distance_rb','X2_p_rb','F_rb','total_classified','total_predicted','MI_frequency','MII_frequency','CoC_1','sum_gc_gain','sum_gc_loss']
    column_map = {'cM_rb' : "Centimorgan (Red-Blue)",
                 'cM_ry': "Centimorgan (Red-Yellow)",
                 'cM_yb': "Centimorgan (Yellow-Blue)",
                 'significance_map_distance_rb': "Significance of Differences (wt:mutant)",
                 'X2_p_rb': "interference p-value",
                 'F_rb': "Papazzian Interfernce (Red-Blue)",
                 'total_classified': "Total Classified",
                  'total_predicted': "Total Predicted",
                  'MI_frequency': "MI NonDisjunction (frequency)",
                  'MII_frequency': "MII NonDisjunction (frequency)",
                  'I1': "Coefficient of Coincidence (Red-Blue)",
                  'sum_gc_gain': "Gene Conversion with Gain of marker (frequency)",
                  'sum_gc_loss': "Gene Conversion with Loss of marker (frequency)",
                 }
    
    # If output_extent is limited, filter the columns and rename them
    if output_extent == 'limited':
    if output_extent == 'limited':
        result_dict = {k: v for k, v in result_dict.items() if k in limited_columns}
        detected_tetrad = {}  # Empty the detected_tetrad dictionary
        sheet = sheet.rename(columns=column_map)
    
    # Check if prefix already exists in the sheet
    hit = sheet.Prefix.astype(str).str.match('^' + prefix + '$')
    if hit.sum() == 0:
        new_row = {'Prefix': prefix}
        new_row.update(result_dict)
        new_row.update(detected_tetrad)
        sheet = pd.concat([sheet, pd.DataFrame([new_row])], ignore_index=True)

    else:
        row_index = sheet.index[hit][0]
        
        for dictionary in [result_dict, detected_tetrad]:
            for col, value in dictionary.items():
                if isinstance(value, bool) and sheet[col].dtype != 'boolean':
                    sheet[col] = sheet[col].astype('boolean')
                sheet.at[row_index, col] = value
            
    return sheet

def output_results(input_df, root_path, a_prefix, detected_tetrad, result_dict, prefix, output_extent, result_type):
    """
    Output results to the specified sheet within the input DataFrame.
    
    Args:
    - input_df: Dictionary with sheet names as keys and DataFrames as values.
    - root_path: Path where results should be saved.
    - a_prefix: Prefix for annotations.
    - detected_tetrad: Dictionary containing tetrad values.
    - result_dict: Dictionary containing result values.
    - prefix: Prefix string to identify rows.
    - output_extent: Extent specification for output.
    - result_type: Type of results ('sum', 'tetrad', or 'triad').
    """
    sheet_name = get_sheet_name(result_type)
    
    # Ensure sheet exists
    if sheet_name not in input_df:
        input_df[sheet_name] = pd.DataFrame(columns=['Prefix'])
    
    # Ensure sheet exists
    sheet = input_df[sheet_name].copy()
    updated_sheet = update_sheet_with_results(sheet, prefix, result_dict, detected_tetrad,output_extent)
    input_df[sheet_name] = updated_sheet


