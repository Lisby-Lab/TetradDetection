import glob
import os
import pickle
import numpy as np
from collections import Counter
import pandas as pd


def output_results(input_df, output_csv, root_path, a_prefix, detected_tetrad, result_dict, prefix, result_type):
    """
    Function to output results of a processing step into a DataFrame.
    
    Parameters:
    - input_df (pd.DataFrame): Input DataFrame to which the results should be appended.
    - output_csv (str): Path to save the CSV file.
    - root_path (str): Root directory path.
    - a_prefix (str): Prefix to identify the data.
    - detected_tetrad (dict): Dictionary containing detected tetrad results.
    - result_dict (dict): Dictionary containing results to be appended.
    - prefix (str): Prefix string for the result.
    - result_type (str): Type of result ("sum", "tetrad", "triad").
    
    Returns:
    - None
    """
    # Determine the sheet name based on the result type
    if result_type == "sum":
        sheet_name = "sum"
    elif result_type == "tetrad":
        sheet_name = "tetrad"
    else: 
        sheet_name = "triad"
   
    # If sheet doesn't exist, create it with given prefix
    if sheet_name not in input_df:
        input_df[sheet_name] = pd.DataFrame(columns=[...])  
        input_df[sheet_name].loc[0, 'Prefix'] = prefix
        
    sheet = input_df[sheet_name]
    # Find the row with the given prefix or append a new one
    if 'Prefix' in sheet.columns:
        hit = sheet.Prefix.astype(str).str.match('^' + prefix + '$')
        if hit.sum() == 0:
            row_index = len(sheet)
            sheet.loc[row_index, 'Prefix'] = prefix
        else:
            row_index = sheet.index[hit][0]
    else:
        row_index = len(sheet)
        sheet.loc[row_index, 'Prefix'] = prefix

    # Add results from the result_dict to the sheet
    new_data = {}
    for col, value in result_dict.items():
        #if col not in sheet.columns:
        if col not in new_data:
            new_data[col] = [np.nan] * len(sheet)
        new_data[col][row_index] = value

    for col in detected_tetrad:
        if col not in new_data:
            new_data[col] = [np.nan] * len(sheet)
        new_data[col][row_index] = detected_tetrad[col]

    # Add the new columns to the sheet
    new_df = pd.DataFrame(new_data)
    sheet = pd.concat([sheet, new_df], axis=1)

    input_df[sheet_name] = sheet

