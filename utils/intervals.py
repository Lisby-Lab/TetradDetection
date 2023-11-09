import numpy as np
import math
from scipy.stats import chisquare
from scipy import stats
from typing import Dict, Union, List, Any
import tkinter as tk
from tkinter import ttk

class TetradCalculatorApp:
    """
    A GUI application for the Tetrad Calculator.
    
    This application allows users to input color counter, tetrad types, and the total number of tetrads
    to calculate the results and display them within the GUI.
    
    Attributes:
        color_counter_entry (tk.Entry): Entry widget for inputting the color counter dictionary.
        tetrad_type_entry (tk.Entry): Entry widget for inputting the tetrad types dictionary.
        total_tetrad_entry (tk.Entry): Entry widget for inputting the total number of tetrads.
        compute_button (tk.Button): Button widget to initiate the computation.
        output_text (tk.Text): Text widget to display the computed results.
    """
    def __init__(self, root):
        """
        Initialize the TetradCalculatorApp with the provided root window.
        
        Args:
            root (tk.Tk): The main application window.
        """
        # Set title and size
        root.title("Tetrad Calculator")
        root.geometry("600x400")
        
        # Create and place labels and input fields for color_counter
        ttk.Label(root, text="Enter Color Counter (e.g. {'red': 5, 'blue': 3}):").grid(row=0, column=0, padx=10, pady=5)
        self.color_counter_entry = ttk.Entry(root, width=50)
        self.color_counter_entry.grid(row=0, column=1, padx=10, pady=5)

        # Create and place labels and input fields for tetrad_type_dict
        ttk.Label(root, text="Enter Tetrad Types (e.g. {'A': 5, 'B': 3}):").grid(row=1, column=0, padx=10, pady=5)
        self.tetrad_type_entry = ttk.Entry(root, width=50)
        self.tetrad_type_entry.grid(row=1, column=1, padx=10, pady=5)

        # Create and place labels and input fields for total_tetrads
        ttk.Label(root, text="Enter Total Tetrad Number (e.g. 8):").grid(row=2, column=0, padx=10, pady=5)
        self.total_tetrad_entry = ttk.Entry(root, width=50)
        self.total_tetrad_entry.grid(row=2, column=1, padx=10, pady=5)

        # Create and place the Compute and Exit buttons
        self.compute_button = ttk.Button(root, text="Compute", command=self.compute)
        self.compute_button.grid(row=3, column=0, padx=10, pady=20)
        ttk.Button(root, text="Exit", command=root.quit).grid(row=3, column=1, padx=10, pady=20)

        # Create and place a text box for output display
        self.output_text = tk.Text(root, height=10, width=60)
        self.output_text.grid(row=4, column=0, columnspan=2, padx=10, pady=10)
        
    def compute(self) -> None:
        """
        Retrieve the input data, compute the results using TetradCalculator,
        and display the results in the output text box.
        """
        # Retrieve input data
        color_counter = eval(self.color_counter_entry.get())
        tetrad_type_dict = eval(self.tetrad_type_entry.get())
        total_tetrads = eval(self.total_tetrad_entry.get())

        # Create an instance of TetradCalculator using the input data and compute the results
        calculator = TetradCalculator(color_counter, tetrad_type_dict, total_tetrads) 
        results = calculator.compute()
        
        # Display the results in the output text box
        self.output_text.insert(tk.END, str(results))

class TetradCalculator:
    CONSTANTS = {
        'wt_X': {'rb': 0.1748, 'ry': 0.1288, 'yb': 0.049}, # Map Distance (X = Morgan) calculated on large scale (>1000) Wild-Type dataset
        'se_wt': {'rb': 0.0046, 'ry': 0.004025, 'yb': 0.0023}, # Standard Error calculated according to Perkins equation on large scale Wild-Type dataset
        'wt_perkins_var': {'rb': 0.0000243, 'ry': 0.0000185, 'yb': 0.00000125}, # Sampling Variance of the Perkins Equation
         'non_disjunction_frequency_' : {'MI':0.05, 'MII': 0.01, 'events': 16063}, # tetrads from 4 independent datasets calculated for mean
        'gene_conversion_' : {'gain_wt_freq': 0.1472, 'loss_wt_freq':0.52, 'events': 16063}
    }
    """
    A class to perform calculations on tetrads based on provided frequencies.
    
    Attributes:
        result_dict (Dict[str, Union[str, float]]): Dictionary containing the calculated results.
        color_counter (Dict[str, int]): Counter of colors.
        tetrad_type_dict (Dict[str, int]): Dictionary of tetrad types and counts.
        total_tetrads (int): Total number of tetrads.
        n_total_all_tetrads (List[str]): Keys excluded for total tetrads computation.
        n_total_excluded_tetrads (List[str]): Keys excluded for genetic distance and interference calculations.
    """
    def __init__(self, color_counter: Dict[str, int],
                tetrad_type_dict: Dict[str, int],
                 total_tetrads: int) -> Dict[str, Union[str, float]]:            
        """
        Initialize the TetradCalculator with provided tetrad_type dictionary.
        
        Args:
            color_counter (Dict[str, int]): Counter of colors.
            tetrad_type_dict (Dict[str, int]): Dictionary containing tetrad types and their respective counts.
            total_tetrads (int): Total number of tetrads.
        """
            
        keys_str = """
        cM_rb cM_ry cM_yb significance_map_distance_rb significance_map_distance_ry significance_map_distance_yb
        X_rb X_ry X_yb significant_overlap_rb significant_overlap_ry significant_overlap_yb
        X2_rb X2_p_rb F_rb F_ry F_yb X2_ry X2_yb X2_p_ry X2_p_yb
        total_classified total_predicted 
        fP_rb fNPD_rb fT_rb
        var_fT_rb var_fNPD_rb cov_fT_fNPD_rb 
        SE_rb SE_yb SE_ry
        R_rb R_ry R_yb 
        E_rb E_ry E_yb
        fDCO_exp fDCO_obs DCO_exp DCO_obs
        fT_exp_rb fNPD_exp_rb fP_exp_rb
        MI_frequency MII_frequency MI_p_value MII_p_value
        I1 I2 X2_p_I1 X2_I1
        fP_ry fNPD_ry fT_ry
        var_fT_ry var_fNPD_ry cov_fT_fNPD_ry
        fP_yb fNPD_yb fT_yb
        var_fT_yb var_fNPD_yb cov_fT_fNPD_yb
        P_obs_rb P_exp_rb T_obs_rb T_exp_rb NPD_obs_rb NPD_exp_rb
        P_obs_yb P_exp_yb T_obs_yb T_exp_yb NPD_obs_yb NPD_exp_yb
        P_obs_ry P_exp_ry T_obs_ry T_exp_ry NPD_obs_ry NPD_exp_ry
        se_wt_mutant_rb se_wt_mutant_ry se_wt_mutant_yb 
        X_rb_perk X_ry_perk X_yb_perk
        perkins_var_rb perkins_var_ry perkins_var_yb 
        one_gain_sum one_loss_sum
        two_gain_sum two_loss_sum
        sum_gc_gain sum_gc_loss gc_gain_p_value gc_loss_p_value
        X_rb_interfering X_ry_interfering X_yb_interfering X_rb_non_interfereing X_ry_non_interfereing X_yb_non_interfereing
        non_interfering_coc_rb non_interfering_coc_ry non_interfering_coc_yb
        interfering_co_percentage_rb interfering_co_percentage_ry interfering_co_percentage_yb
        blue_fraction red_fraction yellow_fraction orange_fraction purple_fraction green_fraction white_fraction empty_fraction
        """
        keys = keys_str.split()
        self.result_dict = {key: 'Z' for key in keys}     
        self.color_counter = color_counter
        self.tetrad_type_dict = tetrad_type_dict
        self.total_tetrads = total_tetrads
        
        # excluded keys for n_total_all_tetrads and n_total computations
        self.n_total_all_tetrads = ['broken_mask', 'bbox_error', 'predictions_below_confidence', 
                                    'less_than_4_spores', 'mask_error', 'less_than_3_spores'] 
        
        self.n_total_excluded_tetrads = self.n_total_all_tetrads + ['Z', 'L', 'M', 'N1G', 'N2G', 'N3G', 'N1L', 'N2L', 'N3L', 'O1G', 'O2G', 'O3G', 'O4GL', 'O5G',  'O6G', 'O1L', 'O2L', 'O3L', 'O4L', 'O5L', 'O6L' ]  # Number of tetrads used for calculating genetic distance and interference - removed gene conversion, nondisjunction and quality control events
        self.n_total_tetrads = self.calculate_total(tetrad_type_dict, self.n_total_all_tetrads)
        self.n_total_viable_tetrads = self.calculate_total(tetrad_type_dict, self.n_total_excluded_tetrads)
        
        self.result_dict['total_predicted'] = self.total_tetrads
        self.result_dict['total_classified'] = self.n_total_tetrads
        
           
    def compute(self) -> Dict[str, Union[str, float]]:
        """
        Calculate and update the results using various calculators.
        
        Returns:
            Dict[str, Union[str, float]]: Dictionary of calculated results.
        """
        
        try:
            self.distance_calculator = ComputeGeneticDistance(self.tetrad_type_dict, self.n_total_viable_tetrads)
            self.result_dict.update(self.distance_calculator.calculate())

            self.interference_calculator = ComputeGeneticInterference(self.tetrad_type_dict, self.result_dict,self.n_total_viable_tetrads) 
            self.result_dict.update(self.interference_calculator.calculate())

            self.color_fraction_calculator = ComputeColorFraction(self.color_counter, self.result_dict)
            self.result_dict.update(self.color_fraction_calculator.calculate())

            self.non_disjunction_calculator = ComputeNonDisjunction(self.tetrad_type_dict,self.result_dict,self.n_total_tetrads)
            self.result_dict.update(self.non_disjunction_calculator.calculate())  

            self.gene_conversion_calculator = ComputeGeneConversion(self.tetrad_type_dict,self.result_dict,self.n_total_tetrads)
            self.result_dict.update(self.gene_conversion_calculator.calculate())
        except ZeroDivisionError:
            # You can print an error message or log it if needed
            print("Error: Division by zero encountered during calculations. Returning unmodified results.")
        except IndexError:
            print("Error: List index out of range.")
        except (TypeError, ValueError):  # Catching multiple exceptions in one except clause
            print("Error: Invalid type or value.")
        except Exception as e:  # This will catch any type of exception
            print(f"An unexpected error occurred: {e}")
            
        return self.result_dict


    def calculate_total(self, tetrad_type_dict: Dict[str, int], excluded_keys: list[str]) -> float:
        """
        Calculate the total count of tetrads excluding certain keys.
        
        Args:
            tetrad_type_dict (dict): Dictionary containing tetrad types and their respective counts.
            excluded_keys (list): List of keys to be excluded from the total count.
        
        Returns:
            float: Total count of tetrads after excluding the given keys.
        """
        total = sum(val for key, val in tetrad_type_dict.items() if key not in excluded_keys)
        return float(total)
        

class ComputeGeneticDistance:
    """
    A class to compute genetic distances based on provided tetrad frequencies.

    Attributes:
        tetrad_type_dict (Dict[str, int]): Dictionary of tetrad types and their respective counts.
        n_total_viable_tetrads (int): Total number of viable tetrads.
        result_dict (Dict[str, Union[float, str]]): Dictionary containing the computed results.
    """
    def __init__(self, tetrad_type_dict: Dict[str, int], n_total_viable_tetrads: int):
        """
        Initialize the ComputeGeneticDistance with provided tetrad_type dictionary and total number of viable tetrads.
        
        Args:
            tetrad_type_dict (Dict[str, int]): Dictionary of tetrad types and their respective counts.
            n_total_viable_tetrads (int): Total number of viable tetrads.
        """
        self.tetrad_type_dict = tetrad_type_dict 
        self.n_total_viable_tetrads = n_total_viable_tetrads
        self.result_dict = {}

    def calculate(self):
        """
        Calculate and return the genetic distances based on the tetrad frequencies.

        Returns:
            Dict[str, Union[float, str]]: Dictionary of computed results.
        """
        self.tetrad_frequency()
        self.perkins_distance_cM()
        self.perkins_distance_assumption_no_interference() 
        self.compute_significant_overlap()
        self.compute_significance_map_distance()
        
        return self.result_dict

    def tetrad_frequency(self) -> None:
        """
        Compute the genetic distances based on the tetrad type dictionary and total count.

        Args:
            tetrad_type_dict (dict): Dictionary containing tetrad types and their respective counts.
            total (int): Total count of tetrads.
            
        Compute the sampling variance of the Perkins equation for each color combination.

        Args:
            perkins_var_rb (float): Perkins variance for Red-Blue combination.
            perkins_var_ry (float): Perkins variance for Red-Yellow combination.
            perkins_var_yb (float): Perkins variance for Yellow-Blue combination.
            var_fT_rb, var_fNPD_rb, cov_fT_fNPD_rb: Variance and covariance for Red-Blue combination.
            var_fT_ry, var_fNPD_ry, cov_fT_fNPD_ry: Variance and covariance for Red-Yellow combination.
            var_fT_yb, var_fNPD_yb, cov_fT_fNPD_yb: Variance and covariance for Yellow-Blue combination.
        
        """
        # Compute fP, fNPD, and fT for each combination (red-blue, red-yellow, yellow-blue)
        self.result_dict['fP_rb'] = self._compute_fP(self.tetrad_type_dict, ['A', 'D', 'K'], self.n_total_viable_tetrads)
        self.result_dict['fNPD_rb'] = self._compute_fNPD(self.tetrad_type_dict, ['G', 'J', 'H'], self.n_total_viable_tetrads)
        self.result_dict['fT_rb'] = self._compute_fT(self.tetrad_type_dict, ['B', 'C', 'E', 'F', 'I'], self.n_total_viable_tetrads)

        self.result_dict['fP_ry'] = self._compute_fP(self.tetrad_type_dict, ['A', 'C'], self.n_total_viable_tetrads)
        self.result_dict['fNPD_ry']  = self._compute_fNPD(self.tetrad_type_dict, ['K', 'J', 'I', 'H'], self.n_total_viable_tetrads)
        self.result_dict['fT_ry'] = self._compute_fT(self.tetrad_type_dict, ['B', 'D', 'E', 'F', 'G'], self.n_total_viable_tetrads)

        self.result_dict['fP_yb'] = self._compute_fP(self.tetrad_type_dict, ['A', 'B', 'H', 'J'], self.n_total_viable_tetrads)
        self.result_dict['fNPD_yb'] = self._compute_fNPD(self.tetrad_type_dict, ['K'], self.n_total_viable_tetrads)
        self.result_dict['fT_yb'] = self._compute_fT(self.tetrad_type_dict, ['C', 'D', 'E', 'F', 'G', 'I'], self.n_total_viable_tetrads)

        self.result_dict['var_fT_rb'] = self.result_dict['fT_rb'] * (1 - self.result_dict['fT_rb']) / self.n_total_viable_tetrads
        self.result_dict['var_fNPD_rb'] = self.result_dict['fNPD_rb'] * (1 - self.result_dict['fNPD_rb']) / self.n_total_viable_tetrads
        self.result_dict['cov_fT_fNPD_rb'] = -(self.result_dict['fT_rb']) * (self.result_dict['fNPD_rb']) / self.n_total_viable_tetrads

        # Red - Yellow
        self.result_dict['var_fT_ry'] = self.result_dict['fT_ry'] * (1 - self.result_dict['fT_ry']) / self.n_total_viable_tetrads
        self.result_dict['var_fNPD_ry'] = self.result_dict['fNPD_ry'] * (1 - self.result_dict['fNPD_ry']) / self.n_total_viable_tetrads
        self.result_dict['cov_fT_fNPD_ry'] = -(self.result_dict['fT_ry']) * (self.result_dict['fNPD_ry']) / self.n_total_viable_tetrads

        # Yellow - Blue
        self.result_dict['var_fT_yb'] = self.result_dict['fT_yb'] * (1 - self.result_dict['fT_yb']) / self.n_total_viable_tetrads
        self.result_dict['var_fNPD_yb'] = self.result_dict['fNPD_yb'] * (1 - self.result_dict['fNPD_yb']) / self.n_total_viable_tetrads
        self.result_dict['cov_fT_fNPD_yb'] = -(self.result_dict['fT_yb']) * (self.result_dict['fNPD_yb']) / self.n_total_viable_tetrads
        
        self.result_dict['perkins_var_rb'] = 0.25 * self.result_dict['var_fT_rb'] + 9 * self.result_dict['var_fNPD_rb'] + 3 * self.result_dict['cov_fT_fNPD_rb']
        self.result_dict['perkins_var_ry'] = 0.25 * self.result_dict['var_fT_ry'] + 9 * self.result_dict['var_fNPD_ry'] + 3 * self.result_dict['cov_fT_fNPD_ry']
        self.result_dict['perkins_var_yb'] = 0.25 * self.result_dict['var_fT_yb'] + 9 * self.result_dict['var_fNPD_yb'] + 3 * self.result_dict['cov_fT_fNPD_yb'] 

        self.result_dict['SE_rb'] = math.sqrt(self.result_dict['perkins_var_rb'])
        self.result_dict['SE_ry'] = math.sqrt(self.result_dict['perkins_var_ry'])
        self.result_dict['SE_yb'] = math.sqrt(self.result_dict['perkins_var_yb'])  
        

    def _compute_fP(self, tetrad_type_dict: Dict[str, int], keys: List[str], total: int) -> float:
        """Helper method to compute fP frequency."""
        return sum(tetrad_type_dict.get(key, 0) for key in keys) / total
    
    def _compute_fNPD(self, tetrad_type_dict: Dict[str, int], keys: List[str], total: int) -> float:
        """Helper method to compute fNPD frequency."""
        return sum(tetrad_type_dict.get(key, 0) for key in keys) / total
    
    def _compute_fT(self, tetrad_type_dict: Dict[str, int], keys: List[str], total: int) -> float:
        """Helper method to compute fT frequency."""
        return sum(tetrad_type_dict.get(key, 0) for key in keys) / total


    def perkins_distance_cM(self) -> None:
        """
        Compute the map distance (cM) using the Perkins equation for each color combination.

        Args:
            tetrad_type_dict (dict): Dictionary containing tetrad types and their respective counts.
        """
        # Red - Blue
        self.result_dict['cM_rb'] = (100 * ((6 * (self.tetrad_type_dict['G'] + self.tetrad_type_dict['J'] + self.tetrad_type_dict['H'] ))
                                            + (self.tetrad_type_dict['B']  + self.tetrad_type_dict['C'] + self.tetrad_type_dict['E'] 
                                               + self.tetrad_type_dict['F'] + self.tetrad_type_dict['I']))) / (2 * self.n_total_viable_tetrads)

        # Red - Yellow
        self.result_dict['cM_ry'] = (100 * ((6 * (self.tetrad_type_dict['K'] + self.tetrad_type_dict['J'] + self.tetrad_type_dict['I']  + self.tetrad_type_dict['H'] ))
                                            + (self.tetrad_type_dict['B'] + self.tetrad_type_dict['D'] + self.tetrad_type_dict['E'] 
                                               + self.tetrad_type_dict['F'] + self.tetrad_type_dict['G']))) / (2 * self.n_total_viable_tetrads)

        # Yellow - Blue
        self.result_dict['cM_yb'] = (100 * ((6 * (self.tetrad_type_dict['K'] ))
                                            + (self.tetrad_type_dict['C'] + self.tetrad_type_dict['D'] + self.tetrad_type_dict['E'] 
                                               + self.tetrad_type_dict['F'] + self.tetrad_type_dict['G'] + self.tetrad_type_dict['I']))) / (2 * self.n_total_viable_tetrads)

    def perkins_distance_assumption_no_interference(self) -> None:
        """Compute X using the Perkins equation without interference assumption."""
        # Computing X using the Perkins equation
        self.result_dict['X_rb_perk'] = (self.result_dict['fT_rb']/2) + (3*self.result_dict['fNPD_rb'])
        self.result_dict['X_ry_perk'] = (self.result_dict['fT_ry']/2) + (3*self.result_dict['fNPD_ry'])
        self.result_dict['X_yb_perk'] = (self.result_dict['fT_yb']/2) + (3*self.result_dict['fNPD_yb'])


    def compute_significant_overlap(self) -> None:
        """Compute significant overlap using method 1."""
        wt_X_rb, se_wt_rb = TetradCalculator.CONSTANTS['wt_X']['rb'], TetradCalculator.CONSTANTS['se_wt']['rb']
        wt_X_ry, se_wt_ry = TetradCalculator.CONSTANTS['wt_X']['ry'], TetradCalculator.CONSTANTS['se_wt']['ry']
        wt_X_yb, se_wt_yb = TetradCalculator.CONSTANTS['wt_X']['yb'], TetradCalculator.CONSTANTS['se_wt']['yb']

        self.result_dict['significant_overlap_rb'] = not (wt_X_rb + se_wt_rb < self.result_dict['X_rb_perk'] - self.result_dict['SE_rb'] or wt_X_rb - se_wt_rb > self.result_dict['X_rb_perk'] + self.result_dict['SE_rb'])
        self.result_dict['significant_overlap_ry'] = not (wt_X_ry + se_wt_ry < self.result_dict['X_ry_perk'] - self.result_dict['SE_ry'] or wt_X_ry - se_wt_ry > self.result_dict['X_ry_perk'] + self.result_dict['SE_ry'])
        self.result_dict['significant_overlap_yb'] = not (wt_X_yb + se_wt_yb < self.result_dict['X_yb_perk'] - self.result_dict['SE_yb'] or wt_X_yb - se_wt_yb > self.result_dict['X_yb_perk'] + self.result_dict['SE_yb'])       

    def compute_significance_map_distance(self) -> None:
        """Compute significance of map distance using method 2."""
        wt_perkins_var_rb = TetradCalculator.CONSTANTS['wt_perkins_var']['rb']
        wt_perkins_var_ry = TetradCalculator.CONSTANTS['wt_perkins_var']['ry']
        wt_perkins_var_yb = TetradCalculator.CONSTANTS['wt_perkins_var']['yb']

        var_wt_mutant_rb = self.result_dict['perkins_var_rb'] + wt_perkins_var_rb
        se_wt_mutant_rb = math.sqrt(var_wt_mutant_rb)

        var_wt_mutant_ry = self.result_dict['perkins_var_ry'] + wt_perkins_var_ry
        se_wt_mutant_ry = math.sqrt(var_wt_mutant_ry)

        var_wt_mutant_yb = self.result_dict['perkins_var_yb'] + wt_perkins_var_yb
        se_wt_mutant_yb = math.sqrt(var_wt_mutant_yb)

        self.result_dict['significance_map_distance_rb'] = 'yes' if (2*se_wt_mutant_rb) < abs(TetradCalculator.CONSTANTS['wt_X']['rb'] - self.result_dict['X_rb_perk']) else 'no'
        self.result_dict['significance_map_distance_ry'] = 'yes' if (2*se_wt_mutant_ry) < abs(TetradCalculator.CONSTANTS['wt_X']['ry'] - self.result_dict['X_ry_perk']) else 'no'
        self.result_dict['significance_map_distance_yb'] = 'yes' if (2*se_wt_mutant_yb) < abs(TetradCalculator.CONSTANTS['wt_X']['yb'] - self.result_dict['X_yb_perk']) else 'no'


class ComputeGeneticInterference:
    def __init__(self, 
                 tetrad_type_dict: Dict[str, int], 
                 result_dict: Dict[str, Any], 
                 n_total_viable_tetrads: int) -> None:
        """
        Initialize the ComputeGeneticInterference class.

        Args:
            tetrad_type_dict (Dict[str, int]): Dictionary containing tetrad types and their respective counts.
            result_dict (Dict[str, Any]): Dictionary containing result values.
            n_total_viable_tetrads (int): Total count of viable tetrads.
        """
        self.tetrad_type_dict = tetrad_type_dict
        self.result_dict = result_dict
        self.n_total_viable_tetrads = n_total_viable_tetrads

    def calculate(self) -> Dict[str, Any]:  
        """
        Perform calculations to compute genetic interference metrics.

        Returns:
            Dict[str, Any]: Dictionary containing the results.
        """
        self.compute_recombination_frequency()
        self.compute_map_length_no_interference()
        self.compute_interference_papazian()
        self.compute_interference_stahl()
        self.compute_map_length_no_interference()
        self.compute_interference_by_CoC()
        self.compute_interference_papazian_mod()
        self.compute_non_interfering_map_stahl()
        
        return self.result_dict


    def compute_interference_papazian(self) -> None:
        """Compute interference using Papazian's method."""
        self.result_dict['E_rb'] = 0.5 * ((1 - self.result_dict['fT_rb']) - (1 - (3 * self.result_dict['fT_rb'] / 2))**(2/3))
        self.result_dict['E_ry'] = 0.5 * ((1 - self.result_dict['fT_ry']) - (1 - (3 * self.result_dict['fT_ry'] / 2))**(2/3))
        self.result_dict['E_yb'] = 0.5 * ((1 - self.result_dict['fT_yb']) - (1 - (3 * self.result_dict['fT_yb'] / 2))**(2/3))

        
    def compute_interference_stahl(self) -> None:
        """Compute interference using Franklin Stahl (2008) method."""
        self.result_dict['P_obs_rb'] = float(self.tetrad_type_dict['A'] + self.tetrad_type_dict['D'] + self.tetrad_type_dict['K'])
        self.result_dict['T_obs_rb'] = float(self.tetrad_type_dict['B'] + self.tetrad_type_dict['C'] + self.tetrad_type_dict['E'] + self.tetrad_type_dict['F'] + self.tetrad_type_dict['I'])
        self.result_dict['NPD_obs_rb'] = float(self.tetrad_type_dict['G'] + self.tetrad_type_dict['J'] + self.tetrad_type_dict['H'])

        self.result_dict['P_obs_ry'] = float(self.tetrad_type_dict['A'] + self.tetrad_type_dict['C'])
        self.result_dict['T_obs_ry'] = float(self.tetrad_type_dict['B'] + self.tetrad_type_dict['D'] + self.tetrad_type_dict['E'] + self.tetrad_type_dict['F'] + self.tetrad_type_dict['G'])
        self.result_dict['NPD_obs_ry'] = float(self.tetrad_type_dict['J'] + self.tetrad_type_dict['K'] + self.tetrad_type_dict['I'] + self.tetrad_type_dict['H'])

        self.result_dict['P_obs_yb'] = float(self.tetrad_type_dict['A'] + self.tetrad_type_dict['B'] + self.tetrad_type_dict['H'] + self.tetrad_type_dict['J'])
        self.result_dict['T_obs_yb'] = float(self.tetrad_type_dict['C'] + self.tetrad_type_dict['D'] + self.tetrad_type_dict['E'] + self.tetrad_type_dict['F'] + self.tetrad_type_dict['G'] + self.tetrad_type_dict['I'])
        self.result_dict['NPD_obs_yb'] = float(self.tetrad_type_dict['K'])

        # Frequency of TTTetrad types expected
        self.result_dict['fT_exp_rb'] = 2/3 * (1 - math.exp(-3 * self.result_dict['X_rb']))
        self.result_dict['fT_exp_ry'] = 2/3 * (1 - math.exp(-3 * self.result_dict['X_ry']))
        self.result_dict['fT_exp_yb'] = 2/3 * (1 - math.exp(-3 * self.result_dict['X_yb']))

        # Frequency of NPDs expected
        self.result_dict['fNPD_exp_rb'] = 0.5 * ((1 - math.exp(-2 * self.result_dict['X_rb'])) - self.result_dict['fT_exp_rb'])
        self.result_dict['fNPD_exp_ry'] = 0.5 * ((1 - math.exp(-2 * self.result_dict['X_ry'])) - self.result_dict['fT_exp_ry'])
        self.result_dict['fNPD_exp_yb'] = 0.5 * ((1 - math.exp(-2 * self.result_dict['X_yb'])) - self.result_dict['fT_exp_yb'])

        # Frequency of Parental type expected
        self.result_dict['fP_exp_rb'] = 1 - self.result_dict['fT_exp_rb'] - self.result_dict['fNPD_exp_rb']
        self.result_dict['fP_exp_ry'] = 1 - self.result_dict['fT_exp_ry'] - self.result_dict['fNPD_exp_ry']
        self.result_dict['fP_exp_yb'] = 1 - self.result_dict['fT_exp_yb'] - self.result_dict['fNPD_exp_yb']

        # Total expected classes
        self.result_dict['P_exp_rb'] = self.result_dict['fP_exp_rb'] * float(self.n_total_viable_tetrads)
        self.result_dict['T_exp_rb'] = self.result_dict['fT_exp_rb'] * float(self.n_total_viable_tetrads)
        self.result_dict['NPD_exp_rb'] = self.result_dict['fNPD_exp_rb'] * float(self.n_total_viable_tetrads)

        self.result_dict['P_exp_ry'] = self.result_dict['fP_exp_ry'] * float(self.n_total_viable_tetrads)
        self.result_dict['T_exp_ry'] = self.result_dict['fT_exp_ry'] * float(self.n_total_viable_tetrads)
        self.result_dict['NPD_exp_ry'] = self.result_dict['fNPD_exp_ry'] * float(self.n_total_viable_tetrads)

        self.result_dict['P_exp_yb'] = self.result_dict['fP_exp_yb'] * float(self.n_total_viable_tetrads)
        self.result_dict['T_exp_yb'] = self.result_dict['fT_exp_yb'] * float(self.n_total_viable_tetrads)
        self.result_dict['NPD_exp_yb'] = self.result_dict['fNPD_exp_yb'] * float(self.n_total_viable_tetrads)

        # Chi-Square with one degree of freedom
        self.result_dict['X2_rb'] = (self.result_dict['P_obs_rb'] - self.result_dict['P_exp_rb'])**2 / self.result_dict['P_exp_rb'] + \
                                    (self.result_dict['T_obs_rb'] - self.result_dict['T_exp_rb'])**2 / self.result_dict['T_exp_rb'] + \
                                    (self.result_dict['NPD_obs_rb'] - self.result_dict['NPD_exp_rb'])**2 / self.result_dict['NPD_exp_rb']
                
        self.result_dict['X2_p_rb'] = 1 - stats.chi2.cdf(self.result_dict['X2_rb'], 1)

        self.result_dict['X2_ry'] = (self.result_dict['P_obs_ry'] - self.result_dict['P_exp_ry'])**2 / self.result_dict['P_exp_ry'] + \
                                    (self.result_dict['T_obs_ry'] - self.result_dict['T_exp_ry'])**2 / self.result_dict['T_exp_ry'] + \
                                    (self.result_dict['NPD_obs_ry'] - self.result_dict['NPD_exp_ry'])**2 / self.result_dict['NPD_exp_ry']
        self.result_dict['X2_p_ry'] = 1 - stats.chi2.cdf(self.result_dict['X2_ry'], 1)

        self.result_dict['X2_yb'] = (self.result_dict['P_obs_yb'] - self.result_dict['P_exp_yb'])**2 / self.result_dict['P_exp_yb'] + \
                                    (self.result_dict['T_obs_yb'] - self.result_dict['T_exp_yb'])**2 / self.result_dict['T_exp_yb'] + \
                                    (self.result_dict['NPD_obs_yb'] - self.result_dict['NPD_exp_yb'])**2 / self.result_dict['NPD_exp_yb']
        self.result_dict['X2_p_yb'] = 1 - stats.chi2.cdf(self.result_dict['X2_yb'], 1)

    def compute_recombination_frequency(self) -> None:
        """Compute recombination frequency for the three color combinations."""
        self.result_dict['R_rb'] = float(self.result_dict['fNPD_rb'] + (self.result_dict['fT_rb'] / 2))
        self.result_dict['R_ry'] = float(self.result_dict['fNPD_ry'] + (self.result_dict['fT_ry'] / 2))
        self.result_dict['R_yb'] = float(self.result_dict['fNPD_yb'] + (self.result_dict['fT_yb'] / 2))

    def compute_map_length_value(self, R_value):
        """Helper function to compute map length value given R_value."""
        value = 1 - 2 * R_value
        if value > 0:
            return -0.5 * np.log(value)
        else:
            return np.nan

    def compute_map_length_no_interference(self) -> None:
        """Compute map length under the assumption of no interference for the three color combinations."""
        self.result_dict['X_rb'] = self.compute_map_length_value(self.result_dict['R_rb'])
        self.result_dict['X_ry'] = self.compute_map_length_value(self.result_dict['R_ry'])
        self.result_dict['X_yb'] = self.compute_map_length_value(self.result_dict['R_yb'])

    def compute_interference_by_CoC(self) -> None:
        """Compute interference using CoC method."""
        # Expected DCO
        self.result_dict['fDCO_exp'] = 100 * (self.result_dict['X_ry'] * self.result_dict['X_yb'])
        self.result_dict['DCO_exp'] = (self.n_total_viable_tetrads / 100) * self.result_dict['fDCO_exp']

        # Frequency of DCOobs
        self.result_dict['DCO_obs'] = float(self.tetrad_type_dict['E'] + self.tetrad_type_dict['F'] + self.tetrad_type_dict['G'] + 
                                            self.tetrad_type_dict['D'] + self.tetrad_type_dict['K'] + self.tetrad_type_dict['I'] + 
                                            self.tetrad_type_dict['J'])
        self.result_dict['fDCO_obs'] = self.result_dict['DCO_obs'] / self.n_total_viable_tetrads

        # Compute I and its significance
        self.result_dict['I1'] = (self.result_dict['R_ry'] + self.result_dict['R_yb'] - self.result_dict['R_rb']) / (2 * self.result_dict['R_ry'] * self.result_dict['R_yb'])
        self.result_dict['X2_I1'] = (self.result_dict['DCO_obs'] - self.result_dict['DCO_exp'])**2 / self.result_dict['DCO_exp']
        self.result_dict['X2_p_I1'] = 1 - stats.chi2.cdf(self.result_dict['X2_I1'], 1)

        self.result_dict['I2'] = self.result_dict['fDCO_obs'] / (self.result_dict['R_ry'] * self.result_dict['R_yb'])

    def compute_interference_papazian_mod(self) -> None:
        """Compute interference using the modified Papazian NPD ratio."""
        self.result_dict['F_rb'] = self.result_dict['fNPD_rb'] / self.result_dict['fNPD_exp_rb']
        self.result_dict['F_ry'] = self.result_dict['fNPD_ry'] / self.result_dict['fNPD_exp_ry']
        self.result_dict['F_yb'] = self.result_dict['fNPD_yb'] / self.result_dict['fNPD_exp_yb']


    def compute_non_interfering_map_stahl(self) -> None:
        """Compute non-interfering map distance and CoC using Stahl's method."""
        self.result_dict['X_rb_interfering'] = self.result_dict['X_rb_perk'] * (1 - self.result_dict['F_rb'])**0.5
        self.result_dict['X_ry_interfering'] = self.result_dict['X_ry_perk'] * (1 - self.result_dict['F_ry'])**0.5
        self.result_dict['X_yb_interfering'] = self.result_dict['X_yb_perk'] * (1 - self.result_dict['F_yb'])**0.5

        self.result_dict['X_rb_non_interfereing'] = self.result_dict['X_rb_perk'] - self.result_dict['X_rb_interfering']
        self.result_dict['X_ry_non_interfereing'] = self.result_dict['X_ry_perk'] - self.result_dict['X_ry_interfering']
        self.result_dict['X_yb_non_interfereing'] = self.result_dict['X_yb_perk'] - self.result_dict['X_yb_interfering']

        self.result_dict['non_interfering_coc_rb'] = 1 - ((self.result_dict['X_rb_perk']**2) / (self.result_dict['X_rb_interfering']**2)) + ((self.result_dict['X_rb_perk']**2) * (1 - self.result_dict['F_rb'])) / (self.result_dict['X_rb_interfering']**2)

        self.result_dict['non_interfering_coc_ry'] = 1 - ((self.result_dict['X_ry_perk']**2) / (self.result_dict['X_ry_interfering']**2)) + ((self.result_dict['X_ry_perk']**2) * (1 - self.result_dict['F_ry'])) / (self.result_dict['X_ry_interfering']**2)

        self.result_dict['non_interfering_coc_yb'] = 1 - ((self.result_dict['X_yb_perk']**2) / (self.result_dict['X_yb_interfering']**2)) + ((self.result_dict['X_yb_perk']**2) * (1 - self.result_dict['F_yb'])) / (self.result_dict['X_yb_interfering']**2)

        self.result_dict['interfering_co_percentage_rb'] = (self.result_dict['X_rb_interfering'] / self.result_dict['X_rb_perk']) * 100
        self.result_dict['interfering_co_percentage_ry'] = (self.result_dict['X_ry_interfering'] / self.result_dict['X_ry_perk']) * 100
        self.result_dict['interfering_co_percentage_yb'] = (self.result_dict['X_yb_interfering'] / self.result_dict['X_yb_perk']) * 100
        if self.result_dict['F_rb'] > 1:
            self.result_dict['X_rb_interfering'], self.result_dict['X_ry_interfering'], self.result_dict['X_yb_interfering'], self.result_dict['X_rb_non_interfereing'], self.result_dict['X_ry_non_interfereing'], self.result_dict['X_yb_non_interfereing'], self.result_dict['non_interfering_coc_rb'], self.result_dict['non_interfering_coc_ry'], self.result_dict['non_interfering_coc_yb'], self.result_dict['interfering_co_percentage_rb'], self.result_dict['interfering_co_percentage_ry'], self.result_dict['interfering_co_percentage_yb'] = [0] * 12

class ComputeNonDisjunction:
    def __init__(self,tetrad_type_dict,result_dict,n_total_tetrads):
        self.tetrad_type_dict = tetrad_type_dict
        self.result_dict = result_dict
        self.n_total_tetrads = n_total_tetrads
        
    def calculate(self):
        self.compute_nondisjunction_events()
        return self.result_dict

    def compute_nondisjunction_events(self) -> None:
        # Extracting the values from the CONSTANTS dictionary
        classified_events_wt = TetradCalculator.CONSTANTS['non_disjunction_frequency_']['events']
        MI_freq_wt = TetradCalculator.CONSTANTS['non_disjunction_frequency_']['MI']
        MII_freq_wt = TetradCalculator.CONSTANTS['non_disjunction_frequency_']['MII']

        # MI Frequency
        MI_proportion = (self.tetrad_type_dict['L'] / self.n_total_tetrads)
        self.result_dict['MI_frequency'] = (MI_proportion * 100)

        # MII Frequency
        MII_proportion = (self.tetrad_type_dict['M'] / self.n_total_tetrads)
        self.result_dict['MII_frequency'] = (MII_proportion * 100)

        # Pooled proportions
        MI_pooled_proportion = ((MI_proportion * self.n_total_tetrads) + (MI_freq_wt * classified_events_wt)) / (classified_events_wt + self.n_total_tetrads)
        MII_pooled_proportion = ((MII_proportion * self.n_total_tetrads) + (MII_freq_wt * classified_events_wt)) / (classified_events_wt + self.n_total_tetrads)

        # Z-scores
        MI_z_score = (MI_proportion - MI_freq_wt) / math.sqrt(MI_pooled_proportion * (1 - MI_pooled_proportion) * (1 / classified_events_wt + 1 / self.n_total_tetrads))
        MII_z_score = (MII_proportion - MII_freq_wt) / math.sqrt(MII_pooled_proportion * (1 - MII_pooled_proportion) * (1 / classified_events_wt + 1 / self.n_total_tetrads))

        # P-values
        self.result_dict['MI_p_value'] = 2 * (1 - stats.norm.cdf(abs(MI_z_score)))
        self.result_dict['MII_p_value'] = 2 * (1 - stats.norm.cdf(abs(MII_z_score)))

class ComputeGeneConversion:  
    def __init__(self,
                 tetrad_type_dict: Dict[str, int],
                 result_dict: Dict[str, Any],
                 n_total_tetrads: int) -> None:
        """
        Initialize the ComputeGeneConversion class.

        Args:
            tetrad_type_dict (Dict[str, int]): Dictionary containing tetrad types and their respective counts.
            result_dict (Dict[str, Any]): Dictionary containing result values.
            n_total_tetrads (int): Total count of tetrads.
        """
        self.tetrad_type_dict = tetrad_type_dict
        self.result_dict = result_dict
        self.n_total_tetrads = n_total_tetrads
        
    def calculate(self) -> Dict[str, Any]:
        """
        Perform calculations to compute gene conversion events.

        Returns:
            Dict[str, Any]: Dictionary containing the results.
        """
        self.compute_gene_conversion_events()
        return self.result_dict
        
    def calculate_gain_loss(self, 
                            event_types: List[str], 
                            tetrad_type_dict: Dict[str, int], 
                            total_tetrads: int) -> Dict[str, float]:
        """
        Calculate gain and loss for given event types.

        Args:
            event_types (List[str]): List of event types to consider.
            tetrad_type_dict (Dict[str, int]): Dictionary containing tetrad types and their respective counts.
            total_tetrads (int): Total count of tetrads.

        Returns:
            Dict[str, float]: Dictionary containing gain and loss percentages for the given event types.
        """
        return {event: (tetrad_type_dict[event] / total_tetrads) * 100 for event in event_types}
    
    def compute_gene_conversion_events(self) -> None: 
        # 1 GC events gain and loss
        events_one_gain = ['N1G', 'N2G', 'N3G']
        events_one_loss = ['N1L', 'N2L', 'N3L']

        # 2 GC events gain and loss
        events_two_gain = ['O1G', 'O2G', 'O3G', 'O4G', 'O5G', 'O6G']
        events_two_loss = ['O1L', 'O2L', 'O3L', 'O4L', 'O5L', 'O6L']

        one_gain_results = self.calculate_gain_loss(events_one_gain, self.tetrad_type_dict, self.n_total_tetrads)
        one_loss_results = self.calculate_gain_loss(events_one_loss, self.tetrad_type_dict, self.n_total_tetrads)

        two_gain_results = self.calculate_gain_loss(events_two_gain, self.tetrad_type_dict, self.n_total_tetrads)
        two_loss_results = self.calculate_gain_loss(events_two_loss, self.tetrad_type_dict, self.n_total_tetrads)

        self.result_dict['one_gain_sum'] = sum(one_gain_results.values())
        self.result_dict['one_loss_sum'] = sum(one_loss_results.values())

        self.result_dict['two_gain_sum'] = sum(two_gain_results.values())
        self.result_dict['two_loss_sum'] = sum(two_loss_results.values())


        # Extracting total gain and loss frequencies
        sum_gc_gain_frequency = (sum(one_gain_results.values()) + sum(two_gain_results.values())) / 100
        sum_gc_loss_frequency = (sum(one_loss_results.values()) + sum(two_loss_results.values())) / 100
        self.result_dict['sum_gc_gain'] = sum_gc_gain_frequency * 100
        self.result_dict['sum_gc_loss'] = sum_gc_loss_frequency * 100

class ComputeColorFraction:
    def __init__(self, 
                 color_fraction: Dict[str, int], 
                 result_dict: Dict[str, Any]) -> None:
        """
        Initialize the ComputeColorFraction class.

        Args:
            color_fraction (Dict[str, int]): Dictionary containing counts of different colors.
            result_dict (Dict[str, Any]): Dictionary containing result values.
        """
        self.color_fraction = color_fraction
        self.result_dict = result_dict
        
    def calculate(self) -> Dict[str, Any]:
        self.calculate_color_fractions()
        return self.result_dict
    
    def calculate_color_fractions(self) -> None:
        """
        Calculate the distribution of colors within the dataset.
        
        Args:
            color_counter (dict): Dictionary containing counts of different colors.
        
        Returns:
            dict: Fractions of each color.
        """
        total_color_count = sum(self.color_fraction.values())
        color_count = {color: count/total_color_count for color, count in self.color_fraction.items()}

        colors = ['blue', 'red', 'yellow', 'orange', 'purple', 'green', 'white', 'empty']
        for color in colors:
            self.result_dict[f'{color}_fraction'] = color_count.get(color, 0)

        
# Run the GUI
#root = tk.Tk()
#app = TetradCalculatorApp(root)
#root.mainloop()        