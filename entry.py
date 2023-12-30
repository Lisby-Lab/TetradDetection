#!/usr/bin/env python3

from utils import predict
from utils import tetrad_types 
from utils import file_utils
import numpy as np
from collections import defaultdict, namedtuple 
from utils import intervals
import argparse
import glob
import pandas as pd
import os
from typing import List, Tuple, Dict, Counter, Any, DefaultDict
import logging
import torch
from torchvision.models.detection.transform import GeneralizedRCNNTransform
from utils import engine

class TetradCounter:
    """
    Represents the primary class for detecting and classifying tetrads.
    This class is responsible for loading models, processing datasets, processing images, and writing to Excel the results.
    """
    Counts = namedtuple("Counts", ["tetrads", "triads"])

    def __init__(self, root_path: str, output_path: str, prefix, excel_file, output_extent, 
                 tetrad_model_path: str = None, color_model_path: str = None) -> None:
        """
        Initializes the tetrad and color detection models and sets the root path.

        Parameters:
        - root_path (str): The root directory containing the datasets.
        - output_path (str): Directory to save the results.
        - prefix (str): Prefix name (e.g. slx4).
        - excel_file (str): Name of the Excel file.
        - output_extent (str): Extent of the output to Excel file ('limited' or 'extended').
        - tetrad_model_path (str): Path to the tetrad/triad model.
        - color_model_path (str): Path to the color model.
        """ 
        self.device = torch.device('cpu')
        self.root_path = root_path
        self.output_path = output_path
        self.excel_file = excel_file
        self.output_extent = output_extent
        self.prefix = prefix
        self.tetrad_model_path = tetrad_model_path
        self.color_model_path = color_model_path
        self.tetrads_well = 0
        self.triads_well = 0
        self._tetrad_model = None
        self._color_model = None
        
    @property
    def tetrad_model(self):
        """Lazy loading of tetrad model."""
        if self._tetrad_model is None:
            self._tetrad_model = self._load_tetrad_prediction_model()
        return self._tetrad_model

    def _load_tetrad_prediction_model(self) -> torch.nn.Module:
        """
        Loads the tetrad prediction model.

        Returns:
            torch.model: The loaded tetrad prediction model.
        """
        model = engine.TetradPredictionArchitecture(num_classes=3)
        model.load_state_dict(torch.load(self.tetrad_model_path, map_location=self.device))
        if not os.path.exists(self.tetrad_model_path):
            logging.error(f"Model not found at {self.tetrad_model_path}")
            raise FileNotFoundError(f"Model not found at {self.tetrad_model_path}")
        grcnn = GeneralizedRCNNTransform(min_size=1080, max_size=1080, 
                                         image_mean=[0.485, 0.456, 0.406], 
                                         image_std=[0.229, 0.224, 0.225])
        model.transform = grcnn
        model.to(self.device)
        model.eval()
        return model    
    
    @property
    def color_model(self):
        """Lazy loading of color model."""
        if self._color_model is None:
            self._color_model = self._load_color_prediction_model()
        return self._color_model

    def _load_color_prediction_model(self) -> torch.nn.Module:
        """
        Loads the color prediction model.

        Returns:
            torch.model: The loaded color prediction model.
        """
        
        model = engine.ColorClassificationArchitecture(num_classes=9)
        model.load_state_dict(torch.load(self.color_model_path, map_location=self.device))
        if not os.path.exists(self.color_model_path):
            logging.error(f"Model not found at {self.color_model_path}")
            raise FileNotFoundError(f"Model not found at {self.color_model_path}")
        model.to(self.device)
        model.eval()
        return model    
    

    def run_tetrad_count_all_process(self) -> None:
        """Main process for running tetrad detection and color classification on all datasets."""
        if not os.path.exists(self.root_path):
            logging.error(f"Error: The directory {self.root_path} does not exist.")
            return
        files = sorted(os.listdir(self.root_path))
         # Filter out non-image files and directories
        image_files = [f for f in files if f.endswith('.tiff') and not os.path.isdir(os.path.join(self.root_path, f))]
        # Check if there are any image files to process
        if not image_files:
            error_message = f"No valid image files (.tiff) found in {self.root_path}"
            logging.error(error_message)
            raise FileNotFoundError(error_message)

        well_prefixes = ["_" + f.split('/')[-1].split('_')[1] + "_" for f in image_files]
        well_prefixes = list(np.unique(well_prefixes))
        self._process_dataset(well_prefixes, self.root_path)

    def _process_dataset(self, dataset: str, root_path: str) -> None:
        """
        Processes a single dataset.

        Parameters:
        - dataset (str): The dataset to be processed.
        - root_path (str): The root directory containing the dataset.
        """
        logging.debug(f'[current_dataset_run]: {dataset}')  
        output_csv = os.path.join(self.output_path, self.excel_file)
        input_csv = os.path.join(self.output_path, self.excel_file)
      #  if not os.path.exists(input_csv):
      #      logging.error(f"Excel file not found at {input_csv}")
      #      return
        
        logging.debug(f'[current_output_excel]: {input_csv}')
        
       # if not os.path.exists(input_csv):
        TetradCounter.create_empty_excel_file(input_csv)
                
        input_df = pd.read_excel(input_csv, sheet_name=None)
        
        tetrad_type_dict = defaultdict(int)
        triad_type_dict = defaultdict(int)  
        sum_type_dict = defaultdict(int)
        color_dict = defaultdict(int)
        combination_count= defaultdict(int)
        for a_prefix in dataset:  
            imgs = sorted(glob.glob(os.path.join(root_path, f'*{a_prefix}*.tiff')))
            counts,detected_tetrads,detected_triads,detected_sum, color_dict, combination_count= self._process_image(imgs,dataset,tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count)
            n_verified_tetrads, n_verified_triad, n_verified_sum = [sum(d.values()) for d in (tetrad_type_dict, triad_type_dict, sum_type_dict)]

            self.tetrads_well += counts.tetrads
            self.triads_well += counts.triads
            
            logging.debug(f'[total_nmb_tetrads / image] {self.tetrads_well}')
            logging.debug(f'[detected_tetrad_types] {tetrad_type_dict}')

            logging.debug(f'[total_nmb_triads / image] {self.triads_well}')        
            logging.debug(f'[detected_triad_types] {triad_type_dict}')

            logging.debug(f'[detected_total_types] {sum_type_dict}')

            logging.debug(f'[detected_colors] {color_dict}')
            logging.debug(f'[combination_patterns] {combination_count}')
 
        tetrad_counter_sum = Counter(detected_tetrads)
        triad_counter_sum = Counter(detected_triads)
        counter_sum = Counter(detected_sum)
        color_counter = Counter(color_dict)
        
        if sum(counter_sum.values()) != 0: 
            # Create an instance of the TetradCalculator
            calculator = intervals.TetradCalculator(color_counter,counter_sum, (self.tetrads_well+self.triads_well))
            # Call the compute method on the calculator instance
            sum_result_dict = calculator.compute()   
            logging.debug(f'[total_classified] {sum_result_dict["total_classified"]}')
            file_utils.output_results(input_df, self.root_path, a_prefix, counter_sum, sum_result_dict,self.prefix,self.output_extent,'sum')
        if sum(tetrad_counter_sum.values()) != 0:
            # Create an instance of the TetradCalculator
            calculator = intervals.TetradCalculator(color_counter,tetrad_counter_sum, self.tetrads_well)
            # Call the compute method on the calculator instance
            tetrad_result_dict = calculator.compute()   
            logging.debug(f'[total_tetrads_classified] {tetrad_result_dict["total_classified"]}')
            file_utils.output_results(input_df, self.root_path, a_prefix, tetrad_counter_sum, tetrad_result_dict,self.prefix,self.output_extent, 'tetrad')
            
        if sum(triad_counter_sum.values()) != 0:
            # Create an instance of the TetradCalculator
            calculator = intervals.TetradCalculator(color_counter,triad_counter_sum, self.triads_well)
            # Call the compute method on the calculator instance
            triad_result_dict = calculator.compute() 
            logging.debug(f'[total_triads_classified] {triad_result_dict["total_classified"]}')
            file_utils.output_results(input_df, self.root_path, a_prefix, triad_counter_sum, triad_result_dict,self.prefix,self.output_extent, 'triad')

        # Write the updated DataFrame back to the Excel file once after all the loop iterations
        with pd.ExcelWriter(output_csv, engine='openpyxl') as writer:
            for name, sheet in input_df.items():
                sheet.to_excel(writer, sheet_name=name, index=False)


    def _process_image(self, image: str,dataset: str, tetrad_type_dict: defaultdict, triad_type_dict: defaultdict, sum_type_dict: defaultdict, color_dict: defaultdict, combination_count: defaultdict):
        """
        Processes a single image.

        Parameters:
        - image (str): The image to be processed.
        - dataset (str): The dataset to which the image belongs.
        - tetrad_type_dict (defaultdict): Dictionary to hold tetrad types and counts.
        - triad_type_dict (defaultdict): Dictionary to hold triad types and counts.
        - sum_type_dict (defaultdict): Dictionary to hold summed types and counts.
        - color_dict (defaultdict): Dictionary to hold color counts.
        - combination_count (defaultdict): Dictionary to hold combination counts.
        
        Returns:
        - Counts: Namedtuple containing counts of detected tetrads and triads.
        - tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count: Updated dictionaries.
        """
        if not all(os.path.exists(img) for img in image):
            logging.error(f"One or more images in {image} not found.")
            return TetradCounter.Counts(0, 0), tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count
        if len(image) < 3:
            logging.debug(f"Error: Not enough images for {image}. Skipping...")
            return TetradCounter.Counts(0, 0), tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count
        image_paths_dict = {
        'blue': next(img for img in image if img.endswith("1.tiff")),
        'red': next(img for img in image if img.endswith("2.tiff")),
        'yellow': next(img for img in image if img.endswith("3.tiff"))
    }
        tetrad_predictor = predict.TetradPredictor(tetrad_model=self.tetrad_model)
        n_predictions, n_tetrads, n_triads,  scaled_transform_list, pred_class_thr   = tetrad_predictor.predict_tetrads(image_paths_dict)

        if n_predictions == 0:
            return TetradCounter.Counts(0, 0) , tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count

        color_classifier = predict.TetradClassifier(color_model=self.color_model)
        detected_tetrad, detected_triad, detected_sum, color_count_dict, tetrad_counts = color_classifier.classify_tetrads(scaled_transform_list,pred_class_thr)
        
        # Update predicted tetrad classes into dictionary
        self._update_counts(detected_tetrad, tetrad_type_dict)
        self._update_counts(detected_triad, triad_type_dict)
        self._update_counts(detected_sum, sum_type_dict)  
        self._update_counts(color_count_dict, color_dict)
        if tetrad_counts is not None:
            self._update_counts(tetrad_counts, combination_count)
        return TetradCounter.Counts(n_tetrads, n_triads), tetrad_type_dict, triad_type_dict, sum_type_dict, color_dict, combination_count
          
    def _update_counts(self, source_dict: Dict[str, int], target_dict: defaultdict) -> None:
        """
        Utility function to update the target dictionary with counts from the source dictionary.

        Parameters:
        - source_dict (Dict[str, int]): Dictionary containing counts to be added.
        - target_dict (defaultdict): Dictionary to which counts are to be added.
        """
        for k, v in source_dict.items():
            target_dict[k] += v
            
    @staticmethod
    def create_empty_excel_file(filename, sheet_names=['sum', 'tetrad', 'triad']):
        """
        Create an empty Excel file with given sheet names.

        Parameters:
        - filename (str): Name of the Excel file to be created.
        - sheet_names (list): List of sheet names to be added to the Excel file.
        """
        try:
            with pd.ExcelWriter(filename) as writer:
                for sheet in sheet_names:
                    pd.DataFrame().to_excel(writer, sheet_name=sheet)
        except ValueError as e:
            # This logs the error message with a timestamp
            logging.error(f"Failed to create an Excel writer for the file '{filename}': {e}")
    
        
def main():
    """Main entry point of the script."""
    parser = argparse.ArgumentParser(description='Run Tetrad Detection')

    parser.add_argument("--tetrad_model", default='./weights/tetrad_triad_prediction_weights.pth', 
                        help="Path to the tetrad/triad model. Default: args.tetrad_model")
    
    parser.add_argument("--color_model", default='./weights/color_classification_weights.pth', 
                        help="Path to the color model. Default: args.color_model")

    parser.add_argument("-p", "--path", required=False, default="./", metavar="input path to fluorescent images",
                        help="Path to images to run TetradDetection on (default: .)")
    
    parser.add_argument("-o", "--output", required=False, default="./outputs/", metavar="path to save excel results",
                    help="Directory to save the results (default: ./)")
    
    parser.add_argument("-v", "--verbose", action="store_true", 
                    help="Enable verbose logging")
    
    parser.add_argument('--excel_file', default="default.xlsx", help="Name of the Excel file.")
    
    parser.add_argument('--prefix', default='sample1', help="Prefix name (e.g. slx4).")
    
    parser.add_argument('--output_extent', choices=['limited', 'extended'], default='limited',
                        help="Extent of the output to Excel file. Choices: 'limited' (default) or 'extended'. Extended version contains all outputs from TetradCalculator.")
    
    args = parser.parse_args()

    # Setting up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    tetrad_counter = TetradCounter(args.path, args.output, args.prefix, args.excel_file,args.output_extent, args.tetrad_model, args.color_model)
    tetrad_counter.run_tetrad_count_all_process()
    
if __name__ == "__main__":
    main()



