import gc
import cv2
from skimage.io import imread, imsave
import torchvision.transforms.functional as transform
import torch
import numpy as np
from PIL import Image
import asyncio
import aiofiles
import io
from typing import List, Dict, Optional, Tuple, DefaultDict
from torchvision.transforms import ToTensor
from torch.utils.data import Dataset, DataLoader
from collections import defaultdict
from skimage import exposure, img_as_ubyte, img_as_float32
from utils import tetrad_types
import logging
from skimage.transform import resize

class TetradPredictor:
    """
    Class for analyzing and predicting tetrads from image data.
    
    Attributes:
        tetrad_model (torch.nn.Module): The loaded tetrad prediction model.
        CLASS_NAMES (list): List of class names used in prediction.
        tetrad_confidence (float): Confidence threshold for predictions.
        device (torch.device): Device to use for torch computations.
    """

    CLASS_NAMES = ['__background__', 'tetrad', 'triad']

    def __init__(self, tetrad_model: torch.nn.Module, tetrad_confidence: float = 0.85, device: Optional[torch.device] = None):
        """
        Initializes the TetradPredictor.

        Args:
            tetrad_model (torch.nn.Module): Pre-loaded tetrad prediction model.
            tetrad_confidence (float): Confidence threshold for predictions.
            device (torch.device, optional): Device to use for torch computations. If None, uses CUDA if available.
        """
        self.tetrad_model = tetrad_model
        self.tetrad_confidence = tetrad_confidence
        self.device =  torch.device('cpu') 

    def resize_to_model_input(self, img: np.ndarray, target_size=(1080, 1080)) -> np.ndarray:
        """
        Resize image to model input size if needed.
        Preserves intensity range and dtype.
        """
        if img.shape == target_size:
            return img
        
        if img.shape != (1080, 1080):
            logging.info(f"Resizing image from {img.shape} to (1080, 1080)")
      
        img_resized = resize(
            img,
            output_shape=target_size,
            order=3,               # bicubic
            preserve_range=True,   # critical for microscopy
            anti_aliasing=False
        )
        
        return img_resized.astype(img.dtype)
    
    def predict_tetrads(self, prefix: Dict[str, str]) -> Tuple[int, int, int, List[np.ndarray], List[np.ndarray], List[str]]:
        """
        Predicts tetrads from the given image paths.

        Args:
            prefix (dict): Dictionary containing paths to images with keys 'blue', 'red', 'yellow'.

        Returns:
            tuple: Contains the number of predictions, number of tetrads, number of triads, list of cropped RGBA images, 
                   list of scaled transformed images, and a list of predicted classes.
        """
        image_paths = [prefix['blue'], prefix['red'], prefix['yellow']]
        image_data = [imread(path) for path in image_paths]
        
        # Resize each channel to model input size
        image_data = [
            self.resize_to_model_input(img, target_size=(1080, 1080))
            for img in image_data
        ]
        print(image_data[0].shape)
        # Rescale images and merge
        blue, red, yellow = [img_as_ubyte(exposure.rescale_intensity(img)) for img in image_data]
        merged_rescaled = cv2.merge((blue, red, yellow))
        # Convert image to tensor and perform prediction
        tensor_image = transform.to_tensor(merged_rescaled)
        outputs = self.tetrad_model([tensor_image])        
        # Process predictions
        pred_score = outputs[0]['scores'].detach().cpu().numpy()
        above_threshold_indices = pred_score > self.tetrad_confidence
        if not any(above_threshold_indices):
            return 0, 0, 0, [], []
        
        pred_classes = outputs[0]['labels'][above_threshold_indices].tolist()
        n_tetrads = pred_classes.count(1)
        n_triads = pred_classes.count(2)
        
        # Extract bounding boxes of predictions
        bboxes = outputs[0]['boxes'][above_threshold_indices].int()
        cropped_rgba_img_list = [merged_rescaled[ymin:ymax, xmin:xmax] for xmin, ymin, xmax, ymax in bboxes.tolist()]
        
        # Rescale intensity of cropped images
        scaled_transform_list = [cv2.merge([img_as_ubyte(exposure.rescale_intensity(channel)) for channel in cv2.split(cropped_img)]) 
                                 for cropped_img in cropped_rgba_img_list]

        # Convert numerical labels to string labels
        pred_class_labels = [self.CLASS_NAMES[label] for label in pred_classes]
        
        if 'triad' in pred_class_labels:
            pred_class_labels = ['triad'] * len(pred_class_labels)
        return len(pred_classes), n_tetrads, n_triads, scaled_transform_list, pred_class_labels

'''''''''
Classify tetrad types based on color predictions of segregating fluorescent markers.
- Cast tetrad list to torch.Dataset and torch.DataLoader
- Predict color combinations with color detection model
- Cast color predictions to lookup table for classification
'''''''''

def custom_collate_fn(batch):
    return batch

class TetradClassifier:
    def __init__(self, color_model: Optional[torch.nn.Module] = None, color_confidence: float = 0.85, 
                 device: Optional[torch.device] = None):
        """
        Initialize TetradClassifier.

        Parameters:
        color_model : Model used for color prediction.
        color_confidence : Threshold confidence for a color prediction to be accepted.
        device : Torch device to be used. If None, the code will use CUDA if available.
        """
        self.CLASS_NAMES = ['__background__', 'blue' , 'orange', 'red', 'yellow', 'empty', 'white', 'purple', 'green']
        self.color_model = color_model
        self.color_confidence = color_confidence
        self.device =  torch.device('cpu') 
        
    def classify_tetrads(self, scaled_transform_list: List[np.ndarray], pred_class_thr)  -> Dict[str, int]: 
        """
        Classify a list of cropped images into tetrads.
        
        Parameters:
        cropped_rgba_img_list : List of cropped images to be classified.
        scaled_transform_list : List of transformed images to be classified.
        
        Returns:
        Dictionary counting the number of each type of tetrad.
        """
        batch_color_outputs = self._load_and_apply_model(scaled_transform_list)
        tetrad_type_dict, triad_type_dict, sum_type_dict, sum_color_dict, tetrad_counts = self._process_predictions(scaled_transform_list, batch_color_outputs, pred_class_thr)
        return tetrad_type_dict, triad_type_dict, sum_type_dict, sum_color_dict, tetrad_counts
    
    def _load_and_apply_model(self, scaled_transform_list: List[np.ndarray]) -> List[Dict[str, torch.Tensor]]:
        """
        Load and apply the color model to a set of transformed images.

        Parameters:
        scaled_transform_list : List of transformed images to be classified.

        Returns:
        List of model outputs for each image.
        """
        color_dataset = TetradClassifier.ColorDataset(scaled_transform_list)
        dataset_iterator = DataLoader(color_dataset, batch_size=10, shuffle=False, num_workers=1, collate_fn=custom_collate_fn)

        batch_color_outputs = []
        for images in dataset_iterator:
            image = [img.to(torch.device('cpu')) for img in images]
            with torch.no_grad():
                color_outputs = self.color_model(image)
                batch_color_outputs.extend(color_outputs)
        
        return batch_color_outputs

    def _process_predictions(
            self, 
            scaled_transform_list: List[np.ndarray], 
            batch_color_outputs: List[Dict[str, torch.Tensor]], 
            tetrad_pred_class_thr: List[str]
    ) -> Tuple[DefaultDict[str, int], DefaultDict[str, int], Dict[str, int]]:
        """
        Process a list of images with their respective model outputs to count tetrads and triads.

        Parameters:
        - scaled_transform_list: List of images to process.
        - batch_color_outputs: Corresponding model outputs for each image.
        - tetrad_pred_class_thr: Predicted classifications for each image.

        Returns:
        - Tuple of three dictionaries: tetrad counts, triad counts, and combined counts.
        """
        # Default dict to collect various tetrad types
        tetrad_type_dict = defaultdict(int)
        triad_type_dict = defaultdict(int)
        sum_type_dict = defaultdict(int)
        
        # Default dict to collect color counts
        tetrad_color_count = defaultdict(int)
        triad_color_count = defaultdict(int)
        sum_color_dict = defaultdict(int)
        combined_tetrad_counts = {}
        # Process each image and its corresponding prediction
        for i, cropped_image in enumerate(scaled_transform_list):
            current_prediction = tetrad_pred_class_thr[i]
            if current_prediction == 'tetrad':
                tetrad_counts = self._type_classification(4, 4, tetrad_type_dict, tetrad_color_count, current_prediction, cropped_image, batch_color_outputs[i])
                triad_counts = {}
            elif current_prediction == 'triad':
                triad_counts = self._type_classification(4, 3, triad_type_dict, triad_color_count, current_prediction, cropped_image, batch_color_outputs[i])
                tetrad_counts = {}
                
        # Combine tetrad and triad counts
        for key in set(tetrad_type_dict) | set(triad_type_dict):
            sum_type_dict[key] = tetrad_type_dict.get(key, 0) + triad_type_dict.get(key, 0)
            
        # Combine tetrad and triad color counts
        for key in set(tetrad_color_count) | set(triad_color_count):
            sum_color_dict[key] = tetrad_color_count.get(key, 0) + triad_color_count.get(key, 0)
          
        if tetrad_counts is None:
            tetrad_counts = {}
        if triad_counts is None:
            triad_counts = {}
        
        for key in set(tetrad_counts) | set(triad_counts):
            combined_tetrad_counts[key] = tetrad_counts.get(key, 0) + triad_counts.get(key, 0)    
        
        return tetrad_type_dict, triad_type_dict, sum_type_dict, sum_color_dict, combined_tetrad_counts


    def _type_classification(self,
        pred_class_limit: int, 
        pred_class_count: int, 
        tetrad_or_triad_dict: DefaultDict[str, int], 
        color_count_dict: DefaultDict[str,int],
        current_prediction: str, 
        cropped_image: np.ndarray, 
        color_outputs: Dict[str, torch.Tensor]
                                    ) -> None:
        """
        Processes individual predictions based on given parameters.

        Parameters:
        - pred_class_limit: The limit on the number of classes to consider.
        - pred_class_count: Expected number of classes for tetrad or triad.
        - tetrad_or_triad_dict: Dictionary to update with counts.
        - color_count_dict: Dictionary to collect color distribution accross images
        - current_prediction: 'tetrad' or 'triad'.
        - cropped_image: Current image being processed.
        - color_outputs: Model outputs for the current image.
        """

        # Convert tensor scores to a list of prediction scores
        pred_score = list(color_outputs['scores'].detach().cpu().numpy())

        # If only one score, update broken mask count
        if len(pred_score) == 1:
            tetrad_or_triad_dict['broken_mask'] += 1
            return

        # Get the threshold of prediction scores that are above the color confidence 
        # Exception in case not a single tetrad or triad is detected within the image e.g. empty image, non-sporulated cells.
        try:
            pred_threshold = [pred_score.index(x) for x in pred_score if x > self.color_confidence][-1]
        except IndexError:
            tetrad_or_triad_dict['predictions_below_confidence'] += 1
            return

        # Extract class names based on labels and prediction threshold
        pred_class = [self.CLASS_NAMES[j] for j in list(color_outputs['labels'][:pred_class_limit])]
        pred_class_thr = pred_class[:pred_threshold+1]
        
        # Check if number of classes is below expected count
        if len(pred_class_thr) < pred_class_count:
            tetrad_or_triad_dict[f'less_than_{pred_class_count}_spores'] += 1
            return
        
        for color in pred_class_thr:
            color_count_dict[color] += 1

        bboxes = []
        # Process individual object data in prediction
        for object_id in range(len(pred_class_thr)):
            xmin, ymin, xmax, ymax = color_outputs['boxes'][object_id]
            bboxes.append((int(xmin), int(ymin), int(xmax), int(ymax)))


        # If bounding boxes don't pass the threshold, update bounding box error count
        if TetradClassifier.bboxes_threshold(bboxes):
            tetrad_or_triad_dict['bbox_error'] += 1
            return

        # Sort predicted classes and look up tetrad type
        sorted_pred_class = pred_class_thr.copy()
        sorted_pred_class.sort()
        current_type, tetrad_counts = tetrad_types.lookup_table(sorted_pred_class)
        tetrad_or_triad_dict[current_type] += 1
        
        if not isinstance(tetrad_counts, dict):
            raise ValueError("Expected tetrad_counts to be a dictionary.")
        
        return tetrad_counts
    
    @staticmethod
    def bboxes_threshold(bboxes: List[Tuple[int]]) -> bool:
        flag = False
        for box1 in bboxes:
            for box2 in bboxes:
                iou = TetradClassifier.intersection_over_union(box1, box2)
                if 1 > iou > 0.5:
                    flag = True
                    return flag
        return flag

    @staticmethod
    def intersection_over_union(box1: Tuple[int], box2: Tuple[int]) -> float:
        # Get coordinates of the intersection 
        x1 = max(box1[0], box2[0])
        y1 = max(box1[1], box2[1])
        x2 = min(box1[2], box2[2])
        y2 = min(box1[3], box2[3])

        # Get the area of intersection rectangle
        intersection = max(0, x2 - x1 + 1) * max(0, y2 - y1 + 1)

        # Get the area of both rectangles
        box1Area = (box1[2] - box1[0] + 1) * (box1[3] - box1[1] + 1)
        box2Area = (box2[2] - box2[0] + 1) * (box2[3] - box2[1] + 1)

        iou = intersection / float(box1Area + box2Area - intersection)
        return iou
            
            

    class ColorDataset(torch.utils.data.Dataset):
        """
        Torch Dataset for color images.
        -Initialize dataset with a list of images.
        """
        def __init__(self, img_list: List[np.ndarray]):
            super(TetradClassifier.ColorDataset, self).__init__()
            self.img_list = img_list
            self.to_tensor = ToTensor()
        """
        Returns the length of the dataset.
        """
        def __len__(self) -> int:
            return len(self.img_list)
        """
        Fetch an image by its index.
        """
        def __getitem__(self, idx: int) -> torch.Tensor:
            img = self.img_list[idx]
            img = self.to_tensor(img) 
            return img


