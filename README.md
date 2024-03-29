# TetradDetection

## Description
High-throughout _S.cerevisiae_ tetrad image analysis pipeline. The pipeline relies on two deeplearning models for accurate detection of tetrads and triads followed by subsequent classification into corresponding tetrad types based on recombination (non-crossover, crossover). The pipeline takes input fluorescent images (mCerulean, tdTomato, YFP) and outputs calculated map distance, interference, non-discjuntion and gene-conversion frequencies in an .xlsx file.

![Pipeline Architecture](./assets/pipeline_architecture.png)

## Installation

### 1. Install Anaconda
[Download an Anaconda distribution](https://www.anaconda.com/download)

### 2. Create a virtual environment and activate it
```python
conda create -n TetradDetection python=3.11.7
conda activate TetradDetection
```

### 3. Install TetradDetection in Windows
```python
git clone https://github.com/Lisby-Lab/TetradDetection.git
```

### 4. Install modules using requirements.txt
```python
cd TetradDetection/
```

```python
pip install -r requirements.txt
```

### 5. Download zip files for the neural network
[Tetrad Prediction weights](https://sid.erda.dk/share_redirect/GY4SXrNM8T)

[Color Classification weights](https://sid.erda.dk/share_redirect/ChqIgVyH5h)

### 6. Unzip the weights and place them in the weights folder

## Usage
### Prerequisites
1. Standard computer with enough RAM to apply the neural networks. 8 GB RAM is enough to predict the 1080 x 1080 px image provided as the test data. The RAM requirements scale with the number of image pixels.
2. Single-channel grayscale fluorescent .tiff/.tif images captured at 60x magnification with a 425 nm laser (mCerulean), 561 nm laser (tdTomato) and 488 nm laser (YFP).
3. Images must be named following these guidelines, where the identifier indicates the sequential order of the image. Specifically, the identifier denotes the position of the image in the sequence: e.g., the first image would have an identifier of 1, the second image an identifier of 2, the third image an identifier of 3, and so on. Image Name can be changed according to preference. The suffix indicates which fluorescent channel the image belongs to.
   
| Name      | Identifier | Suffix | Example Filename       |
|-----------|------------|--------|------------------------|
| mCerulean | \_1\_     | \_1    | mCerulean_1_1.tiff      |
| tdTomato  | \_1\_     | \_2    | tdTomato_1_2.tiff       |
| YFP       | \_1\_     | \_3    | YFP_1_3.tiff            |


<table border="0">
  <tr>
    <td><img src="./assets/cyan.png" alt="CFP"></td>
    <td><img src="./assets/red.png" alt="RFP"></td>
    <td><img src="./assets/yellow.png" alt="YFP"></td>
  </tr>
</table>

4. Within the TetradDetection environment open a new Terminal and write type in:
```python
python3 entry.py --path ./imgs/
```
5. After the program finishes a new default.xlsx file will be created on the root folder containing calculated Map Distance, Interference, Non-Disjunction and Gene Conversion frequency values.

### Side note:
A tetrad detection only model is also available for [download](https://sid.erda.dk/share_redirect/B2zhGjLLiu). Unzip the weights and place them in the weights folder. 
Within the Terminal type in:
```python
python3 entry.py --path ./imgs/ --tetrad_model ./weights/tetrad_only_prediction_weights.pth
```

## Support
Please contact bszucs@sund.ku.dk if you have any questions.


