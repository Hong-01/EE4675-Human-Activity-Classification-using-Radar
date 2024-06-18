# EE4675-Human-Activity-Classification-using-Radar
Project Code of EE4675 Object Classification with Radar

_Author: Alexandros Theocharous | Yanqi Hong_

# File Structure and Organization

This project is organized into three main folders: `Code`, `Documentation`, and `Processed Data`. Each folder contains specific components essential for the project's execution and understanding.

## Table of Contents
    1. [Code]
        1. [Data Processing]
        2. [Neural Network]
    2. [Documentation]
    3. [Processed Data]

## Code 

    The `Code` folder is divided into two main sections: Data Processing and Neural Network.

    ### Data Processing

        This section contains MATLAB scripts and functions necessary for processing the radar data.

        #### Main Script
            - **Data_Processing.m**: This is the main script to process the data. Execute this file to start the data processing workflow. It orchestrates the execution of various functions and handles the data flow.

        #### Functions
            - **CFAR.m**: Implements the Constant False Alarm Rate (CFAR) detection algorithm. This function is crucial for identifying targets in the radar data. Note: This function may take a significant amount of time to execute, depending on the dataset size and complexity.
            - **DetectionMatProcess.m**: Processes the detection matrix obtained from the CFAR detection. This function further refines the data for subsequent use.
            - **Label_extract4.m**: A utility function derived from the dataset example code, used for extracting labels.
            - **oddnumber.m**: Another utility function from the dataset example code, performing specific operations as needed.

        #### CFAR Detection Results
            To expedite the processing, this folder contains precomputed CFAR detection results. Using these files allows you to bypass the CFAR detection step, saving time.
            - **detection_DR.mat**: Detection results for dataset DR. (DR: Doppler Range)
            - **detection_DT.mat**: Detection results for dataset DT. (DT: Doppler Time)
            - **detection_RT.mat**: Detection results for dataset RT. (RT: Range Time)

    ### Neural Network

        This section includes Jupyter notebooks for training and evaluating different neural network models on the processed radar data and the saved best model for CNN.
        #### Notebooks
            - **radar-classification-cnn.ipynb**: Notebook for training a Convolutional Neural Network (CNN) to classify radar data.
            - **radar-classification-mobilenetv2.ipynb**: Notebook for training a MobileNetV2 model for radar data classification.
            - **radar-classification-resent.ipynb**: Notebook for training a ResNet (Residual Network) model, known for its deep architecture and effective training on large datasets.
        #### Saved Model
            - **model_best_cnn.keras**: The best model obtained from training the CNN model. It is also available in the Kaggel online platform. Please follow the link in the **radar-classification-cnn.ipynb**.

## Documentation

    The `Documentation` folder contains the project's overview report.
    - **Report.pdf**: Overview for this project.
    - **Presentation radar classification.pdf**: Presentation Slides

## Processed Data

    The `Processed Data` folder includes the processed data files, ready for use in the neural network training and evaluation steps.
    - **image_data_cfar.zip**: This zip file contains the image data resulting from the CFAR detection process. It is formatted and ready for input into the neural network models provided in the `Neural Network` section.
