# DSP Circle Detection Project

This repository contains the final project for the DSP course (2023), which focuses on circle detection in images using various methods such as correlation and gradient analysis.

## Overview

In many image processing applications, detecting a specific shape is crucial. This project aims to detect circles, which are insensitive to rotation, with two different approaches:

- **Known radius detection**: Detect circles of a fixed size.
- **Unknown radius detection**: Detect circles of varying sizes.

The project is implemented in MATLAB, and several key methods are used:
1. **Direct correlation** for detecting known-size circles.
2. **FFT-based correlation** for faster detection.
3. **Gradient vector pairs** to identify circles with unknown radii.

## Files Structure

The project consists of four parts, each handled by a different function:
- `direct_correlator.m`: Implements direct correlation for detecting circles.
- `fft2_correlator.m`: Uses FFT to speed up the circle detection process.
- `circle_locator.m`: Detects circles using gradient vectors for varying radii.

### Key Directories:
- `Functions/`: Contains the MATLAB functions you'll edit for the project.
- `Subroutines/`: Contains additional scripts to support the main functions.
- `Results/`: This folder stores the output figures generated during execution.
- `input-images/`: Contains the input images for testing the algorithms.

## Getting Started

To run the project:
1. Clone this repository to your local machine.
2. Open MATLAB and navigate to the project folder.
3. Run the `main.m` file using either F5 or the "Run" button.

## Usage

1. **Run the main code**: The `main.m` script controls the overall flow of the project. It clears pre-defined parameters, adds necessary paths, and runs the desired part of the project.
   
2. **Enable different parts**: You can change the parameter `Part_Enable` in the `main.m` file to select the portion of the project you want to run:
   - `'A1'`: Known circle detection in grayscale images using direct correlation.
   - `'A2'`: Known circle detection in grayscale images using FFT-based correlation.
   - `'A3'`: Known circle detection in color images.
   - `'B'`: Unknown radius circle detection in color images using gradient vector pairs.

3. **Run and analyze**: After selecting a part, run the code and analyze the results stored in the `Results/` folder.

## Output

- The results include figures showing detected circles and correlation analysis.
- Execution time for each method is displayed in the MATLAB console.

## Requirements

- MATLAB R2016b or later
