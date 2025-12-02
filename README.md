# BMmod - Markov Chain Spatial Modeling Software

![Fortran](https://img.shields.io/badge/Fortran-f2018-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Build-CLI-lightgrey)
![Platforms](https://img.shields.io/badge/OS-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)
![Maintainer](https://img.shields.io/badge/Maintainer-Bosszz-purple)

## Overview

BMmod is a modern Fortran 2008/2018 codebase for spatial Markov chain modeling and spectral methods. It provides numerical routines for balancing, Hessenberg reduction, QR eigenvalue analysis, and spectral decomposition, plus domain functions to derive transition rates, probabilities, and 3-D models.

## Features

- **Markov Chain Modeling**: Build 1-D and 3-D Markov chain models for spatial analysis
- **Multiple Rate Matrix Options**:
  1. Specified transition rates
  2. Transition probabilities at specified lag
  3. Embedded Markov chain transition probabilities
  4. Embedded Markov chain transition frequencies
  5. Independent/Maximum entropy transition frequencies
  6. Volumetric proportions
  7. Frequency of embedded occurrences
- **Linear Algebra Utilities**: Implements robust EISPACK-style algorithms
  - Matrix balancing (`balanc`)
  - Hessenberg reduction (`elmhes`)
  - QR eigenvalue computation (`hqr`)
  - Spectral decomposition (`spectral`)
- **3-D Model Generation**: Constructs 3-D transition probability models
- **CPT Data Integration**: Includes MATLAB scripts for CPT (Cone Penetration Test) data analysis

## Project Structure

```
Msource_BM/
├── BMmod.f90          # Main Fortran source (entry: program BMmod)
├── exam_site/         # Example site data and outputs
│   ├── BMmod.par      # Parameter file for example
│   ├── data.*         # Various output files
│   └── exam_site.eas  # Example output matrix
└── CPTsoft/           # MATLAB utilities for CPT analysis
    ├── SBT_cloud_fig.m    # Soil behavior type classification
    ├── SCdata.csv         # Soil classification data
    ├── exam_CPT.csv       # Example CPT data
    ├── *.png              # Generated plots
    └── soil_classification_results.txt
```

Note: Compiled artifacts (e.g., `.exe`, `.mod`) are not tracked; build locally.

## Installation

### Prerequisites

- Fortran compiler (gfortran, Intel Fortran, or compatible with Fortran 2008/2018)
- MATLAB for CPT analysis utilities

## Usage

### Basic Usage

1. Prepare a parameter file (see `exam_site/BMmod.par` for example)
2. Run the program:
   ```bash
   ./BMmod.exe
   ```
3. When prompted, enter the name of your parameter file

### Parameter File Format

The parameter file should contain:

```
[# of categories]
[proportions for each category]
[background category index]
[name of debugging file]
[name of 3-D model output file]
[name of determinant file]
[determinant limits for 3-D model; x,y,z direction]
[dx,dy,dz for 3-D model]
[... direction-specific parameters ...]
```

### Example Parameter File

See `exam_site/BMmod.par` for a complete example with 4 categories.

## CPT Analysis (MATLAB)

The `CPTsoft` directory contains MATLAB scripts for soil classification based on CPT data:

### Running CPT Analysis

```matlab
% In MATLAB, navigate to CPTsoft directory
cd CPTsoft

% Run the analysis script
SBT_cloud_fig
```

### Input Files

- `SCdata.csv`: Soil classification statistics (means, std devs, correlations)
- `exam_CPT.csv`: Test points for classification

### Output Files

- `SBT_probability_density.png`: Probability density plots for each soil type
- `SBT_posterior_probability.png`: Posterior probability maps
- `test_points_classification.png`: Classification results for test points
- `soil_classification_results.txt`: Numerical classification results

## Technical Details

### Modules

1. **bmmod_types**: Precision kinds and size constants
2. **bmmod_data**: Global model parameters and shared state
3. **bmmod_linalg**: Linear algebra utilities (balancing, Hessenberg, QR)
4. **bmmod_core**: Domain logic for Markov chain operations

### Key Algorithms

- **Balancing**: Reduces matrix norm for improved numerical stability
- **Hessenberg Reduction**: Converts matrix to upper Hessenberg form
- **QR Algorithm**: Computes eigenvalues using QR decomposition
- **Spectral Decomposition**: Constructs spectral components for transition matrices

## Output Files

- `.eas` files: Transition probability matrices at different lags
- `.bgr` files: 3-D model data (binary format)
- `.dbg` files: Debugging information and diagnostics

## Contributing

This software is provided as-is for educational and research purposes. The code follows modern Fortran best practices and is fully documented.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Bosszz — Created on 2025-12-02

## Citation

If you use BMmod in research, please cite:

```
Bosszz (2025). BMmod: Markov Chain Spatial Modeling Software.
MIT Licensed. https://github.com/your-org/BMmod
```

## Disclaimer

This software is provided "as-is", without warranty of any kind. The author assumes no responsibility for errors or omissions in the software or documentation.
