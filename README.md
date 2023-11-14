# Motion Robust CMR

This repository contains the implementation of reconstruction algorithms discussed in the journal article, **["Motion-robust free-running cardiovascular MRI"](https://arxiv.org/abs/2308.02088)**, currently available on arXiv. The primary focus of this study is to enhance the robustness of Cardiovascular Magnetic Resonance (CMR) imaging against motion artifacts.
<p align="center">
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/c1340d64-1922-45a8-ad7f-f759ef4a93e5" width="200"/>
</p>

## About the Project

The "Motion Robust CMR" project aims to address the challenges posed by motion artifacts during CMR imaging. Employing advanced reconstruction algorithms presented in our research work, this project demonstrates significant improvements in image quality and robustness in scenarios where motion artifacts are unavoidable. Our proposed algorithm, CORe (Compressive Recovery with Outlier Rejection), incorporates outlier rejection to suppress artifacts.

### Key Features

- Implementation of known reconstruction algorithms for CMR, such as Compressed Sensing 'CS' (Lustig et al., 2008), Robust Regression 'RR' (Nikolova et al., 2004), and Sparse Outliers 'SO' (Dong et al., 2012).
- Introduction of CORe (Compressive Recovery with Outlier Rejection), a novel method focusing on outlier rejection to improve motion robustness.
- Multiple studies on the comparison of algorithms:
   - Study I: Static phantom
   - Study II: Dynamic phantom
   - Study III: High-resolution 3D cine
   - Study IV: Rest 4D flow
   - Study V: Stress 4D flow
- ADMM (Alternating Direction Method of Multipliers)/Split Bregman implementation of reconstruction algorithms.

## Getting Started

To use these algorithms, clone the repository and follow the setup instructions provided below.

### Installation

Requires MATLAB 2019 or later.

### Usage
- Study I and II utilize simulated phantom data, which is generated directly within the main MATLAB script.
- Studies III, IV, and V employ in-vivo CMR datasets, available at [Cardiovascular MRI (CMR) 3D Cine, 4D Flow, and Exercise Stress 4D Flow Datasets](https://zenodo.org/records/8105485).

### [Datasets](https://zenodo.org/records/8105485)
The datasets for Studies III, IV, and V consist of Cartesian undersampled acquisitions with self-gating. The provided .mat files contain self-gated expiratory phase k-space data, sorted into 20 cardiac phases.


### Cite as:

Arshad SM, Potter LC, Chen C, Liu Y, Chandrasekaran P, Crabtree C, Han Y, Ahmad R (2023). Motion-robust free-running cardiovascular MRI. arXiv preprint arXiv:2308.02088.
