# Motion Robust CMR

This repository contains the implementation of reconstruction algorithms discussed in the journal article, **["Motion-robust free-running cardiovascular MRI"](https://arxiv.org/abs/2308.02088)**, currently available on [arXiv](https://arxiv.org/abs/2308.02088). The primary focus of this project is our innovative reconstruction algorithm, CORe (Compressive Recovery with Outlier Rejection), which is specifically designed to suppress motion artifacts.
<p align="center">
  <h3 align="center">Exercise Stress Imaging</h3>
</p>
<p align="center">
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/5ca0ff76-b84a-433f-ab53-3df9becaaf6e" height="150"/>
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/1d14c2dc-047f-4a81-96e1-8940af78a3a0" height="150"/>
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/fe8ad0c8-60b3-4d99-bbe3-794dab6009a7" height="150"/>
  <!-- <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/bffd8697-082f-41a7-810b-2556f7321d41"height="100"/>-->
  <br>
  <em>Conventional CS Reconstruction</em>
  <br>
  <br>
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/457984bd-d556-4409-a0dc-e6d967d0fa4f" height="150"/>
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/6873abb3-a6b2-4a3c-b824-a01b100d7404" height="150"/>
  <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/b0a7ce42-0aac-4bd7-8f53-5bf85c7f6f01" height="150"/>
  <br>
  <em>Proposed CORe Reconstruction</em>
<!--  <em>Conventional CS Reconstruction</em></span>-->
  <!-- <img src="https://github.com/OSU-MR/motion-robust-CMR/assets/97550963/db49819d-3aa0-4614-ad37-46f904f9bf22" width="200" height="20"/>-->
 <!-- <em>Proposed CORe Reconstruction</em></span>-->
</p>

## About the Project

The "Motion Robust CMR" study aims to address the challenges posed by motion artifacts during CMR imaging. Our proposed reconstruction algorithm, CORe incorporates outlier rejection to suppress motion artifacts. This work compares the performance of CORe against known reconstruction algorithms in different CMR techniques and scenarios, dicussed in our [research work](https://arxiv.org/abs/2308.02088). The results demonstrate significant improvements in image quality and robustness in scenarios where motion artifacts are unavoidable.

### Key Features
- Introduction of CORe (Compressive Recovery with Outlier Rejection), a novel method equiped with outlier rejection to improve motion robustness.
- Implementation of known CMR reconstruction algorithms for comparison, such as,
  - Compressed Sensing 'CS' (Lustig et al., 2008)
  - Robust Regression 'RR' (Nikolova et al., 2004)
  - Sparse Outliers 'SO' (Dong et al., 2012).
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

Arshad SM, Potter LC, Chen C, Liu Y, Chandrasekaran P, Crabtree C, Han Y, Ahmad R (2023). Motion-robust free-running cardiovascular MRI. arXiv preprint [arXiv:2308.02088](https://arxiv.org/abs/2308.02088).

## Authors 

- **Syed Murtaza Arshad**, Primary Developer, [GitHub Profile](https://github.com/syedmurtazaarshad), arshad.32@osu.edu
- **Rizwan Ahmad**, PI, [GitHub Profile](https://github.com/OSU-CMR), ahmad.46@osu.edu
