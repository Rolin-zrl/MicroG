# MicroG
《Micro-Structures Graph-Based Point Cloud Registration for Balancing Efficiency and Accuracy》[[arxiv](http://arxiv.org/abs/2410.21857)]
## News
- [2024-11-17]: [[Early Access](https://ieeexplore.ieee.org/document/10755047)]
- [2024-10-26]: Our paper has been accpeted by **TGRS** (IEEE Transactions on Geoscience and Remote Sensing).
## Table of Contents
1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Usage](#usage)
4. [Data-preparation](#data-preparation)
5. [Acknowledgments](#acknowledgments)
6. [Citation](#citation)

## Introduction
 This paper introduces a micro-structures graph-based coarse-to-fine global point cloud registration method. This method employs a hierarchical outlier removal strategy based on graph nodes and edges, combined with the GNC-Welsch estimator, to ensure robustness during coarse registration. At finer scales, PA-AA optimization is utilized to further exploit the geometric features of corresponding micro-structures, enhancing accuracy with minimal additional computational cost.
 <img src="/imgs/framework.png" alt="teaser" width="100%">

We have conducted experiments on ETH and 3DMatch dataset.

 <img src="/imgs/teaser.png" alt="teaser" width="50%">
 <img src="/imgs/resultsETH.png" alt="teaser" width="100%">



## Dependencies
We've tested it on CLion 2024 running on Ubuntu 18.04. 
</br>A compiler that supports OpenMP.
- CMake >= 3.10
- PCL >=1.12
- Eigen3 >=3.3.0
## Usage
1. Clone the repository:
   ```bash
   git clone git@github.com:Rolin-zrl/MicroG.git
2. Compile
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ./ MicroG "../data/s1.ply" "../data/s1.ply" "./data/s1-s2.tfm" 0.1 800 2

#### Input six parameters introduced in the paper
1. source_path
2. target_path
3. gt_path *(optional)*
4. resolution: *downsampling resolution*
5. K: *default 800*
6. L: *default 2.0*
## Data-preparation

You can test on the online available point cloud data and registration dataset, such as 
</br>[WHU TLS Registration Dataset](https://3s.whu.edu.cn/ybs/en/benchmark.htm),
</br>[ETH PRS TLS Registration Dataset](https://prs.igp.ethz.ch/research/completed_projects/automatic_registration_of_point_clouds.html),
</br>[ETH ASL Robotics Registration Dataset](https://projects.asl.ethz.ch/datasets/doku.php?id=laserregistration:laserregistration),
</br>[3D Match](https://3dmatch.cs.princeton.edu/),
</br>[Robotic 3D Scan Repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/),
etc.

## Acknowledgments
People who inspired this idea, gave suggestions and reported errors, especially [Pengcheng Wei](https://github.com/WPC-WHU) contributing to the removal of outliers.

## Citation
@ARTICLE{10755047,
  author={Zhang, Rongling and Yan, Li and Wei, Pengcheng and Xie, Hong and Wang, Pinzhuo and Wang, Binbing},
  journal={IEEE Transactions on Geoscience and Remote Sensing}, 
  title={Micro-Structures Graph-Based Point Cloud Registration for Balancing Efficiency and Accuracy}, 
  year={2024},
  volume={},
  number={},
  pages={1-1},
  keywords={Point cloud registration;correspondence graph;robust estimator;planar adjustment;Anderson acceleration},
  doi={10.1109/TGRS.2024.3488502}}


