# clouldremoval
Cloud simulation and removal for remote sensing monochromatic images 

# Main file to run
multitests.m
It sweep through many values of regularisation parameter $\lambda$. 

# Others
- simcloud.m: the simulation code to generate simulated cloud images (pure cloud)
- perlin_noise.m: Perlin noise generator. The main generator for cloud. Only one version works and slow. 
- ATMcouldremovers.m: the main function file for 4 methods for cloud removal, all RPCA based. 
- various m files: supporting files for ATMcouldremovers.m (solvers)

# Question
Please contact me for any questions. 

# Paper
- Yi Guo, Feng Li & Zhuo Wang (2023). **Laplacian convolutional representation for traffic time series imputation**. International Journal of Remote Sensing. [[DOI](https://doi.org/10.1080/01431161.2023.2208710)][[Preprint](https://arxiv.org/abs/2210.01981)] [[Matlab code](https://github.com/yguo-wsu/clouldremoval)]

# Citation
```
@article{doi:10.1080/01431161.2023.2208710,
author = {Guo,Yi and Li,Feng and Wang,Zhuo},
title = {Cloud Removal Using Scattering Model and Evaluation via Semi-realistic simulation},
journal = {International Journal of Remote Sensing},
year = {2023},
doi = {10.1080/01431161.2023.2208710},
}
```
