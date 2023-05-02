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
Please contact me for any questions: y.guo at wsu dot edu dot au. 
