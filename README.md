# Neural-network-for-Wiener-Hopf-kernel-factorization

## Introduction

We propose an approximation method based on the neural network to solve the Wiener-Hopf kernel matrix factorization problem, 
which takes advantage of the automatic differentiation feature of the neural network to simplify the calculation of the Cauchy-Riemann equation, without the need of mesh generation. 
To validate the proposed method, we study a class of parallel-plate acoustic diffraction problems with rigid and soft boundary conditions, which have explicit solutions for the matrix kernel factorization, facilitating the comparison with the neural network optimization results. 
However, despite the potential of the neural network to map complex nonlinear functions, it is still unsatisfactory in fitting poles and infinity points. 
Therefore, an additional loss function based on the boundary conditions is added during the training process, to ensure that the near-field results meet the requirements. 
Compared with the explicit solutions, the results obtained by the neural network decomposition are consistent in both far-field and near-field, 
indicating that it has the potential to become a new choice for solving the matrix kernel factorization problem.

## Environment

Tensorflow 2.4 / Cuda 11.0 / Python 3.8

## Running cases

1. MINN.ipynb: running the notebook to factorize the Wiener-Hopf kernel. Integral path xy_pre.mat and k0_pre.mat can be generated from MATLAB code.

2. Scattering.m: running the MATLAB code to see the parallel plate scattering results.

## Reference

1.  @article{raissi2017physicsI,
  title={Physics Informed Deep Learning (Part I): Data-driven Solutions of Nonlinear Partial Differential Equations},
  author={Raissi, Maziar and Perdikaris, Paris and Karniadakis, George Em},
  journal={arXiv preprint arXiv:1711.10561},
  year={2017}
}
**Github Repo** : https://github.com/maziarraissi/PINNs

2. @article{hurd1981scattering,
  title={Scattering by hard and soft parallel half-planes},
  author={Hurd, RA and L{\"u}neburg, E},
  journal={Canadian Journal of Physics},
  volume={59},
  number={12},
  pages={1879--1885},
  year={1981},
  publisher={NRC Research Press Ottawa, Canada}
}
