This project is an implementation of an _interaction expansion continuos time quantum Monte-carlo algorithm for fermions_. It is designed to perform numerical solution of Anderson impurity model in continuous time domain. Written in mixture of c and c++ at Moscow State University by A. N. Rubtsov et al.

**Corresponding paper - http://arxiv.org/abs/cond-mat/0411344**

The development of this code is stopped. The latest source download (07.11.2011) is available as <a href='http://code.google.com/p/ct-qmc/downloads/detail?name=ct-qmc-1.6.2.tar.gz'>ct-qmc-1.6.2.tar.gz</a>

This code supports MPI-parallelization for calculating the Green's function in Matsubara representation. This version supports multi-orbital calculations. The irreducible vertex parts are provided in a single-cpu version only.

Look into the file default.h for the explanation of the parameters of the code.


---

Other projects:
  * http://pomerol.googlecode.com. An exact-diagonalization code aimed at solving condensed matter second-quantized models of interacting fermions on finite size lattices at finite temperatures.