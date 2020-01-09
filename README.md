# SPQR
Penalized Kernel Quantile Regression for Varying Coefficient Models



## Overview

*SPQR* is a novel and efficient penalized kernel smoothing algorithm that consistently identifies the partially linear structure in the varying coefficient quantile regression model.
Efficient ADMM algorithm is developed for implementation. 


## Main functions

- [main_admm.m]: Perform the proposed ADMM algorithm to obtain SPQR estimate.

- [simulation_example.m]: Perform simulation analysis when n=400, p=7, and tau=0.25.

## Supplementary functions

- [keru.m]:  Epanechnikov kernel function used for constructing kernels.

- [scad.m]: Scad derivation function used to obtain weight parameters in SPQR.

- [rq.m]: Solve conventional quantile regression problem using the dual problem. Created by Roger Koenker (http://www.econ.uiuc.edu/~roger/research/rq/rq.html).



## Author

* Eun Ryung Lee and Seyoung Park
 Department of Statistics, Sungkyunkwan University


## Contact

* ishspsy@skku.edu


## License

This project is licensed under the MIT License.




