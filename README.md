linear/integer programing solver

based on https://github.com/dieram3/competitive-programming-library/blob/master/include/cpl/math/simplex.hpp

- slightly reworked the primal simplex to set artificial variables only if necessary
- added the possibility to use bound for variables
- added dual simplex (add very little value by itself, used for gomory after primal complete so preconditions don't need to be checked)
- added an integer search using naive gomory cuts (fails often)

The sudoku part uses not optimized gaussian/gaussian-jordan algos to transform
the combinatorial pb to linear integer prog pb. in order to test the simplex solver. Env. > 99% of random
sudoku grids are solved and > 75% of very hard grids are solved. It seems that for the most part failure
is due to cuts not really relevant being selected and not sufficient care is taken in manipulating
floating point numbers.

```
$ c++ -O2 -std=c++11 main.cpp -o blpsolver
$ echo .............7..52.85.2...3...4....7..2....3.31...96......83.1.16...2.4.9.8.4.... | ./blpsolver 
.............7..52.85.2...3...4....7..2....3.31...96......83.1.16...2.4.9.8.4....
mat size 209 sol size 25
gaussian_algo                  m 149 n 209
gaussian_jordan_algo           m 149 n 209
get_free_vars_augm_matrix      m 123 n 61
get_free_vars_augm_matrix      m 115 n 61
primal pivot cnt 120
primal pivot cnt 59
dual pivot cnt 17
dual pivot cnt 17
dual pivot cnt 15
dual pivot cnt 72
dual pivot cnt 96
dual pivot cnt 47
dual pivot cnt 66
dual pivot cnt 103
dual pivot cnt 106
dual pivot cnt 15
gomory cnt 11
res 13.000000
0.000000 1.000000 0.000000 -0.000000 1.000000 -0.000000 1.000000 1.000000 -0.000000 1.000000 -0.000000 1.000000 1.000000 -0.000000 -0.000000 -0.000000 -0.000000 -0.000000 -0.000000 -0.000000 1.000000 1.000000 -0.000000 -0.000000 -0.000000 1.000000 0.000000 -0.000000 -0.000000 1.000000 0.000000 -0.000000 -0.000000 -0.000000 -0.000000 1.000000 1.000000 -0.000000 
mat size 209 sol size 81
239516478641378952785924163856431297492867531317259684524683719163792845978145326
solved 1 / 1 100.000% time grid  49363.000 us time total 49363 us
```
  
-------
This code is released under the Boost Software License, Version 1.0.

Copyright (C) 2018 rafirafi.
