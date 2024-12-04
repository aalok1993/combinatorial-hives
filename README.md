# On the randomized Horn problem and the surface tension of hives

[Link for the paper](https://arxiv.org/abs/2410.12619).


## Octahedron Recurrence



## Lozenge Tiling

Here is an illustration of Speyer's Matching for obtaining the lozenge tiling of a point in the double hive.

![Speyer's Matching](https://github.com/aalok1993/combinatorial-hives/blob/main/res/Speyers_Matching.gif?raw=true)

To obtain the lozenge tiling, use the following command
`python lozenge.py --matrix_type GUE --N 100 --V_i 50 --V_j 50 --sig_X 1 --sig_Y 1 --LA_diff 1000`

**matrix_type** denotes the type of matrix to be used for the eigenvalue minor process to obtain a double. **N** denotes the size of the double hive. (**V_i**,**V_j**) denote the x-coordinate and y-coordinate of the base point V in the interior of the double hive for which the lozenge tiling will be generated. **sig_X** and **sig_Y** represents the scaling factor of the first and second matrix, respectively.
**LA_diff** denotes the difference factor in Lambda used for converting a pair of Gelfand-Tsetlin patterns into a double hive. This needs to be set sufficiently high.

## Surface Tension
