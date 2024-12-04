# On the randomized Horn problem and the surface tension of hives

[Link for the paper](https://arxiv.org/abs/2410.12619).

## Hives

Consider the triangle $T_n = \{(i,j) \in \mathbb{Z}^2 \mid 0\leq i \leq j \leq n\}$ for some fixed value of $n$. The figure below on the left illustrates the triangle for $n=5$. The triangle $T_n$ contains three types of rhombi having the following form.

$$\textit{Type 1: } (A,B,C,D) = ((i,j),(i+1,j),(i+2,j+1),(i+1,j+1))$$
$$\textit{Type 2: } (A,B,C,D) = ((i,j),(i+1,j+1),(i+1,j+2),(i,j+1)$$
$$\textit{Type 3: } (A,B,C,D) = ((i,j),(i,j-1),(i+1,j-1),(i+1,j)$$

Here, in all three types of rhombi, $AC$ is the long diagonal, whereas $BD$ is the short diagonal. The figure below on the right illustrates the three types of rhombi. A hive is a function $h_n : T_n \to \mathbb{R}$ on the triangle that is rhombus-concave, that is, for each rhombus $ABCD$ in $T_n$,

$$h_n(A)+h_n(C)\leq h_n(B)+h_n(D).$$

The triangle $T_5$             |  The three types of rhombi
:-------------------------:|:-------------------------:
<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Triangle.png?raw=true" width="360">  |  <img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Rhombi.png?raw=true" width="360">

## Gelfand-Tsetlin pattern

Consider the pattern $\gamma = (\lambda_{j,k})_{1 \leq j \leq k \leq n}$, where $\lambda_{j,k}$ are real numbers. $\gamma$ is called a Gelfand-Tsetlin pattern if it satisfies the interlacing condition  

$$\lambda_{j,k+1} \geq \lambda_{j,k} \geq \lambda_{j+1,k+1}$$

The figure below depicts the Gelfand-Tsetlin pattern, wherein each $\lambda_{j,k}$ is greater than or equal to the number immediately to the southeast or northeast of it.

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/GT_pattern.png?raw=true" width="360">

## Octahedron Recurrence



## Lozenge Tiling

Here is an illustration of Speyer's Matching for obtaining the lozenge tiling of a point in the double hive.

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Speyers_Matching.gif?raw=true" width="360">
<!-- ![Speyer's Matching](https://github.com/aalok1993/combinatorial-hives/blob/main/res/Speyers_Matching.gif?raw=true) -->

To obtain the lozenge tiling, use the following command

`python lozenge.py --matrix_type GUE --N 100 --V_i 50 --V_j 50 --sig_X 1 --sig_Y 1 --LA_diff 1000`

In the above command, **_matrix_type_** denotes the type of matrix to be used for the eigenvalue minor process to obtain a double hive. 
**_N_** denotes the size of the double hive. 
(**_V_i_**,**_V_j_**) denote the x-coordinate and y-coordinate of the base point V in the interior of the double hive for which the lozenge tiling will be generated. 
**_sig_X_** and **_sig_Y_** represent the scaling factor of the first and second matrices, respectively.
**_LA_diff_** denotes the difference factor in Lambda used for converting a pair of Gelfand-Tsetlin patterns into a double hive. This needs to be set sufficiently high.

## Surface Tension
