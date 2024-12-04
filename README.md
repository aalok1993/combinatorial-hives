# On the randomized Horn problem and the surface tension of hives

[Link for the paper](https://arxiv.org/abs/2410.12619).

## Hives

Consider the triangle $`T_n = \{(i,j) \in \mathbb{Z}^2 \mid 0\leq i \leq j \leq n\}`$, for some fixed value of $n$. The figure below on the left illustrates the triangle for $n=5$. The triangle $T_n$ contains three types of rhombi having the following form.

$$\textit{Type 1: } (A,B,C,D) = ((i,j),(i+1,j),(i+2,j+1),(i+1,j+1))$$
$$\textit{Type 2: } (A,B,C,D) = ((i,j),(i+1,j+1),(i+1,j+2),(i,j+1)$$
$$\textit{Type 3: } (A,B,C,D) = ((i,j),(i,j-1),(i+1,j-1),(i+1,j)$$

Here, in all three types of rhombi, $AC$ is the long diagonal, whereas $BD$ is the short diagonal. The figure below on the right illustrates the three types of rhombi. A hive is a function $h_n : T_n \to \mathbb{R}$ on the triangle that is rhombus-concave, that is, for each rhombus $ABCD$ in $T_n$,

$$h_n(A)+h_n(C)\leq h_n(B)+h_n(D)$$

The triangle $T_5$             |  The three types of rhombi
:-------------------------:|:-------------------------:
<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Triangle.png?raw=true" width="360">  |  <img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Rhombi.png?raw=true" width="360">

## Gelfand-Tsetlin pattern

Consider the pattern $`\gamma=(\lambda_{j,k})_{1\leq j\leq k\leq n}`$ where $\lambda_{j,k}$ are real numbers. $\gamma$ is called a Gelfand-Tsetlin pattern if it satisfies the interlacing condition  

$$\lambda_{j,k+1} \geq \lambda_{j,k} \geq \lambda_{j+1,k+1}$$

The figure below depicts the Gelfand-Tsetlin pattern, wherein each $\lambda_{j,k}$ is greater than or equal to the number immediately to the southeast or northeast of it.

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/GT_pattern.png?raw=true" width="360">

## Constructing hives using the minor process and large eigengaps

Given a $n \times n$ Hermitian matrix $A$, let $\lambda_{1,k} \geq \dots \geq \lambda_{k,k}$ be the eigenvalues of the top left $k \times k$ minor of $A$ for $1 \leq k \leq n$. Then, $(\lambda_{j,k})_{1 \leq j \leq k \leq n}$ forms a Gelfand-Tsetlin pattern. A weakly decreasing sequence $\Lambda$ is said to have large gaps if

$$\min_{1 \leq i < n} \Lambda_i - \Lambda_{i+1} > \lambda_{1,n} - \lambda_{n,n}$$

Using $(\lambda_{j,k})_{1 \leq j \leq k \leq n}$ and $\Lambda$ we can construct the hive $h_n : T_n \to \mathbb{R}$ 
as depicted in the figures below.

Minor Process with $n=3$ |  Minor Process with $n=3$
:-------------------------:|:-------------------------:
<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Hive_MP1.png?raw=true" width="360">  |  <img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Hive_MP2.png?raw=true" width="360">

Let $X$ and $Y$ be two Hermitian matrices of size $n \times n$. Let $\mu_{1,k} \geq \dots \geq \mu_{k,k}$ and $\nu_{1,k} \geq \dots \geq \nu_{k,k}$ be the eigenvalues of the top left $k \times k$ minor of $X$ and $Y$, respectively for $1 \leq k \leq n$. Let $\Lambda$ be a weakly decreasing sequence with large gaps. Then, using 
$`(\mu_{j,k})_{1 \leq j \leq k \leq n}`$, $`(\nu_{j,k})_{1 \leq j \leq k \leq n}`$, and $\Lambda$ we obtain a double hive as depicted in the figure below.

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/MinorProcess.png?raw=true" width="360">



## Octahedron Recurrence

Consider the tetrahedron 

$$` \mathtt{tet} = \{ [x,y,z,w] \in \mathbb{Z}^4 \mid x,y,z,w \geq 0; x+y+z+w=n \}`$$

with vertices $R=(n,0,0,0), P=(0,n,0,0), Q=(0,0,n,0), S=(0,0,0,n)$ as illustrated in the figure below. It consists of four hives on the boundary surface: $PQS$, $QRS$, $PQR$, and $PRS$. The hives $PQS$ and $QRS$ form the top panel double hive, whereas $PQR$ and $PRS$ form the lower panel double hive. A point lying on each of these triangles is of the form

$$PQS: (0, j-i, n-j, i)$$
$$QRS: (i-j, 0, n-i, j)$$
$$PQR: (i, j, n-i-j,0)$$
$$PRS: (n-j,n-i,0,i+j-n)$$   

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Tetrahedron.png?raw=true" width="360">

Let $o:\mathtt{tet}\to \mathbb{R}$ be a function that assigns a real value to each vertex of the tetrahdron.
The tetrahedron consists of several unit sized tetrahedra and several octahedra. Each octahedron inside the tetrahedron consists of six vertices having values and indices of the form 

$$a = o([x+1,y,z+1,w])$$
$$b = o([x+1,y,z,w+1])$$
$$c = o([x,y+1,z,w+1])$$
$$d = o([x,y+1,z+1,w])$$
$$e = o([x+1,y+1,z,w])$$
$$f = o([x,y,z+1,w+1])$$

For each such octahedron inside $\mathtt{tet}$, we have the octahedron rule which asserts that $f = \max(a+c,b+d)-e$ as depicted in the figure below.

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/Octahedron.png?raw=true" width="360">

We now explain a method for the exact sampling of special types of augmented hives. Using the minor process explained earlier, we obtain a double hive. A tetrahedron is constructed and the top panel consisting of triangle $PQS$ and $QRS$ is initialized with values of the double hive. We then sequentially excavate the tetrahedron using the octahedron recurrence to obtain the value at each point in the tetrahedron. This gives us the double hive on the bottom panel.

We consider four special types of matrices for which we sample augmented hives: 

  1. Gaussian Unitary Ensembles (GUE)
  2. Random Projections (RP)
  3. Sequential Spectrum Matrix (SSM)
  4. Alternating Sequential Spectrum Matrix (ASSM)

For sampling augmenting hives using octahedron recurrence for the four different types of matrices use the following commands. The command and the resulting augmented hives are shown below.

`python octahedron_recurrence.py --matrix_type GUE --M 1000 --N 250 --sig_X 1 --sig_Y 1 --LA_diff 1000`

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/GUE_1_1.png?raw=true" width="720">

`python octahedron_recurrence.py --matrix_type RP --M 1000 --N 250 --sig_X 1 --sig_Y 1 --LA_diff 1000`

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/RP.png?raw=true" width="720">

`python octahedron_recurrence.py --matrix_type SSM --M 1000 --N 250 --sig_X 1 --sig_Y 1 --LA_diff 1000`

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/SSM.png?raw=true" width="720">

`python octahedron_recurrence.py --matrix_type ASSM --M 1000 --N 250 --sig_X 1 --sig_Y 1 --LA_diff 1000`

<img src="https://github.com/aalok1993/combinatorial-hives/blob/main/res/ASSM.png?raw=true" width="720">

In the above commands, **_matrix_type_** denotes the type of matrix to be used for the eigenvalue minor process to obtain a double hive. 
**_M_** denotes the number of matrices to sample for calculating the average augmented hive.
**_N_** denotes the size of the double hive. 
**_sig_X_** and **_sig_Y_** represent the scaling factor of the first and second matrices, respectively.
**_LA_diff_** denotes the difference factor in Lambda used for converting a pair of Gelfand-Tsetlin patterns into a double hive. This needs to be set sufficiently high.

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
