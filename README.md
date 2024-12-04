# On the randomized Horn problem and the surface tension of hives

[Link for the paper](https://arxiv.org/abs/2410.12619).


## Octahedron Recurrence



## Lozenge Tiling

Here is an illustration of Speyer's Matching for obtaining the lozenge tiling of a point in the double hive.

![Speyer's Matching](https://github.com/aalok1993/combinatorial-hives/blob/main/res/Speyers_Matching.gif?raw=true)

For obtaining the lozenge tiling, use the following command
`python lozenge.py --matrix_type GUE --N 100 --V_i 50 --V_j 50 --sig_X 1 --sig_Y 1 --LA_diff 1000`


`
Usage: lozenge.py [-h] [--matrix_type {GUE,RP,SSM,ASSM}] [--N N] [--V_i V_I] [--V_j V_J] [--sig_X SIG_X] [--sig_Y SIG_Y] [--LA_diff LA_DIFF]

Generate Lozenge Tilings using Speyers Perfect Matching

optional arguments:
  -h, --help            show this help message and exit
  --matrix_type {GUE,RP,SSM,ASSM}
                        Type of matrix to use for minor process
  --N N                 Size of the double hive
  --V_i V_I             x-coordinate of the base point V in the interior of the double hive for which the lozenge tiling is to be generated.
  --V_j V_J             y-coordinate of the base point V in the interior of the double hive for which the lozenge tiling is to be generated.
  --sig_X SIG_X         Scaling factor for X
  --sig_Y SIG_Y         Scaling factor for Y
  --LA_diff LA_DIFF     Difference factor in large lambda
`

## Surface Tension
