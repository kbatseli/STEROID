Symmetric Tensor Eigen-Rank-One Iterative Decomposition (STEROID)
-----------------------------------------------------------------

The Symmetric Tensor Eigen-Rank-One Iterative Decomposition (STEROID) decomposes an arbitrary symmetric tensor A into a real linear combination of unit-norm symmetric rank-1 terms.

1. Functions
------------

* [V,D,lambdas,e,X,tail]=steroid(A,desiredd)

Use this function to compute the STEROID of a symmetric tensor A for the desired order 'desiredd'.

* A=randsymten(d,n)

Creates a random symmetric tensor of order d and dimension n.

* e=symcheck(A)

Checks the symmetry of A for all permutations of its indices.

* demo.m

Small demo that illustrates the use of steroid.m

2. Reference
------------

"Symmetric Tensor Decomposition by an Iterative Eigendecomposition Algorithm"

http://arxiv.org/abs/1409.4926


Authors: Kim Batselier, Ngai Wong