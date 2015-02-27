Symmetric Tensor Eigen-Rank-One Iterative Decomposition (STEROID)
-----------------------------------------------------------------

The Symmetric Tensor Eigen-Rank-One Iterative Decomposition (STEROID) decomposes an arbitrary symmetric tensor A into a real linear combination of unit-norm symmetric rank-1 terms.

1. Functions
------------

* [V,d,lambdas,e,tail]=steroid(A) or [V,d,lambdas,e,tail]=steroid(A,method)

Use this function to compute the STEROID of a symmetric tensor A.

* A=randsymten(d,n)

Creates a random symmetric tensor of order d and dimension n.

* e=symcheck(A)

Checks the symmetry of A for all permutations of its indices.

* B=embed(A)

Embeds a symmetric d-way tensor A with odd d into a symmetric d+1-way tensor B such that B(:,:,...,:,1)=A.

* [lindex index]=exp2ind(mons)

Converts homogeneous monomial exponents into a linear and tensor index of a d-way symmetric tensor of dimension n.

* [U,V,D,e]=svdsteroid(A,method)

Combination of an SVD and STEROID for a tensor A that is symmetric in all but the first mode.

* demo.m

Small demo that illustrates the use of steroid.m

2. Reference
------------

"Symmetric Tensor Decomposition by an Iterative Eigendecomposition Algorithm"

http://arxiv.org/abs/1409.4926


Authors: Kim Batselier, Ngai Wong