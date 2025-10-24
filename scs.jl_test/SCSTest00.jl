using ECOS
using SparseArrays

"""
min y 
s.t. sqrt(x₁² + x₂²) ≤ y x₁ = 1 x₂ = 2
     
This is a classic problem where we are finding the minimum Euclidean norm of a 2D vector 
whose components are fixed. 
The optimal solution is clearly x₁=1, x₂=2, and y = sqrt(1² + 2²) = sqrt(5) ≈ 2.236.
"""
# --- Problem Data ---
# 1. Define the problem in ECOS standard form
# z = [x₁, x₂, y]
n = 3 # Number of variables

# Objective: minimize y=c'z
c = [0.0, 0.0, 1.0]

# Equality constraints: Ax = b
# x₁ = 1
# x₂ = 2
A = sparse([1.0 0.0 0.0; 0.0 1.0 0.0])
b = [1.0, 2.0]
p = length(b)

# Conic constraints: h - Gx ∈ K
# sqrt(x₁² + x₂²) <= y  is equivalent to [y, x₁, x₂] ∈ SecondOrderCone(3)
# We represent [y, x₁, x₂] as h - Gx.
# Let h = 0, then -Gx = [y, x₁, x₂].
G = sparse(-[0.0 0.0 1.0;
    1.0 0.0 0.0;
    0.0 1.0 0.0])
h = [0.0, 0.0, 0.0]
m = length(h) # Total dimension of cone constraints
# Cone dimensions: One second-order cone of dimension 3.
q =Int[3]
ncones=1
# no positvie cone and exponential constraints
l=0
nex=0

# Set the problem
pwork = ECOS.ECOS_setup(
    n,
    m,
    p,
    l,
    ncones,
    q,
    nex,
    G.nzval,
    G.colptr .-1,
    G.rowval .-1,
    A.nzval,
    A.colptr .-1,
    A.rowval .-1,
    c,
    h,
    b,
)
# --- Solve the problem ---
# Call ECOS.solve with the defined parameters
solu = ECOS.ECOS_solve(pwork) 
ECOS.ECOS_cleanup(pwork, 0)
