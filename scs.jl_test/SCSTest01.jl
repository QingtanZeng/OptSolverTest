using ECOS
using SparseArrays

"""
min x₁² + x₂² 
s.t. [1,1   [x₁    [2
     -1,1] * x₂] =  0]
     
The optimal solution is clearly x₁=1, x₂=1, and cost=2.
The equivalent linear cone problem is
min t
s.t. Ax=b
     [t;[x1;x2]] <₌K 0

"""
# --- Problem Data ---
# 1. Define the problem in ECOS standard form
# z = [x₁, x₂, t]
n = 3 # Number of variables

# Objective: minimize y=c'z
c = [0.0, 0.0, 1.0]

# Equality constraints: Ax = b
# x₁ = 1
# x₂ = 2
A = sparse([1.0 1.0 0.0; -1.0 1.0 0.0])
b = [2.0; -2]
p = length(b)

# Conic constraints: h - Gx ∈ K
# (x₁² + x₂²) <= t²  is equivalent to [t, x₁, x₂] ∈ SecondOrderCone(3)
# We represent [t, x₁, x₂] as h - Gx.
# Let h = 0, then -Gx = [t, x₁, x₂].
G = sparse([0.0 0.0 -1.0;
    -1.0 0.0 0.0;
    0.0 -1.0 0.0])
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
