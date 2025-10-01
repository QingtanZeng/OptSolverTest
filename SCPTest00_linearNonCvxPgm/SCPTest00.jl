using ECOS
using LinearAlgebra, SparseArrays

"""
Solve the following linear Objective but nonconvex constraints problem
using standard conic programming with second-order IPM method.
minimize     z₂
subject to   z₂ - z₁² = 0
             z₂ + 0.1*z₁ = 0.6
             z₁ ≤ 0.85 
"""
# define common parameters
iterMax = 20
zref_0 = [10, 100]
tolabs = 1e-2
# define common structure
z_hist = [zeros(7) for k=1:iterMax]
gap = zeros(iterMax)
pcost=zeros(iterMax)

mutable struct ConPgm
    z::Vector{Float64}
    gap::Float64
    pcost::Float64
end
conpgm=ConPgm([zref_0; zeros(5)], -1.0, -1.0)

function main()
    exitcode = Int(11)
    flgFea = false
    flgOpt = false
    itrEnd = 0;

    for itr=1:iterMax
        itrEnd = itr
        println("--- Iteration $itr ---")

        exitcode = ScpSubPrm!(conpgm.z[1:2], conpgm)
        z_hist[itr] = copy(conpgm.z)
        z1 = conpgm.z[1];   z2 = conpgm.z[2];
        gap[itr] = conpgm.gap
        pcost[itr] = conpgm.pcost

        vcabs = maximum([abs.([z2-z1^2, z2+0.1*z1-0.6]); z1-0.85])

        flgFea = (exitcode == ECOS.ECOS_OPTIMAL) && ( vcabs< tolabs)
        # flgOpt = itr>1? (norm(z_hist[itr][1:2] - z_hist[itr-1][1:2])<tolabs) : false
        if itr > 1
            # 计算当前解与上一个解的距离
            zDlt = norm(z_hist[itr][1:2] - z_hist[itr-1][1:2])
            flgOpt = zDlt < tolabs
        else
            # 第一次迭代没有之前的解可以比较
            flgOpt = false
        end

        println("variables:z=$(conpgm.z[1:2]); cost:$(conpgm.pcost); gap:$(conpgm.gap); 
                Vc:$(conpgm.z[3]); Tr:$(conpgm.z[5])")

        if flgFea && flgOpt 
            println("--- END: Feasible and Optimal ---") 
            println()
            break
        elseif flgFea && !flgOpt
            println("---Iteration $itr  Continue: Feasible but Unoptimal ---")
        else
            println("--- Iteration $itr Continue: Unfeasible ---")
        end
        println()
    end

    println("variables: z = $(z_hist[itrEnd][1:2]); cost: $(pcost[itrEnd]); gap: $(gap[itrEnd])")
    println("--- END: itr=$itrEnd, flgFea=$flgFea, flgOpt=$flgOpt ---")
    return itrEnd
end

function ScpSubPrm!(zref::Vector{<:Real}, conpgm::ConPgm)::Int
# define common parameters
    wvc = 1e2
    wtr = 10

# 1. Define the problem in ECOS standard form
#     [z,     ,v',eta_vc, chi_tr ]
# z = [[z₁,z₂],v',eta_vc, [eta_tr, mu_1x2]]
    n = 7 # Number of variables

# Objective function
# z = [[z₁,z₂],v',eta_vc, [eta_tr, mu_1x2]]
# c = [[0, 1], 0, wvc,    wtr[1,   0_1x2] ]
    c = [[0, 1]; [0, wvc]; wtr*[1 ; [0, 0]]]

# Equality constraints: Ax = b
# z =      [[z₁,         z₂], v',eta_vc, [eta_tr, mu_1x2]]
# 1. 线性化的非凸约束: -2*z₁_ref*z₁ + z₂ + v = -z₁_ref²
# 2. 线性约束: 0.1*z₁ + z₂ = 0.6
# 3. 信赖域约束: z_orig - μ = z_ref
    A = sparse([ [-2*zref[1]  1]   1  0       zeros(1,3);
                [0.1          1]  0  0       zeros(1,3);
                -1.0*I(2)        zeros(2,2)  [zeros(2,1) 1.0*I(2)] ])
    b = [-zref[1]^2; 0.6; -zref]
    p = length(b)

# Conic constraints: s = h-Gx ∈ K
# {{eta_vc,v'}_(1+1), {eta_tr, mu_1x2}_(1+2)} ∈ K₂xK₃ 
# represent by s = h-Gx ∈ K 
# Positive orthant: Fz<=d >> d-Fz>=0 >> d-Fz∈K_posi
# SOC: h-Gz <=K 0
# K_pos (l=1): z₁ <= 0.85  => 0.85 - z₁ >= 0
# K_soc1 (q[1]=2): ||v||₂ <= η_vc  => [η_vc, v] ∈ Q²
# K_soc2 (q[2]=3): ||μ||₂ <= η_tr  => [η_tr, μ₁, μ₂] ∈ Q³
    G = sparse([    [1.0 0       zeros(1,2)  zeros(1,3)];
                    -[zeros(1,2)  0  1        zeros(1,3);
                    zeros(1,2)  1  0        zeros(1,3);
                    zeros(3,4)              1.0*I(3)  ]
                ])
    h = [0.9;[0;0];zeros(3)]
    m = length(h) # Total dimension of cone constraints

# positvie cone and exponential constraints
    l=1
    nex=0

# Cone dimensions: One second-order cone of dimension 3.
    q =Int[2,3]
    ncones=2

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
    exitcode = ECOS.ECOS_solve(pwork)

# Copy variables to structure
    if exitcode == ECOS.ECOS_OPTIMAL
        pwork_loaded = unsafe_load(pwork)
        info = unsafe_load(pwork_loaded.info)

        # --- CRITICAL FIX: Safely copy the solution vector ---
        # The original `unsafe_copyto!` is dangerous and causes memory corruption.
        # The correct way is to create a safe, managed copy of the result.
        conpgm.z = copy(unsafe_wrap(Array, pwork_loaded.x, pwork_loaded.n))
        conpgm.gap = info.gap
        conpgm.pcost = info.pcost
    end

    ECOS.ECOS_cleanup(pwork, 0)

    return exitcode

end

main()