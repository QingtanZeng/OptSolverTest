## Overview
The repo serves the following three development objects:
1. implementation of conic optimization algorithms at principle level;
2. test open-source convex solvers' interface and performance in Julia and C/C++;
3. customized solvers based on open-source projects, for embeded system (CPU + Parallel acceleration).

Based on solvers' classification in three dimensions, <b>{First-Order, Second-Order} x {IPM, Exterior Point Method(Penalty)} x {Primal, Primal-dual}</b>, \
only the three types of solvers and responding projects are focused, for real-time trajectory generation's OCP on embeded system:
1. <b>Homogeneous Self-dual Embedding(HSDE) IPM</b> for Linear SOCP :   ECOS
2. <b>Homogeneous Embedding(HE) IPM</b> for Quadratic SOCP:             Clarabel 
3. <b>HE Douglas-Rachford Splitting(DRS)</b> for Quadratic SOCP:             SCS \
   and a similar variant, Proportional-integral projected gradient method (PIPG)

<b>HE DRS</b> is preferred for real-time platform, due to the advantages:
1. <b>Infeasibility Certification</b> concurrent with convergence, by Homogeneous Embedding;
2. Converge rapidly to local <b>modest precision</b> solution, by first-order splitting iteration;
3. <b>Low computational</b> load of single iteration, by conical projection and fixed LDLT decomposition, \
   compare to updated decomposition of each iteration in IPM's newton direction;
4. <b>Warm-start</b>, beneficial for MPC.
   
