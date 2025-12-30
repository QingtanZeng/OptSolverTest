## Overview
The repo serves the following three development objectives:
1. implementation of conic optimization algorithms at principle level;
2. test open-source convex solvers' interface and performance in Julia and C/C++;
3. customized solvers based on open-source projects, for embedded system (CPU + Parallel acceleration).

Based on solvers' classification in three dimensions,

<b>{First-Order, Second-Order} x {IPM, Exterior Point Method(Penalty)} x {Primal, Primal-dual}</b>,

only the three types of solvers and corresponding projects are focused, for real-time trajectory generation's OCP on embedded system:
1. <b>Homogeneous Self-dual Embedding(HSDE) IPM</b> for Linear SOCP :   ECOS
2. <b>Homogeneous Embedding(HE) IPM</b> for Quadratic SOCP:             Clarabel 
3. <b>HE Douglas-Rachford Splitting(DRS)</b> for Quadratic SOCP:             SCS \
   and a similar variant, Proportional-integral projected gradient method (PIPG)

<b>HE DRS</b> is preferred for real-time platform, due to the advantages:
1. <b>Infeasibility Certification</b> concurrent with convergence, by Homogeneous Embedding;
2. Converges rapidly to local <b>modest precision</b> solution, by first-order splitting iteration;
3. <b>Low computational</b> load of single iteration, by conical projection and fixed LDLT decomposition, \
   compared to updated decomposition of each iteration in IPM's newton direction;
4. <b>Warm-start</b>, beneficial for MPC.


## An Implementation of Douglas-Rachford Splitting(DRS) for Quadratic SOCP
An example:

<img width="200" height="90" alt="image" src="https://github.com/user-attachments/assets/713afa10-2484-43a5-94dc-db51dcf95037" />

The optimal solution is clearly x₁=1, x₂=1, and cost=2.

The convergence process is shown below,
where <b>O(1/k) sublinear convergence rate</b> is obvious
that gap, percent of equality error, percent of primal-dual variables' change is linear 
in the logarithmic coordinate. \
Therefore converge rapidly to modest precision solution with <3% error, which is acceptable under most 
real-time trajectory generation cases. <b>Feasibility is far more important than optimality and high precision</b>,
especially requiring >10Hz updates to deal with constantly changing environment, sudden disturbance, or MPC tasks.

<p align="center">
<img alt="Convergence of DRS for Quadratic SOCP"
    title="Convergence of DRS for Quadratic SOCP"
    src="drs_impl/convergence_plot.png"
    width="600px" />
</p>
