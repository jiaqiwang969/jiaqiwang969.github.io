/**
@page step_43 The step-43 tutorial program
This tutorial depends on step-31.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Advectiondominatedtwophaseflowmathematicalmodel">Advection-dominated two-phase flow mathematical model.</a>
        <li><a href="#Adaptiveoperatorsplittingandtimestepping">Adaptive operator splitting and time stepping.</a>
        <li><a href="#Timediscretization">Time discretization.</a>
        <li><a href="#Weakformspacediscretizationforthepressurevelocitypart">Weak form, space discretization for the pressure-velocity part.</a>
        <li><a href="#Stabilizationweakformandspacediscretizationforthesaturationtransportequation">Stabilization, weak form and space discretization for the saturation transport equation.</a>
        <li><a href="#Adaptivemeshrefinement">Adaptive mesh refinement.</a>
        <li><a href="#Linearsystemanditspreconditioning">Linear system and its preconditioning.</a>
        <li><a href="#Thetestcases">The test cases.</a>
        <li><a href="#Listofreferences">List of references</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Boundaryandinitialvalueclasses">Boundary and initial value classes</a>
        <li><a href="#Permeabilitymodels">Permeability models</a>
        <li><a href="#Physicalquantities">Physical quantities</a>
        <li><a href="#Helperclassesforsolversandpreconditioners">Helper classes for solvers and preconditioners</a>
        <li><a href="#TheTwoPhaseFlowProblemclass">The TwoPhaseFlowProblem class</a>
        <li><a href="#TwoPhaseFlowProblemdimTwoPhaseFlowProblem">TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem</a>
        <li><a href="#TwoPhaseFlowProblemdimsetup_dofs">TwoPhaseFlowProblem<dim>::setup_dofs</a>
        <li><a href="#Assemblingmatricesandpreconditioners">Assembling matrices and preconditioners</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemdimassemble_darcy_preconditioner">TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner</a>
        <li><a href="#TwoPhaseFlowProblemdimbuild_darcy_preconditioner">TwoPhaseFlowProblem<dim>::build_darcy_preconditioner</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_darcy_system">TwoPhaseFlowProblem<dim>::assemble_darcy_system</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_system">TwoPhaseFlowProblem<dim>::assemble_saturation_system</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_matrix">TwoPhaseFlowProblem<dim>::assemble_saturation_matrix</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs_cell_term">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs_boundary_term">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term</a>
      </ul>
        <li><a href="#TwoPhaseFlowProblemdimsolve">TwoPhaseFlowProblem<dim>::solve</a>
        <li><a href="#TwoPhaseFlowProblemdimrefine_mesh">TwoPhaseFlowProblem<dim>::refine_mesh</a>
        <li><a href="#TwoPhaseFlowProblemdimoutput_results">TwoPhaseFlowProblem<dim>::output_results</a>
        <li><a href="#Toolfunctions">Tool functions</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemdimdetermine_whether_to_solve_for_pressure_and_velocity">TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity</a>
        <li><a href="#TwoPhaseFlowProblemdimproject_back_saturation">TwoPhaseFlowProblem<dim>::project_back_saturation</a>
        <li><a href="#TwoPhaseFlowProblemdimget_max_u_F_prime">TwoPhaseFlowProblem<dim>::get_max_u_F_prime</a>
        <li><a href="#TwoPhaseFlowProblemdimget_extrapolated_saturation_range">TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range</a>
        <li><a href="#TwoPhaseFlowProblemdimcompute_viscosity">TwoPhaseFlowProblem<dim>::compute_viscosity</a>
      </ul>
        <li><a href="#TwoPhaseFlowProblemdimrun">TwoPhaseFlowProblem<dim>::run</a>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program was contributed by Chih-Che Chueh (University of Victoria) and
Wolfgang Bangerth. Results from this program are used and discussed in the
following publications (in particular in the second one):
- Chih-Che Chueh, Marc Secanell, Wolfgang Bangerth, Ned Djilali. Multi-level
  adaptive simulation of transient two-phase flow in heterogeneous porous
  media. Computers &amp; Fluids, 39:1585-1596, 2010
- Chih-Che Chueh, Ned Djilali, Wolfgang Bangerth. An h-adaptive operator
  splitting method for two-phase flow in 3D heterogeneous porous
  media. SIAM Journal on Scientific Computing, 35:B149-B175, 2013.

The implementation discussed here uses and extends
parts of the step-21 and step-31 tutorial programs.

The work of the Chih-Che Chueh was funded through the Canada Research Chairs
Program and the MITACS Network of Centres of Excellence. Parts of the work by
Wolfgang Bangerth were funded through Award No. KUS-C1-016-04, made by the King
Abdullah University of Science and Technology, and through an Alfred P. Sloan
Research Fellowship.
This material is also in parts based upon work supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology; and in a continuation by the National Science
Foundation under Award No. EAR-0949446 and The University of California
&ndash; Davis. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation, The
California Institute of Technology, or of The University of California
&ndash; Davis.
</i>


<a name="Introduction"></a><a name="Intro"></a> <h1>Introduction</h1>


The simulation of multiphase flow in porous media is a ubiquitous problem, and
we have previously addressed it already in some form in step-20 and
step-21. However, as was easy to see there, it faces two major difficulties:
numerical accuracy and efficiency. The first is easy to see in the stationary
solver step-20: using lowest order Raviart-Thomas elements can not be expected
to yield highly accurate solutions. We need more accurate methods. The second
reason is apparent from the time dependent step-21: that program is
excruciatingly slow, and there is no hope to get highly accurate solutions in
3d within reasonable time frames.

In this
program, in order to overcome these two problems, there are five areas which
we are trying to improve for a high performance simulator:

<ul>
<li> Higher order spatial discretizations
<li> Adaptive mesh refinement
<li> Adaptive time stepping
<li> Operator splitting
<li> Efficient solver and preconditioning
</ul>

Much inspiration for this program comes from step-31 but several of the
techniques discussed here are original.


<a name="Advectiondominatedtwophaseflowmathematicalmodel"></a><h3>Advection-dominated two-phase flow mathematical model.</h3>


We consider the flow of a two-phase immiscible, incompressible
fluid. Capillary and gravity effects are neglected, and viscous
effects are assumed dominant. The governing equations for such a
flow that are identical to those used in step-21 and are
@f{align*}
  \mathbf{u}_t &= - \mathbf{K} \lambda_t \left(S\right) \nabla p, \\
  \nabla \cdot \mathbf{u}_t &= q, \\
  \epsilon \frac{\partial S}{\partial t} + \nabla \cdot \left( \mathbf{u}_t  F\left( S \right) \right)&=0,
@f}
where $S$ is the saturation (volume fraction between zero and one) of the second (wetting) phase, $p$ is the pressure, $\mathbf{K}$ is the permeability tensor, $\lambda_t$ is the total mobility, $\epsilon$ is the porosity, $F$ is the fractional flow of the wetting phase, $q$ is the source term and $\mathbf{u}_t$ is the total velocity. The total mobility, fractional flow of the wetting phase and total velocity are respectively given by
@f{align*}
   \lambda_t(S)&= \lambda_w + \lambda_{nw} = \frac{k_{rw}(S)}{\mu_w} + \frac{k_{rnw}(S)}{\mu_{nw}}, \\
   F(S) &= \frac{\lambda_w}{\lambda_t} = \frac{\lambda_w}{\lambda_w + \lambda_{nw}} = \frac{k_{rw}(S)/\mu_w}{k_{rw}(S)/\mu_w + k_{rnw}(S)/\mu_{nw}}, \\
   \mathbf{u}_t &= \mathbf{u}_w + \mathbf{u}_{nw} = -\lambda_t(S)\mathbf{K} \cdot \nabla p,
@f}
where subscripts $w, nw$ represent the wetting and non-wetting phases,
respectively.

For convenience, the
porosity $\epsilon$ in the saturation equation, which can be considered a
scaling factor for the time variable, is set to
one. Following a commonly used prescription for the dependence of the relative
permeabilities $k_{rw}$ and $k_{rnw}$ on saturation, we use
@f{align*}
   k_{rw}  &= S^2, \qquad&\qquad
   k_{rnw} &= \left( 1-S \right)^2.
@f}

The porous media equations above are
augmented by initial conditions for the saturation and boundary conditions for
the pressure. Since saturation and the gradient of the pressure uniquely
determine the velocity, no boundary conditions are necessary for the velocity.
Since the flow equations do not contain time derivatives, initial conditions for the velocity and pressure
variables are not required. The flow field separates the boundary into inflow or outflow
parts. Specifically,
@f[
   \mathbf{\Gamma}_{in}(t) = \left\{\mathbf{x} \in \partial \Omega:\mathbf{n} \cdot \mathbf{u}_t<0\right\},
@f]
and we arrive at a complete model by also imposing boundary values for the
saturation variable on the inflow boundary $\mathbf{\Gamma}_{in}$.


<a name="Adaptiveoperatorsplittingandtimestepping"></a><h3>Adaptive operator splitting and time stepping.</h3>


As seen in step-21, solving the flow equations for velocity and pressure are
the parts of the program that take far longer than the (explicit) updating
step for the saturation variable once we know the flow variables. On the other
hand,  the pressure and velocity depend only weakly on saturation, so one may
think about only solving for pressure and velocity every few time steps while
updating the saturation in every step. If we can find a criterion for when the
flow variables need to be updated, we call this splitting an "adaptive
operator splitting" scheme.

Here, we use the following a posteriori criterion to decide when to re-compute
pressure and velocity variables
(detailed derivations and descriptions can be found in [Chueh, Djilali
and Bangerth 2011]):
@f{align*}
  \theta(n,n_p)
  =
    \max_{\kappa\in{\mathbb T}}
    \left(
    \left\|
      \frac 1{\lambda_t\left(S^{(n-1)}\right)}
      - \frac 1{\lambda_t\left(S^{(n_p)}\right)} \right\|_{L^\infty(\kappa)}
    \left\|\|\mathbf{K}^{-1}\|_1\right\|_{L^\infty(\kappa)}
    \right).
@f}
where superscripts in parentheses denote the number of the saturation time
step at which any quantity is defined and $n_p<n$ represents the last step
where we actually computed the pressure and velocity. If $\theta(n,n_p)$
exceeds a certain threshold we re-compute the flow variables; otherwise, we
skip this computation in time step $n$ and only move the saturation variable
one time step forward.

In short, the algorithm allows us to perform a number of
saturation time steps of length $\Delta t_c^{(n)}=t^{(n)}_c-t^{(n-1)}_c$ until
the criterion above tells us to re-compute velocity and pressure
variables, leading to a macro time step of length
@f[
   \Delta t_p^{(n)} = \sum_{i=n_p+1}^{n} \Delta t_c^{(i)}.
@f]
We choose the length of (micro) steps subject to the Courant-Friedrichs-Lewy
(CFL) restriction according to the criterion
@f[
  \Delta t_c = \frac{\textrm{min}_{K}h_{K}}{7 \|\mathbf{u}_t\|_{L^{\infty}\left(\Omega\right)}},
@f]
which we have confirmed to be stable for the choice of finite element and time
stepping scheme for the saturation equation discussed below ($h_K$ denotes the
diameter of cell $K$).
The result is a scheme where neither micro nor macro time
steps are of uniform length, and both are chosen adaptively.

<a name="Timediscretization"></a><h3>Time discretization.</h3>

Using this time discretization, we obtain the following set of equations for
each time step from the IMPES approach (see step-21):
@f{align*}
   \mathbf{u}^{(n)}_t + \lambda_t\left(S^{(n-1)}\right) \mathbf{K} \nabla p^{(n)} =0, \\
   \nabla \cdot \mathbf{u}^{(n)}_t = q, \\
   \epsilon \left( \frac{S^{(n-1)}-S^{(n)}}{\Delta t^{(n)}_c} \right) + \mathbf{u}^{(n)}_t \cdot \nabla F\left(S^{(n-1)}\right) + F\left(S^{(n-1)}\right) \nabla \cdot \mathbf{u}^{(n)}_t =0.
@f}


Using the fact that $\nabla \cdot \mathbf{u}_t = q$, the time discrete
saturation equation becomes
@f{align*}
  &\epsilon \left( \frac{S^{(n)}-S^{(n-1)}}{\Delta t^{(n)}_c} \right) + \mathbf{u}^{(n)}_t \cdot \nabla F\left(S^{(n-1)}\right) + F\left(S^{(n-1)}\right)q=0.
@f}

<a name="Weakformspacediscretizationforthepressurevelocitypart"></a><h3>Weak form, space discretization for the pressure-velocity part.</h3>


By multiplying the equations defining the total velocity $\mathbf u_t^{(n)}$ and
the equation that expresses its divergence in terms of source terms, with test
functions $\mathbf{v}$ and $w$
respectively and then integrating terms by parts as necessary, the weak form
of the problem reads: Find $\mathbf u, p$ so that for all test functions
$\mathbf{v}, w$ there holds
@f{gather*}
   \left( \left( \mathbf{K} \lambda_t\left(S^{(n-1)}\right) \right)^{-1} \mathbf{u}^{(n)}_t, \mathbf{v}\right)_{\Omega} - \left(p^{(n)}, \nabla \cdot \mathbf{v}\right)_{\Omega} = -\left(p^{(n)}, \mathbf{n} \cdot \mathbf{v} \right)_{\partial \Omega}, \\
   - \left( \nabla \cdot \mathbf{u}^{(n)}_t,w\right)_{\Omega} = - \big(q,w\big)_{\Omega}.
@f}
Here, $\mathbf{n}$ represents the unit outward normal vector to $\partial
\Omega$ and the pressure $p^{(n)}$ can be prescribed weakly on the open part
of the boundary $\partial \Omega$ whereas on those parts where a velocity is
prescribed (for example impermeable boundaries with $\mathbf n \cdot \mathbf
u=0$ the term disappears altogether because $\mathbf n \cdot \mathbf
v=0$.

We use continuous finite elements to discretize the velocity and pressure
equations. Specifically, we use mixed finite elements to ensure high order approximation
for both vector (e.g. a fluid velocity) and scalar variables (e.g. pressure)
simultaneously. For saddle point problems, it is well established that
the so-called Babuska-Brezzi or Ladyzhenskaya-Babuska-Brezzi (LBB) conditions
[Brezzi 1991, Chen 2005] need to be satisfied to ensure stability of
the pressure-velocity system. These stability conditions are satisfied in the
present work by using elements for velocity that are one order higher than for
the pressure, i.e. $u_h \in Q^d_{p+1}$ and $p_h \in Q_p$, where $p=1$, $d$ is
the space dimension, and $Q_s$ denotes the space of tensor product Lagrange
polynomials of degree $s$ in each variable.

<a name="Stabilizationweakformandspacediscretizationforthesaturationtransportequation"></a><h3>Stabilization, weak form and space discretization for the saturation transport equation.</h3>

The chosen $Q_1$ elements for the saturation equation do not lead to a stable
discretization without upwinding or other kinds of stabilization, and spurious
oscillations will appear in the numerical solution. Adding an artificial
diffusion term is one approach to eliminating these oscillations
[Chen 2005]. On the other hand, adding too much diffusion smears sharp
fronts in the solution and suffers from grid-orientation difficulties
[Chen 2005]. To avoid these effects, we use the artificial diffusion
term proposed by [Guermond and Pasquetti 2008] and
validated in [Chueh, Djilali, Bangerth 2011] and
[Kronbichler, Heister and Bangerth, 2011], as well as in step-31.

This method modifies the (discrete) weak form of the saturation equation
to read
@f{align*}
  \left(\epsilon \frac{\partial S_h}{\partial t},\sigma_h\right)
  -
  \left(\mathbf{u}_t  F\left( S_h \right),
    \nabla \sigma_h\right)
  +
  \left(\mathbf n \cdot \mathbf{u}_t  \hat F\left( S_h \right),
    \sigma_h\right)_{\partial\Omega}
  +
  (\nu(S_h) \nabla S_h, \nabla \sigma_h)
  &=0
  \qquad
  \forall \sigma_h,
@f}
where $\nu$ is the artificial diffusion parameter and $\hat F$ is an
appropriately chosen numerical flux on the boundary of the domain (we choose
the obvious full upwind flux for this).

Following [Guermond and Pasquetti 2008] (and as detailed in
[Chueh, Djilali and Bangerth 2011]), we use
the parameter as a piecewise
constant function set on each cell $K$ with the diameter $h_{K}$ as
@f[
   \nu(S_h)|_{K} = \beta \| \mathbf{u}_t \max\{F'(S_h),1\} \|_{L^{\infty}(K)} \textrm{min} \left\{ h_{K},h^{\alpha}_{K} \frac{\|\textrm{Res}(S_h)\|_{L^{\infty}(K)}}{c(\mathbf{u}_t,S)} \right\}
@f]
where $\alpha$ is a stabilization exponent and $\beta$ is a dimensionless
user-defined stabilization constant. Following [Guermond and Pasquetti 2008]
as well as the implementation in step-31, the velocity and saturation global
normalization constant, $c(\mathbf{u}_t,S)$, and the residual $\textrm{Res}(S)$
are respectively given by
@f[
   c(\mathbf{u}_t,S) = c_R \|\mathbf{u}_t \max\{F'(S),1\}\|_{L^{\infty}(\Omega)} \textrm{var}(S)^\alpha | \textrm{diam} (\Omega) |^{\alpha - 2}
@f]
and
@f[
   \textrm{Res}(S) = \left( \epsilon \frac{\partial S}{\partial t} + \mathbf{u}_t \cdot \nabla F(S) + F(S)q \right) \cdot S^{\alpha - 1}
@f]
where $c_R$ is a second dimensionless user-defined constant,
$\textrm{diam}(\Omega)$ is the diameter of the domain and $\textrm{var}(S) =
\textrm{max}_{\Omega} S - \textrm{min}_{\Omega} S$ is the range of the present
saturation values in the entire computational domain $\Omega$.

This stabilization scheme has a number of advantages over simpler schemes such
as finite volume (or discontinuous Galerkin) methods or streamline upwind
Petrov Galerkin (SUPG) discretizations. In particular, the artificial
diffusion term acts primarily in the vicinity of discontinuities
since the residual is small in areas where the saturation is smooth. It
therefore provides for a higher degree of accuracy. On the other hand, it is
nonlinear since $\nu$ depends on the saturation $S$. We avoid this difficulty
by treating all nonlinear terms explicitly, which leads to the following
fully discrete problem at time step $n$:
@f{align*}
   &\left( \epsilon S_h^{(n)},\sigma_h\right)_{\Omega} - \Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\mathbf{u}^{*}_t,\nabla\sigma_h\Big)_{\Omega} + \Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{*}_t\right),\sigma_h\Big)_{\partial\Omega} \nonumber \\
   & \quad = \left( \epsilon S_h^{(n-1)},\sigma_h\right)_{\Omega} - \Delta t^{(n)}_c \bigg(\nu\left(S_h^{(n-1)}\right)\nabla S_h^{(n-1)},\nabla\sigma_h\bigg)_{\Omega} \nonumber \\
   & \qquad + \Delta t^{(n)}_c \bigg(\mathbf{n}\cdot\nu\left(S_h^{(n-1)}\right)\nabla S^{(n-1)},\sigma_h\bigg)_{\partial\Omega}
@f}
where $\mathbf{u}_t^{*}$ is the velocity linearly extrapolated from
$\mathbf{u}^{(n_p)}_t$ and $\mathbf{u}^{(n_{pp})}_t$ to the current time $t^{(n)}$ if $\theta<\theta^*$ while $\mathbf{u}_t^{*}$ is $\mathbf{u}^{(n_p)}_t$ if $\theta>\theta^*$.
Consequently, the equation is linear in $S_h^{(n)}$ and all that is required
is to solve with a mass matrix on the saturation space.

Since the Dirichlet boundary conditions for saturation are only imposed on the
inflow boundaries, the third term on the left hand side of the equation above
needs to be split further into two parts:
@f{align*}
  &\Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_t\right),\sigma_h\Big)_{\partial\Omega} \nonumber \\
  &\qquad= \Delta t^{(n)}_c \Big(F\left(S^{(n-1)}_{(+)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_{t(+)}\right),\sigma_h\Big)_{\partial\Omega_{(+)}} + \Delta t^{(n)}_c \Big(F\left(S^{(n-1)}_{(-)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_{t(-)}\right),\sigma_h\Big)_{\partial\Omega_{(-)}}
@f}
where $\partial\Omega_{(-)} = \left\{\mathbf{x} \in \partial\Omega : \mathbf{n}
  \cdot \mathbf{u}_t<0\right\}$ and
$\partial\Omega_{(+)} = \left\{\mathbf{x} \in \partial\Omega : \mathbf{n} \cdot
  \mathbf{u}_t>0\right\}$ represent inflow and outflow boundaries,
respectively. We choose values using an
upwind formulation, i.e. $S^{(n-1)}_{(+)}$ and $\mathbf{u}^{(n)}_{t(+)}$
correspond to the values taken from the present cell, while the values of
$S^{(n-1)}_{(-)}$ and $\mathbf{u}^{(n)}_{t(-)}$ are those taken from the
neighboring boundary $\partial\Omega_{(-)}$.


<a name="Adaptivemeshrefinement"></a><h3>Adaptive mesh refinement.</h3>


Choosing meshes adaptively to resolve sharp
saturation fronts is an essential ingredient to achieve efficiency in our
algorithm. Here, we use the same shock-type refinement approach used in
[Chueh, Djilali and Bangerth 2011] to select those cells that should be refined or
coarsened. The refinement indicator for each cell $K$ of the triangulation is
computed by
@f[
   \eta_{K} = |\nabla S_h(\mathbf x_K)|
@f]
where $\nabla S_h(\mathbf x_K)$ is the gradient of the discrete saturation
variable evaluated at the center $\mathbf x_K$ of cell $K$. This approach is
analogous to ones frequently used in compressible flow problems, where density
gradients are used to indicate refinement. That said, as we will
discuss at the end of the <a href="#Results">results section</a>, this turns
out to not be a very useful criterion since it leads to refinement basically
everywhere. We only show it here for illustrative purposes.


<a name="Linearsystemanditspreconditioning"></a><h3>Linear system and its preconditioning.</h3>


Following the discretization of the governing equations
discussed above, we
obtain a linear system of equations in time step $(n)$ of the following form:
@f[
 \left(
  \begin{array}{ccc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T} & \mathbf{0}  \\
   \mathbf{B}           & \mathbf{0}     & \mathbf{0}   \\
   \mathbf{H}           & \mathbf{0}     & \mathbf{M}^{S}
  \end{array}
 \right)
 \left(
  \begin{array}{c}
   \mathbf{U}^{(n)} \\
   \mathbf{P}^{(n)} \\
   \mathbf{S}^{(n)}
  \end{array}
 \right)
 =
 \left(
  \begin{array}{c}
   0 \\
   \mathbf{F}_{2} \\
   \mathbf{F}_{3}
  \end{array}
 \right)
@f]
where the individual matrices and vectors are defined as follows using shape functions $\mathbf{v}_i$ for velocity, and $\phi_i$ for both pressure and saturation:
@f{align*}
  \mathbf{M}^{\mathbf{u}}_{ij}
  &= \left( \left( \mathbf{K} \lambda_t\left(S^{(n-1)}\right) \right)^{-1}
  \mathbf{v}_{i},\mathbf{v}_{j}\right)_{\Omega},
  &
  \mathbf{M}^{S}_{ij}           &= \left(\epsilon \phi_i,\phi_j\right)_{\Omega}
  \\
  \mathbf{B}_{ij}
  &= - \left( \nabla \cdot \mathbf{v}_{j},\phi_{i}\right)_{\Omega},
  &
  \mathbf{H}_{ij}
  &= - \Delta t^{(n)}_c \Big( F\left(S^{(n-1)}\right) \mathbf{v}_i,\nabla\phi_j\Big)_{\Omega}
  \\
  \left(\mathbf{F}_{2}\right)_i
  &= - \big(F\left(S^{(n-1)}\right)q,\phi_i\big)_{\Omega},
@f}
and $\mathbf{F}_{3}$ as given in the definition of the stabilized transport
equation.

The linear system above is of block triangular form if we consider the top
left $2\times 2$ panel of matrices as one block. We can therefore first solve
for the velocity and pressure (unless we decide to use $\mathbf U^{(n_p)}$ in
place of the velocity)
followed by a solve for the saturation variable. The first of these steps
requires us to solve
@f[
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T}  \\
   \mathbf{B}           & \mathbf{0}
  \end{array}
 \right)
 \left(
  \begin{array}{c}
   \mathbf{U}^{(n)} \\
   \mathbf{P}^{(n)}
  \end{array}
 \right)
 =
 \left(
  \begin{array}{c}
   0 \\
   \mathbf{F}_{2}
  \end{array}
 \right)
@f]
We apply the Generalized Minimal Residual (GMRES) method [Saad and Schultz
1986] to this linear system. The ideal preconditioner for the
velocity-pressure system is
@f{align*}
\mathbf{P} =
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} &  \mathbf{0}  \\
   \mathbf{B}           & -\mathbf{S}
  \end{array}
 \right),
 & \qquad
 \mathbf{P}^{-1} =
 \left(
  \begin{array}{cc}
   \left(\mathbf{M}^{\mathbf{u}}\right)^{-1}                              &  \mathbf{0}  \\
   \mathbf{S}^{-1} \mathbf{B} \left(\mathbf{M}^{\mathbf{u}}\right)^{-1}   & -\mathbf{S}^{-1}
  \end{array}
 \right)
 @f}
where
$\mathbf{S}=\mathbf{B}\left(\mathbf{M}^{\mathbf{u}}\right)^{-1}\mathbf{B}^T$ is
the Schur complement [Zhang 2005] of the system. This preconditioner is
optimal since
@f{align*}
 \mathbf{P}^{-1}
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T}  \\
   \mathbf{B}           & \mathbf{0}
  \end{array}
 \right)
 =
  \left(
  \begin{array}{cc}
   \mathbf{I}         &  \left(\mathbf{M}^{\mathbf{u}}\right)^{-1} \mathbf{B}^{T}  \\
   \mathbf{0}         &  \mathbf{I}
  \end{array}
 \right),
@f}
for which it can be shown that GMRES converges in two iterations.

However, we cannot of course expect to use exact inverses of the
velocity mass matrix and the Schur complement. We therefore follow the
approach by [Silvester and Wathen 1994] originally proposed for
the Stokes system. Adapting it to the current set of equations yield the
preconditioner
@f{align*}
 \mathbf{\tilde{P}}^{-1} =
 \left(
  \begin{array}{cc}
   \widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}
                              &  \mathbf{0}  \\
   \widetilde{\mathbf{{S}}^{-1}} \mathbf{B} \widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}   & -\widetilde{\mathbf{{S}}^{-1}}
  \end{array}
 \right)
@f}
where a tilde indicates an approximation of the exact inverse matrix. In
particular, since $\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}=\left( \left(
    \mathbf{K} \lambda_t \right)^{-1}
  \mathbf{v}_{i},\mathbf{v}_{j}\right)_{\Omega}$
is a sparse symmetric and positive definite matrix, we choose for
$\widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}$ a single application of
a sparse incomplete Cholesky decomposition of this matrix
[Golub and Van Loan 1996].
We note that the Schur complement that corresponds to the porous
media flow operator in non-mixed form, $-\nabla \cdot [\mathbf K
\lambda_t(S)]\nabla$ and
$\mathbf{\tilde {S}} = \left( \left( \mathbf{K} \lambda_t \right) \nabla \phi_{i},\nabla \phi_{j}\right)_{\Omega}$
should be a good approximation of the actual Schur complement matrix $\mathbf
S$. Since both of these matrices are again symmetric and positive definite, we
use an incomplete Cholesky decomposition of $\mathbf{\tilde S}$ for $\widetilde
{\mathbf{{S}}^{-1}}$. It is important to note that $\mathbf{\tilde S}$ needs
to be built with Dirichlet boundary conditions to ensure its invertibility.

Once the velocity $\mathbf{U}^{(n)} \equiv \mathbf{u}^*_t$  is available, we
can assemble $\mathbf{H}$ and
$\mathbf{F}_{3}$ and solve for the saturations using
@f{align*}
  \mathbf{M}^{S} \mathbf{S}^{(n)} = \mathbf{F}_{3} - \mathbf{H} \mathbf{U}^{(n)}.
@f}
where the mass matrix $\mathbf{M}^{S}$ is solved by the conjugate gradient
method, using an incomplete Cholesky decomposition as preconditioner once
more.

<a name="Thetestcases"></a><h3>The test cases.</h3>


@note
The implementation discussed here uses and extends
parts of the step-21, step-31 and step-33 tutorial programs of this
library. In particular, if you want to understand how it works, please
consult step-21 for a discussion of the mathematical problem, and
step-31 from which most of the implementation is derived. We will not
discuss aspects of the implementation that have already been discussed
in step-31.

We show numerical results for some two-phase flow equations augmented by
appropriate initial and boundary conditions in conjunction with two different
choices of the permeability model. In the problems considered, there is no
internal source term ($q=0$). As mentioned above, quantitative numerical
results are presented in [Chueh, Djilali and Bangerth 2011].

For simplicity, we choose $\Omega=[0,1]^d,d=2,3$, though all methods (as well
as our implementation) should work equally well on general unstructured meshes.

Initial conditions are only required for the saturation variable, and we
choose $S(\mathbf{x},0)=0.2$, i.e. the porous medium is initially filled by a
mixture of the non-wetting (80%) and wetting (20%) phases. This differs from
the initial condition in step-21 where we had taken $S(\mathbf{x},0)=0$, but
for complicated mathematical reasons that are mentioned there in a longish
remark, the current method using an entropy-based artificial diffusion term
does not converge to the viscosity solution with this initial condition
without additional modifications to the method. We therefore choose this
modified version for the current program.

Furthermore, we prescribe a linear pressure on
the boundaries:
@f[
   p(\mathbf{x},t) = 1 - x \qquad
   \textrm{on} \quad \partial \Omega \times [0,T].
@f]
Pressure and saturation uniquely
determine a velocity, and the velocity determines whether a boundary segment
is an inflow or outflow boundary. On the inflow part of the boundary,
$\mathbf{\Gamma}_{in}(t)$, we impose
@f{align*}
   S(\mathbf{x},t) = 1 \qquad & \textrm{on} \quad \mathbf{\Gamma}_{in}(t) \cap \left\{x = 0\right\}, \\
   S(\mathbf{x},t) = 0 \qquad & \textrm{on} \quad \mathbf{\Gamma}_{in}(t) \backslash \left\{x = 0\right\}.
@f}
In other words, the domain is flooded by the wetting phase from the left.
No boundary conditions for the saturation are required for the outflow parts
of the boundary.

All the numerical and physical parameters used for the 2D/3D
cases are listed in the following table:

<table align="center" class="tutorial" width="50%">
<tr>
    <th>Parameter                           </th><th>Symbol          </th><th>Value               </th><th>units     </th></tr><tr>
    <td>Porosity                            </td><td>$\epsilon$      </td><td>1.0                 </td><td>-                   </td></tr><tr>
    <td>Viscosity (wetting)                 </td><td>$\mu_w$         </td><td>0.2                 </td><td>$kg \cdot m^{-1} \cdot sec^{-1}$   </td></tr><tr>
    <td>Viscosity (nonwetting)              </td><td>$\mu_{nw}$      </td><td>1.0                 </td><td>$kg \cdot m^{-1} \cdot sec^{-1}$      </td></tr><tr>
    <td>Stabilization exponent              </td><td>$\alpha$        </td><td>1.0                 </td><td>-     </td></tr><tr>
    <td>Stabilization constant              </td><td>$\beta$         </td><td>2D: 0.3; 3D: 0.27   </td><td>- </td></tr><tr>
    <td>Normalization constant              </td><td>$c_R$           </td><td>1.0                 </td><td>- </td></tr><tr>
    <td>Number of high-permeability regions </td><td>$N$             </td><td>50; 200             </td><td>- </td></tr><tr>
    <td>Operator splitting threshold        </td><td>$\theta^\ast$   </td><td>5.0              </td><td>- </td></tr>
</table>


<a name="Listofreferences"></a><h3>List of references</h3>



<ol>
<li>
CC Chueh, N Djilali and W Bangerth.
<br> An h-adaptive operator splitting method for two-phase flow in 3D
  heterogeneous porous media.
<br> SIAM Journal on Scientific Computing, vol. 35 (2013), pp. B149-B175

<li>
M. Kronbichler, T. Heister, and W. Bangerth
<br> High Accuracy Mantle Convection Simulation through Modern Numerical
Methods.
<br> Geophysics Journal International, vol. 191 (2012), pp. 12-29

<li>
F Brezzi and M Fortin.
<br> <i>Mixed and Hybrid Finite Element Methods</i>.
<br> Springer-Verlag, 1991.

<li>
Z Chen.
<br> <i>Finite Element Methods and Their Applications</i>.
<br> Springer, 2005.

<li>
JL Guermond and R Pasquetti.
<br> Entropy-based nonlinear viscosity for Fourier approximations of
  conservation laws.
<br> <i>Comptes Rendus Mathematique</i>, 346(13-14):801-806, 2008.

<li>
CC Chueh, M Secanell, W Bangerth, and N Djilali.
<br> Multi-level adaptive simulation of transient two-phase flow in
  heterogeneous porous media.
<br> <i>Computers and Fluids</i>, 39:1585-1596, 2010.

<li>
Y Saad and MH Schultz.
<br> Gmres: A generalized minimal residual algorithm for solving
  nonsymmetric linear systems.
<br> <i>SIAM Journal on Scientific and Statistical Computing</i>,
  7(3):856-869, 1986.

<li>
F Zhang.
<br> <i>The Schur Complement and its Applications</i>.
<br> Springer, 2005.

<li>
D Silvester and A Wathen.
<br> Fast iterative solution of stabilised Stokes systems part ii: Using
  general block preconditioners.
<br> <i>SIAM Journal on Numerical Analysis</i>, 31(5):1352-1367, 1994.

<li>
GH Golub and CF van Loan.
<br> <i>Matrix Computations</i>.
<br> 3rd Edition, Johns Hopkins, 1996.

<li>
SE Buckley and MC Leverett.
<br> Mechanism of fluid displacements in sands.
<br> <i>AIME Trans.</i>, 146:107-116, 1942.

</ol>
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * The first step, as always, is to include the functionality of a number of
 * deal.II and C++ header files.
 * 

 * 
 * The list includes some header files that provide vector, matrix, and
 * preconditioner classes that implement interfaces to the respective Trilinos
 * classes; some more information on these may be found in step-31.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/tensor_function.h>
 * #include <deal.II/base/index_set.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/block_sparsity_pattern.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_vector.h>
 * #include <deal.II/lac/trilinos_parallel_block_vector.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * #include <memory>
 * 
 * 
 * @endcode
 * 
 * At the end of this top-matter, we open a namespace for the current project
 * into which all the following material will go, and then import all deal.II
 * names into this namespace:
 * 
 * @code
 * namespace Step43
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Boundaryandinitialvalueclasses"></a> 
 * <h3>Boundary and initial value classes</h3>
 * 

 * 
 * The following part is taken directly from step-21 so there is no need to
 * repeat the descriptions found there.
 * 
 * @code
 *   template <int dim>
 *   class PressureBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     PressureBoundaryValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double
 *   PressureBoundaryValues<dim>::value(const Point<dim> &p,
 *                                      const unsigned int /*component*/) const
 *   {
 *     return 1 - p[0];
 *   }
 * 
 * 
 *   template <int dim>
 *   class SaturationBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     SaturationBoundaryValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double
 *   SaturationBoundaryValues<dim>::value(const Point<dim> &p,
 *                                        const unsigned int /*component*/) const
 *   {
 *     if (p[0] == 0)
 *       return 1;
 *     else
 *       return 0;
 *   }
 * 
 * 
 *   template <int dim>
 *   class SaturationInitialValues : public Function<dim>
 *   {
 *   public:
 *     SaturationInitialValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  value) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double
 *   SaturationInitialValues<dim>::value(const Point<dim> & /*p*/,
 *                                       const unsigned int /*component*/) const
 *   {
 *     return 0.2;
 *   }
 * 
 * 
 *   template <int dim>
 *   void SaturationInitialValues<dim>::vector_value(const Point<dim> &p,
 *                                                   Vector<double> &values) const
 *   {
 *     for (unsigned int c = 0; c < this->n_components; ++c)
 *       values(c) = SaturationInitialValues<dim>::value(p, c);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Permeabilitymodels"></a> 
 * <h3>Permeability models</h3>
 * 

 * 
 * In this tutorial, we still use the two permeability models previously
 * used in step-21 so we again refrain from commenting in detail about them.
 * 
 * @code
 *   namespace SingleCurvingCrack
 *   {
 *     template <int dim>
 *     class KInverse : public TensorFunction<2, dim>
 *     {
 *     public:
 *       KInverse()
 *         : TensorFunction<2, dim>()
 *       {}
 * 
 *       virtual void
 *       value_list(const std::vector<Point<dim>> &points,
 *                  std::vector<Tensor<2, dim>> &  values) const override;
 *     };
 * 
 * 
 *     template <int dim>
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                    std::vector<Tensor<2, dim>> &  values) const
 *     {
 *       Assert(points.size() == values.size(),
 *              ExcDimensionMismatch(points.size(), values.size()));
 * 
 *       for (unsigned int p = 0; p < points.size(); ++p)
 *         {
 *           values[p].clear();
 * 
 *           const double distance_to_flowline =
 *             std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0]));
 * 
 *           const double permeability =
 *             std::max(std::exp(-(distance_to_flowline * distance_to_flowline) /
 *                               (0.1 * 0.1)),
 *                      0.01);
 * 
 *           for (unsigned int d = 0; d < dim; ++d)
 *             values[p][d][d] = 1. / permeability;
 *         }
 *     }
 *   } // namespace SingleCurvingCrack
 * 
 * 
 *   namespace RandomMedium
 *   {
 *     template <int dim>
 *     class KInverse : public TensorFunction<2, dim>
 *     {
 *     public:
 *       KInverse()
 *         : TensorFunction<2, dim>()
 *       {}
 * 
 *       virtual void
 *       value_list(const std::vector<Point<dim>> &points,
 *                  std::vector<Tensor<2, dim>> &  values) const override;
 * 
 *     private:
 *       static std::vector<Point<dim>> centers;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     std::vector<Point<dim>> KInverse<dim>::centers = []() {
 *       const unsigned int N =
 *         (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented()));
 * 
 *       std::vector<Point<dim>> centers_list(N);
 *       for (unsigned int i = 0; i < N; ++i)
 *         for (unsigned int d = 0; d < dim; ++d)
 *           centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX;
 * 
 *       return centers_list;
 *     }();
 * 
 * 
 * 
 *     template <int dim>
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                    std::vector<Tensor<2, dim>> &  values) const
 *     {
 *       AssertDimension(points.size(), values.size());
 * 
 *       for (unsigned int p = 0; p < points.size(); ++p)
 *         {
 *           values[p].clear();
 * 
 *           double permeability = 0;
 *           for (unsigned int i = 0; i < centers.size(); ++i)
 *             permeability +=
 *               std::exp(-(points[p] - centers[i]).norm_square() / (0.05 * 0.05));
 * 
 *           const double normalized_permeability =
 *             std::min(std::max(permeability, 0.01), 4.);
 * 
 *           for (unsigned int d = 0; d < dim; ++d)
 *             values[p][d][d] = 1. / normalized_permeability;
 *         }
 *     }
 *   } // namespace RandomMedium
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Physicalquantities"></a> 
 * <h3>Physical quantities</h3>
 * 

 * 
 * The implementations of all the physical quantities such as total mobility
 * $\lambda_t$ and fractional flow of water $F$ are taken from step-21 so
 * again we don't have do any comment about them. Compared to step-21 we
 * have added checks that the saturation passed to these functions is in
 * fact within the physically valid range. Furthermore, given that the
 * wetting phase moves at speed $\mathbf u F'(S)$ it is clear that $F'(S)$
 * must be greater or equal to zero, so we assert that as well to make sure
 * that our calculations to get at the formula for the derivative made
 * sense.
 * 
 * @code
 *   double mobility_inverse(const double S, const double viscosity)
 *   {
 *     return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S));
 *   }
 * 
 * 
 *   double fractional_flow(const double S, const double viscosity)
 *   {
 *     Assert((S >= 0) && (S <= 1),
 *            ExcMessage("Saturation is outside its physically valid range."));
 * 
 *     return S * S / (S * S + viscosity * (1 - S) * (1 - S));
 *   }
 * 
 * 
 *   double fractional_flow_derivative(const double S, const double viscosity)
 *   {
 *     Assert((S >= 0) && (S <= 1),
 *            ExcMessage("Saturation is outside its physically valid range."));
 * 
 *     const double temp = (S * S + viscosity * (1 - S) * (1 - S));
 * 
 *     const double numerator =
 *       2.0 * S * temp - S * S * (2.0 * S - 2.0 * viscosity * (1 - S));
 *     const double denominator = std::pow(temp, 2.0);
 * 
 *     const double F_prime = numerator / denominator;
 * 
 *     Assert(F_prime >= 0, ExcInternalError());
 * 
 *     return F_prime;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Helperclassesforsolversandpreconditioners"></a> 
 * <h3>Helper classes for solvers and preconditioners</h3>
 * 

 * 
 * In this first part we define a number of classes that we need in the
 * construction of linear solvers and preconditioners. This part is
 * essentially the same as that used in step-31. The only difference is that
 * the original variable name stokes_matrix is replaced by another name
 * darcy_matrix to match our problem.
 * 
 * @code
 *   namespace LinearSolvers
 *   {
 *     template <class MatrixType, class PreconditionerType>
 *     class InverseMatrix : public Subscriptor
 *     {
 *     public:
 *       InverseMatrix(const MatrixType &        m,
 *                     const PreconditionerType &preconditioner);
 * 
 * 
 *       template <typename VectorType>
 *       void vmult(VectorType &dst, const VectorType &src) const;
 * 
 *     private:
 *       const SmartPointer<const MatrixType> matrix;
 *       const PreconditionerType &           preconditioner;
 *     };
 * 
 * 
 *     template <class MatrixType, class PreconditionerType>
 *     InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(
 *       const MatrixType &        m,
 *       const PreconditionerType &preconditioner)
 *       : matrix(&m)
 *       , preconditioner(preconditioner)
 *     {}
 * 
 * 
 * 
 *     template <class MatrixType, class PreconditionerType>
 *     template <typename VectorType>
 *     void InverseMatrix<MatrixType, PreconditionerType>::vmult(
 *       VectorType &      dst,
 *       const VectorType &src) const
 *     {
 *       SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm());
 *       SolverCG<VectorType> cg(solver_control);
 * 
 *       dst = 0;
 * 
 *       try
 *         {
 *           cg.solve(*matrix, dst, src, preconditioner);
 *         }
 *       catch (std::exception &e)
 *         {
 *           Assert(false, ExcMessage(e.what()));
 *         }
 *     }
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     class BlockSchurPreconditioner : public Subscriptor
 *     {
 *     public:
 *       BlockSchurPreconditioner(
 *         const TrilinosWrappers::BlockSparseMatrix &S,
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                             PreconditionerTypeMp> &Mpinv,
 *         const PreconditionerTypeA &                Apreconditioner);
 * 
 *       void vmult(TrilinosWrappers::MPI::BlockVector &      dst,
 *                  const TrilinosWrappers::MPI::BlockVector &src) const;
 * 
 *     private:
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix>
 *         darcy_matrix;
 *       const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                                              PreconditionerTypeMp>>
 *                                  m_inverse;
 *       const PreconditionerTypeA &a_preconditioner;
 * 
 *       mutable TrilinosWrappers::MPI::Vector tmp;
 *     };
 * 
 * 
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::
 *       BlockSchurPreconditioner(
 *         const TrilinosWrappers::BlockSparseMatrix &S,
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                             PreconditionerTypeMp> &Mpinv,
 *         const PreconditionerTypeA &                Apreconditioner)
 *       : darcy_matrix(&S)
 *       , m_inverse(&Mpinv)
 *       , a_preconditioner(Apreconditioner)
 *       , tmp(complete_index_set(darcy_matrix->block(1, 1).m()))
 *     {}
 * 
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     void
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult(
 *       TrilinosWrappers::MPI::BlockVector &      dst,
 *       const TrilinosWrappers::MPI::BlockVector &src) const
 *     {
 *       a_preconditioner.vmult(dst.block(0), src.block(0));
 *       darcy_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1));
 *       tmp *= -1;
 *       m_inverse->vmult(dst.block(1), tmp);
 *     }
 *   } // namespace LinearSolvers
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheTwoPhaseFlowProblemclass"></a> 
 * <h3>The TwoPhaseFlowProblem class</h3>
 * 

 * 
 * The definition of the class that defines the top-level logic of solving
 * the time-dependent advection-dominated two-phase flow problem (or
 * Buckley-Leverett problem [Buckley 1942]) is mainly based on tutorial
 * programs step-21 and step-33, and in particular on step-31 where we have
 * used basically the same general structure as done here. As in step-31,
 * the key routines to look for in the implementation below are the
 * <code>run()</code> and <code>solve()</code> functions.
 *   

 * 
 * The main difference to step-31 is that, since adaptive operator splitting
 * is considered, we need a couple more member variables to hold the last
 * two computed Darcy (velocity/pressure) solutions in addition to the
 * current one (which is either computed directly, or extrapolated from the
 * previous two), and we need to remember the last two times we computed the
 * Darcy solution. We also need a helper function that figures out whether
 * we do indeed need to recompute the Darcy solution.
 *   

 * 
 * Unlike step-31, this step uses one more AffineConstraints object called
 * darcy_preconditioner_constraints. This constraint object is used only for
 * assembling the matrix for the Darcy preconditioner and includes hanging
 * node constraints as well as Dirichlet boundary value constraints for the
 * pressure variable. We need this because we are building a Laplace matrix
 * for the pressure as an approximation of the Schur complement) which is
 * only positive definite if boundary conditions are applied.
 *   

 * 
 * The collection of member functions and variables thus declared in this
 * class is then rather similar to those in step-31:
 * 
 * @code
 *   template <int dim>
 *   class TwoPhaseFlowProblem
 *   {
 *   public:
 *     TwoPhaseFlowProblem(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void setup_dofs();
 *     void assemble_darcy_preconditioner();
 *     void build_darcy_preconditioner();
 *     void assemble_darcy_system();
 *     void assemble_saturation_system();
 *     void assemble_saturation_matrix();
 *     void assemble_saturation_rhs();
 *     void assemble_saturation_rhs_cell_term(
 *       const FEValues<dim> &                       saturation_fe_values,
 *       const FEValues<dim> &                       darcy_fe_values,
 *       const double                                global_max_u_F_prime,
 *       const double                                global_S_variation,
 *       const std::vector<types::global_dof_index> &local_dof_indices);
 *     void assemble_saturation_rhs_boundary_term(
 *       const FEFaceValues<dim> &                   saturation_fe_face_values,
 *       const FEFaceValues<dim> &                   darcy_fe_face_values,
 *       const std::vector<types::global_dof_index> &local_dof_indices);
 *     void solve();
 *     void refine_mesh(const unsigned int min_grid_level,
 *                      const unsigned int max_grid_level);
 *     void output_results() const;
 * 
 * @endcode
 * 
 * We follow with a number of helper functions that are used in a variety
 * of places throughout the program:
 * 
 * @code
 *     double                    get_max_u_F_prime() const;
 *     std::pair<double, double> get_extrapolated_saturation_range() const;
 *     bool   determine_whether_to_solve_for_pressure_and_velocity() const;
 *     void   project_back_saturation();
 *     double compute_viscosity(
 *       const std::vector<double> &        old_saturation,
 *       const std::vector<double> &        old_old_saturation,
 *       const std::vector<Tensor<1, dim>> &old_saturation_grads,
 *       const std::vector<Tensor<1, dim>> &old_old_saturation_grads,
 *       const std::vector<Vector<double>> &present_darcy_values,
 *       const double                       global_max_u_F_prime,
 *       const double                       global_S_variation,
 *       const double                       cell_diameter) const;
 * 
 * 
 * @endcode
 * 
 * This all is followed by the member variables, most of which are similar
 * to the ones in step-31, with the exception of the ones that pertain to
 * the macro time stepping for the velocity/pressure system:
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     double             global_Omega_diameter;
 * 
 *     const unsigned int degree;
 * 
 *     const unsigned int        darcy_degree;
 *     FESystem<dim>             darcy_fe;
 *     DoFHandler<dim>           darcy_dof_handler;
 *     AffineConstraints<double> darcy_constraints;
 * 
 *     AffineConstraints<double> darcy_preconditioner_constraints;
 * 
 *     TrilinosWrappers::BlockSparseMatrix darcy_matrix;
 *     TrilinosWrappers::BlockSparseMatrix darcy_preconditioner_matrix;
 * 
 *     TrilinosWrappers::MPI::BlockVector darcy_solution;
 *     TrilinosWrappers::MPI::BlockVector darcy_rhs;
 * 
 *     TrilinosWrappers::MPI::BlockVector last_computed_darcy_solution;
 *     TrilinosWrappers::MPI::BlockVector second_last_computed_darcy_solution;
 * 
 * 
 *     const unsigned int        saturation_degree;
 *     FE_Q<dim>                 saturation_fe;
 *     DoFHandler<dim>           saturation_dof_handler;
 *     AffineConstraints<double> saturation_constraints;
 * 
 *     TrilinosWrappers::SparseMatrix saturation_matrix;
 * 
 * 
 *     TrilinosWrappers::MPI::Vector saturation_solution;
 *     TrilinosWrappers::MPI::Vector old_saturation_solution;
 *     TrilinosWrappers::MPI::Vector old_old_saturation_solution;
 *     TrilinosWrappers::MPI::Vector saturation_rhs;
 * 
 *     TrilinosWrappers::MPI::Vector
 *       saturation_matching_last_computed_darcy_solution;
 * 
 *     const double saturation_refinement_threshold;
 * 
 *     double       time;
 *     const double end_time;
 * 
 *     double current_macro_time_step;
 *     double old_macro_time_step;
 * 
 *     double       time_step;
 *     double       old_time_step;
 *     unsigned int timestep_number;
 * 
 *     const double viscosity;
 *     const double porosity;
 *     const double AOS_threshold;
 * 
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC> Amg_preconditioner;
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC> Mp_preconditioner;
 * 
 *     bool rebuild_saturation_matrix;
 * 
 * @endcode
 * 
 * At the very end we declare a variable that denotes the material
 * model. Compared to step-21, we do this here as a member variable since
 * we will want to use it in a variety of places and so having a central
 * place where such a variable is declared will make it simpler to replace
 * one class by another (e.g. replace RandomMedium::KInverse by
 * SingleCurvingCrack::KInverse).
 * 
 * @code
 *     const RandomMedium::KInverse<dim> k_inverse;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimTwoPhaseFlowProblem"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem</h3>
 * 

 * 
 * The constructor of this class is an extension of the constructors in
 * step-21 and step-31. We need to add the various variables that concern
 * the saturation. As discussed in the introduction, we are going to use
 * $Q_2 \times Q_1$ (Taylor-Hood) elements again for the Darcy system, an
 * element combination that fulfills the Ladyzhenskaya-Babuska-Brezzi (LBB)
 * conditions [Brezzi and Fortin 1991, Chen 2005], and $Q_1$ elements for
 * the saturation. However, by using variables that store the polynomial
 * degree of the Darcy and temperature finite elements, it is easy to
 * consistently modify the degree of the elements as well as all quadrature
 * formulas used on them downstream. Moreover, we initialize the time
 * stepping variables related to operator splitting as well as the option
 * for matrix assembly and preconditioning:
 * 
 * @code
 *   template <int dim>
 *   TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree)
 *     : triangulation(Triangulation<dim>::maximum_smoothing)
 *     , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN())
 *     , degree(degree)
 *     , darcy_degree(degree)
 *     , darcy_fe(FE_Q<dim>(darcy_degree + 1), dim, FE_Q<dim>(darcy_degree), 1)
 *     , darcy_dof_handler(triangulation)
 *     ,
 * 
 *     saturation_degree(degree + 1)
 *     , saturation_fe(saturation_degree)
 *     , saturation_dof_handler(triangulation)
 *     ,
 * 
 *     saturation_refinement_threshold(0.5)
 *     ,
 * 
 *     time(0)
 *     , end_time(10)
 *     ,
 * 
 *     current_macro_time_step(0)
 *     , old_macro_time_step(0)
 *     ,
 * 
 *     time_step(0)
 *     , old_time_step(0)
 *     , timestep_number(0)
 *     , viscosity(0.2)
 *     , porosity(1.0)
 *     , AOS_threshold(3.0)
 *     ,
 * 
 *     rebuild_saturation_matrix(true)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimsetup_dofs"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::setup_dofs</h3>
 * 

 * 
 * This is the function that sets up the DoFHandler objects we have here
 * (one for the Darcy part and one for the saturation part) as well as set
 * to the right sizes the various objects required for the linear algebra in
 * this program. Its basic operations are similar to what step-31 did.
 *   

 * 
 * The body of the function first enumerates all degrees of freedom for the
 * Darcy and saturation systems. For the Darcy part, degrees of freedom are
 * then sorted to ensure that velocities precede pressure DoFs so that we
 * can partition the Darcy matrix into a $2 \times 2$ matrix.
 *   

 * 
 * Then, we need to incorporate hanging node constraints and Dirichlet
 * boundary value constraints into darcy_preconditioner_constraints.  The
 * boundary condition constraints are only set on the pressure component
 * since the Schur complement preconditioner that corresponds to the porous
 * media flow operator in non-mixed form, $-\nabla \cdot [\mathbf K
 * \lambda_t(S)]\nabla$, acts only on the pressure variable. Therefore, we
 * use a component_mask that filters out the velocity component, so that the
 * condensation is performed on pressure degrees of freedom only.
 *   

 * 
 * After having done so, we count the number of degrees of freedom in the
 * various blocks. This information is then used to create the sparsity
 * pattern for the Darcy and saturation system matrices as well as the
 * preconditioner matrix from which we build the Darcy preconditioner. As in
 * step-31, we choose to create the pattern using the blocked version of
 * DynamicSparsityPattern. So, for this, we follow the same way as step-31
 * did and we don't have to repeat descriptions again for the rest of the
 * member function.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::setup_dofs()
 *   {
 *     std::vector<unsigned int> darcy_block_component(dim + 1, 0);
 *     darcy_block_component[dim] = 1;
 *     {
 *       darcy_dof_handler.distribute_dofs(darcy_fe);
 *       DoFRenumbering::Cuthill_McKee(darcy_dof_handler);
 *       DoFRenumbering::component_wise(darcy_dof_handler, darcy_block_component);
 * 
 *       darcy_constraints.clear();
 *       DoFTools::make_hanging_node_constraints(darcy_dof_handler,
 *                                               darcy_constraints);
 *       darcy_constraints.close();
 *     }
 *     {
 *       saturation_dof_handler.distribute_dofs(saturation_fe);
 * 
 *       saturation_constraints.clear();
 *       DoFTools::make_hanging_node_constraints(saturation_dof_handler,
 *                                               saturation_constraints);
 *       saturation_constraints.close();
 *     }
 *     {
 *       darcy_preconditioner_constraints.clear();
 * 
 *       FEValuesExtractors::Scalar pressure(dim);
 * 
 *       DoFTools::make_hanging_node_constraints(darcy_dof_handler,
 *                                               darcy_preconditioner_constraints);
 *       DoFTools::make_zero_boundary_constraints(darcy_dof_handler,
 *                                                darcy_preconditioner_constraints,
 *                                                darcy_fe.component_mask(
 *                                                  pressure));
 * 
 *       darcy_preconditioner_constraints.close();
 *     }
 * 
 * 
 *     const std::vector<types::global_dof_index> darcy_dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(darcy_dof_handler,
 *                                         darcy_block_component);
 *     const unsigned int n_u = darcy_dofs_per_block[0],
 *                        n_p = darcy_dofs_per_block[1],
 *                        n_s = saturation_dof_handler.n_dofs();
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << " (on " << triangulation.n_levels() << " levels)" << std::endl
 *               << "Number of degrees of freedom: " << n_u + n_p + n_s << " ("
 *               << n_u << '+' << n_p << '+' << n_s << ')' << std::endl
 *               << std::endl;
 * 
 *     {
 *       darcy_matrix.clear();
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2);
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u);
 *       dsp.block(0, 1).reinit(n_u, n_p);
 *       dsp.block(1, 0).reinit(n_p, n_u);
 *       dsp.block(1, 1).reinit(n_p, n_p);
 * 
 *       dsp.collect_sizes();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 * 
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (!((c == dim) && (d == dim)))
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 * 
 *       DoFTools::make_sparsity_pattern(
 *         darcy_dof_handler, coupling, dsp, darcy_constraints, false);
 * 
 *       darcy_matrix.reinit(dsp);
 *     }
 * 
 *     {
 *       Amg_preconditioner.reset();
 *       Mp_preconditioner.reset();
 *       darcy_preconditioner_matrix.clear();
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2);
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u);
 *       dsp.block(0, 1).reinit(n_u, n_p);
 *       dsp.block(1, 0).reinit(n_p, n_u);
 *       dsp.block(1, 1).reinit(n_p, n_p);
 * 
 *       dsp.collect_sizes();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (c == d)
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       DoFTools::make_sparsity_pattern(
 *         darcy_dof_handler, coupling, dsp, darcy_constraints, false);
 * 
 *       darcy_preconditioner_matrix.reinit(dsp);
 *     }
 * 
 * 
 *     {
 *       saturation_matrix.clear();
 * 
 *       DynamicSparsityPattern dsp(n_s, n_s);
 * 
 *       DoFTools::make_sparsity_pattern(saturation_dof_handler,
 *                                       dsp,
 *                                       saturation_constraints,
 *                                       false);
 * 
 * 
 *       saturation_matrix.reinit(dsp);
 *     }
 * 
 *     std::vector<IndexSet> darcy_partitioning(2);
 *     darcy_partitioning[0] = complete_index_set(n_u);
 *     darcy_partitioning[1] = complete_index_set(n_p);
 *     darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD);
 *     darcy_solution.collect_sizes();
 * 
 *     last_computed_darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD);
 *     last_computed_darcy_solution.collect_sizes();
 * 
 *     second_last_computed_darcy_solution.reinit(darcy_partitioning,
 *                                                MPI_COMM_WORLD);
 *     second_last_computed_darcy_solution.collect_sizes();
 * 
 *     darcy_rhs.reinit(darcy_partitioning, MPI_COMM_WORLD);
 *     darcy_rhs.collect_sizes();
 * 
 *     IndexSet saturation_partitioning = complete_index_set(n_s);
 *     saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD);
 *     old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD);
 *     old_old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD);
 * 
 *     saturation_matching_last_computed_darcy_solution.reinit(
 *       saturation_partitioning, MPI_COMM_WORLD);
 * 
 *     saturation_rhs.reinit(saturation_partitioning, MPI_COMM_WORLD);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Assemblingmatricesandpreconditioners"></a> 
 * <h3>Assembling matrices and preconditioners</h3>
 * 

 * 
 * The next few functions are devoted to setting up the various system and
 * preconditioner matrices and right hand sides that we have to deal with in
 * this program.
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_darcy_preconditioner"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner</h4>
 * 

 * 
 * This function assembles the matrix we use for preconditioning the Darcy
 * system. What we need are a vector mass matrix weighted by
 * $\left(\mathbf{K} \lambda_t\right)^{-1}$ on the velocity components and a
 * mass matrix weighted by $\left(\mathbf{K} \lambda_t\right)$ on the
 * pressure component. We start by generating a quadrature object of
 * appropriate order, the FEValues object that can give values and gradients
 * at the quadrature points (together with quadrature weights). Next we
 * create data structures for the cell matrix and the relation between local
 * and global DoFs. The vectors phi_u and grad_phi_p are going to hold the
 * values of the basis functions in order to faster build up the local
 * matrices, as was already done in step-22. Before we start the loop over
 * all active cells, we have to specify which components are pressure and
 * which are velocity.
 *   

 * 
 * The creation of the local matrix is rather simple. There are only a term
 * weighted by $\left(\mathbf{K} \lambda_t\right)^{-1}$ (on the velocity)
 * and a Laplace matrix weighted by $\left(\mathbf{K} \lambda_t\right)$ to
 * be generated, so the creation of the local matrix is done in essentially
 * two lines. Since the material model functions at the top of this file
 * only provide the inverses of the permeability and mobility, we have to
 * compute $\mathbf K$ and $\lambda_t$ by hand from the given values, once
 * per quadrature point.
 *   

 * 
 * Once the local matrix is ready (loop over rows and columns in the local
 * matrix on each quadrature point), we get the local DoF indices and write
 * the local information into the global matrix. We do this by directly
 * applying the constraints (i.e. darcy_preconditioner_constraints) that
 * takes care of hanging node and zero Dirichlet boundary condition
 * constraints. By doing so, we don't have to do that afterwards, and we
 * later don't have to use AffineConstraints::condense and
 * MatrixTools::apply_boundary_values, both functions that would need to
 * modify matrix and vector entries and so are difficult to write for the
 * Trilinos classes where we don't immediately have access to individual
 * memory locations.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner()
 *   {
 *     std::cout << "   Rebuilding darcy preconditioner..." << std::endl;
 * 
 *     darcy_preconditioner_matrix = 0;
 * 
 *     const QGauss<dim> quadrature_formula(darcy_degree + 2);
 *     FEValues<dim>     darcy_fe_values(darcy_fe,
 *                                   quadrature_formula,
 *                                   update_JxW_values | update_values |
 *                                     update_gradients |
 *                                     update_quadrature_points);
 *     FEValues<dim>     saturation_fe_values(saturation_fe,
 *                                        quadrature_formula,
 *                                        update_values);
 * 
 *     const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 *     std::vector<double> old_saturation_values(n_q_points);
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
 *     std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 *     auto       cell            = darcy_dof_handler.begin_active();
 *     const auto endc            = darcy_dof_handler.end();
 *     auto       saturation_cell = saturation_dof_handler.begin_active();
 * 
 *     for (; cell != endc; ++cell, ++saturation_cell)
 *       {
 *         darcy_fe_values.reinit(cell);
 *         saturation_fe_values.reinit(saturation_cell);
 * 
 *         local_matrix = 0;
 * 
 *         saturation_fe_values.get_function_values(old_saturation_solution,
 *                                                  old_saturation_values);
 * 
 *         k_inverse.value_list(darcy_fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double old_s = old_saturation_values[q];
 * 
 *             const double inverse_mobility = mobility_inverse(old_s, viscosity);
 *             const double mobility         = 1.0 / inverse_mobility;
 *             const Tensor<2, dim> permeability = invert(k_inverse_values[q]);
 * 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 phi_u[k]      = darcy_fe_values[velocities].value(k, q);
 *                 grad_phi_p[k] = darcy_fe_values[pressure].gradient(k, q);
 *               }
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   local_matrix(i, j) +=
 *                     (k_inverse_values[q] * inverse_mobility * phi_u[i] *
 *                        phi_u[j] +
 *                      permeability * mobility * grad_phi_p[i] * grad_phi_p[j]) *
 *                     darcy_fe_values.JxW(q);
 *                 }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         darcy_preconditioner_constraints.distribute_local_to_global(
 *           local_matrix, local_dof_indices, darcy_preconditioner_matrix);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimbuild_darcy_preconditioner"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::build_darcy_preconditioner</h4>
 * 

 * 
 * After calling the above functions to assemble the preconditioner matrix,
 * this function generates the inner preconditioners that are going to be
 * used for the Schur complement block preconditioner. The preconditioners
 * need to be regenerated at every saturation time step since they depend on
 * the saturation $S$ that varies with time.
 *   

 * 
 * In here, we set up the preconditioner for the velocity-velocity matrix
 * $\mathbf{M}^{\mathbf{u}}$ and the Schur complement $\mathbf{S}$. As
 * explained in the introduction, we are going to use an IC preconditioner
 * based on the vector matrix $\mathbf{M}^{\mathbf{u}}$ and another based on
 * the scalar Laplace matrix $\tilde{\mathbf{S}}^p$ (which is spectrally
 * close to the Schur complement of the Darcy matrix). Usually, the
 * TrilinosWrappers::PreconditionIC class can be seen as a good black-box
 * preconditioner which does not need any special knowledge of the matrix
 * structure and/or the operator that's behind it.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::build_darcy_preconditioner()
 *   {
 *     assemble_darcy_preconditioner();
 * 
 *     Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>();
 *     Amg_preconditioner->initialize(darcy_preconditioner_matrix.block(0, 0));
 * 
 *     Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>();
 *     Mp_preconditioner->initialize(darcy_preconditioner_matrix.block(1, 1));
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_darcy_system"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_darcy_system</h4>
 * 

 * 
 * This is the function that assembles the linear system for the Darcy
 * system.
 *   

 * 
 * Regarding the technical details of implementation, the procedures are
 * similar to those in step-22 and step-31. We reset matrix and vector,
 * create a quadrature formula on the cells, and then create the respective
 * FEValues object.
 *   

 * 
 * There is one thing that needs to be commented: since we have a separate
 * finite element and DoFHandler for the saturation, we need to generate a
 * second FEValues object for the proper evaluation of the saturation
 * solution. This isn't too complicated to realize here: just use the
 * saturation structures and set an update flag for the basis function
 * values which we need for evaluation of the saturation solution. The only
 * important part to remember here is that the same quadrature formula is
 * used for both FEValues objects to ensure that we get matching information
 * when we loop over the quadrature points of the two objects.
 *   

 * 
 * The declarations proceed with some shortcuts for array sizes, the
 * creation of the local matrix, right hand side as well as the vector for
 * the indices of the local dofs compared to the global system.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_darcy_system()
 *   {
 *     darcy_matrix = 0;
 *     darcy_rhs    = 0;
 * 
 *     QGauss<dim>     quadrature_formula(darcy_degree + 2);
 *     QGauss<dim - 1> face_quadrature_formula(darcy_degree + 2);
 * 
 *     FEValues<dim> darcy_fe_values(darcy_fe,
 *                                   quadrature_formula,
 *                                   update_values | update_gradients |
 *                                     update_quadrature_points |
 *                                     update_JxW_values);
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe,
 *                                        quadrature_formula,
 *                                        update_values);
 * 
 *     FEFaceValues<dim> darcy_fe_face_values(darcy_fe,
 *                                            face_quadrature_formula,
 *                                            update_values |
 *                                              update_normal_vectors |
 *                                              update_quadrature_points |
 *                                              update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const Functions::ZeroFunction<dim> pressure_right_hand_side;
 *     const PressureBoundaryValues<dim>  pressure_boundary_values;
 * 
 *     std::vector<double>         pressure_rhs_values(n_q_points);
 *     std::vector<double>         boundary_values(n_face_q_points);
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 * @endcode
 * 
 * Next we need a vector that will contain the values of the saturation
 * solution at the previous time level at the quadrature points to
 * assemble the saturation dependent coefficients in the Darcy equations.
 *     

 * 
 * The set of vectors we create next hold the evaluations of the basis
 * functions as well as their gradients that will be used for creating the
 * matrices. Putting these into their own arrays rather than asking the
 * FEValues object for this information each time it is needed is an
 * optimization to accelerate the assembly process, see step-22 for
 * details.
 *     

 * 
 * The last two declarations are used to extract the individual blocks
 * (velocity, pressure, saturation) from the total FE system.
 * 
 * @code
 *     std::vector<double> old_saturation_values(n_q_points);
 * 
 *     std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
 *     std::vector<double>         div_phi_u(dofs_per_cell);
 *     std::vector<double>         phi_p(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 * @endcode
 * 
 * Now start the loop over all cells in the problem. We are working on two
 * different DoFHandlers for this assembly routine, so we must have two
 * different cell iterators for the two objects in use. This might seem a
 * bit peculiar, but since both the Darcy system and the saturation system
 * use the same grid we can assume that the two iterators run in sync over
 * the cells of the two DoFHandler objects.
 *     

 * 
 * The first statements within the loop are again all very familiar, doing
 * the update of the finite element data as specified by the update flags,
 * zeroing out the local arrays and getting the values of the old solution
 * at the quadrature points.  At this point we also have to get the values
 * of the saturation function of the previous time step at the quadrature
 * points. To this end, we can use the FEValues::get_function_values
 * (previously already used in step-9, step-14 and step-15), a function
 * that takes a solution vector and returns a list of function values at
 * the quadrature points of the present cell. In fact, it returns the
 * complete vector-valued solution at each quadrature point, i.e. not only
 * the saturation but also the velocities and pressure.
 *     

 * 
 * Then we are ready to loop over the quadrature points on the cell to do
 * the integration. The formula for this follows in a straightforward way
 * from what has been discussed in the introduction.
 *     

 * 
 * Once this is done, we start the loop over the rows and columns of the
 * local matrix and feed the matrix with the relevant products.
 *     

 * 
 * The last step in the loop over all cells is to enter the local
 * contributions into the global matrix and vector structures to the
 * positions specified in local_dof_indices. Again, we let the
 * AffineConstraints class do the insertion of the cell matrix
 * elements to the global matrix, which already condenses the hanging node
 * constraints.
 * 
 * @code
 *     auto       cell            = darcy_dof_handler.begin_active();
 *     const auto endc            = darcy_dof_handler.end();
 *     auto       saturation_cell = saturation_dof_handler.begin_active();
 * 
 *     for (; cell != endc; ++cell, ++saturation_cell)
 *       {
 *         darcy_fe_values.reinit(cell);
 *         saturation_fe_values.reinit(saturation_cell);
 * 
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 *         saturation_fe_values.get_function_values(old_saturation_solution,
 *                                                  old_saturation_values);
 * 
 *         pressure_right_hand_side.value_list(
 *           darcy_fe_values.get_quadrature_points(), pressure_rhs_values);
 *         k_inverse.value_list(darcy_fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 phi_u[k]     = darcy_fe_values[velocities].value(k, q);
 *                 div_phi_u[k] = darcy_fe_values[velocities].divergence(k, q);
 *                 phi_p[k]     = darcy_fe_values[pressure].value(k, q);
 *               }
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 const double old_s = old_saturation_values[q];
 *                 for (unsigned int j = 0; j <= i; ++j)
 *                   {
 *                     local_matrix(i, j) +=
 *                       (phi_u[i] * k_inverse_values[q] *
 *                          mobility_inverse(old_s, viscosity) * phi_u[j] -
 *                        div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
 *                       darcy_fe_values.JxW(q);
 *                   }
 * 
 *                 local_rhs(i) +=
 *                   (-phi_p[i] * pressure_rhs_values[q]) * darcy_fe_values.JxW(q);
 *               }
 *           }
 * 
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary())
 *             {
 *               darcy_fe_face_values.reinit(cell, face);
 * 
 *               pressure_boundary_values.value_list(
 *                 darcy_fe_face_values.get_quadrature_points(), boundary_values);
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   {
 *                     const Tensor<1, dim> phi_i_u =
 *                       darcy_fe_face_values[velocities].value(i, q);
 * 
 *                     local_rhs(i) +=
 *                       -(phi_i_u * darcy_fe_face_values.normal_vector(q) *
 *                         boundary_values[q] * darcy_fe_face_values.JxW(q));
 *                   }
 *             }
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
 *             local_matrix(i, j) = local_matrix(j, i);
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         darcy_constraints.distribute_local_to_global(
 *           local_matrix, local_rhs, local_dof_indices, darcy_matrix, darcy_rhs);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_system"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_system</h4>
 * 

 * 
 * This function is to assemble the linear system for the saturation
 * transport equation. It calls, if necessary, two other member functions:
 * assemble_saturation_matrix() and assemble_saturation_rhs(). The former
 * function then assembles the saturation matrix that only needs to be
 * changed occasionally. On the other hand, the latter function that
 * assembles the right hand side must be called at every saturation time
 * step.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_system()
 *   {
 *     if (rebuild_saturation_matrix == true)
 *       {
 *         saturation_matrix = 0;
 *         assemble_saturation_matrix();
 *       }
 * 
 *     saturation_rhs = 0;
 *     assemble_saturation_rhs();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_matrix"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_matrix</h4>
 * 

 * 
 * This function is easily understood since it only forms a simple mass
 * matrix for the left hand side of the saturation linear system by basis
 * functions phi_i_s and phi_j_s only. Finally, as usual, we enter the local
 * contribution into the global matrix by specifying the position in
 * local_dof_indices. This is done by letting the AffineConstraints class do
 * the insertion of the cell matrix elements to the global matrix, which
 * already condenses the hanging node constraints.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix()
 *   {
 *     QGauss<dim> quadrature_formula(saturation_degree + 2);
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe,
 *                                        quadrature_formula,
 *                                        update_values | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = saturation_fe.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *       {
 *         saturation_fe_values.reinit(cell);
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const double phi_i_s = saturation_fe_values.shape_value(i, q);
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const double phi_j_s = saturation_fe_values.shape_value(j, q);
 *                   local_matrix(i, j) +=
 *                     porosity * phi_i_s * phi_j_s * saturation_fe_values.JxW(q);
 *                 }
 *             }
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         saturation_constraints.distribute_local_to_global(local_matrix,
 *                                                           local_dof_indices,
 *                                                           saturation_matrix);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs</h4>
 * 

 * 
 * This function is to assemble the right hand side of the saturation
 * transport equation. Before going about it, we have to create two FEValues
 * objects for the Darcy and saturation systems respectively and, in
 * addition, two FEFaceValues objects for the two systems because we have a
 * boundary integral term in the weak form of saturation equation. For the
 * FEFaceValues object of the saturation system, we also require normal
 * vectors, which we request using the update_normal_vectors flag.
 *   

 * 
 * Next, before looping over all the cells, we have to compute some
 * parameters (e.g. global_u_infty, global_S_variation, and
 * global_Omega_diameter) that the artificial viscosity $\nu$ needs. This is
 * largely the same as was done in step-31, so you may see there for more
 * information.
 *   

 * 
 * The real works starts with the loop over all the saturation and Darcy
 * cells to put the local contributions into the global vector. In this
 * loop, in order to simplify the implementation, we split some of the work
 * into two helper functions: assemble_saturation_rhs_cell_term and
 * assemble_saturation_rhs_boundary_term.  We note that we insert cell or
 * boundary contributions into the global vector in the two functions rather
 * than in this present function.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs()
 *   {
 *     QGauss<dim>     quadrature_formula(saturation_degree + 2);
 *     QGauss<dim - 1> face_quadrature_formula(saturation_degree + 2);
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe,
 *                                        quadrature_formula,
 *                                        update_values | update_gradients |
 *                                          update_quadrature_points |
 *                                          update_JxW_values);
 *     FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values);
 *     FEFaceValues<dim> saturation_fe_face_values(saturation_fe,
 *                                                 face_quadrature_formula,
 *                                                 update_values |
 *                                                   update_normal_vectors |
 *                                                   update_quadrature_points |
 *                                                   update_JxW_values);
 *     FEFaceValues<dim> darcy_fe_face_values(darcy_fe,
 *                                            face_quadrature_formula,
 *                                            update_values);
 *     FEFaceValues<dim> saturation_fe_face_values_neighbor(
 *       saturation_fe, face_quadrature_formula, update_values);
 * 
 *     const unsigned int dofs_per_cell =
 *       saturation_dof_handler.get_fe().n_dofs_per_cell();
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const double                    global_max_u_F_prime = get_max_u_F_prime();
 *     const std::pair<double, double> global_S_range =
 *       get_extrapolated_saturation_range();
 *     const double global_S_variation =
 *       global_S_range.second - global_S_range.first;
 * 
 *     auto       cell       = saturation_dof_handler.begin_active();
 *     const auto endc       = saturation_dof_handler.end();
 *     auto       darcy_cell = darcy_dof_handler.begin_active();
 *     for (; cell != endc; ++cell, ++darcy_cell)
 *       {
 *         saturation_fe_values.reinit(cell);
 *         darcy_fe_values.reinit(darcy_cell);
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         assemble_saturation_rhs_cell_term(saturation_fe_values,
 *                                           darcy_fe_values,
 *                                           global_max_u_F_prime,
 *                                           global_S_variation,
 *                                           local_dof_indices);
 * 
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary())
 *             {
 *               darcy_fe_face_values.reinit(darcy_cell, face);
 *               saturation_fe_face_values.reinit(cell, face);
 *               assemble_saturation_rhs_boundary_term(saturation_fe_face_values,
 *                                                     darcy_fe_face_values,
 *                                                     local_dof_indices);
 *             }
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs_cell_term"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term</h4>
 * 

 * 
 * This function takes care of integrating the cell terms of the right hand
 * side of the saturation equation, and then assembling it into the global
 * right hand side vector. Given the discussion in the introduction, the
 * form of these contributions is clear. The only tricky part is getting the
 * artificial viscosity and all that is necessary to compute it. The first
 * half of the function is devoted to this task.
 *   

 * 
 * The last part of the function is copying the local contributions into the
 * global vector with position specified in local_dof_indices.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term(
 *     const FEValues<dim> &                       saturation_fe_values,
 *     const FEValues<dim> &                       darcy_fe_values,
 *     const double                                global_max_u_F_prime,
 *     const double                                global_S_variation,
 *     const std::vector<types::global_dof_index> &local_dof_indices)
 *   {
 *     const unsigned int dofs_per_cell = saturation_fe_values.dofs_per_cell;
 *     const unsigned int n_q_points    = saturation_fe_values.n_quadrature_points;
 * 
 *     std::vector<double>         old_saturation_solution_values(n_q_points);
 *     std::vector<double>         old_old_saturation_solution_values(n_q_points);
 *     std::vector<Tensor<1, dim>> old_grad_saturation_solution_values(n_q_points);
 *     std::vector<Tensor<1, dim>> old_old_grad_saturation_solution_values(
 *       n_q_points);
 *     std::vector<Vector<double>> present_darcy_solution_values(
 *       n_q_points, Vector<double>(dim + 1));
 * 
 *     saturation_fe_values.get_function_values(old_saturation_solution,
 *                                              old_saturation_solution_values);
 *     saturation_fe_values.get_function_values(
 *       old_old_saturation_solution, old_old_saturation_solution_values);
 *     saturation_fe_values.get_function_gradients(
 *       old_saturation_solution, old_grad_saturation_solution_values);
 *     saturation_fe_values.get_function_gradients(
 *       old_old_saturation_solution, old_old_grad_saturation_solution_values);
 *     darcy_fe_values.get_function_values(darcy_solution,
 *                                         present_darcy_solution_values);
 * 
 *     const double nu =
 *       compute_viscosity(old_saturation_solution_values,
 *                         old_old_saturation_solution_values,
 *                         old_grad_saturation_solution_values,
 *                         old_old_grad_saturation_solution_values,
 *                         present_darcy_solution_values,
 *                         global_max_u_F_prime,
 *                         global_S_variation,
 *                         saturation_fe_values.get_cell()->diameter());
 * 
 *     Vector<double> local_rhs(dofs_per_cell);
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 *           const double   old_s = old_saturation_solution_values[q];
 *           Tensor<1, dim> present_u;
 *           for (unsigned int d = 0; d < dim; ++d)
 *             present_u[d] = present_darcy_solution_values[q](d);
 * 
 *           const double         phi_i_s = saturation_fe_values.shape_value(i, q);
 *           const Tensor<1, dim> grad_phi_i_s =
 *             saturation_fe_values.shape_grad(i, q);
 * 
 *           local_rhs(i) +=
 *             (time_step * fractional_flow(old_s, viscosity) * present_u *
 *                grad_phi_i_s -
 *              time_step * nu * old_grad_saturation_solution_values[q] *
 *                grad_phi_i_s +
 *              porosity * old_s * phi_i_s) *
 *             saturation_fe_values.JxW(q);
 *         }
 * 
 *     saturation_constraints.distribute_local_to_global(local_rhs,
 *                                                       local_dof_indices,
 *                                                       saturation_rhs);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs_boundary_term"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term</h4>
 * 

 * 
 * The next function is responsible for the boundary integral terms in the
 * right hand side form of the saturation equation.  For these, we have to
 * compute the upwinding flux on the global boundary faces, i.e. we impose
 * Dirichlet boundary conditions weakly only on inflow parts of the global
 * boundary. As before, this has been described in step-21 so we refrain
 * from giving more descriptions about that.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term(
 *     const FEFaceValues<dim> &                   saturation_fe_face_values,
 *     const FEFaceValues<dim> &                   darcy_fe_face_values,
 *     const std::vector<types::global_dof_index> &local_dof_indices)
 *   {
 *     const unsigned int dofs_per_cell = saturation_fe_face_values.dofs_per_cell;
 *     const unsigned int n_face_q_points =
 *       saturation_fe_face_values.n_quadrature_points;
 * 
 *     Vector<double> local_rhs(dofs_per_cell);
 * 
 *     std::vector<double> old_saturation_solution_values_face(n_face_q_points);
 *     std::vector<Vector<double>> present_darcy_solution_values_face(
 *       n_face_q_points, Vector<double>(dim + 1));
 *     std::vector<double> neighbor_saturation(n_face_q_points);
 * 
 *     saturation_fe_face_values.get_function_values(
 *       old_saturation_solution, old_saturation_solution_values_face);
 *     darcy_fe_face_values.get_function_values(
 *       darcy_solution, present_darcy_solution_values_face);
 * 
 *     SaturationBoundaryValues<dim> saturation_boundary_values;
 *     saturation_boundary_values.value_list(
 *       saturation_fe_face_values.get_quadrature_points(), neighbor_saturation);
 * 
 *     for (unsigned int q = 0; q < n_face_q_points; ++q)
 *       {
 *         Tensor<1, dim> present_u_face;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           present_u_face[d] = present_darcy_solution_values_face[q](d);
 * 
 *         const double normal_flux =
 *           present_u_face * saturation_fe_face_values.normal_vector(q);
 * 
 *         const bool is_outflow_q_point = (normal_flux >= 0);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           local_rhs(i) -=
 *             time_step * normal_flux *
 *             fractional_flow((is_outflow_q_point == true ?
 *                                old_saturation_solution_values_face[q] :
 *                                neighbor_saturation[q]),
 *                             viscosity) *
 *             saturation_fe_face_values.shape_value(i, q) *
 *             saturation_fe_face_values.JxW(q);
 *       }
 *     saturation_constraints.distribute_local_to_global(local_rhs,
 *                                                       local_dof_indices,
 *                                                       saturation_rhs);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimsolve"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::solve</h3>
 * 

 * 
 * This function implements the operator splitting algorithm, i.e. in each
 * time step it either re-computes the solution of the Darcy system or
 * extrapolates velocity/pressure from previous time steps, then determines
 * the size of the time step, and then updates the saturation variable. The
 * implementation largely follows similar code in step-31. It is, next to
 * the run() function, the central one in this program.
 *   

 * 
 * At the beginning of the function, we ask whether to solve the
 * pressure-velocity part by evaluating the a posteriori criterion (see the
 * following function). If necessary, we will solve the pressure-velocity
 * part using the GMRES solver with the Schur complement block
 * preconditioner as is described in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::solve()
 *   {
 *     const bool solve_for_pressure_and_velocity =
 *       determine_whether_to_solve_for_pressure_and_velocity();
 * 
 *     if (solve_for_pressure_and_velocity == true)
 *       {
 *         std::cout << "   Solving Darcy (pressure-velocity) system..."
 *                   << std::endl;
 * 
 *         assemble_darcy_system();
 *         build_darcy_preconditioner();
 * 
 *         {
 *           const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                                              TrilinosWrappers::PreconditionIC>
 *             mp_inverse(darcy_preconditioner_matrix.block(1, 1),
 *                        *Mp_preconditioner);
 * 
 *           const LinearSolvers::BlockSchurPreconditioner<
 *             TrilinosWrappers::PreconditionIC,
 *             TrilinosWrappers::PreconditionIC>
 *             preconditioner(darcy_matrix, mp_inverse, *Amg_preconditioner);
 * 
 *           SolverControl solver_control(darcy_matrix.m(),
 *                                        1e-16 * darcy_rhs.l2_norm());
 * 
 *           SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres(
 *             solver_control,
 *             SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
 *               100));
 * 
 *           for (unsigned int i = 0; i < darcy_solution.size(); ++i)
 *             if (darcy_constraints.is_constrained(i))
 *               darcy_solution(i) = 0;
 * 
 *           gmres.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner);
 * 
 *           darcy_constraints.distribute(darcy_solution);
 * 
 *           std::cout << "        ..." << solver_control.last_step()
 *                     << " GMRES iterations." << std::endl;
 *         }
 * 
 *         {
 *           second_last_computed_darcy_solution = last_computed_darcy_solution;
 *           last_computed_darcy_solution        = darcy_solution;
 * 
 *           saturation_matching_last_computed_darcy_solution =
 *             saturation_solution;
 *         }
 *       }
 * @endcode
 * 
 * On the other hand, if we have decided that we don't want to compute the
 * solution of the Darcy system for the current time step, then we need to
 * simply extrapolate the previous two Darcy solutions to the same time as
 * we would have computed the velocity/pressure at. We do a simple linear
 * extrapolation, i.e. given the current length $dt$ of the macro time
 * step from the time when we last computed the Darcy solution to now
 * (given by <code>current_macro_time_step</code>), and $DT$ the length of
 * the last macro time step (given by <code>old_macro_time_step</code>),
 * then we get $u^\ast = u_p + dt \frac{u_p-u_{pp}}{DT} = (1+dt/DT)u_p -
 * dt/DT u_{pp}$, where $u_p$ and $u_{pp}$ are the last two computed Darcy
 * solutions. We can implement this formula using just two lines of code.
 *     

 * 
 * Note that the algorithm here only works if we have at least two
 * previously computed Darcy solutions from which we can extrapolate to
 * the current time, and this is ensured by requiring re-computation of
 * the Darcy solution for the first 2 time steps.
 * 
 * @code
 *     else
 *       {
 *         darcy_solution = last_computed_darcy_solution;
 *         darcy_solution.sadd(1 + current_macro_time_step / old_macro_time_step,
 *                             -current_macro_time_step / old_macro_time_step,
 *                             second_last_computed_darcy_solution);
 *       }
 * 
 * 
 * @endcode
 * 
 * With the so computed velocity vector, compute the optimal time step
 * based on the CFL criterion discussed in the introduction...
 * 
 * @code
 *     {
 *       old_time_step = time_step;
 * 
 *       const double max_u_F_prime = get_max_u_F_prime();
 *       if (max_u_F_prime > 0)
 *         time_step = porosity * GridTools::minimal_cell_diameter(triangulation) /
 *                     saturation_degree / max_u_F_prime / 50;
 *       else
 *         time_step = end_time - time;
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * ...and then also update the length of the macro time steps we use while
 * we're dealing with time step sizes. In particular, this involves: (i)
 * If we have just recomputed the Darcy solution, then the length of the
 * previous macro time step is now fixed and the length of the current
 * macro time step is, up to now, simply the length of the current (micro)
 * time step. (ii) If we have not recomputed the Darcy solution, then the
 * length of the current macro time step has just grown by
 * <code>time_step</code>.
 * 
 * @code
 *     if (solve_for_pressure_and_velocity == true)
 *       {
 *         old_macro_time_step     = current_macro_time_step;
 *         current_macro_time_step = time_step;
 *       }
 *     else
 *       current_macro_time_step += time_step;
 * 
 * @endcode
 * 
 * The last step in this function is to recompute the saturation solution
 * based on the velocity field we've just obtained. This naturally happens
 * in every time step, and we don't skip any of these computations. At the
 * end of computing the saturation, we project back into the allowed
 * interval $[0,1]$ to make sure our solution remains physical.
 * 
 * @code
 *     {
 *       std::cout << "   Solving saturation transport equation..." << std::endl;
 * 
 *       assemble_saturation_system();
 * 
 *       SolverControl solver_control(saturation_matrix.m(),
 *                                    1e-16 * saturation_rhs.l2_norm());
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);
 * 
 *       TrilinosWrappers::PreconditionIC preconditioner;
 *       preconditioner.initialize(saturation_matrix);
 * 
 *       cg.solve(saturation_matrix,
 *                saturation_solution,
 *                saturation_rhs,
 *                preconditioner);
 * 
 *       saturation_constraints.distribute(saturation_solution);
 *       project_back_saturation();
 * 
 *       std::cout << "        ..." << solver_control.last_step()
 *                 << " CG iterations." << std::endl;
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimrefine_mesh"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::refine_mesh</h3>
 * 

 * 
 * The next function does the refinement and coarsening of the mesh. It does
 * its work in three blocks: (i) Compute refinement indicators by looking at
 * the gradient of a solution vector extrapolated linearly from the previous
 * two using the respective sizes of the time step (or taking the only
 * solution we have if this is the first time step). (ii) Flagging those
 * cells for refinement and coarsening where the gradient is larger or
 * smaller than a certain threshold, preserving minimal and maximal levels
 * of mesh refinement. (iii) Transferring the solution from the old to the
 * new mesh. None of this is particularly difficult.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::refine_mesh(const unsigned int min_grid_level,
 *                                              const unsigned int max_grid_level)
 *   {
 *     Vector<double> refinement_indicators(triangulation.n_active_cells());
 *     {
 *       const QMidpoint<dim>        quadrature_formula;
 *       FEValues<dim>               fe_values(saturation_fe,
 *                               quadrature_formula,
 *                               update_gradients);
 *       std::vector<Tensor<1, dim>> grad_saturation(1);
 * 
 *       TrilinosWrappers::MPI::Vector extrapolated_saturation_solution(
 *         saturation_solution);
 *       if (timestep_number != 0)
 *         extrapolated_saturation_solution.sadd((1. + time_step / old_time_step),
 *                                               time_step / old_time_step,
 *                                               old_saturation_solution);
 * 
 *       for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *         {
 *           const unsigned int cell_no = cell->active_cell_index();
 *           fe_values.reinit(cell);
 *           fe_values.get_function_gradients(extrapolated_saturation_solution,
 *                                            grad_saturation);
 * 
 *           refinement_indicators(cell_no) = grad_saturation[0].norm();
 *         }
 *     }
 * 
 *     {
 *       for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *         {
 *           const unsigned int cell_no = cell->active_cell_index();
 *           cell->clear_coarsen_flag();
 *           cell->clear_refine_flag();
 * 
 *           if ((static_cast<unsigned int>(cell->level()) < max_grid_level) &&
 *               (std::fabs(refinement_indicators(cell_no)) >
 *                saturation_refinement_threshold))
 *             cell->set_refine_flag();
 *           else if ((static_cast<unsigned int>(cell->level()) >
 *                     min_grid_level) &&
 *                    (std::fabs(refinement_indicators(cell_no)) <
 *                     0.5 * saturation_refinement_threshold))
 *             cell->set_coarsen_flag();
 *         }
 *     }
 * 
 *     triangulation.prepare_coarsening_and_refinement();
 * 
 *     {
 *       std::vector<TrilinosWrappers::MPI::Vector> x_saturation(3);
 *       x_saturation[0] = saturation_solution;
 *       x_saturation[1] = old_saturation_solution;
 *       x_saturation[2] = saturation_matching_last_computed_darcy_solution;
 * 
 *       std::vector<TrilinosWrappers::MPI::BlockVector> x_darcy(2);
 *       x_darcy[0] = last_computed_darcy_solution;
 *       x_darcy[1] = second_last_computed_darcy_solution;
 * 
 *       SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> saturation_soltrans(
 *         saturation_dof_handler);
 * 
 *       SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> darcy_soltrans(
 *         darcy_dof_handler);
 * 
 * 
 *       triangulation.prepare_coarsening_and_refinement();
 *       saturation_soltrans.prepare_for_coarsening_and_refinement(x_saturation);
 * 
 *       darcy_soltrans.prepare_for_coarsening_and_refinement(x_darcy);
 * 
 *       triangulation.execute_coarsening_and_refinement();
 *       setup_dofs();
 * 
 *       std::vector<TrilinosWrappers::MPI::Vector> tmp_saturation(3);
 *       tmp_saturation[0].reinit(saturation_solution);
 *       tmp_saturation[1].reinit(saturation_solution);
 *       tmp_saturation[2].reinit(saturation_solution);
 *       saturation_soltrans.interpolate(x_saturation, tmp_saturation);
 * 
 *       saturation_solution                              = tmp_saturation[0];
 *       old_saturation_solution                          = tmp_saturation[1];
 *       saturation_matching_last_computed_darcy_solution = tmp_saturation[2];
 * 
 *       saturation_constraints.distribute(saturation_solution);
 *       saturation_constraints.distribute(old_saturation_solution);
 *       saturation_constraints.distribute(
 *         saturation_matching_last_computed_darcy_solution);
 * 
 *       std::vector<TrilinosWrappers::MPI::BlockVector> tmp_darcy(2);
 *       tmp_darcy[0].reinit(darcy_solution);
 *       tmp_darcy[1].reinit(darcy_solution);
 *       darcy_soltrans.interpolate(x_darcy, tmp_darcy);
 * 
 *       last_computed_darcy_solution        = tmp_darcy[0];
 *       second_last_computed_darcy_solution = tmp_darcy[1];
 * 
 *       darcy_constraints.distribute(last_computed_darcy_solution);
 *       darcy_constraints.distribute(second_last_computed_darcy_solution);
 * 
 *       rebuild_saturation_matrix = true;
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimoutput_results"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::output_results</h3>
 * 

 * 
 * This function generates graphical output. It is in essence a copy of the
 * implementation in step-31.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::output_results() const
 *   {
 *     const FESystem<dim> joint_fe(darcy_fe, 1, saturation_fe, 1);
 *     DoFHandler<dim>     joint_dof_handler(triangulation);
 *     joint_dof_handler.distribute_dofs(joint_fe);
 *     Assert(joint_dof_handler.n_dofs() ==
 *              darcy_dof_handler.n_dofs() + saturation_dof_handler.n_dofs(),
 *            ExcInternalError());
 * 
 *     Vector<double> joint_solution(joint_dof_handler.n_dofs());
 * 
 *     {
 *       std::vector<types::global_dof_index> local_joint_dof_indices(
 *         joint_fe.n_dofs_per_cell());
 *       std::vector<types::global_dof_index> local_darcy_dof_indices(
 *         darcy_fe.n_dofs_per_cell());
 *       std::vector<types::global_dof_index> local_saturation_dof_indices(
 *         saturation_fe.n_dofs_per_cell());
 * 
 *       auto       joint_cell      = joint_dof_handler.begin_active();
 *       const auto joint_endc      = joint_dof_handler.end();
 *       auto       darcy_cell      = darcy_dof_handler.begin_active();
 *       auto       saturation_cell = saturation_dof_handler.begin_active();
 * 
 *       for (; joint_cell != joint_endc;
 *            ++joint_cell, ++darcy_cell, ++saturation_cell)
 *         {
 *           joint_cell->get_dof_indices(local_joint_dof_indices);
 *           darcy_cell->get_dof_indices(local_darcy_dof_indices);
 *           saturation_cell->get_dof_indices(local_saturation_dof_indices);
 * 
 *           for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i)
 *             if (joint_fe.system_to_base_index(i).first.first == 0)
 *               {
 *                 Assert(joint_fe.system_to_base_index(i).second <
 *                          local_darcy_dof_indices.size(),
 *                        ExcInternalError());
 *                 joint_solution(local_joint_dof_indices[i]) = darcy_solution(
 *                   local_darcy_dof_indices[joint_fe.system_to_base_index(i)
 *                                             .second]);
 *               }
 *             else
 *               {
 *                 Assert(joint_fe.system_to_base_index(i).first.first == 1,
 *                        ExcInternalError());
 *                 Assert(joint_fe.system_to_base_index(i).second <
 *                          local_darcy_dof_indices.size(),
 *                        ExcInternalError());
 *                 joint_solution(local_joint_dof_indices[i]) =
 *                   saturation_solution(
 *                     local_saturation_dof_indices
 *                       [joint_fe.system_to_base_index(i).second]);
 *               }
 *         }
 *     }
 *     std::vector<std::string> joint_solution_names(dim, "velocity");
 *     joint_solution_names.emplace_back("pressure");
 *     joint_solution_names.emplace_back("saturation");
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(joint_dof_handler);
 *     data_out.add_data_vector(joint_solution,
 *                              joint_solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 * 
 *     data_out.build_patches();
 * 
 *     std::string filename =
 *       "solution-" + Utilities::int_to_string(timestep_number, 5) + ".vtu";
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Toolfunctions"></a> 
 * <h3>Tool functions</h3>
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemdimdetermine_whether_to_solve_for_pressure_and_velocity"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity</h4>
 * 

 * 
 * This function implements the a posteriori criterion for adaptive operator
 * splitting. The function is relatively straightforward given the way we
 * have implemented other functions above and given the formula for the
 * criterion derived in the paper.
 *   

 * 
 * If one decides that one wants the original IMPES method in which the
 * Darcy equation is solved in every time step, then this can be achieved by
 * setting the threshold value <code>AOS_threshold</code> (with a default of
 * $5.0$) to zero, thereby forcing the function to always return true.
 *   

 * 
 * Finally, note that the function returns true unconditionally for the
 * first two time steps to ensure that we have always solved the Darcy
 * system at least twice when skipping its solution, thereby allowing us to
 * extrapolate the velocity from the last two solutions in
 * <code>solve()</code>.
 * 
 * @code
 *   template <int dim>
 *   bool TwoPhaseFlowProblem<
 *     dim>::determine_whether_to_solve_for_pressure_and_velocity() const
 *   {
 *     if (timestep_number <= 2)
 *       return true;
 * 
 *     const QGauss<dim>  quadrature_formula(saturation_degree + 2);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(saturation_fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points);
 * 
 *     std::vector<double> old_saturation_after_solving_pressure(n_q_points);
 *     std::vector<double> present_saturation(n_q_points);
 * 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 *     double max_global_aop_indicator = 0.0;
 * 
 *     for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *       {
 *         double max_local_mobility_reciprocal_difference = 0.0;
 *         double max_local_permeability_inverse_l1_norm   = 0.0;
 * 
 *         fe_values.reinit(cell);
 *         fe_values.get_function_values(
 *           saturation_matching_last_computed_darcy_solution,
 *           old_saturation_after_solving_pressure);
 *         fe_values.get_function_values(saturation_solution, present_saturation);
 * 
 *         k_inverse.value_list(fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double mobility_reciprocal_difference = std::fabs(
 *               mobility_inverse(present_saturation[q], viscosity) -
 *               mobility_inverse(old_saturation_after_solving_pressure[q],
 *                                viscosity));
 * 
 *             max_local_mobility_reciprocal_difference =
 *               std::max(max_local_mobility_reciprocal_difference,
 *                        mobility_reciprocal_difference);
 * 
 *             max_local_permeability_inverse_l1_norm =
 *               std::max(max_local_permeability_inverse_l1_norm,
 *                        l1_norm(k_inverse_values[q]));
 *           }
 * 
 *         max_global_aop_indicator =
 *           std::max(max_global_aop_indicator,
 *                    (max_local_mobility_reciprocal_difference *
 *                     max_local_permeability_inverse_l1_norm));
 *       }
 * 
 *     return (max_global_aop_indicator > AOS_threshold);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimproject_back_saturation"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::project_back_saturation</h4>
 * 

 * 
 * The next function simply makes sure that the saturation values always
 * remain within the physically reasonable range of $[0,1]$. While the
 * continuous equations guarantee that this is so, the discrete equations
 * don't. However, if we allow the discrete solution to escape this range we
 * get into trouble because terms like $F(S)$ and $F'(S)$ will produce
 * unreasonable results (e.g. $F'(S)<0$ for $S<0$, which would imply that
 * the wetting fluid phase flows <i>against</i> the direction of the bulk
 * fluid velocity)). Consequently, at the end of each time step, we simply
 * project the saturation field back into the physically reasonable region.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::project_back_saturation()
 *   {
 *     for (unsigned int i = 0; i < saturation_solution.size(); ++i)
 *       if (saturation_solution(i) < 0.2)
 *         saturation_solution(i) = 0.2;
 *       else if (saturation_solution(i) > 1)
 *         saturation_solution(i) = 1;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimget_max_u_F_prime"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::get_max_u_F_prime</h4>
 *   

 * 
 * Another simpler helper function: Compute the maximum of the total
 * velocity times the derivative of the fraction flow function, i.e.,
 * compute $\|\mathbf{u} F'(S)\|_{L_\infty(\Omega)}$. This term is used in
 * both the computation of the time step as well as in normalizing the
 * entropy-residual term in the artificial viscosity.
 * 
 * @code
 *   template <int dim>
 *   double TwoPhaseFlowProblem<dim>::get_max_u_F_prime() const
 *   {
 *     const QGauss<dim>  quadrature_formula(darcy_degree + 2);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values);
 *     FEValues<dim> saturation_fe_values(saturation_fe,
 *                                        quadrature_formula,
 *                                        update_values);
 * 
 *     std::vector<Vector<double>> darcy_solution_values(n_q_points,
 *                                                       Vector<double>(dim + 1));
 *     std::vector<double>         saturation_values(n_q_points);
 * 
 *     double max_velocity_times_dF_dS = 0;
 * 
 *     auto       cell            = darcy_dof_handler.begin_active();
 *     const auto endc            = darcy_dof_handler.end();
 *     auto       saturation_cell = saturation_dof_handler.begin_active();
 *     for (; cell != endc; ++cell, ++saturation_cell)
 *       {
 *         darcy_fe_values.reinit(cell);
 *         saturation_fe_values.reinit(saturation_cell);
 * 
 *         darcy_fe_values.get_function_values(darcy_solution,
 *                                             darcy_solution_values);
 *         saturation_fe_values.get_function_values(old_saturation_solution,
 *                                                  saturation_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             Tensor<1, dim> velocity;
 *             for (unsigned int i = 0; i < dim; ++i)
 *               velocity[i] = darcy_solution_values[q](i);
 * 
 *             const double dF_dS =
 *               fractional_flow_derivative(saturation_values[q], viscosity);
 * 
 *             max_velocity_times_dF_dS =
 *               std::max(max_velocity_times_dF_dS, velocity.norm() * dF_dS);
 *           }
 *       }
 * 
 *     return max_velocity_times_dF_dS;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimget_extrapolated_saturation_range"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range</h4>
 *   

 * 
 * For computing the stabilization term, we need to know the range of the
 * saturation variable. Unlike in step-31, this range is trivially bounded
 * by the interval $[0,1]$ but we can do a bit better by looping over a
 * collection of quadrature points and seeing what the values are there. If
 * we can, i.e., if there are at least two timesteps around, we can even
 * take the values extrapolated to the next time step.
 *   

 * 
 * As before, the function is taken with minimal modifications from step-31.
 * 
 * @code
 *   template <int dim>
 *   std::pair<double, double>
 *   TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range() const
 *   {
 *     const QGauss<dim>  quadrature_formula(saturation_degree + 2);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(saturation_fe, quadrature_formula, update_values);
 *     std::vector<double> old_saturation_values(n_q_points);
 *     std::vector<double> old_old_saturation_values(n_q_points);
 * 
 *     if (timestep_number != 0)
 *       {
 *         double min_saturation = std::numeric_limits<double>::max(),
 *                max_saturation = -std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *           {
 *             fe_values.reinit(cell);
 *             fe_values.get_function_values(old_saturation_solution,
 *                                           old_saturation_values);
 *             fe_values.get_function_values(old_old_saturation_solution,
 *                                           old_old_saturation_values);
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q)
 *               {
 *                 const double saturation =
 *                   (1. + time_step / old_time_step) * old_saturation_values[q] -
 *                   time_step / old_time_step * old_old_saturation_values[q];
 * 
 *                 min_saturation = std::min(min_saturation, saturation);
 *                 max_saturation = std::max(max_saturation, saturation);
 *               }
 *           }
 * 
 *         return std::make_pair(min_saturation, max_saturation);
 *       }
 *     else
 *       {
 *         double min_saturation = std::numeric_limits<double>::max(),
 *                max_saturation = -std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : saturation_dof_handler.active_cell_iterators())
 *           {
 *             fe_values.reinit(cell);
 *             fe_values.get_function_values(old_saturation_solution,
 *                                           old_saturation_values);
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q)
 *               {
 *                 const double saturation = old_saturation_values[q];
 * 
 *                 min_saturation = std::min(min_saturation, saturation);
 *                 max_saturation = std::max(max_saturation, saturation);
 *               }
 *           }
 * 
 *         return std::make_pair(min_saturation, max_saturation);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimcompute_viscosity"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::compute_viscosity</h4>
 *   

 * 
 * The final tool function is used to compute the artificial viscosity on a
 * given cell. This isn't particularly complicated if you have the formula
 * for it in front of you, and looking at the implementation in step-31. The
 * major difference to that tutorial program is that the velocity here is
 * not simply $\mathbf u$ but $\mathbf u F'(S)$ and some of the formulas
 * need to be adjusted accordingly.
 * 
 * @code
 *   template <int dim>
 *   double TwoPhaseFlowProblem<dim>::compute_viscosity(
 *     const std::vector<double> &        old_saturation,
 *     const std::vector<double> &        old_old_saturation,
 *     const std::vector<Tensor<1, dim>> &old_saturation_grads,
 *     const std::vector<Tensor<1, dim>> &old_old_saturation_grads,
 *     const std::vector<Vector<double>> &present_darcy_values,
 *     const double                       global_max_u_F_prime,
 *     const double                       global_S_variation,
 *     const double                       cell_diameter) const
 *   {
 *     const double beta  = .4 * dim;
 *     const double alpha = 1;
 * 
 *     if (global_max_u_F_prime == 0)
 *       return 5e-3 * cell_diameter;
 * 
 *     const unsigned int n_q_points = old_saturation.size();
 * 
 *     double max_residual             = 0;
 *     double max_velocity_times_dF_dS = 0;
 * 
 *     const bool use_dF_dS = true;
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         Tensor<1, dim> u;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           u[d] = present_darcy_values[q](d);
 * 
 *         const double dS_dt = porosity *
 *                              (old_saturation[q] - old_old_saturation[q]) /
 *                              old_time_step;
 * 
 *         const double dF_dS = fractional_flow_derivative(
 *           (old_saturation[q] + old_old_saturation[q]) / 2.0, viscosity);
 * 
 *         const double u_grad_S =
 *           u * dF_dS * (old_saturation_grads[q] + old_old_saturation_grads[q]) /
 *           2.0;
 * 
 *         const double residual =
 *           std::abs((dS_dt + u_grad_S) *
 *                    std::pow((old_saturation[q] + old_old_saturation[q]) / 2,
 *                             alpha - 1.));
 * 
 *         max_residual = std::max(residual, max_residual);
 *         max_velocity_times_dF_dS =
 *           std::max(std::sqrt(u * u) * (use_dF_dS ? std::max(dF_dS, 1.) : 1),
 *                    max_velocity_times_dF_dS);
 *       }
 * 
 *     const double c_R            = 1.0;
 *     const double global_scaling = c_R * porosity *
 *                                   (global_max_u_F_prime)*global_S_variation /
 *                                   std::pow(global_Omega_diameter, alpha - 2.);
 * 
 *     return (beta *
 *             (max_velocity_times_dF_dS)*std::min(cell_diameter,
 *                                                 std::pow(cell_diameter, alpha) *
 *                                                   max_residual /
 *                                                   global_scaling));
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimrun"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::run</h3>
 * 

 * 
 * This function is, besides <code>solve()</code>, the primary function of
 * this program as it controls the time iteration as well as when the
 * solution is written into output files and when to do mesh refinement.
 *   

 * 
 * With the exception of the startup code that loops back to the beginning
 * of the function through the <code>goto start_time_iteration</code> label,
 * everything should be relatively straightforward. In any case, it mimics
 * the corresponding function in step-31.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::run()
 *   {
 *     const unsigned int initial_refinement     = (dim == 2 ? 5 : 2);
 *     const unsigned int n_pre_refinement_steps = (dim == 2 ? 3 : 2);
 * 
 * 
 *     GridGenerator::hyper_cube(triangulation, 0, 1);
 *     triangulation.refine_global(initial_refinement);
 *     global_Omega_diameter = GridTools::diameter(triangulation);
 * 
 *     setup_dofs();
 * 
 *     unsigned int pre_refinement_step = 0;
 * 
 *   start_time_iteration:
 * 
 *     VectorTools::project(saturation_dof_handler,
 *                          saturation_constraints,
 *                          QGauss<dim>(saturation_degree + 2),
 *                          SaturationInitialValues<dim>(),
 *                          old_saturation_solution);
 * 
 *     time_step = old_time_step = 0;
 *     current_macro_time_step = old_macro_time_step = 0;
 * 
 *     time = 0;
 * 
 *     do
 *       {
 *         std::cout << "Timestep " << timestep_number << ":  t=" << time
 *                   << ", dt=" << time_step << std::endl;
 * 
 *         solve();
 * 
 *         std::cout << std::endl;
 * 
 *         if (timestep_number % 200 == 0)
 *           output_results();
 * 
 *         if (timestep_number % 25 == 0)
 *           refine_mesh(initial_refinement,
 *                       initial_refinement + n_pre_refinement_steps);
 * 
 *         if ((timestep_number == 0) &&
 *             (pre_refinement_step < n_pre_refinement_steps))
 *           {
 *             ++pre_refinement_step;
 *             goto start_time_iteration;
 *           }
 * 
 *         time += time_step;
 *         ++timestep_number;
 * 
 *         old_old_saturation_solution = old_saturation_solution;
 *         old_saturation_solution     = saturation_solution;
 *       }
 *     while (time <= end_time);
 *   }
 * } // namespace Step43
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * The main function looks almost the same as in all other programs. The need
 * to initialize the MPI subsystem for a program that uses Trilinos -- even
 * for programs that do not actually run in parallel -- is explained in
 * step-31.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step43;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *         argc, argv, numbers::invalid_unsigned_int);
 * 
 * @endcode
 * 
 * This program can only be run in serial. Otherwise, throw an exception.
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
 *                   ExcMessage(
 *                     "This program can only be run in serial, use ./step-43"));
 * 
 *       TwoPhaseFlowProblem<2> two_phase_flow_problem(1);
 *       two_phase_flow_problem.run();
 *     }
 *   catch (std::exception &exc)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Exception on processing: " << std::endl
 *                 << exc.what() << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 * 
 *       return 1;
 *     }
 *   catch (...)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Unknown exception!" << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       return 1;
 *     }
 * 
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>



The output of this program is not really much different from that of
step-21: it solves the same problem, after all. Of more importance are
quantitative metrics such as the accuracy of the solution as well as
the time needed to compute it. These are documented in detail in the
two publications listed at the top of this page and we won't repeat
them here.

That said, no tutorial program is complete without a couple of good
pictures, so here is some output of a run in 3d:

<table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.velocity.png" alt="">
	<p align="center">
        Velocity vectors of flow through the porous medium with random
        permeability model. Streaming paths of high permeability and resulting
        high velocity are clearly visible.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.streamlines.png" alt="">
	<p align="center">
        Streamlines colored by the saturation along the streamline path. Blue
        streamlines indicate low saturations, i.e., the flow along these
	streamlines must be slow or else more fluid would have been
        transported along them. On the other hand, green paths indicate high
        velocities since the fluid front has already reached further into the
        domain.
	</p>
    </td>
  </tr>
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.saturation.png" alt="">
	<p align="center">
        Streamlines with a volume rendering of the saturation, showing how far
        the fluid front has advanced at this time.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.mesh.png" alt="">
	<p align="center">
	Surface of the mesh showing the adaptive refinement along the front.
	</p>
    </td>
  </tr>
</table>


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


The primary objection one may have to this program is that it is still too
slow: 3d computations on reasonably fine meshes are simply too expensive to be
done routinely and with reasonably quick turn-around. This is similar to the
situation we were in when we wrote step-31, from which this program has taken
much inspiration. The solution is similar as it was there as well: We need to
parallelize the program in a way similar to how we derived step-32 out of
step-31. In fact, all of the techniques used in step-32 would be transferable
to this program as well, making the program run on dozens or hundreds of
processors immediately.

A different direction is to make the program more relevant to many other
porous media applications. Specifically, one avenue is to go to the primary
user of porous media flow simulators, namely the oil industry. There,
applications in this area are dominated by multiphase flow (i.e., more than
the two phases we have here), and the reactions they may have with each other
(or any other way phases may exchange mass, such as through dissolution in and
bubbling out of gas from the oil phase). Furthermore, the presence of gas
often leads to compressibility effects of the fluid. Jointly, these effects
are typically formulated in the widely-used "black oil model". True reactions
between multiple phases also play a role in oil reservoir modeling when
considering controlled burns of oil in the reservoir to raise pressure and
temperature. These are much more complex problems, though, and left for future
projects.

Finally, from a mathematical perspective, we have derived the
criterion for re-computing the velocity/pressure solution at a given
time step under the assumption that we want to compare the solution we
would get at the current time step with that computed the last time we
actually solved this system. However, in the program, whenever we did
not re-compute the solution, we didn't just use the previously
computed solution but instead extrapolated from the previous two times
we solved the system. Consequently, the criterion was pessimistically
stated: what we should really compare is the solution we would get at
the current time step with the extrapolated one. Re-stating the
theorem in this regard is left as an exercise.

There are also other ways to extend the mathematical foundation of
this program; for example, one may say that it isn't the velocity we
care about, but in fact the saturation. Thus, one may ask whether the
criterion we use here to decide whether $\mathbf u$ needs to be
recomputed is appropriate; one may, for example, suggest that it is
also important to decide whether (and by how much) a wrong velocity
field in fact affects the solution of the saturation equation. This
would then naturally lead to a sensitivity analysis.

From an algorithmic viewpoint, we have here used a criterion for refinement
that is often used in engineering, namely by looking at the gradient of
the solution. However, if you inspect the solution, you will find that
it quickly leads to refinement almost everywhere, even in regions where it
is clearly not necessary: frequently used therefore does not need to imply
that it is a useful criterion to begin with. On the other hand, replacing
this criterion by a different and better one should not be very difficult.
For example, the KellyErrorEstimator class used in many other programs
should certainly be applicable to the current problem as well.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-43.cc"
*/
