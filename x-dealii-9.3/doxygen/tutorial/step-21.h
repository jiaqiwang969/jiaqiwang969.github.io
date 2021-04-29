/**
@page step_21 The step-21 tutorial program
This tutorial depends on step-20.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetwophaseflowproblem">The two phase flow problem</a>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Spacediscretization">Space discretization</a>
        <li><a href="#Linearsolvers">Linear solvers</a>
        <li><a href="#Choosingatimestep">Choosing a time step</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeTwoPhaseFlowProblemcodeclass">The <code>TwoPhaseFlowProblem</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
      <ul>
        <li><a href="#Pressurerighthandside">Pressure right hand side</a>
        <li><a href="#Pressureboundaryvalues">Pressure boundary values</a>
        <li><a href="#Saturationboundaryvalues">Saturation boundary values</a>
        <li><a href="#Initialdata">Initial data</a>
      </ul>
        <li><a href="#Theinversepermeabilitytensor">The inverse permeability tensor</a>
      <ul>
        <li><a href="#Singlecurvingcrackpermeability">Single curving crack permeability</a>
        <li><a href="#Randommediumpermeability">Random medium permeability</a>
      </ul>
        <li><a href="#Theinversemobilityandsaturationfunctions">The inverse mobility and saturation functions</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
        <li><a href="#codeTwoPhaseFlowProblemcodeclassimplementation"><code>TwoPhaseFlowProblem</code> class implementation</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemTwoPhaseFlowProblem">TwoPhaseFlowProblem::TwoPhaseFlowProblem</a>
        <li><a href="#TwoPhaseFlowProblemmake_grid_and_dofs">TwoPhaseFlowProblem::make_grid_and_dofs</a>
        <li><a href="#TwoPhaseFlowProblemassemble_system">TwoPhaseFlowProblem::assemble_system</a>
        <li><a href="#TwoPhaseFlowProblemassemble_rhs_S">TwoPhaseFlowProblem::assemble_rhs_S</a>
        <li><a href="#TwoPhaseFlowProblemsolve">TwoPhaseFlowProblem::solve</a>
        <li><a href="#TwoPhaseFlowProblemoutput_results">TwoPhaseFlowProblem::output_results</a>
        <li><a href="#TwoPhaseFlowProblemproject_back_saturation">TwoPhaseFlowProblem::project_back_saturation</a>
        <li><a href="#TwoPhaseFlowProblemget_maximal_velocity">TwoPhaseFlowProblem::get_maximal_velocity</a>
        <li><a href="#TwoPhaseFlowProblemrun">TwoPhaseFlowProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Solvers">Solvers</a>
        <li><a href="#Timestepping">Time stepping</a>
        <li><a href="#Adaptivity">Adaptivity</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Introduction"></a><a name="Intro"></a> <h1>Introduction</h1>


This program grew out of a student project by Yan Li at Texas A&amp;M
University. Most of the work for this program is by her.

In this project, we propose a numerical simulation for two phase
flow problems in porous media. This problem includes one
elliptic equation and one nonlinear, time dependent transport
equation. This is therefore also the first time-dependent tutorial
program (besides the somewhat strange time-dependence of @ref step_18
"step-18").

The equations covered here are an extension of the material already covered in
step-20. In particular, they fall into the class of
vector-valued problems. A toplevel overview of this topic can be found in the
@ref vector_valued module.


<a name="Thetwophaseflowproblem"></a><h3>The two phase flow problem</h3>


Modeling of two phase flow in porous media is important for both
environmental remediation and the management of petroleum and groundwater
reservoirs. Practical situations involving two phase flow include the
dispersal of a nonaqueous phase liquid in an aquifer, or the joint
movement of a mixture of fluids such as oil and water in a
reservoir. Simulation models, if they are to provide realistic
predictions, must accurately account for these effects.

To derive the governing equations, consider two phase flow in a
reservoir $\Omega$ under the assumption that the movement of fluids is
dominated by viscous effects; i.e. we neglect the effects of gravity,
compressibility, and capillary pressure. Porosity will be considered
to be constant. We will denote variables referring to either of the two
phases using subscripts $w$ and $o$, short for water and oil. The
derivation of the equations holds for other pairs of fluids as well,
however.

The velocity with which molecules of each of the two phases move is
determined by Darcy's law that states that the velocity is
proportional to the pressure gradient:
@f{eqnarray*}
  \mathbf{u}_{j}
  =
  -\frac{k_{rj}(S)}{\mu_{j}} \mathbf{K} \cdot \nabla p
@f}
where $\mathbf{u}_{j}$ is the velocity of phase $j=o,w$, $K$ is the
permeability tensor, $k_{rj}$ is the relative permeability of phase
$j$, $p$ is the
pressure and $\mu_{j}$ is the viscosity of phase $j$. Finally, $S$ is
the saturation (volume fraction), i.e. a function with values between
0 and 1 indicating the composition of the mixture of fluids. In
general, the coefficients $K, k_{rj}, \mu$ may be spatially dependent
variables, and we will always treat them as non-constant functions in
the following.

We combine Darcy's law with the statement of conservation of mass for
each phase,
@f[
  \textrm{div}\ \mathbf{u}_{j} = q_j,
@f]
with a source term for each phase. By summing over the two phases,
we can express the governing equations in terms of the
so-called pressure equation:
@f{eqnarray*}
- \nabla \cdot (\mathbf{K}\lambda(S) \nabla p)= q.
@f}
Here, $q$ is the sum source term, and
@f[
  \lambda(S) = \frac{k_{rw}(S)}{\mu_{w}}+\frac{k_{ro}(S)}{\mu_{o}}
@f]
is the total mobility.

So far, this looks like an ordinary stationary, Poisson-like equation that we
can solve right away with the techniques of the first few tutorial programs
(take a look at step-6, for example, for something very
similar). However, we have not said anything yet about the saturation, which
of course is going to change as the fluids move around.

The second part of the equations is the description of the
dynamics of the saturation, i.e., how the relative concentration of the
two fluids changes with time. The saturation equation for the displacing
fluid (water) is given by the following conservation law:
@f{eqnarray*}
  S_{t} + \nabla \cdot (F(S) \mathbf{u}) = q_{w},
@f}
which can be rewritten by using the product rule of the divergence operator
in the previous equation:
@f{eqnarray*}
  S_{t} + F(S) \left[\nabla \cdot \mathbf{u}\right]
        + \mathbf{u} \cdot \left[ \nabla F(S)\right]
  = S_{t} + F(S) q + \mathbf{u} \cdot \nabla F(S) = q_{w}.
@f}
Here, $q=\nabla\cdot \mathbf{u}$ is the total influx introduced
above, and $q_{w}$ is the flow rate of the displacing fluid (water).
These two are related to the fractional flow $F(S)$ in the following way:
@f[
  q_{w} = F(S) q,
@f]
where the fractional flow is often parameterized via the (heuristic) expression
@f[
  F(S)
  =
  \frac{k_{rw}(S)/\mu_{w}}{k_{rw}(S)/\mu_{w} + k_{ro}(S)/\mu_{o}}.
@f]
Putting it all together yields the saturation equation in the following,
advected form:
@f{eqnarray*}
  S_{t} + \mathbf{u} \cdot \nabla F(S) = 0,
@f}
where $\mathbf u$ is the total velocity
@f[
  \mathbf{u} =
  \mathbf{u}_{o} + \mathbf{u}_{w} = -\lambda(S) \mathbf{K}\cdot\nabla p.
@f]
Note that the advection equation contains the term $\mathbf{u} \cdot \nabla
F(S)$ rather than $\mathbf{u} \cdot \nabla S$ to indicate that the saturation
is not simply transported along; rather, since the two phases move with
different velocities, the saturation can actually change even in the advected
coordinate system. To see this, rewrite $\mathbf{u} \cdot \nabla F(S)
= \mathbf{u} F'(S) \cdot \nabla S$ to observe that the <i>actual</i>
velocity with which the phase with saturation $S$ is transported is
$\mathbf u F'(S)$ whereas the other phase is transported at velocity
$\mathbf u (1-F'(S))$. $F(S)$ is consequently often referred to as the
<i>fractional flow</i>.

In summary, what we get are the following two equations:
@f{eqnarray*}
  - \nabla \cdot (\mathbf{K}\lambda(S) \nabla p) &=& q
  \qquad \textrm{in}\ \Omega\times[0,T],
  \\
  S_{t} + \mathbf{u} \cdot \nabla F(S) &=& 0
  \qquad \textrm{in}\ \Omega\times[0,T].
@f}
Here, $p=p(\mathbf x, t), S=S(\mathbf x, t)$ are now time dependent
functions: while at every time instant the flow field is in
equilibrium with the pressure (i.e. we neglect dynamic
accelerations), the saturation is transported along with the flow and
therefore changes over time, in turn affected the flow field again
through the dependence of the first equation on $S$.

This set of equations has a peculiar character: one of the two
equations has a time derivative, the other one doesn't. This
corresponds to the character that the pressure and velocities are
coupled through an instantaneous constraint, whereas the saturation
evolves over finite time scales.

Such systems of equations are called Differential Algebraic Equations
(DAEs), since one of the equations is a differential equation, the
other is not (at least not with respect to the time variable) and is
therefore an "algebraic" equation. (The notation comes from the field
of ordinary differential equations, where everything that does not
have derivatives with respect to the time variable is necessarily an
algebraic equation.) This class of equations contains pretty
well-known cases: for example, the time dependent Stokes and
Navier-Stokes equations (where the algebraic constraint is that the
divergence of the flow field, $\textrm{div}\ \mathbf u$, must be zero)
as well as the time dependent Maxwell equations (here, the algebraic
constraint is that the divergence of the electric displacement field
equals the charge density, $\textrm{div}\ \mathbf D = \rho$ and that the
divergence of the magnetic flux density is zero: $\textrm{div}\ \mathbf
B = 0$); even the quasistatic model of step-18 falls into this
category. We will see that the different character of the two equations
will inform our discretization strategy for the two equations.


<a name="Timediscretization"></a><h3>Time discretization</h3>


In the reservoir simulation community, it is common to solve the equations
derived above by going back to the first order, mixed formulation. To this
end, we re-introduce the total velocity $\mathbf u$ and write the equations in
the following form:
@f{eqnarray*}
  \mathbf{u}+\mathbf{K}\lambda(S) \nabla p&=&0 \\
  \nabla \cdot\mathbf{u} &=& q \\
  S_{t} + \mathbf{u} \cdot \nabla F(S) &=& 0.
@f}
This formulation has the additional benefit that we do not have to express the
total velocity $\mathbf u$ appearing in the transport equation as a function
of the pressure, but can rather take the primary variable for it. Given the
saddle point structure of the first two equations and their similarity to the
mixed Laplace formulation we have introduced in step-20, it
will come as no surprise that we will use a mixed discretization again.

But let's postpone this for a moment. The first business we have with these
equations is to think about the time discretization. In reservoir simulation,
there is a rather standard algorithm that we will use here. It first solves
the pressure using an implicit equation, then the saturation using an explicit
time stepping scheme. The algorithm is called IMPES for IMplicit Pressure
Explicit Saturation and was first proposed a long time ago: by Sheldon et
al. in 1959 and Stone and Gardner in 1961 (J. W. Sheldon, B. Zondek and
W. T. Cardwell: <i>One-dimensional, incompressible, non-capillary, two-phase
fluid flow in a porous medium</i>, Trans. SPE AIME, 216 (1959), pp. 290-296; H.
L. Stone and A. O. Gardner Jr: <i>Analysis of gas-cap or dissolved-gas
reservoirs</i>, Trans. SPE AIME, 222 (1961), pp. 92-104).
In a slightly modified form, this algorithm can be
written as follows: for each time step, solve
@f{eqnarray*}
  \mathbf{u}^{n+1}+\mathbf{K}\lambda(S^n) \nabla p^{n+1}&=&0 \\
  \nabla \cdot\mathbf{u}^{n+1} &=& q^{n+1} \\
  \frac {S^{n+1}-S^n}{\triangle t} + \mathbf{u}^{n+1} \cdot \nabla F(S^n) &=& 0,
@f}
where $\triangle t$ is the length of a time step. Note how we solve the
implicit pressure-velocity system that only depends on the previously computed
saturation $S^n$, and then do an explicit time step for $S^{n+1}$ that only
depends on the previously known $S^n$ and the just computed
$\mathbf{u}^{n+1}$. This way, we never have to iterate for the nonlinearities
of the system as we would have if we used a fully implicit method. (In
a more modern perspective, this should be seen as an "operator
splitting" method. step-58 has a long description of the idea behind this.)

We can then state the problem in weak form as follows, by multiplying each
equation with test functions $\mathbf v$, $\phi$, and $\sigma$ and integrating
terms by parts:
@f{eqnarray*}
  \left((\mathbf{K}\lambda(S^n))^{-1} \mathbf{u}^{n+1},\mathbf v\right)_\Omega -
  (p^{n+1}, \nabla\cdot\mathbf v)_\Omega &=&
  - (p^{n+1}, \mathbf v)_{\partial\Omega}
  \\
  (\nabla \cdot\mathbf{u}^{n+1}, \phi)_\Omega &=& (q^{n+1},\phi)_\Omega
@f}
Note that in the first term, we have to prescribe the pressure $p^{n+1}$ on
the boundary $\partial\Omega$ as boundary values for our problem. $\mathbf n$
denotes the unit outward normal vector to $\partial K$, as usual.

For the saturation equation, we obtain after integrating by parts
@f{eqnarray*}
  (S^{n+1}, \sigma)_\Omega
  -
  \triangle t
  \sum_K
  \left\{
  \left(F(S^n), \nabla \cdot (\mathbf{u}^{n+1} \sigma)\right)_K
  -
  \left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}, \sigma\right)_{\partial K}
  \right\}
  &=&
  (S^n,\sigma)_\Omega.
@f}
Using the fact that $\nabla \cdot \mathbf{u}^{n+1}=q^{n+1}$, we can rewrite the
cell term to get an equation as follows:
@f{eqnarray*}
  (S^{n+1}, \sigma)_\Omega
  -
  \triangle t
  \sum_K
  \left\{
  \left(F(S^n) \mathbf{u}^{n+1}, \nabla \sigma\right)_K
  -
  \left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}), \sigma\right)_{\partial K}
  \right\}
  &=&
  (S^n,\sigma)_\Omega +
  \triangle t \sum_K  \left(F(S^n) q^{n+1}, \sigma\right)_K.
@f}
We introduce an object of type DiscreteTime in order to keep track of the
current value of time and time step in the code. This class encapsulates many
complexities regarding adjusting time step size and stopping at a specified
final time.



<a name="Spacediscretization"></a><h3>Space discretization</h3>


In each time step, we then apply the mixed finite method of @ref step_20
"step-20" to the velocity and pressure. To be well-posed, we choose
Raviart-Thomas spaces $RT_{k}$ for $\mathbf{u}$ and discontinuous elements of
class $DGQ_{k}$ for $p$. For the saturation, we will also choose $DGQ_{k}$
spaces.

Since we have discontinuous spaces, we have to think about how to evaluate
terms on the interfaces between cells, since discontinuous functions are not
really defined there. In particular, we have to give a meaning to the last
term on the left hand side of the saturation equation. To this end, let us
define that we want to evaluate it in the following sense:
@f{eqnarray*}
  &&\left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}), \sigma\right)_{\partial K}
  \\
  &&\qquad =
  \left(F(S^n_+) (\mathbf n \cdot \mathbf{u}^{n+1}_+), \sigma\right)_{\partial K_+}
  +
  \left(F(S^n_-) (\mathbf n \cdot \mathbf{u}^{n+1}_-), \sigma\right)_{\partial K_-},
@f}
where $\partial K_{-} \dealcoloneq \{x\in \partial K, \mathbf{u}(x) \cdot \mathbf{n}<0\}$
denotes the inflow boundary and $\partial K_{+} \dealcoloneq \{\partial K \setminus
\partial K_{-}\}$ is the outflow part of the boundary.
The quantities $S_+,\mathbf{u}_+$ then correspond to the values of these
variables on the present cell, whereas $S_-,\mathbf{u}_-$ (needed on the
inflow part of the boundary of $K$) are quantities taken from the neighboring
cell. Some more context on discontinuous element techniques and evaluation of
fluxes can also be found in step-12 and step-12b.


<a name="Linearsolvers"></a><h3>Linear solvers</h3>


The linear solvers used in this program are a straightforward extension of the
ones used in step-20 (but without LinearOperator). Essentially, we simply have
to extend everything from
two to three solution components. If we use the discrete spaces
mentioned above and put shape functions into the bilinear forms, we
arrive at the following linear system to be solved for time step $n+1$:
@f[
\left(
\begin{array}{ccc}
M^u(S^{n}) & B^{T}& 0\\
B &    0 & 0\\
\triangle t\; H &    0& M^S
\end{array}
\right)
\left(
\begin{array}{c}
\mathbf{U}^{n+1} \\ P^{n+1} \\ S^{n+1}
\end{array}
\right)
=
\left(
\begin{array}{c}
0 \\ F_2 \\ F_3
\end{array}
\right)
@f]
where the individual matrices and vectors are defined as follows using
shape functions $\mathbf v_i$ (of type Raviart Thomas $RT_k$) for
velocities and $\phi_i$ (of type $DGQ_k$) for both pressures and saturations:
@f{eqnarray*}
M^u(S^n)_{ij} &=&
\left((\mathbf{K}\lambda(S^n))^{-1} \mathbf{v}_i,\mathbf
v_j\right)_\Omega,
\\
B_{ij} &=&
-(\nabla \cdot \mathbf v_j, \phi_i)_\Omega,
\\
H_{ij} &=&
  -
  \sum_K
  \left\{
  \left(F(S^n) \mathbf v_i, \nabla \phi_j)\right)_K
  -
  \left(F(S^n_+) (\mathbf n \cdot (\mathbf v_i)_+), \phi_j\right)_{\partial K_+}
  -
  \left(F(S^n_-) (\mathbf n \cdot (\mathbf v_i)_-), \phi_j\right)_{\partial K_-},
  \right\}
\\
M^S_{ij} &=&
(\phi_i, \phi_j)_\Omega,
\\
(F_2)_i &=&
-(q^{n+1},\phi_i)_\Omega,
\\
(F_3)_i &=&
(S^n,\phi_i)_\Omega +\triangle t \sum_K  \left(F(S^n) q^{n+1}, \phi_i\right)_K.
@f}

@note Due to historical accidents, the role of matrices $B$ and $B^T$
has been reverted in this program compared to step-20. In other words,
here $B$ refers to the divergence and $B^T$ to the gradient operators
when it was the other way around in step-20.

The system above presents a complication: Since the matrix $H_{ij}$
depends on $\mathbf u^{n+1}$ implicitly (the velocities are needed to
determine which parts of the boundaries $\partial K$ of cells are
influx or outflux parts), we can only assemble this matrix after we
have solved for the velocities.

The solution scheme then involves the following steps:
<ol>
  <li>Solve for the pressure $p^{n+1}$ using the Schur complement
  technique introduced in step-20.

  <li>Solve for the velocity $\mathbf u^{n+1}$ as also discussed in
  step-20.

  <li>Compute the term $F_3-\triangle t\; H \mathbf u^{n+1}$, using
  the just computed velocities.

  <li>Solve for the saturation $S^{n+1}$.
</ol>

In this scheme, we never actually build the matrix $H$, but rather
generate the right hand side of the third equation once we are ready
to do so.

In the program, we use a variable <code>solution</code> to store the
solution of the present time step. At the end of each step, we copy
its content, i.e. all three of its block components, into the variable
<code>old_solution</code> for use in the next time step.


<a name="Choosingatimestep"></a><h3>Choosing a time step</h3>


A general rule of thumb in hyperbolic transport equations like the equation we
have to solve for the saturation equation is that if we use an explicit time
stepping scheme, then we should use a time step such that the distance that a
particle can travel within one time step is no larger than the diameter of a
single cell. In other words, here, we should choose
@f[
  \triangle t_{n+1} \le \frac h{|\mathbf{u}^{n+1}(\mathbf{x})|}.
@f]
Fortunately, we are in a position where we can do that: we only need the
time step when we want to assemble the right hand side of the saturation
equation, which is after we have already solved for $\mathbf{u}^{n+1}$. All we
therefore have to do after solving for the velocity is to loop over all
quadrature points in the domain and determine the maximal magnitude of the
velocity. We can then set the time step for the saturation equation to
@f[
  \triangle t_{n+1} = \frac {\min_K h_K}{\max_{\mathbf{x}}|\mathbf{u}^{n+1}(\mathbf{x})|}.
@f]

Why is it important to do this? If we don't, then we will end up with lots of
places where our saturation is larger than one or less than zero, as can
easily be verified. (Remember that the saturation corresponds to something
like the water fraction in the fluid mixture, and therefore must physically be
between 0 and 1.) On the other hand, if we choose our time step according to
the criterion listed above, this only happens very very infrequently &mdash;
in fact only once for the entire run of the program. However, to be on the
safe side, however, we run a function <code>project_back_saturation</code> at
the end of each time step, that simply projects the saturation back onto the
interval $[0,1]$, should it have gotten out of the physical range. This is
useful since the functions $\lambda(S)$ and $F(S)$ do not represent anything
physical outside this range, and we should not expect the program to do
anything useful once we have negative saturations or ones larger than one.

Note that we will have similar restrictions on the time step also in
step-23 and step-24 where we solve the time dependent
wave equation, another hyperbolic problem. We will also come back to the issue
of time step choice below in the section on <a href="#extensions">possible
extensions to this program</a>.


<a name="Thetestcase"></a><h3>The test case</h3>


For simplicity, this program assumes that there is no source, $q=0$, and that
the heterogeneous porous medium is isotropic $\mathbf{K}(\mathbf{x}) =
k(\mathbf{x}) \mathbf{I}$. The first one of these is a realistic assumption in
oil reservoirs: apart from injection and production wells, there are usually
no mechanisms for fluids to appear or disappear out of the blue. The second
one is harder to justify: on a microscopic level, most rocks are isotropic,
because they consist of a network of interconnected pores. However, this
microscopic scale is out of the range of today's computer simulations, and we
have to be content with simulating things on the scale of meters. On that
scale, however, fluid transport typically happens through a network of cracks
in the rock, rather than through pores. However, cracks often result from
external stress fields in the rock layer (for example from tectonic faulting)
and the cracks are therefore roughly aligned. This leads to a situation where
the permeability is often orders of magnitude larger in the direction parallel
to the cracks than perpendicular to the cracks. A problem typically faces in
reservoir simulation, however, is that the modeler doesn't know the direction
of cracks because oil reservoirs are not accessible to easy inspection. The
only solution in that case is to assume an effective, isotropic permeability.

Whatever the matter, both of these restrictions, no sources and isotropy,
would be easy to lift with a few lines of code in the program.

Next, for simplicity, our numerical simulation will be done on the
unit cell $\Omega = [0,1]\times [0,1]$ for $t\in [0,T]$. Our initial
conditions are $S(\mathbf{x},0)=0$; in the oil reservoir picture, where $S$
would indicate the water saturation, this means that the reservoir contains
pure oil at the beginning. Note that we do not need any initial
conditions for pressure or velocity, since the equations do not contain time
derivatives of these variables. Finally, we impose the following pressure
boundary conditions:
@f[
  p(\mathbf{x},t)=1-x_1 \qquad \textrm{on}\ \partial\Omega.
@f]
Since the pressure and velocity solve a mixed form Poisson equation, the
imposed pressure leads to a resulting flow field for the velocity. On the
other hand, this flow field determines whether a piece of the boundary is of
inflow or outflow type, which is of relevance because we have to impose
boundary conditions for the saturation on the inflow part of the boundary,
@f[
  \Gamma_{in}(t) = \{\mathbf{x}\in\partial\Omega:
                     \mathbf{n} \cdot \mathbf{u}(\mathbf{x},t) < 0\}.
@f]
On this inflow boundary, we impose the following saturation values:
@f{eqnarray}
  S(\mathbf{x},t) = 1 & \textrm{on}\ \Gamma_{in}\cap\{x_1=0\},
  \\
  S(\mathbf{x},t) = 0 & \textrm{on}\ \Gamma_{in}\backslash \{x_1=0\}.
@f}
In other words, we have pure water entering the reservoir at the left, whereas
the other parts of the boundary are in contact with undisturbed parts of the
reservoir and whenever influx occurs on these boundaries, pure oil will enter.

In our simulations, we choose the total mobility as
@f[
  \lambda (S) = \frac{1.0}{\mu} S^2 +(1-S)^2
@f]
where we use $\mu=0.2$ for the viscosity. In addition, the fractional flow of
water is given by
@f[
  F(S)=\frac{S^2}{S^2+\mu (1-S)^2}
@f]

@note Coming back to this testcase in step-43 several years later revealed an
oddity in the setup of this testcase. To this end, consider that we can
rewrite the advection equation for the saturation as $S_{t} + (\mathbf{u}
F'(S)) \cdot \nabla S = 0$. Now, at the initial time, we have $S=0$, and with
the given choice of function $F(S)$, we happen to have $F'(0)=0$. In other
words, at $t=0$, the equation reduces to $S_t=0$ for all $\mathbf x$, so the
saturation is zero everywhere and it is going to stay zero everywhere! This is
despite the fact that $\mathbf u$ is not necessarily zero: the combined fluid
is moving, but we've chosen our partial flux $F(S)$ in such a way that
infinitesimal amounts of wetting fluid also only move at infinitesimal speeds
(i.e., they stick to the medium more than the non-wetting phase in which they
are embedded). That said, how can we square this with the knowledge that
wetting fluid is invading from the left, leading to the flow patterns seen in
the <a href="#Results">results section</a>? That's where we get into
mathematics: Equations like the transport equation we are considering here
have infinitely many solutions, but only one of them is physical: the one that
results from the so-called viscosity limit, called the <a
href="http://en.wikipedia.org/wiki/Viscosity_solution">viscosity
solution</a>. The thing is that with discontinuous elements we arrive at this
viscosity limit because using a numerical flux introduces a finite amount of
artificial viscosity into the numerical scheme. On the other hand, in step-43,
we use an artificial viscosity that is proportional to $\|\mathbf u F'(S)\|$
on every cell, which at the initial time is zero. Thus, the saturation there is
zero and remains zero; the solution we then get is <i>one</i> solution of the
advection equation, but the method does not converge to the viscosity solution
without further changes. We will therefore use a different initial condition in
that program.


Finally, to come back to the description of the testcase, we will show results
for computations with the two permeability
functions introduced at the end of the results section of @ref step_20
"step-20":
<ul>
  <li>A function that models a single, winding crack that snakes through the
  domain. In analogy to step-20, but taking care of the slightly
  different geometry we have here, we describe this by the following function:
  @f[
    k(\mathbf x)
    =
    \max \left\{ e^{-\left(\frac{x_2-\frac 12 - 0.1\sin(10x_1)}{0.1}\right)^2}, 0.01 \right\}.
  @f]
  Taking the maximum is necessary to ensure that the ratio between maximal and
  minimal permeability remains bounded. If we don't do that, permeabilities
  will span many orders of magnitude. On the other hand, the ratio between
  maximal and minimal permeability is a factor in the condition number of the
  Schur complement matrix, and if too large leads to problems for which our
  linear solvers will no longer converge properly.

  <li>A function that models a somewhat random medium. Here, we choose
  @f{eqnarray*}
    k(\mathbf x)
    &=&
    \min \left\{ \max \left\{ \sum_{i=1}^N \sigma_i(\mathbf{x}), 0.01 \right\}, 4\right\},
    \\
    \sigma_i(\mathbf x)
    &=&
    e^{-\left(\frac{|\mathbf{x}-\mathbf{x}_i|}{0.05}\right)^2},
  @f}
  where the centers $\mathbf{x}_i$ are $N$ randomly chosen locations inside
  the domain. This function models a domain in which there are $N$ centers of
  higher permeability (for example where rock has cracked) embedded in a
  matrix of more pristine, unperturbed background rock. Note that here we have
  cut off the permeability function both above and below to ensure a bounded
  condition number.
</ul>
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * This program is an adaptation of step-20 and includes some technique of DG
 * methods from step-12. A good part of the program is therefore very similar
 * to step-20 and we will not comment again on these parts. Only the new stuff
 * will be discussed in more detail.
 * 

 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * All of these include files have been used before:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
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
 * #include <deal.II/fe/fe_raviart_thomas.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * @endcode
 * 
 * In this program, we use a tensor-valued coefficient. Since it may have a
 * spatial dependence, we consider it a tensor-valued function. The following
 * include file provides the <code>TensorFunction</code> class that offers
 * such functionality:
 * 
 * @code
 * #include <deal.II/base/tensor_function.h>
 * 
 * @endcode
 * 
 * Additionally, we use the class <code>DiscreteTime</code> to perform
 * operations related to time incrementation.
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step21
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeTwoPhaseFlowProblemcodeclass"></a> 
 * <h3>The <code>TwoPhaseFlowProblem</code> class</h3>
 * 

 * 
 * This is the main class of the program. It is close to the one of step-20,
 * but with a few additional functions:
 *   

 * 
 * <ul> <li><code>assemble_rhs_S</code> assembles the right hand side of the
 * saturation equation. As explained in the introduction, this can't be
 * integrated into <code>assemble_rhs</code> since it depends on the
 * velocity that is computed in the first part of the time step.
 *   

 * 
 * <li><code>get_maximal_velocity</code> does as its name suggests. This
 * function is used in the computation of the time step size.
 *   

 * 
 * <li><code>project_back_saturation</code> resets all saturation degrees
 * of freedom with values less than zero to zero, and all those with
 * saturations greater than one to one.  </ul>
 *   

 * 
 * The rest of the class should be pretty much obvious. The
 * <code>viscosity</code> variable stores the viscosity $\mu$ that enters
 * several of the formulas in the nonlinear equations. The variable
 * <code>time</code> keeps track of the time information within the
 * simulation.
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
 *     void   make_grid_and_dofs();
 *     void   assemble_system();
 *     void   assemble_rhs_S();
 *     double get_maximal_velocity() const;
 *     void   solve();
 *     void   project_back_saturation();
 *     void   output_results() const;
 * 
 *     const unsigned int degree;
 * 
 *     Triangulation<dim> triangulation;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     BlockSparsityPattern      sparsity_pattern;
 *     BlockSparseMatrix<double> system_matrix;
 * 
 *     const unsigned int n_refinement_steps;
 * 
 *     DiscreteTime time;
 *     double       viscosity;
 * 
 *     BlockVector<double> solution;
 *     BlockVector<double> old_solution;
 *     BlockVector<double> system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 
 * <a name="Pressurerighthandside"></a> 
 * <h4>Pressure right hand side</h4>
 * 

 * 
 * At present, the right hand side of the pressure equation is simply the
 * zero function. However, the rest of the program is fully equipped to deal
 * with anything else, if this is desired:
 * 
 * @code
 *   template <int dim>
 *   class PressureRightHandSide : public Function<dim>
 *   {
 *   public:
 *     PressureRightHandSide()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Pressureboundaryvalues"></a> 
 * <h4>Pressure boundary values</h4>
 * 

 * 
 * The next are pressure boundary values. As mentioned in the introduction,
 * we choose a linear pressure field:
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
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       return 1 - p[0];
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Saturationboundaryvalues"></a> 
 * <h4>Saturation boundary values</h4>
 * 

 * 
 * Then we also need boundary values on the inflow portions of the
 * boundary. The question whether something is an inflow part is decided
 * when assembling the right hand side, we only have to provide a functional
 * description of the boundary values. This is as explained in the
 * introduction:
 * 
 * @code
 *   template <int dim>
 *   class SaturationBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     SaturationBoundaryValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       if (p[0] == 0)
 *         return 1;
 *       else
 *         return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Initialdata"></a> 
 * <h4>Initial data</h4>
 * 

 * 
 * Finally, we need initial data. In reality, we only need initial data for
 * the saturation, but we are lazy, so we will later, before the first time
 * step, simply interpolate the entire solution for the previous time step
 * from a function that contains all vector components.
 *   

 * 
 * We therefore simply create a function that returns zero in all
 * components. We do that by simply forward every function to the
 * Functions::ZeroFunction class. Why not use that right away in the places of
 * this program where we presently use the <code>InitialValues</code> class?
 * Because this way it is simpler to later go back and choose a different
 * function for initial values.
 * 
 * @code
 *   template <int dim>
 *   class InitialValues : public Function<dim>
 *   {
 *   public:
 *     InitialValues()
 *       : Function<dim>(dim + 2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override
 *     {
 *       return Functions::ZeroFunction<dim>(dim + 2).value(p, component);
 *     }
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override
 *     {
 *       Functions::ZeroFunction<dim>(dim + 2).vector_value(p, values);
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversepermeabilitytensor"></a> 
 * <h3>The inverse permeability tensor</h3>
 * 

 * 
 * As announced in the introduction, we implement two different permeability
 * tensor fields. Each of them we put into a namespace of its own, so that
 * it will be easy later to replace use of one by the other in the code.
 * 

 * 
 * 
 * <a name="Singlecurvingcrackpermeability"></a> 
 * <h4>Single curving crack permeability</h4>
 * 

 * 
 * The first function for the permeability was the one that models a single
 * curving crack. It was already used at the end of step-20, and its
 * functional form is given in the introduction of the present tutorial
 * program. As in some previous programs, we have to declare a (seemingly
 * unnecessary) default constructor of the KInverse class to avoid warnings
 * from some compilers:
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
 *                  std::vector<Tensor<2, dim>> &  values) const override
 *       {
 *         Assert(points.size() == values.size(),
 *                ExcDimensionMismatch(points.size(), values.size()));
 * 
 *         for (unsigned int p = 0; p < points.size(); ++p)
 *           {
 *             values[p].clear();
 * 
 *             const double distance_to_flowline =
 *               std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0]));
 * 
 *             const double permeability =
 *               std::max(std::exp(-(distance_to_flowline * distance_to_flowline) /
 *                                 (0.1 * 0.1)),
 *                        0.01);
 * 
 *             for (unsigned int d = 0; d < dim; ++d)
 *               values[p][d][d] = 1. / permeability;
 *           }
 *       }
 *     };
 *   } // namespace SingleCurvingCrack
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Randommediumpermeability"></a> 
 * <h4>Random medium permeability</h4>
 * 

 * 
 * This function does as announced in the introduction, i.e. it creates an
 * overlay of exponentials at random places. There is one thing worth
 * considering for this class. The issue centers around the problem that the
 * class creates the centers of the exponentials using a random function. If
 * we therefore created the centers each time we create an object of the
 * present type, we would get a different list of centers each time. That's
 * not what we expect from classes of this type: they should reliably
 * represent the same function.
 *   

 * 
 * The solution to this problem is to make the list of centers a static
 * member variable of this class, i.e. there exists exactly one such
 * variable for the entire program, rather than for each object of this
 * type. That's exactly what we are going to do.
 *   

 * 
 * The next problem, however, is that we need a way to initialize this
 * variable. Since this variable is initialized at the beginning of the
 * program, we can't use a regular member function for that since there may
 * not be an object of this type around at the time. The C++ standard
 * therefore says that only non-member and static member functions can be
 * used to initialize a static variable. We use the latter possibility by
 * defining a function <code>get_centers</code> that computes the list of
 * center points when called.
 *   

 * 
 * Note that this class works just fine in both 2d and 3d, with the only
 * difference being that we use more points in 3d: by experimenting we find
 * that we need more exponentials in 3d than in 2d (we have more ground to
 * cover, after all, if we want to keep the distance between centers roughly
 * equal), so we choose 40 in 2d and 100 in 3d. For any other dimension, the
 * function does presently not know what to do so simply throws an exception
 * indicating exactly this.
 * 
 * @code
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
 *                  std::vector<Tensor<2, dim>> &  values) const override
 *       {
 *         Assert(points.size() == values.size(),
 *                ExcDimensionMismatch(points.size(), values.size()));
 * 
 *         for (unsigned int p = 0; p < points.size(); ++p)
 *           {
 *             values[p].clear();
 * 
 *             double permeability = 0;
 *             for (unsigned int i = 0; i < centers.size(); ++i)
 *               permeability += std::exp(-(points[p] - centers[i]).norm_square() /
 *                                        (0.05 * 0.05));
 * 
 *             const double normalized_permeability =
 *               std::min(std::max(permeability, 0.01), 4.);
 * 
 *             for (unsigned int d = 0; d < dim; ++d)
 *               values[p][d][d] = 1. / normalized_permeability;
 *           }
 *       }
 * 
 *     private:
 *       static std::vector<Point<dim>> centers;
 * 
 *       static std::vector<Point<dim>> get_centers()
 *       {
 *         const unsigned int N =
 *           (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented()));
 * 
 *         std::vector<Point<dim>> centers_list(N);
 *         for (unsigned int i = 0; i < N; ++i)
 *           for (unsigned int d = 0; d < dim; ++d)
 *             centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX;
 * 
 *         return centers_list;
 *       }
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     std::vector<Point<dim>>
 *       KInverse<dim>::centers = KInverse<dim>::get_centers();
 *   } // namespace RandomMedium
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversemobilityandsaturationfunctions"></a> 
 * <h3>The inverse mobility and saturation functions</h3>
 * 

 * 
 * There are two more pieces of data that we need to describe, namely the
 * inverse mobility function and the saturation curve. Their form is also
 * given in the introduction:
 * 
 * @code
 *   double mobility_inverse(const double S, const double viscosity)
 *   {
 *     return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S));
 *   }
 * 
 *   double fractional_flow(const double S, const double viscosity)
 *   {
 *     return S * S / (S * S + viscosity * (1 - S) * (1 - S));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * The linear solvers we use are also completely analogous to the ones used
 * in step-20. The following classes are therefore copied verbatim from
 * there. Note that the classes here are not only copied from
 * step-20, but also duplicate classes in deal.II. In a future version of this
 * example, they should be replaced by an efficient method, though. There is a
 * single change: if the size of a linear system is small, i.e. when the mesh
 * is very coarse, then it is sometimes not sufficient to set a maximum of
 * <code>src.size()</code> CG iterations before the solver in the
 * <code>vmult()</code> function converges. (This is, of course, a result of
 * numerical round-off, since we know that on paper, the CG method converges
 * in at most <code>src.size()</code> steps.) As a consequence, we set the
 * maximum number of iterations equal to the maximum of the size of the linear
 * system and 200.
 * 
 * @code
 *   template <class MatrixType>
 *   class InverseMatrix : public Subscriptor
 *   {
 *   public:
 *     InverseMatrix(const MatrixType &m)
 *       : matrix(&m)
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       SolverControl solver_control(std::max<unsigned int>(src.size(), 200),
 *                                    1e-8 * src.l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       dst = 0;
 * 
 *       cg.solve(*matrix, dst, src, PreconditionIdentity());
 *     }
 * 
 *   private:
 *     const SmartPointer<const MatrixType> matrix;
 *   };
 * 
 * 
 * 
 *   class SchurComplement : public Subscriptor
 *   {
 *   public:
 *     SchurComplement(const BlockSparseMatrix<double> &          A,
 *                     const InverseMatrix<SparseMatrix<double>> &Minv)
 *       : system_matrix(&A)
 *       , m_inverse(&Minv)
 *       , tmp1(A.block(0, 0).m())
 *       , tmp2(A.block(0, 0).m())
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       system_matrix->block(0, 1).vmult(tmp1, src);
 *       m_inverse->vmult(tmp2, tmp1);
 *       system_matrix->block(1, 0).vmult(dst, tmp2);
 *     }
 * 
 *   private:
 *     const SmartPointer<const BlockSparseMatrix<double>>           system_matrix;
 *     const SmartPointer<const InverseMatrix<SparseMatrix<double>>> m_inverse;
 * 
 *     mutable Vector<double> tmp1, tmp2;
 *   };
 * 
 * 
 * 
 *   class ApproximateSchurComplement : public Subscriptor
 *   {
 *   public:
 *     ApproximateSchurComplement(const BlockSparseMatrix<double> &A)
 *       : system_matrix(&A)
 *       , tmp1(A.block(0, 0).m())
 *       , tmp2(A.block(0, 0).m())
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       system_matrix->block(0, 1).vmult(tmp1, src);
 *       system_matrix->block(0, 0).precondition_Jacobi(tmp2, tmp1);
 *       system_matrix->block(1, 0).vmult(dst, tmp2);
 *     }
 * 
 *   private:
 *     const SmartPointer<const BlockSparseMatrix<double>> system_matrix;
 * 
 *     mutable Vector<double> tmp1, tmp2;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeTwoPhaseFlowProblemcodeclassimplementation"></a> 
 * <h3><code>TwoPhaseFlowProblem</code> class implementation</h3>
 * 

 * 
 * Here now the implementation of the main class. Much of it is actually
 * copied from step-20, so we won't comment on it in much detail. You should
 * try to get familiar with that program first, then most of what is
 * happening here should be mostly clear.
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemTwoPhaseFlowProblem"></a> 
 * <h4>TwoPhaseFlowProblem::TwoPhaseFlowProblem</h4>
 * 

 * 
 * First for the constructor. We use $RT_k \times DQ_k \times DQ_k$
 * spaces. For initializing the DiscreteTime object, we don't set the time
 * step size in the constructor because we don't have its value yet.
 * The time step size is initially set to zero, but it will be computed
 * before it is needed to increment time, as described in a subsection of
 * the introduction. The time object internally prevents itself from being
 * incremented when $dt = 0$, forcing us to set a non-zero desired size for
 * $dt$ before advancing time.
 * 
 * @code
 *   template <int dim>
 *   TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree)
 *     : degree(degree)
 *     , fe(FE_RaviartThomas<dim>(degree),
 *          1,
 *          FE_DGQ<dim>(degree),
 *          1,
 *          FE_DGQ<dim>(degree),
 *          1)
 *     , dof_handler(triangulation)
 *     , n_refinement_steps(5)
 *     , time(/*start time*/ 0., /*end time*/ 1.)
 *     , viscosity(0.2)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemmake_grid_and_dofs"></a> 
 * <h4>TwoPhaseFlowProblem::make_grid_and_dofs</h4>
 * 

 * 
 * This next function starts out with well-known functions calls that create
 * and refine a mesh, and then associate degrees of freedom with it. It does
 * all the same things as in step-20, just now for three components instead
 * of two.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0, 1);
 *     triangulation.refine_global(n_refinement_steps);
 * 
 *     dof_handler.distribute_dofs(fe);
 *     DoFRenumbering::component_wise(dof_handler);
 * 
 *     const std::vector<types::global_dof_index> dofs_per_component =
 *       DoFTools::count_dofs_per_fe_component(dof_handler);
 *     const unsigned int n_u = dofs_per_component[0],
 *                        n_p = dofs_per_component[dim],
 *                        n_s = dofs_per_component[dim + 1];
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (" << n_u << '+' << n_p << '+' << n_s << ')' << std::endl
 *               << std::endl;
 * 
 *     const unsigned int n_couplings = dof_handler.max_couplings_between_dofs();
 * 
 *     sparsity_pattern.reinit(3, 3);
 *     sparsity_pattern.block(0, 0).reinit(n_u, n_u, n_couplings);
 *     sparsity_pattern.block(1, 0).reinit(n_p, n_u, n_couplings);
 *     sparsity_pattern.block(2, 0).reinit(n_s, n_u, n_couplings);
 *     sparsity_pattern.block(0, 1).reinit(n_u, n_p, n_couplings);
 *     sparsity_pattern.block(1, 1).reinit(n_p, n_p, n_couplings);
 *     sparsity_pattern.block(2, 1).reinit(n_s, n_p, n_couplings);
 *     sparsity_pattern.block(0, 2).reinit(n_u, n_s, n_couplings);
 *     sparsity_pattern.block(1, 2).reinit(n_p, n_s, n_couplings);
 *     sparsity_pattern.block(2, 2).reinit(n_s, n_s, n_couplings);
 * 
 *     sparsity_pattern.collect_sizes();
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
 *     sparsity_pattern.compress();
 * 
 * 
 *     system_matrix.reinit(sparsity_pattern);
 * 
 * 
 *     solution.reinit(3);
 *     solution.block(0).reinit(n_u);
 *     solution.block(1).reinit(n_p);
 *     solution.block(2).reinit(n_s);
 *     solution.collect_sizes();
 * 
 *     old_solution.reinit(3);
 *     old_solution.block(0).reinit(n_u);
 *     old_solution.block(1).reinit(n_p);
 *     old_solution.block(2).reinit(n_s);
 *     old_solution.collect_sizes();
 * 
 *     system_rhs.reinit(3);
 *     system_rhs.block(0).reinit(n_u);
 *     system_rhs.block(1).reinit(n_p);
 *     system_rhs.block(2).reinit(n_s);
 *     system_rhs.collect_sizes();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemassemble_system"></a> 
 * <h4>TwoPhaseFlowProblem::assemble_system</h4>
 * 

 * 
 * This is the function that assembles the linear system, or at least
 * everything except the (1,3) block that depends on the still-unknown
 * velocity computed during this time step (we deal with this in
 * <code>assemble_rhs_S</code>). Much of it is again as in step-20, but we
 * have to deal with some nonlinearity this time.  However, the top of the
 * function is pretty much as usual (note that we set matrix and right hand
 * side to zero at the beginning &mdash; something we didn't have to do for
 * stationary problems since there we use each matrix object only once and
 * it is empty at the beginning anyway).
 *   

 * 
 * Note that in its present form, the function uses the permeability
 * implemented in the RandomMedium::KInverse class. Switching to the single
 * curved crack permeability function is as simple as just changing the
 * namespace name.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     QGauss<dim>     quadrature_formula(degree + 2);
 *     QGauss<dim - 1> face_quadrature_formula(degree + 2);
 * 
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const PressureRightHandSide<dim>  pressure_right_hand_side;
 *     const PressureBoundaryValues<dim> pressure_boundary_values;
 *     const RandomMedium::KInverse<dim> k_inverse;
 * 
 *     std::vector<double>         pressure_rhs_values(n_q_points);
 *     std::vector<double>         boundary_values(n_face_q_points);
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 *     std::vector<Vector<double>>              old_solution_values(n_q_points,
 *                                                                  Vector<double>(dim + 2));
 *     std::vector<std::vector<Tensor<1, dim>>> old_solution_grads(
 *       n_q_points, std::vector<Tensor<1, dim>>(dim + 2));
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 *     const FEValuesExtractors::Scalar saturation(dim + 1);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 * @endcode
 * 
 * Here's the first significant difference: We have to get the values
 * of the saturation function of the previous time step at the
 * quadrature points. To this end, we can use the
 * FEValues::get_function_values (previously already used in step-9,
 * step-14 and step-15), a function that takes a solution vector and
 * returns a list of function values at the quadrature points of the
 * present cell. In fact, it returns the complete vector-valued
 * solution at each quadrature point, i.e. not only the saturation but
 * also the velocities and pressure:
 * 
 * @code
 *         fe_values.get_function_values(old_solution, old_solution_values);
 * 
 * @endcode
 * 
 * Then we also have to get the values of the pressure right hand side
 * and of the inverse permeability tensor at the quadrature points:
 * 
 * @code
 *         pressure_right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                             pressure_rhs_values);
 *         k_inverse.value_list(fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 * @endcode
 * 
 * With all this, we can now loop over all the quadrature points and
 * shape functions on this cell and assemble those parts of the matrix
 * and right hand side that we deal with in this function. The
 * individual terms in the contributions should be self-explanatory
 * given the explicit form of the bilinear form stated in the
 * introduction:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const double old_s = old_solution_values[q](dim + 1);
 * 
 *               const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q);
 *               const double div_phi_i_u = fe_values[velocities].divergence(i, q);
 *               const double phi_i_p     = fe_values[pressure].value(i, q);
 *               const double phi_i_s     = fe_values[saturation].value(i, q);
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const Tensor<1, dim> phi_j_u =
 *                     fe_values[velocities].value(j, q);
 *                   const double div_phi_j_u =
 *                     fe_values[velocities].divergence(j, q);
 *                   const double phi_j_p = fe_values[pressure].value(j, q);
 *                   const double phi_j_s = fe_values[saturation].value(j, q);
 * 
 *                   local_matrix(i, j) +=
 *                     (phi_i_u * k_inverse_values[q] *
 *                        mobility_inverse(old_s, viscosity) * phi_j_u -
 *                      div_phi_i_u * phi_j_p - phi_i_p * div_phi_j_u +
 *                      phi_i_s * phi_j_s) *
 *                     fe_values.JxW(q);
 *                 }
 * 
 *               local_rhs(i) +=
 *                 (-phi_i_p * pressure_rhs_values[q]) * fe_values.JxW(q);
 *             }
 * 
 * 
 * @endcode
 * 
 * Next, we also have to deal with the pressure boundary values. This,
 * again is as in step-20:
 * 
 * @code
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary())
 *             {
 *               fe_face_values.reinit(cell, face);
 * 
 *               pressure_boundary_values.value_list(
 *                 fe_face_values.get_quadrature_points(), boundary_values);
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   {
 *                     const Tensor<1, dim> phi_i_u =
 *                       fe_face_values[velocities].value(i, q);
 * 
 *                     local_rhs(i) +=
 *                       -(phi_i_u * fe_face_values.normal_vector(q) *
 *                         boundary_values[q] * fe_face_values.JxW(q));
 *                   }
 *             }
 * 
 * @endcode
 * 
 * The final step in the loop over all cells is to transfer local
 * contributions into the global matrix and right hand side vector:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               local_matrix(i, j));
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           system_rhs(local_dof_indices[i]) += local_rhs(i);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * So much for assembly of matrix and right hand side. Note that we do not
 * have to interpolate and apply boundary values since they have all been
 * taken care of in the weak form already.
 * 

 * 
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemassemble_rhs_S"></a> 
 * <h4>TwoPhaseFlowProblem::assemble_rhs_S</h4>
 * 

 * 
 * As explained in the introduction, we can only evaluate the right hand
 * side of the saturation equation once the velocity has been computed. We
 * therefore have this separate function to this end.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_rhs_S()
 *   {
 *     QGauss<dim>       quadrature_formula(degree + 2);
 *     QGauss<dim - 1>   face_quadrature_formula(degree + 2);
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 *     FEFaceValues<dim> fe_face_values_neighbor(fe,
 *                                               face_quadrature_formula,
 *                                               update_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     Vector<double> local_rhs(dofs_per_cell);
 * 
 *     std::vector<Vector<double>> old_solution_values(n_q_points,
 *                                                     Vector<double>(dim + 2));
 *     std::vector<Vector<double>> old_solution_values_face(n_face_q_points,
 *                                                          Vector<double>(dim +
 *                                                                         2));
 *     std::vector<Vector<double>> old_solution_values_face_neighbor(
 *       n_face_q_points, Vector<double>(dim + 2));
 *     std::vector<Vector<double>> present_solution_values(n_q_points,
 *                                                         Vector<double>(dim +
 *                                                                        2));
 *     std::vector<Vector<double>> present_solution_values_face(
 *       n_face_q_points, Vector<double>(dim + 2));
 * 
 *     std::vector<double>                  neighbor_saturation(n_face_q_points);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     SaturationBoundaryValues<dim> saturation_boundary_values;
 * 
 *     const FEValuesExtractors::Scalar saturation(dim + 1);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         local_rhs = 0;
 *         fe_values.reinit(cell);
 * 
 *         fe_values.get_function_values(old_solution, old_solution_values);
 *         fe_values.get_function_values(solution, present_solution_values);
 * 
 * @endcode
 * 
 * First for the cell terms. These are, following the formulas in the
 * introduction, $(S^n,\sigma)-(F(S^n) \mathbf{v}^{n+1},\nabla
 * \sigma)$, where $\sigma$ is the saturation component of the test
 * function:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const double   old_s = old_solution_values[q](dim + 1);
 *               Tensor<1, dim> present_u;
 *               for (unsigned int d = 0; d < dim; ++d)
 *                 present_u[d] = present_solution_values[q](d);
 * 
 *               const double         phi_i_s = fe_values[saturation].value(i, q);
 *               const Tensor<1, dim> grad_phi_i_s =
 *                 fe_values[saturation].gradient(i, q);
 * 
 *               local_rhs(i) +=
 *                 (time.get_next_step_size() * fractional_flow(old_s, viscosity) *
 *                    present_u * grad_phi_i_s +
 *                  old_s * phi_i_s) *
 *                 fe_values.JxW(q);
 *             }
 * 
 * @endcode
 * 
 * Secondly, we have to deal with the flux parts on the face
 * boundaries. This was a bit more involved because we first have to
 * determine which are the influx and outflux parts of the cell
 * boundary. If we have an influx boundary, we need to evaluate the
 * saturation on the other side of the face (or the boundary values,
 * if we are at the boundary of the domain).
 *         

 * 
 * All this is a bit tricky, but has been explained in some detail
 * already in step-9. Take a look there how this is supposed to work!
 * 
 * @code
 *         for (const auto face_no : cell->face_indices())
 *           {
 *             fe_face_values.reinit(cell, face_no);
 * 
 *             fe_face_values.get_function_values(old_solution,
 *                                                old_solution_values_face);
 *             fe_face_values.get_function_values(solution,
 *                                                present_solution_values_face);
 * 
 *             if (cell->at_boundary(face_no))
 *               saturation_boundary_values.value_list(
 *                 fe_face_values.get_quadrature_points(), neighbor_saturation);
 *             else
 *               {
 *                 const auto         neighbor = cell->neighbor(face_no);
 *                 const unsigned int neighbor_face =
 *                   cell->neighbor_of_neighbor(face_no);
 * 
 *                 fe_face_values_neighbor.reinit(neighbor, neighbor_face);
 * 
 *                 fe_face_values_neighbor.get_function_values(
 *                   old_solution, old_solution_values_face_neighbor);
 * 
 *                 for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                   neighbor_saturation[q] =
 *                     old_solution_values_face_neighbor[q](dim + 1);
 *               }
 * 
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 Tensor<1, dim> present_u_face;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   present_u_face[d] = present_solution_values_face[q](d);
 * 
 *                 const double normal_flux =
 *                   present_u_face * fe_face_values.normal_vector(q);
 * 
 *                 const bool is_outflow_q_point = (normal_flux >= 0);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   local_rhs(i) -=
 *                     time.get_next_step_size() * normal_flux *
 *                     fractional_flow((is_outflow_q_point == true ?
 *                                        old_solution_values_face[q](dim + 1) :
 *                                        neighbor_saturation[q]),
 *                                     viscosity) *
 *                     fe_face_values[saturation].value(i, q) *
 *                     fe_face_values.JxW(q);
 *               }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           system_rhs(local_dof_indices[i]) += local_rhs(i);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemsolve"></a> 
 * <h4>TwoPhaseFlowProblem::solve</h4>
 * 

 * 
 * After all these preparations, we finally solve the linear system for
 * velocity and pressure in the same way as in step-20. After that, we have
 * to deal with the saturation equation (see below):
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::solve()
 *   {
 *     const InverseMatrix<SparseMatrix<double>> m_inverse(
 *       system_matrix.block(0, 0));
 *     Vector<double> tmp(solution.block(0).size());
 *     Vector<double> schur_rhs(solution.block(1).size());
 *     Vector<double> tmp2(solution.block(2).size());
 * 
 * 
 * @endcode
 * 
 * First the pressure, using the pressure Schur complement of the first
 * two equations:
 * 
 * @code
 *     {
 *       m_inverse.vmult(tmp, system_rhs.block(0));
 *       system_matrix.block(1, 0).vmult(schur_rhs, tmp);
 *       schur_rhs -= system_rhs.block(1);
 * 
 * 
 *       SchurComplement schur_complement(system_matrix, m_inverse);
 * 
 *       ApproximateSchurComplement approximate_schur_complement(system_matrix);
 * 
 *       InverseMatrix<ApproximateSchurComplement> preconditioner(
 *         approximate_schur_complement);
 * 
 * 
 *       SolverControl            solver_control(solution.block(1).size(),
 *                                    1e-12 * schur_rhs.l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       cg.solve(schur_complement, solution.block(1), schur_rhs, preconditioner);
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " CG Schur complement iterations for pressure." << std::endl;
 *     }
 * 
 * @endcode
 * 
 * Now the velocity:
 * 
 * @code
 *     {
 *       system_matrix.block(0, 1).vmult(tmp, solution.block(1));
 *       tmp *= -1;
 *       tmp += system_rhs.block(0);
 * 
 *       m_inverse.vmult(solution.block(0), tmp);
 *     }
 * 
 * @endcode
 * 
 * Finally, we have to take care of the saturation equation. The first
 * business we have here is to determine the time step using the formula
 * in the introduction. Knowing the shape of our domain and that we
 * created the mesh by regular subdivision of cells, we can compute the
 * diameter of each of our cells quite easily (in fact we use the linear
 * extensions in coordinate directions of the cells, not the
 * diameter). Note that we will learn a more general way to do this in
 * step-24, where we use the GridTools::minimal_cell_diameter function.
 *     

 * 
 * The maximal velocity we compute using a helper function to compute the
 * maximal velocity defined below, and with all this we can evaluate our
 * new time step length. We use the method
 * DiscreteTime::set_desired_next_time_step() to suggest the new
 * calculated value of the time step to the DiscreteTime object. In most
 * cases, the time object uses the exact provided value to increment time.
 * It some case, the step size may be modified further by the time object.
 * For example, if the calculated time increment overshoots the end time,
 * it is truncated accordingly.
 * 
 * @code
 *     time.set_desired_next_step_size(std::pow(0.5, double(n_refinement_steps)) /
 *                                     get_maximal_velocity());
 * 
 * @endcode
 * 
 * The next step is to assemble the right hand side, and then to pass
 * everything on for solution. At the end, we project back saturations
 * onto the physically reasonable range:
 * 
 * @code
 *     assemble_rhs_S();
 *     {
 *       SolverControl            solver_control(system_matrix.block(2, 2).m(),
 *                                    1e-8 * system_rhs.block(2).l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 *       cg.solve(system_matrix.block(2, 2),
 *                solution.block(2),
 *                system_rhs.block(2),
 *                PreconditionIdentity());
 * 
 *       project_back_saturation();
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " CG iterations for saturation." << std::endl;
 *     }
 * 
 * 
 *     old_solution = solution;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemoutput_results"></a> 
 * <h4>TwoPhaseFlowProblem::output_results</h4>
 * 

 * 
 * There is nothing surprising here. Since the program will do a lot of time
 * steps, we create an output file only every fifth time step and skip all
 * other time steps at the top of the file already.
 *   

 * 
 * When creating file names for output close to the bottom of the function,
 * we convert the number of the time step to a string representation that
 * is padded by leading zeros to four digits. We do this because this way
 * all output file names have the same length, and consequently sort well
 * when creating a directory listing.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::output_results() const
 *   {
 *     if (time.get_step_number() % 5 != 0)
 *       return;
 * 
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 2:
 *           solution_names = {"u", "v", "p", "S"};
 *           break;
 * 
 *         case 3:
 *           solution_names = {"u", "v", "w", "p", "S"};
 *           break;
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, solution_names);
 * 
 *     data_out.build_patches(degree + 1);
 * 
 *     std::ofstream output("solution-" +
 *                          Utilities::int_to_string(time.get_step_number(), 4) +
 *                          ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemproject_back_saturation"></a> 
 * <h4>TwoPhaseFlowProblem::project_back_saturation</h4>
 * 

 * 
 * In this function, we simply run over all saturation degrees of freedom
 * and make sure that if they should have left the physically reasonable
 * range, that they be reset to the interval $[0,1]$. To do this, we only
 * have to loop over all saturation components of the solution vector; these
 * are stored in the block 2 (block 0 are the velocities, block 1 are the
 * pressures).
 *   

 * 
 * It may be instructive to note that this function almost never triggers
 * when the time step is chosen as mentioned in the introduction. However,
 * if we choose the timestep only slightly larger, we get plenty of values
 * outside the proper range. Strictly speaking, the function is therefore
 * unnecessary if we choose the time step small enough. In a sense, the
 * function is therefore only a safety device to avoid situations where our
 * entire solution becomes unphysical because individual degrees of freedom
 * have become unphysical a few time steps earlier.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::project_back_saturation()
 *   {
 *     for (unsigned int i = 0; i < solution.block(2).size(); ++i)
 *       if (solution.block(2)(i) < 0)
 *         solution.block(2)(i) = 0;
 *       else if (solution.block(2)(i) > 1)
 *         solution.block(2)(i) = 1;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemget_maximal_velocity"></a> 
 * <h4>TwoPhaseFlowProblem::get_maximal_velocity</h4>
 * 

 * 
 * The following function is used in determining the maximal allowable time
 * step. What it does is to loop over all quadrature points in the domain
 * and find what the maximal magnitude of the velocity is.
 * 
 * @code
 *   template <int dim>
 *   double TwoPhaseFlowProblem<dim>::get_maximal_velocity() const
 *   {
 *     QGauss<dim>        quadrature_formula(degree + 2);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(fe, quadrature_formula, update_values);
 *     std::vector<Vector<double>> solution_values(n_q_points,
 *                                                 Vector<double>(dim + 2));
 *     double                      max_velocity = 0;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         fe_values.get_function_values(solution, solution_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             Tensor<1, dim> velocity;
 *             for (unsigned int i = 0; i < dim; ++i)
 *               velocity[i] = solution_values[q](i);
 * 
 *             max_velocity = std::max(max_velocity, velocity.norm());
 *           }
 *       }
 * 
 *     return max_velocity;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemrun"></a> 
 * <h4>TwoPhaseFlowProblem::run</h4>
 * 

 * 
 * This is the final function of our main class. Its brevity speaks for
 * itself. There are only two points worth noting: First, the function
 * projects the initial values onto the finite element space at the
 * beginning; the VectorTools::project function doing this requires an
 * argument indicating the hanging node constraints. We have none in this
 * program (we compute on a uniformly refined mesh), but the function
 * requires the argument anyway, of course. So we have to create a
 * constraint object. In its original state, constraint objects are
 * unsorted, and have to be sorted (using the AffineConstraints::close
 * function) before they can be used. This is what we do here, and which is
 * why we can't simply call the VectorTools::project function with an
 * anonymous temporary object <code>AffineConstraints<double>()</code> as the
 * second argument.
 *   

 * 
 * The second point worth mentioning is that we only compute the length of
 * the present time step in the middle of solving the linear system
 * corresponding to each time step. We can therefore output the present
 * time of a time step only at the end of the time step.
 * We increment time by calling the method DiscreteTime::advance_time()
 * inside the loop. Since we are reporting the time and dt after we
 * increment it, we have to call the method
 * DiscreteTime::get_previous_step_size() instead of
 * DiscreteTime::get_next_step_size(). After many steps, when the simulation
 * reaches the end time, the last dt is chosen by the DiscreteTime class in
 * such a way that the last step finishes exactly at the end time.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::run()
 *   {
 *     make_grid_and_dofs();
 * 
 *     {
 *       AffineConstraints<double> constraints;
 *       constraints.close();
 * 
 *       VectorTools::project(dof_handler,
 *                            constraints,
 *                            QGauss<dim>(degree + 2),
 *                            InitialValues<dim>(),
 *                            old_solution);
 *     }
 * 
 *     do
 *       {
 *         std::cout << "Timestep " << time.get_step_number() + 1 << std::endl;
 * 
 *         assemble_system();
 * 
 *         solve();
 * 
 *         output_results();
 * 
 *         time.advance_time();
 *         std::cout << "   Now at t=" << time.get_current_time()
 *                   << ", dt=" << time.get_previous_step_size() << '.'
 *                   << std::endl
 *                   << std::endl;
 *       }
 *     while (time.is_at_end() == false);
 *   }
 * } // namespace Step21
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * That's it. In the main function, we pass the degree of the finite element
 * space to the constructor of the TwoPhaseFlowProblem object.  Here, we use
 * zero-th degree elements, i.e. $RT_0\times DQ_0 \times DQ_0$. The rest is as
 * in all the other programs.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step21;
 * 
 *       TwoPhaseFlowProblem<2> two_phase_flow_problem(0);
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


The code as presented here does not actually compute the results
found on the web page. The reason is, that even on a decent
computer it runs more than a day. If you want to reproduce these
results, modify the end time of the DiscreteTime object to `250` within the
constructor of TwoPhaseFlowProblem.

If we run the program, we get the following kind of output:
@code
Number of active cells: 1024
Number of degrees of freedom: 4160 (2112+1024+1024)

Timestep 1
   22 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0326742, dt=0.0326742.

Timestep 2
   17 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0653816, dt=0.0327074.

Timestep 3
   17 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0980651, dt=0.0326836.

...
@endcode
As we can see, the time step is pretty much constant right from the start,
which indicates that the velocities in the domain are not strongly dependent
on changes in saturation, although they certainly are through the factor
$\lambda(S)$ in the pressure equation.

Our second observation is that the number of CG iterations needed to solve the
pressure Schur complement equation drops from 22 to 17 between the first and
the second time step (in fact, it remains around 17 for the rest of the
computations). The reason is actually simple: Before we solve for the pressure
during a time step, we don't reset the <code>solution</code> variable to
zero. The pressure (and the other variables) therefore have the previous time
step's values at the time we get into the CG solver. Since the velocities and
pressures don't change very much as computations progress, the previous time
step's pressure is actually a good initial guess for this time step's
pressure. Consequently, the number of iterations we need once we have computed
the pressure once is significantly reduced.

The final observation concerns the number of iterations needed to solve for
the saturation, i.e. one. This shouldn't surprise us too much: the matrix we
have to solve with is the mass matrix. However, this is the mass matrix for
the $DGQ_0$ element of piecewise constants where no element couples with the
degrees of freedom on neighboring cells. The matrix is therefore a diagonal
one, and it is clear that we should be able to invert this matrix in a single
CG iteration.


With all this, here are a few movies that show how the saturation progresses
over time. First, this is for the single crack model, as implemented in the
<code>SingleCurvingCrack::KInverse</code> class:

<img src="https://www.dealii.org/images/steps/developer/step-21.centerline.gif" alt="">

As can be seen, the water rich fluid snakes its way mostly along the
high-permeability zone in the middle of the domain, whereas the rest of the
domain is mostly impermeable. This and the next movie are generated using
<code>n_refinement_steps=7</code>, leading to a $128\times 128$ mesh with some
16,000 cells and about 66,000 unknowns in total.


The second movie shows the saturation for the random medium model of class
<code>RandomMedium::KInverse</code>, where we have randomly distributed
centers of high permeability and fluid hops from one of these zones to
the next:

<img src="https://www.dealii.org/images/steps/developer/step-21.random2d.gif" alt="">


Finally, here is the same situation in three space dimensions, on a mesh with
<code>n_refinement_steps=5</code>, which produces a mesh of some 32,000 cells
and 167,000 degrees of freedom:

<img src="https://www.dealii.org/images/steps/developer/step-21.random3d.gif" alt="">

To repeat these computations, all you have to do is to change the line
@code
      TwoPhaseFlowProblem<2> two_phase_flow_problem(0);
@endcode
in the main function to
@code
      TwoPhaseFlowProblem<3> two_phase_flow_problem(0);
@endcode
The visualization uses a cloud technique, where the saturation is indicated by
colored but transparent clouds for each cell. This way, one can also see
somewhat what happens deep inside the domain. A different way of visualizing
would have been to show isosurfaces of the saturation evolving over
time. There are techniques to plot isosurfaces transparently, so that one can
see several of them at the same time like the layers of an onion.

So why don't we show such isosurfaces? The problem lies in the way isosurfaces
are computed: they require that the field to be visualized is continuous, so
that the isosurfaces can be generated by following contours at least across a
single cell. However, our saturation field is piecewise constant and
discontinuous. If we wanted to plot an isosurface for a saturation $S=0.5$,
chances would be that there is no single point in the domain where that
saturation is actually attained. If we had to define isosurfaces in that
context at all, we would have to take the interfaces between cells, where one
of the two adjacent cells has a saturation greater than and the other cell a
saturation less than 0.5. However, it appears that most visualization programs
are not equipped to do this kind of transformation.


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


There are a number of areas where this program can be improved. Three of them
are listed below. All of them are, in fact, addressed in a tutorial program
that forms the continuation of the current one: step-43.


<a name="Solvers"></a><h4>Solvers</h4>


At present, the program is not particularly fast: the 2d random medium
computation took about a day for the 1,000 or so time steps. The corresponding
3d computation took almost two days for 800 time steps. The reason why it
isn't faster than this is twofold. First, we rebuild the entire matrix in
every time step, although some parts such as the $B$, $B^T$, and $M^S$ blocks
never change.

Second, we could do a lot better with the solver and
preconditioners. Presently, we solve the Schur complement $B^TM^u(S)^{-1}B$
with a CG method, using $[B^T (\textrm{diag}(M^u(S)))^{-1} B]^{-1}$ as a
preconditioner. Applying this preconditioner is expensive, since it involves
solving a linear system each time. This may have been appropriate for @ref
step_20 "step-20", where we have to solve the entire problem only
once. However, here we have to solve it hundreds of times, and in such cases
it is worth considering a preconditioner that is more expensive to set up the
first time, but cheaper to apply later on.

One possibility would be to realize that the matrix we use as preconditioner,
$B^T (\textrm{diag}(M^u(S)))^{-1} B$ is still sparse, and symmetric on top of
that. If one looks at the flow field evolve over time, we also see that while
$S$ changes significantly over time, the pressure hardly does and consequently
$B^T (\textrm{diag}(M^u(S)))^{-1} B \approx B^T (\textrm{diag}(M^u(S^0)))^{-1}
B$. In other words, the matrix for the first time step should be a good
preconditioner also for all later time steps.  With a bit of
back-and-forthing, it isn't hard to actually get a representation of it as a
SparseMatrix object. We could then hand it off to the SparseMIC class to form
a sparse incomplete Cholesky decomposition. To form this decomposition is
expensive, but we have to do it only once in the first time step, and can then
use it as a cheap preconditioner in the future. We could do better even by
using the SparseDirectUMFPACK class that produces not only an incomplete, but
a complete decomposition of the matrix, which should yield an even better
preconditioner.

Finally, why use the approximation $B^T (\textrm{diag}(M^u(S)))^{-1} B$ to
precondition $B^T M^u(S)^{-1} B$? The latter matrix, after all, is the mixed
form of the Laplace operator on the pressure space, for which we use linear
elements. We could therefore build a separate matrix $A^p$ on the side that
directly corresponds to the non-mixed formulation of the Laplacian, for
example using the bilinear form $(\mathbf{K}\lambda(S^n) \nabla
\varphi_i,\nabla\varphi_j)$. We could then form an incomplete or complete
decomposition of this non-mixed matrix and use it as a preconditioner of the
mixed form.

Using such techniques, it can reasonably be expected that the solution process
will be faster by at least an order of magnitude.


<a name="Timestepping"></a><h4>Time stepping</h4>


In the introduction we have identified the time step restriction
@f[
  \triangle t_{n+1} \le \frac h{|\mathbf{u}^{n+1}(\mathbf{x})|}
@f]
that has to hold globally, i.e. for all $\mathbf x$. After discretization, we
satisfy it by choosing
@f[
  \triangle t_{n+1} = \frac {\min_K h_K}{\max_{\mathbf{x}}|\mathbf{u}^{n+1}(\mathbf{x})|}.
@f]

This restriction on the time step is somewhat annoying: the finer we make the
mesh the smaller the time step; in other words, we get punished twice: each
time step is more expensive to solve and we have to do more time steps.

This is particularly annoying since the majority of the additional work is
spent solving the implicit part of the equations, i.e. the pressure-velocity
system, whereas it is the hyperbolic transport equation for the saturation
that imposes the time step restriction.

To avoid this bottleneck, people have invented a number of approaches. For
example, they may only re-compute the pressure-velocity field every few time
steps (or, if you want, use different time step sizes for the
pressure/velocity and saturation equations). This keeps the time step
restriction on the cheap explicit part while it makes the solution of the
implicit part less frequent. Experiments in this direction are
certainly worthwhile; one starting point for such an approach is the paper by
Zhangxin Chen, Guanren Huan and Baoyan Li: <i>An improved IMPES method for
two-phase flow in porous media</i>, Transport in Porous Media, 54 (2004),
pp. 361&mdash;376. There are certainly many other papers on this topic as well, but
this one happened to land on our desk a while back.



<a name="Adaptivity"></a><h4>Adaptivity</h4>


Adaptivity would also clearly help. Looking at the movies, one clearly sees
that most of the action is confined to a relatively small part of the domain
(this particularly obvious for the saturation, but also holds for the
velocities and pressures). Adaptivity can therefore be expected to keep the
necessary number of degrees of freedom low, or alternatively increase the
accuracy.

On the other hand, adaptivity for time dependent problems is not a trivial
thing: we would have to change the mesh every few time steps, and we would
have to transport our present solution to the next mesh every time we change
it (something that the SolutionTransfer class can help with). These are not
insurmountable obstacles, but they do require some additional coding and more
than we felt comfortable was worth packing into this tutorial program.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-21.cc"
*/
