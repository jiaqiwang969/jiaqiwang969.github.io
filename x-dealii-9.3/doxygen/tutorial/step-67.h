/**
@page step_67 The step-67 tutorial program
This tutorial depends on step-33, step-48, step-59.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#TheEulerequations">The Euler equations</a>
        <li><a href="#HighorderdiscontinuousGalerkindiscretization">High-order discontinuous Galerkin discretization</a>
        <li><a href="#Explicittimeintegration">Explicit time integration</a>
        <li><a href="#Fastevaluationofintegralsbymatrixfreetechniques">Fast evaluation of integrals by matrix-free techniques</a>
        <li><a href="#Evaluationoftheinversemassmatrixwithmatrixfreetechniques">Evaluation of the inverse mass matrix with matrix-free techniques</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#LowstorageexplicitRungeKuttatimeintegrators">Low-storage explicit Runge&mdash;Kutta time integrators</a>
        <li><a href="#ImplementationofpointwiseoperationsoftheEulerequations">Implementation of point-wise operations of the Euler equations</a>
        <li><a href="#TheEulerOperationclass">The EulerOperation class</a>
      <ul>
        <li><a href="#Localevaluators">Local evaluators</a>
        <li><a href="#Theapplyandrelatedfunctions">The apply() and related functions</a>
      </ul>
        <li><a href="#TheEulerProblemclass">The EulerProblem class</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Convergenceratesfortheanalyticaltestcase">Convergence rates for the analytical test case</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin2D">Results for flow in channel around cylinder in 2D</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin3D">Results for flow in channel around cylinder in 3D</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Moreadvancednumericalfluxfunctionsandskewsymmetricformulations">More advanced numerical flux functions and skew-symmetric formulations</a>
        <li><a href="#Equippingthecodeforsupersoniccalculations">Equipping the code for supersonic calculations</a>
        <li><a href="#ExtensiontothelinearizedEulerequations">Extension to the linearized Euler equations</a>
        <li><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly

<br>

<i>
This program was contributed by Martin Kronbichler. Many ideas presented here
are the result of common code development with Niklas Fehn, Katharina Kormann,
Peter Munch, and Svenja Schoeder.

This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA).
</i>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This tutorial program solves the Euler equations of fluid dynamics using an
explicit time integrator with the matrix-free framework applied to a
high-order discontinuous Galerkin discretization in space. For details about
the Euler system and an alternative implicit approach, we also refer to the
step-33 tutorial program. You might also want to look at step-69 for
an alternative approach to solving these equations.


<a name="TheEulerequations"></a><h3>The Euler equations</h3>


The Euler equations are a conservation law, describing the motion of a
compressible inviscid gas,
@f[
\frac{\partial \mathbf{w}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{w}) =
\mathbf{G}(\mathbf w),
@f]
where the $d+2$ components of the solution vector are $\mathbf{w}=(\rho, \rho
u_1,\ldots,\rho u_d,E)^{\mathrm T}$. Here, $\rho$ denotes the fluid density,
${\mathbf u}=(u_1,\ldots, u_d)^\mathrm T$ the fluid velocity, and $E$ the
energy density of the gas. The velocity is not directly solved for, but rather
the variable $\rho \mathbf{u}$, the linear momentum (since this is the
conserved quantity).

The Euler flux function, a $(d+2)\times d$ matrix, is defined as
@f[
  \mathbf F(\mathbf w)
  =
  \begin{pmatrix}
  \rho \mathbf{u}\\
  \rho \mathbf{u} \otimes \mathbf{u} + \mathbb{I}p\\
  (E+p)\mathbf{u}
  \end{pmatrix}
@f]
with $\mathbb{I}$ the $d\times d$ identity matrix and $\otimes$ the outer
product; its components denote the mass, momentum, and energy fluxes, respectively.
The right hand side forcing is given by
@f[
  \mathbf G(\mathbf w)
  =
  \begin{pmatrix}
  0\\
  \rho\mathbf{g}\\
  \rho \mathbf{u} \cdot \mathbf{g}
  \end{pmatrix},
@f]
where the vector $\mathbf g$ denotes the direction and magnitude of
gravity. It could, however, also denote any other external force per unit mass
that is acting on the fluid. (Think, for example, of the electrostatic
forces exerted by an external electric field on charged particles.)

The three blocks of equations, the second involving $d$ components, describe
the conservation of mass, momentum, and energy. The pressure is not a
solution variable but needs to be expressed through a "closure relationship"
by the other variables; we here choose the relationship appropriate
for a gas with molecules composed of two atoms, which at moderate
temperatures is given by $p=(\gamma - 1) \left(E-\frac 12 \rho
\mathbf{u}\cdot \mathbf{u}\right)$ with the constant $\gamma = 1.4$.


<a name="HighorderdiscontinuousGalerkindiscretization"></a><h3>High-order discontinuous Galerkin discretization</h3>


For spatial discretization, we use a high-order discontinuous Galerkin (DG)
discretization, using a solution expansion of the form
@f[
\mathbf{w}_h(\mathbf{x}, t) =
\sum_{j=1}^{n_\mathbf{dofs}} \boldsymbol{\varphi}_j(\mathbf{x}) {w}_j(t).
@f]
Here, $\boldsymbol{\varphi}_j$ denotes the $j$th basis function, written
in vector form with separate shape functions for the different components and
letting $w_j(t)$ go through the density, momentum, and energy variables,
respectively. In this form, the space dependence is contained in the shape
functions and the time dependence in the unknown coefficients $w_j$. As
opposed to the continuous finite element method where some shape functions
span across element boundaries, the shape functions are local to a single
element in DG methods, with a discontinuity from one element to the next. The
connection of the solution from one cell to its neighbors is instead
imposed by the numerical fluxes
specified below. This allows for some additional flexibility, for example to
introduce directionality in the numerical method by, e.g., upwinding.

DG methods are popular methods for solving problems of transport character
because they combine low dispersion errors with controllable dissipation on
barely resolved scales. This makes them particularly attractive for simulation
in the field of fluid dynamics where a wide range of active scales needs to be
represented and inadequately resolved features are prone to disturb the
important well-resolved features. Furthermore, high-order DG methods are
well-suited for modern hardware with the right implementation. At the same
time, DG methods are no silver bullet. In particular when the solution
develops discontinuities (shocks), as is typical for the Euler equations in
some flow regimes, high-order DG methods tend to oscillatory solutions, like
all high-order methods when not using flux- or slope-limiters. This is a consequence of <a
href="https://en.wikipedia.org/wiki/Godunov%27s_theorem">Godunov's theorem</a>
that states that any total variation limited (TVD) scheme that is linear (like
a basic DG discretization) can at most be first-order accurate. Put
differently, since DG methods aim for higher order accuracy, they cannot be
TVD on solutions that develop shocks. Even though some communities claim that
the numerical flux in DG methods can control dissipation, this is of limited
value unless <b>all</b> shocks in a problem align with cell boundaries. Any
shock that passes through the interior of cells will again produce oscillatory
components due to the high-order polynomials. In the finite element and DG
communities, there exist a number of different approaches to deal with shocks,
for example the introduction of artificial diffusion on troubled cells (using
a troubled-cell indicator based e.g. on a modal decomposition of the
solution), a switch to dissipative low-order finite volume methods on a
subgrid, or the addition of some limiting procedures. Given the ample
possibilities in this context, combined with the considerable implementation
effort, we here refrain from the regime of the Euler equations with pronounced
shocks, and rather concentrate on the regime of subsonic flows with wave-like
phenomena. For a method that works well with shocks (but is more expensive per
unknown), we refer to the step-69 tutorial program.

For the derivation of the DG formulation, we multiply the Euler equations with
test functions $\mathbf{v}$ and integrate over an individual cell $K$, which
gives
@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}
+ \left(\mathbf{v}, \nabla \cdot \mathbf{F}(\mathbf{w})\right)_{K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.
@f]

We then integrate the second term by parts, moving the divergence
from the solution slot to the test function slot, and producing an integral
over the element boundary:
@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}
- \left(\nabla \mathbf{v}, \mathbf{F}(\mathbf{w})\right)_{K}
+ \left<\mathbf{v}, \mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})
\right>_{\partial K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.
@f]
In the surface integral, we have replaced the term $\mathbf{F}(\mathbf w)$ by
the term $\widehat{\mathbf{F}}(\mathbf w)$, the numerical flux. The role of
the numerical flux is to connect the solution on neighboring elements and
weakly impose continuity of the solution. This ensures that the global
coupling of the PDE is reflected in the discretization, despite independent
basis functions on the cells. The connectivity to the neighbor is included by
defining the numerical flux as a function $\widehat{\mathbf{F}}(\mathbf w^-,
\mathbf w^+)$ of the solution from both sides of an interior face, $\mathbf
w^-$ and $\mathbf w^+$. A basic property we require is that the numerical flux
needs to be <b>conservative</b>. That is, we want all information (i.e.,
mass, momentum, and energy) that leaves a cell over
a face to enter the neighboring cell in its entirety and vice versa. This can
be expressed as $\widehat{\mathbf{F}}(\mathbf w^-, \mathbf w^+) =
\widehat{\mathbf{F}}(\mathbf w^+, \mathbf w^-)$, meaning that the numerical
flux evaluates to the same result from either side. Combined with the fact
that the numerical flux is multiplied by the unit outer normal vector on the
face under consideration, which points in opposite direction from the two
sides, we see that the conservation is fulfilled. An alternative point of view
of the numerical flux is as a single-valued intermediate state that links the
solution weakly from both sides.

There is a large number of numerical flux functions available, also called
Riemann solvers. For the Euler equations, there exist so-called exact Riemann
solvers -- meaning that the states from both sides are combined in a way that
is consistent with the Euler equations along a discontinuity -- and
approximate Riemann solvers, which violate some physical properties and rely
on other mechanisms to render the scheme accurate overall. Approxiate Riemann
solvers have the advantage of beging cheaper to compute. Most flux functions
have their origin in the finite volume community, which are similar to DG
methods with polynomial degree 0 within the cells (called volumes). As the
volume integral of the Euler operator $\mathbf{F}$ would disappear for
constant solution and test functions, the numerical flux must fully represent
the physical operator, explaining why there has been a large body of research
in that community. For DG methods, consistency is guaranteed by higher order
polynomials within the cells, making the numerical flux less of an issue and
usually affecting only the convergence rate, e.g., whether the solution
converges as $\mathcal O(h^p)$, $\mathcal O(h^{p+1/2})$ or $\mathcal
O(h^{p+1})$ in the $L_2$ norm for polynomials of degree $p$. The numerical
flux can thus be seen as a mechanism to select more advantageous
dissipation/dispersion properties or regarding the extremal eigenvalue of the
discretized and linearized operator, which affect the maximal admissible time
step size in explicit time integrators.

In this tutorial program, we implement two variants of fluxes that can be
controlled via a switch in the program (of course, it would be easy to make
them a run time parameter controlled via an input file). The first flux is
the local Lax--Friedrichs flux
@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
   \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
   \mathbf{n^-}.
@f]

In the original definition of the Lax--Friedrichs flux, a factor $\lambda =
\max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ is used
(corresponding to the maximal speed at which information is moving on
the two sides of the interface), stating
that the difference between the two states, $[\![\mathbf{w}]\!]$ is penalized
by the largest eigenvalue in the Euler flux, which is $\|\mathbf{u}\|+c$,
where $c=\sqrt{\gamma p / \rho}$ is the speed of sound. In the implementation
below, we modify the penalty term somewhat, given that the penalty is of
approximate nature anyway. We use
@f{align*}{
\lambda
&=
\frac{1}{2}\max\left(\sqrt{\|\mathbf{u^-}\|^2+(c^-)^2},
                     \sqrt{\|\mathbf{u}^+\|^2+(c^+)^2}\right)
\\
&=
\frac{1}{2}\sqrt{\max\left(\|\mathbf{u^-}\|^2+(c^-)^2,
                           \|\mathbf{u}^+\|^2+(c^+)^2\right)}.
@f}
The additional factor $\frac 12$ reduces the penalty strength (which results
in a reduced negative real part of the eigenvalues, and thus increases the
admissible time step size). Using the squares within the sums allows us to
reduce the number of expensive square root operations, which is 4 for the
original Lax--Friedrichs definition, to a single one.
This simplification leads to at most a factor of
2 in the reduction of the parameter $\lambda$, since $\|\mathbf{u}\|^2+c^2 \leq
\|\mathbf{u}\|^2+2 c |\mathbf{u}\| + c^2 = \left(\|\mathbf{u}\|+c\right)^2
\leq 2 \left(\|\mathbf{u}\|^2+c^2\right)$, with the last inequality following
from Young's inequality.

The second numerical flux is one proposed by Harten, Lax and van Leer, called
the HLL flux. It takes the different directions of propagation of the Euler
equations into account, depending on the speed of sound. It utilizes some
intermediate states $\bar{\mathbf{u}}$ and $\bar{c}$ to define the two
branches $s^\mathrm{p} = \max\left(0, \bar{\mathbf{u}}\cdot \mathbf{n} +
\bar{c}\right)$ and $s^\mathrm{n} = \min\left(0, \bar{\mathbf{u}}\cdot
\mathbf{n} - \bar{c}\right)$. From these branches, one then defines the flux
@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{s^\mathrm{p} \mathbf{F}(\mathbf{w}^-)-s^\mathrm{n} \mathbf{F}(\mathbf{w}^+)}
                   {s^\mathrm p - s^\mathrm{n} } +
\frac{s^\mathrm{p} s^\mathrm{n}}{s^\mathrm{p}-s^\mathrm{n}}
\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes \mathbf{n^-}.
@f]
Regarding the definition of the intermediate state $\bar{\mathbf{u}}$ and
$\bar{c}$, several variants have been proposed. The variant originally
proposed uses a density-averaged definition of the velocity, $\bar{\mathbf{u}}
= \frac{\sqrt{\rho^-} \mathbf{u}^- + \sqrt{\rho^+}\mathbf{u}^+}{\sqrt{\rho^-}
+ \sqrt{\rho^+}}$. Since we consider the Euler equations without shocks, we
simply use arithmetic means, $\bar{\mathbf{u}} = \frac{\mathbf{u}^- +
\mathbf{u}^+}{2}$ and $\bar{c} = \frac{c^- + c^+}{2}$, with $c^{\pm} =
\sqrt{\gamma p^{\pm} / \rho^{\pm}}$, in this tutorial program, and leave other
variants to a possible extension. We also note that the HLL flux has been
extended in the literature to the so-called HLLC flux, where C stands for the
ability to represent contact discontinuities.

At the boundaries with no neighboring state $\mathbf{w}^+$ available, it is
common practice to deduce suitable exterior values from the boundary
conditions (see the general literature on DG methods for details). In this
tutorial program, we consider three types of boundary conditions, namely
<b>inflow boundary conditions</b> where all components are prescribed,
@f[
\mathbf{w}^+ = \begin{pmatrix} \rho_\mathrm{D}(t)\\
(\rho \mathbf u)_{\mathrm D}(t) \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(Dirichlet)},
@f]
<b>subsonic outflow boundaries</b>, where we do not prescribe exterior
solutions as the flow field is leaving the domain and use the interior values
instead; we still need to prescribe the energy as there is one incoming
characteristic left in the Euler flux,
@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(mixed Neumann/Dirichlet)},
@f]
and <b>wall boundary condition</b> which describe a no-penetration
configuration:
@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- - 2 [(\rho \mathbf u)^-\cdot \mathbf n] \mathbf{n}
 \\ E^-\end{pmatrix}.
@f]

The polynomial expansion of the solution is finally inserted to the weak form
and test functions are replaced by the basis functions. This gives a discrete
in space, continuous in time nonlinear system with a finite number of unknown
coefficient values $w_j$, $j=1,\ldots,n_\text{dofs}$. Regarding the choice of
the polynomial degree in the DG method, there is no consensus in literature as
of 2019 as to what polynomial degrees are most efficient and the decision is
problem-dependent. Higher order polynomials ensure better convergence rates
and are thus superior for moderate to high accuracy requirements for
<b>smooth</b> solutions. At the same time, the volume-to-surface ratio
of where degrees of freedom are located,
increases with higher degrees, and this makes the effect of the numerical flux
weaker, typically reducing dissipation. However, in most of the cases the
solution is not smooth, at least not compared to the resolution that can be
afforded. This is true for example in incompressible fluid dynamics,
compressible fluid dynamics, and the related topic of wave propagation. In this
pre-asymptotic regime, the error is approximately proportional to the
numerical resolution, and other factors such as dispersion errors or the
dissipative behavior become more important. Very high order methods are often
ruled out because they come with more restrictive CFL conditions measured
against the number of unknowns, and they are also not as flexible when it
comes to representing complex geometries. Therefore, polynomial degrees
between two and six are most popular in practice, see e.g. the efficiency
evaluation in @cite FehnWallKronbichler2019 and references cited therein.

<a name="Explicittimeintegration"></a><h3>Explicit time integration</h3>


To discretize in time, we slightly rearrange the weak form and sum over all
cells:
@f[
\sum_{K \in \mathcal T_h} \left(\boldsymbol{\varphi}_i,
\frac{\partial \mathbf{w}}{\partial t}\right)_{K}
=
\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w})\right)_{K}
-\left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})\right>_{\partial K} +
\left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w)\right)_{K}
\right],
@f]
where $\boldsymbol{\varphi}_i$ runs through all basis functions with from 1 to
$n_\text{dofs}$.

We now denote by $\mathcal M$ the mass matrix with entries $\mathcal M_{ij} =
\sum_{K} \left(\boldsymbol{\varphi}_i,
\boldsymbol{\varphi}_j\right)_K$, and by
@f[
\mathcal L_h(t,\mathbf{w}_h) = \left[\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w}_h)\right)_{K}
- \left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w}_h)\right>_{\partial K}
+ \left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w_h)\right)_{K}
\right]\right]_{i=1,\ldots,n_\text{dofs}}.
@f]
the operator evaluating the right-hand side of the Euler operator, given a
function $\mathbf{w}_h$ associated with a global vector of unknowns
and the finite element in use. This function $\mathcal L_h$ is explicitly time-dependent as the
numerical flux evaluated at the boundary will involve time-dependent data
$\rho_\mathrm{D}$, $(\rho \mathbf{u})_\mathrm{D}$, and $E_\mathbf{D}$ on some
parts of the boundary, depending on the assignment of boundary
conditions. With this notation, we can write the discrete in space, continuous
in time system compactly as
@f[
\mathcal M \frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal L_h(t, \mathbf{w}_h),
@f]
where we have taken the liberty to also denote the global solution
vector by $\mathbf{w}_h$ (in addition to the the corresponding finite
element function). Equivalently, the system above has the form
@f[
\frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal M^{-1} \mathcal L_h(t, \mathbf{w}_h).
@f]

For hyperbolic systems discretized by high-order discontinuous Galerkin
methods, explicit time integration of this system is very popular. This is due
to the fact that the mass matrix $\mathcal M$ is block-diagonal (with each
block corresponding to only variables of the same kind defined on the same
cell) and thus easily inverted. In each time step -- or stage of a
Runge--Kutta scheme -- one only needs to evaluate the differential operator
once using the given data and subsequently apply the inverse of the mass
matrix. For implicit time stepping, on the other hand, one would first have to
linearize the equations and then iteratively solve the linear system, which
involves several residual evaluations and at least a dozen applications of
the linearized operator, as has been demonstrated in the step-33 tutorial
program.

Of course, the simplicity of explicit time stepping comes with a price, namely
conditional stability due to the so-called Courant--Friedrichs--Lewy (CFL)
condition. It states that the time step cannot be larger than the fastest
propagation of information by the discretized differential operator. In more
modern terms, the speed of propagation corresponds to the largest eigenvalue
in the discretized operator, and in turn depends on the mesh size, the
polynomial degree $p$ and the physics of the Euler operator, i.e., the
eigenvalues of the linearization of $\mathbf F(\mathbf w)$ with respect to
$\mathbf{w}$. In this program, we set the time step as follows:
@f[
\Delta t = \frac{\mathrm{Cr}}{p^{1.5}}\left(\frac{1}
           {\max\left[\frac{\|\mathbf{u}\|}{h_u} + \frac{c}{h_c}\right]}\right),
@f]

with the maximum taken over all quadrature points and all cells. The
dimensionless number $\mathrm{Cr}$ denotes the Courant number and can be
chosen up to a maximally stable number $\mathrm{Cr}_\text{max}$, whose value
depends on the selected time stepping method and its stability properties. The
power $p^{1.5}$ used for the polynomial scaling is heuristic and represents
the closest fit for polynomial degrees between 1 and 8, see e.g.
@cite SchoederKormann2018. In the limit of higher degrees, $p>10$, a scaling of
$p^2$ is more accurate, related to the inverse estimates typically used for
interior penalty methods. Regarding the <i>effective</i> mesh sizes $h_u$ and
$h_c$ used in the formula, we note that the convective transport is
directional. Thus an appropriate scaling is to use the element length in the
direction of the velocity $\mathbf u$. The code below derives this scaling
from the inverse of the Jacobian from the reference to real cell, i.e., we
approximate $\frac{\|\mathbf{u}\|}{h_u} \approx \|J^{-1} \mathbf
u\|_{\infty}$. The acoustic waves, instead, are isotropic in character, which
is why we use the smallest feature size, represented by the smallest singular
value of $J$, for the acoustic scaling $h_c$. Finally, we need to add the
convective and acoustic limits, as the Euler equations can transport
information with speed $\|\mathbf{u}\|+c$.

In this tutorial program, we use a specific variant of <a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">explicit
Runge--Kutta methods</a>, which in general use the following update procedure
from the state $\mathbf{w}_h^{n}$ at time $t^n$ to the new time $t^{n+1}$ with
$\Delta t = t^{n+1}-t^n$:
@f[
\begin{aligned}
\mathbf{k}_1 &= \mathcal M^{-1} \mathcal L_h\left(t^n, \mathbf{w}_h^n\right),
\\
\mathbf{k}_2 &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_2\Delta t,
                       \mathbf{w}_h^n + a_{21} \Delta t \mathbf{k}_1\right),
\\
&\vdots \\
\mathbf{k}_s &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_s\Delta t,
  \mathbf{w}_h^n + \sum_{j=1}^{s-1} a_{sj} \Delta t \mathbf{k}_j\right),
\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^n + \Delta t\left(b_1 \mathbf{k}_1 +
b_2 \mathbf{k}_2 + \ldots + b_s \mathbf{k}_s\right).
\end{aligned}
@f]
The vectors $\mathbf{k}_i$, $i=1,\ldots,s$, in an $s$-stage scheme are
evaluations of the operator at some intermediate state and used to define the
end-of-step value $\mathbf{w}_h^{n+1}$ via some linear combination. The scalar
coefficients in this scheme, $c_i$, $a_{ij}$, and $b_j$, are defined such that
certain conditions are satisfied for higher order schemes, the most basic one
being $c_i = \sum_{j=1}^{i-1}a_{ij}$. The parameters are typically collected in
the form of a so-called <a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge%E2%80%93Kutta_methods">Butcher
tableau</a> that collects all of the coefficients that define the
scheme. For a five-stage scheme, it would look like this:
@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_{21} \\
c_3 & a_{31} & a_{32} \\
c_4 & a_{41} & a_{42} & a_{43} \\
c_5 & a_{51} & a_{52} & a_{53} & a_{54} \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}
@f]

In this tutorial program, we use a subset of explicit Runge--Kutta methods,
so-called low-storage Runge--Kutta methods (LSRK), which assume additional
structure in the coefficients. In the variant used by reference
@cite KennedyCarpenterLewis2000, the assumption is to use Butcher tableaus of
the form
@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_1 \\
c_3 & b_1 & a_2 \\
c_4 & b_1 & b_2 & a_3 \\
c_5 & b_1 & b_2 & b_3 & a_4 \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}
@f]
With such a definition, the update to $\mathbf{w}_h^n$ shares the storage with
the information for the intermediate values $\mathbf{k}_i$. Starting with
$\mathbf{w}^{n+1}=\mathbf{w}^n$ and $\mathbf{r}_1 = \mathbf{w}^n$, the update
in each of the $s$ stages simplifies to
@f[
\begin{aligned}
\mathbf{k}_i &=
\mathcal M^{-1} \mathcal L_h\left(t^n+c_i\Delta t, \mathbf{r}_{i} \right),\\
\mathbf{r}_{i+1} &= \mathbf{w}_h^{n+1} + \Delta t \, a_i \mathbf{k}_i,\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^{n+1} + \Delta t \, b_i \mathbf{k}_i.
\end{aligned}
@f]
Besides the vector $\mathbf w_h^{n+1}$ that is successively updated, this scheme
only needs two auxiliary vectors, namely the vector $\mathbf{k}_i$ to hold the
evaluation of the differential operator, and the vector $\mathbf{r}_i$ that
holds the right-hand side for the differential operator application. In
subsequent stages $i$, the values $\mathbf{k}_i$ and $\mathbf{r}_i$ can use
the same storage.

The main advantages of low-storage variants are the reduced memory consumption
on the one hand (if a very large number of unknowns must be fit in memory,
holding all $\mathbf{k}_i$ to compute subsequent updates can be a limit
already for $s$ in between five and eight -- recall that we are using
an explicit scheme, so we do not need to store any matrices that are
typically much larger than a few vectors), and the reduced memory access on
the other. In this program, we are particularly interested in the latter
aspect. Since cost of operator evaluation is only a small multiple of the cost
of simply streaming the input and output vector from memory with the optimized
matrix-free methods of deal.II, we must consider the cost of vector updates,
and low-storage variants can deliver up to twice the throughput of
conventional explicit Runge--Kutta methods for this reason, see e.g. the
analysis in @cite SchoederKormann2018.

Besides three variants for third, fourth and fifth order accuracy from the
reference @cite KennedyCarpenterLewis2000, we also use a fourth-order accurate
variant with seven stages that was optimized for acoustics setups from
@cite TseliosSimos2007. Acoustic problems are one of the interesting aspects of
the subsonic regime of the Euler equations where compressibility leads to the
transmission of sound waves; often, one uses further simplifications of the
linearized Euler equations around a background state or the acoustic wave
equation around a fixed frame.


<a name="Fastevaluationofintegralsbymatrixfreetechniques"></a><h3>Fast evaluation of integrals by matrix-free techniques</h3>


The major ingredients used in this program are the fast matrix-free techniques
we use to evaluate the operator $\mathcal L_h$ and the inverse mass matrix
$\mathcal M$. Actually, the term <i>matrix-free</i> is a slight misnomer,
because we are working with a nonlinear operator and do not linearize the
operator that in turn could be represented by a matrix. However, fast
evaluation of integrals has become popular as a replacement of sparse
matrix-vector products, as shown in step-37 and step-59, and we have coined
this infrastructure <i>matrix-free functionality</i> in deal.II for this
reason. Furthermore, the inverse mass matrix is indeed applied in a
matrix-free way, as detailed below.

The matrix-free infrastructure allows us to quickly evaluate the integrals in
the weak forms. The ingredients are the fast interpolation from solution
coefficients into values and derivatives at quadrature points, point-wise
operations at quadrature points (where we implement the differential operator
as derived above), as well as multiplication by all test functions and
summation over quadrature points. The first and third component make use of
sum factorization and have been extensively discussed in the step-37 tutorial
program for the cell integrals and step-59 for the face integrals. The only
difference is that we now deal with a system of $d+2$ components, rather than
the scalar systems in previous tutorial programs. In the code, all that
changes is a template argument of the FEEvaluation and FEFaceEvaluation
classes, the one to set the number of components. The access to the vector is
the same as before, all handled transparently by the evaluator. We also note
that the variant with a single evaluator chosen in the code below is not the
only choice -- we could also have used separate evalators for the separate
components $\rho$, $\rho \mathbf u$, and $E$; given that we treat all
components similarly (also reflected in the way we state the equation as a
vector system), this would be more complicated here. As before, the
FEEvaluation class provides explicit vectorization by combining the operations
on several cells (and faces), involving data types called
VectorizedArray. Since the arithmetic operations are overloaded for this type,
we do not have to bother with it all that much, except for the evaluation of
functions through the Function interface, where we need to provide particular
<i>vectorized</i> evaluations for several quadrature point locations at once.

A more substantial change in this program is the operation at quadrature
points: Here, the multi-component evaluators provide us with return types not
discussed before. Whereas FEEvaluation::get_value() would return a scalar
(more precisely, a VectorizedArray type due to vectorization across cells) for
the Laplacian of step-37, it now returns a type that is
`Tensor<1,dim+2,VectorizedArray<Number>>`. Likewise, the gradient type is now
`Tensor<1,dim+2,Tensor<1,dim,VectorizedArray<Number>>>`, where the outer
tensor collects the `dim+2` components of the Euler system, and the inner
tensor the partial derivatives in the various directions. For example, the
flux $\mathbf{F}(\mathbf{w})$ of the Euler system is of this type. In order to reduce the amount of
code we have to write for spelling out these types, we use the C++ `auto`
keyword where possible.

From an implementation point of view, the nonlinearity is not a big
difficulty: It is introduced naturally as we express the terms of the Euler
weak form, for example in the form of the momentum term $\rho \mathbf{u}
\otimes \mathbf{u}$. To obtain this expression, we first deduce the velocity
$\mathbf{u}$ from the momentum variable $\rho \mathbf{u}$. Given that $\rho
\mathbf{u}$ is represented as a $p$-degree polynomial, as is $\rho$, the
velocity $\mathbf{u}$ is a rational expression in terms of the reference
coordinates $\hat{\mathbf{x}}$. As we perform the multiplication $(\rho
\mathbf{u})\otimes \mathbf{u}$, we obtain an expression that is the
ratio of two polynomials, with polynomial degree $2p$ in the
numerator and polynomial degree $p$ in the denominator. Combined with the
gradient of the test function, the integrand is of degree $3p$ in the
numerator and $p$ in the denominator already for affine cells, i.e.,
for parallelograms/ parallelepipeds.
For curved cells, additional polynomial and rational expressions
appear when multiplying the integrand by the determinant of the Jacobian of
the mapping. At this point, one usually needs to give up on insisting on exact
integration, and take whatever accuracy the Gaussian (more precisely,
Gauss--Legrende) quadrature provides. The situation is then similar to the one
for the Laplace equation, where the integrand contains rational expressions on
non-affince cells and is also only integrated approximately. As these formulas
only integrate polynomials exactly, we have to live with the <a
href="https://mathoverflow.net/questions/26018/what-are-variational-crimes-and-who-coined-the-term">variational
crime</a> in the form of an integration error.

While inaccurate integration is usually tolerable for elliptic problems, for
hyperbolic problems inexact integration causes some headache in the form of an
effect called <b>aliasing</b>. The term comes from signal processing and
expresses the situation of inappropriate, too coarse sampling. In terms of
quadrature, the inappropriate sampling means that we use too few quadrature
points compared to what would be required to accurately sample the
variable-coefficient integrand. It has been shown in the DG literature that
aliasing errors can introduce unphysical oscillations in the numerical
solution for <i>barely</i> resolved simulations. The fact that aliasing mostly
affects coarse resolutions -- whereas finer meshes with the same scheme
work fine -- is not surprising because well-resolved simulations
tend to be smooth on length-scales of a cell (i.e., they have
small coefficients in the higher polynomial degrees that are missed by
too few quadrature points, whereas the main solution contribution in the lower
polynomial degrees is still well-captured -- this is simply a consequence of Taylor's
theorem). To address this topic, various approaches have been proposed in the
DG literature. One technique is filtering which damps the solution components
pertaining to higher polynomial degrees. As the chosen nodal basis is not
hierarchical, this would mean to transform from the nodal basis into a
hierarchical one (e.g., a modal one based on Legendre polynomials) where the
contributions within a cell are split by polynomial degrees. In that basis,
one could then multiply the solution coefficients associated with higher
degrees by a small number, keep the lower ones intact (to not destroy consistency), and
then transform back to the nodal basis. However, filters reduce the accuracy of the
method. Another, in some sense simpler, strategy is to use more quadrature
points to capture non-linear terms more accurately. Using more than $p+1$
quadrature points per coordinate directions is sometimes called
over-integration or consistent integration. The latter name is most common in
the context of the incompressible Navier-Stokes equations, where the
$\mathbf{u}\otimes \mathbf{u}$ nonlinearity results in polynomial integrands
of degree $3p$ (when also considering the test function), which can be
integrated exactly with $\textrm{floor}\left(\frac{3p}{2}\right)+1$ quadrature
points per direction as long as the element geometry is affine. In the context
of the Euler equations with non-polynomial integrands, the choice is less
clear. Depending on the variation in the various variables both
$\textrm{floor}\left(\frac{3p}{2}\right)+1$ or $2p+1$ points (integrating
exactly polynomials of degree $3p$ or $4p$, respectively) are common.

To reflect this variability in the choice of quadrature in the program, we
keep the number of quadrature points a variable to be specified just as the
polynomial degree, and note that one would make different choices depending
also on the flow configuration. The default choice is $p+2$ points -- a bit
more than the minimum possible of $p+1$ points. The FEEvaluation and
FEFaceEvaluation classes allow to seamlessly change the number of points by a
template parameter, such that the program does not get more complicated
because of that.


<a name="Evaluationoftheinversemassmatrixwithmatrixfreetechniques"></a><h3>Evaluation of the inverse mass matrix with matrix-free techniques</h3>


The last ingredient is the evaluation of the inverse mass matrix $\mathcal
M^{-1}$. In DG methods with explicit time integration, mass matrices are
block-diagonal and thus easily inverted -- one only needs to invert the
diagonal blocks. However, given the fact that matrix-free evaluation of
integrals is closer in cost to the access of the vectors only, even the
application of a block-diagonal matrix (e.g. via an array of LU factors) would
be several times more expensive than evaluation of $\mathcal L_h$
simply because just storing and loading matrices of size
`dofs_per_cell` times `dofs_per_cell` for higher order finite elements
repeatedly is expensive. As this is
clearly undesirable, part of the community has moved to bases where the mass
matrix is diagonal, for example the <i>L<sub>2</sub></i>-orthogonal Legendre basis using
hierarchical polynomials or Lagrange polynomials on the points of the Gaussian
quadrature (which is just another way of utilizing Legendre
information). While the diagonal property breaks down for deformed elements,
the error made by taking a diagonal mass matrix and ignoring the rest (a
variant of mass lumping, though not the one with an additional integration
error as utilized in step-48) has been shown to not alter discretization
accuracy. The Lagrange basis in the points of Gaussian quadrature is sometimes
also referred to as a collocation setup, as the nodal points of the
polynomials coincide (= are "co-located") with the points of quadrature, obviating some
interpolation operations. Given the fact that we want to use more quadrature
points for nonlinear terms in $\mathcal L_h$, however, the collocation
property is lost. (More precisely, it is still used in FEEvaluation and
FEFaceEvaluation after a change of basis, see the matrix-free paper
@cite KronbichlerKormann2019.)

In this tutorial program, we use the collocation idea for the application of
the inverse mass matrix, but with a slight twist. Rather than using the
collocation via Lagrange polynomials in the points of Gaussian quadrature, we
prefer a conventional Lagrange basis in Gauss-Lobatto points as those make the
evaluation of face integrals cheap. This is because for Gauss-Lobatto
points, some of the node points are located on the faces of the cell
and it is not difficult to show that on any given face, the only shape
functions with non-zero values are exactly the ones whose node points
are in fact located on that face. One could of course also use the
Gauss-Lobatto quadrature (with some additional integration error) as was done
in step-48, but we do not want to sacrifice accuracy as these
quadrature formulas are generally of lower order than the general
Gauss quadrature formulas. Instead, we use an idea described in the reference
@cite KronbichlerSchoeder2016 where it was proposed to change the basis for the
sake of applying the inverse mass matrix. Let us denote by $S$ the matrix of
shape functions evaluated at quadrature points, with shape functions in the row
of the matrix and quadrature points in columns. Then, the mass matrix on a cell
$K$ is given by
@f[
\mathcal M^K = S J^K S^\mathrm T.
@f]
Here, $J^K$ is the diagonal matrix with the determinant of the Jacobian times
the quadrature weight (JxW) as entries. The matrix $S$ is constructed as the
Kronecker product (tensor product) of one-dimensional matrices, e.g. in 3D as
@f[
S = S_{\text{1D}}\otimes S_{\text{1D}}\otimes S_{\text{1D}},
@f]
which is the result of the basis functions being a tensor product of
one-dimensional shape functions and the quadrature formula being the tensor
product of 1D quadrature formulas. For the case that the number of polynomials
equals the number of quadrature points, all matrices in $S J^K S^\mathrm T$
are square, and also the ingredients to $S$ in the Kronecker product are
square. Thus, one can invert each matrix to form the overall inverse,
@f[
\left(\mathcal M^K\right)^{-1} = S_{\text{1D}}^{-\mathrm T}\otimes
S_{\text{1D}}^{-\mathrm T}\otimes S_{\text{1D}}^{-\mathrm T}
\left(J^K\right)^{-1}
S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}.
@f]
This formula is of exactly the same structure as the steps in the forward
evaluation of integrals with sum factorization techniques (i.e., the
FEEvaluation and MatrixFree framework of deal.II). Hence, we can utilize the
same code paths with a different interpolation matrix,
$S_{\mathrm{1D}}^{-\mathrm{T}}$ rather than $S_{\mathrm{1D}}$.

The class MatrixFreeOperators::CellwiseInverseMassMatrix implements this
operation: It changes from the basis contained in the finite element (in this
case, FE_DGQ) to the Lagrange basis in Gaussian quadrature points. Here, the
inverse of a diagonal mass matrix can be evaluated, which is simply the inverse
of the `JxW` factors (i.e., the quadrature weight times the determinant of the
Jacobian from reference to real coordinates). Once this is done, we can change
back to the standard nodal Gauss-Lobatto basis.

The advantage of this particular way of applying the inverse mass matrix is
a cost similar to the forward application of a mass matrix, which is cheaper
than the evaluation of the spatial operator $\mathcal L_h$
with over-integration and face integrals. (We
will demonstrate this with detailed timing information in the
<a href="#Results">results section</a>.) In fact, it
is so cheap that it is limited by the bandwidth of reading the source vector,
reading the diagonal, and writing into the destination vector on most modern
architectures. The hardware used for the result section allows to do the
computations at least twice as fast as the streaming of the vectors from
memory.


<a name="Thetestcase"></a><h3>The test case</h3>


In this tutorial program, we implement two test cases. The first case is a
convergence test limited to two space dimensions. It runs a so-called
isentropic vortex which is transported via a background flow field. The second
case uses a more exciting setup: We start with a cylinder immersed in a
channel, using the GridGenerator::channel_with_cylinder() function. Here, we
impose a subsonic initial field at Mach number of $\mathrm{Ma}=0.307$ with a
constant velocity in $x$ direction. At the top and bottom walls as well as at
the cylinder, we impose a no-penetration (i.e., tangential flow)
condition. This setup forces the flow to re-orient as compared to the initial
condition, which results in a big sound wave propagating away from the
cylinder. In upstream direction, the wave travels more slowly (as it
has to move against the oncoming gas), including a
discontinuity in density and pressure. In downstream direction, the transport
is faster as sound propagation and fluid flow go in the same direction, which smears
out the discontinuity somewhat. Once the sound wave hits the upper and lower
walls, the sound is reflected back, creating some nice shapes as illustrated
in the <a href="#Results">results section</a> below.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The include files are similar to the previous matrix-free tutorial programs
 * step-37, step-48, and step-59
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/time_stepping.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/vectorization.h>
 * 
 * #include <deal.II/distributed/tria.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * 
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/tria.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * 
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * #include <deal.II/matrix_free/matrix_free.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iomanip>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The following file includes the CellwiseInverseMassMatrix data structure
 * that we will use for the mass matrix inversion, the only new include
 * file for this tutorial program:
 * 
 * @code
 * #include <deal.II/matrix_free/operators.h>
 * 
 * 
 * 
 * namespace Euler_DG
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Similarly to the other matrix-free tutorial programs, we collect all
 * parameters that control the execution of the program at the top of the
 * file. Besides the dimension and polynomial degree we want to run with, we
 * also specify a number of points in the Gaussian quadrature formula we
 * want to use for the nonlinear terms in the Euler equations. Furthermore,
 * we specify the time interval for the time-dependent problem, and
 * implement two different test cases. The first one is an analytical
 * solution in 2D, whereas the second is a channel flow around a cylinder as
 * described in the introduction. Depending on the test case, we also change
 * the final time up to which we run the simulation, and a variable
 * `output_tick` that specifies in which intervals we want to write output
 * (assuming that the tick is larger than the time step size).
 * 
 * @code
 *   constexpr unsigned int testcase             = 0;
 *   constexpr unsigned int dimension            = 2;
 *   constexpr unsigned int n_global_refinements = 3;
 *   constexpr unsigned int fe_degree            = 5;
 *   constexpr unsigned int n_q_points_1d        = fe_degree + 2;
 * 
 *   using Number = double;
 * 
 *   constexpr double gamma       = 1.4;
 *   constexpr double final_time  = testcase == 0 ? 10 : 2.0;
 *   constexpr double output_tick = testcase == 0 ? 1 : 0.05;
 * 
 * @endcode
 * 
 * Next off are some details of the time integrator, namely a Courant number
 * that scales the time step size in terms of the formula $\Delta t =
 * \text{Cr} n_\text{stages} \frac{h}{(p+1)^{1.5} (\|\mathbf{u} +
 * c)_\text{max}}$, as well as a selection of a few low-storage Runge--Kutta
 * methods. We specify the Courant number per stage of the Runge--Kutta
 * scheme, as this gives a more realistic expression of the numerical cost
 * for schemes of various numbers of stages.
 * 
 * @code
 *   const double courant_number = 0.15 / std::pow(fe_degree, 1.5);
 *   enum LowStorageRungeKuttaScheme
 *   {
 *     stage_3_order_3, /* Kennedy, Carpenter, Lewis, 2000 */
 *     stage_5_order_4, /* Kennedy, Carpenter, Lewis, 2000 */
 *     stage_7_order_4, /* Tselios, Simos, 2007 */
 *     stage_9_order_5, /* Kennedy, Carpenter, Lewis, 2000 */
 *   };
 *   constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4;
 * 
 * @endcode
 * 
 * Eventually, we select a detail of the spatial discretization, namely the
 * numerical flux (Riemann solver) at the faces between cells. For this
 * program, we have implemented a modified variant of the Lax--Friedrichs
 * flux and the Harten--Lax--van Leer (HLL) flux.
 * 
 * @code
 *   enum EulerNumericalFlux
 *   {
 *     lax_friedrichs_modified,
 *     harten_lax_vanleer,
 *   };
 *   constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified;
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * We now define a class with the exact solution for the test case 0 and one
 * with a background flow field for test case 1 of the channel. Given that
 * the Euler equations are a problem with $d+2$ equations in $d$ dimensions,
 * we need to tell the Function base class about the correct number of
 * components.
 * 
 * @code
 *   template <int dim>
 *   class ExactSolution : public Function<dim>
 *   {
 *   public:
 *     ExactSolution(const double time)
 *       : Function<dim>(dim + 2, time)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * As far as the actual function implemented is concerned, the analytical
 * test case is an isentropic vortex case (see e.g. the book by Hesthaven
 * and Warburton, Example 6.1 in Section 6.6 on page 209) which fulfills the
 * Euler equations with zero force term on the right hand side. Given that
 * definition, we return either the density, the momentum, or the energy
 * depending on which component is requested. Note that the original
 * definition of the density involves the $\frac{1}{\gamma -1}$-th power of
 * some expression. Since `std::pow()` has pretty slow implementations on
 * some systems, we replace it by logarithm followed by exponentiation (of
 * base 2), which is mathematically equivalent but usually much better
 * optimized. This formula might lose accuracy in the last digits
 * for very small numbers compared to `std::pow()`, but we are happy with
 * it anyway, since small numbers map to data close to 1.
 *   

 * 
 * For the channel test case, we simply select a density of 1, a velocity of
 * 0.4 in $x$ direction and zero in the other directions, and an energy that
 * corresponds to a speed of sound of 1.3 measured against the background
 * velocity field, computed from the relation $E = \frac{c^2}{\gamma (\gamma
 * -1)} + \frac 12 \rho \|u\|^2$.
 * 
 * @code
 *   template <int dim>
 *   double ExactSolution<dim>::value(const Point<dim> & x,
 *                                    const unsigned int component) const
 *   {
 *     const double t = this->get_time();
 * 
 *     switch (testcase)
 *       {
 *         case 0:
 *           {
 *             Assert(dim == 2, ExcNotImplemented());
 *             const double beta = 5;
 * 
 *             Point<dim> x0;
 *             x0[0] = 5.;
 *             const double radius_sqr =
 *               (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t;
 *             const double factor =
 *               beta / (numbers::PI * 2) * std::exp(1. - radius_sqr);
 *             const double density_log = std::log2(
 *               std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor));
 *             const double density = std::exp2(density_log * (1. / (gamma - 1.)));
 *             const double u       = 1. - factor * (x[1] - x0[1]);
 *             const double v       = factor * (x[0] - t - x0[0]);
 * 
 *             if (component == 0)
 *               return density;
 *             else if (component == 1)
 *               return density * u;
 *             else if (component == 2)
 *               return density * v;
 *             else
 *               {
 *                 const double pressure =
 *                   std::exp2(density_log * (gamma / (gamma - 1.)));
 *                 return pressure / (gamma - 1.) +
 *                        0.5 * (density * u * u + density * v * v);
 *               }
 *           }
 * 
 *         case 1:
 *           {
 *             if (component == 0)
 *               return 1.;
 *             else if (component == 1)
 *               return 0.4;
 *             else if (component == dim + 1)
 *               return 3.097857142857143;
 *             else
 *               return 0.;
 *           }
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *           return 0.;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LowstorageexplicitRungeKuttatimeintegrators"></a> 
 * <h3>Low-storage explicit Runge--Kutta time integrators</h3>
 * 

 * 
 * The next few lines implement a few low-storage variants of Runge--Kutta
 * methods. These methods have specific Butcher tableaux with coefficients
 * $b_i$ and $a_i$ as shown in the introduction. As usual in Runge--Kutta
 * method, we can deduce time steps, $c_i = \sum_{j=1}^{i-2} b_i + a_{i-1}$
 * from those coefficients. The main advantage of this kind of scheme is the
 * fact that only two vectors are needed per stage, namely the accumulated
 * part of the solution $\mathbf{w}$ (that will hold the solution
 * $\mathbf{w}^{n+1}$ at the new time $t^{n+1}$ after the last stage), the
 * update vector $\mathbf{r}_i$ that gets evaluated during the stages, plus
 * one vector $\mathbf{k}_i$ to hold the evaluation of the operator. Such a
 * Runge--Kutta setup reduces the memory storage and memory access. As the
 * memory bandwidth is often the performance-limiting factor on modern
 * hardware when the evaluation of the differential operator is
 * well-optimized, performance can be improved over standard time
 * integrators. This is true also when taking into account that a
 * conventional Runge--Kutta scheme might allow for slightly larger time
 * steps as more free parameters allow for better stability properties.
 *   

 * 
 * In this tutorial programs, we concentrate on a few variants of
 * low-storage schemes defined in the article by Kennedy, Carpenter, and
 * Lewis (2000), as well as one variant described by Tselios and Simos
 * (2007). There is a large series of other schemes available, which could
 * be addressed by additional sets of coefficients or slightly different
 * update formulas.
 *   

 * 
 * We define a single class for the four integrators, distinguished by the
 * enum described above. To each scheme, we then fill the vectors for the
 * $b_i$ and $a_i$ to the given variables in the class.
 * 
 * @code
 *   class LowStorageRungeKuttaIntegrator
 *   {
 *   public:
 *     LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme)
 *     {
 *       TimeStepping::runge_kutta_method lsrk;
 * @endcode
 * 
 * First comes the three-stage scheme of order three by Kennedy et al.
 * (2000). While its stability region is significantly smaller than for
 * the other schemes, it only involves three stages, so it is very
 * competitive in terms of the work per stage.
 * 
 * @code
 *       switch (scheme)
 *         {
 *           case stage_3_order_3:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The next scheme is a five-stage scheme of order four, again
 * defined in the paper by Kennedy et al. (2000).
 * 
 * @code
 *           case stage_5_order_4:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The following scheme of seven stages and order four has been
 * explicitly derived for acoustics problems. It is a balance of
 * accuracy for imaginary eigenvalues among fourth order schemes,
 * combined with a large stability region. Since DG schemes are
 * dissipative among the highest frequencies, this does not
 * necessarily translate to the highest possible time step per
 * stage. In the context of the present tutorial program, the
 * numerical flux plays a crucial role in the dissipation and thus
 * also the maximal stable time step size. For the modified
 * Lax--Friedrichs flux, this scheme is similar to the
 * `stage_5_order_4` scheme in terms of step size per stage if only
 * stability is considered, but somewhat less efficient for the HLL
 * flux.
 * 
 * @code
 *           case stage_7_order_4:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The last scheme included here is the nine-stage scheme of order
 * five from Kennedy et al. (2000). It is the most accurate among
 * the schemes used here, but the higher order of accuracy
 * sacrifices some stability, so the step length normalized per
 * stage is less than for the fourth order schemes.
 * 
 * @code
 *           case stage_9_order_5:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5;
 *               break;
 *             }
 * 
 *           default:
 *             AssertThrow(false, ExcNotImplemented());
 *         }
 *       TimeStepping::LowStorageRungeKutta<
 *         LinearAlgebra::distributed::Vector<Number>>
 *         rk_integrator(lsrk);
 *       rk_integrator.get_coefficients(ai, bi, ci);
 *     }
 * 
 *     unsigned int n_stages() const
 *     {
 *       return bi.size();
 *     }
 * 
 * @endcode
 * 
 * The main function of the time integrator is to go through the stages,
 * evaluate the operator, prepare the $\mathbf{r}_i$ vector for the next
 * evaluation, and update the solution vector $\mathbf{w}$. We hand off
 * the work to the `pde_operator` involved in order to be able to merge
 * the vector operations of the Runge--Kutta setup with the evaluation of
 * the differential operator for better performance, so all we do here is
 * to delegate the vectors and coefficients.
 *     

 * 
 * We separately call the operator for the first stage because we need
 * slightly modified arguments there: We evaluate the solution from
 * the old solution $\mathbf{w}^n$ rather than a $\mathbf r_i$ vector, so
 * the first argument is `solution`. We here let the stage vector
 * $\mathbf{r}_i$ also hold the temporary result of the evaluation, as it
 * is not used otherwise. For all subsequent stages, we use the vector
 * `vec_ki` as the second vector argument to store the result of the
 * operator evaluation. Finally, when we are at the last stage, we must
 * skip the computation of the vector $\mathbf{r}_{s+1}$ as there is no
 * coefficient $a_s$ available (nor will it be used).
 * 
 * @code
 *     template <typename VectorType, typename Operator>
 *     void perform_time_step(const Operator &pde_operator,
 *                            const double    current_time,
 *                            const double    time_step,
 *                            VectorType &    solution,
 *                            VectorType &    vec_ri,
 *                            VectorType &    vec_ki) const
 *     {
 *       AssertDimension(ai.size() + 1, bi.size());
 * 
 *       pde_operator.perform_stage(current_time,
 *                                  bi[0] * time_step,
 *                                  ai[0] * time_step,
 *                                  solution,
 *                                  vec_ri,
 *                                  solution,
 *                                  vec_ri);
 * 
 *       for (unsigned int stage = 1; stage < bi.size(); ++stage)
 *         {
 *           const double c_i = ci[stage];
 *           pde_operator.perform_stage(current_time + c_i * time_step,
 *                                      bi[stage] * time_step,
 *                                      (stage == bi.size() - 1 ?
 *                                         0 :
 *                                         ai[stage] * time_step),
 *                                      vec_ri,
 *                                      vec_ki,
 *                                      solution,
 *                                      vec_ri);
 *         }
 *     }
 * 
 *   private:
 *     std::vector<double> bi;
 *     std::vector<double> ai;
 *     std::vector<double> ci;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofpointwiseoperationsoftheEulerequations"></a> 
 * <h3>Implementation of point-wise operations of the Euler equations</h3>
 * 

 * 
 * In the following functions, we implement the various problem-specific
 * operators pertaining to the Euler equations. Each function acts on the
 * vector of conserved variables $[\rho, \rho\mathbf{u}, E]$ that we hold in
 * the solution vectors, and computes various derived quantities.
 *   

 * 
 * First out is the computation of the velocity, that we derive from the
 * momentum variable $\rho \mathbf{u}$ by division by $\rho$. One thing to
 * note here is that we decorate all those functions with the keyword
 * `DEAL_II_ALWAYS_INLINE`. This is a special macro that maps to a
 * compiler-specific keyword that tells the compiler to never create a
 * function call for any of those functions, and instead move the
 * implementation <a
 * href="https://en.wikipedia.org/wiki/Inline_function">inline</a> to where
 * they are called. This is critical for performance because we call into some
 * of those functions millions or billions of times: For example, we both use
 * the velocity for the computation of the flux further down, but also for the
 * computation of the pressure, and both of these places are evaluated at
 * every quadrature point of every cell. Making sure these functions are
 * inlined ensures not only that the processor does not have to execute a jump
 * instruction into the function (and the corresponding return jump), but also
 * that the compiler can re-use intermediate information from one function's
 * context in code that comes after the place where the function was called.
 * (We note that compilers are generally quite good at figuring out which
 * functions to inline by themselves. Here is a place where compilers may or
 * may not have figured it out by themselves but where we know for sure that
 * inlining is a win.)
 *   

 * 
 * Another trick we apply is a separate variable for the inverse density
 * $\frac{1}{\rho}$. This enables the compiler to only perform a single
 * division for the flux, despite the division being used at several
 * places. As divisions are around ten to twenty times as expensive as
 * multiplications or additions, avoiding redundant divisions is crucial for
 * performance. We note that taking the inverse first and later multiplying
 * with it is not equivalent to a division in floating point arithmetic due
 * to roundoff effects, so the compiler is not allowed to exchange one way by
 * the other with standard optimization flags. However, it is also not
 * particularly difficult to write the code in the right way.
 *   

 * 
 * To summarize, the chosen strategy of always inlining and careful
 * definition of expensive arithmetic operations allows us to write compact
 * code without passing all intermediate results around, despite making sure
 * that the code maps to excellent machine code.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim, Number>
 *     euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Number inverse_density = Number(1.) / conserved_variables[0];
 * 
 *     Tensor<1, dim, Number> velocity;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       velocity[d] = conserved_variables[1 + d] * inverse_density;
 * 
 *     return velocity;
 *   }
 * 
 * @endcode
 * 
 * The next function computes the pressure from the vector of conserved
 * variables, using the formula $p = (\gamma - 1) \left(E - \frac 12 \rho
 * \mathbf{u}\cdot \mathbf{u}\right)$. As explained above, we use the
 * velocity from the `euler_velocity()` function. Note that we need to
 * specify the first template argument `dim` here because the compiler is
 * not able to deduce it from the arguments of the tensor, whereas the
 * second argument (number type) can be automatically deduced.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Number
 *     euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Tensor<1, dim, Number> velocity =
 *       euler_velocity<dim>(conserved_variables);
 * 
 *     Number rho_u_dot_u = conserved_variables[1] * velocity[0];
 *     for (unsigned int d = 1; d < dim; ++d)
 *       rho_u_dot_u += conserved_variables[1 + d] * velocity[d];
 * 
 *     return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u);
 *   }
 * 
 * @endcode
 * 
 * Here is the definition of the Euler flux function, i.e., the definition
 * of the actual equation. Given the velocity and pressure (that the
 * compiler optimization will make sure are done only once), this is
 * straight-forward given the equation stated in the introduction.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>>
 *     euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Tensor<1, dim, Number> velocity =
 *       euler_velocity<dim>(conserved_variables);
 *     const Number pressure = euler_pressure<dim>(conserved_variables);
 * 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         flux[0][d] = conserved_variables[1 + d];
 *         for (unsigned int e = 0; e < dim; ++e)
 *           flux[e + 1][d] = conserved_variables[e + 1] * velocity[d];
 *         flux[d + 1][d] += pressure;
 *         flux[dim + 1][d] =
 *           velocity[d] * (conserved_variables[dim + 1] + pressure);
 *       }
 * 
 *     return flux;
 *   }
 * 
 * @endcode
 * 
 * This next function is a helper to simplify the implementation of the
 * numerical flux, implementing the action of a tensor of tensors (with
 * non-standard outer dimension of size `dim + 2`, so the standard overloads
 * provided by deal.II's tensor classes do not apply here) with another
 * tensor of the same inner dimension, i.e., a matrix-vector product.
 * 
 * @code
 *   template <int n_components, int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, n_components, Number>
 *     operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
 *               const Tensor<1, dim, Number> &                         vector)
 *   {
 *     Tensor<1, n_components, Number> result;
 *     for (unsigned int d = 0; d < n_components; ++d)
 *       result[d] = matrix[d] * vector;
 *     return result;
 *   }
 * 
 * @endcode
 * 
 * This function implements the numerical flux (Riemann solver). It gets the
 * state from the two sides of an interface and the normal vector, oriented
 * from the side of the solution $\mathbf{w}^-$ towards the solution
 * $\mathbf{w}^+$. In finite volume methods which rely on piece-wise
 * constant data, the numerical flux is the central ingredient as it is the
 * only place where the physical information is entered. In DG methods, the
 * numerical flux is less central due to the polynomials within the elements
 * and the physical flux used there. As a result of higher-degree
 * interpolation with consistent values from both sides in the limit of a
 * continuous solution, the numerical flux can be seen as a control of the
 * jump of the solution from both sides to weakly impose continuity. It is
 * important to realize that a numerical flux alone cannot stabilize a
 * high-order DG method in the presence of shocks, and thus any DG method
 * must be combined with further shock-capturing techniques to handle those
 * cases. In this tutorial, we focus on wave-like solutions of the Euler
 * equations in the subsonic regime without strong discontinuities where our
 * basic scheme is sufficient.
 *   

 * 
 * Nonetheless, the numerical flux is decisive in terms of the numerical
 * dissipation of the overall scheme and influences the admissible time step
 * size with explicit Runge--Kutta methods. We consider two choices, a
 * modified Lax--Friedrichs scheme and the widely used Harten--Lax--van Leer
 * (HLL) flux. For both variants, we first need to get the velocities and
 * pressures from both sides of the interface and evaluate the physical
 * Euler flux.
 *   

 * 
 * For the local Lax--Friedrichs flux, the definition is $\hat{\mathbf{F}}
 * =\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
 * \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
 * \mathbf{n^-}$, where the factor $\lambda =
 * \max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ gives the
 * maximal wave speed and $c = \sqrt{\gamma p / \rho}$ is the speed of
 * sound. Here, we choose two modifications of that expression for reasons
 * of computational efficiency, given the small impact of the flux on the
 * solution. For the above definition of the factor $\lambda$, we would need
 * to take four square roots, two for the two velocity norms and two for the
 * speed of sound on either side. The first modification is hence to rather
 * use $\sqrt{\|\mathbf{u}\|^2+c^2}$ as an estimate of the maximal speed
 * (which is at most a factor of 2 away from the actual maximum, as shown in
 * the introduction). This allows us to pull the square root out of the
 * maximum and get away with a single square root computation. The second
 * modification is to further relax on the parameter $\lambda$---the smaller
 * it is, the smaller the dissipation factor (which is multiplied by the
 * jump in $\mathbf{w}$, which might result in a smaller or bigger
 * dissipation in the end). This allows us to fit the spectrum into the
 * stability region of the explicit Runge--Kutta integrator with bigger time
 * steps. However, we cannot make dissipation too small because otherwise
 * imaginary eigenvalues grow larger. Finally, the current conservative
 * formulation is not energy-stable in the limit of $\lambda\to 0$ as it is
 * not skew-symmetric, and would need additional measures such as split-form
 * DG schemes in that case.
 *   

 * 
 * For the HLL flux, we follow the formula from literature, introducing an
 * additional weighting of the two states from Lax--Friedrichs by a
 * parameter $s$. It is derived from the physical transport directions of
 * the Euler equations in terms of the current direction of velocity and
 * sound speed. For the velocity, we here choose a simple arithmetic average
 * which is sufficient for DG scenarios and moderate jumps in material
 * parameters.
 *   

 * 
 * Since the numerical flux is multiplied by the normal vector in the weak
 * form, we multiply by the result by the normal vector for all terms in the
 * equation. In these multiplications, the `operator*` defined above enables
 * a compact notation similar to the mathematical definition.
 *   

 * 
 * In this and the following functions, we use variable suffixes `_m` and
 * `_p` to indicate quantities derived from $\mathbf{w}^-$ and $\mathbf{w}^+$,
 * i.e., values "here" and "there" relative to the current cell when looking
 * at a neighbor cell.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim + 2, Number>
 *     euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
 *                          const Tensor<1, dim + 2, Number> &u_p,
 *                          const Tensor<1, dim, Number> &    normal)
 *   {
 *     const auto velocity_m = euler_velocity<dim>(u_m);
 *     const auto velocity_p = euler_velocity<dim>(u_p);
 * 
 *     const auto pressure_m = euler_pressure<dim>(u_m);
 *     const auto pressure_p = euler_pressure<dim>(u_p);
 * 
 *     const auto flux_m = euler_flux<dim>(u_m);
 *     const auto flux_p = euler_flux<dim>(u_p);
 * 
 *     switch (numerical_flux_type)
 *       {
 *         case lax_friedrichs_modified:
 *           {
 *             const auto lambda =
 *               0.5 * std::sqrt(std::max(velocity_p.norm_square() +
 *                                          gamma * pressure_p * (1. / u_p[0]),
 *                                        velocity_m.norm_square() +
 *                                          gamma * pressure_m * (1. / u_m[0])));
 * 
 *             return 0.5 * (flux_m * normal + flux_p * normal) +
 *                    0.5 * lambda * (u_m - u_p);
 *           }
 * 
 *         case harten_lax_vanleer:
 *           {
 *             const auto avg_velocity_normal =
 *               0.5 * ((velocity_m + velocity_p) * normal);
 *             const auto   avg_c = std::sqrt(std::abs(
 *               0.5 * gamma *
 *               (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0]))));
 *             const Number s_pos =
 *               std::max(Number(), avg_velocity_normal + avg_c);
 *             const Number s_neg =
 *               std::min(Number(), avg_velocity_normal - avg_c);
 *             const Number inverse_s = Number(1.) / (s_pos - s_neg);
 * 
 *             return inverse_s *
 *                    ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) -
 *                     s_pos * s_neg * (u_m - u_p));
 *           }
 * 
 *         default:
 *           {
 *             Assert(false, ExcNotImplemented());
 *             return {};
 *           }
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This and the next function are helper functions to provide compact
 * evaluation calls as multiple points get batched together via a
 * VectorizedArray argument (see the step-37 tutorial for details). This
 * function is used for the subsonic outflow boundary conditions where we
 * need to set the energy component to a prescribed value. The next one
 * requests the solution on all components and is used for inflow boundaries
 * where all components of the solution are set.
 * 
 * @code
 *   template <int dim, typename Number>
 *   VectorizedArray<Number>
 *   evaluate_function(const Function<dim> &                      function,
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized,
 *                     const unsigned int                         component)
 *   {
 *     VectorizedArray<Number> result;
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *       {
 *         Point<dim> p;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           p[d] = p_vectorized[d][v];
 *         result[v] = function.value(p, component);
 *       }
 *     return result;
 *   }
 * 
 * 
 *   template <int dim, typename Number, int n_components = dim + 2>
 *   Tensor<1, n_components, VectorizedArray<Number>>
 *   evaluate_function(const Function<dim> &                      function,
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized)
 *   {
 *     AssertDimension(function.n_components, n_components);
 *     Tensor<1, n_components, VectorizedArray<Number>> result;
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *       {
 *         Point<dim> p;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           p[d] = p_vectorized[d][v];
 *         for (unsigned int d = 0; d < n_components; ++d)
 *           result[d][v] = function.value(p, d);
 *       }
 *     return result;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerOperationclass"></a> 
 * <h3>The EulerOperation class</h3>
 * 

 * 
 * This class implements the evaluators for the Euler problem, in analogy to
 * the `LaplaceOperator` class of step-37 or step-59. Since the present
 * operator is non-linear and does not require a matrix interface (to be
 * handed over to preconditioners), we skip the various `vmult` functions
 * otherwise present in matrix-free operators and only implement an `apply`
 * function as well as the combination of `apply` with the required vector
 * updates for the low-storage Runge--Kutta time integrator mentioned above
 * (called `perform_stage`). Furthermore, we have added three additional
 * functions involving matrix-free routines, namely one to compute an
 * estimate of the time step scaling (that is combined with the Courant
 * number for the actual time step size) based on the velocity and speed of
 * sound in the elements, one for the projection of solutions (specializing
 * VectorTools::project() for the DG case), and one to compute the errors
 * against a possible analytical solution or norms against some background
 * state.
 *   

 * 
 * The rest of the class is similar to other matrix-free tutorials. As
 * discussed in the introduction, we provide a few functions to allow a user
 * to pass in various forms of boundary conditions on different parts of the
 * domain boundary marked by types::boundary_id variables, as well as
 * possible body forces.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   class EulerOperator
 *   {
 *   public:
 *     static constexpr unsigned int n_quadrature_points_1d = n_points_1d;
 * 
 *     EulerOperator(TimerOutput &timer_output);
 * 
 *     void reinit(const Mapping<dim> &   mapping,
 *                 const DoFHandler<dim> &dof_handler);
 * 
 *     void set_inflow_boundary(const types::boundary_id       boundary_id,
 *                              std::unique_ptr<Function<dim>> inflow_function);
 * 
 *     void set_subsonic_outflow_boundary(
 *       const types::boundary_id       boundary_id,
 *       std::unique_ptr<Function<dim>> outflow_energy);
 * 
 *     void set_wall_boundary(const types::boundary_id boundary_id);
 * 
 *     void set_body_force(std::unique_ptr<Function<dim>> body_force);
 * 
 *     void apply(const double                                      current_time,
 *                const LinearAlgebra::distributed::Vector<Number> &src,
 *                LinearAlgebra::distributed::Vector<Number> &      dst) const;
 * 
 *     void
 *     perform_stage(const Number cur_time,
 *                   const Number factor_solution,
 *                   const Number factor_ai,
 *                   const LinearAlgebra::distributed::Vector<Number> &current_ri,
 *                   LinearAlgebra::distributed::Vector<Number> &      vec_ki,
 *                   LinearAlgebra::distributed::Vector<Number> &      solution,
 *                   LinearAlgebra::distributed::Vector<Number> &next_ri) const;
 * 
 *     void project(const Function<dim> &                       function,
 *                  LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     std::array<double, 3> compute_errors(
 *       const Function<dim> &                             function,
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     double compute_cell_transport_speed(
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     void
 *     initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const;
 * 
 *   private:
 *     MatrixFree<dim, Number> data;
 * 
 *     TimerOutput &timer;
 * 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
 *       inflow_boundaries;
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
 *                                    subsonic_outflow_boundaries;
 *     std::set<types::boundary_id>   wall_boundaries;
 *     std::unique_ptr<Function<dim>> body_force;
 * 
 *     void local_apply_inverse_mass_matrix(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const;
 * 
 *     void local_apply_cell(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const;
 * 
 *     void local_apply_face(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     face_range) const;
 * 
 *     void local_apply_boundary_face(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     face_range) const;
 *   };
 * 
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer)
 *     : timer(timer)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * For the initialization of the Euler operator, we set up the MatrixFree
 * variable contained in the class. This can be done given a mapping to
 * describe possible curved boundaries as well as a DoFHandler object
 * describing the degrees of freedom. Since we use a discontinuous Galerkin
 * discretization in this tutorial program where no constraints are imposed
 * strongly on the solution field, we do not need to pass in an
 * AffineConstraints object and rather use a dummy for the
 * construction. With respect to quadrature, we want to select two different
 * ways of computing the underlying integrals: The first is a flexible one,
 * based on a template parameter `n_points_1d` (that will be assigned the
 * `n_q_points_1d` value specified at the top of this file). More accurate
 * integration is necessary to avoid the aliasing problem due to the
 * variable coefficients in the Euler operator. The second less accurate
 * quadrature formula is a tight one based on `fe_degree+1` and needed for
 * the inverse mass matrix. While that formula provides an exact inverse
 * only on affine element shapes and not on deformed elements, it enables
 * the fast inversion of the mass matrix by tensor product techniques,
 * necessary to ensure optimal computational efficiency overall.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::reinit(
 *     const Mapping<dim> &   mapping,
 *     const DoFHandler<dim> &dof_handler)
 *   {
 *     const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler};
 *     const AffineConstraints<double>            dummy;
 *     const std::vector<const AffineConstraints<double> *> constraints = {&dummy};
 *     const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
 *                                                     QGauss<1>(fe_degree + 1)};
 * 
 *     typename MatrixFree<dim, Number>::AdditionalData additional_data;
 *     additional_data.mapping_update_flags =
 *       (update_gradients | update_JxW_values | update_quadrature_points |
 *        update_values);
 *     additional_data.mapping_update_flags_inner_faces =
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors |
 *        update_values);
 *     additional_data.mapping_update_flags_boundary_faces =
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors |
 *        update_values);
 *     additional_data.tasks_parallel_scheme =
 *       MatrixFree<dim, Number>::AdditionalData::none;
 * 
 *     data.reinit(
 *       mapping, dof_handlers, constraints, quadratures, additional_data);
 *   }
 * 
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::initialize_vector(
 *     LinearAlgebra::distributed::Vector<Number> &vector) const
 *   {
 *     data.initialize_dof_vector(vector);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The subsequent four member functions are the ones that must be called from
 * outside to specify the various types of boundaries. For an inflow boundary,
 * we must specify all components in terms of density $\rho$, momentum $\rho
 * \mathbf{u}$ and energy $E$. Given this information, we then store the
 * function alongside the respective boundary id in a map member variable of
 * this class. Likewise, we proceed for the subsonic outflow boundaries (where
 * we request a function as well, which we use to retrieve the energy) and for
 * wall (no-penetration) boundaries where we impose zero normal velocity (no
 * function necessary, so we only request the boundary id). For the present
 * DG code where boundary conditions are solely applied as part of the weak
 * form (during time integration), the call to set the boundary conditions
 * can appear both before or after the `reinit()` call to this class. This
 * is different from continuous finite element codes where the boundary
 * conditions determine the content of the AffineConstraints object that is
 * sent into MatrixFree for initialization, thus requiring to be set before
 * the initialization of the matrix-free data structures.
 *   

 * 
 * The checks added in each of the four function are used to
 * ensure that boundary conditions are mutually exclusive on the various
 * parts of the boundary, i.e., that a user does not accidentally designate a
 * boundary as both an inflow and say a subsonic outflow boundary.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary(
 *     const types::boundary_id       boundary_id,
 *     std::unique_ptr<Function<dim>> inflow_function)
 *   {
 *     AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
 *                     subsonic_outflow_boundaries.end() &&
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as inflow"));
 *     AssertThrow(inflow_function->n_components == dim + 2,
 *                 ExcMessage("Expected function with dim+2 components"));
 * 
 *     inflow_boundaries[boundary_id] = std::move(inflow_function);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary(
 *     const types::boundary_id       boundary_id,
 *     std::unique_ptr<Function<dim>> outflow_function)
 *   {
 *     AssertThrow(inflow_boundaries.find(boundary_id) ==
 *                     inflow_boundaries.end() &&
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as subsonic outflow"));
 *     AssertThrow(outflow_function->n_components == dim + 2,
 *                 ExcMessage("Expected function with dim+2 components"));
 * 
 *     subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary(
 *     const types::boundary_id boundary_id)
 *   {
 *     AssertThrow(inflow_boundaries.find(boundary_id) ==
 *                     inflow_boundaries.end() &&
 *                   subsonic_outflow_boundaries.find(boundary_id) ==
 *                     subsonic_outflow_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as wall boundary"));
 * 
 *     wall_boundaries.insert(boundary_id);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_body_force(
 *     std::unique_ptr<Function<dim>> body_force)
 *   {
 *     AssertDimension(body_force->n_components, dim);
 * 
 *     this->body_force = std::move(body_force);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Localevaluators"></a> 
 * <h4>Local evaluators</h4>
 * 

 * 
 * Now we proceed to the local evaluators for the Euler problem. The
 * evaluators are relatively simple and follow what has been presented in
 * step-37, step-48, or step-59. The first notable difference is the fact
 * that we use an FEEvaluation with a non-standard number of quadrature
 * points. Whereas we previously always set the number of quadrature points
 * to equal the polynomial degree plus one (ensuring exact integration on
 * affine element shapes), we now set the number quadrature points as a
 * separate variable (e.g. the polynomial degree plus two or three halves of
 * the polynomial degree) to more accurately handle nonlinear terms. Since
 * the evaluator is fed with the appropriate loop lengths via the template
 * argument and keeps the number of quadrature points in the whole cell in
 * the variable FEEvaluation::n_q_points, we now automatically operate on
 * the more accurate formula without further changes.
 *   

 * 
 * The second difference is due to the fact that we are now evaluating a
 * multi-component system, as opposed to the scalar systems considered
 * previously. The matrix-free framework provides several ways to handle the
 * multi-component case. The variant shown here utilizes an FEEvaluation
 * object with multiple components embedded into it, specified by the fourth
 * template argument `dim + 2` for the components in the Euler system. As a
 * consequence, the return type of FEEvaluation::get_value() is not a scalar
 * any more (that would return a VectorizedArray type, collecting data from
 * several elements), but a Tensor of `dim+2` components. The functionality
 * is otherwise similar to the scalar case; it is handled by a template
 * specialization of a base class, called FEEvaluationAccess. An alternative
 * variant would have been to use several FEEvaluation objects, a scalar one
 * for the density, a vector-valued one with `dim` components for the
 * momentum, and another scalar evaluator for the energy. To ensure that
 * those components point to the correct part of the solution, the
 * constructor of FEEvaluation takes three optional integer arguments after
 * the required MatrixFree field, namely the number of the DoFHandler for
 * multi-DoFHandler systems (taking the first by default), the number of the
 * quadrature point in case there are multiple Quadrature objects (see more
 * below), and as a third argument the component within a vector system. As
 * we have a single vector for all components, we would go with the third
 * argument, and set it to `0` for the density, `1` for the vector-valued
 * momentum, and `dim+1` for the energy slot. FEEvaluation then picks the
 * appropriate subrange of the solution vector during
 * FEEvaluationBase::read_dof_values() and
 * FEEvaluation::distributed_local_to_global() or the more compact
 * FEEvaluation::gather_evaluate() and FEEvaluation::integrate_scatter()
 * calls.
 *   

 * 
 * When it comes to the evaluation of the body force vector, we distinguish
 * between two cases for efficiency reasons: In case we have a constant
 * function (derived from Functions::ConstantFunction), we can precompute
 * the value outside the loop over quadrature points and simply use the
 * value everywhere. For a more general function, we instead need to call
 * the `evaluate_function()` method we provided above; this path is more
 * expensive because we need to access the memory associated with the
 * quadrature point data.
 *   

 * 
 * The rest follows the other tutorial programs. Since we have implemented
 * all physics for the Euler equations in the separate `euler_flux()`
 * function, all we have to do here is to call this function
 * given the current solution evaluated at quadrature points, returned by
 * `phi.get_value(q)`, and tell the FEEvaluation object to queue the flux
 * for testing it by the gradients of the shape functions (which is a Tensor
 * of outer `dim+2` components, each holding a tensor of `dim` components
 * for the $x,y,z$ component of the Euler flux). One final thing worth
 * mentioning is the order in which we queue the data for testing by the
 * value of the test function, `phi.submit_value()`, in case we are given an
 * external function: We must do this after calling `phi.get_value(q)`,
 * because `get_value()` (reading the solution) and `submit_value()`
 * (queuing the value for multiplication by the test function and summation
 * over quadrature points) access the same underlying data field. Here it
 * would be easy to achieve also without temporary variable `w_q` since
 * there is no mixing between values and gradients. For more complicated
 * setups, one has to first copy out e.g. both the value and gradient at a
 * quadrature point and then queue results again by
 * FEEvaluationBase::submit_value() and FEEvaluationBase::submit_gradient().
 *   

 * 
 * As a final note, we mention that we do not use the first MatrixFree
 * argument of this function, which is a call-back from MatrixFree::loop().
 * The interfaces imposes the present list of arguments, but since we are in
 * a member function where the MatrixFree object is already available as the
 * `data` variable, we stick with that to avoid confusion.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_cell(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data);
 * 
 *     Tensor<1, dim, VectorizedArray<Number>> constant_body_force;
 *     const Functions::ConstantFunction<dim> *constant_function =
 *       dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());
 * 
 *     if (constant_function)
 *       constant_body_force = evaluate_function<dim, Number, dim>(
 *         *constant_function, Point<dim, VectorizedArray<Number>>());
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto w_q = phi.get_value(q);
 *             phi.submit_gradient(euler_flux<dim>(w_q), q);
 *             if (body_force.get() != nullptr)
 *               {
 *                 const Tensor<1, dim, VectorizedArray<Number>> force =
 *                   constant_function ? constant_body_force :
 *                                       evaluate_function<dim, Number, dim>(
 *                                         *body_force, phi.quadrature_point(q));
 * 
 *                 Tensor<1, dim + 2, VectorizedArray<Number>> forcing;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   forcing[d + 1] = w_q[0] * force[d];
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   forcing[dim + 1] += force[d] * w_q[d + 1];
 * 
 *                 phi.submit_value(forcing, q);
 *               }
 *           }
 * 
 *         phi.integrate_scatter(((body_force.get() != nullptr) ?
 *                                  EvaluationFlags::values :
 *                                  EvaluationFlags::nothing) |
 *                                 EvaluationFlags::gradients,
 *                               dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function concerns the computation of integrals on interior
 * faces, where we need evaluators from both cells adjacent to the face. We
 * associate the variable `phi_m` with the solution component $\mathbf{w}^-$
 * and the variable `phi_p` with the solution component $\mathbf{w}^+$. We
 * distinguish the two sides in the constructor of FEFaceEvaluation by the
 * second argument, with `true` for the interior side and `false` for the
 * exterior side, with interior and exterior denoting the orientation with
 * respect to the normal vector.
 *   

 * 
 * Note that the calls FEFaceEvaluation::gather_evaluate() and
 * FEFaceEvaluation::integrate_scatter() combine the access to the vectors
 * and the sum factorization parts. This combined operation not only saves a
 * line of code, but also contains an important optimization: Given that we
 * use a nodal basis in terms of the Lagrange polynomials in the points of
 * the Gauss-Lobatto quadrature formula, only $(p+1)^{d-1}$ out of the
 * $(p+1)^d$ basis functions evaluate to non-zero on each face. Thus, the
 * evaluator only accesses the necessary data in the vector and skips the
 * parts which are multiplied by zero. If we had first read the vector, we
 * would have needed to load all data from the vector, as the call in
 * isolation would not know what data is required in subsequent
 * operations. If the subsequent FEFaceEvaluation::evaluate() call requests
 * values and derivatives, indeed all $(p+1)^d$ vector entries for each
 * component are needed, as the normal derivative is nonzero for all basis
 * functions.
 *   

 * 
 * The arguments to the evaluators as well as the procedure is similar to
 * the cell evaluation. We again use the more accurate (over-)integration
 * scheme due to the nonlinear terms, specified as the third template
 * argument in the list. At the quadrature points, we then go to our
 * free-standing function for the numerical flux. It receives the solution
 * evaluated at quadrature points from both sides (i.e., $\mathbf{w}^-$ and
 * $\mathbf{w}^+$), as well as the normal vector onto the minus side. As
 * explained above, the numerical flux is already multiplied by the normal
 * vector from the minus side. We need to switch the sign because the
 * boundary term comes with a minus sign in the weak form derived in the
 * introduction. The flux is then queued for testing both on the minus sign
 * and on the plus sign, with switched sign as the normal vector from the
 * plus side is exactly opposed to the one from the minus side.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_face(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data,
 *                                                                       true);
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data,
 *                                                                       false);
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 *         phi_p.reinit(face);
 *         phi_p.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         phi_m.reinit(face);
 *         phi_m.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
 *           {
 *             const auto numerical_flux =
 *               euler_numerical_flux<dim>(phi_m.get_value(q),
 *                                         phi_p.get_value(q),
 *                                         phi_m.get_normal_vector(q));
 *             phi_m.submit_value(-numerical_flux, q);
 *             phi_p.submit_value(numerical_flux, q);
 *           }
 * 
 *         phi_p.integrate_scatter(EvaluationFlags::values, dst);
 *         phi_m.integrate_scatter(EvaluationFlags::values, dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For faces located at the boundary, we need to impose the appropriate
 * boundary conditions. In this tutorial program, we implement four cases as
 * mentioned above. (A fifth case, for supersonic outflow conditions is
 * discussed in the "Results" section below.) The discontinuous Galerkin
 * method imposes boundary conditions not as constraints, but only
 * weakly. Thus, the various conditions are imposed by finding an appropriate
 * <i>exterior</i> quantity $\mathbf{w}^+$ that is then handed to the
 * numerical flux function also used for the interior faces. In essence,
 * we "pretend" a state on the outside of the domain in such a way that
 * if that were reality, the solution of the PDE would satisfy the boundary
 * conditions we want.
 *   

 * 
 * For wall boundaries, we need to impose a no-normal-flux condition on the
 * momentum variable, whereas we use a Neumann condition for the density and
 * energy with $\rho^+ = \rho^-$ and $E^+ = E^-$. To achieve the no-normal
 * flux condition, we set the exterior values to the interior values and
 * subtract two times the velocity in wall-normal direction, i.e., in the
 * direction of the normal vector.
 *   

 * 
 * For inflow boundaries, we simply set the given Dirichlet data
 * $\mathbf{w}_\mathrm{D}$ as a boundary value. An alternative would have been
 * to use $\mathbf{w}^+ = -\mathbf{w}^- + 2 \mathbf{w}_\mathrm{D}$, the
 * so-called mirror principle.
 *   

 * 
 * The imposition of outflow is essentially a Neumann condition, i.e.,
 * setting $\mathbf{w}^+ = \mathbf{w}^-$. For the case of subsonic outflow,
 * we still need to impose a value for the energy, which we derive from the
 * respective function. A special step is needed for the case of
 * <i>backflow</i>, i.e., the case where there is a momentum flux into the
 * domain on the Neumann portion. According to the literature (a fact that can
 * be derived by appropriate energy arguments), we must switch to another
 * variant of the flux on inflow parts, see Gravemeier, Comerford,
 * Yoshihara, Ismail, Wall, "A novel formulation for Neumann inflow
 * conditions in biomechanics", Int. J. Numer. Meth. Biomed. Eng., vol. 28
 * (2012). Here, the momentum term needs to be added once again, which
 * corresponds to removing the flux contribution on the momentum
 * variables. We do this in a post-processing step, and only for the case
 * when we both are at an outflow boundary and the dot product between the
 * normal vector and the momentum (or, equivalently, velocity) is
 * negative. As we work on data of several quadrature points at once for
 * SIMD vectorizations, we here need to explicitly loop over the array
 * entries of the SIMD array.
 *   

 * 
 * In the implementation below, we check for the various types
 * of boundaries at the level of quadrature points. Of course, we could also
 * have moved the decision out of the quadrature point loop and treat entire
 * faces as of the same kind, which avoids some map/set lookups in the inner
 * loop over quadrature points. However, the loss of efficiency is hardly
 * noticeable, so we opt for the simpler code here. Also note that the final
 * `else` clause will catch the case when some part of the boundary was not
 * assigned any boundary condition via `EulerOperator::set_..._boundary(...)`.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_boundary_face(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true);
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 *         phi.reinit(face);
 *         phi.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto w_m    = phi.get_value(q);
 *             const auto normal = phi.get_normal_vector(q);
 * 
 *             auto rho_u_dot_n = w_m[1] * normal[0];
 *             for (unsigned int d = 1; d < dim; ++d)
 *               rho_u_dot_n += w_m[1 + d] * normal[d];
 * 
 *             bool at_outflow = false;
 * 
 *             Tensor<1, dim + 2, VectorizedArray<Number>> w_p;
 *             const auto boundary_id = data.get_boundary_id(face);
 *             if (wall_boundaries.find(boundary_id) != wall_boundaries.end())
 *               {
 *                 w_p[0] = w_m[0];
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
 *                 w_p[dim + 1] = w_m[dim + 1];
 *               }
 *             else if (inflow_boundaries.find(boundary_id) !=
 *                      inflow_boundaries.end())
 *               w_p =
 *                 evaluate_function(*inflow_boundaries.find(boundary_id)->second,
 *                                   phi.quadrature_point(q));
 *             else if (subsonic_outflow_boundaries.find(boundary_id) !=
 *                      subsonic_outflow_boundaries.end())
 *               {
 *                 w_p          = w_m;
 *                 w_p[dim + 1] = evaluate_function(
 *                   *subsonic_outflow_boundaries.find(boundary_id)->second,
 *                   phi.quadrature_point(q),
 *                   dim + 1);
 *                 at_outflow = true;
 *               }
 *             else
 *               AssertThrow(false,
 *                           ExcMessage("Unknown boundary id, did "
 *                                      "you set a boundary condition for "
 *                                      "this part of the domain boundary?"));
 * 
 *             auto flux = euler_numerical_flux<dim>(w_m, w_p, normal);
 * 
 *             if (at_outflow)
 *               for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *                 {
 *                   if (rho_u_dot_n[v] < -1e-12)
 *                     for (unsigned int d = 0; d < dim; ++d)
 *                       flux[d + 1][v] = 0.;
 *                 }
 * 
 *             phi.submit_value(-flux, q);
 *           }
 * 
 *         phi.integrate_scatter(EvaluationFlags::values, dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function implements the inverse mass matrix operation. The
 * algorithms and rationale have been discussed extensively in the
 * introduction, so we here limit ourselves to the technicalities of the
 * MatrixFreeOperators::CellwiseInverseMassMatrix class. It does similar
 * operations as the forward evaluation of the mass matrix, except with a
 * different interpolation matrix, representing the inverse $S^{-1}$
 * factors. These represent a change of basis from the specified basis (in
 * this case, the Lagrange basis in the points of the Gauss--Lobatto
 * quadrature formula) to the Lagrange basis in the points of the Gauss
 * quadrature formula. In the latter basis, we can apply the inverse of the
 * point-wise `JxW` factor, i.e., the quadrature weight times the
 * determinant of the Jacobian of the mapping from reference to real
 * coordinates. Once this is done, the basis is changed back to the nodal
 * Gauss-Lobatto basis again. All of these operations are done by the
 * `apply()` function below. What we need to provide is the local fields to
 * operate on (which we extract from the global vector by an FEEvaluation
 * object) and write the results back to the destination vector of the mass
 * matrix operation.
 *   

 * 
 * One thing to note is that we added two integer arguments (that are
 * optional) to the constructor of FEEvaluation, the first being 0
 * (selecting among the DoFHandler in multi-DoFHandler systems; here, we
 * only have one) and the second being 1 to make the quadrature formula
 * selection. As we use the quadrature formula 0 for the over-integration of
 * nonlinear terms, we use the formula 1 with the default $p+1$ (or
 * `fe_degree+1` in terms of the variable name) points for the mass
 * matrix. This leads to square contributions to the mass matrix and ensures
 * exact integration, as explained in the introduction.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_inverse_mass_matrix(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
 *       inverse(phi);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.read_dof_values(src);
 * 
 *         inverse.apply(phi.begin_dof_values(), phi.begin_dof_values());
 * 
 *         phi.set_dof_values(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theapplyandrelatedfunctions"></a> 
 * <h4>The apply() and related functions</h4>
 * 

 * 
 * We now come to the function which implements the evaluation of the Euler
 * operator as a whole, i.e., $\mathcal M^{-1} \mathcal L(t, \mathbf{w})$,
 * calling into the local evaluators presented above. The steps should be
 * clear from the previous code. One thing to note is that we need to adjust
 * the time in the functions we have associated with the various parts of
 * the boundary, in order to be consistent with the equation in case the
 * boundary data is time-dependent. Then, we call MatrixFree::loop() to
 * perform the cell and face integrals, including the necessary ghost data
 * exchange in the `src` vector. The seventh argument to the function,
 * `true`, specifies that we want to zero the `dst` vector as part of the
 * loop, before we start accumulating integrals into it. This variant is
 * preferred over explicitly calling `dst = 0.;` before the loop as the
 * zeroing operation is done on a subrange of the vector in parts that are
 * written by the integrals nearby. This enhances data locality and allows
 * for caching, saving one roundtrip of vector data to main memory and
 * enhancing performance. The last two arguments to the loop determine which
 * data is exchanged: Since we only access the values of the shape functions
 * one faces, typical of first-order hyperbolic problems, and since we have
 * a nodal basis with nodes at the reference element surface, we only need
 * to exchange those parts. This again saves precious memory bandwidth.
 *   

 * 
 * Once the spatial operator $\mathcal L$ is applied, we need to make a
 * second round and apply the inverse mass matrix. Here, we call
 * MatrixFree::cell_loop() since only cell integrals appear. The cell loop
 * is cheaper than the full loop as access only goes to the degrees of
 * freedom associated with the locally owned cells, which is simply the
 * locally owned degrees of freedom for DG discretizations. Thus, no ghost
 * exchange is needed here.
 *   

 * 
 * Around all these functions, we put timer scopes to record the
 * computational time for statistics about the contributions of the various
 * parts.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::apply(
 *     const double                                      current_time,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     LinearAlgebra::distributed::Vector<Number> &      dst) const
 *   {
 *     {
 *       TimerOutput::Scope t(timer, "apply - integrals");
 * 
 *       for (auto &i : inflow_boundaries)
 *         i.second->set_time(current_time);
 *       for (auto &i : subsonic_outflow_boundaries)
 *         i.second->set_time(current_time);
 * 
 *       data.loop(&EulerOperator::local_apply_cell,
 *                 &EulerOperator::local_apply_face,
 *                 &EulerOperator::local_apply_boundary_face,
 *                 this,
 *                 dst,
 *                 src,
 *                 true,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values);
 *     }
 * 
 *     {
 *       TimerOutput::Scope t(timer, "apply - inverse mass");
 * 
 *       data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix,
 *                      this,
 *                      dst,
 *                      dst);
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Let us move to the function that does an entire stage of a Runge--Kutta
 * update. It calls EulerOperator::apply() followed by some updates
 * to the vectors, namely `next_ri = solution + factor_ai * k_i` and
 * `solution += factor_solution * k_i`. Rather than performing these
 * steps through the vector interfaces, we here present an alternative
 * strategy that is faster on cache-based architectures. As the memory
 * consumed by the vectors is often much larger than what fits into caches,
 * the data has to effectively come from the slow RAM memory. The situation
 * can be improved by loop fusion, i.e., performing both the updates to
 * `next_ki` and `solution` within a single sweep. In that case, we would
 * read the two vectors `rhs` and `solution` and write into `next_ki` and
 * `solution`, compared to at least 4 reads and two writes in the baseline
 * case. Here, we go one step further and perform the loop immediately when
 * the mass matrix inversion has finished on a part of the
 * vector. MatrixFree::cell_loop() provides a mechanism to attach an
 * `std::function` both before the loop over cells first touches a vector
 * entry (which we do not use here, but is e.g. used for zeroing the vector)
 * and a second `std::function` to be called after the loop last touches
 * an entry. The callback is in form of a range over the given vector (in
 * terms of the local index numbering in the MPI universe) that can be
 * addressed by `local_element()` functions.
 *   

 * 
 * For this second callback, we create a lambda that works on a range and
 * write the respective update on this range. Ideally, we would add the
 * `DEAL_II_OPENMP_SIMD_PRAGMA` before the local loop to suggest to the
 * compiler to SIMD parallelize this loop (which means in practice that we
 * ensure that there is no overlap, also called aliasing, between the index
 * ranges of the pointers we use inside the loops). It turns out that at the
 * time of this writing, GCC 7.2 fails to compile an OpenMP pragma inside a
 * lambda function, so we comment this pragma out below. If your compiler is
 * newer, you should be able to uncomment these lines again.
 *   

 * 
 * Note that we select a different code path for the last
 * Runge--Kutta stage when we do not need to update the `next_ri`
 * vector. This strategy gives a considerable speedup. Whereas the inverse
 * mass matrix and vector updates take more than 60% of the computational
 * time with default vector updates on a 40-core machine, the percentage is
 * around 35% with the more optimized variant. In other words, this is a
 * speedup of around a third.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::perform_stage(
 *     const Number                                      current_time,
 *     const Number                                      factor_solution,
 *     const Number                                      factor_ai,
 *     const LinearAlgebra::distributed::Vector<Number> &current_ri,
 *     LinearAlgebra::distributed::Vector<Number> &      vec_ki,
 *     LinearAlgebra::distributed::Vector<Number> &      solution,
 *     LinearAlgebra::distributed::Vector<Number> &      next_ri) const
 *   {
 *     {
 *       TimerOutput::Scope t(timer, "rk_stage - integrals L_h");
 * 
 *       for (auto &i : inflow_boundaries)
 *         i.second->set_time(current_time);
 *       for (auto &i : subsonic_outflow_boundaries)
 *         i.second->set_time(current_time);
 * 
 *       data.loop(&EulerOperator::local_apply_cell,
 *                 &EulerOperator::local_apply_face,
 *                 &EulerOperator::local_apply_boundary_face,
 *                 this,
 *                 vec_ki,
 *                 current_ri,
 *                 true,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values);
 *     }
 * 
 * 
 *     {
 *       TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd");
 *       data.cell_loop(
 *         &EulerOperator::local_apply_inverse_mass_matrix,
 *         this,
 *         next_ri,
 *         vec_ki,
 *         std::function<void(const unsigned int, const unsigned int)>(),
 *         [&](const unsigned int start_range, const unsigned int end_range) {
 *           const Number ai = factor_ai;
 *           const Number bi = factor_solution;
 *           if (ai == Number())
 *             {
 *               /* DEAL_II_OPENMP_SIMD_PRAGMA */
 *               for (unsigned int i = start_range; i < end_range; ++i)
 *                 {
 *                   const Number k_i          = next_ri.local_element(i);
 *                   const Number sol_i        = solution.local_element(i);
 *                   solution.local_element(i) = sol_i + bi * k_i;
 *                 }
 *             }
 *           else
 *             {
 *               /* DEAL_II_OPENMP_SIMD_PRAGMA */
 *               for (unsigned int i = start_range; i < end_range; ++i)
 *                 {
 *                   const Number k_i          = next_ri.local_element(i);
 *                   const Number sol_i        = solution.local_element(i);
 *                   solution.local_element(i) = sol_i + bi * k_i;
 *                   next_ri.local_element(i)  = sol_i + ai * k_i;
 *                 }
 *             }
 *         });
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Having discussed the implementation of the functions that deal with
 * advancing the solution by one time step, let us now move to functions
 * that implement other, ancillary operations. Specifically, these are
 * functions that compute projections, evaluate errors, and compute the speed
 * of information transport on a cell.
 *   

 * 
 * The first of these functions is essentially equivalent to
 * VectorTools::project(), just much faster because it is specialized for DG
 * elements where there is no need to set up and solve a linear system, as
 * each element has independent basis functions. The reason why we show the
 * code here, besides a small speedup of this non-critical operation, is that
 * it shows additional functionality provided by
 * MatrixFreeOperators::CellwiseInverseMassMatrix.
 *   

 * 
 * The projection operation works as follows: If we denote the matrix of
 * shape functions evaluated at quadrature points by $S$, the projection on
 * cell $K$ is an operation of the form $\underbrace{S J^K S^\mathrm
 * T}_{\mathcal M^K} \mathbf{w}^K = S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$, where $J^K$ is the diagonal
 * matrix containing the determinant of the Jacobian times the quadrature
 * weight (JxW), $\mathcal M^K$ is the cell-wise mass matrix, and
 * $\tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ is the evaluation of the
 * field to be projected onto quadrature points. (In reality the matrix $S$
 * has additional structure through the tensor product, as explained in the
 * introduction.) This system can now equivalently be written as
 * $\mathbf{w}^K = \left(S J^K S^\mathrm T\right)^{-1} S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q} = S^{-\mathrm T}
 * \left(J^K\right)^{-1} S^{-1} S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. Now, the term $S^{-1} S$ and
 * then $\left(J^K\right)^{-1} J^K$ cancel, resulting in the final
 * expression $\mathbf{w}^K = S^{-\mathrm T}
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. This operation is
 * implemented by
 * MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis().
 * The name is derived from the fact that this projection is simply
 * the multiplication by $S^{-\mathrm T}$, a basis change from the
 * nodal basis in the points of the Gaussian quadrature to the given finite
 * element basis. Note that we call FEEvaluation::set_dof_values() to write
 * the result into the vector, overwriting previous content, rather than
 * accumulating the results as typical in integration tasks -- we can do
 * this because every vector entry has contributions from only a single
 * cell for discontinuous Galerkin discretizations.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::project(
 *     const Function<dim> &                       function,
 *     LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
 *       inverse(phi);
 *     solution.zero_out_ghost_values();
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_dof_value(evaluate_function(function,
 *                                                  phi.quadrature_point(q)),
 *                                q);
 *         inverse.transform_from_q_points_to_basis(dim + 2,
 *                                                  phi.begin_dof_values(),
 *                                                  phi.begin_dof_values());
 *         phi.set_dof_values(solution);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function again repeats functionality also provided by the
 * deal.II library, namely VectorTools::integrate_difference(). We here show
 * the explicit code to highlight how the vectorization across several cells
 * works and how to accumulate results via that interface: Recall that each
 * <i>lane</i> of the vectorized array holds data from a different cell. By
 * the loop over all cell batches that are owned by the current MPI process,
 * we could then fill a VectorizedArray of results; to obtain a global sum,
 * we would need to further go on and sum across the entries in the SIMD
 * array. However, such a procedure is not stable as the SIMD array could in
 * fact not hold valid data for all its lanes. This happens when the number
 * of locally owned cells is not a multiple of the SIMD width. To avoid
 * invalid data, we must explicitly skip those invalid lanes when accessing
 * the data. While one could imagine that we could make it work by simply
 * setting the empty lanes to zero (and thus, not contribute to a sum), the
 * situation is more complicated than that: What if we were to compute a
 * velocity out of the momentum? Then, we would need to divide by the
 * density, which is zero -- the result would consequently be NaN and
 * contaminate the result. This trap is avoided by accumulating the results
 * from the valid SIMD range as we loop through the cell batches, using the
 * function MatrixFree::n_active_entries_per_cell_batch() to give us the
 * number of lanes with valid data. It equals VectorizedArray::size() on
 * most cells, but can be less on the last cell batch if the number of cells
 * has a remainder compared to the SIMD width.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors(
 *     const Function<dim> &                             function,
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     TimerOutput::Scope t(timer, "compute errors");
 *     double             errors_squared[3] = {};
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0);
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(solution, EvaluationFlags::values);
 *         VectorizedArray<Number> local_errors_squared[3] = {};
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto error =
 *               evaluate_function(function, phi.quadrature_point(q)) -
 *               phi.get_value(q);
 *             const auto JxW = phi.JxW(q);
 * 
 *             local_errors_squared[0] += error[0] * error[0] * JxW;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW;
 *             local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW;
 *           }
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
 *              ++v)
 *           for (unsigned int d = 0; d < 3; ++d)
 *             errors_squared[d] += local_errors_squared[d][v];
 *       }
 * 
 *     Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);
 * 
 *     std::array<double, 3> errors;
 *     for (unsigned int d = 0; d < 3; ++d)
 *       errors[d] = std::sqrt(errors_squared[d]);
 * 
 *     return errors;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This final function of the EulerOperator class is used to estimate the
 * transport speed, scaled by the mesh size, that is relevant for setting
 * the time step size in the explicit time integrator. In the Euler
 * equations, there are two speeds of transport, namely the convective
 * velocity $\mathbf{u}$ and the propagation of sound waves with sound
 * speed $c = \sqrt{\gamma p/\rho}$ relative to the medium moving at
 * velocity $\mathbf u$.
 *   

 * 
 * In the formula for the time step size, we are interested not by
 * these absolute speeds, but by the amount of time it takes for
 * information to cross a single cell. For information transported along with
 * the medium, $\mathbf u$ is scaled by the mesh size,
 * so an estimate of the maximal velocity can be obtained by computing
 * $\|J^{-\mathrm T} \mathbf{u}\|_\infty$, where $J$ is the Jacobian of the
 * transformation from real to the reference domain. Note that
 * FEEvaluationBase::inverse_jacobian() returns the inverse and transpose
 * Jacobian, representing the metric term from real to reference
 * coordinates, so we do not need to transpose it again. We store this limit
 * in the variable `convective_limit` in the code below.
 *   

 * 
 * The sound propagation is isotropic, so we need to take mesh sizes in any
 * direction into account. The appropriate mesh size scaling is then given
 * by the minimal singular value of $J$ or, equivalently, the maximal
 * singular value of $J^{-1}$. Note that one could approximate this quantity
 * by the minimal distance between vertices of a cell when ignoring curved
 * cells. To get the maximal singular value of the Jacobian, the general
 * strategy would be some LAPACK function. Since all we need here is an
 * estimate, we can avoid the hassle of decomposing a tensor of
 * VectorizedArray numbers into several matrices and go into an (expensive)
 * eigenvalue function without vectorization, and instead use a few
 * iterations (five in the code below) of the power method applied to
 * $J^{-1}J^{-\mathrm T}$. The speed of convergence of this method depends
 * on the ratio of the largest to the next largest eigenvalue and the
 * initial guess, which is the vector of all ones. This might suggest that
 * we get slow convergence on cells close to a cube shape where all
 * lengths are almost the same. However, this slow convergence means that
 * the result will sit between the two largest singular values, which both
 * are close to the maximal value anyway. In all other cases, convergence
 * will be quick. Thus, we can merely hardcode 5 iterations here and be
 * confident that the result is good.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed(
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     TimerOutput::Scope t(timer, "compute transport speed");
 *     Number             max_transport = 0;
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(solution, EvaluationFlags::values);
 *         VectorizedArray<Number> local_max = 0.;
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto solution = phi.get_value(q);
 *             const auto velocity = euler_velocity<dim>(solution);
 *             const auto pressure = euler_pressure<dim>(solution);
 * 
 *             const auto inverse_jacobian = phi.inverse_jacobian(q);
 *             const auto convective_speed = inverse_jacobian * velocity;
 *             VectorizedArray<Number> convective_limit = 0.;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               convective_limit =
 *                 std::max(convective_limit, std::abs(convective_speed[d]));
 * 
 *             const auto speed_of_sound =
 *               std::sqrt(gamma * pressure * (1. / solution[0]));
 * 
 *             Tensor<1, dim, VectorizedArray<Number>> eigenvector;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               eigenvector[d] = 1.;
 *             for (unsigned int i = 0; i < 5; ++i)
 *               {
 *                 eigenvector = transpose(inverse_jacobian) *
 *                               (inverse_jacobian * eigenvector);
 *                 VectorizedArray<Number> eigenvector_norm = 0.;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   eigenvector_norm =
 *                     std::max(eigenvector_norm, std::abs(eigenvector[d]));
 *                 eigenvector /= eigenvector_norm;
 *               }
 *             const auto jac_times_ev   = inverse_jacobian * eigenvector;
 *             const auto max_eigenvalue = std::sqrt(
 *               (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector));
 *             local_max =
 *               std::max(local_max,
 *                        max_eigenvalue * speed_of_sound + convective_limit);
 *           }
 * 
 * @endcode
 * 
 * Similarly to the previous function, we must make sure to accumulate
 * speed only on the valid cells of a cell batch.
 * 
 * @code
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
 *              ++v)
 *           for (unsigned int d = 0; d < 3; ++d)
 *             max_transport = std::max(max_transport, local_max[v]);
 *       }
 * 
 *     max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);
 * 
 *     return max_transport;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerProblemclass"></a> 
 * <h3>The EulerProblem class</h3>
 * 

 * 
 * This class combines the EulerOperator class with the time integrator and
 * the usual global data structures such as FiniteElement and DoFHandler, to
 * actually run the simulations of the Euler problem.
 *   

 * 
 * The member variables are a triangulation, a finite element, a mapping (to
 * create high-order curved surfaces, see e.g. step-10), and a DoFHandler to
 * describe the degrees of freedom. In addition, we keep an instance of the
 * EulerOperator described above around, which will do all heavy lifting in
 * terms of integrals, and some parameters for time integration like the
 * current time or the time step size.
 *   

 * 
 * Furthermore, we use a PostProcessor instance to write some additional
 * information to the output file, in similarity to what was done in
 * step-33. The interface of the DataPostprocessor class is intuitive,
 * requiring us to provide information about what needs to be evaluated
 * (typically only the values of the solution, except for the Schlieren plot
 * that we only enable in 2D where it makes sense), and the names of what
 * gets evaluated. Note that it would also be possible to extract most
 * information by calculator tools within visualization programs such as
 * ParaView, but it is so much more convenient to do it already when writing
 * the output.
 * 
 * @code
 *   template <int dim>
 *   class EulerProblem
 *   {
 *   public:
 *     EulerProblem();
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid_and_dofs();
 * 
 *     void output_results(const unsigned int result_number);
 * 
 *     LinearAlgebra::distributed::Vector<Number> solution;
 * 
 *     ConditionalOStream pcout;
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 * 
 *     FESystem<dim>        fe;
 *     MappingQGeneric<dim> mapping;
 *     DoFHandler<dim>      dof_handler;
 * 
 *     TimerOutput timer;
 * 
 *     EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator;
 * 
 *     double time, time_step;
 * 
 *     class Postprocessor : public DataPostprocessor<dim>
 *     {
 *     public:
 *       Postprocessor();
 * 
 *       virtual void evaluate_vector_field(
 *         const DataPostprocessorInputs::Vector<dim> &inputs,
 *         std::vector<Vector<double>> &computed_quantities) const override;
 * 
 *       virtual std::vector<std::string> get_names() const override;
 * 
 *       virtual std::vector<
 *         DataComponentInterpretation::DataComponentInterpretation>
 *       get_data_component_interpretation() const override;
 * 
 *       virtual UpdateFlags get_needed_update_flags() const override;
 * 
 *     private:
 *       const bool do_schlieren_plot;
 *     };
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   EulerProblem<dim>::Postprocessor::Postprocessor()
 *     : do_schlieren_plot(dim == 2)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * For the main evaluation of the field variables, we first check that the
 * lengths of the arrays equal the expected values (the lengths `2*dim+4` or
 * `2*dim+5` are derived from the sizes of the names we specify in the
 * get_names() function below). Then we loop over all evaluation points and
 * fill the respective information: First we fill the primal solution
 * variables of density $\rho$, momentum $\rho \mathbf{u}$ and energy $E$,
 * then we compute the derived velocity $\mathbf u$, the pressure $p$, the
 * speed of sound $c=\sqrt{\gamma p / \rho}$, as well as the Schlieren plot
 * showing $s = |\nabla \rho|^2$ in case it is enabled. (See step-69 for
 * another example where we create a Schlieren plot.)
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::Postprocessor::evaluate_vector_field(
 *     const DataPostprocessorInputs::Vector<dim> &inputs,
 *     std::vector<Vector<double>> &               computed_quantities) const
 *   {
 *     const unsigned int n_evaluation_points = inputs.solution_values.size();
 * 
 *     if (do_schlieren_plot == true)
 *       Assert(inputs.solution_gradients.size() == n_evaluation_points,
 *              ExcInternalError());
 * 
 *     Assert(computed_quantities.size() == n_evaluation_points,
 *            ExcInternalError());
 *     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
 *     Assert(computed_quantities[0].size() ==
 *              dim + 2 + (do_schlieren_plot == true ? 1 : 0),
 *            ExcInternalError());
 * 
 *     for (unsigned int q = 0; q < n_evaluation_points; ++q)
 *       {
 *         Tensor<1, dim + 2> solution;
 *         for (unsigned int d = 0; d < dim + 2; ++d)
 *           solution[d] = inputs.solution_values[q](d);
 * 
 *         const double         density  = solution[0];
 *         const Tensor<1, dim> velocity = euler_velocity<dim>(solution);
 *         const double         pressure = euler_pressure<dim>(solution);
 * 
 *         for (unsigned int d = 0; d < dim; ++d)
 *           computed_quantities[q](d) = velocity[d];
 *         computed_quantities[q](dim)     = pressure;
 *         computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density);
 * 
 *         if (do_schlieren_plot == true)
 *           computed_quantities[q](dim + 2) =
 *             inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0];
 *       }
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const
 *   {
 *     std::vector<std::string> names;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       names.emplace_back("velocity");
 *     names.emplace_back("pressure");
 *     names.emplace_back("speed_of_sound");
 * 
 *     if (do_schlieren_plot == true)
 *       names.emplace_back("schlieren_plot");
 * 
 *     return names;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For the interpretation of quantities, we have scalar density, energy,
 * pressure, speed of sound, and the Schlieren plot, and vectors for the
 * momentum and the velocity.
 * 
 * @code
 *   template <int dim>
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *   EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const
 *   {
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       interpretation;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       interpretation.push_back(
 *         DataComponentInterpretation::component_is_part_of_vector);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 * 
 *     if (do_schlieren_plot == true)
 *       interpretation.push_back(
 *         DataComponentInterpretation::component_is_scalar);
 * 
 *     return interpretation;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * With respect to the necessary update flags, we only need the values for
 * all quantities but the Schlieren plot, which is based on the density
 * gradient.
 * 
 * @code
 *   template <int dim>
 *   UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const
 *   {
 *     if (do_schlieren_plot == true)
 *       return update_values | update_gradients;
 *     else
 *       return update_values;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The constructor for this class is unsurprising: We set up a parallel
 * triangulation based on the `MPI_COMM_WORLD` communicator, a vector finite
 * element with `dim+2` components for density, momentum, and energy, a
 * high-order mapping of the same degree as the underlying finite element,
 * and initialize the time and time step to zero.
 * 
 * @code
 *   template <int dim>
 *   EulerProblem<dim>::EulerProblem()
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 * #ifdef DEAL_II_WITH_P4EST
 *     , triangulation(MPI_COMM_WORLD)
 * #endif
 *     , fe(FE_DGQ<dim>(fe_degree), dim + 2)
 *     , mapping(fe_degree)
 *     , dof_handler(triangulation)
 *     , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
 *     , euler_operator(timer)
 *     , time(0)
 *     , time_step(0)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * As a mesh, this tutorial program implements two options, depending on the
 * global variable `testcase`: For the analytical variant (`testcase==0`),
 * the domain is $(0, 10) \times (-5, 5)$, with Dirichlet boundary
 * conditions (inflow) all around the domain. For `testcase==1`, we set the
 * domain to a cylinder in a rectangular box, derived from the flow past
 * cylinder testcase for incompressible viscous flow by Sch&auml;fer and
 * Turek (1996). Here, we have a larger variety of boundaries. The inflow
 * part at the left of the channel is given the inflow type, for which we
 * choose a constant inflow profile, whereas we set a subsonic outflow at
 * the right. For the boundary around the cylinder (boundary id equal to 2)
 * as well as the channel walls (boundary id equal to 3) we use the wall
 * boundary type, which is no-normal flow. Furthermore, for the 3D cylinder
 * we also add a gravity force in vertical direction. Having the base mesh
 * in place (including the manifolds set by
 * GridGenerator::channel_with_cylinder()), we can then perform the
 * specified number of global refinements, create the unknown numbering from
 * the DoFHandler, and hand the DoFHandler and Mapping objects to the
 * initialization of the EulerOperator.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::make_grid_and_dofs()
 *   {
 *     switch (testcase)
 *       {
 *         case 0:
 *           {
 *             Point<dim> lower_left;
 *             for (unsigned int d = 1; d < dim; ++d)
 *               lower_left[d] = -5;
 * 
 *             Point<dim> upper_right;
 *             upper_right[0] = 10;
 *             for (unsigned int d = 1; d < dim; ++d)
 *               upper_right[d] = 5;
 * 
 *             GridGenerator::hyper_rectangle(triangulation,
 *                                            lower_left,
 *                                            upper_right);
 *             triangulation.refine_global(2);
 * 
 *             euler_operator.set_inflow_boundary(
 *               0, std::make_unique<ExactSolution<dim>>(0));
 * 
 *             break;
 *           }
 * 
 *         case 1:
 *           {
 *             GridGenerator::channel_with_cylinder(
 *               triangulation, 0.03, 1, 0, true);
 * 
 *             euler_operator.set_inflow_boundary(
 *               0, std::make_unique<ExactSolution<dim>>(0));
 *             euler_operator.set_subsonic_outflow_boundary(
 *               1, std::make_unique<ExactSolution<dim>>(0));
 * 
 *             euler_operator.set_wall_boundary(2);
 *             euler_operator.set_wall_boundary(3);
 * 
 *             if (dim == 3)
 *               euler_operator.set_body_force(
 *                 std::make_unique<Functions::ConstantFunction<dim>>(
 *                   std::vector<double>({0., 0., -0.2})));
 * 
 *             break;
 *           }
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     triangulation.refine_global(n_global_refinements);
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     euler_operator.reinit(mapping, dof_handler);
 *     euler_operator.initialize_vector(solution);
 * 
 * @endcode
 * 
 * In the following, we output some statistics about the problem. Because we
 * often end up with quite large numbers of cells or degrees of freedom, we
 * would like to print them with a comma to separate each set of three
 * digits. This can be done via "locales", although the way this works is
 * not particularly intuitive. step-32 explains this in slightly more
 * detail.
 * 
 * @code
 *     std::locale s = pcout.get_stream().getloc();
 *     pcout.get_stream().imbue(std::locale(""));
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << " ( = " << (dim + 2) << " [vars] x "
 *           << triangulation.n_global_active_cells() << " [cells] x "
 *           << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )"
 *           << std::endl;
 *     pcout.get_stream().imbue(s);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For output, we first let the Euler operator compute the errors of the
 * numerical results. More precisely, we compute the error against the
 * analytical result for the analytical solution case, whereas we compute
 * the deviation against the background field with constant density and
 * energy and constant velocity in $x$ direction for the second test case.
 *   

 * 
 * The next step is to create output. This is similar to what is done in
 * step-33: We let the postprocessor defined above control most of the
 * output, except for the primal field that we write directly. For the
 * analytical solution test case, we also perform another projection of the
 * analytical solution and print the difference between that field and the
 * numerical solution. Once we have defined all quantities to be written, we
 * build the patches for output. Similarly to step-65, we create a
 * high-order VTK output by setting the appropriate flag, which enables us
 * to visualize fields of high polynomial degrees. Finally, we call the
 * `DataOutInterface::write_vtu_in_parallel()` function to write the result
 * to the given file name. This function uses special MPI parallel write
 * facilities, which are typically more optimized for parallel file systems
 * than the standard library's `std::ofstream` variants used in most other
 * tutorial programs. A particularly nice feature of the
 * `write_vtu_in_parallel()` function is the fact that it can combine output
 * from all MPI ranks into a single file, making it unnecessary to have a
 * central record of all such files (namely, the "pvtu" file).
 *   

 * 
 * For parallel programs, it is often instructive to look at the partitioning
 * of cells among processors. To this end, one can pass a vector of numbers
 * to DataOut::add_data_vector() that contains as many entries as the
 * current processor has active cells; these numbers should then be the
 * rank of the processor that owns each of these cells. Such a vector
 * could, for example, be obtained from
 * GridTools::get_subdomain_association(). On the other hand, on each MPI
 * process, DataOut will only read those entries that correspond to locally
 * owned cells, and these of course all have the same value: namely, the rank
 * of the current process. What is in the remaining entries of the vector
 * doesn't actually matter, and so we can just get away with a cheap trick: We
 * just fill *all* values of the vector we give to DataOut::add_data_vector()
 * with the rank of the current MPI process. The key is that on each process,
 * only the entries corresponding to the locally owned cells will be read,
 * ignoring the (wrong) values in other entries. The fact that every process
 * submits a vector in which the correct subset of entries is correct is all
 * that is necessary.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::output_results(const unsigned int result_number)
 *   {
 *     const std::array<double, 3> errors =
 *       euler_operator.compute_errors(ExactSolution<dim>(time), solution);
 *     const std::string quantity_name = testcase == 0 ? "error" : "norm";
 * 
 *     pcout << "Time:" << std::setw(8) << std::setprecision(3) << time
 *           << ", dt: " << std::setw(8) << std::setprecision(2) << time_step
 *           << ", " << quantity_name << " rho: " << std::setprecision(4)
 *           << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4)
 *           << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4)
 *           << std::setw(10) << errors[2] << std::endl;
 * 
 *     {
 *       TimerOutput::Scope t(timer, "output");
 * 
 *       Postprocessor postprocessor;
 *       DataOut<dim>  data_out;
 * 
 *       DataOutBase::VtkFlags flags;
 *       flags.write_higher_order_cells = true;
 *       data_out.set_flags(flags);
 * 
 *       data_out.attach_dof_handler(dof_handler);
 *       {
 *         std::vector<std::string> names;
 *         names.emplace_back("density");
 *         for (unsigned int d = 0; d < dim; ++d)
 *           names.emplace_back("momentum");
 *         names.emplace_back("energy");
 * 
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *           interpretation;
 *         interpretation.push_back(
 *           DataComponentInterpretation::component_is_scalar);
 *         for (unsigned int d = 0; d < dim; ++d)
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_part_of_vector);
 *         interpretation.push_back(
 *           DataComponentInterpretation::component_is_scalar);
 * 
 *         data_out.add_data_vector(dof_handler, solution, names, interpretation);
 *       }
 *       data_out.add_data_vector(solution, postprocessor);
 * 
 *       LinearAlgebra::distributed::Vector<Number> reference;
 *       if (testcase == 0 && dim == 2)
 *         {
 *           reference.reinit(solution);
 *           euler_operator.project(ExactSolution<dim>(time), reference);
 *           reference.sadd(-1., 1, solution);
 *           std::vector<std::string> names;
 *           names.emplace_back("error_density");
 *           for (unsigned int d = 0; d < dim; ++d)
 *             names.emplace_back("error_momentum");
 *           names.emplace_back("error_energy");
 * 
 *           std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *             interpretation;
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_scalar);
 *           for (unsigned int d = 0; d < dim; ++d)
 *             interpretation.push_back(
 *               DataComponentInterpretation::component_is_part_of_vector);
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_scalar);
 * 
 *           data_out.add_data_vector(dof_handler,
 *                                    reference,
 *                                    names,
 *                                    interpretation);
 *         }
 * 
 *       Vector<double> mpi_owner(triangulation.n_active_cells());
 *       mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
 *       data_out.add_data_vector(mpi_owner, "owner");
 * 
 *       data_out.build_patches(mapping,
 *                              fe.degree,
 *                              DataOut<dim>::curved_inner_cells);
 * 
 *       const std::string filename =
 *         "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu";
 *       data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The EulerProblem::run() function puts all pieces together. It starts off
 * by calling the function that creates the mesh and sets up data structures,
 * and then initializing the time integrator and the two temporary vectors of
 * the low-storage integrator. We call these vectors `rk_register_1` and
 * `rk_register_2`, and use the first vector to represent the quantity
 * $\mathbf{r}_i$ and the second one for $\mathbf{k}_i$ in the formulas for
 * the Runge--Kutta scheme outlined in the introduction. Before we start the
 * time loop, we compute the time step size by the
 * `EulerOperator::compute_cell_transport_speed()` function. For reasons of
 * comparison, we compare the result obtained there with the minimal mesh
 * size and print them to screen. For velocities and speeds of sound close
 * to unity as in this tutorial program, the predicted effective mesh size
 * will be close, but they could vary if scaling were different.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::run()
 *   {
 *     {
 *       const unsigned int n_vect_number = VectorizedArray<Number>::size();
 *       const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number;
 * 
 *       pcout << "Running with "
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
 *             << " MPI processes" << std::endl;
 *       pcout << "Vectorization over " << n_vect_number << " "
 *             << (std::is_same<Number, double>::value ? "doubles" : "floats")
 *             << " = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl;
 *     }
 * 
 *     make_grid_and_dofs();
 * 
 *     const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme);
 * 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_1;
 *     LinearAlgebra::distributed::Vector<Number> rk_register_2;
 *     rk_register_1.reinit(solution);
 *     rk_register_2.reinit(solution);
 * 
 *     euler_operator.project(ExactSolution<dim>(time), solution);
 * 
 *     double min_vertex_distance = std::numeric_limits<double>::max();
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         min_vertex_distance =
 *           std::min(min_vertex_distance, cell->minimum_vertex_distance());
 *     min_vertex_distance =
 *       Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);
 * 
 *     time_step = courant_number * integrator.n_stages() /
 *                 euler_operator.compute_cell_transport_speed(solution);
 *     pcout << "Time step size: " << time_step
 *           << ", minimal h: " << min_vertex_distance
 *           << ", initial transport scaling: "
 *           << 1. / euler_operator.compute_cell_transport_speed(solution)
 *           << std::endl
 *           << std::endl;
 * 
 *     output_results(0);
 * 
 * @endcode
 * 
 * Now we are ready to start the time loop, which we run until the time
 * has reached the desired end time. Every 5 time steps, we compute a new
 * estimate for the time step -- since the solution is nonlinear, it is
 * most effective to adapt the value during the course of the
 * simulation. In case the Courant number was chosen too aggressively, the
 * simulation will typically blow up with time step NaN, so that is easy
 * to detect here. One thing to note is that roundoff errors might
 * propagate to the leading digits due to an interaction of slightly
 * different time step selections that in turn lead to slightly different
 * solutions. To decrease this sensitivity, it is common practice to round
 * or truncate the time step size to a few digits, e.g. 3 in this case. In
 * case the current time is near the prescribed 'tick' value for output
 * (e.g. 0.02), we also write the output. After the end of the time loop,
 * we summarize the computation by printing some statistics, which is
 * mostly done by the TimerOutput::print_wall_time_statistics() function.
 * 
 * @code
 *     unsigned int timestep_number = 0;
 * 
 *     while (time < final_time - 1e-12)
 *       {
 *         ++timestep_number;
 *         if (timestep_number % 5 == 0)
 *           time_step =
 *             courant_number * integrator.n_stages() /
 *             Utilities::truncate_to_n_digits(
 *               euler_operator.compute_cell_transport_speed(solution), 3);
 * 
 *         {
 *           TimerOutput::Scope t(timer, "rk time stepping total");
 *           integrator.perform_time_step(euler_operator,
 *                                        time,
 *                                        time_step,
 *                                        solution,
 *                                        rk_register_1,
 *                                        rk_register_2);
 *         }
 * 
 *         time += time_step;
 * 
 *         if (static_cast<int>(time / output_tick) !=
 *               static_cast<int>((time - time_step) / output_tick) ||
 *             time >= final_time - 1e-12)
 *           output_results(
 *             static_cast<unsigned int>(std::round(time / output_tick)));
 *       }
 * 
 *     timer.print_wall_time_statistics(MPI_COMM_WORLD);
 *     pcout << std::endl;
 *   }
 * 
 * } // namespace Euler_DG
 * 
 * 
 * 
 * @endcode
 * 
 * The main() function is not surprising and follows what was done in all
 * previous MPI programs: As we run an MPI program, we need to call `MPI_Init()`
 * and `MPI_Finalize()`, which we do through the
 * Utilities::MPI::MPI_InitFinalize data structure. Note that we run the program
 * only with MPI, and set the thread count to 1.
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   using namespace Euler_DG;
 *   using namespace dealii;
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *   try
 *     {
 *       deallog.depth_console(0);
 * 
 *       EulerProblem<dimension> euler_problem;
 *       euler_problem.run();
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


<a name="Programoutput"></a><h3>Program output</h3>


Running the program with the default settings on a machine with 40 processes
produces the following output:
@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 147,456 ( = 4 [vars] x 1,024 [cells] x 36 [dofs/cell/var] )
Time step size: 0.00689325, minimal h: 0.3125, initial transport scaling: 0.102759

Time:       0, dt:   0.0069, error rho:   2.76e-07, rho * u:  1.259e-06, energy: 2.987e-06
Time:    1.01, dt:   0.0069, error rho:   1.37e-06, rho * u:  2.252e-06, energy: 4.153e-06
Time:    2.01, dt:   0.0069, error rho:  1.561e-06, rho * u:   2.43e-06, energy: 4.493e-06
Time:    3.01, dt:   0.0069, error rho:  1.714e-06, rho * u:  2.591e-06, energy: 4.762e-06
Time:    4.01, dt:   0.0069, error rho:  1.843e-06, rho * u:  2.625e-06, energy: 4.985e-06
Time:    5.01, dt:   0.0069, error rho:  1.496e-06, rho * u:  1.961e-06, energy: 4.142e-06
Time:       6, dt:   0.0083, error rho:  1.007e-06, rho * u:  7.119e-07, energy: 2.972e-06
Time:       7, dt:   0.0095, error rho:  9.096e-07, rho * u:  3.786e-07, energy: 2.626e-06
Time:       8, dt:   0.0096, error rho:  8.439e-07, rho * u:  3.338e-07, energy:  2.43e-06
Time:       9, dt:   0.0096, error rho:  7.822e-07, rho * u:  2.984e-07, energy: 2.248e-06
Time:      10, dt:   0.0096, error rho:  7.231e-07, rho * u:  2.666e-07, energy: 2.074e-06

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     2.249s    30 |     2.249s |     2.249s     8 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |  0.008066s    13 |   0.00952s |   0.01041s    20 |
| compute transport speed       |       258 |   0.01012s    13 |   0.05392s |   0.08574s    25 |
| output                        |        11 |    0.9597s    13 |    0.9613s |    0.9623s     6 |
| rk time stepping total        |      1283 |    0.9827s    25 |     1.015s |      1.06s    13 |
| rk_stage - integrals L_h      |      6415 |    0.8803s    26 |    0.9198s |    0.9619s    14 |
| rk_stage - inv mass + vec upd |      6415 |   0.05677s    15 |   0.06487s |   0.07597s    13 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

The program output shows that all errors are small. This is due to the fact
that we use a relatively fine mesh of $32^2$ cells with polynomials of degree
5 for a solution that is smooth. An interesting pattern shows for the time
step size: whereas it is 0.0069 up to time 5, it increases to 0.0096 for later
times. The step size increases once the vortex with some motion on top of the
speed of sound (and thus faster propagation) leaves the computational domain
between times 5 and 6.5. After that point, the flow is simply uniform
in the same direction, and the maximum velocity of the gas is reduced
compared to the previous state where the uniform velocity was overlaid
by the vortex. Our time step formula recognizes this effect.

The final block of output shows detailed information about the timing
of individual parts of the programs; it breaks this down by showing
the time taken by the fastest and the slowest processor, and the
average time -- this is often useful in very large computations to
find whether there are processors that are consistently overheated
(and consequently are throttling their clock speed) or consistently
slow for other reasons.
The summary shows that 1283 time steps have been performed
in 1.02 seconds (looking at the average time among all MPI processes), while
the output of 11 files has taken additional 0.96 seconds. Broken down per time
step and into the five Runge--Kutta stages, the compute time per evaluation is
0.16 milliseconds. This high performance is typical of matrix-free evaluators
and a reason why explicit time integration is very competitive against
implicit solvers, especially for large-scale simulations. The breakdown of
computational times at the end of the program run shows that the evaluation of
integrals in $\mathcal L_h$ contributes with around 0.92 seconds and the
application of the inverse mass matrix with 0.06 seconds. Furthermore, the
estimation of the transport speed for the time step size computation
contributes with another 0.05 seconds of compute time.

If we use three more levels of global refinement and 9.4 million DoFs in total,
the final statistics are as follows (for the modified Lax--Friedrichs flux,
$p=5$, and the same system of 40 cores of dual-socket Intel Xeon Gold 6230):
@code
+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     244.9s    12 |     244.9s |     244.9s    34 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |    0.4239s    12 |    0.4318s |    0.4408s     9 |
| compute transport speed       |      2053 |     3.962s    12 |     6.727s |     10.12s     7 |
| output                        |        11 |     30.35s    12 |     30.36s |     30.37s     9 |
| rk time stepping total        |     10258 |     201.7s     7 |     205.1s |     207.8s    12 |
| rk_stage - integrals L_h      |     51290 |     121.3s     6 |     126.6s |     136.3s    16 |
| rk_stage - inv mass + vec upd |     51290 |     66.19s    16 |     77.52s |     81.84s    10 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

Per time step, the solver now takes 0.02 seconds, about 25 times as long as
for the small problem with 147k unknowns. Given that the problem involves 64
times as many unknowns, the increase in computing time is not
surprising. Since we also do 8 times as many time steps, the compute time
should in theory increase by a factor of 512. The actual increase is 205 s /
1.02 s = 202. This is because the small problem size cannot fully utilize the
40 cores due to communication overhead. This becomes clear if we look into the
details of the operations done per time step. The evaluation of the
differential operator $\mathcal L_h$ with nearest neighbor communication goes
from 0.92 seconds to 127 seconds, i.e., it increases with a factor of 138. On
the other hand, the cost for application of the inverse mass matrix and the
vector updates, which do not need to communicate between the MPI processes at
all, has increased by a factor of 1195. The increase is more than the
theoretical factor of 512 because the operation is limited by the bandwidth
from RAM memory for the larger size while for the smaller size, all vectors
fit into the caches of the CPU. The numbers show that the mass matrix
evaluation and vector update part consume almost 40% of the time spent by the
Runge--Kutta stages -- despite using a low-storage Runge--Kutta integrator and
merging of vector operations! And despite using over-integration for the
$\mathcal L_h$ operator. For simpler differential operators and more expensive
time integrators, the proportion spent in the mass matrix and vector update
part can also reach 70%. If we compute a throughput number in terms of DoFs
processed per second and Runge--Kutta stage, we obtain @f[ \text{throughput} =
\frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} = \frac{10258 \cdot 5 \cdot
9.4\,\text{MDoFs}}{205s} = 2360\, \text{MDoFs/s} @f] This throughput number is
very high, given that simply copying one vector to another one runs at
only around 10,000 MDoFs/s.

If we go to the next-larger size with 37.7 million DoFs, the overall
simulation time is 2196 seconds, with 1978 seconds spent in the time
stepping. The increase in run time is a factor of 9.3 for the L_h operator
(1179 versus 127 seconds) and a factor of 10.3 for the inverse mass matrix and
vector updates (797 vs 77.5 seconds). The reason for this non-optimal increase
in run time can be traced back to cache effects on the given hardware (with 40
MB of L2 cache and 55 MB of L3 cache): While not all of the relevant data fits
into caches for 9.4 million DoFs (one vector takes 75 MB and we have three
vectors plus some additional data in MatrixFree), there is capacity for one and
a half vector nonetheless. Given that modern caches are more sophisticated than
the naive least-recently-used strategy (where we would have little re-use as
the data is used in a streaming-like fashion), we can assume that a sizeable
fraction of data can indeed be delivered from caches for the 9.4 million DoFs
case. For the larger case, even with optimal caching less than 10 percent of
data would fit into caches, with an associated loss in performance.


<a name="Convergenceratesfortheanalyticaltestcase"></a><h3>Convergence rates for the analytical test case</h3>


For the modified Lax--Friedrichs flux and measuring the error in the momentum
variable, we obtain the following convergence table (the rates are very
similar for the density and energy variables):

<table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.373e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.130e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.899e-03</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">5.577e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">7.381e-03</td>
    <td>3.64</td>
    <td align="right">36,864</td>
    <td align="center">2.082e-04</td>
    <td>5.42</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.724e-03</td>
    <td>3.56</td>
    <td align="right">65,536</td>
    <td align="center">3.072e-04</td>
    <td>4.59</td>
    <td align="right">147,456</td>
    <td align="center">2.625e-06</td>
    <td>6.31</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">6.205e-04</td>
    <td>2.92</td>
    <td align="right">262,144</td>
    <td align="center">1.880e-05</td>
    <td>4.03</td>
    <td align="right">589,824</td>
    <td align="center">3.268e-08</td>
    <td>6.33</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">8.279e-05</td>
    <td>2.91</td>
    <td align="right">1,048,576</td>
    <td align="center">1.224e-06</td>
    <td>3.94</td>
    <td align="right">2,359,296</td>
    <td align="center">9.252e-10</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.105e-05</td>
    <td>2.91</td>
    <td align="right">4,194,304</td>
    <td align="center">7.871e-08</td>
    <td>3.96</td>
    <td align="right">9,437,184</td>
    <td align="center">1.369e-10</td>
    <td>2.77</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">1.615e-06</td>
    <td>2.77</td>
    <td align="right">16,777,216</td>
    <td align="center">4.961e-09</td>
    <td>3.99</td>
    <td align="right">37,748,736</td>
    <td align="center">7.091e-11</td>
    <td>0.95</td>
  </tr>
</table>

If we switch to the Harten-Lax-van Leer flux, the results are as follows:
<table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.339e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.037e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.849e-03</td>
    <td>3.92</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">4.204e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">9.143e-03</td>
    <td>3.31</td>
    <td align="right">36,864</td>
    <td align="center">2.501e-04</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.913e-03</td>
    <td>3.09</td>
    <td align="right">65,536</td>
    <td align="center">3.257e-04</td>
    <td>4.81</td>
    <td align="right">147,456</td>
    <td align="center">3.260e-06</td>
    <td>6.26</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">7.862e-04</td>
    <td>2.64</td>
    <td align="right">262,144</td>
    <td align="center">1.588e-05</td>
    <td>4.36</td>
    <td align="right">589,824</td>
    <td align="center">2.953e-08</td>
    <td>6.79</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">1.137e-04</td>
    <td>2.79</td>
    <td align="right">1,048,576</td>
    <td align="center">9.400e-07</td>
    <td>4.08</td>
    <td align="right">2,359,296</td>
    <td align="center">4.286e-10</td>
    <td>6.11</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.476e-05</td>
    <td>2.95</td>
    <td align="right">4,194,304</td>
    <td align="center">5.799e-08</td>
    <td>4.02</td>
    <td align="right">9,437,184</td>
    <td align="center">2.789e-11</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">2.038e-06</td>
    <td>2.86</td>
    <td align="right">16,777,216</td>
    <td align="center">3.609e-09</td>
    <td>4.01</td>
    <td align="right">37,748,736</td>
    <td align="center">5.730e-11</td>
    <td>-1.04</td>
  </tr>
</table>

The tables show that we get optimal $\mathcal O\left(h^{p+1}\right)$
convergence rates for both numerical fluxes. The errors are slightly smaller
for the Lax--Friedrichs flux for $p=2$, but the picture is reversed for
$p=3$; in any case, the differences on this testcase are relatively
small.

For $p=5$, we reach the roundoff accuracy of $10^{-11}$ with both
fluxes on the finest grids. Also note that the errors are absolute with a
domain length of $10^2$, so relative errors are below $10^{-12}$. The HLL flux
is somewhat better for the highest degree, which is due to a slight inaccuracy
of the Lax--Friedrichs flux: The Lax--Friedrichs flux sets a Dirichlet
condition on the solution that leaves the domain, which results in a small
artificial reflection, which is accentuated for the Lax--Friedrichs
flux. Apart from that, we see that the influence of the numerical flux is
minor, as the polynomial part inside elements is the main driver of the
accucary. The limited influence of the flux also has consequences when trying
to approach more challenging setups with the higher-order DG setup: Taking for
example the parameters and grid of step-33, we get oscillations (which in turn
make density negative and make the solution explode) with both fluxes once the
high-mass part comes near the boundary, as opposed to the low-order finite
volume case ($p=0$). Thus, any case that leads to shocks in the solution
necessitates some form of limiting or artificial dissipation. For another
alternative, see the step-69 tutorial program.


<a name="Resultsforflowinchannelaroundcylinderin2D"></a><h3>Results for flow in channel around cylinder in 2D</h3>


For the test case of the flow around a cylinder in a channel, we need to
change the first code line to
@code
  constexpr unsigned int testcase = 1;
@endcode
This test case starts with a background field of a constant velocity
of Mach number 0.31 and a constant initial density; the flow will have
to go around an obstacle in the form of a cylinder. Since we impose a
no-penetration condition on the cylinder walls, the flow that
initially impinges head-on onto to cylinder has to rearrange,
which creates a big sound wave. The following pictures show the pressure at
times 0.1, 0.25, 0.5, and 1.0 (top left to bottom right) for the 2D case with
5 levels of global refinement, using 102,400 cells with polynomial degree of
5 and 14.7 million degrees of freedom over all 4 solution variables.
We clearly see the discontinuity that
propagates slowly in the upstream direction and more quickly in downstream
direction in the first snapshot at time 0.1. At time 0.25, the sound wave has
reached the top and bottom walls and reflected back to the interior. From the
different distances of the reflected waves from lower and upper walls we can
see the slight asymmetry of the Sch&auml;fer-Turek test case represented by
GridGenerator::channel_with_cylinder() with somewhat more space above the
cylinder compared to below. At later times, the picture is more chaotic with
many sound waves all over the place.

<table align="center" class="doxtable" style="width:85%">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_010.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_025.png" alt="" width="100%">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_050.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_100.png" alt="" width="100%">
    </td>
  </tr>
</table>

The next picture shows an elevation plot of the pressure at time 1.0 looking
from the channel inlet towards the outlet at the same resolution -- here,
we can see the large number
of reflections. In the figure, two types of waves are visible. The
larger-amplitude waves correspond to various reflections that happened as the
initial discontinuity hit the walls, whereas the small-amplitude waves of
size similar to the elements correspond to numerical artifacts. They have their
origin in the finite resolution of the scheme and appear as the discontinuity
travels through elements with high-order polynomials. This effect can be cured
by increasing resolution. Apart from this effect, the rich wave structure is
the result of the transport accuracy of the high-order DG method.

<img src="https://www.dealii.org/images/steps/developer/step-67.pressure_elevated.jpg" alt="" width="40%">

With 2 levels of global refinement with 1,600 cells, the mesh and its
partitioning on 40 MPI processes looks as follows:

<img src="https://www.dealii.org/images/steps/developer/step-67.grid-owner.png" alt="" width="70%">

When we run the code with 4 levels of global refinements on 40 cores, we get
the following output:
@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 3,686,400 ( = 4 [vars] x 25,600 [cells] x 36 [dofs/cell/var] )
Time step size: 7.39876e-05, minimal h: 0.001875, initial transport scaling: 0.00110294

Time:       0, dt:  7.4e-05, norm rho:   4.17e-16, rho * u:  1.629e-16, energy: 1.381e-15
Time:    0.05, dt:  6.3e-05, norm rho:    0.02075, rho * u:    0.03801, energy:   0.08772
Time:     0.1, dt:  5.9e-05, norm rho:    0.02211, rho * u:    0.04515, energy:   0.08953
Time:    0.15, dt:  5.7e-05, norm rho:    0.02261, rho * u:    0.04592, energy:   0.08967
Time:     0.2, dt:  5.8e-05, norm rho:    0.02058, rho * u:    0.04361, energy:   0.08222
Time:    0.25, dt:  5.9e-05, norm rho:    0.01695, rho * u:    0.04203, energy:   0.06873
Time:     0.3, dt:  5.9e-05, norm rho:    0.01653, rho * u:     0.0401, energy:   0.06604
Time:    0.35, dt:  5.7e-05, norm rho:    0.01774, rho * u:    0.04264, energy:    0.0706

...

Time:    1.95, dt:  5.8e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05185
Time:       2, dt:  5.7e-05, norm rho:    0.01432, rho * u:    0.03969, energy:   0.04889

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     273.6s    13 |     273.6s |     273.6s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.01112s    35 |    0.0672s |    0.1337s     0 |
| compute transport speed       |      6914 |     5.422s    35 |     15.96s |     29.99s     1 |
| output                        |        41 |     37.24s    35 |      37.3s |     37.37s     0 |
| rk time stepping total        |     34564 |     205.4s     1 |     219.5s |     230.1s    35 |
| rk_stage - integrals L_h      |    172820 |     153.6s     1 |     164.9s |     175.6s    27 |
| rk_stage - inv mass + vec upd |    172820 |     47.13s    13 |     53.09s |     64.05s    33 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

The norms shown here for the various quantities are the deviations
$\rho'$, $(\rho u)'$, and $E'$ against the background field (namely, the
initial condition). The distribution of run time is overall similar as in the
previous test case. The only slight difference is the larger proportion of
time spent in $\mathcal L_h$ as compared to the inverse mass matrix and vector
updates. This is because the geometry is deformed and the matrix-free
framework needs to load additional arrays for the geometry from memory that
are compressed in the affine mesh case.

Increasing the number of global refinements to 5, the output becomes:
@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 14,745,600 ( = 4 [vars] x 102,400 [cells] x 36 [dofs/cell/var] )

...

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |      2693s    32 |      2693s |      2693s    23 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.04537s    32 |     0.173s |    0.3489s     0 |
| compute transport speed       |     13858 |     40.75s    32 |     85.99s |     149.8s     0 |
| output                        |        41 |     153.8s    32 |     153.9s |     154.1s     0 |
| rk time stepping total        |     69284 |      2386s     0 |      2450s |      2496s    32 |
| rk_stage - integrals L_h      |    346420 |      1365s    32 |      1574s |      1718s    19 |
| rk_stage - inv mass + vec upd |    346420 |     722.5s    10 |     870.7s |      1125s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

The effect on performance is similar to the analytical test case -- in
theory, computation times should increase by a factor of 8, but we actually
see an increase by a factor of 11 for the time steps (219.5 seconds versus
2450 seconds). This can be traced back to caches, with the small case mostly
fitting in caches. An interesting effect, typical of programs with a mix of
local communication (integrals $\mathcal L_h$) and global communication (computation of
transport speed) with some load imbalance, can be observed by looking at the
MPI ranks that encounter the minimal and maximal time of different phases,
respectively. Rank 0 reports the fastest throughput for the "rk time stepping
total" part. At the same time, it appears to be slowest for the "compute
transport speed" part, almost a factor of 2 slower than the
average and almost a factor of 4 compared to the faster rank.
Since the latter involves global communication, we can attribute the
slowness in this part to the fact that the local Runge--Kutta stages have
advanced more quickly on this rank and need to wait until the other processors
catch up. At this point, one can wonder about the reason for this imbalance:
The number of cells is almost the same on all MPI processes.
However, the matrix-free framework is faster on affine and Cartesian
cells located towards the outlet of the channel, to which the lower MPI ranks
are assigned. On the other hand, rank 32, which reports the highest run time
for the Runga--Kutta stages, owns the curved cells near the cylinder, for
which no data compression is possible. To improve throughput, we could assign
different weights to different cell types when partitioning the
parallel::distributed::Triangulation object, or even measure the run time for a
few time steps and try to rebalance then.

The throughput per Runge--Kutta stage can be computed to 2085 MDoFs/s for the
14.7 million DoFs test case over the 346,000 Runge--Kutta stages, slightly slower
than the Cartesian mesh throughput of 2360 MDoFs/s reported above.

Finally, if we add one additional refinement, we record the following output:
@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 58,982,400 ( = 4 [vars] x 409,600 [cells] x 36 [dofs/cell/var] )

...

Time:    1.95, dt:  1.4e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05183
Time:       2, dt:  1.4e-05, norm rho:    0.01431, rho * u:    0.03969, energy:   0.04887

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 2.166e+04s    26 | 2.166e+04s | 2.166e+04s    24 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.1758s    30 |     0.672s |     1.376s     1 |
| compute transport speed       |     27748 |     321.3s    34 |     678.8s |      1202s     1 |
| output                        |        41 |     616.3s    32 |     616.4s |     616.4s    34 |
| rk time stepping total        |    138733 | 1.983e+04s     1 | 2.036e+04s | 2.072e+04s    34 |
| rk_stage - integrals L_h      |    693665 | 1.052e+04s    32 | 1.248e+04s | 1.387e+04s    19 |
| rk_stage - inv mass + vec upd |    693665 |      6404s    10 |      7868s | 1.018e+04s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

The "rk time stepping total" part corresponds to a throughput of 2010 MDoFs/s. The
overall run time to perform 139k time steps is 20k seconds (5.7 hours) or 7
time steps per second -- not so bad for having nearly 60 million
unknowns. More throughput can be achieved by adding more cores to
the computation.


<a name="Resultsforflowinchannelaroundcylinderin3D"></a><h3>Results for flow in channel around cylinder in 3D</h3>


Switching the channel test case to 3D with 3 global refinements, the output is
@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 221,184,000 ( = 5 [vars] x 204,800 [cells] x 216 [dofs/cell/var] )

...

Time:    1.95, dt:  0.00011, norm rho:    0.01131, rho * u:    0.03056, energy:   0.04091
Time:       2, dt:  0.00011, norm rho:     0.0119, rho * u:    0.03142, energy:   0.04425

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 1.734e+04s     4 | 1.734e+04s | 1.734e+04s    38 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.6551s    34 |     3.216s |     7.281s     0 |
| compute transport speed       |      3546 |       160s    34 |     393.2s |     776.9s     0 |
| output                        |        41 |      1350s    34 |      1353s |      1357s     0 |
| rk time stepping total        |     17723 | 1.519e+04s     0 | 1.558e+04s | 1.582e+04s    34 |
| rk_stage - integrals L_h      |     88615 | 1.005e+04s    32 | 1.126e+04s |  1.23e+04s    11 |
| rk_stage - inv mass + vec upd |     88615 |      3056s    11 |      4322s |      5759s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

The physics are similar to the 2D case, with a slight motion in the z
direction due to the gravitational force. The throughput per Runge--Kutta
stage in this case is
@f[
\text{throughput} = \frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} =
\frac{17723 \cdot 5 \cdot 221.2\,\text{M}}{15580s} = 1258\, \text{MDoFs/s}.
@f]

The throughput is lower than in 2D because the computation of the $\mathcal L_h$ term
is more expensive. This is due to over-integration with `degree+2` points and
the larger fraction of face integrals (worse volume-to-surface ratio) with
more expensive flux computations. If we only consider the inverse mass matrix
and vector update part, we record a throughput of 4857 MDoFs/s for the 2D case
of the isentropic vortex with 37.7 million unknowns, whereas the 3D case
runs with 4535 MDoFs/s. The performance is similar because both cases are in
fact limited by the memory bandwidth.

If we go to four levels of global refinement, we need to increase the number
of processes to fit everything in memory -- the computation needs around 350
GB of RAM memory in this case. Also, the time it takes to complete 35k time
steps becomes more tolerable by adding additional resources. We therefore use
6 nodes with 40 cores each, resulting in a computation with 240 MPI processes:
@code
Running with 240 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 1,769,472,000 ( = 5 [vars] x 1,638,400 [cells] x 216 [dofs/cell/var] )

...

Time:    1.95, dt:  5.6e-05, norm rho:    0.01129, rho * u:     0.0306, energy:   0.04086
Time:       2, dt:  5.6e-05, norm rho:    0.01189, rho * u:    0.03145, energy:   0.04417

+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 5.396e+04s   151 | 5.396e+04s | 5.396e+04s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |     2.632s   178 |     7.221s |     16.56s     0 |
| compute transport speed       |      7072 |       714s   193 |      1553s |      3351s     0 |
| output                        |        41 |      8065s   176 |      8070s |      8079s     0 |
| rk time stepping total        |     35350 |  4.25e+04s     0 |  4.43e+04s | 4.515e+04s   193 |
| rk_stage - integrals L_h      |    176750 | 2.936e+04s   134 | 3.222e+04s |  3.67e+04s    99 |
| rk_stage - inv mass + vec upd |    176750 |      7004s    99 | 1.207e+04s |  1.55e+04s   132 |
+-------------------------------------------+------------------+------------+------------------+
@endcode
This simulation had nearly 2 billion unknowns -- quite a large
computation indeed, and still only needed around 1.5 seconds per time
step.


<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


The code presented here straight-forwardly extends to adaptive meshes, given
appropriate indicators for setting the refinement flags. Large-scale
adaptivity of a similar solver in the context of the acoustic wave equation
has been achieved by the <a href="https://github.com/kronbichler/exwave">exwave
project</a>. However, in the present context, the benefits of adaptivity are often
limited to early times and effects close to the origin of sound waves, as the
waves eventually reflect and diffract. This leads to steep gradients all over
the place, similar to turbulent flow, and a more or less globally
refined mesh.

Another topic that we did not discuss in the results section is a comparison
of different time integration schemes. The program provides four variants of
low-storage Runga--Kutta integrators that each have slightly different
accuracy and stability behavior. Among the schemes implemented here, the
higher-order ones provide additional accuracy but come with slightly lower
efficiency in terms of step size per stage before they violate the CFL
condition. An interesting extension would be to compare the low-storage
variants proposed here with standard Runge--Kutta integrators or to use vector
operations that are run separate from the mass matrix operation and compare
performance.


<a name="Moreadvancednumericalfluxfunctionsandskewsymmetricformulations"></a><h4>More advanced numerical flux functions and skew-symmetric formulations</h4>


As mentioned in the introduction, the modified Lax--Friedrichs flux and the
HLL flux employed in this program are only two variants of a large body of
numerical fluxes available in the literature on the Euler equations. One
example is the HLLC flux (Harten-Lax-van Leer-Contact) flux which adds the
effect of rarefaction waves missing in the HLL flux, or the Roe flux. As
mentioned in the introduction, the effect of numerical fluxes on high-order DG
schemes is debatable (unlike for the case of low-order discretizations).

A related improvement to increase the stability of the solver is to also
consider the spatial integral terms. A shortcoming in the rather naive
implementation used above is the fact that the energy conservation of the
original Euler equations (in the absence of shocks) only holds up to a
discretization error. If the solution is under-resolved, the discretization
error can give rise to an increase in the numerical energy and eventually
render the discretization unstable. This is because of the inexact numerical
integration of the terms in the Euler equations, which both contain rational
nonlinearities and higher-degree content from curved cells. A way out of this
dilemma are so-called skew-symmetric formulations, see @cite Gassner2013 for a
simple variant. Skew symmetry means that switching the role of the solution
$\mathbf{w}$ and test functions $\mathbf{v}$ in the weak form produces the
exact negative of the original quantity, apart from some boundary terms. In
the discrete setting, the challenge is to keep this skew symmetry also when
the integrals are only computed approximately (in the continuous case,
skew-symmetry is a consequence of integration by parts). Skew-symmetric
numerical schemes balance spatial derivatives in the conservative form
$(\nabla \mathbf v, \mathbf{F}(\mathbf w))_{K}$ with contributions in the
convective form $(\mathbf v, \tilde{\mathbf{F}}(\mathbf w)\nabla
\mathbf{w})_{K}$ for some $\tilde{\mathbf{F}}$. The precise terms depend on
the equation and the integration formula, and can in some cases by understood
by special skew-symmetric finite difference schemes.

To get started, interested readers could take a look at
https://github.com/kronbichler/advection_miniapp, where a
skew-symmetric DG formulation is implemented with deal.II for a simple advection
equation.

<a name="Equippingthecodeforsupersoniccalculations"></a><h4>Equipping the code for supersonic calculations</h4>


As mentioned in the introduction, the solution to the Euler equations develops
shocks as the Mach number increases, which require additional mechanisms to
stabilize the scheme, e.g. in the form of limiters. The main challenge besides
actually implementing the limiter or artificial viscosity approach would be to
load-balance the computations, as the additional computations involved for
limiting the oscillations in troubled cells would make them more expensive than the
plain DG cells without limiting. Furthermore, additional numerical fluxes that
better cope with the discontinuities would also be an option.

One ingredient also necessary for supersonic flows are appropriate boundary
conditions. As opposed to the subsonic outflow boundaries discussed in the
introduction and implemented in the program, all characteristics are outgoing
for supersonic outflow boundaries, so we do not want to prescribe any external
data,
@f[
\mathbf{w}^+ = \mathbf{w}^- = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E^-\end{pmatrix} \quad
 \text{(Neumann)}.
@f]

In the code, we would simply add the additional statement
@code
            else if (supersonic_outflow_boundaries.find(boundary_id) !=
                     supersonic_outflow_boundaries.end())
              {
                w_p        = w_m;
                at_outflow = true;
              }
@endcode
in the `local_apply_boundary_face()` function.

<a name="ExtensiontothelinearizedEulerequations"></a><h4>Extension to the linearized Euler equations</h4>


When the interest with an Euler solution is mostly in the propagation of sound
waves, it often makes sense to linearize the Euler equations around a
background state, i.e., a given density, velocity and energy (or pressure)
field, and only compute the change against these fields. This is the setting
of the wide field of aeroacoustics. Even though the resolution requirements
are sometimes considerably reduced, implementation gets somewhat more
complicated as the linearization gives rise to additional terms. From a code
perspective, in the operator evaluation we also need to equip the code with
the state to linearize against. This information can be provided either by
analytical functions (that are evaluated in terms of the position of the
quadrature points) or by a vector similar to the solution. Based on that
vector, we would create an additional FEEvaluation object to read from it and
provide the values of the field at quadrature points. If the background
velocity is zero and the density is constant, the linearized Euler equations
further simplify and can equivalently be written in the form of the
acoustic wave equation.

A challenge in the context of sound propagation is often the definition of
boundary conditions, as the computational domain needs to be of finite size,
whereas the actual simulation often spans an infinite (or at least much
larger) physical domain. Conventional Dirichlet or Neumann boundary conditions
give rise to reflections of the sound waves that eventually propagate back to
the region of interest and spoil the solution. Therefore, various variants of
non-reflecting boundary conditions or sponge layers, often in the form of
<a
href="https://en.wikipedia.org/wiki/Perfectly_matched_layer">perfectly
matched layers</a> -- where the solution is damped without reflection
-- are common.


<a name="ExtensiontothecompressibleNavierStokesequations"></a><h4>Extension to the compressible Navier-Stokes equations</h4>


The solver presented in this tutorial program can also be extended to the
compressible Navier--Stokes equations by adding viscous terms, as described in
@cite FehnWallKronbichler2019. To keep as much of the performance obtained
here despite the additional cost of elliptic terms, e.g. via an interior
penalty method, one can switch the basis from FE_DGQ to FE_DGQHermite like in
the step-59 tutorial program.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-67.cc"
*/
