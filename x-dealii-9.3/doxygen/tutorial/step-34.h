/**
@page step_34 The step-34 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Irrotationalflow"> Irrotational flow </a>
        <li><a href="#Thenumericalapproximation">The numerical approximation</a>
        <li><a href="#Collocationboundaryelementmethod"> Collocation boundary element method </a>
        <li><a href="#Treatingthesingularintegrals"> Treating the singular integrals. </a>
        <li><a href="#Implementation">Implementation</a>
        <li><a href="#Testcase">Testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Singleanddoublelayeroperatorkernels">Single and double layer operator kernels</a>
        <li><a href="#TheBEMProblemclass">The BEMProblem class</a>
      <ul>
        <li><a href="#BEMProblemBEMProblemandBEMProblemread_parameters">BEMProblem::BEMProblem and BEMProblem::read_parameters</a>
        <li><a href="#BEMProblemread_domain">BEMProblem::read_domain</a>
        <li><a href="#BEMProblemrefine_and_resize">BEMProblem::refine_and_resize</a>
        <li><a href="#BEMProblemassemble_system">BEMProblem::assemble_system</a>
        <li><a href="#BEMProblemsolve_system">BEMProblem::solve_system</a>
        <li><a href="#BEMProblemcompute_errors">BEMProblem::compute_errors</a>
        <li><a href="#BEMProblemcompute_exterior_solution">BEMProblem::compute_exterior_solution</a>
        <li><a href="#BEMProblemoutput_results">BEMProblem::output_results</a>
        <li><a href="#BEMProblemrun">BEMProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
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

<i>This program was contributed by Luca Heltai (thanks to Michael
Gratton for pointing out what the exact solution should have been in
the three dimensional case).  </i>

@dealiiTutorialDOI{10.5281/zenodo.495473,https://zenodo.org/badge/DOI/10.5281/zenodo.495473.svg}

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Irrotationalflow"></a><h3> Irrotational flow </h3>

The incompressible motion of an inviscid fluid past a body (for
example air past an airplane wing, or air or water past a propeller) is
usually modeled by the Euler equations of fluid dynamics:

\f{align*}
  \frac{\partial }{\partial t}\mathbf{v} + (\mathbf{v}\cdot\nabla)\mathbf{v}
  &=
  -\frac{1}{\rho}\nabla p + \mathbf{g}
  \qquad &\text{in } \mathbb{R}^n \backslash \Omega
  \\
  \nabla \cdot \mathbf{v}&=0
  &\text{in } \mathbb{R}^n\backslash\Omega
\f}
where the fluid density $\rho$ and the acceleration $\mathbf{g}$ due
to external forces are given and the velocity $\mathbf{v}$ and the
pressure $p$ are the unknowns. Here $\Omega$ is a closed bounded
region representing the body around which the fluid moves.

The above equations can be derived from Navier-Stokes equations
assuming that the effects due to viscosity are negligible compared to
those due to the pressure gradient, inertial forces and the external
forces. This is the opposite case of the Stokes equations discussed in
step-22 which are the limit case of dominant viscosity,
i.e. where the velocity is so small that inertia forces can be
neglected. On the other hand, owing to the assumed incompressibility,
the equations are not suited for very high speed gas flows where
compressibility and the equation of state of the gas have to be taken
into account, leading to the Euler equations of gas dynamics, a
hyperbolic system.

For the purpose of this tutorial program, we will consider only stationary
flow without external forces:
\f{align*}
  (\mathbf{v}\cdot\nabla)\mathbf{v}
  &=
  -\frac{1}{\rho}\nabla p
  \qquad &\text{in } \mathbb{R}^n \backslash \Omega
  \\
  \nabla \cdot \mathbf{v}&=0
  &\text{in } \mathbb{R}^n\backslash\Omega
\f}


Uniqueness of the solution of the Euler equations is ensured by adding the
boundary conditions
\f[
  \label{eq:boundary-conditions}
  \begin{aligned}
    \mathbf{n}\cdot\mathbf{v}& = 0 \qquad && \text{ on } \partial\Omega \\
    \mathbf{v}& = \mathbf{v}_\infty && \text{ when } |\mathbf{x}| \to \infty,
  \end{aligned}
\f]

which is to say that the body is at rest in our coordinate systems and
is not permeable, and that the fluid has (constant) velocity
$\mathbf{v}_\infty$ at infinity. An alternative viewpoint is that our
coordinate system moves along with the body whereas the background
fluid is at rest at infinity. Notice that we define the normal
$\mathbf{n}$ as the <i>outer</i> normal to the domain $\Omega$, which
is the opposite of the outer normal to the integration domain.

For both stationary and non stationary flow, the solution process
starts by solving for the velocity in the second equation and
substituting in the first equation in order to find the pressure.
The solution of the stationary Euler equations is typically performed
in order to understand the behavior of the given (possibly complex)
geometry when a prescribed motion is enforced on the system.

The first step in this process is to change the frame of reference from a
coordinate system moving along with the body to one in which the body moves
through a fluid that is at rest at infinity. This can be expressed by
introducing a new velocity $\mathbf{\tilde{v}}=\mathbf{v}-\mathbf{v}_\infty$ for
which we find that the same equations hold (because $\nabla\cdot
\mathbf{v}_\infty=0$) and we have boundary conditions
\f[
  \label{eq:boundary-conditions-tilde}
  \begin{aligned}
    \mathbf{n}\cdot\mathbf{\tilde{v}}& = -\mathbf{n}\cdot\mathbf{v}_\infty \qquad && \text{ on } \partial\Omega \\
    \mathbf{\tilde{v}}& = 0 && \text{ when } |\mathbf{x}| \to \infty,
  \end{aligned}
\f]

If we assume that the fluid is irrotational, i.e., $\nabla \times
\mathbf{v}=0$ in $\mathbb{R}^n\backslash\Omega$, we can represent the
velocity, and consequently also the perturbation velocity, as the
gradient of a scalar function:
\f[
  \mathbf{\tilde{v}}=\nabla\phi,
\f]
and so the second part of Euler equations above can be rewritten
as the homogeneous Laplace equation for the unknown $\phi$:
\f{align*}
\label{laplace}
\Delta\phi &= 0 \qquad &&\text{in}\ \mathbb{R}^n\backslash\Omega,
	   \\
	   \mathbf{n}\cdot\nabla\phi &= -\mathbf{n}\cdot\mathbf{v}_\infty
	   && \text{on}\ \partial\Omega
\f}
while the momentum equation reduces to Bernoulli's equation that expresses the
pressure $p$ as a function of the potential $\phi$:
\f[
\frac{p}{\rho} +\frac{1}{2} | \nabla \phi |^2 = 0 \in \Omega.
\f]

So we can solve the problem by solving the Laplace equation for the
potential.  We recall that the following functions, called fundamental
solutions of the Laplace equation,

\f[ \begin{aligned}
\label{eq:3} G(\mathbf{y}-\mathbf{x}) = &
-\frac{1}{2\pi}\ln|\mathbf{y}-\mathbf{x}| \qquad && \text{for } n=2 \\
G(\mathbf{y}-\mathbf{x}) = &
\frac{1}{4\pi}\frac{1}{|\mathbf{y}-\mathbf{x}|}&& \text{for } n=3,
\end{aligned}
\f]

satisfy in a distributional sense the equation:

\f[
-\Delta_y G(\mathbf{y}-\mathbf{x}) = \delta(\mathbf{y}-\mathbf{x}),
\f]

where the derivative is done in the variable $\mathbf{y}$. By using
the usual Green identities, our problem can be written on the boundary
$\partial\Omega = \Gamma$ only. We recall the general definition of
the second Green %identity:

\f[\label{green}
  \int_{\omega}
  (-\Delta u)v\,dx + \int_{\partial\omega} \frac{\partial u}{\partial \tilde{\mathbf{n}} }v \,ds
  =
  \int_{\omega}
  (-\Delta v)u\,dx + \int_{\partial\omega} u\frac{\partial v}{\partial \tilde{\mathbf{n}}} \,ds,
\f]

where $\tilde{\mathbf{n}}$ is the normal to the surface of $\omega$ pointing
outwards from the domain of integration $\omega$.

In our case the domain of integration is the domain
$\mathbb{R}^n\backslash\Omega$, whose boundary is $ \Gamma_\infty \cup
\Gamma$, where the "boundary" at infinity is defined as

\f[
\Gamma_\infty \dealcoloneq \lim_{r\to\infty} \partial B_r(0).
\f]

In our program the normals are defined as <i>outer</i> to the domain
$\Omega$, that is, they are in fact <i>inner</i> to the integration
domain, and some care is required in defining the various integrals
with the correct signs for the normals, i.e. replacing $\tilde{\mathbf{n}}$
by $-\mathbf{n}$.

If we substitute $u$ and $v$ in the Green %identity with the solution
$\phi$ and with the fundamental solution of the Laplace equation
respectively, as long as $\mathbf{x}$ is chosen in the region
$\mathbb{R}^n\backslash\Omega$, we obtain:
\f[
  \phi(\mathbf{x}) -
  \int_{\Gamma\cup\Gamma_\infty}\frac{\partial G(\mathbf{y}-\mathbf{x})}{\partial \mathbf{n}_y}\phi(\mathbf{y})\,ds_y
  =
  -\int_{\Gamma\cup\Gamma_\infty}G(\mathbf{y}-\mathbf{x})\frac{\partial \phi}{\partial \mathbf{n}_y}(\mathbf{y})\,ds_y
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega
\f]

where the normals are now pointing <i>inward</i> the domain of
integration.

Notice that in the above equation, we also have the integrals on the
portion of the boundary at $\Gamma_\infty$. Using the boundary
conditions of our problem, we have that $\nabla \phi$ is zero at
infinity (which simplifies the integral on $\Gamma_\infty$ on the
right hand side).

The integral on $\Gamma_\infty$ that appears on the left hand side can
be treated by observing that $\nabla\phi=0$ implies that $\phi$ at
infinity is necessarily constant. We define its value to be
$\phi_\infty$.  It is an easy exercise to prove that

\f[
-\int_{\Gamma_\infty} \frac{\partial G(\mathbf{y}-\mathbf{x})}
{\partial \mathbf{n}_y}\phi_\infty \,ds_y =
\lim_{r\to\infty} \int_{\partial B_r(0)} \frac{\mathbf{r}}{r} \cdot \nabla G(\mathbf{y}-\mathbf{x})
\phi_\infty \,ds_y = -\phi_\infty.
\f]

Using this result, we can reduce the above equation only on the
boundary $\Gamma$ using the so-called Single and Double Layer
Potential operators:

\f[\label{integral}
  \phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty
  -\left(S \frac{\partial \phi}{\partial n_y}\right)(\mathbf{x})
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega.
\f]

(The name of these operators comes from the fact that they describe the
electric potential in $\mathbb{R}^n$ due to a single thin sheet of charges
along a surface, and due to a double sheet of charges and anti-charges along
the surface, respectively.)

In our case, we know the Neumann values of $\phi$ on the boundary:
$\mathbf{n}\cdot\nabla\phi = -\mathbf{n}\cdot\mathbf{v}_\infty$.
Consequently,
\f[
  \phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty +
   \left(S[\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
   \qquad \forall\mathbf{x} \in \mathbb{R}^n\backslash\Omega.
\f]
If we take the limit for $\mathbf{x}$ tending to $\Gamma$ of
the above equation, using well known properties of the single and double layer
operators, we obtain an equation for $\phi$ just on the boundary $\Gamma$ of
$\Omega$:

\f[\label{SD}
  \alpha(\mathbf{x})\phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty +
  \left(S [\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
  \quad \mathbf{x}\in \partial\Omega,
\f]

which is the Boundary Integral Equation (BIE) we were looking for,
where the quantity $\alpha(\mathbf{x})$ is the fraction of angle or
solid angle by which the point $\mathbf{x}$ sees the domain of
integration $\mathbb{R}^n\backslash\Omega$.

In particular, at points $\mathbf{x}$ where the boundary
$\partial\Omega$ is differentiable (i.e. smooth) we have
$\alpha(\mathbf{x})=\frac 12$, but the value may be smaller or larger
at points where the boundary has a corner or an edge.

Substituting the single and double layer operators we get:
\f[
  \alpha(\mathbf{x}) \phi(\mathbf{x})
  + \frac{1}{2\pi}\int_{\partial \Omega}  \frac{
  (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }{ |\mathbf{y}-\mathbf{x}|^2 }
  \phi(\mathbf{y}) \,ds_y
  = \phi_\infty
    -\frac{1}{2\pi}\int_{\partial \Omega}  \ln|\mathbf{y}-\mathbf{x}| \, \mathbf{n}\cdot\mathbf{v_\infty}\,ds_y
\f]
for two dimensional flows and
\f[
  \alpha(\mathbf{x}) \phi(\mathbf{x})
   + \frac{1}{4\pi}\int_{\partial \Omega} \frac{ (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }{ |\mathbf{y}-\mathbf{x}|^3 }\phi(\mathbf{y})\,ds_y
  = \phi_\infty +
  \frac{1}{4\pi}\int_{\partial \Omega} \frac{1}{|\mathbf{y}-\mathbf{x}|} \, \mathbf{n}\cdot\mathbf{v_\infty}\,ds_y
\f]
for three dimensional flows, where the normal derivatives of the fundamental
solutions have been written in a form that makes computation easier. In either
case, $\phi$ is the solution of an integral equation posed entirely on the
boundary since both $\mathbf{x},\mathbf{y}\in\partial\Omega$.

Notice that the fraction of angle (in 2d) or solid angle (in 3d)
$\alpha(\mathbf{x})$ by which the point $\mathbf{x}$ sees the domain
$\Omega$ can be defined using the double layer potential itself:
\f[
\alpha(\mathbf{x}) \dealcoloneq 1 -
\frac{1}{2(n-1)\pi}\int_{\partial \Omega} \frac{ (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }
{ |\mathbf{y}-\mathbf{x}|^{n} }\phi(\mathbf{y})\,ds_y = 1+
\int_{\partial \Omega} \frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y.
\f]

The reason why this is possible can be understood if we consider the
fact that the solution of a pure Neumann problem is known up to an
arbitrary constant $c$, which means that, if we set the Neumann data
to be zero, then any constant $\phi = \phi_\infty$ will be a solution.
Inserting the constant solution and the Neumann boundary condition in the
boundary integral equation, we have
@f{align*}
\alpha\left(\mathbf{x}\right)\phi\left(\mathbf{x}\right)
&=\int_{\Omega}\phi\left(\mathbf{y}\right)\delta\left(\mathbf{y}-\mathbf{x}\right)\, dy\\
\Rightarrow
\alpha\left(\mathbf{x}\right)\phi_\infty
&=\phi_\infty\int_{\Gamma\cup\Gamma_\infty}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
=\phi_\infty\left[\int_{\Gamma_\infty}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
+\int_{\Gamma}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
\right]
@f}
The integral on $\Gamma_\infty$ is unity, see above, so division by the constant $\phi_\infty$ gives us the explicit
expression above for $\alpha(\mathbf{x})$.

While this example program is really only focused on the solution of the
boundary integral equation, in a realistic setup one would still need to solve
for the velocities. To this end, note that we have just computed
$\phi(\mathbf{x})$ for all $\mathbf{x}\in\partial\Omega$. In the next step, we
can compute (analytically, if we want) the solution $\phi(\mathbf{x})$ in all
of $\mathbb{R}^n\backslash\Omega$. To this end, recall that we had
\f[
  \phi(\mathbf{x})
  =
  \phi_\infty +
  (D\phi)(\mathbf{x})
  +
  \left(S[\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega.
\f]
where now we have everything that is on the right hand side ($S$ and $D$ are
integrals we can evaluate, the normal velocity on the boundary is given, and
$\phi$ on the boundary we have just computed). Finally, we can then recover
the velocity as $\mathbf{\tilde v}=\nabla \phi$.

Notice that the evaluation of the above formula for $\mathbf{x} \in
\Omega$ should yield zero as a result, since the integration of the
Dirac delta $\delta(\mathbf{x})$ in the domain
$\mathbb{R}^n\backslash\Omega$ is always zero by definition.

As a final test, let us verify that this velocity indeed satisfies the
momentum balance equation for a stationary flow field, i.e., whether
$\mathbf{v}\cdot\nabla\mathbf{v} = -\frac 1\rho \nabla p$ where
$\mathbf{v}=\mathbf{\tilde
v}+\mathbf{v}_\infty=\nabla\phi+\mathbf{v}_\infty$ for some (unknown) pressure
$p$ and a given constant $\rho$. In other words, we would like to verify that
Bernoulli's law as stated above indeed holds. To show this, we use that
the left hand side of this equation equates to
@f{align*}
  \mathbf{v}\cdot\nabla\mathbf{v}
  &=
  [(\nabla\phi+\mathbf{v}_\infty)\cdot\nabla] (\nabla\phi+\mathbf{v}_\infty)
  \\
  &=
  [(\nabla\phi+\mathbf{v}_\infty)\cdot\nabla] (\nabla\phi)
@f}
where we have used that $\mathbf{v}_\infty$ is constant. We would like to
write this expression as the gradient of something (remember that $\rho$ is a
constant). The next step is more
convenient if we consider the components of the equation individually
(summation over indices that appear twice is implied):
@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  (\partial_j\phi+v_{\infty,j}) \partial_j \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]
  -
  \partial_j [(\partial_j\phi+v_{\infty,j})] \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]
@f}
because $\partial_j \partial_j\phi = \Delta \phi = 0$ and $\textrm{div}
\ \mathbf{v}_\infty=0$. Next,
@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  \partial_j [v_{\infty,j} \partial_i\phi]
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  \partial_j [v_{\infty,j}] \partial_i\phi
  +
  v_{\infty,j} \partial_j \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  v_{\infty,j} \partial_j \partial_i\phi
  \\
  &=
  \partial_i \partial_j [(\partial_j\phi) \phi]
  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  +
  \partial_i [v_{\infty,j} \partial_j \phi]
  -
  \partial_i [v_{\infty,j}] \partial_j \phi
@f}
Again, the last term disappears because $\mathbf{v}_\infty$ is constant and we
can merge the first and third term into one:
@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_i (\partial_j [(\partial_j\phi) \phi + v_{\infty,j} \partial_j \phi])
  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  \\
  &=
  \partial_i [(\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi]
  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
@f}

We now only need to massage that last term a bit more. Using the product rule,
we get
@f{align*}
  \partial_j [\partial_i (\partial_j\phi) \phi]
  &=
  \partial_i [\partial_j \partial_j\phi] \phi
  +
  \partial_i [\partial_j \phi] (\partial_j \phi).
@f}
The first of these terms is zero (because, again, the summation over $j$ gives
$\Delta\phi$, which is zero). The last term can be written as $\frac 12
\partial_i [(\partial_j\phi)(\partial_j\phi)]$ which is in the desired gradient
form. As a consequence, we can now finally state that
@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_i (\partial_j [(\partial_j\phi) \phi + v_{\infty,j} \partial_j \phi])
  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  \\
  &=
  \partial_i
  \left[
    (\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi
    -
    \frac 12 (\partial_j\phi)(\partial_j\phi)
  \right],
  \\
  &=
  \partial_i
  \left[
    \frac 12 (\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi
  \right],
@f}
or in vector form:
@f[
  \mathbf{v}\cdot\nabla\mathbf{v}
  =
  \nabla
  \left[
    \frac 12 \mathbf{\tilde v}^2
    + \mathbf{v}_{\infty} \cdot \mathbf{\tilde v}
  \right],
@f]
or in other words:
@f[
  p
  =
  -\rho
  \left[
    \frac 12 \mathbf{\tilde v}^2
    + \mathbf{v}_{\infty} \cdot \mathbf{\tilde v}
  \right]
  =
  -\rho
  \left[
    \frac 12 \mathbf{v}^2
    -
    \frac 12 \mathbf{v}_{\infty}^2
  \right]
  .
@f]
Because the pressure is only determined up to a constant (it appears only with
a gradient in the equations), an equally valid definition is
@f[
  p
  =
  -\frac 12 \rho \mathbf{v}^2
  .
@f]
This is exactly Bernoulli's law mentioned above.


<a name="Thenumericalapproximation"></a><h3>The numerical approximation</h3>


Numerical approximations of Boundary Integral Equations (BIE) are commonly
referred to as the boundary element method or panel method (the latter
expression being used mostly in the computational fluid dynamics community).
The goal of the following test problem is to solve the integral
formulation of the Laplace equation with Neumann boundary conditions,
using a circle and a sphere respectively in two and three space
dimensions, illustrating along the way the features that allow one to
treat boundary element problems almost as easily as finite element
problems using the deal.II library.

To this end, let $\mathcal{T}_h = \bigcup_i K_i$ be a subdivision of the
manifold $\Gamma = \partial \Omega$ into $M$ line segments if $n=2$, or $M$
quadrilaterals if $n=3$. We will call each individual segment or
quadrilateral an <i>element</i> or <i>cell</i>, independently of the
dimension $n$ of the surrounding space $\mathbb{R}^n$.
We define the finite dimensional space $V_h$ as
\f[
  \label{eq:definition-Vh}
  V_h \dealcoloneq \{ v \in C^0(\Gamma) \text{ s.t. } v|_{K_i} \in \mathcal{Q}^1(K_i),
  \forall i\},
\f]
with basis functions $\psi_i(\mathbf{x})$ for which we will use the usual FE_Q
finite element, with the catch that this time it is defined on a manifold of
codimension one (which we do by using the second template argument that is
usually defaulted to equal the first; here, we will create objects
<code>FE_Q@<dim-1,dim@></code> to indicate that we have <code>dim-1</code>
dimensional cells in a <code>dim</code> dimensional space).
An element $\phi_h$ of $V_h$ is uniquely
identified by the vector $\boldsymbol{\phi}$ of its coefficients
$\phi_i$, that is:
\f[
  \label{eq:definition-of-element}
  \phi_h(\mathbf{x}) \dealcoloneq \phi_i \psi_i(\mathbf{x}), \qquad
  \boldsymbol{\phi} \dealcoloneq \{ \phi_i \},
\f]
where summation  is implied over repeated indexes. Note that we could use
discontinuous elements here &mdash; in fact, there is no real reason to use
continuous ones since the integral formulation does not
imply any derivatives on our trial functions so continuity is unnecessary,
and often in the literature only piecewise constant elements are used.

<a name="Collocationboundaryelementmethod"></a><h3> Collocation boundary element method </h3>


By far, the most common approximation of boundary integral equations
is by use of the collocation based boundary element method.

This method requires the evaluation of the boundary integral equation
at a number of collocation points which is equal to the number of
unknowns of the system. The choice of these points is a delicate
matter, that requires a careful study. Assume that these points are
known for the moment, and call them $\mathbf x_i$ with $i=0...n\_dofs$.

The problem then becomes:
Given the datum $\mathbf{v}_\infty$, find a function $\phi_h$ in $V_h$
such that the following $n\_dofs$ equations are satisfied:

\f{align*}
    \alpha(\mathbf{x}_i) \phi_h(\mathbf{x}_i)
    - \int_{\Gamma_y} \frac{ \partial G(\mathbf{y}-\mathbf{x}_i)}{\partial\mathbf{n}_y }
    \phi_h(\mathbf{y}) \,ds_y =
    \int_{\Gamma_y} G(\mathbf{y}-\mathbf{x}_i) \,
    \mathbf{n}_y\cdot\mathbf{v_\infty} \,ds_y
    ,
\f}

where the quantity $\alpha(\mathbf{x}_i)$ is the fraction of (solid)
angle by which the point $\mathbf{x}_i$ sees the domain $\Omega$, as
explained above, and we set $\phi_\infty$ to be zero.  If the support
points $\mathbf{x}_i$ are chosen appropriately, then the problem can
be written as the following linear system:

\f[
\label{eq:linear-system}
(\mathbf{A}+\mathbf{N})\boldsymbol\phi = \mathbf{b},
\f]

where

\f[
\begin{aligned}
\mathbf{A}_{ij}&=
\alpha(\mathbf{x}_i) \psi_j(\mathbf{x}_i)
= 1+\int_\Gamma
\frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}\,ds_y
\psi_j(\mathbf{x}_i)
\\
\mathbf{N}_{ij}&= - \int_\Gamma
  \frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}
  \psi_j(\mathbf{y}) \,ds_y
\\
\mathbf{b}_i&= \int_\Gamma
   G(\mathbf{y}-\mathbf{x}_i)  \, \mathbf{n}_y\cdot\mathbf{v_\infty}
   ds_y.
\end{aligned}
\f]

From a linear algebra point of view, the best possible choice of the
collocation points is the one that renders the matrix
$\mathbf{A}+\mathbf{N}$ the most diagonally dominant. A natural choice
is then to select the $\mathbf{x}_i$ collocation points to be the
support points of the nodal basis functions $\psi_i(\mathbf{x})$. In that
case, $\psi_j(\mathbf{x}_i)=\delta_{ij}$, and as a consequence the matrix
$\mathbf{A}$ is diagonal with entries
\f[
  \mathbf{A}_{ii}
  =
  1+\int_\Gamma
  \frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}\,ds_y
  =
  1-\sum_j N_{ij},
\f]
where we have used that $\sum_j \psi_j(\mathbf{y})=1$ for the usual Lagrange
elements.
With this choice of collocation points, the computation of the entries
of the matrices $\mathbf{A}$, $\mathbf{N}$ and of the right hand side
$\mathbf{b}$ requires the evaluation of singular integrals on the
elements $K_i$ of the triangulation $\mathcal{T}_h$.
As usual in these cases, all integrations are performed on a reference
simple domain, i.e., we assume that each element $K_i$ of
$\mathcal{T}_h$ can be expressed as a linear (in two dimensions) or
bi-linear (in three dimensions) transformation of the reference
boundary element $\hat K \dealcoloneq [0,1]^{n-1}$, and we perform the integrations after a
change of variables from the real element $K_i$ to the reference
element $\hat K$.

<a name="Treatingthesingularintegrals"></a><h3> Treating the singular integrals. </h3>


In two dimensions it is not necessary to compute the diagonal elements
$\mathbf{N}_{ii}$ of the system matrix, since, even if the denominator
goes to zero when $\mathbf{x}=\mathbf{y}$, the numerator is always
zero because $\mathbf{n}_y$ and $(\mathbf{y}-\mathbf{x})$ are
orthogonal (on our polygonal approximation of the boundary of $\Omega$), and
the only singular integral arises in the computation
of $\mathbf{b}_i$ on the i-th element of $\mathcal{T}_h$:
\f[
  \frac{1}{\pi}
  \int_{K_i}
  \ln|\mathbf{y}-\mathbf{x}_i| \, \mathbf{n}_y\cdot\mathbf{v_\infty} \,ds_y.
\f]

This can be easily treated by the QGaussLogR quadrature
formula.

Similarly, it is possible to use the QGaussOneOverR quadrature formula
to perform the singular integrations in three dimensions. The
interested reader will find detailed explanations on how these
quadrature rules work in their documentation.

The resulting matrix $\mathbf{A}+\mathbf{N}$ is full. Depending on its
size, it might be convenient to use a direct solver or an iterative
one. For the purpose of this example code, we chose to use only an
iterative solver, without providing any preconditioner.

If this were a production code rather than a demonstration of principles,
there are techniques that are available to not store full matrices but instead
store only those entries that are large and/or relevant. In the literature on
boundary element methods, a plethora of methods is available that allows to
determine which elements are important and which are not, leading to a
significantly sparser representation of these matrices that also facilitates
rapid evaluations of the scalar product between vectors and matrices. This not
being the goal of this program, we leave this for more sophisticated
implementations.


<a name="Implementation"></a><h3>Implementation</h3>


The implementation is rather straight forward. The main point that hasn't been
used in any of the previous tutorial programs is that most classes in deal.II
are not only templated on the dimension, but in fact on the dimension of the
manifold on which we pose the differential equation as well as the dimension
of the space into which this manifold is embedded. By default, the second
template argument equals the first, meaning for example that we want to solve
on a two-dimensional region of two-dimensional space. The triangulation class
to use in this case would be <code>Triangulation@<2@></code>, which is an
equivalent way of writing <code>Triangulation@<2,2@></code>.

However, this doesn't have to be so: in the current example, we will for
example want to solve on the surface of a sphere, which is a two-dimensional
manifold embedded in a three-dimensional space. Consequently, the right class
will be <code>Triangulation@<2,3@></code>, and correspondingly we will use
<code>DoFHandler@<2,3@></code> as the DoF handler class and
<code>FE_Q@<2,3@></code> for finite elements.

Some further details on what one can do with things that live on
curved manifolds can be found in the report
<a target="_top"
href="http://www.dealii.org/reports/codimension-one/desimone-heltai-manigrasso.pdf"><i>Tools
for the Solution of PDEs Defined on Curved Manifolds with the deal.II
Library</i> by A. DeSimone, L. Heltai, C. Manigrasso</a>. In addition, the
step-38 tutorial program extends what we show here to cases where the equation
posed on the manifold is not an integral operator but in fact involves
derivatives.


<a name="Testcase"></a><h3>Testcase</h3>


The testcase we will be solving is for a circular (in 2d) or spherical
(in 3d) obstacle. Meshes for these geometries will be read in from
files in the current directory and an object of type SphericalManifold
will then be attached to the triangulation to allow mesh refinement
that respects the continuous geometry behind the discrete initial
mesh.

For a sphere of radius $a$ translating at a velocity of $U$ in the $x$ direction, the potential reads
@f{align*}
\phi = -\frac{1}{2}U \left(\frac{a}{r}\right)3 r \cos\theta
@f}
see, e.g. J. N. Newman, <i>Marine Hydrodynamics</i>, 1977,
pp. 127. For unit speed and radius, and restricting $(x,y,z)$ to lie
on the surface of the sphere,
$\phi = -x/2$. In the test problem,
the flow is $(1,1,1)$, so the appropriate exact solution on the
surface of the sphere is the superposition of the above solution with
the analogous solution along the $y$ and $z$ axes, or $\phi =
\frac{1}{2}(x + y + z)$.
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
 * The program starts with including a bunch of include files that we will use
 * in the various parts of the program. Most of them have been discussed in
 * previous tutorials already:
 * 
 * @code
 * #include <deal.II/base/smartpointer.h>
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/quadrature_selector.h>
 * #include <deal.II/base/parsed_function.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/solver_control.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_in.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/grid/manifold_lib.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * And here are a few C++ standard header files that we will need:
 * 
 * @code
 * #include <cmath>
 * #include <iostream>
 * #include <fstream>
 * #include <string>
 * 
 * @endcode
 * 
 * The last part of this preamble is to import everything in the dealii
 * namespace into the one into which everything in this program will go:
 * 
 * @code
 * namespace Step34
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Singleanddoublelayeroperatorkernels"></a> 
 * <h3>Single and double layer operator kernels</h3>
 * 

 * 
 * First, let us define a bit of the boundary integral equation machinery.
 * 

 * 
 * The following two functions are the actual calculations of the single and
 * double layer potential kernels, that is $G$ and $\nabla G$. They are well
 * defined only if the vector $R = \mathbf{y}-\mathbf{x}$ is different from
 * zero.
 * 
 * @code
 *   namespace LaplaceKernel
 *   {
 *     template <int dim>
 *     double single_layer(const Tensor<1, dim> &R)
 *     {
 *       switch (dim)
 *         {
 *           case 2:
 *             return (-std::log(R.norm()) / (2 * numbers::PI));
 * 
 *           case 3:
 *             return (1. / (R.norm() * 4 * numbers::PI));
 * 
 *           default:
 *             Assert(false, ExcInternalError());
 *             return 0.;
 *         }
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> double_layer(const Tensor<1, dim> &R)
 *     {
 *       switch (dim)
 *         {
 *           case 2:
 *             return R / (-2 * numbers::PI * R.norm_square());
 *           case 3:
 *             return R / (-4 * numbers::PI * R.norm_square() * R.norm());
 * 
 *           default:
 *             Assert(false, ExcInternalError());
 *             return Tensor<1, dim>();
 *         }
 *     }
 *   } // namespace LaplaceKernel
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheBEMProblemclass"></a> 
 * <h3>The BEMProblem class</h3>
 * 

 * 
 * The structure of a boundary element method code is very similar to the
 * structure of a finite element code, and so the member functions of this
 * class are like those of most of the other tutorial programs. In
 * particular, by now you should be familiar with reading parameters from an
 * external file, and with the splitting of the different tasks into
 * different modules. The same applies to boundary element methods, and we
 * won't comment too much on them, except on the differences.
 * 
 * @code
 *   template <int dim>
 *   class BEMProblem
 *   {
 *   public:
 *     BEMProblem(const unsigned int fe_degree      = 1,
 *                const unsigned int mapping_degree = 1);
 * 
 *     void run();
 * 
 *   private:
 *     void read_parameters(const std::string &filename);
 * 
 *     void read_domain();
 * 
 *     void refine_and_resize();
 * 
 * @endcode
 * 
 * The only really different function that we find here is the assembly
 * routine. We wrote this function in the most possible general way, in
 * order to allow for easy generalization to higher order methods and to
 * different fundamental solutions (e.g., Stokes or Maxwell).
 *     

 * 
 * The most noticeable difference is the fact that the final matrix is
 * full, and that we have a nested loop inside the usual loop on cells
 * that visits all support points of the degrees of freedom.  Moreover,
 * when the support point lies inside the cell which we are visiting, then
 * the integral we perform becomes singular.
 *     

 * 
 * The practical consequence is that we have two sets of quadrature
 * formulas, finite element values and temporary storage, one for standard
 * integration and one for the singular integration, which are used where
 * necessary.
 * 
 * @code
 *     void assemble_system();
 * 
 * @endcode
 * 
 * There are two options for the solution of this problem. The first is to
 * use a direct solver, and the second is to use an iterative solver. We
 * opt for the second option.
 *     

 * 
 * The matrix that we assemble is not symmetric, and we opt to use the
 * GMRES method; however the construction of an efficient preconditioner
 * for boundary element methods is not a trivial issue. Here we use a non
 * preconditioned GMRES solver. The options for the iterative solver, such
 * as the tolerance, the maximum number of iterations, are selected
 * through the parameter file.
 * 
 * @code
 *     void solve_system();
 * 
 * @endcode
 * 
 * Once we obtained the solution, we compute the $L^2$ error of the
 * computed potential as well as the $L^\infty$ error of the approximation
 * of the solid angle. The mesh we are using is an approximation of a
 * smooth curve, therefore the computed diagonal matrix of fraction of
 * angles or solid angles $\alpha(\mathbf{x})$ should be constantly equal
 * to $\frac 12$. In this routine we output the error on the potential and
 * the error in the approximation of the computed angle. Notice that the
 * latter error is actually not the error in the computation of the angle,
 * but a measure of how well we are approximating the sphere and the
 * circle.
 *     

 * 
 * Experimenting a little with the computation of the angles gives very
 * accurate results for simpler geometries. To verify this you can comment
 * out, in the read_domain() method, the tria.set_manifold(1, manifold)
 * line, and check the alpha that is generated by the program. By removing
 * this call, whenever the mesh is refined new nodes will be placed along
 * the straight lines that made up the coarse mesh, rather than be pulled
 * onto the surface that we really want to approximate. In the three
 * dimensional case, the coarse grid of the sphere is obtained starting
 * from a cube, and the obtained values of alphas are exactly $\frac 12$
 * on the nodes of the faces, $\frac 34$ on the nodes of the edges and
 * $\frac 78$ on the 8 nodes of the vertices.
 * 
 * @code
 *     void compute_errors(const unsigned int cycle);
 * 
 * @endcode
 * 
 * Once we obtained a solution on the codimension one domain, we want to
 * interpolate it to the rest of the space. This is done by performing
 * again the convolution of the solution with the kernel in the
 * compute_exterior_solution() function.
 *     

 * 
 * We would like to plot the velocity variable which is the gradient of
 * the potential solution. The potential solution is only known on the
 * boundary, but we use the convolution with the fundamental solution to
 * interpolate it on a standard dim dimensional continuous finite element
 * space. The plot of the gradient of the extrapolated solution will give
 * us the velocity we want.
 *     

 * 
 * In addition to the solution on the exterior domain, we also output the
 * solution on the domain's boundary in the output_results() function, of
 * course.
 * 
 * @code
 *     void compute_exterior_solution();
 * 
 *     void output_results(const unsigned int cycle);
 * 
 * @endcode
 * 
 * To allow for dimension independent programming, we specialize this
 * single function to extract the singular quadrature formula needed to
 * integrate the singular kernels in the interior of the cells.
 * 
 * @code
 *     const Quadrature<dim - 1> &get_singular_quadrature(
 *       const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell,
 *       const unsigned int index) const;
 * 
 * 
 * @endcode
 * 
 * The usual deal.II classes can be used for boundary element methods by
 * specifying the "codimension" of the problem. This is done by setting
 * the optional second template arguments to Triangulation, FiniteElement
 * and DoFHandler to the dimension of the embedding space. In our case we
 * generate either 1 or 2 dimensional meshes embedded in 2 or 3
 * dimensional spaces.
 *     

 * 
 * The optional argument by default is equal to the first argument, and
 * produces the usual finite element classes that we saw in all previous
 * examples.
 *     

 * 
 * The class is constructed in a way to allow for arbitrary order of
 * approximation of both the domain (through high order mapping) and the
 * finite element space. The order of the finite element space and of the
 * mapping can be selected in the constructor of the class.
 * 

 * 
 * 
 * @code
 *     Triangulation<dim - 1, dim> tria;
 *     FE_Q<dim - 1, dim>          fe;
 *     DoFHandler<dim - 1, dim>    dof_handler;
 *     MappingQ<dim - 1, dim>      mapping;
 * 
 * @endcode
 * 
 * In BEM methods, the matrix that is generated is dense. Depending on the
 * size of the problem, the final system might be solved by direct LU
 * decomposition, or by iterative methods. In this example we use an
 * unpreconditioned GMRES method. Building a preconditioner for BEM method
 * is non trivial, and we don't treat this subject here.
 * 

 * 
 * 
 * @code
 *     FullMatrix<double> system_matrix;
 *     Vector<double>     system_rhs;
 * 
 * @endcode
 * 
 * The next two variables will denote the solution $\phi$ as well as a
 * vector that will hold the values of $\alpha(\mathbf x)$ (the fraction
 * of $\Omega$ visible from a point $\mathbf x$) at the support points of
 * our shape functions.
 * 

 * 
 * 
 * @code
 *     Vector<double> phi;
 *     Vector<double> alpha;
 * 
 * @endcode
 * 
 * The convergence table is used to output errors in the exact solution
 * and in the computed alphas.
 * 

 * 
 * 
 * @code
 *     ConvergenceTable convergence_table;
 * 
 * @endcode
 * 
 * The following variables are the ones that we fill through a parameter
 * file.  The new objects that we use in this example are the
 * Functions::ParsedFunction object and the QuadratureSelector object.
 *     

 * 
 * The Functions::ParsedFunction class allows us to easily and quickly
 * define new function objects via parameter files, with custom
 * definitions which can be very complex (see the documentation of that
 * class for all the available options).
 *     

 * 
 * We will allocate the quadrature object using the QuadratureSelector
 * class that allows us to generate quadrature formulas based on an
 * identifying string and on the possible degree of the formula itself. We
 * used this to allow custom selection of the quadrature formulas for the
 * standard integration, and to define the order of the singular
 * quadrature rule.
 *     

 * 
 * We also define a couple of parameters which are used in case we wanted
 * to extend the solution to the entire domain.
 * 

 * 
 * 
 * @code
 *     Functions::ParsedFunction<dim> wind;
 *     Functions::ParsedFunction<dim> exact_solution;
 * 
 *     unsigned int                         singular_quadrature_order;
 *     std::shared_ptr<Quadrature<dim - 1>> quadrature;
 * 
 *     SolverControl solver_control;
 * 
 *     unsigned int n_cycles;
 *     unsigned int external_refinement;
 * 
 *     bool run_in_this_dimension;
 *     bool extend_solution;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemBEMProblemandBEMProblemread_parameters"></a> 
 * <h4>BEMProblem::BEMProblem and BEMProblem::read_parameters</h4>
 * 

 * 
 * The constructor initializes the various object in much the same way as
 * done in the finite element programs such as step-4 or step-6. The only
 * new ingredient here is the ParsedFunction object, which needs, at
 * construction time, the specification of the number of components.
 *   

 * 
 * For the exact solution the number of vector components is one, and no
 * action is required since one is the default value for a ParsedFunction
 * object. The wind, however, requires dim components to be
 * specified. Notice that when declaring entries in a parameter file for the
 * expression of the Functions::ParsedFunction, we need to specify the
 * number of components explicitly, since the function
 * Functions::ParsedFunction::declare_parameters is static, and has no
 * knowledge of the number of components.
 * 
 * @code
 *   template <int dim>
 *   BEMProblem<dim>::BEMProblem(const unsigned int fe_degree,
 *                               const unsigned int mapping_degree)
 *     : fe(fe_degree)
 *     , dof_handler(tria)
 *     , mapping(mapping_degree, true)
 *     , wind(dim)
 *     , singular_quadrature_order(5)
 *     , n_cycles(4)
 *     , external_refinement(5)
 *     , run_in_this_dimension(true)
 *     , extend_solution(true)
 *   {}
 * 
 * 
 *   template <int dim>
 *   void BEMProblem<dim>::read_parameters(const std::string &filename)
 *   {
 *     deallog << std::endl
 *             << "Parsing parameter file " << filename << std::endl
 *             << "for a " << dim << " dimensional simulation. " << std::endl;
 * 
 *     ParameterHandler prm;
 * 
 *     prm.declare_entry("Number of cycles", "4", Patterns::Integer());
 *     prm.declare_entry("External refinement", "5", Patterns::Integer());
 *     prm.declare_entry("Extend solution on the -2,2 box",
 *                       "true",
 *                       Patterns::Bool());
 *     prm.declare_entry("Run 2d simulation", "true", Patterns::Bool());
 *     prm.declare_entry("Run 3d simulation", "true", Patterns::Bool());
 * 
 *     prm.enter_subsection("Quadrature rules");
 *     {
 *       prm.declare_entry(
 *         "Quadrature type",
 *         "gauss",
 *         Patterns::Selection(
 *           QuadratureSelector<(dim - 1)>::get_quadrature_names()));
 *       prm.declare_entry("Quadrature order", "4", Patterns::Integer());
 *       prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
 *     }
 *     prm.leave_subsection();
 * 
 * @endcode
 * 
 * For both two and three dimensions, we set the default input data to be
 * such that the solution is $x+y$ or $x+y+z$. The actually computed
 * solution will have value zero at infinity. In this case, this coincide
 * with the exact solution, and no additional corrections are needed, but
 * you should be aware of the fact that we arbitrarily set $\phi_\infty$,
 * and the exact solution we pass to the program needs to have the same
 * value at infinity for the error to be computed correctly.
 *     

 * 
 * The use of the Functions::ParsedFunction object is pretty straight
 * forward. The Functions::ParsedFunction::declare_parameters function
 * takes an additional integer argument that specifies the number of
 * components of the given function. Its default value is one. When the
 * corresponding Functions::ParsedFunction::parse_parameters method is
 * called, the calling object has to have the same number of components
 * defined here, otherwise an exception is thrown.
 *     

 * 
 * When declaring entries, we declare both 2 and three dimensional
 * functions. However only the dim-dimensional one is ultimately
 * parsed. This allows us to have only one parameter file for both 2 and 3
 * dimensional problems.
 *     

 * 
 * Notice that from a mathematical point of view, the wind function on the
 * boundary should satisfy the condition $\int_{\partial\Omega}
 * \mathbf{v}\cdot \mathbf{n} d \Gamma = 0$, for the problem to have a
 * solution. If this condition is not satisfied, then no solution can be
 * found, and the solver will not converge.
 * 
 * @code
 *     prm.enter_subsection("Wind function 2d");
 *     {
 *       Functions::ParsedFunction<2>::declare_parameters(prm, 2);
 *       prm.set("Function expression", "1; 1");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Wind function 3d");
 *     {
 *       Functions::ParsedFunction<3>::declare_parameters(prm, 3);
 *       prm.set("Function expression", "1; 1; 1");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Exact solution 2d");
 *     {
 *       Functions::ParsedFunction<2>::declare_parameters(prm);
 *       prm.set("Function expression", "x+y");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Exact solution 3d");
 *     {
 *       Functions::ParsedFunction<3>::declare_parameters(prm);
 *       prm.set("Function expression", "x+y+z");
 *     }
 *     prm.leave_subsection();
 * 
 * 
 * @endcode
 * 
 * In the solver section, we set all SolverControl parameters. The object
 * will then be fed to the GMRES solver in the solve_system() function.
 * 
 * @code
 *     prm.enter_subsection("Solver");
 *     SolverControl::declare_parameters(prm);
 *     prm.leave_subsection();
 * 
 * @endcode
 * 
 * After declaring all these parameters to the ParameterHandler object,
 * let's read an input file that will give the parameters their values. We
 * then proceed to extract these values from the ParameterHandler object:
 * 
 * @code
 *     prm.parse_input(filename);
 * 
 *     n_cycles            = prm.get_integer("Number of cycles");
 *     external_refinement = prm.get_integer("External refinement");
 *     extend_solution     = prm.get_bool("Extend solution on the -2,2 box");
 * 
 *     prm.enter_subsection("Quadrature rules");
 *     {
 *       quadrature = std::shared_ptr<Quadrature<dim - 1>>(
 *         new QuadratureSelector<dim - 1>(prm.get("Quadrature type"),
 *                                         prm.get_integer("Quadrature order")));
 *       singular_quadrature_order = prm.get_integer("Singular quadrature order");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Wind function " + std::to_string(dim) + "d");
 *     {
 *       wind.parse_parameters(prm);
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Exact solution " + std::to_string(dim) + "d");
 *     {
 *       exact_solution.parse_parameters(prm);
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Solver");
 *     solver_control.parse_parameters(prm);
 *     prm.leave_subsection();
 * 
 * 
 * @endcode
 * 
 * Finally, here's another example of how to use parameter files in
 * dimension independent programming.  If we wanted to switch off one of
 * the two simulations, we could do this by setting the corresponding "Run
 * 2d simulation" or "Run 3d simulation" flag to false:
 * 
 * @code
 *     run_in_this_dimension =
 *       prm.get_bool("Run " + std::to_string(dim) + "d simulation");
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemread_domain"></a> 
 * <h4>BEMProblem::read_domain</h4>
 * 

 * 
 * A boundary element method triangulation is basically the same as a
 * (dim-1) dimensional triangulation, with the difference that the vertices
 * belong to a (dim) dimensional space.
 *   

 * 
 * Some of the mesh formats supported in deal.II use by default three
 * dimensional points to describe meshes. These are the formats which are
 * compatible with the boundary element method capabilities of deal.II. In
 * particular we can use either UCD or GMSH formats. In both cases, we have
 * to be particularly careful with the orientation of the mesh, because,
 * unlike in the standard finite element case, no reordering or
 * compatibility check is performed here.  All meshes are considered as
 * oriented, because they are embedded in a higher dimensional space. (See
 * the documentation of the GridIn and of the Triangulation for further
 * details on orientation of cells in a triangulation.) In our case, the
 * normals to the mesh are external to both the circle in 2d or the sphere
 * in 3d.
 *   

 * 
 * The other detail that is required for appropriate refinement of
 * the boundary element mesh is an accurate description of the
 * manifold that the mesh approximates. We already saw this
 * several times for the boundary of standard finite element meshes
 * (for example in step-5 and step-6), and here the principle and
 * usage is the same, except that the SphericalManifold class takes
 * an additional template parameter that specifies the embedding
 * space dimension.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::read_domain()
 *   {
 *     const Point<dim>                      center = Point<dim>();
 *     const SphericalManifold<dim - 1, dim> manifold(center);
 * 
 *     std::ifstream in;
 *     switch (dim)
 *       {
 *         case 2:
 *           in.open("coarse_circle.inp");
 *           break;
 * 
 *         case 3:
 *           in.open("coarse_sphere.inp");
 *           break;
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     GridIn<dim - 1, dim> gi;
 *     gi.attach_triangulation(tria);
 *     gi.read_ucd(in);
 * 
 *     tria.set_all_manifold_ids(1);
 * @endcode
 * 
 * The call to Triangulation::set_manifold copies the manifold (via
 * Manifold::clone()), so we do not need to worry about invalid pointers
 * to <code>manifold</code>:
 * 
 * @code
 *     tria.set_manifold(1, manifold);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemrefine_and_resize"></a> 
 * <h4>BEMProblem::refine_and_resize</h4>
 * 

 * 
 * This function globally refines the mesh, distributes degrees of freedom,
 * and resizes matrices and vectors.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::refine_and_resize()
 *   {
 *     tria.refine_global(1);
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     const unsigned int n_dofs = dof_handler.n_dofs();
 * 
 *     system_matrix.reinit(n_dofs, n_dofs);
 * 
 *     system_rhs.reinit(n_dofs);
 *     phi.reinit(n_dofs);
 *     alpha.reinit(n_dofs);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemassemble_system"></a> 
 * <h4>BEMProblem::assemble_system</h4>
 * 

 * 
 * The following is the main function of this program, assembling the matrix
 * that corresponds to the boundary integral equation.
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::assemble_system()
 *   {
 * @endcode
 * 
 * First we initialize an FEValues object with the quadrature formula for
 * the integration of the kernel in non singular cells. This quadrature is
 * selected with the parameter file, and needs to be quite precise, since
 * the functions we are integrating are not polynomial functions.
 * 
 * @code
 *     FEValues<dim - 1, dim> fe_v(mapping,
 *                                 fe,
 *                                 *quadrature,
 *                                 update_values | update_normal_vectors |
 *                                   update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int n_q_points = fe_v.n_quadrature_points;
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(
 *       fe.n_dofs_per_cell());
 * 
 *     std::vector<Vector<double>> cell_wind(n_q_points, Vector<double>(dim));
 *     double                      normal_wind;
 * 
 * @endcode
 * 
 * Unlike in finite element methods, if we use a collocation boundary
 * element method, then in each assembly loop we only assemble the
 * information that refers to the coupling between one degree of freedom
 * (the degree associated with support point $i$) and the current
 * cell. This is done using a vector of fe.dofs_per_cell elements, which
 * will then be distributed to the matrix in the global row $i$. The
 * following object will hold this information:
 * 
 * @code
 *     Vector<double> local_matrix_row_i(fe.n_dofs_per_cell());
 * 
 * @endcode
 * 
 * The index $i$ runs on the collocation points, which are the support
 * points of the $i$th basis function, while $j$ runs on inner integration
 * points.
 * 

 * 
 * We construct a vector of support points which will be used in the local
 * integrations:
 * 
 * @code
 *     std::vector<Point<dim>> support_points(dof_handler.n_dofs());
 *     DoFTools::map_dofs_to_support_points<dim - 1, dim>(mapping,
 *                                                        dof_handler,
 *                                                        support_points);
 * 
 * 
 * @endcode
 * 
 * After doing so, we can start the integration loop over all cells, where
 * we first initialize the FEValues object and get the values of
 * $\mathbf{\tilde v}$ at the quadrature points (this vector field should
 * be constant, but it doesn't hurt to be more general):
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_v.reinit(cell);
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
 *         const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();
 *         wind.vector_value_list(q_points, cell_wind);
 * 
 * @endcode
 * 
 * We then form the integral over the current cell for all degrees of
 * freedom (note that this includes degrees of freedom not located on
 * the current cell, a deviation from the usual finite element
 * integrals). The integral that we need to perform is singular if one
 * of the local degrees of freedom is the same as the support point
 * $i$. A the beginning of the loop we therefore check whether this is
 * the case, and we store which one is the singular index:
 * 
 * @code
 *         for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *           {
 *             local_matrix_row_i = 0;
 * 
 *             bool         is_singular    = false;
 *             unsigned int singular_index = numbers::invalid_unsigned_int;
 * 
 *             for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
 *               if (local_dof_indices[j] == i)
 *                 {
 *                   singular_index = j;
 *                   is_singular    = true;
 *                   break;
 *                 }
 * 
 * @endcode
 * 
 * We then perform the integral. If the index $i$ is not one of
 * the local degrees of freedom, we simply have to add the single
 * layer terms to the right hand side, and the double layer terms
 * to the matrix:
 * 
 * @code
 *             if (is_singular == false)
 *               {
 *                 for (unsigned int q = 0; q < n_q_points; ++q)
 *                   {
 *                     normal_wind = 0;
 *                     for (unsigned int d = 0; d < dim; ++d)
 *                       normal_wind += normals[q][d] * cell_wind[q](d);
 * 
 *                     const Tensor<1, dim> R = q_points[q] - support_points[i];
 * 
 *                     system_rhs(i) += (LaplaceKernel::single_layer(R) *
 *                                       normal_wind * fe_v.JxW(q));
 * 
 *                     for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
 * 
 *                       local_matrix_row_i(j) -=
 *                         ((LaplaceKernel::double_layer(R) * normals[q]) *
 *                          fe_v.shape_value(j, q) * fe_v.JxW(q));
 *                   }
 *               }
 *             else
 *               {
 * @endcode
 * 
 * Now we treat the more delicate case. If we are here, this
 * means that the cell that runs on the $j$ index contains
 * support_point[i]. In this case both the single and the
 * double layer potential are singular, and they require
 * special treatment.
 *                 

 * 
 * Whenever the integration is performed with the singularity
 * inside the given cell, then a special quadrature formula is
 * used that allows one to integrate arbitrary functions
 * against a singular weight on the reference cell.
 *                 

 * 
 * The correct quadrature formula is selected by the
 * get_singular_quadrature function, which is explained in
 * detail below.
 * 
 * @code
 *                 Assert(singular_index != numbers::invalid_unsigned_int,
 *                        ExcInternalError());
 * 
 *                 const Quadrature<dim - 1> &singular_quadrature =
 *                   get_singular_quadrature(cell, singular_index);
 * 
 *                 FEValues<dim - 1, dim> fe_v_singular(
 *                   mapping,
 *                   fe,
 *                   singular_quadrature,
 *                   update_jacobians | update_values | update_normal_vectors |
 *                     update_quadrature_points);
 * 
 *                 fe_v_singular.reinit(cell);
 * 
 *                 std::vector<Vector<double>> singular_cell_wind(
 *                   singular_quadrature.size(), Vector<double>(dim));
 * 
 *                 const std::vector<Tensor<1, dim>> &singular_normals =
 *                   fe_v_singular.get_normal_vectors();
 *                 const std::vector<Point<dim>> &singular_q_points =
 *                   fe_v_singular.get_quadrature_points();
 * 
 *                 wind.vector_value_list(singular_q_points, singular_cell_wind);
 * 
 *                 for (unsigned int q = 0; q < singular_quadrature.size(); ++q)
 *                   {
 *                     const Tensor<1, dim> R =
 *                       singular_q_points[q] - support_points[i];
 *                     double normal_wind = 0;
 *                     for (unsigned int d = 0; d < dim; ++d)
 *                       normal_wind +=
 *                         (singular_cell_wind[q](d) * singular_normals[q][d]);
 * 
 *                     system_rhs(i) += (LaplaceKernel::single_layer(R) *
 *                                       normal_wind * fe_v_singular.JxW(q));
 * 
 *                     for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
 *                       {
 *                         local_matrix_row_i(j) -=
 *                           ((LaplaceKernel::double_layer(R) *
 *                             singular_normals[q]) *
 *                            fe_v_singular.shape_value(j, q) *
 *                            fe_v_singular.JxW(q));
 *                       }
 *                   }
 *               }
 * 
 * @endcode
 * 
 * Finally, we need to add the contributions of the current cell
 * to the global matrix.
 * 
 * @code
 *             for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
 *               system_matrix(i, local_dof_indices[j]) += local_matrix_row_i(j);
 *           }
 *       }
 * 
 * @endcode
 * 
 * The second part of the integral operator is the term
 * $\alpha(\mathbf{x}_i) \phi_j(\mathbf{x}_i)$. Since we use a collocation
 * scheme, $\phi_j(\mathbf{x}_i)=\delta_{ij}$ and the corresponding matrix
 * is a diagonal one with entries equal to $\alpha(\mathbf{x}_i)$.
 * 

 * 
 * One quick way to compute this diagonal matrix of the solid angles, is
 * to use the Neumann matrix itself. It is enough to multiply the matrix
 * with a vector of elements all equal to -1, to get the diagonal matrix
 * of the alpha angles, or solid angles (see the formula in the
 * introduction for this). The result is then added back onto the system
 * matrix object to yield the final form of the matrix:
 * 
 * @code
 *     Vector<double> ones(dof_handler.n_dofs());
 *     ones.add(-1.);
 * 
 *     system_matrix.vmult(alpha, ones);
 *     alpha.add(1);
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *       system_matrix(i, i) += alpha(i);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemsolve_system"></a> 
 * <h4>BEMProblem::solve_system</h4>
 * 

 * 
 * The next function simply solves the linear system.
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::solve_system()
 *   {
 *     SolverGMRES<Vector<double>> solver(solver_control);
 *     solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemcompute_errors"></a> 
 * <h4>BEMProblem::compute_errors</h4>
 * 

 * 
 * The computation of the errors is exactly the same in all other example
 * programs, and we won't comment too much. Notice how the same methods that
 * are used in the finite element methods can be used here.
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::compute_errors(const unsigned int cycle)
 *   {
 *     Vector<float> difference_per_cell(tria.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       phi,
 *                                       exact_solution,
 *                                       difference_per_cell,
 *                                       QGauss<(dim - 1)>(2 * fe.degree + 1),
 *                                       VectorTools::L2_norm);
 *     const double L2_error =
 *       VectorTools::compute_global_error(tria,
 *                                         difference_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 * @endcode
 * 
 * The error in the alpha vector can be computed directly using the
 * Vector::linfty_norm() function, since on each node, the value should be
 * $\frac 12$. All errors are then output and appended to our
 * ConvergenceTable object for later computation of convergence rates:
 * 
 * @code
 *     Vector<double> difference_per_node(alpha);
 *     difference_per_node.add(-.5);
 * 
 *     const double       alpha_error    = difference_per_node.linfty_norm();
 *     const unsigned int n_active_cells = tria.n_active_cells();
 *     const unsigned int n_dofs         = dof_handler.n_dofs();
 * 
 *     deallog << "Cycle " << cycle << ':' << std::endl
 *             << "   Number of active cells:       " << n_active_cells
 *             << std::endl
 *             << "   Number of degrees of freedom: " << n_dofs << std::endl;
 * 
 *     convergence_table.add_value("cycle", cycle);
 *     convergence_table.add_value("cells", n_active_cells);
 *     convergence_table.add_value("dofs", n_dofs);
 *     convergence_table.add_value("L2(phi)", L2_error);
 *     convergence_table.add_value("Linfty(alpha)", alpha_error);
 *   }
 * 
 * 
 * @endcode
 * 
 * Singular integration requires a careful selection of the quadrature
 * rules. In particular the deal.II library provides quadrature rules which
 * are tailored for logarithmic singularities (QGaussLog, QGaussLogR), as
 * well as for 1/R singularities (QGaussOneOverR).
 *   

 * 
 * Singular integration is typically obtained by constructing weighted
 * quadrature formulas with singular weights, so that it is possible to
 * write
 *   

 * 
 * \f[ \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i) \f]
 *   

 * 
 * where $s(x)$ is a given singularity, and the weights and quadrature
 * points $w_i,q_i$ are carefully selected to make the formula above an
 * equality for a certain class of functions $f(x)$.
 *   

 * 
 * In all the finite element examples we have seen so far, the weight of the
 * quadrature itself (namely, the function $s(x)$), was always constantly
 * equal to 1.  For singular integration, we have two choices: we can use
 * the definition above, factoring out the singularity from the integrand
 * (i.e., integrating $f(x)$ with the special quadrature rule), or we can
 * ask the quadrature rule to "normalize" the weights $w_i$ with $s(q_i)$:
 *   

 * 
 * \f[ \int_K f(x) s(x) dx = \int_K g(x) dx = \sum_{i=1}^N
 * \frac{w_i}{s(q_i)} g(q_i) \f]
 *   

 * 
 * We use this second option, through the @p factor_out_singularity
 * parameter of both QGaussLogR and QGaussOneOverR.
 *   

 * 
 * These integrals are somewhat delicate, especially in two dimensions, due
 * to the transformation from the real to the reference cell, where the
 * variable of integration is scaled with the determinant of the
 * transformation.
 *   

 * 
 * In two dimensions this process does not result only in a factor appearing
 * as a constant factor on the entire integral, but also on an additional
 * integral altogether that needs to be evaluated:
 *   

 * 
 * \f[ \int_0^1 f(x)\ln(x/\alpha) dx = \int_0^1 f(x)\ln(x) dx - \int_0^1
 * f(x) \ln(\alpha) dx.  \f]
 *   

 * 
 * This process is taken care of by the constructor of the QGaussLogR class,
 * which adds additional quadrature points and weights to take into
 * consideration also the second part of the integral.
 *   

 * 
 * A similar reasoning should be done in the three dimensional case, since
 * the singular quadrature is tailored on the inverse of the radius $r$ in
 * the reference cell, while our singular function lives in real space,
 * however in the three dimensional case everything is simpler because the
 * singularity scales linearly with the determinant of the
 * transformation. This allows us to build the singular two dimensional
 * quadrature rules only once and, reuse them over all cells.
 *   

 * 
 * In the one dimensional singular integration this is not possible, since
 * we need to know the scaling parameter for the quadrature, which is not
 * known a priori. Here, the quadrature rule itself depends also on the size
 * of the current cell. For this reason, it is necessary to create a new
 * quadrature for each singular integration.
 *   

 * 
 * The different quadrature rules are built inside the
 * get_singular_quadrature, which is specialized for dim=2 and dim=3, and
 * they are retrieved inside the assemble_system function. The index given
 * as an argument is the index of the unit support point where the
 * singularity is located.
 * 

 * 
 * 
 * @code
 *   template <>
 *   const Quadrature<2> &BEMProblem<3>::get_singular_quadrature(
 *     const DoFHandler<2, 3>::active_cell_iterator &,
 *     const unsigned int index) const
 *   {
 *     Assert(index < fe.n_dofs_per_cell(),
 *            ExcIndexRange(0, fe.n_dofs_per_cell(), index));
 * 
 *     static std::vector<QGaussOneOverR<2>> quadratures;
 *     if (quadratures.size() == 0)
 *       for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
 *         quadratures.emplace_back(singular_quadrature_order,
 *                                  fe.get_unit_support_points()[i],
 *                                  true);
 *     return quadratures[index];
 *   }
 * 
 * 
 *   template <>
 *   const Quadrature<1> &BEMProblem<2>::get_singular_quadrature(
 *     const DoFHandler<1, 2>::active_cell_iterator &cell,
 *     const unsigned int                            index) const
 *   {
 *     Assert(index < fe.n_dofs_per_cell(),
 *            ExcIndexRange(0, fe.n_dofs_per_cell(), index));
 * 
 *     static Quadrature<1> *q_pointer = nullptr;
 *     if (q_pointer)
 *       delete q_pointer;
 * 
 *     q_pointer = new QGaussLogR<1>(singular_quadrature_order,
 *                                   fe.get_unit_support_points()[index],
 *                                   1. / cell->measure(),
 *                                   true);
 *     return (*q_pointer);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemcompute_exterior_solution"></a> 
 * <h4>BEMProblem::compute_exterior_solution</h4>
 * 

 * 
 * We'd like to also know something about the value of the potential $\phi$
 * in the exterior domain: after all our motivation to consider the boundary
 * integral problem was that we wanted to know the velocity in the exterior
 * domain!
 *   

 * 
 * To this end, let us assume here that the boundary element domain is
 * contained in the box $[-2,2]^{\text{dim}}$, and we extrapolate the actual
 * solution inside this box using the convolution with the fundamental
 * solution. The formula for this is given in the introduction.
 *   

 * 
 * The reconstruction of the solution in the entire space is done on a
 * continuous finite element grid of dimension dim. These are the usual
 * ones, and we don't comment any further on them. At the end of the
 * function, we output this exterior solution in, again, much the usual way.
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::compute_exterior_solution()
 *   {
 *     Triangulation<dim> external_tria;
 *     GridGenerator::hyper_cube(external_tria, -2, 2);
 * 
 *     FE_Q<dim>       external_fe(1);
 *     DoFHandler<dim> external_dh(external_tria);
 *     Vector<double>  external_phi;
 * 
 *     external_tria.refine_global(external_refinement);
 *     external_dh.distribute_dofs(external_fe);
 *     external_phi.reinit(external_dh.n_dofs());
 * 
 *     FEValues<dim - 1, dim> fe_v(mapping,
 *                                 fe,
 *                                 *quadrature,
 *                                 update_values | update_normal_vectors |
 *                                   update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int n_q_points = fe_v.n_quadrature_points;
 * 
 *     std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell());
 * 
 *     std::vector<double>         local_phi(n_q_points);
 *     std::vector<double>         normal_wind(n_q_points);
 *     std::vector<Vector<double>> local_wind(n_q_points, Vector<double>(dim));
 * 
 *     std::vector<Point<dim>> external_support_points(external_dh.n_dofs());
 *     DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping,
 *                                               external_dh,
 *                                               external_support_points);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_v.reinit(cell);
 * 
 *         const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
 *         const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();
 * 
 *         cell->get_dof_indices(dofs);
 *         fe_v.get_function_values(phi, local_phi);
 * 
 *         wind.vector_value_list(q_points, local_wind);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             normal_wind[q] = 0;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               normal_wind[q] += normals[q][d] * local_wind[q](d);
 *           }
 * 
 *         for (unsigned int i = 0; i < external_dh.n_dofs(); ++i)
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             {
 *               const Tensor<1, dim> R = q_points[q] - external_support_points[i];
 * 
 *               external_phi(i) +=
 *                 ((LaplaceKernel::single_layer(R) * normal_wind[q] +
 *                   (LaplaceKernel::double_layer(R) * normals[q]) *
 *                     local_phi[q]) *
 *                  fe_v.JxW(q));
 *             }
 *       }
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(external_dh);
 *     data_out.add_data_vector(external_phi, "external_phi");
 *     data_out.build_patches();
 * 
 *     const std::string filename = std::to_string(dim) + "d_external.vtk";
 *     std::ofstream     file(filename);
 * 
 *     data_out.write_vtk(file);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemoutput_results"></a> 
 * <h4>BEMProblem::output_results</h4>
 * 

 * 
 * Outputting the results of our computations is a rather mechanical
 * tasks. All the components of this function have been discussed before.
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::output_results(const unsigned int cycle)
 *   {
 *     DataOut<dim - 1, DoFHandler<dim - 1, dim>> dataout;
 * 
 *     dataout.attach_dof_handler(dof_handler);
 *     dataout.add_data_vector(
 *       phi, "phi", DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);
 *     dataout.add_data_vector(
 *       alpha,
 *       "alpha",
 *       DataOut<dim - 1, DoFHandler<dim - 1, dim>>::type_dof_data);
 *     dataout.build_patches(
 *       mapping,
 *       mapping.get_degree(),
 *       DataOut<dim - 1, DoFHandler<dim - 1, dim>>::curved_inner_cells);
 * 
 *     const std::string filename = std::to_string(dim) + "d_boundary_solution_" +
 *                                  std::to_string(cycle) + ".vtk";
 *     std::ofstream file(filename);
 * 
 *     dataout.write_vtk(file);
 * 
 *     if (cycle == n_cycles - 1)
 *       {
 *         convergence_table.set_precision("L2(phi)", 3);
 *         convergence_table.set_precision("Linfty(alpha)", 3);
 * 
 *         convergence_table.set_scientific("L2(phi)", true);
 *         convergence_table.set_scientific("Linfty(alpha)", true);
 * 
 *         convergence_table.evaluate_convergence_rates(
 *           "L2(phi)", ConvergenceTable::reduction_rate_log2);
 *         convergence_table.evaluate_convergence_rates(
 *           "Linfty(alpha)", ConvergenceTable::reduction_rate_log2);
 *         deallog << std::endl;
 *         convergence_table.write_text(std::cout);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemrun"></a> 
 * <h4>BEMProblem::run</h4>
 * 

 * 
 * This is the main function. It should be self explanatory in its
 * briefness:
 * 
 * @code
 *   template <int dim>
 *   void BEMProblem<dim>::run()
 *   {
 *     read_parameters("parameters.prm");
 * 
 *     if (run_in_this_dimension == false)
 *       {
 *         deallog << "Run in dimension " << dim
 *                 << " explicitly disabled in parameter file. " << std::endl;
 *         return;
 *       }
 * 
 *     read_domain();
 * 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         refine_and_resize();
 *         assemble_system();
 *         solve_system();
 *         compute_errors(cycle);
 *         output_results(cycle);
 *       }
 * 
 *     if (extend_solution == true)
 *       compute_exterior_solution();
 *   }
 * } // namespace Step34
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * This is the main function of this program. It is exactly like all previous
 * tutorial programs:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step34;
 * 
 *       const unsigned int degree         = 1;
 *       const unsigned int mapping_degree = 1;
 * 
 *       deallog.depth_console(3);
 *       BEMProblem<2> laplace_problem_2d(degree, mapping_degree);
 *       laplace_problem_2d.run();
 * 
 *       BEMProblem<3> laplace_problem_3d(degree, mapping_degree);
 *       laplace_problem_3d.run();
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


We ran the program using the following <code>parameters.prm</code> file (which
can also be found in the directory in which all the other source files are):
@verbatim
# Listing of Parameters
# ---------------------
set Extend solution on the -2,2 box = true
set External refinement             = 5
set Number of cycles                = 4
set Run 2d simulation               = true
set Run 3d simulation               = true


subsection Exact solution 2d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =

  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = x+y   # default: 0

  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,t
end


subsection Exact solution 3d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =

  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = .5*(x+y+z)   # default: 0

  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,z,t
end


subsection Quadrature rules
  set Quadrature order          = 4
  set Quadrature type           = gauss
  set Singular quadrature order = 5
end


subsection Solver
  set Log frequency = 1
  set Log history   = false
  set Log result    = true
  set Max steps     = 100
  set Tolerance     = 1.e-10
end


subsection Wind function 2d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =

  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = 1; 1  # default: 0; 0

  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,t
end


subsection Wind function 3d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =

  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = 1; 1; 1 # default: 0; 0; 0

  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,z,t
end
@endverbatim

When we run the program, the following is printed on screen:
@verbatim
DEAL::
DEAL::Parsing parameter file parameters.prm
DEAL::for a 2 dimensional simulation.
DEAL:GMRES::Starting value 2.21576
DEAL:GMRES::Convergence step 1 value 2.37635e-13
DEAL::Cycle 0:
DEAL::   Number of active cells:       20
DEAL::   Number of degrees of freedom: 20
DEAL:GMRES::Starting value 3.15543
DEAL:GMRES::Convergence step 1 value 2.89310e-13
DEAL::Cycle 1:
DEAL::   Number of active cells:       40
DEAL::   Number of degrees of freedom: 40
DEAL:GMRES::Starting value 4.46977
DEAL:GMRES::Convergence step 1 value 3.11815e-13
DEAL::Cycle 2:
DEAL::   Number of active cells:       80
DEAL::   Number of degrees of freedom: 80
DEAL:GMRES::Starting value 6.32373
DEAL:GMRES::Convergence step 1 value 3.22474e-13
DEAL::Cycle 3:
DEAL::   Number of active cells:       160
DEAL::   Number of degrees of freedom: 160
DEAL::
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    20   20 4.465e-02    - 5.000e-02    -
    1    40   40 1.081e-02 2.05 2.500e-02 1.00
    2    80   80 2.644e-03 2.03 1.250e-02 1.00
    3   160  160 6.529e-04 2.02 6.250e-03 1.00
DEAL::
DEAL::Parsing parameter file parameters.prm
DEAL::for a 3 dimensional simulation.
DEAL:GMRES::Starting value 2.84666
DEAL:GMRES::Convergence step 3 value 8.68638e-18
DEAL::Cycle 0:
DEAL::   Number of active cells:       24
DEAL::   Number of degrees of freedom: 26
DEAL:GMRES::Starting value 6.34288
DEAL:GMRES::Convergence step 5 value 1.38740e-11
DEAL::Cycle 1:
DEAL::   Number of active cells:       96
DEAL::   Number of degrees of freedom: 98
DEAL:GMRES::Starting value 12.9780
DEAL:GMRES::Convergence step 5 value 3.29225e-11
DEAL::Cycle 2:
DEAL::   Number of active cells:       384
DEAL::   Number of degrees of freedom: 386
DEAL:GMRES::Starting value 26.0874
DEAL:GMRES::Convergence step 6 value 1.47271e-12
DEAL::Cycle 3:
DEAL::   Number of active cells:       1536
DEAL::   Number of degrees of freedom: 1538
DEAL::
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    24   26 3.437e-01    - 2.327e-01    -
    1    96   98 9.794e-02 1.81 1.239e-01 0.91
    2   384  386 2.417e-02 2.02 6.319e-02 0.97
    3  1536 1538 5.876e-03 2.04 3.176e-02 0.99
@endverbatim

As we can see from the convergence table in 2d, if we choose
quadrature formulas which are accurate enough, then the error we
obtain for $\alpha(\mathbf{x})$ should be exactly the inverse of the
number of elements. The approximation of the circle with N segments of
equal size generates a regular polygon with N faces, whose angles are
exactly $\pi-\frac {2\pi}{N}$, therefore the error we commit should be
exactly $\frac 12 - (\frac 12 -\frac 1N) = \frac 1N$. In fact this is
a very good indicator that we are performing the singular integrals in
an appropriate manner.

The error in the approximation of the potential $\phi$ is largely due
to approximation of the domain. A much better approximation could be
obtained by using higher order mappings.

If we modify the main() function, setting fe_degree and mapping_degree
to two, and raise the order of the quadrature formulas  in
the parameter file, we obtain the following convergence table for the
two dimensional simulation

@verbatim
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    20   40 5.414e-05    - 2.306e-04    -
    1    40   80 3.623e-06 3.90 1.737e-05 3.73
    2    80  160 2.690e-07 3.75 1.253e-05 0.47
    3   160  320 2.916e-08 3.21 7.670e-06 0.71
@endverbatim

and

@verbatim
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    24   98 3.770e-03    - 8.956e-03    -
    1    96  386 1.804e-04 4.39 1.182e-03 2.92
    2   384 1538 9.557e-06 4.24 1.499e-04 2.98
    3  1536 6146 6.617e-07 3.85 1.892e-05 2.99
@endverbatim

for the three dimensional case. As we can see, convergence results are
much better with higher order mapping, mainly due to a better
resolution of the curved geometry. Notice that, given the same number
of degrees of freedom, for example in step 3 of the Q1 case and step 2
of Q2 case in the three dimensional simulation, the error is roughly
three orders of magnitude lower.

The result of running these computations is a bunch of output files that we
can pass to our visualization program of choice.
The output files are of two kind: the potential on the boundary
element surface, and the potential extended to the outer and inner
domain. The combination of the two for the two dimensional case looks
like

<img src="https://www.dealii.org/images/steps/developer/step-34_2d.png" alt="">

while in three dimensions we show first the potential on the surface,
together with a contour plot,

<img src="https://www.dealii.org/images/steps/developer/step-34_3d.png" alt="">

and then the external contour plot of the potential, with opacity set to 25%:

<img src="https://www.dealii.org/images/steps/developer/step-34_3d-2.png" alt="">


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


This is the first tutorial program that considers solving equations defined on
surfaces embedded in higher dimensional spaces. But the equation discussed
here was relatively simple because it only involved an integral operator, not
derivatives which are more difficult to define on the surface. The step-38
tutorial program considers such problems and provides the necessary tools.

From a practical perspective, the Boundary Element Method (BEM) used
here suffers from two bottlenecks. The first is that assembling the
matrix has a cost that is *quadratic* in the number of unknowns, that
is ${\cal O}(N^2)$ where $N$ is the total number of unknowns. This can
be seen by looking at the `assemble_system()` function, which has this
structure:
@code
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        ...

        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
          ...
@endcode
Here, the first loop walks over all cells (one factor of $N$) whereas
the inner loop contributes another factor of $N$.

This has to be contrasted with the finite element method for *local*
differential operators: There, we loop over all cells (one factor of
$N$) and on each cell do an amount of work that is independent of how
many cells or unknowns there are. This clearly presents a
bottleneck.

The second bottleneck is that the system matrix is dense (i.e., is of
type FullMatrix) because every degree of freedom couples with every
other degree of freedom. As pointed out above, just *computing* this
matrix with its $N^2$ nonzero entries necessarily requires at least
${\cal O}(N^2)$ operations, but it's worth pointing out that it also
costs this many operations to just do one matrix-vector product. If
the GMRES method used to solve the linear system requires a number of
iterations that grows with the size of the problem, as is typically
the case, then solving the linear system will require a number of
operations that grows even faster than just ${\cal O}(N^2)$.

"Real" boundary element methods address these issues by strategies
that determine which entries of the matrix will be small and can
consequently be neglected (at the cost of introducing an additional
error, of course). This is possible by recognizing that the matrix
entries decay with the (physical) distance between the locations where
degrees of freedom $i$ and $j$ are defined. This can be exploited in
methods such as the Fast Multipole Method (FMM) that control which
matrix entries must be stored and computed to achieve a certain
accuracy, and -- if done right -- result in methods in which both
assembly and solution of the linear system requires less than
${\cal O}(N^2)$ operations.

Implementing these methods clearly presents opportunities to extend
the current program.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-34.cc"
*/
