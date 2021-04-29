/**
@page step_23 The step-23 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Timediscretization">Time discretization</a>
      <ul>
        <li><a href="#MethodoflinesorRothesmethod">Method of lines or Rothe's method?</a>
        <li><a href="#Rothesmethod">Rothe's method!</a>
      </ul>
        <li><a href="#Spacediscretization">Space discretization</a>
        <li><a href="#Energyconservation">Energy conservation</a>
        <li><a href="#WhoareCourantFriedrichsandLewy">Who are Courant, Friedrichs, and Lewy?</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeWaveEquationcodeclass">The <code>WaveEquation</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeWaveEquationcodeclass">Implementation of the <code>WaveEquation</code> class</a>
      <ul>
        <li><a href="#WaveEquationsetup_system">WaveEquation::setup_system</a>
        <li><a href="#WaveEquationsolve_uandWaveEquationsolve_v">WaveEquation::solve_u and WaveEquation::solve_v</a>
        <li><a href="#WaveEquationoutput_results">WaveEquation::output_results</a>
        <li><a href="#WaveEquationrun">WaveEquation::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{28}

This is the first of a number of tutorial programs that will finally
cover "real" time-dependent problems, not the slightly odd form of time
dependence found in step-18 or the DAE model of step-21. In particular, this program introduces
the wave equation in a bounded domain. Later, step-24
will consider an example of absorbing boundary conditions, and @ref
step_25 "step-25" a kind of nonlinear wave equation producing
solutions called solitons.

The wave equation in its prototypical form reads as follows: find
$u(x,t), x\in\Omega, t\in[0,T]$ that satisfies
@f{eqnarray*}
	\frac{\partial^2 u}{\partial t^2}
	-
	\Delta u &=& f
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	u(x,t) &=& g
	\qquad
	\textrm{on}\ \partial\Omega\times [0,T],
\\
	u(x,0) &=& u_0(x)
	\qquad
	\textrm{in}\ \Omega,
\\
	\frac{\partial u(x,0)}{\partial t} &=& u_1(x)
	\qquad
	\textrm{in}\ \Omega.
@f}
Note that since this is an equation with second-order time
derivatives, we need to pose two initial conditions, one for the value
and one for the time derivative of the solution.

Physically, the equation describes the motion of an elastic medium. In
2-d, one can think of how a membrane moves if subjected to a
force. The Dirichlet boundary conditions above indicate that the
membrane is clamped at the boundary at a height $g(x,t)$ (this height
might be moving as well &mdash; think of people holding a blanket and
shaking it up and down). The first initial condition equals the
initial deflection of the membrane, whereas the second one gives its
velocity. For example, one could think of pushing the membrane down
with a finger and then letting it go at $t=0$ (nonzero deflection but
zero initial velocity), or hitting it with a hammer at $t=0$ (zero
deflection but nonzero velocity). Both cases would induce motion in
the membrane.


<a name="Timediscretization"></a><h3>Time discretization</h3>


<a name="MethodoflinesorRothesmethod"></a><h4>Method of lines or Rothe's method?</h4>

There is a long-standing debate in the numerical analysis community
over whether a discretization of time dependent equations should
involve first discretizing the time variable leading to a stationary
PDE at each time step that is then solved using standard finite
element techniques (this is called the Rothe method), or whether
one should first discretize the spatial variables, leading to a large
system of ordinary differential equations that can then be handled by
one of the usual ODE solvers (this is called the method of lines).

Both of these methods have advantages and disadvantages.
Traditionally, people have preferred the method of lines, since it
allows to use the very well developed machinery of high-order ODE
solvers available for the rather stiff ODEs resulting from this
approach, including step length control and estimation of the temporal
error.

On the other hand, Rothe's method becomes awkward when using
higher-order time stepping method, since one then has to write down a
PDE that couples the solution of the present time step not only with
that at the previous time step, but possibly also even earlier
solutions, leading to a significant number of terms.

For these reasons, the method of lines was the method of choice for a
long time. However, it has one big drawback: if we discretize the
spatial variable first, leading to a large ODE system, we have to
choose a mesh once and for all. If we are willing to do this, then
this is a legitimate and probably superior approach.

If, on the other hand, we are looking at the wave equation and many
other time dependent problems, we find that the character of a
solution changes as time progresses. For example, for the wave
equation, we may have a single wave travelling through the domain,
where the solution is smooth or even constant in front of and behind
the wave &mdash; adaptivity would be really useful for such cases, but the
key is that the area where we need to refine the mesh changes from
time step to time step!

If we intend to go that way, i.e. choose a different mesh for each
time step (or set of time steps), then the method of lines is not
appropriate any more: instead of getting one ODE system with a number
of variables equal to the number of unknowns in the finite element
mesh, our number of unknowns now changes all the time, a fact that
standard ODE solvers are certainly not prepared to deal with at
all. On the other hand, for the Rothe method, we just get a PDE for
each time step that we may choose to discretize independently of the
mesh used for the previous time step; this approach is not without
perils and difficulties, but at least is a sensible and well-defined
procedure.

For all these reasons, for the present program, we choose to use the
Rothe method for discretization, i.e. we first discretize in time and
then in space. We will not actually use adaptive meshes at all, since
this involves a large amount of additional code, but we will comment
on this some more in the <a href="#Results">results section below</a>.


<a name="Rothesmethod"></a><h4>Rothe's method!</h4>


Given these considerations, here is how we will proceed: let us first
define a simple time stepping method for this second order problem,
and then in a second step do the spatial discretization, i.e. we will
follow Rothe's approach.

For the first step, let us take a little detour first: in order to
discretize a second time derivative, we can either discretize it
directly, or we can introduce an additional variable and transform the
system into a first order system. In many cases, this turns out to be
equivalent, but dealing with first order systems is often simpler. To
this end, let us introduce
@f[
	v = \frac{\partial u}{\partial t},
@f]
and call this variable the <i>velocity</i> for obvious reasons. We can
then reformulate the original wave equation as follows:
@f{eqnarray*}
	\frac{\partial u}{\partial t}
	-
	v
	&=& 0
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	\frac{\partial v}{\partial t}
	-
	\Delta u &=& f
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	u(x,t) &=& g
	\qquad
	\textrm{on}\ \partial\Omega\times [0,T],
\\
	u(x,0) &=& u_0(x)
	\qquad
	\textrm{in}\ \Omega,
\\
	v(x,0) &=& u_1(x)
	\qquad
	\textrm{in}\ \Omega.
@f}
The advantage of this formulation is that it now only contains first
time derivatives for both variables, for which it is simple to write
down time stepping schemes. Note that we do not have boundary
conditions for $v$ at first. However, we could enforce $v=\frac{\partial
g}{\partial t}$ on the boundary. It turns out in numerical examples that this
is actually necessary: without doing so the solution doesn't look particularly
wrong, but the Crank-Nicolson scheme does not conserve energy if one doesn't
enforce these boundary conditions.

With this formulation, let us introduce the following time
discretization where a superscript $n$ indicates the number of a time
step and $k=t_n-t_{n-1}$ is the length of the present time step:
\f{eqnarray*}
  \frac{u^n - u^{n-1}}{k}
  - \left[\theta v^n + (1-\theta) v^{n-1}\right] &=& 0,
  \\
  \frac{v^n - v^{n-1}}{k}
  - \Delta\left[\theta u^n + (1-\theta) u^{n-1}\right]
  &=& \theta f^n + (1-\theta) f^{n-1}.
\f}
Note how we introduced a parameter $\theta$ here. If we chose
$\theta=0$, for example, the first equation would reduce to
$\frac{u^n - u^{n-1}}{k}  - v^{n-1} = 0$, which is well-known as the
forward or explicit Euler method. On the other hand, if we set
$\theta=1$, then we would get
$\frac{u^n - u^{n-1}}{k}  - v^n = 0$, which corresponds to the
backward or implicit Euler method. Both these methods are first order
accurate methods. They are simple to implement, but they are not
really very accurate.

The third case would be to choose $\theta=\frac 12$. The first of the
equations above would then read $\frac{u^n - u^{n-1}}{k}
- \frac 12 \left[v^n + v^{n-1}\right] = 0$. This method is known as
the Crank-Nicolson method and has the advantage that it is second
order accurate. In addition, it has the nice property that it
preserves the energy in the solution (physically, the energy is the
sum of the kinetic energy of the particles in the membrane plus the
potential energy present due to the fact that it is locally stretched;
this quantity is a conserved one in the continuous equation, but most
time stepping schemes do not conserve it after time
discretization). Since $v^n$ also appears in the equation for $u^n$,
the Crank-Nicolson scheme is also implicit.

In the program, we will leave $\theta$ as a parameter, so that it will
be easy to play with it. The results section will show some numerical
evidence comparing the different schemes.

The equations above (called the <i>semidiscretized</i> equations
because we have only discretized the time, but not space), can be
simplified a bit by eliminating $v^n$ from the first equation and
rearranging terms. We then get
\f{eqnarray*}
  \left[ 1-k^2\theta^2\Delta \right] u^n &=&
  	 \left[ 1+k^2\theta(1-\theta)\Delta\right] u^{n-1} + k v^{n-1}
   	 + k^2\theta\left[\theta f^n + (1-\theta) f^{n-1}\right],\\
   v^n &=& v^{n-1} + k\Delta\left[ \theta u^n + (1-\theta) u^{n-1}\right]
   + k\left[\theta f^n + (1-\theta) f^{n-1}\right].
\f}
In this form, we see that if we are given the solution
$u^{n-1},v^{n-1}$ of the previous timestep, that we can then solve for
the variables $u^n,v^n$ separately, i.e. one at a time. This is
convenient. In addition, we recognize that the operator in the first
equation is positive definite, and the second equation looks
particularly simple.


<a name="Spacediscretization"></a><h3>Space discretization</h3>


We have now derived equations that relate the approximate
(semi-discrete) solution $u^n(x)$ and its time derivative $v^n(x)$ at
time $t_n$ with the solutions $u^{n-1}(x),v^{n-1}(x)$ of the previous
time step at $t_{n-1}$. The next step is to also discretize the
spatial variable using the usual finite element methodology. To this
end, we multiply each equation with a test function, integrate over
the entire domain, and integrate by parts where necessary. This leads
to
\f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla u^{n-1},\nabla \varphi)
  +
  k(v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right],
  \\
  (v^n,\varphi)
   &=&
   (v^{n-1},\varphi)
    -
    k\left[ \theta (\nabla u^n,\nabla\varphi) +
    (1-\theta) (\nabla u^{n-1},\nabla \varphi)\right]
  + k
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right].
\f}

It is then customary to approximate $u^n(x) \approx u^n_h(x) = \sum_i
U_i^n\phi_i^n(x)$, where $\phi_i^n(x)$ are the shape functions used
for the discretization of the $n$-th time step and $U_i^n$ are the
unknown nodal values of the solution. Similarly, $v^n(x) \approx
v^n_h(x) = \sum_i V_i^n\phi_i^n(x)$. Finally, we have the solutions of
the previous time step, $u^{n-1}(x) \approx u^{n-1}_h(x) = \sum_i
U_i^{n-1}\phi_i^{n-1}(x)$ and $v^{n-1}(x) \approx v^{n-1}_h(x) = \sum_i
V_i^{n-1}\phi_i^{n-1}(x)$. Note that since the solution of the previous
time step has already been computed by the time we get to time step
$n$, $U^{n-1},V^{n-1}$ are known. Furthermore, note that the solutions
of the previous step may have been computed on a different mesh, so
we have to use shape functions $\phi^{n-1}_i(x)$.

If we plug these expansions into above equations and test with the
test functions from the present mesh, we get the following linear
system:
\f{eqnarray*}
  (M^n + k^2\theta^2 A^n)U^n &=&
  M^{n,n-1}U^{n-1} - k^2\theta(1-\theta) A^{n,n-1}U^{n-1}
  +
  kM^{n,n-1}V^{n-1}
  + k^2\theta
  \left[
  \theta F^n + (1-\theta) F^{n-1}
  \right],
  \\
  M^nV^n
   &=&
   M^{n,n-1}V^{n-1}
    -
    k\left[ \theta A^n U^n +
    (1-\theta) A^{n,n-1} U^{n-1}\right]
   + k
  \left[
  \theta F^n + (1-\theta) F^{n-1}
  \right],
\f}
where
@f{eqnarray*}
	M^n_{ij} &=& (\phi_i^n, \phi_j^n),
	\\
	A^n_{ij} &=& (\nabla\phi_i^n, \nabla\phi_j^n),
	\\
	M^{n,n-1}_{ij} &=& (\phi_i^n, \phi_j^{n-1}),
	\\
	A^{n,n-1}_{ij} &=& (\nabla\phi_i^n, \nabla\phi_j^{n-1}),
	\\
	F^n_{i} &=& (f^n,\phi_i^n),
	\\
	F^{n-1}_{i} &=& (f^{n-1},\phi_i^n).
@f}

If we solve these two equations, we can move the solution one step
forward and go on to the next time step.

It is worth noting that if we choose the same mesh on each time step
(as we will in fact do in the program below), then we have the same
shape functions on time step $n$ and $n-1$,
i.e. $\phi^n_i=\phi_i^{n-1}=\phi_i$. Consequently, we get
$M^n=M^{n,n-1}=M$ and $A^n=A^{n,n-1}=A$. On the other hand, if we had
used different shape functions, then we would have to compute
integrals that contain shape functions defined on two meshes. This is a
somewhat messy process that we omit here, but that is treated in some
detail in step-28.

Under these conditions (i.e. a mesh that doesn't change), one can optimize the
solution procedure a bit by basically eliminating the solution of the second
linear system. We will discuss this in the introduction of the @ref step_25
"step-25" program.

<a name="Energyconservation"></a><h3>Energy conservation</h3>


One way to compare the quality of a time stepping scheme is to see whether the
numerical approximation preserves conservation properties of the continuous
equation. For the wave equation, the natural quantity to look at is the
energy. By multiplying the wave equation by $u_t$, integrating over $\Omega$,
and integrating by parts where necessary, we find that
@f[
	\frac{d}{d t}
	\left[\frac 12 \int_\Omega \left(\frac{\partial u}{\partial
	t}\right)^2 + (\nabla u)^2 \; dx\right]
	=
	\int_\Omega f \frac{\partial u}{\partial t} \; dx
	+
	\int_{\partial\Omega} n\cdot\nabla u
	\frac{\partial g}{\partial t} \; dx.
@f]
By consequence, in absence of body forces and constant boundary values, we get
that
@f[
	E(t) = \frac 12 \int_\Omega \left(\frac{\partial u}{\partial
	t}\right)^2 + (\nabla u)^2 \; dx
@f]
is a conserved quantity, i.e. one that doesn't change with time. We
will compute this quantity after each time
step. It is straightforward to see that if we replace $u$ by its finite
element approximation, and $\frac{\partial u}{\partial t}$ by the finite
element approximation of the velocity $v$, then
@f[
	E(t_n) = \frac 12 \left<V^n, M^n V^n\right>
	+
	\frac 12 \left<U^n, A^n U^n\right>.
@f]
As we will see in the results section, the Crank-Nicolson scheme does indeed
conserve the energy, whereas neither the forward nor the backward Euler scheme
do.


<a name="WhoareCourantFriedrichsandLewy"></a><h3>Who are Courant, Friedrichs, and Lewy?</h3>


One of the reasons why the wave equation is nasty to solve numerically is that
explicit time discretizations are only stable if the time step is small
enough. In particular, it is coupled to the spatial mesh width $h$. For the
lowest order discretization we use here, the relationship reads
@f[
	k\le \frac hc
@f]
where $c$ is the wave speed, which in our formulation of the wave equation has
been normalized to one. Consequently, unless we use the implicit schemes with
$\theta>0$, our solutions will not be numerically stable if we violate this
restriction. Implicit schemes do not have this restriction for stability, but
they become inaccurate if the time step is too large.

This condition was first recognized by Courant, Friedrichs, and Lewy &mdash;
in 1928, long before computers became available for numerical
computations! (This result appeared in the German language article
R. Courant, K. Friedrichs and H. Lewy: <i>&Uuml;ber die partiellen
Differenzengleichungen der mathematischen Physik</i>, Mathematische
Annalen, vol. 100, no. 1, pages 32-74, 1928.)
This condition on the time step is most frequently just referred
to as the <i>CFL</i> condition. Intuitively, the CFL condition says
that the time step must not be larger than the time it takes a wave to
cross a single cell.

In the program, we will refine the square
$[-1,1]^2$ seven times uniformly, giving a mesh size of $h=\frac 1{64}$, which
is what we set the time step to. The fact that we set the time step and mesh
size individually in two different places is error prone: it is too easy to
refine the mesh once more but forget to also adjust the time step. @ref
step_24 "step-24" shows a better way how to keep these things in sync.


<a name="Thetestcase"></a><h3>The test case</h3>


Although the program has all the hooks to deal with nonzero initial and
boundary conditions and body forces, we take a simple case where the domain is
a square $[-1,1]^2$ and
@f{eqnarray*}
	f &=& 0,
	\\
	u_0 &=& 0,
	\\
	u_1 &=& 0,
	\\
	g &=& \left\{\begin{matrix}\sin (4\pi t)
	&\qquad& \text{for }\ t\le \frac 12, x=-1, -\frac 13<y<\frac 13
	\\
	 0
	&&\text{otherwise}
	\end{matrix}
	\right.
@f}
This corresponds to a membrane initially at rest and clamped all around, where
someone is waving a part of the clamped boundary once up and down, thereby
shooting a wave into the domain.
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
 * We start with the usual assortment of include files that we've seen in so
 * many of the previous tests:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * Here are the only three include files of some new interest: The first one
 * is already used, for example, for the
 * VectorTools::interpolate_boundary_values and
 * MatrixTools::apply_boundary_values functions. However, we here use another
 * function in that class, VectorTools::project to compute our initial values
 * as the $L^2$ projection of the continuous initial values. Furthermore, we
 * use VectorTools::create_right_hand_side to generate the integrals
 * $(f^n,\phi^n_i)$. These were previously always generated by hand in
 * <code>assemble_system</code> or similar functions in application
 * code. However, we're too lazy to do that here, so simply use a library
 * function:
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * In a very similar vein, we are also too lazy to write the code to assemble
 * mass and Laplace matrices, although it would have only taken copying the
 * relevant code from any number of previous tutorial programs. Rather, we
 * want to focus on the things that are truly new to this program and
 * therefore use the MatrixCreator::create_mass_matrix and
 * MatrixCreator::create_laplace_matrix functions. They are declared here:
 * 
 * @code
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * @endcode
 * 
 * Finally, here is an include file that contains all sorts of tool functions
 * that one sometimes needs. In particular, we need the
 * Utilities::int_to_string class that, given an integer argument, returns a
 * string representation of it. It is particularly useful since it allows for
 * a second parameter indicating the number of digits to which we want the
 * result padded with leading zeros. We will use this to write output files
 * that have the form <code>solution-XXX.vtu</code> where <code>XXX</code>
 * denotes the number of the time step and always consists of three digits
 * even if we are still in the single or double digit time steps.
 * 
 * @code
 * #include <deal.II/base/utilities.h>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step23
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeWaveEquationcodeclass"></a> 
 * <h3>The <code>WaveEquation</code> class</h3>
 * 

 * 
 * Next comes the declaration of the main class. It's public interface of
 * functions is like in most of the other tutorial programs. Worth
 * mentioning is that we now have to store four matrices instead of one: the
 * mass matrix $M$, the Laplace matrix $A$, the matrix $M+k^2\theta^2A$ used
 * for solving for $U^n$, and a copy of the mass matrix with boundary
 * conditions applied used for solving for $V^n$. Note that it is a bit
 * wasteful to have an additional copy of the mass matrix around. We will
 * discuss strategies for how to avoid this in the section on possible
 * improvements.
 *   

 * 
 * Likewise, we need solution vectors for $U^n,V^n$ as well as for the
 * corresponding vectors at the previous time step, $U^{n-1},V^{n-1}$. The
 * <code>system_rhs</code> will be used for whatever right hand side vector
 * we have when solving one of the two linear systems in each time
 * step. These will be solved in the two functions <code>solve_u</code> and
 * <code>solve_v</code>.
 *   

 * 
 * Finally, the variable <code>theta</code> is used to indicate the
 * parameter $\theta$ that is used to define which time stepping scheme to
 * use, as explained in the introduction. The rest is self-explanatory.
 * 
 * @code
 *   template <int dim>
 *   class WaveEquation
 *   {
 *   public:
 *     WaveEquation();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void solve_u();
 *     void solve_v();
 *     void output_results() const;
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> mass_matrix;
 *     SparseMatrix<double> laplace_matrix;
 *     SparseMatrix<double> matrix_u;
 *     SparseMatrix<double> matrix_v;
 * 
 *     Vector<double> solution_u, solution_v;
 *     Vector<double> old_solution_u, old_solution_v;
 *     Vector<double> system_rhs;
 * 
 *     double       time_step;
 *     double       time;
 *     unsigned int timestep_number;
 *     const double theta;
 *   };
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
 * Before we go on filling in the details of the main class, let us define
 * the equation data corresponding to the problem, i.e. initial and boundary
 * values for both the solution $u$ and its time derivative $v$, as well as
 * a right hand side class. We do so using classes derived from the Function
 * class template that has been used many times before, so the following
 * should not be a surprise.
 *   

 * 
 * Let's start with initial values and choose zero for both the value $u$ as
 * well as its time derivative, the velocity $v$:
 * 
 * @code
 *   template <int dim>
 *   class InitialValuesU : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   class InitialValuesV : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * Secondly, we have the right hand side forcing term. Boring as we are, we
 * choose zero here as well:
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * Finally, we have boundary values for $u$ and $v$. They are as described
 * in the introduction, one being the time derivative of the other:
 * 
 * @code
 *   template <int dim>
 *   class BoundaryValuesU : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 * 
 *       if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) &&
 *           (p[1] > -1. / 3))
 *         return std::sin(this->get_time() * 4 * numbers::PI);
 *       else
 *         return 0;
 *     }
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   class BoundaryValuesV : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 * 
 *       if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) &&
 *           (p[1] > -1. / 3))
 *         return (std::cos(this->get_time() * 4 * numbers::PI) * 4 * numbers::PI);
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
 * <a name="ImplementationofthecodeWaveEquationcodeclass"></a> 
 * <h3>Implementation of the <code>WaveEquation</code> class</h3>
 * 

 * 
 * The implementation of the actual logic is actually fairly short, since we
 * relegate things like assembling the matrices and right hand side vectors
 * to the library. The rest boils down to not much more than 130 lines of
 * actual code, a significant fraction of which is boilerplate code that can
 * be taken from previous example programs (e.g. the functions that solve
 * linear systems, or that generate output).
 *   

 * 
 * Let's start with the constructor (for an explanation of the choice of
 * time step, see the section on Courant, Friedrichs, and Lewy in the
 * introduction):
 * 
 * @code
 *   template <int dim>
 *   WaveEquation<dim>::WaveEquation()
 *     : fe(1)
 *     , dof_handler(triangulation)
 *     , time_step(1. / 64)
 *     , time(time_step)
 *     , timestep_number(1)
 *     , theta(0.5)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationsetup_system"></a> 
 * <h4>WaveEquation::setup_system</h4>
 * 

 * 
 * The next function is the one that sets up the mesh, DoFHandler, and
 * matrices and vectors at the beginning of the program, i.e. before the
 * first time step. The first few lines are pretty much standard if you've
 * read through the tutorial programs at least up to step-6:
 * 
 * @code
 *   template <int dim>
 *   void WaveEquation<dim>::setup_system()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -1, 1);
 *     triangulation.refine_global(7);
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl;
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl
 *               << std::endl;
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *     sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * Then comes a block where we have to initialize the 3 matrices we need
 * in the course of the program: the mass matrix, the Laplace matrix, and
 * the matrix $M+k^2\theta^2A$ used when solving for $U^n$ in each time
 * step.
 *     

 * 
 * When setting up these matrices, note that they all make use of the same
 * sparsity pattern object. Finally, the reason why matrices and sparsity
 * patterns are separate objects in deal.II (unlike in many other finite
 * element or linear algebra classes) becomes clear: in a significant
 * fraction of applications, one has to hold several matrices that happen
 * to have the same sparsity pattern, and there is no reason for them not
 * to share this information, rather than re-building and wasting memory
 * on it several times.
 *     

 * 
 * After initializing all of these matrices, we call library functions
 * that build the Laplace and mass matrices. All they need is a DoFHandler
 * object and a quadrature formula object that is to be used for numerical
 * integration. Note that in many respects these functions are better than
 * what we would usually do in application programs, for example because
 * they automatically parallelize building the matrices if multiple
 * processors are available in a machine: for more information see the
 * documentation of WorkStream or the
 * @ref threads "Parallel computing with multiple processors"
 * module. The matrices for solving linear systems will be filled in the
 * run() method because we need to re-apply boundary conditions every time
 * step.
 * 
 * @code
 *     mass_matrix.reinit(sparsity_pattern);
 *     laplace_matrix.reinit(sparsity_pattern);
 *     matrix_u.reinit(sparsity_pattern);
 *     matrix_v.reinit(sparsity_pattern);
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler,
 *                                       QGauss<dim>(fe.degree + 1),
 *                                       mass_matrix);
 *     MatrixCreator::create_laplace_matrix(dof_handler,
 *                                          QGauss<dim>(fe.degree + 1),
 *                                          laplace_matrix);
 * 
 * @endcode
 * 
 * The rest of the function is spent on setting vector sizes to the
 * correct value. The final line closes the hanging node constraints
 * object. Since we work on a uniformly refined mesh, no constraints exist
 * or have been computed (i.e. there was no need to call
 * DoFTools::make_hanging_node_constraints as in other programs), but we
 * need a constraints object in one place further down below anyway.
 * 
 * @code
 *     solution_u.reinit(dof_handler.n_dofs());
 *     solution_v.reinit(dof_handler.n_dofs());
 *     old_solution_u.reinit(dof_handler.n_dofs());
 *     old_solution_v.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.close();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationsolve_uandWaveEquationsolve_v"></a> 
 * <h4>WaveEquation::solve_u and WaveEquation::solve_v</h4>
 * 

 * 
 * The next two functions deal with solving the linear systems associated
 * with the equations for $U^n$ and $V^n$. Both are not particularly
 * interesting as they pretty much follow the scheme used in all the
 * previous tutorial programs.
 *   

 * 
 * One can make little experiments with preconditioners for the two matrices
 * we have to invert. As it turns out, however, for the matrices at hand
 * here, using Jacobi or SSOR preconditioners reduces the number of
 * iterations necessary to solve the linear system slightly, but due to the
 * cost of applying the preconditioner it is no win in terms of run-time. It
 * is not much of a loss either, but let's keep it simple and just do
 * without:
 * 
 * @code
 *   template <int dim>
 *   void WaveEquation<dim>::solve_u()
 *   {
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity());
 * 
 *     std::cout << "   u-equation: " << solver_control.last_step()
 *               << " CG iterations." << std::endl;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void WaveEquation<dim>::solve_v()
 *   {
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity());
 * 
 *     std::cout << "   v-equation: " << solver_control.last_step()
 *               << " CG iterations." << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationoutput_results"></a> 
 * <h4>WaveEquation::output_results</h4>
 * 

 * 
 * Likewise, the following function is pretty much what we've done
 * before. The only thing worth mentioning is how here we generate a string
 * representation of the time step number padded with leading zeros to 3
 * character length using the Utilities::int_to_string function's second
 * argument.
 * 
 * @code
 *   template <int dim>
 *   void WaveEquation<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution_u, "U");
 *     data_out.add_data_vector(solution_v, "V");
 * 
 *     data_out.build_patches();
 * 
 *     const std::string filename =
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
 * @endcode
 * 
 * Like step-15, since we write output at every time step (and the system
 * we have to solve is relatively easy), we instruct DataOut to use the
 * zlib compression algorithm that is optimized for speed instead of disk
 * usage since otherwise plotting the output becomes a bottleneck:
 * 
 * @code
 *     DataOutBase::VtkFlags vtk_flags;
 *     vtk_flags.compression_level =
 *       DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
 *     data_out.set_flags(vtk_flags);
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationrun"></a> 
 * <h4>WaveEquation::run</h4>
 * 

 * 
 * The following is really the only interesting function of the program. It
 * contains the loop over all time steps, but before we get to that we have
 * to set up the grid, DoFHandler, and matrices. In addition, we have to
 * somehow get started with initial values. To this end, we use the
 * VectorTools::project function that takes an object that describes a
 * continuous function and computes the $L^2$ projection of this function
 * onto the finite element space described by the DoFHandler object. Can't
 * be any simpler than that:
 * 
 * @code
 *   template <int dim>
 *   void WaveEquation<dim>::run()
 *   {
 *     setup_system();
 * 
 *     VectorTools::project(dof_handler,
 *                          constraints,
 *                          QGauss<dim>(fe.degree + 1),
 *                          InitialValuesU<dim>(),
 *                          old_solution_u);
 *     VectorTools::project(dof_handler,
 *                          constraints,
 *                          QGauss<dim>(fe.degree + 1),
 *                          InitialValuesV<dim>(),
 *                          old_solution_v);
 * 
 * @endcode
 * 
 * The next thing is to loop over all the time steps until we reach the
 * end time ($T=5$ in this case). In each time step, we first have to
 * solve for $U^n$, using the equation $(M^n + k^2\theta^2 A^n)U^n =$
 * $(M^{n,n-1} - k^2\theta(1-\theta) A^{n,n-1})U^{n-1} + kM^{n,n-1}V^{n-1}
 * +$ $k\theta \left[k \theta F^n + k(1-\theta) F^{n-1} \right]$. Note
 * that we use the same mesh for all time steps, so that $M^n=M^{n,n-1}=M$
 * and $A^n=A^{n,n-1}=A$. What we therefore have to do first is to add up
 * $MU^{n-1} - k^2\theta(1-\theta) AU^{n-1} + kMV^{n-1}$ and the forcing
 * terms, and put the result into the <code>system_rhs</code> vector. (For
 * these additions, we need a temporary vector that we declare before the
 * loop to avoid repeated memory allocations in each time step.)
 *     

 * 
 * The one thing to realize here is how we communicate the time variable
 * to the object describing the right hand side: each object derived from
 * the Function class has a time field that can be set using the
 * Function::set_time and read by Function::get_time. In essence, using
 * this mechanism, all functions of space and time are therefore
 * considered functions of space evaluated at a particular time. This
 * matches well what we typically need in finite element programs, where
 * we almost always work on a single time step at a time, and where it
 * never happens that, for example, one would like to evaluate a
 * space-time function for all times at any given spatial location.
 * 
 * @code
 *     Vector<double> tmp(solution_u.size());
 *     Vector<double> forcing_terms(solution_u.size());
 * 
 *     for (; time <= 5; time += time_step, ++timestep_number)
 *       {
 *         std::cout << "Time step " << timestep_number << " at t=" << time
 *                   << std::endl;
 * 
 *         mass_matrix.vmult(system_rhs, old_solution_u);
 * 
 *         mass_matrix.vmult(tmp, old_solution_v);
 *         system_rhs.add(time_step, tmp);
 * 
 *         laplace_matrix.vmult(tmp, old_solution_u);
 *         system_rhs.add(-theta * (1 - theta) * time_step * time_step, tmp);
 * 
 *         RightHandSide<dim> rhs_function;
 *         rhs_function.set_time(time);
 *         VectorTools::create_right_hand_side(dof_handler,
 *                                             QGauss<dim>(fe.degree + 1),
 *                                             rhs_function,
 *                                             tmp);
 *         forcing_terms = tmp;
 *         forcing_terms *= theta * time_step;
 * 
 *         rhs_function.set_time(time - time_step);
 *         VectorTools::create_right_hand_side(dof_handler,
 *                                             QGauss<dim>(fe.degree + 1),
 *                                             rhs_function,
 *                                             tmp);
 * 
 *         forcing_terms.add((1 - theta) * time_step, tmp);
 * 
 *         system_rhs.add(theta * time_step, forcing_terms);
 * 
 * @endcode
 * 
 * After so constructing the right hand side vector of the first
 * equation, all we have to do is apply the correct boundary
 * values. As for the right hand side, this is a space-time function
 * evaluated at a particular time, which we interpolate at boundary
 * nodes and then use the result to apply boundary values as we
 * usually do. The result is then handed off to the solve_u()
 * function:
 * 
 * @code
 *         {
 *           BoundaryValuesU<dim> boundary_values_u_function;
 *           boundary_values_u_function.set_time(time);
 * 
 *           std::map<types::global_dof_index, double> boundary_values;
 *           VectorTools::interpolate_boundary_values(dof_handler,
 *                                                    0,
 *                                                    boundary_values_u_function,
 *                                                    boundary_values);
 * 
 * @endcode
 * 
 * The matrix for solve_u() is the same in every time steps, so one
 * could think that it is enough to do this only once at the
 * beginning of the simulation. However, since we need to apply
 * boundary values to the linear system (which eliminate some matrix
 * rows and columns and give contributions to the right hand side),
 * we have to refill the matrix in every time steps before we
 * actually apply boundary data. The actual content is very simple:
 * it is the sum of the mass matrix and a weighted Laplace matrix:
 * 
 * @code
 *           matrix_u.copy_from(mass_matrix);
 *           matrix_u.add(theta * theta * time_step * time_step, laplace_matrix);
 *           MatrixTools::apply_boundary_values(boundary_values,
 *                                              matrix_u,
 *                                              solution_u,
 *                                              system_rhs);
 *         }
 *         solve_u();
 * 
 * 
 * @endcode
 * 
 * The second step, i.e. solving for $V^n$, works similarly, except
 * that this time the matrix on the left is the mass matrix (which we
 * copy again in order to be able to apply boundary conditions, and
 * the right hand side is $MV^{n-1} - k\left[ \theta A U^n +
 * (1-\theta) AU^{n-1}\right]$ plus forcing terms. Boundary values
 * are applied in the same way as before, except that now we have to
 * use the BoundaryValuesV class:
 * 
 * @code
 *         laplace_matrix.vmult(system_rhs, solution_u);
 *         system_rhs *= -theta * time_step;
 * 
 *         mass_matrix.vmult(tmp, old_solution_v);
 *         system_rhs += tmp;
 * 
 *         laplace_matrix.vmult(tmp, old_solution_u);
 *         system_rhs.add(-time_step * (1 - theta), tmp);
 * 
 *         system_rhs += forcing_terms;
 * 
 *         {
 *           BoundaryValuesV<dim> boundary_values_v_function;
 *           boundary_values_v_function.set_time(time);
 * 
 *           std::map<types::global_dof_index, double> boundary_values;
 *           VectorTools::interpolate_boundary_values(dof_handler,
 *                                                    0,
 *                                                    boundary_values_v_function,
 *                                                    boundary_values);
 *           matrix_v.copy_from(mass_matrix);
 *           MatrixTools::apply_boundary_values(boundary_values,
 *                                              matrix_v,
 *                                              solution_v,
 *                                              system_rhs);
 *         }
 *         solve_v();
 * 
 * @endcode
 * 
 * Finally, after both solution components have been computed, we
 * output the result, compute the energy in the solution, and go on to
 * the next time step after shifting the present solution into the
 * vectors that hold the solution at the previous time step. Note the
 * function SparseMatrix::matrix_norm_square that can compute
 * $\left<V^n,MV^n\right>$ and $\left<U^n,AU^n\right>$ in one step,
 * saving us the expense of a temporary vector and several lines of
 * code:
 * 
 * @code
 *         output_results();
 * 
 *         std::cout << "   Total energy: "
 *                   << (mass_matrix.matrix_norm_square(solution_v) +
 *                       laplace_matrix.matrix_norm_square(solution_u)) /
 *                        2
 *                   << std::endl;
 * 
 *         old_solution_u = solution_u;
 *         old_solution_v = solution_v;
 *       }
 *   }
 * } // namespace Step23
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * What remains is the main function of the program. There is nothing here
 * that hasn't been shown in several of the previous programs:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step23;
 * 
 *       WaveEquation<2> wave_equation_solver;
 *       wave_equation_solver.run();
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


When the program is run, it produces the following output:
@code
Number of active cells: 16384
Number of degrees of freedom: 16641

Time step 1 at t=0.015625
   u-equation: 8 CG iterations.
   v-equation: 22 CG iterations.
   Total energy: 1.17887
Time step 2 at t=0.03125
   u-equation: 8 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 2.9655
Time step 3 at t=0.046875
   u-equation: 8 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 4.33761
Time step 4 at t=0.0625
   u-equation: 7 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 5.35499
Time step 5 at t=0.078125
   u-equation: 7 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 6.18652
Time step 6 at t=0.09375
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 6.6799

...

Time step 31 at t=0.484375
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 21.9068
Time step 32 at t=0.5
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.3394
Time step 33 at t=0.515625
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019

...

Time step 319 at t=4.98438
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019
Time step 320 at t=5
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019
@endcode

What we see immediately is that the energy is a constant at least after
$t=\frac 12$ (until which the boundary source term $g$ is nonzero, injecting
energy into the system).

In addition to the screen output, the program writes the solution of each time
step to an output file. If we process them adequately and paste them into a
movie, we get the following:

<img src="https://www.dealii.org/images/steps/developer/step-23.movie.gif" alt="Animation of the solution of step 23.">

The movie shows the generated wave nice traveling through the domain and back,
being reflected at the clamped boundary. Some numerical noise is trailing the
wave, an artifact of a too-large mesh size that can be reduced by reducing the
mesh width and the time step.


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


If you want to explore a bit, try out some of the following things:
<ul>
  <li>Varying $\theta$. This gives different time stepping schemes, some of
  which are stable while others are not. Take a look at how the energy
  evolves.

  <li>Different initial and boundary conditions, right hand sides.

  <li>More complicated domains or more refined meshes. Remember that the time
  step needs to be bounded by the mesh width, so changing the mesh should
  always involve also changing the time step. We will come back to this issue
  in step-24.

  <li>Variable coefficients: In real media, the wave speed is often
  variable. In particular, the "real" wave equation in realistic media would
  read
  @f[
     \rho(x) \frac{\partial^2 u}{\partial t^2}
     -
     \nabla \cdot
     a(x) \nabla u = f,
  @f]
  where $\rho(x)$ is the density of the material, and $a(x)$ is related to the
  stiffness coefficient. The wave speed is then $c=\sqrt{a/\rho}$.

  To make such a change, we would have to compute the mass and Laplace
  matrices with a variable coefficient. Fortunately, this isn't too hard: the
  functions MatrixCreator::create_laplace_matrix and
  MatrixCreator::create_mass_matrix have additional default parameters that can
  be used to pass non-constant coefficient functions to them. The required
  changes are therefore relatively small. On the other hand, care must be
  taken again to make sure the time step is within the allowed range.

  <li>In the in-code comments, we discussed the fact that the matrices for
  solving for $U^n$ and $V^n$ need to be reset in every time because of
  boundary conditions, even though the actual content does not change. It is
  possible to avoid copying by not eliminating columns in the linear systems,
  which is implemented by appending a @p false argument to the call:
  @code
    MatrixTools::apply_boundary_values(boundary_values,
                                       matrix_u,
                                       solution_u,
                                       system_rhs,
                                       false);
  @endcode

  <li>deal.II being a library that supports adaptive meshes it would of course be
  nice if this program supported change the mesh every few time steps. Given the
  structure of the solution &mdash; a wave that travels through the domain &mdash;
  it would seem appropriate if we only refined the mesh where the wave currently is,
  and not simply everywhere. It is intuitively clear that we should be able to
  save a significant amount of cells this way. (Though upon further thought one
  realizes that this is really only the case in the initial stages of the simulation.
  After some time, for wave phenomena, the domain is filled with reflections of
  the initial wave going in every direction and filling every corner of the domain.
  At this point, there is in general little one can gain using local mesh
  refinement.)

  To make adaptively changing meshes possible, there are basically two routes.
  The "correct" way would be to go back to the weak form we get using Rothe's
  method. For example, the first of the two equations to be solved in each time
  step looked like this:
  \f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla u^{n-1},\nabla \varphi)
  +
  k(v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right].
  \f}
  Now, note that we solve for $u^n$ on mesh ${\mathbb T}^n$, and
  consequently the test functions $\varphi$ have to be from the space
  $V_h^n$ as well. As discussed in the introduction, terms like
  $(u^{n-1},\varphi)$ then require us to integrate the solution of the
  previous step (which may have been computed on a different mesh
  ${\mathbb T}^{n-1}$) against the test functions of the current mesh,
  leading to a matrix $M^{n,n-1}$. This process of integrating shape
  functions from different meshes is, at best, awkward. It can be done
  but because it is difficult to ensure that ${\mathbb T}^{n-1}$ and
  ${\mathbb T}^{n}$ differ by at most one level of refinement, one
  has to recursively match cells from both meshes. It is feasible to
  do this, but it leads to lengthy and not entirely obvious code.

  The second approach is the following: whenever we change the mesh,
  we simply interpolate the solution from the last time step on the old
  mesh to the new mesh, using the SolutionTransfer class. In other words,
  instead of the equation above, we would solve
  \f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (I^n u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla I^n u^{n-1},\nabla \varphi)
  +
  k(I^n v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right],
  \f}
  where $I^n$ interpolates a given function onto mesh ${\mathbb T}^n$.
  This is a much simpler approach because, in each time step, we no
  longer have to worry whether $u^{n-1},v^{n-1}$ were computed on the
  same mesh as we are using now or on a different mesh. Consequently,
  the only changes to the code necessary are the addition of a function
  that computes the error, marks cells for refinement, sets up a
  SolutionTransfer object, transfers the solution to the new mesh, and
  rebuilds matrices and right hand side vectors on the new mesh. Neither
  the functions building the matrices and right hand sides, nor the
  solvers need to be changed.

  While this second approach is, strictly speaking,
  not quite correct in the Rothe framework (it introduces an addition source
  of error, namely the interpolation), it is nevertheless what
  almost everyone solving time dependent equations does. We will use this
  method in step-31, for example.
</ul>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-23.cc"
*/
