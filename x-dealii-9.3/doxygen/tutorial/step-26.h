/**
@page step_26 The step-26 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Adaptingmeshesfortimedependentproblems"> Adapting meshes for time dependent problems </a>
        <li><a href="#WhatcouldpossiblygowrongVerifyingwhetherthecodeiscorrect"> What could possibly go wrong? Verifying whether the code is correct </a>
        <li><a href="#Thetestcase"> The testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ThecodeHeatEquationcodeclass">The <code>HeatEquation</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ThecodeHeatEquationcodeimplementation">The <code>HeatEquation</code> implementation</a>
      <ul>
        <li><a href="#codeHeatEquationsetup_systemcode"><code>HeatEquation::setup_system</code></a>
        <li><a href="#codeHeatEquationsolve_time_stepcode"><code>HeatEquation::solve_time_step</code></a>
        <li><a href="#codeHeatEquationoutput_resultscode"><code>HeatEquation::output_results</code></a>
        <li><a href="#codeHeatEquationrefine_meshcode"><code>HeatEquation::refine_mesh</code></a>
        <li><a href="#codeHeatEquationruncode"><code>HeatEquation::run</code></a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Adaptivetimestepping">Adaptive time stepping</a>
        <li><a href="#Bettertimesteppingmethods">Better time stepping methods</a>
        <li><a href="#Betterrefinementcriteria">Better refinement criteria</a>
        <li><a href="#Positivitypreservation">Positivity preservation</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{29,30}
(@dealiiVideoLectureSeeAlso{31.7})


This program implements the heat equation
@f{align*}
  \frac{\partial u(\mathbf x, t)}{\partial t}
  -
  \Delta u(\mathbf x, t)
  &=
  f(\mathbf x, t),
  \qquad\qquad &&
  \forall \mathbf x \in \Omega, t\in (0,T),
  \\
  u(\mathbf x, 0) &= u_0(\mathbf x) &&
  \forall \mathbf x \in \Omega, \\
  \\
  u(\mathbf x, t) &= g(\mathbf x,t) &&
  \forall \mathbf x \in \partial\Omega, t \in (0,T).
@f}
In some sense, this equation is simpler than the ones we have discussed in the
preceding programs step-23, step-24, step-25, namely the wave equation. This
is due to the fact that the heat equation smoothes out the solution over time,
and is consequently more forgiving in many regards. For example, when using
implicit time stepping methods, we can actually take large time steps, we have
less trouble with the small disturbances we introduce through adapting the
mesh every few time steps, etc.

Our goal here will be to solve the equations above using the theta-scheme that
discretizes the equation in time using the following approach, where we would
like $u^n(\mathbf x)$ to approximate $u(\mathbf x, t_n)$ at some time $t_n$:
@f{align*}
  \frac{u^n(\mathbf x)-u^{n-1}(\mathbf x)}{k_n}
  -
  \left[
  (1-\theta)\Delta u^{n-1}(\mathbf x)
  +
  \theta\Delta u^n(\mathbf x)
  \right]
  &=
  \left[
  (1-\theta)f(\mathbf x, t_{n-1})
  +
  \theta f(\mathbf x, t_n)
  \right].
@f}
Here, $k_n=t_n-t_{n-1}$ is the time step size. The theta-scheme generalizes
the explicit Euler ($\theta=0$), implicit Euler ($\theta=1$) and
Crank-Nicolson ($\theta=\frac 12$) time discretizations. Since the latter has
the highest convergence order, we will choose $\theta=\frac 12$ in the program
below, but make it so that playing with this parameter remains simple. (If you
are interested in playing with higher order methods, take a look at step-52.)

Given this time discretization, space discretization happens as it always
does, by multiplying with test functions, integrating by parts, and then
restricting everything to a finite dimensional subspace. This yields the
following set of fully discrete equations after multiplying through with
$k_n$:
@f{align*}
  M U^n-MU^{n-1}
  +
  k_n \left[
  (1-\theta)A U^{n-1}
  +
  \theta A U^n
  \right]
  &=
  k_n
  \left[
  (1-\theta)F^{n-1}
  +
  \theta F^n
  \right],
@f}
where $M$ is the mass matrix and $A$ is the stiffness matrix that results from
discretizing the Laplacian. Bringing all known quantities to the right hand
side yields the linear system we have to solve in every step:
@f{align*}
  (M
  +
  k_n \theta A) U^n
  &=
  MU^{n-1}
  -
  k_n
  (1-\theta)A U^{n-1}
  +
  k_n
  \left[
  (1-\theta)F^{n-1}
  +
  \theta F^n
  \right].
@f}
The linear system on the left hand side is symmetric and positive definite, so
we should have no trouble solving it with the Conjugate Gradient method.

We can start the iteration above if we have the set of nodal coefficients
$U^0$ at the initial time. Here, we take the ones we get by interpolating the
initial values $u_0(\mathbf x)$ onto the mesh used for the first time step. We
will also need to choose a time step; we will here just choose it as fixed,
but clearly advanced simulators will want to choose it adaptively. We will
briefly come back to this in the <a href="#Results">results section
below</a>.


<a name="Adaptingmeshesfortimedependentproblems"></a><h3> Adapting meshes for time dependent problems </h3>


When solving the wave equation and its variants in the previous few programs,
we kept the mesh fixed. Just as for stationary equations, one can make a good
case that this is not the smartest approach and that significant savings can
be had by adapting the mesh. There are, however, significant difficulties
compared to the stationary case. Let us go through them in turn:

<ul>
  <li><i>Time step size and minimal mesh size</i>: For stationary problems, the
  general approach is "make the mesh as fine as it is necessary". For problems
  with singularities, this often leads to situations where we get many levels
  of refinement into corners or along interfaces. The very first tutorial to
  use adaptive meshes, step-6, is a point in case already.

  However, for time dependent problems, we typically need to choose the time
  step related to the mesh size. For explicit time discretizations, this is
  obvious, since we need to respect a CFL condition that ties the time step
  size to the smallest mesh size. For implicit time discretizations, no such
  hard restriction exists, but in practice we still want to make the time step
  smaller if we make the mesh size smaller since we typically have error
  estimates of the form $\|e\| \le {\cal O}(k^p + h^q)$ where $p,q$ are the
  convergence orders of the time and space discretization, respectively. We
  can only make the error small if we decrease both terms. Ideally, an
  estimate like this would suggest to choose $k \propto h^{q/p}$. Because, at
  least for problems with non-smooth solutions, the error is typically
  localized in the cells with the smallest mesh size, we have to indeed choose
  $k \propto h_{\text{min}}^{q/p}$, using the <i>smallest</i> mesh size.

  The consequence is that refining the mesh further in one place implies not
  only the moderate additional effort of increasing the number of degrees of
  freedom slightly, but also the much larger effort of having the solve the
  <i>global</i> linear system more often because of the smaller time step.

  In practice, one typically deals with this by acknowledging that we can not
  make the time step arbitrarily small, and consequently can not make the
  local mesh size arbitrarily small. Rather, we set a maximal level of
  refinement and when we flag cells for refinement, we simply do not refine
  those cells whose children would exceed this maximal level of refinement.

  There is a similar problem in that we will choose a right hand side that
  will switch on in different parts of the domain at different times. To avoid
  being caught flat footed with too coarse a mesh in areas where we suddenly
  need a finer mesh, we will also enforce in our program a <i>minimal</i> mesh
  refinement level.

  <li><i>Test functions from different meshes</i>: Let us consider again the
  semi-discrete equations we have written down above:
  @f{align*}
    \frac{u^n(\mathbf x)-u^{n-1}(\mathbf x)}{k_n}
    -
    \left[
    (1-\theta)\Delta u^{n-1}(\mathbf x)
    +
    \theta\Delta u^n(\mathbf x)
    \right]
    &=
    \left[
    (1-\theta)f(\mathbf x, t_{n-1})
    +
    \theta f(\mathbf x, t_n)
    \right].
  @f}
  We can here consider $u^{n-1}$ as data since it has presumably been computed
  before. Now, let us replace
  @f{align*}
    u^n(\mathbf x)\approx u_h^n(\mathbf x)
    =
    \sum_j U^n \varphi_j(\mathbf x),
  @f}
  multiply with test functions $\varphi_i(\mathbf x)$ and integrate by parts
  where necessary. In a process as outlined above, this would yield
  @f{align*}
    \sum_j
    (M
    +
    k_n \theta A)_{ij} U^n_j
    &=
    (\varphi_i, u_h^{n-1})
    -
    k_n
    (1-\theta)(\nabla \varphi_i, \nabla u_h^{n-1})
    +
    k_n
    \left[
    (1-\theta)F^{n-1}
    +
    \theta F^n
    \right].
  @f}
  Now imagine that we have changed the mesh between time steps $n-1$ and
  $n$. Then the problem is that the basis functions we use for $u_h^n$ and
  $u^{n-1}$ are different! This pertains to the terms on the right hand side,
  the first of which we could more clearly write as (the second follows the
  same pattern)
  @f{align*}
    (\varphi_i, u_h^{n-1})
    =
    (\varphi_i^n, u_h^{n-1})
    =
    \sum_{j=1}^{N_{n-1}}
    (\varphi_i^n, \varphi_j^{n-1}) U^{n-1}_j,
    \qquad\qquad
    i=1\ldots N_n.
  @f}
  If the meshes used in these two time steps are the same, then
  $(\varphi_i^n, \varphi_j^{n-1})$ forms a square mass matrix
  $M_{ij}$. However, if the meshes are not the same, then in general the matrix
  is rectangular. Worse, it is difficult to even compute these integrals
  because if we loop over the cells of the mesh at time step $n$, then we need
  to evaluate $\varphi_j^{n-1}$ at the quadrature points of these cells, but
  they do not necessarily correspond to the cells of the mesh at time step
  $n-1$ and $\varphi_j^{n-1}$ is not defined via these cells; the same of
  course applies if we wanted to compute the integrals via integration on the
  cells of mesh $n-1$.

  In any case, what we have to face is a situation where we need to integrate
  shape functions defined on two different meshes. This can be done, and is in
  fact demonstrated in step-28, but the process is at best described by the
  word "awkward".

  In practice, one does not typically want to do this. Rather, we avoid the
  whole situation by interpolating the solution from the old to the new mesh
  every time we adapt the mesh. In other words, rather than solving the
  equations above, we instead solve the problem
  @f{align*}
    \sum_j
    (M
    +
    k_n \theta A)_{ij} U^n_j
    &=
    (\varphi_i, I_h^n u_h^{n-1})
    -
    k_n
    (1-\theta)(\nabla \varphi_i, \nabla I_h^n u_h^{n-1})
    +
    k_n
    \left[
    (1-\theta)F^{n-1}
    +
    \theta F^n
    \right],
  @f}
  where $I_h^n$ is the interpolation operator onto the finite element space
  used in time step $n$. This is not the optimal approach since it introduces
  an additional error besides time and space discretization, but it is a
  pragmatic one that makes it feasible to do time adapting meshes.
</ul>



<a name="WhatcouldpossiblygowrongVerifyingwhetherthecodeiscorrect"></a><h3> What could possibly go wrong? Verifying whether the code is correct </h3>


There are a number of things one can typically get wrong when implementing a
finite element code. In particular, for time dependent problems, the following
are common sources of bugs:
- The time integration, for example by getting the coefficients in front of
  the terms involving the current and previous time steps wrong (e.g., mixing
  up a factor $\theta$ for $1-\theta$).
- Handling the right hand side, for example forgetting a factor of $k_n$ or
  $\theta$.
- Mishandling the boundary values, again for example forgetting a factor of
  $k_n$ or $\theta$, or forgetting to apply nonzero boundary values not only
  to the right hand side but also to the system matrix.

A less common problem is getting the initial conditions wrong because one can
typically see that it is wrong by just outputting the first time step. In any
case, in order to verify the correctness of the code, it is helpful to have a
testing protocol that allows us to verify each of these components
separately. This means:
- Testing the code with nonzero initial conditions but zero right hand side
  and boundary values and verifying that the time evolution is correct.
- Then testing with zero initial conditions and boundary values but nonzero
  right hand side and again ensuring correctness.
- Finally, testing with zero initial conditions and right hand side but
  nonzero boundary values.

This sounds complicated, but fortunately, for linear partial differential
equations without coefficients (or constant coefficients) like the one here,
there is a fairly standard protocol that rests on the following observation:
if you choose as your domain a square $[0,1]^2$ (or, with slight
modifications, a rectangle), then the exact solution can be written as
@f{align*}
  u(x,y,t) = a(t) \sin(n_x \pi x) \sin(n_y \pi y)
@f}
(with integer constants $n_x,n_y$)
if only the initial condition, right hand side and boundary values are all
of the form $\sin(n_x \pi x) \sin(n_y \pi y)$ as well. This is due to the fact
that the function $\sin(n_x \pi x) \sin(n_y \pi y)$ is an eigenfunction of the
Laplace operator and allows us to compute things like the time factor $a(t)$
analytically and, consequently, compare with what we get numerically.

As an example, let us consider the situation where we have
$u_0(x,y)=\sin(n_x \pi x) \sin(n_x \pi y)$ and
$f(x,y,t)=0$. With the claim (ansatz) of the form for
$u(x,y,t)$ above, we get that
@f{align*}
  \left(\frac{\partial}{\partial t} -\Delta\right)
  u(x,y,t)
  &=
  \left(\frac{\partial}{\partial t} -\Delta\right)
  a(t) \sin(n_x \pi x) \sin(n_y \pi y)
  \\
  &=
  \left(a'(t) + (n_x^2+n_y^2)\pi^2 a(t) \right) \sin(n_x \pi x) \sin(n_y \pi y).
@f}
For this to be equal to $f(x,y,t)=0$, we need that
@f{align*}
  a'(t) + (n_x^2+n_y^2)\pi^2 a(t) = 0
@f}
and due to the initial conditions, $a(0)=1$. This differential equation can be
integrated to yield
@f{align*}
  a(t) = - e^{-(n_x^2+n_y^2)\pi^2 t}.
@f}
In other words, if the initial condition is a product of sines, then the
solution has exactly the same shape of a product of sines that decays to zero
with a known time dependence. This is something that is easy to test if you
have a sufficiently fine mesh and sufficiently small time step.

What is typically going to happen if you get the time integration scheme wrong
(e.g., by having the wrong factors of $\theta$ or $k$ in front of the various
terms) is that you don't get the right temporal behavior of the
solution. Double check the various factors until you get the right
behavior. You may also want to verify that the temporal decay rate (as
determined, for example, by plotting the value of the solution at a fixed
point) does not double or halve each time you double or halve the time step or
mesh size. You know that it's not the handling of the
boundary conditions or right hand side because these were both zero.

If you have so verified that the time integrator is correct, take the
situation where the right hand side is nonzero but the initial conditions are
zero: $u_0(x,y)=0$ and
$f(x,y,t)=\sin(n_x \pi x) \sin(n_x \pi y)$. Again,
@f{align*}
  \left(\frac{\partial}{\partial t} -\Delta\right)
  u(x,y,t)
  &=
  \left(\frac{\partial}{\partial t} -\Delta\right)
  a(t) \sin(n_x \pi x) \sin(n_y \pi y)
  \\
  &=
  \left(a'(t) + (n_x^2+n_y^2)\pi^2 a(t) \right) \sin(n_x \pi x) \sin(n_y \pi y),
@f}
and for this to be equal to $f(x,y,t)$, we need that
@f{align*}
  a'(t) + (n_x^2+n_y^2)\pi^2 a(t) = 1
@f}
and due to the initial conditions, $a(0)=0$. Integrating this equation in time
yields
@f{align*}
  a(t) = \frac{1}{(n_x^2+n_y^2)\pi^2} \left[ 1 - e^{-(n_x^2+n_y^2)\pi^2 t} \right].
@f}

Again, if you have the wrong factors of $\theta$ or $k$ in front of the right
hand side terms you will either not get the right temporal behavior of the
solution, or it will converge to a maximum value other than
$\frac{1}{(n_x^2+n_y^2)\pi^2}$.

Once we have verified that the time integration and right hand side handling
are correct using this scheme, we can go on to verifying that we have the
boundary values correct, using a very similar approach.



<a name="Thetestcase"></a><h3> The testcase </h3>


Solving the heat equation on a simple domain with a simple right hand side
almost always leads to solutions that are exceedingly boring, since they
become very smooth very quickly and then do not move very much any
more. Rather, we here solve the equation on the L-shaped domain with zero
Dirichlet boundary values and zero initial conditions, but as right hand side
we choose
@f{align*}
  f(\mathbf x, t)
  =
  \left\{
  \begin{array}{ll}
    \chi_1(\mathbf x)
    & \text{if \(0\le t \le 0.2\tau\) or \(\tau\le t \le 1.2\tau\) or \(2\tau\le t
    \le 2.2\tau\), etc}
    \\
    \chi_2(\mathbf x)
    & \text{if \(0.5\le t \le 0.7\tau\) or \(1.5\tau\le t \le 1.7\tau\) or \(2.5\tau\le t
    \le 2.7\tau\), etc}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.
@f}
Here,
@f{align*}
  \chi_1(\mathbf x) &=
  \left\{
  \begin{array}{ll}
    1
    & \text{if \(x>0.5\) and \(y>-0.5\)}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.
  \\
  \chi_2(\mathbf x) &=
  \left\{
  \begin{array}{ll}
    1
    & \text{if \(x>-0.5\) and \(y>0.5\)}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.
@f}
In other words, in every period of length $\tau$, the right hand side first
flashes on in domain 1, then off completely, then on in domain 2, then off
completely again. This pattern is probably best observed via the little
animation of the solution shown in the <a href="#Results">results
section</a>.

If you interpret the heat equation as finding the spatially and temporally
variable temperature distribution of a conducting solid, then the test case
above corresponds to an L-shaped body where we keep the boundary at zero
temperature, and heat alternatingly in two parts of the domain. While heating
is in effect, the temperature rises in these places, after which it diffuses
and diminishes again. The point of these initial conditions is that they
provide us with a solution that has singularities both in time (when sources
switch on and off) as well as time (at the reentrant corner as well as at the
edges and corners of the regions where the source acts).
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The program starts with the usual include files, all of which you should
 * have seen before by now:
 * 
 * @code
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/solution_transfer.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * Then the usual placing of all content of this program into a namespace and
 * the importation of the deal.II namespace into the one we will work in:
 * 
 * @code
 * namespace Step26
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeHeatEquationcodeclass"></a> 
 * <h3>The <code>HeatEquation</code> class</h3>
 *   

 * 
 * The next piece is the declaration of the main class of this program. It
 * follows the well trodden path of previous examples. If you have looked at
 * step-6, for example, the only thing worth noting here is that we need to
 * build two matrices (the mass and Laplace matrix) and keep the current and
 * previous time step's solution. We then also need to store the current
 * time, the size of the time step, and the number of the current time
 * step. The last of the member variables denotes the theta parameter
 * discussed in the introduction that allows us to treat the explicit and
 * implicit Euler methods as well as the Crank-Nicolson method and other
 * generalizations all in one program.
 *   

 * 
 * As far as member functions are concerned, the only possible surprise is
 * that the <code>refine_mesh</code> function takes arguments for the
 * minimal and maximal mesh refinement level. The purpose of this is
 * discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   class HeatEquation
 *   {
 *   public:
 *     HeatEquation();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void solve_time_step();
 *     void output_results() const;
 *     void refine_mesh(const unsigned int min_grid_level,
 *                      const unsigned int max_grid_level);
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
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> old_solution;
 *     Vector<double> system_rhs;
 * 
 *     double       time;
 *     double       time_step;
 *     unsigned int timestep_number;
 * 
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
 * In the following classes and functions, we implement the various pieces
 * of data that define this problem (right hand side and boundary values)
 * that are used in this program and for which we need function objects. The
 * right hand side is chosen as discussed at the end of the
 * introduction. For boundary values, we choose zero values, but this is
 * easily changed below.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     RightHandSide()
 *       : Function<dim>()
 *       , period(0.2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 *     const double period;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 *     (void)component;
 *     AssertIndexRange(component, 1);
 *     Assert(dim == 2, ExcNotImplemented());
 * 
 *     const double time = this->get_time();
 *     const double point_within_period =
 *       (time / period - std::floor(time / period));
 * 
 *     if ((point_within_period >= 0.0) && (point_within_period <= 0.2))
 *       {
 *         if ((p[0] > 0.5) && (p[1] > -0.5))
 *           return 1;
 *         else
 *           return 0;
 *       }
 *     else if ((point_within_period >= 0.5) && (point_within_period <= 0.7))
 *       {
 *         if ((p[0] > -0.5) && (p[1] > 0.5))
 *           return 1;
 *         else
 *           return 0;
 *       }
 *     else
 *       return 0;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/,
 *                                     const unsigned int component) const
 *   {
 *     (void)component;
 *     Assert(component == 0, ExcIndexRange(component, 0, 1));
 *     return 0;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeHeatEquationcodeimplementation"></a> 
 * <h3>The <code>HeatEquation</code> implementation</h3>
 *   

 * 
 * It is time now for the implementation of the main class. Let's
 * start with the constructor which selects a linear element, a time
 * step constant at 1/500 (remember that one period of the source
 * on the right hand side was set to 0.2 above, so we resolve each
 * period with 100 time steps) and chooses the Crank Nicolson method
 * by setting $\theta=1/2$.
 * 
 * @code
 *   template <int dim>
 *   HeatEquation<dim>::HeatEquation()
 *     : fe(1)
 *     , dof_handler(triangulation)
 *     , time_step(1. / 500)
 *     , theta(0.5)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationsetup_systemcode"></a> 
 * <h4><code>HeatEquation::setup_system</code></h4>
 *   

 * 
 * The next function is the one that sets up the DoFHandler object,
 * computes the constraints, and sets the linear algebra objects
 * to their correct sizes. We also compute the mass and Laplace
 * matrix here by simply calling two functions in the library.
 *   

 * 
 * Note that we do not take the hanging node constraints into account when
 * assembling the matrices (both functions have an AffineConstraints argument
 * that defaults to an empty object). This is because we are going to
 * condense the constraints in run() after combining the matrices for the
 * current time-step.
 * 
 * @code
 *   template <int dim>
 *   void HeatEquation<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::cout << std::endl
 *               << "===========================================" << std::endl
 *               << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl
 *               << std::endl;
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ true);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     mass_matrix.reinit(sparsity_pattern);
 *     laplace_matrix.reinit(sparsity_pattern);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler,
 *                                       QGauss<dim>(fe.degree + 1),
 *                                       mass_matrix);
 *     MatrixCreator::create_laplace_matrix(dof_handler,
 *                                          QGauss<dim>(fe.degree + 1),
 *                                          laplace_matrix);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     old_solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationsolve_time_stepcode"></a> 
 * <h4><code>HeatEquation::solve_time_step</code></h4>
 *   

 * 
 * The next function is the one that solves the actual linear system
 * for a single time step. There is nothing surprising here:
 * 
 * @code
 *   template <int dim>
 *   void HeatEquation<dim>::solve_time_step()
 *   {
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.0);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 * 
 *     std::cout << "     " << solver_control.last_step() << " CG iterations."
 *               << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationoutput_resultscode"></a> 
 * <h4><code>HeatEquation::output_results</code></h4>
 *   

 * 
 * Neither is there anything new in generating graphical output other than the
 * fact that we tell the DataOut object what the current time and time step
 * number is, so that this can be written into the output file:
 * 
 * @code
 *   template <int dim>
 *   void HeatEquation<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "U");
 * 
 *     data_out.build_patches();
 * 
 *     data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));
 * 
 *     const std::string filename =
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk";
 *     std::ofstream output(filename);
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationrefine_meshcode"></a> 
 * <h4><code>HeatEquation::refine_mesh</code></h4>
 *   

 * 
 * This function is the interesting part of the program. It takes care of
 * the adaptive mesh refinement. The three tasks
 * this function performs is to first find out which cells to
 * refine/coarsen, then to actually do the refinement and eventually
 * transfer the solution vectors between the two different grids. The first
 * task is simply achieved by using the well-established Kelly error
 * estimator on the solution. The second task is to actually do the
 * remeshing. That involves only basic functions as well, such as the
 * <code>refine_and_coarsen_fixed_fraction</code> that refines those cells
 * with the largest estimated error that together make up 60 per cent of the
 * error, and coarsens those cells with the smallest error that make up for
 * a combined 40 per cent of the error. Note that for problems such as the
 * current one where the areas where something is going on are shifting
 * around, we want to aggressively coarsen so that we can move cells
 * around to where it is necessary.
 *   

 * 
 * As already discussed in the introduction, too small a mesh leads to
 * too small a time step, whereas too large a mesh leads to too little
 * resolution. Consequently, after the first two steps, we have two
 * loops that limit refinement and coarsening to an allowable range of
 * cells:
 * 
 * @code
 *   template <int dim>
 *   void HeatEquation<dim>::refine_mesh(const unsigned int min_grid_level,
 *                                       const unsigned int max_grid_level)
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       solution,
 *       estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.6,
 *                                                       0.4);
 * 
 *     if (triangulation.n_levels() > max_grid_level)
 *       for (const auto &cell :
 *            triangulation.active_cell_iterators_on_level(max_grid_level))
 *         cell->clear_refine_flag();
 *     for (const auto &cell :
 *          triangulation.active_cell_iterators_on_level(min_grid_level))
 *       cell->clear_coarsen_flag();
 * @endcode
 * 
 * These two loops above are slightly different but this is easily
 * explained. In the first loop, instead of calling
 * <code>triangulation.end()</code> we may as well have called
 * <code>triangulation.end_active(max_grid_level)</code>. The two
 * calls should yield the same iterator since iterators are sorted
 * by level and there should not be any cells on levels higher than
 * on level <code>max_grid_level</code>. In fact, this very piece
 * of code makes sure that this is the case.
 * 

 * 
 * As part of mesh refinement we need to transfer the solution vectors
 * from the old mesh to the new one. To this end we use the
 * SolutionTransfer class and we have to prepare the solution vectors that
 * should be transferred to the new grid (we will lose the old grid once
 * we have done the refinement so the transfer has to happen concurrently
 * with refinement). At the point where we call this function, we will
 * have just computed the solution, so we no longer need the old_solution
 * variable (it will be overwritten by the solution just after the mesh
 * may have been refined, i.e., at the end of the time step; see below).
 * In other words, we only need the one solution vector, and we copy it
 * to a temporary object where it is safe from being reset when we further
 * down below call <code>setup_system()</code>.
 *     

 * 
 * Consequently, we initialize a SolutionTransfer object by attaching
 * it to the old DoF handler. We then prepare the triangulation and the
 * data vector for refinement (in this order).
 * 
 * @code
 *     SolutionTransfer<dim> solution_trans(dof_handler);
 * 
 *     Vector<double> previous_solution;
 *     previous_solution = solution;
 *     triangulation.prepare_coarsening_and_refinement();
 *     solution_trans.prepare_for_coarsening_and_refinement(previous_solution);
 * 
 * @endcode
 * 
 * Now everything is ready, so do the refinement and recreate the DoF
 * structure on the new grid, and finally initialize the matrix structures
 * and the new vectors in the <code>setup_system</code> function. Next, we
 * actually perform the interpolation of the solution from old to new
 * grid. The final step is to apply the hanging node constraints to the
 * solution vector, i.e., to make sure that the values of degrees of
 * freedom located on hanging nodes are so that the solution is
 * continuous. This is necessary since SolutionTransfer only operates on
 * cells locally, without regard to the neighborhoof.
 * 
 * @code
 *     triangulation.execute_coarsening_and_refinement();
 *     setup_system();
 * 
 *     solution_trans.interpolate(previous_solution, solution);
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationruncode"></a> 
 * <h4><code>HeatEquation::run</code></h4>
 *   

 * 
 * This is the main driver of the program, where we loop over all
 * time steps. At the top of the function, we set the number of
 * initial global mesh refinements and the number of initial cycles of
 * adaptive mesh refinement by repeating the first time step a few
 * times. Then we create a mesh, initialize the various objects we will
 * work with, set a label for where we should start when re-running
 * the first time step, and interpolate the initial solution onto
 * out mesh (we choose the zero function here, which of course we could
 * do in a simpler way by just setting the solution vector to zero). We
 * also output the initial time step once.
 *   

 * 
 * @note If you're an experienced programmer, you may be surprised
 * that we use a <code>goto</code> statement in this piece of code!
 * <code>goto</code> statements are not particularly well liked any
 * more since Edsgar Dijkstra, one of the greats of computer science,
 * wrote a letter in 1968 called "Go To Statement considered harmful"
 * (see <a href="http://en.wikipedia.org/wiki/Considered_harmful">here</a>).
 * The author of this code subscribes to this notion whole-heartedly:
 * <code>goto</code> is hard to understand. In fact, deal.II contains
 * virtually no occurrences: excluding code that was essentially
 * transcribed from books and not counting duplicated code pieces,
 * there are 3 locations in about 600,000 lines of code at the time
 * this note is written; we also use it in 4 tutorial programs, in
 * exactly the same context as here. Instead of trying to justify
 * the occurrence here, let's first look at the code and we'll come
 * back to the issue at the end of function.
 * 
 * @code
 *   template <int dim>
 *   void HeatEquation<dim>::run()
 *   {
 *     const unsigned int initial_global_refinement       = 2;
 *     const unsigned int n_adaptive_pre_refinement_steps = 4;
 * 
 *     GridGenerator::hyper_L(triangulation);
 *     triangulation.refine_global(initial_global_refinement);
 * 
 *     setup_system();
 * 
 *     unsigned int pre_refinement_step = 0;
 * 
 *     Vector<double> tmp;
 *     Vector<double> forcing_terms;
 * 
 *   start_time_iteration:
 * 
 *     time            = 0.0;
 *     timestep_number = 0;
 * 
 *     tmp.reinit(solution.size());
 *     forcing_terms.reinit(solution.size());
 * 
 * 
 *     VectorTools::interpolate(dof_handler,
 *                              Functions::ZeroFunction<dim>(),
 *                              old_solution);
 *     solution = old_solution;
 * 
 *     output_results();
 * 
 * @endcode
 * 
 * Then we start the main loop until the computed time exceeds our
 * end time of 0.5. The first task is to build the right hand
 * side of the linear system we need to solve in each time step.
 * Recall that it contains the term $MU^{n-1}-(1-\theta)k_n AU^{n-1}$.
 * We put these terms into the variable system_rhs, with the
 * help of a temporary vector:
 * 
 * @code
 *     while (time <= 0.5)
 *       {
 *         time += time_step;
 *         ++timestep_number;
 * 
 *         std::cout << "Time step " << timestep_number << " at t=" << time
 *                   << std::endl;
 * 
 *         mass_matrix.vmult(system_rhs, old_solution);
 * 
 *         laplace_matrix.vmult(tmp, old_solution);
 *         system_rhs.add(-(1 - theta) * time_step, tmp);
 * 
 * @endcode
 * 
 * The second piece is to compute the contributions of the source
 * terms. This corresponds to the term $k_n
 * \left[ (1-\theta)F^{n-1} + \theta F^n \right]$. The following
 * code calls VectorTools::create_right_hand_side to compute the
 * vectors $F$, where we set the time of the right hand side
 * (source) function before we evaluate it. The result of this
 * all ends up in the forcing_terms variable:
 * 
 * @code
 *         RightHandSide<dim> rhs_function;
 *         rhs_function.set_time(time);
 *         VectorTools::create_right_hand_side(dof_handler,
 *                                             QGauss<dim>(fe.degree + 1),
 *                                             rhs_function,
 *                                             tmp);
 *         forcing_terms = tmp;
 *         forcing_terms *= time_step * theta;
 * 
 *         rhs_function.set_time(time - time_step);
 *         VectorTools::create_right_hand_side(dof_handler,
 *                                             QGauss<dim>(fe.degree + 1),
 *                                             rhs_function,
 *                                             tmp);
 * 
 *         forcing_terms.add(time_step * (1 - theta), tmp);
 * 
 * @endcode
 * 
 * Next, we add the forcing terms to the ones that
 * come from the time stepping, and also build the matrix
 * $M+k_n\theta A$ that we have to invert in each time step.
 * The final piece of these operations is to eliminate
 * hanging node constrained degrees of freedom from the
 * linear system:
 * 
 * @code
 *         system_rhs += forcing_terms;
 * 
 *         system_matrix.copy_from(mass_matrix);
 *         system_matrix.add(theta * time_step, laplace_matrix);
 * 
 *         constraints.condense(system_matrix, system_rhs);
 * 
 * @endcode
 * 
 * There is one more operation we need to do before we
 * can solve it: boundary values. To this end, we create
 * a boundary value object, set the proper time to the one
 * of the current time step, and evaluate it as we have
 * done many times before. The result is used to also
 * set the correct boundary values in the linear system:
 * 
 * @code
 *         {
 *           BoundaryValues<dim> boundary_values_function;
 *           boundary_values_function.set_time(time);
 * 
 *           std::map<types::global_dof_index, double> boundary_values;
 *           VectorTools::interpolate_boundary_values(dof_handler,
 *                                                    0,
 *                                                    boundary_values_function,
 *                                                    boundary_values);
 * 
 *           MatrixTools::apply_boundary_values(boundary_values,
 *                                              system_matrix,
 *                                              solution,
 *                                              system_rhs);
 *         }
 * 
 * @endcode
 * 
 * With this out of the way, all we have to do is solve the
 * system, generate graphical data, and...
 * 
 * @code
 *         solve_time_step();
 * 
 *         output_results();
 * 
 * @endcode
 * 
 * ...take care of mesh refinement. Here, what we want to do is
 * (i) refine the requested number of times at the very beginning
 * of the solution procedure, after which we jump to the top to
 * restart the time iteration, (ii) refine every fifth time
 * step after that.
 *         

 * 
 * The time loop and, indeed, the main part of the program ends
 * with starting into the next time step by setting old_solution
 * to the solution we have just computed.
 * 
 * @code
 *         if ((timestep_number == 1) &&
 *             (pre_refinement_step < n_adaptive_pre_refinement_steps))
 *           {
 *             refine_mesh(initial_global_refinement,
 *                         initial_global_refinement +
 *                           n_adaptive_pre_refinement_steps);
 *             ++pre_refinement_step;
 * 
 *             tmp.reinit(solution.size());
 *             forcing_terms.reinit(solution.size());
 * 
 *             std::cout << std::endl;
 * 
 *             goto start_time_iteration;
 *           }
 *         else if ((timestep_number > 0) && (timestep_number % 5 == 0))
 *           {
 *             refine_mesh(initial_global_refinement,
 *                         initial_global_refinement +
 *                           n_adaptive_pre_refinement_steps);
 *             tmp.reinit(solution.size());
 *             forcing_terms.reinit(solution.size());
 *           }
 * 
 *         old_solution = solution;
 *       }
 *   }
 * } // namespace Step26
 * @endcode
 * 
 * Now that you have seen what the function does, let us come back to the issue
 * of the <code>goto</code>. In essence, what the code does is
 * something like this:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   void run ()
 *   {
 *     initialize;
 *   start_time_iteration:
 *     for (timestep=1...)
 *     {
 *        solve timestep;
 *        if (timestep==1 && not happy with the result)
 *        {
 *          adjust some data structures;
 *          goto start_time_iteration; // simply try again
 *        }
 *        postprocess;
 *     }
 *   }
 * @endcode
 * </div>
 * Here, the condition "happy with the result" is whether we'd like to keep
 * the current mesh or would rather refine the mesh and start over on the
 * new mesh. We could of course replace the use of the <code>goto</code>
 * by the following:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   void run ()
 *   {
 *     initialize;
 *     while (true)
 *     {
 *        solve timestep;
 *        if (not happy with the result)
 *           adjust some data structures;
 *        else
 *           break;
 *     }
 *     postprocess;
 *

 *     for (timestep=2...)
 *     {
 *        solve timestep;
 *        postprocess;
 *     }
 *   }
 * @endcode
 * </div>
 * This has the advantage of getting rid of the <code>goto</code>
 * but the disadvantage of having to duplicate the code that implements
 * the "solve timestep" and "postprocess" operations in two different
 * places. This could be countered by putting these parts of the code
 * (sizable chunks in the actual implementation above) into their
 * own functions, but a <code>while(true)</code> loop with a
 * <code>break</code> statement is not really all that much easier
 * to read or understand than a <code>goto</code>.
 * 

 * 
 * In the end, one might simply agree that <i>in general</i>
 * <code>goto</code> statements are a bad idea but be pragmatic and
 * state that there may be occasions where they can help avoid code
 * duplication and awkward control flow. This may be one of these
 * places, and it matches the position Steve McConnell takes in his
 * excellent book "Code Complete" @cite CodeComplete about good
 * programming practices (see the mention of this book in the
 * introduction of step-1) that spends a surprising ten pages on the
 * question of <code>goto</code> in general.
 * 

 * 
 * 

 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * Having made it this far,  there is, again, nothing
 * much to discuss for the main function of this
 * program: it looks like all such functions since step-6.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step26;
 * 
 *       HeatEquation<2> heat_equation_solver;
 *       heat_equation_solver.run();
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


As in many of the tutorials, the actual output of the program matters less
than how we arrived there. Nonetheless, here it is:
@code
===========================================
Number of active cells: 48
Number of degrees of freedom: 65

Time step 1 at t=0.002
     7 CG iterations.

===========================================
Number of active cells: 60
Number of degrees of freedom: 81


Time step 1 at t=0.002
     7 CG iterations.

===========================================
Number of active cells: 105
Number of degrees of freedom: 136


Time step 1 at t=0.002
     7 CG iterations.

[...]

Time step 249 at t=0.498
     13 CG iterations.
Time step 250 at t=0.5
     14 CG iterations.

===========================================
Number of active cells: 1803
Number of degrees of freedom: 2109
@endcode

Maybe of more interest is a visualization of the solution and the mesh on which
it was computed:

<img src="https://www.dealii.org/images/steps/developer/step-26.movie.gif" alt="Animation of the solution of step 26.">

The movie shows how the two sources switch on and off and how the mesh reacts
to this. It is quite obvious that the mesh as is is probably not the best we
could come up with. We'll get back to this in the next section.


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


There are at least two areas where one can improve this program significantly:
adaptive time stepping and a better choice of the mesh.

<a name="Adaptivetimestepping"></a><h4>Adaptive time stepping</h4>


Having chosen an implicit time stepping scheme, we are not bound by any
CFL-like condition on the time step. Furthermore, because the time scales on
which change happens on a given cell in the heat equation are not bound to the
cells diameter (unlike the case with the wave equation, where we had a fixed
speed of information transport that couples the temporal and spatial scales),
we can choose the time step as we please. Or, better, choose it as we deem
necessary for accuracy.

Looking at the solution, it is clear that the action does not happen uniformly
over time: a lot is changing around the time we switch on a source, things
become less dramatic once a source is on for a little while, and we enter a
long phase of decline when both sources are off. During these times, we could
surely get away with a larger time step than before without sacrificing too
much accuracy.

The literature has many suggestions on how to choose the time step size
adaptively. Much can be learned, for example, from the way ODE solvers choose
their time steps. One can also be inspired by a posteriori error estimators
that can, ideally, be written in a way that the consist of a temporal and a
spatial contribution to the overall error. If the temporal one is too large,
we should choose a smaller time step. Ideas in this direction can be found,
for example, in the PhD thesis of a former principal developer of deal.II,
Ralf Hartmann, published by the University of Heidelberg, Germany, in 2002.


<a name="Bettertimesteppingmethods"></a><h4>Better time stepping methods</h4>


We here use one of the simpler time stepping methods, namely the second order
in time Crank-Nicolson method. However, more accurate methods such as
Runge-Kutta methods are available and should be used as they do not represent
much additional effort. It is not difficult to implement this for the current
program, but a more systematic treatment is also given in step-52.


<a name="Betterrefinementcriteria"></a><h4>Better refinement criteria</h4>


If you look at the meshes in the movie above, it is clear that they are not
particularly well suited to the task at hand. In fact, they look rather
random.

There are two factors at play. First, there are some islands where cells
have been refined but that are surrounded by non-refined cells (and there
are probably also a few occasional coarsened islands). These are not terrible,
as they most of the time do not affect the approximation quality of the mesh,
but they also don't help because so many of their additional degrees of
freedom are in fact constrained by hanging node constraints. That said,
this is easy to fix: the Triangulation class takes an argument to its
constructor indicating a level of "mesh smoothing". Passing one of many
possible flags, this instructs the triangulation to refine some additional
cells, or not to refine some cells, so that the resulting mesh does not have
these artifacts.

The second problem is more severe: the mesh appears to lag the solution.
The underlying reason is that we only adapt the mesh once every fifth
time step, and only allow for a single refinement in these cases. Whenever a
source switches on, the solution had been very smooth in this area before and
the mesh was consequently rather coarse. This implies that the next time step
when we refine the mesh, we will get one refinement level more in this area,
and five time steps later another level, etc. But this is not enough: first,
we should refine immediately when a source switches on (after all, in the
current context we at least know what the right hand side is), and we should
allow for more than one refinement level. Of course, all of this can be done
using deal.II, it just requires a bit of algorithmic thinking in how to make
this work!


<a name="Positivitypreservation"></a><h4>Positivity preservation</h4>


To increase the accuracy and resolution of your simulation in time, one
typically decreases the time step size $k_n$. If you start playing around
with the time step in this particular example, you will notice that the
solution becomes partly negative, if $k_n$ is below a certain threshold.
This is not what we would expect to happen (in nature).

To get an idea of this behavior mathematically, let us consider a general,
fully discrete problem:
@f{align*}
  A u^{n} = B u^{n-1}.
@f}
The general form of the $i$th equation then reads:
@f{align*}
  a_{ii} u^{n}_i &= b_{ii} u^{n-1}_i +
  \sum\limits_{j \in S_i} \left( b_{ij} u^{n-1}_j - a_{ij} u^{n}_j \right),
@f}
where $S_i$ is the set of degrees of freedom that DoF $i$ couples with (i.e.,
for which either the matrix $A$ or matrix $B$ has a nonzero entry at position
$(i,j)$). If all coefficients
fulfill the following conditions:
@f{align*}
  a_{ii} &> 0, & b_{ii} &\geq 0, & a_{ij} &\leq 0, & b_{ij} &\geq 0,
  &
  \forall j &\in S_i,
@f}
all solutions $u^{n}$ keep their sign from the previous ones $u^{n-1}$, and
consequently from the initial values $u^0$. See e.g.
<a href="http://bookstore.siam.org/cs14/">Kuzmin, H&auml;m&auml;l&auml;inen</a>
for more information on positivity preservation.

Depending on the PDE to solve and the time integration scheme used, one is
able to deduce conditions for the time step $k_n$. For the heat equation with
the Crank-Nicolson scheme,
<a href="https://doi.org/10.2478/cmam-2010-0025">Schatz et. al.</a> have
translated it to the following ones:
@f{align*}
  (1 - \theta) k a_{ii} &\leq m_{ii},\qquad \forall i,
  &
  \theta k \left| a_{ij} \right| &\geq m_{ij},\qquad j \neq i,
@f}
where $M = m_{ij}$ denotes the mass matrix and $A = a_{ij}$ the stiffness
matrix with $a_{ij} \leq 0$ for $j \neq i$, respectively. With
$a_{ij} \leq 0$, we can formulate bounds for the global time step $k$ as
follows:
@f{align*}
  k_{\text{max}} &= \frac{ 1 }{ 1 - \theta }
  \min\left( \frac{ m_{ii} }{ a_{ii} } \right),~ \forall i,
  &
  k_{\text{min}} &= \frac{ 1 }{ \theta  }
  \max\left( \frac{ m_{ij} }{ \left|a_{ij}\right| } \right),~ j \neq i.
@f}
In other words, the time step is constrained by <i>both a lower
and upper bound</i> in case of a Crank-Nicolson scheme. These bounds should be
considered along with the CFL condition to ensure significance of the performed
simulations.

Being unable to make the time step as small as we want to get more
accuracy without losing the positivity property is annoying. It raises
the question of whether we can at least <i>compute</i> the minimal time step
we can choose  to ensure positivity preservation in this particular tutorial.
Indeed, we can use
the SparseMatrix objects for both mass and stiffness that are created via
the MatrixCreator functions. Iterating through each entry via SparseMatrixIterators
lets us check for diagonal and off-diagonal entries to set a proper time step
dynamically. For quadratic matrices, the diagonal element is stored as the
first member of a row (see SparseMatrix documentation). An exemplary code
snippet on how to grab the entries of interest from the <code>mass_matrix</code>
is shown below.

@code
Assert (mass_matrix.m() == mass_matrix.n(), ExcNotQuadratic());
const unsigned int num_rows = mass_matrix.m();
double mass_matrix_min_diag    = std::numeric_limits<double>::max(),
       mass_matrix_max_offdiag = 0.;

SparseMatrixIterators::Iterator<double,true> row_it (&mass_matrix, 0);

for(unsigned int m = 0; m<num_rows; ++m)
{
  // check the diagonal element
  row_it = mass_matrix.begin(m);
  mass_matrix_min_diag = std::min(row_it->value(), mass_matrix_min_diag);
  ++row_it;

  // check the off-diagonal elements
  for(; row_it != mass_matrix.end(m); ++row_it)
    mass_matrix_max_offdiag = std::max(row_it->value(), mass_matrix_max_offdiag);
}
@endcode

Using the information so computed, we can bound the time step via the formulas
above.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-26.cc"
*/
