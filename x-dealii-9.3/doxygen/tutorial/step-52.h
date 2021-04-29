/**
@page step_52 The step-52 tutorial program
This tutorial depends on step-26.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemstatement">Problem statement</a>
        <li><a href="#RungeKuttamethods">Runge-Kutta methods</a>
      <ul>
        <li><a href="#ExplicitRungeKuttamethods">Explicit Runge-Kutta methods</a>
        <li><a href="#EmbeddedRungeKuttamethods">Embedded Runge-Kutta methods</a>
        <li><a href="#ImplicitRungeKuttamethods">Implicit Runge-Kutta methods</a>
      </ul>
        <li><a href="#Spatiallydiscreteformulation">Spatially discrete formulation</a>
        <li><a href="#Notesonthetestcase">Notes on the testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeDiffusioncodeclass">The <code>Diffusion</code> class</a>
      <ul>
        <li><a href="#codeDiffusionsetup_systemcode"><code>Diffusion::setup_system</code></a>
        <li><a href="#codeDiffusionassemble_systemcode"><code>Diffusion::assemble_system</code></a>
        <li><a href="#codeDiffusionget_sourcecode"><code>Diffusion::get_source</code></a>
        <li><a href="#codeDiffusionevaluate_diffusioncode"><code>Diffusion::evaluate_diffusion</code></a>
        <li><a href="#codeDiffusionid_minus_tau_J_inversecode"><code>Diffusion::id_minus_tau_J_inverse</code></a>
        <li><a href="#codeDiffusionoutput_resultscode"><code>Diffusion::output_results</code></a>
        <li><a href="#codeDiffusionexplicit_methodcode"><code>Diffusion::explicit_method</code></a>
        <li><a href="#codeDiffusionimplicit_methodcode"><code>Diffusion::implicit_method</code></a>
        <li><a href="#codeDiffusionembedded_explicit_methodcode"><code>Diffusion::embedded_explicit_method</code></a>
        <li><a href="#codeDiffusionruncode"><code>Diffusion::run</code></a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Bruno Turcksin and Damien Lebrun-Grandie.</i>

@note In order to run this program, deal.II must be configured to use
the UMFPACK sparse direct solver. Refer to the <a
href="../../readme.html#umfpack">ReadMe</a> for instructions how to do this.

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This program shows how to use Runge-Kutta methods to solve a time-dependent
problem. It solves a small variation of the heat equation discussed first in
step-26 but, since the purpose of this program is only to demonstrate using
more advanced ways to interface with deal.II's time stepping algorithms, only
solves a simple problem on a uniformly refined mesh.


<a name="Problemstatement"></a><h3>Problem statement</h3>


In this example, we solve the one-group time-dependent diffusion
approximation of the neutron transport equation (see step-28 for the
time-independent multigroup diffusion). This is a model for how neutrons move
around highly scattering media, and consequently it is a variant of the
time-dependent diffusion equation -- which is just a different name for the
heat equation discussed in step-26, plus some extra terms.
We assume that the medium is not
fissible and therefore, the neutron flux satisfies the following equation:
@f{eqnarray*}
\frac{1}{v}\frac{\partial \phi(x,t)}{\partial t} = \nabla \cdot D(x) \nabla \phi(x,t)
- \Sigma_a(x) \phi(x,t) + S(x,t)
@f}
augmented by appropriate boundary conditions. Here, $v$ is the velocity of
neutrons (for simplicity we assume it is equal to 1 which can be achieved by
simply scaling the time variable), $D$ is the diffusion coefficient,
$\Sigma_a$ is the absorption cross section, and $S$ is a source. Because we are
only interested in the time dependence, we assume that $D$ and $\Sigma_a$ are
constant.

Since this program only intends to demonstrate how to use advanced time
stepping algorithms, we will only look for the solutions of relatively simple
problems. Specifically, we are looking for a solution on a square domain
$[0,b]\times[0,b]$ of the form
@f{eqnarray*}
\phi(x,t) = A\sin(\omega t)(bx-x^2).
@f}
By using quadratic finite elements, we can represent this function exactly at
any particular time, and all the error will be due to the time
discretization. We do this because it is then easy to observe the order of
convergence of the various time stepping schemes we will consider, without
having to separate spatial and temporal errors.

We impose the following boundary conditions: homogeneous Dirichlet for $x=0$ and
$x=b$ and homogeneous Neumann conditions for $y=0$ and $y=b$. We choose the
source term so that the corresponding solution is
in fact of the form stated above:
@f{eqnarray*}
S=A\left(\frac{1}{v}\omega \cos(\omega t)(bx -x^2) + \sin(\omega t)
\left(\Sigma_a (bx-x^2)+2D\right) \right).
@f}
Because the solution is a sine in time, we know that the exact solution
satisfies $\phi\left(x,\frac{\pi}{\omega}\right) = 0$.
Therefore, the error at time $t=\frac{\pi}{\omega}$ is simply the norm of the numerical
solution, i.e., $\|e(\cdot,t=\frac{\pi}{\omega})\|_{L_2} = \|\phi_h(\cdot,t=\frac{\pi}{\omega})\|_{L_2}$,
and is particularly easily evaluated. In the code, we evaluate the $l_2$ norm
of the vector of nodal values of $\phi_h$ instead of the $L_2$ norm of the
associated spatial function, since the former is simpler to compute; however,
on uniform meshes, the two are just related by a constant and we can
consequently observe the temporal convergence order with either.


<a name="RungeKuttamethods"></a><h3>Runge-Kutta methods</h3>


The Runge-Kutta methods implemented in deal.II assume that the equation to be
solved can be written as:
@f{eqnarray*}
\frac{dy}{dt} = g(t,y).
@f}
On the other hand, when using finite elements, discretized time derivatives always result in the
presence of a mass matrix on the left hand side. This can easily be seen by
considering that if the solution vector $y(t)$ in the equation above is in fact the vector
of nodal coefficients $U(t)$ for a variable of the form
@f{eqnarray*}
  u_h(x,t) = \sum_j U_j(t) \varphi_j(x)
@f}
with spatial shape functions $\varphi_j(x)$, then multiplying an equation of
the form
@f{eqnarray*}
  \frac{\partial u(x,t)}{\partial t} = q(t,u(x,t))
@f}
by test functions, integrating over $\Omega$, substituting $u\rightarrow u_h$
and restricting the test functions to the $\varphi_i(x)$ from above, then this
spatially discretized equation has the form
@f{eqnarray*}
M\frac{dU}{dt} = f(t,U),
@f}
where $M$ is the mass matrix and $f(t,U)$ is the spatially discretized version
of $q(t,u(x,t))$ (where $q$ is typically the place where spatial
derivatives appear, but this is not of much concern for the moment given that
we only consider time derivatives). In other words, this form fits the general
scheme above if we write
@f{eqnarray*}
\frac{dy}{dt} = g(t,y) = M^{-1}f(t,y).
@f}

Runke-Kutta methods are time stepping schemes that approximate $y(t_n)\approx
y_{n}$ through a particular one-step approach. They are typically written in the form
@f{eqnarray*}
y_{n+1} = y_n + \sum_{i=1}^s b_i k_i
@f}
where for the form of the right hand side above
@f{eqnarray*}
k_i = h M^{-1} f\left(t_n+c_ih,y_n+\sum_{j=1}^sa_{ij}k_j\right).
@f}
Here $a_{ij}$, $b_i$, and $c_i$ are known coefficients that identify which
particular Runge-Kutta scheme you want to use, and $h=t_{n+1}-t_n$ is the time step
used. Different time stepping methods of the Runge-Kutta class differ in the
number of stages $s$ and the values they use for the coefficients $a_{ij}$,
$b_i$, and $c_i$ but are otherwise easy to implement since one can look up
tabulated values for these coefficients. (These tables are often called
Butcher tableaus.)

At the time of the writing of this tutorial, the methods implemented in
deal.II can be divided in three categories:
<ol>
<li> Explicit Runge-Kutta; in order for a method to be explicit, it is
necessary that in the formula above defining $k_i$, $k_i$ does not appear
on the right hand side. In other words, these methods have to satisfy
$a_{ii}=0, i=1,\ldots,s$.
<li> Embedded (or adaptive) Runge-Kutta; we will discuss their properties below.
<li> Implicit Runge-Kutta; this class of methods require the solution of a
possibly nonlinear system the stages $k_i$ above, i.e., they have
$a_{ii}\neq 0$ for at least one of the stages $i=1,\ldots,s$.
</ol>
Many well known time stepping schemes that one does not typically associate
with the names Runge or Kutta can in fact be written in a way so that they,
too, can be expressed in these categories. They oftentimes represent the
lowest-order members of these families.


<a name="ExplicitRungeKuttamethods"></a><h4>Explicit Runge-Kutta methods</h4>


These methods, only require a function to evaluate $M^{-1}f(t,y)$ but not
(as implicit methods) to solve an equation that involves
$f(t,y)$ for $y$. As all explicit time stepping methods, they become unstable
when the time step chosen is too large.

Well known methods in this class include forward Euler, third order
Runge-Kutta, and fourth order Runge-Kutta (often abbreviated as RK4).


<a name="EmbeddedRungeKuttamethods"></a><h4>Embedded Runge-Kutta methods</h4>


These methods use both a lower and a higher order method to
estimate the error and decide if the time step needs to be shortened or can be
increased. The term "embedded" refers to the fact that the lower-order method
does not require additional evaluates of the function $M^{-1}f(\cdot,\cdot)$
but reuses data that has to be computed for the high order method anyway. It
is, in other words, essentially free, and we get the error estimate as a side
product of using the higher order method.

This class of methods include Heun-Euler, Bogacki-Shampine, Dormand-Prince (ode45 in
Matlab and often abbreviated as RK45 to indicate that the lower and higher order methods
used here are 4th and 5th order Runge-Kutta methods, respectively), Fehlberg,
and Cash-Karp.

At the time of the writing, only embedded explicit methods have been implemented.


<a name="ImplicitRungeKuttamethods"></a><h4>Implicit Runge-Kutta methods</h4>


Implicit methods require the solution of (possibly nonlinear) systems of the
form $\alpha y = f(t,y)$
for $y$ in each (sub-)timestep. Internally, this is
done using a Newton-type method and, consequently, they require that the user
provide functions that can evaluate $M^{-1}f(t,y)$ and
$\left(I-\tau M^{-1} \frac{\partial f}{\partial y}\right)^{-1}$ or equivalently
$\left(M - \tau \frac{\partial f}{\partial y}\right)^{-1} M$.

The particular form of this operator results from the fact that each Newton
step requires the solution of an equation of the form
@f{align*}
  \left(M - \tau \frac{\partial f}{\partial y}\right) \Delta y
  = -M h(t,y)
@f}
for some (given) $h(t,y)$. Implicit methods are
always stable, regardless of the time step size, but too large time steps of
course affect the <i>accuracy</i> of the solution, even if the numerical
solution remains stable and bounded.

Methods in this class include backward Euler, implicit midpoint,
Crank-Nicolson, and the two stage SDIRK method (short for "singly diagonally
implicit Runge-Kutta", a term coined to indicate that the diagonal elements
$a_{ii}$ defining the time stepping method are all equal; this property
allows for the Newton matrix $I-\tau M^{-1}\frac{\partial f}{\partial y}$ to
be re-used between stages because $\tau$ is the same every time).


<a name="Spatiallydiscreteformulation"></a><h3>Spatially discrete formulation</h3>


By expanding the solution of our model problem
as always using shape functions $\psi_j$ and writing
@f{eqnarray*}
\phi_h(x,t) = \sum_j U_j(t) \psi_j(x),
@f}
we immediately get the spatially discretized version of the diffusion equation as
@f{eqnarray*}
  M \frac{dU(t)}{dt}
  = -{\cal D} U(t) - {\cal A} U(t) + {\cal S}(t)
@f}
where
@f{eqnarray*}
  M_{ij}  &=& (\psi_i,\psi_j), \\
  {\cal D}_{ij}  &=& (D\nabla\psi_i,\nabla\psi_j)_\Omega, \\
  {\cal A}_{ij}  &=& (\Sigma_a\psi_i,\psi_j)_\Omega, \\
  {\cal S}_{i}(t)  &=& (\psi_i,S(x,t))_\Omega.
@f}
See also step-24 and step-26 to understand how we arrive here.
Boundary terms are not necessary due to the chosen boundary conditions for
the current problem. To use the Runge-Kutta methods, we recast this
as follows:
@f{eqnarray*}
f(y) = -{\cal D}y - {\cal A}y + {\cal S}.
@f}
In the code, we will need to be able to evaluate this function $f(U)$ along
with its derivative,
@f{eqnarray*}
\frac{\partial f}{\partial y} = -{\cal D} - {\cal A}.
@f}


<a name="Notesonthetestcase"></a><h3>Notes on the testcase</h3>


To simplify the problem, the domain is two dimensional and the mesh is
uniformly refined (there is no need to adapt the mesh since we use quadratic
finite elements and the exact solution is quadratic). Going from a two
dimensional domain to a three dimensional domain is not very
challenging. However if you intend to solve more complex problems where the
mesh must be adapted (as is done, for example, in step-26), then it is
important to remember the following issues:

<ol>
<li> You will need to project the solution to the new mesh when the mesh is changed. Of course,
     the mesh
     used should be the same from the beginning to the end of each time step,
     a question that arises because Runge-Kutta methods use multiple
     evaluations of the equations within each time step.
<li> You will need to update the mass matrix and its inverse every time the
     mesh is changed.
</ol>
The techniques for these steps are readily available by looking at step-26.
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
 * The first task as usual is to include the functionality of these well-known
 * deal.II library files and some C++ header files.
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_out.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparse_direct.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * #include <cmath>
 * #include <map>
 * 
 * @endcode
 * 
 * This is the only include file that is new: It includes all the Runge-Kutta
 * methods.
 * 
 * @code
 * #include <deal.II/base/time_stepping.h>
 * 
 * 
 * @endcode
 * 
 * The next step is like in all previous tutorial programs: We put everything
 * into a namespace of its own and then import the deal.II classes and functions
 * into it.
 * 
 * @code
 * namespace Step52
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeDiffusioncodeclass"></a> 
 * <h3>The <code>Diffusion</code> class</h3>
 * 

 * 
 * The next piece is the declaration of the main class. Most of the
 * functions in this class are not new and have been explained in previous
 * tutorials. The only interesting functions are
 * <code>evaluate_diffusion()</code> and
 * <code>id_minus_tau_J_inverse()</code>. <code>evaluate_diffusion()</code>
 * evaluates the diffusion equation, $M^{-1}(f(t,y))$, at a given time and a
 * given $y$. <code>id_minus_tau_J_inverse()</code> evaluates $\left(I-\tau
 * M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ or equivalently
 * $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$ at a given
 * time, for a given $\tau$ and $y$. This function is needed when an
 * implicit method is used.
 * 
 * @code
 *   class Diffusion
 *   {
 *   public:
 *     Diffusion();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 * 
 *     void assemble_system();
 * 
 *     double get_source(const double time, const Point<2> &point) const;
 * 
 *     Vector<double> evaluate_diffusion(const double          time,
 *                                       const Vector<double> &y) const;
 * 
 *     Vector<double> id_minus_tau_J_inverse(const double          time,
 *                                           const double          tau,
 *                                           const Vector<double> &y);
 * 
 *     void output_results(const double                     time,
 *                         const unsigned int               time_step,
 *                         TimeStepping::runge_kutta_method method) const;
 * 
 * @endcode
 * 
 * The next three functions are the drivers for the explicit methods, the
 * implicit methods, and the embedded explicit methods respectively. The
 * driver function for embedded explicit methods returns the number of
 * steps executed given that it only takes the number of time steps passed
 * as an argument as a hint, but internally computed the optimal time step
 * itself.
 * 
 * @code
 *     void explicit_method(const TimeStepping::runge_kutta_method method,
 *                          const unsigned int                     n_time_steps,
 *                          const double                           initial_time,
 *                          const double                           final_time);
 * 
 *     void implicit_method(const TimeStepping::runge_kutta_method method,
 *                          const unsigned int                     n_time_steps,
 *                          const double                           initial_time,
 *                          const double                           final_time);
 * 
 *     unsigned int
 *     embedded_explicit_method(const TimeStepping::runge_kutta_method method,
 *                              const unsigned int n_time_steps,
 *                              const double       initial_time,
 *                              const double       final_time);
 * 
 * 
 *     const unsigned int fe_degree;
 * 
 *     const double diffusion_coefficient;
 *     const double absorption_cross_section;
 * 
 *     Triangulation<2> triangulation;
 * 
 *     const FE_Q<2> fe;
 * 
 *     DoFHandler<2> dof_handler;
 * 
 *     AffineConstraints<double> constraint_matrix;
 * 
 *     SparsityPattern sparsity_pattern;
 * 
 *     SparseMatrix<double> system_matrix;
 *     SparseMatrix<double> mass_matrix;
 *     SparseMatrix<double> mass_minus_tau_Jacobian;
 * 
 *     SparseDirectUMFPACK inverse_mass_matrix;
 * 
 *     Vector<double> solution;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * We choose quadratic finite elements and we initialize the parameters.
 * 
 * @code
 *   Diffusion::Diffusion()
 *     : fe_degree(2)
 *     , diffusion_coefficient(1. / 30.)
 *     , absorption_cross_section(1.)
 *     , fe(fe_degree)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionsetup_systemcode"></a> 
 * <h4><code>Diffusion::setup_system</code></h4>
 * Now, we create the constraint matrix and the sparsity pattern. Then, we
 * initialize the matrices and the solution vector.
 * 
 * @code
 *   void Diffusion::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              1,
 *                                              Functions::ZeroFunction<2>(),
 *                                              constraint_matrix);
 *     constraint_matrix.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraint_matrix);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *     mass_matrix.reinit(sparsity_pattern);
 *     mass_minus_tau_Jacobian.reinit(sparsity_pattern);
 *     solution.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionassemble_systemcode"></a> 
 * <h4><code>Diffusion::assemble_system</code></h4>
 * In this function, we compute $-\int D \nabla b_i \cdot \nabla b_j
 * d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$ and the mass
 * matrix $\int b_i b_j d\boldsymbol{r}$. The mass matrix is then
 * inverted using a direct solver; the <code>inverse_mass_matrix</code>
 * variable will then store the inverse of the mass matrix so that
 * $M^{-1}$ can be applied to a vector using the <code>vmult()</code>
 * function of that object. (Internally, UMFPACK does not really store
 * the inverse of the matrix, but its LU factors; applying the inverse
 * matrix is then equivalent to doing one forward and one backward solves
 * with these two factors, which has the same complexity as applying an
 * explicit inverse of the matrix).
 * 
 * @code
 *   void Diffusion::assemble_system()
 *   {
 *     system_matrix = 0.;
 *     mass_matrix   = 0.;
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1);
 * 
 *     FEValues<2> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients | update_JxW_values);
 * 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix      = 0.;
 *         cell_mass_matrix = 0.;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               {
 *                 cell_matrix(i, j) +=
 *                   ((-diffusion_coefficient *                // (-D
 *                       fe_values.shape_grad(i, q_point) *    //  * grad phi_i
 *                       fe_values.shape_grad(j, q_point)      //  * grad phi_j
 *                     - absorption_cross_section *            //  -Sigma
 *                         fe_values.shape_value(i, q_point) * //  * phi_i
 *                         fe_values.shape_value(j, q_point))  //  * phi_j)
 *                    * fe_values.JxW(q_point));               // * dx
 *                 cell_mass_matrix(i, j) += fe_values.shape_value(i, q_point) *
 *                                           fe_values.shape_value(j, q_point) *
 *                                           fe_values.JxW(q_point);
 *               }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraint_matrix.distribute_local_to_global(cell_matrix,
 *                                                      local_dof_indices,
 *                                                      system_matrix);
 *         constraint_matrix.distribute_local_to_global(cell_mass_matrix,
 *                                                      local_dof_indices,
 *                                                      mass_matrix);
 *       }
 * 
 *     inverse_mass_matrix.initialize(mass_matrix);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionget_sourcecode"></a> 
 * <h4><code>Diffusion::get_source</code></h4>
 *   

 * 
 * In this function, the source term of the equation for a given time and a
 * given point is computed.
 * 
 * @code
 *   double Diffusion::get_source(const double time, const Point<2> &point) const
 *   {
 *     const double intensity = 10.;
 *     const double frequency = numbers::PI / 10.;
 *     const double b         = 5.;
 *     const double x         = point(0);
 * 
 *     return intensity *
 *            (frequency * std::cos(frequency * time) * (b * x - x * x) +
 *             std::sin(frequency * time) *
 *               (absorption_cross_section * (b * x - x * x) +
 *                2. * diffusion_coefficient));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionevaluate_diffusioncode"></a> 
 * <h4><code>Diffusion::evaluate_diffusion</code></h4>
 *   

 * 
 * Next, we evaluate the weak form of the diffusion equation at a given time
 * $t$ and for a given vector $y$. In other words, as outlined in the
 * introduction, we evaluate $M^{-1}(-{\cal D}y - {\cal A}y + {\cal
 * S})$. For this, we have to apply the matrix $-{\cal D} - {\cal A}$
 * (previously computed and stored in the variable
 * <code>system_matrix</code>) to $y$ and then add the source term which we
 * integrate as we usually do. (Integrating up the solution could be done
 * using VectorTools::create_right_hand_side() if you wanted to save a few
 * lines of code, or wanted to take advantage of doing the integration in
 * parallel.) The result is then multiplied by $M^{-1}$.
 * 
 * @code
 *   Vector<double> Diffusion::evaluate_diffusion(const double          time,
 *                                                const Vector<double> &y) const
 *   {
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     tmp = 0.;
 *     system_matrix.vmult(tmp, y);
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1);
 * 
 *     FEValues<2> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_quadrature_points |
 *                             update_JxW_values);
 * 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double> cell_source(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_source = 0.;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           {
 *             const double source =
 *               get_source(time, fe_values.quadrature_point(q_point));
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_source(i) += fe_values.shape_value(i, q_point) * // phi_i(x)
 *                                 source *                            // * S(x)
 *                                 fe_values.JxW(q_point);             // * dx
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraint_matrix.distribute_local_to_global(cell_source,
 *                                                      local_dof_indices,
 *                                                      tmp);
 *       }
 * 
 *     Vector<double> value(dof_handler.n_dofs());
 *     inverse_mass_matrix.vmult(value, tmp);
 * 
 *     return value;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionid_minus_tau_J_inversecode"></a> 
 * <h4><code>Diffusion::id_minus_tau_J_inverse</code></h4>
 *   

 * 
 * We compute $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$. This
 * is done in several steps:
 * - compute $M-\tau \frac{\partial f}{\partial y}$
 * - invert the matrix to get $\left(M-\tau \frac{\partial f}
 * {\partial y}\right)^{-1}$
 * - compute $tmp=My$
 * - compute $z=\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} tmp =
 * \left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} My$
 * - return z.
 * 
 * @code
 *   Vector<double> Diffusion::id_minus_tau_J_inverse(const double /*time*/,
 *                                                    const double          tau,
 *                                                    const Vector<double> &y)
 *   {
 *     SparseDirectUMFPACK inverse_mass_minus_tau_Jacobian;
 * 
 *     mass_minus_tau_Jacobian.copy_from(mass_matrix);
 *     mass_minus_tau_Jacobian.add(-tau, system_matrix);
 * 
 *     inverse_mass_minus_tau_Jacobian.initialize(mass_minus_tau_Jacobian);
 * 
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     mass_matrix.vmult(tmp, y);
 * 
 *     Vector<double> result(y);
 *     inverse_mass_minus_tau_Jacobian.vmult(result, tmp);
 * 
 *     return result;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionoutput_resultscode"></a> 
 * <h4><code>Diffusion::output_results</code></h4>
 *   

 * 
 * The following function then outputs the solution in vtu files indexed by
 * the number of the time step and the name of the time stepping method. Of
 * course, the (exact) result should really be the same for all time
 * stepping method, but the output here at least allows us to compare them.
 * 
 * @code
 *   void Diffusion::output_results(const double                     time,
 *                                  const unsigned int               time_step,
 *                                  TimeStepping::runge_kutta_method method) const
 *   {
 *     std::string method_name;
 * 
 *     switch (method)
 *       {
 *         case TimeStepping::FORWARD_EULER:
 *           {
 *             method_name = "forward_euler";
 *             break;
 *           }
 *         case TimeStepping::RK_THIRD_ORDER:
 *           {
 *             method_name = "rk3";
 *             break;
 *           }
 *         case TimeStepping::RK_CLASSIC_FOURTH_ORDER:
 *           {
 *             method_name = "rk4";
 *             break;
 *           }
 *         case TimeStepping::BACKWARD_EULER:
 *           {
 *             method_name = "backward_euler";
 *             break;
 *           }
 *         case TimeStepping::IMPLICIT_MIDPOINT:
 *           {
 *             method_name = "implicit_midpoint";
 *             break;
 *           }
 *         case TimeStepping::SDIRK_TWO_STAGES:
 *           {
 *             method_name = "sdirk";
 *             break;
 *           }
 *         case TimeStepping::HEUN_EULER:
 *           {
 *             method_name = "heun_euler";
 *             break;
 *           }
 *         case TimeStepping::BOGACKI_SHAMPINE:
 *           {
 *             method_name = "bocacki_shampine";
 *             break;
 *           }
 *         case TimeStepping::DOPRI:
 *           {
 *             method_name = "dopri";
 *             break;
 *           }
 *         case TimeStepping::FEHLBERG:
 *           {
 *             method_name = "fehlberg";
 *             break;
 *           }
 *         case TimeStepping::CASH_KARP:
 *           {
 *             method_name = "cash_karp";
 *             break;
 *           }
 *         default:
 *           {
 *             break;
 *           }
 *       }
 * 
 *     DataOut<2> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 * 
 *     data_out.build_patches();
 * 
 *     data_out.set_flags(DataOutBase::VtkFlags(time, time_step));
 * 
 *     const std::string filename = "solution_" + method_name + "-" +
 *                                  Utilities::int_to_string(time_step, 3) +
 *                                  ".vtu";
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 * 
 *     static std::vector<std::pair<double, std::string>> times_and_names;
 * 
 *     static std::string method_name_prev = "";
 *     static std::string pvd_filename;
 *     if (method_name_prev != method_name)
 *       {
 *         times_and_names.clear();
 *         method_name_prev = method_name;
 *         pvd_filename     = "solution_" + method_name + ".pvd";
 *       }
 *     times_and_names.emplace_back(time, filename);
 *     std::ofstream pvd_output(pvd_filename);
 *     DataOutBase::write_pvd_record(pvd_output, times_and_names);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionexplicit_methodcode"></a> 
 * <h4><code>Diffusion::explicit_method</code></h4>
 *   

 * 
 * This function is the driver for all the explicit methods. At the
 * top it initializes the time stepping and the solution (by setting
 * it to zero and then ensuring that boundary value and hanging node
 * constraints are respected; of course, with the mesh we use here,
 * hanging node constraints are not in fact an issue). It then calls
 * <code>evolve_one_time_step</code> which performs one time step.
 * Time is stored and incremented through a DiscreteTime object.
 *   

 * 
 * For explicit methods, <code>evolve_one_time_step</code> needs to
 * evaluate $M^{-1}(f(t,y))$, i.e, it needs
 * <code>evaluate_diffusion</code>. Because
 * <code>evaluate_diffusion</code> is a member function, it needs to
 * be bound to <code>this</code>. After each evolution step, we
 * again apply the correct boundary values and hanging node
 * constraints.
 *   

 * 
 * Finally, the solution is output
 * every 10 time steps.
 * 
 * @code
 *   void Diffusion::explicit_method(const TimeStepping::runge_kutta_method method,
 *                                   const unsigned int n_time_steps,
 *                                   const double       initial_time,
 *                                   const double       final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::ExplicitRungeKutta<Vector<double>> explicit_runge_kutta(
 *       method);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         explicit_runge_kutta.evolve_one_time_step(
 *           [this](const double time, const Vector<double> &y) {
 *             return this->evaluate_diffusion(time, y);
 *           },
 *           time.get_current_time(),
 *           time.get_next_step_size(),
 *           solution);
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionimplicit_methodcode"></a> 
 * <h4><code>Diffusion::implicit_method</code></h4>
 * This function is equivalent to <code>explicit_method</code> but for
 * implicit methods. When using implicit methods, we need to evaluate
 * $M^{-1}(f(t,y))$ and $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial
 * y}\right)^{-1}$ for which we use the two member functions previously
 * introduced.
 * 
 * @code
 *   void Diffusion::implicit_method(const TimeStepping::runge_kutta_method method,
 *                                   const unsigned int n_time_steps,
 *                                   const double       initial_time,
 *                                   const double       final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::ImplicitRungeKutta<Vector<double>> implicit_runge_kutta(
 *       method);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         implicit_runge_kutta.evolve_one_time_step(
 *           [this](const double time, const Vector<double> &y) {
 *             return this->evaluate_diffusion(time, y);
 *           },
 *           [this](const double time, const double tau, const Vector<double> &y) {
 *             return this->id_minus_tau_J_inverse(time, tau, y);
 *           },
 *           time.get_current_time(),
 *           time.get_next_step_size(),
 *           solution);
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionembedded_explicit_methodcode"></a> 
 * <h4><code>Diffusion::embedded_explicit_method</code></h4>
 * This function is the driver for the embedded explicit methods. It requires
 * more parameters:
 * - coarsen_param: factor multiplying the current time step when the error
 * is below the threshold.
 * - refine_param: factor multiplying the current time step when the error
 * is above the threshold.
 * - min_delta: smallest time step acceptable.
 * - max_delta: largest time step acceptable.
 * - refine_tol: threshold above which the time step is refined.
 * - coarsen_tol: threshold below which the time step is coarsen.
 *   

 * 
 * Embedded methods use a guessed time step. If the error using this time step
 * is too large, the time step will be reduced. If the error is below the
 * threshold, a larger time step will be tried for the next time step.
 * <code>delta_t_guess</code> is the guessed time step produced by the
 * embedded method. In summary, time step size is potentially modified in
 * three ways:
 * - Reducing or increasing time step size within
 * TimeStepping::EmbeddedExplicitRungeKutta::evolve_one_time_step().
 * - Using the calculated <code>delta_t_guess</code>.
 * - Automatically adjusting the step size of the last time step to ensure
 * simulation ends precisely at <code>final_time</code>. This adjustment
 * is handled inside the DiscreteTime instance.
 * 
 * @code
 *   unsigned int Diffusion::embedded_explicit_method(
 *     const TimeStepping::runge_kutta_method method,
 *     const unsigned int                     n_time_steps,
 *     const double                           initial_time,
 *     const double                           final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 *     const double coarsen_param = 1.2;
 *     const double refine_param  = 0.8;
 *     const double min_delta     = 1e-8;
 *     const double max_delta     = 10 * time_step;
 *     const double refine_tol    = 1e-1;
 *     const double coarsen_tol   = 1e-5;
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>>
 *       embedded_explicit_runge_kutta(method,
 *                                     coarsen_param,
 *                                     refine_param,
 *                                     min_delta,
 *                                     max_delta,
 *                                     refine_tol,
 *                                     coarsen_tol);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         const double new_time =
 *           embedded_explicit_runge_kutta.evolve_one_time_step(
 *             [this](const double time, const Vector<double> &y) {
 *               return this->evaluate_diffusion(time, y);
 *             },
 *             time.get_current_time(),
 *             time.get_next_step_size(),
 *             solution);
 *         time.set_next_step_size(new_time - time.get_current_time());
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 * 
 *         time.set_desired_next_step_size(
 *           embedded_explicit_runge_kutta.get_status().delta_t_guess);
 *       }
 * 
 *     return time.get_step_number();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionruncode"></a> 
 * <h4><code>Diffusion::run</code></h4>
 *   

 * 
 * The following is the main function of the program. At the top, we create
 * the grid (a [0,5]x[0,5] square) and refine it four times to get a mesh
 * that has 16 by 16 cells, for a total of 256.  We then set the boundary
 * indicator to 1 for those parts of the boundary where $x=0$ and $x=5$.
 * 
 * @code
 *   void Diffusion::run()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0., 5.);
 *     triangulation.refine_global(4);
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             if ((face->center()[0] == 0.) || (face->center()[0] == 5.))
 *               face->set_boundary_id(1);
 *             else
 *               face->set_boundary_id(0);
 *           }
 * 
 * @endcode
 * 
 * Next, we set up the linear systems and fill them with content so that
 * they can be used throughout the time stepping process:
 * 
 * @code
 *     setup_system();
 * 
 *     assemble_system();
 * 
 * @endcode
 * 
 * Finally, we solve the diffusion problem using several of the
 * Runge-Kutta methods implemented in namespace TimeStepping, each time
 * outputting the error at the end time. (As explained in the
 * introduction, since the exact solution is zero at the final time, the
 * error equals the numerical solution and can be computed by just taking
 * the $l_2$ norm of the solution vector.)
 * 
 * @code
 *     unsigned int       n_steps      = 0;
 *     const unsigned int n_time_steps = 200;
 *     const double       initial_time = 0.;
 *     const double       final_time   = 10.;
 * 
 *     std::cout << "Explicit methods:" << std::endl;
 *     explicit_method(TimeStepping::FORWARD_EULER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Forward Euler:            error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     explicit_method(TimeStepping::RK_THIRD_ORDER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Third order Runge-Kutta:  error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     explicit_method(TimeStepping::RK_CLASSIC_FOURTH_ORDER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Fourth order Runge-Kutta: error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << std::endl;
 * 
 * 
 *     std::cout << "Implicit methods:" << std::endl;
 *     implicit_method(TimeStepping::BACKWARD_EULER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Backward Euler:           error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::IMPLICIT_MIDPOINT,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Implicit Midpoint:        error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::CRANK_NICOLSON,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Crank-Nicolson:           error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::SDIRK_TWO_STAGES,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   SDIRK:                    error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << std::endl;
 * 
 * 
 *     std::cout << "Embedded explicit methods:" << std::endl;
 *     n_steps = embedded_explicit_method(TimeStepping::HEUN_EULER,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Heun-Euler:               error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::BOGACKI_SHAMPINE,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Bogacki-Shampine:         error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::DOPRI,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Dopri:                    error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::FEHLBERG,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Fehlberg:                 error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::CASH_KARP,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Cash-Karp:                error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 *   }
 * } // namespace Step52
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
 * The following <code>main</code> function is similar to previous examples
 * and need not be commented on.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step52::Diffusion diffusion;
 *       diffusion.run();
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
 *     };
 * 
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


The point of this program is less to show particular results, but instead to
show how it is done. This we have already demonstrated simply by discussing
the code above. Consequently, the output the program yields is relatively
sparse and consists only of the console output and the solutions given in VTU
format for visualization.

The console output contains both errors and, for some of the methods, the
number of steps they performed:
@code
Explicit methods:
   Forward Euler:            error=1.00883
   Third order Runge-Kutta:  error=0.000227982
   Fourth order Runge-Kutta: error=1.90541e-06

Implicit methods:
   Backward Euler:           error=1.03428
   Implicit Midpoint:        error=0.00862702
   Crank-Nicolson:           error=0.00862675
   SDIRK:                    error=0.0042349

Embedded explicit methods:
   Heun-Euler:               error=0.0073012
                   steps performed=284
   Bogacki-Shampine:         error=0.000408407
                   steps performed=181
   Dopri:                    error=0.000836695
                   steps performed=120
   Fehlberg:                 error=0.00248922
                   steps performed=106
   Cash-Karp:                error=0.0787735
                   steps performed=106
@endcode

As expected the higher order methods give (much) more accurate solutions. We
also see that the (rather inaccurate) Heun-Euler method increased the number of
time steps in order to satisfy the tolerance. On the other hand, the other
embedded methods used a lot less time steps than what was prescribed.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-52.cc"
*/
