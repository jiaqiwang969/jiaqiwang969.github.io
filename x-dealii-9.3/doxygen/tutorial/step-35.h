/**
@page step_35 The step-35 tutorial program
This tutorial depends on step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Motivation"> Motivation </a>
        <li><a href="#Projectionmethods"> Projection methods </a>
        <li><a href="#TheFullyDiscreteSetting"> The Fully Discrete Setting </a>
        <li><a href="#Implementation"> Implementation </a>
        <li><a href="#TheTestcase"> The Testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Runtimeparameters">Run time parameters</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ThecodeNavierStokesProjectioncodeclass">The <code>NavierStokesProjection</code> class</a>
      <ul>
        <li><a href="#codeNavierStokesProjectionNavierStokesProjectioncode"> <code>NavierStokesProjection::NavierStokesProjection</code> </a>
        <li><a href="#codeNavierStokesProjectioncreate_triangulation_and_dofscode"><code>NavierStokesProjection::create_triangulation_and_dofs</code></a>
        <li><a href="#codeNavierStokesProjectioninitializecode"> <code>NavierStokesProjection::initialize</code> </a>
        <li><a href="#codeNavierStokesProjectioninitialize__matricescode"> <code>NavierStokesProjection::initialize_*_matrices</code> </a>
        <li><a href="#codeNavierStokesProjectionruncode"> <code>NavierStokesProjection::run</code> </a>
        <li><a href="#codeNavierStokesProjectiondiffusion_stepcode"><code>NavierStokesProjection::diffusion_step</code></a>
        <li><a href="#codeNavierStokesProjectionassemble_advection_termcode"> <code>NavierStokesProjection::assemble_advection_term</code> </a>
        <li><a href="#codeNavierStokesProjectionprojection_stepcode"><code>NavierStokesProjection::projection_step</code></a>
        <li><a href="#codeNavierStokesProjectionupdate_pressurecode"> <code>NavierStokesProjection::update_pressure</code> </a>
        <li><a href="#codeNavierStokesProjectionoutput_resultscode"> <code>NavierStokesProjection::output_results</code> </a>
      </ul>
        <li><a href="#Themainfunction"> The main function </a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Re100"> Re = 100 </a>
        <li><a href="#Re500"> Re = 500 </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program grew out of a student project by Abner Salgado at Texas A&M
University. Most of the work for this program is by him.
</i>

<a name="Intro"></a>
<a name="Introduction"></a><h1> Introduction </h1>


<a name="Motivation"></a>
<a name="Motivation"></a><h3> Motivation </h3>

The purpose of this program is to show how to effectively solve the incompressible time-dependent
Navier-Stokes equations. These equations describe the flow of a viscous incompressible fluid and read
@f{align*}
  u_t + u \cdot \nabla u - \nu \Delta u + \nabla p = f, \\
  \nabla \cdot u = 0,
@f}
where $u$ represents the velocity of the flow and $p$ the pressure. This system of equations is supplemented by
the initial condition
@f[
  u |_{t=0} = u_0,
@f]
with $u_0$ sufficiently smooth and solenoidal, and suitable boundary conditions. For instance, an admissible boundary
condition, is
@f[
  u|_{\partial\Omega} = u_b.
@f]
It is possible to prescribe other boundary conditions as well. In the test case that we solve here the boundary
is partitioned into two disjoint subsets $\partial\Omega = \Gamma_1 \cup \Gamma_2$ and we have
@f[
  u|_{\Gamma_1} = u_b,
@f]
and
@f[
 u\times n|_{\Gamma_2} = 0, \quad p|_{\Gamma_2} = 0
@f]
where $n$ is the outer unit normal. The boundary conditions on $\Gamma_2$ are often
used to model outflow conditions.

In previous tutorial programs (see for instance step-20 and
step-22) we have seen
how to solve the time-independent Stokes equations using a Schur complement approach. For the
time-dependent case, after time discretization, we would arrive at a system like
@f{align*}
  \frac1\tau u^k - \nu \Delta u^k + \nabla p^k = F^k, \\
  \nabla \cdot u^k = 0,
@f}
where $\tau$ is the time-step. Although the structure of this system is similar to the Stokes system and thus
it could be solved using a Schur complement approach, it turns out that the condition number of the
Schur complement is proportional to $\tau^{-2}$. This makes the system very
difficult to solve, and means that for the Navier-Stokes equations, this is
not a useful avenue to the solution.

<a name="Projection"></a>
<a name="Projectionmethods"></a><h3> Projection methods </h3>


Rather, we need to come up with a different approach to solve the time-dependent Navier-Stokes
equations. The difficulty in their solution comes from the fact that the velocity and the pressure are coupled
through the constraint
@f[
  \nabla \cdot u = 0,
@f]
for which the pressure is the Lagrange multiplier.
Projection methods aim at decoupling this constraint from the diffusion (Laplace) operator.

Let us shortly describe how the projection methods look like in a semi-discrete setting. The objective is to
obtain a sequence of velocities $\{u^k\}$ and pressures $\{p^k\}$. We will
also obtain a sequence $\{\phi^k\}$ of auxiliary variables.
Suppose that from the initial conditions, and an application of a first order method we have found
$(u^0,p^0,\phi^0=0)$ and $(u^1,p^1,\phi^1=p^1-p^0)$. Then the projection method consists of the following steps:
<ul>
  <li> <b>Step 0</b>: Extrapolation. Define:
  @f[
    u^\star = 2u^k - u^{k-1}, \quad p^\sharp = p^k + \frac43 \phi^k - \frac13 \phi^{k-1}.
  @f]
  <li> <b>Step 1</b>: Diffusion step. We find $u^{k+1}$ that solves the single
  linear equation
  @f[
    \frac1{2\tau}\left( 3u^{k+1} - 4u^k + u^{k-1} \right)
    + u^\star \cdot\nabla u^{k+1} + \frac12 \left( \nabla \cdot u^\star \right) u^{k+1}
    -\nu \Delta u^{k+1} + \nabla p^\sharp
    = f^{k+1},
    \quad
    u^{k+1}|_{\Gamma_1} = u_b,
    \quad
    u^{k+1} \times n|_{\Gamma_2} = 0.
  @f]

  <li> <b>Step 2</b>: Projection. Find $\phi^{k+1}$ that solves
  @f[
    \Delta \phi^{k+1} = \frac3{2\tau} \nabla \cdot u^{k+1},
    \quad
    \partial_n \phi^{k+1}|_{\Gamma_1} = 0,
    \quad
    \phi^{k+1}|_{\Gamma_2} = 0
  @f]
  <li> <b>Step 3</b>: Pressure correction. Here we have two options:
    <ul>
      <li> <i>Incremental Method in Standard Form</i>. The pressure is updated by:
      @f[
        p^{k+1} = p^k + \phi^{k+1}.
      @f]
      <li> <i>Incremental Method in Rotational Form</i>. In this case
      @f[
        p^{k+1} = p^k + \phi^{k+1} - \nu \nabla \cdot u^{k+1}.
      @f]
    </ul>
</ul>

Without going into details, let us remark a few things about the projection methods that we have just described:
<ul>
  <li> The advection term $u\cdot\nabla u$ is replaced by its <i>skew symmetric form</i>
  @f[
    u \cdot \nabla u + \frac12 \left( \nabla\cdot u \right) u.
  @f]
  This is consistent with the continuous equation (because $\nabla\cdot u = 0$,
  though this is not true pointwise for the discrete solution) and it is needed to
  guarantee unconditional stability of the
  time-stepping scheme. Moreover, to linearize the term we use the second order extrapolation $u^\star$ of
  $u^{k+1}$.
  <li> The projection step is a realization of the Helmholtz decomposition
  @f[
    L^2(\Omega)^d = H \oplus \nabla H^1_{\Gamma_2}(\Omega),
  @f]
  where
  @f[
    H = \left\{ v \in L^2(\Omega)^d:\  \nabla\cdot v =0, \  v\cdot n|_{\Gamma_1} = 0 \right\},
  @f]
  and
  @f[
    H^1_{\Gamma_2}(\Omega) = \left\{ q \in H^1(\Omega):\ q|_{\Gamma_2} = 0 \right\}.
  @f]
  Indeed, if we use this decomposition on $u^{k+1}$ we obtain
  @f[
    u^{k+1} = v^{k+1} + \nabla \left( \frac{2\tau}{3}  \phi^{k+1} \right),
  @f]
  with $v^{k+1}\in H$. Taking the divergence of this equation we arrive at the projection equation.
  <li> The more accurate of the two variants outlined above is the rotational
  one. However, the program below implements both variants. Moreover, in the author's experience,
  the standard form is the one that should be used if, for instance, the viscosity $\nu$ is variable.
</ul>


<p>
The standard incremental scheme and the rotational incremental scheme were first considered by van Kan in
<ul>
  <li> J. van Kan, "A second-order accurate pressure-correction scheme for viscous incompressible flow",
       SIAM Journal on Scientific and Statistical Computing, vol. 7, no. 3, pp. 870–891, 1986
</ul>
and is analyzed by Guermond in
<ul>
  <li> J.-L. Guermond, "Un résultat de convergence d’ordre deux en temps pour
                        l’approximation des équations de Navier–Stokes par une technique de projection incrémentale",
       ESAIM: Mathematical Modelling and Numerical Analysis, vol. 33, no. 1, pp. 169–189, 1999
</ul>
for the case $\nu = 1$.
It turns out that this technique suffers from unphysical boundary conditions for the kinematic pressure that
lead to reduced rates of convergence. To prevent this, Timmermans et al. proposed in
<ul>
  <li> L. Timmermans, P. Minev, and F. Van De Vosse,
       "An approximate projection scheme for incompressible flow using spectral elements",
       International Journal for Numerical Methods in Fluids, vol. 22, no. 7, pp. 673–688, 1996
</ul>
the rotational pressure-correction projection method that uses a divergence correction for the kinematic pressure.
A thorough analysis for scheme has first been performed in
<ul>
  <li> J.-L. Guermond and J. Shen, "On the error estimates for the rotational pressure-correction projection methods",
       Mathematics of Computation, vol. 73, no. 248, pp. 1719–1737, 2004
</ul>
for the Stokes problem.
</p>

<a name ="fullydiscrete"></a>
<a name="TheFullyDiscreteSetting"></a><h3> The Fully Discrete Setting </h3>

To obtain a fully discrete setting of the method we, as always, need a variational formulation. There is one
subtle issue here given the nature of the boundary conditions. When we multiply the equation by a suitable test
function one of the term that arises is
@f[
  -\nu \int_\Omega \Delta u \cdot v.
@f]
If we, say, had Dirichlet boundary conditions on the whole boundary then after integration by parts we would
obtain
@f[
  -\nu \int_\Omega \Delta u \cdot v = \nu \int_\Omega \nabla u : \nabla v
                                    - \int_{\partial\Omega} \partial_n u \cdot v
                                    = \nu \int_\Omega \nabla u : \nabla v.
@f]
One of the advantages of this formulation is that it fully decouples the components of the velocity. Moreover,
they all share the same system matrix. This can be exploited in the program.

However, given the nonstandard boundary conditions, to be able to take them into account we need to use
the following %identity
@f[
  \Delta u = \nabla\nabla\cdot u - \nabla\times\nabla\times u,
@f]
so that when we integrate by parts and take into account the boundary conditions we obtain
@f[
  -\nu \int_\Omega \Delta u \cdot v = \nu \int_\Omega \left[ \nabla \cdot u \nabla \cdot v
                                    + \nabla \times u \nabla \times v \right],
@f]
which is the form that we would have to use. Notice that this couples the components of the velocity.
Moreover, to enforce the boundary condition on the pressure, we need to rewrite
@f[
  \int_\Omega \nabla p \cdot v = -\int_\Omega p \nabla \cdot v + \int_{\Gamma_1} p v\cdot n
                                + \int_{\Gamma_2} p v\cdot n
                               = -\int_\Omega p \nabla \cdot v,
@f]
where the boundary integral in $\Gamma_1$ equals zero given the boundary conditions for the velocity,
and the one in $\Gamma_2$ given the boundary conditions for the pressure.

In the simplified case where the boundary $\Gamma_2$ is %parallel to a coordinate axis, which holds for
the testcase that we carry out below, it can actually be shown that
@f[
  \nu \int_\Omega \nabla u : \nabla v = \nu \int_\Omega \left[ \nabla \cdot u \nabla \cdot v
                                    + \nabla \times u \nabla \times v \right].
@f]
This issue is not very often addressed in the literature. For more information the reader can consult, for
instance,
<ul>
  <li> J.-L. GUERMOND, L. QUARTAPELLE, On the approximation of the unsteady Navier-Stokes equations by
  finite element projection methods, Numer. Math., 80  (1998) 207-238
  <li> J.-L. GUERMOND, P. MINEV, J. SHEN, Error analysis of pressure-correction schemes for the
  Navier-Stokes equations with open boundary conditions, SIAM J. Numer. Anal., 43  1 (2005) 239--258.
</ul>



<a name = "implementation"></a>
<a name="Implementation"></a><h3> Implementation </h3>


Our implementation of the projection methods follows <i>verbatim</i> the description given above. We must note,
however, that as opposed to most other problems that have several solution components, we do not use
vector-valued finite elements. Instead, we use separate finite elements the components of the velocity
and the pressure, respectively, and use different <code>DoFHandler</code>'s for those as well. The main
reason for doing this is that, as we see from the description of the scheme, the <code>dim</code> components
of the velocity and the pressure are decoupled. As a consequence, the equations for all the velocity components
look all the same, have the same system matrix, and can be solved in %parallel. Obviously, this approach
has also its disadvantages. For instance, we need to keep several <code>DoFHandler</code>s and iterators
synchronized when assembling matrices and right hand sides; obtaining quantities that are inherent to
vector-valued functions (e.g. divergences) becomes a little awkward, and others.

<a name ="testcase"></a>
<a name="TheTestcase"></a><h3> The Testcase </h3>


The testcase that we use for this program consists of the flow around a square obstacle. The geometry is
as follows:

<img src="https://www.dealii.org/images/steps/developer/step-35.geometry.png" alt="">

with $H=4.1$, making the geometry slightly non-symmetric.

We impose no-slip boundary conditions on both the top and bottom walls and the obstacle. On the left side we
have the inflow boundary condition
@f[
  u =
  \left( \begin{array}{c} 4 U_m y (H-y)/H^2 \\ 0 \end{array} \right),
@f]
with $U_m = 1.5$, i.e. the inflow boundary conditions correspond to Poiseuille flow for this configuration.
Finally, on the right vertical wall we impose the condition that the vertical component of the velocity
and the pressure should both be zero.
The final time $T=10$.
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
 * We start by including all the necessary deal.II header files and some C++
 * related ones. Each one of them has been discussed in previous tutorial
 * programs, so we will not get into details here.
 * 
 * @code
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/base/point.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/multithread_info.h>
 * #include <deal.II/base/thread_management.h>
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/base/parallel.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/conditional_ostream.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/sparse_ilu.h>
 * #include <deal.II/lac/sparse_direct.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/grid_in.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_tools.h>
 * #include <deal.II/fe/fe_system.h>
 * 
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <cmath>
 * #include <iostream>
 * 
 * @endcode
 * 
 * Finally this is as in all previous programs:
 * 
 * @code
 * namespace Step35
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameters"></a> 
 * <h3>Run time parameters</h3>
 *   

 * 
 * Since our method has several parameters that can be fine-tuned we put them
 * into an external file, so that they can be determined at run-time.
 *   

 * 
 * This includes, in particular, the formulation of the equation for the
 * auxiliary variable $\phi$, for which we declare an <code>enum</code>. Next,
 * we declare a class that is going to read and store all the parameters that
 * our program needs to run.
 * 
 * @code
 *   namespace RunTimeParameters
 *   {
 *     enum class Method
 *     {
 *       standard,
 *       rotational
 *     };
 * 
 *     class Data_Storage
 *     {
 *     public:
 *       Data_Storage();
 * 
 *       void read_data(const std::string &filename);
 * 
 *       Method form;
 * 
 *       double dt;
 *       double initial_time;
 *       double final_time;
 * 
 *       double Reynolds;
 * 
 *       unsigned int n_global_refines;
 * 
 *       unsigned int pressure_degree;
 * 
 *       unsigned int vel_max_iterations;
 *       unsigned int vel_Krylov_size;
 *       unsigned int vel_off_diagonals;
 *       unsigned int vel_update_prec;
 *       double       vel_eps;
 *       double       vel_diag_strength;
 * 
 *       bool         verbose;
 *       unsigned int output_interval;
 * 
 *     protected:
 *       ParameterHandler prm;
 *     };
 * 
 * @endcode
 * 
 * In the constructor of this class we declare all the parameters. The
 * details of how this works have been discussed elsewhere, for example in
 * step-29.
 * 
 * @code
 *     Data_Storage::Data_Storage()
 *       : form(Method::rotational)
 *       , dt(5e-4)
 *       , initial_time(0.)
 *       , final_time(1.)
 *       , Reynolds(1.)
 *       , n_global_refines(0)
 *       , pressure_degree(1)
 *       , vel_max_iterations(1000)
 *       , vel_Krylov_size(30)
 *       , vel_off_diagonals(60)
 *       , vel_update_prec(15)
 *       , vel_eps(1e-12)
 *       , vel_diag_strength(0.01)
 *       , verbose(true)
 *       , output_interval(15)
 *     {
 *       prm.declare_entry("Method_Form",
 *                         "rotational",
 *                         Patterns::Selection("rotational|standard"),
 *                         " Used to select the type of method that we are going "
 *                         "to use. ");
 *       prm.enter_subsection("Physical data");
 *       {
 *         prm.declare_entry("initial_time",
 *                           "0.",
 *                           Patterns::Double(0.),
 *                           " The initial time of the simulation. ");
 *         prm.declare_entry("final_time",
 *                           "1.",
 *                           Patterns::Double(0.),
 *                           " The final time of the simulation. ");
 *         prm.declare_entry("Reynolds",
 *                           "1.",
 *                           Patterns::Double(0.),
 *                           " The Reynolds number. ");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Time step data");
 *       {
 *         prm.declare_entry("dt",
 *                           "5e-4",
 *                           Patterns::Double(0.),
 *                           " The time step size. ");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Space discretization");
 *       {
 *         prm.declare_entry("n_of_refines",
 *                           "0",
 *                           Patterns::Integer(0, 15),
 *                           " The number of global refines we do on the mesh. ");
 *         prm.declare_entry("pressure_fe_degree",
 *                           "1",
 *                           Patterns::Integer(1, 5),
 *                           " The polynomial degree for the pressure space. ");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Data solve velocity");
 *       {
 *         prm.declare_entry(
 *           "max_iterations",
 *           "1000",
 *           Patterns::Integer(1, 1000),
 *           " The maximal number of iterations GMRES must make. ");
 *         prm.declare_entry("eps",
 *                           "1e-12",
 *                           Patterns::Double(0.),
 *                           " The stopping criterion. ");
 *         prm.declare_entry("Krylov_size",
 *                           "30",
 *                           Patterns::Integer(1),
 *                           " The size of the Krylov subspace to be used. ");
 *         prm.declare_entry("off_diagonals",
 *                           "60",
 *                           Patterns::Integer(0),
 *                           " The number of off-diagonal elements ILU must "
 *                           "compute. ");
 *         prm.declare_entry("diag_strength",
 *                           "0.01",
 *                           Patterns::Double(0.),
 *                           " Diagonal strengthening coefficient. ");
 *         prm.declare_entry("update_prec",
 *                           "15",
 *                           Patterns::Integer(1),
 *                           " This number indicates how often we need to "
 *                           "update the preconditioner");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.declare_entry("verbose",
 *                         "true",
 *                         Patterns::Bool(),
 *                         " This indicates whether the output of the solution "
 *                         "process should be verbose. ");
 * 
 *       prm.declare_entry("output_interval",
 *                         "1",
 *                         Patterns::Integer(1),
 *                         " This indicates between how many time steps we print "
 *                         "the solution. ");
 *     }
 * 
 * 
 * 
 *     void Data_Storage::read_data(const std::string &filename)
 *     {
 *       std::ifstream file(filename);
 *       AssertThrow(file, ExcFileNotOpen(filename));
 * 
 *       prm.parse_input(file);
 * 
 *       if (prm.get("Method_Form") == std::string("rotational"))
 *         form = Method::rotational;
 *       else
 *         form = Method::standard;
 * 
 *       prm.enter_subsection("Physical data");
 *       {
 *         initial_time = prm.get_double("initial_time");
 *         final_time   = prm.get_double("final_time");
 *         Reynolds     = prm.get_double("Reynolds");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Time step data");
 *       {
 *         dt = prm.get_double("dt");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Space discretization");
 *       {
 *         n_global_refines = prm.get_integer("n_of_refines");
 *         pressure_degree  = prm.get_integer("pressure_fe_degree");
 *       }
 *       prm.leave_subsection();
 * 
 *       prm.enter_subsection("Data solve velocity");
 *       {
 *         vel_max_iterations = prm.get_integer("max_iterations");
 *         vel_eps            = prm.get_double("eps");
 *         vel_Krylov_size    = prm.get_integer("Krylov_size");
 *         vel_off_diagonals  = prm.get_integer("off_diagonals");
 *         vel_diag_strength  = prm.get_double("diag_strength");
 *         vel_update_prec    = prm.get_integer("update_prec");
 *       }
 *       prm.leave_subsection();
 * 
 *       verbose = prm.get_bool("verbose");
 * 
 *       output_interval = prm.get_integer("output_interval");
 *     }
 *   } // namespace RunTimeParameters
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
 * In the next namespace, we declare the initial and boundary conditions:
 * 
 * @code
 *   namespace EquationData
 *   {
 * @endcode
 * 
 * As we have chosen a completely decoupled formulation, we will not take
 * advantage of deal.II's capabilities to handle vector valued problems. We
 * do, however, want to use an interface for the equation data that is
 * somehow dimension independent. To be able to do that, our functions
 * should be able to know on which spatial component we are currently
 * working, and we should be able to have a common interface to do that. The
 * following class is an attempt in that direction.
 * 
 * @code
 *     template <int dim>
 *     class MultiComponentFunction : public Function<dim>
 *     {
 *     public:
 *       MultiComponentFunction(const double initial_time = 0.);
 *       void set_component(const unsigned int d);
 * 
 *     protected:
 *       unsigned int comp;
 *     };
 * 
 *     template <int dim>
 *     MultiComponentFunction<dim>::MultiComponentFunction(
 *       const double initial_time)
 *       : Function<dim>(1, initial_time)
 *       , comp(0)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void MultiComponentFunction<dim>::set_component(const unsigned int d)
 *     {
 *       Assert(d < dim, ExcIndexRange(d, 0, dim));
 *       comp = d;
 *     }
 * 
 * 
 * @endcode
 * 
 * With this class defined, we declare classes that describe the boundary
 * conditions for velocity and pressure:
 * 
 * @code
 *     template <int dim>
 *     class Velocity : public MultiComponentFunction<dim>
 *     {
 *     public:
 *       Velocity(const double initial_time = 0.0);
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void value_list(const std::vector<Point<dim>> &points,
 *                               std::vector<double> &          values,
 *                               const unsigned int component = 0) const override;
 *     };
 * 
 * 
 *     template <int dim>
 *     Velocity<dim>::Velocity(const double initial_time)
 *       : MultiComponentFunction<dim>(initial_time)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void Velocity<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                    std::vector<double> &          values,
 *                                    const unsigned int) const
 *     {
 *       const unsigned int n_points = points.size();
 *       Assert(values.size() == n_points,
 *              ExcDimensionMismatch(values.size(), n_points));
 *       for (unsigned int i = 0; i < n_points; ++i)
 *         values[i] = Velocity<dim>::value(points[i]);
 *     }
 * 
 * 
 *     template <int dim>
 *     double Velocity<dim>::value(const Point<dim> &p, const unsigned int) const
 *     {
 *       if (this->comp == 0)
 *         {
 *           const double Um = 1.5;
 *           const double H  = 4.1;
 *           return 4. * Um * p(1) * (H - p(1)) / (H * H);
 *         }
 *       else
 *         return 0.;
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     class Pressure : public Function<dim>
 *     {
 *     public:
 *       Pressure(const double initial_time = 0.0);
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void value_list(const std::vector<Point<dim>> &points,
 *                               std::vector<double> &          values,
 *                               const unsigned int component = 0) const override;
 *     };
 * 
 *     template <int dim>
 *     Pressure<dim>::Pressure(const double initial_time)
 *       : Function<dim>(1, initial_time)
 *     {}
 * 
 * 
 *     template <int dim>
 *     double Pressure<dim>::value(const Point<dim> & p,
 *                                 const unsigned int component) const
 *     {
 *       (void)component;
 *       AssertIndexRange(component, 1);
 *       return 25. - p(0);
 *     }
 * 
 *     template <int dim>
 *     void Pressure<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                    std::vector<double> &          values,
 *                                    const unsigned int component) const
 *     {
 *       (void)component;
 *       AssertIndexRange(component, 1);
 *       const unsigned int n_points = points.size();
 *       Assert(values.size() == n_points,
 *              ExcDimensionMismatch(values.size(), n_points));
 *       for (unsigned int i = 0; i < n_points; ++i)
 *         values[i] = Pressure<dim>::value(points[i]);
 *     }
 *   } // namespace EquationData
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeNavierStokesProjectioncodeclass"></a> 
 * <h3>The <code>NavierStokesProjection</code> class</h3>
 * 

 * 
 * Now for the main class of the program. It implements the various versions
 * of the projection method for Navier-Stokes equations. The names for all the
 * methods and member variables should be self-explanatory, taking into
 * account the implementation details given in the introduction.
 * 
 * @code
 *   template <int dim>
 *   class NavierStokesProjection
 *   {
 *   public:
 *     NavierStokesProjection(const RunTimeParameters::Data_Storage &data);
 * 
 *     void run(const bool verbose = false, const unsigned int n_plots = 10);
 * 
 *   protected:
 *     RunTimeParameters::Method type;
 * 
 *     const unsigned int deg;
 *     const double       dt;
 *     const double       t_0;
 *     const double       T;
 *     const double       Re;
 * 
 *     EquationData::Velocity<dim>               vel_exact;
 *     std::map<types::global_dof_index, double> boundary_values;
 *     std::vector<types::boundary_id>           boundary_ids;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     FE_Q<dim> fe_velocity;
 *     FE_Q<dim> fe_pressure;
 * 
 *     DoFHandler<dim> dof_handler_velocity;
 *     DoFHandler<dim> dof_handler_pressure;
 * 
 *     QGauss<dim> quadrature_pressure;
 *     QGauss<dim> quadrature_velocity;
 * 
 *     SparsityPattern sparsity_pattern_velocity;
 *     SparsityPattern sparsity_pattern_pressure;
 *     SparsityPattern sparsity_pattern_pres_vel;
 * 
 *     SparseMatrix<double> vel_Laplace_plus_Mass;
 *     SparseMatrix<double> vel_it_matrix[dim];
 *     SparseMatrix<double> vel_Mass;
 *     SparseMatrix<double> vel_Laplace;
 *     SparseMatrix<double> vel_Advection;
 *     SparseMatrix<double> pres_Laplace;
 *     SparseMatrix<double> pres_Mass;
 *     SparseMatrix<double> pres_Diff[dim];
 *     SparseMatrix<double> pres_iterative;
 * 
 *     Vector<double> pres_n;
 *     Vector<double> pres_n_minus_1;
 *     Vector<double> phi_n;
 *     Vector<double> phi_n_minus_1;
 *     Vector<double> u_n[dim];
 *     Vector<double> u_n_minus_1[dim];
 *     Vector<double> u_star[dim];
 *     Vector<double> force[dim];
 *     Vector<double> v_tmp;
 *     Vector<double> pres_tmp;
 *     Vector<double> rot_u;
 * 
 *     SparseILU<double>   prec_velocity[dim];
 *     SparseILU<double>   prec_pres_Laplace;
 *     SparseDirectUMFPACK prec_mass;
 *     SparseDirectUMFPACK prec_vel_mass;
 * 
 *     DeclException2(ExcInvalidTimeStep,
 *                    double,
 *                    double,
 *                    << " The time step " << arg1 << " is out of range."
 *                    << std::endl
 *                    << " The permitted range is (0," << arg2 << "]");
 * 
 *     void create_triangulation_and_dofs(const unsigned int n_refines);
 * 
 *     void initialize();
 * 
 *     void interpolate_velocity();
 * 
 *     void diffusion_step(const bool reinit_prec);
 * 
 *     void projection_step(const bool reinit_prec);
 * 
 *     void update_pressure(const bool reinit_prec);
 * 
 *   private:
 *     unsigned int vel_max_its;
 *     unsigned int vel_Krylov_size;
 *     unsigned int vel_off_diagonals;
 *     unsigned int vel_update_prec;
 *     double       vel_eps;
 *     double       vel_diag_strength;
 * 
 *     void initialize_velocity_matrices();
 * 
 *     void initialize_pressure_matrices();
 * 
 * @endcode
 * 
 * The next few structures and functions are for doing various things in
 * parallel. They follow the scheme laid out in @ref threads, using the
 * WorkStream class. As explained there, this requires us to declare two
 * structures for each of the assemblers, a per-task data and a scratch data
 * structure. These are then handed over to functions that assemble local
 * contributions and that copy these local contributions to the global
 * objects.
 *     

 * 
 * One of the things that are specific to this program is that we don't just
 * have a single DoFHandler object that represents both the velocities and
 * the pressure, but we use individual DoFHandler objects for these two
 * kinds of variables. We pay for this optimization when we want to assemble
 * terms that involve both variables, such as the divergence of the velocity
 * and the gradient of the pressure, times the respective test functions.
 * When doing so, we can't just anymore use a single FEValues object, but
 * rather we need two, and they need to be initialized with cell iterators
 * that point to the same cell in the triangulation but different
 * DoFHandlers.
 *     

 * 
 * To do this in practice, we declare a "synchronous" iterator -- an object
 * that internally consists of several (in our case two) iterators, and each
 * time the synchronous iteration is moved forward one step, each of the
 * iterators stored internally is moved forward one step as well, thereby
 * always staying in sync. As it so happens, there is a deal.II class that
 * facilitates this sort of thing. (What is important here is to know that
 * two DoFHandler objects built on the same triangulation will walk over the
 * cells of the triangulation in the same order.)
 * 
 * @code
 *     using IteratorTuple =
 *       std::tuple<typename DoFHandler<dim>::active_cell_iterator,
 *                  typename DoFHandler<dim>::active_cell_iterator>;
 * 
 *     using IteratorPair = SynchronousIterators<IteratorTuple>;
 * 
 *     void initialize_gradient_operator();
 * 
 *     struct InitGradPerTaskData
 *     {
 *       unsigned int                         d;
 *       unsigned int                         vel_dpc;
 *       unsigned int                         pres_dpc;
 *       FullMatrix<double>                   local_grad;
 *       std::vector<types::global_dof_index> vel_local_dof_indices;
 *       std::vector<types::global_dof_index> pres_local_dof_indices;
 * 
 *       InitGradPerTaskData(const unsigned int dd,
 *                           const unsigned int vdpc,
 *                           const unsigned int pdpc)
 *         : d(dd)
 *         , vel_dpc(vdpc)
 *         , pres_dpc(pdpc)
 *         , local_grad(vdpc, pdpc)
 *         , vel_local_dof_indices(vdpc)
 *         , pres_local_dof_indices(pdpc)
 *       {}
 *     };
 * 
 *     struct InitGradScratchData
 *     {
 *       unsigned int  nqp;
 *       FEValues<dim> fe_val_vel;
 *       FEValues<dim> fe_val_pres;
 *       InitGradScratchData(const FE_Q<dim> &  fe_v,
 *                           const FE_Q<dim> &  fe_p,
 *                           const QGauss<dim> &quad,
 *                           const UpdateFlags  flags_v,
 *                           const UpdateFlags  flags_p)
 *         : nqp(quad.size())
 *         , fe_val_vel(fe_v, quad, flags_v)
 *         , fe_val_pres(fe_p, quad, flags_p)
 *       {}
 *       InitGradScratchData(const InitGradScratchData &data)
 *         : nqp(data.nqp)
 *         , fe_val_vel(data.fe_val_vel.get_fe(),
 *                      data.fe_val_vel.get_quadrature(),
 *                      data.fe_val_vel.get_update_flags())
 *         , fe_val_pres(data.fe_val_pres.get_fe(),
 *                       data.fe_val_pres.get_quadrature(),
 *                       data.fe_val_pres.get_update_flags())
 *       {}
 *     };
 * 
 *     void assemble_one_cell_of_gradient(const IteratorPair & SI,
 *                                        InitGradScratchData &scratch,
 *                                        InitGradPerTaskData &data);
 * 
 *     void copy_gradient_local_to_global(const InitGradPerTaskData &data);
 * 
 * @endcode
 * 
 * The same general layout also applies to the following classes and
 * functions implementing the assembly of the advection term:
 * 
 * @code
 *     void assemble_advection_term();
 * 
 *     struct AdvectionPerTaskData
 *     {
 *       FullMatrix<double>                   local_advection;
 *       std::vector<types::global_dof_index> local_dof_indices;
 *       AdvectionPerTaskData(const unsigned int dpc)
 *         : local_advection(dpc, dpc)
 *         , local_dof_indices(dpc)
 *       {}
 *     };
 * 
 *     struct AdvectionScratchData
 *     {
 *       unsigned int                nqp;
 *       unsigned int                dpc;
 *       std::vector<Point<dim>>     u_star_local;
 *       std::vector<Tensor<1, dim>> grad_u_star;
 *       std::vector<double>         u_star_tmp;
 *       FEValues<dim>               fe_val;
 *       AdvectionScratchData(const FE_Q<dim> &  fe,
 *                            const QGauss<dim> &quad,
 *                            const UpdateFlags  flags)
 *         : nqp(quad.size())
 *         , dpc(fe.n_dofs_per_cell())
 *         , u_star_local(nqp)
 *         , grad_u_star(nqp)
 *         , u_star_tmp(nqp)
 *         , fe_val(fe, quad, flags)
 *       {}
 * 
 *       AdvectionScratchData(const AdvectionScratchData &data)
 *         : nqp(data.nqp)
 *         , dpc(data.dpc)
 *         , u_star_local(nqp)
 *         , grad_u_star(nqp)
 *         , u_star_tmp(nqp)
 *         , fe_val(data.fe_val.get_fe(),
 *                  data.fe_val.get_quadrature(),
 *                  data.fe_val.get_update_flags())
 *       {}
 *     };
 * 
 *     void assemble_one_cell_of_advection(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       AdvectionScratchData &                                scratch,
 *       AdvectionPerTaskData &                                data);
 * 
 *     void copy_advection_local_to_global(const AdvectionPerTaskData &data);
 * 
 * @endcode
 * 
 * The final few functions implement the diffusion solve as well as
 * postprocessing the output, including computing the curl of the velocity:
 * 
 * @code
 *     void diffusion_component_solve(const unsigned int d);
 * 
 *     void output_results(const unsigned int step);
 * 
 *     void assemble_vorticity(const bool reinit_prec);
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionNavierStokesProjectioncode"></a> 
 * <h4> <code>NavierStokesProjection::NavierStokesProjection</code> </h4>
 * 

 * 
 * In the constructor, we just read all the data from the
 * <code>Data_Storage</code> object that is passed as an argument, verify that
 * the data we read is reasonable and, finally, create the triangulation and
 * load the initial data.
 * 
 * @code
 *   template <int dim>
 *   NavierStokesProjection<dim>::NavierStokesProjection(
 *     const RunTimeParameters::Data_Storage &data)
 *     : type(data.form)
 *     , deg(data.pressure_degree)
 *     , dt(data.dt)
 *     , t_0(data.initial_time)
 *     , T(data.final_time)
 *     , Re(data.Reynolds)
 *     , vel_exact(data.initial_time)
 *     , fe_velocity(deg + 1)
 *     , fe_pressure(deg)
 *     , dof_handler_velocity(triangulation)
 *     , dof_handler_pressure(triangulation)
 *     , quadrature_pressure(deg + 1)
 *     , quadrature_velocity(deg + 2)
 *     , vel_max_its(data.vel_max_iterations)
 *     , vel_Krylov_size(data.vel_Krylov_size)
 *     , vel_off_diagonals(data.vel_off_diagonals)
 *     , vel_update_prec(data.vel_update_prec)
 *     , vel_eps(data.vel_eps)
 *     , vel_diag_strength(data.vel_diag_strength)
 *   {
 *     if (deg < 1)
 *       std::cout
 *         << " WARNING: The chosen pair of finite element spaces is not stable."
 *         << std::endl
 *         << " The obtained results will be nonsense" << std::endl;
 * 
 *     AssertThrow(!((dt <= 0.) || (dt > .5 * T)), ExcInvalidTimeStep(dt, .5 * T));
 * 
 *     create_triangulation_and_dofs(data.n_global_refines);
 *     initialize();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectioncreate_triangulation_and_dofscode"></a> 
 * <h4><code>NavierStokesProjection::create_triangulation_and_dofs</code></h4>
 * 

 * 
 * The method that creates the triangulation and refines it the needed number
 * of times. After creating the triangulation, it creates the mesh dependent
 * data, i.e. it distributes degrees of freedom and renumbers them, and
 * initializes the matrices and vectors that we will use.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::create_triangulation_and_dofs(
 *     const unsigned int n_refines)
 *   {
 *     GridIn<dim> grid_in;
 *     grid_in.attach_triangulation(triangulation);
 * 
 *     {
 *       std::string   filename = "nsbench2.inp";
 *       std::ifstream file(filename);
 *       Assert(file, ExcFileNotOpen(filename.c_str()));
 *       grid_in.read_ucd(file);
 *     }
 * 
 *     std::cout << "Number of refines = " << n_refines << std::endl;
 *     triangulation.refine_global(n_refines);
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl;
 * 
 *     boundary_ids = triangulation.get_boundary_ids();
 * 
 *     dof_handler_velocity.distribute_dofs(fe_velocity);
 *     DoFRenumbering::boost::Cuthill_McKee(dof_handler_velocity);
 *     dof_handler_pressure.distribute_dofs(fe_pressure);
 *     DoFRenumbering::boost::Cuthill_McKee(dof_handler_pressure);
 * 
 *     initialize_velocity_matrices();
 *     initialize_pressure_matrices();
 *     initialize_gradient_operator();
 * 
 *     pres_n.reinit(dof_handler_pressure.n_dofs());
 *     pres_n_minus_1.reinit(dof_handler_pressure.n_dofs());
 *     phi_n.reinit(dof_handler_pressure.n_dofs());
 *     phi_n_minus_1.reinit(dof_handler_pressure.n_dofs());
 *     pres_tmp.reinit(dof_handler_pressure.n_dofs());
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         u_n[d].reinit(dof_handler_velocity.n_dofs());
 *         u_n_minus_1[d].reinit(dof_handler_velocity.n_dofs());
 *         u_star[d].reinit(dof_handler_velocity.n_dofs());
 *         force[d].reinit(dof_handler_velocity.n_dofs());
 *       }
 *     v_tmp.reinit(dof_handler_velocity.n_dofs());
 *     rot_u.reinit(dof_handler_velocity.n_dofs());
 * 
 *     std::cout << "dim (X_h) = " << (dof_handler_velocity.n_dofs() * dim) 
 *               << std::endl                                               
 *               << "dim (M_h) = " << dof_handler_pressure.n_dofs()         
 *               << std::endl                                               
 *               << "Re        = " << Re << std::endl                       
 *               << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectioninitializecode"></a> 
 * <h4> <code>NavierStokesProjection::initialize</code> </h4>
 * 

 * 
 * This method creates the constant matrices and loads the initial data
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::initialize()
 *   {
 *     vel_Laplace_plus_Mass = 0.;
 *     vel_Laplace_plus_Mass.add(1. / Re, vel_Laplace);
 *     vel_Laplace_plus_Mass.add(1.5 / dt, vel_Mass);
 * 
 *     EquationData::Pressure<dim> pres(t_0);
 *     VectorTools::interpolate(dof_handler_pressure, pres, pres_n_minus_1);
 *     pres.advance_time(dt);
 *     VectorTools::interpolate(dof_handler_pressure, pres, pres_n);
 *     phi_n         = 0.;
 *     phi_n_minus_1 = 0.;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         vel_exact.set_time(t_0);
 *         vel_exact.set_component(d);
 *         VectorTools::interpolate(dof_handler_velocity,
 *                                  vel_exact,
 *                                  u_n_minus_1[d]);
 *         vel_exact.advance_time(dt);
 *         VectorTools::interpolate(dof_handler_velocity, vel_exact, u_n[d]);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectioninitialize__matricescode"></a> 
 * <h4> <code>NavierStokesProjection::initialize_*_matrices</code> </h4>
 * 

 * 
 * In this set of methods we initialize the sparsity patterns, the constraints
 * (if any) and assemble the matrices that do not depend on the timestep
 * <code>dt</code>. Note that for the Laplace and mass matrices, we can use
 * functions in the library that do this. Because the expensive operations of
 * this function -- creating the two matrices -- are entirely independent, we
 * could in principle mark them as tasks that can be worked on in %parallel
 * using the Threads::new_task functions. We won't do that here since these
 * functions internally already are parallelized, and in particular because
 * the current function is only called once per program run and so does not
 * incur a cost in each time step. The necessary modifications would be quite
 * straightforward, however.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::initialize_velocity_matrices()
 *   {
 *     {
 *       DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(),
 *                                  dof_handler_velocity.n_dofs());
 *       DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp);
 *       sparsity_pattern_velocity.copy_from(dsp);
 *     }
 *     vel_Laplace_plus_Mass.reinit(sparsity_pattern_velocity);
 *     for (unsigned int d = 0; d < dim; ++d)
 *       vel_it_matrix[d].reinit(sparsity_pattern_velocity);
 *     vel_Mass.reinit(sparsity_pattern_velocity);
 *     vel_Laplace.reinit(sparsity_pattern_velocity);
 *     vel_Advection.reinit(sparsity_pattern_velocity);
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler_velocity,
 *                                       quadrature_velocity,
 *                                       vel_Mass);
 *     MatrixCreator::create_laplace_matrix(dof_handler_velocity,
 *                                          quadrature_velocity,
 *                                          vel_Laplace);
 *   }
 * 
 * @endcode
 * 
 * The initialization of the matrices that act on the pressure space is
 * similar to the ones that act on the velocity space.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::initialize_pressure_matrices()
 *   {
 *     {
 *       DynamicSparsityPattern dsp(dof_handler_pressure.n_dofs(),
 *                                  dof_handler_pressure.n_dofs());
 *       DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp);
 *       sparsity_pattern_pressure.copy_from(dsp);
 *     }
 * 
 *     pres_Laplace.reinit(sparsity_pattern_pressure);
 *     pres_iterative.reinit(sparsity_pattern_pressure);
 *     pres_Mass.reinit(sparsity_pattern_pressure);
 * 
 *     MatrixCreator::create_laplace_matrix(dof_handler_pressure,
 *                                          quadrature_pressure,
 *                                          pres_Laplace);
 *     MatrixCreator::create_mass_matrix(dof_handler_pressure,
 *                                       quadrature_pressure,
 *                                       pres_Mass);
 *   }
 * 
 * 
 * @endcode
 * 
 * For the gradient operator, we start by initializing the sparsity pattern
 * and compressing it. It is important to notice here that the gradient
 * operator acts from the pressure space into the velocity space, so we have
 * to deal with two different finite element spaces. To keep the loops
 * synchronized, we use the alias that we have defined before, namely
 * <code>PairedIterators</code> and <code>IteratorPair</code>.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::initialize_gradient_operator()
 *   {
 *     {
 *       DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(),
 *                                  dof_handler_pressure.n_dofs());
 *       DoFTools::make_sparsity_pattern(dof_handler_velocity,
 *                                       dof_handler_pressure,
 *                                       dsp);
 *       sparsity_pattern_pres_vel.copy_from(dsp);
 *     }
 * 
 *     InitGradPerTaskData per_task_data(0,
 *                                       fe_velocity.n_dofs_per_cell(),
 *                                       fe_pressure.n_dofs_per_cell());
 *     InitGradScratchData scratch_data(fe_velocity,
 *                                      fe_pressure,
 *                                      quadrature_velocity,
 *                                      update_gradients | update_JxW_values,
 *                                      update_values);
 * 
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         pres_Diff[d].reinit(sparsity_pattern_pres_vel);
 *         per_task_data.d = d;
 *         WorkStream::run(
 *           IteratorPair(IteratorTuple(dof_handler_velocity.begin_active(),
 *                                      dof_handler_pressure.begin_active())),
 *           IteratorPair(IteratorTuple(dof_handler_velocity.end(),
 *                                      dof_handler_pressure.end())),
 *           *this,
 *           &NavierStokesProjection<dim>::assemble_one_cell_of_gradient,
 *           &NavierStokesProjection<dim>::copy_gradient_local_to_global,
 *           scratch_data,
 *           per_task_data);
 *       }
 *   }
 * 
 *   template <int dim>
 *   void NavierStokesProjection<dim>::assemble_one_cell_of_gradient(
 *     const IteratorPair & SI,
 *     InitGradScratchData &scratch,
 *     InitGradPerTaskData &data)
 *   {
 *     scratch.fe_val_vel.reinit(std::get<0>(*SI));
 *     scratch.fe_val_pres.reinit(std::get<1>(*SI));
 * 
 *     std::get<0>(*SI)->get_dof_indices(data.vel_local_dof_indices);
 *     std::get<1>(*SI)->get_dof_indices(data.pres_local_dof_indices);
 * 
 *     data.local_grad = 0.;
 *     for (unsigned int q = 0; q < scratch.nqp; ++q)
 *       {
 *         for (unsigned int i = 0; i < data.vel_dpc; ++i)
 *           for (unsigned int j = 0; j < data.pres_dpc; ++j)
 *             data.local_grad(i, j) +=
 *               -scratch.fe_val_vel.JxW(q) *
 *               scratch.fe_val_vel.shape_grad(i, q)[data.d] *
 *               scratch.fe_val_pres.shape_value(j, q);
 *       }
 *   }
 * 
 * 
 *   template <int dim>
 *   void NavierStokesProjection<dim>::copy_gradient_local_to_global(
 *     const InitGradPerTaskData &data)
 *   {
 *     for (unsigned int i = 0; i < data.vel_dpc; ++i)
 *       for (unsigned int j = 0; j < data.pres_dpc; ++j)
 *         pres_Diff[data.d].add(data.vel_local_dof_indices[i],
 *                               data.pres_local_dof_indices[j],
 *                               data.local_grad(i, j));
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionruncode"></a> 
 * <h4> <code>NavierStokesProjection::run</code> </h4>
 * 

 * 
 * This is the time marching function, which starting at <code>t_0</code>
 * advances in time using the projection method with time step <code>dt</code>
 * until <code>T</code>.
 *   

 * 
 * Its second parameter, <code>verbose</code> indicates whether the function
 * should output information what it is doing at any given moment: for
 * example, it will say whether we are working on the diffusion, projection
 * substep; updating preconditioners etc. Rather than implementing this
 * output using code like
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   if (verbose) std::cout << "something";
 * @endcode
 * </div>
 * we use the ConditionalOStream class to do that for us. That
 * class takes an output stream and a condition that indicates whether the
 * things you pass to it should be passed through to the given output
 * stream, or should just be ignored. This way, above code simply becomes
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   verbose_cout << "something";
 * @endcode
 * </div>
 * and does the right thing in either case.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::run(const bool         verbose,
 *                                         const unsigned int output_interval)
 *   {
 *     ConditionalOStream verbose_cout(std::cout, verbose);
 * 
 *     const auto n_steps = static_cast<unsigned int>((T - t_0) / dt);
 *     vel_exact.set_time(2. * dt);
 *     output_results(1);
 *     for (unsigned int n = 2; n <= n_steps; ++n)
 *       {
 *         if (n % output_interval == 0)
 *           {
 *             verbose_cout << "Plotting Solution" << std::endl;
 *             output_results(n);
 *           }
 *         std::cout << "Step = " << n << " Time = " << (n * dt) << std::endl;
 *         verbose_cout << "  Interpolating the velocity " << std::endl;
 * 
 *         interpolate_velocity();
 *         verbose_cout << "  Diffusion Step" << std::endl;
 *         if (n % vel_update_prec == 0)
 *           verbose_cout << "    With reinitialization of the preconditioner"
 *                        << std::endl;
 *         diffusion_step((n % vel_update_prec == 0) || (n == 2));
 *         verbose_cout << "  Projection Step" << std::endl;
 *         projection_step((n == 2));
 *         verbose_cout << "  Updating the Pressure" << std::endl;
 *         update_pressure((n == 2));
 *         vel_exact.advance_time(dt);
 *       }
 *     output_results(n_steps);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void NavierStokesProjection<dim>::interpolate_velocity()
 *   {
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         u_star[d].equ(2., u_n[d]);
 *         u_star[d] -= u_n_minus_1[d];
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectiondiffusion_stepcode"></a> 
 * <h4><code>NavierStokesProjection::diffusion_step</code></h4>
 * 

 * 
 * The implementation of a diffusion step. Note that the expensive operation
 * is the diffusion solve at the end of the function, which we have to do once
 * for each velocity component. To accelerate things a bit, we allow to do
 * this in %parallel, using the Threads::new_task function which makes sure
 * that the <code>dim</code> solves are all taken care of and are scheduled to
 * available processors: if your machine has more than one processor core and
 * no other parts of this program are using resources currently, then the
 * diffusion solves will run in %parallel. On the other hand, if your system
 * has only one processor core then running things in %parallel would be
 * inefficient (since it leads, for example, to cache congestion) and things
 * will be executed sequentially.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::diffusion_step(const bool reinit_prec)
 *   {
 *     pres_tmp.equ(-1., pres_n);
 *     pres_tmp.add(-4. / 3., phi_n, 1. / 3., phi_n_minus_1);
 * 
 *     assemble_advection_term();
 * 
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         force[d] = 0.;
 *         v_tmp.equ(2. / dt, u_n[d]);
 *         v_tmp.add(-.5 / dt, u_n_minus_1[d]);
 *         vel_Mass.vmult_add(force[d], v_tmp);
 * 
 *         pres_Diff[d].vmult_add(force[d], pres_tmp);
 *         u_n_minus_1[d] = u_n[d];
 * 
 *         vel_it_matrix[d].copy_from(vel_Laplace_plus_Mass);
 *         vel_it_matrix[d].add(1., vel_Advection);
 * 
 *         vel_exact.set_component(d);
 *         boundary_values.clear();
 *         for (const auto &boundary_id : boundary_ids)
 *           {
 *             switch (boundary_id)
 *               {
 *                 case 1:
 *                   VectorTools::interpolate_boundary_values(
 *                     dof_handler_velocity,
 *                     boundary_id,
 *                     Functions::ZeroFunction<dim>(),
 *                     boundary_values);
 *                   break;
 *                 case 2:
 *                   VectorTools::interpolate_boundary_values(dof_handler_velocity,
 *                                                            boundary_id,
 *                                                            vel_exact,
 *                                                            boundary_values);
 *                   break;
 *                 case 3:
 *                   if (d != 0)
 *                     VectorTools::interpolate_boundary_values(
 *                       dof_handler_velocity,
 *                       boundary_id,
 *                       Functions::ZeroFunction<dim>(),
 *                       boundary_values);
 *                   break;
 *                 case 4:
 *                   VectorTools::interpolate_boundary_values(
 *                     dof_handler_velocity,
 *                     boundary_id,
 *                     Functions::ZeroFunction<dim>(),
 *                     boundary_values);
 *                   break;
 *                 default:
 *                   Assert(false, ExcNotImplemented());
 *               }
 *           }
 *         MatrixTools::apply_boundary_values(boundary_values,
 *                                            vel_it_matrix[d],
 *                                            u_n[d],
 *                                            force[d]);
 *       }
 * 
 * 
 *     Threads::TaskGroup<void> tasks;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         if (reinit_prec)
 *           prec_velocity[d].initialize(vel_it_matrix[d],
 *                                       SparseILU<double>::AdditionalData(
 *                                         vel_diag_strength, vel_off_diagonals));
 *         tasks += Threads::new_task(
 *           &NavierStokesProjection<dim>::diffusion_component_solve, *this, d);
 *       }
 *     tasks.join_all();
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void
 *   NavierStokesProjection<dim>::diffusion_component_solve(const unsigned int d)
 *   {
 *     SolverControl solver_control(vel_max_its, vel_eps * force[d].l2_norm());
 *     SolverGMRES<Vector<double>> gmres(
 *       solver_control,
 *       SolverGMRES<Vector<double>>::AdditionalData(vel_Krylov_size));
 *     gmres.solve(vel_it_matrix[d], u_n[d], force[d], prec_velocity[d]);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionassemble_advection_termcode"></a> 
 * <h4> <code>NavierStokesProjection::assemble_advection_term</code> </h4>
 * 

 * 
 * The following few functions deal with assembling the advection terms, which
 * is the part of the system matrix for the diffusion step that changes at
 * every time step. As mentioned above, we will run the assembly loop over all
 * cells in %parallel, using the WorkStream class and other
 * facilities as described in the documentation module on @ref threads.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::assemble_advection_term()
 *   {
 *     vel_Advection = 0.;
 *     AdvectionPerTaskData data(fe_velocity.n_dofs_per_cell());
 *     AdvectionScratchData scratch(fe_velocity,
 *                                  quadrature_velocity,
 *                                  update_values | update_JxW_values |
 *                                    update_gradients);
 *     WorkStream::run(
 *       dof_handler_velocity.begin_active(),
 *       dof_handler_velocity.end(),
 *       *this,
 *       &NavierStokesProjection<dim>::assemble_one_cell_of_advection,
 *       &NavierStokesProjection<dim>::copy_advection_local_to_global,
 *       scratch,
 *       data);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void NavierStokesProjection<dim>::assemble_one_cell_of_advection(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     AdvectionScratchData &                                scratch,
 *     AdvectionPerTaskData &                                data)
 *   {
 *     scratch.fe_val.reinit(cell);
 *     cell->get_dof_indices(data.local_dof_indices);
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         scratch.fe_val.get_function_values(u_star[d], scratch.u_star_tmp);
 *         for (unsigned int q = 0; q < scratch.nqp; ++q)
 *           scratch.u_star_local[q](d) = scratch.u_star_tmp[q];
 *       }
 * 
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         scratch.fe_val.get_function_gradients(u_star[d], scratch.grad_u_star);
 *         for (unsigned int q = 0; q < scratch.nqp; ++q)
 *           {
 *             if (d == 0)
 *               scratch.u_star_tmp[q] = 0.;
 *             scratch.u_star_tmp[q] += scratch.grad_u_star[q][d];
 *           }
 *       }
 * 
 *     data.local_advection = 0.;
 *     for (unsigned int q = 0; q < scratch.nqp; ++q)
 *       for (unsigned int i = 0; i < scratch.dpc; ++i)
 *         for (unsigned int j = 0; j < scratch.dpc; ++j)
 *           data.local_advection(i, j) += (scratch.u_star_local[q] *            
 *                                            scratch.fe_val.shape_grad(j, q) *  
 *                                            scratch.fe_val.shape_value(i, q)   
 *                                          +                                    
 *                                          0.5 *                                
 *                                            scratch.u_star_tmp[q] *            
 *                                            scratch.fe_val.shape_value(i, q) * 
 *                                            scratch.fe_val.shape_value(j, q))  
 *                                         * scratch.fe_val.JxW(q);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void NavierStokesProjection<dim>::copy_advection_local_to_global(
 *     const AdvectionPerTaskData &data)
 *   {
 *     for (unsigned int i = 0; i < fe_velocity.n_dofs_per_cell(); ++i)
 *       for (unsigned int j = 0; j < fe_velocity.n_dofs_per_cell(); ++j)
 *         vel_Advection.add(data.local_dof_indices[i],
 *                           data.local_dof_indices[j],
 *                           data.local_advection(i, j));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionprojection_stepcode"></a> 
 * <h4><code>NavierStokesProjection::projection_step</code></h4>
 * 

 * 
 * This implements the projection step:
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::projection_step(const bool reinit_prec)
 *   {
 *     pres_iterative.copy_from(pres_Laplace);
 * 
 *     pres_tmp = 0.;
 *     for (unsigned d = 0; d < dim; ++d)
 *       pres_Diff[d].Tvmult_add(pres_tmp, u_n[d]);
 * 
 *     phi_n_minus_1 = phi_n;
 * 
 *     static std::map<types::global_dof_index, double> bval;
 *     if (reinit_prec)
 *       VectorTools::interpolate_boundary_values(dof_handler_pressure,
 *                                                3,
 *                                                Functions::ZeroFunction<dim>(),
 *                                                bval);
 * 
 *     MatrixTools::apply_boundary_values(bval, pres_iterative, phi_n, pres_tmp);
 * 
 *     if (reinit_prec)
 *       prec_pres_Laplace.initialize(pres_iterative,
 *                                    SparseILU<double>::AdditionalData(
 *                                      vel_diag_strength, vel_off_diagonals));
 * 
 *     SolverControl solvercontrol(vel_max_its, vel_eps * pres_tmp.l2_norm());
 *     SolverCG<Vector<double>> cg(solvercontrol);
 *     cg.solve(pres_iterative, phi_n, pres_tmp, prec_pres_Laplace);
 * 
 *     phi_n *= 1.5 / dt;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionupdate_pressurecode"></a> 
 * <h4> <code>NavierStokesProjection::update_pressure</code> </h4>
 * 

 * 
 * This is the pressure update step of the projection method. It implements
 * the standard formulation of the method, that is @f[ p^{n+1} = p^n +
 * \phi^{n+1}, @f] or the rotational form, which is @f[ p^{n+1} = p^n +
 * \phi^{n+1} - \frac{1}{Re} \nabla\cdot u^{n+1}. @f]
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::update_pressure(const bool reinit_prec)
 *   {
 *     pres_n_minus_1 = pres_n;
 *     switch (type)
 *       {
 *         case RunTimeParameters::Method::standard:
 *           pres_n += phi_n;
 *           break;
 *         case RunTimeParameters::Method::rotational:
 *           if (reinit_prec)
 *             prec_mass.initialize(pres_Mass);
 *           pres_n = pres_tmp;
 *           prec_mass.solve(pres_n);
 *           pres_n.sadd(1. / Re, 1., pres_n_minus_1);
 *           pres_n += phi_n;
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       };
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNavierStokesProjectionoutput_resultscode"></a> 
 * <h4> <code>NavierStokesProjection::output_results</code> </h4>
 * 

 * 
 * This method plots the current solution. The main difficulty is that we want
 * to create a single output file that contains the data for all velocity
 * components, the pressure, and also the vorticity of the flow. On the other
 * hand, velocities and the pressure live on separate DoFHandler objects, and
 * so can't be written to the same file using a single DataOut object. As a
 * consequence, we have to work a bit harder to get the various pieces of data
 * into a single DoFHandler object, and then use that to drive graphical
 * output.
 *   

 * 
 * We will not elaborate on this process here, but rather refer to step-32,
 * where a similar procedure is used (and is documented) to create a joint
 * DoFHandler object for all variables.
 *   

 * 
 * Let us also note that we here compute the vorticity as a scalar quantity in
 * a separate function, using the $L^2$ projection of the quantity
 * $\text{curl} u$ onto the finite element space used for the components of
 * the velocity. In principle, however, we could also have computed as a
 * pointwise quantity from the velocity, and do so through the
 * DataPostprocessor mechanism discussed in step-29 and step-33.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::output_results(const unsigned int step)
 *   {
 *     assemble_vorticity((step == 1));
 *     const FESystem<dim> joint_fe(
 *       fe_velocity, dim, fe_pressure, 1, fe_velocity, 1);
 *     DoFHandler<dim> joint_dof_handler(triangulation);
 *     joint_dof_handler.distribute_dofs(joint_fe);
 *     Assert(joint_dof_handler.n_dofs() ==
 *              ((dim + 1) * dof_handler_velocity.n_dofs() +
 *               dof_handler_pressure.n_dofs()),
 *            ExcInternalError());
 *     Vector<double> joint_solution(joint_dof_handler.n_dofs());
 *     std::vector<types::global_dof_index> loc_joint_dof_indices(
 *       joint_fe.n_dofs_per_cell()),
 *       loc_vel_dof_indices(fe_velocity.n_dofs_per_cell()),
 *       loc_pres_dof_indices(fe_pressure.n_dofs_per_cell());
 *     typename DoFHandler<dim>::active_cell_iterator
 *       joint_cell = joint_dof_handler.begin_active(),
 *       joint_endc = joint_dof_handler.end(),
 *       vel_cell   = dof_handler_velocity.begin_active(),
 *       pres_cell  = dof_handler_pressure.begin_active();
 *     for (; joint_cell != joint_endc; ++joint_cell, ++vel_cell, ++pres_cell)
 *       {
 *         joint_cell->get_dof_indices(loc_joint_dof_indices);
 *         vel_cell->get_dof_indices(loc_vel_dof_indices);
 *         pres_cell->get_dof_indices(loc_pres_dof_indices);
 *         for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i)
 *           switch (joint_fe.system_to_base_index(i).first.first)
 *             {
 *               case 0:
 *                 Assert(joint_fe.system_to_base_index(i).first.second < dim,
 *                        ExcInternalError());
 *                 joint_solution(loc_joint_dof_indices[i]) =
 *                   u_n[joint_fe.system_to_base_index(i).first.second](
 *                     loc_vel_dof_indices[joint_fe.system_to_base_index(i)
 *                                           .second]);
 *                 break;
 *               case 1:
 *                 Assert(joint_fe.system_to_base_index(i).first.second == 0,
 *                        ExcInternalError());
 *                 joint_solution(loc_joint_dof_indices[i]) =
 *                   pres_n(loc_pres_dof_indices[joint_fe.system_to_base_index(i)
 *                                                 .second]);
 *                 break;
 *               case 2:
 *                 Assert(joint_fe.system_to_base_index(i).first.second == 0,
 *                        ExcInternalError());
 *                 joint_solution(loc_joint_dof_indices[i]) = rot_u(
 *                   loc_vel_dof_indices[joint_fe.system_to_base_index(i).second]);
 *                 break;
 *               default:
 *                 Assert(false, ExcInternalError());
 *             }
 *       }
 *     std::vector<std::string> joint_solution_names(dim, "v");
 *     joint_solution_names.emplace_back("p");
 *     joint_solution_names.emplace_back("rot_u");
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(joint_dof_handler);
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       component_interpretation(
 *         dim + 2, DataComponentInterpretation::component_is_part_of_vector);
 *     component_interpretation[dim] =
 *       DataComponentInterpretation::component_is_scalar;
 *     component_interpretation[dim + 1] =
 *       DataComponentInterpretation::component_is_scalar;
 *     data_out.add_data_vector(joint_solution,
 *                              joint_solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              component_interpretation);
 *     data_out.build_patches(deg + 1);
 *     std::ofstream output("solution-" + Utilities::int_to_string(step, 5) +
 *                          ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Following is the helper function that computes the vorticity by projecting
 * the term $\text{curl} u$ onto the finite element space used for the
 * components of the velocity. The function is only called whenever we
 * generate graphical output, so not very often, and as a consequence we
 * didn't bother parallelizing it using the WorkStream concept as we do for
 * the other assembly functions. That should not be overly complicated,
 * however, if needed. Moreover, the implementation that we have here only
 * works for 2d, so we bail if that is not the case.
 * 
 * @code
 *   template <int dim>
 *   void NavierStokesProjection<dim>::assemble_vorticity(const bool reinit_prec)
 *   {
 *     Assert(dim == 2, ExcNotImplemented());
 *     if (reinit_prec)
 *       prec_vel_mass.initialize(vel_Mass);
 * 
 *     FEValues<dim>      fe_val_vel(fe_velocity,
 *                              quadrature_velocity,
 *                              update_gradients | update_JxW_values |
 *                                update_values);
 *     const unsigned int dpc = fe_velocity.n_dofs_per_cell(),
 *                        nqp = quadrature_velocity.size();
 *     std::vector<types::global_dof_index> ldi(dpc);
 *     Vector<double>                       loc_rot(dpc);
 * 
 *     std::vector<Tensor<1, dim>> grad_u1(nqp), grad_u2(nqp);
 *     rot_u = 0.;
 * 
 *     for (const auto &cell : dof_handler_velocity.active_cell_iterators())
 *       {
 *         fe_val_vel.reinit(cell);
 *         cell->get_dof_indices(ldi);
 *         fe_val_vel.get_function_gradients(u_n[0], grad_u1);
 *         fe_val_vel.get_function_gradients(u_n[1], grad_u2);
 *         loc_rot = 0.;
 *         for (unsigned int q = 0; q < nqp; ++q)
 *           for (unsigned int i = 0; i < dpc; ++i)
 *             loc_rot(i) += (grad_u2[q][0] - grad_u1[q][1]) * 
 *                           fe_val_vel.shape_value(i, q) *    
 *                           fe_val_vel.JxW(q);
 * 
 *         for (unsigned int i = 0; i < dpc; ++i)
 *           rot_u(ldi[i]) += loc_rot(i);
 *       }
 * 
 *     prec_vel_mass.solve(rot_u);
 *   }
 * } // namespace Step35
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3> The main function </h3>
 * 

 * 
 * The main function looks very much like in all the other tutorial programs, so
 * there is little to comment on here:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step35;
 * 
 *       RunTimeParameters::Data_Storage data;
 *       data.read_data("parameter-file.prm");
 * 
 *       deallog.depth_console(data.verbose ? 2 : 0);
 * 
 *       NavierStokesProjection<2> test(data);
 *       test.run(data.verbose, data.output_interval);
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
 *     }
 *   std::cout << "----------------------------------------------------"
 *             << std::endl
 *             << "Apparently everything went fine!" << std::endl
 *             << "Don't forget to brush your teeth :-)" << std::endl
 *             << std::endl;
 *   return 0;
 * }
 * @endcode
<a name="results"></a>
<a name="Results"></a><h1>Results</h1>


<a name="Re100"></a>
<a name="Re100"></a><h3> Re = 100 </h3>


We run the code with the following <code>parameter-file.prm</code>, which can be found in the
same directory as the source:
@verbatim
  # First a global definition
  # the type of method we want to use
  set Method_Form = rotational

  subsection Physical data
    # In this subsection we declare the physical data
    # The initial and final time, and the Reynolds number
    set initial_time = 0.
    set final_time   = 25.
    set Reynolds     = 100
  end

  subsection Time step data
    # In this subsection we declare the data that is to be used for time discretization,
    # i.e. the time step dt
    set dt = 5e-3
  end

  subsection Space discretization
    # In this subsection we declare the data that is relevant to the space discretization
    # we set the number of global refines the triangulation must have
    # and the degree k of the pair Q_(k+1)--Q_k of velocity--pressure finite element spaces
    set n_of_refines = 3
    set pressure_fe_degree = 1
  end

  subsection Data solve velocity
    # In this section we declare the parameters that are going to control the solution process
    # for the velocity.
    set max_iterations = 1000  # maximal number of iterations that GMRES must make
    set eps            = 1e-6  # stopping criterion
    set Krylov_size    = 30    # size of the Krylov subspace to be used in GMRES
    set off_diagonals  = 60    # number of off diagonals that ILU must compute
    set diag_strength  = 0.01  # diagonal strengthening value
    set update_prec    = 10    # this number indicates how often the preconditioner must be updated
  end

  #The output frequency
  set output = 50

  #Finally we set the verbosity level
  set verbose = false
@endverbatim

Since the <code>verbose</code> parameter is set to <code>false</code>,
we do not get any kind of output besides the number of the time step
the program is currently working on.
If we were to set it to <code>true</code> we would get information on what the program is doing and
how many steps each iterative process had to make to converge, etc.

Let us plot the obtained results for $t=1,5,12,20,25$ (i.e. time steps
200, 1000, 2400, 4000, and 5000), where in the left column we show the
vorticity and in the right the velocity field:

<table>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.vorticity.0.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.velocity.0.png" alt=""> </td>
  </tr>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.vorticity.1.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.velocity.1.png" alt=""> </td>
  </tr>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.vorticity.2.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.velocity.2.png" alt=""> </td>
  </tr>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.vorticity.3.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.velocity.3.png" alt=""> </td>
  </tr>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.vorticity.4.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_100.velocity.4.png" alt=""> </td>
  </tr>
</table>

The images show nicely the development and extension of a vortex chain
behind the obstacles, with the sign of the vorticity indicating
whether this is a left or right turning vortex.


<a name="Re500"></a>
<a name="Re500"></a><h3> Re = 500 </h3>


We can change the Reynolds number, $Re$, in the parameter file to a
value of $500$. Doing so, and reducing the time step somewhat as well,
yields the following images at times $t=20,40$:

<table>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.vorticity.0.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.velocity.0.png" alt=""> </td>
  </tr>
  <tr>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.vorticity.1.png" alt=""> </td>
    <td> <img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.velocity.1.png" alt=""> </td>
  </tr>
</table>

For this larger Reynolds number, we observe unphysical oscillations, especially
for the vorticity. The discretization scheme has now difficulties in correctly
resolving the flow, which should still be laminar and well-organized.
These phenomena are typical of discretization schemes that lack robustness
in under-resolved scenarios, where under-resolved means that the Reynolds
number computed with the mesh size instead of the physical dimensions of
the geometry is large. We look at a zoom at the region behind the obstacle, and
the mesh size we have there:


<img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.zoom.png" alt="">

We can easily test our hypothesis by re-running the simulation with one more
mesh refinement set in the parameter file:

<img src="https://www.dealii.org/images/steps/developer/step-35.Re_500.zoom_2.png" alt="">

Indeed, the vorticity field now looks much smoother. While we can expect that
further refining the mesh will suppress the remaining oscillations as well,
one should take measures to obtain a robust scheme in the limit of coarse
resolutions, as described below.


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


This program can be extended in the following directions:
<ul>
  <li> Adaptive mesh refinement: As we have seen, we computed everything on a single fixed mesh.
  Using adaptive mesh refinement can lead to increased accuracy while not significantly increasing the
  computational time.

  <li> Adaptive time-stepping: Although there apparently is currently no theory about
  projection methods with variable time step,
  practice shows that they perform very well.

  <li> High Reynolds %numbers: As we can see from the results, increasing the Reynolds number changes significantly
  the behavior of the discretization scheme. Using well-known stabilization techniques we could be able to
  compute the flow in this, or many other problems, when the Reynolds number is very large and where computational
  costs demand spatial resolutions for which the flow is only marginally resolved, especially for 3D turbulent
  flows.

  <li> Variable density incompressible flows: There are projection-like methods for the case of incompressible
  flows with variable density. Such flows play a role if fluids of different
  density mix, for example fresh water and salt water, or alcohol and water.

  <li> Compressible Navier-Stokes equations: These equations are relevant for
  cases where
  velocities are high enough so that the fluid becomes compressible, but not
  fast enough that we get into a regime where viscosity becomes negligible
  and the Navier-Stokes equations need to be replaced by the hyperbolic Euler
  equations of gas dynamics. Compressibility starts to become a factor if the
  velocity becomes greater than about one third of the speed of sound, so it
  is not a factor for almost all terrestrial vehicles. On the other hand,
  commercial jetliners fly at about 85 per cent of the speed of sound, and
  flow over the wings becomes significantly supersonic, a regime in which the
  compressible Navier-Stokes equations are not applicable any more
  either. There are significant applications for the range in between,
  however, such as for small aircraft or the fast trains in many European and
  East Asian countries.
</ul>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-35.cc"
*/
