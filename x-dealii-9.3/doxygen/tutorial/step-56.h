/**
@page step_56 The step-56 tutorial program
This tutorial depends on step-16, step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#StokesProblem"> Stokes Problem </a>
        <li><a href="#LinearSolverandPreconditioningIssues"> Linear Solver and Preconditioning Issues </a>
        <li><a href="#ReferenceSolution"> Reference Solution </a>
        <li><a href="#ComputingErrors"> Computing Errors </a>
        <li><a href="#DoFHandlers"> DoF Handlers </a>
        <li><a href="#DifferencesfromtheStep22tutorial"> Differences from the Step 22 tutorial </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#FunctionsforSolutionandRighthandside">Functions for Solution and Righthand side</a>
        <li><a href="#ASPECTBlockSchurPreconditioner">ASPECT BlockSchurPreconditioner</a>
        <li><a href="#TheStokesProblemclass">The StokesProblem class</a>
      <ul>
        <li><a href="#StokesProblemsetup_dofs">StokesProblem::setup_dofs</a>
        <li><a href="#StokesProblemassemble_system">StokesProblem::assemble_system</a>
        <li><a href="#StokesProblemassemble_multigrid">StokesProblem::assemble_multigrid</a>
        <li><a href="#StokesProblemsolve">StokesProblem::solve</a>
        <li><a href="#StokesProblemprocess_solution">StokesProblem::process_solution</a>
        <li><a href="#StokesProblemoutput_results">StokesProblem::output_results</a>
        <li><a href="#StokesProblemrun">StokesProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Errors"> Errors </a>
        <li><a href="#TimingResults"> Timing Results </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Checkhigherorderdiscretizations"> Check higher order discretizations </a>
        <li><a href="#Comparewithcheappreconditioner"> Compare with cheap preconditioner </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<i>This program was contributed by Ryan Grove and Timo Heister.

This material is based upon work partially supported by National Science
Foundation grant DMS1522191 and the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-0949446 and The University of California-Davis.

The authors would like to thank the Isaac Newton Institute for
Mathematical Sciences, Cambridge, for support and hospitality during
the programme Melt in the Mantle where work on this tutorial was
undertaken. This work was supported by EPSRC grant no EP/K032208/1.
</i>

@dealiiTutorialDOI{10.5281/zenodo.400995,https://zenodo.org/badge/DOI/10.5281/zenodo.400995.svg}

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


<a name="StokesProblem"></a><h3> Stokes Problem </h3>


The purpose of this tutorial is to create an efficient linear solver for the
Stokes equation and compare it to alternative approaches.  Here, we will use
FGMRES with geometric multigrid as a preconditioner velocity block, and we
will show in the results section that this is a fundamentally better approach
than the linear solvers used in step-22 (including the scheme described in
"Possible Extensions").  Fundamentally, this is because only with multigrid it
is possible to get $O(n)$ solve time, where $n$ is the number of unknowns of
the linear system. Using the Timer class, we collect some statistics to
compare setup times, solve times, and number of iterations. We also compute
errors to make sure that what we have implemented is correct.

Let $u \in H_0^1 = \{ u \in H^1(\Omega), u|_{\partial \Omega} = 0 \}$
and $p \in L_*^2 = \{ p \in L^2(\Omega), \int_\Omega p = 0
\}$. The Stokes equations read as follows in non-dimensionalized form:

@f{eqnarray*}
 - 2 \text{div} \frac {1}{2} \left[ (\nabla \textbf{u})
 + (\nabla \textbf{u})^T\right] + \nabla p & =& f \\
 - \nabla \cdot u &=& 0
@f}

Note that we are using the deformation tensor instead of $\Delta u$ (a
detailed description of the difference between the two can be found in
step-22, but in summary, the deformation tensor is more physical as
well as more expensive).

<a name="LinearSolverandPreconditioningIssues"></a><h3> Linear %Solver and Preconditioning Issues </h3>


The weak form of
the discrete equations naturally leads to the following linear system
for the nodal values of the velocity and pressure fields:
@f{eqnarray*}
\left(\begin{array}{cc} A & B^T \\ B & 0
\end{array}\right) \left(\begin{array}{c} U \\ P \end{array}\right) =
\left(\begin{array}{c} F \\ 0 \end{array}\right).
@f}

Our goal is to compare several solution approaches.  While step-22
solves the linear system using a "Schur complement approach" in two
separate steps, we instead attack the
block system at once using FMGRES with an efficient
preconditioner, in the spirit of the approach outlined in the "Results"
section of step-22. The idea is as follows: if we find a block
preconditioner $P$ such that the matrix

@f{eqnarray*}
\left(\begin{array}{cc} A & B^T \\ B & 0 \end{array}\right) P^{-1}
@f}

is simple, then an iterative solver with that preconditioner will
converge in a few iterations. Notice that we are doing right
preconditioning here.  Using the Schur complement $S=BA^{-1}B^T$,
we find that

@f{eqnarray*}
P^{-1} = \left(\begin{array}{cc} A & B^T \\ 0 &
 S \end{array}\right)^{-1}
@f}

is a good choice. Let $\widetilde{A^{-1}}$ be an approximation of $A^{-1}$
and $\widetilde{S^{-1}}$ of $S^{-1}$, we see
@f{eqnarray*}
P^{-1} =
\left(\begin{array}{cc} A^{-1} & 0 \\ 0 & I \end{array}\right)
\left(\begin{array}{cc} I & B^T \\ 0 & -I \end{array}\right)
\left(\begin{array}{cc} I & 0 \\ 0 & S^{-1} \end{array}\right)
\approx
\left(\begin{array}{cc} \widetilde{A^{-1}} & 0 \\ 0 & I \end{array}\right)
\left(\begin{array}{cc} I & B^T \\ 0 & -I \end{array}\right)
\left(\begin{array}{cc} I & 0 \\ 0 & \widetilde{S^{-1}} \end{array}\right).
  @f}

Since $P$ is aimed to be a preconditioner only, we shall use
the approximations on the right in the equation above.

As discussed in step-22, $-M_p^{-1}=:\widetilde{S^{-1}} \approx
S^{-1}$, where $M_p$ is the pressure mass matrix and is solved approximately by using CG
with ILU as a preconditioner, and $\widetilde{A^{-1}}$ is obtained by one of
multiple methods: solving a linear system with CG and ILU as
preconditioner, just using one application of an ILU, solving a linear
system with CG and GMG (Geometric
Multigrid as described in step-16) as a preconditioner, or just performing a single V-cycle
of GMG.

As a comparison, instead of FGMRES, we also use the direct solver
UMFPACK on the whole system to compare our results with.  If you want to use
a direct solver (like UMFPACK), the system needs to be invertible. To avoid
the one dimensional null space given by the constant pressures, we fix the first pressure unknown
 to zero. This is not necessary for the iterative solvers.


<a name="ReferenceSolution"></a><h3> Reference Solution </h3>


The test problem is a "Manufactured Solution" (see step-7 for
details), and we choose $u=(u_1,u_2,u_3)=(2\sin (\pi x), - \pi y \cos
(\pi x),- \pi z \cos (\pi x))$ and $p = \sin (\pi x)\cos (\pi y)\sin
(\pi z)$.
We apply Dirichlet boundary conditions for the velocity on the whole
boundary of the domain $\Omega=[0,1]\times[0,1]\times[0,1]$.
To enforce the boundary conditions we can just use our reference solution.

If you look up in the deal.II manual what is needed to create a class
derived from <code>Function@<dim@></code>, you will find that this
class has numerous @p virtual functions, including
Function::value(), Function::vector_value(), Function::value_list(),
etc., all of which can be overloaded.  Different parts of deal.II
will require different ones of these particular
functions. This can be confusing at first, but luckily the only thing
you actually have to implement is @p value().  The other virtual
functions in the Function class have default
implementations inside that will call your implementation of @p value
by default.

Notice that our reference solution fulfills $\nabla \cdot u = 0$. In
addition, the pressure is chosen to have a mean value of zero.  For
the "Method of Manufactured Solutions" of step-7, we need to find $\bf
f$ such that:

@f{align*}
{\bf f} =   - 2 \text{div} \frac {1}{2} \left[ (\nabla \textbf{u}) + (\nabla \textbf{u})^T\right] + \nabla p.
@f}

Using the reference solution above, we obtain:

@f{eqnarray*}
{\bf f} &=& (2 \pi^2 \sin (\pi x),- \pi^3 y \cos(\pi
x),- \pi^3 z \cos(\pi x))\\ & & + (\pi \cos(\pi x) \cos(\pi y)
\sin(\pi z) ,- \pi \sin(\pi y) \sin(\pi x) \sin(\pi z), \pi \cos(\pi
z) \sin(\pi x) \cos(\pi y)) @f}

<a name="ComputingErrors"></a><h3> Computing Errors </h3>


Because we do not enforce the mean
pressure to be zero for our numerical solution in the linear system,
we need to post process the solution after solving. To do this we use
the VectorTools::compute_mean_value() function to compute the mean value
of the pressure to subtract it from the pressure.


<a name="DoFHandlers"></a><h3> DoF Handlers </h3>


The way we implement geometric multigrid here only executes it on the
velocity variables (i.e., the $A$ matrix described above) but not the
pressure. One could implement this in different ways, including one in
which one considers all coarse grid operations as acting on $2\times
2$ block systems where we only consider the top left
block. Alternatively, we can implement things by really only
considering a linear system on the velocity part of the overall finite
element discretization. The latter is the way we want to use here.

To implement this, one would need to be able to ask questions such as
"May I have just part of a DoFHandler?". This is not possible at the
time when this program was written, so in order to answer this request
for our needs, we simply create a separate, second DoFHandler for just the
velocities. We then build linear systems for the multigrid
preconditioner based on only this second DoFHandler, and simply
transfer the first block of (overall) vectors into corresponding
vectors for the entire second DoFHandler. To make this work, we have
to assure that the <i>order</i> in which the (velocity) degrees of freedom are
ordered in the two DoFHandler objects is the same. This is in fact the
case by first distributing degrees of freedom on both, and then using
the same sequence of DoFRenumbering operations on both.


<a name="DifferencesfromtheStep22tutorial"></a><h3> Differences from the Step 22 tutorial </h3>


The main difference between step-56 and step-22 is that we use block
solvers instead of the Schur Complement approach used in
step-22. Details of this approach can be found under the "Block Schur
complement preconditioner" subsection of the "Possible Extensions"
section of step-22. For the preconditioner of the velocity block, we
borrow a class from <a href="https://aspect.geodynamics.org">ASPECT</a>
called @p BlockSchurPreconditioner that has the option to solve for
the inverse of $A$ or just apply one preconditioner sweep for it
instead, which provides us with an expensive and cheap approach,
respectively.
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
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/block_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_gmres.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/grid_refinement.h>
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
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * #include <deal.II/lac/sparse_direct.h>
 * 
 * #include <deal.II/lac/sparse_ilu.h>
 * #include <deal.II/grid/grid_out.h>
 * 
 * @endcode
 * 
 * We need to include the following file to do timings:
 * 
 * @code
 * #include <deal.II/base/timer.h>
 * 
 * @endcode
 * 
 * This includes the files necessary for us to use geometric Multigrid
 * 
 * @code
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * namespace Step56
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * In order to make it easy to switch between the different solvers that are
 * being used, we declare an enum that can be passed as an argument to the
 * constructor of the main class.
 * 
 * @code
 *   enum class SolverType
 *   {
 *     FGMRES_ILU,
 *     FGMRES_GMG,
 *     UMFPACK
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionsforSolutionandRighthandside"></a> 
 * <h3>Functions for Solution and Righthand side</h3>
 *   

 * 
 * The class Solution is used to define the boundary conditions and to
 * compute errors of the numerical solution. Note that we need to define the
 * values and gradients in order to compute L2 and H1 errors. Here we
 * decided to separate the implementations for 2d and 3d using template
 * specialization.
 *   

 * 
 * Note that the first dim components are the velocity components
 * and the last is the pressure.
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
 *   {
 *   public:
 *     Solution()
 *       : Function<dim>(dim + 1)
 *     {}
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *     virtual Tensor<1, dim>
 *     gradient(const Point<dim> & p,
 *              const unsigned int component = 0) const override;
 *   };
 * 
 *   template <>
 *   double Solution<2>::value(const Point<2> &   p,
 *                             const unsigned int component) const
 *   {
 *     Assert(component <= 2 + 1, ExcIndexRange(component, 0, 2 + 1));
 * 
 *     using numbers::PI;
 *     const double x = p(0);
 *     const double y = p(1);
 * 
 *     if (component == 0)
 *       return sin(PI * x);
 *     if (component == 1)
 *       return -PI * y * cos(PI * x);
 *     if (component == 2)
 *       return sin(PI * x) * cos(PI * y);
 * 
 *     return 0;
 *   }
 * 
 *   template <>
 *   double Solution<3>::value(const Point<3> &   p,
 *                             const unsigned int component) const
 *   {
 *     Assert(component <= 3 + 1, ExcIndexRange(component, 0, 3 + 1));
 * 
 *     using numbers::PI;
 *     const double x = p(0);
 *     const double y = p(1);
 *     const double z = p(2);
 * 
 *     if (component == 0)
 *       return 2.0 * sin(PI * x);
 *     if (component == 1)
 *       return -PI * y * cos(PI * x);
 *     if (component == 2)
 *       return -PI * z * cos(PI * x);
 *     if (component == 3)
 *       return sin(PI * x) * cos(PI * y) * sin(PI * z);
 * 
 *     return 0;
 *   }
 * 
 * @endcode
 * 
 * Note that for the gradient we need to return a Tensor<1,dim>
 * 
 * @code
 *   template <>
 *   Tensor<1, 2> Solution<2>::gradient(const Point<2> &   p,
 *                                      const unsigned int component) const
 *   {
 *     Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));
 * 
 *     using numbers::PI;
 *     const double x = p(0);
 *     const double y = p(1);
 * 
 *     Tensor<1, 2> return_value;
 *     if (component == 0)
 *       {
 *         return_value[0] = PI * cos(PI * x);
 *         return_value[1] = 0.0;
 *       }
 *     else if (component == 1)
 *       {
 *         return_value[0] = y * PI * PI * sin(PI * x);
 *         return_value[1] = -PI * cos(PI * x);
 *       }
 *     else if (component == 2)
 *       {
 *         return_value[0] = PI * cos(PI * x) * cos(PI * y);
 *         return_value[1] = -PI * sin(PI * x) * sin(PI * y);
 *       }
 * 
 *     return return_value;
 *   }
 * 
 *   template <>
 *   Tensor<1, 3> Solution<3>::gradient(const Point<3> &   p,
 *                                      const unsigned int component) const
 *   {
 *     Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));
 * 
 *     using numbers::PI;
 *     const double x = p(0);
 *     const double y = p(1);
 *     const double z = p(2);
 * 
 *     Tensor<1, 3> return_value;
 *     if (component == 0)
 *       {
 *         return_value[0] = 2 * PI * cos(PI * x);
 *         return_value[1] = 0.0;
 *         return_value[2] = 0.0;
 *       }
 *     else if (component == 1)
 *       {
 *         return_value[0] = y * PI * PI * sin(PI * x);
 *         return_value[1] = -PI * cos(PI * x);
 *         return_value[2] = 0.0;
 *       }
 *     else if (component == 2)
 *       {
 *         return_value[0] = z * PI * PI * sin(PI * x);
 *         return_value[1] = 0.0;
 *         return_value[2] = -PI * cos(PI * x);
 *       }
 *     else if (component == 3)
 *       {
 *         return_value[0] = PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
 *         return_value[1] = -PI * sin(PI * x) * sin(PI * y) * sin(PI * z);
 *         return_value[2] = PI * sin(PI * x) * cos(PI * y) * cos(PI * z);
 *       }
 * 
 *     return return_value;
 *   }
 * 
 * @endcode
 * 
 * Implementation of $f$. See the introduction for more information.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     RightHandSide()
 *       : Function<dim>(dim + 1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 *   template <>
 *   double RightHandSide<2>::value(const Point<2> &   p,
 *                                  const unsigned int component) const
 *   {
 *     Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));
 * 
 *     using numbers::PI;
 *     double x = p(0);
 *     double y = p(1);
 *     if (component == 0)
 *       return PI * PI * sin(PI * x) + PI * cos(PI * x) * cos(PI * y);
 *     if (component == 1)
 *       return -PI * PI * PI * y * cos(PI * x) - PI * sin(PI * y) * sin(PI * x);
 *     if (component == 2)
 *       return 0;
 * 
 *     return 0;
 *   }
 * 
 *   template <>
 *   double RightHandSide<3>::value(const Point<3> &   p,
 *                                  const unsigned int component) const
 *   {
 *     Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));
 * 
 *     using numbers::PI;
 *     double x = p(0);
 *     double y = p(1);
 *     double z = p(2);
 *     if (component == 0)
 *       return 2 * PI * PI * sin(PI * x) +
 *              PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
 *     if (component == 1)
 *       return -PI * PI * PI * y * cos(PI * x) +
 *              PI * (-1) * sin(PI * y) * sin(PI * x) * sin(PI * z);
 *     if (component == 2)
 *       return -PI * PI * PI * z * cos(PI * x) +
 *              PI * cos(PI * z) * sin(PI * x) * cos(PI * y);
 *     if (component == 3)
 *       return 0;
 * 
 *     return 0;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ASPECTBlockSchurPreconditioner"></a> 
 * <h3>ASPECT BlockSchurPreconditioner</h3>
 * 

 * 
 * In the following, we will implement a preconditioner that expands
 * on the ideas discussed in the Results section of step-22.
 * Specifically, we
 * 1. use an upper block-triangular preconditioner because we want to
 * use right preconditioning.
 * 2. optionally allow using an inner solver for the velocity block instead
 * of a single preconditioner application.
 * 3. do not use InverseMatrix but explicitly call SolverCG.
 * This approach is also used in the ASPECT code
 * (see https://aspect.geodynamics.org) that solves the Stokes equations in
 * the context of simulating convection in the earth mantle, and which
 * has been used to solve problems on many thousands of processors.
 *   

 * 
 * The bool flag @p do_solve_A in the constructor allows us to either
 * apply the preconditioner for the velocity block once or use an inner
 * iterative solver for a more accurate approximation instead.
 *   

 * 
 * Notice how we keep track of the sum of the inner iterations
 * (preconditioner applications).
 * 
 * @code
 *   template <class PreconditionerAType, class PreconditionerSType>
 *   class BlockSchurPreconditioner : public Subscriptor
 *   {
 *   public:
 *     BlockSchurPreconditioner(
 *       const BlockSparseMatrix<double> &system_matrix,
 *       const SparseMatrix<double> &     schur_complement_matrix,
 *       const PreconditionerAType &      preconditioner_A,
 *       const PreconditionerSType &      preconditioner_S,
 *       const bool                       do_solve_A);
 * 
 *     void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;
 * 
 *     mutable unsigned int n_iterations_A;
 *     mutable unsigned int n_iterations_S;
 * 
 *   private:
 *     const BlockSparseMatrix<double> &system_matrix;
 *     const SparseMatrix<double> &     schur_complement_matrix;
 *     const PreconditionerAType &      preconditioner_A;
 *     const PreconditionerSType &      preconditioner_S;
 * 
 *     const bool do_solve_A;
 *   };
 * 
 *   template <class PreconditionerAType, class PreconditionerSType>
 *   BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::
 *     BlockSchurPreconditioner(
 *       const BlockSparseMatrix<double> &system_matrix,
 *       const SparseMatrix<double> &     schur_complement_matrix,
 *       const PreconditionerAType &      preconditioner_A,
 *       const PreconditionerSType &      preconditioner_S,
 *       const bool                       do_solve_A)
 *     : n_iterations_A(0)
 *     , n_iterations_S(0)
 *     , system_matrix(system_matrix)
 *     , schur_complement_matrix(schur_complement_matrix)
 *     , preconditioner_A(preconditioner_A)
 *     , preconditioner_S(preconditioner_S)
 *     , do_solve_A(do_solve_A)
 *   {}
 * 
 * 
 * 
 *   template <class PreconditionerAType, class PreconditionerSType>
 *   void
 *   BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::vmult(
 *     BlockVector<double> &      dst,
 *     const BlockVector<double> &src) const
 *   {
 *     Vector<double> utmp(src.block(0));
 * 
 * @endcode
 * 
 * First solve with the approximation for S
 * 
 * @code
 *     {
 *       SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       dst.block(1) = 0.0;
 *       cg.solve(schur_complement_matrix,
 *                dst.block(1),
 *                src.block(1),
 *                preconditioner_S);
 * 
 *       n_iterations_S += solver_control.last_step();
 *       dst.block(1) *= -1.0;
 *     }
 * 
 * @endcode
 * 
 * Second, apply the top right block (B^T)
 * 
 * @code
 *     {
 *       system_matrix.block(0, 1).vmult(utmp, dst.block(1));
 *       utmp *= -1.0;
 *       utmp += src.block(0);
 *     }
 * 
 * @endcode
 * 
 * Finally, either solve with the top left block
 * or just apply one preconditioner sweep
 * 
 * @code
 *     if (do_solve_A == true)
 *       {
 *         SolverControl            solver_control(10000, utmp.l2_norm() * 1e-4);
 *         SolverCG<Vector<double>> cg(solver_control);
 * 
 *         dst.block(0) = 0.0;
 *         cg.solve(system_matrix.block(0, 0),
 *                  dst.block(0),
 *                  utmp,
 *                  preconditioner_A);
 * 
 *         n_iterations_A += solver_control.last_step();
 *       }
 *     else
 *       {
 *         preconditioner_A.vmult(dst.block(0), utmp);
 *         n_iterations_A += 1;
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="TheStokesProblemclass"></a> 
 * <h3>The StokesProblem class</h3>
 *   

 * 
 * This is the main class of the problem.
 * 
 * @code
 *   template <int dim>
 *   class StokesProblem
 *   {
 *   public:
 *     StokesProblem(const unsigned int pressure_degree,
 *                   const SolverType   solver_type);
 *     void run();
 * 
 *   private:
 *     void setup_dofs();
 *     void assemble_system();
 *     void assemble_multigrid();
 *     void solve();
 *     void compute_errors();
 *     void output_results(const unsigned int refinement_cycle) const;
 * 
 *     const unsigned int pressure_degree;
 *     const SolverType   solver_type;
 * 
 *     Triangulation<dim> triangulation;
 *     FESystem<dim>      velocity_fe;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 *     DoFHandler<dim>    velocity_dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     BlockSparsityPattern      sparsity_pattern;
 *     BlockSparseMatrix<double> system_matrix;
 *     SparseMatrix<double>      pressure_mass_matrix;
 * 
 *     BlockVector<double> solution;
 *     BlockVector<double> system_rhs;
 * 
 *     MGLevelObject<SparsityPattern>      mg_sparsity_patterns;
 *     MGLevelObject<SparseMatrix<double>> mg_matrices;
 *     MGLevelObject<SparseMatrix<double>> mg_interface_matrices;
 *     MGConstrainedDoFs                   mg_constrained_dofs;
 * 
 *     TimerOutput computing_timer;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   StokesProblem<dim>::StokesProblem(const unsigned int pressure_degree,
 *                                     const SolverType   solver_type)
 * 
 *     : pressure_degree(pressure_degree)
 *     , solver_type(solver_type)
 *     , triangulation(Triangulation<dim>::maximum_smoothing)
 *     ,
 * @endcode
 * 
 * Finite element for the velocity only:
 * 
 * @code
 *     velocity_fe(FE_Q<dim>(pressure_degree + 1), dim)
 *     ,
 * @endcode
 * 
 * Finite element for the whole system:
 * 
 * @code
 *     fe(velocity_fe, 1, FE_Q<dim>(pressure_degree), 1)
 *     , dof_handler(triangulation)
 *     , velocity_dof_handler(triangulation)
 *     , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsetup_dofs"></a> 
 * <h4>StokesProblem::setup_dofs</h4>
 * 

 * 
 * This function sets up the DoFHandler, matrices, vectors, and Multigrid
 * structures (if needed).
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::setup_dofs()
 *   {
 *     TimerOutput::Scope scope(computing_timer, "Setup");
 * 
 *     system_matrix.clear();
 *     pressure_mass_matrix.clear();
 * 
 * @endcode
 * 
 * The main DoFHandler only needs active DoFs, so we are not calling
 * distribute_mg_dofs() here
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * This block structure separates the dim velocity components from
 * the pressure component (used for reordering). Note that we have
 * 2 instead of dim+1 blocks like in step-22, because our FESystem
 * is nested and the dim velocity components appear as one block.
 * 
 * @code
 *     std::vector<unsigned int> block_component(2);
 *     block_component[0] = 0;
 *     block_component[1] = 1;
 * 
 * @endcode
 * 
 * Velocities start at component 0:
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 * @endcode
 * 
 * ILU behaves better if we apply a reordering to reduce fillin. There
 * is no advantage in doing this for the other solvers.
 * 
 * @code
 *     if (solver_type == SolverType::FGMRES_ILU)
 *       {
 *         TimerOutput::Scope ilu_specific(computing_timer, "(ILU specific)");
 *         DoFRenumbering::Cuthill_McKee(dof_handler);
 *       }
 * 
 * @endcode
 * 
 * This ensures that all velocities DoFs are enumerated before the
 * pressure unknowns. This allows us to use blocks for vectors and
 * matrices and allows us to get the same DoF numbering for
 * dof_handler and velocity_dof_handler.
 * 
 * @code
 *     DoFRenumbering::block_wise(dof_handler);
 * 
 *     if (solver_type == SolverType::FGMRES_GMG)
 *       {
 *         TimerOutput::Scope multigrid_specific(computing_timer,
 *                                               "(Multigrid specific)");
 *         TimerOutput::Scope setup_multigrid(computing_timer,
 *                                            "Setup - Multigrid");
 * 
 * @endcode
 * 
 * This distributes the active dofs and multigrid dofs for the
 * velocity space in a separate DoFHandler as described in the
 * introduction.
 * 
 * @code
 *         velocity_dof_handler.distribute_dofs(velocity_fe);
 *         velocity_dof_handler.distribute_mg_dofs();
 * 
 * @endcode
 * 
 * The following block of code initializes the MGConstrainedDofs
 * (using the boundary conditions for the velocity), and the
 * sparsity patterns and matrices for each level. The resize()
 * function of MGLevelObject<T> will destroy all existing contained
 * objects.
 * 
 * @code
 *         std::set<types::boundary_id> zero_boundary_ids;
 *         zero_boundary_ids.insert(0);
 * 
 *         mg_constrained_dofs.clear();
 *         mg_constrained_dofs.initialize(velocity_dof_handler);
 *         mg_constrained_dofs.make_zero_boundary_constraints(velocity_dof_handler,
 *                                                            zero_boundary_ids);
 *         const unsigned int n_levels = triangulation.n_levels();
 * 
 *         mg_interface_matrices.resize(0, n_levels - 1);
 *         mg_matrices.resize(0, n_levels - 1);
 *         mg_sparsity_patterns.resize(0, n_levels - 1);
 * 
 *         for (unsigned int level = 0; level < n_levels; ++level)
 *           {
 *             DynamicSparsityPattern csp(velocity_dof_handler.n_dofs(level),
 *                                        velocity_dof_handler.n_dofs(level));
 *             MGTools::make_sparsity_pattern(velocity_dof_handler, csp, level);
 *             mg_sparsity_patterns[level].copy_from(csp);
 * 
 *             mg_matrices[level].reinit(mg_sparsity_patterns[level]);
 *             mg_interface_matrices[level].reinit(mg_sparsity_patterns[level]);
 *           }
 *       }
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
 *     const unsigned int n_u = dofs_per_block[0];
 *     const unsigned int n_p = dofs_per_block[1];
 * 
 *     {
 *       constraints.clear();
 * @endcode
 * 
 * The following makes use of a component mask for interpolation of the
 * boundary values for the velocity only, which is further explained in
 * the vector valued dealii step-20 tutorial.
 * 
 * @code
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                Solution<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 * 
 * @endcode
 * 
 * As discussed in the introduction, we need to fix one degree of freedom
 * of the pressure variable to ensure solvability of the problem. We do
 * this here by marking the first pressure dof, which has index n_u as a
 * constrained dof.
 * 
 * @code
 *       if (solver_type == SolverType::UMFPACK)
 *         constraints.add_line(n_u);
 * 
 *       constraints.close();
 *     }
 * 
 *     std::cout << "\tNumber of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "\tNumber of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (" << n_u << '+' << n_p << ')' << std::endl;
 * 
 *     {
 *       BlockDynamicSparsityPattern csp(dofs_per_block, dofs_per_block);
 *       DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
 *       sparsity_pattern.copy_from(csp);
 *     }
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     solution.reinit(dofs_per_block);
 *     system_rhs.reinit(dofs_per_block);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemassemble_system"></a> 
 * <h4>StokesProblem::assemble_system</h4>
 * 

 * 
 * In this function, the system matrix is assembled. We assemble the pressure
 * mass matrix in the (1,1) block (if needed) and move it out of this location
 * at the end of this function.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::assemble_system()
 *   {
 *     TimerOutput::Scope assemble(computing_timer, "Assemble");
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 * @endcode
 * 
 * If true, we will assemble the pressure mass matrix in the (1,1) block:
 * 
 * @code
 *     const bool assemble_pressure_mass_matrix =
 *       (solver_type == SolverType::UMFPACK) ? false : true;
 * 
 *     QGauss<dim> quadrature_formula(pressure_degree + 2);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values | update_gradients);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const RightHandSide<dim>    right_hand_side;
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
 *     std::vector<double>                  div_phi_u(dofs_per_cell);
 *     std::vector<double>                  phi_p(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 *         right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
 *                                           rhs_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 symgrad_phi_u[k] =
 *                   fe_values[velocities].symmetric_gradient(k, q);
 *                 div_phi_u[k] = fe_values[velocities].divergence(k, q);
 *                 phi_p[k]     = fe_values[pressure].value(k, q);
 *               }
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 for (unsigned int j = 0; j <= i; ++j)
 *                   {
 *                     local_matrix(i, j) +=
 *                       (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) -
 *                        div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
 *                        (assemble_pressure_mass_matrix ? phi_p[i] * phi_p[j] :
 *                                                         0)) *
 *                       fe_values.JxW(q);
 *                   }
 * 
 *                 const unsigned int component_i =
 *                   fe.system_to_component_index(i).first;
 *                 local_rhs(i) += fe_values.shape_value(i, q) *
 *                                 rhs_values[q](component_i) * fe_values.JxW(q);
 *               }
 *           }
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
 *             local_matrix(i, j) = local_matrix(j, i);
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(local_matrix,
 *                                                local_rhs,
 *                                                local_dof_indices,
 *                                                system_matrix,
 *                                                system_rhs);
 *       }
 * 
 *     if (solver_type != SolverType::UMFPACK)
 *       {
 *         pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
 *         pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
 *         system_matrix.block(1, 1) = 0;
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemassemble_multigrid"></a> 
 * <h4>StokesProblem::assemble_multigrid</h4>
 * 

 * 
 * Here, like in step-16, we have a function that assembles the level
 * and interface matrices necessary for the multigrid preconditioner.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::assemble_multigrid()
 *   {
 *     TimerOutput::Scope multigrid_specific(computing_timer,
 *                                           "(Multigrid specific)");
 *     TimerOutput::Scope assemble_multigrid(computing_timer,
 *                                           "Assemble Multigrid");
 * 
 *     mg_matrices = 0.;
 * 
 *     QGauss<dim> quadrature_formula(pressure_degree + 2);
 * 
 *     FEValues<dim> fe_values(velocity_fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values | update_gradients);
 * 
 *     const unsigned int dofs_per_cell = velocity_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
 * 
 *     std::vector<AffineConstraints<double>> boundary_constraints(
 *       triangulation.n_levels());
 *     std::vector<AffineConstraints<double>> boundary_interface_constraints(
 *       triangulation.n_levels());
 *     for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
 *       {
 *         boundary_constraints[level].add_lines(
 *           mg_constrained_dofs.get_refinement_edge_indices(level));
 *         boundary_constraints[level].add_lines(
 *           mg_constrained_dofs.get_boundary_indices(level));
 *         boundary_constraints[level].close();
 * 
 *         IndexSet idx = mg_constrained_dofs.get_refinement_edge_indices(level) &
 *                        mg_constrained_dofs.get_boundary_indices(level);
 * 
 *         boundary_interface_constraints[level].add_lines(idx);
 *         boundary_interface_constraints[level].close();
 *       }
 * 
 * @endcode
 * 
 * This iterator goes over all cells (not just active)
 * 
 * @code
 *     for (const auto &cell : velocity_dof_handler.cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         cell_matrix = 0;
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient(k, q);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int j = 0; j <= i; ++j)
 *                 {
 *                   cell_matrix(i, j) +=
 *                     (symgrad_phi_u[i] * symgrad_phi_u[j]) * fe_values.JxW(q);
 *                 }
 *           }
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
 *             cell_matrix(i, j) = cell_matrix(j, i);
 * 
 *         cell->get_mg_dof_indices(local_dof_indices);
 * 
 *         boundary_constraints[cell->level()].distribute_local_to_global(
 *           cell_matrix, local_dof_indices, mg_matrices[cell->level()]);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             if (!mg_constrained_dofs.at_refinement_edge(cell->level(),
 *                                                         local_dof_indices[i]) ||
 *                 mg_constrained_dofs.at_refinement_edge(cell->level(),
 *                                                        local_dof_indices[j]))
 *               cell_matrix(i, j) = 0;
 * 
 *         boundary_interface_constraints[cell->level()]
 *           .distribute_local_to_global(cell_matrix,
 *                                       local_dof_indices,
 *                                       mg_interface_matrices[cell->level()]);
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsolve"></a> 
 * <h4>StokesProblem::solve</h4>
 * 

 * 
 * This function sets up things differently based on if you want to use ILU
 * or GMG as a preconditioner.  Both methods share the same solver (FGMRES)
 * but require a different preconditioner to be initialized. Here we time not
 * only the entire solve function, but we separately time the setup of the
 * preconditioner as well as the solve itself.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::solve()
 *   {
 *     TimerOutput::Scope solve(computing_timer, "Solve");
 *     constraints.set_zero(solution);
 * 
 *     if (solver_type == SolverType::UMFPACK)
 *       {
 *         computing_timer.enter_subsection("(UMFPACK specific)");
 *         computing_timer.enter_subsection("Solve - Initialize");
 * 
 *         SparseDirectUMFPACK A_direct;
 *         A_direct.initialize(system_matrix);
 * 
 *         computing_timer.leave_subsection();
 *         computing_timer.leave_subsection();
 * 
 *         {
 *           TimerOutput::Scope solve_backslash(computing_timer,
 *                                              "Solve - Backslash");
 *           A_direct.vmult(solution, system_rhs);
 *         }
 * 
 *         constraints.distribute(solution);
 *         return;
 *       }
 * 
 * @endcode
 * 
 * Here we must make sure to solve for the residual with "good enough"
 * accuracy
 * 
 * @code
 *     SolverControl solver_control(system_matrix.m(),
 *                                  1e-10 * system_rhs.l2_norm());
 *     unsigned int  n_iterations_A;
 *     unsigned int  n_iterations_S;
 * 
 * @endcode
 * 
 * This is used to pass whether or not we want to solve for A inside
 * the preconditioner.  One could change this to false to see if
 * there is still convergence and if so does the program then run
 * faster or slower
 * 
 * @code
 *     const bool use_expensive = true;
 * 
 *     SolverFGMRES<BlockVector<double>> solver(solver_control);
 * 
 *     if (solver_type == SolverType::FGMRES_ILU)
 *       {
 *         computing_timer.enter_subsection("(ILU specific)");
 *         computing_timer.enter_subsection("Solve - Set-up Preconditioner");
 * 
 *         std::cout << "   Computing preconditioner..." << std::endl
 *                   << std::flush;
 * 
 *         SparseILU<double> A_preconditioner;
 *         A_preconditioner.initialize(system_matrix.block(0, 0));
 * 
 *         SparseILU<double> S_preconditioner;
 *         S_preconditioner.initialize(pressure_mass_matrix);
 * 
 *         const BlockSchurPreconditioner<SparseILU<double>, SparseILU<double>>
 *           preconditioner(system_matrix,
 *                          pressure_mass_matrix,
 *                          A_preconditioner,
 *                          S_preconditioner,
 *                          use_expensive);
 * 
 *         computing_timer.leave_subsection();
 *         computing_timer.leave_subsection();
 * 
 *         {
 *           TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES");
 * 
 *           solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *           n_iterations_A = preconditioner.n_iterations_A;
 *           n_iterations_S = preconditioner.n_iterations_S;
 *         }
 *       }
 *     else
 *       {
 *         computing_timer.enter_subsection("(Multigrid specific)");
 *         computing_timer.enter_subsection("Solve - Set-up Preconditioner");
 * 
 * @endcode
 * 
 * Transfer operators between levels
 * 
 * @code
 *         MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs);
 *         mg_transfer.build(velocity_dof_handler);
 * 
 * @endcode
 * 
 * Setup coarse grid solver
 * 
 * @code
 *         FullMatrix<double> coarse_matrix;
 *         coarse_matrix.copy_from(mg_matrices[0]);
 *         MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
 *         coarse_grid_solver.initialize(coarse_matrix);
 * 
 *         using Smoother = PreconditionSOR<SparseMatrix<double>>;
 *         mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother;
 *         mg_smoother.initialize(mg_matrices);
 *         mg_smoother.set_steps(2);
 * 
 * @endcode
 * 
 * Multigrid, when used as a preconditioner for CG, needs to be a
 * symmetric operator, so the smoother must be symmetric
 * 
 * @code
 *         mg_smoother.set_symmetric(true);
 * 
 *         mg::Matrix<Vector<double>> mg_matrix(mg_matrices);
 *         mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices);
 *         mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices);
 * 
 * @endcode
 * 
 * Now, we are ready to set up the V-cycle operator and the multilevel
 * preconditioner.
 * 
 * @code
 *         Multigrid<Vector<double>> mg(
 *           mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother);
 *         mg.set_edge_matrices(mg_interface_down, mg_interface_up);
 * 
 *         PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
 *           A_Multigrid(velocity_dof_handler, mg, mg_transfer);
 * 
 *         SparseILU<double> S_preconditioner;
 *         S_preconditioner.initialize(pressure_mass_matrix,
 *                                     SparseILU<double>::AdditionalData());
 * 
 *         const BlockSchurPreconditioner<
 *           PreconditionMG<dim,
 *                          Vector<double>,
 *                          MGTransferPrebuilt<Vector<double>>>,
 *           SparseILU<double>>
 *           preconditioner(system_matrix,
 *                          pressure_mass_matrix,
 *                          A_Multigrid,
 *                          S_preconditioner,
 *                          use_expensive);
 * 
 *         computing_timer.leave_subsection();
 *         computing_timer.leave_subsection();
 * 
 *         {
 *           TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES");
 *           solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *           n_iterations_A = preconditioner.n_iterations_A;
 *           n_iterations_S = preconditioner.n_iterations_S;
 *         }
 *       }
 * 
 *     constraints.distribute(solution);
 * 
 *     std::cout
 *       << std::endl
 *       << "\tNumber of FGMRES iterations: " << solver_control.last_step()
 *       << std::endl
 *       << "\tTotal number of iterations used for approximation of A inverse: "
 *       << n_iterations_A << std::endl
 *       << "\tTotal number of iterations used for approximation of S inverse: "
 *       << n_iterations_S << std::endl
 *       << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemprocess_solution"></a> 
 * <h4>StokesProblem::process_solution</h4>
 * 

 * 
 * This function computes the L2 and H1 errors of the solution. For this,
 * we need to make sure the pressure has mean zero.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::compute_errors()
 *   {
 * @endcode
 * 
 * Compute the mean pressure $\frac{1}{\Omega} \int_{\Omega} p(x) dx $
 * and then subtract it from each pressure coefficient. This will result
 * in a pressure with mean value zero. Here we make use of the fact that
 * the pressure is component $dim$ and that the finite element space
 * is nodal.
 * 
 * @code
 *     const double mean_pressure = VectorTools::compute_mean_value(
 *       dof_handler, QGauss<dim>(pressure_degree + 2), solution, dim);
 *     solution.block(1).add(-mean_pressure);
 *     std::cout << "   Note: The mean value was adjusted by " << -mean_pressure
 *               << std::endl;
 * 
 *     const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
 *     const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
 *                                                      dim + 1);
 * 
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(pressure_degree + 2),
 *                                       VectorTools::L2_norm,
 *                                       &velocity_mask);
 * 
 *     const double Velocity_L2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(pressure_degree + 2),
 *                                       VectorTools::L2_norm,
 *                                       &pressure_mask);
 * 
 *     const double Pressure_L2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(pressure_degree + 2),
 *                                       VectorTools::H1_norm,
 *                                       &velocity_mask);
 * 
 *     const double Velocity_H1_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::H1_norm);
 * 
 *     std::cout << std::endl
 *               << "   Velocity L2 Error: " << Velocity_L2_error << std::endl
 *               << "   Pressure L2 Error: " << Pressure_L2_error << std::endl
 *               << "   Velocity H1 Error: " << Velocity_H1_error << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemoutput_results"></a> 
 * <h4>StokesProblem::output_results</h4>
 * 

 * 
 * This function generates graphical output like it is done in step-22.
 * 
 * @code
 *   template <int dim>
 *   void
 *   StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const
 *   {
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     solution_names.emplace_back("pressure");
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution,
 *                              solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 *     data_out.build_patches();
 * 
 *     std::ofstream output(
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemrun"></a> 
 * <h4>StokesProblem::run</h4>
 * 

 * 
 * The last step in the Stokes class is, as usual, the function that
 * generates the initial grid and calls the other functions in the
 * respective order.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_cube(triangulation);
 *     triangulation.refine_global(6 - dim);
 * 
 *     if (solver_type == SolverType::FGMRES_ILU)
 *       std::cout << "Now running with ILU" << std::endl;
 *     else if (solver_type == SolverType::FGMRES_GMG)
 *       std::cout << "Now running with Multigrid" << std::endl;
 *     else
 *       std::cout << "Now running with UMFPACK" << std::endl;
 * 
 * 
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 3;
 *          ++refinement_cycle)
 *       {
 *         std::cout << "Refinement cycle " << refinement_cycle << std::endl;
 * 
 *         if (refinement_cycle > 0)
 *           triangulation.refine_global(1);
 * 
 *         std::cout << "   Set-up..." << std::endl;
 *         setup_dofs();
 * 
 *         std::cout << "   Assembling..." << std::endl;
 *         assemble_system();
 * 
 *         if (solver_type == SolverType::FGMRES_GMG)
 *           {
 *             std::cout << "   Assembling Multigrid..." << std::endl;
 * 
 *             assemble_multigrid();
 *           }
 * 
 *         std::cout << "   Solving..." << std::flush;
 *         solve();
 * 
 *         compute_errors();
 * 
 *         output_results(refinement_cycle);
 * 
 *         Utilities::System::MemoryStats mem;
 *         Utilities::System::get_memory_stats(mem);
 *         std::cout << "   VM Peak: " << mem.VmPeak << std::endl;
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 *       }
 *   }
 * } // namespace Step56
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step56;
 * 
 *       const int degree = 1;
 *       const int dim    = 3;
 * @endcode
 * 
 * options for SolverType: UMFPACK FGMRES_ILU FGMRES_GMG
 * 
 * @code
 *       StokesProblem<dim> flow_problem(degree, SolverType::FGMRES_GMG);
 * 
 *       flow_problem.run();
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


<a name="Errors"></a><h3> Errors </h3>


We first run the code and confirm that the finite element solution converges
with the correct rates as predicted by the error analysis of mixed finite
element problems. Given sufficiently smooth exact solutions $u$ and $p$,
the errors of the Taylor-Hood element $Q_k \times Q_{k-1}$ should be

@f[
\| u -u_h \|_0 + h ( \| u- u_h\|_1 + \|p - p_h \|_0)
\leq C h^{k+1} ( \|u \|_{k+1} + \| p \|_k )
@f]

see for example Ern/Guermond "Theory and Practice of Finite Elements", Section
4.2.5 p195. This is indeed what we observe, using the $Q_2 \times Q_1$
element as an example (this is what is done in the code, but is easily
changed in <code>main()</code>):

<table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th>L2 Velocity</th>
    <th>Reduction</th>
    <th>L2 Pressure</th>
    <th>Reduction</th>
    <th>H1 Velocity</th>
    <th>Reduction</th>
  </tr>
  <tr>
    <td>3D, 3 global refinements</td>
    <td>0.000670888</td>
    <td align="center">-</td>
    <td>0.0036533</td>
    <td align="center">-</td>
    <td>0.0414704</td>
    <td align="center">-</td>
  </tr>
  <tr>
    <td>3D, 4 global refinements</td>
    <td>8.38E-005</td>
    <td>8.0</td>
    <td>0.00088494</td>
    <td>4.1</td>
    <td>0.0103781</td>
    <td>4.0</td>
  </tr>
  <tr>
    <td>3D, 5 global refinements</td>
    <td>1.05E-005</td>
    <td>8.0</td>
    <td>0.000220253</td>
    <td>4.0</td>
    <td>0.00259519</td>
    <td>4.0</td>
</th>
  </tr>
</table>

<a name="TimingResults"></a><h3> Timing Results </h3>


Let us compare the direct solver approach using UMFPACK to the two
methods in which we choose $\widetilde {A^{-1}}=A^{-1}$ and
$\widetilde{S^{-1}}=S^{-1}$ by solving linear systems with $A,S$ using
CG. The preconditioner for CG is then either ILU or GMG.
The following table summarizes solver iterations, timings, and virtual
memory (VM) peak usage:

<table align="center" class="doxtable">
<tr>
  <th></th>
  <th colspan="3">General</th>
  <th colspan="6">GMG</th>
  <th colspan="6">ILU</th>
  <th colspan="3">UMFPACK</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th colspan="2">Timings</th>
  <th colspan="2">Timings</th>
  <th colspan="3">Iterations</th>
  <th></th>
  <th colspan="2">Timings</th>
  <th colspan="3">Iterations</th>
  <th></th>
  <th colspan="2">Timings</th>
  <th></th>
</tr>
<tr>
  <th>Cycle</th>
  <th>DoFs</th>
  <th>Setup</th>
  <th>Assembly</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>Outer</th>
  <th>Inner (A)</th>
  <th>Inner (S)</th>
  <th>VM Peak</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>Outer</th>
  <th>Inner (A)</th>
  <th>Inner (S)</th>
  <th>VM Peak</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>VM Peak</th>
</tr>
<tr>
  <td>0</td>
  <td>15468</td>
  <td>0.1s</td>
  <td>0.3s</td>
  <td>0.3s</td>
  <td>1.3s</td>
  <td>21</td>
  <td>67</td>
  <td>22</td>
  <td>4805</td>
  <td>0.3s</td>
  <td>0.6s</td>
  <td>21</td>
  <td>180</td>
  <td>22</td>
  <td>4783</td>
  <td>2.65s</td>
  <td>2.8s</td>
  <td>5054</td>
</tr>
<tr>
  <td>1</td>
  <td>112724</td>
  <td>1.0s</td>
  <td>2.4s</td>
  <td>2.6s</td>
  <td>14s</td>
  <td>21</td>
  <td>67</td>
  <td>22</td>
  <td>5441</td>
  <td>2.8s</td>
  <td>15.8s</td>
  <td>21</td>
  <td>320</td>
  <td>22</td>
  <td>5125</td>
  <td>236s</td>
  <td>237s</td>
  <td>11288</td>
</tr>
<tr>
  <td>2</td>
  <td>859812</td>
  <td>9.0s</td>
  <td>20s</td>
  <td>20s</td>
  <td>101s</td>
  <td>20</td>
  <td>65</td>
  <td>21</td>
  <td>10641</td>
  <td>27s</td>
  <td>268s</td>
  <td>21</td>
  <td>592</td>
  <td>22</td>
  <td>8307</td>
  <td>-</td>
  <td>-</td>
  <td>-</td>
</tr>
</table>

As can be seen from the table:

1. UMFPACK uses large amounts of memory, especially in 3d. Also, UMFPACK
timings do not scale favorably with problem size.

2. Because we are using inner solvers for $A$ and $S$, ILU and GMG require the
same number of outer iterations.

3. The number of (inner) iterations for $A$ increases for ILU with refinement, leading
to worse than linear scaling in solve time. In contrast, the number of inner
iterations for $A$ stays constant with GMG leading to nearly perfect scaling in
solve time.

4. GMG needs slightly more memory than ILU to store the level and interface
matrices.

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="Checkhigherorderdiscretizations"></a><h4> Check higher order discretizations </h4>


Experiment with higher order stable FE pairs and check that you observe the
correct convergence rates.

<a name="Comparewithcheappreconditioner"></a><h4> Compare with cheap preconditioner </h4>


The introduction also outlined another option to precondition the
overall system, namely one in which we do not choose $\widetilde
{A^{-1}}=A^{-1}$ as in the table above, but in which
$\widetilde{A^{-1}}$ is only a single preconditioner application with
GMG or ILU, respectively.

This is in fact implemented in the code: Currently, the boolean
<code>use_expensive</code> in <code>solve()</code> is set to @p true. The
option mentioned above is obtained by setting it to @p false.

What you will find is that the number of FGMRES iterations stays
constant under refinement if you use GMG this way. This means that the
Multigrid is optimal and independent of $h$.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-56.cc"
*/
