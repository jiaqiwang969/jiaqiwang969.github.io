/**
@page step_55 The step-55 tutorial program
This tutorial depends on step-40, step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Optimalpreconditioners">Optimal preconditioners</a>
        <li><a href="#Thesolverandpreconditioner">The solver and preconditioner</a>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
        <li><a href="#Problemsetup">Problem setup</a>
        <li><a href="#Themainprogram">The main program</a>
        <li><a href="#SystemSetup">System Setup</a>
        <li><a href="#Assembly">Assembly</a>
        <li><a href="#Solving">Solving</a>
        <li><a href="#Therest">The rest</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#InvestigateTrilinositerations">Investigate Trilinos iterations</a>
        <li><a href="#SolvetheOseenprobleminsteadoftheStokessystem">Solve the Oseen problem instead of the Stokes system</a>
        <li><a href="#Adaptiverefinement">Adaptive refinement</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Timo Heister. Special thanks to Sander
Rhebergen for the inspiration to finally write this tutorial.

This material is based upon work partially supported by National Science
Foundation grant DMS1522191 and the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-0949446 and The University of California-Davis.

The authors would like to thank the Isaac Newton Institute for
Mathematical Sciences, Cambridge, for support and hospitality during
the programme Melt in the Mantle where work on this tutorial was
undertaken. This work was supported by EPSRC grant no EP/K032208/1.
</i>


@note As a prerequisite of this program, you need to have PETSc or Trilinos
and the p4est library installed. The installation of deal.II together with
these additional libraries is described in the <a href="../../readme.html"
target="body">README</a> file.

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


Building on step-40, this tutorial shows how to solve linear PDEs with several
components in parallel using MPI with PETSc or Trilinos for the linear
algebra. For this, we return to the Stokes equations as discussed in
step-22. The motivation for writing this tutorial is to provide an
intermediate step (pun intended) between step-40 (parallel Laplace) and
step-32 (parallel coupled Stokes with Boussinesq for a time dependent
problem).

The learning outcomes for this tutorial are:

- You are able to solve PDEs with several variables in parallel and can
  apply this to different problems.

- You understand the concept of optimal preconditioners and are able to check
  this for a particular problem.

- You are able to construct manufactured solutions using the free computer
  algreba system SymPy (https://sympy.org).

- You can implement various other tasks for parallel programs: error
  computation, writing graphical output, etc.

- You can visualize vector fields, stream lines, and contours of vector
  quantities.

We are solving for a velocity $\textbf{u}$ and pressure $p$ that satisfy the
Stokes equation, which reads
@f{eqnarray*}
  - \triangle \textbf{u} + \nabla p &=& \textbf{f}, \\
  -\textrm{div}\; \textbf{u} &=& 0.
@f}


<a name="Optimalpreconditioners"></a><h3>Optimal preconditioners</h3>


Make sure that you read (even better: try) what is described in "Block Schur
complement preconditioner" in the "Possible Extensions" section in step-22.
Like described there, we are going to solve the block system using a Krylov
method and a block preconditioner.

Our goal here is to construct a very simple (maybe the simplest?) optimal
preconditioner for the linear system. A preconditioner is called "optimal" or
"of optimal complexity", if the number of iterations of the preconditioned
system is independent of the mesh size $h$. You can extend that definition to
also require indepence of the number of processors used (we will discuss that
in the results section), the computational domain and the mesh quality, the
test case itself, the polynomial degree of the finite element space, and more.

Why is a constant number of iterations considered to be "optimal"? Assume the
discretized PDE gives a linear system with N unknowns. Because the matrix
coming from the FEM discretization is sparse, a matrix-vector product can be
done in O(N) time. A preconditioner application can also only be O(N) at best
(for example doable with multigrid methods). If the number of iterations
required to solve the linear system is independent of $h$ (and therefore N),
the total cost of solving the system will be O(N). It is not possible to beat
this complexity, because even looking at all the entries of the right-hand
side already takes O(N) time. For more information see @cite elman2005,
Chapter 2.5 (Multigrid).

The preconditioner described here is even simpler than the one described in
step-22 and will typically require more iterations and consequently time to
solve. When considering preconditioners, optimality is not the only important
metric. But an optimal and expensive preconditioner is typically more
desirable than a cheaper, non-optimal one. This is because, eventually, as the
mesh size becomes smaller and smaller and linear problems become bigger and
bigger, the former will eventually beat the latter.

<a name="Thesolverandpreconditioner"></a><h3>The solver and preconditioner</h3>


We precondition the linear system
@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{c}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{c}
    F \\ 0
  \end{array}\right),
@f}

with the block diagonal preconditioner
@f{eqnarray*}
  P^{-1}
  =
  \left(\begin{array}{cc}
    A & 0 \\ 0 & S
  \end{array}\right) ^{-1},
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ 0 & S^{-1}
  \end{array}\right),
@f}
where $S=-BA^{-1} B^T$ is the Schur complement.

With this choice of $P$, assuming that we handle $A^{-1}$ and $S^{-1}$ exactly
(which is an "idealized" situation), the preconditioned linear system has
three distinct eigenvalues independent of $h$ and is therefore "optimal".  See
section 6.2.1 (especially p. 292) in @cite elman2005. For comparison,
using the ideal version of the upper block-triangular preconditioner in
step-22 (also used in step-56) would have all eigenvalues be equal to one.

We will use approximations of the inverse operations in $P^{-1}$ that are
(nearly) independent of $h$. In this situation, one can again show, that the
eigenvalues are independent of $h$. For the Krylov method we choose MINRES,
which is attractive for the analysis (iteration count is proven to be
independent of $h$, see the remainder of the chapter 6.2.1 in the book
mentioned above), great from the computational standpoint (simpler and cheaper
than GMRES for example), and applicable (matrix and preconditioner are
symmetric).

For the approximations we will use a CG solve with the mass matrix in the
pressure space for approximating the action of $S^{-1}$. Note that the mass
matrix is spectrally equivalent to $S$. We can expect the number of CG
iterations to be independent of $h$, even with a simple preconditioner like
ILU.

For the approximation of the velocity block $A$ we will perform a single AMG
V-cycle. In practice this choice is not exactly independent of $h$, which can
explain the slight increase in iteration numbers. A possible explanation is
that the coarsest level will be solved exactly and the number of levels and
size of the coarsest matrix is not predictable.


<a name="Thetestcase"></a><h3>The testcase</h3>


We will construct a manufactured solution based on the classical Kovasznay problem,
see @cite kovasznay1948laminar. Here
is an image of the solution colored by the x velocity including
streamlines of the velocity:

 <img src="https://www.dealii.org/images/steps/developer/step-55.solution.png" alt="">

We have to cheat here, though, because we are not solving the non-linear
Navier-Stokes equations, but the linear Stokes system without convective
term. Therefore, to recreate the exact same solution, we use the method of
manufactured solutions with the solution of the Kovasznay problem. This will
effectively move the convective term into the right-hand side $f$.

The right-hand side is computed using the script "reference.py" and we use
the exact solution for boundary conditions and error computation.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * 
 * @endcode
 * 
 * The following chunk out code is identical to step-40 and allows
 * switching between PETSc and Trilinos:
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/generic_linear_algebra.h>
 * 
 * /* #define FORCE_USE_OF_TRILINOS */
 * 
 * namespace LA
 * {
 * #if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
 *   !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
 *   using namespace dealii::LinearAlgebraPETSc;
 * #  define USE_PETSC_LA
 * #elif defined(DEAL_II_WITH_TRILINOS)
 *   using namespace dealii::LinearAlgebraTrilinos;
 * #else
 * #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
 * #endif
 * } // namespace LA
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/solver_minres.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * #include <deal.II/lac/petsc_sparse_matrix.h>
 * #include <deal.II/lac/petsc_vector.h>
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/petsc_precondition.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/lac/sparsity_tools.h>
 * #include <deal.II/distributed/tria.h>
 * #include <deal.II/distributed/grid_refinement.h>
 * 
 * #include <cmath>
 * #include <fstream>
 * #include <iostream>
 * 
 * namespace Step55
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * We need a few helper classes to represent our solver strategy
 * described in the introduction.
 * 

 * 
 * 
 * @code
 *   namespace LinearSolvers
 *   {
 * @endcode
 * 
 * This class exposes the action of applying the inverse of a giving
 * matrix via the function InverseMatrix::vmult(). Internally, the
 * inverse is not formed explicitly. Instead, a linear solver with CG
 * is performed. This class extends the InverseMatrix class in step-22
 * with an option to specify a preconditioner, and to allow for different
 * vector types in the vmult function.
 * 
 * @code
 *     template <class Matrix, class Preconditioner>
 *     class InverseMatrix : public Subscriptor
 *     {
 *     public:
 *       InverseMatrix(const Matrix &m, const Preconditioner &preconditioner);
 * 
 *       template <typename VectorType>
 *       void vmult(VectorType &dst, const VectorType &src) const;
 * 
 *     private:
 *       const SmartPointer<const Matrix> matrix;
 *       const Preconditioner &           preconditioner;
 *     };
 * 
 * 
 *     template <class Matrix, class Preconditioner>
 *     InverseMatrix<Matrix, Preconditioner>::InverseMatrix(
 *       const Matrix &        m,
 *       const Preconditioner &preconditioner)
 *       : matrix(&m)
 *       , preconditioner(preconditioner)
 *     {}
 * 
 * 
 * 
 *     template <class Matrix, class Preconditioner>
 *     template <typename VectorType>
 *     void
 *     InverseMatrix<Matrix, Preconditioner>::vmult(VectorType &      dst,
 *                                                  const VectorType &src) const
 *     {
 *       SolverControl solver_control(src.size(), 1e-8 * src.l2_norm());
 *       SolverCG<LA::MPI::Vector> cg(solver_control);
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
 * 
 * @endcode
 * 
 * The class A template class for a simple block diagonal preconditioner
 * for 2x2 matrices.
 * 
 * @code
 *     template <class PreconditionerA, class PreconditionerS>
 *     class BlockDiagonalPreconditioner : public Subscriptor
 *     {
 *     public:
 *       BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A,
 *                                   const PreconditionerS &preconditioner_S);
 * 
 *       void vmult(LA::MPI::BlockVector &      dst,
 *                  const LA::MPI::BlockVector &src) const;
 * 
 *     private:
 *       const PreconditionerA &preconditioner_A;
 *       const PreconditionerS &preconditioner_S;
 *     };
 * 
 *     template <class PreconditionerA, class PreconditionerS>
 *     BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::
 *       BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A,
 *                                   const PreconditionerS &preconditioner_S)
 *       : preconditioner_A(preconditioner_A)
 *       , preconditioner_S(preconditioner_S)
 *     {}
 * 
 * 
 *     template <class PreconditionerA, class PreconditionerS>
 *     void BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::vmult(
 *       LA::MPI::BlockVector &      dst,
 *       const LA::MPI::BlockVector &src) const
 *     {
 *       preconditioner_A.vmult(dst.block(0), src.block(0));
 *       preconditioner_S.vmult(dst.block(1), src.block(1));
 *     }
 * 
 *   } // namespace LinearSolvers
 * 
 * @endcode
 * 
 * 
 * <a name="Problemsetup"></a> 
 * <h3>Problem setup</h3>
 * 

 * 
 * The following classes represent the right hand side and the exact
 * solution for the test problem.
 * 

 * 
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
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  value) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   void RightHandSide<dim>::vector_value(const Point<dim> &p,
 *                                         Vector<double> &  values) const
 *   {
 *     const double R_x = p[0];
 *     const double R_y = p[1];
 * 
 *     const double pi  = numbers::PI;
 *     const double pi2 = pi * pi;
 *     values[0] =
 *       -1.0L / 2.0L * (-2 * sqrt(25.0 + 4 * pi2) + 10.0) *
 *         exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
 *       0.4 * pi2 * exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) +
 *       0.1 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 2) *
 *         exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi);
 *     values[1] = 0.2 * pi * (-sqrt(25.0 + 4 * pi2) + 5.0) *
 *                   exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) -
 *                 0.05 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 3) *
 *                   exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
 *                   pi;
 *     values[2] = 0;
 *   }
 * 
 * 
 *   template <int dim>
 *   class ExactSolution : public Function<dim>
 *   {
 *   public:
 *     ExactSolution()
 *       : Function<dim>(dim + 1)
 *     {}
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  value) const override;
 *   };
 * 
 *   template <int dim>
 *   void ExactSolution<dim>::vector_value(const Point<dim> &p,
 *                                         Vector<double> &  values) const
 *   {
 *     const double R_x = p[0];
 *     const double R_y = p[1];
 * 
 *     const double pi  = numbers::PI;
 *     const double pi2 = pi * pi;
 *     values[0] =
 *       -exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 1;
 *     values[1] = (1.0L / 2.0L) * (-sqrt(25.0 + 4 * pi2) + 5.0) *
 *                 exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
 *                 pi;
 *     values[2] =
 *       -1.0L / 2.0L * exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
 *       2.0 *
 *         (-6538034.74494422 +
 *          0.0134758939981709 * exp(4 * sqrt(25.0 + 4 * pi2))) /
 *         (-80.0 * exp(3 * sqrt(25.0 + 4 * pi2)) +
 *          16.0 * sqrt(25.0 + 4 * pi2) * exp(3 * sqrt(25.0 + 4 * pi2))) -
 *       1634508.68623606 * exp(-3.0 * sqrt(25.0 + 4 * pi2)) /
 *         (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)) +
 *       (-0.00673794699908547 * exp(sqrt(25.0 + 4 * pi2)) +
 *        3269017.37247211 * exp(-3 * sqrt(25.0 + 4 * pi2))) /
 *         (-8 * sqrt(25.0 + 4 * pi2) + 40.0) +
 *       0.00336897349954273 * exp(1.0 * sqrt(25.0 + 4 * pi2)) /
 *         (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainprogram"></a> 
 * <h3>The main program</h3>
 *   

 * 
 * The main class is very similar to step-40, except that matrices and
 * vectors are now block versions, and we store a std::vector<IndexSet>
 * for owned and relevant DoFs instead of a single IndexSet. We have
 * exactly two IndexSets, one for all velocity unknowns and one for all
 * pressure unknowns.
 * 
 * @code
 *   template <int dim>
 *   class StokesProblem
 *   {
 *   public:
 *     StokesProblem(unsigned int velocity_degree);
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     unsigned int velocity_degree;
 *     double       viscosity;
 *     MPI_Comm     mpi_communicator;
 * 
 *     FESystem<dim>                             fe;
 *     parallel::distributed::Triangulation<dim> triangulation;
 *     DoFHandler<dim>                           dof_handler;
 * 
 *     std::vector<IndexSet> owned_partitioning;
 *     std::vector<IndexSet> relevant_partitioning;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     LA::MPI::BlockSparseMatrix system_matrix;
 *     LA::MPI::BlockSparseMatrix preconditioner_matrix;
 *     LA::MPI::BlockVector       locally_relevant_solution;
 *     LA::MPI::BlockVector       system_rhs;
 * 
 *     ConditionalOStream pcout;
 *     TimerOutput        computing_timer;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   StokesProblem<dim>::StokesProblem(unsigned int velocity_degree)
 *     : velocity_degree(velocity_degree)
 *     , viscosity(0.1)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree - 1), 1)
 *     , triangulation(mpi_communicator,
 *                     typename Triangulation<dim>::MeshSmoothing(
 *                       Triangulation<dim>::smoothing_on_refinement |
 *                       Triangulation<dim>::smoothing_on_coarsening))
 *     , dof_handler(triangulation)
 *     , pcout(std::cout,
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *     , computing_timer(mpi_communicator,
 *                       pcout,
 *                       TimerOutput::summary,
 *                       TimerOutput::wall_times)
 *   {}
 * 
 * 
 * @endcode
 * 
 * The Kovasnay flow is defined on the domain [-0.5, 1.5]^2, which we
 * create by passing the min and max values to GridGenerator::hyper_cube.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -0.5, 1.5);
 *     triangulation.refine_global(3);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="SystemSetup"></a> 
 * <h3>System Setup</h3>
 *   

 * 
 * The construction of the block matrices and vectors is new compared to
 * step-40 and is different compared to serial codes like step-22, because
 * we need to supply the set of rows that belong to our processor.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::setup_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "setup");
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * Put all dim velocities into block 0 and the pressure into block 1,
 * then reorder the unknowns by block. Finally count how many unknowns
 * we have per block.
 * 
 * @code
 *     std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
 *     stokes_sub_blocks[dim] = 1;
 *     DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks);
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks);
 * 
 *     const unsigned int n_u = dofs_per_block[0];
 *     const unsigned int n_p = dofs_per_block[1];
 * 
 *     pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
 *           << n_u << '+' << n_p << ')' << std::endl;
 * 
 * @endcode
 * 
 * We split up the IndexSet for locally owned and locally relevant DoFs
 * into two IndexSets based on how we want to create the block matrices
 * and vectors.
 * 
 * @code
 *     owned_partitioning.resize(2);
 *     owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u);
 *     owned_partitioning[1] =
 *       dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p);
 * 
 *     IndexSet locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *     relevant_partitioning.resize(2);
 *     relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
 *     relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);
 * 
 * @endcode
 * 
 * Setting up the constraints for boundary conditions and hanging nodes
 * is identical to step-40. Rven though we don't have any hanging nodes
 * because we only perform global refinement, it is still a good idea
 * to put this function call in, in case adaptive refinement gets
 * introduced later.
 * 
 * @code
 *     {
 *       constraints.reinit(locally_relevant_dofs);
 * 
 *       FEValuesExtractors::Vector velocities(0);
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                ExactSolution<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 *       constraints.close();
 *     }
 * 
 * @endcode
 * 
 * Now we create the system matrix based on a BlockDynamicSparsityPattern.
 * We know that we won't have coupling between different velocity
 * components (because we use the laplace and not the deformation tensor)
 * and no coupling between pressure with its test functions, so we use
 * a Table to communicate this coupling information to
 * DoFTools::make_sparsity_pattern.
 * 
 * @code
 *     {
 *       system_matrix.clear();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (c == dim && d == dim)
 *             coupling[c][d] = DoFTools::none;
 *           else if (c == dim || d == dim || c == d)
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
 * 
 *       DoFTools::make_sparsity_pattern(
 *         dof_handler, coupling, dsp, constraints, false);
 * 
 *       SparsityTools::distribute_sparsity_pattern(
 *         dsp,
 *         dof_handler.locally_owned_dofs(),
 *         mpi_communicator,
 *         locally_relevant_dofs);
 * 
 *       system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
 *     }
 * 
 * @endcode
 * 
 * The preconditioner matrix has a different coupling (we only fill in
 * the 1,1 block with the mass matrix), otherwise this code is identical
 * to the construction of the system_matrix above.
 * 
 * @code
 *     {
 *       preconditioner_matrix.clear();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (c == dim && d == dim)
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
 * 
 *       DoFTools::make_sparsity_pattern(
 *         dof_handler, coupling, dsp, constraints, false);
 *       SparsityTools::distribute_sparsity_pattern(
 *         dsp,
 *         Utilities::MPI::all_gather(mpi_communicator,
 *                                    dof_handler.locally_owned_dofs()),
 *         mpi_communicator,
 *         locally_relevant_dofs);
 *       preconditioner_matrix.reinit(owned_partitioning,
 * @endcode
 * 
 * owned_partitioning,
 * 
 * @code
 *                                    dsp,
 *                                    mpi_communicator);
 *     }
 * 
 * @endcode
 * 
 * Finally, we construct the block vectors with the right sizes. The
 * function call with two std::vector<IndexSet> will create a ghosted
 * vector.
 * 
 * @code
 *     locally_relevant_solution.reinit(owned_partitioning,
 *                                      relevant_partitioning,
 *                                      mpi_communicator);
 *     system_rhs.reinit(owned_partitioning, mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Assembly"></a> 
 * <h3>Assembly</h3>
 *   

 * 
 * This function assembles the system matrix, the preconditioner matrix,
 * and the right hand side. The code is pretty standard.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::assemble_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "assembly");
 * 
 *     system_matrix         = 0;
 *     preconditioner_matrix = 0;
 *     system_rhs            = 0;
 * 
 *     const QGauss<dim> quadrature_formula(velocity_degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     const RightHandSide<dim>    right_hand_side;
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));
 * 
 *     std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
 *     std::vector<double>         div_phi_u(dofs_per_cell);
 *     std::vector<double>         phi_p(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 *     const FEValuesExtractors::Vector     velocities(0);
 *     const FEValuesExtractors::Scalar     pressure(dim);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix  = 0;
 *           cell_matrix2 = 0;
 *           cell_rhs     = 0;
 * 
 *           fe_values.reinit(cell);
 *           right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
 *                                             rhs_values);
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             {
 *               for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *                 {
 *                   grad_phi_u[k] = fe_values[velocities].gradient(k, q);
 *                   div_phi_u[k]  = fe_values[velocities].divergence(k, q);
 *                   phi_p[k]      = fe_values[pressure].value(k, q);
 *                 }
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     {
 *                       cell_matrix(i, j) +=
 *                         (viscosity *
 *                            scalar_product(grad_phi_u[i], grad_phi_u[j]) -
 *                          div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
 *                         fe_values.JxW(q);
 * 
 *                       cell_matrix2(i, j) += 1.0 / viscosity * phi_p[i] *
 *                                             phi_p[j] * fe_values.JxW(q);
 *                     }
 * 
 *                   const unsigned int component_i =
 *                     fe.system_to_component_index(i).first;
 *                   cell_rhs(i) += fe_values.shape_value(i, q) *
 *                                  rhs_values[q](component_i) * fe_values.JxW(q);
 *                 }
 *             }
 * 
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints.distribute_local_to_global(cell_matrix,
 *                                                  cell_rhs,
 *                                                  local_dof_indices,
 *                                                  system_matrix,
 *                                                  system_rhs);
 * 
 *           constraints.distribute_local_to_global(cell_matrix2,
 *                                                  local_dof_indices,
 *                                                  preconditioner_matrix);
 *         }
 * 
 *     system_matrix.compress(VectorOperation::add);
 *     preconditioner_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solving"></a> 
 * <h3>Solving</h3>
 *   

 * 
 * This function solves the linear system with MINRES with a block diagonal
 * preconditioner and AMG for the two diagonal blocks as described in the
 * introduction. The preconditioner applies a v cycle to the 0,0 block
 * and a CG with the mass matrix for the 1,1 block (the Schur complement).
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::solve()
 *   {
 *     TimerOutput::Scope t(computing_timer, "solve");
 * 
 *     LA::MPI::PreconditionAMG prec_A;
 *     {
 *       LA::MPI::PreconditionAMG::AdditionalData data;
 * 
 * #ifdef USE_PETSC_LA
 *       data.symmetric_operator = true;
 * #endif
 *       prec_A.initialize(system_matrix.block(0, 0), data);
 *     }
 * 
 *     LA::MPI::PreconditionAMG prec_S;
 *     {
 *       LA::MPI::PreconditionAMG::AdditionalData data;
 * 
 * #ifdef USE_PETSC_LA
 *       data.symmetric_operator = true;
 * #endif
 *       prec_S.initialize(preconditioner_matrix.block(1, 1), data);
 *     }
 * 
 * @endcode
 * 
 * The InverseMatrix is used to solve for the mass matrix:
 * 
 * @code
 *     using mp_inverse_t = LinearSolvers::InverseMatrix<LA::MPI::SparseMatrix,
 *                                                       LA::MPI::PreconditionAMG>;
 *     const mp_inverse_t mp_inverse(preconditioner_matrix.block(1, 1), prec_S);
 * 
 * @endcode
 * 
 * This constructs the block preconditioner based on the preconditioners
 * for the individual blocks defined above.
 * 
 * @code
 *     const LinearSolvers::BlockDiagonalPreconditioner<LA::MPI::PreconditionAMG,
 *                                                      mp_inverse_t>
 *       preconditioner(prec_A, mp_inverse);
 * 
 * @endcode
 * 
 * With that, we can finally set up the linear solver and solve the system:
 * 
 * @code
 *     SolverControl solver_control(system_matrix.m(),
 *                                  1e-10 * system_rhs.l2_norm());
 * 
 *     SolverMinRes<LA::MPI::BlockVector> solver(solver_control);
 * 
 *     LA::MPI::BlockVector distributed_solution(owned_partitioning,
 *                                               mpi_communicator);
 * 
 *     constraints.set_zero(distributed_solution);
 * 
 *     solver.solve(system_matrix,
 *                  distributed_solution,
 *                  system_rhs,
 *                  preconditioner);
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *           << std::endl;
 * 
 *     constraints.distribute(distributed_solution);
 * 
 * @endcode
 * 
 * Like in step-56, we subtract the mean pressure to allow error
 * computations against our reference solution, which has a mean value
 * of zero.
 * 
 * @code
 *     locally_relevant_solution = distributed_solution;
 *     const double mean_pressure =
 *       VectorTools::compute_mean_value(dof_handler,
 *                                       QGauss<dim>(velocity_degree + 2),
 *                                       locally_relevant_solution,
 *                                       dim);
 *     distributed_solution.block(1).add(-mean_pressure);
 *     locally_relevant_solution.block(1) = distributed_solution.block(1);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Therest"></a> 
 * <h3>The rest</h3>
 *   

 * 
 * The remainder of the code that deals with mesh refinement, output, and
 * the main loop is pretty standard.
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::refine_grid()
 *   {
 *     TimerOutput::Scope t(computing_timer, "refine");
 * 
 *     triangulation.refine_global();
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void StokesProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     {
 *       const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
 *       const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
 *                                                        dim + 1);
 * 
 *       Vector<double> cellwise_errors(triangulation.n_active_cells());
 *       QGauss<dim>    quadrature(velocity_degree + 2);
 * 
 *       VectorTools::integrate_difference(dof_handler,
 *                                         locally_relevant_solution,
 *                                         ExactSolution<dim>(),
 *                                         cellwise_errors,
 *                                         quadrature,
 *                                         VectorTools::L2_norm,
 *                                         &velocity_mask);
 * 
 *       const double error_u_l2 =
 *         VectorTools::compute_global_error(triangulation,
 *                                           cellwise_errors,
 *                                           VectorTools::L2_norm);
 * 
 *       VectorTools::integrate_difference(dof_handler,
 *                                         locally_relevant_solution,
 *                                         ExactSolution<dim>(),
 *                                         cellwise_errors,
 *                                         quadrature,
 *                                         VectorTools::L2_norm,
 *                                         &pressure_mask);
 * 
 *       const double error_p_l2 =
 *         VectorTools::compute_global_error(triangulation,
 *                                           cellwise_errors,
 *                                           VectorTools::L2_norm);
 * 
 *       pcout << "error: u_0: " << error_u_l2 << " p_0: " << error_p_l2
 *             << std::endl;
 *     }
 * 
 * 
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     solution_names.emplace_back("pressure");
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(locally_relevant_solution,
 *                              solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 * 
 *     LA::MPI::BlockVector interpolated;
 *     interpolated.reinit(owned_partitioning, MPI_COMM_WORLD);
 *     VectorTools::interpolate(dof_handler, ExactSolution<dim>(), interpolated);
 * 
 *     LA::MPI::BlockVector interpolated_relevant(owned_partitioning,
 *                                                relevant_partitioning,
 *                                                MPI_COMM_WORLD);
 *     interpolated_relevant = interpolated;
 *     {
 *       std::vector<std::string> solution_names(dim, "ref_u");
 *       solution_names.emplace_back("ref_p");
 *       data_out.add_data_vector(interpolated_relevant,
 *                                solution_names,
 *                                DataOut<dim>::type_dof_data,
 *                                data_component_interpretation);
 *     }
 * 
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells());
 *     for (unsigned int i = 0; i < subdomain.size(); ++i)
 *       subdomain(i) = triangulation.locally_owned_subdomain();
 *     data_out.add_data_vector(subdomain, "subdomain");
 * 
 *     data_out.build_patches();
 * 
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", cycle, mpi_communicator, 2);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void StokesProblem<dim>::run()
 *   {
 * #ifdef USE_PETSC_LA
 *     pcout << "Running using PETSc." << std::endl;
 * #else
 *     pcout << "Running using Trilinos." << std::endl;
 * #endif
 *     const unsigned int n_cycles = 5;
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           make_grid();
 *         else
 *           refine_grid();
 * 
 *         setup_system();
 * 
 *         assemble_system();
 *         solve();
 * 
 *         if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
 *           {
 *             TimerOutput::Scope t(computing_timer, "output");
 *             output_results(cycle);
 *           }
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 * 
 *         pcout << std::endl;
 *       }
 *   }
 * } // namespace Step55
 * 
 * 
 * 
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step55;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       StokesProblem<2> problem(2);
 *       problem.run();
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


As expected from the discussion above, the number of iterations is independent
of the number of processors and only very slightly dependent on $h$:

<table>
<tr>
  <th colspan="2">PETSc</th>
  <th colspan="8">number of processors</th>
</tr>
<tr>
  <th>cycle</th>
  <th>dofs</th>
  <th>1</th>
  <th>2</th>
  <th>4</th>
  <th>8</th>
  <th>16</th>
  <th>32</th>
  <th>64</th>
  <th>128</th>
</tr>
<tr>
  <td>0</td>
  <td>659</td>
  <td>49</td>
  <td>49</td>
  <td>49</td>
  <td>51</td>
  <td>51</td>
  <td>51</td>
  <td>49</td>
  <td>49</td>
</tr>
<tr>
  <td>1</td>
  <td>2467</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>54</td>
  <td>54</td>
  <td>53</td>
</tr>
<tr>
  <td>2</td>
  <td>9539</td>
  <td>56</td>
  <td>56</td>
  <td>56</td>
  <td>54</td>
  <td>56</td>
  <td>56</td>
  <td>54</td>
  <td>56</td>
</tr>
<tr>
  <td>3</td>
  <td>37507</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>56</td>
  <td>57</td>
  <td>56</td>
</tr>
<tr>
  <td>4</td>
  <td>148739</td>
  <td>58</td>
  <td>59</td>
  <td>57</td>
  <td>59</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
</tr>
<tr>
  <td>5</td>
  <td>592387</td>
  <td>60</td>
  <td>60</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
</tr>
<tr>
  <td>6</td>
  <td>2364419</td>
  <td>62</td>
  <td>62</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
</tr>
</table>

<table>
<tr>
  <th colspan="2">Trilinos</th>
  <th colspan="8">number of processors</th>
</tr>
<tr>
  <th>cycle</th>
  <th>dofs</th>
  <th>1</th>
  <th>2</th>
  <th>4</th>
  <th>8</th>
  <th>16</th>
  <th>32</th>
  <th>64</th>
  <th>128</th>
</tr>
<tr>
  <td>0</td>
  <td>659</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
</tr>
<tr>
  <td>1</td>
  <td>2467</td>
  <td>92</td>
  <td>89</td>
  <td>89</td>
  <td>82</td>
  <td>86</td>
  <td>81</td>
  <td>78</td>
  <td>78</td>
</tr>
<tr>
  <td>2</td>
  <td>9539</td>
  <td>102</td>
  <td>99</td>
  <td>96</td>
  <td>95</td>
  <td>95</td>
  <td>88</td>
  <td>83</td>
  <td>95</td>
</tr>
<tr>
  <td>3</td>
  <td>37507</td>
  <td>107</td>
  <td>105</td>
  <td>104</td>
  <td>99</td>
  <td>100</td>
  <td>96</td>
  <td>96</td>
  <td>90</td>
</tr>
<tr>
  <td>4</td>
  <td>148739</td>
  <td>112</td>
  <td>112</td>
  <td>111</td>
  <td>111</td>
  <td>127</td>
  <td>126</td>
  <td>115</td>
  <td>117</td>
</tr>
<tr>
  <td>5</td>
  <td>592387</td>
  <td>116</td>
  <td>115</td>
  <td>114</td>
  <td>112</td>
  <td>118</td>
  <td>120</td>
  <td>131</td>
  <td>130</td>
</tr>
<tr>
  <td>6</td>
  <td>2364419</td>
  <td>130</td>
  <td>126</td>
  <td>120</td>
  <td>120</td>
  <td>121</td>
  <td>122</td>
  <td>121</td>
  <td>123</td>
</tr>
</table>

While the PETSc results show a constant number of iterations, the iterations
increase when using Trilinos. This is likely because of the different settings
used for the AMG preconditioner. For performance reasons we do not allow
coarsening below a couple thousand unknowns. As the coarse solver is an exact
solve (we are using LU by default), a change in number of levels will
influence the quality of a V-cycle. Therefore, a V-cycle is closer to an exact
solver for smaller problem sizes.

<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="InvestigateTrilinositerations"></a><h4>Investigate Trilinos iterations</h4>


Play with the smoothers, smoothing steps, and other properties for the
Trilinos AMG to achieve an optimal preconditioner.

<a name="SolvetheOseenprobleminsteadoftheStokessystem"></a><h4>Solve the Oseen problem instead of the Stokes system</h4>


This change requires changing the outer solver to GMRES or BiCGStab, because
the system is no longer symmetric.

You can prescribe the exact flow solution as $b$ in the convective term $b
\cdot \nabla u$. This should give the same solution as the original problem,
if you set the right hand side to zero.

<a name="Adaptiverefinement"></a><h4>Adaptive refinement</h4>


So far, this tutorial program refines the mesh globally in each step.
Replacing the code in StokesProblem::refine_grid() by something like
@code
Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

FEValuesExtractors::Vector velocities(0);
KellyErrorEstimator<dim>::estimate(
  dof_handler,
  QGauss<dim - 1>(fe.degree + 1),
  std::map<types::boundary_id, const Function<dim> *>(),
  locally_relevant_solution,
  estimated_error_per_cell,
  fe.component_mask(velocities));
parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
  triangulation, estimated_error_per_cell, 0.3, 0.0);
triangulation.execute_coarsening_and_refinement();
@endcode
makes it simple to explore adaptive mesh refinement.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-55.cc"
*/
