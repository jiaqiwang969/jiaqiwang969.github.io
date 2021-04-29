/**
@page step_16b The step-16b tutorial program
This tutorial depends on step-16.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Theintegratoroneachcell">The integrator on each cell</a>
        <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The <code>LaplaceProblem</code> class template</a>
        <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The <code>LaplaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemassemble_multigrid">LaplaceProblem::assemble_multigrid</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#Postprocessing">Postprocessing</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This is a variant of step-16 with the only change that we are using the
MeshWorker framework with the pre-made LocalIntegrator helper classes instead
of manually assembling the matrices.

The details of this framework on how it is used in practice will be explained
as part of this tutorial program.

<a name="Thetestcase"></a><h3>The testcase</h3>


The problem we solve here is the same as the one in step-16.
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
 * Again, the first few include files are already known, so we won't comment
 * on them:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * These, now, are the include necessary for the multilevel methods. The first
 * one declares how to handle Dirichlet boundary conditions on each of the
 * levels of the multigrid method. For the actual description of the degrees
 * of freedom, we do not need any new include file because DoFHandler already
 * has all necessary methods implemented. We will only need to distribute the
 * DoFs for the levels further down.
 * 

 * 
 * The rest of the include files deals with the mechanics of multigrid as a
 * linear operator (solver or preconditioner).
 * 
 * @code
 * #include <deal.II/multigrid/mg_constrained_dofs.h>
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * @endcode
 * 
 * Finally we include the MeshWorker framework. This framework through its
 * function loop() and integration_loop(), automates loops over cells and
 * assembling of data into vectors, matrices, etc. It obeys constraints
 * automatically. Since we have to build several matrices and have to be aware
 * of several sets of constraints, this will save us a lot of headache.
 * 
 * @code
 * #include <deal.II/meshworker/dof_info.h>
 * #include <deal.II/meshworker/integration_info.h>
 * #include <deal.II/meshworker/simple.h>
 * #include <deal.II/meshworker/output.h>
 * #include <deal.II/meshworker/loop.h>
 * 
 * @endcode
 * 
 * In order to save effort, we use the pre-implemented Laplacian found in
 * 
 * @code
 * #include <deal.II/integrators/laplace.h>
 * #include <deal.II/integrators/l2.h>
 * 
 * @endcode
 * 
 * This is C++:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * 
 * using namespace dealii;
 * 
 * namespace Step16
 * {
 * @endcode
 * 
 * 
 * <a name="Theintegratoroneachcell"></a> 
 * <h3>The integrator on each cell</h3>
 * 

 * 
 * The MeshWorker::integration_loop() expects a class that provides functions
 * for integration on cells and boundary and interior faces. This is done by
 * the following class. In the constructor, we tell the loop that cell
 * integrals should be computed (the 'true'), but integrals should not be
 * computed on boundary and interior faces (the two 'false'). Accordingly, we
 * only need a cell function, but none for the faces.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceIntegrator : public MeshWorker::LocalIntegrator<dim>
 *   {
 *   public:
 *     LaplaceIntegrator();
 *     virtual void cell(MeshWorker::DoFInfo<dim> &        dinfo,
 *                       MeshWorker::IntegrationInfo<dim> &info) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   LaplaceIntegrator<dim>::LaplaceIntegrator()
 *     : MeshWorker::LocalIntegrator<dim>(true, false, false)
 *   {}
 * 
 * 
 * @endcode
 * 
 * Next the actual integrator on each cell. We solve a Poisson problem with a
 * coefficient one in the right half plane and one tenth in the left
 * half plane.
 * 

 * 
 * The MeshWorker::LocalResults base class of MeshWorker::DoFInfo contains
 * objects that can be filled in this local integrator. How many objects are
 * created is determined inside the MeshWorker framework by the assembler
 * class. Here, we test for instance that one matrix is required
 * (MeshWorker::LocalResults::n_matrices()). The matrices are accessed through
 * MeshWorker::LocalResults::matrix(), which takes the number of the matrix as
 * its first argument. The second argument is only used for integrals over
 * faces when there are two matrices for each test function used. Then, a
 * second matrix with indicator 'true' would exist with the same index.
 * 

 * 
 * MeshWorker::IntegrationInfo provides one or several FEValues objects, which
 * below are used by LocalIntegrators::Laplace::cell_matrix() or
 * LocalIntegrators::L2::L2(). Since we are assembling only a single PDE,
 * there is also only one of these objects with index zero.
 * 

 * 
 * In addition, we note that this integrator serves to compute the matrices
 * for the multilevel preconditioner as well as the matrix and the right hand
 * side for the global system. Since the assembler for a system requires an
 * additional vector, MeshWorker::LocalResults::n_vectors() is returning a
 * nonzero value. Accordingly, we fill a right hand side vector at the end of
 * this function. Since LocalResults can deal with several BlockVector
 * objects, but we are again in the simplest case here, we enter the
 * information into block zero of vector zero.
 * 
 * @code
 *   template <int dim>
 *   void
 *   LaplaceIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &        dinfo,
 *                                MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     AssertDimension(dinfo.n_matrices(), 1);
 *     const double coefficient = (dinfo.cell->center()(0) > 0.) ? .1 : 1.;
 * 
 *     LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix,
 *                                            info.fe_values(0),
 *                                            coefficient);
 * 
 *     if (dinfo.n_vectors() > 0)
 *       {
 *         std::vector<double> rhs(info.fe_values(0).n_quadrature_points, 1.);
 *         LocalIntegrators::L2::L2(dinfo.vector(0).block(0),
 *                                  info.fe_values(0),
 *                                  rhs);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceProblem</code> class template</h3>
 * 

 * 
 * This main class is basically the same class as in step-6. As far as
 * member functions is concerned, the only addition is the
 * <code>assemble_multigrid</code> function that assembles the matrices that
 * correspond to the discrete operators on intermediate levels:
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void assemble_multigrid();
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     const unsigned int degree;
 * 
 * @endcode
 * 
 * The following members are the essential data structures for the multigrid
 * method. The first two represent the sparsity patterns and the matrices on
 * individual levels of the multilevel hierarchy, very much like the objects
 * for the global mesh above.
 *     

 * 
 * Then we have two new matrices only needed for multigrid methods with
 * local smoothing on adaptive meshes. They convey data between the interior
 * part of the refined region and the refinement edge, as outlined in detail
 * in the @ref mg_paper "multigrid paper".
 *     

 * 
 * The last object stores information about the boundary indices on each
 * level and information about indices lying on a refinement edge between
 * two different refinement levels. It thus serves a similar purpose as
 * AffineConstraints, but on each level.
 * 
 * @code
 *     MGLevelObject<SparsityPattern>      mg_sparsity_patterns;
 *     MGLevelObject<SparseMatrix<double>> mg_matrices;
 *     MGLevelObject<SparseMatrix<double>> mg_interface_in;
 *     MGLevelObject<SparseMatrix<double>> mg_interface_out;
 *     MGConstrainedDoFs                   mg_constrained_dofs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>LaplaceProblem</code> class implementation</h3>
 * 

 * 
 * Just one short remark about the constructor of the Triangulation:
 * by convention, all adaptively refined triangulations in deal.II never
 * change by more than one level across a face between cells. For our
 * multigrid algorithms, however, we need a slightly stricter guarantee,
 * namely that the mesh also does not change by more than refinement level
 * across vertices that might connect two cells. In other words, we must
 * prevent the following situation:
 *   

 * 
 * @image html limit_level_difference_at_vertices.png ""
 *   

 * 
 * This is achieved by passing the
 * Triangulation::limit_level_difference_at_vertices flag to the constructor
 * of the triangulation class.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree)
 *     : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     , fe(degree)
 *     , dof_handler(triangulation)
 *     , degree(degree)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * In addition to just distributing the degrees of freedom in
 * the DoFHandler, we do the same on each level. Then, we follow the
 * same procedure as before to set up the system on the leaf mesh.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler.distribute_mg_dofs();
 * 
 *     deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << " (by level: ";
 *     for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
 *       deallog << dof_handler.n_dofs(level)
 *               << (level == triangulation.n_levels() - 1 ? ")" : ", ");
 *     deallog << std::endl;
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     std::set<types::boundary_id> dirichlet_boundary_ids = {0};
 *     Functions::ZeroFunction<dim> homogeneous_dirichlet_bc;
 *     const std::map<types::boundary_id, const Function<dim> *>
 *       dirichlet_boundary_functions = {
 *         {types::boundary_id(0), &homogeneous_dirichlet_bc}};
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              dirichlet_boundary_functions,
 *                                              constraints);
 *     constraints.close();
 *     constraints.condense(dsp);
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 * @endcode
 * 
 * The multigrid constraints have to be initialized. They need to know
 * about the boundary values as well, so we pass the
 * <code>dirichlet_boundary</code> here as well.
 * 
 * @code
 *     mg_constrained_dofs.clear();
 *     mg_constrained_dofs.initialize(dof_handler);
 *     mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
 *                                                        dirichlet_boundary_ids);
 * 
 * 
 * @endcode
 * 
 * Now for the things that concern the multigrid data structures. First, we
 * resize the multilevel objects to hold matrices and sparsity patterns for
 * every level. The coarse level is zero (this is mandatory right now but
 * may change in a future revision). Note that these functions take a
 * complete, inclusive range here (not a starting index and size), so the
 * finest level is <code>n_levels-1</code>. We first have to resize the
 * container holding the SparseMatrix classes, since they have to release
 * their SparsityPattern before the can be destroyed upon resizing.
 * 
 * @code
 *     const unsigned int n_levels = triangulation.n_levels();
 * 
 *     mg_interface_in.resize(0, n_levels - 1);
 *     mg_interface_in.clear_elements();
 *     mg_interface_out.resize(0, n_levels - 1);
 *     mg_interface_out.clear_elements();
 *     mg_matrices.resize(0, n_levels - 1);
 *     mg_matrices.clear_elements();
 *     mg_sparsity_patterns.resize(0, n_levels - 1);
 * 
 * @endcode
 * 
 * Now, we have to provide a matrix on each level. To this end, we first use
 * the MGTools::make_sparsity_pattern function to generate a preliminary
 * compressed sparsity pattern on each level (see the @ref Sparsity module
 * for more information on this topic) and then copy it over to the one we
 * really want. The next step is to initialize both kinds of level matrices
 * with these sparsity patterns.
 *     

 * 
 * It may be worth pointing out that the interface matrices only have
 * entries for degrees of freedom that sit at or next to the interface
 * between coarser and finer levels of the mesh. They are therefore even
 * sparser than the matrices on the individual levels of our multigrid
 * hierarchy. If we were more concerned about memory usage (and possibly the
 * speed with which we can multiply with these matrices), we should use
 * separate and different sparsity patterns for these two kinds of matrices.
 * 
 * @code
 *     for (unsigned int level = 0; level < n_levels; ++level)
 *       {
 *         DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
 *                                    dof_handler.n_dofs(level));
 *         MGTools::make_sparsity_pattern(dof_handler, dsp, level);
 * 
 *         mg_sparsity_patterns[level].copy_from(dsp);
 * 
 *         mg_matrices[level].reinit(mg_sparsity_patterns[level]);
 *         mg_interface_in[level].reinit(mg_sparsity_patterns[level]);
 *         mg_interface_out[level].reinit(mg_sparsity_patterns[level]);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system</h4>
 * 

 * 
 * The following function assembles the linear system on the finest level of
 * the mesh. Since we want to reuse the code here for the level assembling
 * below, we use the local integrator class LaplaceIntegrator and leave the
 * loops to the MeshWorker framework. Thus, this function first sets up the
 * objects necessary for this framework, namely
 * - a MeshWorker::IntegrationInfoBox object, which will provide all the
 * required data in quadrature points on the cell. This object can be seen
 * as an extension of FEValues, providing a lot more useful information,
 * - a MeshWorker::DoFInfo object, which on the one hand side extends the
 * functionality of cell iterators, but also provides space for return
 * values in its base class LocalResults,
 * - an assembler, in this case for the whole system. The term 'simple' here
 * refers to the fact that the global system does not have a block
 * structure,
 * - the local integrator, which implements the actual forms.
 *   

 * 
 * After the loop has combined all of these into a matrix and a right hand
 * side, there is one thing left to do: the assemblers leave matrix rows and
 * columns of constrained degrees of freedom untouched. Therefore, we put a
 * one on the diagonal to make the whole system well posed. The value one, or
 * any fixed value has the advantage, that its effect on the spectrum of the
 * matrix is easily understood. Since the corresponding eigenvectors form an
 * invariant subspace, the value chosen does not affect the convergence of
 * Krylov space solvers.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_system()
 *   {
 *     MappingQ1<dim>                      mapping;
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     UpdateFlags                         update_flags =
 *       update_values | update_gradients | update_hessians;
 *     info_box.add_update_flags_all(update_flags);
 *     info_box.initialize(fe, mapping);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 *     MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>>
 *       assembler;
 *     assembler.initialize(constraints);
 *     assembler.initialize(system_matrix, system_rhs);
 * 
 *     LaplaceIntegrator<dim> matrix_integrator;
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
 *                                            dof_handler.end(),
 *                                            dof_info,
 *                                            info_box,
 *                                            matrix_integrator,
 *                                            assembler);
 * 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *       if (constraints.is_constrained(i))
 *         system_matrix.set(i, i, 1.);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_multigrid"></a> 
 * <h4>LaplaceProblem::assemble_multigrid</h4>
 * 

 * 
 * The next function is the one that builds the linear operators (matrices)
 * that define the multigrid method on each level of the mesh. The integration
 * core is the same as above, but the loop below will go over all existing
 * cells instead of just the active ones, and the results must be entered into
 * the correct level matrices. Fortunately, MeshWorker hides most of that from
 * us, and thus the difference between this function and the previous lies
 * only in the setup of the assembler and the different iterators in the loop.
 * Also, fixing up the matrices in the end is a little more complicated.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_multigrid()
 *   {
 *     MappingQ1<dim>                      mapping;
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     UpdateFlags                         update_flags =
 *       update_values | update_gradients | update_hessians;
 *     info_box.add_update_flags_all(update_flags);
 *     info_box.initialize(fe, mapping);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 *     MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler;
 *     assembler.initialize(mg_constrained_dofs);
 *     assembler.initialize(mg_matrices);
 *     assembler.initialize_interfaces(mg_interface_in, mg_interface_out);
 * 
 *     LaplaceIntegrator<dim> matrix_integrator;
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(),
 *                                            dof_handler.end_mg(),
 *                                            dof_info,
 *                                            info_box,
 *                                            matrix_integrator,
 *                                            assembler);
 * 
 *     const unsigned int nlevels = triangulation.n_levels();
 *     for (unsigned int level = 0; level < nlevels; ++level)
 *       {
 *         for (unsigned int i = 0; i < dof_handler.n_dofs(level); ++i)
 *           if (mg_constrained_dofs.is_boundary_index(level, i) ||
 *               mg_constrained_dofs.at_refinement_edge(level, i))
 *             mg_matrices[level].set(i, i, 1.);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * This is the other function that is significantly different in support of
 * the multigrid solver (or, in fact, the preconditioner for which we use
 * the multigrid method).
 *   

 * 
 * Let us start out by setting up two of the components of multilevel
 * methods: transfer operators between levels, and a solver on the coarsest
 * level. In finite element methods, the transfer operators are derived from
 * the finite element function spaces involved and can often be computed in
 * a generic way independent of the problem under consideration. In that
 * case, we can use the MGTransferPrebuilt class that, given the constraints
 * of the final linear system and the MGConstrainedDoFs object that knows
 * about the boundary conditions on the each level and the degrees of
 * freedom on interfaces between different refinement level can build the
 * matrices for those transfer operations from a DoFHandler object with
 * level degrees of freedom.
 *   

 * 
 * The second part of the following lines deals with the coarse grid
 * solver. Since our coarse grid is very coarse indeed, we decide for a
 * direct solver (a Householder decomposition of the coarsest level matrix),
 * even if its implementation is not particularly sophisticated. If our
 * coarse mesh had many more cells than the five we have here, something
 * better suited would obviously be necessary here.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs);
 *     mg_transfer.build(dof_handler);
 * 
 *     FullMatrix<double> coarse_matrix;
 *     coarse_matrix.copy_from(mg_matrices[0]);
 *     MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
 *     coarse_grid_solver.initialize(coarse_matrix);
 * 
 * @endcode
 * 
 * The next component of a multilevel solver or preconditioner is that we
 * need a smoother on each level. A common choice for this is to use the
 * application of a relaxation method (such as the SOR, Jacobi or Richardson
 * method) or a small number of iterations of a solver method (such as CG or
 * GMRES). The mg::SmootherRelaxation and MGSmootherPrecondition classes
 * provide support for these two kinds of smoothers. Here, we opt for the
 * application of a single SOR iteration. To this end, we define an
 * appropriate alias and then setup a smoother object.
 *     

 * 
 * The last step is to initialize the smoother object with our level
 * matrices and to set some smoothing parameters. The
 * <code>initialize()</code> function can optionally take additional
 * arguments that will be passed to the smoother object on each level. In
 * the current case for the SOR smoother, this could, for example, include
 * a relaxation parameter. However, we here leave these at their default
 * values. The call to <code>set_steps()</code> indicates that we will use
 * two pre- and two post-smoothing steps on each level; to use a variable
 * number of smoother steps on different levels, more options can be set
 * in the constructor call to the <code>mg_smoother</code> object.
 *     

 * 
 * The last step results from the fact that we use the SOR method as a
 * smoother - which is not symmetric - but we use the conjugate gradient
 * iteration (which requires a symmetric preconditioner) below, we need to
 * let the multilevel preconditioner make sure that we get a symmetric
 * operator even for nonsymmetric smoothers:
 * 
 * @code
 *     using Smoother = PreconditionSOR<SparseMatrix<double>>;
 *     mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother;
 *     mg_smoother.initialize(mg_matrices);
 *     mg_smoother.set_steps(2);
 *     mg_smoother.set_symmetric(true);
 * 
 * @endcode
 * 
 * The next preparatory step is that we must wrap our level and interface
 * matrices in an object having the required multiplication functions. We
 * will create two objects for the interface objects going from coarse to
 * fine and the other way around; the multigrid algorithm will later use
 * the transpose operator for the latter operation, allowing us to
 * initialize both up and down versions of the operator with the matrices
 * we already built:
 * 
 * @code
 *     mg::Matrix<Vector<double>> mg_matrix(mg_matrices);
 *     mg::Matrix<Vector<double>> mg_interface_up(mg_interface_in);
 *     mg::Matrix<Vector<double>> mg_interface_down(mg_interface_out);
 * 
 * @endcode
 * 
 * Now, we are ready to set up the V-cycle operator and the multilevel
 * preconditioner.
 * 
 * @code
 *     Multigrid<Vector<double>> mg(
 *       mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother);
 *     mg.set_edge_matrices(mg_interface_down, mg_interface_up);
 * 
 *     PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
 *       preconditioner(dof_handler, mg, mg_transfer);
 * 
 * @endcode
 * 
 * With all this together, we can finally get about solving the linear
 * system in the usual way:
 * 
 * @code
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     solution = 0;
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Postprocessing"></a> 
 * <h4>Postprocessing</h4>
 * 

 * 
 * The following two functions postprocess a solution once it is
 * computed. In particular, the first one refines the mesh at the beginning
 * of each cycle while the second one outputs results at the end of each
 * such cycle. The functions are almost unchanged from those in step-6, with
 * the exception of one minor difference: we generate output in VTK
 * format, to use the more modern visualization programs available today
 * compared to those that were available when step-6 was written.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::refine_grid()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       solution,
 *       estimated_error_per_cell);
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches();
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * Like several of the functions above, this is almost exactly a copy of
 * the corresponding function in step-6. The only difference is the call to
 * <code>assemble_multigrid</code> that takes care of forming the matrices
 * on every level that we need in the multigrid method.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 8; ++cycle)
 *       {
 *         deallog << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_ball(triangulation);
 *             triangulation.refine_global(1);
 *           }
 *         else
 *           refine_grid();
 * 
 *         deallog << "   Number of active cells:       "
 *                 << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         assemble_system();
 *         assemble_multigrid();
 * 
 *         solve();
 *         output_results(cycle);
 *       }
 *   }
 * } // namespace Step16
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * This is again the same function as in step-6:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step16;
 * 
 *       deallog.depth_console(2);
 * 
 *       LaplaceProblem<2> laplace_problem(1);
 *       laplace_problem.run();
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


As in step-16, the solution looks like this on the finest mesh:

<p align="center">
  <img src="https://www.dealii.org/images/steps/developer/step-16.solution.png" alt="">
</p>

The output is formatted in a slightly different way compared to step-16 but is
functionally the same and shows the same convergence properties:
<pre>
DEAL::Cycle 0
DEAL::   Number of active cells:       20
DEAL::   Number of degrees of freedom: 25 (by level: 8, 25)
DEAL:cg::Starting value 0.510691
DEAL:cg::Convergence step 6 value 4.59193e-14
DEAL::Cycle 1
DEAL::   Number of active cells:       44
DEAL::   Number of degrees of freedom: 55 (by level: 8, 25, 45)
DEAL:cg::Starting value 0.440678
DEAL:cg::Convergence step 8 value 1.99419e-13
DEAL::Cycle 2
DEAL::   Number of active cells:       86
DEAL::   Number of degrees of freedom: 105 (by level: 8, 25, 69, 49)
DEAL:cg::Starting value 0.371855
DEAL:cg::Convergence step 9 value 1.13984e-13
DEAL::Cycle 3
DEAL::   Number of active cells:       170
DEAL::   Number of degrees of freedom: 200 (by level: 8, 25, 77, 174)
DEAL:cg::Starting value 0.318967
DEAL:cg::Convergence step 9 value 2.62112e-13
DEAL::Cycle 4
DEAL::   Number of active cells:       332
DEAL::   Number of degrees of freedom: 388 (by level: 8, 25, 86, 231, 204)
DEAL:cg::Starting value 0.276534
DEAL:cg::Convergence step 10 value 1.69562e-13
DEAL::Cycle 5
DEAL::   Number of active cells:       632
DEAL::   Number of degrees of freedom: 714 (by level: 8, 25, 89, 231, 514, 141)
DEAL:cg::Starting value 0.215300
DEAL:cg::Convergence step 10 value 6.47463e-13
DEAL::Cycle 6
DEAL::   Number of active cells:       1202
DEAL::   Number of degrees of freedom: 1332 (by level: 8, 25, 89, 282, 771, 435, 257)
DEAL:cg::Starting value 0.175848
DEAL:cg::Convergence step 10 value 1.80664e-13
DEAL::Cycle 7
DEAL::   Number of active cells:       2288
DEAL::   Number of degrees of freedom: 2511 (by level: 8, 25, 89, 318, 779, 1420, 829, 30)
DEAL:cg::Starting value 0.136724
DEAL:cg::Convergence step 11 value 9.73331e-14
</pre>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-16b.cc"
*/
