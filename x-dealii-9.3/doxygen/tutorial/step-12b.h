/**
@page step_12b The step-12b tutorial program
This tutorial depends on step-16, step-7, step-39.

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
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#TheAdvectionProblemclass">The AdvectionProblem class</a>
      <ul>
        <li><a href="#Theassemble_systemfunction">The assemble_system function</a>
        <li><a href="#Thelocalintegrators">The local integrators</a>
      </ul>
        <li><a href="#Alltherest">All the rest</a>
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


This is a variant of step-16 with the only change that we are using the
MeshWorker framework with the pre-made LocalIntegrator helper classes instead
of assembling the face terms using FEInterfaceValues.

The details of this framework on how it is used in practice will be explained
as part of this tutorial program.

<a name="Thetestcase"></a><h3>The testcase</h3>


The problem we solve here is the same as the one in step-12.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The first few files have already been covered in previous examples and will
 * thus not be further commented on:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/fe/mapping_q1.h>
 * @endcode
 * 
 * Here the discontinuous finite elements are defined. They are used in the same
 * way as all other finite elements, though -- as you have seen in previous
 * tutorial programs -- there isn't much user interaction with finite element
 * classes at all: they are passed to <code>DoFHandler</code> and
 * <code>FEValues</code> objects, and that is about it.
 * 
 * @code
 * #include <deal.II/fe/fe_dgq.h>
 * @endcode
 * 
 * We are going to use the simplest possible solver, called Richardson
 * iteration, that represents a simple defect correction. This, in combination
 * with a block SSOR preconditioner (defined in precondition_block.h), that
 * uses the special block matrix structure of system matrices arising from DG
 * discretizations.
 * 
 * @code
 * #include <deal.II/lac/solver_richardson.h>
 * #include <deal.II/lac/precondition_block.h>
 * @endcode
 * 
 * We are going to use gradients as refinement indicator.
 * 
 * @code
 * #include <deal.II/numerics/derivative_approximation.h>
 * 
 * @endcode
 * 
 * Here come the new include files for using the MeshWorker framework. The first
 * contains the class MeshWorker::DoFInfo, which provides local integrators with
 * a mapping between local and global degrees of freedom. It stores the results
 * of local integrals as well in its base class MeshWorker::LocalResults.
 * In the second of these files, we find an object of type
 * MeshWorker::IntegrationInfo, which is mostly a wrapper around a group of
 * FEValues objects. The file <tt>meshworker/simple.h</tt> contains classes
 * assembling locally integrated data into a global system containing only a
 * single matrix. Finally, we will need the file that runs the loop over all
 * mesh cells and faces.
 * 
 * @code
 * #include <deal.II/meshworker/dof_info.h>
 * #include <deal.II/meshworker/integration_info.h>
 * #include <deal.II/meshworker/simple.h>
 * #include <deal.II/meshworker/loop.h>
 * 
 * @endcode
 * 
 * Like in all programs, we finish this section by including the needed C++
 * headers and declaring we want to use objects in the dealii namespace without
 * prefix.
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * 
 * 
 * namespace Step12
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 *   

 * 
 * First, we define a class describing the inhomogeneous boundary data. Since
 * only its values are used, we implement value_list(), but leave all other
 * functions of Function undefined.
 * 
 * @code
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     BoundaryValues() = default;
 *     virtual void value_list(const std::vector<Point<dim>> &points,
 *                             std::vector<double> &          values,
 *                             const unsigned int component = 0) const override;
 *   };
 * 
 * @endcode
 * 
 * Given the flow direction, the inflow boundary of the unit square $[0,1]^2$
 * are the right and the lower boundaries. We prescribe discontinuous boundary
 * values 1 and 0 on the x-axis and value 0 on the right boundary. The values
 * of this function on the outflow boundaries will not be used within the DG
 * scheme.
 * 
 * @code
 *   template <int dim>
 *   void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                        std::vector<double> &          values,
 *                                        const unsigned int component) const
 *   {
 *     (void)component;
 *     AssertIndexRange(component, 1);
 *     Assert(values.size() == points.size(),
 *            ExcDimensionMismatch(values.size(), points.size()));
 * 
 *     for (unsigned int i = 0; i < values.size(); ++i)
 *       {
 *         if (points[i](0) < 0.5)
 *           values[i] = 1.;
 *         else
 *           values[i] = 0.;
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally, a function that computes and returns the wind field
 * $\beta=\beta(\mathbf x)$. As explained in the introduction, we will use a
 * rotational field around the origin in 2d. In 3d, we simply leave the
 * $z$-component unset (i.e., at zero), whereas the function can not be used
 * in 1d in its current implementation:
 * 
 * @code
 *   template <int dim>
 *   Tensor<1, dim> beta(const Point<dim> &p)
 *   {
 *     Assert(dim >= 2, ExcNotImplemented());
 * 
 *     Point<dim> wind_field;
 *     wind_field(0) = -p(1);
 *     wind_field(1) = p(0);
 *     wind_field /= wind_field.norm();
 * 
 *     return wind_field;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheAdvectionProblemclass"></a> 
 * <h3>The AdvectionProblem class</h3>
 *   

 * 
 * After this preparations, we proceed with the main class of this program,
 * called AdvectionProblem. It is basically the main class of step-6. We do
 * not have an AffineConstraints object, because there are no hanging node
 * constraints in DG discretizations.
 * 

 * 
 * Major differences will only come up in the implementation of the assemble
 * functions, since here, we not only need to cover the flux integrals over
 * faces, we also use the MeshWorker interface to simplify the loops
 * involved.
 * 
 * @code
 *   template <int dim>
 *   class AdvectionProblem
 *   {
 *   public:
 *     AdvectionProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve(Vector<double> &solution);
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim>   triangulation;
 *     const MappingQ1<dim> mapping;
 * 
 * @endcode
 * 
 * Furthermore we want to use DG elements of degree 1 (but this is only
 * specified in the constructor). If you want to use a DG method of a
 * different degree the whole program stays the same, only replace 1 in
 * the constructor by the desired polynomial degree.
 * 
 * @code
 *     FE_DGQ<dim>     fe;
 *     DoFHandler<dim> dof_handler;
 * 
 * @endcode
 * 
 * The next four members represent the linear system to be solved.
 * <code>system_matrix</code> and <code>right_hand_side</code> are generated
 * by <code>assemble_system()</code>, the <code>solution</code> is computed
 * in <code>solve()</code>. The <code>sparsity_pattern</code> is used to
 * determine the location of nonzero elements in <code>system_matrix</code>.
 * 
 * @code
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> right_hand_side;
 * 
 * @endcode
 * 
 * Finally, we have to provide functions that assemble the cell, boundary,
 * and inner face terms. Within the MeshWorker framework, the loop over all
 * cells and much of the setup of operations will be done outside this
 * class, so all we have to provide are these three operations. They will
 * then work on intermediate objects for which first, we here define
 * alias to the info objects handed to the local integration functions
 * in order to make our life easier below.
 * 
 * @code
 *     using DoFInfo  = MeshWorker::DoFInfo<dim>;
 *     using CellInfo = MeshWorker::IntegrationInfo<dim>;
 * 
 * @endcode
 * 
 * The following three functions are then the ones that get called inside
 * the generic loop over all cells and faces. They are the ones doing the
 * actual integration.
 *     

 * 
 * In our code below, these functions do not access member variables of the
 * current class, so we can mark them as <code>static</code> and simply pass
 * pointers to these functions to the MeshWorker framework. If, however,
 * these functions would want to access member variables (or needed
 * additional arguments beyond the ones specified below), we could use the
 * facilities of lambda functions to provide the
 * MeshWorker framework with objects that act as if they had the required
 * number and types of arguments, but have in fact other arguments already
 * bound.
 * 
 * @code
 *     static void integrate_cell_term(DoFInfo &dinfo, CellInfo &info);
 *     static void integrate_boundary_term(DoFInfo &dinfo, CellInfo &info);
 *     static void integrate_face_term(DoFInfo & dinfo1,
 *                                     DoFInfo & dinfo2,
 *                                     CellInfo &info1,
 *                                     CellInfo &info2);
 *   };
 * 
 * 
 * @endcode
 * 
 * We start with the constructor. The 1 in the constructor call of
 * <code>fe</code> is the polynomial degree.
 * 
 * @code
 *   template <int dim>
 *   AdvectionProblem<dim>::AdvectionProblem()
 *     : mapping()
 *     , fe(1)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 *   template <int dim>
 *   void AdvectionProblem<dim>::setup_system()
 *   {
 * @endcode
 * 
 * In the function that sets up the usual finite element data structures, we
 * first need to distribute the DoFs.
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * We start by generating the sparsity pattern. To this end, we first fill
 * an intermediate object of type DynamicSparsityPattern with the couplings
 * appearing in the system. After building the pattern, this object is
 * copied to <code>sparsity_pattern</code> and can be discarded.
 * 

 * 
 * To build the sparsity pattern for DG discretizations, we can call the
 * function analogue to DoFTools::make_sparsity_pattern, which is called
 * DoFTools::make_flux_sparsity_pattern:
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
 *     sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * Finally, we set up the structure of all components of the linear system.
 * 
 * @code
 *     system_matrix.reinit(sparsity_pattern);
 *     solution.reinit(dof_handler.n_dofs());
 *     right_hand_side.reinit(dof_handler.n_dofs());
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Theassemble_systemfunction"></a> 
 * <h4>The assemble_system function</h4>
 * 

 * 
 * Here we see the major difference to assembling by hand. Instead of writing
 * loops over cells and faces, we leave all this to the MeshWorker framework.
 * In order to do so, we just have to define local integration functions and
 * use one of the classes in namespace MeshWorker::Assembler to build the
 * global system.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::assemble_system()
 *   {
 * @endcode
 * 
 * This is the magic object, which knows everything about the data
 * structures and local integration.  This is the object doing the work in
 * the function MeshWorker::loop(), which is implicitly called by
 * MeshWorker::integration_loop() below. After the functions to which we
 * provide pointers did the local integration, the
 * MeshWorker::Assembler::SystemSimple object distributes these into the
 * global sparse matrix and the right hand side vector.
 * 
 * @code
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 * 
 * @endcode
 * 
 * First, we initialize the quadrature formulae and the update flags in the
 * worker base class. For quadrature, we play safe and use a QGauss formula
 * with number of points one higher than the polynomial degree used. Since
 * the quadratures for cells, boundary and interior faces can be selected
 * independently, we have to hand over this value three times.
 * 
 * @code
 *     const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
 *     info_box.initialize_gauss_quadrature(n_gauss_points,
 *                                          n_gauss_points,
 *                                          n_gauss_points);
 * 
 * @endcode
 * 
 * These are the types of values we need for integrating our system. They
 * are added to the flags used on cells, boundary and interior faces, as
 * well as interior neighbor faces, which is forced by the four @p true
 * values.
 * 
 * @code
 *     info_box.initialize_update_flags();
 *     UpdateFlags update_flags =
 *       update_quadrature_points | update_values | update_gradients;
 *     info_box.add_update_flags(update_flags, true, true, true, true);
 * 
 * @endcode
 * 
 * After preparing all data in <tt>info_box</tt>, we initialize the FEValues
 * objects in there.
 * 
 * @code
 *     info_box.initialize(fe, mapping);
 * 
 * @endcode
 * 
 * The object created so far helps us do the local integration on each cell
 * and face. Now, we need an object which receives the integrated (local)
 * data and forwards them to the assembler.
 * 
 * @code
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 * @endcode
 * 
 * Now, we have to create the assembler object and tell it, where to put the
 * local data. These will be our system matrix and the right hand side.
 * 
 * @code
 *     MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>>
 *       assembler;
 *     assembler.initialize(system_matrix, right_hand_side);
 * 
 * @endcode
 * 
 * Finally, the integration loop over all active cells (determined by the
 * first argument, which is an active iterator).
 *     

 * 
 * As noted in the discussion when declaring the local integration functions
 * in the class declaration, the arguments expected by the assembling
 * integrator class are not actually function pointers. Rather, they are
 * objects that can be called like functions with a certain number of
 * arguments. Consequently, we could also pass objects with appropriate
 * operator() implementations here, or lambda functions if the local
 * integrators were, for example, non-static member functions.
 * 
 * @code
 *     MeshWorker::loop<dim,
 *                      dim,
 *                      MeshWorker::DoFInfo<dim>,
 *                      MeshWorker::IntegrationInfoBox<dim>>(
 *       dof_handler.begin_active(),
 *       dof_handler.end(),
 *       dof_info,
 *       info_box,
 *       &AdvectionProblem<dim>::integrate_cell_term,
 *       &AdvectionProblem<dim>::integrate_boundary_term,
 *       &AdvectionProblem<dim>::integrate_face_term,
 *       assembler);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thelocalintegrators"></a> 
 * <h4>The local integrators</h4>
 * 

 * 
 * These are the functions given to the MeshWorker::integration_loop() called
 * just above. They compute the local contributions to the system matrix and
 * right hand side on cells and faces.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::integrate_cell_term(DoFInfo & dinfo,
 *                                                   CellInfo &info)
 *   {
 * @endcode
 * 
 * First, let us retrieve some of the objects used here from @p info. Note
 * that these objects can handle much more complex structures, thus the
 * access here looks more complicated than might seem necessary.
 * 
 * @code
 *     const FEValuesBase<dim> &  fe_values    = info.fe_values();
 *     FullMatrix<double> &       local_matrix = dinfo.matrix(0).matrix;
 *     const std::vector<double> &JxW          = fe_values.get_JxW_values();
 * 
 * @endcode
 * 
 * With these objects, we continue local integration like always. First, we
 * loop over the quadrature points and compute the advection vector in the
 * current point.
 * 
 * @code
 *     for (unsigned int point = 0; point < fe_values.n_quadrature_points; ++point)
 *       {
 *         const Tensor<1, dim> beta_at_q_point =
 *           beta(fe_values.quadrature_point(point));
 * 
 * @endcode
 * 
 * We solve a homogeneous equation, thus no right hand side shows up in
 * the cell term.  What's left is integrating the matrix entries.
 * 
 * @code
 *         for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < fe_values.dofs_per_cell; ++j)
 *             local_matrix(i, j) += -beta_at_q_point *                
 *                                   fe_values.shape_grad(i, point) *  
 *                                   fe_values.shape_value(j, point) * 
 *                                   JxW[point];
 *       }
 *   }
 * 
 * @endcode
 * 
 * Now the same for the boundary terms. Note that now we use FEValuesBase, the
 * base class for both FEFaceValues and FESubfaceValues, in order to get
 * access to normal vectors.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::integrate_boundary_term(DoFInfo & dinfo,
 *                                                       CellInfo &info)
 *   {
 *     const FEValuesBase<dim> &fe_face_values = info.fe_values();
 *     FullMatrix<double> &     local_matrix   = dinfo.matrix(0).matrix;
 *     Vector<double> &         local_vector   = dinfo.vector(0).block(0);
 * 
 *     const std::vector<double> &        JxW = fe_face_values.get_JxW_values();
 *     const std::vector<Tensor<1, dim>> &normals =
 *       fe_face_values.get_normal_vectors();
 * 
 *     std::vector<double> g(fe_face_values.n_quadrature_points);
 * 
 *     static BoundaryValues<dim> boundary_function;
 *     boundary_function.value_list(fe_face_values.get_quadrature_points(), g);
 * 
 *     for (unsigned int point = 0; point < fe_face_values.n_quadrature_points;
 *          ++point)
 *       {
 *         const double beta_dot_n =
 *           beta(fe_face_values.quadrature_point(point)) * normals[point];
 *         if (beta_dot_n > 0)
 *           for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < fe_face_values.dofs_per_cell; ++j)
 *               local_matrix(i, j) += beta_dot_n *                           
 *                                     fe_face_values.shape_value(j, point) * 
 *                                     fe_face_values.shape_value(i, point) * 
 *                                     JxW[point];
 *         else
 *           for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i)
 *             local_vector(i) += -beta_dot_n *                          
 *                                g[point] *                             
 *                                fe_face_values.shape_value(i, point) * 
 *                                JxW[point];
 *       }
 *   }
 * 
 * @endcode
 * 
 * Finally, the interior face terms. The difference here is that we receive
 * two info objects, one for each cell adjacent to the face and we assemble
 * four matrices, one for each cell and two for coupling back and forth.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::integrate_face_term(DoFInfo & dinfo1,
 *                                                   DoFInfo & dinfo2,
 *                                                   CellInfo &info1,
 *                                                   CellInfo &info2)
 *   {
 * @endcode
 * 
 * For quadrature points, weights, etc., we use the FEValuesBase object of
 * the first argument.
 * 
 * @code
 *     const FEValuesBase<dim> &fe_face_values = info1.fe_values();
 *     const unsigned int       dofs_per_cell  = fe_face_values.dofs_per_cell;
 * 
 * @endcode
 * 
 * For additional shape functions, we have to ask the neighbors
 * FEValuesBase.
 * 
 * @code
 *     const FEValuesBase<dim> &fe_face_values_neighbor = info2.fe_values();
 *     const unsigned int       neighbor_dofs_per_cell =
 *       fe_face_values_neighbor.dofs_per_cell;
 * 
 * @endcode
 * 
 * Then we get references to the four local matrices. The letters u and v
 * refer to trial and test functions, respectively. The %numbers indicate
 * the cells provided by info1 and info2. By convention, the two matrices
 * in each info object refer to the test functions on the respective cell.
 * The first matrix contains the interior couplings of that cell, while the
 * second contains the couplings between cells.
 * 
 * @code
 *     FullMatrix<double> &u1_v1_matrix = dinfo1.matrix(0, false).matrix;
 *     FullMatrix<double> &u2_v1_matrix = dinfo1.matrix(0, true).matrix;
 *     FullMatrix<double> &u1_v2_matrix = dinfo2.matrix(0, true).matrix;
 *     FullMatrix<double> &u2_v2_matrix = dinfo2.matrix(0, false).matrix;
 * 
 * @endcode
 * 
 * Here, following the previous functions, we would have the local right
 * hand side vectors. Fortunately, the interface terms only involve the
 * solution and the right hand side does not receive any contributions.
 * 

 * 
 * 
 * @code
 *     const std::vector<double> &        JxW = fe_face_values.get_JxW_values();
 *     const std::vector<Tensor<1, dim>> &normals =
 *       fe_face_values.get_normal_vectors();
 * 
 *     for (unsigned int point = 0; point < fe_face_values.n_quadrature_points;
 *          ++point)
 *       {
 *         const double beta_dot_n =
 *           beta(fe_face_values.quadrature_point(point)) * normals[point];
 *         if (beta_dot_n > 0)
 *           {
 * @endcode
 * 
 * This term we've already seen:
 * 
 * @code
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 u1_v1_matrix(i, j) += beta_dot_n *                           
 *                                       fe_face_values.shape_value(j, point) * 
 *                                       fe_face_values.shape_value(i, point) * 
 *                                       JxW[point];
 * 
 * @endcode
 * 
 * We additionally assemble the term $(\beta\cdot n u,\hat
 * v)_{\partial \kappa_+}$,
 * 
 * @code
 *             for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k)
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 u1_v2_matrix(k, j) +=
 *                   -beta_dot_n *                                   
 *                   fe_face_values.shape_value(j, point) *          
 *                   fe_face_values_neighbor.shape_value(k, point) * 
 *                   JxW[point];
 *           }
 *         else
 *           {
 * @endcode
 * 
 * This one we've already seen, too:
 * 
 * @code
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l)
 *                 u2_v1_matrix(i, l) +=
 *                   beta_dot_n *                                    
 *                   fe_face_values_neighbor.shape_value(l, point) * 
 *                   fe_face_values.shape_value(i, point) *          
 *                   JxW[point];
 * 
 * @endcode
 * 
 * And this is another new one: $(\beta\cdot n \hat u,\hat
 * v)_{\partial \kappa_-}$:
 * 
 * @code
 *             for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k)
 *               for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l)
 *                 u2_v2_matrix(k, l) +=
 *                   -beta_dot_n *                                   
 *                   fe_face_values_neighbor.shape_value(l, point) * 
 *                   fe_face_values_neighbor.shape_value(k, point) * 
 *                   JxW[point];
 *           }
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Alltherest"></a> 
 * <h3>All the rest</h3>
 *   

 * 
 * For this simple problem we use the simplest possible solver, called
 * Richardson iteration, that represents a simple defect correction. This, in
 * combination with a block SSOR preconditioner, that uses the special block
 * matrix structure of system matrices arising from DG discretizations. The
 * size of these blocks are the number of DoFs per cell. Here, we use a SSOR
 * preconditioning as we have not renumbered the DoFs according to the flow
 * field. If the DoFs are renumbered in the downstream direction of the flow,
 * then a block Gauss-Seidel preconditioner (see the PreconditionBlockSOR
 * class with relaxation=1) does a much better job.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::solve(Vector<double> &solution)
 *   {
 *     SolverControl                    solver_control(1000, 1e-12);
 *     SolverRichardson<Vector<double>> solver(solver_control);
 * 
 * @endcode
 * 
 * Here we create the preconditioner,
 * 
 * @code
 *     PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
 * 
 * @endcode
 * 
 * then assign the matrix to it and set the right block size:
 * 
 * @code
 *     preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());
 * 
 * @endcode
 * 
 * After these preparations we are ready to start the linear solver.
 * 
 * @code
 *     solver.solve(system_matrix, solution, right_hand_side, preconditioner);
 *   }
 * 
 * 
 * @endcode
 * 
 * We refine the grid according to a very simple refinement criterion, namely
 * an approximation to the gradient of the solution. As here we consider the
 * DG(1) method (i.e. we use piecewise bilinear shape functions) we could
 * simply compute the gradients on each cell. But we do not want to base our
 * refinement indicator on the gradients on each cell only, but want to base
 * them also on jumps of the discontinuous solution function over faces
 * between neighboring cells. The simplest way of doing that is to compute
 * approximative gradients by difference quotients including the cell under
 * consideration and its neighbors. This is done by the
 * <code>DerivativeApproximation</code> class that computes the approximate
 * gradients in a way similar to the <code>GradientEstimation</code> described
 * in step-9 of this tutorial. In fact, the
 * <code>DerivativeApproximation</code> class was developed following the
 * <code>GradientEstimation</code> class of step-9. Relating to the discussion
 * in step-9, here we consider $h^{1+d/2}|\nabla_h u_h|$. Furthermore we note
 * that we do not consider approximate second derivatives because solutions to
 * the linear advection equation are in general not in $H^2$ but only in $H^1$
 * (or, to be more precise: in $H^1_\beta$, i.e., the space of functions whose
 * derivatives in direction $\beta$ are square integrable).
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::refine_grid()
 *   {
 * @endcode
 * 
 * The <code>DerivativeApproximation</code> class computes the gradients to
 * float precision. This is sufficient as they are approximate and serve as
 * refinement indicators only.
 * 
 * @code
 *     Vector<float> gradient_indicator(triangulation.n_active_cells());
 * 
 * @endcode
 * 
 * Now the approximate gradients are computed
 * 
 * @code
 *     DerivativeApproximation::approximate_gradient(mapping,
 *                                                   dof_handler,
 *                                                   solution,
 *                                                   gradient_indicator);
 * 
 * @endcode
 * 
 * and they are cell-wise scaled by the factor $h^{1+d/2}$
 * 
 * @code
 *     unsigned int cell_no = 0;
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       gradient_indicator(cell_no++) *=
 *         std::pow(cell->diameter(), 1 + 1.0 * dim / 2);
 * 
 * @endcode
 * 
 * Finally they serve as refinement indicator.
 * 
 * @code
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     gradient_indicator,
 *                                                     0.3,
 *                                                     0.1);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * The output of this program consists of eps-files of the adaptively refined
 * grids and the numerical solutions given in gnuplot format.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 * @endcode
 * 
 * First write the grid in eps format.
 * 
 * @code
 *     {
 *       const std::string filename = "grid-" + std::to_string(cycle) + ".eps";
 *       deallog << "Writing grid to <" << filename << ">" << std::endl;
 *       std::ofstream eps_output(filename);
 * 
 *       GridOut grid_out;
 *       grid_out.write_eps(triangulation, eps_output);
 *     }
 * 
 * @endcode
 * 
 * Then output the solution in gnuplot format.
 * 
 * @code
 *     {
 *       const std::string filename = "sol-" + std::to_string(cycle) + ".gnuplot";
 *       deallog << "Writing solution to <" << filename << ">" << std::endl;
 *       std::ofstream gnuplot_output(filename);
 * 
 *       DataOut<dim> data_out;
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "u");
 * 
 *       data_out.build_patches();
 * 
 *       data_out.write_gnuplot(gnuplot_output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * The following <code>run</code> function is similar to previous examples.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle)
 *       {
 *         deallog << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation);
 * 
 *             triangulation.refine_global(3);
 *           }
 *         else
 *           refine_grid();
 * 
 * 
 *         deallog << "Number of active cells:       "
 *                 << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *                 << std::endl;
 * 
 *         assemble_system();
 *         solve(solution);
 * 
 *         output_results(cycle);
 *       }
 *   }
 * } // namespace Step12
 * 
 * 
 * @endcode
 * 
 * The following <code>main</code> function is similar to previous examples as
 * well, and need not be commented on.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       dealii::deallog.depth_console(5);
 * 
 *       Step12::AdvectionProblem<2> dgmethod;
 *       dgmethod.run();
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
 * 
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>



The output of this program is very similar to step-16 and we are not repeating the output here.

We show the solutions on the initial mesh, the mesh after two
and after five adaptive refinement steps.

<img src="https://www.dealii.org/images/steps/developer/step-12.sol-0.png" alt="">
<img src="https://www.dealii.org/images/steps/developer/step-12.sol-2.png" alt="">
<img src="https://www.dealii.org/images/steps/developer/step-12.sol-5.png" alt="">


Then we show the final grid (after 5 refinement steps) and the solution again,
this time with a nicer 3d rendering (obtained using the DataOutBase::write_vtk
function and the VTK-based VisIt visualization program) that better shows the
sharpness of the jump on the refined mesh and the over- and undershoots of the
solution along the interface:

<img src="https://www.dealii.org/images/steps/developer/step-12.grid-5.png" alt="">
<img src="https://www.dealii.org/images/steps/developer/step-12.3d-solution.png" alt="">


And finally we show a plot of a 3d computation.

<img src="https://www.dealii.org/images/steps/developer/step-12.sol-5-3d.png" alt="">


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


For ideas for further extensions, please see see step-12.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-12b.cc"
*/
