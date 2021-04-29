/**
@page step_11 The step-11 tutorial program
This tutorial depends on step-10.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


The problem we will be considering is the solution of Laplace's problem with
Neumann boundary conditions only:
@f{eqnarray*}
  -\Delta u &=& f \qquad \mathrm{in}\ \Omega,
  \\
  \partial_n u &=& g \qquad \mathrm{on}\ \partial\Omega.
@f}
It is well known that if this problem is to have a solution, then the forces
need to satisfy the compatibility condition
@f[
  \int_\Omega f\; dx + \int_{\partial\Omega} g\; ds = 0.
@f]
We will consider the special case that $\Omega$ is the circle of radius 1
around the origin, and $f=-2$, $g=1$. This choice satisfies the compatibility
condition.

The compatibility condition allows a solution of the above equation, but it
nevertheless retains an ambiguity: since only derivatives of the solution
appear in the equations, the solution is only determined up to a constant. For
this reason, we have to pose another condition for the numerical solution,
which fixes this constant.

For this, there are various possibilities:
<ol>
<li> Fix one node of the discretization to zero or any other fixed value.
  This amounts to an additional condition $u_h(x_0)=0$. Although this is
  common practice, it is not necessarily a good idea, since we know that the
  solutions of Laplace's equation are only in $H^1$, which does not allow for
  the definition of point values because it is not a subset of the continuous
  functions. Therefore, even though fixing one node is allowed for
  discretized functions, it is not for continuous functions, and one can
  often see this in a resulting error spike at this point in the numerical
  solution.

<li> Fixing the mean value over the domain to zero or any other value. This
  is allowed on the continuous level, since $H^1(\Omega)\subset L^1(\Omega)$
  by Sobolev's inequality, and thus also on the discrete level since we
  there only consider subsets of $H^1$.

<li> Fixing the mean value over the boundary of the domain to zero or any
  other value. This is also allowed on the continuous level, since
  $H^{1/2}(\partial\Omega)\subset L^1(\partial\Omega)$, again by Sobolev's
  inequality.
</ol>
We will choose the last possibility, since we want to demonstrate another
technique with it.

While this describes the problem to be solved, we still have to figure out how
to implement it. Basically, except for the additional mean value constraint,
we have solved this problem several times, using Dirichlet boundary values,
and we only need to drop the treatment of Dirichlet boundary nodes. The use of
higher order mappings is also rather trivial and will be explained at the
various places where we use it; in almost all conceivable cases, you will only
consider the objects describing mappings as a black box which you need not
worry about, because their only uses seem to be to be passed to places deep
inside the library where functions know how to handle them (i.e. in the
<code>FEValues</code> classes and their descendants).

The tricky point in this program is the use of the mean value
constraint. Fortunately, there is a class in the library which knows how to
handle such constraints, and we have used it quite often already, without
mentioning its generality. Note that if we assume that the boundary nodes are
spaced equally along the boundary, then the mean value constraint
@f[
  \int_{\partial \Omega} u(x) \; ds = 0
@f]
can be written as
@f[
  \sum_{i\in\partial\Omega_h} u_i = 0,
@f]
where the sum shall run over all degree of freedom indices which are located
on the boundary of the computational domain. Let us denote by $i_0$ that index
on the boundary with the lowest number (or any other conveniently chosen
index), then the constraint can also be represented by
@f[
  u_{i_0} = \sum_{i\in\partial\Omega_h\backslash i_0} -u_i.
@f]
This, luckily, is exactly the form of constraints for which the
AffineConstraints class was designed. Note that we have used this
class in several previous examples for the representation of hanging nodes
constraints, which also have this form: there, the middle vertex shall have
the mean of the values of the adjacent vertices. In general, the
AffineConstraints class is designed to handle affine constraints
of the form
@f[
  CU = b
@f]
where $C$ denotes a matrix, $b$ denotes a vector, and $U$ the vector of nodal
values. In this case, since $C$ represents one homogeneous constraint, $b$ is
the zero vector.

In this example, the mean value along the boundary allows just such a
representation, with $C$ being a matrix with just one row (i.e. there is only
one constraint). In the implementation, we will create an AffineConstraints
object, add one constraint (i.e. add another row to the matrix) referring to the
first boundary node $i_0$, and insert the weights with which all the other nodes
contribute, which in this example happens to be just $-1$.

Later, we will use this object to eliminate the first boundary node from the
linear system of equations, reducing it to one which has a solution without
the ambiguity of the constant shift value. One of the problems of the
implementation will be that the explicit elimination of this node results in a
number of additional elements in the matrix, of which we do not know in
advance where they are located and how many additional entries will be in each
of the rows of the matrix. We will show how we can use an intermediate object
to work around this problem.

But now on to the implementation of the program solving this problem...
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * As usual, the program starts with a rather long list of include files which
 * you are probably already used to by now:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/table_handler.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * Just this one is new: it declares a class
 * DynamicSparsityPattern, which we will use and explain
 * further down below.
 * 
 * @code
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * @endcode
 * 
 * We will make use of the std::find algorithm of the C++ standard library, so
 * we have to include the following file for its declaration:
 * 
 * @code
 * #include <algorithm>
 * #include <iostream>
 * #include <iomanip>
 * #include <cmath>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step11
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Then we declare a class which represents the solution of a Laplace
 * problem. As this example program is based on step-5, the class looks
 * rather the same, with the sole structural difference that the functions
 * <code>assemble_system</code> now calls <code>solve</code> itself, and is
 * thus called <code>assemble_and_solve</code>, and that the output function
 * was dropped since the solution function is so boring that it is not worth
 * being viewed.
 *   

 * 
 * The only other noteworthy change is that the constructor takes a value
 * representing the polynomial degree of the mapping to be used later on,
 * and that it has another member variable representing exactly this
 * mapping. In general, this variable will occur in real applications at the
 * same places where the finite element is declared or used.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem(const unsigned int mapping_degree);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_and_solve();
 *     void solve();
 *     void write_high_order_mesh(const unsigned cycle);
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 *     MappingQ<dim>      mapping;
 * 
 *     SparsityPattern           sparsity_pattern;
 *     SparseMatrix<double>      system_matrix;
 *     AffineConstraints<double> mean_value_constraints;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     TableHandler output_table;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * Construct such an object, by initializing the variables. Here, we use
 * linear finite elements (the argument to the <code>fe</code> variable
 * denotes the polynomial degree), and mappings of given order. Print to
 * screen what we are about to do.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem(const unsigned int mapping_degree)
 *     : fe(1)
 *     , dof_handler(triangulation)
 *     , mapping(mapping_degree)
 *   {
 *     std::cout << "Using mapping with degree " << mapping_degree << ":"
 *               << std::endl
 *               << "============================" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The first task is to set up the variables for this problem. This includes
 * generating a valid <code>DoFHandler</code> object, as well as the
 * sparsity patterns for the matrix, and the object representing the
 * constraints that the mean value of the degrees of freedom on the boundary
 * be zero.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 * @endcode
 * 
 * The first task is trivial: generate an enumeration of the degrees of
 * freedom, and initialize solution and right hand side vector to their
 * correct sizes:
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * The next task is to construct the object representing the constraint that
 * the mean value of the degrees of freedom on the boundary shall be
 * zero. For this, we first want a list of those nodes that are actually
 * at the boundary. The <code>DoFTools</code> namespace has a function
 * that returns an IndexSet object that contains the indices of all those
 * degrees of freedom that are at the boundary.
 *     

 * 
 * Once we have this index set, we wanted to know which is the first
 * index corresponding to a degree of freedom on the boundary. We need
 * this because we wanted to constrain one of the nodes on the boundary by
 * the values of all other DoFs on the boundary. To get the index of this
 * "first" degree of freedom is easy enough using the IndexSet class:
 * 
 * @code
 *     const IndexSet boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler);
 * 
 *     const types::global_dof_index first_boundary_dof =
 *       boundary_dofs.nth_index_in_set(0);
 * 
 * @endcode
 * 
 * Then generate a constraints object with just this one constraint. First
 * clear all previous content (which might reside there from the previous
 * computation on a once coarser grid), then add this one line
 * constraining the <code>first_boundary_dof</code> to the sum of other
 * boundary DoFs each with weight -1. Finally, close the constraints
 * object, i.e. do some internal bookkeeping on it for faster processing
 * of what is to come later:
 * 
 * @code
 *     mean_value_constraints.clear();
 *     mean_value_constraints.add_line(first_boundary_dof);
 *     for (types::global_dof_index i : boundary_dofs)
 *       if (i != first_boundary_dof)
 *         mean_value_constraints.add_entry(first_boundary_dof, i, -1);
 *     mean_value_constraints.close();
 * 
 * @endcode
 * 
 * Next task is to generate a sparsity pattern. This is indeed a tricky
 * task here. Usually, we just call
 * <code>DoFTools::make_sparsity_pattern</code> and condense the result
 * using the hanging node constraints. We have no hanging node constraints
 * here (since we only refine globally in this example), but we have this
 * global constraint on the boundary. This poses one severe problem in
 * this context: the <code>SparsityPattern</code> class wants us to state
 * beforehand the maximal number of entries per row, either for all rows
 * or for each row separately. There are functions in the library which
 * can tell you this number in case you just have hanging node constraints
 * (namely DoFHandler::max_couplings_between_dofs), but how is
 * this for the present case? The difficulty arises because the
 * elimination of the constrained degree of freedom requires a number of
 * additional entries in the matrix at places that are not so simple to
 * determine. We would therefore have a problem had we to give a maximal
 * number of entries per row here.
 *     

 * 
 * Since this can be so difficult that no reasonable answer can be given
 * that allows allocation of only a reasonable amount of memory, there is
 * a class DynamicSparsityPattern, that can help us out
 * here. It does not require that we know in advance how many entries rows
 * could have, but allows just about any length. It is thus significantly
 * more flexible in case you do not have good estimates of row lengths,
 * however at the price that building up such a pattern is also
 * significantly more expensive than building up a pattern for which you
 * had information in advance. Nevertheless, as we have no other choice
 * here, we'll just build such an object by initializing it with the
 * dimensions of the matrix and calling another function
 * <code>DoFTools::make_sparsity_pattern</code> to get the sparsity
 * pattern due to the differential operator, then condense it with the
 * constraints object which adds those positions in the sparsity pattern
 * that are required for the elimination of the constraint.
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *     mean_value_constraints.condense(dsp);
 * 
 * @endcode
 * 
 * Finally, once we have the full pattern, we can initialize an object of
 * type <code>SparsityPattern</code> from it and in turn initialize the
 * matrix with it. Note that this is actually necessary, since the
 * DynamicSparsityPattern is so inefficient compared to
 * the <code>SparsityPattern</code> class due to the more flexible data
 * structures it has to use, that we can impossibly base the sparse matrix
 * class on it, but rather need an object of type
 * <code>SparsityPattern</code>, which we generate by copying from the
 * intermediate object.
 *     

 * 
 * As a further sidenote, you will notice that we do not explicitly have
 * to <code>compress</code> the sparsity pattern here. This, of course, is
 * due to the fact that the <code>copy_from</code> function generates a
 * compressed object right from the start, to which you cannot add new
 * entries anymore. The <code>compress</code> call is therefore implicit
 * in the <code>copy_from</code> call.
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function then assembles the linear system of equations, solves
 * it, and evaluates the solution. This then makes three actions, and we
 * will put them into eight true statements (excluding declaration of
 * variables, and handling of temporary vectors). Thus, this function is
 * something for the very lazy. Nevertheless, the functions called are
 * rather powerful, and through them this function uses a good deal of the
 * whole library. But let's look at each of the steps.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_and_solve()
 *   {
 * @endcode
 * 
 * First, we have to assemble the matrix and the right hand side. In all
 * previous examples, we have investigated various ways how to do this
 * manually. However, since the Laplace matrix and simple right hand sides
 * appear so frequently in applications, the library provides functions
 * for actually doing this for you, i.e. they perform the loop over all
 * cells, setting up the local matrices and vectors, and putting them
 * together for the end result.
 *     

 * 
 * The following are the two most commonly used ones: creation of the
 * Laplace matrix and creation of a right hand side vector from body or
 * boundary forces. They take the mapping object, the
 * <code>DoFHandler</code> object representing the degrees of freedom and
 * the finite element in use, a quadrature formula to be used, and the
 * output object. The function that creates a right hand side vector also
 * has to take a function object describing the (continuous) right hand
 * side function.
 *     

 * 
 * Let us look at the way the matrix and body forces are integrated:
 * 
 * @code
 *     const unsigned int gauss_degree =
 *       std::max(static_cast<unsigned int>(
 *                  std::ceil(1. * (mapping.get_degree() + 1) / 2)),
 *                2U);
 *     MatrixTools::create_laplace_matrix(mapping,
 *                                        dof_handler,
 *                                        QGauss<dim>(gauss_degree),
 *                                        system_matrix);
 *     VectorTools::create_right_hand_side(mapping,
 *                                         dof_handler,
 *                                         QGauss<dim>(gauss_degree),
 *                                         Functions::ConstantFunction<dim>(-2),
 *                                         system_rhs);
 * @endcode
 * 
 * That's quite simple, right?
 *     

 * 
 * Two remarks are in order, though: First, these functions are used in a
 * lot of contexts. Maybe you want to create a Laplace or mass matrix for
 * a vector values finite element; or you want to use the default Q1
 * mapping; or you want to assembled the matrix with a coefficient in the
 * Laplace operator. For this reason, there are quite a large number of
 * variants of these functions in the <code>MatrixCreator</code> and
 * <code>MatrixTools</code> namespaces. Whenever you need a slightly
 * different version of these functions than the ones called above, it is
 * certainly worthwhile to take a look at the documentation and to check
 * whether something fits your needs.
 *     

 * 
 * The second remark concerns the quadrature formula we use: we want to
 * integrate over bilinear shape functions, so we know that we have to use
 * at least an order two Gauss quadrature formula. On the other hand, we
 * want the quadrature rule to have at least the order of the boundary
 * approximation. Since the order of Gauss rule with $r$ points is $2r -
 * 1$, and the order of the boundary approximation using polynomials of
 * degree $p$ is $p+1$, we know that $2r \geq p$. Since r has to be an
 * integer and (as mentioned above) has to be at least $2$, this makes up
 * for the formula above computing <code>gauss_degree</code>.
 *     

 * 
 * Since the generation of the body force contributions to the right hand
 * side vector was so simple, we do that all over again for the boundary
 * forces as well: allocate a vector of the right size and call the right
 * function. The boundary function has constant values, so we can generate
 * an object from the library on the fly, and we use the same quadrature
 * formula as above, but this time of lower dimension since we integrate
 * over faces now instead of cells:
 * 
 * @code
 *     Vector<double> tmp(system_rhs.size());
 *     VectorTools::create_boundary_right_hand_side(
 *       mapping,
 *       dof_handler,
 *       QGauss<dim - 1>(gauss_degree),
 *       Functions::ConstantFunction<dim>(1),
 *       tmp);
 * @endcode
 * 
 * Then add the contributions from the boundary to those from the interior
 * of the domain:
 * 
 * @code
 *     system_rhs += tmp;
 * @endcode
 * 
 * For assembling the right hand side, we had to use two different vector
 * objects, and later add them together. The reason we had to do so is
 * that the <code>VectorTools::create_right_hand_side</code> and
 * <code>VectorTools::create_boundary_right_hand_side</code> functions
 * first clear the output vector, rather than adding up their results to
 * previous contents. This can reasonably be called a design flaw in the
 * library made in its infancy, but unfortunately things are as they are
 * for some time now and it is difficult to change such things that
 * silently break existing code, so we have to live with that.
 * 

 * 
 * Now, the linear system is set up, so we can eliminate the one degree of
 * freedom which we constrained to the other DoFs on the boundary for the
 * mean value constraint from matrix and right hand side vector, and solve
 * the system. After that, distribute the constraints again, which in this
 * case means setting the constrained degree of freedom to its proper
 * value
 * 
 * @code
 *     mean_value_constraints.condense(system_matrix);
 *     mean_value_constraints.condense(system_rhs);
 * 
 *     solve();
 *     mean_value_constraints.distribute(solution);
 * 
 * @endcode
 * 
 * Finally, evaluate what we got as solution. As stated in the
 * introduction, we are interested in the H1 semi-norm of the
 * solution. Here, as well, we have a function in the library that does
 * this, although in a slightly non-obvious way: the
 * <code>VectorTools::integrate_difference</code> function integrates the
 * norm of the difference between a finite element function and a
 * continuous function. If we therefore want the norm of a finite element
 * field, we just put the continuous function to zero. Note that this
 * function, just as so many other ones in the library as well, has at
 * least two versions, one which takes a mapping as argument (which we
 * make us of here), and the one which we have used in previous examples
 * which implicitly uses <code>MappingQ1</code>.  Also note that we take a
 * quadrature formula of one degree higher, in order to avoid
 * superconvergence effects where the solution happens to be especially
 * close to the exact solution at certain points (we don't know whether
 * this might be the case here, but there are cases known of this, and we
 * just want to make sure):
 * 
 * @code
 *     Vector<float> norm_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Functions::ZeroFunction<dim>(),
 *                                       norm_per_cell,
 *                                       QGauss<dim>(gauss_degree + 1),
 *                                       VectorTools::H1_seminorm);
 * @endcode
 * 
 * Then, the function just called returns its results as a vector of
 * values each of which denotes the norm on one cell. To get the global
 * norm, we do the following:
 * 
 * @code
 *     const double norm =
 *       VectorTools::compute_global_error(triangulation,
 *                                         norm_per_cell,
 *                                         VectorTools::H1_seminorm);
 * 
 * @endcode
 * 
 * Last task -- generate output:
 * 
 * @code
 *     output_table.add_value("cells", triangulation.n_active_cells());
 *     output_table.add_value("|u|_1", norm);
 *     output_table.add_value("error",
 *                            std::fabs(norm - std::sqrt(3.14159265358 / 2)));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The following function solving the linear system of equations is copied
 * from step-5 and is explained there in some detail:
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Next, we write the solution as well as the
 * material ids to a VTU file. This is similar to what was done in many
 * other tutorial programs. The new ingredient presented in this tutorial
 * program is that we want to ensure that the data written to the file
 * used for visualization is actually a faithful representation of what
 * is used internally by deal.II. That is because most of the visualization
 * data formats only represent cells by their vertex coordinates, but
 * have no way of representing the curved boundaries that are used
 * in deal.II when using higher order mappings -- in other words, what
 * you see in the visualization tool is not actually what you are computing
 * on. (The same, incidentally, is true when using higher order shape
 * functions: Most visualization tools only render bilinear/trilinear
 * representations. This is discussed in detail in DataOut::build_patches().)
 *   

 * 
 * So we need to ensure that a high-order representation is written
 * to the file. We need to consider two particular topics. Firstly, we tell
 * the DataOut object via the DataOutBase::VtkFlags that we intend to
 * interpret the subdivisions of the elements as a high-order Lagrange
 * polynomial rather than a collection of bilinear patches.
 * Recent visualization programs, like ParaView version 5.5
 * or newer, can then render a high-order solution (see a <a
 * href="https://github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
 * page</a> for more details). Secondly, we need to make sure that the mapping
 * is passed to the DataOut::build_patches() method. Finally, the DataOut
 * class only prints curved faces for <i>boundary</i> cells by default, so we
 * need to ensure that also inner cells are printed in a curved representation
 * via the mapping.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::write_high_order_mesh(const unsigned cycle)
 *   {
 *     DataOut<dim> data_out;
 * 
 *     DataOutBase::VtkFlags flags;
 *     flags.write_higher_order_cells = true;
 *     data_out.set_flags(flags);
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 * 
 *     data_out.build_patches(mapping,
 *                            mapping.get_degree(),
 *                            DataOut<dim>::curved_inner_cells);
 * 
 *     std::ofstream file("solution-c=" + std::to_string(cycle) +
 *                        ".p=" + std::to_string(mapping.get_degree()) + ".vtu");
 * 
 *     data_out.write_vtu(file);
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally the main function controlling the different steps to be
 * performed. Its content is rather straightforward, generating a
 * triangulation of a circle, associating a boundary to it, and then doing
 * several cycles on subsequently finer grids. Note again that we have put
 * mesh refinement into the loop header; this may be something for a test
 * program, but for real applications you should consider that this implies
 * that the mesh is refined after the loop is executed the last time since
 * the increment clause (the last part of the three-parted loop header) is
 * executed before the comparison part (the second one), which may be rather
 * costly if the mesh is already quite refined. In that case, you should
 * arrange code such that the mesh is not further refined after the last
 * loop run (or you should do it at the beginning of each run except for the
 * first one).
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_ball(triangulation);
 * 
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle)
 *       {
 *         setup_system();
 *         assemble_and_solve();
 *         write_high_order_mesh(cycle);
 * 
 *         triangulation.refine_global();
 *       }
 * 
 * @endcode
 * 
 * After all the data is generated, write a table of results to the
 * screen:
 * 
 * @code
 *     output_table.set_precision("|u|_1", 6);
 *     output_table.set_precision("error", 6);
 *     output_table.write_text(std::cout);
 *     std::cout << std::endl;
 *   }
 * } // namespace Step11
 * 
 * 
 * 
 * @endcode
 * 
 * Finally the main function. It's structure is the same as that used in
 * several of the previous examples, so probably needs no more explanation.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       std::cout.precision(5);
 * 
 * @endcode
 * 
 * This is the main loop, doing the computations with mappings of linear
 * through cubic mappings. Note that since we need the object of type
 * <code>LaplaceProblem@<2@></code> only once, we do not even name it,
 * but create an unnamed such object and call the <code>run</code>
 * function of it, subsequent to which it is immediately destroyed
 * again.
 * 
 * @code
 *       for (unsigned int mapping_degree = 1; mapping_degree <= 3;
 *            ++mapping_degree)
 *         Step11::LaplaceProblem<2>(mapping_degree).run();
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


This is what the program outputs:
@code
Using mapping with degree 1:
============================
cells  |u|_1    error
    5 0.680402 0.572912
   20 1.088141 0.165173
   80 1.210142 0.043172
  320 1.242375 0.010939
 1280 1.250569 0.002745
 5120 1.252627 0.000687

Using mapping with degree 2:
============================
cells  |u|_1    error
    5 1.177062 0.076252
   20 1.228978 0.024336
   80 1.245175 0.008139
  320 1.250881 0.002433
 1280 1.252646 0.000668
 5120 1.253139 0.000175

Using mapping with degree 3:
============================
cells  |u|_1    error
    5 1.193493 0.059821
   20 1.229825 0.023489
   80 1.245221 0.008094
  320 1.250884 0.002430
 1280 1.252646 0.000668
 5120 1.253139 0.000175
@endcode
As we expected, the convergence order for each of the different
mappings is clearly quadratic in the mesh size. What <em>is</em>
interesting, though, is that the error for a bilinear mapping
(i.e. degree 1) is more than three times larger than that for the
higher order mappings; it is therefore clearly advantageous in this
case to use a higher order mapping, not because it improves the order
of convergence but just to reduce the constant before the convergence
order. On the other hand, using a cubic mapping only improves the
result further by insignificant amounts, except on very coarse
grids.

We can also visualize the underlying meshes by using, for instance,
ParaView. The image below shows initial meshes for different mapping
degrees:

<img src="https://www.dealii.org/images/steps/developer/step-11.cycle_0.png" alt="">

Clearly, the effect is most pronounced when we go from the linear to
quadratic mapping. This is also reflected in the error values given
in the table above. The effect of going from quadratic to cubic degree
is less dramatic, but still tangible owing to a more accurate
description of the circular boundary.

Next, let's look at the meshes after three global refinements

<img src="https://www.dealii.org/images/steps/developer/step-11.cycle_3.png" alt="">

Here, the differences are much less visible, especially for higher order
mappings. Indeed, at this refinement level the error values reported
in the table are essentially identical between mappings of degrees two
and three.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-11.cc"
*/
