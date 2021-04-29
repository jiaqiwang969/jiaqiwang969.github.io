/**
@page step_48 The step-48 tutorial program
This tutorial depends on step-25, step-37.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemstatementanddiscretization"> Problem statement and discretization </a>
        <li><a href="#Implementationofconstraints">Implementation of constraints</a>
        <li><a href="#Parallelization"> Parallelization </a>
        <li><a href="#Thetestcase"> The test case </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#SineGordonOperation">SineGordonOperation</a>
      <ul>
        <li><a href="#SineGordonOperationSineGordonOperation">SineGordonOperation::SineGordonOperation</a>
        <li><a href="#SineGordonOperationlocal_apply">SineGordonOperation::local_apply</a>
        <li><a href="#SineGordonOperationapply">SineGordonOperation::apply</a>
      </ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#SineGordonProblemclass">SineGordonProblem class</a>
      <ul>
        <li><a href="#SineGordonProblemSineGordonProblem">SineGordonProblem::SineGordonProblem</a>
        <li><a href="#SineGordonProblemmake_grid_and_dofs">SineGordonProblem::make_grid_and_dofs</a>
        <li><a href="#SineGordonProblemoutput_results">SineGordonProblem::output_results</a>
        <li><a href="#SineGordonProblemrun">SineGordonProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Comparisonwithasparsematrix">Comparison with a sparse matrix</a>
        <li><a href="#Parallelrunin2Dand3D">Parallel run in 2D and 3D</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly

<i>
This program was contributed by Katharina Kormann and Martin
Kronbichler.

The algorithm for the matrix-vector product is based on the article <a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic
interface for parallel cell-based finite element operator
application</a> by Martin Kronbichler and Katharina Kormann, Computers
and Fluids 63:135&ndash;147, 2012, and the paper &quot;Parallel finite element operator
application: Graph partitioning and coloring&quot; by Katharina
Kormann and Martin Kronbichler in: Proceedings of the 7th IEEE
International Conference on e-Science, 2011.  </i>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This program demonstrates how to use the cell-based implementation of finite
element operators with the MatrixFree class, first introduced in step-37, to
solve nonlinear partial differential equations. Moreover, we have another look
at the handling of constraints within the matrix-free framework.
Finally, we will use an explicit time-stepping
method to solve the problem and introduce Gauss-Lobatto finite elements that
are very convenient in this case since their mass matrix can be accurately
approximated by a diagonal, and thus trivially invertible, matrix. The two
ingredients to this property are firstly a distribution of the nodal points of
Lagrange polynomials according to the point distribution of the Gauss-Lobatto
quadrature rule. Secondly, the quadrature is done with the same Gauss-Lobatto
quadrature rule. In this formula, the integrals $\int_K \varphi_i \varphi_j
dx\approx \sum_q \varphi_i \varphi_j \mathrm{det}(J) \big |_{x_q}$ become
zero whenever $i\neq j$, because exactly one function $\varphi_j$ is one and
all others zero in the points defining the Lagrange polynomials.
Moreover, the Gauss-Lobatto distribution of nodes of Lagrange
polynomials clusters the nodes towards the element boundaries. This results in
a well-conditioned polynomial basis for high-order discretization
methods. Indeed, the condition number of an FE_Q elements with equidistant
nodes grows exponentially with the degree, which destroys any benefit for
orders of about five and higher. For this reason, Gauss-Lobatto points are the
default distribution for the FE_Q element (but at degrees one and two, those
are equivalent to the equidistant points).

<a name="Problemstatementanddiscretization"></a><h3> Problem statement and discretization </h3>


As an example, we choose to solve the sine-Gordon soliton equation
\f{eqnarray*}
u_{tt} &=& \Delta u -\sin(u) \quad\mbox{for}\quad (x,t) \in
\Omega \times (t_0,t_f],\\
{\mathbf n} \cdot \nabla u &=& 0
\quad\mbox{for}\quad (x,t) \in \partial\Omega \times (t_0,t_f],\\
u(x,t_0) &=& u_0(x).
\f}

that was already introduced in step-25. As a simple explicit time
integration method, we choose leap frog scheme using the second-order
formulation of the equation. With this time stepping, the scheme reads in
weak form

\f{eqnarray*}
(v,u^{n+1}) = (v,2 u^n-u^{n-1} -
(\Delta t)^2 \sin(u^n)) - (\nabla v, (\Delta t)^2 \nabla u^n),
\f}
where <i> v</i> denotes a test function and the index <i>n</i> stands for
the time step number.

For the spatial discretization, we choose FE_Q elements
with basis functions defined to interpolate the support points of the
Gauss-Lobatto quadrature rule. Moreover, when we compute the integrals
over the basis functions to form the mass matrix and the operator on
the right hand side of the equation above, we use the
Gauss-Lobatto quadrature rule with the same support points as the
node points of the finite element to evaluate the integrals. Since the
finite element is Lagrangian, this will yield a diagonal mass matrix
on the left hand side of the equation, making the solution of the
linear system in each time step trivial.

Using this quadrature rule, for a <i>p</i>th order finite element, we use a
<i>(2p-1)</i>th order accurate formula to evaluate the integrals. Since the
product of two <i>p</i>th order basis functions when computing a mass matrix
gives a function with polynomial degree <i>2p</i> in each direction, the
integrals are not computed exactly.  However, the overall convergence
properties are not disturbed by the quadrature error on meshes with affine
element shapes with L2 errors proportional to <i>h<sup>p+1</sup></i>. Note
however that order reduction with sub-optimal convergence rates of the L2
error of <i>O(h<sup>p</sup>)</i> or even <i>O(h<sup>p-1</sup>)</i> for some 3D
setups has been reported <a href="https://dx.doi.org/10.1002/num.20353">in
literature</a> on deformed (non-affine) element shapes for wave equations
when the integrand is not a polynomial any more.

Apart from the fact that we avoid solving linear systems with this
type of elements when using explicit time-stepping, they come with two
other advantages. When we are using the sum-factorization approach to
evaluate the finite element operator (cf. step-37), we have to
evaluate the function at the quadrature points. In the case of
Gauss-Lobatto elements, where quadrature points and node points of the
finite element coincide, this operation is trivial since the value
of the function at the quadrature points is given by its one-dimensional
coefficients. In this way, the arithmetic work for the finite element operator
evaluation is reduced by approximately a factor of two compared to the generic
Gaussian quadrature.

To sum up the discussion, by using the right finite element and
quadrature rule combination, we end up with a scheme where we
only need to compute the right hand side vector corresponding
to the formulation above and then multiply it by the inverse of the
diagonal mass matrix in each time step. In practice, of course, we extract
the diagonal elements and invert them only once at the beginning of the
program.

<a name="Implementationofconstraints"></a><h3>Implementation of constraints</h3>


The usual way to handle constraints in <code>deal.II</code> is to use
the AffineConstraints class that builds a sparse matrix storing
information about which degrees of freedom (DoF) are constrained and
how they are constrained. This format uses an unnecessarily large
amount of memory since there are not so many different types of
constraints: for example, in the case of hanging nodes when using
linear finite element on every cell, most constraints have the form
$x_k = \frac 12 x_i + \frac 12 x_j$ where the coefficients $\frac 12$
are always the same and only $i,j,k$ are different. While storing this
redundant information is not a problem in general because it is only
needed once during matrix and right hand side assembly, it becomes a
bottleneck in the matrix-free approach since there this
information has to be accessed every time we apply the operator, and the
remaining components of the operator evaluation are so fast. Thus,
instead of an AffineConstraints object, MatrixFree uses a variable that
we call <code>constraint_pool</code> that collects the weights of the
different constraints. Then, only an identifier of each constraint in the
mesh instead of all the weights have to be stored. Moreover,
the constraints are not applied in a pre- and postprocessing step
but rather as we evaluate the finite element
operator. Therefore, the constraint information is embedded into the
variable <code>indices_local_to_global</code> that is used to extract
the cell information from the global vector. If a DoF is constrained,
the <code>indices_local_to_global</code> variable contains the global
indices of the DoFs that it is constrained to. Then, we have another
variable <code>constraint_indicator</code> at hand that holds, for
each cell, the local indices of DoFs that are constrained as well as
the identifier of the type of constraint. Fortunately, you will not see
these data structures in the example program since the class
<code>FEEvaluation</code> takes care of the constraints without user
interaction.

In the presence of hanging nodes, the diagonal mass matrix obtained on the
element level via the Gauss-Lobatto quadrature/node point procedure does not
directly translate to a diagonal global mass matrix, as following the
constraints on rows and columns would also add off-diagonal entries. As
explained in <a href="https://dx.doi.org/10.4208/cicp.101214.021015a">Kormann
(2016)</a>, interpolating constraints on a vector, which maintains the
diagonal shape of the mass matrix, is consistent with the equations up to an
error of the same magnitude as the quadrature error. In the program below, we
will simply assemble the diagonal of the mass matrix as if it were a vector to
enable this approximation.


<a name="Parallelization"></a><h3> Parallelization </h3>


The MatrixFree class comes with the option to be parallelized on three levels:
MPI parallelization on clusters of distributed nodes, thread parallelization
scheduled by the Threading Building Blocks library, and finally with a
vectorization by working on a batch of two (or more) cells via SIMD data type
(sometimes called cross-element or external vectorization).
As we have already discussed in step-37, you will
get best performance by using an instruction set specific to your system,
e.g. with the cmake variable <tt>-DCMAKE_CXX_FLAGS="-march=native"</tt>. The
MPI parallelization was already exploited in step-37. Here, we additionally
consider thread parallelization with TBB. This is fairly simple, as all we
need to do is to tell the initialization of the MatrixFree object about the
fact that we want to use a thread parallel scheme through the variable
MatrixFree::AdditionalData::thread_parallel_scheme. During setup, a dependency
graph is set up similar to the one described in the @ref workstream_paper ,
which allows to schedule the work of the @p local_apply function on chunks of
cells without several threads accessing the same vector indices. As opposed to
the WorkStream loops, some additional clever tricks to avoid global
synchronizations as described in <a
href="https://dx.doi.org/10.1109/eScience.2011.53">Kormann and Kronbichler
(2011)</a> are also applied.

Note that this program is designed to be run with a distributed triangulation
(parallel::distributed::Triangulation), which requires deal.II to be
configured with <a href="http://www.p4est.org/">p4est</a> as described
in the <a href="../../readme.html">deal.II ReadMe</a> file. However, a
non-distributed triangulation is also supported, in which case the
computation will be run in serial.

<a name="Thetestcase"></a><h3> The test case </h3>


In our example, we choose the initial value to be \f{eqnarray*} u(x,t) =
\prod_{i=1}^{d} -4 \arctan \left(
\frac{m}{\sqrt{1-m^2}}\frac{\sin\left(\sqrt{1-m^2} t +c_2\right)}{\cosh(mx_i+c_1)}\right)
\f} and solve the equation over the time interval [-10,10]. The
constants are chosen to be $c_1=c_1=0$ and <i> m=0.5</i>. As mentioned
in step-25, in one dimension <i>u</i> as a function of <i>t</i> is the exact
solution of the sine-Gordon equation. For higher dimension, this is however
not the case.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The necessary files from the deal.II library.
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/distributed/tria.h>
 * 
 * @endcode
 * 
 * This includes the data structures for the efficient implementation of
 * matrix-free methods.
 * 
 * @code
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * #include <iomanip>
 * 
 * 
 * namespace Step48
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * We start by defining two global variables to collect all parameters
 * subject to changes at one place: One for the dimension and one for the
 * finite element degree. The dimension is used in the main function as a
 * template argument for the actual classes (like in all other deal.II
 * programs), whereas the degree of the finite element is more crucial, as
 * it is passed as a template argument to the implementation of the
 * Sine-Gordon operator. Therefore, it needs to be a compile-time constant.
 * 
 * @code
 *   const unsigned int dimension = 2;
 *   const unsigned int fe_degree = 4;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperation"></a> 
 * <h3>SineGordonOperation</h3>
 * 

 * 
 * The <code>SineGordonOperation</code> class implements the cell-based
 * operation that is needed in each time step. This nonlinear operation can
 * be implemented straight-forwardly based on the <code>MatrixFree</code>
 * class, in the same way as a linear operation would be treated by this
 * implementation of the finite element operator application. We apply two
 * template arguments to the class, one for the dimension and one for the
 * degree of the finite element. This is a difference to other functions in
 * deal.II where only the dimension is a template argument. This is
 * necessary to provide the inner loops in @p FEEvaluation with information
 * about loop lengths etc., which is essential for efficiency. On the other
 * hand, it makes it more challenging to implement the degree as a run-time
 * parameter.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class SineGordonOperation
 *   {
 *   public:
 *     SineGordonOperation(const MatrixFree<dim, double> &data_in,
 *                         const double                   time_step);
 * 
 *     void apply(LinearAlgebra::distributed::Vector<double> &dst,
 *                const std::vector<LinearAlgebra::distributed::Vector<double> *>
 *                  &src) const;
 * 
 *   private:
 *     const MatrixFree<dim, double> &            data;
 *     const VectorizedArray<double>              delta_t_sqr;
 *     LinearAlgebra::distributed::Vector<double> inv_mass_matrix;
 * 
 *     void local_apply(
 *       const MatrixFree<dim, double> &                                  data,
 *       LinearAlgebra::distributed::Vector<double> &                     dst,
 *       const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
 *       const std::pair<unsigned int, unsigned int> &cell_range) const;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationSineGordonOperation"></a> 
 * <h4>SineGordonOperation::SineGordonOperation</h4>
 * 

 * 
 * This is the constructor of the SineGordonOperation class. It receives a
 * reference to the MatrixFree holding the problem information and the time
 * step size as input parameters. The initialization routine sets up the
 * mass matrix. Since we use Gauss-Lobatto elements, the mass matrix is a
 * diagonal matrix and can be stored as a vector. The computation of the
 * mass matrix diagonal is simple to achieve with the data structures
 * provided by FEEvaluation: Just loop over all cell batches, i.e.,
 * collections of cells due to SIMD vectorization, and integrate over the
 * function that is constant one on all quadrature points by using the
 * <code>integrate</code> function with @p true argument at the slot for
 * values. Finally, we invert the diagonal entries to have the inverse mass
 * matrix directly available in each time step.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   SineGordonOperation<dim, fe_degree>::SineGordonOperation(
 *     const MatrixFree<dim, double> &data_in,
 *     const double                   time_step)
 *     : data(data_in)
 *     , delta_t_sqr(make_vectorized_array(time_step * time_step))
 *   {
 *     data.initialize_dof_vector(inv_mass_matrix);
 * 
 *     FEEvaluation<dim, fe_degree> fe_eval(data);
 *     const unsigned int           n_q_points = fe_eval.n_q_points;
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         fe_eval.reinit(cell);
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           fe_eval.submit_value(make_vectorized_array(1.), q);
 *         fe_eval.integrate(EvaluationFlags::values);
 *         fe_eval.distribute_local_to_global(inv_mass_matrix);
 *       }
 * 
 *     inv_mass_matrix.compress(VectorOperation::add);
 *     for (unsigned int k = 0; k < inv_mass_matrix.locally_owned_size(); ++k)
 *       if (inv_mass_matrix.local_element(k) > 1e-15)
 *         inv_mass_matrix.local_element(k) =
 *           1. / inv_mass_matrix.local_element(k);
 *       else
 *         inv_mass_matrix.local_element(k) = 1;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationlocal_apply"></a> 
 * <h4>SineGordonOperation::local_apply</h4>
 * 

 * 
 * This operator implements the core operation of the program, the
 * integration over a range of cells for the nonlinear operator of the
 * Sine-Gordon problem. The implementation is based on the FEEvaluation
 * class as in step-37. Due to the special structure in Gauss-Lobatto
 * elements, certain operations become simpler, in particular the evaluation
 * of shape function values on quadrature points which is simply the
 * injection of the values of cell degrees of freedom. The MatrixFree class
 * detects possible structure of the finite element at quadrature points
 * when initializing, which is then automatically used by FEEvaluation for
 * selecting the most appropriate numerical kernel.
 * 

 * 
 * The nonlinear function that we have to evaluate for the time stepping
 * routine includes the value of the function at the present time @p current
 * as well as the value at the previous time step @p old. Both values are
 * passed to the operator in the collection of source vectors @p src, which
 * is simply a <tt>std::vector</tt> of pointers to the actual solution
 * vectors. This construct of collecting several source vectors into one is
 * necessary as the cell loop in @p MatrixFree takes exactly one source and
 * one destination vector, even if we happen to use many vectors like the
 * two in this case. Note that the cell loop accepts any valid class for
 * input and output, which does not only include vectors but general data
 * types.  However, only in case it encounters a
 * LinearAlgebra::distributed::Vector<Number> or a <tt>std::vector</tt>
 * collecting these vectors, it calls functions that exchange ghost data due
 * to MPI at the beginning and the end of the loop. In the loop over the
 * cells, we first have to read in the values in the vectors related to the
 * local values.  Then, we evaluate the value and the gradient of the
 * current solution vector and the values of the old vector at the
 * quadrature points. Next, we combine the terms in the scheme in the loop
 * over the quadrature points. Finally, we integrate the result against the
 * test function and accumulate the result to the global solution vector @p
 * dst.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void SineGordonOperation<dim, fe_degree>::local_apply(
 *     const MatrixFree<dim> &                                          data,
 *     LinearAlgebra::distributed::Vector<double> &                     dst,
 *     const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
 *     const std::pair<unsigned int, unsigned int> &cell_range) const
 *   {
 *     AssertDimension(src.size(), 2);
 *     FEEvaluation<dim, fe_degree> current(data), old(data);
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         current.reinit(cell);
 *         old.reinit(cell);
 * 
 *         current.read_dof_values(*src[0]);
 *         old.read_dof_values(*src[1]);
 * 
 *         current.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
 *         old.evaluate(EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < current.n_q_points; ++q)
 *           {
 *             const VectorizedArray<double> current_value = current.get_value(q);
 *             const VectorizedArray<double> old_value     = old.get_value(q);
 * 
 *             current.submit_value(2. * current_value - old_value -
 *                                    delta_t_sqr * std::sin(current_value),
 *                                  q);
 *             current.submit_gradient(-delta_t_sqr * current.get_gradient(q), q);
 *           }
 * 
 *         current.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
 *         current.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationapply"></a> 
 * <h4>SineGordonOperation::apply</h4>
 * 

 * 
 * This function performs the time stepping routine based on the cell-local
 * strategy. Note that we need to set the destination vector to zero before
 * we add the integral contributions of the current time step (via the
 * FEEvaluation::distribute_local_to_global() call). In this tutorial, we
 * let the cell-loop do the zero operation via the fifth `true` argument
 * passed to MatrixFree::cell_loop. The loop can schedule the zero operation
 * closer to the operations on vector entries for supported vector entries,
 * thereby possibly increasing data locality (the vector entries that first
 * get zeroed are later re-used in the `distribute_local_to_global()`
 * call). The structure of the cell loop is implemented in the cell finite
 * element operator class. On each cell it applies the routine defined as
 * the <code>local_apply()</code> method of the class
 * <code>SineGordonOperation</code>, i.e., <code>this</code>. One could also
 * provide a function with the same signature that is not part of a
 * class. Finally, the result of the integration is multiplied by the
 * inverse mass matrix.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void SineGordonOperation<dim, fe_degree>::apply(
 *     LinearAlgebra::distributed::Vector<double> &                     dst,
 *     const std::vector<LinearAlgebra::distributed::Vector<double> *> &src) const
 *   {
 *     data.cell_loop(
 *       &SineGordonOperation<dim, fe_degree>::local_apply, this, dst, src, true);
 *     dst.scale(inv_mass_matrix);
 *   }
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
 * We define a time-dependent function that is used as initial
 * value. Different solutions can be obtained by varying the starting
 * time. This function, taken from step-25, would represent an analytic
 * solution in 1D for all times, but is merely used for setting some
 * starting solution of interest here. More elaborate choices that could
 * test the convergence of this program are given in step-25.
 * 
 * @code
 *   template <int dim>
 *   class InitialCondition : public Function<dim>
 *   {
 *   public:
 *     InitialCondition(const unsigned int n_components = 1,
 *                      const double       time         = 0.)
 *       : Function<dim>(n_components, time)
 *     {}
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/) const override
 *     {
 *       double t = this->get_time();
 * 
 *       const double m  = 0.5;
 *       const double c1 = 0.;
 *       const double c2 = 0.;
 *       const double factor =
 *         (m / std::sqrt(1. - m * m) * std::sin(std::sqrt(1. - m * m) * t + c2));
 *       double result = 1.;
 *       for (unsigned int d = 0; d < dim; ++d)
 *         result *= -4. * std::atan(factor / std::cosh(m * p[d] + c1));
 *       return result;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemclass"></a> 
 * <h3>SineGordonProblem class</h3>
 * 

 * 
 * This is the main class that builds on the class in step-25.  However, we
 * replaced the SparseMatrix<double> class by the MatrixFree class to store
 * the geometry data. Also, we use a distributed triangulation in this
 * example.
 * 
 * @code
 *   template <int dim>
 *   class SineGordonProblem
 *   {
 *   public:
 *     SineGordonProblem();
 *     void run();
 * 
 *   private:
 *     ConditionalOStream pcout;
 * 
 *     void make_grid_and_dofs();
 *     void output_results(const unsigned int timestep_number);
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     MappingQ1<dim> mapping;
 * 
 *     AffineConstraints<double> constraints;
 *     IndexSet                  locally_relevant_dofs;
 * 
 *     MatrixFree<dim, double> matrix_free_data;
 * 
 *     LinearAlgebra::distributed::Vector<double> solution, old_solution,
 *       old_old_solution;
 * 
 *     const unsigned int n_global_refinements;
 *     double             time, time_step;
 *     const double       final_time;
 *     const double       cfl_number;
 *     const unsigned int output_timestep_skip;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemSineGordonProblem"></a> 
 * <h4>SineGordonProblem::SineGordonProblem</h4>
 * 

 * 
 * This is the constructor of the SineGordonProblem class. The time interval
 * and time step size are defined here. Moreover, we use the degree of the
 * finite element that we defined at the top of the program to initialize a
 * FE_Q finite element based on Gauss-Lobatto support points. These points
 * are convenient because in conjunction with a QGaussLobatto quadrature
 * rule of the same order they give a diagonal mass matrix without
 * compromising accuracy too much (note that the integration is inexact,
 * though), see also the discussion in the introduction. Note that FE_Q
 * selects the Gauss-Lobatto nodal points by default due to their improved
 * conditioning versus equidistant points. To make things more explicit, we
 * state the selection of the nodal points nonetheless.
 * 
 * @code
 *   template <int dim>
 *   SineGordonProblem<dim>::SineGordonProblem()
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *     ,
 * #ifdef DEAL_II_WITH_P4EST
 *     triangulation(MPI_COMM_WORLD)
 *     ,
 * #endif
 *     fe(QGaussLobatto<1>(fe_degree + 1))
 *     , dof_handler(triangulation)
 *     , n_global_refinements(10 - 2 * dim)
 *     , time(-10)
 *     , time_step(10.)
 *     , final_time(10.)
 *     , cfl_number(.1 / fe_degree)
 *     , output_timestep_skip(200)
 *   {}
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemmake_grid_and_dofs"></a> 
 * <h4>SineGordonProblem::make_grid_and_dofs</h4>
 * 

 * 
 * As in step-25 this functions sets up a cube grid in <code>dim</code>
 * dimensions of extent $[-15,15]$. We refine the mesh more in the center of
 * the domain since the solution is concentrated there. We first refine all
 * cells whose center is within a radius of 11, and then refine once more
 * for a radius 6.  This simple ad hoc refinement could be done better by
 * adapting the mesh to the solution using error estimators during the time
 * stepping as done in other example programs, and using
 * parallel::distributed::SolutionTransfer to transfer the solution to the
 * new mesh.
 * 
 * @code
 *   template <int dim>
 *   void SineGordonProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -15, 15);
 *     triangulation.refine_global(n_global_refinements);
 *     {
 *       typename Triangulation<dim>::active_cell_iterator
 *         cell     = triangulation.begin_active(),
 *         end_cell = triangulation.end();
 *       for (; cell != end_cell; ++cell)
 *         if (cell->is_locally_owned())
 *           if (cell->center().norm() < 11)
 *             cell->set_refine_flag();
 *       triangulation.execute_coarsening_and_refinement();
 * 
 *       cell     = triangulation.begin_active();
 *       end_cell = triangulation.end();
 *       for (; cell != end_cell; ++cell)
 *         if (cell->is_locally_owned())
 *           if (cell->center().norm() < 6)
 *             cell->set_refine_flag();
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * 
 *     pcout << "   Number of global active cells: "
 * #ifdef DEAL_II_WITH_P4EST
 *           << triangulation.n_global_active_cells()
 * #else
 *           << triangulation.n_active_cells()
 * #endif
 *           << std::endl;
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl;
 * 
 * 
 * @endcode
 * 
 * We generate hanging node constraints for ensuring continuity of the
 * solution. As in step-40, we need to equip the constraint matrix with
 * the IndexSet of locally relevant degrees of freedom to avoid it to
 * consume too much memory for big problems. Next, the <code> MatrixFree
 * </code> object for the problem is set up. Note that we specify a
 * particular scheme for shared-memory parallelization (hence one would
 * use multithreading for intra-node parallelism and not MPI; we here
 * choose the standard option &mdash; if we wanted to disable shared
 * memory parallelization even in case where there is more than one TBB
 * thread available in the program, we would choose
 * MatrixFree::AdditionalData::TasksParallelScheme::none). Also note that,
 * instead of using the default QGauss quadrature argument, we supply a
 * QGaussLobatto quadrature formula to enable the desired
 * behavior. Finally, three solution vectors are initialized. MatrixFree
 * expects a particular layout of ghost indices (as it handles index
 * access in MPI-local numbers that need to match between the vector and
 * MatrixFree), so we just ask it to initialize the vectors to be sure the
 * ghost exchange is properly handled.
 * 
 * @code
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     constraints.close();
 * 
 *     typename MatrixFree<dim>::AdditionalData additional_data;
 *     additional_data.tasks_parallel_scheme =
 *       MatrixFree<dim>::AdditionalData::TasksParallelScheme::partition_partition;
 * 
 *     matrix_free_data.reinit(mapping,
 *                             dof_handler,
 *                             constraints,
 *                             QGaussLobatto<1>(fe_degree + 1),
 *                             additional_data);
 * 
 *     matrix_free_data.initialize_dof_vector(solution);
 *     old_solution.reinit(solution);
 *     old_old_solution.reinit(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemoutput_results"></a> 
 * <h4>SineGordonProblem::output_results</h4>
 * 

 * 
 * This function prints the norm of the solution and writes the solution
 * vector to a file. The norm is standard (except for the fact that we need
 * to accumulate the norms over all processors for the parallel grid which
 * we do via the VectorTools::compute_global_error() function), and the
 * second is similar to what we did in step-40 or step-37. Note that we can
 * use the same vector for output as the one used during computations: The
 * vectors in the matrix-free framework always provide full information on
 * all locally owned cells (this is what is needed in the local evaluations,
 * too), including ghost vector entries on these cells. This is the only
 * data that is needed in the VectorTools::integrate_difference() function
 * as well as in DataOut. The only action to take at this point is to make
 * sure that the vector updates its ghost values before we read from
 * them, and to reset ghost values once done. This is a feature present only
 * in the LinearAlgebra::distributed::Vector class. Distributed vectors with
 * PETSc and Trilinos, on the other hand, need to be copied to special
 * vectors including ghost values (see the relevant section in step-40). If
 * we also wanted to access all degrees of freedom on ghost cells (e.g. when
 * computing error estimators that use the jump of solution over cell
 * boundaries), we would need more information and create a vector
 * initialized with locally relevant dofs just as in step-40. Observe also
 * that we need to distribute constraints for output - they are not filled
 * during computations (rather, they are interpolated on the fly in the
 * matrix-free method FEEvaluation::read_dof_values()).
 * 
 * @code
 *   template <int dim>
 *   void
 *   SineGordonProblem<dim>::output_results(const unsigned int timestep_number)
 *   {
 *     constraints.distribute(solution);
 * 
 *     Vector<float> norm_per_cell(triangulation.n_active_cells());
 *     solution.update_ghost_values();
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Functions::ZeroFunction<dim>(),
 *                                       norm_per_cell,
 *                                       QGauss<dim>(fe_degree + 1),
 *                                       VectorTools::L2_norm);
 *     const double solution_norm =
 *       VectorTools::compute_global_error(triangulation,
 *                                         norm_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 *     pcout << "   Time:" << std::setw(8) << std::setprecision(3) << time
 *           << ", solution norm: " << std::setprecision(5) << std::setw(7)
 *           << solution_norm << std::endl;
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches(mapping);
 * 
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", timestep_number, MPI_COMM_WORLD, 3);
 * 
 *     solution.zero_out_ghost_values();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemrun"></a> 
 * <h4>SineGordonProblem::run</h4>
 * 

 * 
 * This function is called by the main function and steps into the
 * subroutines of the class.
 *   

 * 
 * After printing some information about the parallel setup, the first
 * action is to set up the grid and the cell operator. Then, the time step
 * is computed from the CFL number given in the constructor and the finest
 * mesh size. The finest mesh size is computed as the diameter of the last
 * cell in the triangulation, which is the last cell on the finest level of
 * the mesh. This is only possible for meshes where all elements on a level
 * have the same size, otherwise, one needs to loop over all cells. Note
 * that we need to query all the processors for their finest cell since
 * not all processors might hold a region where the mesh is at the finest
 * level. Then, we readjust the time step a little to hit the final time
 * exactly.
 * 
 * @code
 *   template <int dim>
 *   void SineGordonProblem<dim>::run()
 *   {
 *     {
 *       pcout << "Number of MPI ranks:            "
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
 *       pcout << "Number of threads on each rank: "
 *             << MultithreadInfo::n_threads() << std::endl;
 *       const unsigned int n_vect_doubles = VectorizedArray<double>::size();
 *       const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;
 *       pcout << "Vectorization over " << n_vect_doubles
 *             << " doubles = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl
 *             << std::endl;
 *     }
 *     make_grid_and_dofs();
 * 
 *     const double local_min_cell_diameter =
 *       triangulation.last()->diameter() / std::sqrt(dim);
 *     const double global_min_cell_diameter =
 *       -Utilities::MPI::max(-local_min_cell_diameter, MPI_COMM_WORLD);
 *     time_step = cfl_number * global_min_cell_diameter;
 *     time_step = (final_time - time) / (int((final_time - time) / time_step));
 *     pcout << "   Time step size: " << time_step
 *           << ", finest cell: " << global_min_cell_diameter << std::endl
 *           << std::endl;
 * 
 * @endcode
 * 
 * Next the initial value is set. Since we have a two-step time stepping
 * method, we also need a value of the solution at time-time_step. For
 * accurate results, one would need to compute this from the time
 * derivative of the solution at initial time, but here we ignore this
 * difficulty and just set it to the initial value function at that
 * artificial time.
 * 

 * 
 * We then go on by writing the initial state to file and collecting
 * the two starting solutions in a <tt>std::vector</tt> of pointers that
 * get later consumed by the SineGordonOperation::apply() function. Next,
 * an instance of the <code> SineGordonOperation class </code> based on
 * the finite element degree specified at the top of this file is set up.
 * 
 * @code
 *     VectorTools::interpolate(mapping,
 *                              dof_handler,
 *                              InitialCondition<dim>(1, time),
 *                              solution);
 *     VectorTools::interpolate(mapping,
 *                              dof_handler,
 *                              InitialCondition<dim>(1, time - time_step),
 *                              old_solution);
 *     output_results(0);
 * 
 *     std::vector<LinearAlgebra::distributed::Vector<double> *>
 *       previous_solutions({&old_solution, &old_old_solution});
 * 
 *     SineGordonOperation<dim, fe_degree> sine_gordon_op(matrix_free_data,
 *                                                        time_step);
 * 
 * @endcode
 * 
 * Now loop over the time steps. In each iteration, we shift the solution
 * vectors by one and call the `apply` function of the
 * `SineGordonOperator` class. Then, we write the solution to a file. We
 * clock the wall times for the computational time needed as wall as the
 * time needed to create the output and report the numbers when the time
 * stepping is finished.
 *     

 * 
 * Note how this shift is implemented: We simply call the swap method on
 * the two vectors which swaps only some pointers without the need to copy
 * data around, a relatively expensive operation within an explicit time
 * stepping method. Let us see what happens in more detail: First, we
 * exchange <code>old_solution</code> with <code>old_old_solution</code>,
 * which means that <code>old_old_solution</code> gets
 * <code>old_solution</code>, which is what we expect. Similarly,
 * <code>old_solution</code> gets the content from <code>solution</code>
 * in the next step. After this, <code>solution</code> holds
 * <code>old_old_solution</code>, but that will be overwritten during this
 * step.
 * 
 * @code
 *     unsigned int timestep_number = 1;
 * 
 *     Timer  timer;
 *     double wtime       = 0;
 *     double output_time = 0;
 *     for (time += time_step; time <= final_time;
 *          time += time_step, ++timestep_number)
 *       {
 *         timer.restart();
 *         old_old_solution.swap(old_solution);
 *         old_solution.swap(solution);
 *         sine_gordon_op.apply(solution, previous_solutions);
 *         wtime += timer.wall_time();
 * 
 *         timer.restart();
 *         if (timestep_number % output_timestep_skip == 0)
 *           output_results(timestep_number / output_timestep_skip);
 * 
 *         output_time += timer.wall_time();
 *       }
 *     timer.restart();
 *     output_results(timestep_number / output_timestep_skip + 1);
 *     output_time += timer.wall_time();
 * 
 *     pcout << std::endl
 *           << "   Performed " << timestep_number << " time steps." << std::endl;
 * 
 *     pcout << "   Average wallclock time per time step: "
 *           << wtime / timestep_number << "s" << std::endl;
 * 
 *     pcout << "   Spent " << output_time << "s on output and " << wtime
 *           << "s on computations." << std::endl;
 *   }
 * } // namespace Step48
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * As in step-40, we initialize MPI at the start of the program. Since we will
 * in general mix MPI parallelization with threads, we also set the third
 * argument in MPI_InitFinalize that controls the number of threads to an
 * invalid number, which means that the TBB library chooses the number of
 * threads automatically, typically to the number of available cores in the
 * system. As an alternative, you can also set this number manually if you
 * want to set a specific number of threads (e.g. when MPI-only is required).
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   using namespace Step48;
 *   using namespace dealii;
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *     argc, argv, numbers::invalid_unsigned_int);
 * 
 *   try
 *     {
 *       SineGordonProblem<dimension> sg_problem;
 *       sg_problem.run();
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


<a name="Comparisonwithasparsematrix"></a><h3>Comparison with a sparse matrix</h3>


In order to demonstrate the gain in using the MatrixFree class instead of
the standard <code>deal.II</code> assembly routines for evaluating the
information from old time steps, we study a simple serial run of the code on a
nonadaptive mesh. Since much time is spent on evaluating the sine function, we
do not only show the numbers of the full sine-Gordon equation but also for the
wave equation (the sine-term skipped from the sine-Gordon equation). We use
both second and fourth order elements. The results are summarized in the
following table.

<table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3">wave equation</th>
    <th colspan="2">sine-Gordon</th>
  </tr>
  <tr>
    <th>&nbsp;</th>
    <th>MF</th>
    <th>SpMV</th>
    <th>dealii</th>
    <th>MF</th>
    <th>dealii</th>
  </tr>
  <tr>
    <td>2D, $\mathcal{Q}_2$</td>
    <td align="right"> 0.0106</td>
    <td align="right"> 0.00971</td>
    <td align="right"> 0.109</td>
    <td align="right"> 0.0243</td>
    <td align="right"> 0.124</td>
  </tr>
  <tr>
    <td>2D, $\mathcal{Q}_4$</td>
    <td align="right"> 0.0328</td>
    <td align="right"> 0.0706</td>
    <td align="right"> 0.528</td>
    <td align="right"> 0.0714</td>
    <td align="right"> 0.502</td>
   </tr>
   <tr>
    <td>3D, $\mathcal{Q}_2$</td>
    <td align="right"> 0.0151</td>
    <td align="right"> 0.0320</td>
    <td align="right"> 0.331</td>
    <td align="right"> 0.0376</td>
    <td align="right"> 0.364</td>
   </tr>
   <tr>
    <td>3D, $\mathcal{Q}_4$</td>
    <td align="right"> 0.0918</td>
    <td align="right"> 0.844</td>
    <td align="right"> 6.83</td>
    <td align="right"> 0.194</td>
    <td align="right"> 6.95</td>
   </tr>
</table>

It is apparent that the matrix-free code outperforms the standard assembly
routines in deal.II by far. In 3D and for fourth order elements, one operator
evaluation is also almost ten times as fast as a sparse matrix-vector
product.

<a name="Parallelrunin2Dand3D"></a><h3>Parallel run in 2D and 3D</h3>


We start with the program output obtained on a workstation with 12 cores / 24
threads (one Intel Xeon E5-2687W v4 CPU running at 3.2 GHz, hyperthreading
enabled), running the program in release mode:
@code
\$ make run
Number of MPI ranks:            1
Number of threads on each rank: 24
Vectorization over 4 doubles = 256 bits (AVX)

   Number of global active cells: 15412
   Number of degrees of freedom: 249065
   Time step size: 0.00292997, finest cell: 0.117188

   Time:     -10, solution norm:  9.5599
   Time:   -9.41, solution norm:  17.678
   Time:   -8.83, solution norm:  23.504
   Time:   -8.24, solution norm:    27.5
   Time:   -7.66, solution norm:  29.513
   Time:   -7.07, solution norm:  29.364
   Time:   -6.48, solution norm:   27.23
   Time:    -5.9, solution norm:  23.527
   Time:   -5.31, solution norm:  18.439
   Time:   -4.73, solution norm:  11.935
   Time:   -4.14, solution norm:  5.5284
   Time:   -3.55, solution norm:  8.0354
   Time:   -2.97, solution norm:  14.707
   Time:   -2.38, solution norm:      20
   Time:    -1.8, solution norm:  22.834
   Time:   -1.21, solution norm:  22.771
   Time:  -0.624, solution norm:  20.488
   Time: -0.0381, solution norm:  16.697
   Time:   0.548, solution norm:  11.221
   Time:    1.13, solution norm:  5.3912
   Time:    1.72, solution norm:  8.4528
   Time:    2.31, solution norm:  14.335
   Time:    2.89, solution norm:  18.555
   Time:    3.48, solution norm:  20.894
   Time:    4.06, solution norm:  21.305
   Time:    4.65, solution norm:  19.903
   Time:    5.24, solution norm:  16.864
   Time:    5.82, solution norm:  12.223
   Time:    6.41, solution norm:   6.758
   Time:    6.99, solution norm:  7.2423
   Time:    7.58, solution norm:  12.888
   Time:    8.17, solution norm:  17.273
   Time:    8.75, solution norm:  19.654
   Time:    9.34, solution norm:  19.838
   Time:    9.92, solution norm:  17.964
   Time:      10, solution norm:  17.595

   Performed 6826 time steps.
   Average wallclock time per time step: 0.0013453s
   Spent 14.976s on output and 9.1831s on computations.
@endcode

In 3D, the respective output looks like
@code
\$ make run
Number of MPI ranks:            1
Number of threads on each rank: 24
Vectorization over 4 doubles = 256 bits (AVX)

   Number of global active cells: 17592
   Number of degrees of freedom: 1193881
   Time step size: 0.0117233, finest cell: 0.46875

   Time:     -10, solution norm:  29.558
   Time:   -7.66, solution norm:  129.13
   Time:   -5.31, solution norm:  67.753
   Time:   -2.97, solution norm:  79.245
   Time:  -0.621, solution norm:  123.52
   Time:    1.72, solution norm:  43.525
   Time:    4.07, solution norm:  93.285
   Time:    6.41, solution norm:  97.722
   Time:    8.76, solution norm:  36.734
   Time:      10, solution norm:  94.115

   Performed 1706 time steps.
   Average wallclock time per time step: 0.0084542s
   Spent 16.766s on output and 14.423s on computations.
@endcode

It takes 0.008 seconds for one time step with more than a million
degrees of freedom (note that we would need many processors to reach such
numbers when solving linear systems).

If we replace the thread-parallelization by a pure MPI parallelization, the
timings change into:
@code
\$ mpirun -n 24 ./step-48
Number of MPI ranks:            24
Number of threads on each rank: 1
Vectorization over 4 doubles = 256 bits (AVX)
...
   Performed 1706 time steps.
   Average wallclock time per time step: 0.0051747s
   Spent 2.0535s on output and 8.828s on computations.
@endcode

We observe a dramatic speedup for the output (which makes sense, given that
most code of the output is not parallelized via threads, whereas it is for
MPI), but less than the theoretical factor of 12 we would expect from the
parallelism. More interestingly, the computations also get faster when
switching from the threads-only variant to the MPI-only variant. This is a
general observation for the MatrixFree framework (as of updating this data in
2019). The main reason is that the decisions regarding work on conflicting
cell batches made to enable execution in parallel are overly pessimistic:
While they ensure that no work on neighboring cells is done on different
threads at the same time, this conservative setting implies that data from
neighboring cells is also evicted from caches by the time neighbors get
touched. Furthermore, the current scheme is not able to provide a constant
load for all 24 threads for the given mesh with 17,592 cells.

The current program allows to also mix MPI parallelization with thread
parallelization. This is most beneficial when running programs on clusters
with multiple nodes, using MPI for the inter-node parallelization and threads
for the intra-node parallelization. On the workstation used above, we can run
threads in the hyperthreading region (i.e., using 2 threads for each of the 12
MPI ranks). An important setting for mixing MPI with threads is to ensure
proper binning of tasks to CPUs. On many clusters the placing is either
automatically via the `mpirun/mpiexec` environment, or there can be manual
settings. Here, we simply report the run times the plain version of the
program (noting that things could be improved towards the timings of the
MPI-only program when proper pinning is done):
@code
\$ mpirun -n 12 ./step-48
Number of MPI ranks:            12
Number of threads on each rank: 2
Vectorization over 4 doubles = 256 bits (AVX)
...
   Performed 1706 time steps.
   Average wallclock time per time step: 0.0056651s
   Spent 2.5175s on output and 9.6646s on computations.
@endcode



<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


There are several things in this program that could be improved to make it
even more efficient (besides improved boundary conditions and physical
stuff as discussed in step-25):

<ul> <li> <b>Faster evaluation of sine terms:</b> As becomes obvious
  from the comparison of the plain wave equation and the sine-Gordon
  equation above, the evaluation of the sine terms dominates the total
  time for the finite element operator application. There are a few
  reasons for this: Firstly, the deal.II sine computation of a
  VectorizedArray field is not vectorized (as opposed to the rest of
  the operator application). This could be cured by handing the sine
  computation to a library with vectorized sine computations like
  Intel's math kernel library (MKL). By using the function
  <code>vdSin</code> in MKL, the program uses half the computing time
  in 2D and 40 percent less time in 3D. On the other hand, the sine
  computation is structurally much more complicated than the simple
  arithmetic operations like additions and multiplications in the rest
  of the local operation.

  <li> <b>Higher order time stepping:</b> While the implementation allows for
  arbitrary order in the spatial part (by adjusting the degree of the finite
  element), the time stepping scheme is a standard second-order leap-frog
  scheme. Since solutions in wave propagation problems are usually very
  smooth, the error is likely dominated by the time stepping part. Of course,
  this could be cured by using smaller time steps (at a fixed spatial
  resolution), but it would be more efficient to use higher order time
  stepping as well. While it would be straight-forward to do so for a
  first-order system (use some Runge&ndash;Kutta scheme of higher order,
  probably combined with adaptive time step selection like the <a
  href="http://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method">Dormand&ndash;Prince
  method</a>), it is more challenging for the second-order formulation. At
  least in the finite difference community, people usually use the PDE to find
  spatial correction terms that improve the temporal error.

</ul>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-48.cc"
*/
