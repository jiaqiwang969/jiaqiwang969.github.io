/**
@page step_7 The step-7 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Verificationofcorrectness">Verification of correctness</a>
        <li><a href="#NonhomogeneousNeumannboundaryconditions">Non-homogeneous Neumann boundary conditions</a>
        <li><a href="#Anoteongoodprogrammingpractice">A note on good programming practice</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#TheHelmholtzsolverclass">The Helmholtz solver class</a>
        <li><a href="#TheHelmholtzProblemclassimplementation">The HelmholtzProblem class implementation</a>
      <ul>
        <li><a href="#HelmholtzProblemHelmholtzProblemconstructor">HelmholtzProblem::HelmholtzProblem constructor</a>
        <li><a href="#HelmholtzProblemsetup_system">HelmholtzProblem::setup_system</a>
        <li><a href="#HelmholtzProblemassemble_system">HelmholtzProblem::assemble_system</a>
        <li><a href="#HelmholtzProblemsolve">HelmholtzProblem::solve</a>
        <li><a href="#HelmholtzProblemrefine_grid">HelmholtzProblem::refine_grid</a>
        <li><a href="#HelmholtzProblemprocess_solution">HelmholtzProblem::process_solution</a>
        <li><a href="#HelmholtzProblemrun">HelmholtzProblem::run</a>
      <ul>
        <li><a href="#Outputofgraphicaldata">Output of graphical data</a>
        <li><a href="#Outputofconvergencetables">Output of convergence tables</a>
        <li><a href="#Furthertablemanipulations">Further table manipulations</a>
      </ul>
      </ul>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
      <ul>
        <li><a href="#Whenistheerrorsmall"> When is the error "small"? </a>
      </ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#HigherOrderElements"> Higher Order Elements </a>
        <li><a href="#ConvergenceComparison"> Convergence Comparison </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


In this program, we will mainly consider two aspects:
<ol>
  <li> Verification of correctness of the program and generation of convergence
  tables;
  <li> Non-homogeneous Neumann boundary conditions for the Helmholtz equation.
</ol>
Besides these topics, again a variety of improvements and tricks will be
shown.


<a name="Verificationofcorrectness"></a><h3>Verification of correctness</h3>


There has probably never been a
non-trivial finite element program that worked right from the start. It is
therefore necessary to find ways to verify whether a computed solution is
correct or not. Usually, this is done by choosing the set-up of a simulation
in such a way that we know the exact continuous solution and evaluate the difference
between continuous and computed discrete solution. If this difference
converges to zero with the right order of convergence, this is already a good
indication of correctness, although there may be other sources of error
persisting which have only a small contribution to the total error or are of
higher order. In the context of finite element simulations, this technique
of picking the solution by choosing appropriate right hand sides and
boundary conditions
is often called the <i>Method of Manufactured Solution</i>.

In this example, we will not go into the theories of systematic software
verification which is a very complicated problem. Rather we will demonstrate
the tools which deal.II can offer in this respect. This is basically centered
around the functionality of a single function, VectorTools::integrate_difference().
This function computes the difference between a given continuous function and
a finite element field in various norms on each cell.
Of course, like with any other integral, we can only evaluate these norms using quadrature formulas;
the choice of the right quadrature formula is therefore crucial to the
accurate evaluation of the error. This holds in particular for the $L_\infty$
norm, where we evaluate the maximal deviation of numerical and exact solution
only at the quadrature points; one should then not try to use a quadrature
rule whose evaluation occurs only at points where
[super-convergence](https://en.wikipedia.org/wiki/Superconvergence) might occur, such as
the Gauss points of the lowest-order Gauss quadrature formula for which the
integrals in the assembly of the matrix is correct (e.g., for linear elements,
do not use the QGauss(2) quadrature formula). In fact, this is generally good
advice also for the other norms: if your quadrature points are fortuitously
chosen at locations where the error happens to be particularly small due to
superconvergence, the computed error will look like it is much smaller than
it really is and may even suggest a higher convergence order. Consequently,
we will choose a different quadrature formula for the integration of these
error norms than for the assembly of the linear system.

The function VectorTools::integrate_difference() evaluates the desired norm on each
cell $K$ of the triangulation and returns a vector which holds these
values for each cell. From the local values, we can then obtain the global error. For
example, if the vector $\mathbf e$ with element $e_K$ for all cells
$K$ contains the local $L_2$ norms $\|u-u_h\|_K$, then
@f[
  E = \| {\mathbf e} \| = \left( \sum_K e_K^2 \right)^{1/2}
@f]
is the global $L_2$ error $E=\|u-u_h\|_\Omega$.

In the program, we will show how to evaluate and use these quantities, and we
will monitor their values under mesh refinement. Of course, we have to choose
the problem at hand such that we can explicitly state the solution and its
derivatives, but since we want to evaluate the correctness of the program,
this is only reasonable. If we know that the program produces the correct
solution for one (or, if one wants to be really sure: many) specifically
chosen right hand sides, we can be rather confident that it will also compute
the correct solution for problems where we don't know the exact values.

In addition to simply computing these quantities, we will show how to generate
nicely formatted tables from the data generated by this program that
automatically computes convergence rates etc. In addition, we will compare
different strategies for mesh refinement.


<a name="NonhomogeneousNeumannboundaryconditions"></a><h3>Non-homogeneous Neumann boundary conditions</h3>


The second, totally
unrelated, subject of this example program is the use of non-homogeneous
boundary conditions. These are included into the variational form using
boundary integrals which we have to evaluate numerically when assembling the
right hand side vector.

Before we go into programming, let's have a brief look at the mathematical
formulation. The equation that we want to solve here is the Helmholtz equation
"with the nice sign":
@f[
  -\Delta u + \alpha u = f,
@f]
on the square $[-1,1]^2$ with $\alpha=1$, augmented by Dirichlet boundary conditions
@f[
  u = g_1
@f]
on some part $\Gamma_1$ of the boundary $\Gamma$, and Neumann conditions
@f[
  {\mathbf n}\cdot \nabla u = g_2
@f]
on the rest $\Gamma_2 = \Gamma \backslash \Gamma_1$.
In our particular testcase, we will use $\Gamma_1=\Gamma \cap\{\{x=1\}
\cup \{y=1\}\}$.
(We say that this equation has the "nice sign" because the operator
$-\Delta + \alpha I$ with the identity $I$ and $\alpha>0$ is a positive definite
operator; the <a
href="https://en.wikipedia.org/wiki/Helmholtz_equation">equation with
the "bad sign"</a> is $-\Delta u - \alpha u$ and results from modeling
time-harmonic processes. The operator is not positive
definite if $\alpha>0$ is large, and this leads to all sorts of issues
we need not discuss here. The operator may also not be invertible --
i.e., the equation does not have a unique solution -- if $\alpha$
happens to be one of the eigenvalues of $-\Delta$.)

Because we want to verify the convergence of our numerical solution $u_h$,
we want a setup so that we know the exact solution $u$. This is where
the Method of Manufactured Solutions comes in. To this end, let us
choose a function
@f[
  \bar u(x) = \sum_{i=1}^3 \exp\left(-\frac{|x-x_i|^2}{\sigma^2}\right)
@f]
where the centers $x_i$ of the exponentials are
  $x_1=(-\frac 12,\frac 12)$,
  $x_2=(-\frac 12,-\frac 12)$, and
  $x_3=(\frac 12,-\frac 12)$,
and the half width is set to $\sigma=\frac {1}{8}$. The method of manufactured
solution then says: choose
@f{align*}
  f &= -\Delta \bar u + \bar u, \\
  g_1 &= \bar u|_{\Gamma_1}, \\
  g_2 &= {\mathbf n}\cdot \nabla\bar u|_{\Gamma_2}.
@f}
With this particular choice, we infer that of course the solution of the
original problem happens to be $u=\bar u$. In other words, by choosing
the right hand sides of the equation and the boundary conditions in a
particular way, we have manufactured ourselves a problem to which we
know the solution. This allows us then to compute the error of our
numerical solution. In the code below, we represent $\bar u$ by the
<code>Solution</code> class, and other classes will be used to
denote $\bar u|_{\Gamma_1}=g_1$ and ${\mathbf n}\cdot \nabla\bar u|_{\Gamma_2}=g_2$.

Using the above definitions, we can state the weak formulation of the
equation, which reads: find $u\in H^1_g=\{v\in H^1: v|_{\Gamma_1}=g_1\}$ such
that
@f[
  {(\nabla v, \nabla u)}_\Omega + {(v,u)}_\Omega
  =
  {(v,f)}_\Omega + {(v,g_2)}_{\Gamma_2}
@f]
for all test functions $v\in H^1_0=\{v\in H^1: v|_{\Gamma_1}=0\}$. The
boundary term ${(v,g_2)}_{\Gamma_2}$ has appeared by integration by parts and
using $\partial_n u=g_2$ on $\Gamma_2$ and $v=0$ on $\Gamma_1$. The cell
matrices and vectors which we use to build the global matrices and right hand
side vectors in the discrete formulation therefore look like this:
@f{eqnarray*}
  A_{ij}^K &=& \left(\nabla \varphi_i, \nabla \varphi_j\right)_K
              +\left(\varphi_i, \varphi_j\right)_K,
  \\
  F_i^K &=& \left(\varphi_i, f\right)_K
           +\left(\varphi_i, g_2\right)_{\partial K\cap \Gamma_2}.
@f}
Since the generation of the domain integrals has been shown in previous
examples several times, only the generation of the contour integral is of
interest here. It basically works along the following lines: for domain
integrals we have the <code>FEValues</code> class that provides values and
gradients of the shape values, as well as Jacobian determinants and other
information and specified quadrature points in the cell; likewise, there is a
class <code>FEFaceValues</code> that performs these tasks for integrations on
faces of cells. One provides it with a quadrature formula for a manifold with
dimension one less than the dimension of the domain is, and the cell and the
number of its face on which we want to perform the integration. The class will
then compute the values, gradients, normal vectors, weights, etc. at the
quadrature points on this face, which we can then use in the same way as for
the domain integrals. The details of how this is done are shown in the
following program.


<a name="Anoteongoodprogrammingpractice"></a><h3>A note on good programming practice</h3>


Besides the mathematical topics outlined above, we also want to use this
program to illustrate one aspect of good programming practice, namely the use
of namespaces. In programming the deal.II library, we have take great care not
to use names for classes and global functions that are overly generic, say
<code>f(), sz(), rhs()</code> etc. Furthermore, we have put everything into
namespace <code>dealii</code>. But when one writes application programs that
aren't meant for others to use, one doesn't always pay this much attention. If
you follow the programming style of step-1 through step-6, these functions
then end up in the global namespace where, unfortunately, a lot of other stuff
also lives (basically everything the C language provides, along with
everything you get from the operating system through header files). To make
things a bit worse, the designers of the C language were also not always
careful in avoiding generic names; for example, the symbols <code>j1,
jn</code> are defined in C header files (they denote Bessel functions).

To avoid the problems that result if names of different functions or variables
collide (often with confusing error messages), it is good practice to put
everything you do into a <a
href="http://en.wikipedia.org/wiki/Namespace_(computer_science)">namespace</a>. Following
this style, we will open a namespace <code>Step7</code> at the top of the
program, import the deal.II namespace into it, put everything that's specific
to this program (with the exception of <code>main()</code>, which must be in
the global namespace) into it, and only close it at the bottom of the file. In
other words, the structure of the program is of the kind
@code
  #includes ...

  namespace Step7
  {
    using namespace dealii;

    ...everything to do with the program...
  }

  int main ()
  {
    ...do whatever main() does...
  }
@endcode
We will follow this scheme throughout the remainder of the deal.II tutorial.
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
 * These first include files have all been treated in previous examples, so we
 * won't explain what is in them again.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * In this example, we will not use the numeration scheme which is used per
 * default by the DoFHandler class, but will renumber them using the
 * Cuthill-McKee algorithm. As has already been explained in step-2, the
 * necessary functions are declared in the following file:
 * 
 * @code
 * #include <deal.II/dofs/dof_renumbering.h>
 * @endcode
 * 
 * Then we will show a little trick how we can make sure that objects are not
 * deleted while they are still in use. For this purpose, deal.II has the
 * SmartPointer helper class, which is declared in this file:
 * 
 * @code
 * #include <deal.II/base/smartpointer.h>
 * @endcode
 * 
 * Next, we will want to use the function VectorTools::integrate_difference()
 * mentioned in the introduction, and we are going to use a ConvergenceTable
 * that collects all important data during a run and prints it at the end as a
 * table. These comes from the following two files:
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/base/convergence_table.h>
 * @endcode
 * 
 * And finally, we need to use the FEFaceValues class, which is declared in
 * the same file as the FEValues class:
 * 
 * @code
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <array>
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The last step before we go on with the actual implementation is to open a
 * namespace <code>Step7</code> into which we will put everything, as
 * discussed at the end of the introduction, and to import the members of
 * namespace <code>dealii</code> into it:
 * 
 * @code
 * namespace Step7
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
 * Before implementing the classes that actually solve something, we first
 * declare and define some function classes that represent right hand side
 * and solution classes. Since we want to compare the numerically obtained
 * solution to the exact continuous one, we need a function object that
 * represents the continuous solution. On the other hand, we need the right
 * hand side function, and that one of course shares some characteristics
 * with the solution. In order to reduce dependencies which arise if we have
 * to change something in both classes at the same time, we move the common
 * characteristics of both functions into a base class.
 *   

 * 
 * The common characteristics for solution (as explained in the
 * introduction, we choose a sum of three exponentials) and right hand side,
 * are these: the number of exponentials, their centers, and their half
 * width. We declare them in the following class. Since the number of
 * exponentials is a compile-time constant we use a fixed-length
 * <code>std::array</code> to store the center points:
 * 
 * @code
 *   template <int dim>
 *   class SolutionBase
 *   {
 *   protected:
 *     static const std::array<Point<dim>, 3> source_centers;
 *     static const double                    width;
 *   };
 * 
 * 
 * @endcode
 * 
 * The variables which denote the centers and the width of the exponentials
 * have just been declared, now we still need to assign values to
 * them. Here, we can show another small piece of template sorcery, namely
 * how we can assign different values to these variables depending on the
 * dimension. We will only use the 2d case in the program, but we show the
 * 1d case for exposition of a useful technique.
 *   

 * 
 * First we assign values to the centers for the 1d case, where we place the
 * centers equidistantly at -1/3, 0, and 1/3. The <code>template
 * &lt;&gt;</code> header for this definition indicates an explicit
 * specialization. This means, that the variable belongs to a template, but
 * that instead of providing the compiler with a template from which it can
 * specialize a concrete variable by substituting <code>dim</code> with some
 * concrete value, we provide a specialization ourselves, in this case for
 * <code>dim=1</code>. If the compiler then sees a reference to this
 * variable in a place where the template argument equals one, it knows that
 * it doesn't have to generate the variable from a template by substituting
 * <code>dim</code>, but can immediately use the following definition:
 * 
 * @code
 *   template <>
 *   const std::array<Point<1>, 3> SolutionBase<1>::source_centers = {
 *     {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}};
 * 
 * @endcode
 * 
 * Likewise, we can provide an explicit specialization for
 * <code>dim=2</code>. We place the centers for the 2d case as follows:
 * 
 * @code
 *   template <>
 *   const std::array<Point<2>, 3> SolutionBase<2>::source_centers = {
 *     {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}};
 * 
 * @endcode
 * 
 * There remains to assign a value to the half-width of the exponentials. We
 * would like to use the same value for all dimensions. In this case, we
 * simply provide the compiler with a template from which it can generate a
 * concrete instantiation by substituting <code>dim</code> with a concrete
 * value:
 * 
 * @code
 *   template <int dim>
 *   const double SolutionBase<dim>::width = 1. / 8.;
 * 
 * 
 * 
 * @endcode
 * 
 * After declaring and defining the characteristics of solution and right
 * hand side, we can declare the classes representing these two. They both
 * represent continuous functions, so they are derived from the
 * Function&lt;dim&gt; base class, and they also inherit the characteristics
 * defined in the SolutionBase class.
 *   

 * 
 * The actual classes are declared in the following. Note that in order to
 * compute the error of the numerical solution against the continuous one in
 * the L2 and H1 (semi-)norms, we have to provide value and gradient of the
 * exact solution. This is more than we have done in previous examples, where
 * all we provided was the value at one or a list of points. Fortunately, the
 * Function class also has virtual functions for the gradient, so we can
 * simply overload the respective virtual member functions in the Function
 * base class. Note that the gradient of a function in <code>dim</code>
 * space dimensions is a vector of size <code>dim</code>, i.e. a tensor of
 * rank 1 and dimension <code>dim</code>. As for so many other things, the
 * library provides a suitable class for this. One new thing about this
 * class is that it explicitly uses the Tensor objects, which previously
 * appeared as intermediate terms in step-3 and step-4. A tensor is a
 * generalization of scalars (rank zero tensors), vectors (rank one
 * tensors), and matrices (rank two tensors), as well as higher dimensional
 * objects. The Tensor class requires two template arguments: the tensor
 * rank and tensor dimension. For example, here we use tensors of rank one
 * (vectors) with dimension <code>dim</code> (so they have <code>dim</code>
 * entries.) While this is a bit less flexible than using Vector, the
 * compiler can generate faster code when the length of the vector is known
 * at compile time. Additionally, specifying a Tensor of rank one and
 * dimension <code>dim</code> guarantees that the tensor will have the right
 * shape (since it is built into the type of the object itself), so the
 * compiler can catch most size-related mistakes for us.
 *   

 * 
 * Like in step-4, for compatibility with some compilers we explicitly
 * declare the default constructor:
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>, protected SolutionBase<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual Tensor<1, dim>
 *     gradient(const Point<dim> & p,
 *              const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * The actual definition of the values and gradients of the exact solution
 * class is according to their mathematical definition and does not need
 * much explanation.
 *   

 * 
 * The only thing that is worth mentioning is that if we access
 * elements of a base class that is template dependent (in this case
 * the elements of SolutionBase&lt;dim&gt;), then the C++ language
 * forces us to write <code>this-&gt;source_centers</code>, and
 * similarly for other members of the base class. C++ does not
 * require the <code>this-&gt;</code> qualification if the base
 * class is not template dependent. The reason why this is necessary
 * is complicated; C++ books will explain under the phrase
 * <i>two-stage (name) lookup</i>, and there is also a lengthy
 * description in the deal.II FAQs.
 * 
 * @code
 *   template <int dim>
 *   double Solution<dim>::value(const Point<dim> &p, const unsigned int) const
 *   {
 *     double return_value = 0;
 *     for (const auto &center : this->source_centers)
 *       {
 *         const Tensor<1, dim> x_minus_xi = p - center;
 *         return_value +=
 *           std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
 *       }
 * 
 *     return return_value;
 *   }
 * 
 * 
 * @endcode
 * 
 * Likewise, this is the computation of the gradient of the solution.  In
 * order to accumulate the gradient from the contributions of the
 * exponentials, we allocate an object <code>return_value</code> that
 * denotes the mathematical quantity of a tensor of rank <code>1</code> and
 * dimension <code>dim</code>. Its default constructor sets it to the vector
 * containing only zeroes, so we need not explicitly care for its
 * initialization.
 *   

 * 
 * Note that we could as well have taken the type of the object to be
 * Point&lt;dim&gt; instead of Tensor&lt;1,dim&gt;. Tensors of rank 1 and
 * points are almost exchangeable, and have only very slightly different
 * mathematical meanings. In fact, the Point&lt;dim&gt; class is derived
 * from the Tensor&lt;1,dim&gt; class, which makes up for their mutual
 * exchange ability. Their main difference is in what they logically mean:
 * points are points in space, such as the location at which we want to
 * evaluate a function (see the type of the first argument of this function
 * for example). On the other hand, tensors of rank 1 share the same
 * transformation properties, for example that they need to be rotated in a
 * certain way when we change the coordinate system; however, they do not
 * share the same connotation that points have and are only objects in a
 * more abstract space than the one spanned by the coordinate
 * directions. (In fact, gradients live in `reciprocal' space, since the
 * dimension of their components is not that of a length, but of one over
 * length).
 * 
 * @code
 *   template <int dim>
 *   Tensor<1, dim> Solution<dim>::gradient(const Point<dim> &p,
 *                                          const unsigned int) const
 *   {
 *     Tensor<1, dim> return_value;
 * 
 *     for (const auto &center : this->source_centers)
 *       {
 *         const Tensor<1, dim> x_minus_xi = p - center;
 * 
 * @endcode
 * 
 * For the gradient, note that its direction is along (x-x_i), so we
 * add up multiples of this distance vector, where the factor is given
 * by the exponentials.
 * 
 * @code
 *         return_value +=
 *           (-2. / (this->width * this->width) *
 *            std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) *
 *            x_minus_xi);
 *       }
 * 
 *     return return_value;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Besides the function that represents the exact solution, we also need a
 * function which we can use as right hand side when assembling the linear
 * system of discretized equations. This is accomplished using the following
 * class and the following definition of its function. Note that here we
 * only need the value of the function, not its gradients or higher
 * derivatives.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>, protected SolutionBase<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * The value of the right hand side is given by the negative Laplacian of
 * the solution plus the solution itself, since we wanted to solve
 * Helmholtz's equation:
 * 
 * @code
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> &p,
 *                                    const unsigned int) const
 *   {
 *     double return_value = 0;
 *     for (const auto &center : this->source_centers)
 *       {
 *         const Tensor<1, dim> x_minus_xi = p - center;
 * 
 * @endcode
 * 
 * The first contribution is the Laplacian:
 * 
 * @code
 *         return_value +=
 *           ((2. * dim -
 *             4. * x_minus_xi.norm_square() / (this->width * this->width)) /
 *            (this->width * this->width) *
 *            std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
 * @endcode
 * 
 * And the second is the solution itself:
 * 
 * @code
 *         return_value +=
 *           std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
 *       }
 * 
 *     return return_value;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheHelmholtzsolverclass"></a> 
 * <h3>The Helmholtz solver class</h3>
 * 

 * 
 * Then we need the class that does all the work. Except for its name, its
 * interface is mostly the same as in previous examples.
 *   

 * 
 * One of the differences is that we will use this class in several modes:
 * for different finite elements, as well as for adaptive and global
 * refinement. The decision whether global or adaptive refinement shall be
 * used is communicated to the constructor of this class through an
 * enumeration type declared at the top of the class. The constructor then
 * takes a finite element object and the refinement mode as arguments.
 *   

 * 
 * The rest of the member functions are as before except for the
 * <code>process_solution</code> function: After the solution has been
 * computed, we perform some analysis on it, such as computing the error in
 * various norms. To enable some output, it requires the number of the
 * refinement cycle, and consequently gets it as an argument.
 * 
 * @code
 *   template <int dim>
 *   class HelmholtzProblem
 *   {
 *   public:
 *     enum RefinementMode
 *     {
 *       global_refinement,
 *       adaptive_refinement
 *     };
 * 
 *     HelmholtzProblem(const FiniteElement<dim> &fe,
 *                      const RefinementMode      refinement_mode);
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void refine_grid();
 *     void process_solution(const unsigned int cycle);
 * 
 * @endcode
 * 
 * Now for the data elements of this class. Among the variables that we
 * have already used in previous examples, only the finite element object
 * differs: The finite elements which the objects of this class operate on
 * are passed to the constructor of this class. It has to store a pointer
 * to the finite element for the member functions to use. Now, for the
 * present class there is no big deal in that, but since we want to show
 * techniques rather than solutions in these programs, we will here point
 * out a problem that often occurs -- and of course the right solution as
 * well.
 *     

 * 
 * Consider the following situation that occurs in all the example
 * programs: we have a triangulation object, and we have a finite element
 * object, and we also have an object of type DoFHandler that uses both of
 * the first two. These three objects all have a lifetime that is rather
 * long compared to most other objects: they are basically set at the
 * beginning of the program or an outer loop, and they are destroyed at
 * the very end. The question is: can we guarantee that the two objects
 * which the DoFHandler uses, live at least as long as they are in use?
 * This means that the DoFHandler must have some kind of knowledge on the
 * destruction of the other objects.
 *     

 * 
 * We will show here how the library managed to find out that there are
 * still active references to an object and the object is still alive
 * from the point of view of a using object. Basically, the method is along
 * the following line: all objects that are subject to such potentially
 * dangerous pointers are derived from a class called Subscriptor. For
 * example, the Triangulation, DoFHandler, and a base class of the
 * FiniteElement class are derived from Subscriptor. This latter class
 * does not offer much functionality, but it has a built-in counter which
 * we can subscribe to, thus the name of the class. Whenever we initialize
 * a pointer to that object, we can increase its use counter, and when we
 * move away our pointer or do not need it any more, we decrease the
 * counter again. This way, we can always check how many objects still use
 * that object. Additionally, the class requires to know about a pointer
 * that it can use to tell the subscribing object about its invalidation.
 *     

 * 
 * If an object of a class that is derived from the Subscriptor class is
 * destroyed, it also has to call the destructor of the Subscriptor class.
 * In this destructor, we tell all the subscribing objects about the
 * invalidation of the object using the stored pointers. The same happens
 * when the object appears on the right hand side of a move expression,
 * i.e., it will no longer contain valid content after the operation. The
 * subscribing class is expected to check the value stored in its
 * corresponding pointer before trying to access the object subscribed to.
 *     

 * 
 * This is exactly what the SmartPointer class is doing. It basically acts
 * just like a pointer, i.e. it can be dereferenced, can be assigned to and
 * from other pointers, and so on. On top of that it uses the mechanism
 * described above to find out if the pointer this class is representing is
 * dangling when we try to dereference it. In that case an exception is
 * thrown.
 *     

 * 
 * In the present example program, we want to protect the finite element
 * object from the situation that for some reason the finite element
 * pointed to is destroyed while still in use. We therefore use a
 * SmartPointer to the finite element object; since the finite element
 * object is actually never changed in our computations, we pass a const
 * FiniteElement&lt;dim&gt; as template argument to the SmartPointer
 * class. Note that the pointer so declared is assigned at construction
 * time of the solve object, and destroyed upon destruction, so the lock
 * on the destruction of the finite element object extends throughout the
 * lifetime of this HelmholtzProblem object.
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     SmartPointer<const FiniteElement<dim>> fe;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 * @endcode
 * 
 * The second to last variable stores the refinement mode passed to the
 * constructor. Since it is only set in the constructor, we can declare
 * this variable constant, to avoid that someone sets it involuntarily
 * (e.g. in an `if'-statement where == was written as = by chance).
 * 
 * @code
 *     const RefinementMode refinement_mode;
 * 
 * @endcode
 * 
 * For each refinement level some data (like the number of cells, or the
 * L2 error of the numerical solution) will be generated and later
 * printed. The TableHandler can be used to collect all this data and to
 * output it at the end of the run as a table in a simple text or in LaTeX
 * format. Here we don't only use the TableHandler but we use the derived
 * class ConvergenceTable that additionally evaluates rates of
 * convergence:
 * 
 * @code
 *     ConvergenceTable convergence_table;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheHelmholtzProblemclassimplementation"></a> 
 * <h3>The HelmholtzProblem class implementation</h3>
 * 

 * 
 * 
 * <a name="HelmholtzProblemHelmholtzProblemconstructor"></a> 
 * <h4>HelmholtzProblem::HelmholtzProblem constructor</h4>
 * 

 * 
 * In the constructor of this class, we only set the variables passed as
 * arguments, and associate the DoF handler object with the triangulation
 * (which is empty at present, however).
 * 
 * @code
 *   template <int dim>
 *   HelmholtzProblem<dim>::HelmholtzProblem(const FiniteElement<dim> &fe,
 *                                           const RefinementMode refinement_mode)
 *     : dof_handler(triangulation)
 *     , fe(&fe)
 *     , refinement_mode(refinement_mode)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemsetup_system"></a> 
 * <h4>HelmholtzProblem::setup_system</h4>
 * 

 * 
 * The following function sets up the degrees of freedom, sizes of matrices
 * and vectors, etc. Most of its functionality has been showed in previous
 * examples, the only difference being the renumbering step immediately
 * after first distributing degrees of freedom.
 *   

 * 
 * Renumbering the degrees of freedom is not overly difficult, as long as
 * you use one of the algorithms included in the library. It requires only a
 * single line of code. Some more information on this can be found in
 * step-2.
 *   

 * 
 * Note, however, that when you renumber the degrees of freedom, you must do
 * so immediately after distributing them, since such things as hanging
 * nodes, the sparsity pattern etc. depend on the absolute numbers which are
 * altered by renumbering.
 *   

 * 
 * The reason why we introduce renumbering here is that it is a relatively
 * cheap operation but often has a beneficial effect: While the CG iteration
 * itself is independent of the actual ordering of degrees of freedom, we
 * will use SSOR as a preconditioner. SSOR goes through all degrees of
 * freedom and does some operations that depend on what happened before; the
 * SSOR operation is therefore not independent of the numbering of degrees
 * of freedom, and it is known that its performance improves by using
 * renumbering techniques. A little experiment shows that indeed, for
 * example, the number of CG iterations for the fifth refinement cycle of
 * adaptive refinement with the Q1 program used here is 40 without, but 36
 * with renumbering. Similar savings can generally be observed for all the
 * computations in this program.
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(*fe);
 *     DoFRenumbering::Cuthill_McKee(dof_handler);
 * 
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *     hanging_node_constraints.condense(dsp);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemassemble_system"></a> 
 * <h4>HelmholtzProblem::assemble_system</h4>
 * 

 * 
 * Assembling the system of equations for the problem at hand is mostly as
 * for the example programs before. However, some things have changed
 * anyway, so we comment on this function fairly extensively.
 *   

 * 
 * At the top of the function you will find the usual assortment of variable
 * declarations. Compared to previous programs, of importance is only that
 * we expect to solve problems also with bi-quadratic elements and therefore
 * have to use sufficiently accurate quadrature formula. In addition, we
 * need to compute integrals over faces, i.e. <code>dim-1</code> dimensional
 * objects. The declaration of a face quadrature formula is then
 * straightforward:
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim>     quadrature_formula(fe->degree + 1);
 *     QGauss<dim - 1> face_quadrature_formula(fe->degree + 1);
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Then we need objects which can evaluate the values, gradients, etc of
 * the shape functions at the quadrature points. While it seems that it
 * should be feasible to do it with one object for both domain and face
 * integrals, there is a subtle difference since the weights in the domain
 * integrals include the measure of the cell in the domain, while the face
 * integral quadrature requires the measure of the face in a
 * lower-dimensional manifold. Internally these two classes are rooted in
 * a common base class which does most of the work and offers the same
 * interface to both domain and interface integrals.
 *     

 * 
 * For the domain integrals in the bilinear form for Helmholtz's equation,
 * we need to compute the values and gradients, as well as the weights at
 * the quadrature points. Furthermore, we need the quadrature points on
 * the real cell (rather than on the unit cell) to evaluate the right hand
 * side function. The object we use to get at this information is the
 * FEValues class discussed previously.
 *     

 * 
 * For the face integrals, we only need the values of the shape functions,
 * as well as the weights. We also need the normal vectors and quadrature
 * points on the real cell since we want to determine the Neumann values
 * from the exact solution object (see below). The class that gives us
 * this information is called FEFaceValues:
 * 
 * @code
 *     FEValues<dim> fe_values(*fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values(*fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_quadrature_points |
 *                                        update_normal_vectors |
 *                                        update_JxW_values);
 * 
 * @endcode
 * 
 * Then we need some objects already known from previous examples: An
 * object denoting the right hand side function, its values at the
 * quadrature points on a cell, the cell matrix and right hand side, and
 * the indices of the degrees of freedom on a cell.
 *     

 * 
 * Note that the operations we will do with the right hand side object are
 * only querying data, never changing the object. We can therefore declare
 * it <code>const</code>:
 * 
 * @code
 *     const RightHandSide<dim> right_hand_side;
 *     std::vector<double>      rhs_values(n_q_points);
 * 
 * @endcode
 * 
 * Finally we define an object denoting the exact solution function. We
 * will use it to compute the Neumann values at the boundary from
 * it. Usually, one would of course do so using a separate object, in
 * particular since the exact solution is generally unknown while the
 * Neumann values are prescribed. We will, however, be a little bit lazy
 * and use what we already have in information. Real-life programs would
 * to go other ways here, of course.
 * 
 * @code
 *     Solution<dim> exact_solution;
 * 
 * @endcode
 * 
 * Now for the main loop over all cells. This is mostly unchanged from
 * previous examples, so we only comment on the things that have changed.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0.;
 *         cell_rhs    = 0.;
 * 
 *         fe_values.reinit(cell);
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                    rhs_values);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 * @endcode
 * 
 * The first thing that has changed is the bilinear form. It
 * now contains the additional term from the Helmholtz
 * equation:
 * 
 * @code
 *                 cell_matrix(i, j) +=
 *                   ((fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
 *                       fe_values.shape_grad(j, q_point)     // grad phi_j(x_q)
 *                     +                                      
 *                     fe_values.shape_value(i, q_point) *    // phi_i(x_q)
 *                       fe_values.shape_value(j, q_point)) * // phi_j(x_q)
 *                    fe_values.JxW(q_point));                // dx
 * 
 * 
 *               cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
 *                               rhs_values[q_point] *               // f(x_q)
 *                               fe_values.JxW(q_point));            // dx
 *             }
 * 
 * @endcode
 * 
 * Then there is that second term on the right hand side, the contour
 * integral. First we have to find out whether the intersection of the
 * faces of this cell with the boundary part Gamma2 is nonzero. To
 * this end, we loop over all faces and check whether its boundary
 * indicator equals <code>1</code>, which is the value that we have
 * assigned to that portions of the boundary composing Gamma2 in the
 * <code>run()</code> function further below. (The default value of
 * boundary indicators is <code>0</code>, so faces can only have an
 * indicator equal to <code>1</code> if we have explicitly set it.)
 * 
 * @code
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary() && (face->boundary_id() == 1))
 *             {
 * @endcode
 * 
 * If we came into here, then we have found an external face
 * belonging to Gamma2. Next, we have to compute the values of
 * the shape functions and the other quantities which we will
 * need for the computation of the contour integral. This is
 * done using the <code>reinit</code> function which we already
 * know from the FEValue class:
 * 
 * @code
 *               fe_face_values.reinit(cell, face);
 * 
 * @endcode
 * 
 * And we can then perform the integration by using a loop over
 * all quadrature points.
 *               

 * 
 * On each quadrature point, we first compute the value of the
 * normal derivative. We do so using the gradient of the exact
 * solution and the normal vector to the face at the present
 * quadrature point obtained from the
 * <code>fe_face_values</code> object. This is then used to
 * compute the additional contribution of this face to the right
 * hand side:
 * 
 * @code
 *               for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                    ++q_point)
 *                 {
 *                   const double neumann_value =
 *                     (exact_solution.gradient(
 *                        fe_face_values.quadrature_point(q_point)) *
 *                      fe_face_values.normal_vector(q_point));
 * 
 *                   for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                     cell_rhs(i) +=
 *                       (fe_face_values.shape_value(i, q_point) * // phi_i(x_q)
 *                        neumann_value *                          // g(x_q)
 *                        fe_face_values.JxW(q_point));            // dx
 *                 }
 *             }
 * 
 * @endcode
 * 
 * Now that we have the contributions of the present cell, we can
 * transfer it to the global matrix and right hand side vector, as in
 * the examples before:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               system_matrix.add(local_dof_indices[i],
 *                                 local_dof_indices[j],
 *                                 cell_matrix(i, j));
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *           }
 *       }
 * 
 * @endcode
 * 
 * Likewise, elimination and treatment of boundary values has been shown
 * previously.
 *     

 * 
 * We note, however that now the boundary indicator for which we
 * interpolate boundary values (denoted by the second parameter to
 * <code>interpolate_boundary_values</code>) does not represent the whole
 * boundary any more. Rather, it is that portion of the boundary which we
 * have not assigned another indicator (see below). The degrees of freedom
 * at the boundary that do not belong to Gamma1 are therefore excluded
 * from the interpolation of boundary values, just as we want.
 * 
 * @code
 *     hanging_node_constraints.condense(system_matrix);
 *     hanging_node_constraints.condense(system_rhs);
 * 
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Solution<dim>(),
 *                                              boundary_values);
 *     MatrixTools::apply_boundary_values(boundary_values,
 *                                        system_matrix,
 *                                        solution,
 *                                        system_rhs);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemsolve"></a> 
 * <h4>HelmholtzProblem::solve</h4>
 * 

 * 
 * Solving the system of equations is done in the same way as before:
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     hanging_node_constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemrefine_grid"></a> 
 * <h4>HelmholtzProblem::refine_grid</h4>
 * 

 * 
 * Now for the function doing grid refinement. Depending on the refinement
 * mode passed to the constructor, we do global or adaptive refinement.
 *   

 * 
 * Global refinement is simple, so there is not much to comment on.  In case
 * of adaptive refinement, we use the same functions and classes as in the
 * previous example program. Note that one could treat Neumann boundaries
 * differently than Dirichlet boundaries, and one should in fact do so here
 * since we have Neumann boundary conditions on part of the boundaries, but
 * since we don't have a function here that describes the Neumann values (we
 * only construct these values from the exact solution when assembling the
 * matrix), we omit this detail even though doing this in a strictly correct
 * way would not be hard to add.
 *   

 * 
 * At the end of the switch, we have a default case that looks slightly
 * strange: an <code>Assert</code> statement with a <code>false</code>
 * condition. Since the <code>Assert</code> macro raises an error whenever
 * the condition is false, this means that whenever we hit this statement
 * the program will be aborted. This in intentional: Right now we have only
 * implemented two refinement strategies (global and adaptive), but someone
 * might want to add a third strategy (for example adaptivity with a
 * different refinement criterion) and add a third member to the enumeration
 * that determines the refinement mode. If it weren't for the default case
 * of the switch statement, this function would simply run to its end
 * without doing anything. This is most likely not what was intended. One of
 * the defensive programming techniques that you will find all over the
 * deal.II library is therefore to always have default cases that abort, to
 * make sure that values not considered when listing the cases in the switch
 * statement are eventually caught, and forcing programmers to add code to
 * handle them. We will use this same technique in other places further down
 * as well.
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::refine_grid()
 *   {
 *     switch (refinement_mode)
 *       {
 *         case global_refinement:
 *           {
 *             triangulation.refine_global(1);
 *             break;
 *           }
 * 
 *         case adaptive_refinement:
 *           {
 *             Vector<float> estimated_error_per_cell(
 *               triangulation.n_active_cells());
 * 
 *             KellyErrorEstimator<dim>::estimate(
 *               dof_handler,
 *               QGauss<dim - 1>(fe->degree + 1),
 *               std::map<types::boundary_id, const Function<dim> *>(),
 *               solution,
 *               estimated_error_per_cell);
 * 
 *             GridRefinement::refine_and_coarsen_fixed_number(
 *               triangulation, estimated_error_per_cell, 0.3, 0.03);
 * 
 *             triangulation.execute_coarsening_and_refinement();
 * 
 *             break;
 *           }
 * 
 *         default:
 *           {
 *             Assert(false, ExcNotImplemented());
 *           }
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemprocess_solution"></a> 
 * <h4>HelmholtzProblem::process_solution</h4>
 * 

 * 
 * Finally we want to process the solution after it has been computed. For
 * this, we integrate the error in various (semi-)norms, and we generate
 * tables that will later be used to display the convergence against the
 * continuous solution in a nice format.
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::process_solution(const unsigned int cycle)
 *   {
 * @endcode
 * 
 * Our first task is to compute error norms. In order to integrate the
 * difference between computed numerical solution and the continuous
 * solution (described by the Solution class defined at the top of this
 * file), we first need a vector that will hold the norm of the error on
 * each cell. Since accuracy with 16 digits is not so important for these
 * quantities, we save some memory by using <code>float</code> instead of
 * <code>double</code> values.
 *     

 * 
 * The next step is to use a function from the library which computes the
 * error in the L2 norm on each cell.  We have to pass it the DoF handler
 * object, the vector holding the nodal values of the numerical solution,
 * the continuous solution as a function object, the vector into which it
 * shall place the norm of the error on each cell, a quadrature rule by
 * which this norm shall be computed, and the type of norm to be
 * used. Here, we use a Gauss formula with three points in each space
 * direction, and compute the L2 norm.
 *     

 * 
 * Finally, we want to get the global L2 norm. This can of course be
 * obtained by summing the squares of the norms on each cell, and taking
 * the square root of that value. This is equivalent to taking the l2
 * (lower case <code>l</code>) norm of the vector of norms on each cell:
 * 
 * @code
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe->degree + 1),
 *                                       VectorTools::L2_norm);
 *     const double L2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 * @endcode
 * 
 * By same procedure we get the H1 semi-norm. We re-use the
 * <code>difference_per_cell</code> vector since it is no longer used
 * after computing the <code>L2_error</code> variable above. The global
 * $H^1$ semi-norm error is then computed by taking the sum of squares
 * of the errors on each individual cell, and then the square root of
 * it -- an operation that is conveniently performed by
 * VectorTools::compute_global_error.
 * 
 * @code
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe->degree + 1),
 *                                       VectorTools::H1_seminorm);
 *     const double H1_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::H1_seminorm);
 * 
 * @endcode
 * 
 * Finally, we compute the maximum norm. Of course, we can't actually
 * compute the true maximum of the error over *all* points in the domain,
 * but only the maximum over a finite set of evaluation points that, for
 * convenience, we will still call "quadrature points" and represent by
 * an object of type Quadrature even though we do not actually perform any
 * integration.
 *     

 * 
 * There is then the question of what points precisely we want to evaluate
 * at. It turns out that the result we get depends quite sensitively on the
 * "quadrature" points being used. There is also the issue of
 * superconvergence: Finite element solutions are, on some meshes and for
 * polynomial degrees $k\ge 2$, particularly accurate at the node points as
 * well as at Gauss-Lobatto points, much more accurate than at randomly
 * chosen points. (See
 * @cite Li2019 and the discussion and references in Section 1.2 for more
 * information on this.) In other words, if we are interested in finding
 * the largest difference $u(\mathbf x)-u_h(\mathbf x)$, then we ought to
 * look at points $\mathbf x$ that are specifically not of this "special"
 * kind of points and we should specifically not use
 * `QGauss(fe->degree+1)` to define where we evaluate. Rather, we use a
 * special quadrature rule that is obtained by iterating the trapezoidal
 * rule by the degree of the finite element times two plus one in each space
 * direction. Note that the constructor of the QIterated class takes a
 * one-dimensional quadrature rule and a number that tells it how often it
 * shall repeat this rule in each space direction.
 *     

 * 
 * Using this special quadrature rule, we can then try to find the maximal
 * error on each cell. Finally, we compute the global L infinity error
 * from the L infinity errors on each cell with a call to
 * VectorTools::compute_global_error.
 * 
 * @code
 *     const QTrapezoid<1>  q_trapez;
 *     const QIterated<dim> q_iterated(q_trapez, fe->degree * 2 + 1);
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       difference_per_cell,
 *                                       q_iterated,
 *                                       VectorTools::Linfty_norm);
 *     const double Linfty_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::Linfty_norm);
 * 
 * @endcode
 * 
 * After all these errors have been computed, we finally write some
 * output. In addition, we add the important data to the TableHandler by
 * specifying the key of the column and the value.  Note that it is not
 * necessary to define column keys beforehand -- it is sufficient to just
 * add values, and columns will be introduced into the table in the order
 * values are added the first time.
 * 
 * @code
 *     const unsigned int n_active_cells = triangulation.n_active_cells();
 *     const unsigned int n_dofs         = dof_handler.n_dofs();
 * 
 *     std::cout << "Cycle " << cycle << ':' << std::endl
 *               << "   Number of active cells:       " << n_active_cells
 *               << std::endl
 *               << "   Number of degrees of freedom: " << n_dofs << std::endl;
 * 
 *     convergence_table.add_value("cycle", cycle);
 *     convergence_table.add_value("cells", n_active_cells);
 *     convergence_table.add_value("dofs", n_dofs);
 *     convergence_table.add_value("L2", L2_error);
 *     convergence_table.add_value("H1", H1_error);
 *     convergence_table.add_value("Linfty", Linfty_error);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemrun"></a> 
 * <h4>HelmholtzProblem::run</h4>
 * 

 * 
 * As in previous example programs, the <code>run</code> function controls
 * the flow of execution. The basic layout is as in previous examples: an
 * outer loop over successively refined grids, and in this loop first
 * problem setup, assembling the linear system, solution, and
 * post-processing.
 *   

 * 
 * The first task in the main loop is creation and refinement of grids. This
 * is as in previous examples, with the only difference that we want to have
 * part of the boundary marked as Neumann type, rather than Dirichlet.
 *   

 * 
 * For this, we will use the following convention: Faces belonging to Gamma1
 * will have the boundary indicator <code>0</code> (which is the default, so
 * we don't have to set it explicitly), and faces belonging to Gamma2 will
 * use <code>1</code> as boundary indicator.  To set these values, we loop
 * over all cells, then over all faces of a given cell, check whether it is
 * part of the boundary that we want to denote by Gamma2, and if so set its
 * boundary indicator to <code>1</code>. For the present program, we
 * consider the left and bottom boundaries as Gamma2. We determine whether a
 * face is part of that boundary by asking whether the x or y coordinates
 * (i.e. vector components 0 and 1) of the midpoint of a face equals -1, up
 * to some small wiggle room that we have to give since it is instable to
 * compare floating point numbers that are subject to round off in
 * intermediate computations.
 *   

 * 
 * It is worth noting that we have to loop over all cells here, not only the
 * active ones. The reason is that upon refinement, newly created faces
 * inherit the boundary indicator of their parent face. If we now only set
 * the boundary indicator for active faces, coarsen some cells and refine
 * them later on, they will again have the boundary indicator of the parent
 * cell which we have not modified, instead of the one we
 * intended. Consequently, we have to change the boundary indicators of
 * faces of all cells on Gamma2, whether they are active or not.
 * Alternatively, we could of course have done this job on the coarsest mesh
 * (i.e. before the first refinement step) and refined the mesh only after
 * that.
 * 
 * @code
 *   template <int dim>
 *   void HelmholtzProblem<dim>::run()
 *   {
 *     const unsigned int n_cycles =
 *       (refinement_mode == global_refinement) ? 5 : 9;
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, -1., 1.);
 *             triangulation.refine_global(3);
 * 
 *             for (const auto &cell : triangulation.cell_iterators())
 *               for (const auto &face : cell->face_iterators())
 *                 {
 *                   const auto center = face->center();
 *                   if ((std::fabs(center(0) - (-1.0)) < 1e-12) ||
 *                       (std::fabs(center(1) - (-1.0)) < 1e-12))
 *                     face->set_boundary_id(1);
 *                 }
 *           }
 *         else
 *           refine_grid();
 * 
 * 
 * @endcode
 * 
 * The next steps are already known from previous examples. This is
 * mostly the basic set-up of every finite element program:
 * 
 * @code
 *         setup_system();
 * 
 *         assemble_system();
 *         solve();
 * 
 * @endcode
 * 
 * The last step in this chain of function calls is usually the
 * evaluation of the computed solution for the quantities one is
 * interested in. This is done in the following function. Since the
 * function generates output that indicates the number of the present
 * refinement step, we pass this number as an argument.
 * 
 * @code
 *         process_solution(cycle);
 *       }
 * 
 * @endcode
 * 
 * 
 * <a name="Outputofgraphicaldata"></a> 
 * <h5>Output of graphical data</h5>
 * 

 * 
 * After the last iteration we output the solution on the finest
 * grid. This is done using the following sequence of statements which we
 * have already discussed in previous examples. The first step is to
 * generate a suitable filename (called <code>vtk_filename</code> here,
 * since we want to output data in VTK format; we add the prefix to
 * distinguish the filename from that used for other output files further
 * down below). Here, we augment the name by the mesh refinement
 * algorithm, and as above we make sure that we abort the program if
 * another refinement method is added and not handled by the following
 * switch statement:
 * 
 * @code
 *     std::string vtk_filename;
 *     switch (refinement_mode)
 *       {
 *         case global_refinement:
 *           vtk_filename = "solution-global";
 *           break;
 *         case adaptive_refinement:
 *           vtk_filename = "solution-adaptive";
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 * @endcode
 * 
 * We augment the filename by a postfix denoting the finite element which
 * we have used in the computation. To this end, the finite element base
 * class stores the maximal polynomial degree of shape functions in each
 * coordinate variable as a variable <code>degree</code>, and we use for
 * the switch statement (note that the polynomial degree of bilinear shape
 * functions is really 2, since they contain the term <code>x*y</code>;
 * however, the polynomial degree in each coordinate variable is still
 * only 1). We again use the same defensive programming technique to
 * safeguard against the case that the polynomial degree has an unexpected
 * value, using the <code>Assert (false, ExcNotImplemented())</code> idiom
 * in the default branch of the switch statement:
 * 
 * @code
 *     switch (fe->degree)
 *       {
 *         case 1:
 *           vtk_filename += "-q1";
 *           break;
 *         case 2:
 *           vtk_filename += "-q2";
 *           break;
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 * @endcode
 * 
 * Once we have the base name for the output file, we add an extension
 * appropriate for VTK output, open a file, and add the solution vector to
 * the object that will do the actual output:
 * 
 * @code
 *     vtk_filename += ".vtk";
 *     std::ofstream output(vtk_filename);
 * 
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 * 
 * @endcode
 * 
 * Now building the intermediate format as before is the next step. We
 * introduce one more feature of deal.II here. The background is the
 * following: in some of the runs of this function, we have used
 * biquadratic finite elements. However, since almost all output formats
 * only support bilinear data, the data is written only bilinear, and
 * information is consequently lost.  Of course, we can't change the
 * format in which graphic programs accept their inputs, but we can write
 * the data differently such that we more closely resemble the information
 * available in the quadratic approximation. We can, for example, write
 * each cell as four sub-cells with bilinear data each, such that we have
 * nine data points for each cell in the triangulation. The graphic
 * programs will, of course, display this data still only bilinear, but at
 * least we have given some more of the information we have.
 *     

 * 
 * In order to allow writing more than one sub-cell per actual cell, the
 * <code>build_patches</code> function accepts a parameter (the default is
 * <code>1</code>, which is why you haven't seen this parameter in
 * previous examples). This parameter denotes into how many sub-cells per
 * space direction each cell shall be subdivided for output. For example,
 * if you give <code>2</code>, this leads to 4 cells in 2D and 8 cells in
 * 3D. For quadratic elements, two sub-cells per space direction is
 * obviously the right choice, so this is what we choose. In general, for
 * elements of polynomial order <code>q</code>, we use <code>q</code>
 * subdivisions, and the order of the elements is determined in the same
 * way as above.
 *     

 * 
 * With the intermediate format so generated, we can then actually write
 * the graphical output:
 * 
 * @code
 *     data_out.build_patches(fe->degree);
 *     data_out.write_vtk(output);
 * 
 * @endcode
 * 
 * 
 * <a name="Outputofconvergencetables"></a> 
 * <h5>Output of convergence tables</h5>
 * 

 * 
 * After graphical output, we would also like to generate tables from the
 * error computations we have done in
 * <code>process_solution</code>. There, we have filled a table object
 * with the number of cells for each refinement step as well as the errors
 * in different norms.
 * 

 * 
 * For a nicer textual output of this data, one may want to set the
 * precision with which the values will be written upon output. We use 3
 * digits for this, which is usually sufficient for error norms. By
 * default, data is written in fixed point notation. However, for columns
 * one would like to see in scientific notation another function call sets
 * the <code>scientific_flag</code> to <code>true</code>, leading to
 * floating point representation of numbers.
 * 
 * @code
 *     convergence_table.set_precision("L2", 3);
 *     convergence_table.set_precision("H1", 3);
 *     convergence_table.set_precision("Linfty", 3);
 * 
 *     convergence_table.set_scientific("L2", true);
 *     convergence_table.set_scientific("H1", true);
 *     convergence_table.set_scientific("Linfty", true);
 * 
 * @endcode
 * 
 * For the output of a table into a LaTeX file, the default captions of
 * the columns are the keys given as argument to the
 * <code>add_value</code> functions. To have TeX captions that differ from
 * the default ones you can specify them by the following function calls.
 * Note, that `\\' is reduced to `\' by the compiler such that the real
 * TeX caption is, e.g., `$L^\infty$-error'.
 * 
 * @code
 *     convergence_table.set_tex_caption("cells", "\\# cells");
 *     convergence_table.set_tex_caption("dofs", "\\# dofs");
 *     convergence_table.set_tex_caption("L2", "L^2-error");
 *     convergence_table.set_tex_caption("H1", "H^1-error");
 *     convergence_table.set_tex_caption("Linfty", "L^\\infty-error");
 * 
 * @endcode
 * 
 * Finally, the default LaTeX format for each column of the table is `c'
 * (centered). To specify a different (e.g. `right') one, the following
 * function may be used:
 * 
 * @code
 *     convergence_table.set_tex_format("cells", "r");
 *     convergence_table.set_tex_format("dofs", "r");
 * 
 * @endcode
 * 
 * After this, we can finally write the table to the standard output
 * stream <code>std::cout</code> (after one extra empty line, to make
 * things look prettier). Note, that the output in text format is quite
 * simple and that captions may not be printed directly above the specific
 * columns.
 * 
 * @code
 *     std::cout << std::endl;
 *     convergence_table.write_text(std::cout);
 * 
 * @endcode
 * 
 * The table can also be written into a LaTeX file.  The (nicely)
 * formatted table can be viewed at after calling `latex filename' and
 * e.g. `xdvi filename', where filename is the name of the file to which
 * we will write output now. We construct the file name in the same way as
 * before, but with a different prefix "error":
 * 
 * @code
 *     std::string error_filename = "error";
 *     switch (refinement_mode)
 *       {
 *         case global_refinement:
 *           error_filename += "-global";
 *           break;
 *         case adaptive_refinement:
 *           error_filename += "-adaptive";
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     switch (fe->degree)
 *       {
 *         case 1:
 *           error_filename += "-q1";
 *           break;
 *         case 2:
 *           error_filename += "-q2";
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     error_filename += ".tex";
 *     std::ofstream error_table_file(error_filename);
 * 
 *     convergence_table.write_tex(error_table_file);
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Furthertablemanipulations"></a> 
 * <h5>Further table manipulations</h5>
 * 

 * 
 * In case of global refinement, it might be of interest to also output
 * the convergence rates. This may be done by the functionality the
 * ConvergenceTable offers over the regular TableHandler. However, we do
 * it only for global refinement, since for adaptive refinement the
 * determination of something like an order of convergence is somewhat
 * more involved. While we are at it, we also show a few other things that
 * can be done with tables.
 * 
 * @code
 *     if (refinement_mode == global_refinement)
 *       {
 * @endcode
 * 
 * The first thing is that one can group individual columns together
 * to form so-called super columns. Essentially, the columns remain
 * the same, but the ones that were grouped together will get a
 * caption running across all columns in a group. For example, let's
 * merge the "cycle" and "cells" columns into a super column named "n
 * cells":
 * 
 * @code
 *         convergence_table.add_column_to_supercolumn("cycle", "n cells");
 *         convergence_table.add_column_to_supercolumn("cells", "n cells");
 * 
 * @endcode
 * 
 * Next, it isn't necessary to always output all columns, or in the
 * order in which they were originally added during the run.
 * Selecting and re-ordering the columns works as follows (note that
 * this includes super columns):
 * 
 * @code
 *         std::vector<std::string> new_order;
 *         new_order.emplace_back("n cells");
 *         new_order.emplace_back("H1");
 *         new_order.emplace_back("L2");
 *         convergence_table.set_column_order(new_order);
 * 
 * @endcode
 * 
 * For everything that happened to the ConvergenceTable until this
 * point, it would have been sufficient to use a simple
 * TableHandler. Indeed, the ConvergenceTable is derived from the
 * TableHandler but it offers the additional functionality of
 * automatically evaluating convergence rates. For example, here is
 * how we can let the table compute reduction and convergence rates
 * (convergence rates are the binary logarithm of the reduction rate):
 * 
 * @code
 *         convergence_table.evaluate_convergence_rates(
 *           "L2", ConvergenceTable::reduction_rate);
 *         convergence_table.evaluate_convergence_rates(
 *           "L2", ConvergenceTable::reduction_rate_log2);
 *         convergence_table.evaluate_convergence_rates(
 *           "H1", ConvergenceTable::reduction_rate);
 *         convergence_table.evaluate_convergence_rates(
 *           "H1", ConvergenceTable::reduction_rate_log2);
 * @endcode
 * 
 * Each of these function calls produces an additional column that is
 * merged with the original column (in our example the `L2' and the
 * `H1' column) to a supercolumn.
 * 

 * 
 * Finally, we want to write this convergence chart again, first to
 * the screen and then, in LaTeX format, to disk. The filename is
 * again constructed as above.
 * 
 * @code
 *         std::cout << std::endl;
 *         convergence_table.write_text(std::cout);
 * 
 *         std::string conv_filename = "convergence";
 *         switch (refinement_mode)
 *           {
 *             case global_refinement:
 *               conv_filename += "-global";
 *               break;
 *             case adaptive_refinement:
 *               conv_filename += "-adaptive";
 *               break;
 *             default:
 *               Assert(false, ExcNotImplemented());
 *           }
 *         switch (fe->degree)
 *           {
 *             case 1:
 *               conv_filename += "-q1";
 *               break;
 *             case 2:
 *               conv_filename += "-q2";
 *               break;
 *             default:
 *               Assert(false, ExcNotImplemented());
 *           }
 *         conv_filename += ".tex";
 * 
 *         std::ofstream table_file(conv_filename);
 *         convergence_table.write_tex(table_file);
 *       }
 *   }
 * 
 * @endcode
 * 
 * The final step before going to <code>main()</code> is then to close the
 * namespace <code>Step7</code> into which we have put everything we needed
 * for this program:
 * 
 * @code
 * } // namespace Step7
 * 
 * @endcode
 * 
 * 
 * <a name="Mainfunction"></a> 
 * <h3>Main function</h3>
 * 

 * 
 * The main function is mostly as before. The only difference is that we solve
 * three times, once for Q1 and adaptive refinement, once for Q1 elements and
 * global refinement, and once for Q2 elements and global refinement.
 * 

 * 
 * Since we instantiate several template classes below for two space
 * dimensions, we make this more generic by declaring a constant at the
 * beginning of the function denoting the number of space dimensions. If you
 * want to run the program in 1d or 2d, you will then only have to change this
 * one instance, rather than all uses below:
 * 
 * @code
 * int main()
 * {
 *   const unsigned int dim = 2;
 * 
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step7;
 * 
 * @endcode
 * 
 * Now for the three calls to the main class. Each call is blocked into
 * curly braces in order to destroy the respective objects (i.e. the
 * finite element and the HelmholtzProblem object) at the end of the
 * block and before we go to the next run. This avoids conflicts with
 * variable names, and also makes sure that memory is released
 * immediately after one of the three runs has finished, and not only at
 * the end of the <code>try</code> block.
 * 
 * @code
 *       {
 *         std::cout << "Solving with Q1 elements, adaptive refinement"
 *                   << std::endl
 *                   << "============================================="
 *                   << std::endl
 *                   << std::endl;
 * 
 *         FE_Q<dim>             fe(1);
 *         HelmholtzProblem<dim> helmholtz_problem_2d(
 *           fe, HelmholtzProblem<dim>::adaptive_refinement);
 * 
 *         helmholtz_problem_2d.run();
 * 
 *         std::cout << std::endl;
 *       }
 * 
 *       {
 *         std::cout << "Solving with Q1 elements, global refinement" << std::endl
 *                   << "===========================================" << std::endl
 *                   << std::endl;
 * 
 *         FE_Q<dim>             fe(1);
 *         HelmholtzProblem<dim> helmholtz_problem_2d(
 *           fe, HelmholtzProblem<dim>::global_refinement);
 * 
 *         helmholtz_problem_2d.run();
 * 
 *         std::cout << std::endl;
 *       }
 * 
 *       {
 *         std::cout << "Solving with Q2 elements, global refinement" << std::endl
 *                   << "===========================================" << std::endl
 *                   << std::endl;
 * 
 *         FE_Q<dim>             fe(2);
 *         HelmholtzProblem<dim> helmholtz_problem_2d(
 *           fe, HelmholtzProblem<dim>::global_refinement);
 * 
 *         helmholtz_problem_2d.run();
 * 
 *         std::cout << std::endl;
 *       }
 *       {
 *         std::cout << "Solving with Q2 elements, adaptive refinement"
 *                   << std::endl
 *                   << "===========================================" << std::endl
 *                   << std::endl;
 * 
 *         FE_Q<dim>             fe(2);
 *         HelmholtzProblem<dim> helmholtz_problem_2d(
 *           fe, HelmholtzProblem<dim>::adaptive_refinement);
 * 
 *         helmholtz_problem_2d.run();
 * 
 *         std::cout << std::endl;
 *       }
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



The program generates two kinds of output. The first are the output
files <code>solution-adaptive-q1.vtk</code>,
<code>solution-global-q1.vtk</code>, and
<code>solution-global-q2.vtk</code>. We show the latter in a 3d view
here:


<img src="https://www.dealii.org/images/steps/developer/step-7.solution.png" alt="">




Secondly, the program writes tables not only to disk, but also to the
screen while running. The output looks like the following (recall that
columns labeled as "<code>H1</code>" actually show the $H^1$ <i>semi-</i>norm
of the error, not the full $H^1$ norm):


@code
examples/\step-7> make run
Solving with Q1 elements, adaptive refinement
=============================================

Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 81
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 157
Cycle 2:
   Number of active cells:       280
   Number of degrees of freedom: 341
Cycle 3:
   Number of active cells:       577
   Number of degrees of freedom: 690
Cycle 4:
   Number of active cells:       1099
   Number of degrees of freedom: 1264
Cycle 5:
   Number of active cells:       2191
   Number of degrees of freedom: 2452
Cycle 6:
   Number of active cells:       4165
   Number of degrees of freedom: 4510
Cycle 7:
   Number of active cells:       7915
   Number of degrees of freedom: 8440
Cycle 8:
   Number of active cells:       15196
   Number of degrees of freedom: 15912

cycle cells dofs     L2        H1      Linfty
    0    64    81 1.840e+00 2.858e+00 1.835e+00
    1   124   157 5.190e-02 1.200e+00 1.344e-01
    2   280   341 1.439e-02 7.892e-01 7.554e-02
    3   577   690 8.627e-03 5.061e-01 2.805e-02
    4  1099  1264 3.217e-03 3.030e-01 1.073e-02
    5  2191  2452 1.445e-03 2.097e-01 5.073e-03
    6  4165  4510 8.387e-04 1.460e-01 2.013e-03
    7  7915  8440 7.051e-04 1.053e-01 1.804e-03
    8 15196 15912 2.774e-04 7.463e-02 6.911e-04

Solving with Q1 elements, global refinement
===========================================

Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 81
Cycle 1:
   Number of active cells:       256
   Number of degrees of freedom: 289
Cycle 2:
   Number of active cells:       1024
   Number of degrees of freedom: 1089
Cycle 3:
   Number of active cells:       4096
   Number of degrees of freedom: 4225
Cycle 4:
   Number of active cells:       16384
   Number of degrees of freedom: 16641

cycle cells dofs     L2        H1      Linfty
    0    64    81 1.840e+00 2.858e+00 1.835e+00
    1   256   289 3.570e-02 1.199e+00 1.307e-01
    2  1024  1089 1.192e-02 7.565e-01 7.168e-02
    3  4096  4225 3.047e-03 3.823e-01 2.128e-02
    4 16384 16641 7.660e-04 1.917e-01 5.554e-03

n cells         H1                   L2
0    64 2.858e+00    -    - 1.840e+00     -    -
1   256 1.199e+00 2.38 1.25 3.570e-02 51.54 5.69
2  1024 7.565e-01 1.58 0.66 1.192e-02  2.99 1.58
3  4096 3.823e-01 1.98 0.98 3.047e-03  3.91 1.97
4 16384 1.917e-01 1.99 1.00 7.660e-04  3.98 1.99

Solving with Q2 elements, global refinement
===========================================

Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 289
Cycle 1:
   Number of active cells:       256
   Number of degrees of freedom: 1089
Cycle 2:
   Number of active cells:       1024
   Number of degrees of freedom: 4225
Cycle 3:
   Number of active cells:       4096
   Number of degrees of freedom: 16641
Cycle 4:
   Number of active cells:       16384
   Number of degrees of freedom: 66049

cycle cells dofs     L2        H1      Linfty
    0    64   289 1.606e-01 1.278e+00 3.029e-01
    1   256  1089 7.638e-03 5.248e-01 4.816e-02
    2  1024  4225 8.601e-04 1.086e-01 4.827e-03
    3  4096 16641 1.107e-04 2.756e-02 7.802e-04
    4 16384 66049 1.393e-05 6.915e-03 9.971e-05

n cells         H1                   L2
0    64 1.278e+00    -    - 1.606e-01     -    -
1   256 5.248e-01 2.43 1.28 7.638e-03 21.03 4.39
2  1024 1.086e-01 4.83 2.27 8.601e-04  8.88 3.15
3  4096 2.756e-02 3.94 1.98 1.107e-04  7.77 2.96
4 16384 6.915e-03 3.99 1.99 1.393e-05  7.94 2.99

Solving with Q2 elements, adaptive refinement
===========================================

Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 289
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 577
Cycle 2:
   Number of active cells:       289
   Number of degrees of freedom: 1353
Cycle 3:
   Number of active cells:       547
   Number of degrees of freedom: 2531
Cycle 4:
   Number of active cells:       1057
   Number of degrees of freedom: 4919
Cycle 5:
   Number of active cells:       2059
   Number of degrees of freedom: 9223
Cycle 6:
   Number of active cells:       3913
   Number of degrees of freedom: 17887
Cycle 7:
   Number of active cells:       7441
   Number of degrees of freedom: 33807
Cycle 8:
   Number of active cells:       14212
   Number of degrees of freedom: 64731

cycle cells dofs     L2        H1      Linfty
    0    64   289 1.606e-01 1.278e+00 3.029e-01
    1   124   577 7.891e-03 5.256e-01 4.852e-02
    2   289  1353 1.070e-03 1.155e-01 4.868e-03
    3   547  2531 5.962e-04 5.101e-02 1.876e-03
    4  1057  4919 1.977e-04 3.094e-02 7.923e-04
    5  2059  9223 7.738e-05 1.974e-02 7.270e-04
    6  3913 17887 2.925e-05 8.772e-03 1.463e-04
    7  7441 33807 1.024e-05 4.121e-03 8.567e-05
    8 14212 64731 3.761e-06 2.108e-03 2.167e-05
@endcode


One can see the error reduction upon grid refinement, and for the
cases where global refinement was performed, also the convergence
rates can be seen. The linear and quadratic convergence rates of Q1
and Q2 elements in the $H^1$ semi-norm can clearly be seen, as
are the quadratic and cubic rates in the $L_2$ norm.




Finally, the program also generated LaTeX versions of the tables (not shown
here) that is written into a file in a way so that it could be
copy-pasted into a LaTeX document.


<a name="Whenistheerrorsmall"></a><h4> When is the error "small"? </h4>


What we showed above is how to determine the size of the error
$\|u-u_h\|$ in a number of different norms. We did this primarily
because we were interested in testing that our solutions *converge*.
But from an engineering perspective, the question is often more
practical: How fine do I have to make my mesh so that the error is
"small enough"? In other words, if in the table above the $H^1$
semi-norm has been reduced to `4.121e-03`, is this good enough for me
to sign the blueprint and declare that our numerical simulation showed
that the bridge is strong enough?

In practice, we are rarely in this situation because I can not
typically compare the numerical solution $u_h$ against the exact
solution $u$ in situations that matter -- if I knew $u$, I would not
have to compute $u_h$. But even if I could, the question to ask in
general is then: `4.121e-03` *what*? The solution will have physical
units, say kg-times-meter-squared, and I'm integrating a function with
units square of the above over the domain, and then take the square
root. So if the domain is two-dimensional, the units of
$\|u-u_h\|_{L_2}$ are kg-times-meter-cubed. The question is then: Is
$4.121\times 10^{-3}$ kg-times-meter-cubed small? That depends on what
you're trying to simulate: If you're an astronomer used to masses
measured in solar masses and distances in light years, then yes, this
is a fantastically small number. But if you're doing atomic physics,
then no: That's not small, and your error is most certainly not
sufficiently small; you need a finer mesh.

In other words, when we look at these sorts of numbers, we generally
need to compare against a "scale". One way to do that is to not look
at the *absolute* error $\|u-u_h\|$ in whatever norm, but at the
*relative* error $\|u-u_h\|/\|u\|$. If this ratio is $10^{-5}$, then
you know that *on average*, the difference between $u$ and $u_h$ is
0.001 per cent -- probably small enough for engineering purposes.

How do we compute $\|u\|$? We just need to do an integration loop over
all cells, quadrature points on these cells, and then sum things up
and take the square root at the end. But there is a simpler way often
used: You can call
@code
    Vector<double> zero_vector (dof_handler.n_dofs());
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      zero_vector,
                                      Solution<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
@endcode
which computes $\|u-0\|_{L_2}$. Alternatively, if you're particularly
lazy and don't feel like creating the `zero_vector`, you could use
that if the mesh is not too coarse, then $\|u\| \approx \|u_h\|$, and
we can compute $\|u\| \approx \|u_h\|=\|0-u_h\|$ by calling
@code
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      ZeroFunction<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
@endcode
In both cases, one then only has to combine the vector of cellwise
norms into one global norm as we already do in the program, by calling
@code
    const double L2_norm =
      VectorTools::compute_global_error(triangulation,
                                        norm_per_cell,
                                        VectorTools::L2_norm);
@endcode



<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="HigherOrderElements"></a><h4> Higher Order Elements </h4>


Go ahead and run the program with higher order elements ($Q_3$, $Q_4$, ...). You
will notice that assertions in several parts of the code will trigger (for
example in the generation of the filename for the data output). You might have to address these,
but it should not be very hard to get the program to work!

<a name="ConvergenceComparison"></a><h4> Convergence Comparison </h4>


Is Q1 or Q2 better? What about adaptive versus global refinement? A (somewhat
unfair but typical) metric to compare them, is to look at the error as a
function of the number of unknowns.

To see this, create a plot in log-log style with the number of unknowns on the
$x$ axis and the $L_2$ error on the $y$ axis. You can add reference lines for
$h^2=N^{-1}$ and $h^3=N^{-3/2}$ and check that global and adaptive refinement
follow those. If one makes the (not completely unreasonable)
assumption that with a good linear solver, the computational effort is
proportional to the number of unknowns $N$, then it is clear that an
error reduction of ${\cal O}(N^{-3/2})$ is substantially better than a
reduction of the form ${\cal O}(N^{-1})$: That is, that adaptive
refinement gives us the desired error level with less computational
work than if we used global refinement. This is not a particularly
surprising conclusion, but it's worth checking these sorts of
assumptions in practice.

Of course, a fairer comparison would be to plot runtime (switch to release
mode first!) instead of number of unknowns on the $x$ axis. If you
plotted run time against the number of unknowns by timing each
refinement step (e.g., using the Timer class), you will notice that
the linear solver is not perfect -- its run time grows faster than
proportional to the linear system size -- and picking a better
linear solver might be appropriate for this kind of comparison.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-7.cc"
*/
