/**
@page step_14 The step-14 tutorial program
This tutorial depends on step-13.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Themaths">The maths</a>
        <li><a href="#Thesoftware">The software</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Evaluatingthesolution">Evaluating the solution</a>
      <ul>
        <li><a href="#TheEvaluationBaseclass">The EvaluationBase class</a>
        <li><a href="#ThePointValueEvaluationclass">The PointValueEvaluation class</a>
        <li><a href="#ThePointXDerivativeEvaluationclass">The PointXDerivativeEvaluation class</a>
        <li><a href="#TheGridOutputclass">The GridOutput class</a>
      </ul>
        <li><a href="#TheLaplacesolverclasses">The Laplace solver classes</a>
      <ul>
        <li><a href="#TheLaplacesolverbaseclass">The Laplace solver base class</a>
        <li><a href="#TheLaplaceSolverclass">The Laplace Solver class</a>
        <li><a href="#ThePrimalSolverclass">The PrimalSolver class</a>
        <li><a href="#TheRefinementGlobalandRefinementKellyclasses">The RefinementGlobal and RefinementKelly classes</a>
        <li><a href="#TheRefinementWeightedKellyclass">The RefinementWeightedKelly class</a>
      </ul>
        <li><a href="#Equationdata">Equation data</a>
      <ul>
        <li><a href="#TheSetUpBaseandSetUpclasses">The SetUpBase and SetUp classes</a>
        <li><a href="#TheCurvedRidgesclass">The CurvedRidges class</a>
        <li><a href="#TheExercise_2_3class">The Exercise_2_3 class</a>
        <li><a href="#Discussion">Discussion</a>
      </ul>
        <li><a href="#Dualfunctionals">Dual functionals</a>
      <ul>
        <li><a href="#TheDualFunctionalBaseclass">The DualFunctionalBase class</a>
        <li><a href="#ThedualfunctionalPointValueEvaluationclass">The dual functional PointValueEvaluation class</a>
        <li><a href="#ThedualfunctionalPointXDerivativeEvaluationclass">The dual functional PointXDerivativeEvaluation class</a>
      </ul>
        <li><a href="#ExtendingtheLaplaceSolvernamespace">Extending the LaplaceSolver namespace</a>
      <ul>
        <li><a href="#TheDualSolverclass">The DualSolver class</a>
        <li><a href="#TheWeightedResidualclass">The WeightedResidual class</a>
      </ul>
        <li><a href="#Estimatingerrors">Estimating errors</a>
      <ul>
        <li><a href="#Errorestimationdriverfunctions">Error estimation driver functions</a>
        <li><a href="#Estimatingonasinglecell">Estimating on a single cell</a>
        <li><a href="#Computingcelltermerrorcontributions">Computing cell term error contributions</a>
        <li><a href="#Computingedgetermerrorcontributions1">Computing edge term error contributions &mdash; 1</a>
        <li><a href="#Computingedgetermerrorcontributions2">Computing edge term error contributions &mdash; 2</a>
      </ul>
        <li><a href="#Asimulationframework">A simulation framework</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Pointvalues">Point values</a>
        <li><a href="#Comparingrefinementcriteria">Comparing refinement criteria</a>
        <li><a href="#Evaluationofpointstresses">Evaluation of point stresses</a>
        <li><a href="#step13revisited">step-13 revisited</a>
        <li><a href="#Conclusionsandoutlook">Conclusions and outlook</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


<a name="Themaths"></a><h3>The maths</h3>


The Heidelberg group of Professor Rolf Rannacher, to which the three initial
authors of the deal.II library belonged during their PhD time and partly also
afterwards, has been involved with adaptivity and error estimation for finite
element discretizations since the mid-1990ies. The main achievement is the
development of error estimates for arbitrary functionals of the solution, and
of optimal mesh refinement for its computation.

We will not discuss the derivation of these concepts in too great detail, but
will implement the main ideas in the present example program. For a thorough
introduction into the general idea, we refer to the seminal work of Becker and
Rannacher @cite BR95, @cite BR96r, and the overview article of the same authors in
Acta Numerica @cite BR01; the first introduces the concept of error
estimation and adaptivity for general functional output for the Laplace
equation, while the second gives many examples of applications of these
concepts to a large number of other, more complicated equations. For
applications to individual types of equations, see also the publications by
Becker @cite Bec95, @cite Bec98, Kanschat @cite Kan96, @cite FK97, Suttmeier
@cite Sut96, @cite RS97, @cite RS98c, @cite RS99, Bangerth @cite BR99b,
@cite Ban00w, @cite BR01a, @cite Ban02, and Hartmann @cite Har02, @cite HH01,
@cite HH01b. All of these works, from the original introduction by Becker and
Rannacher to individual contributions to particular equations, have later been
summarized in a book by Bangerth and Rannacher that covers all of these topics,
see @cite BR03.


The basic idea is the following: in applications, one is not usually
interested in the solution per se, but rather in certain aspects of it. For
example, in simulations of flow problems, one may want to know the lift or
drag of a body immersed in the fluid; it is this quantity that we want to know
to best accuracy, and whether the rest of the solution of the describing
equations is well resolved is not of primary interest. Likewise, in elasticity
one might want to know about values of the stress at certain points to guess
whether maximal load values of joints are safe, for example. Or, in radiative
transfer problems, mean flux intensities are of interest.

In all the cases just listed, it is the evaluation of a functional $J(u)$ of
the solution which we are interested in, rather than the values of $u$
everywhere. Since the exact solution $u$ is not available, but only its
numerical approximation $u_h$, it is sensible to ask whether the computed
value $J(u_h)$ is within certain limits of the exact value $J(u)$, i.e. we
want to bound the error with respect to this functional, $J(u)-J(u_h)$.

For simplicity of exposition, we henceforth assume that both the quantity of
interest $J$, as well as the equation are linear, and we will in particular
show the derivation for the Laplace equation with homogeneous Dirichlet
boundary conditions, although the concept is much more general. For this
general case, we refer to the references listed above.  The goal is to obtain
bounds on the error, $J(e)=J(u)-J(u_h)$. For this, let us denote by $z$ the
solution of a dual problem, defined as follows:
@f[
  a(\varphi,z) = J(\varphi) \qquad \forall \varphi,
@f]
where $a(\cdot,\cdot)$ is the bilinear form associated with the differential
equation, and the test functions are chosen from the corresponding solution
space. Then, taking as special test function $\varphi=e$ the error, we have
that
@f[
  J(e) = a(e,z)
@f]
and we can, by Galerkin orthogonality, rewrite this as
@f[
  J(e) = a(e,z-\varphi_h)
@f]
where $\varphi_h$ can be chosen from the discrete test space in
whatever way we find convenient.

Concretely, for Laplace's equation, the error identity reads
@f[
  J(e) = (\nabla e, \nabla(z-\varphi_h)).
@f]
Because we want to use this formula not only to compute error, but
also to refine the mesh, we need to rewrite the expression above as a
sum over cells where each cell's contribution can then be used as an
error indicator for this cell.
Thus, we split the scalar products into terms for each cell, and
integrate by parts on each of them:
@f{eqnarray*}
  J(e)
  &=&
  \sum_K (\nabla (u-u_h), \nabla (z-\varphi_h))_K
  \\
  &=&
  \sum_K (-\Delta (u-u_h), z-\varphi_h)_K
  + (\partial_n (u-u_h), z-z_h)_{\partial K}.
@f}
Next we use that $-\Delta u=f$, and that for solutions of the Laplace
equation, the solution is smooth enough that $\partial_n u$ is
continuous almost everywhere -- so the terms involving $\partial_n u$ on one
cell cancels with that on its neighbor, where the normal vector has the
opposite sign. (The same is not true for $\partial_n u_h$, though.)
At the boundary of the domain, where there is no neighbor cell
with which this term could cancel, the weight $z-\varphi_h$ can be chosen as
zero, and the whole term disappears.

Thus, we have
@f{eqnarray*}
  J(e)
  &=&
  \sum_K (f+\Delta u_h, z-\varphi_h)_K
  - (\partial_n u_h, z-\varphi_h)_{\partial K\backslash \partial\Omega}.
@f}
In a final step, note that when taking the normal derivative of $u_h$, we mean
the value of this quantity as taken from this side of the cell (for the usual
Lagrange elements, derivatives are not continuous across edges). We then
rewrite the above formula by exchanging half of the edge integral of cell $K$
with the neighbor cell $K'$, to obtain
@f{eqnarray*}
  J(e)
  &=&
  \sum_K (f+\Delta u_h, z-\varphi_h)_K
  - \frac 12 (\partial_n u_h|_K + \partial_{n'} u_h|_{K'},
              z-\varphi_h)_{\partial K\backslash \partial\Omega}.
@f}
Using that for the normal vectors on adjacent cells we have $n'=-n$, we define the jump of the
normal derivative by
@f[
  [\partial_n u_h] \dealcoloneq \partial_n u_h|_K + \partial_{n'} u_h|_{K'}
  =
  \partial_n u_h|_K - \partial_n u_h|_{K'},
@f]
and get the final form after setting the discrete function $\varphi_h$, which
is by now still arbitrary, to the point interpolation of the dual solution,
$\varphi_h=I_h z$:
@f{eqnarray*}
  J(e)
  &=&
  \sum_K (f+\Delta u_h, z-I_h z)_K
  - \frac 12 ([\partial_n u_h],
              z-I_h z)_{\partial K\backslash \partial\Omega}.
@f}

With this, we have obtained an exact representation of the error of the finite
element discretization with respect to arbitrary (linear) functionals
$J(\cdot)$. Its structure is a weighted form of a residual estimator, as both
$f+\Delta u_h$ and $[\partial_n u_h]$ are cell and edge residuals that vanish
on the exact solution, and $z-I_h z$ are weights indicating how important the
residuals on a certain cell is for the evaluation of the given functional.
Furthermore, it is a cell-wise quantity, so we can use it as a mesh refinement
criterion. The question, is: how to evaluate it? After all, the evaluation
requires knowledge of the dual solution $z$, which carries the information
about the quantity we want to know to best accuracy.

In some, very special cases, this dual solution is known. For example, if the
functional $J(\cdot)$ is the point evaluation, $J(\varphi)=\varphi(x_0)$, then
the dual solution has to satisfy
@f[
  -\Delta z = \delta(x-x_0),
@f]
with the Dirac delta function on the right hand side, and the dual solution is
the Green's function with respect to the point $x_0$. For simple geometries,
this function is analytically known, and we could insert it into the error
representation formula.

However, we do not want to restrict ourselves to such special cases. Rather,
we will compute the dual solution numerically, and approximate $z$ by some
numerically obtained $\tilde z$. We note that it is not sufficient to compute
this approximation $\tilde z$ using the same method as used for the primal
solution $u_h$, since then $\tilde z-I_h \tilde z=0$, and the overall error
estimate would be zero. Rather, the approximation $\tilde z$ has to be from a
larger space than the primal finite element space. There are various ways to
obtain such an approximation (see the cited literature), and we will choose to
compute it with a higher order finite element space. While this is certainly
not the most efficient way, it is simple since we already have all we need to
do that in place, and it also allows for simple experimenting. For more
efficient methods, again refer to the given literature, in particular
@cite BR95, @cite BR03.

With this, we end the discussion of the mathematical side of this program and
turn to the actual implementation.


@note There are two steps above that do not seem necessary if all you
care about is computing the error: namely, (i) the subtraction of
$\phi_h$ from $z$, and (ii) splitting the integral into a sum of cells
and integrating by parts on each. Indeed, neither of these two steps
change $J(e)$ at all, as we only ever consider identities above until
the substitution of $z$ by $\tilde z$. In other words, if you care
only about <i>estimating the global error</i> $J(e)$, then these steps
are not necessary. On the other hand, if you want to use the error
estimate also as a refinement criterion for each cell of the mesh,
then it is necessary to (i) break the estimate into a sum of cells,
and (ii) massage the formulas in such a way that each cell's
contributions have something to do with the local error. (While the
contortions above do not change the value of the <i>sum</i> $J(e)$,
they change the values we compute for each cell $K$.) To this end, we
want to write everything in the form "residual times dual weight"
where a "residual" is something that goes to zero as the approximation
becomes $u_h$ better and better. For example, the quantity $\partial_n
u_h$ is not a residual, since it simply converges to the (normal
component of) the gradient of the exact solution. On the other hand,
$[\partial_n u_h]$ is a residual because it converges to $[\partial_n
u]=0$. All of the steps we have taken above in developing the final
form of $J(e)$ have indeed had the goal of bringing the final formula
into a form where each term converges to zero as the discrete solution
$u_h$ converges to $u$. This then allows considering each cell's
contribution as an "error indicator" that also converges to zero -- as
it should as the mesh is refined.



<a name="Thesoftware"></a><h3>The software</h3>


The step-14 example program builds heavily on the techniques already used in
the step-13 program. Its implementation of the dual weighted residual error
estimator explained above is done by deriving a second class, properly called
<code>DualSolver</code>, from the <code>Solver</code> base class, and having a class
(<code>WeightedResidual</code>) that joins the two again and controls the solution
of the primal and dual problem, and then uses both to compute the error
indicator for mesh refinement.

The program continues the modular concept of the previous example, by
implementing the dual functional, describing quantity of interest, by an
abstract base class, and providing two different functionals which implement
this interface. Adding a different quantity of interest is thus simple.

One of the more fundamental differences is the handling of data. A common case
is that you develop a program that solves a certain equation, and test it with
different right hand sides, different domains, different coefficients and
boundary values, etc. Usually, these have to match, so that exact solutions
are known, or that their combination makes sense at all.

We demonstrate a way how this can be achieved in a simple, yet very flexible
way. We will put everything that belongs to a certain setup into one class,
and provide a little C++ mortar around it, so that entire setups (domains,
coefficients, right hand sides, etc.) can be exchanged by only changing
something in <em>one</em> place.

Going this way a little further, we have also centralized all the other
parameters that describe how the program is to work in one place, such as the
order of the finite element, the maximal number of degrees of freedom, the
evaluation objects that shall be executed on the computed solutions, and so
on. This allows for simpler configuration of the program, and we will show in
a later program how to use a library class that can handle setting these
parameters by reading an input file. The general aim is to reduce the places
within a program where one may have to look when wanting to change some
parameter, as it has turned out in practice that one forgets where they are as
programs grow. Furthermore, putting all options describing what the program
does in a certain run into a file (that can be stored with the results) helps
repeatability of results more than if the various flags were set somewhere in
the program, where their exact values are forgotten after the next change to
this place.

Unfortunately, the program has become rather long. While this admittedly
reduces its usefulness as an example program, we think that it is a very good
starting point for development of a program for other kinds of problems,
involving different equations than the Laplace equation treated here.
Furthermore, it shows everything that we can show you about our way of a
posteriori error estimation, and its structure should make it simple for you
to adjust this method to other problems, other functionals, other geometries,
coefficients, etc.

The author believes that the present program is his masterpiece among the
example programs, regarding the mathematical complexity, as well as the
simplicity to add extensions. If you use this program as a basis for your own
programs, we would kindly like to ask you to state this fact and the name of
the author of the example program, Wolfgang Bangerth, in publications that
arise from that, of your program consists in a considerable part of the
example program.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * Start out with well known things...
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/thread_management.h>
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_tools.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * #include <algorithm>
 * #include <fstream>
 * #include <iostream>
 * #include <list>
 * #include <memory>
 * #include <numeric>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step14
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Evaluatingthesolution"></a> 
 * <h3>Evaluating the solution</h3>
 * 

 * 
 * As mentioned in the introduction, significant parts of the program have
 * simply been taken over from the step-13 example program. We therefore
 * only comment on those things that are new.
 *   

 * 
 * First, the framework for evaluation of solutions is unchanged, i.e. the
 * base class is the same, and the class to evaluate the solution at a grid
 * point is unchanged:
 * 
 * @code
 *   namespace Evaluation
 *   {
 * @endcode
 * 
 * 
 * <a name="TheEvaluationBaseclass"></a> 
 * <h4>The EvaluationBase class</h4>
 * 
 * @code
 *     template <int dim>
 *     class EvaluationBase
 *     {
 *     public:
 *       virtual ~EvaluationBase() = default;
 * 
 *       void set_refinement_cycle(const unsigned int refinement_cycle);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const = 0;
 * 
 *     protected:
 *       unsigned int refinement_cycle;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step)
 *     {
 *       refinement_cycle = step;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePointValueEvaluationclass"></a> 
 * <h4>The PointValueEvaluation class</h4>
 * 
 * @code
 *     template <int dim>
 *     class PointValueEvaluation : public EvaluationBase<dim>
 *     {
 *     public:
 *       PointValueEvaluation(const Point<dim> &evaluation_point);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const override;
 * 
 *       DeclException1(
 *         ExcEvaluationPointNotFound,
 *         Point<dim>,
 *         << "The evaluation point " << arg1
 *         << " was not found among the vertices of the present grid.");
 * 
 *     private:
 *       const Point<dim> evaluation_point;
 *     };
 * 
 * 
 *     template <int dim>
 *     PointValueEvaluation<dim>::PointValueEvaluation(
 *       const Point<dim> &evaluation_point)
 *       : evaluation_point(evaluation_point)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void PointValueEvaluation<dim>::
 *          operator()(const DoFHandler<dim> &dof_handler,
 *                const Vector<double> & solution) const
 *     {
 *       double point_value = 1e20;
 * 
 *       bool evaluation_point_found = false;
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         if (!evaluation_point_found)
 *           for (const auto vertex : cell->vertex_indices())
 *             if (cell->vertex(vertex).distance(evaluation_point) <
 *                 cell->diameter() * 1e-8)
 *               {
 *                 point_value = solution(cell->vertex_dof_index(vertex, 0));
 * 
 *                 evaluation_point_found = true;
 *                 break;
 *               }
 * 
 *       AssertThrow(evaluation_point_found,
 *                   ExcEvaluationPointNotFound(evaluation_point));
 * 
 *       std::cout << "   Point value=" << point_value << std::endl;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePointXDerivativeEvaluationclass"></a> 
 * <h4>The PointXDerivativeEvaluation class</h4>
 * 

 * 
 * Besides the class implementing the evaluation of the solution at one
 * point, we here provide one which evaluates the gradient at a grid
 * point. Since in general the gradient of a finite element function is
 * not continuous at a vertex, we have to be a little bit more careful
 * here. What we do is to loop over all cells, even if we have found the
 * point already on one cell, and use the mean value of the gradient at
 * the vertex taken from all adjacent cells.
 *     

 * 
 * Given the interface of the <code>PointValueEvaluation</code> class, the
 * declaration of this class provides little surprise, and neither does
 * the constructor:
 * 
 * @code
 *     template <int dim>
 *     class PointXDerivativeEvaluation : public EvaluationBase<dim>
 *     {
 *     public:
 *       PointXDerivativeEvaluation(const Point<dim> &evaluation_point);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const;
 * 
 *       DeclException1(
 *         ExcEvaluationPointNotFound,
 *         Point<dim>,
 *         << "The evaluation point " << arg1
 *         << " was not found among the vertices of the present grid.");
 * 
 *     private:
 *       const Point<dim> evaluation_point;
 *     };
 * 
 * 
 *     template <int dim>
 *     PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation(
 *       const Point<dim> &evaluation_point)
 *       : evaluation_point(evaluation_point)
 *     {}
 * 
 * 
 * @endcode
 * 
 * The more interesting things happen inside the function doing the actual
 * evaluation:
 * 
 * @code
 *     template <int dim>
 *     void PointXDerivativeEvaluation<dim>::
 *          operator()(const DoFHandler<dim> &dof_handler,
 *                const Vector<double> & solution) const
 *     {
 * @endcode
 * 
 * This time initialize the return value with something useful, since we
 * will have to add up a number of contributions and take the mean value
 * afterwards...
 * 
 * @code
 *       double point_derivative = 0;
 * 
 * @endcode
 * 
 * ...then have some objects of which the meaning will become clear
 * below...
 * 
 * @code
 *       QTrapezoid<dim>             vertex_quadrature;
 *       FEValues<dim>               fe_values(dof_handler.get_fe(),
 *                               vertex_quadrature,
 *                               update_gradients | update_quadrature_points);
 *       std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size());
 * 
 * @endcode
 * 
 * ...and next loop over all cells and their vertices, and count how
 * often the vertex has been found:
 * 
 * @code
 *       unsigned int evaluation_point_hits = 0;
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         for (const auto vertex : cell->vertex_indices())
 *           if (cell->vertex(vertex) == evaluation_point)
 *             {
 * @endcode
 * 
 * Things are now no more as simple, since we can't get the
 * gradient of the finite element field as before, where we
 * simply had to pick one degree of freedom at a vertex.
 *               

 * 
 * Rather, we have to evaluate the finite element field on this
 * cell, and at a certain point. As you know, evaluating finite
 * element fields at certain points is done through the
 * <code>FEValues</code> class, so we use that. The question is:
 * the <code>FEValues</code> object needs to be a given a
 * quadrature formula and can then compute the values of finite
 * element quantities at the quadrature points. Here, we don't
 * want to do quadrature, we simply want to specify some points!
 *               

 * 
 * Nevertheless, the same way is chosen: use a special
 * quadrature rule with points at the vertices, since these are
 * what we are interested in. The appropriate rule is the
 * trapezoidal rule, so that is the reason why we used that one
 * above.
 *               

 * 
 * Thus: initialize the <code>FEValues</code> object on this
 * cell,
 * 
 * @code
 *               fe_values.reinit(cell);
 * @endcode
 * 
 * and extract the gradients of the solution vector at the
 * vertices:
 * 
 * @code
 *               fe_values.get_function_gradients(solution, solution_gradients);
 * 
 * @endcode
 * 
 * Now we have the gradients at all vertices, so pick out that
 * one which belongs to the evaluation point (note that the
 * order of vertices is not necessarily the same as that of the
 * quadrature points):
 * 
 * @code
 *               unsigned int q_point = 0;
 *               for (; q_point < solution_gradients.size(); ++q_point)
 *                 if (fe_values.quadrature_point(q_point) == evaluation_point)
 *                   break;
 * 
 * @endcode
 * 
 * Check that the evaluation point was indeed found,
 * 
 * @code
 *               Assert(q_point < solution_gradients.size(), ExcInternalError());
 * @endcode
 * 
 * and if so take the x-derivative of the gradient there as the
 * value which we are interested in, and increase the counter
 * indicating how often we have added to that variable:
 * 
 * @code
 *               point_derivative += solution_gradients[q_point][0];
 *               ++evaluation_point_hits;
 * 
 * @endcode
 * 
 * Finally break out of the innermost loop iterating over the
 * vertices of the present cell, since if we have found the
 * evaluation point at one vertex it cannot be at a following
 * vertex as well:
 * 
 * @code
 *               break;
 *             }
 * 
 * @endcode
 * 
 * Now we have looped over all cells and vertices, so check whether the
 * point was found:
 * 
 * @code
 *       AssertThrow(evaluation_point_hits > 0,
 *                   ExcEvaluationPointNotFound(evaluation_point));
 * 
 * @endcode
 * 
 * We have simply summed up the contributions of all adjacent cells, so
 * we still have to compute the mean value. Once this is done, report
 * the status:
 * 
 * @code
 *       point_derivative /= evaluation_point_hits;
 *       std::cout << "   Point x-derivative=" << point_derivative << std::endl;
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheGridOutputclass"></a> 
 * <h4>The GridOutput class</h4>
 * 

 * 
 * Since this program has a more difficult structure (it computed a dual
 * solution in addition to a primal one), writing out the solution is no
 * more done by an evaluation object since we want to write both solutions
 * at once into one file, and that requires some more information than
 * available to the evaluation classes.
 *     

 * 
 * However, we also want to look at the grids generated. This again can be
 * done with one such class. Its structure is analog to the
 * <code>SolutionOutput</code> class of the previous example program, so
 * we do not discuss it here in more detail. Furthermore, everything that
 * is used here has already been used in previous example programs.
 * 
 * @code
 *     template <int dim>
 *     class GridOutput : public EvaluationBase<dim>
 *     {
 *     public:
 *       GridOutput(const std::string &output_name_base);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const override;
 * 
 *     private:
 *       const std::string output_name_base;
 *     };
 * 
 * 
 *     template <int dim>
 *     GridOutput<dim>::GridOutput(const std::string &output_name_base)
 *       : output_name_base(output_name_base)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void GridOutput<dim>::operator()(const DoFHandler<dim> &dof_handler,
 *                                      const Vector<double> & /*solution*/) const
 *     {
 *       std::ofstream out(output_name_base + "-" +
 *                         std::to_string(this->refinement_cycle) + ".svg");
 *       GridOut().write_svg(dof_handler.get_triangulation(), out);
 *     }
 *   } // namespace Evaluation
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheLaplacesolverclasses"></a> 
 * <h3>The Laplace solver classes</h3>
 * 

 * 
 * Next are the actual solver classes. Again, we discuss only the
 * differences to the previous program.
 * 
 * @code
 *   namespace LaplaceSolver
 *   {
 * @endcode
 * 
 * 
 * <a name="TheLaplacesolverbaseclass"></a> 
 * <h4>The Laplace solver base class</h4>
 * 

 * 
 * This class is almost unchanged, with the exception that it declares two
 * more functions: <code>output_solution</code> will be used to generate
 * output files from the actual solutions computed by derived classes, and
 * the <code>set_refinement_cycle</code> function by which the testing
 * framework sets the number of the refinement cycle to a local variable
 * in this class; this number is later used to generate filenames for the
 * solution output.
 * 
 * @code
 *     template <int dim>
 *     class Base
 *     {
 *     public:
 *       Base(Triangulation<dim> &coarse_grid);
 *       virtual ~Base() = default;
 * 
 *       virtual void solve_problem() = 0;
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
 *       virtual void         refine_grid()                            = 0;
 *       virtual unsigned int n_dofs() const                           = 0;
 * 
 *       virtual void set_refinement_cycle(const unsigned int cycle);
 * 
 *       virtual void output_solution() const = 0;
 * 
 *     protected:
 *       const SmartPointer<Triangulation<dim>> triangulation;
 * 
 *       unsigned int refinement_cycle;
 *     };
 * 
 * 
 *     template <int dim>
 *     Base<dim>::Base(Triangulation<dim> &coarse_grid)
 *       : triangulation(&coarse_grid)
 *       , refinement_cycle(numbers::invalid_unsigned_int)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void Base<dim>::set_refinement_cycle(const unsigned int cycle)
 *     {
 *       refinement_cycle = cycle;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheLaplaceSolverclass"></a> 
 * <h4>The Laplace Solver class</h4>
 * 

 * 
 * Likewise, the <code>Solver</code> class is entirely unchanged and will
 * thus not be discussed.
 * 
 * @code
 *     template <int dim>
 *     class Solver : public virtual Base<dim>
 *     {
 *     public:
 *       Solver(Triangulation<dim> &       triangulation,
 *              const FiniteElement<dim> & fe,
 *              const Quadrature<dim> &    quadrature,
 *              const Quadrature<dim - 1> &face_quadrature,
 *              const Function<dim> &      boundary_values);
 *       virtual ~Solver() override;
 * 
 *       virtual void solve_problem() override;
 * 
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const override;
 * 
 *       virtual unsigned int n_dofs() const override;
 * 
 *     protected:
 *       const SmartPointer<const FiniteElement<dim>>  fe;
 *       const SmartPointer<const Quadrature<dim>>     quadrature;
 *       const SmartPointer<const Quadrature<dim - 1>> face_quadrature;
 *       DoFHandler<dim>                               dof_handler;
 *       Vector<double>                                solution;
 *       const SmartPointer<const Function<dim>>       boundary_values;
 * 
 *       virtual void assemble_rhs(Vector<double> &rhs) const = 0;
 * 
 *     private:
 *       struct LinearSystem
 *       {
 *         LinearSystem(const DoFHandler<dim> &dof_handler);
 * 
 *         void solve(Vector<double> &solution) const;
 * 
 *         AffineConstraints<double> hanging_node_constraints;
 *         SparsityPattern           sparsity_pattern;
 *         SparseMatrix<double>      matrix;
 *         Vector<double>            rhs;
 *       };
 * 
 * 
 * @endcode
 * 
 * The remainder of the class is essentially a copy of step-13
 * as well, including the data structures and functions
 * necessary to compute the linear system in parallel using the
 * WorkStream framework:
 * 
 * @code
 *       struct AssemblyScratchData
 *       {
 *         AssemblyScratchData(const FiniteElement<dim> &fe,
 *                             const Quadrature<dim> &   quadrature);
 *         AssemblyScratchData(const AssemblyScratchData &scratch_data);
 * 
 *         FEValues<dim> fe_values;
 *       };
 * 
 *       struct AssemblyCopyData
 *       {
 *         FullMatrix<double>                   cell_matrix;
 *         std::vector<types::global_dof_index> local_dof_indices;
 *       };
 * 
 * 
 *       void assemble_linear_system(LinearSystem &linear_system);
 * 
 *       void local_assemble_matrix(
 *         const typename DoFHandler<dim>::active_cell_iterator &cell,
 *         AssemblyScratchData &                                 scratch_data,
 *         AssemblyCopyData &                                    copy_data) const;
 * 
 * 
 *       void copy_local_to_global(const AssemblyCopyData &copy_data,
 *                                 LinearSystem &          linear_system) const;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     Solver<dim>::Solver(Triangulation<dim> &       triangulation,
 *                         const FiniteElement<dim> & fe,
 *                         const Quadrature<dim> &    quadrature,
 *                         const Quadrature<dim - 1> &face_quadrature,
 *                         const Function<dim> &      boundary_values)
 *       : Base<dim>(triangulation)
 *       , fe(&fe)
 *       , quadrature(&quadrature)
 *       , face_quadrature(&face_quadrature)
 *       , dof_handler(triangulation)
 *       , boundary_values(&boundary_values)
 *     {}
 * 
 * 
 *     template <int dim>
 *     Solver<dim>::~Solver()
 *     {
 *       dof_handler.clear();
 *     }
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::solve_problem()
 *     {
 *       dof_handler.distribute_dofs(*fe);
 *       solution.reinit(dof_handler.n_dofs());
 * 
 *       LinearSystem linear_system(dof_handler);
 *       assemble_linear_system(linear_system);
 *       linear_system.solve(solution);
 *     }
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::postprocess(
 *       const Evaluation::EvaluationBase<dim> &postprocessor) const
 *     {
 *       postprocessor(dof_handler, solution);
 *     }
 * 
 * 
 *     template <int dim>
 *     unsigned int Solver<dim>::n_dofs() const
 *     {
 *       return dof_handler.n_dofs();
 *     }
 * 
 * 
 * @endcode
 * 
 * The following few functions and constructors are verbatim
 * copies taken from step-13:
 * 
 * @code
 *     template <int dim>
 *     void Solver<dim>::assemble_linear_system(LinearSystem &linear_system)
 *     {
 *       Threads::Task<void> rhs_task =
 *         Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs);
 * 
 *       auto worker =
 *         [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                AssemblyScratchData &scratch_data,
 *                AssemblyCopyData &   copy_data) {
 *           this->local_assemble_matrix(cell, scratch_data, copy_data);
 *         };
 * 
 *       auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) {
 *         this->copy_local_to_global(copy_data, linear_system);
 *       };
 * 
 *       WorkStream::run(dof_handler.begin_active(),
 *                       dof_handler.end(),
 *                       worker,
 *                       copier,
 *                       AssemblyScratchData(*fe, *quadrature),
 *                       AssemblyCopyData());
 *       linear_system.hanging_node_constraints.condense(linear_system.matrix);
 * 
 *       std::map<types::global_dof_index, double> boundary_value_map;
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                *boundary_values,
 *                                                boundary_value_map);
 * 
 *       rhs_task.join();
 *       linear_system.hanging_node_constraints.condense(linear_system.rhs);
 * 
 *       MatrixTools::apply_boundary_values(boundary_value_map,
 *                                          linear_system.matrix,
 *                                          solution,
 *                                          linear_system.rhs);
 *     }
 * 
 * 
 *     template <int dim>
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData(
 *       const FiniteElement<dim> &fe,
 *       const Quadrature<dim> &   quadrature)
 *       : fe_values(fe, quadrature, update_gradients | update_JxW_values)
 *     {}
 * 
 * 
 *     template <int dim>
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData(
 *       const AssemblyScratchData &scratch_data)
 *       : fe_values(scratch_data.fe_values.get_fe(),
 *                   scratch_data.fe_values.get_quadrature(),
 *                   update_gradients | update_JxW_values)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::local_assemble_matrix(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       AssemblyScratchData &                                 scratch_data,
 *       AssemblyCopyData &                                    copy_data) const
 *     {
 *       const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
 *       const unsigned int n_q_points    = quadrature->size();
 * 
 *       copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 * 
 *       copy_data.local_dof_indices.resize(dofs_per_cell);
 * 
 *       scratch_data.fe_values.reinit(cell);
 * 
 *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             copy_data.cell_matrix(i, j) +=
 *               (scratch_data.fe_values.shape_grad(i, q_point) *
 *                scratch_data.fe_values.shape_grad(j, q_point) *
 *                scratch_data.fe_values.JxW(q_point));
 * 
 *       cell->get_dof_indices(copy_data.local_dof_indices);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
 *                                            LinearSystem &linear_system) const
 *     {
 *       for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
 *         for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j)
 *           linear_system.matrix.add(copy_data.local_dof_indices[i],
 *                                    copy_data.local_dof_indices[j],
 *                                    copy_data.cell_matrix(i, j));
 *     }
 * 
 * 
 * @endcode
 * 
 * Now for the functions that implement actions in the linear
 * system class. First, the constructor initializes all data
 * elements to their correct sizes, and sets up a number of
 * additional data structures, such as constraints due to hanging
 * nodes. Since setting up the hanging nodes and finding out about
 * the nonzero elements of the matrix is independent, we do that
 * in parallel (if the library was configured to use concurrency,
 * at least; otherwise, the actions are performed
 * sequentially). Note that we start only one thread, and do the
 * second action in the main thread. Since only one thread is
 * generated, we don't use the <code>Threads::ThreadGroup</code>
 * class here, but rather use the one created thread object
 * directly to wait for this particular thread's exit. The
 * approach is generally the same as the one we have used in
 * <code>Solver::assemble_linear_system()</code> above.
 *     

 * 
 * Note that taking the address of the
 * <code>DoFTools::make_hanging_node_constraints</code> function
 * is a little tricky, since there are actually three functions of
 * this name, one for each supported space dimension. Taking
 * addresses of overloaded functions is somewhat complicated in
 * C++, since the address-of operator <code>&</code> in that case
 * returns a set of values (the addresses of all
 * functions with that name), and selecting the right one is then
 * the next step. If the context dictates which one to take (for
 * example by assigning to a function pointer of known type), then
 * the compiler can do that by itself, but if this set of pointers
 * shall be given as the argument to a function that takes a
 * template, the compiler could choose all without having a
 * preference for one. We therefore have to make it clear to the
 * compiler which one we would like to have; for this, we could
 * use a cast, but for more clarity, we assign it to a temporary
 * <code>mhnc_p</code> (short for <code>pointer to
 * make_hanging_node_constraints</code>) with the right type, and
 * using this pointer instead.
 * 
 * @code
 *     template <int dim>
 *     Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler)
 *     {
 *       hanging_node_constraints.clear();
 * 
 *       void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
 *         &DoFTools::make_hanging_node_constraints;
 * 
 * @endcode
 * 
 * Start a side task then continue on the main thread
 * 
 * @code
 *       Threads::Task<void> side_task =
 *         Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints);
 * 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 * 
 * 
 * @endcode
 * 
 * Wait for the side task to be done before going further
 * 
 * @code
 *       side_task.join();
 * 
 *       hanging_node_constraints.close();
 *       hanging_node_constraints.condense(dsp);
 *       sparsity_pattern.copy_from(dsp);
 * 
 *       matrix.reinit(sparsity_pattern);
 *       rhs.reinit(dof_handler.n_dofs());
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const
 *     {
 *       SolverControl            solver_control(5000, 1e-12);
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *       preconditioner.initialize(matrix, 1.2);
 * 
 *       cg.solve(matrix, solution, rhs, preconditioner);
 * 
 *       hanging_node_constraints.distribute(solution);
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePrimalSolverclass"></a> 
 * <h4>The PrimalSolver class</h4>
 * 

 * 
 * The <code>PrimalSolver</code> class is also mostly unchanged except for
 * implementing the <code>output_solution</code> function. We keep the
 * <code>GlobalRefinement</code> and <code>RefinementKelly</code> classes
 * in this program, and they can then rely on the default implementation
 * of this function which simply outputs the primal solution. The class
 * implementing dual weighted error estimators will overload this function
 * itself, to also output the dual solution.
 * 
 * @code
 *     template <int dim>
 *     class PrimalSolver : public Solver<dim>
 *     {
 *     public:
 *       PrimalSolver(Triangulation<dim> &       triangulation,
 *                    const FiniteElement<dim> & fe,
 *                    const Quadrature<dim> &    quadrature,
 *                    const Quadrature<dim - 1> &face_quadrature,
 *                    const Function<dim> &      rhs_function,
 *                    const Function<dim> &      boundary_values);
 * 
 *       virtual void output_solution() const override;
 * 
 *     protected:
 *       const SmartPointer<const Function<dim>> rhs_function;
 *       virtual void assemble_rhs(Vector<double> &rhs) const override;
 *     };
 * 
 * 
 *     template <int dim>
 *     PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation,
 *                                     const FiniteElement<dim> & fe,
 *                                     const Quadrature<dim> &    quadrature,
 *                                     const Quadrature<dim - 1> &face_quadrature,
 *                                     const Function<dim> &      rhs_function,
 *                                     const Function<dim> &      boundary_values)
 *       : Base<dim>(triangulation)
 *       , Solver<dim>(triangulation,
 *                     fe,
 *                     quadrature,
 *                     face_quadrature,
 *                     boundary_values)
 *       , rhs_function(&rhs_function)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void PrimalSolver<dim>::output_solution() const
 *     {
 *       DataOut<dim> data_out;
 *       data_out.attach_dof_handler(this->dof_handler);
 *       data_out.add_data_vector(this->solution, "solution");
 *       data_out.build_patches();
 * 
 *       std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
 *                         ".vtu");
 *       data_out.write(out, DataOutBase::vtu);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
 *     {
 *       FEValues<dim> fe_values(*this->fe,
 *                               *this->quadrature,
 *                               update_values | update_quadrature_points |
 *                                 update_JxW_values);
 * 
 *       const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
 *       const unsigned int n_q_points    = this->quadrature->size();
 * 
 *       Vector<double>                       cell_rhs(dofs_per_cell);
 *       std::vector<double>                  rhs_values(n_q_points);
 *       std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *       for (const auto &cell : this->dof_handler.active_cell_iterators())
 *         {
 *           cell_rhs = 0;
 * 
 *           fe_values.reinit(cell);
 * 
 *           rhs_function->value_list(fe_values.get_quadrature_points(),
 *                                    rhs_values);
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
 *                               rhs_values[q_point] *               // f((x_q)
 *                               fe_values.JxW(q_point));            // dx
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRefinementGlobalandRefinementKellyclasses"></a> 
 * <h4>The RefinementGlobal and RefinementKelly classes</h4>
 * 

 * 
 * For the following two classes, the same applies as for most of the
 * above: the class is taken from the previous example as-is:
 * 
 * @code
 *     template <int dim>
 *     class RefinementGlobal : public PrimalSolver<dim>
 *     {
 *     public:
 *       RefinementGlobal(Triangulation<dim> &       coarse_grid,
 *                        const FiniteElement<dim> & fe,
 *                        const Quadrature<dim> &    quadrature,
 *                        const Quadrature<dim - 1> &face_quadrature,
 *                        const Function<dim> &      rhs_function,
 *                        const Function<dim> &      boundary_values);
 * 
 *       virtual void refine_grid() override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     RefinementGlobal<dim>::RefinementGlobal(
 *       Triangulation<dim> &       coarse_grid,
 *       const FiniteElement<dim> & fe,
 *       const Quadrature<dim> &    quadrature,
 *       const Quadrature<dim - 1> &face_quadrature,
 *       const Function<dim> &      rhs_function,
 *       const Function<dim> &      boundary_values)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           fe,
 *                           quadrature,
 *                           face_quadrature,
 *                           rhs_function,
 *                           boundary_values)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void RefinementGlobal<dim>::refine_grid()
 *     {
 *       this->triangulation->refine_global(1);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     class RefinementKelly : public PrimalSolver<dim>
 *     {
 *     public:
 *       RefinementKelly(Triangulation<dim> &       coarse_grid,
 *                       const FiniteElement<dim> & fe,
 *                       const Quadrature<dim> &    quadrature,
 *                       const Quadrature<dim - 1> &face_quadrature,
 *                       const Function<dim> &      rhs_function,
 *                       const Function<dim> &      boundary_values);
 * 
 *       virtual void refine_grid() override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     RefinementKelly<dim>::RefinementKelly(
 *       Triangulation<dim> &       coarse_grid,
 *       const FiniteElement<dim> & fe,
 *       const Quadrature<dim> &    quadrature,
 *       const Quadrature<dim - 1> &face_quadrature,
 *       const Function<dim> &      rhs_function,
 *       const Function<dim> &      boundary_values)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           fe,
 *                           quadrature,
 *                           face_quadrature,
 *                           rhs_function,
 *                           boundary_values)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void RefinementKelly<dim>::refine_grid()
 *     {
 *       Vector<float> estimated_error_per_cell(
 *         this->triangulation->n_active_cells());
 *       KellyErrorEstimator<dim>::estimate(
 *         this->dof_handler,
 *         QGauss<dim - 1>(this->fe->degree + 1),
 *         std::map<types::boundary_id, const Function<dim> *>(),
 *         this->solution,
 *         estimated_error_per_cell);
 *       GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.3,
 *                                                       0.03);
 *       this->triangulation->execute_coarsening_and_refinement();
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRefinementWeightedKellyclass"></a> 
 * <h4>The RefinementWeightedKelly class</h4>
 * 

 * 
 * This class is a variant of the previous one, in that it allows to
 * weight the refinement indicators we get from the library's Kelly
 * indicator by some function. We include this class since the goal of
 * this example program is to demonstrate automatic refinement criteria
 * even for complex output quantities such as point values or stresses. If
 * we did not solve a dual problem and compute the weights thereof, we
 * would probably be tempted to give a hand-crafted weighting to the
 * indicators to account for the fact that we are going to evaluate these
 * quantities. This class accepts such a weighting function as argument to
 * its constructor:
 * 
 * @code
 *     template <int dim>
 *     class RefinementWeightedKelly : public PrimalSolver<dim>
 *     {
 *     public:
 *       RefinementWeightedKelly(Triangulation<dim> &       coarse_grid,
 *                               const FiniteElement<dim> & fe,
 *                               const Quadrature<dim> &    quadrature,
 *                               const Quadrature<dim - 1> &face_quadrature,
 *                               const Function<dim> &      rhs_function,
 *                               const Function<dim> &      boundary_values,
 *                               const Function<dim> &      weighting_function);
 * 
 *       virtual void refine_grid() override;
 * 
 *     private:
 *       const SmartPointer<const Function<dim>> weighting_function;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     RefinementWeightedKelly<dim>::RefinementWeightedKelly(
 *       Triangulation<dim> &       coarse_grid,
 *       const FiniteElement<dim> & fe,
 *       const Quadrature<dim> &    quadrature,
 *       const Quadrature<dim - 1> &face_quadrature,
 *       const Function<dim> &      rhs_function,
 *       const Function<dim> &      boundary_values,
 *       const Function<dim> &      weighting_function)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           fe,
 *                           quadrature,
 *                           face_quadrature,
 *                           rhs_function,
 *                           boundary_values)
 *       , weighting_function(&weighting_function)
 *     {}
 * 
 * 
 * 
 * @endcode
 * 
 * Now, here comes the main function, including the weighting:
 * 
 * @code
 *     template <int dim>
 *     void RefinementWeightedKelly<dim>::refine_grid()
 *     {
 * @endcode
 * 
 * First compute some residual based error indicators for all cells by a
 * method already implemented in the library. What exactly we compute
 * here is described in more detail in the documentation of that class.
 * 
 * @code
 *       Vector<float> estimated_error_per_cell(
 *         this->triangulation->n_active_cells());
 *       std::map<types::boundary_id, const Function<dim> *> dummy_function_map;
 *       KellyErrorEstimator<dim>::estimate(this->dof_handler,
 *                                          *this->face_quadrature,
 *                                          dummy_function_map,
 *                                          this->solution,
 *                                          estimated_error_per_cell);
 * 
 * @endcode
 * 
 * Next weigh each entry in the vector of indicators by the value of the
 * function given to the constructor, evaluated at the cell center. We
 * need to write the result into the vector entry that corresponds to the
 * current cell, which we can obtain by asking the cell what its index
 * among all active cells is using CellAccessor::active_cell_index(). (In
 * reality, this index is zero for the first cell we handle in the loop,
 * one for the second cell, etc., and we could as well just keep track of
 * this index using an integer counter; but using
 * CellAccessor::active_cell_index() makes this more explicit.)
 * 
 * @code
 *       for (const auto &cell : this->dof_handler.active_cell_iterators())
 *         estimated_error_per_cell(cell->active_cell_index()) *=
 *           weighting_function->value(cell->center());
 * 
 *       GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.3,
 *                                                       0.03);
 *       this->triangulation->execute_coarsening_and_refinement();
 *     }
 * 
 *   } // namespace LaplaceSolver
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 *   

 * 
 * In this example program, we work with the same data sets as in the
 * previous one, but as it may so happen that someone wants to run the
 * program with different boundary values and right hand side functions, or
 * on a different grid, we show a simple technique to do exactly that. For
 * more clarity, we furthermore pack everything that has to do with equation
 * data into a namespace of its own.
 *   

 * 
 * The underlying assumption is that this is a research program, and that
 * there we often have a number of test cases that consist of a domain, a
 * right hand side, boundary values, possibly a specified coefficient, and a
 * number of other parameters. They often vary all at the same time when
 * shifting from one example to another. To make handling such sets of
 * problem description parameters simple is the goal of the following.
 *   

 * 
 * Basically, the idea is this: let us have a structure for each set of
 * data, in which we pack everything that describes a test case: here, these
 * are two subclasses, one called <code>BoundaryValues</code> for the
 * boundary values of the exact solution, and one called
 * <code>RightHandSide</code>, and then a way to generate the coarse
 * grid. Since the solution of the previous example program looked like
 * curved ridges, we use this name here for the enclosing class. Note that
 * the names of the two inner classes have to be the same for all enclosing
 * test case classes, and also that we have attached the dimension template
 * argument to the enclosing class rather than to the inner ones, to make
 * further processing simpler.  (From a language viewpoint, a namespace
 * would be better to encapsulate these inner classes, rather than a
 * structure. However, namespaces cannot be given as template arguments, so
 * we use a structure to allow a second object to select from within its
 * given argument. The enclosing structure, of course, has no member
 * variables apart from the classes it declares, and a static function to
 * generate the coarse mesh; it will in general never be instantiated.)
 *   

 * 
 * The idea is then the following (this is the right time to also take a
 * brief look at the code below): we can generate objects for boundary
 * values and right hand side by simply giving the name of the outer class
 * as a template argument to a class which we call here
 * <code>Data::SetUp</code>, and it then creates objects for the inner
 * classes. In this case, to get all that characterizes the curved ridge
 * solution, we would simply generate an instance of
 * <code>Data::SetUp@<Data::CurvedRidge@></code>, and everything we need to
 * know about the solution would be static member variables and functions of
 * that object.
 *   

 * 
 * This approach might seem like overkill in this case, but will become very
 * handy once a certain set up is not only characterized by Dirichlet
 * boundary values and a right hand side function, but in addition by
 * material properties, Neumann values, different boundary descriptors,
 * etc. In that case, the <code>SetUp</code> class might consist of a dozen
 * or more objects, and each descriptor class (like the
 * <code>CurvedRidges</code> class below) would have to provide them. Then,
 * you will be happy to be able to change from one set of data to another by
 * only changing the template argument to the <code>SetUp</code> class at
 * one place, rather than at many.
 *   

 * 
 * With this framework for different test cases, we are almost finished, but
 * one thing remains: by now we can select statically, by changing one
 * template argument, which data set to choose. In order to be able to do
 * that dynamically, i.e. at run time, we need a base class. This we provide
 * in the obvious way, see below, with virtual abstract functions. It forces
 * us to introduce a second template parameter <code>dim</code> which we
 * need for the base class (which could be avoided using some template
 * magic, but we omit that), but that's all.
 *   

 * 
 * Adding new testcases is now simple, you don't have to touch the framework
 * classes, only a structure like the <code>CurvedRidges</code> one is
 * needed.
 * 
 * @code
 *   namespace Data
 *   {
 * @endcode
 * 
 * 
 * <a name="TheSetUpBaseandSetUpclasses"></a> 
 * <h4>The SetUpBase and SetUp classes</h4>
 * 

 * 
 * Based on the above description, the <code>SetUpBase</code> class then
 * looks as follows. To allow using the <code>SmartPointer</code> class
 * with this class, we derived from the <code>Subscriptor</code> class.
 * 
 * @code
 *     template <int dim>
 *     struct SetUpBase : public Subscriptor
 *     {
 *       virtual const Function<dim> &get_boundary_values() const = 0;
 * 
 *       virtual const Function<dim> &get_right_hand_side() const = 0;
 * 
 *       virtual void
 *       create_coarse_grid(Triangulation<dim> &coarse_grid) const = 0;
 *     };
 * 
 * 
 * @endcode
 * 
 * And now for the derived class that takes the template argument as
 * explained above.
 *     

 * 
 * Here we pack the data elements into private variables, and allow access
 * to them through the methods of the base class.
 * 
 * @code
 *     template <class Traits, int dim>
 *     struct SetUp : public SetUpBase<dim>
 *     {
 *       virtual const Function<dim> &get_boundary_values() const override;
 * 
 *       virtual const Function<dim> &get_right_hand_side() const override;
 * 
 * 
 *       virtual void
 *       create_coarse_grid(Triangulation<dim> &coarse_grid) const override;
 * 
 *     private:
 *       static const typename Traits::BoundaryValues boundary_values;
 *       static const typename Traits::RightHandSide  right_hand_side;
 *     };
 * 
 * @endcode
 * 
 * We have to provide definitions for the static member variables of the
 * above class:
 * 
 * @code
 *     template <class Traits, int dim>
 *     const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values;
 *     template <class Traits, int dim>
 *     const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side;
 * 
 * @endcode
 * 
 * And definitions of the member functions:
 * 
 * @code
 *     template <class Traits, int dim>
 *     const Function<dim> &SetUp<Traits, dim>::get_boundary_values() const
 *     {
 *       return boundary_values;
 *     }
 * 
 * 
 *     template <class Traits, int dim>
 *     const Function<dim> &SetUp<Traits, dim>::get_right_hand_side() const
 *     {
 *       return right_hand_side;
 *     }
 * 
 * 
 *     template <class Traits, int dim>
 *     void SetUp<Traits, dim>::create_coarse_grid(
 *       Triangulation<dim> &coarse_grid) const
 *     {
 *       Traits::create_coarse_grid(coarse_grid);
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheCurvedRidgesclass"></a> 
 * <h4>The CurvedRidges class</h4>
 * 

 * 
 * The class that is used to describe the boundary values and right hand
 * side of the <code>curved ridge</code> problem already used in the
 * step-13 example program is then like so:
 * 
 * @code
 *     template <int dim>
 *     struct CurvedRidges
 *     {
 *       class BoundaryValues : public Function<dim>
 *       {
 *       public:
 *         virtual double value(const Point<dim> & p,
 *                              const unsigned int component) const;
 *       };
 * 
 * 
 *       class RightHandSide : public Function<dim>
 *       {
 *       public:
 *         virtual double value(const Point<dim> & p,
 *                              const unsigned int component) const;
 *       };
 * 
 *       static void create_coarse_grid(Triangulation<dim> &coarse_grid);
 *     };
 * 
 * 
 *     template <int dim>
 *     double CurvedRidges<dim>::BoundaryValues::value(
 *       const Point<dim> &p,
 *       const unsigned int /*component*/) const
 *     {
 *       double q = p(0);
 *       for (unsigned int i = 1; i < dim; ++i)
 *         q += std::sin(10 * p(i) + 5 * p(0) * p(0));
 *       const double exponential = std::exp(q);
 *       return exponential;
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     double CurvedRidges<dim>::RightHandSide::value(
 *       const Point<dim> &p,
 *       const unsigned int /*component*/) const
 *     {
 *       double q = p(0);
 *       for (unsigned int i = 1; i < dim; ++i)
 *         q += std::sin(10 * p(i) + 5 * p(0) * p(0));
 *       const double u  = std::exp(q);
 *       double       t1 = 1, t2 = 0, t3 = 0;
 *       for (unsigned int i = 1; i < dim; ++i)
 *         {
 *           t1 += std::cos(10 * p(i) + 5 * p(0) * p(0)) * 10 * p(0);
 *           t2 += 10 * std::cos(10 * p(i) + 5 * p(0) * p(0)) -
 *                 100 * std::sin(10 * p(i) + 5 * p(0) * p(0)) * p(0) * p(0);
 *           t3 += 100 * std::cos(10 * p(i) + 5 * p(0) * p(0)) *
 *                   std::cos(10 * p(i) + 5 * p(0) * p(0)) -
 *                 100 * std::sin(10 * p(i) + 5 * p(0) * p(0));
 *         }
 *       t1 = t1 * t1;
 * 
 *       return -u * (t1 + t2 + t3);
 *     }
 * 
 * 
 *     template <int dim>
 *     void CurvedRidges<dim>::create_coarse_grid(Triangulation<dim> &coarse_grid)
 *     {
 *       GridGenerator::hyper_cube(coarse_grid, -1, 1);
 *       coarse_grid.refine_global(2);
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheExercise_2_3class"></a> 
 * <h4>The Exercise_2_3 class</h4>
 * 

 * 
 * This example program was written while giving practical courses for a
 * lecture on adaptive finite element methods and duality based error
 * estimates. For these courses, we had one exercise, which required to
 * solve the Laplace equation with constant right hand side on a square
 * domain with a square hole in the center, and zero boundary
 * values. Since the implementation of the properties of this problem is
 * so particularly simple here, lets do it. As the number of the exercise
 * was 2.3, we take the liberty to retain this name for the class as well.
 * 
 * @code
 *     template <int dim>
 *     struct Exercise_2_3
 *     {
 * @endcode
 * 
 * We need a class to denote the boundary values of the problem. In this
 * case, this is simple: it's the zero function, so don't even declare a
 * class, just an alias:
 * 
 * @code
 *       using BoundaryValues = Functions::ZeroFunction<dim>;
 * 
 * @endcode
 * 
 * Second, a class that denotes the right hand side. Since they are
 * constant, just subclass the corresponding class of the library and be
 * done:
 * 
 * @code
 *       class RightHandSide : public Functions::ConstantFunction<dim>
 *       {
 *       public:
 *         RightHandSide()
 *           : Functions::ConstantFunction<dim>(1.)
 *         {}
 *       };
 * 
 * @endcode
 * 
 * Finally a function to generate the coarse grid. This is somewhat more
 * complicated here, see immediately below.
 * 
 * @code
 *       static void create_coarse_grid(Triangulation<dim> &coarse_grid);
 *     };
 * 
 * 
 * @endcode
 * 
 * As stated above, the grid for this example is the square [-1,1]^2 with
 * the square [-1/2,1/2]^2 as hole in it. We create the coarse grid as 4
 * times 4 cells with the middle four ones missing. To understand how
 * exactly the mesh is going to look, it may be simplest to just look
 * at the "Results" section of this tutorial program first. In general,
 * if you'd like to understand more about creating meshes either from
 * scratch by hand, as we do here, or using other techniques, you
 * should take a look at step-49.
 *     

 * 
 * Of course, the example has an extension to 3d, but since this function
 * cannot be written in a dimension independent way we choose not to
 * implement this here, but rather only specialize the template for
 * dim=2. If you compile the program for 3d, you'll get a message from the
 * linker that this function is not implemented for 3d, and needs to be
 * provided.
 *     

 * 
 * For the creation of this geometry, the library has no predefined
 * method. In this case, the geometry is still simple enough to do the
 * creation by hand, rather than using a mesh generator.
 * 
 * @code
 *     template <>
 *     void Exercise_2_3<2>::create_coarse_grid(Triangulation<2> &coarse_grid)
 *     {
 * @endcode
 * 
 * We first define the space dimension, to allow those parts of the
 * function that are actually dimension independent to use this
 * variable. That makes it simpler if you later take this as a starting
 * point to implement a 3d version of this mesh. The next step is then
 * to have a list of vertices. Here, they are 24 (5 times 5, with the
 * middle one omitted). It is probably best to draw a sketch here.
 * 
 * @code
 *       const unsigned int dim = 2;
 * 
 *       const std::vector<Point<2>> vertices = {
 *         {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, 
 *         {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, 
 *         {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               
 *         {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, 
 *         {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};
 * 
 * @endcode
 * 
 * Next, we have to define the cells and the vertices they contain.
 * 
 * @code
 *       const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
 *         cell_vertices = {{{0, 1, 5, 6}},
 *                          {{1, 2, 6, 7}},
 *                          {{2, 3, 7, 8}},
 *                          {{3, 4, 8, 9}},
 *                          {{5, 6, 10, 11}},
 *                          {{8, 9, 12, 13}},
 *                          {{10, 11, 14, 15}},
 *                          {{12, 13, 17, 18}},
 *                          {{14, 15, 19, 20}},
 *                          {{15, 16, 20, 21}},
 *                          {{16, 17, 21, 22}},
 *                          {{17, 18, 22, 23}}};
 * 
 *       const unsigned int n_cells = cell_vertices.size();
 * 
 * @endcode
 * 
 * Again, we generate a C++ vector type from this, but this time by
 * looping over the cells (yes, this is boring). Additionally, we set
 * the material indicator to zero for all the cells:
 * 
 * @code
 *       std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
 *       for (unsigned int i = 0; i < n_cells; ++i)
 *         {
 *           for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
 *             cells[i].vertices[j] = cell_vertices[i][j];
 *           cells[i].material_id = 0;
 *         }
 * 
 * @endcode
 * 
 * Finally pass all this information to the library to generate a
 * triangulation. The last parameter may be used to pass information
 * about non-zero boundary indicators at certain faces of the
 * triangulation to the library, but we don't want that here, so we give
 * an empty object:
 * 
 * @code
 *       coarse_grid.create_triangulation(vertices, cells, SubCellData());
 * 
 * @endcode
 * 
 * And since we want that the evaluation point (3/4,3/4) in this example
 * is a grid point, we refine once globally:
 * 
 * @code
 *       coarse_grid.refine_global(1);
 *     }
 *   } // namespace Data
 * 
 * @endcode
 * 
 * 
 * <a name="Discussion"></a> 
 * <h4>Discussion</h4>
 *   

 * 
 * As you have now read through this framework, you may be wondering why we
 * have not chosen to implement the classes implementing a certain setup
 * (like the <code>CurvedRidges</code> class) directly as classes derived
 * from <code>Data::SetUpBase</code>. Indeed, we could have done very well
 * so. The only reason is that then we would have to have member variables
 * for the solution and right hand side classes in the
 * <code>CurvedRidges</code> class, as well as member functions overloading
 * the abstract functions of the base class giving access to these member
 * variables. The <code>SetUp</code> class has the sole reason to relieve us
 * from the need to reiterate these member variables and functions that
 * would be necessary in all such classes. In some way, the template
 * mechanism here only provides a way to have default implementations for a
 * number of functions that depend on external quantities and can thus not
 * be provided using normal virtual functions, at least not without the help
 * of templates.
 *   

 * 
 * However, there might be good reasons to actually implement classes
 * derived from <code>Data::SetUpBase</code>, for example if the solution or
 * right hand side classes require constructors that take arguments, which
 * the <code>Data::SetUpBase</code> class cannot provide. In that case,
 * subclassing is a worthwhile strategy. Other possibilities for special
 * cases are to derive from <code>Data::SetUp@<SomeSetUp@></code> where
 * <code>SomeSetUp</code> denotes a class, or even to explicitly specialize
 * <code>Data::SetUp@<SomeSetUp@></code>. The latter allows to transparently
 * use the way the <code>SetUp</code> class is used for other set-ups, but
 * with special actions taken for special arguments.
 *   

 * 
 * A final observation favoring the approach taken here is the following: we
 * have found numerous times that when starting a project, the number of
 * parameters (usually boundary values, right hand side, coarse grid, just
 * as here) was small, and the number of test cases was small as well. One
 * then starts out by handcoding them into a number of <code>switch</code>
 * statements. Over time, projects grow, and so does the number of test
 * cases. The number of <code>switch</code> statements grows with that, and
 * their length as well, and one starts to find ways to consider impossible
 * examples where domains, boundary values, and right hand sides do not fit
 * together any more, and starts losing the overview over the whole
 * structure. Encapsulating everything belonging to a certain test case into
 * a structure of its own has proven worthwhile for this, as it keeps
 * everything that belongs to one test case in one place. Furthermore, it
 * allows to put these things all in one or more files that are only devoted
 * to test cases and their data, without having to bring their actual
 * implementation into contact with the rest of the program.
 * 

 * 
 * 

 * 
 * 
 * <a name="Dualfunctionals"></a> 
 * <h3>Dual functionals</h3>
 * 

 * 
 * As with the other components of the program, we put everything we need to
 * describe dual functionals into a namespace of its own, and define an
 * abstract base class that provides the interface the class solving the
 * dual problem needs for its work.
 *   

 * 
 * We will then implement two such classes, for the evaluation of a point
 * value and of the derivative of the solution at that point. For these
 * functionals we already have the corresponding evaluation objects, so they
 * are complementary.
 * 
 * @code
 *   namespace DualFunctional
 *   {
 * @endcode
 * 
 * 
 * <a name="TheDualFunctionalBaseclass"></a> 
 * <h4>The DualFunctionalBase class</h4>
 * 

 * 
 * First start with the base class for dual functionals. Since for linear
 * problems the characteristics of the dual problem play a role only in
 * the right hand side, we only need to provide for a function that
 * assembles the right hand side for a given discretization:
 * 
 * @code
 *     template <int dim>
 *     class DualFunctionalBase : public Subscriptor
 *     {
 *     public:
 *       virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
 *                                 Vector<double> &       rhs) const = 0;
 *     };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThedualfunctionalPointValueEvaluationclass"></a> 
 * <h4>The dual functional PointValueEvaluation class</h4>
 * 

 * 
 * As a first application, we consider the functional corresponding to the
 * evaluation of the solution's value at a given point which again we
 * assume to be a vertex. Apart from the constructor that takes and stores
 * the evaluation point, this class consists only of the function that
 * implements assembling the right hand side.
 * 
 * @code
 *     template <int dim>
 *     class PointValueEvaluation : public DualFunctionalBase<dim>
 *     {
 *     public:
 *       PointValueEvaluation(const Point<dim> &evaluation_point);
 * 
 *       virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
 *                                 Vector<double> &       rhs) const override;
 * 
 *       DeclException1(
 *         ExcEvaluationPointNotFound,
 *         Point<dim>,
 *         << "The evaluation point " << arg1
 *         << " was not found among the vertices of the present grid.");
 * 
 *     protected:
 *       const Point<dim> evaluation_point;
 *     };
 * 
 * 
 *     template <int dim>
 *     PointValueEvaluation<dim>::PointValueEvaluation(
 *       const Point<dim> &evaluation_point)
 *       : evaluation_point(evaluation_point)
 *     {}
 * 
 * 
 * @endcode
 * 
 * As for doing the main purpose of the class, assembling the right hand
 * side, let us first consider what is necessary: The right hand side of
 * the dual problem is a vector of values J(phi_i), where J is the error
 * functional, and phi_i is the i-th shape function. Here, J is the
 * evaluation at the point x0, i.e. J(phi_i)=phi_i(x0).
 *     

 * 
 * Now, we have assumed that the evaluation point is a vertex. Thus, for
 * the usual finite elements we might be using in this program, we can
 * take for granted that at such a point exactly one shape function is
 * nonzero, and in particular has the value one. Thus, we set the right
 * hand side vector to all-zeros, then seek for the shape function
 * associated with that point and set the corresponding value of the right
 * hand side vector to one:
 * 
 * @code
 *     template <int dim>
 *     void
 *     PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler,
 *                                             Vector<double> &       rhs) const
 *     {
 * @endcode
 * 
 * So, first set everything to zeros...
 * 
 * @code
 *       rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * ...then loop over cells and find the evaluation point among the
 * vertices (or very close to a vertex, which may happen due to floating
 * point round-off):
 * 
 * @code
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         for (const auto vertex : cell->vertex_indices())
 *           if (cell->vertex(vertex).distance(evaluation_point) <
 *               cell->diameter() * 1e-8)
 *             {
 * @endcode
 * 
 * Ok, found, so set corresponding entry, and leave function
 * since we are finished:
 * 
 * @code
 *               rhs(cell->vertex_dof_index(vertex, 0)) = 1;
 *               return;
 *             }
 * 
 * @endcode
 * 
 * Finally, a sanity check: if we somehow got here, then we must have
 * missed the evaluation point, so raise an exception unconditionally:
 * 
 * @code
 *       AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point));
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThedualfunctionalPointXDerivativeEvaluationclass"></a> 
 * <h4>The dual functional PointXDerivativeEvaluation class</h4>
 * 

 * 
 * As second application, we again consider the evaluation of the
 * x-derivative of the solution at one point. Again, the declaration of
 * the class, and the implementation of its constructor is not too
 * interesting:
 * 
 * @code
 *     template <int dim>
 *     class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
 *     {
 *     public:
 *       PointXDerivativeEvaluation(const Point<dim> &evaluation_point);
 * 
 *       virtual void assemble_rhs(const DoFHandler<dim> &dof_handler,
 *                                 Vector<double> &       rhs) const;
 * 
 *       DeclException1(
 *         ExcEvaluationPointNotFound,
 *         Point<dim>,
 *         << "The evaluation point " << arg1
 *         << " was not found among the vertices of the present grid.");
 * 
 *     protected:
 *       const Point<dim> evaluation_point;
 *     };
 * 
 * 
 *     template <int dim>
 *     PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation(
 *       const Point<dim> &evaluation_point)
 *       : evaluation_point(evaluation_point)
 *     {}
 * 
 * 
 * @endcode
 * 
 * What is interesting is the implementation of this functional: here,
 * J(phi_i)=d/dx phi_i(x0).
 *     

 * 
 * We could, as in the implementation of the respective evaluation object
 * take the average of the gradients of each shape function phi_i at this
 * evaluation point. However, we take a slightly different approach: we
 * simply take the average over all cells that surround this point. The
 * question which cells <code>surrounds</code> the evaluation point is
 * made dependent on the mesh width by including those cells for which the
 * distance of the cell's midpoint to the evaluation point is less than
 * the cell's diameter.
 *     

 * 
 * Taking the average of the gradient over the area/volume of these cells
 * leads to a dual solution which is very close to the one which would
 * result from the point evaluation of the gradient. It is simple to
 * justify theoretically that this does not change the method
 * significantly.
 * 
 * @code
 *     template <int dim>
 *     void PointXDerivativeEvaluation<dim>::assemble_rhs(
 *       const DoFHandler<dim> &dof_handler,
 *       Vector<double> &       rhs) const
 *     {
 * @endcode
 * 
 * Again, first set all entries to zero:
 * 
 * @code
 *       rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * Initialize a <code>FEValues</code> object with a quadrature formula,
 * have abbreviations for the number of quadrature points and shape
 * functions...
 * 
 * @code
 *       QGauss<dim>        quadrature(dof_handler.get_fe().degree + 1);
 *       FEValues<dim>      fe_values(dof_handler.get_fe(),
 *                               quadrature,
 *                               update_gradients | update_quadrature_points |
 *                                 update_JxW_values);
 *       const unsigned int n_q_points    = fe_values.n_quadrature_points;
 *       const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
 * 
 * @endcode
 * 
 * ...and have two objects that are used to store the global indices of
 * the degrees of freedom on a cell, and the values of the gradients of
 * the shape functions at the quadrature points:
 * 
 * @code
 *       Vector<double>            cell_rhs(dofs_per_cell);
 *       std::vector<unsigned int> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Finally have a variable in which we will sum up the area/volume of
 * the cells over which we integrate, by integrating the unit functions
 * on these cells:
 * 
 * @code
 *       double total_volume = 0;
 * 
 * @endcode
 * 
 * Then start the loop over all cells, and select those cells which are
 * close enough to the evaluation point:
 * 
 * @code
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         if (cell->center().distance(evaluation_point) <= cell->diameter())
 *           {
 * @endcode
 * 
 * If we have found such a cell, then initialize the
 * <code>FEValues</code> object and integrate the x-component of
 * the gradient of each shape function, as well as the unit
 * function for the total area/volume.
 * 
 * @code
 *             fe_values.reinit(cell);
 *             cell_rhs = 0;
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q)
 *               {
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   cell_rhs(i) +=
 *                     fe_values.shape_grad(i, q)[0] // (d/dx phi_i(x_q))
 *                     * fe_values.JxW(q);           // * dx
 *                 total_volume += fe_values.JxW(q);
 *               }
 * 
 * @endcode
 * 
 * If we have the local contributions, distribute them to the
 * global vector:
 * 
 * @code
 *             cell->get_dof_indices(local_dof_indices);
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               rhs(local_dof_indices[i]) += cell_rhs(i);
 *           }
 * 
 * @endcode
 * 
 * After we have looped over all cells, check whether we have found any
 * at all, by making sure that their volume is non-zero. If not, then
 * the results will be botched, as the right hand side should then still
 * be zero, so throw an exception:
 * 
 * @code
 *       AssertThrow(total_volume > 0,
 *                   ExcEvaluationPointNotFound(evaluation_point));
 * 
 * @endcode
 * 
 * Finally, we have by now only integrated the gradients of the shape
 * functions, not taking their mean value. We fix this by dividing by
 * the measure of the volume over which we have integrated:
 * 
 * @code
 *       rhs /= total_volume;
 *     }
 * 
 * 
 *   } // namespace DualFunctional
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ExtendingtheLaplaceSolvernamespace"></a> 
 * <h3>Extending the LaplaceSolver namespace</h3>
 * 
 * @code
 *   namespace LaplaceSolver
 *   {
 * @endcode
 * 
 * 
 * <a name="TheDualSolverclass"></a> 
 * <h4>The DualSolver class</h4>
 * 

 * 
 * In the same way as the <code>PrimalSolver</code> class above, we now
 * implement a <code>DualSolver</code>. It has all the same features, the
 * only difference is that it does not take a function object denoting a
 * right hand side object, but now takes a <code>DualFunctionalBase</code>
 * object that will assemble the right hand side vector of the dual
 * problem. The rest of the class is rather trivial.
 *     

 * 
 * Since both primal and dual solver will use the same triangulation, but
 * different discretizations, it now becomes clear why we have made the
 * <code>Base</code> class a virtual one: since the final class will be
 * derived from both <code>PrimalSolver</code> as well as
 * <code>DualSolver</code>, it would have two <code>Base</code> instances,
 * would we not have marked the inheritance as virtual. Since in many
 * applications the base class would store much more information than just
 * the triangulation which needs to be shared between primal and dual
 * solvers, we do not usually want to use two such base classes.
 * 
 * @code
 *     template <int dim>
 *     class DualSolver : public Solver<dim>
 *     {
 *     public:
 *       DualSolver(
 *         Triangulation<dim> &                           triangulation,
 *         const FiniteElement<dim> &                     fe,
 *         const Quadrature<dim> &                        quadrature,
 *         const Quadrature<dim - 1> &                    face_quadrature,
 *         const DualFunctional::DualFunctionalBase<dim> &dual_functional);
 * 
 *     protected:
 *       const SmartPointer<const DualFunctional::DualFunctionalBase<dim>>
 *                    dual_functional;
 *       virtual void assemble_rhs(Vector<double> &rhs) const override;
 * 
 *       static const Functions::ZeroFunction<dim> boundary_values;
 *     };
 * 
 *     template <int dim>
 *     const Functions::ZeroFunction<dim> DualSolver<dim>::boundary_values;
 * 
 *     template <int dim>
 *     DualSolver<dim>::DualSolver(
 *       Triangulation<dim> &                           triangulation,
 *       const FiniteElement<dim> &                     fe,
 *       const Quadrature<dim> &                        quadrature,
 *       const Quadrature<dim - 1> &                    face_quadrature,
 *       const DualFunctional::DualFunctionalBase<dim> &dual_functional)
 *       : Base<dim>(triangulation)
 *       , Solver<dim>(triangulation,
 *                     fe,
 *                     quadrature,
 *                     face_quadrature,
 *                     boundary_values)
 *       , dual_functional(&dual_functional)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const
 *     {
 *       dual_functional->assemble_rhs(this->dof_handler, rhs);
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheWeightedResidualclass"></a> 
 * <h4>The WeightedResidual class</h4>
 * 

 * 
 * Here finally comes the main class of this program, the one that
 * implements the dual weighted residual error estimator. It joins the
 * primal and dual solver classes to use them for the computation of
 * primal and dual solutions, and implements the error representation
 * formula for use as error estimate and mesh refinement.
 *     

 * 
 * The first few of the functions of this class are mostly overriders of
 * the respective functions of the base class:
 * 
 * @code
 *     template <int dim>
 *     class WeightedResidual : public PrimalSolver<dim>, public DualSolver<dim>
 *     {
 *     public:
 *       WeightedResidual(
 *         Triangulation<dim> &                           coarse_grid,
 *         const FiniteElement<dim> &                     primal_fe,
 *         const FiniteElement<dim> &                     dual_fe,
 *         const Quadrature<dim> &                        quadrature,
 *         const Quadrature<dim - 1> &                    face_quadrature,
 *         const Function<dim> &                          rhs_function,
 *         const Function<dim> &                          boundary_values,
 *         const DualFunctional::DualFunctionalBase<dim> &dual_functional);
 * 
 *       virtual void solve_problem() override;
 * 
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const override;
 * 
 *       virtual unsigned int n_dofs() const override;
 * 
 *       virtual void refine_grid() override;
 * 
 *       virtual void output_solution() const override;
 * 
 *     private:
 * @endcode
 * 
 * In the private section, we have two functions that are used to call
 * the <code>solve_problem</code> functions of the primal and dual base
 * classes. These two functions will be called in parallel by the
 * <code>solve_problem</code> function of this class.
 * 
 * @code
 *       void solve_primal_problem();
 *       void solve_dual_problem();
 * @endcode
 * 
 * Then declare abbreviations for active cell iterators, to avoid that
 * we have to write this lengthy name over and over again:
 * 

 * 
 * 
 * @code
 *       using active_cell_iterator =
 *         typename DoFHandler<dim>::active_cell_iterator;
 * 
 * @endcode
 * 
 * Next, declare a data type that we will us to store the contribution
 * of faces to the error estimator. The idea is that we can compute the
 * face terms from each of the two cells to this face, as they are the
 * same when viewed from both sides. What we will do is to compute them
 * only once, based on some rules explained below which of the two
 * adjacent cells will be in charge to do so. We then store the
 * contribution of each face in a map mapping faces to their values, and
 * only collect the contributions for each cell by looping over the
 * cells a second time and grabbing the values from the map.
 *       

 * 
 * The data type of this map is declared here:
 * 
 * @code
 *       using FaceIntegrals =
 *         typename std::map<typename DoFHandler<dim>::face_iterator, double>;
 * 
 * @endcode
 * 
 * In the computation of the error estimates on cells and faces, we need
 * a number of helper objects, such as <code>FEValues</code> and
 * <code>FEFaceValues</code> functions, but also temporary objects
 * storing the values and gradients of primal and dual solutions, for
 * example. These fields are needed in the three functions that do the
 * integration on cells, and regular and irregular faces, respectively.
 *       

 * 
 * There are three reasonable ways to provide these fields: first, as
 * local variables in the function that needs them; second, as member
 * variables of this class; third, as arguments passed to that function.
 *       

 * 
 * These three alternatives all have drawbacks: the third that their
 * number is not negligible and would make calling these functions a
 * lengthy enterprise. The second has the drawback that it disallows
 * parallelization, since the threads that will compute the error
 * estimate have to have their own copies of these variables each, so
 * member variables of the enclosing class will not work. The first
 * approach, although straightforward, has a subtle but important
 * drawback: we will call these functions over and over again, many
 * thousands of times maybe; it now turns out that allocating
 * vectors and other objects that need memory from the heap is an
 * expensive business in terms of run-time, since memory allocation is
 * expensive when several threads are involved. It is thus
 * significantly better to allocate the memory only once, and recycle
 * the objects as often as possible.
 *       

 * 
 * What to do? Our answer is to use a variant of the third strategy.
 * In fact, this is exactly what the WorkStream concept is supposed to
 * do (we have already introduced it above, but see also @ref threads).
 * To avoid that we have to give these functions a dozen or so
 * arguments, we pack all these variables into two structures, one which
 * is used for the computations on cells, the other doing them on the
 * faces. Both are then joined into the WeightedResidualScratchData class
 * that will serve as the "scratch data" class of the WorkStream concept:
 * 
 * @code
 *       struct CellData
 *       {
 *         FEValues<dim>                           fe_values;
 *         const SmartPointer<const Function<dim>> right_hand_side;
 * 
 *         std::vector<double> cell_residual;
 *         std::vector<double> rhs_values;
 *         std::vector<double> dual_weights;
 *         std::vector<double> cell_laplacians;
 *         CellData(const FiniteElement<dim> &fe,
 *                  const Quadrature<dim> &   quadrature,
 *                  const Function<dim> &     right_hand_side);
 *         CellData(const CellData &cell_data);
 *       };
 * 
 *       struct FaceData
 *       {
 *         FEFaceValues<dim>    fe_face_values_cell;
 *         FEFaceValues<dim>    fe_face_values_neighbor;
 *         FESubfaceValues<dim> fe_subface_values_cell;
 * 
 *         std::vector<double>                  jump_residual;
 *         std::vector<double>                  dual_weights;
 *         typename std::vector<Tensor<1, dim>> cell_grads;
 *         typename std::vector<Tensor<1, dim>> neighbor_grads;
 *         FaceData(const FiniteElement<dim> & fe,
 *                  const Quadrature<dim - 1> &face_quadrature);
 *         FaceData(const FaceData &face_data);
 *       };
 * 
 *       struct WeightedResidualScratchData
 *       {
 *         WeightedResidualScratchData(
 *           const FiniteElement<dim> & primal_fe,
 *           const Quadrature<dim> &    primal_quadrature,
 *           const Quadrature<dim - 1> &primal_face_quadrature,
 *           const Function<dim> &      rhs_function,
 *           const Vector<double> &     primal_solution,
 *           const Vector<double> &     dual_weights);
 * 
 *         WeightedResidualScratchData(
 *           const WeightedResidualScratchData &scratch_data);
 * 
 *         CellData       cell_data;
 *         FaceData       face_data;
 *         Vector<double> primal_solution;
 *         Vector<double> dual_weights;
 *       };
 * 
 * 
 * @endcode
 * 
 * WorkStream::run generally wants both a scratch object and a copy
 * object. Here, for reasons similar to what we had in step-9 when
 * discussing the computation of an approximation of the gradient, we
 * don't actually need a "copy data" structure. Since WorkStream insists
 * on having one of these, we just declare an empty structure that does
 * nothing other than being there.
 * 
 * @code
 *       struct WeightedResidualCopyData
 *       {};
 * 
 * 
 * 
 * @endcode
 * 
 * Regarding the evaluation of the error estimator, we have one driver
 * function that uses WorkStream::run() to call the second function on
 * every cell:
 * 
 * @code
 *       void estimate_error(Vector<float> &error_indicators) const;
 * 
 *       void estimate_on_one_cell(const active_cell_iterator & cell,
 *                                 WeightedResidualScratchData &scratch_data,
 *                                 WeightedResidualCopyData &   copy_data,
 *                                 Vector<float> &              error_indicators,
 *                                 FaceIntegrals &face_integrals) const;
 * 
 * @endcode
 * 
 * Then we have functions that do the actual integration of the error
 * representation formula. They will treat the terms on the cell
 * interiors, on those faces that have no hanging nodes, and on those
 * faces with hanging nodes, respectively:
 * 
 * @code
 *       void integrate_over_cell(const active_cell_iterator &cell,
 *                                const Vector<double> &      primal_solution,
 *                                const Vector<double> &      dual_weights,
 *                                CellData &                  cell_data,
 *                                Vector<float> &error_indicators) const;
 * 
 *       void integrate_over_regular_face(const active_cell_iterator &cell,
 *                                        const unsigned int          face_no,
 *                                        const Vector<double> &primal_solution,
 *                                        const Vector<double> &dual_weights,
 *                                        FaceData &            face_data,
 *                                        FaceIntegrals &face_integrals) const;
 *       void integrate_over_irregular_face(const active_cell_iterator &cell,
 *                                          const unsigned int          face_no,
 *                                          const Vector<double> &primal_solution,
 *                                          const Vector<double> &dual_weights,
 *                                          FaceData &            face_data,
 *                                          FaceIntegrals &face_integrals) const;
 *     };
 * 
 * 
 * 
 * @endcode
 * 
 * In the implementation of this class, we first have the constructors of
 * the <code>CellData</code> and <code>FaceData</code> member classes, and
 * the <code>WeightedResidual</code> constructor. They only initialize
 * fields to their correct lengths, so we do not have to discuss them in
 * too much detail:
 * 
 * @code
 *     template <int dim>
 *     WeightedResidual<dim>::CellData::CellData(
 *       const FiniteElement<dim> &fe,
 *       const Quadrature<dim> &   quadrature,
 *       const Function<dim> &     right_hand_side)
 *       : fe_values(fe,
 *                   quadrature,
 *                   update_values | update_hessians | update_quadrature_points |
 *                     update_JxW_values)
 *       , right_hand_side(&right_hand_side)
 *       , cell_residual(quadrature.size())
 *       , rhs_values(quadrature.size())
 *       , dual_weights(quadrature.size())
 *       , cell_laplacians(quadrature.size())
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::CellData::CellData(const CellData &cell_data)
 *       : fe_values(cell_data.fe_values.get_fe(),
 *                   cell_data.fe_values.get_quadrature(),
 *                   update_values | update_hessians | update_quadrature_points |
 *                     update_JxW_values)
 *       , right_hand_side(cell_data.right_hand_side)
 *       , cell_residual(cell_data.cell_residual)
 *       , rhs_values(cell_data.rhs_values)
 *       , dual_weights(cell_data.dual_weights)
 *       , cell_laplacians(cell_data.cell_laplacians)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::FaceData::FaceData(
 *       const FiniteElement<dim> & fe,
 *       const Quadrature<dim - 1> &face_quadrature)
 *       : fe_face_values_cell(fe,
 *                             face_quadrature,
 *                             update_values | update_gradients |
 *                               update_JxW_values | update_normal_vectors)
 *       , fe_face_values_neighbor(fe,
 *                                 face_quadrature,
 *                                 update_values | update_gradients |
 *                                   update_JxW_values | update_normal_vectors)
 *       , fe_subface_values_cell(fe, face_quadrature, update_gradients)
 *     {
 *       const unsigned int n_face_q_points = face_quadrature.size();
 * 
 *       jump_residual.resize(n_face_q_points);
 *       dual_weights.resize(n_face_q_points);
 *       cell_grads.resize(n_face_q_points);
 *       neighbor_grads.resize(n_face_q_points);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::FaceData::FaceData(const FaceData &face_data)
 *       : fe_face_values_cell(face_data.fe_face_values_cell.get_fe(),
 *                             face_data.fe_face_values_cell.get_quadrature(),
 *                             update_values | update_gradients |
 *                               update_JxW_values | update_normal_vectors)
 *       , fe_face_values_neighbor(
 *           face_data.fe_face_values_neighbor.get_fe(),
 *           face_data.fe_face_values_neighbor.get_quadrature(),
 *           update_values | update_gradients | update_JxW_values |
 *             update_normal_vectors)
 *       , fe_subface_values_cell(
 *           face_data.fe_subface_values_cell.get_fe(),
 *           face_data.fe_subface_values_cell.get_quadrature(),
 *           update_gradients)
 *       , jump_residual(face_data.jump_residual)
 *       , dual_weights(face_data.dual_weights)
 *       , cell_grads(face_data.cell_grads)
 *       , neighbor_grads(face_data.neighbor_grads)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::WeightedResidualScratchData::
 *       WeightedResidualScratchData(
 *         const FiniteElement<dim> & primal_fe,
 *         const Quadrature<dim> &    primal_quadrature,
 *         const Quadrature<dim - 1> &primal_face_quadrature,
 *         const Function<dim> &      rhs_function,
 *         const Vector<double> &     primal_solution,
 *         const Vector<double> &     dual_weights)
 *       : cell_data(primal_fe, primal_quadrature, rhs_function)
 *       , face_data(primal_fe, primal_face_quadrature)
 *       , primal_solution(primal_solution)
 *       , dual_weights(dual_weights)
 *     {}
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::WeightedResidualScratchData::
 *       WeightedResidualScratchData(
 *         const WeightedResidualScratchData &scratch_data)
 *       : cell_data(scratch_data.cell_data)
 *       , face_data(scratch_data.face_data)
 *       , primal_solution(scratch_data.primal_solution)
 *       , dual_weights(scratch_data.dual_weights)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     WeightedResidual<dim>::WeightedResidual(
 *       Triangulation<dim> &                           coarse_grid,
 *       const FiniteElement<dim> &                     primal_fe,
 *       const FiniteElement<dim> &                     dual_fe,
 *       const Quadrature<dim> &                        quadrature,
 *       const Quadrature<dim - 1> &                    face_quadrature,
 *       const Function<dim> &                          rhs_function,
 *       const Function<dim> &                          bv,
 *       const DualFunctional::DualFunctionalBase<dim> &dual_functional)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           primal_fe,
 *                           quadrature,
 *                           face_quadrature,
 *                           rhs_function,
 *                           bv)
 *       , DualSolver<dim>(coarse_grid,
 *                         dual_fe,
 *                         quadrature,
 *                         face_quadrature,
 *                         dual_functional)
 *     {}
 * 
 * 
 * @endcode
 * 
 * The next five functions are boring, as they simply relay their work to
 * the base classes. The first calls the primal and dual solvers in
 * parallel, while postprocessing the solution and retrieving the number
 * of degrees of freedom is done by the primal class.
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::solve_problem()
 *     {
 *       Threads::TaskGroup<void> tasks;
 *       tasks +=
 *         Threads::new_task(&WeightedResidual<dim>::solve_primal_problem, *this);
 *       tasks +=
 *         Threads::new_task(&WeightedResidual<dim>::solve_dual_problem, *this);
 *       tasks.join_all();
 *     }
 * 
 * 
 *     template <int dim>
 *     void WeightedResidual<dim>::solve_primal_problem()
 *     {
 *       PrimalSolver<dim>::solve_problem();
 *     }
 * 
 *     template <int dim>
 *     void WeightedResidual<dim>::solve_dual_problem()
 *     {
 *       DualSolver<dim>::solve_problem();
 *     }
 * 
 * 
 *     template <int dim>
 *     void WeightedResidual<dim>::postprocess(
 *       const Evaluation::EvaluationBase<dim> &postprocessor) const
 *     {
 *       PrimalSolver<dim>::postprocess(postprocessor);
 *     }
 * 
 * 
 *     template <int dim>
 *     unsigned int WeightedResidual<dim>::n_dofs() const
 *     {
 *       return PrimalSolver<dim>::n_dofs();
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * Now, it is becoming more interesting: the <code>refine_grid()</code>
 * function asks the error estimator to compute the cell-wise error
 * indicators, then uses their absolute values for mesh refinement.
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::refine_grid()
 *     {
 * @endcode
 * 
 * First call the function that computes the cell-wise and global error:
 * 
 * @code
 *       Vector<float> error_indicators(this->triangulation->n_active_cells());
 *       estimate_error(error_indicators);
 * 
 * @endcode
 * 
 * Then note that marking cells for refinement or coarsening only works
 * if all indicators are positive, to allow their comparison. Thus, drop
 * the signs on all these indicators:
 * 
 * @code
 *       for (float &error_indicator : error_indicators)
 *         error_indicator = std::fabs(error_indicator);
 * 
 * @endcode
 * 
 * Finally, we can select between different strategies for
 * refinement. The default here is to refine those cells with the
 * largest error indicators that make up for a total of 80 per cent of
 * the error, while we coarsen those with the smallest indicators that
 * make up for the bottom 2 per cent of the error.
 * 
 * @code
 *       GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation,
 *                                                         error_indicators,
 *                                                         0.8,
 *                                                         0.02);
 *       this->triangulation->execute_coarsening_and_refinement();
 *     }
 * 
 * 
 * @endcode
 * 
 * Since we want to output both the primal and the dual solution, we
 * overload the <code>output_solution</code> function. The only
 * interesting feature of this function is that the primal and dual
 * solutions are defined on different finite element spaces, which is not
 * the format the <code>DataOut</code> class expects. Thus, we have to
 * transfer them to a common finite element space. Since we want the
 * solutions only to see them qualitatively, we contend ourselves with
 * interpolating the dual solution to the (smaller) primal space. For the
 * interpolation, there is a library function, that takes a
 * AffineConstraints object including the hanging node
 * constraints. The rest is standard.
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::output_solution() const
 *     {
 *       AffineConstraints<double> primal_hanging_node_constraints;
 *       DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
 *                                               primal_hanging_node_constraints);
 *       primal_hanging_node_constraints.close();
 *       Vector<double> dual_solution(PrimalSolver<dim>::dof_handler.n_dofs());
 *       FETools::interpolate(DualSolver<dim>::dof_handler,
 *                            DualSolver<dim>::solution,
 *                            PrimalSolver<dim>::dof_handler,
 *                            primal_hanging_node_constraints,
 *                            dual_solution);
 * 
 *       DataOut<dim> data_out;
 *       data_out.attach_dof_handler(PrimalSolver<dim>::dof_handler);
 * 
 * @endcode
 * 
 * Add the data vectors for which we want output. Add them both, the
 * <code>DataOut</code> functions can handle as many data vectors as you
 * wish to write to output:
 * 
 * @code
 *       data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution");
 *       data_out.add_data_vector(dual_solution, "dual_solution");
 * 
 *       data_out.build_patches();
 * 
 *       std::ofstream out("solution-" + std::to_string(this->refinement_cycle) +
 *                         ".vtu");
 *       data_out.write(out, DataOutBase::vtu);
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Estimatingerrors"></a> 
 * <h3>Estimating errors</h3>
 * 

 * 
 * 
 * <a name="Errorestimationdriverfunctions"></a> 
 * <h4>Error estimation driver functions</h4>
 *     

 * 
 * As for the actual computation of error estimates, let's start with the
 * function that drives all this, i.e. calls those functions that actually
 * do the work, and finally collects the results.
 * 
 * @code
 *     template <int dim>
 *     void
 *     WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const
 *     {
 * @endcode
 * 
 * The first task in computing the error is to set up vectors that
 * denote the primal solution, and the weights (z-z_h)=(z-I_hz), both in
 * the finite element space for which we have computed the dual
 * solution. For this, we have to interpolate the primal solution to the
 * dual finite element space, and to subtract the interpolation of the
 * computed dual solution to the primal finite element
 * space. Fortunately, the library provides functions for the
 * interpolation into larger or smaller finite element spaces, so this
 * is mostly obvious.
 *       

 * 
 * First, let's do that for the primal solution: it is cell-wise
 * interpolated into the finite element space in which we have solved
 * the dual problem: But, again as in the
 * <code>WeightedResidual::output_solution</code> function we first need
 * to create an AffineConstraints object including the hanging node
 * constraints, but this time of the dual finite element space.
 * 
 * @code
 *       AffineConstraints<double> dual_hanging_node_constraints;
 *       DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler,
 *                                               dual_hanging_node_constraints);
 *       dual_hanging_node_constraints.close();
 *       Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs());
 *       FETools::interpolate(PrimalSolver<dim>::dof_handler,
 *                            PrimalSolver<dim>::solution,
 *                            DualSolver<dim>::dof_handler,
 *                            dual_hanging_node_constraints,
 *                            primal_solution);
 * 
 * @endcode
 * 
 * Then for computing the interpolation of the numerically approximated
 * dual solution z into the finite element space of the primal solution
 * and subtracting it from z: use the
 * <code>interpolate_difference</code> function, that gives (z-I_hz) in
 * the element space of the dual solution.
 * 
 * @code
 *       AffineConstraints<double> primal_hanging_node_constraints;
 *       DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler,
 *                                               primal_hanging_node_constraints);
 *       primal_hanging_node_constraints.close();
 *       Vector<double> dual_weights(DualSolver<dim>::dof_handler.n_dofs());
 *       FETools::interpolation_difference(DualSolver<dim>::dof_handler,
 *                                         dual_hanging_node_constraints,
 *                                         DualSolver<dim>::solution,
 *                                         PrimalSolver<dim>::dof_handler,
 *                                         primal_hanging_node_constraints,
 *                                         dual_weights);
 * 
 * @endcode
 * 
 * Note that this could probably have been more efficient since those
 * constraints have been used previously when assembling matrix and
 * right hand side for the primal problem and writing out the dual
 * solution. We leave the optimization of the program in this respect as
 * an exercise.
 * 

 * 
 * Having computed the dual weights we now proceed with computing the
 * cell and face residuals of the primal solution. First we set up a map
 * between face iterators and their jump term contributions of faces to
 * the error estimator. The reason is that we compute the jump terms
 * only once, from one side of the face, and want to collect them only
 * afterwards when looping over all cells a second time.
 *       

 * 
 * We initialize this map already with a value of -1e20 for all faces,
 * since this value will stand out in the results if something should go
 * wrong and we fail to compute the value for a face for some
 * reason. Secondly, this initialization already makes the std::map
 * object allocate all objects it may possibly need. This is important
 * since we will write into this structure from parallel threads,
 * and doing so would not be thread-safe if the map needed to allocate
 * memory and thereby reshape its data structures. In other words, the
 * initial initialization relieves us from the necessity to synchronize
 * the threads through a mutex each time they write to (and modify the
 * structure of) this map.
 * 
 * @code
 *       FaceIntegrals face_integrals;
 *       for (const auto &cell :
 *            DualSolver<dim>::dof_handler.active_cell_iterators())
 *         for (const auto &face : cell->face_iterators())
 *           face_integrals[face] = -1e20;
 * 
 *       auto worker = [this,
 *                      &error_indicators,
 *                      &face_integrals](const active_cell_iterator & cell,
 *                                       WeightedResidualScratchData &scratch_data,
 *                                       WeightedResidualCopyData &   copy_data) {
 *         this->estimate_on_one_cell(
 *           cell, scratch_data, copy_data, error_indicators, face_integrals);
 *       };
 * 
 *       auto do_nothing_copier =
 *         std::function<void(const WeightedResidualCopyData &)>();
 * 
 * @endcode
 * 
 * Then hand it all off to WorkStream::run() to compute the
 * estimators for all cells in parallel:
 * 
 * @code
 *       WorkStream::run(
 *         DualSolver<dim>::dof_handler.begin_active(),
 *         DualSolver<dim>::dof_handler.end(),
 *         worker,
 *         do_nothing_copier,
 *         WeightedResidualScratchData(*DualSolver<dim>::fe,
 *                                     *DualSolver<dim>::quadrature,
 *                                     *DualSolver<dim>::face_quadrature,
 *                                     *this->rhs_function,
 *                                     primal_solution,
 *                                     dual_weights),
 *         WeightedResidualCopyData());
 * 
 * @endcode
 * 
 * Once the error contributions are computed, sum them up. For this,
 * note that the cell terms are already set, and that only the edge
 * terms need to be collected. Thus, loop over all cells and their
 * faces, make sure that the contributions of each of the faces are
 * there, and add them up. Only take minus one half of the jump term,
 * since the other half will be taken by the neighboring cell.
 * 
 * @code
 *       unsigned int present_cell = 0;
 *       for (const auto &cell :
 *            DualSolver<dim>::dof_handler.active_cell_iterators())
 *         {
 *           for (const auto &face : cell->face_iterators())
 *             {
 *               Assert(face_integrals.find(face) != face_integrals.end(),
 *                      ExcInternalError());
 *               error_indicators(present_cell) -= 0.5 * face_integrals[face];
 *             }
 *           ++present_cell;
 *         }
 *       std::cout << "   Estimated error="
 *                 << std::accumulate(error_indicators.begin(),
 *                                    error_indicators.end(),
 *                                    0.)
 *                 << std::endl;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Estimatingonasinglecell"></a> 
 * <h4>Estimating on a single cell</h4>
 * 

 * 
 * Next we have the function that is called to estimate the error on a
 * single cell. The function may be called multiple times if the library was
 * configured to use multithreading. Here it goes:
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::estimate_on_one_cell(
 *       const active_cell_iterator & cell,
 *       WeightedResidualScratchData &scratch_data,
 *       WeightedResidualCopyData &   copy_data,
 *       Vector<float> &              error_indicators,
 *       FaceIntegrals &              face_integrals) const
 *     {
 * @endcode
 * 
 * Because of WorkStream, estimate_on_one_cell requires a CopyData object
 * even if it is no used. The next line silences a warning about this
 * unused variable.
 * 
 * @code
 *       (void)copy_data;
 * 
 * @endcode
 * 
 * First task on each cell is to compute the cell residual
 * contributions of this cell, and put them into the
 * <code>error_indicators</code> variable:
 * 
 * @code
 *       integrate_over_cell(cell,
 *                           scratch_data.primal_solution,
 *                           scratch_data.dual_weights,
 *                           scratch_data.cell_data,
 *                           error_indicators);
 * 
 * @endcode
 * 
 * After computing the cell terms, turn to the face terms. For this,
 * loop over all faces of the present cell, and see whether
 * something needs to be computed on it:
 * 
 * @code
 *       for (const auto face_no : cell->face_indices())
 *         {
 * @endcode
 * 
 * First, if this face is part of the boundary, then there is
 * nothing to do. However, to make things easier when summing up
 * the contributions of the faces of cells, we enter this face
 * into the list of faces with a zero contribution to the error.
 * 
 * @code
 *           if (cell->face(face_no)->at_boundary())
 *             {
 *               face_integrals[cell->face(face_no)] = 0;
 *               continue;
 *             }
 * 
 * @endcode
 * 
 * Next, note that since we want to compute the jump terms on
 * each face only once although we access it twice (if it is not
 * at the boundary), we have to define some rules who is
 * responsible for computing on a face:
 *           

 * 
 * First, if the neighboring cell is on the same level as this
 * one, i.e. neither further refined not coarser, then the one
 * with the lower index within this level does the work. In
 * other words: if the other one has a lower index, then skip
 * work on this face:
 * 
 * @code
 *           if ((cell->neighbor(face_no)->has_children() == false) &&
 *               (cell->neighbor(face_no)->level() == cell->level()) &&
 *               (cell->neighbor(face_no)->index() < cell->index()))
 *             continue;
 * 
 * @endcode
 * 
 * Likewise, we always work from the coarser cell if this and
 * its neighbor differ in refinement. Thus, if the neighboring
 * cell is less refined than the present one, then do nothing
 * since we integrate over the subfaces when we visit the coarse
 * cell.
 * 
 * @code
 *           if (cell->at_boundary(face_no) == false)
 *             if (cell->neighbor(face_no)->level() < cell->level())
 *               continue;
 * 
 * 
 * @endcode
 * 
 * Now we know that we are in charge here, so actually compute
 * the face jump terms. If the face is a regular one, i.e.  the
 * other side's cell is neither coarser not finer than this
 * cell, then call one function, and if the cell on the other
 * side is further refined, then use another function. Note that
 * the case that the cell on the other side is coarser cannot
 * happen since we have decided above that we handle this case
 * when we pass over that other cell.
 * 
 * @code
 *           if (cell->face(face_no)->has_children() == false)
 *             integrate_over_regular_face(cell,
 *                                         face_no,
 *                                         scratch_data.primal_solution,
 *                                         scratch_data.dual_weights,
 *                                         scratch_data.face_data,
 *                                         face_integrals);
 *           else
 *             integrate_over_irregular_face(cell,
 *                                           face_no,
 *                                           scratch_data.primal_solution,
 *                                           scratch_data.dual_weights,
 *                                           scratch_data.face_data,
 *                                           face_integrals);
 *         }
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingcelltermerrorcontributions"></a> 
 * <h4>Computing cell term error contributions</h4>
 * 

 * 
 * As for the actual computation of the error contributions, first turn to
 * the cell terms:
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::integrate_over_cell(
 *       const active_cell_iterator &cell,
 *       const Vector<double> &      primal_solution,
 *       const Vector<double> &      dual_weights,
 *       CellData &                  cell_data,
 *       Vector<float> &             error_indicators) const
 *     {
 * @endcode
 * 
 * The tasks to be done are what appears natural from looking at the
 * error estimation formula: first get the right hand side and Laplacian
 * of the numerical solution at the quadrature points for the cell
 * residual,
 * 
 * @code
 *       cell_data.fe_values.reinit(cell);
 *       cell_data.right_hand_side->value_list(
 *         cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values);
 *       cell_data.fe_values.get_function_laplacians(primal_solution,
 *                                                   cell_data.cell_laplacians);
 * 
 * @endcode
 * 
 * ...then get the dual weights...
 * 
 * @code
 *       cell_data.fe_values.get_function_values(dual_weights,
 *                                               cell_data.dual_weights);
 * 
 * @endcode
 * 
 * ...and finally build the sum over all quadrature points and store it
 * with the present cell:
 * 
 * @code
 *       double sum = 0;
 *       for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p)
 *         sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) *
 *                 cell_data.dual_weights[p] * cell_data.fe_values.JxW(p));
 *       error_indicators(cell->active_cell_index()) += sum;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingedgetermerrorcontributions1"></a> 
 * <h4>Computing edge term error contributions -- 1</h4>
 * 

 * 
 * On the other hand, computation of the edge terms for the error estimate
 * is not so simple. First, we have to distinguish between faces with and
 * without hanging nodes. Because it is the simple case, we first consider
 * the case without hanging nodes on a face (let's call this the `regular'
 * case):
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::integrate_over_regular_face(
 *       const active_cell_iterator &cell,
 *       const unsigned int          face_no,
 *       const Vector<double> &      primal_solution,
 *       const Vector<double> &      dual_weights,
 *       FaceData &                  face_data,
 *       FaceIntegrals &             face_integrals) const
 *     {
 *       const unsigned int n_q_points =
 *         face_data.fe_face_values_cell.n_quadrature_points;
 * 
 * @endcode
 * 
 * The first step is to get the values of the gradients at the
 * quadrature points of the finite element field on the present
 * cell. For this, initialize the <code>FEFaceValues</code> object
 * corresponding to this side of the face, and extract the gradients
 * using that object.
 * 
 * @code
 *       face_data.fe_face_values_cell.reinit(cell, face_no);
 *       face_data.fe_face_values_cell.get_function_gradients(
 *         primal_solution, face_data.cell_grads);
 * 
 * @endcode
 * 
 * The second step is then to extract the gradients of the finite
 * element solution at the quadrature points on the other side of the
 * face, i.e. from the neighboring cell.
 *       

 * 
 * For this, do a sanity check before: make sure that the neighbor
 * actually exists (yes, we should not have come here if the neighbor
 * did not exist, but in complicated software there are bugs, so better
 * check this), and if this is not the case throw an error.
 * 
 * @code
 *       Assert(cell->neighbor(face_no).state() == IteratorState::valid,
 *              ExcInternalError());
 * @endcode
 * 
 * If we have that, then we need to find out with which face of the
 * neighboring cell we have to work, i.e. the <code>how-many'th</code> the
 * neighbor the present cell is of the cell behind the present face. For
 * this, there is a function, and we put the result into a variable with
 * the name <code>neighbor_neighbor</code>:
 * 
 * @code
 *       const unsigned int neighbor_neighbor =
 *         cell->neighbor_of_neighbor(face_no);
 * @endcode
 * 
 * Then define an abbreviation for the neighbor cell, initialize the
 * <code>FEFaceValues</code> object on that cell, and extract the
 * gradients on that cell:
 * 
 * @code
 *       const active_cell_iterator neighbor = cell->neighbor(face_no);
 *       face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor);
 *       face_data.fe_face_values_neighbor.get_function_gradients(
 *         primal_solution, face_data.neighbor_grads);
 * 
 * @endcode
 * 
 * Now that we have the gradients on this and the neighboring cell,
 * compute the jump residual by multiplying the jump in the gradient
 * with the normal vector:
 * 
 * @code
 *       for (unsigned int p = 0; p < n_q_points; ++p)
 *         face_data.jump_residual[p] =
 *           ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
 *            face_data.fe_face_values_cell.normal_vector(p));
 * 
 * @endcode
 * 
 * Next get the dual weights for this face:
 * 
 * @code
 *       face_data.fe_face_values_cell.get_function_values(dual_weights,
 *                                                         face_data.dual_weights);
 * 
 * @endcode
 * 
 * Finally, we have to compute the sum over jump residuals, dual
 * weights, and quadrature weights, to get the result for this face:
 * 
 * @code
 *       double face_integral = 0;
 *       for (unsigned int p = 0; p < n_q_points; ++p)
 *         face_integral +=
 *           (face_data.jump_residual[p] * face_data.dual_weights[p] *
 *            face_data.fe_face_values_cell.JxW(p));
 * 
 * @endcode
 * 
 * Double check that the element already exists and that it was not
 * already written to...
 * 
 * @code
 *       Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
 *              ExcInternalError());
 *       Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError());
 * 
 * @endcode
 * 
 * ...then store computed value at assigned location. Note that the
 * stored value does not contain the factor 1/2 that appears in the
 * error representation. The reason is that the term actually does not
 * have this factor if we loop over all faces in the triangulation, but
 * only appears if we write it as a sum over all cells and all faces of
 * each cell; we thus visit the same face twice. We take account of this
 * by using this factor -1/2 later, when we sum up the contributions for
 * each cell individually.
 * 
 * @code
 *       face_integrals[cell->face(face_no)] = face_integral;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingedgetermerrorcontributions2"></a> 
 * <h4>Computing edge term error contributions -- 2</h4>
 * 

 * 
 * We are still missing the case of faces with hanging nodes. This is what
 * is covered in this function:
 * 
 * @code
 *     template <int dim>
 *     void WeightedResidual<dim>::integrate_over_irregular_face(
 *       const active_cell_iterator &cell,
 *       const unsigned int          face_no,
 *       const Vector<double> &      primal_solution,
 *       const Vector<double> &      dual_weights,
 *       FaceData &                  face_data,
 *       FaceIntegrals &             face_integrals) const
 *     {
 * @endcode
 * 
 * First again two abbreviations, and some consistency checks whether
 * the function is called only on faces for which it is supposed to be
 * called:
 * 
 * @code
 *       const unsigned int n_q_points =
 *         face_data.fe_face_values_cell.n_quadrature_points;
 * 
 *       const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
 *       const typename DoFHandler<dim>::cell_iterator neighbor =
 *         cell->neighbor(face_no);
 *       Assert(neighbor.state() == IteratorState::valid, ExcInternalError());
 *       Assert(neighbor->has_children(), ExcInternalError());
 *       (void)neighbor;
 * 
 * @endcode
 * 
 * Then find out which neighbor the present cell is of the adjacent
 * cell. Note that we will operate on the children of this adjacent
 * cell, but that their orientation is the same as that of their mother,
 * i.e. the neighbor direction is the same.
 * 
 * @code
 *       const unsigned int neighbor_neighbor =
 *         cell->neighbor_of_neighbor(face_no);
 * 
 * @endcode
 * 
 * Then simply do everything we did in the previous function for one
 * face for all the sub-faces now:
 * 
 * @code
 *       for (unsigned int subface_no = 0; subface_no < face->n_children();
 *            ++subface_no)
 *         {
 * @endcode
 * 
 * Start with some checks again: get an iterator pointing to the
 * cell behind the present subface and check whether its face is a
 * subface of the one we are considering. If that were not the case,
 * then there would be either a bug in the
 * <code>neighbor_neighbor</code> function called above, or -- worse
 * -- some function in the library did not keep to some underlying
 * assumptions about cells, their children, and their faces. In any
 * case, even though this assertion should not be triggered, it does
 * not harm to be cautious, and in optimized mode computations the
 * assertion will be removed anyway.
 * 
 * @code
 *           const active_cell_iterator neighbor_child =
 *             cell->neighbor_child_on_subface(face_no, subface_no);
 *           Assert(neighbor_child->face(neighbor_neighbor) ==
 *                    cell->face(face_no)->child(subface_no),
 *                  ExcInternalError());
 * 
 * @endcode
 * 
 * Now start the work by again getting the gradient of the solution
 * first at this side of the interface,
 * 
 * @code
 *           face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no);
 *           face_data.fe_subface_values_cell.get_function_gradients(
 *             primal_solution, face_data.cell_grads);
 * @endcode
 * 
 * then at the other side,
 * 
 * @code
 *           face_data.fe_face_values_neighbor.reinit(neighbor_child,
 *                                                    neighbor_neighbor);
 *           face_data.fe_face_values_neighbor.get_function_gradients(
 *             primal_solution, face_data.neighbor_grads);
 * 
 * @endcode
 * 
 * and finally building the jump residuals. Since we take the normal
 * vector from the other cell this time, revert the sign of the
 * first term compared to the other function:
 * 
 * @code
 *           for (unsigned int p = 0; p < n_q_points; ++p)
 *             face_data.jump_residual[p] =
 *               ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
 *                face_data.fe_face_values_neighbor.normal_vector(p));
 * 
 * @endcode
 * 
 * Then get dual weights:
 * 
 * @code
 *           face_data.fe_face_values_neighbor.get_function_values(
 *             dual_weights, face_data.dual_weights);
 * 
 * @endcode
 * 
 * At last, sum up the contribution of this sub-face, and set it in
 * the global map:
 * 
 * @code
 *           double face_integral = 0;
 *           for (unsigned int p = 0; p < n_q_points; ++p)
 *             face_integral +=
 *               (face_data.jump_residual[p] * face_data.dual_weights[p] *
 *                face_data.fe_face_values_neighbor.JxW(p));
 *           face_integrals[neighbor_child->face(neighbor_neighbor)] =
 *             face_integral;
 *         }
 * 
 * @endcode
 * 
 * Once the contributions of all sub-faces are computed, loop over all
 * sub-faces to collect and store them with the mother face for simple
 * use when later collecting the error terms of cells. Again make safety
 * checks that the entries for the sub-faces have been computed and do
 * not carry an invalid value.
 * 
 * @code
 *       double sum = 0;
 *       for (unsigned int subface_no = 0; subface_no < face->n_children();
 *            ++subface_no)
 *         {
 *           Assert(face_integrals.find(face->child(subface_no)) !=
 *                    face_integrals.end(),
 *                  ExcInternalError());
 *           Assert(face_integrals[face->child(subface_no)] != -1e20,
 *                  ExcInternalError());
 * 
 *           sum += face_integrals[face->child(subface_no)];
 *         }
 * @endcode
 * 
 * Finally store the value with the parent face.
 * 
 * @code
 *       face_integrals[face] = sum;
 *     }
 * 
 *   } // namespace LaplaceSolver
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Asimulationframework"></a> 
 * <h3>A simulation framework</h3>
 * 

 * 
 * In the previous example program, we have had two functions that were used
 * to drive the process of solving on subsequently finer grids. We extend
 * this here to allow for a number of parameters to be passed to these
 * functions, and put all of that into framework class.
 *   

 * 
 * You will have noted that this program is built up of a number of small
 * parts (evaluation functions, solver classes implementing various
 * refinement methods, different dual functionals, different problem and
 * data descriptions), which makes the program relatively simple to extend,
 * but also allows to solve a large number of different problems by
 * replacing one part by another. We reflect this flexibility by declaring a
 * structure in the following framework class that holds a number of
 * parameters that may be set to test various combinations of the parts of
 * this program, and which can be used to test it at various problems and
 * discretizations in a simple way.
 * 
 * @code
 *   template <int dim>
 *   struct Framework
 *   {
 *   public:
 * @endcode
 * 
 * First, we declare two abbreviations for simple use of the respective
 * data types:
 * 
 * @code
 *     using Evaluator     = Evaluation::EvaluationBase<dim>;
 *     using EvaluatorList = std::list<Evaluator *>;
 * 
 * 
 * @endcode
 * 
 * Then we have the structure which declares all the parameters that may
 * be set. In the default constructor of the structure, these values are
 * all set to default values, for simple use.
 * 
 * @code
 *     struct ProblemDescription
 *     {
 * @endcode
 * 
 * First allow for the degrees of the piecewise polynomials by which the
 * primal and dual problems will be discretized. They default to (bi-,
 * tri-)linear ansatz functions for the primal, and (bi-, tri-)quadratic
 * ones for the dual problem. If a refinement criterion is chosen that
 * does not need the solution of a dual problem, the value of the dual
 * finite element degree is of course ignored.
 * 
 * @code
 *       unsigned int primal_fe_degree;
 *       unsigned int dual_fe_degree;
 * 
 * @endcode
 * 
 * Then have an object that describes the problem type, i.e. right hand
 * side, domain, boundary values, etc. The pointer needed here defaults
 * to the Null pointer, i.e. you will have to set it in actual instances
 * of this object to make it useful.
 * 
 * @code
 *       std::unique_ptr<const Data::SetUpBase<dim>> data;
 * 
 * @endcode
 * 
 * Since we allow to use different refinement criteria (global
 * refinement, refinement by the Kelly error indicator, possibly with a
 * weight, and using the dual estimator), define a number of enumeration
 * values, and subsequently a variable of that type. It will default to
 * <code>dual_weighted_error_estimator</code>.
 * 
 * @code
 *       enum RefinementCriterion
 *       {
 *         dual_weighted_error_estimator,
 *         global_refinement,
 *         kelly_indicator,
 *         weighted_kelly_indicator
 *       };
 * 
 *       RefinementCriterion refinement_criterion;
 * 
 * @endcode
 * 
 * Next, an object that describes the dual functional. It is only needed
 * if the dual weighted residual refinement is chosen, and also defaults
 * to a Null pointer.
 * 
 * @code
 *       std::unique_ptr<const DualFunctional::DualFunctionalBase<dim>>
 *         dual_functional;
 * 
 * @endcode
 * 
 * Then a list of evaluation objects. Its default value is empty,
 * i.e. no evaluation objects.
 * 
 * @code
 *       EvaluatorList evaluator_list;
 * 
 * @endcode
 * 
 * Next to last, a function that is used as a weight to the
 * <code>RefinementWeightedKelly</code> class. The default value of this
 * pointer is zero, but you have to set it to some other value if you
 * want to use the <code>weighted_kelly_indicator</code> refinement
 * criterion.
 * 
 * @code
 *       std::unique_ptr<const Function<dim>> kelly_weight;
 * 
 * @endcode
 * 
 * Finally, we have a variable that denotes the maximum number of
 * degrees of freedom we allow for the (primal) discretization. If it is
 * exceeded, we stop the process of solving and intermittent mesh
 * refinement. Its default value is 20,000.
 * 
 * @code
 *       unsigned int max_degrees_of_freedom;
 * 
 * @endcode
 * 
 * Finally the default constructor of this class:
 * 
 * @code
 *       ProblemDescription();
 *     };
 * 
 * @endcode
 * 
 * The driver framework class only has one method which calls solver and
 * mesh refinement intermittently, and does some other small tasks in
 * between. Since it does not need data besides the parameters given to
 * it, we make it static:
 * 
 * @code
 *     static void run(const ProblemDescription &descriptor);
 *   };
 * 
 * 
 * @endcode
 * 
 * As for the implementation, first the constructor of the parameter object,
 * setting all values to their defaults:
 * 
 * @code
 *   template <int dim>
 *   Framework<dim>::ProblemDescription::ProblemDescription()
 *     : primal_fe_degree(1)
 *     , dual_fe_degree(2)
 *     , refinement_criterion(dual_weighted_error_estimator)
 *     , max_degrees_of_freedom(20000)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * Then the function which drives the whole process:
 * 
 * @code
 *   template <int dim>
 *   void Framework<dim>::run(const ProblemDescription &descriptor)
 *   {
 * @endcode
 * 
 * First create a triangulation from the given data object,
 * 
 * @code
 *     Triangulation<dim> triangulation(
 *       Triangulation<dim>::smoothing_on_refinement);
 *     descriptor.data->create_coarse_grid(triangulation);
 * 
 * @endcode
 * 
 * then a set of finite elements and appropriate quadrature formula:
 * 
 * @code
 *     const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree);
 *     const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree);
 *     const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1);
 *     const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1);
 * 
 * @endcode
 * 
 * Next, select one of the classes implementing different refinement
 * criteria.
 * 
 * @code
 *     std::unique_ptr<LaplaceSolver::Base<dim>> solver;
 *     switch (descriptor.refinement_criterion)
 *       {
 *         case ProblemDescription::dual_weighted_error_estimator:
 *           {
 *             solver = std::make_unique<LaplaceSolver::WeightedResidual<dim>>(
 *               triangulation,
 *               primal_fe,
 *               dual_fe,
 *               quadrature,
 *               face_quadrature,
 *               descriptor.data->get_right_hand_side(),
 *               descriptor.data->get_boundary_values(),
 *               *descriptor.dual_functional);
 *             break;
 *           }
 * 
 *         case ProblemDescription::global_refinement:
 *           {
 *             solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>(
 *               triangulation,
 *               primal_fe,
 *               quadrature,
 *               face_quadrature,
 *               descriptor.data->get_right_hand_side(),
 *               descriptor.data->get_boundary_values());
 *             break;
 *           }
 * 
 *         case ProblemDescription::kelly_indicator:
 *           {
 *             solver = std::make_unique<LaplaceSolver::RefinementKelly<dim>>(
 *               triangulation,
 *               primal_fe,
 *               quadrature,
 *               face_quadrature,
 *               descriptor.data->get_right_hand_side(),
 *               descriptor.data->get_boundary_values());
 *             break;
 *           }
 * 
 *         case ProblemDescription::weighted_kelly_indicator:
 *           {
 *             solver =
 *               std::make_unique<LaplaceSolver::RefinementWeightedKelly<dim>>(
 *                 triangulation,
 *                 primal_fe,
 *                 quadrature,
 *                 face_quadrature,
 *                 descriptor.data->get_right_hand_side(),
 *                 descriptor.data->get_boundary_values(),
 *                 *descriptor.kelly_weight);
 *             break;
 *           }
 * 
 *         default:
 *           AssertThrow(false, ExcInternalError());
 *       }
 * 
 * @endcode
 * 
 * Now that all objects are in place, run the main loop. The stopping
 * criterion is implemented at the bottom of the loop.
 *     

 * 
 * In the loop, first set the new cycle number, then solve the problem,
 * output its solution(s), apply the evaluation objects to it, then decide
 * whether we want to refine the mesh further and solve again on this
 * mesh, or jump out of the loop.
 * 
 * @code
 *     for (unsigned int step = 0; true; ++step)
 *       {
 *         std::cout << "Refinement cycle: " << step << std::endl;
 * 
 *         solver->set_refinement_cycle(step);
 *         solver->solve_problem();
 *         solver->output_solution();
 * 
 *         std::cout << "   Number of degrees of freedom=" << solver->n_dofs()
 *                   << std::endl;
 * 
 *         for (const auto &evaluator : descriptor.evaluator_list)
 *           {
 *             evaluator->set_refinement_cycle(step);
 *             solver->postprocess(*evaluator);
 *           }
 * 
 * 
 *         if (solver->n_dofs() < descriptor.max_degrees_of_freedom)
 *           solver->refine_grid();
 *         else
 *           break;
 *       }
 * 
 * @endcode
 * 
 * Clean up the screen after the loop has run:
 * 
 * @code
 *     std::cout << std::endl;
 *   }
 * 
 * } // namespace Step14
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * Here finally comes the main function. It drives the whole process by
 * specifying a set of parameters to be used for the simulation (polynomial
 * degrees, evaluation and dual functionals, etc), and passes them packed into
 * a structure to the frame work class above.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step14;
 * 
 * @endcode
 * 
 * Describe the problem we want to solve here by passing a descriptor
 * object to the function doing the rest of the work:
 * 
 * @code
 *       const unsigned int                 dim = 2;
 *       Framework<dim>::ProblemDescription descriptor;
 * 
 * @endcode
 * 
 * First set the refinement criterion we wish to use:
 * 
 * @code
 *       descriptor.refinement_criterion =
 *         Framework<dim>::ProblemDescription::dual_weighted_error_estimator;
 * @endcode
 * 
 * Here, we could as well have used <code>global_refinement</code> or
 * <code>weighted_kelly_indicator</code>. Note that the information
 * given about dual finite elements, dual functional, etc is only
 * important for the given choice of refinement criterion, and is
 * ignored otherwise.
 * 

 * 
 * Then set the polynomial degrees of primal and dual problem. We choose
 * here bi-linear and bi-quadratic ones:
 * 
 * @code
 *       descriptor.primal_fe_degree = 1;
 *       descriptor.dual_fe_degree   = 2;
 * 
 * @endcode
 * 
 * Then set the description of the test case, i.e. domain, boundary
 * values, and right hand side. These are prepackaged in classes. We
 * take here the description of <code>Exercise_2_3</code>, but you can
 * also use <code>CurvedRidges@<dim@></code>:
 * 
 * @code
 *       descriptor.data =
 *         std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>();
 * 
 * @endcode
 * 
 * Next set first a dual functional, then a list of evaluation
 * objects. We choose as default the evaluation of the value at an
 * evaluation point, represented by the classes
 * <code>PointValueEvaluation</code> in the namespaces of evaluation and
 * dual functional classes. You can also set the
 * <code>PointXDerivativeEvaluation</code> classes for the x-derivative
 * instead of the value at the evaluation point.
 *       

 * 
 * Note that dual functional and evaluation objects should
 * match. However, you can give as many evaluation functionals as you
 * want, so you can have both point value and derivative evaluated after
 * each step.  One such additional evaluation is to output the grid in
 * each step.
 * 
 * @code
 *       const Point<dim> evaluation_point(0.75, 0.75);
 *       descriptor.dual_functional =
 *         std::make_unique<DualFunctional::PointValueEvaluation<dim>>(
 *           evaluation_point);
 * 
 *       Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point);
 *       Evaluation::GridOutput<dim>           postprocessor2("grid");
 * 
 *       descriptor.evaluator_list.push_back(&postprocessor1);
 *       descriptor.evaluator_list.push_back(&postprocessor2);
 * 
 * @endcode
 * 
 * Set the maximal number of degrees of freedom after which we want the
 * program to stop refining the mesh further:
 * 
 * @code
 *       descriptor.max_degrees_of_freedom = 20000;
 * 
 * @endcode
 * 
 * Finally pass the descriptor object to a function that runs the entire
 * solution with it:
 * 
 * @code
 *       Framework<dim>::run(descriptor);
 *     }
 * 
 * @endcode
 * 
 * Catch exceptions to give information about things that failed:
 * 
 * @code
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


<a name="Pointvalues"></a><h3>Point values</h3>



This program offers a lot of possibilities to play around. We can thus
only show a small part of all possible results that can be obtained
with the help of this program. However, you are encouraged to just try
it out, by changing the settings in the main program. Here, we start
by simply letting it run, unmodified:
@code
Refinement cycle: 0
   Number of degrees of freedom=72
   Point value=0.03243
   Estimated error=0.000702385
Refinement cycle: 1
   Number of degrees of freedom=67
   Point value=0.0324827
   Estimated error=0.000888953
Refinement cycle: 2
   Number of degrees of freedom=130
   Point value=0.0329619
   Estimated error=0.000454606
Refinement cycle: 3
   Number of degrees of freedom=307
   Point value=0.0331934
   Estimated error=0.000241254
Refinement cycle: 4
   Number of degrees of freedom=718
   Point value=0.0333675
   Estimated error=7.4912e-05
Refinement cycle: 5
   Number of degrees of freedom=1665
   Point value=0.0334083
   Estimated error=3.69111e-05
Refinement cycle: 6
   Number of degrees of freedom=3975
   Point value=0.033431
   Estimated error=1.54218e-05
Refinement cycle: 7
   Number of degrees of freedom=8934
   Point value=0.0334406
   Estimated error=6.28359e-06
Refinement cycle: 8
   Number of degrees of freedom=21799
   Point value=0.0334444
@endcode


First let's look what the program actually computed. On the seventh
grid, primal and dual numerical solutions look like this (using a
color scheme intended to evoke the snow-capped mountains of
Colorado that the original author of this program now calls
home):
<table align="center">
  <tr>
    <td width="50%">
      <img src="https://www.dealii.org/images/steps/developer/step-14.point-value.solution-7.9.2.png" alt="">
    </td>
    <td width="50%">
      <img src="https://www.dealii.org/images/steps/developer/step-14.point-value.solution-7-dual.9.2.png" alt="">
    </td>
  </tr>
</table>
Apparently, the region at the bottom left is so unimportant for the
point value evaluation at the top right that the grid is left entirely
unrefined there, even though the solution has singularities at the inner
corner of that cell! Due
to the symmetry in right hand side and domain, the solution should
actually look like at the top right in all four corners, but the mesh
refinement criterion involving the dual solution chose to refine them
differently -- because we said that we really only care about a single
function value somewhere at the top right.



Here are some of the meshes that are produced in refinement cycles 0,
2, 4 (top row), and 5, 7, and 8 (bottom row):

<table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-0.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-2.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-4.9.2.png" alt="" width="100%"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-5.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-7.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.grid-8.9.2.png" alt="" width="100%"></td>
  </tr>
</table>

Note the subtle interplay between resolving the corner singularities,
and resolving around the point of evaluation. It will be rather
difficult to generate such a mesh by hand, as this would involve to
judge quantitatively how much which of the four corner singularities
should be resolved, and to set the weight compared to the vicinity of
the evaluation point.



The program prints the point value and the estimated error in this
quantity. From extrapolating it, we can guess that the exact value is
somewhere close to 0.0334473, plus or minus 0.0000001 (note that we get
almost 6 valid digits from only 22,000 (primal) degrees of
freedom. This number cannot be obtained from the value of the
functional alone, but I have used the assumption that the error
estimator is mostly exact, and extrapolated the computed value plus
the estimated error, to get an approximation of the true
value. Computing with more degrees of freedom shows that this
assumption is indeed valid.



From the computed results, we can generate two graphs: one that shows
the convergence of the error $J(u)-J(u_h)$ (taking the
extrapolated value as correct) in the point value, and the value that
we get by adding up computed value $J(u_h)$ and estimated
error eta (if the error estimator $eta$ were exact, then the value
$J(u_h)+\eta$ would equal the exact point value, and the error
in this quantity would always be zero; however, since the error
estimator is only a - good - approximation to the true error, we can
by this only reduce the size of the error). In this graph, we also
indicate the complexity ${\cal O}(1/N)$ to show that mesh refinement
acts optimal in this case. The second chart compares
true and estimated error, and shows that the two are actually very
close to each other, even for such a complicated quantity as the point
value:


<table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.error.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-value.error-estimation.png" alt="" width="100%"></td>
  </tr>
</table>


<a name="Comparingrefinementcriteria"></a><h3>Comparing refinement criteria</h3>



Since we have accepted quite some effort when using the mesh
refinement driven by the dual weighted error estimator (for solving
the dual problem, and for evaluating the error representation), it is
worth while asking whether that effort was successful. To this end, we
first compare the achieved error levels for different mesh refinement
criteria. To generate this data, simply change the value of the mesh
refinement criterion variable in the main program. The results are
thus (for the weight in the Kelly indicator, we have chosen the
function $1/(r^2+0.1^2)$, where $r$
is the distance to the evaluation point; it can be shown that this is
the optimal weight if we neglect the effects of boundaries):

<img src="https://www.dealii.org/images/steps/developer/step-14.point-value.error-comparison.png" alt="">



Checking these numbers, we see that for global refinement, the error
is proportional to $O(1/(sqrt(N) log(N)))$, and for the dual
estimator $O(1/N)$. Generally speaking, we see that the dual
weighted error estimator is better than the other refinement
indicators, at least when compared with those that have a similarly
regular behavior. The Kelly indicator produces smaller errors, but
jumps about the picture rather irregularly, with the error also
changing signs sometimes. Therefore, its behavior does not allow to
extrapolate the results to larger values of N. Furthermore, if we
trust the error estimates of the dual weighted error estimator, the
results can be improved by adding the estimated error to the computed
values. In terms of reliability, the weighted estimator is thus better
than the Kelly indicator, although the latter sometimes produces
smaller errors.



<a name="Evaluationofpointstresses"></a><h3>Evaluation of point stresses</h3>



Besides evaluating the values of the solution at a certain point, the
program also offers the possibility to evaluate the x-derivatives at a
certain point, and also to tailor mesh refinement for this. To let the
program compute these quantities, simply replace the two occurrences of
<code>PointValueEvaluation</code> in the main function by
<code>PointXDerivativeEvaluation</code>, and let the program run:
@code
Refinement cycle: 0
   Number of degrees of freedom=72
   Point x-derivative=-0.0719397
   Estimated error=-0.0126173
Refinement cycle: 1
   Number of degrees of freedom=61
   Point x-derivative=-0.0707956
   Estimated error=-0.00774316
Refinement cycle: 2
   Number of degrees of freedom=131
   Point x-derivative=-0.0568671
   Estimated error=-0.00313426
Refinement cycle: 3
   Number of degrees of freedom=247
   Point x-derivative=-0.053033
   Estimated error=-0.00136114
Refinement cycle: 4
   Number of degrees of freedom=532
   Point x-derivative=-0.0526429
   Estimated error=-0.000558868
Refinement cycle: 5
   Number of degrees of freedom=1267
   Point x-derivative=-0.0526955
   Estimated error=-0.000220116
Refinement cycle: 6
   Number of degrees of freedom=2864
   Point x-derivative=-0.0527495
   Estimated error=-9.46731e-05
Refinement cycle: 7
   Number of degrees of freedom=6409
   Point x-derivative=-0.052785
   Estimated error=-4.21543e-05
Refinement cycle: 8
   Number of degrees of freedom=14183
   Point x-derivative=-0.0528028
   Estimated error=-2.04241e-05
Refinement cycle: 9
   Number of degrees of freedom=29902
   Point x-derivative=-0.052814
@endcode



The solution looks roughly the same as before (the exact solution of
course <em>is</em> the same, only the grid changed a little), but the
dual solution is now different. A close-up around the point of
evaluation shows this:
<table align="center">
  <tr>
    <td width="50%">
      <img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.solution-7-dual.png" alt="">
    </td>
    <td width="50%">
      <img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.solution-7-dual-close-up.png" alt="">
    </td>
</table>
This time, the grids in refinement cycles 0, 5, 6, 7, 8, and 9 look
like this:

<table align="center" width="80%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-0.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-5.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-6.9.2.png" alt="" width="100%"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-7.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-8.9.2.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.grid-9.9.2.png" alt="" width="100%"></td>
  </tr>
</table>

Note the asymmetry of the grids compared with those we obtained for
the point evaluation. This is due to the fact that the domain and the primal
solution may be symmetric about the diagonal, but the $x$-derivative is
not, and the latter enters the refinement criterion.



Then, it is interesting to compare actually computed values of the
quantity of interest (i.e. the x-derivative of the solution at one
point) with a reference value of -0.0528223... plus or minus
0.0000005. We get this reference value by computing on finer grid after
some more mesh refinements, with approximately 130,000 cells.
Recall that if the error is $O(1/N)$ in the optimal case, then
taking a mesh with ten times more cells gives us one additional digit
in the result.



In the left part of the following chart, you again see the convergence
of the error towards this extrapolated value, while on the right you
see a comparison of true and estimated error:

<table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.error.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.point-derivative.error-estimation.png" alt="" width="100%"></td>
  </tr>
</table>

After an initial phase where the true error changes its sign, the
estimated error matches it quite well, again. Also note the dramatic
improvement in the error when using the estimated error to correct the
computed value of $J(u_h)$.



<a name="step13revisited"></a><h3>step-13 revisited</h3>



If instead of the <code>Exercise_2_3</code> data set, we choose
<code>CurvedRidges</code> in the main function, and choose $(0.5,0.5)$
as the evaluation point, then we can redo the
computations of the previous example program, to compare whether the
results obtained with the help of the dual weighted error estimator
are better than those we had previously.



First, the meshes after 9 adaptive refinement cycles obtained with
the point evaluation and derivative evaluation refinement
criteria, respectively, look like this:

<table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.step-13.point-value.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.step-13.point-derivative.png" alt="" width="100%"></td>
  </tr>
</table>

The features of the solution can still be seen in the mesh, but since the
solution is smooth, the singularities of the dual solution entirely
dominate the mesh refinement criterion, and lead to strongly
concentrated meshes. The solution after the seventh refinement step looks
like the following:

<table width="40%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-14.step-13.solution-7.9.2.png" alt="" width="100%"></td>
  </tr>
</table>

Obviously, the solution is worse at some places, but the mesh
refinement process should have taken care that these places are not
important for computing the point value.




The next point is to compare the new (duality based) mesh refinement
criterion with the old ones. These are the results:

<img src="https://www.dealii.org/images/steps/developer/step-14.step-13.error-comparison.png" alt="">



The results are, well, somewhat mixed. First, the Kelly indicator
disqualifies itself by its unsteady behavior, changing the sign of the
error several times, and with increasing errors under mesh
refinement. The dual weighted error estimator has a monotone decrease
in the error, and is better than the weighted Kelly and global
refinement, but the margin is not as large as expected. This is, here,
due to the fact the global refinement can exploit the regular
structure of the meshes around the point of evaluation, which leads to
a better order of convergence for the point error. However, if we had
a mesh that is not locally rectangular, for example because we had to
approximate curved boundaries, or if the coefficients were not
constant, then this advantage of globally refinement meshes would
vanish, while the good performance of the duality based estimator
would remain.




<a name="Conclusionsandoutlook"></a><h3>Conclusions and outlook</h3>



The results here are not too clearly indicating the superiority of the
dual weighted error estimation approach for mesh refinement over other
mesh refinement criteria, such as the Kelly indicator. This is due to
the relative simplicity of the shown applications. If you are not
convinced yet that this approach is indeed superior, you are invited
to browse through the literature indicated in the introduction, where
plenty of examples are provided where the dual weighted approach can
reduce the necessary numerical work by orders of magnitude, making
this the only way to compute certain quantities to reasonable
accuracies at all.



Besides the objections you may raise against its use as a mesh
refinement criterion, consider that accurate knowledge of the error in
the quantity one might want to compute is of great use, since we can
stop computations when we are satisfied with the accuracy. Using more
traditional approaches, it is very difficult to get accurate estimates
for arbitrary quantities, except for, maybe, the error in the energy
norm, and we will then have no guarantee that the result we computed
satisfies any requirements on its accuracy. Also, as was shown for the
evaluation of point values and derivatives, the error estimate can be
used to extrapolate the results, yielding much higher accuracy in the
quantity we want to know.



Leaving these mathematical considerations, we tried to write the
program in a modular way, such that implementing another test case, or
another evaluation and dual functional is simple. You are encouraged
to take the program as a basis for your own experiments, and to play a
little.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-14.cc"
*/
