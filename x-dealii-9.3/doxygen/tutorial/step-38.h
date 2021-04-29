/**
@page step_38 The step-38 tutorial program
This tutorial depends on step-34.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Testcase">Testcase</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeLaplaceBeltramiProblemcodeclasstemplate">The <code>LaplaceBeltramiProblem</code> class template</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeLaplaceBeltramiProblemcodeclass">Implementation of the <code>LaplaceBeltramiProblem</code> class</a>
      <ul>
        <li><a href="#LaplaceBeltramiProblemmake_grid_and_dofs">LaplaceBeltramiProblem::make_grid_and_dofs</a>
        <li><a href="#LaplaceBeltramiProblemassemble_system">LaplaceBeltramiProblem::assemble_system</a>
        <li><a href="#LaplaceBeltramiProblemsolve">LaplaceBeltramiProblem::solve</a>
        <li><a href="#LaplaceBeltramiProblemoutput_result">LaplaceBeltramiProblem::output_result</a>
        <li><a href="#LaplaceBeltramiProblemcompute_error">LaplaceBeltramiProblem::compute_error</a>
        <li><a href="#LaplaceBeltramiProblemrun">LaplaceBeltramiProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Andrea Bonito and M. Sebastian Pauletti,
with editing and writing by Wolfgang Bangerth.
<br>
This material is based upon work supported by the National Science
Foundation under Grant No. DMS-0914977. Any opinions, findings and conclusions
or recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the National Science Foundation
(NSF).
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


In this example, we show how to solve a partial differential equation (PDE)
on a codimension one surface $\Gamma \subset \mathbb R^3$
made of quadrilaterals, i.e. on a surface in 3d or a line in 2d.
We focus on the following elliptic second order PDE
@f{align*}
-\Delta_\Gamma u &= f \qquad \text{on } \qquad \Gamma,\\
u  &= g \qquad \text{on} \qquad \partial \Gamma,
@f}
which generalized the Laplace equation we have previously solved in several of
the early tutorial programs. Our implementation is based on step-4. step-34
also solves problems on lower dimensional surfaces; however, there we only
consider integral equations that do not involve derivatives on the solution
variable, while here we actually have to investigate what it means to take
derivatives of a function only defined on a (possibly curved) surface.

In order to define the above operator, we start by introducing some notations.
Let $\mathbf x_S:\hat S \rightarrow S$ be a parameterization of
a surface $S$ from a reference element $\hat S \subset \mathbb R^2$,
i.e. each point $\hat{\mathbf x}\in\hat S$ induces a point ${\mathbf
  x}_S(\hat{\mathbf x}) \in S$. Then let
@f[
G_S\dealcoloneq (D \mathbf{x}_S)^T \ D \mathbf{x}_S
@f]
denotes the corresponding first fundamental form, where $D
\mathbf{x}_S=\left(\frac{\partial x_{S,i}(\hat{\mathbf x})}{\partial \hat x_j}\right)_{ij}$ is the
derivative (Jacobian) of the mapping.
In the following, $S$ will be either the entire surface $\Gamma$ or,
more convenient for the finite element method, any face $S \in
{\mathbb T}$, where ${\mathbb T}$ is a partition (triangulation) of $\Gamma$
constituted of quadrilaterals.
We are now in position to define the tangential gradient of a function $v : S \rightarrow \mathbb
R$ by
@f[
(\nabla_S v)\circ \mathbf x_S \dealcoloneq  D \mathbf x_S \ G_S^{-1} \ \nabla (v \circ \mathbf x_S).
@f]
The surface Laplacian (also called the Laplace-Beltrami operator) is then
defined as  $\Delta_S \dealcoloneq \nabla_S \cdot \nabla_S$.
Note that an alternate way to compute the surface gradient on smooth surfaces $\Gamma$ is
@f[
\nabla_S v = \nabla \tilde v - \mathbf n (\mathbf n \cdot \nabla \tilde v),
@f]
where $\tilde v$ is a "smooth" extension of $v$ in a tubular neighborhood of $\Gamma$ and
$\mathbf n$ is the normal of $\Gamma$.
Since $\Delta_S = \nabla_S \cdot \nabla_S$, we deduce
@f[
\Delta_S v = \Delta \tilde v - \mathbf n^T \ D^2 \tilde v \ \mathbf n - (\mathbf n \cdot \nabla \tilde v) (\nabla \cdot \mathbf n - \mathbf n^T \ D \mathbf n \ \mathbf n ).
@f]
Worth mentioning, the term $\nabla \cdot \mathbf n - \mathbf n \ D \mathbf n \ \mathbf n$ appearing in the above expression is the total curvature of the surface (sum of principal curvatures).

As usual, we are only interested in weak solutions for which we can use $C^0$
finite elements (rather than requiring $C^1$ continuity as for strong
solutions). We therefore resort to the weak formulation
@f[
\int_\Gamma \nabla_\Gamma u \cdot
\nabla_\Gamma v = \int_\Gamma f \ v  \qquad \forall v \in H^1_0(\Gamma)
@f]
and take advantage of the partition ${\mathbb T}$ to further write
@f[
\sum_{K\in  {\mathbb T}}\int_K \nabla_{K} u \cdot \nabla_{K} v = \sum_{K\in
  {\mathbb T}} \int_K f \ v  \qquad \forall v \in H^1_0(\Gamma).
@f]
Moreover, each integral in the above expression is computed in the reference
element $\hat K \dealcoloneq [0,1]^2$
so that
@f{align*}
\int_{K} \nabla_{K} u \cdot \nabla_{K} v
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} (D \mathbf
  x_K)^T D \mathbf x_K G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}
\\
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}
@f}
and
@f[
\int_{K} f \ v = \int_{\hat K} (f \circ \mathbf x_K) (v \circ \mathbf
x_K)  \sqrt{\det
    (G_K)}.
@f]
Finally, we use a quadrature formula defined by points $\{p_l\}_{l=1}^N\subset
\hat K$ and weights $\{w_l\}_{l=1}^N \subset \mathbb R^+_*$ to
evaluate the above integrals and
obtain
@f[\int_{K} \nabla_{K} u \cdot \nabla_{K} v \approx \sum_{l=1}^N
 (\nabla (u \circ \mathbf x_K)(p_l))^T G^{-1}(p_l)  \nabla (v \circ \mathbf x_K)
(p_l) \sqrt{\det (G(p_l))} \ w_l
@f]
and
@f[
\int_{K} f \ v \approx \sum_{l=1}^N (f \circ \mathbf x_K)(p_l) \ (v \circ \mathbf x_K)(p_l) \sqrt{\det (G(p_l))} \ w_l.
@f]


Fortunately, deal.II has already all the tools to compute the above
expressions.
In fact, they barely differ from the ways in which we solve the usual
Laplacian, only requiring the surface coordinate mapping to be provided in the
constructor of the FEValues class.
This surface description given, in the codimension one surface case, the two
routines FEValues::shape_grad and FEValues::JxW
return
@f{align*}
\text{FEValues::shape\_grad}(i,l)&=D \mathbf x_K(p_l) G^{-1}(p_l)D(\varphi_i \circ \mathbf x_K)
  (p_l)
\\
\text{FEValues::JxW}(l) &=  \sqrt{\det (G(p_l))} \ w_l.
@f}
This provides exactly the terms we need for our computations.

On a more general note, details for the finite element approximation on
surfaces can be found for instance in
[Dziuk, in Partial differential equations and calculus of
variations 1357, Lecture Notes in Math., 1988],
[Demlow, SIAM J. Numer. Anal.  47(2), 2009]
and
[Bonito, Nochetto, and Pauletti, SIAM J. Numer. Anal. 48(5), 2010].



<a name="Testcase"></a><h3>Testcase</h3>


In general when you want to test numerically the accuracy and/or order of
convergence of an algorithm you need to provide an exact solution. The usual
trick is to pick a function that we want to be the solution, then apply the
differential operator to it that defines a forcing term for the right hand
side. This is what we do in this example. In the current case, the form of the
domain is obviously also essential.

We produce one test case for a 2d problem and another one for 3d:

<ul>
<li>
  In 2d, let's choose as domain a half circle. On this domain, we choose the
  function $u(\mathbf x)=-2x_1x_2$ as the solution. To compute the right hand
  side, we have to compute the surface Laplacian of the
  solution function. There are (at least) two ways to do that. The first one
  is to project away the normal derivative as described above using the natural extension of $u(\mathbf x)$ (still denoted by $u$) over $\mathbb R^d$, i.e. to compute
  @f[
    -\Delta_\Gamma u =  \Delta u - \mathbf n^T \ D^2 u \ \mathbf n - (\mathbf n \cdot \nabla u)\ \kappa,
  @f]
  where $\kappa$ is the total curvature of $\Gamma$.
  Since we are on the unit circle, $\mathbf n=\mathbf x$ and $\kappa = 1$ so that
  @f[
    -\Delta_\Gamma u = -8 x_1x_2.
  @f]

  A somewhat simpler way, at least for the current case of a curve in
  two-dimensional space, is to note that we can map the interval $t \in
  [0,\pi]$ onto the domain $\Omega$ using the transformation
  $\mathbf x(t)= \left(\begin{array}{c} \cos t \\ \sin t \end{array}\right)$.
  At position $\mathbf x=\mathbf x(t)$, the value of the solution is then
  $u(\mathbf x(t)) = -2\cos t \sin t$.
  Taking into account that the transformation is length preserving, i.e. a
  segment of length $dt$ is mapped onto a piece of curve of exactly the same
  length, the tangential Laplacian then satisfies
  @f{align*}
    \Delta_\Gamma u
    &= \frac{d^2}{dt^2}(-2\cos t \sin t)
    = -2 \frac{d}{dt}(-\sin^2 t + \cos^2 t)
    = -2 (-2 \sin t \cos t - 2 \cos t \sin t)
    \\
    &= 8 \sin t \cos t
    \\
    &= 8 x_1x_2,
  @f}
  which is of course the same result as we had above.
</li>
<li>
  In 3d, the domain is again half of the surface of the unit ball, i.e. a half
  sphere or dome. We choose $u(\mathbf x)=-2\sin(\pi x_1)\cos(\pi x_2)e^z$ as
  the solution. We can compute the right hand side of the
  equation, $f=-\Delta_\Gamma u$, in the same way as the method above (with $\kappa = 2$), yielding an
  awkward and lengthy expression. You can find the full expression in the
  source code.
</li>
</ul>

In the program, we will also compute the $H^1$ seminorm error of the
solution. Since the solution function and its numerical approximation are only
defined on the manifold, the obvious definition of this error functional is
$| e |_{H^1(\Gamma)}
  = | \nabla_\Gamma e |_{L_2(\Gamma)}
  = \left( \int_\Gamma | \nabla_\Gamma (u-u_h) |^2 \right)^{1/2}$. This requires us to provide the
<i>tangential</i> gradient $\nabla_\Gamma u$ to the function VectorTools::integrate_difference
(first introduced in step-7), which we
will do by implementing the function <code>Solution::gradient</code> in the
program below.


<a name="Implementation"></a><h3>Implementation</h3>


If you've read through step-4 and understand the discussion above of how
solution and right hand side correspond to each other, you will be immediately
familiar with this program as well. In fact, there are only two things that
are of significance:

- The way we generate the mesh that triangulates the computational domain.

- The way we use Mapping objects to describe that the domain on which we solve
  the partial differential equation is not planar but in fact curved.

Mapping objects were already introduced in step-10 and step-11 and as
explained there, there is usually not a whole lot you have to know about how
they work as long as you have a working description of how the boundary
looks. In essence, we will simply declare an appropriate object of type
MappingQ that will automatically obtain the boundary description from the
Triangulation. The mapping object will then be passed to the appropriate
functions, and we will get a boundary description for half circles or half
spheres that is predefined in the library.

The rest of the program follows closely step-4 and, as far as computing the
error, step-7. Some aspects of this program, in particular the use of two
template arguments on the classes Triangulation, DoFHandler, and similar, are
already described in detail in step-34; you may wish to read through this
tutorial program as well.
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
 * If you've read through step-4 and step-7, you will recognize that we have
 * used all of the following include files there already. Consequently, we
 * will not explain their meaning here again.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/solver_control.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * namespace Step38
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceBeltramiProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceBeltramiProblem</code> class template</h3>
 * 

 * 
 * This class is almost exactly similar to the <code>LaplaceProblem</code>
 * class in step-4.
 * 

 * 
 * The essential differences are these:
 *   

 * 
 * - The template parameter now denotes the dimensionality of the embedding
 * space, which is no longer the same as the dimensionality of the domain
 * and the triangulation on which we compute. We indicate this by calling
 * the parameter @p spacedim, and introducing a constant @p dim equal to
 * the dimensionality of the domain -- here equal to
 * <code>spacedim-1</code>.
 * - All member variables that have geometric aspects now need to know about
 * both their own dimensionality as well as that of the embedding
 * space. Consequently, we need to specify both of their template
 * parameters one for the dimension of the mesh @p dim, and the other for
 * the dimension of the embedding space, @p spacedim. This is exactly what
 * we did in step-34, take a look there for a deeper explanation.
 * - We need an object that describes which kind of mapping to use from the
 * reference cell to the cells that the triangulation is composed of. The
 * classes derived from the Mapping base class do exactly this. Throughout
 * most of deal.II, if you don't do anything at all, the library assumes
 * that you want an object of kind MappingQ1 that uses a (bi-, tri-)linear
 * mapping. In many cases, this is quite sufficient, which is why the use
 * of these objects is mostly optional: for example, if you have a
 * polygonal two-dimensional domain in two-dimensional space, a bilinear
 * mapping of the reference cell to the cells of the triangulation yields
 * an exact representation of the domain. If you have a curved domain, one
 * may want to use a higher order mapping for those cells that lie at the
 * boundary of the domain -- this is what we did in step-11, for
 * example. However, here we have a curved domain, not just a curved
 * boundary, and while we can approximate it with bilinearly mapped cells,
 * it is really only prudent to use a higher order mapping for all
 * cells. Consequently, this class has a member variable of type MappingQ;
 * we will choose the polynomial degree of the mapping equal to the
 * polynomial degree of the finite element used in the computations to
 * ensure optimal approximation, though this iso-parametricity is not
 * required.
 * 
 * @code
 *   template <int spacedim>
 *   class LaplaceBeltramiProblem
 *   {
 *   public:
 *     LaplaceBeltramiProblem(const unsigned degree = 2);
 *     void run();
 * 
 *   private:
 *     static constexpr unsigned int dim = spacedim - 1;
 * 
 *     void make_grid_and_dofs();
 *     void assemble_system();
 *     void solve();
 *     void output_results() const;
 *     void compute_error() const;
 * 
 * 
 *     Triangulation<dim, spacedim> triangulation;
 *     FE_Q<dim, spacedim>          fe;
 *     DoFHandler<dim, spacedim>    dof_handler;
 *     MappingQ<dim, spacedim>      mapping;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * Next, let us define the classes that describe the exact solution and the
 * right hand sides of the problem. This is in analogy to step-4 and step-7
 * where we also defined such objects. Given the discussion in the
 * introduction, the actual formulas should be self-explanatory. A point of
 * interest may be how we define the value and gradient functions for the 2d
 * and 3d cases separately, using explicit specializations of the general
 * template. An alternative to doing it this way might have been to define
 * the general template and have a <code>switch</code> statement (or a
 * sequence of <code>if</code>s) for each possible value of the spatial
 * dimension.
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
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
 *   template <>
 *   double Solution<2>::value(const Point<2> &p, const unsigned int) const
 *   {
 *     return (-2. * p(0) * p(1));
 *   }
 * 
 * 
 *   template <>
 *   Tensor<1, 2> Solution<2>::gradient(const Point<2> &p,
 *                                      const unsigned int) const
 *   {
 *     Tensor<1, 2> return_value;
 *     return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0));
 *     return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1));
 * 
 *     return return_value;
 *   }
 * 
 * 
 *   template <>
 *   double Solution<3>::value(const Point<3> &p, const unsigned int) const
 *   {
 *     return (std::sin(numbers::PI * p(0)) * std::cos(numbers::PI * p(1)) *
 *             exp(p(2)));
 *   }
 * 
 * 
 *   template <>
 *   Tensor<1, 3> Solution<3>::gradient(const Point<3> &p,
 *                                      const unsigned int) const
 *   {
 *     using numbers::PI;
 * 
 *     Tensor<1, 3> return_value;
 * 
 *     return_value[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     return_value[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     return_value[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     return return_value;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 *   template <>
 *   double RightHandSide<2>::value(const Point<2> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     return (-8. * p(0) * p(1));
 *   }
 * 
 * 
 *   template <>
 *   double RightHandSide<3>::value(const Point<3> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     using numbers::PI;
 * 
 *     Tensor<2, 3> hessian;
 * 
 *     hessian[0][0] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[1][1] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[2][2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     hessian[0][1] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     hessian[1][0] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 * 
 *     hessian[0][2] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[2][0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     hessian[1][2] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     hessian[2][1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 * 
 *     Tensor<1, 3> gradient;
 *     gradient[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     gradient[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     gradient[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     Point<3> normal = p;
 *     normal /= p.norm();
 * 
 *     return (-trace(hessian) + 2 * (gradient * normal) +
 *             (hessian * normal) * normal);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeLaplaceBeltramiProblemcodeclass"></a> 
 * <h3>Implementation of the <code>LaplaceBeltramiProblem</code> class</h3>
 * 

 * 
 * The rest of the program is actually quite unspectacular if you know
 * step-4. Our first step is to define the constructor, setting the
 * polynomial degree of the finite element and mapping, and associating the
 * DoF handler to the triangulation:
 * 
 * @code
 *   template <int spacedim>
 *   LaplaceBeltramiProblem<spacedim>::LaplaceBeltramiProblem(
 *     const unsigned degree)
 *     : fe(degree)
 *     , dof_handler(triangulation)
 *     , mapping(degree)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemmake_grid_and_dofs"></a> 
 * <h4>LaplaceBeltramiProblem::make_grid_and_dofs</h4>
 * 

 * 
 * The next step is to create the mesh, distribute degrees of freedom, and
 * set up the various variables that describe the linear system. All of
 * these steps are standard with the exception of how to create a mesh that
 * describes a surface. We could generate a mesh for the domain we are
 * interested in, generate a triangulation using a mesh generator, and read
 * it in using the GridIn class. Or, as we do here, we generate the mesh
 * using the facilities in the GridGenerator namespace.
 *   

 * 
 * In particular, what we're going to do is this (enclosed between the set
 * of braces below): we generate a <code>spacedim</code> dimensional mesh
 * for the half disk (in 2d) or half ball (in 3d), using the
 * GridGenerator::half_hyper_ball function. This function sets the boundary
 * indicators of all faces on the outside of the boundary to zero for the
 * ones located on the perimeter of the disk/ball, and one on the straight
 * part that splits the full disk/ball into two halves. The next step is the
 * main point: The GridGenerator::extract_boundary_mesh function creates a
 * mesh that consists of those cells that are the faces of the previous mesh,
 * i.e. it describes the <i>surface</i> cells of the original (volume)
 * mesh. However, we do not want all faces: only those on the perimeter of
 * the disk or ball which carry boundary indicator zero; we can select these
 * cells using a set of boundary indicators that we pass to
 * GridGenerator::extract_boundary_mesh.
 *   

 * 
 * There is one point that needs to be mentioned. In order to refine a
 * surface mesh appropriately if the manifold is curved (similarly to
 * refining the faces of cells that are adjacent to a curved boundary), the
 * triangulation has to have an object attached to it that describes where
 * new vertices should be located. If you don't attach such a boundary
 * object, they will be located halfway between existing vertices; this is
 * appropriate if you have a domain with straight boundaries (e.g. a
 * polygon) but not when, as here, the manifold has curvature. So for things
 * to work properly, we need to attach a manifold object to our (surface)
 * triangulation, in much the same way as we've already done in 1d for the
 * boundary. We create such an object and attach it to the triangulation.
 *   

 * 
 * The final step in creating the mesh is to refine it a number of
 * times. The rest of the function is the same as in previous tutorial
 * programs.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs()
 *   {
 *     {
 *       Triangulation<spacedim> volume_mesh;
 *       GridGenerator::half_hyper_ball(volume_mesh);
 * 
 *       std::set<types::boundary_id> boundary_ids;
 *       boundary_ids.insert(0);
 * 
 *       GridGenerator::extract_boundary_mesh(volume_mesh,
 *                                            triangulation,
 *                                            boundary_ids);
 *     }
 *     triangulation.set_all_manifold_ids(0);
 *     triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());
 * 
 *     triangulation.refine_global(4);
 * 
 *     std::cout << "Surface mesh has " << triangulation.n_active_cells()
 *               << " cells." << std::endl;
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::cout << "Surface mesh has " << dof_handler.n_dofs()
 *               << " degrees of freedom." << std::endl;
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
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
 * <a name="LaplaceBeltramiProblemassemble_system"></a> 
 * <h4>LaplaceBeltramiProblem::assemble_system</h4>
 * 

 * 
 * The following is the central function of this program, assembling the
 * matrix that corresponds to the surface Laplacian (Laplace-Beltrami
 * operator). Maybe surprisingly, it actually looks exactly the same as for
 * the regular Laplace operator discussed in, for example, step-4. The key
 * is that the FEValues::shape_grad() function does the magic: It returns
 * the surface gradient $\nabla_K \phi_i(x_q)$ of the $i$th shape function
 * at the $q$th quadrature point. The rest then does not need any changes
 * either:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const QGauss<dim>       quadrature_formula(2 * fe.degree);
 *     FEValues<dim, spacedim> fe_values(mapping,
 *                                       fe,
 *                                       quadrature_formula,
 *                                       update_values | update_gradients |
 *                                         update_quadrature_points |
 *                                         update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<double>                  rhs_values(n_q_points);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     RightHandSide<spacedim> rhs;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 *         rhs.value_list(fe_values.get_quadrature_points(), rhs_values);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *               cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
 *                                    fe_values.shape_grad(j, q_point) *
 *                                    fe_values.JxW(q_point);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             cell_rhs(i) += fe_values.shape_value(i, q_point) *
 *                            rhs_values[q_point] * fe_values.JxW(q_point);
 * 
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
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 0, Solution<spacedim>(), boundary_values);
 * 
 *     MatrixTools::apply_boundary_values(
 *       boundary_values, system_matrix, solution, system_rhs, false);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemsolve"></a> 
 * <h4>LaplaceBeltramiProblem::solve</h4>
 * 

 * 
 * The next function is the one that solves the linear system. Here, too, no
 * changes are necessary:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::solve()
 *   {
 *     SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm());
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
 * 
 * <a name="LaplaceBeltramiProblemoutput_result"></a> 
 * <h4>LaplaceBeltramiProblem::output_result</h4>
 * 

 * 
 * This is the function that generates graphical output from the
 * solution. Most of it is boilerplate code, but there are two points worth
 * pointing out:
 *   

 * 
 * - The DataOut::add_data_vector() function can take two kinds of vectors:
 * Either vectors that have one value per degree of freedom defined by the
 * DoFHandler object previously attached via DataOut::attach_dof_handler();
 * and vectors that have one value for each cell of the triangulation, for
 * example to output estimated errors for each cell. Typically, the
 * DataOut class knows to tell these two kinds of vectors apart: there are
 * almost always more degrees of freedom than cells, so we can
 * differentiate by the two kinds looking at the length of a vector. We
 * could do the same here, but only because we got lucky: we use a half
 * sphere. If we had used the whole sphere as domain and $Q_1$ elements,
 * we would have the same number of cells as vertices and consequently the
 * two kinds of vectors would have the same number of elements. To avoid
 * the resulting confusion, we have to tell the DataOut::add_data_vector()
 * function which kind of vector we have: DoF data. This is what the third
 * argument to the function does.
 * - The DataOut::build_patches() function can generate output that subdivides
 * each cell so that visualization programs can resolve curved manifolds
 * or higher polynomial degree shape functions better. We here subdivide
 * each element in each coordinate direction as many times as the
 * polynomial degree of the finite element in use.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::output_results() const
 *   {
 *     DataOut<dim, DoFHandler<dim, spacedim>> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(
 *       solution,
 *       "solution",
 *       DataOut<dim, DoFHandler<dim, spacedim>>::type_dof_data);
 *     data_out.build_patches(mapping, mapping.get_degree());
 * 
 *     const std::string filename =
 *       "solution-" + std::to_string(spacedim) + "d.vtk";
 *     std::ofstream output(filename);
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemcompute_error"></a> 
 * <h4>LaplaceBeltramiProblem::compute_error</h4>
 * 

 * 
 * This is the last piece of functionality: we want to compute the error in
 * the numerical solution. It is a verbatim copy of the code previously
 * shown and discussed in step-7. As mentioned in the introduction, the
 * <code>Solution</code> class provides the (tangential) gradient of the
 * solution. To avoid evaluating the error only a superconvergence points,
 * we choose a quadrature rule of sufficiently high order.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::compute_error() const
 *   {
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Solution<spacedim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(2 * fe.degree + 1),
 *                                       VectorTools::H1_norm);
 * 
 *     double h1_error = VectorTools::compute_global_error(triangulation,
 *                                                         difference_per_cell,
 *                                                         VectorTools::H1_norm);
 *     std::cout << "H1 error = " << h1_error << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemrun"></a> 
 * <h4>LaplaceBeltramiProblem::run</h4>
 * 

 * 
 * The last function provides the top-level logic. Its contents are
 * self-explanatory:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::run()
 *   {
 *     make_grid_and_dofs();
 *     assemble_system();
 *     solve();
 *     output_results();
 *     compute_error();
 *   }
 * } // namespace Step38
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * The remainder of the program is taken up by the <code>main()</code>
 * function. It follows exactly the general layout first introduced in step-6
 * and used in all following tutorial programs:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step38;
 * 
 *       LaplaceBeltramiProblem<3> laplace_beltrami;
 *       laplace_beltrami.run();
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


When you run the program, the following output should be printed on screen:

@verbatim
Surface mesh has 1280 cells.
Surface mesh has 5185 degrees of freedom.
H1 error = 0.0217136
@endverbatim


By playing around with the number of global refinements in the
<code>LaplaceBeltrami::make_grid_and_dofs</code> function you increase or decrease mesh
refinement. For example, doing one more refinement and only running the 3d surface
problem yields the following
output:

@verbatim
Surface mesh has 5120 cells.
Surface mesh has 20609 degrees of freedom.
H1 error = 0.00543481
@endverbatim

This is what we expect: make the mesh size smaller by a factor of two and the
error goes down by a factor of four (remember that we use bi-quadratic
elements). The full sequence of errors from one to five refinements looks like
this, neatly following the theoretically predicted pattern:
@verbatim
0.339438
0.0864385
0.0217136
0.00543481
0.00135913
@endverbatim

Finally, the program produces graphical output that we can visualize. Here is
a plot of the results:

<img src="https://www.dealii.org/images/steps/developer/step-38.solution-3d.png" alt="">

The program also works for 1d curves in 2d, not just 2d surfaces in 3d. You
can test this by changing the template argument in <code>main()</code> like
so:
@code
      LaplaceBeltramiProblem<2> laplace_beltrami;
@endcode
The domain is a curve in 2d, and we can visualize the solution by using the
third dimension (and color) to denote the value of the function $u(x)$. This
then looks like so (the white curve is the domain, the colored curve is the
solution extruded into the third dimension, clearly showing the change in sign
as the curve moves from one quadrant of the domain into the adjacent one):

<img src="https://www.dealii.org/images/steps/developer/step-38.solution-2d.png" alt="">


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


Computing on surfaces only becomes interesting if the surface is more
interesting than just a half sphere. To achieve this, deal.II can read
meshes that describe surfaces through the usual GridIn class. Or, in case you
have an analytic description, a simple mesh can sometimes be stretched and
bent into a shape we are interested in.

Let us consider a relatively simple example: we take the half sphere we used
before, we stretch it by a factor of 10 in the z-direction, and then we jumble
the x- and y-coordinates a bit. Let's show the computational domain and the
solution first before we go into details of the implementation below:

<img src="https://www.dealii.org/images/steps/developer/step-38.warp-1.png" alt="">

<img src="https://www.dealii.org/images/steps/developer/step-38.warp-2.png" alt="">

The way to produce such a mesh is by using the GridTools::transform()
function. It needs a way to transform each individual mesh point to a
different position. Let us here use the following, rather simple function
(remember: stretch in one direction, jumble in the other two):

@code
template <int spacedim>
Point<spacedim> warp(const Point<spacedim> &p)
{
  Point<spacedim> q = p;
  q[spacedim-1] *= 10;

  if (spacedim >= 2)
    q[0] += 2*std::sin(q[spacedim-1]);
  if (spacedim >= 3)
    q[1] += 2*std::cos(q[spacedim-1]);

  return q;
}
@endcode

If we followed the <code>LaplaceBeltrami::make_grid_and_dofs</code> function, we would
extract the half spherical surface mesh as before, warp it into the shape we
want, and refine as often as necessary. This is not quite as simple as we'd
like here, though: refining requires that we have an appropriate manifold
object attached to the triangulation that describes where new vertices of the
mesh should be located upon refinement. I'm sure it's possible to describe
this manifold in a not-too-complicated way by simply undoing the
transformation above (yielding the spherical surface again), finding the
location of a new point on the sphere, and then re-warping the result. But I'm
a lazy person, and since doing this is not really the point here, let's just
make our lives a bit easier: we'll extract the half sphere, refine it as
often as necessary, get rid of the object that describes the manifold since we
now no longer need it, and then finally warp the mesh. With the function
above, this would look as follows:

@code
template <int spacedim>
void LaplaceBeltrami<spacedim>::make_grid_and_dofs()
{
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::half_hyper_ball(volume_mesh);

    volume_mesh.refine_global(4);

    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);

    GridGenerator::extract_boundary_mesh(volume_mesh, triangulation,
                                         boundary_ids);
    GridTools::transform(&warp<spacedim>, triangulation);       /* ** */
    std::ofstream x("x"), y("y");
    GridOut().write_gnuplot(volume_mesh, x);
    GridOut().write_gnuplot(triangulation, y);
  }

  std::cout << "Surface mesh has " << triangulation.n_active_cells()
            << " cells."
            << std::endl;
  ...
}
@endcode

Note that the only essential addition is the line marked with
asterisks. It is worth pointing out one other thing here, though: because we
detach the manifold description from the surface mesh, whenever we use a
mapping object in the rest of the program, it has no curves boundary
description to go on any more. Rather, it will have to use the implicit,
FlatManifold class that is used on all parts of the domain not
explicitly assigned a different manifold object. Consequently, whether we use
MappingQ(2), MappingQ(15) or MappingQ1, each cell of our mesh will be mapped
using a bilinear approximation.

All these drawbacks aside, the resulting pictures are still pretty. The only
other differences to what's in step-38 is that we changed the right hand side
to $f(\mathbf x)=\sin x_3$ and the boundary values (through the
<code>Solution</code> class) to $u(\mathbf x)|_{\partial\Omega}=\cos x_3$. Of
course, we now no longer know the exact solution, so the computation of the
error at the end of <code>LaplaceBeltrami::run</code> will yield a meaningless
number.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-38.cc"
*/
