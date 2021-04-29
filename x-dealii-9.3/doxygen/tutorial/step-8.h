/**
@page step_8 The step-8 tutorial program
This tutorial depends on step-6.

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
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeElasticProblemcodeclasstemplate">The <code>ElasticProblem</code> class template</a>
        <li><a href="#Righthandsidevalues">Right hand side values</a>
        <li><a href="#ThecodeElasticProblemcodeclassimplementation">The <code>ElasticProblem</code> class implementation</a>
      <ul>
        <li><a href="#ElasticProblemElasticProblemconstructor">ElasticProblem::ElasticProblem constructor</a>
        <li><a href="#ElasticProblemsetup_system">ElasticProblem::setup_system</a>
        <li><a href="#ElasticProblemassemble_system">ElasticProblem::assemble_system</a>
        <li><a href="#ElasticProblemsolve">ElasticProblem::solve</a>
        <li><a href="#ElasticProblemrefine_grid">ElasticProblem::refine_grid</a>
        <li><a href="#ElasticProblemoutput_results">ElasticProblem::output_results</a>
        <li><a href="#ElasticProblemrun">ElasticProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
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



In real life, most partial differential equations are really systems
of equations. Accordingly, the solutions are usually
vector-valued. The deal.II library supports such problems (see the
extensive documentation in the @ref vector_valued module), and we will show
that that is mostly rather simple. The only more complicated problems
are in assembling matrix and right hand side, but these are easily
understood as well.

@dealiiVideoLecture{19}

In this tutorial program we will want to solve the
<a href="https://en.wikipedia.org/wiki/Linear_elasticity">elastic equations</a>.
They are an extension to Laplace's equation with a vector-valued solution that
describes the displacement in each space direction of a rigid body
which is subject to a force. Of course, the force is also
vector-valued, meaning that in each point it has a direction and an
absolute value.

One can write the elasticity equations in a number of ways. The one that shows
the symmetry with the Laplace equation in the most obvious way is to write it
as
@f[
  -
  \text{div}\,
  ({\mathbf C} \nabla \mathbf{u})
  =
  \mathbf f,
@f]
where $\mathbf u$ is the vector-valued displacement at each point,
$\mathbf f$ the force, and ${\mathbf C}$ is a rank-4 tensor (i.e., it has four
indices) that encodes the stress-strain relationship -- in essence,
it represents the
<a href="https://en.wikipedia.org/wiki/Hooke%27s_law">"spring constant"</a> in
Hookes law that relates the displacement to the forces. ${\mathbf C}$ will, in many
cases, depend on $\mathbf x$ if the body whose deformation we want to
simulate is composed of different materials.

While the form of the equations above is correct, it is not the way
they are usually derived. In truth, the gradient of the displacement
$\nabla\mathbf u$ (a matrix) has no physical meaning whereas its
symmetrized version,
@f[
\varepsilon(\mathbf u)_{kl} =\frac{1}{2}(\partial_k u_l + \partial_l u_k),
@f]
does and is typically called the "strain". (Here and in the following,
$\partial_k=\frac{\partial}{\partial x_k}$. We will also use the
<a href="https://en.wikipedia.org/wiki/Einstein_notation">Einstein summation
convention</a> that whenever the same index appears twice in an equation,
summation over this index is implied; we will, however, not distinguish
between upper and lower indices.)
With this definition of the strain, the elasticity equations
then read as
@f[
  -
  \text{div}\,
  ({\mathbf C} \varepsilon(\mathbf u))
  =
  \mathbf f,
@f]
which you can think of as the more natural generalization of the Laplace
equation to vector-valued problems. (The form shown first is equivalent to
this form because the tensor ${\mathbf C}$ has certain symmetries, namely that
$C_{ijkl}=C_{ijlk}$, and consequently ${\mathbf C} \varepsilon(\mathbf u)_{kl}
= {\mathbf C} \nabla\mathbf u$.)

One can of course alternatively write these equations in component form:
@f[
  -
  \partial_j (c_{ijkl} \varepsilon_{kl})
  =
  f_i,
  \qquad
  i=1\ldots d.
@f]

In many cases, one knows that the material under consideration is
isotropic, in which case by introduction of the two coefficients
$\lambda$ and $\mu$ the coefficient tensor reduces to
@f[
  c_{ijkl}
  =
  \lambda \delta_{ij} \delta_{kl} +
  \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).
@f]

The elastic equations can then be rewritten in much simpler a form:
@f[
   -
   \nabla \lambda (\nabla\cdot {\mathbf u})
   -
   (\nabla \cdot \mu \nabla) {\mathbf u}
   -
   \nabla\cdot \mu (\nabla {\mathbf u})^T
   =
   {\mathbf f},
@f]
and the respective bilinear form is then
@f[
  a({\mathbf u}, {\mathbf v}) =
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega,
@f]
or also writing the first term a sum over components:
@f[
  a({\mathbf u}, {\mathbf v}) =
  \sum_{k,l}
  \left(
    \lambda \partial_l u_l, \partial_k v_k
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega.
@f]

@note As written, the equations above are generally considered to be the right
description for the displacement of three-dimensional objects if the
displacement is small and we can assume that <a
href="http://en.wikipedia.org/wiki/Hookes_law">Hooke's law</a> is valid. In
that case, the indices $i,j,k,l$ above all run over the set $\{1,2,3\}$ (or,
in the C++ source, over $\{0,1,2\}$). However, as is, the program runs in 2d,
and while the equations above also make mathematical sense in that case, they
would only describe a truly two-dimensional solid. In particular, they are not
the appropriate description of an $x-y$ cross-section of a body infinite in
the $z$ direction; this is in contrast to many other two-dimensional equations
that can be obtained by assuming that the body has infinite extent in
$z$-direction and that the solution function does not depend on the $z$
coordinate. On the other hand, there are equations for two-dimensional models
of elasticity; see for example the Wikipedia article on <a
href="http://en.wikipedia.org/wiki/Infinitesimal_strain_theory#Special_cases">plane
strain</a>, <a
href="http://en.wikipedia.org/wiki/Antiplane_shear">antiplane shear</a> and <a
href="http://en.wikipedia.org/wiki/Plane_stress#Plane_stress">plan stress</a>.

But let's get back to the original problem.
How do we assemble the matrix for such an equation? A very long answer
with a number of different alternatives is given in the documentation of the
@ref vector_valued module. Historically, the solution shown below was the only
one available in the early years of the library. It turns out to also be the
fastest. On the other hand, if a few per cent of compute time do not matter,
there are simpler and probably more intuitive ways to assemble the linear
system than the one discussed below but that weren't available until several
years after this tutorial program was first written; if you are interested in
them, take a look at the @ref vector_valued module.

Let us go back to the question of how to assemble the linear system. The first
thing we need is some knowledge about how the shape functions work in the case
of vector-valued finite elements. Basically, this comes down to the following:
let $n$ be the number of shape functions for the scalar finite element of
which we build the vector element (for example, we will use bilinear functions
for each component of the vector-valued finite element, so the scalar finite
element is the <code>FE_Q(1)</code> element which we have used in previous
examples already, and $n=4$ in two space dimensions). Further, let $N$ be the
number of shape functions for the vector element; in two space dimensions, we
need $n$ shape functions for each component of the vector, so $N=2n$. Then,
the $i$th shape function of the vector element has the form
@f[
  \Phi_i({\mathbf x}) = \varphi_{\text{base}(i)}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)},
@f]
where $e_l$ is the $l$th unit vector, $\text{comp}(i)$ is the function that tells
us which component of $\Phi_i$ is the one that is nonzero (for
each vector shape function, only one component is nonzero, and all others are
zero). $\varphi_{\text{base}(i)}(x)$ describes the space dependence of the shape
function, which is taken to be the $\text{base}(i)$-th shape function of the scalar
element. Of course, while $i$ is in the range $0,\ldots,N-1$, the functions
$\text{comp}(i)$ and $\text{base}(i)$ have the ranges $0,1$ (in 2D) and $0,\ldots,n-1$,
respectively.

For example (though this sequence of shape functions is not
guaranteed, and you should not rely on it),
the following layout could be used by the library:
@f{eqnarray*}
  \Phi_0({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_0({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_1({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_0({\mathbf x})
  \end{array}\right),
  \\
  \Phi_2({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_1({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_3({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_1({\mathbf x})
  \end{array}\right),
  \ldots
@f}
where here
@f[
  \text{comp}(0)=0, \quad  \text{comp}(1)=1, \quad  \text{comp}(2)=0, \quad  \text{comp}(3)=1, \quad  \ldots
@f]
@f[
  \text{base}(0)=0, \quad  \text{base}(1)=0, \quad  \text{base}(2)=1, \quad  \text{base}(3)=1, \quad  \ldots
@f]

In all but very rare cases, you will not need to know which shape function
$\varphi_{\text{base}(i)}$ of the scalar element belongs to a shape function $\Phi_i$
of the vector element. Let us therefore define
@f[
  \phi_i = \varphi_{\text{base}(i)}
@f]
by which we can write the vector shape function as
@f[
  \Phi_i({\mathbf x}) = \phi_{i}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)}.
@f]
You can now safely forget about the function $\text{base}(i)$, at least for the rest
of this example program.

Now using this vector shape functions, we can write the discrete finite
element solution as
@f[
  {\mathbf u}_h({\mathbf x}) =
  \sum_i \Phi_i({\mathbf x})\ U_i
@f]
with scalar coefficients $U_i$. If we define an analog function ${\mathbf v}_h$ as
test function, we can write the discrete problem as follows: Find coefficients
$U_i$ such that
@f[
  a({\mathbf u}_h, {\mathbf v}_h) = ({\mathbf f}, {\mathbf v}_h)
  \qquad
  \forall {\mathbf v}_h.
@f]

If we insert the definition of the bilinear form and the representation of
${\mathbf u}_h$ and ${\mathbf v}_h$ into this formula:
@f{eqnarray*}
  \sum_{i,j}
    U_i V_j
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_\Omega
  \right\}
\\
=
  \sum_j V_j
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_\Omega.
@f}
We note that here and in the following, the indices $k,l$ run over spatial
directions, i.e. $0\le k,l < d$, and that indices $i,j$ run over degrees
of freedom.

The local stiffness matrix on cell $K$ therefore has the following entries:
@f[
  A^K_{ij}
  =
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_K
  \right\},
@f]
where $i,j$ now are local degrees of freedom and therefore $0\le i,j < N$.
In these formulas, we always take some component of the vector shape functions
$\Phi_i$, which are of course given as follows (see their definition):
@f[
  (\Phi_i)_l = \phi_i \delta_{l,\text{comp}(i)},
@f]
with the Kronecker symbol $\delta_{nm}$. Due to this, we can delete some of
the sums over $k$ and $l$:
@f{eqnarray*}
  A^K_{ij}
  &=&
  \sum_{k,l}
  \Bigl\{
  \left(
    \lambda \partial_l \phi_i\ \delta_{l,\text{comp}(i)},
            \partial_k \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
\\
  &\qquad\qquad& +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_l \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
  +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_k \phi_j\ \delta_{l,\text{comp}(j)}
  \right)_K
  \Bigr\}
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \sum_l
  \left(
    \mu \partial_l \phi_i,
        \partial_l \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \left(
    \mu \nabla \phi_i,
        \nabla \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K.
@f}

Likewise, the contribution of cell $K$ to the right hand side vector is
@f{eqnarray*}
  f^K_j
  &=&
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_K
\\
  &=&
  \sum_l
  \left(
    f_l,
    \phi_j \delta_{l,\text{comp}(j)}
  \right)_K
\\
  &=&
  \left(
    f_{\text{comp}(j)},
    \phi_j
  \right)_K.
@f}

This is the form in which we will implement the local stiffness matrix and
right hand side vectors.

As a final note: in the step-17 example program, we will
revisit the elastic problem laid out here, and will show how to solve it in
%parallel on a cluster of computers. The resulting program will thus be able to
solve this problem to significantly higher accuracy, and more efficiently if
this is required. In addition, in step-20, @ref step_21
"step-21", as well as a few other of the later tutorial programs, we will
revisit some vector-valued problems and show a few techniques that may make it
simpler to actually go through all the stuff shown above, with
FiniteElement::system_to_component_index etc.
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
 * As usual, the first few include files are already known, so we will not
 * comment on them further.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/tensor.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * In this example, we need vector-valued finite elements. The support for
 * these can be found in the following include file:
 * 
 * @code
 * #include <deal.II/fe/fe_system.h>
 * @endcode
 * 
 * We will compose the vector-valued finite elements from regular Q1 elements
 * which can be found here, as usual:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * 
 * @endcode
 * 
 * This again is C++:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The last step is as in previous programs. In particular, just like in
 * step-7, we pack everything that's specific to this program into a namespace
 * of its own.
 * 
 * @code
 * namespace Step8
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclasstemplate"></a> 
 * <h3>The <code>ElasticProblem</code> class template</h3>
 * 

 * 
 * The main class is, except for its name, almost unchanged with respect to
 * the step-6 example.
 *   

 * 
 * The only change is the use of a different class for the <code>fe</code>
 * variable: Instead of a concrete finite element class such as FE_Q, we now
 * use a more generic one, FESystem. In fact, FESystem is not really a
 * finite element itself in that it does not implement shape functions of
 * its own. Rather, it is a class that can be used to stack several other
 * elements together to form one vector-valued finite element. In our case,
 * we will compose the vector-valued element of <code>FE_Q(1)</code>
 * objects, as shown below in the constructor of this class.
 * 
 * @code
 *   template <int dim>
 *   class ElasticProblem
 *   {
 *   public:
 *     ElasticProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     FESystem<dim> fe;
 * 
 *     AffineConstraints<double> constraints;
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
 * <a name="Righthandsidevalues"></a> 
 * <h3>Right hand side values</h3>
 * 

 * 
 * Before going over to the implementation of the main class, we declare and
 * define the function which describes the right hand side. This time, the
 * right hand side is vector-valued, as is the solution, so we will describe
 * the changes required for this in some more detail.
 *   

 * 
 * To prevent cases where the return vector has not previously been set to
 * the right size we test for this case and otherwise throw an exception at
 * the beginning of the function. Note that enforcing that output arguments
 * already have the correct size is a convention in deal.II, and enforced
 * almost everywhere. The reason is that we would otherwise have to check at
 * the beginning of the function and possibly change the size of the output
 * vector. This is expensive, and would almost always be unnecessary (the
 * first call to the function would set the vector to the right size, and
 * subsequent calls would only have to do redundant checks). In addition,
 * checking and possibly resizing the vector is an operation that can not be
 * removed if we can't rely on the assumption that the vector already has
 * the correct size; this is in contract to the Assert call that is
 * completely removed if the program is compiled in optimized mode.
 *   

 * 
 * Likewise, if by some accident someone tried to compile and run the
 * program in only one space dimension (in which the elastic equations do
 * not make much sense since they reduce to the ordinary Laplace equation),
 * we terminate the program in the second assertion. The program will work
 * just fine in 3d, however.
 * 
 * @code
 *   template <int dim>
 *   void right_hand_side(const std::vector<Point<dim>> &points,
 *                        std::vector<Tensor<1, dim>> &  values)
 *   {
 *     Assert(values.size() == points.size(),
 *            ExcDimensionMismatch(values.size(), points.size()));
 *     Assert(dim >= 2, ExcNotImplemented());
 * 
 * @endcode
 * 
 * The rest of the function implements computing force values. We will use
 * a constant (unit) force in x-direction located in two little circles
 * (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
 * an area around the origin; in 3d, the z-component of these centers is
 * zero as well.
 *     

 * 
 * For this, let us first define two objects that denote the centers of
 * these areas. Note that upon construction of the Point objects, all
 * components are set to zero.
 * 
 * @code
 *     Point<dim> point_1, point_2;
 *     point_1(0) = 0.5;
 *     point_2(0) = -0.5;
 * 
 *     for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
 *       {
 * @endcode
 * 
 * If <code>points[point_n]</code> is in a circle (sphere) of radius
 * 0.2 around one of these points, then set the force in x-direction
 * to one, otherwise to zero:
 * 
 * @code
 *         if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
 *             ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
 *           values[point_n][0] = 1.0;
 *         else
 *           values[point_n][0] = 0.0;
 * 
 * @endcode
 * 
 * Likewise, if <code>points[point_n]</code> is in the vicinity of the
 * origin, then set the y-force to one, otherwise to zero:
 * 
 * @code
 *         if (points[point_n].norm_square() < 0.2 * 0.2)
 *           values[point_n][1] = 1.0;
 *         else
 *           values[point_n][1] = 0.0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclassimplementation"></a> 
 * <h3>The <code>ElasticProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="ElasticProblemElasticProblemconstructor"></a> 
 * <h4>ElasticProblem::ElasticProblem constructor</h4>
 * 

 * 
 * Following is the constructor of the main class. As said before, we would
 * like to construct a vector-valued finite element that is composed of
 * several scalar finite elements (i.e., we want to build the vector-valued
 * element so that each of its vector components consists of the shape
 * functions of a scalar element). Of course, the number of scalar finite
 * elements we would like to stack together equals the number of components
 * the solution function has, which is <code>dim</code> since we consider
 * displacement in each space direction. The FESystem class can handle this:
 * we pass it the finite element of which we would like to compose the
 * system of, and how often it shall be repeated:
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   ElasticProblem<dim>::ElasticProblem()
 *     : dof_handler(triangulation)
 *     , fe(FE_Q<dim>(1), dim)
 *   {}
 * @endcode
 * 
 * In fact, the FESystem class has several more constructors which can
 * perform more complex operations than just stacking together several
 * scalar finite elements of the same type into one; we will get to know
 * these possibilities in later examples.
 * 

 * 
 * 

 * 
 * 
 * <a name="ElasticProblemsetup_system"></a> 
 * <h4>ElasticProblem::setup_system</h4>
 * 

 * 
 * Setting up the system of equations is identical to the function used in
 * the step-6 example. The DoFHandler class and all other classes used here
 * are fully aware that the finite element we want to use is vector-valued,
 * and take care of the vector-valuedness of the finite element
 * themselves. (In fact, they do not, but this does not need to bother you:
 * since they only need to know how many degrees of freedom there are per
 * vertex, line and cell, and they do not ask what they represent,
 * i.e. whether the finite element under consideration is vector-valued or
 * whether it is, for example, a scalar Hermite element with several degrees
 * of freedom on each vertex).
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(dim),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ false);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemassemble_system"></a> 
 * <h4>ElasticProblem::assemble_system</h4>
 * 

 * 
 * The big changes in this program are in the creation of matrix and right
 * hand side, since they are problem-dependent. We will go through that
 * process step-by-step, since it is a bit more complicated than in previous
 * examples.
 *   

 * 
 * The first parts of this function are the same as before, however: setting
 * up a suitable quadrature formula, initializing an FEValues object for the
 * (vector-valued) finite element we use as well as the quadrature object,
 * and declaring a number of auxiliary arrays. In addition, we declare the
 * ever same two abbreviations: <code>n_q_points</code> and
 * <code>dofs_per_cell</code>. The number of degrees of freedom per cell we
 * now obviously ask from the composed finite element rather than from the
 * underlying scalar Q1 element. Here, it is <code>dim</code> times the
 * number of degrees of freedom per cell of the Q1 element, though this is
 * not explicit knowledge we need to care about:
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim> quadrature_formula(fe.degree + 1);
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
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * As was shown in previous examples as well, we need a place where to
 * store the values of the coefficients at all the quadrature points on a
 * cell. In the present situation, we have two coefficients, lambda and
 * mu.
 * 
 * @code
 *     std::vector<double> lambda_values(n_q_points);
 *     std::vector<double> mu_values(n_q_points);
 * 
 * @endcode
 * 
 * Well, we could as well have omitted the above two arrays since we will
 * use constant coefficients for both lambda and mu, which can be declared
 * like this. They both represent functions always returning the constant
 * value 1.0. Although we could omit the respective factors in the
 * assemblage of the matrix, we use them here for purpose of
 * demonstration.
 * 
 * @code
 *     Functions::ConstantFunction<dim> lambda(1.), mu(1.);
 * 
 * @endcode
 * 
 * Like the two constant functions above, we will call the function
 * right_hand_side just once per cell to make things simpler.
 * 
 * @code
 *     std::vector<Tensor<1, dim>> rhs_values(n_q_points);
 * 
 * @endcode
 * 
 * Now we can begin with the loop over all cells:
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Next we get the values of the coefficients at the quadrature
 * points. Likewise for the right hand side:
 * 
 * @code
 *         lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
 *         mu.value_list(fe_values.get_quadrature_points(), mu_values);
 *         right_hand_side(fe_values.get_quadrature_points(), rhs_values);
 * 
 * @endcode
 * 
 * Then assemble the entries of the local stiffness matrix and right
 * hand side vector. This follows almost one-to-one the pattern
 * described in the introduction of this example.  One of the few
 * comments in place is that we can compute the number
 * <code>comp(i)</code>, i.e. the index of the only nonzero vector
 * component of shape function <code>i</code> using the
 * <code>fe.system_to_component_index(i).first</code> function call
 * below.
 *         

 * 
 * (By accessing the <code>first</code> variable of the return value
 * of the <code>system_to_component_index</code> function, you might
 * already have guessed that there is more in it. In fact, the
 * function returns a <code>std::pair@<unsigned int, unsigned
 * int@></code>, of which the first element is <code>comp(i)</code>
 * and the second is the value <code>base(i)</code> also noted in the
 * introduction, i.e.  the index of this shape function within all the
 * shape functions that are nonzero in this component,
 * i.e. <code>base(i)</code> in the diction of the introduction. This
 * is not a number that we are usually interested in, however.)
 *         

 * 
 * With this knowledge, we can assemble the local matrix
 * contributions:
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             const unsigned int component_i =
 *               fe.system_to_component_index(i).first;
 * 
 *             for (const unsigned int j : fe_values.dof_indices())
 *               {
 *                 const unsigned int component_j =
 *                   fe.system_to_component_index(j).first;
 * 
 *                 for (const unsigned int q_point :
 *                      fe_values.quadrature_point_indices())
 *                   {
 *                     cell_matrix(i, j) +=
 * @endcode
 * 
 * The first term is $\lambda \partial_i u_i, \partial_j
 * v_j) + (\mu \partial_i u_j, \partial_j v_i)$. Note
 * that <code>shape_grad(i,q_point)</code> returns the
 * gradient of the only nonzero component of the i-th
 * shape function at quadrature point q_point. The
 * component <code>comp(i)</code> of the gradient, which
 * is the derivative of this only nonzero vector
 * component of the i-th shape function with respect to
 * the comp(i)th coordinate is accessed by the appended
 * brackets.
 * 
 * @code
 *                       (                                                  
 *                         (fe_values.shape_grad(i, q_point)[component_i] * 
 *                          fe_values.shape_grad(j, q_point)[component_j] * 
 *                          lambda_values[q_point])                         
 *                         +                                                
 *                         (fe_values.shape_grad(i, q_point)[component_j] * 
 *                          fe_values.shape_grad(j, q_point)[component_i] * 
 *                          mu_values[q_point])                             
 *                         +                                                
 * @endcode
 * 
 * The second term is $(\mu \nabla u_i, \nabla
 * v_j)$. We need not access a specific component of
 * the gradient, since we only have to compute the
 * scalar product of the two gradients, of which an
 * overloaded version of <tt>operator*</tt> takes
 * care, as in previous examples.
 *                         

 * 
 * Note that by using the <tt>?:</tt> operator, we only
 * do this if <tt>component_i</tt> equals
 * <tt>component_j</tt>, otherwise a zero is added
 * (which will be optimized away by the compiler).
 * 
 * @code
 *                         ((component_i == component_j) ?        
 *                            (fe_values.shape_grad(i, q_point) * 
 *                             fe_values.shape_grad(j, q_point) * 
 *                             mu_values[q_point]) :              
 *                            0)                                  
 *                         ) *                                    
 *                       fe_values.JxW(q_point);                  
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * Assembling the right hand side is also just as discussed in the
 * introduction:
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             const unsigned int component_i =
 *               fe.system_to_component_index(i).first;
 * 
 *             for (const unsigned int q_point :
 *                  fe_values.quadrature_point_indices())
 *               cell_rhs(i) += fe_values.shape_value(i, q_point) *
 *                              rhs_values[q_point][component_i] *
 *                              fe_values.JxW(q_point);
 *           }
 * 
 * @endcode
 * 
 * The transfer from local degrees of freedom into the global matrix
 * and right hand side vector does not depend on the equation under
 * consideration, and is thus the same as in all previous
 * examples.
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsolve"></a> 
 * <h4>ElasticProblem::solve</h4>
 * 

 * 
 * The solver does not care about where the system of equations comes, as
 * long as it stays positive definite and symmetric (which are the
 * requirements for the use of the CG solver), which the system indeed
 * is. Therefore, we need not change anything.
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrefine_grid"></a> 
 * <h4>ElasticProblem::refine_grid</h4>
 * 

 * 
 * The function that does the refinement of the grid is the same as in the
 * step-6 example. The quadrature formula is adapted to the linear elements
 * again. Note that the error estimator by default adds up the estimated
 * obtained from all components of the finite element solution, i.e., it
 * uses the displacement in all directions with the same weight. If we would
 * like the grid to be adapted to the x-displacement only, we could pass the
 * function an additional parameter which tells it to do so and do not
 * consider the displacements in all other directions for the error
 * indicators. However, for the current problem, it seems appropriate to
 * consider all displacement components with equal weight.
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::refine_grid()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                        QGauss<dim - 1>(fe.degree + 1),
 *                                        {},
 *                                        solution,
 *                                        estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemoutput_results"></a> 
 * <h4>ElasticProblem::output_results</h4>
 * 

 * 
 * The output happens mostly as has been shown in previous examples
 * already. The only difference is that the solution function is vector
 * valued. The DataOut class takes care of this automatically, but we have
 * to give each component of the solution vector a different name.
 *   

 * 
 * To do this, the DataOut::add_vector() function wants a vector of
 * strings. Since the number of components is the same as the number
 * of dimensions we are working in, we use the <code>switch</code>
 * statement below.
 *   

 * 
 * We note that some graphics programs have restriction on what
 * characters are allowed in the names of variables. deal.II therefore
 * supports only the minimal subset of these characters that is supported
 * by all programs. Basically, these are letters, numbers, underscores,
 * and some other characters, but in particular no whitespace and
 * minus/hyphen. The library will throw an exception otherwise, at least
 * if in debug mode.
 *   

 * 
 * After listing the 1d, 2d, and 3d case, it is good style to let the
 * program die if we run upon a case which we did not consider. Remember
 * that the Assert macro generates an exception if the condition in the
 * first parameter is not satisfied. Of course, the condition
 * <code>false</code> can never be satisfied, so the program will always
 * abort whenever it gets to the default statement:
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 * 
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 1:
 *           solution_names.emplace_back("displacement");
 *           break;
 *         case 2:
 *           solution_names.emplace_back("x_displacement");
 *           solution_names.emplace_back("y_displacement");
 *           break;
 *         case 3:
 *           solution_names.emplace_back("x_displacement");
 *           solution_names.emplace_back("y_displacement");
 *           solution_names.emplace_back("z_displacement");
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 * @endcode
 * 
 * After setting up the names for the different components of the
 * solution vector, we can add the solution vector to the list of
 * data vectors scheduled for output. Note that the following
 * function takes a vector of strings as second argument, whereas
 * the one which we have used in all previous examples accepted a
 * string there. (In fact, the function we had used before would
 * convert the single string into a vector with only one element
 * and forwards that to the other function.)
 * 
 * @code
 *     data_out.add_data_vector(solution, solution_names);
 *     data_out.build_patches();
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrun"></a> 
 * <h4>ElasticProblem::run</h4>
 * 

 * 
 * The <code>run</code> function does the same things as in step-6, for
 * example. This time, we use the square [-1,1]^d as domain, and we refine
 * it globally four times before starting the first iteration.
 *   

 * 
 * The reason for refining is a bit accidental: we use the QGauss
 * quadrature formula with two points in each direction for integration of the
 * right hand side; that means that there are four quadrature points on each
 * cell (in 2D). If we only refine the initial grid once globally, then there
 * will be only four quadrature points in each direction on the
 * domain. However, the right hand side function was chosen to be rather
 * localized and in that case, by pure chance, it happens that all quadrature
 * points lie at points where the right hand side function is zero (in
 * mathematical terms, the quadrature points happen to be at points outside
 * the <i>support</i> of the right hand side function). The right hand side
 * vector computed with quadrature will then contain only zeroes (even though
 * it would of course be nonzero if we had computed the right hand side vector
 * exactly using the integral) and the solution of the system of
 * equations is the zero vector, i.e., a finite element function that is zero
 * everywhere. In a sense, we
 * should not be surprised that this is happening since we have chosen
 * an initial grid that is totally unsuitable for the problem at hand.
 *   

 * 
 * The unfortunate thing is that if the discrete solution is constant, then
 * the error indicators computed by the KellyErrorEstimator class are zero
 * for each cell as well, and the call to
 * Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
 * for refinement (why should it if the indicated error is zero for each
 * cell?). The grid in the next iteration will therefore consist of four
 * cells only as well, and the same problem occurs again.
 *   

 * 
 * The conclusion needs to be: while of course we will not choose the
 * initial grid to be well-suited for the accurate solution of the problem,
 * we must at least choose it such that it has the chance to capture the
 * important features of the solution. In this case, it needs to be able to
 * see the right hand side. Thus, we refine globally four times. (Any larger
 * number of global refinement steps would of course also work.)
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 8; ++cycle)
 *       {
 *         std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, -1, 1);
 *             triangulation.refine_global(4);
 *           }
 *         else
 *           refine_grid();
 * 
 *         std::cout << "   Number of active cells:       "
 *                   << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *                   << std::endl;
 * 
 *         assemble_system();
 *         solve();
 *         output_results(cycle);
 *       }
 *   }
 * } // namespace Step8
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * After closing the <code>Step8</code> namespace in the last line above, the
 * following is the main function of the program and is again exactly like in
 * step-6 (apart from the changed class names, of course).
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step8::ElasticProblem<2> elastic_problem_2d;
 *       elastic_problem_2d.run();
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



There is not much to be said about the results of this program, other than
that they look nice. All images were made using VisIt from the
output files that the program wrote to disk. The first two pictures show
the $x$- and $y$-displacements as scalar components:

<table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.x.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.y.png" alt="">
</td>
</tr>
</table>


You can clearly see the sources of $x$-displacement around $x=0.5$ and
$x=-0.5$, and of $y$-displacement at the origin.

What one frequently would like to do is to show the displacement as a vector
field, i.e., vectors that for each point illustrate the direction and magnitude
of displacement. Unfortunately, that's a bit more involved. To understand why
this is so, remember that we have just defined our finite element as a
collection of two  components (in <code>dim=2</code> dimensions). Nowhere have
we said that this is not just a pressure and a concentration (two scalar
quantities) but that the two components actually are the parts of a
vector-valued quantity, namely the displacement. Absent this knowledge, the
DataOut class assumes that all individual variables we print are separate
scalars, and VisIt and Paraview then faithfully assume that this is indeed what it is. In
other words, once we have written the data as scalars, there is nothing in
these programs that allows us to paste these two scalar fields back together as a
vector field. Where we would have to attack this problem is at the root,
namely in <code>ElasticProblem::output_results()</code>. We won't do so here but
instead refer the reader to the step-22 program where we show how to do this
for a more general situation. That said, we couldn't help generating the data
anyway that would show how this would look if implemented as discussed in
step-22. The vector field then looks like this (VisIt and Paraview
randomly select a few
hundred vertices from which to draw the vectors; drawing them from each
individual vertex would make the picture unreadable):

<img src="https://www.dealii.org/images/steps/developer/step-8.vectors.png" alt="">


We note that one may have intuitively expected the
solution to be symmetric about the $x$- and $y$-axes since the $x$- and
$y$-forces are symmetric with respect to these axes. However, the force
considered as a vector is not symmetric and consequently neither is
the solution.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-8.cc"
*/
