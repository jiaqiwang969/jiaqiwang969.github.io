/**
@page step_20 The step-20 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Theequations">The equations</a>
        <li><a href="#Formulationweakformanddiscreteproblem">Formulation, weak form, and discrete problem</a>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
      <ul>
        <li><a href="#SolvingusingtheSchurcomplement">Solving using the Schur complement</a>
        <li><a href="#TheLinearOperatorframeworkindealII">The LinearOperator framework in deal.II</a>
        <li><a href="#ApreconditionerfortheSchurcomplement">A preconditioner for the Schur complement</a>
      </ul>
        <li><a href="#Definitionofthetestcase">Definition of the test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMixedLaplaceProblemcodeclasstemplate">The <code>MixedLaplaceProblem</code> class template</a>
        <li><a href="#Righthandsideboundaryvaluesandexactsolution">Right hand side, boundary values, and exact solution</a>
        <li><a href="#Theinversepermeabilitytensor">The inverse permeability tensor</a>
        <li><a href="#MixedLaplaceProblemclassimplementation">MixedLaplaceProblem class implementation</a>
      <ul>
        <li><a href="#MixedLaplaceProblemMixedLaplaceProblem">MixedLaplaceProblem::MixedLaplaceProblem</a>
        <li><a href="#MixedLaplaceProblemmake_grid_and_dofs">MixedLaplaceProblem::make_grid_and_dofs</a>
        <li><a href="#MixedLaplaceProblemassemble_system">MixedLaplaceProblem::assemble_system</a>
      </ul>
        <li><a href="#Implementationoflinearsolversandpreconditioners">Implementation of linear solvers and preconditioners</a>
      <ul>
        <li><a href="#MixedLaplacesolve">MixedLaplace::solve</a>
      </ul>
        <li><a href="#MixedLaplaceProblemclassimplementationcontinued">MixedLaplaceProblem class implementation (continued)</a>
      <ul>
        <li><a href="#MixedLaplacecompute_errors">MixedLaplace::compute_errors</a>
        <li><a href="#MixedLaplaceoutput_results">MixedLaplace::output_results</a>
        <li><a href="#MixedLaplacerun">MixedLaplace::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Outputoftheprogramandgraphicalvisualization">Output of the program and graphical visualization</a>
        <li><a href="#Convergence">Convergence</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Morerealisticpermeabilityfields">More realistic permeability fields</a>
        <li><a href="#Betterlinearsolvers">Better linear solvers</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{19,20,21}

This program is devoted to two aspects: the use of mixed finite elements -- in
particular Raviart-Thomas elements -- and using block matrices to define
solvers, preconditioners, and nested versions of those that use the
substructure of the system matrix. The equation we are going to solve is again
the Poisson equation, though with a matrix-valued coefficient:
@f{eqnarray*}
  -\nabla \cdot K({\mathbf x}) \nabla p &=& f \qquad {\textrm{in}\ } \Omega, \\
  p &=& g \qquad {\textrm{on}\ }\partial\Omega.
@f}
$K({\mathbf x})$ is assumed to be uniformly positive definite, i.e., there is
$\alpha>0$ such that the eigenvalues $\lambda_i({\mathbf x})$ of $K(x)$ satisfy
$\lambda_i({\mathbf x})\ge \alpha$. The use of the symbol $p$ instead of the usual
$u$ for the solution variable will become clear in the next section.

After discussing the equation and the formulation we are going to use to solve
it, this introduction will cover the use of block matrices and vectors, the
definition of solvers and preconditioners, and finally the actual test case we
are going to solve.

We are going to extend this tutorial program in step-21 to
solve not only the mixed Laplace equation, but add another equation that
describes the transport of a mixture of two fluids.

The equations covered here fall into the class of vector-valued problems. A
toplevel overview of this topic can be found in the @ref vector_valued module.


<a name="Theequations"></a><h3>The equations</h3>


In the form above, the Poisson equation (i.e., the Laplace equation with a nonzero
right hand side) is generally considered a good model equation
for fluid flow in porous media. Of course, one typically models fluid flow through
the <a href="https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations">Navier-Stokes
equations</a> or, if fluid velocities are slow or the viscosity is large, the
<a href="https://en.wikipedia.org/wiki/Stokes_flow">Stokes
equations</a>
(which we cover in step-22).
In the first of these two models, the forces that act are inertia and
viscous friction, whereas in the second it is only viscous friction -- i.e.,
forces that one fluid particle exerts on a nearby one. This is appropriate
if you have free flow in a large domain, say a pipe, a river, or in the air.
On the other hand, if the fluid is confined in pores, then friction forces
exerted by the pore walls on the fluid become more and more important and
internal viscous friction becomes less and less important. Modeling this
then first leads to the
<a href="https://en.wikipedia.org/wiki/Darcy%27s_law#Brinkman_form_of_Darcy's_law">Brinkman
model</a> if both effects are important, and in the limit of very small pores
to the <a href="https://en.wikipedia.org/wiki/Darcy%27s_law">Darcy equations</a>.
The latter is just a different name for the Poisson or Laplace equation,
connotating it with the area to which one wants to apply it: slow flow
in a porous medium. In essence it says that the velocity is proportional
to the negative pressure gradient that drives the fluid through the
porous medium.

The Darcy equation models this pressure that drives the flow. (Because the
solution variable is a pressure, we here use the name $p$ instead of the
name $u$ more commonly used for the solution of partial differential equations.)
Typical applications of this view of the Laplace equation are then modeling
groundwater flow, or the flow of hydrocarbons in oil reservoirs. In these
applications, $K$ is the permeability tensor, i.e., a measure for how much
resistance the soil or rock matrix asserts on the fluid flow.

In the applications named above, a desirable feature for a numerical
scheme is that it should be locally conservative, i.e., that whatever
flows into a cell also flows out of it (or the difference is equal to
the integral over the source terms over each cell, if the sources are
nonzero). However, as it turns out, the usual discretizations of the
Laplace equation (such as those used in step-3, step-4, or step-6) do
not satisfy this property. But, one can achieve this by choosing a
different formulation of the problem and a particular combination of
finite element spaces.


<a name="Formulationweakformanddiscreteproblem"></a><h3>Formulation, weak form, and discrete problem</h3>


To this end, one first introduces a second variable, called the velocity,
${\mathbf u}=-K\nabla p$. By its definition, the velocity is a vector in the
negative
direction of the pressure gradient, multiplied by the permeability tensor. If
the permeability tensor is proportional to the unit matrix, this equation is
easy to understand and intuitive: the higher the permeability, the higher the
velocity; and the velocity is proportional to the gradient of the pressure, going from
areas of high pressure to areas of low pressure (thus the negative sign).

With this second variable, one then finds an alternative version of the
Laplace equation, called the <i>mixed formulation</i>:
@f{eqnarray*}
  K^{-1} {\mathbf u} + \nabla p &=& 0 \qquad {\textrm{in}\ } \Omega, \\
  -{\textrm{div}}\ {\mathbf u} &=& -f \qquad {\textrm{in}\ }\Omega, \\
  p &=& g \qquad {\textrm{on}\ } \partial\Omega.
@f}
Here, we have multiplied the equation defining the velocity ${\mathbf
u}$ by $K^{-1}$ because this makes the set of equations symmetric: one
of the equations has the gradient, the second the negative divergence,
and these two are of course adjoints of each other, resulting in a
symmetric bilinear form and a consequently symmetric system matrix
under the common assumption that $K$ is a symmetric tensor.

The weak formulation of this problem is found by multiplying the two
equations with test functions and integrating some terms by parts:
@f{eqnarray*}
  A(\{{\mathbf u},p\},\{{\mathbf v},q\}) = F(\{{\mathbf v},q\}),
@f}
where
@f{eqnarray*}
  A(\{{\mathbf u},p\},\{{\mathbf v},q\})
  &=&
  ({\mathbf v}, K^{-1}{\mathbf u})_\Omega - ({\textrm{div}}\ {\mathbf v}, p)_\Omega
  - (q,{\textrm{div}}\ {\mathbf u})_\Omega
  \\
  F(\{{\mathbf v},q\}) &=& -(g,{\mathbf v}\cdot {\mathbf n})_{\partial\Omega} - (f,q)_\Omega.
@f}
Here, ${\mathbf n}$ is the outward normal vector at the boundary. Note how in this
formulation, Dirichlet boundary values of the original problem are
incorporated in the weak form.

To be well-posed, we have to look for solutions and test functions in the
space $H({\textrm{div}})=\{{\mathbf w}\in L^2(\Omega)^d:\ {\textrm{div}}\ {\mathbf w}\in L^2\}$
for $\mathbf u$,$\mathbf v$, and $L^2$ for $p,q$. It is a well-known fact stated in
almost every book on finite element theory that if one chooses discrete finite
element spaces for the approximation of ${\mathbf u},p$ inappropriately, then the
resulting discrete problem is instable and the discrete solution
will not converge to the exact solution. (Some details on the problem
considered here -- which falls in the class of "saddle-point problems"
-- can be found on the Wikipedia page on the <a
href="https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition">Ladyzhenskaya-Babuska-Brezzi
(LBB) condition</a>.)

To overcome this, a number of different finite element pairs for ${\mathbf u},p$
have been developed that lead to a stable discrete problem. One such pair is
to use the Raviart-Thomas spaces $RT(k)$ for the velocity ${\mathbf u}$ and
discontinuous elements of class $DQ(k)$ for the pressure $p$. For details
about these spaces, we refer in particular to the book on mixed finite element
methods by Brezzi and Fortin, but many other books on the theory of finite
elements, for example the classic book by Brenner and Scott, also state the
relevant results. In any case, with appropriate choices of function
spaces, the discrete formulation reads as follows: Find ${\mathbf
u}_h,p_h$ so that
@f{eqnarray*}
  A(\{{\mathbf u}_h,p_h\},\{{\mathbf v}_h,q_h\}) = F(\{{\mathbf v}_h,q_h\})
  \qquad\qquad \forall {\mathbf v}_h,q_h.
@f}


Before continuing, let us briefly pause and show that the choice of
function spaces above provides us with the desired local conservation
property. In particular, because the pressure space consists of
discontinuous piecewise polynomials, we can choose the test function
$q$ as the function that is equal to one on any given cell $K$ and
zero everywhere else. If we also choose ${\mathbf v}=0$ everywhere
(remember that the weak form above has to hold for <i>all</i> discrete
test functions $q,v$), then putting these choices of test functions
into the weak formulation above implies in particular that
@f{eqnarray*}
  - (1,{\textrm{div}}\ {\mathbf u}_h)_K
  =
  -(1,f)_K,
@f}
which we can of course write in more explicit form as
@f{eqnarray*}
  \int_K {\textrm{div}}\ {\mathbf u}_h
 =
  \int_K f.
@f}
Applying the divergence theorem results in the fact that ${\mathbf
u}_h$ has to satisfy, for every choice of cell $K$, the relationship
@f{eqnarray*}
  \int_{\partial K} {\mathbf u}_h\cdot{\mathbf n}
  =
  \int_K f.
@f}
If you now recall that ${\mathbf u}$ was the velocity, then the
integral on the left is exactly the (discrete) flux across the
boundary of the cell $K$. The statement is then that the flux must be equal
to the integral over the sources within $K$. In particular, if there
are no sources (i.e., $f=0$ in $K$), then the statement is that
<i>total</i> flux is zero, i.e., whatever flows into a cell must flow out
of it through some other part of the cell boundary. This is what we call
<i>local conservation</i> because it holds for every cell.

On the other hand, the usual continuous $Q_k$ elements would not result in
this kind of property when used for the pressure (as, for example, we
do in step-43) because one can not choose a discrete test function
$q_h$ that is one on a cell $K$ and zero everywhere else: It would be
discontinuous and consequently not in the finite element
space. (Strictly speaking, all we can say is that the proof above
would not work for continuous elements. Whether these elements might
still result in local conservation is a different question as one
could think that a different kind of proof might still work; in
reality, however, the property really does not hold.)



<a name="Assemblingthelinearsystem"></a><h3>Assembling the linear system</h3>


The deal.II library (of course) implements Raviart-Thomas elements $RT(k)$ of
arbitrary order $k$, as well as discontinuous elements $DG(k)$. If we forget
about their particular properties for a second, we then have to solve a
discrete problem
@f{eqnarray*}
  A(x_h,w_h) = F(w_h),
@f}
with the bilinear form and right hand side as stated above, and $x_h=\{{\mathbf u}_h,p_h\}$, $w_h=\{{\mathbf v}_h,q_h\}$. Both $x_h$ and $w_h$ are from the space
$X_h=RT(k)\times DQ(k)$, where $RT(k)$ is itself a space of $dim$-dimensional
functions to accommodate for the fact that the flow velocity is vector-valued.
The necessary question then is: how do we do this in a program?

Vector-valued elements have already been discussed in previous tutorial
programs, the first time and in detail in step-8. The main difference there
was that the vector-valued space $V_h$ is uniform in all its components: the
$dim$ components of the displacement vector are all equal and from the same
function space. What we could therefore do was to build $V_h$ as the outer
product of the $dim$ times the usual $Q(1)$ finite element space, and by this
make sure that all our shape functions have only a single non-zero vector
component. Instead of dealing with vector-valued shape functions, all we did
in step-8 was therefore to look at the (scalar) only non-zero component and
use the <code>fe.system_to_component_index(i).first</code> call to figure out
which component this actually is.

This doesn't work with Raviart-Thomas elements: following from their
construction to satisfy certain regularity properties of the space
$H({\textrm{div}})$, the shape functions of $RT(k)$ are usually nonzero in all
their vector components at once. For this reason, were
<code>fe.system_to_component_index(i).first</code> applied to determine the only
nonzero component of shape function $i$, an exception would be generated. What
we really need to do is to get at <em>all</em> vector components of a shape
function. In deal.II diction, we call such finite elements
<em>non-primitive</em>, whereas finite elements that are either scalar or for
which every vector-valued shape function is nonzero only in a single vector
component are called <em>primitive</em>.

So what do we have to do for non-primitive elements? To figure this out, let
us go back in the tutorial programs, almost to the very beginnings. There, we
learned that we use the <code>FEValues</code> class to determine the values and
gradients of shape functions at quadrature points. For example, we would call
<code>fe_values.shape_value(i,q_point)</code> to obtain the value of the
<code>i</code>th shape function on the quadrature point with number
<code>q_point</code>. Later, in step-8 and other tutorial programs, we learned
that this function call also works for vector-valued shape functions (of
primitive finite elements), and that it returned the value of the only
non-zero component of shape function <code>i</code> at quadrature point
<code>q_point</code>.

For non-primitive shape functions, this is clearly not going to work: there is
no single non-zero vector component of shape function <code>i</code>, and the call
to <code>fe_values.shape_value(i,q_point)</code> would consequently not make
much sense. However, deal.II offers a second function call,
<code>fe_values.shape_value_component(i,q_point,comp)</code> that returns the
value of the <code>comp</code>th vector component of shape function  <code>i</code> at
quadrature point <code>q_point</code>, where <code>comp</code> is an index between
zero and the number of vector components of the present finite element; for
example, the element we will use to describe velocities and pressures is going
to have $dim+1$ components. It is worth noting that this function call can
also be used for primitive shape functions: it will simply return zero for all
components except one; for non-primitive shape functions, it will in general
return a non-zero value for more than just one component.

We could now attempt to rewrite the bilinear form above in terms of vector
components. For example, in 2d, the first term could be rewritten like this
(note that $u_0=x_0, u_1=x_1, p=x_2$):
@f{eqnarray*}
  ({\mathbf u}_h^i, K^{-1}{\mathbf u}_h^j)
  =
  &\left((x_h^i)_0, K^{-1}_{00} (x_h^j)_0\right) +
   \left((x_h^i)_0, K^{-1}_{01} (x_h^j)_1\right) + \\
  &\left((x_h^i)_1, K^{-1}_{10} (x_h^j)_0\right) +
   \left((x_h^i)_1, K^{-1}_{11} (x_h^j)_1\right).
@f}
If we implemented this, we would get code like this:

@code
  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        local_matrix(i,j) += (k_inverse_values[q][0][0] *
                              fe_values.shape_value_component(i,q,0) *
                              fe_values.shape_value_component(j,q,0)
                              +
                              k_inverse_values[q][0][1] *
                              fe_values.shape_value_component(i,q,0) *
                              fe_values.shape_value_component(j,q,1)
                              +
                              k_inverse_values[q][1][0] *
                              fe_values.shape_value_component(i,q,1) *
                              fe_values.shape_value_component(j,q,0)
                              +
                              k_inverse_values[q][1][1] *
                              fe_values.shape_value_component(i,q,1) *
                              fe_values.shape_value_component(j,q,1)
                             ) *
                             fe_values.JxW(q);
@endcode

This is, at best, tedious, error prone, and not dimension independent. There
are obvious ways to make things dimension independent, but in the end, the
code is simply not pretty. What would be much nicer is if we could simply
extract the ${\mathbf u}$ and $p$ components of a shape function $x_h^i$. In the
program we do that in the following way:

@code
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  ...

  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        local_matrix(i,j) += (fe_values[velocities].value (i, q) *
                              k_inverse_values[q] *
                              fe_values[velocities].value (j, q)
                              -
                              fe_values[velocities].divergence (i, q) *
                              fe_values[pressure].value (j, q)
                              -
                              fe_values[pressure].value (i, q) *
                              fe_values[velocities].divergence (j, q)) *
                              fe_values.JxW(q);
@endcode

This is, in fact, not only the first term of the bilinear form, but the
whole thing (sans boundary contributions).

What this piece of code does is, given an <code>fe_values</code> object, to extract
the values of the first $dim$ components of shape function <code>i</code> at
quadrature points <code>q</code>, that is the velocity components of that shape
function. Put differently, if we write shape functions $x_h^i$ as the tuple
$\{{\mathbf u}_h^i,p_h^i\}$, then the function returns the velocity part of this
tuple. Note that the velocity is of course a <code>dim</code>-dimensional tensor, and
that the function returns a corresponding object. Similarly, where we
subscript with the pressure extractor, we extract the scalar pressure
component. The whole mechanism is described in more detail in the
@ref vector_valued module.

In practice, it turns out that we can do a bit better if we evaluate the shape
functions, their gradients and divergences only once per outermost loop, and
store the result, as this saves us a few otherwise repeated computations (it is
possible to save even more repeated operations by calculating all relevant
quantities in advance and then only inserting the results in the actual loop,
see step-22 for a realization of that approach).
The final result then looks like this, working in every space dimension:

@code
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      k_inverse.value_list (fe_values.get_quadrature_points(),
                            k_inverse_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim> phi_i_u     = fe_values[velocities].value (i, q);
            const double        div_phi_i_u = fe_values[velocities].divergence (i, q);
            const double        phi_i_p     = fe_values[pressure].value (i, q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u     = fe_values[velocities].value (j, q);
                const double        div_phi_j_u = fe_values[velocities].divergence (j, q);
                const double        phi_j_p     = fe_values[pressure].value (j, q);

                local_matrix(i,j) += (phi_i_u * k_inverse_values[q] * phi_j_u
                                      - div_phi_i_u * phi_j_p
                                      - phi_i_p * div_phi_j_u) *
                                     fe_values.JxW(q);
              }

            local_rhs(i) += -phi_i_p *
                            rhs_values[q] *
                            fe_values.JxW(q);
          }
@endcode

This very closely resembles the form in which we have originally written down
the bilinear form and right hand side.

There is one final term that we have to take care of: the right hand side
contained the term $(g,{\mathbf v}\cdot {\mathbf n})_{\partial\Omega}$, constituting the
weak enforcement of pressure boundary conditions. We have already seen in
step-7 how to deal with face integrals: essentially exactly the same as with
domain integrals, except that we have to use the FEFaceValues class
instead of <code>FEValues</code>. To compute the boundary term we then simply have
to loop over all boundary faces and integrate there. The mechanism works in
the same way as above, i.e. the extractor classes also work on FEFaceValues objects:

@code
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);

              pressure_boundary_values.value_list(
                fe_face_values.get_quadrature_points(), boundary_values);

              for (unsigned int q = 0; q < n_face_q_points; ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  local_rhs(i) += -(fe_face_values[velocities].value(i, q) *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
@endcode

You will find the exact same code as above in the sources for the present
program. We will therefore not comment much on it below.


<a name="Linearsolversandpreconditioners"></a><h3>Linear solvers and preconditioners</h3>


After assembling the linear system we are faced with the task of solving
it. The problem here is that the matrix possesses two undesirable properties:
- It is <a href="https://en.wikipedia.org/wiki/Definiteness_of_a_matrix">indefinite</a>,
  i.e., it has both positive and negative eigenvalues.
  We don't want to prove this property here, but note that this is true
  for all matrices of the form
  $\left(\begin{array}{cc} M & B \\ B^T & 0 \end{array}\right)$
  such as the one here where $M$ is positive definite.
- The matrix has a zero block at the bottom right (there is no term in
  the bilinear form that couples the pressure $p$ with the
  pressure test function $q$).

At least it is symmetric, but the first issue above still means that
the Conjugate Gradient method is not going to work since it is only
applicable to problems in which the matrix is symmetric and positive definite.
We would have to resort to other iterative solvers instead, such as
MinRes, SymmLQ, or GMRES, that can deal with indefinite systems. However, then
the next problem immediately surfaces: Due to the zero block, there are zeros
on the diagonal and none of the usual, "simple" preconditioners (Jacobi, SSOR)
will work as they require division by diagonal elements.

For the matrix sizes we expect to run with this program, the by far simplest
approach would be to just use a direct solver (in particular, the
SparseDirectUMFPACK class that is bundled with deal.II). step-29 goes this
route and shows that solving <i>any</i> linear system can be done in just
3 or 4 lines of code.

But then, this is a tutorial: We teach how to do things. Consequently,
in the following, we will introduce some techniques that can be used in cases
like these. Namely, we will consider the linear system as not consisting of one
large matrix and vectors, but we will want to decompose matrices
into <i>blocks</i> that correspond to the individual operators that appear in
the system. We note that the resulting solver is not optimal -- there are
much better ways to efficiently compute the system, for example those
explained in the results section of step-22 or the one we use in step-43
for a problem similar to the current one. Here, our goal is simply to
introduce new solution techniques and how they can be implemented in
deal.II.


<a name="SolvingusingtheSchurcomplement"></a><h4>Solving using the Schur complement</h4>


In view of the difficulties using standard solvers and preconditioners
mentioned above, let us take another look at the matrix. If we sort our
degrees of freedom so that all velocity come before all pressure variables,
then we can subdivide the linear system $Ax=b$ into the following blocks:
@f{eqnarray*}
  \left(\begin{array}{cc}
    M & B \\ B^T & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),
@f}
where $U,P$ are the values of velocity and pressure degrees of freedom,
respectively, $M$ is the mass matrix on the velocity space, $B^T$ corresponds to
the negative divergence operator, and $B$ is its transpose and corresponds
to the gradient.

By block elimination, we can then re-order this system in the following way
(multiply the first row of the system by $B^TM^{-1}$ and then subtract the
second row from it):
@f{eqnarray*}
  B^TM^{-1}B P &=& B^TM^{-1} F - G, \\
  MU &=& F - BP.
@f}
Here, the matrix $S=B^TM^{-1}B$ (called the
<a href="https://en.wikipedia.org/wiki/Schur_complement">Schur complement</a>
of $A$)
is obviously symmetric and, owing to the positive definiteness of $M$ and the
fact that $B$ has full column rank, $S$ is also positive
definite.

Consequently, if we could compute $S$, we could apply the Conjugate Gradient
method to it. However, computing $S$ is expensive because it requires us
to compute the inverse of the (possibly large) matrix $M$; and $S$ is in fact
also a full matrix because even though $M$ is sparse, its inverse $M^{-1}$
will generally be a dense matrix.
On the other hand, the CG algorithm doesn't require
us to actually have a representation of $S$: It is sufficient to form
matrix-vector products with it. We can do so in steps, using the fact that
matrix products are associative (i.e., we can set parentheses in such a
way that the product is more convenient to compute):
To compute $Sv=(B^TM^{-1}B)v=B^T(M^{-1}(Bv))$, we
<ol>
 <li> compute $w = B v$;
 <li> solve $My = w$ for $y=M^{-1}w$, using the CG method applied to the
  positive definite and symmetric mass matrix $M$;
 <li> compute $z=B^Ty$ to obtain $z=Sv$.
</ol>
Note how we evaluate the expression $B^TM^{-1}Bv$ right to left to
avoid matrix-matrix products; this way, all we have to do is evaluate
matrix-vector products.

In the following, we will then have to come up with ways to represent the
matrix $S$ so that it can be used in a Conjugate Gradient solver,
as well as to define ways in which we can precondition the solution
of the linear system involving $S$, and deal with solving linear systems
with the matrix $M$ (the second step above).

@note The key point in this consideration is to recognize that to implement
an iterative solver such as CG or GMRES, we never actually need the actual
<i>elements</i> of a matrix! All that is required is that we can form
matrix-vector products. The same is true for preconditioners. In deal.II we
encode this requirement by only requiring that matrices and preconditioners
given to solver classes have a <code>vmult()</code> member function that
does the matrix-vector product. How a class chooses to implement this
function is not important to the solver. Consequently, classes can
implement it by, for example, doing a sequence of products and linear
solves as discussed above.


<a name="TheLinearOperatorframeworkindealII"></a><h4>The LinearOperator framework in deal.II</h4>


deal.II includes support for describing such linear operations in a very
general way. This is done with the LinearOperator class that, like
@ref ConceptMatrixType "the MatrixType concept",
defines a minimal interface for <i>applying</i> a linear operation to a
vector:
@code
    std::function<void(Range &, const Domain &)> vmult;
    std::function<void(Range &, const Domain &)> vmult_add;
    std::function<void(Domain &, const Range &)> Tvmult;
    std::function<void(Domain &, const Range &)> Tvmult_add;
@endcode
The key difference between a LinearOperator and an ordinary matrix is
however that a LinearOperator does not allow any further access to the
underlying object. All you can do with a LinearOperator is to apply its
"action" to a vector! We take the opportunity to introduce the
LinearOperator concept at this point because it is a very useful tool that
allows you to construct complex solvers and preconditioners in a very
intuitive manner.

As a first example let us construct a LinearOperator object that represents
$M^{-1}$. This means that whenever the <code>vmult()</code> function of
this operator is called it has to solve a linear system. This requires us
to specify a solver (and corresponding) preconditioner. Assuming that
<code>M</code> is a reference to the upper left block of the system matrix
we can write:
@code
    const auto op_M = linear_operator(M);

    PreconditionJacobi<SparseMatrix<double>> preconditioner_M;
    preconditioner_M.initialize(M);

    ReductionControl reduction_control_M(2000, 1.0e-18, 1.0e-10);
    SolverCG<Vector<double>>    solver_M(reduction_control_M);

    const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
@endcode
Rather than using a SolverControl we use the ReductionControl class here
that stops iterations when either an absolute tolerance is reached (for
which we choose $10^{-18}$) or when the residual is reduced by a certain
factor (here, $10^{-10}$). In contrast the SolverControl class only checks
for absolute tolerances. We have to use ReductionControl in our case to
work around a minor issue: The right hand sides that we  will feed to
<code>op_M_inv</code> are essentially formed by residuals that naturally
decrease vastly in norm as the outer iterations progress. This makes
control by an absolute tolerance very error prone.

We now have a LinearOperator <code>op_M_inv</code> that we can use to
construct more complicated operators such as the Schur complement $S$.
Assuming that <code>B</code> is a reference to the upper right block
constructing a LinearOperator <code>op_S</code> is a matter of two lines:
@code
    const auto op_B = linear_operator(B);
    const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
@endcode
Here, the multiplication of three LinearOperator objects yields a composite
object <code>op_S</code> whose <code>vmult()</code> function first applies
$B$, then $M^{-1}$ (i.e. solving an equation with $M$), and finally $B^T$
to any given input vector. In that sense <code>op_S.vmult()</code> is
similar to the following code:
@code
    B.vmult (tmp1, src); // multiply with the top right block: B
    solver_M(M, tmp2, tmp1, preconditioner_M); // multiply with M^-1
    B.Tvmult (dst, tmp2); // multiply with the bottom left block: B^T
@endcode
(<code>tmp1</code> and <code>tmp2</code> are two temporary vectors). The
key point behind this approach is the fact that we never actually create an
inner product of matrices. Instead, whenever we have to perform a matrix
vector multiplication with <code>op_S</code> we simply run all individual
<code>vmult</code> operations in above sequence.

@note We could have achieved the same goal of creating a "matrix like"
object by implementing a specialized class <code>SchurComplement</code>
that provides a suitable <code>vmult()</code> function. Skipping over some
details this might have looked like the following:
@code
class SchurComplement
{
  public:

  // ...

  void SchurComplement::vmult (Vector<double>       &dst,
                               const Vector<double> &src) const
  {
    B.vmult (tmp1, src);
    solver_M(M, tmp2, tmp1, preconditioner_M);
    B.Tvmult (dst, tmp2);
  }
};
@endcode
Even though both approaches are exactly equivalent, the LinearOperator
class has a big advantage over this manual approach.
It provides so-called
<i><a href="https://en.wikipedia.org/wiki/Syntactic_sugar">syntactic sugar</a></i>:
Mathematically, we think about $S$ as being the composite matrix
$S=B^TM^{-1}B$ and the LinearOperator class allows you to write this out
more or less verbatim,
@code
const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
@endcode
The manual approach on the other hand obscures this fact.

All that is left for us to do now is to form the right hand sides of the
two equations defining $P$ and $U$, and then solve them with the Schur
complement matrix and the mass matrix, respectively. For example the right
hand side of the first equation reads $B^TM^{-1}F-G$. This could be
implemented as follows:
@code
    Vector<double> schur_rhs (P.size());
    Vector<double> tmp (U.size());
    op_M_inv.vmult (tmp, F);
    transpose_operator(op_B).vmult (schur_rhs, tmp);
    schur_rhs -= G;
@endcode
Again, this is a perfectly valid approach, but the fact that deal.II
requires us to manually resize the final and temporary vector, and that
every operation takes up a new line makes this hard to read. This is the
point where a second class in the linear operator framework can will help
us. Similarly in spirit to LinearOperator, a PackagedOperation stores a
"computation":
@code
    std::function<void(Range &)> apply;
    std::function<void(Range &)> apply_add;
@endcode
The class allows
<a href="https://en.wikipedia.org/wiki/Lazy_evaluation">lazy evaluation</a>
of expressions involving vectors and linear operators. This is done by
storing the computational expression and only performing the computation
when either the object is converted to a vector object, or
PackagedOperation::apply() (or PackagedOperation::apply_add()) is invoked
by hand. Assuming that <code>F</code> and <code>G</code> are the two
vectors of the right hand side we can simply write:
@code
    const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G;
@endcode
Here, <code>schur_rhs</code> is a PackagedOperation that <i>records</i> the
computation we specified. It does not create a vector with the actual
result immediately.

With these prerequisites at hand, solving for $P$ and $U$ is a matter of
creating another solver and inverse:
@code
    SolverControl solver_control_S(2000, 1.e-12);
    SolverCG<Vector<double>>    solver_S(solver_control_S);
    PreconditionIdentity preconditioner_S;

    const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S);

    P = op_S_inv * schur_rhs;
    U = op_M_inv * (F - op_B * P);
@endcode

@note The functionality that we developed in this example step by hand is
already readily available in the library. Have a look at
schur_complement(), condense_schur_rhs(), and postprocess_schur_solution().


<a name="ApreconditionerfortheSchurcomplement"></a><h4>A preconditioner for the Schur complement</h4>


One may ask whether it would help if we had a preconditioner for the Schur
complement $S=B^TM^{-1}B$. The general answer, as usual, is: of course. The
problem is only, we don't know anything about this Schur complement matrix. We
do not know its entries, all we know is its action. On the other hand, we have
to realize that our solver is expensive since in each iteration we have to do
one matrix-vector product with the Schur complement, which means that we have
to do invert the mass matrix once in each iteration.

There are different approaches to preconditioning such a matrix. On the one
extreme is to use something that is cheap to apply and therefore has no real
impact on the work done in each iteration. The other extreme is a
preconditioner that is itself very expensive, but in return really brings down
the number of iterations required to solve with $S$.

We will try something along the second approach, as much to improve the
performance of the program as to demonstrate some techniques. To this end, let
us recall that the ideal preconditioner is, of course, $S^{-1}$, but that is
unattainable. However, how about
@f{eqnarray*}
  \tilde S^{-1} = [B^T ({\textrm{diag}\ }M)^{-1}B]^{-1}
@f}
as a preconditioner? That would mean that every time we have to do one
preconditioning step, we actually have to solve with $\tilde S$. At first,
this looks almost as expensive as solving with $S$ right away. However, note
that in the inner iteration, we do not have to calculate $M^{-1}$, but only
the inverse of its diagonal, which is cheap.

Thankfully, the LinearOperator framework makes this very easy to write out.
We already used a Jacobi preconditioner (<code>preconditioner_M</code>) for
the $M$ matrix earlier. So all that is left to do is to write out how the
approximate Schur complement should look like:
@code
    const auto op_aS =
      transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B;
@endcode
Note how this operator differs in simply doing one Jacobi sweep
(i.e. multiplying with the inverses of the diagonal) instead of multiplying
with the full $M^{-1}$. (This is how a single Jacobi preconditioner step
with $M$ is defined: it is the multiplication with the inverse of the
diagonal of $M$; in other words, the operation $({\textrm{diag}\ }M)^{-1}x$
on a vector $x$ is exactly what PreconditionJacobi does.)

With all this we almost have the preconditioner completed: it should be the
inverse of the approximate Schur complement. We implement this again by
creating a linear operator with inverse_operator() function. This time
however we would like to choose a relatively modest tolerance for the CG
solver (that inverts <code>op_aS</code>). The reasoning is that
<code>op_aS</code> is only coarse approximation to <code>op_S</code>, so we
actually do not need to invert it exactly. This, however creates a subtle
problem: <code>preconditioner_S</code> will be used in the final outer CG
iteration to create an orthogonal basis. But for this to work, it must be
precisely the same linear operation for every invocation. We ensure this by
using an IterationNumberControl that allows us to fix the number of CG
iterations that are performed to a fixed small number (in our case 30):
@code
    IterationNumberControl iteration_number_control_aS(30, 1.e-18);
    SolverCG<Vector<double>>           solver_aS(iteration_number_control_aS);
    PreconditionIdentity preconditioner_aS;
    const auto preconditioner_S =
      inverse_operator(op_aS, solver_aS, preconditioner_aS);
@endcode

That's all!

Obviously, applying this inverse of the approximate Schur complement is a very
expensive preconditioner, almost as expensive as inverting the Schur
complement itself. We can expect it to significantly reduce the number of
outer iterations required for the Schur complement. In fact it does: in a
typical run on 7 times refined meshes using elements of order 0, the number of
outer iterations drops from 592 to 39. On the other hand, we now have to apply
a very expensive preconditioner 25 times. A better measure is therefore simply
the run-time of the program: on a current laptop (as of January 2019), it
drops from 3.57 to 2.05 seconds for this test case. That doesn't seem too
impressive, but the savings become more pronounced on finer meshes and with
elements of higher order. For example, an seven times refined mesh and
using elements of order 2 (which amounts to about 0.4 million degrees of
freedom) yields an improvement of 1134 to 83 outer iterations, at a runtime
of 168 seconds to 40 seconds. Not earth shattering, but significant.


<a name="Definitionofthetestcase"></a><h3>Definition of the test case</h3>


In this tutorial program, we will solve the Laplace equation in mixed
formulation as stated above. Since we want to monitor convergence of the
solution inside the program, we choose right hand side, boundary conditions,
and the coefficient so that we recover a solution function known to us. In
particular, we choose the pressure solution
@f{eqnarray*}
  p = -\left(\frac \alpha 2 xy^2 + \beta x - \frac \alpha 6 x^3\right),
@f}
and for the coefficient we choose the unit matrix $K_{ij}=\delta_{ij}$ for
simplicity. Consequently, the exact velocity satisfies
@f{eqnarray*}
  {\mathbf u} =
  \left(\begin{array}{cc}
    \frac \alpha 2 y^2 + \beta - \frac \alpha 2 x^2 \\
    \alpha xy
  \end{array}\right).
@f}
This solution was chosen since it is exactly divergence free, making it a
realistic test case for incompressible fluid flow. By consequence, the right
hand side equals $f=0$, and as boundary values we have to choose
$g=p|_{\partial\Omega}$.

For the computations in this program, we choose $\alpha=0.3,\beta=1$. You can
find the resulting solution in the <a name="#Results">results section
below</a>, after the commented program.
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
 * Since this program is only an adaptation of step-4, there is not much new
 * stuff in terms of header files. In deal.II, we usually list include files
 * in the order base-lac-grid-dofs-fe-numerics, followed by C++ standard
 * include files:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * @endcode
 * 
 * The only two new header files that deserve some attention are those for
 * the LinearOperator and PackagedOperation classes:
 * 
 * @code
 * #include <deal.II/lac/linear_operator.h>
 * #include <deal.II/lac/packaged_operation.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * This is the only significant new header, namely the one in which the
 * Raviart-Thomas finite element is declared:
 * 
 * @code
 * #include <deal.II/fe/fe_raviart_thomas.h>
 * 
 * @endcode
 * 
 * Finally, as a bonus in this program, we will use a tensorial
 * coefficient. Since it may have a spatial dependence, we consider it a
 * tensor-valued function. The following include file provides the
 * <code>TensorFunction</code> class that offers such functionality:
 * 
 * @code
 * #include <deal.II/base/tensor_function.h>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs: We put all of the code relevant
 * to this program into a namespace. (This idea was first introduced in step-7.)
 * 
 * @code
 * namespace Step20
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMixedLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MixedLaplaceProblem</code> class template</h3>
 * 

 * 
 * Again, since this is an adaptation of step-6, the main class is almost
 * the same as the one in that tutorial program. In terms of member
 * functions, the main differences are that the constructor takes the degree
 * of the Raviart-Thomas element as an argument (and that there is a
 * corresponding member variable to store this value) and the addition of
 * the <code>compute_error</code> function in which, no surprise, we will
 * compute the difference between the exact and the numerical solution to
 * determine convergence of our computations:
 * 
 * @code
 *   template <int dim>
 *   class MixedLaplaceProblem
 *   {
 *   public:
 *     MixedLaplaceProblem(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void make_grid_and_dofs();
 *     void assemble_system();
 *     void solve();
 *     void compute_errors() const;
 *     void output_results() const;
 * 
 *     const unsigned int degree;
 * 
 *     Triangulation<dim> triangulation;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 * @endcode
 * 
 * The second difference is that the sparsity pattern, the system matrix,
 * and solution and right hand side vectors are now blocked. What this
 * means and what one can do with such objects is explained in the
 * introduction to this program as well as further down below when we
 * explain the linear solvers and preconditioners for this problem:
 * 
 * @code
 *     BlockSparsityPattern      sparsity_pattern;
 *     BlockSparseMatrix<double> system_matrix;
 * 
 *     BlockVector<double> solution;
 *     BlockVector<double> system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandexactsolution"></a> 
 * <h3>Right hand side, boundary values, and exact solution</h3>
 * 

 * 
 * Our next task is to define the right hand side of our problem (i.e., the
 * scalar right hand side for the pressure in the original Laplace
 * equation), boundary values for the pressure, and a function that
 * describes both the pressure and the velocity of the exact solution for
 * later computations of the error. Note that these functions have one, one,
 * and <code>dim+1</code> components, respectively, and that we pass the
 * number of components down to the <code>Function@<dim@></code> base
 * class. For the exact solution, we only declare the function that actually
 * returns the entire solution vector (i.e. all components of it) at
 * once. Here are the respective declarations:
 * 
 * @code
 *   namespace PrescribedSolution
 *   {
 *     constexpr double alpha = 0.3;
 *     constexpr double beta  = 1;
 * 
 * 
 *     template <int dim>
 *     class RightHandSide : public Function<dim>
 *     {
 *     public:
 *       RightHandSide()
 *         : Function<dim>(1)
 *       {}
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     class PressureBoundaryValues : public Function<dim>
 *     {
 *     public:
 *       PressureBoundaryValues()
 *         : Function<dim>(1)
 *       {}
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 *     };
 * 
 * 
 *     template <int dim>
 *     class ExactSolution : public Function<dim>
 *     {
 *     public:
 *       ExactSolution()
 *         : Function<dim>(dim + 1)
 *       {}
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  value) const override;
 *     };
 * 
 * 
 * @endcode
 * 
 * And then we also have to define these respective functions, of
 * course. Given our discussion in the introduction of how the solution
 * should look, the following computations should be straightforward:
 * 
 * @code
 *     template <int dim>
 *     double RightHandSide<dim>::value(const Point<dim> & /*p*/,
 *                                      const unsigned int /*component*/) const
 *     {
 *       return 0;
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     double
 *     PressureBoundaryValues<dim>::value(const Point<dim> &p,
 *                                        const unsigned int /*component*/) const
 *     {
 *       return -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] -
 *                alpha * p[0] * p[0] * p[0] / 6);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     void ExactSolution<dim>::vector_value(const Point<dim> &p,
 *                                           Vector<double> &  values) const
 *     {
 *       Assert(values.size() == dim + 1,
 *              ExcDimensionMismatch(values.size(), dim + 1));
 * 
 *       values(0) = alpha * p[1] * p[1] / 2 + beta - alpha * p[0] * p[0] / 2;
 *       values(1) = alpha * p[0] * p[1];
 *       values(2) = -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] -
 *                     alpha * p[0] * p[0] * p[0] / 6);
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversepermeabilitytensor"></a> 
 * <h3>The inverse permeability tensor</h3>
 * 

 * 
 * In addition to the other equation data, we also want to use a
 * permeability tensor, or better -- because this is all that appears in the
 * weak form -- the inverse of the permeability tensor,
 * <code>KInverse</code>. For the purpose of verifying the exactness of the
 * solution and determining convergence orders, this tensor is more in the
 * way than helpful. We will therefore simply set it to the identity matrix.
 *     

 * 
 * However, a spatially varying permeability tensor is indispensable in
 * real-life porous media flow simulations, and we would like to use the
 * opportunity to demonstrate the technique to use tensor valued functions.
 *     

 * 
 * Possibly unsurprisingly, deal.II also has a base class not only for
 * scalar and generally vector-valued functions (the <code>Function</code>
 * base class) but also for functions that return tensors of fixed dimension
 * and rank, the <code>TensorFunction</code> template. Here, the function
 * under consideration returns a dim-by-dim matrix, i.e. a tensor of rank 2
 * and dimension <code>dim</code>. We then choose the template arguments of
 * the base class appropriately.
 *     

 * 
 * The interface that the <code>TensorFunction</code> class provides is
 * essentially equivalent to the <code>Function</code> class. In particular,
 * there exists a <code>value_list</code> function that takes a list of
 * points at which to evaluate the function, and returns the values of the
 * function in the second argument, a list of tensors:
 * 
 * @code
 *     template <int dim>
 *     class KInverse : public TensorFunction<2, dim>
 *     {
 *     public:
 *       KInverse()
 *         : TensorFunction<2, dim>()
 *       {}
 * 
 *       virtual void
 *       value_list(const std::vector<Point<dim>> &points,
 *                  std::vector<Tensor<2, dim>> &  values) const override;
 *     };
 * 
 * 
 * @endcode
 * 
 * The implementation is less interesting. As in previous examples, we add a
 * check to the beginning of the class to make sure that the sizes of input
 * and output parameters are the same (see step-5 for a discussion of this
 * technique). Then we loop over all evaluation points, and for each one
 * set the output tensor to the identity matrix.
 *     

 * 
 * There is an oddity at the top of the function (the
 * `(void)points;` statement) that is worth discussing. The values
 * we put into the output `values` array does not actually depend
 * on the `points` arrays of coordinates at which the function is
 * evaluated. In other words, the `points` argument is in fact
 * unused, and we could have just not given it a name if we had
 * wanted. But we want to use the `points` object for checking
 * that the `values` object has the correct size. The problem is
 * that in release mode, `AssertDimension` is defined as a macro
 * that expands to nothing; the compiler will then complain that
 * the `points` object is unused. The idiomatic approach to
 * silencing this warning is to have a statement that evaluates
 * (reads) variable but doesn't actually do anything: That's what
 * `(void)points;` does: It reads from `points`, and then casts
 * the result of the read to `void`, i.e., nothing. This statement
 * is, in other words, completely pointless and implies no actual
 * action except to explain to the compiler that yes, this
 * variable is in fact used even in release mode. (In debug mode,
 * the `AssertDimension` macro expands to something that reads
 * from the variable, and so the funny statement would not be
 * necessary in debug mode.)
 * 
 * @code
 *     template <int dim>
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                    std::vector<Tensor<2, dim>> &  values) const
 *     {
 *       (void)points;
 *       AssertDimension(points.size(), values.size());
 * 
 *       for (auto &value : values)
 *         value = unit_symmetric_tensor<dim>();
 *     }
 *   } // namespace PrescribedSolution
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemclassimplementation"></a> 
 * <h3>MixedLaplaceProblem class implementation</h3>
 * 

 * 
 * 
 * <a name="MixedLaplaceProblemMixedLaplaceProblem"></a> 
 * <h4>MixedLaplaceProblem::MixedLaplaceProblem</h4>
 * 

 * 
 * In the constructor of this class, we first store the value that was
 * passed in concerning the degree of the finite elements we shall use (a
 * degree of zero, for example, means to use RT(0) and DG(0)), and then
 * construct the vector valued element belonging to the space $X_h$ described
 * in the introduction. The rest of the constructor is as in the early
 * tutorial programs.
 *   

 * 
 * The only thing worth describing here is the constructor call of the
 * <code>fe</code> variable. The <code>FESystem</code> class to which this
 * variable belongs has a number of different constructors that all refer to
 * binding simpler elements together into one larger element. In the present
 * case, we want to couple a single RT(degree) element with a single
 * DQ(degree) element. The constructor to <code>FESystem</code> that does
 * this requires us to specify first the first base element (the
 * <code>FE_RaviartThomas</code> object of given degree) and then the number
 * of copies for this base element, and then similarly the kind and number
 * of <code>FE_DGQ</code> elements. Note that the Raviart-Thomas element
 * already has <code>dim</code> vector components, so that the coupled
 * element will have <code>dim+1</code> vector components, the first
 * <code>dim</code> of which correspond to the velocity variable whereas the
 * last one corresponds to the pressure.
 *   

 * 
 * It is also worth comparing the way we constructed this element from its
 * base elements, with the way we have done so in step-8: there, we have
 * built it as <code>fe (FE_Q@<dim@>(1), dim)</code>, i.e. we have simply
 * used <code>dim</code> copies of the <code>FE_Q(1)</code> element, one
 * copy for the displacement in each coordinate direction.
 * 
 * @code
 *   template <int dim>
 *   MixedLaplaceProblem<dim>::MixedLaplaceProblem(const unsigned int degree)
 *     : degree(degree)
 *     , fe(FE_RaviartThomas<dim>(degree), 1, FE_DGQ<dim>(degree), 1)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemmake_grid_and_dofs"></a> 
 * <h4>MixedLaplaceProblem::make_grid_and_dofs</h4>
 * 

 * 
 * This next function starts out with well-known functions calls that create
 * and refine a mesh, and then associate degrees of freedom with it:
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -1, 1);
 *     triangulation.refine_global(5);
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * However, then things become different. As mentioned in the
 * introduction, we want to subdivide the matrix into blocks corresponding
 * to the two different kinds of variables, velocity and pressure. To this
 * end, we first have to make sure that the indices corresponding to
 * velocities and pressures are not intermingled: First all velocity
 * degrees of freedom, then all pressure DoFs. This way, the global matrix
 * separates nicely into a $2 \times 2$ system. To achieve this, we have to
 * renumber degrees of freedom based on their vector component, an
 * operation that conveniently is already implemented:
 * 
 * @code
 *     DoFRenumbering::component_wise(dof_handler);
 * 
 * @endcode
 * 
 * The next thing is that we want to figure out the sizes of these blocks
 * so that we can allocate an appropriate amount of space. To this end, we
 * call the DoFTools::count_dofs_per_fe_component() function that
 * counts how many shape functions are non-zero for a particular vector
 * component. We have <code>dim+1</code> vector components, and
 * DoFTools::count_dofs_per_fe_component() will count how many shape
 * functions belong to each of these components.
 *     

 * 
 * There is one problem here. As described in the documentation of that
 * function, it <i>wants</i> to put the number of $x$-velocity shape
 * functions into <code>dofs_per_component[0]</code>, the number of
 * $y$-velocity shape functions into <code>dofs_per_component[1]</code>
 * (and similar in 3d), and the number of pressure shape functions into
 * <code>dofs_per_component[dim]</code>. But, the Raviart-Thomas element
 * is special in that it is non-@ref GlossPrimitive "primitive", i.e.,
 * for Raviart-Thomas elements all velocity shape functions
 * are nonzero in all components. In other words, the function cannot
 * distinguish between $x$ and $y$ velocity functions because there
 * <i>is</i> no such distinction. It therefore puts the overall number
 * of velocity into each of <code>dofs_per_component[c]</code>,
 * $0\le c\le \text{dim}$. On the other hand, the number
 * of pressure variables equals the number of shape functions that are
 * nonzero in the dim-th component.
 *     

 * 
 * Using this knowledge, we can get the number of velocity shape
 * functions from any of the first <code>dim</code> elements of
 * <code>dofs_per_component</code>, and then use this below to initialize
 * the vector and matrix block sizes, as well as create output.
 *     

 * 
 * @note If you find this concept difficult to understand, you may
 * want to consider using the function DoFTools::count_dofs_per_fe_block()
 * instead, as we do in the corresponding piece of code in step-22.
 * You might also want to read up on the difference between
 * @ref GlossBlock "blocks" and @ref GlossComponent "components"
 * in the glossary.
 * 
 * @code
 *     const std::vector<types::global_dof_index> dofs_per_component =
 *       DoFTools::count_dofs_per_fe_component(dof_handler);
 *     const unsigned int n_u = dofs_per_component[0],
 *                        n_p = dofs_per_component[dim];
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Total number of cells: " << triangulation.n_cells()
 *               << std::endl
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (" << n_u << '+' << n_p << ')' << std::endl;
 * 
 * @endcode
 * 
 * The next task is to allocate a sparsity pattern for the matrix that we
 * will create. We use a compressed sparsity pattern like in the previous
 * steps, but as <code>system_matrix</code> is a block matrix we use the
 * class <code>BlockDynamicSparsityPattern</code> instead of just
 * <code>DynamicSparsityPattern</code>. This block sparsity pattern has
 * four blocks in a $2 \times 2$ pattern. The blocks' sizes depend on
 * <code>n_u</code> and <code>n_p</code>, which hold the number of velocity
 * and pressure variables. In the second step we have to instruct the block
 * system to update its knowledge about the sizes of the blocks it manages;
 * this happens with the <code>dsp.collect_sizes ()</code> call.
 * 
 * @code
 *     BlockDynamicSparsityPattern dsp(2, 2);
 *     dsp.block(0, 0).reinit(n_u, n_u);
 *     dsp.block(1, 0).reinit(n_p, n_u);
 *     dsp.block(0, 1).reinit(n_u, n_p);
 *     dsp.block(1, 1).reinit(n_p, n_p);
 *     dsp.collect_sizes();
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 * @endcode
 * 
 * We use the compressed block sparsity pattern in the same way as the
 * non-block version to create the sparsity pattern and then the system
 * matrix:
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 * @endcode
 * 
 * Then we have to resize the solution and right hand side vectors in
 * exactly the same way as the block compressed sparsity pattern:
 * 
 * @code
 *     solution.reinit(2);
 *     solution.block(0).reinit(n_u);
 *     solution.block(1).reinit(n_p);
 *     solution.collect_sizes();
 * 
 *     system_rhs.reinit(2);
 *     system_rhs.block(0).reinit(n_u);
 *     system_rhs.block(1).reinit(n_p);
 *     system_rhs.collect_sizes();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemassemble_system"></a> 
 * <h4>MixedLaplaceProblem::assemble_system</h4>
 * 

 * 
 * Similarly, the function that assembles the linear system has mostly been
 * discussed already in the introduction to this example. At its top, what
 * happens are all the usual steps, with the addition that we do not only
 * allocate quadrature and <code>FEValues</code> objects for the cell terms,
 * but also for face terms. After that, we define the usual abbreviations
 * for variables, and the allocate space for the local matrix and right hand
 * side contributions, and the array that holds the global numbers of the
 * degrees of freedom local to the present cell.
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim>     quadrature_formula(degree + 2);
 *     QGauss<dim - 1> face_quadrature_formula(degree + 2);
 * 
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * The next step is to declare objects that represent the source term,
 * pressure boundary value, and coefficient in the equation. In addition
 * to these objects that represent continuous functions, we also need
 * arrays to hold their values at the quadrature points of individual
 * cells (or faces, for the boundary values). Note that in the case of the
 * coefficient, the array has to be one of matrices.
 * 
 * @code
 *     const PrescribedSolution::RightHandSide<dim> right_hand_side;
 *     const PrescribedSolution::PressureBoundaryValues<dim>
 *                                             pressure_boundary_values;
 *     const PrescribedSolution::KInverse<dim> k_inverse;
 * 
 *     std::vector<double>         rhs_values(n_q_points);
 *     std::vector<double>         boundary_values(n_face_q_points);
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 * @endcode
 * 
 * Finally, we need a couple of extractors that we will use to get at the
 * velocity and pressure components of vector-valued shape
 * functions. Their function and use is described in detail in the @ref
 * vector_valued report. Essentially, we will use them as subscripts on
 * the FEValues objects below: the FEValues object describes all vector
 * components of shape functions, while after subscription, it will only
 * refer to the velocities (a set of <code>dim</code> components starting
 * at component zero) or the pressure (a scalar component located at
 * position <code>dim</code>):
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 * @endcode
 * 
 * With all this in place, we can go on with the loop over all cells. The
 * body of this loop has been discussed in the introduction, and will not
 * be commented any further here:
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                    rhs_values);
 *         k_inverse.value_list(fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q);
 *               const double div_phi_i_u = fe_values[velocities].divergence(i, q);
 *               const double phi_i_p     = fe_values[pressure].value(i, q);
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const Tensor<1, dim> phi_j_u =
 *                     fe_values[velocities].value(j, q);
 *                   const double div_phi_j_u =
 *                     fe_values[velocities].divergence(j, q);
 *                   const double phi_j_p = fe_values[pressure].value(j, q);
 * 
 *                   local_matrix(i, j) +=
 *                     (phi_i_u * k_inverse_values[q] * phi_j_u 
 *                      - phi_i_p * div_phi_j_u                 
 *                      - div_phi_i_u * phi_j_p)                
 *                     * fe_values.JxW(q);
 *                 }
 * 
 *               local_rhs(i) += -phi_i_p * rhs_values[q] * fe_values.JxW(q);
 *             }
 * 
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary())
 *             {
 *               fe_face_values.reinit(cell, face);
 * 
 *               pressure_boundary_values.value_list(
 *                 fe_face_values.get_quadrature_points(), boundary_values);
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   local_rhs(i) += -(fe_face_values[velocities].value(i, q) * 
 *                                     fe_face_values.normal_vector(q) *        
 *                                     boundary_values[q] *                     
 *                                     fe_face_values.JxW(q));
 *             }
 * 
 * @endcode
 * 
 * The final step in the loop over all cells is to transfer local
 * contributions into the global matrix and right hand side
 * vector. Note that we use exactly the same interface as in previous
 * examples, although we now use block matrices and vectors instead of
 * the regular ones. In other words, to the outside world, block
 * objects have the same interface as matrices and vectors, but they
 * additionally allow to access individual blocks.
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               local_matrix(i, j));
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           system_rhs(local_dof_indices[i]) += local_rhs(i);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Implementationoflinearsolversandpreconditioners"></a> 
 * <h3>Implementation of linear solvers and preconditioners</h3>
 * 

 * 
 * The linear solvers and preconditioners we use in this example have
 * been discussed in significant detail already in the introduction. We
 * will therefore not discuss the rationale for our approach here any
 * more, but rather only comment on some remaining implementational
 * aspects.
 * 

 * 
 * 
 * <a name="MixedLaplacesolve"></a> 
 * <h4>MixedLaplace::solve</h4>
 * 

 * 
 * As already outlined in the introduction, the solve function consists
 * essentially of two steps. First, we have to form the first equation
 * involving the Schur complement and solve for the pressure (component 1
 * of the solution). Then, we can reconstruct the velocities from the
 * second equation (component 0 of the solution).
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::solve()
 *   {
 * @endcode
 * 
 * As a first step we declare references to all block components of the
 * matrix, the right hand side and the solution vector that we will
 * need.
 * 
 * @code
 *     const auto &M = system_matrix.block(0, 0);
 *     const auto &B = system_matrix.block(0, 1);
 * 
 *     const auto &F = system_rhs.block(0);
 *     const auto &G = system_rhs.block(1);
 * 
 *     auto &U = solution.block(0);
 *     auto &P = solution.block(1);
 * 
 * @endcode
 * 
 * Then, we will create corresponding LinearOperator objects and create
 * the <code>op_M_inv</code> operator:
 * 
 * @code
 *     const auto op_M = linear_operator(M);
 *     const auto op_B = linear_operator(B);
 * 
 *     ReductionControl         reduction_control_M(2000, 1.0e-18, 1.0e-10);
 *     SolverCG<Vector<double>> solver_M(reduction_control_M);
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner_M;
 * 
 *     preconditioner_M.initialize(M);
 * 
 *     const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
 * 
 * @endcode
 * 
 * This allows us to declare the Schur complement <code>op_S</code> and
 * the approximate Schur complement <code>op_aS</code>:
 * 
 * @code
 *     const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
 *     const auto op_aS =
 *       transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B;
 * 
 * @endcode
 * 
 * We now create a preconditioner out of <code>op_aS</code> that
 * applies a fixed number of 30 (inexpensive) CG iterations:
 * 
 * @code
 *     IterationNumberControl   iteration_number_control_aS(30, 1.e-18);
 *     SolverCG<Vector<double>> solver_aS(iteration_number_control_aS);
 * 
 *     const auto preconditioner_S =
 *       inverse_operator(op_aS, solver_aS, PreconditionIdentity());
 * 
 * @endcode
 * 
 * Now on to the first equation. The right hand side of it is
 * $B^TM^{-1}F-G$, which is what we compute in the first few lines. We
 * then solve the first equation with a CG solver and the
 * preconditioner we just declared.
 * 
 * @code
 *     const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G;
 * 
 *     SolverControl            solver_control_S(2000, 1.e-12);
 *     SolverCG<Vector<double>> solver_S(solver_control_S);
 * 
 *     const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S);
 * 
 *     P = op_S_inv * schur_rhs;
 * 
 *     std::cout << solver_control_S.last_step()
 *               << " CG Schur complement iterations to obtain convergence."
 *               << std::endl;
 * 
 * @endcode
 * 
 * After we have the pressure, we can compute the velocity. The equation
 * reads $MU=-BP+F$, and we solve it by first computing the right hand
 * side, and then multiplying it with the object that represents the
 * inverse of the mass matrix:
 * 
 * @code
 *     U = op_M_inv * (F - op_B * P);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemclassimplementationcontinued"></a> 
 * <h3>MixedLaplaceProblem class implementation (continued)</h3>
 * 

 * 
 * 
 * <a name="MixedLaplacecompute_errors"></a> 
 * <h4>MixedLaplace::compute_errors</h4>
 * 

 * 
 * After we have dealt with the linear solver and preconditioners, we
 * continue with the implementation of our main class. In particular, the
 * next task is to compute the errors in our numerical solution, in both the
 * pressures as well as velocities.
 *   

 * 
 * To compute errors in the solution, we have already introduced the
 * <code>VectorTools::integrate_difference</code> function in step-7 and
 * step-11. However, there we only dealt with scalar solutions, whereas here
 * we have a vector-valued solution with components that even denote
 * different quantities and may have different orders of convergence (this
 * isn't the case here, by choice of the used finite elements, but is
 * frequently the case in mixed finite element applications). What we
 * therefore have to do is to `mask' the components that we are interested
 * in. This is easily done: the
 * <code>VectorTools::integrate_difference</code> function takes as one of its
 * arguments a pointer to a weight function (the parameter defaults to the
 * null pointer, meaning unit weights). What we have to do is to pass
 * a function object that equals one in the components we are interested in,
 * and zero in the other ones. For example, to compute the pressure error,
 * we should pass a function that represents the constant vector with a unit
 * value in component <code>dim</code>, whereas for the velocity the
 * constant vector should be one in the first <code>dim</code> components,
 * and zero in the location of the pressure.
 *   

 * 
 * In deal.II, the <code>ComponentSelectFunction</code> does exactly this:
 * it wants to know how many vector components the function it is to
 * represent should have (in our case this would be <code>dim+1</code>, for
 * the joint velocity-pressure space) and which individual or range of
 * components should be equal to one. We therefore define two such masks at
 * the beginning of the function, following by an object representing the
 * exact solution and a vector in which we will store the cellwise errors as
 * computed by <code>integrate_difference</code>:
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::compute_errors() const
 *   {
 *     const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
 *     const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
 *                                                      dim + 1);
 * 
 *     PrescribedSolution::ExactSolution<dim> exact_solution;
 *     Vector<double> cellwise_errors(triangulation.n_active_cells());
 * 
 * @endcode
 * 
 * As already discussed in step-7, we have to realize that it is
 * impossible to integrate the errors exactly. All we can do is
 * approximate this integral using quadrature. This actually presents a
 * slight twist here: if we naively chose an object of type
 * <code>QGauss@<dim@>(degree+1)</code> as one may be inclined to do (this
 * is what we used for integrating the linear system), one realizes that
 * the error is very small and does not follow the expected convergence
 * curves at all. What is happening is that for the mixed finite elements
 * used here, the Gauss points happen to be superconvergence points in
 * which the pointwise error is much smaller (and converges with higher
 * order) than anywhere else. These are therefore not particularly good
 * points for integration. To avoid this problem, we simply use a
 * trapezoidal rule and iterate it <code>degree+2</code> times in each
 * coordinate direction (again as explained in step-7):
 * 
 * @code
 *     QTrapezoid<1>  q_trapez;
 *     QIterated<dim> quadrature(q_trapez, degree + 2);
 * 
 * @endcode
 * 
 * With this, we can then let the library compute the errors and output
 * them to the screen:
 * 
 * @code
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       exact_solution,
 *                                       cellwise_errors,
 *                                       quadrature,
 *                                       VectorTools::L2_norm,
 *                                       &pressure_mask);
 *     const double p_l2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         cellwise_errors,
 *                                         VectorTools::L2_norm);
 * 
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       exact_solution,
 *                                       cellwise_errors,
 *                                       quadrature,
 *                                       VectorTools::L2_norm,
 *                                       &velocity_mask);
 *     const double u_l2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         cellwise_errors,
 *                                         VectorTools::L2_norm);
 * 
 *     std::cout << "Errors: ||e_p||_L2 = " << p_l2_error
 *               << ",   ||e_u||_L2 = " << u_l2_error << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceoutput_results"></a> 
 * <h4>MixedLaplace::output_results</h4>
 * 

 * 
 * The last interesting function is the one in which we generate graphical
 * output. Note that all velocity components get the same solution name
 * "u". Together with using
 * DataComponentInterpretation::component_is_part_of_vector this will
 * cause DataOut<dim>::write_vtu() to generate a vector representation of
 * the individual velocity components, see step-22 or the
 * @ref VVOutput "Generating graphical output"
 * section of the
 * @ref vector_valued
 * module for more information. Finally, it seems inappropriate for higher
 * order elements to only show a single bilinear quadrilateral per cell in
 * the graphical output. We therefore generate patches of size
 * (degree+1)x(degree+1) to capture the full information content of the
 * solution. See the step-7 tutorial program for more information on this.
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::output_results() const
 *   {
 *     std::vector<std::string> solution_names(dim, "u");
 *     solution_names.emplace_back("p");
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       interpretation(dim,
 *                      DataComponentInterpretation::component_is_part_of_vector);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 * 
 *     DataOut<dim> data_out;
 *     data_out.add_data_vector(dof_handler,
 *                              solution,
 *                              solution_names,
 *                              interpretation);
 * 
 *     data_out.build_patches(degree + 1);
 * 
 *     std::ofstream output("solution.vtu");
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplacerun"></a> 
 * <h4>MixedLaplace::run</h4>
 * 

 * 
 * This is the final function of our main class. It's only job is to call
 * the other functions in their natural order:
 * 
 * @code
 *   template <int dim>
 *   void MixedLaplaceProblem<dim>::run()
 *   {
 *     make_grid_and_dofs();
 *     assemble_system();
 *     solve();
 *     compute_errors();
 *     output_results();
 *   }
 * } // namespace Step20
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The main function we stole from step-6 instead of step-4. It is almost
 * equal to the one in step-6 (apart from the changed class names, of course),
 * the only exception is that we pass the degree of the finite element space
 * to the constructor of the mixed Laplace problem (here, we use zero-th order
 * elements).
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step20;
 * 
 *       const unsigned int     fe_degree = 0;
 *       MixedLaplaceProblem<2> mixed_laplace_problem(fe_degree);
 *       mixed_laplace_problem.run();
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


<a name="Outputoftheprogramandgraphicalvisualization"></a><h3>Output of the program and graphical visualization</h3>



If we run the program as is, we get this output for the $32\times 32$
mesh we use (for a total of 1024 cells with 1024 pressure degrees of
freedom since we use piecewise constants, and 2112 velocities because
the Raviart-Thomas element defines one degree per freedom per face and
there are $1024 + 32 = 1056$ faces parallel to the $x$-axis and the same
number parallel to the $y$-axis):
@verbatim
\$ make run
[ 66%] Built target \step-20
Scanning dependencies of target run
[100%] Run \step-20 with Release configuration
Number of active cells: 1024
Total number of cells: 1365
Number of degrees of freedom: 3136 (2112+1024)
24 CG Schur complement iterations to obtain convergence.
Errors: ||e_p||_L2 = 0.0445032,   ||e_u||_L2 = 0.010826
[100%] Built target run
@endverbatim

The fact that the number of iterations is so small, of course, is due to
the good (but expensive!) preconditioner we have developed. To get
confidence in the solution, let us take a look at it. The following three
images show (from left to right) the x-velocity, the y-velocity, and the
pressure:

<table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u_new.jpg" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.v_new.jpg" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.p_new.jpg" width="400" alt=""></td>
  </tr>
</table>



Let us start with the pressure: it is highest at the left and lowest at the
right, so flow will be from left to right. In addition, though hardly visible
in the graph, we have chosen the pressure field such that the flow left-right
flow first channels towards the center and then outward again. Consequently,
the x-velocity has to increase to get the flow through the narrow part,
something that can easily be seen in the left image. The middle image
represents inward flow in y-direction at the left end of the domain, and
outward flow in y-direction at the right end of the domain.



As an additional remark, note how the x-velocity in the left image is only
continuous in x-direction, whereas the y-velocity is continuous in
y-direction. The flow fields are discontinuous in the other directions. This
very obviously reflects the continuity properties of the Raviart-Thomas
elements, which are, in fact, only in the space H(div) and not in the space
$H^1$. Finally, the pressure field is completely discontinuous, but
that should not surprise given that we have chosen <code>FE_DGQ(0)</code> as
the finite element for that solution component.



<a name="Convergence"></a><h3>Convergence</h3>



The program offers two obvious places where playing and observing convergence
is in order: the degree of the finite elements used (passed to the constructor
of the <code>MixedLaplaceProblem</code> class from <code>main()</code>), and
the refinement level (determined in
<code>MixedLaplaceProblem::make_grid_and_dofs</code>). What one can do is to
change these values and observe the errors computed later on in the course of
the program run.



If one does this, one finds the following pattern for the $L_2$ error
in the pressure variable:
<table align="center" class="doxtable">
  <tr>
    <th></th>
    <th colspan="3" align="center">Finite element order</th>
  </tr>
  <tr>
    <th>Refinement level</th>
    <th>0</th>
    <th>1</th>
    <th>2</th>
  </tr>
  <tr>
    <th>0</th>  <td>1.45344</td>  <td>0.0831743</td>  <td>0.0235186</td>
  </tr>
  <tr>
    <th>1</th>  <td>0.715099</td>  <td>0.0245341</td>  <td>0.00293983</td>
  </tr>
  <tr>
    <th>2</th>  <td>0.356383</td>  <td>0.0063458</td>  <td>0.000367478</td>
  </tr>
  <tr>
    <th>3</th>  <td>0.178055</td>  <td>0.00159944</td>  <td>4.59349e-05</td>
  </tr>
  <tr>
    <th>4</th>  <td>0.0890105</td>  <td>0.000400669</td>  <td>5.74184e-06</td>
  </tr>
  <tr>
    <th>5</th>  <td>0.0445032</td>  <td>0.000100218</td>  <td>7.17799e-07</td>
  </tr>
  <tr>
    <th>6</th>  <td>0.0222513</td>  <td>2.50576e-05</td>  <td>9.0164e-08</td>
  </tr>
  <tr>
    <th></th>  <th>$O(h)$</th>  <th>$O(h^2)$</th>  <th>$O(h^3)$</th>
  </tr>
</table>

The theoretically expected convergence orders are very nicely reflected by the
experimentally observed ones indicated in the last row of the table.



One can make the same experiment with the $L_2$ error
in the velocity variables:
<table align="center" class="doxtable">
  <tr>
    <th></th>
    <th colspan="3" align="center">Finite element order</th>
  </tr>
  <tr>
    <th>Refinement level</th>
    <th>0</th>
    <th>1</th>
    <th>2</th>
  </tr>
  <tr>
    <th>0</th> <td>0.367423</td> <td>0.127657</td> <td>5.10388e-14</td>
  </tr>
  <tr>
    <th>1</th> <td>0.175891</td> <td>0.0319142</td> <td>9.04414e-15</td>
  </tr>
  <tr>
    <th>2</th> <td>0.0869402</td> <td>0.00797856</td> <td>1.23723e-14</td>
  </tr>
  <tr>
    <th>3</th> <td>0.0433435</td> <td>0.00199464</td> <td>1.86345e-07</td>
  </tr>
  <tr>
    <th>4</th> <td>0.0216559</td> <td>0.00049866</td> <td>2.72566e-07</td>
  </tr>
  <tr>
    <th>5</th> <td>0.010826</td> <td>0.000124664</td> <td>3.57141e-07</td>
  </tr>
  <tr>
    <th>6</th> <td>0.00541274</td> <td>3.1166e-05</td> <td>4.46124e-07</td>
  </tr>
  <tr>
    <th></th>  <td>$O(h)$</td>  <td>$O(h^2)$</td>  <td>$O(h^3)$</td>
  </tr>
</table>
The result concerning the convergence order is the same here.



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Morerealisticpermeabilityfields"></a><h4>More realistic permeability fields</h4>


Realistic flow computations for ground water or oil reservoir simulations will
not use a constant permeability. Here's a first, rather simple way to change
this situation: we use a permeability that decays very rapidly away from a
central flowline until it hits a background value of 0.001. This is to mimic
the behavior of fluids in sandstone: in most of the domain, the sandstone is
homogeneous and, while permeable to fluids, not overly so; on the other stone,
the stone has cracked, or faulted, along one line, and the fluids flow much
easier along this large crack. Here is how we could implement something like
this:
@code
template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();

      const double distance_to_flowline
        = std::fabs(points[p][1]-0.2*std::sin(10*points[p][0]));

      const double permeability = std::max(std::exp(-(distance_to_flowline*
                                                      distance_to_flowline)
                                                    / (0.1 * 0.1)),
                                           0.001);

      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 1./permeability;
    }
}
@endcode
Remember that the function returns the inverse of the permeability tensor.



With a significantly higher mesh resolution, we can visualize this, here with
x- and y-velocity:

<table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u-wiggle.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.v-wiggle.png" alt=""></td>
  </tr>
</table>

It is obvious how fluids flow essentially only along the middle line, and not
anywhere else.



Another possibility would be to use a random permeability field. A simple way
to achieve this would be to scatter a number of centers around the domain and
then use a permeability field that is the sum of (negative) exponentials for
each of these centers. Flow would then try to hop from one center of high
permeability to the next one. This is an entirely unscientific attempt at
describing a random medium, but one possibility to implement this behavior
would look like this:
@code
template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    KInverse ();

    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<Tensor<2,dim> >    &values) const;

  private:
    std::vector<Point<dim> > centers;
};


template <int dim>
KInverse<dim>::KInverse ()
{
  const unsigned int N = 40;
  centers.resize (N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int d=0; d<dim; ++d)
      centers[i][d] = 2.*rand()/RAND_MAX-1;
}


template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();

      double permeability = 0;
      for (unsigned int i=0; i<centers.size(); ++i)
        permeability += std::exp(-(points[p] - centers[i]).norm_square() / (0.1 * 0.1));

      const double normalized_permeability
        = std::max(permeability, 0.005);

      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 1./normalized_permeability;
    }
}
@endcode

A piecewise constant interpolation of the diagonal elements of the
inverse of this tensor (i.e., of <code>normalized_permeability</code>)
looks as follows:

<img src="https://www.dealii.org/images/steps/developer/step-20.k-random.png" alt="">


With a permeability field like this, we would get x-velocities and pressures as
follows:

<table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u-random.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.p-random.png" alt=""></td>
  </tr>
</table>

We will use these permeability fields again in step-21 and step-43.


<a name="Betterlinearsolvers"></a><h4>Better linear solvers</h4>


As mentioned in the introduction, the Schur complement solver used here is not
the best one conceivable (nor is it intended to be a particularly good
one). Better ones can be found in the literature and can be built using the
same block matrix techniques that were introduced here. We pick up on this
theme again in step-22, where we first build a Schur complement solver for the
Stokes equation as we did here, and then in the <a
href="step_22.html#improved-solver">Improved Solvers</a> section discuss better
ways based on solving the system as a whole but preconditioning based on
individual blocks. We will also come back to this in step-43.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-20.cc"
*/
