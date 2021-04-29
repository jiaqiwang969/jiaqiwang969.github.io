/**
@page step_61 The step-61 tutorial program
This tutorial depends on step-51.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#WeakGalerkinfiniteelementmethods"> Weak Galerkin finite element methods </a>
        <li><a href="#Theequationtosolve"> The equation to solve </a>
        <li><a href="#WeakGalerkinscheme"> Weak Galerkin scheme </a>
        <li><a href="#Representingtheweakgradient"> Representing the weak gradient </a>
        <li><a href="#Assemblingthelinearsystem"> Assembling the linear system </a>
        <li><a href="#PostprocessingandiLsub2subierrors"> Post-processing and <i>L<sub>2</sub></i>-errors </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#TheWGDarcyEquationclasstemplate">The WGDarcyEquation class template</a>
        <li><a href="#Righthandsideboundaryvaluesandexactsolution">Right hand side, boundary values, and exact solution</a>
        <li><a href="#WGDarcyEquationclassimplementation">WGDarcyEquation class implementation</a>
      <ul>
        <li><a href="#WGDarcyEquationWGDarcyEquation">WGDarcyEquation::WGDarcyEquation</a>
        <li><a href="#WGDarcyEquationmake_grid">WGDarcyEquation::make_grid</a>
        <li><a href="#WGDarcyEquationsetup_system">WGDarcyEquation::setup_system</a>
        <li><a href="#WGDarcyEquationassemble_system">WGDarcyEquation::assemble_system</a>
        <li><a href="#WGDarcyEquationdimsolve">WGDarcyEquation<dim>::solve</a>
        <li><a href="#WGDarcyEquationdimcompute_postprocessed_velocity">WGDarcyEquation<dim>::compute_postprocessed_velocity</a>
        <li><a href="#WGDarcyEquationdimcompute_pressure_error">WGDarcyEquation<dim>::compute_pressure_error</a>
        <li><a href="#WGDarcyEquationdimcompute_velocity_error">WGDarcyEquation<dim>::compute_velocity_error</a>
        <li><a href="#WGDarcyEquationoutput_results">WGDarcyEquation::output_results</a>
        <li><a href="#WGDarcyEquationrun">WGDarcyEquation::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestresultsoniWGQsub0subQsub0subRTsub0subi">Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik0i">Convergence table for <i>k=0</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub1subQsub1subRTsub1subi">Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik1i">Convergence table for <i>k=1</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub2subQsub2subRTsub2subi">Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik2i">Convergence table for <i>k=2</i></a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program was contributed by Zhuoran Wang.
Some more information about this program, as well as more numerical
results, are presented in @cite Wang2019 .
</i>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This tutorial program presents an implementation of the "weak Galerkin"
finite element method for the Poisson equation. In some sense, the motivation for
considering this method starts from the same point as in step-51: We would like to
consider discontinuous shape functions, but then need to address the fact that
the resulting problem has a much larger number of degrees of freedom compared to
the usual continuous Galerkin method (because, for
example, each vertex carries as many degrees of freedom as there are adjacent cells).
We also have to address the fact that, unlike in the continuous
Galerkin method, <i>every</i> degree of freedom
on one cell couples with all of the degrees of freedom on each of its face neighbor
cells. Consequently, the matrix one gets from the "traditional" discontinuous
Galerkin methods are both large and relatively dense.

Both the hybridized discontinuous Galerkin method (HDG) in step-51 and the weak
Galerkin (WG) method in this tutorial address the issue of coupling by introducing
additional degrees of freedom whose shape functions only live on a face between
cells (i.e., on the "skeleton" of the mesh), and which therefore "insulate" the
degrees of freedom on the adjacent cells from each other: cell degrees of freedom
only couple with other cell degrees of freedom on the same cell, as well as face
degrees of freedom, but not with cell degrees of freedom on neighboring cells.
Consequently, the coupling of shape functions for these cell degrees of freedom
indeed couple on exactly one cell and the degrees of freedom defined on its
faces.

For a given equation, say the second order Poisson equation,
the difference between the HDG and the WG method is how precisely one formulates
the problem that connects all of these different shape functions. (Indeed,
for some WG and HDG formulation, it is possible to show that they are equivalent.)
The HDG does things by reformulating second order problems in terms of a system of first
order equations and then conceptually considers the face degrees of freedom
to be "fluxes" of this first order system. In contrast, the WG method keeps things
in second order form and considers the face degrees of freedom as of the same
type as the primary solution variable, just restricted to the lower-dimensional
faces. For the purposes of the equation, one then needs to somehow "extend"
these shape functions into the interior of the cell when defining what it means
to apply a differential operator to them. Compared to the HDG, the method
has the advantage that it does not lead to a proliferation of unknowns due
to rewriting the equation as a first-order system, but it is also not quite
as easy to implement. However, as we will see in the following, this
additional effort is not prohibitive.


<a name="WeakGalerkinfiniteelementmethods"></a><h3> Weak Galerkin finite element methods </h3>


Weak Galerkin Finite Element Methods (WGFEMs) use discrete weak functions
to approximate scalar unknowns, and discrete weak gradients to
approximate classical gradients.
The method was originally introduced by Junping Wang and Xiu Ye
in the paper
<a href="https://doi.org/10.1016/j.cam.2012.10.003">
<i>A weak Galerkin finite element method for second order elliptic problems</i>,
J. Comput. Appl. Math., 103-115, 2013</a>.
Compared to the continuous Galerkin method,
the weak Galerkin method satisfies important physical properties, namely
local mass conservation and bulk normal flux continuity.
It results in a SPD linear system, and optimal convergence rates can
be obtained with mesh refinement.


<a name="Theequationtosolve"></a><h3> The equation to solve </h3>

This program solves the Poisson equation
using the weak Galerkin finite element method:
@f{align*}{
  \nabla \cdot \left( -\mathbf{K} \nabla p \right)
    &= f,
    \qquad \mathbf{x} \in \Omega, \\
  p &=  p_D,\qquad \mathbf{x} \in \Gamma^D, \\
  \mathbf{u} \cdot \mathbf{n} &= u_N,
  \qquad \mathbf{x} \in \Gamma^N,
@f}
where $\Omega \subset \mathbb{R}^n (n=2,3)$ is a bounded domain.
In the context of the flow of a fluid through a porous medium,
$p$ is the pressure, $\mathbf{K}$ is a permeability tensor,
$f$ is the source term, and
$p_D, u_N$ represent Dirichlet and Neumann boundary conditions.
We can introduce a flux, $\mathbf{u} = -\mathbf{K} \nabla p$, that corresponds
to the Darcy velocity (in the way we did in step-20) and this variable will
be important in the considerations below.

In this program, we will consider a test case where the exact pressure
is $p = \sin \left( \pi x\right)\sin\left(\pi y \right)$ on the unit square domain,
with homogeneous Dirichelet boundary conditions and $\mathbf{K}$ the identity matrix.
Then we will calculate $L_2$ errors of pressure, velocity, and flux.


<a name="WeakGalerkinscheme"></a><h3> Weak Galerkin scheme </h3>


The Poisson equation above has a solution $p$ that needs to satisfy the weak
formulation of the problem,
@f{equation*}
\mathcal{A}\left(p,q \right) = \mathcal{F} \left(q \right),
@f}
for all test functions $q$, where
@f{equation*}
\mathcal{A}\left(p,q\right)
  \dealcoloneq \int_\Omega \left(\mathbf{K} \nabla p\right) \cdot \nabla q \;\mathrm{d}x,
@f}
and
@f{equation*}
\mathcal{F}\left(q\right)
  \dealcoloneq \int_\Omega f \, q \;\mathrm{d}x
  - \int_{\Gamma^N} u_N q \; \mathrm{d}x.
@f}
Here, we have integrated by parts in the bilinear form, and we are evaluating
the gradient of $p,p$ in the interior and the values of $q$ on the boundary
of the domain. All of this is well defined because we assume that the solution
is in $H^1$ for which taking the gradient and evaluating boundary values
are valid operations.

The idea of the weak Galerkin method is now to approximate the exact $p$
solution with a <i>discontinuous function</i> $p_h$. This function may only be
discontinuous along interfaces between cells, and because we will want to
evaluate this function also along interfaces, we have to
prescribe not only what values it is supposed to have in the cell interiors
but also its values along interfaces. We do this by saying that $p_h$ is
actually a tuple, $p_h=(p^\circ,p^\partial)$, though it's really just
a single function that is either equal to $p^\circ(x)$ or $p^\partial(x)$,
depending on whether it is evaluated at a point $x$ that lies in the cell
interior or on cell interfaces.

We would then like to simply stick this approximation into the bilinear
form above. This works for the case where we have to evaluate the
test function $q_h$ on the boundary (where we would simply take its interface
part $q_h^\partial$) but we have to be careful with the gradient because
that is only defined in cell interiors. Consequently,
the weak Galerkin scheme for the Poisson equation is defined by
@f{equation*}
\mathcal{A}_h\left(p_h,q \right) = \mathcal{F} \left(q_h \right),
@f}
for all discrete test functions $q_h$, where
@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,
@f}
and
@f{equation*}
\mathcal{F}\left(q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}} \int_K f \, q_h^\circ \;\mathrm{d}x
  - \sum_{\gamma \in \Gamma_h^N} \int_\gamma u_N q_h^\partial \;\mathrm{d}x,
@f}
The key point is that here, we have replaced the gradient $\nabla p_h$ by the
<i>discrete weak gradient</i> operator
$\nabla_{w,d} p_h$ that makes sense for our peculiarly defined approximation $p_h$.

The question is then how that operator works. For this, let us first say how we
think of the discrete approximation $p_h$ of the pressure. As mentioned above,
the "function" $p_h$ actually consists of two parts: the values $p_h^\circ$ in
the interior of cells, and $p_h^\partial$ on the interfaces. We have to define
discrete (finite-dimensional) function spaces for both of these; in this
program, we will use FE_DGQ for $p_h^\circ$ as the space in the interior of
cells (defined on each cell, but in general discontinuous along interfaces),
and FE_FaceQ for $p_h^\partial$ as the space on the interfaces.

Then let us consider just a single cell (because the integrals above are all
defined cell-wise, and because the weak discrete gradient is defined cell-by-cell).
The restriction of $p_h$ to cell $K$, $p_h|_K$ then consists
of the pair $(p_h^\circ|_K,p_h^\partial|_{\partial K})$. In essence, we can
think of $\nabla_{w,d} p_h$ of some function defined on $K$ that approximates
the gradient; in particular, if $p_h|_K$ was the restriction of a differentiable
function (to the interior and boundary of $K$ -- which would make it continuous
between the interior and boundary), then
$\nabla_{w,d} p_h$ would simply be the exact gradient $\nabla p_h$. But, since
$p_h|_K$ is not continuous between interior and boundary of $K$, we need a more
general definition; furthermore, we can not deal with arbitrary functions, and
so require that $\nabla_{w,d} p_h$ is also in a finite element space (which, since
the gradient is a vector, has to be vector-valued, and because the weak gradient
is defined on each cell separately, will also be discontinuous between cells).

The way this is done is to define this weak gradient operator $\nabla_{w,d}|_K :
DGQ_k(K) \times DGQ_r(\partial K) \rightarrow RT_s(K)$ (where $RT_s(K)$ is the
vector-valued Raviart-Thomas space of order $s$ on cell $K$) in the following way:
@f{equation*}{
  \int_K \mathbf v_h \cdot (\nabla_{w,d} p_h)
  =
  -\int_K (\nabla \cdot \mathbf v_h) p_h^\circ
  +\int_{\partial K} (\mathbf v_h \cdot \mathbf n) p_h^\partial,
@f}
for all test functions $\mathbf v_h \in RT_s(K)$.
This is, in essence, simply an application of the integration-by-parts
formula. In other words, for a given $p_h=(p^\circ_h,p^\partial_h)$,
we need to think of $\nabla_{w,d} p_h|_K$ as that
Raviart-Thomas function of degree $s$ for which the left hand side and right hand side
are equal for all test functions.

A key point to make is then the following: While the usual gradient $\nabla$ is
a *local* operator that computes derivatives based simply on the value of
a function at a point and its (infinitesimal) neighborhood, the weak discrete gradient
$\nabla_{w,d}$ does not have this property: It depends on the values of the function
it is applied to on the entire cell, including the cell's boundary. Both are,
however, linear operators as is clear from the definition of $\nabla_{w,d}$
above, and that will allow us to represent $\nabla_{w,d}$ via a matrix
in the discussion below.

@note It may be worth pointing out that while the weak discrete
  gradient is an element of the Raviart-Thomas space $RT_s(K)$ on each
  cell $K$, it is discontinuous between cells. On the other hand, the
  Raviart-Thomas space $RT_s=RT_s({\mathbb T})$ defined on the entire
  mesh and implemented by the FE_RaviartThomas class represents
  functions that have continuous normal components at interfaces
  between cells. This means that <i>globally</i>, $\nabla_{w,d} p_h$
  is not in $RT_s$, even though it is on every cell $K$ in $RT_s(K)$.
  Rather, it is in a "broken" Raviart-Thomas space that below we will
  represent by the symbol $DGRT_s$. (The term "broken" here refers to
  the process of "breaking something apart", and not to the synonym to
  the expression "not functional".) One might therefore (rightfully) argue that
  the notation used in the weak Galerkin literature is a bit misleading,
  but as so often it all depends on the context in which a certain
  notation is used -- in the current context, references to the
  Raviart-Thomas space or element are always understood to be to the
  "broken" spaces.

@note deal.II happens to have an implementation of this broken Raviart-Thomas
  space: The FE_DGRT class. As a consequence, in this tutorial we will simply
  always use the FE_DGRT class, even though in all of those places where
  we have to compute cell-local matrices and vectors, it makes no difference.


<a name="Representingtheweakgradient"></a><h3> Representing the weak gradient </h3>


Since $p_h$ is an element of a finite element space, we can expand it in a basis
as we always do, i.e., we can write
@f{equation*}{
  p_h(\mathbf x) = \sum_j P_j \varphi_j(\mathbf x).
@f}
Here, since $p_h$ has two components (the interior and the interface components),
the same must hold true for the basis functions $\varphi_j(\mathbf x)$, which we
can write as $\varphi_j = (\varphi_j^\circ,\varphi_j^\partial)$. If you've
followed the descriptions in step-8, step-20, and the
@ref vector_valued "documentation module on vector-valued problems",
it will be no surprise that for some values of $j$, $\varphi_j^\circ$ will be
zero, whereas for other values of $j$, $\varphi_j^\partial$ will be zero -- i.e.,
shape functions will be of either one or the other kind. That is not important,
here, however. What is important is that we need to wonder how we can represent
$\nabla_{w,d} \varphi_j$ because that is clearly what will appear in the
problem when we want to implement the bilinear form
@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  = \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,
@f}

The key point is that $\nabla_{w,d} \varphi_j$ is known to be a member of the
"broken" Raviart-Thomas space $DGRT_s$. What this means is that we can
represent (on each cell $K$ separately)
@f{equation*}
\nabla_{w,d} \varphi_j|_K
  = \sum_k C_{jk}^K \mathbf v_k|_K
@f}
where the functions $\mathbf v_k \in DGRT_s$, and where $C^K$ is a matrix of
dimension
@f{align*}{
 \text{dim}\left(DGQ_k(K) \times DGQ_r(K)\right) &\times \text{dim}\left(RT_s(K)\right)
  \\
 &=
 \left(\text{dim}(DGQ_k(K)) + \text{dim}(DGQ_r(K))\right) \times \text{dim}\left(RT_s(K)\right).
@f}
(That the weak discrete gradient can be represented as a matrix should not come
as a surprise: It is a linear operator from one finite dimensional
space to another finite dimensional space. If one chooses bases
for both of these spaces, then <i>every linear operator</i> can
of course be written as a matrix mapping the vector of expansion coefficients
with regards to the basis of the domain space of the operator, to
the vector of expansion coefficients with regards to the basis in the image
space.)

Using this expansion, we can easily use the definition of the weak
discrete gradient above to define what the matrix is going to be:
@f{equation*}{
  \int_K \mathbf v_i \cdot \left(\sum_k C_{jk}^K \mathbf v_k\right)
  =
  -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
  +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,
@f}
for all test functions $\mathbf v_i \in DGRT_s$.

This clearly leads to a linear system of the form
@f{equation*}{
  \sum_k M_{ik}^K C_{jk}^K
  =
  G_{ij}^K
@f}
with
@f{equation*}{
  M_{ik}^K = \int_K \mathbf v_i \cdot \mathbf v_k,
  \qquad\qquad
  G_{ij}^K = -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
             +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,
@f}
and consequently
@f{equation*}{
  \left(C^K\right)^T = \left(M^K\right)^{-1} G^K.
@f}
(In this last step, we have assumed that the indices $i,j,k$ only range
over those degrees of freedom active on cell $K$,
thereby ensuring that the mass matrix on the space $RT_s(K)$ is invertible.)
Equivalently, using the symmetry of the matrix $M$, we have that
@f{equation*}{
  C^K = \left(G^K\right)^{T} \left(M^K\right)^{-1}.
@f}
Also worth pointing out is that the
matrices $C^K$ and $G^K$ are of course not square but rectangular.


<a name="Assemblingthelinearsystem"></a><h3> Assembling the linear system </h3>


Having explained how the weak discrete gradient is defined, we can now
come back to the question of how the linear system for the equation in question
should be assembled. Specifically, using the definition of the bilinear
form ${\cal A}_h$ shown above, we then need to compute the elements of the
local contribution to the global matrix,
@f{equation*}{
  A^K_{ij} = \int_K \left({\mathbf K} \nabla_{w,d} \varphi_i\right) \cdot \nabla_{w,d} \varphi_j.
@f}
As explained above, we can expand $\nabla_{w,d} \varphi_i$ in terms of the
Raviart-Thomas basis on each cell, and similarly for $\nabla_{w,d} \varphi_j$:
@f{equation*}{
  A^K_{ij} = \int_K
    \left(
      {\mathbf K}
      \sum_k C_{ik}^K \mathbf v_k|_K
    \right)
    \cdot
    \sum_l C_{jl}^K \mathbf v_l|_K.
@f}
By re-arranging sums, this yields the following expression:
@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
     \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K.
@f}
So, if we have the matrix $C^K$ for each cell $K$, then we can easily compute
the contribution $A^K$ for cell $K$ to the matrix $A$ as follows:
@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
    H^K_{kl}
    =
    \sum_k \sum_l C_{ik}^K H^K_{kl} C_{jl}^K
    =
    \left(C^K H^K (C^K)^T \right)_{ij}.
@f}
Here,
@f{equation*}{
  H^K_{kl} =
  \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K,
@f}
which is really just the mass matrix on cell $K$ using the Raviart-Thomas
basis and weighting by the permeability tensor $\mathbf K$. The derivation
here then shows that the weak Galerkin method really just requires us
to compute these $C^K$ and $H^K$ matrices on each cell $K$, and then
$A^K = C^K H^K (C^K)^T$, which is easily computed. The code to be shown
below does exactly this.

Having so computed the contribution $A^K$ of cell $K$ to the global
matrix, all we have to do is to "distribute" these local contributions
into the global matrix. How this is done is first shown in step-3 and
step-4. In the current program, this will be facilitated by calling
AffineConstraints::distribute_local_to_global().

A linear system of course also needs a right hand side. There is no difficulty
associated with computing the right hand side here other than the fact
that we only need to use the cell-interior part $\varphi_i^\circ$ for
each shape function $\varphi_i$.


<a name="PostprocessingandiLsub2subierrors"></a><h3> Post-processing and <i>L<sub>2</sub></i>-errors </h3>


The discussions in the previous sections have given us a linear
system that we can solve for the numerical pressure $p_h$. We can use
this to compute an approximation to the variable $\mathbf u = -{\mathbf K}\nabla p$
that corresponds to the velocity with which the medium flows in a porous
medium if this is the model we are trying to solve. This kind of
step -- computing a derived quantity from the solution of the discrete
problem -- is typically called "post-processing".

Here, instead of using the exact gradient of $p_h$, let us instead
use the discrete weak gradient of $p_h$ to calculate the velocity on each element.
As discussed above,
on each element the gradient of the numerical pressure $\nabla p$ can be
approximated by discrete weak gradients  $ \nabla_{w,d}\phi_i$:
@f{equation*}
\nabla_{w,d} p_h
= \nabla_{w,d} \left(\sum_{i} P_i \phi_i\right)
= \sum_{i} P_i \nabla_{w,d}\phi_i.
@f}

On cell $K$,
the numerical velocity $ \mathbf{u}_h = -\mathbf{K} \nabla_{w,d}p_h$ can be written as
@f{align*}{
  \mathbf{u}_h
  &= -\mathbf{K} \nabla_{w,d} p_h
   = -\mathbf{K}\sum_{i} \sum_{j} P_i C^K_{ij}\mathbf{v}_j,
@f}
where $C^K$ is the expansion matrix from above, and
$\mathbf{v}_j$ is the basis function of the $RT$ space on a cell.

Unfortunately, $\mathbf{K} \mathbf{v}_j$ may not be in the $RT$ space
(unless, of course, if $\mathbf K$ is constant times the identity matrix).
So, in order to represent it in a finite element program, we need to
project it back into a finite dimensional space we can work with. Here,
we will use the $L_2$-projection to project it back to the (broken) $RT$
space.

We define the projection as
$ \mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right) =
\sum_{k} d_{jk}\mathbf{v}_k$ on each cell $K$.
For any $j$,
$\left( \mathbf{Q}_h \left( \mathbf{Kv}_j \right),\mathbf{v}_k \right)_K =
\left( \mathbf{Kv}_j,\mathbf{v}_k \right)_K.$
So, rather than the formula shown above, the numerical velocity on cell $K$
instead becomes
@f{equation*}
\mathbf{u}_h = \mathbf{Q}_h \left( -\mathbf{K}\nabla_{w,d}p_h \right) =
-\sum_i \sum_j P_i B^K_{ij}\mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right),
@f}
and we have the following system to solve for the coefficients $d_{jk}$:
@f{equation*}
 \sum_j
  \left(\mathbf{v}_i,\mathbf{v}_j\right)
   d_{jk}
   =
    \left( \mathbf{Kv}_j,\mathbf{v}_k \right).
@f}
In the implementation below, the matrix with elements
$
   d_{jk}
$
is called <code>cell_matrix_D</code>,
whereas the matrix with elements
$
      \left( \mathbf{Kv}_j,\mathbf{v}_k \right)
$
is called <code>cell_matrix_E</code>.

Then the elementwise velocity is
@f{equation*}
\mathbf{u}_h = -\sum_{i} \sum_{j}P_ic_{ij}\sum_{k}d_{jk}\mathbf{v}_k =
\sum_{k}- \left(\sum_{j} \sum_{i} P_ic_{ij}d_{jk} \right)\mathbf{v}_k,
@f}
where $-\sum_{j} \sum_{i} P_ic_{ij}d_{jk}$ is called
`cell_velocity` in the code.

Using this velocity obtained by "postprocessing" the solution, we can
define the $L_2$-errors of pressure, velocity, and flux
by the following formulas:
@f{align*}{
\|p-p_h^\circ\|^2
  &= \sum_{K \in \mathbb{T}} \|p-p_h^\circ\|_{L_2(K)}^2, \\
 \|\mathbf{u}-\mathbf{u}_h\|^2
  &= \sum_{K \in \mathbb{T}} \|\mathbf{u}-\mathbf{u}_h\|_{L_2(K)^2}^d,\\
\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|^2
  &= \sum_{K \in \mathbb{T}} \sum_{\gamma \subset \partial K}
    \frac{|K|}{|\gamma|} \|\mathbf{u} \cdot \mathbf{n} - \mathbf{u}_h \cdot \mathbf{n}\|_{L_2(\gamma)}^2,
@f}
where $| K |$ is the area of the element,
$\gamma$ are faces of the element,
$\mathbf{n}$ are unit normal vectors of each face. The last of these
norms measures the accuracy of the normal component of the velocity
vectors over the interfaces between the cells of the mesh. The scaling
factor $|K|/|\gamma|$ is chosen so as to scale out the difference in
the length (or area) of the collection of interfaces as the mesh size
changes.

The first of these errors above is easily computed using
VectorTools::integrate_difference. The others require a bit more work
and are implemented in the code below.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * This program is based on step-7, step-20 and step-51,
 * so most of the following header files are familiar. We
 * need the following, of which only the one that
 * imports the FE_DGRaviartThomas class (namely, `deal.II/fe/fe_dg_vector.h`)
 * is really new; the FE_DGRaviartThomas implements the "broken" Raviart-Thomas
 * space discussed in the introduction:
 * 
 * @code
 * #include <deal.II/base/quadrature.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/tensor_function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/point.h>
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_raviart_thomas.h>
 * #include <deal.II/fe/fe_dg_vector.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_face.h>
 * #include <deal.II/fe/component_mask.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/data_out_faces.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * Our first step, as always, is to put everything related to this tutorial
 * program into its own namespace:
 * 
 * @code
 * namespace Step61
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="TheWGDarcyEquationclasstemplate"></a> 
 * <h3>The WGDarcyEquation class template</h3>
 * 

 * 
 * This is the main class of this program. We will solve for the numerical
 * pressure in the interior and on faces using the weak Galerkin (WG) method,
 * and calculate the $L_2$ error of pressure. In the post-processing step, we
 * will also calculate $L_2$-errors of the velocity and flux.
 *   

 * 
 * The structure of the class is not fundamentally different from that of
 * previous tutorial programs, so there is little need to comment on the
 * details with one exception: The class has a member variable `fe_dgrt`
 * that corresponds to the "broken" Raviart-Thomas space mentioned in the
 * introduction. There is a matching `dof_handler_dgrt` that represents a
 * global enumeration of a finite element field created from this element, and
 * a vector `darcy_velocity` that holds nodal values for this field. We will
 * use these three variables after solving for the pressure to compute a
 * postprocessed velocity field for which we can then evaluate the error
 * and which we can output for visualization.
 * 
 * @code
 *   template <int dim>
 *   class WGDarcyEquation
 *   {
 *   public:
 *     WGDarcyEquation(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void compute_postprocessed_velocity();
 *     void compute_velocity_errors();
 *     void compute_pressure_error();
 *     void output_results() const;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     FESystem<dim>   fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     FE_DGRaviartThomas<dim> fe_dgrt;
 *     DoFHandler<dim>         dof_handler_dgrt;
 *     Vector<double>          darcy_velocity;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandexactsolution"></a> 
 * <h3>Right hand side, boundary values, and exact solution</h3>
 * 

 * 
 * Next, we define the coefficient matrix $\mathbf{K}$ (here, the
 * identity matrix), Dirichlet boundary conditions, the right-hand
 * side $f = 2\pi^2 \sin(\pi x) \sin(\pi y)$, and the exact solution
 * that corresponds to these choices for $K$ and $f$, namely $p =
 * \sin(\pi x) \sin(\pi y)$.
 * 
 * @code
 *   template <int dim>
 *   class Coefficient : public TensorFunction<2, dim>
 *   {
 *   public:
 *     Coefficient()
 *       : TensorFunction<2, dim>()
 *     {}
 * 
 *     virtual void value_list(const std::vector<Point<dim>> &points,
 *                             std::vector<Tensor<2, dim>> &values) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                     std::vector<Tensor<2, dim>> &  values) const
 *   {
 *     Assert(points.size() == values.size(),
 *            ExcDimensionMismatch(points.size(), values.size()));
 *     for (unsigned int p = 0; p < points.size(); ++p)
 *       values[p] = unit_symmetric_tensor<dim>();
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     BoundaryValues()
 *       : Function<dim>(2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/,
 *                                     const unsigned int /*component*/) const
 *   {
 *     return 0;
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
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> &p,
 *                                    const unsigned int /*component*/) const
 *   {
 *     return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) *
 *             std::sin(numbers::PI * p[1]));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The class that implements the exact pressure solution has an
 * oddity in that we implement it as a vector-valued one with two
 * components. (We say that it has two components in the constructor
 * where we call the constructor of the base Function class.) In the
 * `value()` function, we do not test for the value of the
 * `component` argument, which implies that we return the same value
 * for both components of the vector-valued function. We do this
 * because we describe the finite element in use in this program as
 * a vector-valued system that contains the interior and the
 * interface pressures, and when we compute errors, we will want to
 * use the same pressure solution to test both of these components.
 * 
 * @code
 *   template <int dim>
 *   class ExactPressure : public Function<dim>
 *   {
 *   public:
 *     ExactPressure()
 *       : Function<dim>(2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double ExactPressure<dim>::value(const Point<dim> &p,
 *                                    const unsigned int /*component*/) const
 *   {
 *     return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class ExactVelocity : public TensorFunction<1, dim>
 *   {
 *   public:
 *     ExactVelocity()
 *       : TensorFunction<1, dim>()
 *     {}
 * 
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const
 *   {
 *     Tensor<1, dim> return_value;
 *     return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) *
 *                       std::sin(numbers::PI * p[1]);
 *     return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) *
 *                       std::cos(numbers::PI * p[1]);
 *     return return_value;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationclassimplementation"></a> 
 * <h3>WGDarcyEquation class implementation</h3>
 * 

 * 
 * 
 * <a name="WGDarcyEquationWGDarcyEquation"></a> 
 * <h4>WGDarcyEquation::WGDarcyEquation</h4>
 * 

 * 
 * In this constructor, we create a finite element space for vector valued
 * functions, which will here include the ones used for the interior and
 * interface pressures, $p^\circ$ and $p^\partial$.
 * 
 * @code
 *   template <int dim>
 *   WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree)
 *     : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1)
 *     , dof_handler(triangulation)
 *     , fe_dgrt(degree)
 *     , dof_handler_dgrt(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationmake_grid"></a> 
 * <h4>WGDarcyEquation::make_grid</h4>
 * 

 * 
 * We generate a mesh on the unit square domain and refine it.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0, 1);
 *     triangulation.refine_global(5);
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "   Total number of cells: " << triangulation.n_cells()
 *               << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationsetup_system"></a> 
 * <h4>WGDarcyEquation::setup_system</h4>
 * 

 * 
 * After we have created the mesh above, we distribute degrees of
 * freedom and resize matrices and vectors. The only piece of
 * interest in this function is how we interpolate the boundary
 * values for the pressure. Since the pressure consists of interior
 * and interface components, we need to make sure that we only
 * interpolate onto that component of the vector-valued solution
 * space that corresponds to the interface pressures (as these are
 * the only ones that are defined on the boundary of the domain). We
 * do this via a component mask object for only the interface
 * pressures.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler_dgrt.distribute_dofs(fe_dgrt);
 * 
 *     std::cout << "   Number of pressure degrees of freedom: "
 *               << dof_handler.n_dofs() << std::endl;
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 * 
 *     {
 *       constraints.clear();
 *       const FEValuesExtractors::Scalar interface_pressure(1);
 *       const ComponentMask              interface_pressure_mask =
 *         fe.component_mask(interface_pressure);
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                interface_pressure_mask);
 *       constraints.close();
 *     }
 * 
 * 
 * @endcode
 * 
 * In the bilinear form, there is no integration term over faces
 * between two neighboring cells, so we can just use
 * <code>DoFTools::make_sparsity_pattern</code> to calculate the sparse
 * matrix.
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationassemble_system"></a> 
 * <h4>WGDarcyEquation::assemble_system</h4>
 * 

 * 
 * This function is more interesting. As detailed in the
 * introduction, the assembly of the linear system requires us to
 * evaluate the weak gradient of the shape functions, which is an
 * element in the Raviart-Thomas space. As a consequence, we need to
 * define a Raviart-Thomas finite element object, and have FEValues
 * objects that evaluate it at quadrature points. We then need to
 * compute the matrix $C^K$ on every cell $K$, for which we need the
 * matrices $M^K$ and $G^K$ mentioned in the introduction.
 *   

 * 
 * A point that may not be obvious is that in all previous tutorial
 * programs, we have always called FEValues::reinit() with a cell
 * iterator from a DoFHandler. This is so that one can call
 * functions such as FEValuesBase::get_function_values() that
 * extract the values of a finite element function (represented by a
 * vector of DoF values) on the quadrature points of a cell. For
 * this operation to work, one needs to know which vector elements
 * correspond to the degrees of freedom on a given cell -- i.e.,
 * exactly the kind of information and operation provided by the
 * DoFHandler class.
 *   

 * 
 * We could create a DoFHandler object for the "broken" Raviart-Thomas space
 * (using the FE_DGRT class), but we really don't want to here: At
 * least in the current function, we have no need for any globally defined
 * degrees of freedom associated with this broken space, but really only
 * need to reference the shape functions of such a space on the current
 * cell. As a consequence, we use the fact that one can call
 * FEValues::reinit() also with cell iterators into Triangulation
 * objects (rather than DoFHandler objects). In this case, FEValues
 * can of course only provide us with information that only
 * references information about cells, rather than degrees of freedom
 * enumerated on these cells. So we can't use
 * FEValuesBase::get_function_values(), but we can use
 * FEValues::shape_value() to obtain the values of shape functions
 * at quadrature points on the current cell. It is this kind of
 * functionality we will make use of below. The variable that will
 * give us this information about the Raviart-Thomas functions below
 * is then the `fe_values_rt` (and corresponding `fe_face_values_rt`)
 * object.
 *   

 * 
 * Given this introduction, the following declarations should be
 * pretty obvious:
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::assemble_system()
 *   {
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     FEValues<dim>     fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell();
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size();
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();
 * 
 *     RightHandSide<dim>  right_hand_side;
 *     std::vector<double> right_hand_side_values(n_q_points);
 * 
 *     const Coefficient<dim>      coefficient;
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * 
 * @endcode
 * 
 * Next, let us declare the various cell matrices discussed in the
 * introduction:
 * 
 * @code
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell);
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt);
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 *     Vector<double>     cell_solution(dofs_per_cell);
 * 
 * @endcode
 * 
 * We need <code>FEValuesExtractors</code> to access the @p interior and
 * @p face component of the shape functions.
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure_interior(0);
 *     const FEValuesExtractors::Scalar pressure_face(1);
 * 
 * @endcode
 * 
 * This finally gets us in position to loop over all cells. On
 * each cell, we will first calculate the various cell matrices
 * used to construct the local matrix -- as they depend on the
 * cell in question, they need to be re-computed on each cell. We
 * need shape functions for the Raviart-Thomas space as well, for
 * which we need to create first an iterator to the cell of the
 * triangulation, which we can obtain by assignment from the cell
 * pointing into the DoFHandler.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 * 
 *         const typename Triangulation<dim>::active_cell_iterator cell_dgrt =
 *           cell;
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                    right_hand_side_values);
 *         coefficient.value_list(fe_values.get_quadrature_points(),
 *                                coefficient_values);
 * 
 * @endcode
 * 
 * The first cell matrix we will compute is the mass matrix
 * for the Raviart-Thomas space.  Hence, we need to loop over
 * all the quadrature points for the velocity FEValues object.
 * 
 * @code
 *         cell_matrix_M = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q);
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *                 {
 *                   const Tensor<1, dim> v_k =
 *                     fe_values_dgrt[velocities].value(k, q);
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q));
 *                 }
 *             }
 * @endcode
 * 
 * Next we take the inverse of this matrix by using
 * FullMatrix::gauss_jordan(). It will be used to calculate
 * the coefficient matrix $C^K$ later. It is worth recalling
 * later that `cell_matrix_M` actually contains the *inverse*
 * of $M^K$ after this call.
 * 
 * @code
 *         cell_matrix_M.gauss_jordan();
 * 
 * @endcode
 * 
 * From the introduction, we know that the right hand side
 * $G^K$ of the equation that defines $C^K$ is the difference
 * between a face integral and a cell integral. Here, we
 * approximate the negative of the contribution in the
 * interior. Each component of this matrix is the integral of
 * a product between a basis function of the polynomial space
 * and the divergence of a basis function of the
 * Raviart-Thomas space. These basis functions are defined in
 * the interior.
 * 
 * @code
 *         cell_matrix_G = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const double div_v_i =
 *                 fe_values_dgrt[velocities].divergence(i, q);
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const double phi_j_interior =
 *                     fe_values[pressure_interior].value(j, q);
 * 
 *                   cell_matrix_G(i, j) -=
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q));
 *                 }
 *             }
 * 
 * 
 * @endcode
 * 
 * Next, we approximate the integral on faces by quadrature.
 * Each component is the integral of a product between a basis function
 * of the polynomial space and the dot product of a basis function of
 * the Raviart-Thomas space and the normal vector. So we loop over all
 * the faces of the element and obtain the normal vector.
 * 
 * @code
 *         for (const auto &face : cell->face_iterators())
 *           {
 *             fe_face_values.reinit(cell, face);
 *             fe_face_values_dgrt.reinit(cell_dgrt, face);
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 const Tensor<1, dim> normal = fe_face_values.normal_vector(q);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *                   {
 *                     const Tensor<1, dim> v_i =
 *                       fe_face_values_dgrt[velocities].value(i, q);
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                       {
 *                         const double phi_j_face =
 *                           fe_face_values[pressure_face].value(j, q);
 * 
 *                         cell_matrix_G(i, j) +=
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
 *                       }
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * @p cell_matrix_C is then the matrix product between the
 * transpose of $G^K$ and the inverse of the mass matrix
 * (where this inverse is stored in @p cell_matrix_M):
 * 
 * @code
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);
 * 
 * @endcode
 * 
 * Finally we can compute the local matrix $A^K$.  Element
 * $A^K_{ij}$ is given by $\int_{E} \sum_{k,l} C_{ik} C_{jl}
 * (\mathbf{K} \mathbf{v}_k) \cdot \mathbf{v}_l
 * \mathrm{d}x$. We have calculated the coefficients $C$ in
 * the previous step, and so obtain the following after
 * suitably re-arranging the loops:
 * 
 * @code
 *         local_matrix = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *               {
 *                 const Tensor<1, dim> v_k =
 *                   fe_values_dgrt[velocities].value(k, q);
 *                 for (unsigned int l = 0; l < dofs_per_cell_dgrt; ++l)
 *                   {
 *                     const Tensor<1, dim> v_l =
 *                       fe_values_dgrt[velocities].value(l, q);
 * 
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                       for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                         local_matrix(i, j) +=
 *                           (coefficient_values[q] * cell_matrix_C[i][k] * v_k) *
 *                           cell_matrix_C[j][l] * v_l * fe_values_dgrt.JxW(q);
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * Next, we calculate the right hand side, $\int_{K} f q \mathrm{d}x$:
 * 
 * @code
 *         cell_rhs = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               cell_rhs(i) += (fe_values[pressure_interior].value(i, q) *
 *                               right_hand_side_values[q] * fe_values.JxW(q));
 *             }
 * 
 * @endcode
 * 
 * The last step is to distribute components of the local
 * matrix into the system matrix and transfer components of
 * the cell right hand side into the system right hand side:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(
 *           local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimsolve"></a> 
 * <h4>WGDarcyEquation<dim>::solve</h4>
 * 

 * 
 * This step is rather trivial and the same as in many previous
 * tutorial programs:
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> solver(solver_control);
 *     solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_postprocessed_velocity"></a> 
 * <h4>WGDarcyEquation<dim>::compute_postprocessed_velocity</h4>
 * 

 * 
 * In this function, compute the velocity field from the pressure
 * solution previously computed. The
 * velocity is defined as $\mathbf{u}_h = \mathbf{Q}_h \left(
 * -\mathbf{K}\nabla_{w,d}p_h \right)$, which requires us to compute
 * many of the same terms as in the assembly of the system matrix.
 * There are also the matrices $E^K,D^K$ we need to assemble (see
 * the introduction) but they really just follow the same kind of
 * pattern.
 *   

 * 
 * Computing the same matrices here as we have already done in the
 * `assemble_system()` function is of course wasteful in terms of
 * CPU time. Likewise, we copy some of the code from there to this
 * function, and this is also generally a poor idea. A better
 * implementation might provide for a function that encapsulates
 * this duplicated code. One could also think of using the classic
 * trade-off between computing efficiency and memory efficiency to
 * only compute the $C^K$ matrices once per cell during the
 * assembly, storing them somewhere on the side, and re-using them
 * here. (This is what step-51 does, for example, where the
 * `assemble_system()` function takes an argument that determines
 * whether the local matrices are recomputed, and a similar approach
 * -- maybe with storing local matrices elsewhere -- could be
 * adapted for the current program.)
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_postprocessed_velocity()
 *   {
 *     darcy_velocity.reinit(dof_handler_dgrt.n_dofs());
 * 
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell();
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size();
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();
 * 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices_dgrt(
 *       dofs_per_cell_dgrt);
 * 
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell);
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_D(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_E(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 * 
 *     Vector<double> cell_solution(dofs_per_cell);
 *     Vector<double> cell_velocity(dofs_per_cell_dgrt);
 * 
 *     const Coefficient<dim>      coefficient;
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points_dgrt);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure_interior(0);
 *     const FEValuesExtractors::Scalar pressure_face(1);
 * 
 * @endcode
 * 
 * In the introduction, we explained how to calculate the numerical velocity
 * on the cell. We need the pressure solution values on each cell,
 * coefficients of the Gram matrix and coefficients of the $L_2$ projection.
 * We have already calculated the global solution, so we will extract the
 * cell solution from the global solution. The coefficients of the Gram
 * matrix have been calculated when we assembled the system matrix for the
 * pressures. We will do the same way here. For the coefficients of the
 * projection, we do matrix multiplication, i.e., the inverse of the Gram
 * matrix times the matrix with $(\mathbf{K} \mathbf{w}, \mathbf{w})$ as
 * components. Then, we multiply all these coefficients and call them beta.
 * The numerical velocity is the product of beta and the basis functions of
 * the Raviart-Thomas space.
 * 
 * @code
 *     typename DoFHandler<dim>::active_cell_iterator
 *       cell = dof_handler.begin_active(),
 *       endc = dof_handler.end(), cell_dgrt = dof_handler_dgrt.begin_active();
 *     for (; cell != endc; ++cell, ++cell_dgrt)
 *       {
 *         fe_values.reinit(cell);
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 *         coefficient.value_list(fe_values_dgrt.get_quadrature_points(),
 *                                coefficient_values);
 * 
 * @endcode
 * 
 * The component of this <code>cell_matrix_E</code> is the integral of
 * $(\mathbf{K} \mathbf{w}, \mathbf{w})$. <code>cell_matrix_M</code> is
 * the Gram matrix.
 * 
 * @code
 *         cell_matrix_M = 0;
 *         cell_matrix_E = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q);
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *                 {
 *                   const Tensor<1, dim> v_k =
 *                     fe_values_dgrt[velocities].value(k, q);
 * 
 *                   cell_matrix_E(i, k) +=
 *                     (coefficient_values[q] * v_i * v_k * fe_values_dgrt.JxW(q));
 * 
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q));
 *                 }
 *             }
 * 
 * @endcode
 * 
 * To compute the matrix $D$ mentioned in the introduction, we
 * then need to evaluate $D=M^{-1}E$ as explained in the
 * introduction:
 * 
 * @code
 *         cell_matrix_M.gauss_jordan();
 *         cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E);
 * 
 * @endcode
 * 
 * Then we also need, again, to compute the matrix $C$ that is
 * used to evaluate the weak discrete gradient. This is the
 * exact same code as used in the assembly of the system
 * matrix, so we just copy it from there:
 * 
 * @code
 *         cell_matrix_G = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const double div_v_i =
 *                 fe_values_dgrt[velocities].divergence(i, q);
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const double phi_j_interior =
 *                     fe_values[pressure_interior].value(j, q);
 * 
 *                   cell_matrix_G(i, j) -=
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q));
 *                 }
 *             }
 * 
 *         for (const auto &face : cell->face_iterators())
 *           {
 *             fe_face_values.reinit(cell, face);
 *             fe_face_values_dgrt.reinit(cell_dgrt, face);
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 const Tensor<1, dim> normal = fe_face_values.normal_vector(q);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *                   {
 *                     const Tensor<1, dim> v_i =
 *                       fe_face_values_dgrt[velocities].value(i, q);
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                       {
 *                         const double phi_j_face =
 *                           fe_face_values[pressure_face].value(j, q);
 * 
 *                         cell_matrix_G(i, j) +=
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
 *                       }
 *                   }
 *               }
 *           }
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);
 * 
 * @endcode
 * 
 * Finally, we need to extract the pressure unknowns that
 * correspond to the current cell:
 * 
 * @code
 *         cell->get_dof_values(solution, cell_solution);
 * 
 * @endcode
 * 
 * We are now in a position to compute the local velocity
 * unknowns (with respect to the Raviart-Thomas space we are
 * projecting the term $-\mathbf K \nabla_{w,d} p_h$ into):
 * 
 * @code
 *         cell_velocity = 0;
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_velocity(k) +=
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j));
 * 
 * @endcode
 * 
 * We compute Darcy velocity.
 * This is same as cell_velocity but used to graph Darcy velocity.
 * 
 * @code
 *         cell_dgrt->get_dof_indices(local_dof_indices_dgrt);
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               darcy_velocity(local_dof_indices_dgrt[k]) +=
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j));
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_pressure_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_pressure_error</h4>
 * 

 * 
 * This part is to calculate the $L_2$ error of the pressure.  We
 * define a vector that holds the norm of the error on each cell.
 * Next, we use VectorTool::integrate_difference() to compute the
 * error in the $L_2$ norm on each cell. However, we really only
 * care about the error in the interior component of the solution
 * vector (we can't even evaluate the interface pressures at the
 * quadrature points because these are all located in the interior
 * of cells) and consequently have to use a weight function that
 * ensures that the interface component of the solution variable is
 * ignored. This is done by using the ComponentSelectFunction whose
 * arguments indicate which component we want to select (component
 * zero, i.e., the interior pressures) and how many components there
 * are in total (two).
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_pressure_error()
 *   {
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     const ComponentSelectFunction<dim> select_interior_pressure(0, 2);
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       ExactPressure<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe.degree + 2),
 *                                       VectorTools::L2_norm,
 *                                       &select_interior_pressure);
 * 
 *     const double L2_error = difference_per_cell.l2_norm();
 *     std::cout << "L2_error_pressure " << L2_error << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_velocity_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_velocity_error</h4>
 * 

 * 
 * In this function, we evaluate $L_2$ errors for the velocity on
 * each cell, and $L_2$ errors for the flux on faces. The function
 * relies on the `compute_postprocessed_velocity()` function having
 * previous computed, which computes the velocity field based on the
 * pressure solution that has previously been computed.
 *   

 * 
 * We are going to evaluate velocities on each cell and calculate
 * the difference between numerical and exact velocities.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_velocity_errors()
 *   {
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 *     const unsigned int n_face_q_points_dgrt =
 *       fe_face_values_dgrt.get_quadrature().size();
 * 
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points_dgrt);
 *     std::vector<Tensor<1, dim>> velocity_face_values(n_face_q_points_dgrt);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     const ExactVelocity<dim> exact_velocity;
 * 
 *     double L2_err_velocity_cell_sqr_global = 0;
 *     double L2_err_flux_sqr                 = 0;
 * 
 * @endcode
 * 
 * Having previously computed the postprocessed velocity, we here
 * only have to extract the corresponding values on each cell and
 * face and compare it to the exact values.
 * 
 * @code
 *     for (const auto &cell_dgrt : dof_handler_dgrt.active_cell_iterators())
 *       {
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 * @endcode
 * 
 * First compute the $L_2$ error between the postprocessed velocity
 * field and the exact one:
 * 
 * @code
 *         fe_values_dgrt[velocities].get_function_values(darcy_velocity,
 *                                                        velocity_values);
 *         double L2_err_velocity_cell_sqr_local = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           {
 *             const Tensor<1, dim> velocity = velocity_values[q];
 *             const Tensor<1, dim> true_velocity =
 *               exact_velocity.value(fe_values_dgrt.quadrature_point(q));
 * 
 *             L2_err_velocity_cell_sqr_local +=
 *               ((velocity - true_velocity) * (velocity - true_velocity) *
 *                fe_values_dgrt.JxW(q));
 *           }
 *         L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local;
 * 
 * @endcode
 * 
 * For reconstructing the flux we need the size of cells and
 * faces. Since fluxes are calculated on faces, we have the
 * loop over all four faces of each cell. To calculate the
 * face velocity, we extract values at the quadrature points from the
 * `darcy_velocity` which we have computed previously. Then, we
 * calculate the squared velocity error in normal direction. Finally, we
 * calculate the $L_2$ flux error on the cell by appropriately scaling
 * with face and cell areas and add it to the global error.
 * 
 * @code
 *         const double cell_area = cell_dgrt->measure();
 *         for (const auto &face_dgrt : cell_dgrt->face_iterators())
 *           {
 *             const double face_length = face_dgrt->measure();
 *             fe_face_values_dgrt.reinit(cell_dgrt, face_dgrt);
 *             fe_face_values_dgrt[velocities].get_function_values(
 *               darcy_velocity, velocity_face_values);
 * 
 *             double L2_err_flux_face_sqr_local = 0;
 *             for (unsigned int q = 0; q < n_face_q_points_dgrt; ++q)
 *               {
 *                 const Tensor<1, dim> velocity = velocity_face_values[q];
 *                 const Tensor<1, dim> true_velocity =
 *                   exact_velocity.value(fe_face_values_dgrt.quadrature_point(q));
 * 
 *                 const Tensor<1, dim> normal =
 *                   fe_face_values_dgrt.normal_vector(q);
 * 
 *                 L2_err_flux_face_sqr_local +=
 *                   ((velocity * normal - true_velocity * normal) *
 *                    (velocity * normal - true_velocity * normal) *
 *                    fe_face_values_dgrt.JxW(q));
 *               }
 *             const double err_flux_each_face =
 *               L2_err_flux_face_sqr_local / face_length * cell_area;
 *             L2_err_flux_sqr += err_flux_each_face;
 *           }
 *       }
 * 
 * @endcode
 * 
 * After adding up errors over all cells and faces, we take the
 * square root and get the $L_2$ errors of velocity and
 * flux. These we output to screen.
 * 
 * @code
 *     const double L2_err_velocity_cell =
 *       std::sqrt(L2_err_velocity_cell_sqr_global);
 *     const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr);
 * 
 *     std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl
 *               << "L2_error_flux: " << L2_err_flux_face << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationoutput_results"></a> 
 * <h4>WGDarcyEquation::output_results</h4>
 * 

 * 
 * We have two sets of results to output: the interior solution and
 * the skeleton solution. We use <code>DataOut</code> to visualize
 * interior results. The graphical output for the skeleton results
 * is done by using the DataOutFaces class.
 *   

 * 
 * In both of the output files, both the interior and the face
 * variables are stored. For the interface output, the output file
 * simply contains the interpolation of the interior pressures onto
 * the faces, but because it is undefined which of the two interior
 * pressure variables you get from the two adjacent cells, it is
 * best to ignore the interior pressure in the interface output
 * file. Conversely, for the cell interior output file, it is of
 * course impossible to show any interface pressures $p^\partial$,
 * because these are only available on interfaces and not cell
 * interiors. Consequently, you will see them shown as an invalid
 * value (such as an infinity).
 *   

 * 
 * For the cell interior output, we also want to output the velocity
 * variables. This is a bit tricky since it lives on the same mesh
 * but uses a different DoFHandler object (the pressure variables live
 * on the `dof_handler` object, the Darcy velocity on the `dof_handler_dgrt`
 * object). Fortunately, there are variations of the
 * DataOut::add_data_vector() function that allow specifying which
 * DoFHandler a vector corresponds to, and consequently we can visualize
 * the data from both DoFHandler objects within the same file.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::output_results() const
 *   {
 *     {
 *       DataOut<dim> data_out;
 * 
 * @endcode
 * 
 * First attach the pressure solution to the DataOut object:
 * 
 * @code
 *       const std::vector<std::string> solution_names = {"interior_pressure",
 *                                                        "interface_pressure"};
 *       data_out.add_data_vector(dof_handler, solution, solution_names);
 * 
 * @endcode
 * 
 * Then do the same with the Darcy velocity field, and continue
 * with writing everything out into a file.
 * 
 * @code
 *       const std::vector<std::string> velocity_names(dim, "velocity");
 *       const std::vector<
 *         DataComponentInterpretation::DataComponentInterpretation>
 *         velocity_component_interpretation(
 *           dim, DataComponentInterpretation::component_is_part_of_vector);
 *       data_out.add_data_vector(dof_handler_dgrt,
 *                                darcy_velocity,
 *                                velocity_names,
 *                                velocity_component_interpretation);
 * 
 *       data_out.build_patches(fe.degree);
 *       std::ofstream output("solution_interior.vtu");
 *       data_out.write_vtu(output);
 *     }
 * 
 *     {
 *       DataOutFaces<dim> data_out_faces(false);
 *       data_out_faces.attach_dof_handler(dof_handler);
 *       data_out_faces.add_data_vector(solution, "Pressure_Face");
 *       data_out_faces.build_patches(fe.degree);
 *       std::ofstream face_output("solution_interface.vtu");
 *       data_out_faces.write_vtu(face_output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationrun"></a> 
 * <h4>WGDarcyEquation::run</h4>
 * 

 * 
 * This is the final function of the main class. It calls the other functions
 * of our class.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::run()
 *   {
 *     std::cout << "Solving problem in " << dim << " space dimensions."
 *               << std::endl;
 *     make_grid();
 *     setup_system();
 *     assemble_system();
 *     solve();
 *     compute_postprocessed_velocity();
 *     compute_pressure_error();
 *     compute_velocity_errors();
 *     output_results();
 *   }
 * 
 * } // namespace Step61
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * This is the main function. We can change the dimension here to run in 3d.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step61::WGDarcyEquation<2> wg_darcy(0);
 *       wg_darcy.run();
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


We run the program with a right hand side that will produce the
solution $p = \sin(\pi x) \sin(\pi y)$ and with homogeneous Dirichlet
boundary conditions in the domain $\Omega = (0,1)^2$. In addition, we
choose the coefficient matrix in the differential operator
$\mathbf{K}$ as the identity matrix. We test this setup using
$\mbox{WG}(Q_0,Q_0;RT_{[0]})$, $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ and
$\mbox{WG}(Q_2,Q_2;RT_{[2]})$ element combinations, which one can
select by using the appropriate constructor argument for the
`WGDarcyEquation` object in `main()`. We will then visualize pressure
values in interiors of cells and on faces. We want to see that the
pressure maximum is around 1 and the minimum is around 0. With mesh
refinement, the convergence rates of pressure, velocity and flux
should then be around 1 for $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ , 2 for
$\mbox{WG}(Q_1,Q_1;RT_{[1]})$, and 3 for
$\mbox{WG}(Q_2,Q_2;RT_{[2]})$.


<a name="TestresultsoniWGQsub0subQsub0subRTsub0subi"></a><h3>Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></h3>


The following figures show interior pressures and face pressures using the
$\mbox{WG}(Q_0,Q_0;RT_{[0]})$ element. The mesh is refined 2 times (top)
and 4 times (bottom), respectively. (This number can be adjusted in the
`make_grid()` function.) When the mesh is coarse, one can see
the face pressures $p^\partial$ neatly between the values of the interior
pressures $p^\circ$ on the two adjacent cells.

<table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_2.png" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_4.png" alt=""></td>
  </tr>
</table>

From the figures, we can see that with the mesh refinement, the maximum and
minimum pressure values are approaching the values we expect.
Since the mesh is a rectangular mesh and numbers of cells in each direction is even, we
have symmetric solutions. From the 3d figures on the right,
we can see that on $\mbox{WG}(Q_0,Q_0;RT_{[0]})$, the pressure is a constant
in the interior of the cell, as expected.

<a name="Convergencetableforik0i"></a><h4>Convergence table for <i>k=0</i></h4>


We run the code with differently refined meshes (chosen in the `make_grid()` function)
and get the following convergence rates of pressure,
velocity, and flux (as defined in the introduction).

<table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    1.587e-01        </td><td>        5.113e-01               </td><td>   7.062e-01 </td>
  </tr>
  <tr>
   <td>   3                  </td><td>    8.000e-02        </td><td>        2.529e-01               </td><td>   3.554e-01 </td>
  </tr>
  <tr>
   <td>   4                  </td><td>    4.006e-02        </td><td>        1.260e-01               </td><td>   1.780e-01 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>    2.004e-02        </td><td>        6.297e-02               </td><td>   8.902e-02 </td>
  </tr>
  <tr>
   <th>Conv.rate             </th><th>      1.00           </th><th>          1.00                  </th><th>      1.00   </th>
  </tr>
</table>

We can see that the convergence rates of $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ are around 1.
This, of course, matches our theoretical expectations.


<a name="TestresultsoniWGQsub1subQsub1subRTsub1subi"></a><h3>Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></h3>


We can repeat the experiment from above using the next higher polynomial
degree:
The following figures are interior pressures and face pressures implemented using
$\mbox{WG}(Q_1,Q_1;RT_{[1]})$. The mesh is refined 4 times.  Compared to the
previous figures using
$\mbox{WG}(Q_0,Q_0;RT_{[0]})$, on each cell, the solution is no longer constant
on each cell, as we now use bilinear polynomials to do the approximation.
Consequently, there are 4 pressure values in one interior, 2 pressure values on
each face.

<table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_3d_4.png" alt=""></td>
  </tr>
</table>

Compared to the corresponding image for the $\mbox{WG}(Q_0,Q_0;RT_{[0]})$
combination, the solution is now substantially more accurate and, in
particular so close to being continuous at the interfaces that we can
no longer distinguish the interface pressures $p^\partial$ from the
interior pressures $p^\circ$ on the adjacent cells.

<a name="Convergencetableforik1i"></a><h4>Convergence table for <i>k=1</i></h4>


The following are the convergence rates of pressure, velocity, and flux
we obtain from using the $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ element combination:

<table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
    <td>  2           </td><td>           1.613e-02      </td><td>          5.093e-02     </td><td>             7.167e-02   </td>
  </tr>
  <tr>
    <td>  3           </td><td>           4.056e-03      </td><td>          1.276e-02     </td><td>             1.802e-02    </td>
  </tr>
  <tr>
    <td>  4           </td><td>           1.015e-03      </td><td>          3.191e-03     </td><td>             4.512e-03  </td>
  </tr>
  <tr>
    <td>  5           </td><td>           2.540e-04      </td><td>          7.979e-04     </td><td>             1.128e-03  </td>
  </tr>
  <tr>
    <th>Conv.rate     </th><th>              2.00        </th><th>             2.00       </th><th>                 2.00    </th>
  </tr>
</table>

The convergence rates of $WG(Q_1,Q_1;RT_{[1]})$ are around 2, as expected.



<a name="TestresultsoniWGQsub2subQsub2subRTsub2subi"></a><h3>Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></h3>


Let us go one polynomial degree higher.
The following are interior pressures and face pressures implemented using
$WG(Q_2,Q_2;RT_{[2]})$, with mesh size $h = 1/32$ (i.e., 5 global mesh
refinement steps). In the program, we use
`data_out_face.build_patches(fe.degree)` when generating graphical output
(see the documentation of DataOut::build_patches()), which here implies that
we divide each 2d cell interior into 4 subcells in order to provide a better
visualization of the quadratic polynomials.
<table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_2d_5.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_3d_5.png" alt=""></td>
  </tr>
</table>


<a name="Convergencetableforik2i"></a><h4>Convergence table for <i>k=2</i></h4>


As before, we can generate convergence data for the
$L_2$ errors of pressure, velocity, and flux
using the $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ combination:

<table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
     <td>  2               </td><td>       1.072e-03       </td><td>         3.375e-03       </td><td>           4.762e-03   </td>
  </tr>
  <tr>
    <td>   3               </td><td>       1.347e-04       </td><td>         4.233e-04       </td><td>           5.982e-04    </td>
  </tr>
  <tr>
    <td>   4               </td><td>       1.685e-05      </td><td>          5.295e-05       </td><td>           7.487e-05  </td>
  </tr>
  <tr>
    <td>   5               </td><td>       2.107e-06      </td><td>          6.620e-06       </td><td>           9.362e-06  </td>
  </tr>
  <tr>
    <th>Conv.rate          </th><th>         3.00         </th><th>            3.00          </th><th>              3.00    </th>
  </tr>
</table>

Once more, the convergence rates of $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ is
as expected, with values around 3.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-61.cc"
*/
