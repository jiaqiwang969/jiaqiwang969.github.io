/**
@page step_60 The step-60 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#NonmatchinggridconstraintsthroughdistributedLagrangemultipliers">Non-matching grid constraints through distributed Lagrange multipliers</a>
        <li><a href="#Thetestcase">The testcase</a>
        <li><a href="#References">References</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#DistributedLagrangeProblem">DistributedLagrangeProblem</a>
        <li><a href="#DistributedLagrangeProblemParameters">DistributedLagrangeProblem::Parameters</a>
        <li><a href="#Setup">Set up</a>
        <li><a href="#Assembly">Assembly</a>
        <li><a href="#Solve">Solve</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Testcase1"> Test case 1: </a>
        <li><a href="#Testcase2and3"> Test case 2 and 3: </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Runningwithspacedimequaltothree"> Running with `spacedim` equal to three</a>
        <li><a href="#Moregeneraldomains"> More general domains </a>
        <li><a href="#Preconditioner"> Preconditioner</a>
        <li><a href="#ParallelCode"> Parallel Code </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Luca Heltai and Giovanni Alzetta, SISSA, Trieste.
</i>

@dealiiTutorialDOI{10.5281/zenodo.1243280,https://zenodo.org/badge/DOI/10.5281/zenodo.1243280.svg}


<a name="Introduction"></a><h1>Introduction</h1>


<a name="NonmatchinggridconstraintsthroughdistributedLagrangemultipliers"></a><h3>Non-matching grid constraints through distributed Lagrange multipliers</h3>



In this tutorial we consider the case of two domains, $\Omega$ in
$R^{\text{spacedim}}$ and $\Gamma$ in $R^{\text{dim}}$, where $\Gamma$ is
embedded in $\Omega$ ($\Gamma \subseteq \Omega$). We want to solve a partial
differential equation on $\Omega$, enforcing some conditions on the solution of
the problem *on the embedded domain* $\Gamma$.

There are two interesting scenarios:

- the geometrical dimension `dim` of the embedded domain $\Gamma$ is the same of
the domain $\Omega$ (`spacedim`), that is, the spacedim-dimensional measure of
$\Gamma$ is not zero, or

- the embedded domain $\Gamma$ has an intrinsic dimension `dim` which is smaller
than that of $\Omega$ (`spacedim`), thus its spacedim-dimensional measure is
zero; for example it is a curve embedded in a two dimensional domain, or a
surface embedded in a three-dimensional domain.

In both cases define the restriction operator $\gamma$ as the operator that,
given a continuous function on $\Omega$, returns its (continuous) restriction on
$\Gamma$, i.e.,

\f[
\gamma : C^0(\Omega) \mapsto C^0(\Gamma), \quad \text{ s.t. } \gamma u = u|_{\Gamma} \in C^0(\Gamma),
\quad \forall u \in C^0(\Omega).
\f]

It is well known that the operator $\gamma$ can be extended to a continuous
operator on $H^1(\Omega)$, mapping functions in $H^1(\Omega)$ to functions in
$H^1(\Gamma)$ when the intrinsic dimension of $\Gamma$ is the same of $\Omega$.

The same is true, with a less regular range space (namely $H^{1/2}(\Gamma)$),
when the dimension of $\Gamma$ is one less with respect to $\Omega$, and
$\Gamma$ does not have a boundary. In this second case, the operator $\gamma$ is
also known as the *trace* operator, and it is well defined for Lipschitz
co-dimension one curves and surfaces $\Gamma$ embedded in $\Omega$ (read  <a
href="https://en.wikipedia.org/wiki/Trace_operator">this wikipedia article</a>
for further details on the trace operator).

The co-dimension two case is a little more complicated, and in general it is not
possible to construct a continuous trace operator, not even from $H^1(\Omega)$ to
$L^2(\Gamma)$, when the dimension of $\Gamma$ is zero or one respectively in two
and three dimensions.

In this tutorial program we're not interested in further details on $\gamma$: we
take the extension $\gamma$ for granted, assuming that the dimension of the
embedded domain (`dim`) is always smaller by one or equal with respect to the
dimension of the embedding domain $\Omega$ (`spacedim`).

We are going to solve the following differential problem: given a sufficiently
regular function $g$ on $\Gamma$, find the solution $u$ to

@f{eqnarray*}{
- \Delta u + \gamma^T \lambda &=& 0  \text{ in } \Omega\\
\gamma u &=& g  \text{ in } \Gamma \\
u & = & 0 \text{ on } \partial\Omega.
@f}

This is a constrained problem, where we are looking for a harmonic function $u$
that satisfies homogeneous boundary conditions on $\partial\Omega$, subject to
the constraint $\gamma u = g$ using a Lagrange multiplier.

This problem has a physical interpretation: harmonic functions, i.e., functions
that satisfy the Laplace equation, can be thought of as the displacements of a
membrane whose boundary values are prescribed. The current situation then
corresponds to finding the shape of a membrane for which not only the
displacement at the boundary, but also on $\Gamma$ is prescribed. For example,
if $\Gamma$ is a closed curve in 2d space, then that would model a soap film
that is held in place by a wire loop along $\partial \Omega$ as well as a second
loop along $\Gamma$. In cases where $\Gamma$ is a whole area, you can think of
this as a membrane that is stretched over an obstacle where $\Gamma$ is the
contact area. (If the contact area is not known we have a different problem --
called the "obstacle problem" -- which is modeled in step-41.)

As a first example we study the zero Dirichlet boundary condition on
$\partial\Omega$. The same equations apply if we apply zero Neumann boundary
conditions on $\partial\Omega$ or a mix of the two.

The variational formulation can be derived by introducing two infinite
dimensional spaces $V(\Omega)$ and $Q^*(\Gamma)$, respectively for the solution
$u$ and for the Lagrange multiplier $\lambda$.

Multiplying the first equation by $v \in V(\Omega)$ and the second by $q \in
Q(\Gamma)$, integrating by parts when possible, and exploiting the boundary
conditions on $\partial\Omega$, we obtain the following variational problem:

Given a sufficiently regular function $g$ on $\Gamma$, find the solution $u$ to
@f{eqnarray*}{
(\nabla u, \nabla v)_{\Omega} + (\lambda, \gamma v)_{\Gamma} &=& 0 \qquad \forall v \in V(\Omega) \\
(\gamma u, q)_{\Gamma} &=& (g,q)_{\Gamma} \qquad \forall q \in Q(\Gamma),
@f}

where $(\cdot, \cdot)_{\Omega}$ and $(\cdot, \cdot)_{\Gamma}$ represent,
respectively, $L^2$ scalar products in $\Omega$ and in $\Gamma$.

Inspection of the variational formulation tells us that the space $V(\Omega)$
can be taken to be $H^1_0(\Omega)$. The space $Q(\Gamma)$, in the co-dimension
zero case, should be taken as $H^1(\Gamma)$, while in the co-dimension one case
should be taken as $H^{1/2}(\Gamma)$.

The function $g$ should therefore be either in $H^1(\Gamma)$ (for the
co-dimension zero case) or $H^{1/2}(\Gamma)$ (for the co-dimension one case).
This leaves us with a Lagrange multiplier $\lambda$ in $Q^*(\Gamma)$, which is
either $H^{-1}(\Gamma)$ or $H^{-1/2}(\Gamma)$.

There are two options for the discretization of the problem above. One could choose
matching discretizations, where the Triangulation for $\Gamma$ is aligned with the
Triangulation for $\Omega$, or one could choose to discretize the two domains in
a completely independent way.

The first option is clearly more indicated for the simple problem we
proposed above: it is sufficient to use a single Triangulation for $\Omega$ and
then impose certain constraints depending $\Gamma$. An example of this approach
is studied in step-40, where the solution has to stay above an obstacle and this
is achieved imposing constraints on $\Omega$.

To solve more complex problems, for example one where the domain $\Gamma$ is time
dependent, the second option could be a more viable solution. Handling
non aligned meshes is complex by itself: to illustrate how is done we study a
simple problem.

The technique we describe here is presented in the literature using one of many names:
the <b>immersed finite element method</b>, the <b>fictitious boundary method</b>, the
<b>distributed Lagrange multiplier method</b>, and others. The main principle is
that the discretization of the two grids and of the two finite element spaces
are kept completely independent. This technique is particularly efficient for
the simulation of fluid-structure interaction problems, where the configuration
of the embedded structure is part of the problem itself, and one solves a
(possibly non-linear) elastic problem to determine the (time dependent)
configuration of $\Gamma$, and a (possibly non-linear) flow problem in $\Omega
\setminus \Gamma$, plus coupling conditions on the interface between the fluid
and the solid.

In this tutorial program we keep things a little simpler, and we assume that the
configuration of the embedded domain is given in one of two possible ways:

- as a deformation mapping $\psi: \Gamma_0 \mapsto \Gamma \subseteq \Omega$,
defined on a continuous finite dimensional space on $\Gamma_0$ and representing,
for any point $x \in \Gamma_0$, its coordinate $\psi(x)$ in $\Omega$;

- as a displacement mapping $\delta \psi(x) = \psi(x)-x$ for $x\in \Gamma_0$,
representing for any point $x$ the displacement vector applied in order to
deform $x$ to its actual configuration $\psi(x) = x +\delta\psi(x)$.

We define the embedded reference domain $\Gamma_0$ `embedded_grid`: on
this triangulation we construct a finite dimensional space (`embedded_configuration_dh`)
to describe either the deformation or the displacement through a FiniteElement
system of FE_Q objects (`embedded_configuration_fe`). This finite dimensional
space is used only to interpolate a user supplied function
(`embedded_configuration_function`) representing either $\psi$ (if the
parameter `use_displacement` is set to @p false) or $\delta\psi$ (if the
parameter `use_displacement` is set to @p true).

The Lagrange multiplier $\lambda$ and the user supplied function $g$ are
defined through another finite dimensional space `embedded_dh`, and through
another FiniteElement `embedded_fe`, using the same reference domain. In
order to take into account the deformation of the domain, either a MappingFEField
or a MappingQEulerian object are initialized with the `embedded_configuration`
vector.

In the embedding space, a standard finite dimensional space `space_dh` is
constructed on the embedding grid `space_grid`, using the
FiniteElement `space_fe`, following almost verbatim the approach taken in step-6.

We represent the discretizations of the spaces $V$ and $Q$ with
\f[
V_h(\Omega) = \text{span} \{v_i\}_{i=1}^n
\f]
and
\f[
Q_h(\Gamma) = \text{span} \{q_i\}_{i=1}^m
\f]
respectively, where $n$ is the dimension of `space_dh`, and $m$
the dimension of `embedded_dh`.

Once all the finite dimensional spaces are defined, the variational formulation
of the problem above leaves us with the following finite dimensional system
of equations:

\f[
\begin{pmatrix}
K & C^T \\
C & 0
\end{pmatrix}
\begin{pmatrix}
u \\
\lambda
\end{pmatrix}
=
\begin{pmatrix}
0 \\
G
\end{pmatrix}
\f]

where

@f{eqnarray*}{
K_{ij} &\dealcoloneq& (\nabla v_j, \nabla v_i)_\Omega   \qquad i,j=1,\dots,n \\
C_{\alpha j} &\dealcoloneq& (v_j, q_\alpha)_\Gamma  \qquad j=1,\dots,n, \alpha = 1,\dots, m \\\\
G_{\alpha} &\dealcoloneq& (g, q_\alpha)_\Gamma \qquad \alpha = 1,\dots, m.
@f}

While the matrix $K$ is the standard stiffness matrix for the Poisson problem on
$\Omega$, and the vector $G$ is a standard right-hand-side vector for a finite
element problem with forcing term $g$ on $\Gamma$, (see, for example, step-3),
the matrix $C$ or its transpose $C^T$ are non-standard since they couple
information on two non-matching grids.

In particular, the integral that appears in the computation of a single entry of
$C$, is computed on $\Gamma$. As usual in finite elements we split this
integral into contributions from all cells of the triangulation used to
discretize $\Gamma$, we transform the integral on $K$ to an integral on the
reference element $\hat K$, where $F_{K}$ is the mapping from $\hat K$ to $K$,
and compute the integral on $\hat K$ using a quadrature formula:

\f[
C_{\alpha j} \dealcoloneq (v_j, q_\alpha)_\Gamma  = \sum_{K\in \Gamma} \int_{\hat K}
\hat q_\alpha(\hat x) (v_j \circ F_{K}) (\hat x) J_K (\hat x) \mathrm{d} \hat x =
\sum_{K\in \Gamma} \sum_{i=1}^{n_q}  \big(\hat q_\alpha(\hat x_i)  (v_j \circ F_{K}) (\hat x_i) J_K (\hat x_i) w_i \big)
\f]

Computing this sum is non-trivial because we have to evaluate $(v_j \circ F_{K})
(\hat x_i)$. In general, if $\Gamma$ and $\Omega$ are not aligned, the point
$F_{K}(\hat x_i)$ is completely arbitrary with respect to $\Omega$, and unless
we figure out a way to interpolate all basis functions of $V_h(\Omega)$ on an
arbitrary point on $\Omega$, we cannot compute the integral needed for an entry
of the matrix $C$.

To evaluate $(v_j \circ F_{K}) (\hat x_i)$ the following steps needs to be
taken (as shown in the picture below):

- For a given cell $K$ in $\Gamma$ compute the real point $y_i \dealcoloneq F_{K} (\hat
x_i)$, where $x_i$ is one of the quadrature points used for the integral on $K
\subseteq \Gamma$.

- Find the cell of $\Omega$ in which $y_i$ lies. We shall call this element $T$.

- To evaluate the basis function use the inverse of the mapping $G_T$ that
transforms the reference element $\hat T$ into the element $T$: $v_j(y_i) = \hat
v_j \circ G^{-1}_{T} (y_i)$.

<p align="center"> <img
  src="https://www.dealii.org/images/steps/developer/step-60.C_interpolation.png"
  alt=""> </p>

The three steps above can be computed by calling, in turn,

- GridTools::find_active_cell_around_point(), followed by

- Mapping::transform_real_to_unit_cell(). We then

- construct a custom Quadrature formula, containing the point in the reference
 cell and then

- construct an FEValues object, with the given quadrature formula, and
 initialized with the cell obtained in the first step.

This is what the deal.II function VectorTools::point_value() does when
evaluating a finite element field (not just a single shape function) at an
arbitrary point; but this would be inefficient in this case.

A better solution is to use a convenient wrapper to perform the first three
steps on a collection of points: GridTools::compute_point_locations(). If one is
actually interested in computing the full coupling matrix, then it is possible
to call the method NonMatching::create_coupling_mass_matrix(), that performs the
above steps in an efficient way, reusing all possible data structures, and
gathering expensive steps together. This is the function we'll be using later in
this tutorial.

We solve the final saddle point problem by an iterative solver, applied to the
Schur complement $S$ (whose construction is described, for example, in step-20),
and we construct $S$ using LinearOperator classes.


<a name="Thetestcase"></a><h3>The testcase</h3>


The problem we solve here is identical to step-4, with the difference that we
impose some constraints on an embedded domain $\Gamma$. The tutorial is written
in a dimension independent way, and in the results section we show how to vary
both `dim` and `spacedim`.

The tutorial is compiled for `dim` equal to one and `spacedim` equal to two. If
you want to run the program in embedding dimension `spacedim` equal to three,
you will most likely want to change the reference domain for $\Gamma$ to be, for
example, something you read from file, or a closed sphere that you later deform
to something more interesting.

In the default scenario, $\Gamma$ has co-dimension one, and this tutorial
program implements the Fictitious Boundary Method. As it turns out, the same
techniques are used in the Variational Immersed Finite Element Method, and
the coupling operator $C$ defined above is the same in almost all of these
non-matching methods.

The embedded domain is assumed to be included in $\Omega$, which we take as the
unit square $[0,1]^2$. The definition of the fictitious domain $\Gamma$ can be
modified through the parameter file, and can be given as a mapping from the
reference interval $[0,1]$ to a curve in $\Omega$.

If the curve is closed, then the results will be similar to running the same
problem on a grid whose boundary is $\Gamma$. The program will happily run also
with a non-closed $\Gamma$, although in those cases the mathematical
formulation of the problem is more difficult, since $\Gamma$ will have a
boundary by itself that has co-dimension two with respect to the domain
$\Omega$.


<a name="References"></a><h3>References</h3>


<ul>
<li> Glowinski, R., T.-W. Pan, T.I. Hesla, and D.D. Joseph. 1999. “A Distributed
  Lagrange Multiplier/fictitious Domain Method for Particulate Flows.”
  International Journal of Multiphase Flow 25 (5). Pergamon: 755–94.

<li> Boffi, D., L. Gastaldi, L. Heltai, and C.S. Peskin. 2008. “On the
  Hyper-Elastic Formulation of the Immersed Boundary Method.” Computer Methods
  in Applied Mechanics and Engineering 197 (25–28).

<li> Heltai, L., and F. Costanzo. 2012. “Variational Implementation of Immersed
  Finite Element Methods.” Computer Methods in Applied Mechanics and Engineering
  229–232.
</ul>
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * Most of these have been introduced elsewhere, we'll comment only on the new
 * ones.
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/timer.h>
 * 
 * @endcode
 * 
 * The parameter acceptor class is the first novelty of this tutorial program:
 * in general parameter files are used to steer the execution of a program at
 * run time. While even a simple approach saves compile time, as the same
 * executable can be run with different parameter settings, it can become
 * difficult to handle hundreds of parameters simultaneously while maintaining
 * compatibility between different programs. This is where the class
 * ParameterAcceptor proves useful.
 * 

 * 
 * This class is used to define a public interface for classes that want to use
 * a single global ParameterHandler to handle parameters. The class provides a
 * static ParameterHandler member, namely ParameterAcceptor::prm, and
 * implements the "Command design pattern" (see, for example, E. Gamma, R. Helm,
 * R. Johnson, J. Vlissides, Design Patterns: Elements of Reusable
 * Object-Oriented Software, Addison-Wesley Professional, 1994.
 * https://goo.gl/FNYByc).
 * 

 * 
 * ParameterAcceptor provides a global subscription mechanism. Whenever an
 * object of a class derived from ParameterAcceptor is constructed, a pointer
 * to that object-of-derived-type is registered, together with a section entry
 * in the parameter file. Such registry is traversed upon invocation of the
 * single function ParameterAcceptor::initialize("file.prm") which in turn makes
 * sure that all classes stored in the global registry declare the parameters
 * they will be using, and after having declared them, it reads the content of
 * `file.prm` to parse the actual parameters.
 * 

 * 
 * If you call the method ParameterHandler::add_parameter for each of the
 * parameters you want to use in your code, there is nothing else you need to
 * do. If you are using an already existing class that provides the two
 * functions `declare_parameters` and `parse_parameters`, you can still use
 * ParameterAcceptor, by encapsulating the existing class into a
 * ParameterAcceptorProxy class.
 * 

 * 
 * In this example, we'll use both strategies, using ParameterAcceptorProxy for
 * deal.II classes, and deriving our own parameter classes directly from
 * ParameterAcceptor.
 * 
 * @code
 * #include <deal.II/base/parameter_acceptor.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * @endcode
 * 
 * The other new include file is the one that contains the GridTools::Cache
 * class. The structure of deal.II, as many modern numerical libraries, is
 * organized following a Directed Acyclic Graph (DAG). A DAG is a directed graph
 * with topological ordering: each node structurally represents an object, and
 * is connected to non-root nodes by one (or more) oriented edges, from the
 * parent to the child. The most significant example of this structure is the
 * Triangulation and its Triangulation::cell_iterator structure. From a
 * Triangulation (the main node), we can access each cell (children nodes of the
 * triangulation). From the cells themselves we can access over all vertices of
 * the cell. In this simple example, the DAG structure can be represented as
 * three node types (the triangulation, the cell iterator, and the vertex)
 * connected by oriented edges from the triangulation to the cell iterators, and
 * from the cell iterator to the vertices. This has several advantages, but it
 * intrinsically creates “asymmetries”, making certain operations fast and their
 * inverse very slow: finding the vertices of a cell has low computational cost,
 * and can be done by simply traversing the DAG, while finding all the cells
 * that share a vertex requires a non-trivial computation unless a new DAG data
 * structure is added that represents the inverse search.
 * 

 * 
 * Since inverse operations are usually not needed in a finite element code,
 * these are implemented in GridTools without the use of extra data structures
 * related to the Triangulation which would make them much faster. One such data
 * structure, for example, is a map from the vertices of a Triangulation to all
 * cells that share those vertices, which would reduce the computations needed
 * to answer to the previous question.
 * 

 * 
 * Some methods, for example GridTools::find_active_cell_around_point, make
 * heavy usage of these non-standard operations. If you need to call these
 * methods more than once, it becomes convenient to store those data structures
 * somewhere. GridTools::Cache does exactly this, giving you access to
 * previously computed objects, or computing them on the fly (and then storing
 * them inside the class for later use), and making sure that whenever the
 * Triangulation is updated, also the relevant data structures are recomputed.
 * 
 * @code
 * #include <deal.II/grid/grid_tools_cache.h>
 * 
 * #include <deal.II/fe/fe.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * 
 * @endcode
 * 
 * In this example, we will be using a reference domain to describe an embedded
 * Triangulation, deformed through a finite element vector field.
 * 

 * 
 * The next two include files contain the definition of two classes that can be
 * used in these cases. MappingQEulerian allows one to describe a domain through
 * a *displacement* field, based on a FESystem[FE_Q(p)^spacedim] finite element
 * space. The second is a little more generic, and allows you to use arbitrary
 * vector FiniteElement spaces, as long as they provide a *continuous*
 * description of your domain. In this case, the description is done through the
 * actual *deformation* field, rather than a *displacement* field.
 * 

 * 
 * Which one is used depends on how the user wants to specify the reference
 * domain, and/or the actual configuration. We'll provide both options, and
 * experiment a little in the results section of this tutorial program.
 * 
 * @code
 * #include <deal.II/fe/mapping_q_eulerian.h>
 * #include <deal.II/fe/mapping_fe_field.h>
 * 
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * The parsed function class is another new entry. It allows one to create a
 * Function object, starting from a string in a parameter file which is parsed
 * into an object that you can use anywhere deal.II accepts a Function (for
 * example, for interpolation, boundary conditions, etc.).
 * 
 * @code
 * #include <deal.II/base/parsed_function.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * @endcode
 * 
 * This is the last new entry for this tutorial program. The namespace
 * NonMatching contains a few methods that are useful when performing
 * computations on non-matching grids, or on curves that are not aligned with
 * the underlying mesh.
 * 

 * 
 * We'll discuss its use in detail later on in the `setup_coupling` method.
 * 
 * @code
 * #include <deal.II/non_matching/coupling.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_direct.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/linear_operator.h>
 * #include <deal.II/lac/linear_operator_tools.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * namespace Step60
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="DistributedLagrangeProblem"></a> 
 * <h3>DistributedLagrangeProblem</h3>
 *   

 * 
 * In the DistributedLagrangeProblem, we need two parameters describing the
 * dimensions of the domain $\Gamma$ (`dim`) and of the domain $\Omega$
 * (`spacedim`).
 *   

 * 
 * These will be used to initialize a Triangulation<dim,spacedim> (for
 * $\Gamma$) and a Triangulation<spacedim,spacedim> (for $\Omega$).
 *   

 * 
 * A novelty with respect to other tutorial programs is the heavy use of
 * std::unique_ptr. These behave like classical pointers, with the advantage
 * of doing automatic house-keeping: the contained object is automatically
 * destroyed as soon as the unique_ptr goes out of scope, even if it is inside
 * a container or there's an exception. Moreover it does not allow for
 * duplicate pointers, which prevents ownership problems. We do this, because
 * we want to be able to i) construct the problem, ii) read the parameters,
 * and iii) initialize all objects according to what is specified in a
 * parameter file.
 *   

 * 
 * We construct the parameters of our problem in the internal class
 * `Parameters`, derived from ParameterAcceptor. The
 * `DistributedLagrangeProblem` class takes a const reference to a
 * `Parameters` object, so that it is not possible
 * to modify the parameters from within the DistributedLagrangeProblem class
 * itself.
 *   

 * 
 * We could have initialized the parameters first, and then pass the
 * parameters to the DistributedLagrangeProblem assuming all entries are set
 * to the desired values, but this has two disadvantages:
 *   

 * 
 * - We should not make assumptions on how the user initializes a class that
 * is not under our direct control. If the user fails to initialize the
 * class, we should notice and throw an exception;
 *   

 * 
 * - Not all objects that need to read parameters from a parameter file may
 * be available when we construct the Parameters;
 * this is often the case for complex programs, with multiple physics, or
 * where we reuse existing code in some external classes. We simulate this by
 * keeping some "complex" objects, like ParsedFunction objects, inside the
 * `DistributedLagrangeProblem` instead of inside the
 * `Parameters`.
 *   

 * 
 * Here we assume that upon construction, the classes that build up our
 * problem are not usable yet. Parsing the parameter file is what ensures we
 * have all ingredients to build up our classes, and we design them so that if
 * parsing fails, or is not executed, the run is aborted.
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim = dim>
 *   class DistributedLagrangeProblem
 *   {
 *   public:
 * @endcode
 * 
 * The `Parameters` class is derived from ParameterAcceptor. This allows us
 * to use the ParameterAcceptor::add_parameter() method in its constructor.
 *     

 * 
 * The members of this function are all non-const, but the
 * `DistributedLagrangeProblem` class takes a const reference to a
 * `Parameters` object: this ensures that
 * parameters are not modified from within the `DistributedLagrangeProblem`
 * class.
 * 
 * @code
 *     class Parameters : public ParameterAcceptor
 *     {
 *     public:
 *       Parameters();
 * 
 * @endcode
 * 
 * The parameters now described can all be set externally using a
 * parameter file: if no parameter file is present when running the
 * executable, the program will create a "parameters.prm" file with the
 * default values defined here, and then abort to give the user a chance
 * to modify the parameters.prm file.
 * 

 * 
 * Initial refinement for the embedding grid, corresponding to the domain
 * $\Omega$.
 * 
 * @code
 *       unsigned int initial_refinement = 4;
 * 
 * @endcode
 * 
 * The interaction between the embedded grid $\Omega$ and the embedding
 * grid $\Gamma$ is handled through the computation of $C$, which
 * involves all cells of $\Omega$ overlapping with parts of $\Gamma$:
 * a higher refinement of such cells might improve quality of our
 * computations.
 * For this reason we define `delta_refinement`: if it is greater
 * than zero, then we mark each cell of the space grid that contains
 * a vertex of the embedded grid and its neighbors, execute the
 * refinement, and repeat this process `delta_refinement` times.
 * 
 * @code
 *       unsigned int delta_refinement = 3;
 * 
 * @endcode
 * 
 * Starting refinement of the embedded grid, corresponding to the domain
 * $\Gamma$.
 * 
 * @code
 *       unsigned int initial_embedded_refinement = 8;
 * 
 * @endcode
 * 
 * The list of boundary ids where we impose homogeneous Dirichlet boundary
 * conditions. On the remaining boundary ids (if any), we impose
 * homogeneous Neumann boundary conditions.
 * As a default problem we have zero Dirichlet boundary conditions on
 * $\partial \Omega$
 * 
 * @code
 *       std::list<types::boundary_id> homogeneous_dirichlet_ids{0, 1, 2, 3};
 * 
 * @endcode
 * 
 * FiniteElement degree of the embedding space: $V_h(\Omega)$
 * 
 * @code
 *       unsigned int embedding_space_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * FiniteElement degree of the embedded space: $Q_h(\Gamma)$
 * 
 * @code
 *       unsigned int embedded_space_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * FiniteElement degree of the space used to describe the deformation
 * of the embedded domain
 * 
 * @code
 *       unsigned int embedded_configuration_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * Order of the quadrature formula used to integrate the coupling
 * 
 * @code
 *       unsigned int coupling_quadrature_order = 3;
 * 
 * @endcode
 * 
 * If set to true, then the embedded configuration function is
 * interpreted as a displacement function
 * 
 * @code
 *       bool use_displacement = false;
 * 
 * @endcode
 * 
 * Level of verbosity to use in the output
 * 
 * @code
 *       unsigned int verbosity_level = 10;
 * 
 * @endcode
 * 
 * A flag to keep track if we were initialized or not
 * 
 * @code
 *       bool initialized = false;
 *     };
 * 
 *     DistributedLagrangeProblem(const Parameters &parameters);
 * 
 * @endcode
 * 
 * Entry point for the DistributedLagrangeProblem
 * 
 * @code
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * Object containing the actual parameters
 * 
 * @code
 *     const Parameters &parameters;
 * 
 * @endcode
 * 
 * The following functions are similar to all other tutorial programs, with
 * the exception that we now need to set up things for two different
 * families of objects, namely the ones related to the *embedding* grids,
 * and the ones related to the *embedded* one.
 * 

 * 
 * 
 * @code
 *     void setup_grids_and_dofs();
 * 
 *     void setup_embedding_dofs();
 * 
 *     void setup_embedded_dofs();
 * 
 * @endcode
 * 
 * The only unconventional function we have here is the `setup_coupling()`
 * method, used to generate the sparsity patter for the coupling matrix $C$.
 * 

 * 
 * 
 * @code
 *     void setup_coupling();
 * 
 *     void assemble_system();
 * 
 *     void solve();
 * 
 *     void output_results();
 * 
 * 
 * @endcode
 * 
 * first we gather all the objects related to the embedding space geometry
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<Triangulation<spacedim>> space_grid;
 *     std::unique_ptr<GridTools::Cache<spacedim, spacedim>>
 *                                              space_grid_tools_cache;
 *     std::unique_ptr<FiniteElement<spacedim>> space_fe;
 *     std::unique_ptr<DoFHandler<spacedim>>    space_dh;
 * 
 * @endcode
 * 
 * Then the ones related to the embedded grid, with the DoFHandler
 * associated to the Lagrange multiplier `lambda`
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<Triangulation<dim, spacedim>> embedded_grid;
 *     std::unique_ptr<FiniteElement<dim, spacedim>> embedded_fe;
 *     std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_dh;
 * 
 * @endcode
 * 
 * And finally, everything that is needed to *deform* the embedded
 * triangulation
 * 
 * @code
 *     std::unique_ptr<FiniteElement<dim, spacedim>> embedded_configuration_fe;
 *     std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_configuration_dh;
 *     Vector<double>                                embedded_configuration;
 * 
 * @endcode
 * 
 * The ParameterAcceptorProxy class is a "transparent" wrapper derived
 * from both ParameterAcceptor and the type passed as its template
 * parameter. At construction, the arguments are split into two parts: the
 * first argument is an std::string, forwarded to the ParameterAcceptor
 * class, and containing the name of the section that should be used for
 * this class, while all the remaining arguments are forwarded to the
 * constructor of the templated type, in this case, to the
 * Functions::ParsedFunction constructor.
 *     

 * 
 * This class allows you to use existing classes in conjunction with the
 * ParameterAcceptor registration mechanism, provided that those classes
 * have the members `declare_parameters()` and `parse_parameters()`.
 *     

 * 
 * This is the case here, making it fairly easy to exploit the
 * Functions::ParsedFunction class: instead of requiring users to create new
 * Function objects in their code for the RHS, boundary functions, etc.,
 * (like it is done in most of the other tutorials), here we allow the user
 * to use deal.II interface to muParser (http://muparser.beltoforion.de),
 * where the specification of the function is not done at compile time, but
 * at run time, using a string that is parsed into an actual Function
 * object.
 *     

 * 
 * In this case, the `embedded_configuration_function` is a vector valued
 * Function that can be interpreted as either a *deformation* or a
 * *displacement* according to the boolean value of
 * `parameters.use_displacement`. The number of components is specified
 * later on in the construction.
 * 

 * 
 * 
 * @code
 *     ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
 *       embedded_configuration_function;
 * 
 *     std::unique_ptr<Mapping<dim, spacedim>> embedded_mapping;
 * 
 * @endcode
 * 
 * We do the same thing to specify the value of the function $g$,
 * which is what we want our solution to be in the embedded space.
 * In this case the Function is a scalar one.
 * 
 * @code
 *     ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
 *       embedded_value_function;
 * 
 * @endcode
 * 
 * Similarly to what we have done with the Functions::ParsedFunction class,
 * we repeat the same for the ReductionControl class, allowing us to
 * specify all possible stopping criteria for the Schur complement
 * iterative solver we'll use later on.
 * 
 * @code
 *     ParameterAcceptorProxy<ReductionControl> schur_solver_control;
 * 
 * @endcode
 * 
 * Next we gather all SparsityPattern, SparseMatrix, and Vector objects
 * we'll need
 * 
 * @code
 *     SparsityPattern stiffness_sparsity;
 *     SparsityPattern coupling_sparsity;
 * 
 *     SparseMatrix<double> stiffness_matrix;
 *     SparseMatrix<double> coupling_matrix;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     Vector<double> solution;
 *     Vector<double> rhs;
 * 
 *     Vector<double> lambda;
 *     Vector<double> embedded_rhs;
 *     Vector<double> embedded_value;
 * 
 * @endcode
 * 
 * The TimerOutput class is used to provide some statistics on
 * the performance of our program.
 * 
 * @code
 *     TimerOutput monitor;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="DistributedLagrangeProblemParameters"></a> 
 * <h3>DistributedLagrangeProblem::Parameters</h3>
 *   

 * 
 * At construction time, we initialize also the ParameterAcceptor class, with
 * the section name we want our problem to use when parsing the parameter
 * file.
 *   

 * 
 * Parameter files can be organized into section/subsection/etc.:
 * this has the advantage that defined objects share parameters when
 * sharing the same section/subsection/etc. ParameterAcceptor allows
 * to specify the section name using Unix conventions on paths.
 * If the section name starts with a slash ("/"), then the section is
 * interpreted as an *absolute path*, ParameterAcceptor enters a subsection
 * for each directory in the path, using the last name it encountered as
 * the landing subsection for the current class.
 *   

 * 
 * For example, if you construct your class using
 * `ParameterAcceptor("/first/second/third/My Class")`, the parameters will be
 * organized as follows:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *        ... # all the parameters
 *       end
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * Internally, the *current path* stored in ParameterAcceptor is now
 * considered to be "/first/second/third/", i.e. when you specify an
 * absolute path, ParameterAcceptor *changes* the current section to the
 * current path, i.e. to the path of the section name until the *last* "/".
 *   

 * 
 * You can now construct another class derived from ParameterAcceptor using a
 * relative path (e.g., `ParameterAcceptor("My Other Class")`) instead of the
 * absolute one (e.g. `ParameterAcceptor("/first/second/third/My Other
 * Class")`), obtaining:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *         ... # all the parameters
 *       end
 *       subsection My Other Class
 *         ... # all the parameters of MyOtherClass
 *       end
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * If the section name *ends* with a slash then subsequent classes will
 * interpret this as a full path: for example, similar to the one above, if
 * we have two classes, one initialized with
 * `ParameterAcceptor("/first/second/third/My Class/")`
 * and the other with `ParameterAcceptor("My Other Class")`, then the
 * resulting parameter file will look like:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *         ... # all the parameters of MyClass
 *         ... # notice My Class subsection does not end here
 *         subsection My Other Class
 *           ... # all the parameters of MyOtherClass
 *         end # of subsection My Other Class
 *       end # of subsection My Class
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * We are going to exploit this, by making our
 * `Parameters` the *parent* of all subsequently
 * constructed classes. Since most of the other classes are members of
 * `DistributedLagrangeProblem` this allows, for example, to construct two
 * `DistributedLagrangeProblem` for two different dimensions, without having
 * conflicts in the parameters for the two problems.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   DistributedLagrangeProblem<dim, spacedim>::Parameters::Parameters()
 *     : ParameterAcceptor("/Distributed Lagrange<" +
 *                         Utilities::int_to_string(dim) + "," +
 *                         Utilities::int_to_string(spacedim) + ">/")
 *   {
 * @endcode
 * 
 * The ParameterAcceptor::add_parameter() function does a few things:
 *     

 * 
 * - enters the subsection specified at construction time to
 * ParameterAcceptor
 *     

 * 
 * - calls the ParameterAcceptor::prm.add_parameter() function
 *     

 * 
 * - calls any signal you may have attached to
 * ParameterAcceptor::declare_parameters_call_back
 *     

 * 
 * - leaves the subsection
 *     

 * 
 * In turn, ParameterAcceptor::prm.add_parameter
 *     

 * 
 * - declares an entry in the parameter handler for the given variable;
 *     

 * 
 * - takes the current value of the variable
 *     

 * 
 * - transforms it to a string, used as the default value for the parameter
 * file
 *     

 * 
 * - attaches an *action* to ParameterAcceptor::prm that monitors when a
 * file is parsed, or when an entry is set, and when this happens, it
 * updates the value of the variable passed to `add_parameter()` by setting
 * it to whatever was specified in the input file (of course, after the
 * input file has been parsed and the text representation converted to the
 * type of the variable).
 * 
 * @code
 *     add_parameter("Initial embedding space refinement", initial_refinement);
 * 
 *     add_parameter("Initial embedded space refinement",
 *                   initial_embedded_refinement);
 * 
 *     add_parameter("Local refinements steps near embedded domain",
 *                   delta_refinement);
 * 
 *     add_parameter("Homogeneous Dirichlet boundary ids",
 *                   homogeneous_dirichlet_ids);
 * 
 *     add_parameter("Use displacement in embedded interface", use_displacement);
 * 
 *     add_parameter("Embedding space finite element degree",
 *                   embedding_space_finite_element_degree);
 * 
 *     add_parameter("Embedded space finite element degree",
 *                   embedded_space_finite_element_degree);
 * 
 *     add_parameter("Embedded configuration finite element degree",
 *                   embedded_configuration_finite_element_degree);
 * 
 *     add_parameter("Coupling quadrature order", coupling_quadrature_order);
 * 
 *     add_parameter("Verbosity level", verbosity_level);
 * 
 * @endcode
 * 
 * Once the parameter file has been parsed, then the parameters are good to
 * go. Set the internal variable `initialized` to true.
 * 
 * @code
 *     parse_parameters_call_back.connect([&]() -> void { initialized = true; });
 *   }
 * 
 * @endcode
 * 
 * The constructor is pretty standard, with the exception of the
 * `ParameterAcceptorProxy` objects, as explained earlier.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   DistributedLagrangeProblem<dim, spacedim>::DistributedLagrangeProblem(
 *     const Parameters &parameters)
 *     : parameters(parameters)
 *     , embedded_configuration_function("Embedded configuration", spacedim)
 *     , embedded_value_function("Embedded value")
 *     , schur_solver_control("Schur solver control")
 *     , monitor(std::cout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
 *   {
 * @endcode
 * 
 * Here is a way to set default values for a ParameterAcceptor class
 * that was constructed using ParameterAcceptorProxy.
 *     

 * 
 * In this case, we set the default deformation of the embedded grid to be a
 * circle with radius $R$ and center $(Cx, Cy)$, we set the default value
 * for the embedded_value_function to be the constant one, and specify some
 * sensible values for the SolverControl object.
 *     

 * 
 * It is fundamental for $\Gamma$ to be embedded: from the definition of
 * $C_{\alpha j}$ is clear that, if $\Gamma \not\subseteq \Omega$, certain
 * rows of the matrix $C$ will be zero. This would be a problem, as the
 * Schur complement method requires $C$ to have full column rank.
 * 
 * @code
 *     embedded_configuration_function.declare_parameters_call_back.connect(
 *       []() -> void {
 *         ParameterAcceptor::prm.set("Function constants", "R=.3, Cx=.4, Cy=.4");
 * 
 * 
 *         ParameterAcceptor::prm.set("Function expression",
 *                                    "R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy");
 *       });
 * 
 *     embedded_value_function.declare_parameters_call_back.connect(
 *       []() -> void { ParameterAcceptor::prm.set("Function expression", "1"); });
 * 
 *     schur_solver_control.declare_parameters_call_back.connect([]() -> void {
 *       ParameterAcceptor::prm.set("Max steps", "1000");
 *       ParameterAcceptor::prm.set("Reduction", "1.e-12");
 *       ParameterAcceptor::prm.set("Tolerance", "1.e-12");
 *     });
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Setup"></a> 
 * <h3>Set up</h3>
 *   

 * 
 * The function `DistributedLagrangeProblem::setup_grids_and_dofs()` is used
 * to set up the finite element spaces. Notice how `std::make_unique` is
 * used to create objects wrapped inside `std::unique_ptr` objects.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_grids_and_dofs()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Setup grids and dofs");
 * 
 * @endcode
 * 
 * Initializing $\Omega$: constructing the Triangulation and wrapping it
 * into a `std::unique_ptr` object
 * 
 * @code
 *     space_grid = std::make_unique<Triangulation<spacedim>>();
 * 
 * @endcode
 * 
 * Next, we actually create the triangulation using
 * GridGenerator::hyper_cube(). The last argument is set to true: this
 * activates colorization (i.e., assigning different boundary indicators to
 * different parts of the boundary), which we use to assign the Dirichlet
 * and Neumann conditions.
 * 
 * @code
 *     GridGenerator::hyper_cube(*space_grid, 0, 1, true);
 * 
 * @endcode
 * 
 * Once we constructed a Triangulation, we refine it globally according to
 * the specifications in the parameter file, and construct a
 * GridTools::Cache with it.
 * 
 * @code
 *     space_grid->refine_global(parameters.initial_refinement);
 *     space_grid_tools_cache =
 *       std::make_unique<GridTools::Cache<spacedim, spacedim>>(*space_grid);
 * 
 * @endcode
 * 
 * The same is done with the embedded grid. Since the embedded grid is
 * deformed, we first need to setup the deformation mapping. We do so in the
 * following few lines:
 * 
 * @code
 *     embedded_grid = std::make_unique<Triangulation<dim, spacedim>>();
 *     GridGenerator::hyper_cube(*embedded_grid);
 *     embedded_grid->refine_global(parameters.initial_embedded_refinement);
 * 
 *     embedded_configuration_fe = std::make_unique<FESystem<dim, spacedim>>(
 *       FE_Q<dim, spacedim>(
 *         parameters.embedded_configuration_finite_element_degree),
 *       spacedim);
 * 
 *     embedded_configuration_dh =
 *       std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);
 * 
 *     embedded_configuration_dh->distribute_dofs(*embedded_configuration_fe);
 *     embedded_configuration.reinit(embedded_configuration_dh->n_dofs());
 * 
 * @endcode
 * 
 * Once we have defined a finite dimensional space for the deformation, we
 * interpolate the `embedded_configuration_function` defined in the
 * parameter file:
 * 
 * @code
 *     VectorTools::interpolate(*embedded_configuration_dh,
 *                              embedded_configuration_function,
 *                              embedded_configuration);
 * 
 * @endcode
 * 
 * Now we can interpret it according to what the user has specified in the
 * parameter file: as a displacement, in which case we construct a mapping
 * that *displaces* the position of each support point of our configuration
 * finite element space by the specified amount on the corresponding
 * configuration vector, or as an absolution position.
 *     

 * 
 * In the first case, the class MappingQEulerian offers its services, while
 * in the second one, we'll use the class MappingFEField. They are in fact
 * very similar. MappingQEulerian will only work for systems of FE_Q finite
 * element spaces, where the displacement vector is stored in the first
 * `spacedim` components of the FESystem, and the degree given as a
 * parameter at construction time, must match the degree of the first
 * `spacedim` components.
 *     

 * 
 * The class MappingFEField is slightly more general, in that it allows you
 * to select arbitrary FiniteElement types when constructing your
 * approximation. Naturally some choices may (or may not) make sense,
 * according to the type of FiniteElement you choose. MappingFEField
 * implements the pure iso-parametric concept, and can be used, for example,
 * to implement iso-geometric analysis codes in deal.II, by combining it
 * with the FE_Bernstein finite element class. In this example, we'll use
 * the two interchangeably, by taking into account the fact that one
 * configuration will be a `displacement`, while the other will be an
 * absolute `deformation` field.
 * 

 * 
 * 
 * @code
 *     if (parameters.use_displacement == true)
 *       embedded_mapping =
 *         std::make_unique<MappingQEulerian<dim, Vector<double>, spacedim>>(
 *           parameters.embedded_configuration_finite_element_degree,
 *           *embedded_configuration_dh,
 *           embedded_configuration);
 *     else
 *       embedded_mapping =
 *         std::make_unique<MappingFEField<dim, spacedim, Vector<double>>>(
 *           *embedded_configuration_dh, embedded_configuration);
 * 
 *     setup_embedded_dofs();
 * 
 * @endcode
 * 
 * In this tutorial program we not only refine $\Omega$ globally,
 * but also allow a local refinement depending on the position of $\Gamma$,
 * according to the value of `parameters.delta_refinement`, that we use to
 * decide how many rounds of local refinement we should do on $\Omega$,
 * corresponding to the position of $\Gamma$.
 *     

 * 
 * With the mapping in place, it is now possible to query what is the
 * location of all support points associated with the `embedded_dh`, by
 * calling the method DoFTools::map_dofs_to_support_points.
 *     

 * 
 * This method has two variants. One that does *not* take a Mapping, and
 * one that takes a Mapping. If you use the second type, like we are doing
 * in this case, the support points are computed through the specified
 * mapping, which can manipulate them accordingly.
 *     

 * 
 * This is precisely what the `embedded_mapping` is there for.
 * 
 * @code
 *     std::vector<Point<spacedim>> support_points(embedded_dh->n_dofs());
 *     if (parameters.delta_refinement != 0)
 *       DoFTools::map_dofs_to_support_points(*embedded_mapping,
 *                                            *embedded_dh,
 *                                            support_points);
 * 
 * @endcode
 * 
 * Once we have the support points of the embedded finite element space, we
 * would like to identify what cells of the embedding space contain what
 * support point, to get a chance at refining the embedding grid where it is
 * necessary, i.e., where the embedded grid is. This can be done manually,
 * by looping over each support point, and then calling the method
 * Mapping::transform_real_to_unit_cell for each cell of the embedding
 * space, until we find one that returns points in the unit reference cell,
 * or it can be done in a more intelligent way.
 *     

 * 
 * The GridTools::find_active_cell_around_point is a possible option that
 * performs the above task in a cheaper way, by first identifying the
 * closest vertex of the embedding Triangulation to the target point, and
 * then by calling Mapping::transform_real_to_unit_cell only for those cells
 * that share the found vertex.
 *     

 * 
 * In fact, there are algorithms in the GridTools namespace that exploit a
 * GridTools::Cache object, and possibly a KDTree object to speed up these
 * operations as much as possible.
 *     

 * 
 * The simplest way to exploit the maximum speed is by calling a
 * specialized method, GridTools::compute_point_locations, that will store a
 * lot of useful information and data structures during the first point
 * search, and then reuse all of this for subsequent points.
 *     

 * 
 * GridTools::compute_point_locations returns a tuple where the first
 * element is a vector of cells containing the input points, in this
 * case support_points. For refinement, this is the only information we
 * need, and this is exactly what happens now.
 *     

 * 
 * When we need to assemble a coupling matrix, however, we'll also need the
 * reference location of each point to evaluate the basis functions of the
 * embedding space. The other elements of the tuple returned by
 * GridTools::compute_point_locations allow you to reconstruct, for each
 * point, what cell contains it, and what is the location in the reference
 * cell of the given point. Since this information is better grouped into
 * cells, then this is what the algorithm returns: a tuple, containing a
 * vector of all cells that have at least one point in them, together with a
 * list of all reference points and their corresponding index in the
 * original vector.
 *     

 * 
 * In the following loop, we will be ignoring all returned objects except
 * the first, identifying all cells contain at least one support point of
 * the embedded space. This allows for a simple adaptive refinement
 * strategy: refining these cells and their neighbors.
 *     

 * 
 * Notice that we need to do some sanity checks, in the sense that we want
 * to have an embedding grid which is well refined around the embedded grid,
 * but where two consecutive support points lie either in the same cell, or
 * in neighbor embedding cells.
 *     

 * 
 * This is only possible if we ensure that the smallest cell size of the
 * embedding grid is nonetheless bigger than the largest cell size of the
 * embedded grid. Since users can modify both levels of refinements, as well
 * as the amount of local refinement they want around the embedded grid, we
 * make sure that the resulting meshes satisfy our requirements, and if this
 * is not the case, we bail out with an exception.
 * 
 * @code
 *     for (unsigned int i = 0; i < parameters.delta_refinement; ++i)
 *       {
 *         const auto point_locations =
 *           GridTools::compute_point_locations(*space_grid_tools_cache,
 *                                              support_points);
 *         const auto &cells = std::get<0>(point_locations);
 *         for (auto &cell : cells)
 *           {
 *             cell->set_refine_flag();
 *             for (const auto face_no : cell->face_indices())
 *               if (!cell->at_boundary(face_no))
 *                 cell->neighbor(face_no)->set_refine_flag();
 *           }
 *         space_grid->execute_coarsening_and_refinement();
 *       }
 * 
 * @endcode
 * 
 * In order to construct a well posed coupling interpolation operator $C$,
 * there are some constraints on the relative dimension of the grids between
 * the embedding and the embedded domains. The coupling operator $C$ and the
 * spaces $V$ and $Q$ have to satisfy an inf-sup condition in order for the
 * problem to have a solution. It turns out that the non-matching $L^2$
 * projection satisfies such inf-sup, provided that the spaces $V$ and $Q$
 * are compatible between each other (for example, provided that they are
 * chosen to be the ones described in the introduction).
 *     

 * 
 * However, the *discrete* inf-sup condition must also hold. No
 * complications arise here, but it turns out that the discrete inf-sup
 * constant deteriorates when the non-matching grids have local diameters
 * that are too far away from each other. In particular, it turns out that
 * if you choose an embedding grid which is *finer* with respect to the
 * embedded grid, the inf-sup constant deteriorates much more than if you
 * let the embedded grid be finer.
 *     

 * 
 * In order to avoid issues, in this tutorial we will throw an exception if
 * the parameters chosen by the user are such that the maximal diameter of
 * the embedded grid is greater than the minimal diameter of the embedding
 * grid.
 *     

 * 
 * This choice guarantees that almost every cell of the embedded grid spans
 * no more than two cells of the embedding grid, with some rare exceptions,
 * that are negligible in terms of the resulting inf-sup.
 * 
 * @code
 *     const double embedded_space_maximal_diameter =
 *       GridTools::maximal_cell_diameter(*embedded_grid, *embedded_mapping);
 *     double embedding_space_minimal_diameter =
 *       GridTools::minimal_cell_diameter(*space_grid);
 * 
 *     deallog << "Embedding minimal diameter: "
 *             << embedding_space_minimal_diameter
 *             << ", embedded maximal diameter: "
 *             << embedded_space_maximal_diameter << ", ratio: "
 *             << embedded_space_maximal_diameter /
 *                  embedding_space_minimal_diameter
 *             << std::endl;
 * 
 *     AssertThrow(embedded_space_maximal_diameter <
 *                   embedding_space_minimal_diameter,
 *                 ExcMessage(
 *                   "The embedding grid is too refined (or the embedded grid "
 *                   "is too coarse). Adjust the parameters so that the minimal "
 *                   "grid size of the embedding grid is larger "
 *                   "than the maximal grid size of the embedded grid."));
 * 
 * @endcode
 * 
 * $\Omega$ has been refined and we can now set up its DoFs
 * 
 * @code
 *     setup_embedding_dofs();
 *   }
 * 
 * @endcode
 * 
 * We now set up the DoFs of $\Omega$ and $\Gamma$: since they are
 * fundamentally independent (except for the fact that $\Omega$'s mesh is more
 * refined "around"
 * $\Gamma$) the procedure is standard.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_embedding_dofs()
 *   {
 *     space_dh = std::make_unique<DoFHandler<spacedim>>(*space_grid);
 *     space_fe = std::make_unique<FE_Q<spacedim>>(
 *       parameters.embedding_space_finite_element_degree);
 *     space_dh->distribute_dofs(*space_fe);
 * 
 *     DoFTools::make_hanging_node_constraints(*space_dh, constraints);
 *     for (auto id : parameters.homogeneous_dirichlet_ids)
 *       {
 *         VectorTools::interpolate_boundary_values(
 *           *space_dh, id, Functions::ZeroFunction<spacedim>(), constraints);
 *       }
 *     constraints.close();
 * 
 * @endcode
 * 
 * By definition the stiffness matrix involves only $\Omega$'s DoFs
 * 
 * @code
 *     DynamicSparsityPattern dsp(space_dh->n_dofs(), space_dh->n_dofs());
 *     DoFTools::make_sparsity_pattern(*space_dh, dsp, constraints);
 *     stiffness_sparsity.copy_from(dsp);
 *     stiffness_matrix.reinit(stiffness_sparsity);
 *     solution.reinit(space_dh->n_dofs());
 *     rhs.reinit(space_dh->n_dofs());
 * 
 *     deallog << "Embedding dofs: " << space_dh->n_dofs() << std::endl;
 *   }
 * 
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_embedded_dofs()
 *   {
 *     embedded_dh = std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);
 *     embedded_fe = std::make_unique<FE_Q<dim, spacedim>>(
 *       parameters.embedded_space_finite_element_degree);
 *     embedded_dh->distribute_dofs(*embedded_fe);
 * 
 * @endcode
 * 
 * By definition the rhs of the system we're solving involves only a zero
 * vector and $G$, which is computed using only $\Gamma$'s DoFs
 * 
 * @code
 *     lambda.reinit(embedded_dh->n_dofs());
 *     embedded_rhs.reinit(embedded_dh->n_dofs());
 *     embedded_value.reinit(embedded_dh->n_dofs());
 * 
 *     deallog << "Embedded dofs: " << embedded_dh->n_dofs() << std::endl;
 *   }
 * 
 * @endcode
 * 
 * Creating the coupling sparsity pattern is a complex operation,
 * but it can be easily done using the
 * NonMatching::create_coupling_sparsity_pattern, which requires the
 * two DoFHandler objects, the quadrature points for the coupling,
 * a DynamicSparsityPattern (which then needs to be copied into the
 * sparsity one, as usual), the component mask for the embedding and
 * embedded Triangulation (which we leave empty) and the mappings
 * for both the embedding and the embedded Triangulation.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_coupling()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Setup coupling");
 * 
 *     QGauss<dim> quad(parameters.coupling_quadrature_order);
 * 
 *     DynamicSparsityPattern dsp(space_dh->n_dofs(), embedded_dh->n_dofs());
 * 
 *     NonMatching::create_coupling_sparsity_pattern(*space_grid_tools_cache,
 *                                                   *space_dh,
 *                                                   *embedded_dh,
 *                                                   quad,
 *                                                   dsp,
 *                                                   AffineConstraints<double>(),
 *                                                   ComponentMask(),
 *                                                   ComponentMask(),
 *                                                   *embedded_mapping);
 *     coupling_sparsity.copy_from(dsp);
 *     coupling_matrix.reinit(coupling_sparsity);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Assembly"></a> 
 * <h3>Assembly</h3>
 *   

 * 
 * The following function creates the matrices: as noted before computing the
 * stiffness matrix and the rhs is a standard procedure.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::assemble_system()
 *   {
 *     {
 *       TimerOutput::Scope timer_section(monitor, "Assemble system");
 * 
 * @endcode
 * 
 * Embedding stiffness matrix $K$, and the right hand side $G$.
 * 
 * @code
 *       MatrixTools::create_laplace_matrix(
 *         *space_dh,
 *         QGauss<spacedim>(2 * space_fe->degree + 1),
 *         stiffness_matrix,
 *         static_cast<const Function<spacedim> *>(nullptr),
 *         constraints);
 * 
 *       VectorTools::create_right_hand_side(*embedded_mapping,
 *                                           *embedded_dh,
 *                                           QGauss<dim>(2 * embedded_fe->degree +
 *                                                       1),
 *                                           embedded_value_function,
 *                                           embedded_rhs);
 *     }
 *     {
 *       TimerOutput::Scope timer_section(monitor, "Assemble coupling system");
 * 
 * @endcode
 * 
 * To compute the coupling matrix we use the
 * NonMatching::create_coupling_mass_matrix tool, which works similarly to
 * NonMatching::create_coupling_sparsity_pattern.
 * 
 * @code
 *       QGauss<dim> quad(parameters.coupling_quadrature_order);
 *       NonMatching::create_coupling_mass_matrix(*space_grid_tools_cache,
 *                                                *space_dh,
 *                                                *embedded_dh,
 *                                                quad,
 *                                                coupling_matrix,
 *                                                AffineConstraints<double>(),
 *                                                ComponentMask(),
 *                                                ComponentMask(),
 *                                                *embedded_mapping);
 * 
 *       VectorTools::interpolate(*embedded_mapping,
 *                                *embedded_dh,
 *                                embedded_value_function,
 *                                embedded_value);
 *     }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solve"></a> 
 * <h3>Solve</h3>
 *   

 * 
 * All parts have been assembled: we solve the system
 * using the Schur complement method
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::solve()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Solve system");
 * 
 * @endcode
 * 
 * Start by creating the inverse stiffness matrix
 * 
 * @code
 *     SparseDirectUMFPACK K_inv_umfpack;
 *     K_inv_umfpack.initialize(stiffness_matrix);
 * 
 * @endcode
 * 
 * Initializing the operators, as described in the introduction
 * 
 * @code
 *     auto K  = linear_operator(stiffness_matrix);
 *     auto Ct = linear_operator(coupling_matrix);
 *     auto C  = transpose_operator(Ct);
 * 
 *     auto K_inv = linear_operator(K, K_inv_umfpack);
 * 
 * @endcode
 * 
 * Using the Schur complement method
 * 
 * @code
 *     auto                     S = C * K_inv * Ct;
 *     SolverCG<Vector<double>> solver_cg(schur_solver_control);
 *     auto S_inv = inverse_operator(S, solver_cg, PreconditionIdentity());
 * 
 *     lambda = S_inv * embedded_rhs;
 * 
 *     solution = K_inv * Ct * lambda;
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * @endcode
 * 
 * The following function simply generates standard result output on two
 * separate files, one for each mesh.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::output_results()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Output results");
 * 
 *     DataOut<spacedim> embedding_out;
 * 
 *     std::ofstream embedding_out_file("embedding.vtu");
 * 
 *     embedding_out.attach_dof_handler(*space_dh);
 *     embedding_out.add_data_vector(solution, "solution");
 *     embedding_out.build_patches(
 *       parameters.embedding_space_finite_element_degree);
 *     embedding_out.write_vtu(embedding_out_file);
 * 
 * @endcode
 * 
 * The only difference between the two output routines is that in the
 * second case, we want to output the data on the current configuration, and
 * not on the reference one. This is possible by passing the actual
 * embedded_mapping to the DataOut::build_patches function. The mapping will
 * take care of outputting the result on the actual deformed configuration.
 * 

 * 
 * 
 * @code
 *     DataOut<dim, DoFHandler<dim, spacedim>> embedded_out;
 * 
 *     std::ofstream embedded_out_file("embedded.vtu");
 * 
 *     embedded_out.attach_dof_handler(*embedded_dh);
 *     embedded_out.add_data_vector(lambda, "lambda");
 *     embedded_out.add_data_vector(embedded_value, "g");
 *     embedded_out.build_patches(*embedded_mapping,
 *                                parameters.embedded_space_finite_element_degree);
 *     embedded_out.write_vtu(embedded_out_file);
 *   }
 * 
 * @endcode
 * 
 * Similar to all other tutorial programs, the `run()` function simply calls
 * all other methods in the correct order. Nothing special to note, except
 * that we check if parsing was done before we actually attempt to run our
 * program.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::run()
 *   {
 *     AssertThrow(parameters.initialized, ExcNotInitialized());
 *     deallog.depth_console(parameters.verbosity_level);
 * 
 *     setup_grids_and_dofs();
 *     setup_coupling();
 *     assemble_system();
 *     solve();
 *     output_results();
 *   }
 * } // namespace Step60
 * 
 * 
 * 
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step60;
 * 
 *       const unsigned int dim = 1, spacedim = 2;
 * 
 * @endcode
 * 
 * Differently to what happens in other tutorial programs, here we use
 * ParameterAcceptor style of initialization, i.e., all objects are first
 * constructed, and then a single call to the static method
 * ParameterAcceptor::initialize is issued to fill all parameters of the
 * classes that are derived from ParameterAcceptor.
 *       

 * 
 * We check if the user has specified a parameter file name to use when
 * the program was launched. If so, try to read that parameter file,
 * otherwise, try to read the file "parameters.prm".
 *       

 * 
 * If the parameter file that was specified (implicitly or explicitly)
 * does not exist, ParameterAcceptor::initialize will create one for you,
 * and exit the program.
 * 

 * 
 * 
 * @code
 *       DistributedLagrangeProblem<dim, spacedim>::Parameters parameters;
 *       DistributedLagrangeProblem<dim, spacedim>             problem(parameters);
 * 
 *       std::string parameter_file;
 *       if (argc > 1)
 *         parameter_file = argv[1];
 *       else
 *         parameter_file = "parameters.prm";
 * 
 *       ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");
 *       problem.run();
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
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


The directory in which this program is run does not contain a parameter file by
default. On the other hand, this program wants to read its parameters from a
file called parameters.prm -- and so, when you execute it the first time, you
will get an exception that no such file can be found:

@code
----------------------------------------------------
Exception on processing:

--------------------------------------------------------
An error occurred in line <74> of file <../source/base/parameter_acceptor.cc> in function
    static void dealii::ParameterAcceptor::initialize(const std::string &, const std::string &, const ParameterHandler::OutputStyle, dealii::ParameterHandler &)
The violated condition was:
    false
Additional information:
    You specified <parameters.prm> as input parameter file, but it does not exist. We created it for you.
--------------------------------------------------------

Aborting!
----------------------------------------------------
@endcode

However, as the error message already states, the code that triggers the
exception will also generate a parameters.prm file that simply contains the
default values for all parameters this program cares about. By inspection of the
parameter file, we see the following:

@code
# Listing of Parameters
# ---------------------
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0, 1, 2, 3
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10


  subsection Embedded configuration
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  = R=.3, Cx=.4, Cy=.4                 # default:

    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy # default: 0

    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end

  subsection Embedded value
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  =

    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = 1     # default: 0

    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end

  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 1000   # default: 100
    set Reduction     = 1.e-12 # default: 1.e-2
    set Tolerance     = 1.e-12 # default: 1.e-10
  end

end
@endcode

If you now run the program, you will get a file called `used_parameters.prm`,
containing a shorter version of the above parameters (without comments and
documentation), documenting all parameters that were used to run your program:
@code
# Parameter file generated with
# DEAL_II_PACKAGE_VERSION = 9.0.0
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0, 1, 2, 3
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.4, Cy=.4
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = 1
    set Variable names      = x,y,t
  end
  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 1000
    set Reduction     = 1.e-12
    set Tolerance     = 1.e-12
  end
end
@endcode

The rationale behind creating first `parameters.prm` file (the first time the
program is run) and then a `used_parameters.prm` (every other times you run the
program), is because you may want to leave most parameters to their default
values, and only modify a handful of them.

For example, you could use the following (perfectly valid) parameter file with
this tutorial program:
@code
subsection Distributed Lagrange<1,2>
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.4, Cy=.4
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = 1
    set Variable names      = x,y,t
  end
end
@endcode

and you would obtain exactly the same results as in test case 1 below.

<a name="Testcase1"></a><h3> Test case 1: </h3>


For the default problem the value of $u$ on $\Gamma$ is set to the constant $1$:
this is like imposing a constant Dirichlet boundary condition on $\Gamma$, seen
as boundary of the portion of $\Omega$ inside $\Gamma$. Similarly on $\partial
\Omega$ we have zero Dirichlet boundary conditions.


<div class="twocolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-60.1_no_grid.png"
           alt = ""
           width="500">
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-60.1_grid.png"
           alt = ""
           width="500">
    </div>
  </div>
</div>

The output of the program will look like the following:

@code
DEAL::Embedded dofs: 129
DEAL::Embedding minimal diameter: 0.0110485, embedded maximal diameter: 0.00781250, ratio: 0.707107
DEAL::Embedding dofs: 2429
DEAL:cg::Starting value 0.166266
DEAL:cg::Convergence step 108 value 7.65958e-13


+---------------------------------------------+------------+------------+
| Total CPU time elapsed since start          |     0.586s |            |
|                                             |            |            |
| Section                         | no. calls |  CPU time  | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble coupling system        |         1 |     0.132s |        23% |
| Assemble system                 |         1 |    0.0733s |        12% |
| Output results                  |         1 |     0.087s |        15% |
| Setup coupling                  |         1 |    0.0244s |       4.2% |
| Setup grids and dofs            |         1 |    0.0907s |        15% |
| Solve system                    |         1 |     0.178s |        30% |
+---------------------------------+-----------+------------+------------+



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.301s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble coupling system        |         1 |    0.0385s |        13% |
| Assemble system                 |         1 |    0.0131s |       4.3% |
| Output results                  |         1 |    0.0736s |        24% |
| Setup coupling                  |         1 |    0.0234s |       7.7% |
| Setup grids and dofs            |         1 |    0.0679s |        23% |
| Solve system                    |         1 |    0.0832s |        28% |
+---------------------------------+-----------+------------+------------+

@endcode

You may notice that, in terms of CPU time, assembling the coupling system is
twice as expensive as assembling the standard Poisson system, even though the
matrix is smaller. This is due to the non-matching nature of the discretization.
Whether this is acceptable or not, depends on the applications.

If the problem was set in a three-dimensional setting, and the immersed mesh was
time dependent, it would be much more expensive to recreate the mesh at each
step rather than use the technique we present here. Moreover, you may be able to
create a very fast and optimized solver on a uniformly refined square or cubic
grid, and embed the domain where you want to perform your computation using the
technique presented here. This would require you to only have a surface
representatio of your domain (a much cheaper and easier mesh to produce).

To play around a little bit, we are going to complicate a little the fictitious
domain as well as the boundary conditions we impose on it.

<a name="Testcase2and3"></a><h3> Test case 2 and 3: </h3>


If we use the following parameter file:
@code
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0,1,2,3
  set Initial embedded space refinement            = 8
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 4
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.5, Cy=.5, r=.1, w=12
    set Function expression = (R+r*cos(w*pi*x))*cos(2*pi*x)+Cx; (R+r*cos(w*pi*x))*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = x-.5
    set Variable names      = x,y,t
  end
  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 100000
    set Reduction     = 1.e-12
    set Tolerance     = 1.e-12
  end
end
@endcode

We get a "flowery" looking domain, where we impose a linear boundary condition
$g=x-.5$. This test shows that the method is actually quite accurate in
recovering an exactly linear function from its boundary conditions, and even
though the meshes are not aligned, we obtain a pretty good result.

Replacing $x-.5$ with $2(x-.5)^2-2(y-.5)^2$, i.e., modifying the parameter file
such that we have
@code
  ...
  subsection Embedded value
    set Function constants  =
    set Function expression = 2*(x-.5)^2-2*(y-.5)^2
    set Variable names      = x,y,t
  end
@endcode
produces the saddle on the right.

<div class="twocolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-60.3_grid.png"
           alt = ""
           width="500">
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-60.4_grid.png"
           alt = ""
           width="500">
    </div>
  </div>
</div>

<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Runningwithspacedimequaltothree"></a><h4> Running with `spacedim` equal to three</h4>


While the current tutorial program is written for `spacedim` equal to two, there
are only minor changes you have to do in order for the program to run in
different combinations of dimensions.

If you want to run with `spacedim` equal to three and `dim` equal to two, then
you will almost certainly want to perform the following changes:

- use a different reference domain for the embedded grid, maybe reading it from
  a file. It is not possible to construct a smooth closed surface with one
  single parametrization of a square domain, therefore you'll most likely want
  to use a reference domain that is topologically equivalent to a the boundary
  of a sphere.

- use a displacement instead of the deformation to map $\Gamma_0$ into $\Gamma$

<a name="Moregeneraldomains"></a><h4> More general domains </h4>


We have seen in other tutorials (for example in step-5 and step-54) how to read
grids from input files. A nice generalization for this tutorial program would be
to allow the user to select a grid to read from the parameter file itself,
instead of hardcoding the mesh type in the tutorial program itself.

<a name="Preconditioner"></a><h4> Preconditioner</h4>


At the moment, we have no preconditioner on the Schur complement. This is ok for
two dimensional problems, where a few hundred iterations bring the residual down
to the machine precision, but it's not going to work in three dimensions.

It is not obvious what a good preconditioner would be here. The physical problem
we are solving with the Schur complement, is to associate to the Dirichlet data
$g$, the value of the Lagrange multiplier $\lambda$. $\lambda$ can be
interpreted as the *jump* in the normal gradient that needs to be imposed on $u$
across $\Gamma$, in order to obtain the Dirichlet data $g$.

So $S$ is some sort of Neumann to Dirichlet map, and we would like to have a
good approximation for the Dirichlet to Neumann map. A possibility would be to
use a Boundary Element approximation of the problem on $\Gamma$, and construct a
rough approximation of the hyper-singular operator for the Poisson problem
associated to $\Gamma$, which is precisely a Dirichlet to Neumann map.

<a name="ParallelCode"></a><h4> Parallel Code </h4>


The simple code proposed here can serve as a starting point for more
complex problems which, to be solved, need to be run on parallel
code, possibly using distributed meshes (see step-17, step-40, and the
documentation for parallel::shared::Triangulation and
parallel::distributed::Triangulation).

When using non-matching grids in parallel a problem arises: to compute the
matrix $C$ a process needs information about both meshes on the same portion of
real space but, when working with distributed meshes, this information may not
be available, because the locally owned part of the $\Omega$ triangulation
stored on a given processor may not be physically co-located with the locally
owned part of the $\Gamma$ triangulation stored on the same processor.

Various strategies can be implemented to tackle this problem:

- distribute the two meshes so that this constraint is satisfied;

- use communication for the parts of real space where the constraint is not
  satisfied;

- use a distributed triangulation for the embedding space, and a shared
  triangulation for the emdedded configuration.

The latter strategy is clearly the easiest to implement, as most of the
functions used in this tutorial program will work unchanged also in the parallel
case. Of course one could use the reversal strategy (that is, have a distributed
embedded Triangulation and a shared embedding Triangulation).

However, this strategy is most likely going to be more expensive, since by
definition the embedding grid is larger than the embedded grid, and it makes
more sense to distribute the largest of the two grids, maintaining the smallest
one shared among all processors.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-60.cc"
*/
