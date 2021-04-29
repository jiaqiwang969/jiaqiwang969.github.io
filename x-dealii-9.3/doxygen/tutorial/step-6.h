/**
@page step_6 The step-6 tutorial program
This tutorial depends on step-5.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Whatadaptivelyrefinedmesheslooklike"> What adaptively refined meshes look like </a>
        <li><a href="#Whyadapativelyrefinedmeshes"> Why adapatively refined meshes? </a>
        <li><a href="#Howtodealwithhangingnodesintheory"> How to deal with hanging nodes in theory </a>
        <li><a href="#Howtodealwithhangingnodesinpractice"> How to deal with hanging nodes in practice </a>
        <li><a href="#Howweobtainlocallyrefinedmeshes"> How we obtain locally refined meshes </a>
        <li><a href="#Boundaryconditions"> Boundary conditions </a>
        <li><a href="#Otherthingsthisprogramshows"> Other things this program shows </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeStep6codeclasstemplate">The <code>Step6</code> class template</a>
        <li><a href="#Nonconstantcoefficients">Nonconstant coefficients</a>
        <li><a href="#ThecodeStep6codeclassimplementation">The <code>Step6</code> class implementation</a>
      <ul>
        <li><a href="#Step6Step6">Step6::Step6</a>
        <li><a href="#Step6setup_system">Step6::setup_system</a>
        <li><a href="#Step6assemble_system">Step6::assemble_system</a>
        <li><a href="#Step6solve">Step6::solve</a>
        <li><a href="#Step6refine_grid">Step6::refine_grid</a>
        <li><a href="#Step6output_results">Step6::output_results</a>
        <li><a href="#Step6run">Step6::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Solversandpreconditioners">Solvers and preconditioners</a>
        <li><a href="#Abettermesh">A better mesh</a>
        <li><a href="#Playingwiththeregularityofthesolution">Playing with the regularity of the solution</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{15,16,17,17.25,17.5,17.75}

This program is finally about one of the main features of deal.II:
the use of adaptively (locally) refined meshes. The program is still
based on step-4 and step-5, and, as you will see, it does not actually
take very much code to enable adaptivity. Indeed, while we do a great
deal of explaining, adaptive meshes can be added to an existing program
with barely a dozen lines of additional code. The program shows what
these lines are, as well as another important ingredient of adaptive
mesh refinement (AMR): a criterion that can be used to determine whether
it is necessary to refine a cell because the error is large on it,
whether the cell can be coarsened because the error is particularly
small on it, or whether we should just leave the cell as it is. We
will discuss all of these issues in the following.


<a name="Whatadaptivelyrefinedmesheslooklike"></a><h3> What adaptively refined meshes look like </h3>


There are a number of ways how one can adaptively refine meshes. The
basic structure of the overall algorithm is always the same and consists
of a loop over the following steps:
- Solve the PDE on the current mesh;
- Estimate the error on each cell using some criterion that is indicative
  of the error;
- Mark those cells that have large errors for refinement, mark those that have
  particularly small errors for coarsening, and leave the rest alone;
- Refine and coarsen the cells so marked to obtain a new mesh;
- Repeat the steps above on the new mesh until the overall error is
  sufficiently small.

For reasons that are probably lost to history (maybe that these functions
used to be implemented in FORTRAN, a language that does not care about
whether something is spelled in lower or UPPER case letters, with programmers
often choosing upper case letters habitually), the loop above is often
referenced in publications about mesh adaptivity as the
SOLVE-ESTIMATE-MARK-REFINE loop (with this spelling).

Beyond this structure, however, there are a variety of ways to achieve
this. Fundamentally, they differ in how exactly one generates one mesh
from the previous one.

If one were to use triangles (which deal.II does not do), then there are
two essential possibilities:
- Longest-edge refinement: In this strategy, a triangle marked for refinement
  is cut into two by introducing one new edge from the midpoint of the longest
  edge to the opposite vertex. Of course, the midpoint from the longest edge
  has to somehow be balanced by *also* refining the cell on the other side of
  that edge (if there is one). If the edge in question is also the longest
  edge of the neighboring cell, then we can just run a new edge through the
  neighbor to the opposite vertex; otherwise a slightly more involved
  construction is necessary that adds more new vertices on at least one
  other edge of the neighboring cell, and then may propagate to the neighbors
  of the neighbor until the algorithm terminates. This is hard to describe
  in words, and because deal.II does not use triangles not worth the time here.
  But if you're curious, you can always watch video lecture 15 at the link
  shown at the top of this introduction.
- Red-green refinement: An alternative is what is called "red-green refinement".
  This strategy is even more difficult to describe (but also discussed in the
  video lecture) and has the advantage that the refinement does not propagate
  beyond the immediate neighbors of the cell that we want to refine. It is,
  however, substantially more difficult to implement.

There are other variations of these approaches, but the important point is
that they always generate a mesh where the lines where two cells touch
are entire edges of both adjacent cells. With a bit of work, this strategy
is readily adapted to three-dimensional meshes made from tetrahedra.

Neither of these methods works for quadrilaterals in 2d and hexahedra in 3d,
or at least not easily. The reason is that the transition elements created
out of the quadrilateral neighbors of a quadrilateral cell that is to be refined
would be triangles, and we don't want this. Consequently,
the approach to adaptivity chosen in deal.II is to use grids in which
neighboring cells may differ in refinement level by one. This then
results in nodes on the interfaces of cells which belong to one
side, but are unbalanced on the other. The common term for these is
&ldquo;hanging nodes&rdquo;, and these meshes then look like this in a very
simple situation:

@image html hanging_nodes.png "A simple mesh with hanging nodes"

A more complicated two-dimensional mesh would look like this (and is
discussed in the "Results" section below):

<img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg"
     alt="Fifth adaptively refined Ladutenko grid: the cells are clustered
          along the inner circle."
     width="300" height="300">

Finally, a three-dimensional mesh (from step-43) with such hanging nodes is shown here:

<img src="https://www.dealii.org/images/steps/developer/step-43.3d.mesh.png" alt=""
     width="300" height="300">

The first and third mesh are of course based on a square and a cube, but as the
second mesh shows, this is not necessary. The important point is simply that we
can refine a mesh independently of its neighbors (subject to the constraint
that a cell can be only refined once more than its neighbors), but that we end
up with these &ldquo;hanging nodes&rdquo; if we do this.


<a name="Whyadapativelyrefinedmeshes"></a><h3> Why adapatively refined meshes? </h3>


Now that you have seen what these adaptively refined meshes look like,
you should ask <i>why</i> we would want to do this. After all, we know from
theory that if we refine the mesh globally, the error will go down to zero
as
@f{align*}{
  \|\nabla(u-u_h)\|_{\Omega} \le C h_\text{max}^p \| \nabla^{p+1} u \|_{\Omega},
@f}
where $C$ is some constant independent of $h$ and $u$,
$p$ is the polynomial degree of the finite element in use, and
$h_\text{max}$ is the diameter of the largest cell. So if the
<i>largest</i> cell is important, then why would we want to make
the mesh fine in some parts of the domain but not all?

The answer lies in the observation that the formula above is not
optimal. In fact, some more work shows that the following
is a better estimate (which you should compare to the square of
the estimate above):
@f{align*}{
  \|\nabla(u-u_h)\|_{\Omega}^2 \le C \sum_K h_K^{2p} \| \nabla^{p+1} u \|^2_K.
@f}
(Because $h_K\le h_\text{max}$, this formula immediately implies the
previous one if you just pull the mesh size out of the sum.)
What this formula suggests is that it is not necessary to make
the <i>largest</i> cell small, but that the cells really only
need to be small <i>where $\| \nabla^{p+1} u \|_K$ is large</i>!
In other words: The mesh really only has to be fine where the
solution has large variations, as indicated by the $p+1$st derivative.
This makes intuitive sense: if, for example, we use a linear element
$p=1$, then places where the solution is nearly linear (as indicated
by $\nabla^2 u$ being small) will be well resolved even if the mesh
is coarse. Only those places where the second derivative is large
will be poorly resolved by large elements, and consequently
that's where we should make the mesh small.

Of course, this <i>a priori estimate</i> is not very useful
in practice since we don't know the exact solution $u$ of the
problem, and consequently, we cannot compute $\nabla^{p+1}u$.
But, and that is the approach commonly taken, we can compute
numerical approximations of $\nabla^{p+1}u$ based only on
the discrete solution $u_h$ that we have computed before. We
will discuss this in slightly more detail below. This will then
help us determine which cells have a large $p+1$st derivative,
and these are then candidates for refining the mesh.


<a name="Howtodealwithhangingnodesintheory"></a><h3> How to deal with hanging nodes in theory </h3>


The methods using triangular meshes mentioned above go to great
lengths to make sure that each vertex is a vertex of all adjacent
cells -- i.e., that there are no hanging nodes. This then
automatically makes sure that we can define shape functions in such a
way that they are globally continuous (if we use the common $Q_p$
Lagrange finite element methods we have been using so far in the
tutorial programs, as represented by the FE_Q class).

On the other hand, if we define shape functions on meshes with hanging
nodes, we may end up with shape functions that are not continuous. To
see this, think about the situation above where the top right cell is
not refined, and consider for a moment the use of a bilinear finite
element. In that case, the shape functions associated with the hanging
nodes are defined in the obvious way on the two small cells adjacent
to each of the hanging nodes. But how do we extend them to the big
adjacent cells? Clearly, the function's extension to the big cell
cannot be bilinear because then it needs to be linear along each edge
of the large cell, and that means that it needs to be zero on the
entire edge because it needs to be zero on the two vertices of the
large cell on that edge. But it is not zero at the hanging node itself
when seen from the small cells' side -- so it is not continuous. The
following three figures show three of the shape functions along the
edges in question that turn out to not be continuous when defined in
the usual way simply based on the cells they are adjacent to:

<div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_1.png "A discontinuous shape function adjacent to a hanging node"
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_2.png "A discontinuous shape function at a hanging node"
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_3.png "A discontinuous shape function adjacent to a hanging node"
    </div>
  </div>
</div>


But we do want the finite element solution to be continuous so that we
have a &ldquo;conforming finite element method&rdquo; where the
discrete finite element space is a proper subset of the $H^1$ function
space in which we seek the solution of the Laplace equation.
To guarantee that the global solution is continuous at these nodes as well, we
have to state some additional constraints on the values of the solution at
these nodes. The trick is to realize that while the shape functions shown
above are discontinuous (and consequently an <i>arbitrary</i> linear combination
of them is also discontinuous), that linear combinations in which the shape
functions are added up as $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$
can be continuous <i>if the coefficients $U_j$ satisfy certain relationships</i>.
In other words, the coefficients $U_j$ can not be chosen arbitrarily
but have to satisfy certain constraints so that the function $u_h$ is in fact
continuous.
What these constraints have to look is relatively easy to
understand conceptually, but the implementation in software is
complicated and takes several thousand lines of code. On the other
hand, in user code, it is only about half a dozen lines you have to
add when dealing with hanging nodes.

In the program below, we will show how we can get these
constraints from deal.II, and how to use them in the solution of the
linear system of equations. Before going over the details of the program
below, you may want to take a look at the @ref constraints documentation
module that explains how these constraints can be computed and what classes in
deal.II work on them.


<a name="Howtodealwithhangingnodesinpractice"></a><h3> How to deal with hanging nodes in practice </h3>


The practice of hanging node constraints is rather simpler than the
theory we have outlined above. In reality, you will really only have to
add about half a dozen lines of additional code to a program like step-4
to make it work with adaptive meshes that have hanging nodes. The
interesting part about this is that it is entirely independent of the
equation you are solving: The algebraic nature of these constraints has nothing
to do with the equation and only depends on the choice of finite element.
As a consequence, the code to deal with these constraints is entirely
contained in the deal.II library itself, and you do not need to worry
about the details.

The steps you need to make this work are essentially like this:
- You have to create an AffineConstraints object, which (as the name
  suggests) will store all constraints on the finite element space. In
  the current context, these are the constraints due to our desire to
  keep the solution space continuous even in the presence of hanging
  nodes. (Below we will also briefly mention that we will also put
  boundary values into this same object, but that is a separate matter.)
- You have to fill this object using the function
  DoFTools::make_hanging_node_constraints() to ensure continuity of
  the elements of the finite element space.
- You have to use this object when you copy the local contributions to
  the matrix and right hand side into the global objects, by using
  AffineConstraints::distribute_local_to_global(). Up until
  now, we have done this ourselves, but now with constraints, this
  is where the magic happens and we apply the constraints to the
  linear system. What this function does is make sure that the
  degrees of freedom located at hanging nodes are not, in fact,
  really free. Rather, they are factually eliminated from the
  linear system by setting their rows and columns to zero and putting
  something on the diagonal to ensure the matrix remains invertible.
  The matrix resulting from this process remains symmetric and
  positive definite for the Laplace equation we solve here, so we can
  continue to use the Conjugate Gradient method for it.
- You then solve the linear system as usual, but at the end of this
  step, you need to make sure that the degrees of "freedom" located
  on hanging nodes get their correct (constrained) value so that the
  solution you then visualize or evaluate in other ways is in
  fact continuous. This is done by calling
  AffineConstraints::distribute() immediately after solving.

These four steps are really all that is necessary -- it's that simple
from a user perspective. The fact that, in the function calls mentioned
above, you will run through several thousand lines of not-so-trivial
code is entirely immaterial to this: In user code, there are really
only four additional steps.


<a name="Howweobtainlocallyrefinedmeshes"></a><h3> How we obtain locally refined meshes </h3>


The next question, now that we know how to <i>deal</i> with meshes that
have these hanging nodes is how we <i>obtain</i> them.

A simple way has already been shown in step-1: If you <i>know</i> where
it is necessary to refine the mesh, then you can create one by hand. But
in reality, we don't know this: We don't know the solution of the PDE
up front (because, if we did, we wouldn't have to use the finite element
method), and consequently we do not know where it is necessary to
add local mesh refinement to better resolve areas where the solution
has strong variations. But the discussion above shows that maybe we
can get away with using the discrete solution $u_h$ on one mesh to
estimate the derivatives $\nabla^{p+1} u$, and then use this to determine
which cells are too large and which already small enough. We can then
generate a new mesh from the current one using local mesh refinement.
If necessary, this step is then repeated until we are happy with our
numerical solution -- or, more commonly, until we run out of computational
resources or patience.

So that's exactly what we will do.
The locally refined grids are produced using an <i>error estimator</i>
which estimates the energy error for numerical solutions of the Laplace
operator. Since it was developed by Kelly and
co-workers, we often refer to it as the &ldquo;Kelly refinement
indicator&rdquo; in the library, documentation, and mailing list. The
class that implements it is called
KellyErrorEstimator, and there is a great deal of information to
be found in the documentation of that class that need not be repeated
here. The summary, however, is that the class computes a vector with
as many entries as there are @ref GlossActive "active cells", and
where each entry contains an estimate of the error on that cell.
This estimate is then used to refine the cells of the mesh: those
cells that have a large error will be marked for refinement, those
that have a particularly small estimate will be marked for
coarsening. We don't have to do this by hand: The functions in
namespace GridRefinement will do all of this for us once we have
obtained the vector of error estimates.

It is worth noting that while the Kelly error estimator was developed
for Laplace's equation, it has proven to be a suitable tool to generate
locally refined meshes for a wide range of equations, not even restricted
to elliptic only problems. Although it will create non-optimal meshes for other
equations, it is often a good way to quickly produce meshes that are
well adapted to the features of solutions, such as regions of great
variation or discontinuities.



<a name="Boundaryconditions"></a><h3> Boundary conditions </h3>


It turns out that one can see Dirichlet boundary conditions as just another
constraint on the degrees of freedom. It's a particularly simple one,
indeed: If $j$ is a degree of freedom on the boundary, with position
$\mathbf x_j$, then imposing the boundary condition $u=g$ on $\partial\Omega$
simply yields the constraint $U_j=g({\mathbf x}_j)$.

The AffineConstraints class can handle such constraints as well, which makes it
convenient to let the same object we use for hanging node constraints
also deal with these Dirichlet boundary conditions.
This way, we don't need to apply the boundary conditions after assembly
(like we did in the earlier steps).
All that is necessary is that we call the variant of
VectorTools::interpolate_boundary_values() that returns its information
in an AffineConstraints object, rather than the `std::map` we have used
in previous tutorial programs.


<a name="Otherthingsthisprogramshows"></a> <h3> Other things this program shows </h3>



Since the concepts used for locally refined grids are so important,
we do not show much other material in this example. The most
important exception is that we show how to use biquadratic elements
instead of the bilinear ones which we have used in all previous
examples. In fact, the use of higher order elements is accomplished by
only replacing three lines of the program, namely the initialization of
the <code>fe</code> member variable in the constructor of the main
class of this program, and the use of an appropriate quadrature formula
in two places. The rest of the program is unchanged.

The only other new thing is a method to catch exceptions in the
<code>main</code> function in order to output some information in case the
program crashes for some reason. This is discussed below in more detail.
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
 * The first few files have already been covered in previous examples and will
 * thus not be further commented on.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * #include <fstream>
 * 
 * @endcode
 * 
 * From the following include file we will import the declaration of
 * H1-conforming finite element shape functions. This family of finite
 * elements is called <code>FE_Q</code>, and was used in all examples before
 * already to define the usual bi- or tri-linear elements, but we will now use
 * it for bi-quadratic elements:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * @endcode
 * 
 * We will not read the grid from a file as in the previous example, but
 * generate it using a function of the library. However, we will want to write
 * out the locally refined grids (just the grid, not the solution) in each
 * step, so we need the following include file instead of
 * <code>grid_in.h</code>:
 * 
 * @code
 * #include <deal.II/grid/grid_out.h>
 * 
 * 
 * @endcode
 * 
 * When using locally refined grids, we will get so-called <code>hanging
 * nodes</code>. However, the standard finite element methods assumes that the
 * discrete solution spaces be continuous, so we need to make sure that the
 * degrees of freedom on hanging nodes conform to some constraints such that
 * the global solution is continuous. We are also going to store the boundary
 * conditions in this object. The following file contains a class which is
 * used to handle these constraints:
 * 
 * @code
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * @endcode
 * 
 * In order to refine our grids locally, we need a function from the library
 * that decides which cells to flag for refinement or coarsening based on the
 * error indicators we have computed. This function is defined here:
 * 
 * @code
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * @endcode
 * 
 * Finally, we need a simple way to actually compute the refinement indicators
 * based on some error estimate. While in general, adaptivity is very
 * problem-specific, the error indicator in the following file often yields
 * quite nicely adapted grids for a wide class of problems.
 * 
 * @code
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * Finally, this is as in previous programs:
 * 
 * @code
 * using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep6codeclasstemplate"></a> 
 * <h3>The <code>Step6</code> class template</h3>
 * 

 * 
 * The main class is again almost unchanged. Two additions, however, are made:
 * we have added the <code>refine_grid</code> function, which is used to
 * adaptively refine the grid (instead of the global refinement in the
 * previous examples), and a variable which will hold the constraints.
 * 
 * @code
 * template <int dim>
 * class Step6
 * {
 * public:
 *   Step6();
 * 
 *   void run();
 * 
 * private:
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void refine_grid();
 *   void output_results(const unsigned int cycle) const;
 * 
 *   Triangulation<dim> triangulation;
 * 
 *   FE_Q<dim>       fe;
 *   DoFHandler<dim> dof_handler;
 * 
 * 
 * @endcode
 * 
 * This is the new variable in the main class. We need an object which holds
 * a list of constraints to hold the hanging nodes and the boundary
 * conditions.
 * 
 * @code
 *   AffineConstraints<double> constraints;
 * 
 *   SparseMatrix<double> system_matrix;
 *   SparsityPattern      sparsity_pattern;
 * 
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Nonconstantcoefficients"></a> 
 * <h3>Nonconstant coefficients</h3>
 * 

 * 
 * The implementation of nonconstant coefficients is copied verbatim from
 * step-5:
 * 
 * @code
 * template <int dim>
 * double coefficient(const Point<dim> &p)
 * {
 *   if (p.square() < 0.5 * 0.5)
 *     return 20;
 *   else
 *     return 1;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep6codeclassimplementation"></a> 
 * <h3>The <code>Step6</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Step6Step6"></a> 
 * <h4>Step6::Step6</h4>
 * 

 * 
 * The constructor of this class is mostly the same as before, but this time
 * we want to use the quadratic element. To do so, we only have to replace the
 * constructor argument (which was <code>1</code> in all previous examples) by
 * the desired polynomial degree (here <code>2</code>):
 * 
 * @code
 * template <int dim>
 * Step6<dim>::Step6()
 *   : fe(2)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6setup_system"></a> 
 * <h4>Step6::setup_system</h4>
 * 

 * 
 * The next function sets up all the variables that describe the linear
 * finite element problem, such as the DoFHandler, matrices, and
 * vectors. The difference to what we did in step-5 is only that we now also
 * have to take care of hanging node constraints. These constraints are
 * handled almost exclusively by the library, i.e. you only need to know
 * that they exist and how to get them, but you do not have to know how they
 * are formed or what exactly is done with them.
 * 

 * 
 * At the beginning of the function, you find all the things that are the same
 * as in step-5: setting up the degrees of freedom (this time we have
 * quadratic elements, but there is no difference from a user code perspective
 * to the linear -- or any other degree, for that matter -- case), generating
 * the sparsity pattern, and initializing the solution and right hand side
 * vectors. Note that the sparsity pattern will have significantly more
 * entries per row now, since there are now 9 degrees of freedom per cell
 * (rather than only four), that can couple with each other.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 * 
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * We may now populate the AffineConstraints object with the hanging node
 * constraints. Since we will call this function in a loop we first clear
 * the current set of constraints from the last system and then compute new
 * ones:
 * 
 * @code
 *   constraints.clear();
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 * 
 * @endcode
 * 
 * Now we are ready to interpolate the boundary values with indicator 0 (the
 * whole boundary) and store the resulting constraints in our
 * <code>constraints</code> object. Note that we do not to apply the
 * boundary conditions after assembly, like we did in earlier steps: instead
 * we put all constraints on our function space in the AffineConstraints
 * object. We can add constraints to the AffineConstraints object in either
 * order: if two constraints conflict then the constraint matrix either abort
 * or throw an exception via the Assert macro.
 * 
 * @code
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<dim>(),
 *                                            constraints);
 * 
 * @endcode
 * 
 * After all constraints have been added, they need to be sorted and
 * rearranged to perform some actions more efficiently. This postprocessing
 * is done using the <code>close()</code> function, after which no further
 * constraints may be added any more:
 * 
 * @code
 *   constraints.close();
 * 
 * @endcode
 * 
 * Now we first build our compressed sparsity pattern like we did in the
 * previous examples. Nevertheless, we do not copy it to the final sparsity
 * pattern immediately.  Note that we call a variant of
 * make_sparsity_pattern that takes the AffineConstraints object as the third
 * argument. We are letting the routine know that we will never write into
 * the locations given by <code>constraints</code> by setting the argument
 * <code>keep_constrained_dofs</code> to false (in other words, that we will
 * never write into entries of the matrix that correspond to constrained
 * degrees of freedom). If we were to condense the
 * constraints after assembling, we would have to pass <code>true</code>
 * instead because then we would first write into these locations only to
 * later set them to zero again during condensation.
 * 
 * @code
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler,
 *                                   dsp,
 *                                   constraints,
 *                                   /*keep_constrained_dofs = */ false);
 * 
 * @endcode
 * 
 * Now all non-zero entries of the matrix are known (i.e. those from
 * regularly assembling the matrix and those that were introduced by
 * eliminating constraints). We may copy our intermediate object to the
 * sparsity pattern:
 * 
 * @code
 *   sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * We may now, finally, initialize the sparse matrix:
 * 
 * @code
 *   system_matrix.reinit(sparsity_pattern);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6assemble_system"></a> 
 * <h4>Step6::assemble_system</h4>
 * 

 * 
 * Next, we have to assemble the matrix. However, to copy the local matrix and
 * vector on each cell into the global system, we are no longer using a
 * hand-written loop. Instead, we use
 * AffineConstraints::distribute_local_to_global() that internally executes
 * this loop while performing Gaussian elimination on rows and columns
 * corresponding to constrained degrees on freedom.
 * 

 * 
 * The rest of the code that forms the local contributions remains
 * unchanged. It is worth noting, however, that under the hood several things
 * are different than before. First, the variable <code>dofs_per_cell</code>
 * and return value of <code>quadrature_formula.size()</code> now are 9 each,
 * where they were 4 before. Introducing such variables as abbreviations is a
 * good strategy to make code work with different elements without having to
 * change too much code. Secondly, the <code>fe_values</code> object of course
 * needs to do other things as well, since the shape functions are now
 * quadratic, rather than linear, in each coordinate variable. Again, however,
 * this is something that is completely handled by the library.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::assemble_system()
 * {
 *   const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *   FEValues<dim> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients |
 *                             update_quadrature_points | update_JxW_values);
 * 
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 *       fe_values.reinit(cell);
 * 
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         {
 *           const double current_coefficient =
 *             coefficient(fe_values.quadrature_point(q_index));
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               for (const unsigned int j : fe_values.dof_indices())
 *                 cell_matrix(i, j) +=
 *                   (current_coefficient *              // a(x_q)
 *                    fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                    fe_values.JxW(q_index));           // dx
 * 
 *               cell_rhs(i) += (1.0 *                               // f(x)
 *                               fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                               fe_values.JxW(q_index));            // dx
 *             }
 *         }
 * 
 * @endcode
 * 
 * Finally, transfer the contributions from @p cell_matrix and
 * @p cell_rhs into the global objects.
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices);
 *       constraints.distribute_local_to_global(
 *         cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *     }
 * @endcode
 * 
 * Now we are done assembling the linear system. The constraint matrix took
 * care of applying the boundary conditions and also eliminated hanging node
 * constraints. The constrained nodes are still in the linear system (there
 * is a nonzero entry, chosen in a way that the matrix is well conditioned,
 * on the diagonal of the matrix and all other entries for this line are set
 * to zero) but the computed values are invalid (i.e., the corresponding
 * entries in <code>system_rhs</code> are currently meaningless). We compute
 * the correct values for these nodes at the end of the <code>solve</code>
 * function.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6solve"></a> 
 * <h4>Step6::solve</h4>
 * 

 * 
 * We continue with gradual improvements. The function that solves the linear
 * system again uses the SSOR preconditioner, and is again unchanged except
 * that we have to incorporate hanging node constraints. As mentioned above,
 * the degrees of freedom from the AffineConstraints object corresponding to
 * hanging node constraints and boundary values have been removed from the
 * linear system by giving the rows and columns of the matrix a special
 * treatment. This way, the values for these degrees of freedom have wrong,
 * but well-defined values after solving the linear system. What we then have
 * to do is to use the constraints to assign to them the values that they
 * should have. This process, called <code>distributing</code> constraints,
 * computes the values of constrained nodes from the values of the
 * unconstrained ones, and requires only a single additional function call
 * that you find at the end of this function:
 * 

 * 
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::solve()
 * {
 *   SolverControl            solver_control(1000, 1e-12);
 *   SolverCG<Vector<double>> solver(solver_control);
 * 
 *   PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *   preconditioner.initialize(system_matrix, 1.2);
 * 
 *   solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *   constraints.distribute(solution);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6refine_grid"></a> 
 * <h4>Step6::refine_grid</h4>
 * 

 * 
 * We use a sophisticated error estimation scheme to refine the mesh instead
 * of global refinement. We will use the KellyErrorEstimator class which
 * implements an error estimator for the Laplace equation; it can in principle
 * handle variable coefficients, but we will not use these advanced features,
 * but rather use its most simple form since we are not interested in
 * quantitative results but only in a quick way to generate locally refined
 * grids.
 * 

 * 
 * Although the error estimator derived by Kelly et al. was originally
 * developed for the Laplace equation, we have found that it is also well
 * suited to quickly generate locally refined grids for a wide class of
 * problems. This error estimator uses the solution gradient's jump at
 * cell faces (which is a measure for the second derivatives) and
 * scales it by the size of the cell. It is therefore a measure for the local
 * smoothness of the solution at the place of each cell and it is thus
 * understandable that it yields reasonable grids also for hyperbolic
 * transport problems or the wave equation as well, although these grids are
 * certainly suboptimal compared to approaches specially tailored to the
 * problem. This error estimator may therefore be understood as a quick way to
 * test an adaptive program.
 * 

 * 
 * The way the estimator works is to take a <code>DoFHandler</code> object
 * describing the degrees of freedom and a vector of values for each degree of
 * freedom as input and compute a single indicator value for each active cell
 * of the triangulation (i.e. one value for each of the active cells). To do
 * so, it needs two additional pieces of information: a face quadrature formula,
 * i.e., a quadrature formula on <code>dim-1</code> dimensional objects. We use
 * a 3-point Gauss rule again, a choice that is consistent and appropriate with
 * the bi-quadratic finite element shape functions in this program.
 * (What constitutes a suitable quadrature rule here of course depends on
 * knowledge of the way the error estimator evaluates the solution field. As
 * said above, the jump of the gradient is integrated over each face, which
 * would be a quadratic function on each face for the quadratic elements in
 * use in this example. In fact, however, it is the square of the jump of the
 * gradient, as explained in the documentation of that class, and that is a
 * quartic function, for which a 3 point Gauss formula is sufficient since it
 * integrates polynomials up to order 5 exactly.)
 * 

 * 
 * Secondly, the function wants a list of boundary indicators for those
 * boundaries where we have imposed Neumann values of the kind
 * $\partial_n u(\mathbf x) = h(\mathbf x)$, along with a function $h(\mathbf
 * x)$ for each such boundary. This information is represented by a map from
 * boundary indicators to function objects describing the Neumann boundary
 * values. In the present example program, we do not use Neumann boundary
 * values, so this map is empty, and in fact constructed using the default
 * constructor of the map in the place where the function call expects the
 * respective function argument.
 * 

 * 
 * The output is a vector of values for all active cells. While it may
 * make sense to compute the <b>value</b> of a solution degree of freedom
 * very accurately, it is usually not necessary to compute the <b>error
 * indicator</b> corresponding to the solution on a cell particularly
 * accurately. We therefore typically use a vector of floats instead of a vector
 * of doubles to represent error indicators.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::refine_grid()
 * {
 *   Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *   KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                      QGauss<dim - 1>(fe.degree + 1),
 *                                      {},
 *                                      solution,
 *                                      estimated_error_per_cell);
 * 
 * @endcode
 * 
 * The above function returned one error indicator value for each cell in
 * the <code>estimated_error_per_cell</code> array. Refinement is now done
 * as follows: refine those 30 per cent of the cells with the highest error
 * values, and coarsen the 3 per cent of cells with the lowest values.
 *   

 * 
 * One can easily verify that if the second number were zero, this would
 * approximately result in a doubling of cells in each step in two space
 * dimensions, since for each of the 30 per cent of cells, four new would be
 * replaced, while the remaining 70 per cent of cells remain untouched. In
 * practice, some more cells are usually produced since it is disallowed
 * that a cell is refined twice while the neighbor cell is not refined; in
 * that case, the neighbor cell would be refined as well.
 *   

 * 
 * In many applications, the number of cells to be coarsened would be set to
 * something larger than only three per cent. A non-zero value is useful
 * especially if for some reason the initial (coarse) grid is already rather
 * refined. In that case, it might be necessary to refine it in some
 * regions, while coarsening in some other regions is useful. In our case
 * here, the initial grid is very coarse, so coarsening is only necessary in
 * a few regions where over-refinement may have taken place. Thus a small,
 * non-zero value is appropriate here.
 *   

 * 
 * The following function now takes these refinement indicators and flags
 * some cells of the triangulation for refinement or coarsening using the
 * method described above. It is from a class that implements several
 * different algorithms to refine a triangulation based on cell-wise error
 * indicators.
 * 
 * @code
 *   GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                   estimated_error_per_cell,
 *                                                   0.3,
 *                                                   0.03);
 * 
 * @endcode
 * 
 * After the previous function has exited, some cells are flagged for
 * refinement, and some other for coarsening. The refinement or coarsening
 * itself is not performed by now, however, since there are cases where
 * further modifications of these flags is useful. Here, we don't want to do
 * any such thing, so we can tell the triangulation to perform the actions
 * for which the cells are flagged:
 * 
 * @code
 *   triangulation.execute_coarsening_and_refinement();
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6output_results"></a> 
 * <h4>Step6::output_results</h4>
 * 

 * 
 * At the end of computations on each grid, and just before we continue the
 * next cycle with mesh refinement, we want to output the results from this
 * cycle.
 * 

 * 
 * We have already seen in step-1 how this can be achieved for the
 * mesh itself. Here, we change a few things:
 * <ol>
 * <li>We use two different formats: gnuplot and VTU.</li>
 * <li>We embed the cycle number in the output file name.</li>
 * <li>For gnuplot output, we set up a GridOutFlags::Gnuplot object to
 * provide a few extra visualization arguments so that edges appear
 * curved. This is explained in further detail in step-10.</li>
 * </ol>
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::output_results(const unsigned int cycle) const
 * {
 *   {
 *     GridOut               grid_out;
 *     std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot");
 *     GridOutFlags::Gnuplot gnuplot_flags(false, 5);
 *     grid_out.set_flags(gnuplot_flags);
 *     MappingQGeneric<dim> mapping(3);
 *     grid_out.write_gnuplot(triangulation, output, &mapping);
 *   }
 * 
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches();
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
 *     data_out.write_vtu(output);
 *   }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6run"></a> 
 * <h4>Step6::run</h4>
 * 

 * 
 * The final function before <code>main()</code> is again the main driver of
 * the class, <code>run()</code>. It is similar to the one of step-5, except
 * that we generate a file in the program again instead of reading it from
 * disk, in that we adaptively instead of globally refine the mesh, and that
 * we output the solution on the final mesh in the present function.
 * 

 * 
 * The first block in the main loop of the function deals with mesh generation.
 * If this is the first cycle of the program, instead of reading the grid from
 * a file on disk as in the previous example, we now again create it using a
 * library function. The domain is again a circle with center at the origin and
 * a radius of one (these are the two hidden arguments to the function, which
 * have default values).
 * 

 * 
 * You will notice by looking at the coarse grid that it is of inferior
 * quality than the one which we read from the file in the previous example:
 * the cells are less equally formed. However, using the library function this
 * program works in any space dimension, which was not the case before.
 * 

 * 
 * In case we find that this is not the first cycle, we want to refine the
 * grid. Unlike the global refinement employed in the last example program, we
 * now use the adaptive procedure described above.
 * 

 * 
 * The rest of the loop looks as before:
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::run()
 * {
 *   for (unsigned int cycle = 0; cycle < 8; ++cycle)
 *     {
 *       std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *       if (cycle == 0)
 *         {
 *           GridGenerator::hyper_ball(triangulation);
 *           triangulation.refine_global(1);
 *         }
 *       else
 *         refine_grid();
 * 
 * 
 *       std::cout << "   Number of active cells:       "
 *                 << triangulation.n_active_cells() << std::endl;
 * 
 *       setup_system();
 * 
 *       std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *                 << std::endl;
 * 
 *       assemble_system();
 *       solve();
 *       output_results(cycle);
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The main function is unaltered in its functionality from the previous
 * example, but we have taken a step of additional caution. Sometimes,
 * something goes wrong (such as insufficient disk space upon writing an
 * output file, not enough memory when trying to allocate a vector or a
 * matrix, or if we can't read from or write to a file for whatever reason),
 * and in these cases the library will throw exceptions. Since these are
 * run-time problems, not programming errors that can be fixed once and for
 * all, this kind of exceptions is not switched off in optimized mode, in
 * contrast to the <code>Assert</code> macro which we have used to test
 * against programming errors. If uncaught, these exceptions propagate the
 * call tree up to the <code>main</code> function, and if they are not caught
 * there either, the program is aborted. In many cases, like if there is not
 * enough memory or disk space, we can't do anything but we can at least print
 * some text trying to explain the reason why the program failed. A way to do
 * so is shown in the following. It is certainly useful to write any larger
 * program in this way, and you can do so by more or less copying this
 * function except for the <code>try</code> block that actually encodes the
 * functionality particular to the present application.
 * 
 * @code
 * int main()
 * {
 * @endcode
 * 
 * The general idea behind the layout of this function is as follows: let's
 * try to run the program as we did before...
 * 
 * @code
 *   try
 *     {
 *       Step6<2> laplace_problem_2d;
 *       laplace_problem_2d.run();
 *     }
 * @endcode
 * 
 * ...and if this should fail, try to gather as much information as
 * possible. Specifically, if the exception that was thrown is an object of
 * a class that is derived from the C++ standard class
 * <code>exception</code>, then we can use the <code>what</code> member
 * function to get a string which describes the reason why the exception was
 * thrown.
 *   

 * 
 * The deal.II exception classes are all derived from the standard class,
 * and in particular, the <code>exc.what()</code> function will return
 * approximately the same string as would be generated if the exception was
 * thrown using the <code>Assert</code> macro. You have seen the output of
 * such an exception in the previous example, and you then know that it
 * contains the file and line number of where the exception occurred, and
 * some other information. This is also what the following statements would
 * print.
 *   

 * 
 * Apart from this, there isn't much that we can do except exiting the
 * program with an error code (this is what the <code>return 1;</code>
 * does):
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
 * 
 *       return 1;
 *     }
 * @endcode
 * 
 * If the exception that was thrown somewhere was not an object of a class
 * derived from the standard <code>exception</code> class, then we can't do
 * anything at all. We then simply print an error message and exit.
 * 
 * @code
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
 * @endcode
 * 
 * If we got to this point, there was no exception which propagated up to
 * the main function (there may have been exceptions, but they were caught
 * somewhere in the program or the library). Therefore, the program
 * performed as was expected and we can return without error.
 * 
 * @code
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>



The output of the program looks as follows:
@code
Cycle 0:
   Number of active cells:       20
   Number of degrees of freedom: 89
Cycle 1:
   Number of active cells:       44
   Number of degrees of freedom: 209
Cycle 2:
   Number of active cells:       92
   Number of degrees of freedom: 449
Cycle 3:
   Number of active cells:       200
   Number of degrees of freedom: 921
Cycle 4:
   Number of active cells:       440
   Number of degrees of freedom: 2017
Cycle 5:
   Number of active cells:       956
   Number of degrees of freedom: 4425
Cycle 6:
   Number of active cells:       1916
   Number of degrees of freedom: 8993
Cycle 7:
   Number of active cells:       3860
   Number of degrees of freedom: 18353
@endcode



As intended, the number of cells roughly doubles in each cycle. The
number of degrees is slightly more than four times the number of
cells; one would expect a factor of exactly four in two spatial
dimensions on an infinite grid (since the spacing between the degrees
of freedom is half the cell width: one additional degree of freedom on
each edge and one in the middle of each cell), but it is larger than
that factor due to the finite size of the mesh and due to additional
degrees of freedom which are introduced by hanging nodes and local
refinement.



The program outputs the solution and mesh in each cycle of the
refinement loop. The solution looks as follows:

<img src="https://www.dealii.org/images/steps/developer/step-6.solution.9.2.png" alt="">

It is interesting to follow how the program arrives at the final mesh:

<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0.svg"
         alt="Initial grid: the five-cell circle grid with one global refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1.svg"
         alt="First grid: the five-cell circle grid with two global refinements."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2.svg"
         alt="Second grid: the five-cell circle grid with one adaptive refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3.svg"
         alt="Third grid: the five-cell circle grid with two adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4.svg"
         alt="Fourth grid: the five-cell circle grid with three adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5.svg"
         alt="Fifth grid: the five-cell circle grid with four adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_6.svg"
         alt="Sixth grid: the five-cell circle grid with five adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_7.svg"
         alt="Last grid: the five-cell circle grid with six adaptive
         refinements, showing that most cells are clustered around the inner circle."
         width="300" height="300">
  </div>
</div>


It is clearly visible that the region where the solution has a kink,
i.e. the circle at radial distance 0.5 from the center, is
refined most. Furthermore, the central region where the solution is
very smooth and almost flat, is almost not refined at all, but this
results from the fact that we did not take into account that the
coefficient is large there. The region outside is refined rather
arbitrarily, since the second derivative is constant there and refinement
is therefore mostly based on the size of the cells and their deviation
from the optimal square.



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Solversandpreconditioners"></a><h4>Solvers and preconditioners</h4>



One thing that is always worth playing around with if one solves
problems of appreciable size (much bigger than the one we have here)
is to try different solvers or preconditioners. In the current case,
the linear system is symmetric and positive definite, which makes the
CG algorithm pretty much the canonical choice for solving. However,
the SSOR preconditioner we use in the <code>solve()</code> function is
up for grabs.

In deal.II, it is relatively simple to change the preconditioner. For
example, by changing the existing lines of code
@code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
@endcode
into
@code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);
@endcode
we can try out different relaxation parameters for SSOR. By using
@code
  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode
we can use Jacobi as a preconditioner. And by using
@code
  SparseILU<double> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode
we can use a simple incomplete LU decomposition without any thresholding or
strengthening of the diagonal (to use this preconditioner, you have to also
add the header file <code>deal.II/lac/sparse_ilu.h</code> to the include list
at the top of the file).

Using these various different preconditioners, we can compare the
number of CG iterations needed (available through the
<code>solver_control.last_step()</code> call, see
step-4) as well as CPU time needed (using the Timer class,
discussed, for example, in step-28) and get the
following results (left: iterations; right: CPU time):

<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table>

As we can see, all preconditioners behave pretty much the same on this
simple problem, with the number of iterations growing like ${\cal
O}(N^{1/2})$ and because each iteration requires around ${\cal
O}(N)$ operations the total CPU time grows like ${\cal
O}(N^{3/2})$ (for the few smallest meshes, the CPU time is so small
that it doesn't record). Note that even though it is the simplest
method, Jacobi is the fastest for this problem.

The situation changes slightly when the finite element is not a
bi-quadratic one as set in the constructor of this program, but a
bi-linear one. If one makes this change, the results are as follows:

<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table>

In other words, while the increase in iterations and CPU time is as
before, Jacobi is now the method that requires the most iterations; it
is still the fastest one, however, owing to the simplicity of the
operations it has to perform. This is not to say that Jacobi
is actually a good preconditioner -- for problems of appreciable size, it is
definitely not, and other methods will be substantially better -- but really
only that it is fast because its implementation is so simple that it can
compensate for a larger number of iterations.

The message to take away from this is not that simplicity in
preconditioners is always best. While this may be true for the current
problem, it definitely is not once we move to more complicated
problems (elasticity or Stokes, for examples step-8 or
step-22). Secondly, all of these preconditioners still
lead to an increase in the number of iterations as the number $N$ of
degrees of freedom grows, for example ${\cal O}(N^\alpha)$; this, in
turn, leads to a total growth in effort as ${\cal O}(N^{1+\alpha})$
since each iteration takes ${\cal O}(N)$ work. This behavior is
undesirable: we would really like to solve linear systems with $N$
unknowns in a total of ${\cal O}(N)$ work; there is a class
of preconditioners that can achieve this, namely geometric (step-16,
step-37, step-39)
or algebraic multigrid (step-31, step-40, and several others)
preconditioners. They are, however, significantly more complex than
the preconditioners outlined above.

Finally, the last message to take
home is that when the data shown above was generated (in 2018), linear
systems with 100,000 unknowns are
easily solved on a desktop machine in about a second, making
the solution of relatively simple 2d problems even to very high
accuracy not that big a task as it used to be even in the
past. At the time, the situation for 3d problems was entirely different,
but even that has changed substantially in the intervening time -- though
solving problems in 3d to high accuracy remains a challenge.


<a name="Abettermesh"></a><h4>A better mesh</h4>


If you look at the meshes above, you will see even though the domain is the
unit disk, and the jump in the coefficient lies along a circle, the cells
that make up the mesh do not track this geometry well. The reason, already hinted
at in step-1, is that in the absence of other information,
the Triangulation class only sees a bunch of
coarse grid cells but has, of course, no real idea what kind of geometry they
might represent when looked at together. For this reason, we need to tell
the Triangulation what to do when a cell is refined: where should the new
vertices at the edge midpoints and the cell midpoint be located so that the
child cells better represent the desired geometry than the parent cell.

To visualize what the triangulation actually knows about the geometry,
it is not enough to just output the location of vertices and draw a
straight line for each edge; instead, we have to output both interior
and boundary lines as multiple segments so that they look
curved. We can do this by making one change to the gnuplot part of
<code>output_results</code>:
@code
{
  GridOut       grid_out;
  std::ofstream output("grid-" + std::to_string(cycle) + ".gnuplot");
  GridOutFlags::Gnuplot gnuplot_flags(false, 5, /*curved_interior_cells*/true);
  grid_out.set_flags(gnuplot_flags);
  MappingQGeneric<dim> mapping(3);
  grid_out.write_gnuplot(triangulation, output, &mapping);
}
@endcode

In the code above, we already do this for faces that sit at the boundary: this
happens automatically since we use GridGenerator::hyper_ball, which attaches a
SphericalManifold to the boundary of the domain. To make the mesh
<i>interior</i> also track a circular domain, we need to work a bit harder,
though. First, recall that our coarse mesh consists of a central square
cell and four cells around it. Now first consider what would happen if we
also attached the SphericalManifold object not only to the four exterior faces
but also the four cells at the perimeter as well as all of their faces. We can
do this by adding the following snippet (testing that the center of a cell is
larger than a small multiple, say one tenth, of the cell diameter away from
center of the mesh only fails for the central square of the mesh):
@code
GridGenerator::hyper_ball(triangulation);
// after GridGenerator::hyper_ball is called the Triangulation has
// a SphericalManifold with id 0. We can use it again on the interior.
const Point<dim> mesh_center;
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance (cell->center()) > cell->diameter()/10)
    cell->set_all_manifold_ids(0);

triangulation.refine_global(1);
@endcode

After a few global refinement steps, this would lead to a mesh of the following
kind:


  <div class="onecolumn" style="width: 80%">
    <div>
      <img src="https://www.dealii.org/images/steps/developer/step_6_bad_grid_4.svg"
           alt="Grid where some central cells are nearly triangular."
           width="300" height="300">
    </div>
  </div>

This is not a good mesh: the central cell has been refined in such a way that
the children located in the four corners of the original central cell
<i>degenerate</i>: they all tend towards triangles as mesh refinement
continues. This means that the Jacobian matrix of the transformation from
reference cell to actual cell degenerates for these cells, and because
all error estimates for finite element solutions contain the norm of the
inverse of the Jacobian matrix, you will get very large errors on these
cells and, in the limit as mesh refinement, a loss of convergence order because
the cells in these corners become worse and worse under mesh refinement.

So we need something smarter. To this end, consider the following solution
originally developed by Konstantin Ladutenko. We will use the following code:
@code
GridGenerator::hyper_ball(triangulation);

const Point<dim> mesh_center;
const double core_radius  = 1.0/5.0,
             inner_radius = 1.0/3.0;

// Step 1: Shrink the inner cell
//
// We cannot get a circle out of the inner cell because of
// the degeneration problem mentioned above. Rather, shrink
// the inner cell to a core radius of 1/5 that stays
// sufficiently far away from the place where the
// coefficient will have a discontinuity and where we want
// to have cell interfaces that actually lie on a circle.
// We do this shrinking by just scaling the location of each
// of the vertices, given that the center of the circle is
// simply the origin of the coordinate system.
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) < 1e-5)
    {
      for (const auto v : cell->vertex_indices())
        cell->vertex(v) *= core_radius/mesh_center.distance(cell->vertex(v));
    }

// Step 2: Refine all cells except the central one
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) >= 1e-5)
    cell->set_refine_flag();
triangulation.execute_coarsening_and_refinement();

// Step 3: Resize the inner children of the outer cells
//
// The previous step replaced each of the four outer cells
// by its four children, but the radial distance at which we
// have intersected is not what we want to later refinement
// steps. Consequently, move the vertices that were just
// created in radial direction to a place where we need
// them.
for (const auto &cell : triangulation.active_cell_iterators())
  for (const auto v : cell->vertex_indices())
    {
      const double dist = mesh_center.distance(cell->vertex(v));
      if (dist > core_radius*1.0001 && dist < 0.9999)
        cell->vertex(v) *= inner_radius/dist;
    }

// Step 4: Apply curved manifold description
//
// As discussed above, we can not expect to subdivide the
// inner four cells (or their faces) onto concentric rings,
// but we can do so for all other cells that are located
// outside the inner radius. To this end, we loop over all
// cells and determine whether it is in this zone. If it
// isn't, then we set the manifold description of the cell
// and all of its bounding faces to the one that describes
// the spherical manifold already introduced above and that
// will be used for all further mesh refinement.
for (const auto &cell : triangulation.active_cell_iterators())
  {
    bool is_in_inner_circle = false;
    for (const auto v : cell->vertex_indices())
      if (mesh_center.distance(cell->vertex(v)) < inner_radius)
        {
          is_in_inner_circle = true;
          break;
        }

    if (is_in_inner_circle == false)
    // The Triangulation already has a SphericalManifold with
    // manifold id 0 (see the documentation of
    // GridGenerator::hyper_ball) so we just attach it to the outer
    // ring here:
      cell->set_all_manifold_ids(0);
  }
@endcode

This code then generates the following, much better sequence of meshes:

<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0_ladutenko.svg"
         alt="Initial grid: the Ladutenko grid with one global refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1_ladutenko.svg"
         alt="First adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2_ladutenko.svg"
         alt="Second adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3_ladutenko.svg"
         alt="Third adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4_ladutenko.svg"
         alt="Fourth adaptively refined Ladutenko grid. The cells are clustered
         along the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg"
         alt="Fifth adaptively refined Ladutenko grid: the cells are clustered
         along the inner circle."
         width="300" height="300">
  </div>
</div>

Creating good meshes, and in particular making them fit the geometry you
want, is a complex topic in itself. You can find much more on this in
step-49, step-53, and step-54, among other tutorial programs that cover
the issue. step-65 shows another, less manual way to achieve a mesh
well fit to the problem here.
Information on curved domains can also be found in the
documentation module on @ref manifold "Manifold descriptions".

Why does it make sense to choose a mesh that tracks the internal
interface? There are a number of reasons, but the most essential one
comes down to what we actually integrate in our bilinear
form. Conceptually, we want to integrate the term $A_{ij}^K=\int_K
a(\mathbf x) \nabla \varphi_i(\mathbf x) \nabla \varphi_j(\mathbf x) ; dx$ as the
contribution of cell $K$ to the matrix entry $A_{ij}$. We can not
compute it exactly and have to resort to quadrature. We know that
quadrature is accurate if the integrand is smooth. That is because
quadrature in essence computes a polynomial approximation to the
integrand that coincides with the integrand in the quadrature points,
and then computes the volume under this polynomial as an approximation
to the volume under the original integrand. This polynomial
interpolant is accurate if the integrand is smooth on a cell, but it
is usually rather inaccurate if the integrand is discontinuous on a
cell.

Consequently, it is worthwhile to align cells in such a way that the
interfaces across which the coefficient is discontinuous are aligned
with cell interfaces. This way, the coefficient is constant on each
cell, following which the integrand will be smooth, and its polynomial
approximation and the quadrature approximation of the integral will
both be accurate. Note that such an alignment is common in many
practical cases, so deal.II provides a number of functions (such as
@ref GlossMaterialId "material_id") to help manage such a scenario.
Refer to step-28 and step-46 for examples of how material ids can be
applied.

Finally, let us consider the case of a coefficient that has a smooth
and non-uniform distribution in space. We can repeat once again all of
the above discussion on the representation of such a function with the
quadrature. So, to simulate it accurately there are a few readily
available options: you could reduce the cell size, increase the order
of the polynomial used in the quadrature formula, select a more
appropriate quadrature formula, or perform a combination of these
steps. The key is that providing the best fit of the coefficient's
spatial dependence with the quadrature polynomial will lead to a more
accurate finite element solution of the PDE.

As a final note: The discussion in the previous paragraphs shows, we here
have a very concrete way of stating what we think of a good mesh -- it should
be aligned with the jump in the coefficient. But one could also have asked
this kind of question in a more general setting: Given some equation with
a smooth solution and smooth coefficients, can we say what a good mesh
would look like? This is a question for which the answer is easier to state
in intuitive terms than mathematically: A good mesh has cells that all,
by and large, look like squares (or cubes, in 3d). A bad mesh would contain
cells that are very elongated in some directions or, more generally, for which
there are cells that have both short and long edges. There are many ways
in which one could assign a numerical quality index to each cell that measures
whether the cell is "good" or "bad"; some of these are often chosen because
they are cheap and easy to compute, whereas others are based on what enters
into proofs of convergence. An example of the former would be the ratio of
the longest to the shortest edge of a cell: In the ideal case, that ratio
would be one; bad cells have values much larger than one. An example of the
latter kind would consider the gradient (the "Jacobian") of the mapping
from the reference cell $\hat K=[0,1]^d$ to the real cell $K$; this
gradient is a matrix, and a quantity that enters into error estimates
is the maximum over all points on the reference cell of the ratio of the
largest to the smallest eigenvalue of this matrix. It is again not difficult
to see that this ratio is constant if the cell $K$ is an affine image of
$\hat K$, and that it is one for squares and cubes.

In practice, it might be interesting to visualize such quality measures.
The function GridTools::compute_aspect_ratio_of_cells() provides one
way to get this kind of information. Even better, visualization tools
such as VisIt often allow you to visualize this sort of information
for a variety of measures from within the visualization software
itself; in the case of VisIt, just add a "pseudo-color" plot and select
one of the mesh quality measures instead of the solution field.


<a name="Playingwiththeregularityofthesolution"></a><h4>Playing with the regularity of the solution</h4>


From a mathematical perspective, solutions of the Laplace equation
@f[
  -\Delta u = f
@f]
on smoothly bounded, convex domains are known to be smooth themselves. The exact degree
of smoothness, i.e., the function space in which the solution lives, depends
on how smooth exactly the boundary of the domain is, and how smooth the right
hand side is. Some regularity of the solution may be lost at the boundary, but
one generally has that the solution is twice more differentiable in
compact subsets of the domain than the right hand side.
If, in particular, the right hand side satisfies $f\in C^\infty(\Omega)$, then
$u \in C^\infty(\Omega_i)$ where $\Omega_i$ is any compact subset of $\Omega$
($\Omega$ is an open domain, so a compact subset needs to keep a positive distance
from $\partial\Omega$).

The situation we chose for the current example is different, however: we look
at an equation with a non-constant coefficient $a(\mathbf x)$:
@f[
  -\nabla \cdot (a \nabla u) = f.
@f]
Here, if $a$ is not smooth, then the solution will not be smooth either,
regardless of $f$. In particular, we expect that wherever $a$ is discontinuous
along a line (or along a plane in 3d),
the solution will have a kink. This is easy to see: if for example $f$
is continuous, then $f=-\nabla \cdot (a \nabla u)$ needs to be
continuous. This means that $a \nabla u$ must be continuously differentiable
(not have a kink). Consequently, if $a$ has a discontinuity, then $\nabla u$
must have an opposite discontinuity so that the two exactly cancel and their
product yields a function without a discontinuity. But for $\nabla u$ to have
a discontinuity, $u$ must have a kink. This is of course exactly what is
happening in the current example, and easy to observe in the pictures of the
solution.

In general, if the coefficient $a(\mathbf x)$ is discontinuous along a line in 2d,
or a plane in 3d, then the solution may have a kink, but the gradient of the
solution will not go to infinity. That means, that the solution is at least
still in the <a href="https://en.wikipedia.org/wiki/Sobolev_space">Sobolev space</a>
$W^{1,\infty}$ (i.e., roughly speaking, in the
space of functions whose derivatives are bounded). On the other hand,
we know that in the most
extreme cases -- i.e., where the domain has reentrant corners, the
right hand side only satisfies $f\in H^{-1}$, or the coefficient $a$ is only in
$L^\infty$ -- all we can expect is that $u\in H^1$ (i.e., the
<a
href="https://en.wikipedia.org/wiki/Sobolev_space#Sobolev_spaces_with_integer_k">Sobolev
space</a> of functions whose derivative is square integrable), a much larger space than
$W^{1,\infty}$. It is not very difficult to create cases where
the solution is in a space $H^{1+s}$ where we can get $s$ to become as small
as we want. Such cases are often used to test adaptive finite element
methods because the mesh will have to resolve the singularity that causes
the solution to not be in $W^{1,\infty}$ any more.

The typical example one uses for this is called the <i>Kellogg problem</i>
(referring to @cite Kel74), which in the commonly used form has a coefficient
$a(\mathbf x)$ that has different values in the four quadrants of the plane
(or eight different values in the octants of ${\mathbb R}^3$). The exact degree
of regularity (the $s$ in the index of the Sobolev space above) depends on the
values of $a(\mathbf x)$ coming together at the origin, and by choosing the
jumps large enough, the regularity of the solution can be made as close as
desired to $H^1$.

To implement something like this, one could replace the coefficient
function by the following (shown here only for the 2d case):
@code
template <int dim>
double coefficient (const Point<dim> &p)
{
  if ((p[0] < 0) && (p[1] < 0))           // lower left quadrant
    return 1;
  else if ((p[0] >= 0) && (p[1] < 0))     // lower right quadrant
    return 10;
  else if ((p[0] < 0) && (p[1] >= 0))     // upper left quadrant
    return 100;
  else if ((p[0] >= 0) && (p[1] >= 0))    // upper right quadrant
    return 1000;
  else
    {
      Assert(false, ExcInternalError());
      return 0;
    }
}
@endcode
(Adding the <code>Assert</code> at the end ensures that either an exception
is thrown or that the program aborts if we ever get to that point
-- which of course we shouldn't,
but this is a good way to insure yourself: we all make mistakes by
sometimes not thinking of all cases, for example by checking
for <code>p[0]</code> to be less than and greater than zero,
rather than greater-or-equal to zero, and thereby forgetting
some cases that would otherwise lead to bugs that are awkward
to find. The <code>return 0;</code> at the end is only there to
avoid compiler warnings that the function does not end in a
<code>return</code> statement -- the compiler cannot see that the
function would never actually get to that point because of the
preceding <code>Assert</code> statement.)

By playing with such cases where four or more sectors come
together and on which the coefficient has different values, one can
construct cases where the solution has singularities at the
origin. One can also see how the meshes are refined in such cases.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-6.cc"
*/
