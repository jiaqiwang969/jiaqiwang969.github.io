/**
@page step_2 The step-2 tutorial program
This tutorial depends on step-1.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Sparsity"> Sparsity </a>
        <li><a href="#Howdegreesoffreedomareenumerated"> How degrees of freedom are enumerated </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Meshgeneration">Mesh generation</a>
        <li><a href="#CreationofaDoFHandler">Creation of a DoFHandler</a>
        <li><a href="#RenumberingofDoFs">Renumbering of DoFs</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{9}

After we have created a grid in the previous example, we now show how
to define degrees of freedom on this mesh. For this example, we
will use the lowest order ($Q_1$) finite elements, for which the degrees
of freedom are associated with the vertices of the mesh. Later
examples will demonstrate higher order elements where degrees of freedom are
not necessarily associated with vertices any more, but can be associated
with edges, faces, or cells.

The term "degree of freedom" is commonly used in the finite element community
to indicate two slightly different, but related things. The first is that we'd
like to represent the finite element solution as a linear combination of shape
functions, in the form $u_h(\mathbf x) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf
x)$. Here, $U_j$ is a vector of expansion coefficients. Because we don't know
their values yet (we will compute them as the solution of a linear or
nonlinear system), they are called "unknowns" or "degrees of freedom". The
second meaning of the term can be explained as follows: A mathematical
description of finite element problems is often to say that we are looking for
a finite dimensional function $u_h \in V_h$ that satisfies some set of equations
(e.g. $a(u_h,\varphi_h)=(f,\varphi_h)$ for all test functions $\varphi_h\in
V_h$). In other words, all we say here that the solution needs to lie in some
space $V_h$. However, to actually solve this problem on a computer we need to
choose a basis of this space; this is the set of shape functions
$\varphi_j(\mathbf x)$ we have used above in the expansion of $u_h(\mathbf x)$
with coefficients $U_j$. There are of course many bases of the space $V_h$,
but we will specifically choose the one that is described by the finite
element functions that are traditionally defined locally on the cells of the
mesh. Describing "degrees of freedom" in this context requires us to simply
<i>enumerate</i> the basis functions of the space $V_h$. For $Q_1$ elements
this means simply enumerating the vertices of the mesh in some way, but for
higher order elements, one also has to enumerate the shape functions that are
associated with edges, faces, or cell interiors of the mesh. In other words,
the enumeration of degrees of freedom is an entirely separate thing from the
indices we use for vertices. The class that
provides this enumeration of the basis functions of $V_h$ is called DoFHandler.

Defining degrees of freedom ("DoF"s in short) on a mesh is a rather
simple task, since the library does all the work for you. Essentially,
all you have to do is create a finite element object (from one of the
many finite element classes deal.II already has, see for example the
@ref fe documentation) and give it to a DoFHandler object through the
DoFHandler::distribute_dofs function ("distributing DoFs" is the term we use
to describe the process of <i>enumerating</i> the basis functions as discussed
above). The DoFHandler is a class that
knows which degrees of freedom live where, i.e., it can answer
questions like "how many degrees of freedom are there globally" and
"on this cell, give me the global indices of the shape functions that
live here". This is the sort of information you need when determining
how big your system matrix should be, and when copying the
contributions of a single cell into the global matrix.

<a name="Sparsity"></a><h3> Sparsity </h3>


The next step would then be to compute a matrix and right hand side
corresponding to a particular differential equation using this finite element
and mesh. We will keep this step for the step-3 program and rather talk about
one practical aspect of a finite element program, namely that finite element
matrices are always very sparse: almost all entries in these
matrices are zero.

To be more precise, we say that a matrix is sparse
if the number of nonzero entries <i>per row</i> in the matrix is
bounded by a number that is independent of the overall number of degrees of
freedom. For example, the simple 5-point stencil of a finite difference
approximation of the Laplace equation leads to a sparse matrix since the
number of nonzero entries per row is five, and therefore independent of the
total size of the matrix. For more complicated problems -- say, the Stokes
problem of step-22 -- and in particular in 3d, the number of entries per row
may be several hundred. But the important point is that this number is
independent of the overall size of the problem: If you refine the mesh, the
maximal number of unknowns per row remains the same.

Sparsity is one of the distinguishing feature of
the finite element method compared to, say, approximating the solution of a
partial differential equation using a Taylor expansion and matching
coefficients, or using a Fourier basis.

In practical terms, it is the sparsity of matrices that enables us to solve
problems with millions or billions of unknowns. To understand this, note that
a matrix with $N$ rows, each with a fixed upper bound for the number of
nonzero entries, requires ${\cal O}(N)$ memory locations for storage, and a
matrix-vector multiplication also requires only ${\cal O}(N)$
operations. Consequently, if we had a linear solver that requires only a fixed
number of matrix-vector multiplications to come up with the solution of a
linear system with this matrix, then we would have a solver that can find the
values of all $N$ unknowns with optimal complexity, i.e., with a total of
${\cal O}(N)$ operations. It is clear that this wouldn't be possible if the
matrix were not sparse (because then the number of entries in the matrix would
have to be ${\cal O}(N^s)$ with some $s>1$, and doing a fixed number of
matrix-vector products would take ${\cal O}(N^s)$ operations),
but it also requires very specialized solvers such as
multigrid methods to satisfy the requirement that the solution requires only a
fixed number of matrix-vector multiplications. We will frequently look at the
question of what solver to use in the remaining programs of this tutorial.

The sparsity is generated by the fact that finite element shape
functions are defined <i>locally</i> on individual cells, rather than
globally, and that the local differential operators in the bilinear
form only couple shape functions whose support overlaps. (The "support" of
a function is the area where it is nonzero. For the finite element method,
the support of a shape function is generally the cells adjacent to the vertex,
edge, or face it is defined on.) In other words, degrees of freedom $i$ and $j$
that are not defined on the same cell do not overlap, and consequently
the matrix entry $A_{ij}$ will be zero.  (In some cases such
as the Discontinuous Galerkin method, shape functions may also connect
to neighboring cells through face integrals. But finite element
methods do not generally couple shape functions beyond the immediate
neighbors of a cell on which the function is defined.)


<a name="Howdegreesoffreedomareenumerated"></a><h3> How degrees of freedom are enumerated </h3>


By default, the DoFHandler class enumerates degrees of freedom on a mesh in a
rather random way; consequently, the sparsity pattern is also not
optimized for any particular purpose. To show this, the code below will
demonstrate a simple way to output the "sparsity pattern" that corresponds to
a DoFHandler, i.e., an object that represents all of the potentially nonzero
elements of a matrix one may build when discretizing a partial differential
equation on a mesh and its DoFHandler. This lack of structure in the sparsity
pattern will be apparent from the pictures we show below.

For most applications and algorithms, the exact way in which degrees of freedom
are numbered does not matter. For example, the Conjugate Gradient method we
use to solve linear systems does not care. On the other hand,
some algorithms do care: in particular, some preconditioners such as SSOR
will work better if they can walk through degrees of freedom in a particular
order, and it would be nice if we could just sort them in such a way that
SSOR can iterate through them from zero to $N$ in this order. Other examples
include computing incomplete LU or Cholesky factorizations, or if we care
about the block structure of matrices (see step-20 for an example).
deal.II therefore has algorithms that can re-enumerate degrees of freedom
in particular ways in namespace DoFRenumbering. Renumbering can be thought
of as choosing a different, permuted basis of the finite element space. The
sparsity pattern and matrices that result from this renumbering are therefore
also simply a permutation of rows and columns compared to the ones we would
get without explicit renumbering.

In the program below, we will use the algorithm of Cuthill and McKee to do
so. We will show the sparsity pattern for both the original enumeration of
degrees of freedom and of the renumbered version below,
in the <a href="#Results">results section</a>.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The first few includes are just like in the previous program, so do not
 * require additional comments:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * @endcode
 * 
 * However, the next file is new. We need this include file for the
 * association of degrees of freedom ("DoF"s) to vertices, lines, and cells:
 * 
 * @code
 * #include <deal.II/dofs/dof_handler.h>
 * 
 * @endcode
 * 
 * The following include contains the description of the bilinear finite
 * element, including the facts that it has one degree of freedom on each
 * vertex of the triangulation, but none on faces and none in the interior of
 * the cells.
 * 

 * 
 * (In fact, the file contains the description of Lagrange elements in
 * general, i.e. also the quadratic, cubic, etc versions, and not only for 2d
 * but also 1d and 3d.)
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * @endcode
 * 
 * In the following file, several tools for manipulating degrees of freedom
 * can be found:
 * 
 * @code
 * #include <deal.II/dofs/dof_tools.h>
 * @endcode
 * 
 * We will use a sparse matrix to visualize the pattern of nonzero entries
 * resulting from the distribution of degrees of freedom on the grid. That
 * class can be found here:
 * 
 * @code
 * #include <deal.II/lac/sparse_matrix.h>
 * @endcode
 * 
 * We will also need to use an intermediate sparsity pattern structure, which
 * is found in this file:
 * 
 * @code
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * @endcode
 * 
 * We will want to use a special algorithm to renumber degrees of freedom. It
 * is declared here:
 * 
 * @code
 * #include <deal.II/dofs/dof_renumbering.h>
 * 
 * @endcode
 * 
 * And this is again needed for C++ output:
 * 
 * @code
 * #include <fstream>
 * 
 * @endcode
 * 
 * Finally, as in step-1, we import the deal.II namespace into the global
 * scope:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Meshgeneration"></a> 
 * <h3>Mesh generation</h3>
 * 

 * 
 * This is the function that produced the circular grid in the previous step-1
 * example program with fewer refinements steps. The sole difference is that it
 * returns the grid it produces via its argument.
 * 
 * @code
 * void make_grid(Triangulation<2> &triangulation)
 * {
 *   const Point<2> center(1, 0);
 *   const double   inner_radius = 0.5, outer_radius = 1.0;
 *   GridGenerator::hyper_shell(
 *     triangulation, center, inner_radius, outer_radius, 5);
 * 
 *   for (unsigned int step = 0; step < 3; ++step)
 *     {
 *       for (auto &cell : triangulation.active_cell_iterators())
 *         for (const auto v : cell->vertex_indices())
 *           {
 *             const double distance_from_center =
 *               center.distance(cell->vertex(v));
 * 
 *             if (std::fabs(distance_from_center - inner_radius) <=
 *                 1e-6 * inner_radius)
 *               {
 *                 cell->set_refine_flag();
 *                 break;
 *               }
 *           }
 * 
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="CreationofaDoFHandler"></a> 
 * <h3>Creation of a DoFHandler</h3>
 * 

 * 
 * Up to now, we only have a grid, i.e. some geometrical (the position of the
 * vertices) and some topological information (how vertices are connected to
 * lines, and lines to cells, as well as which cells neighbor which other
 * cells). To use numerical algorithms, one needs some logic information in
 * addition to that: we would like to associate degree of freedom numbers to
 * each vertex (or line, or cell, in case we were using higher order elements)
 * to later generate matrices and vectors which describe a finite element
 * field on the triangulation.
 * 

 * 
 * This function shows how to do this. The object to consider is the
 * <code>DoFHandler</code> class template.  Before we do so, however, we first
 * need something that describes how many degrees of freedom are to be
 * associated to each of these objects. Since this is one aspect of the
 * definition of a finite element space, the finite element base class stores
 * this information. In the present context, we therefore create an object of
 * the derived class <code>FE_Q</code> that describes Lagrange elements. Its
 * constructor takes one argument that states the polynomial degree of the
 * element, which here is one (indicating a bi-linear element); this then
 * corresponds to one degree of freedom for each vertex, while there are none
 * on lines and inside the quadrilateral. A value of, say, three given to the
 * constructor would instead give us a bi-cubic element with one degree of
 * freedom per vertex, two per line, and four inside the cell. In general,
 * <code>FE_Q</code> denotes the family of continuous elements with complete
 * polynomials (i.e. tensor-product polynomials) up to the specified order.
 * 

 * 
 * We first need to create an object of this class and then pass it on to the
 * <code>DoFHandler</code> object to allocate storage for the degrees of
 * freedom (in deal.II lingo: we <i>distribute degrees of
 * freedom</i>).
 * 
 * @code
 * void distribute_dofs(DoFHandler<2> &dof_handler)
 * {
 *   const FE_Q<2> finite_element(1);
 *   dof_handler.distribute_dofs(finite_element);
 * 
 * @endcode
 * 
 * Now that we have associated a degree of freedom with a global number to
 * each vertex, we wonder how to visualize this?  There is no simple way to
 * directly visualize the DoF number associated with each vertex. However,
 * such information would hardly ever be truly important, since the
 * numbering itself is more or less arbitrary. There are more important
 * factors, of which we will demonstrate one in the following.
 *   

 * 
 * Associated with each vertex of the triangulation is a shape
 * function. Assume we want to solve something like Laplace's equation, then
 * the different matrix entries will be the integrals over the gradient of
 * each pair of such shape functions. Obviously, since the shape functions
 * are nonzero only on the cells adjacent to the vertex they are associated
 * with, matrix entries will be nonzero only if the supports of the shape
 * functions associated to that column and row %numbers intersect. This is
 * only the case for adjacent shape functions, and therefore only for
 * adjacent vertices. Now, since the vertices are numbered more or less
 * randomly by the above function (DoFHandler::distribute_dofs), the pattern
 * of nonzero entries in the matrix will be somewhat ragged, and we will
 * take a look at it now.
 *   

 * 
 * First we have to create a structure which we use to store the places of
 * nonzero elements. This can then later be used by one or more sparse
 * matrix objects that store the values of the entries in the locations
 * stored by this sparsity pattern. The class that stores the locations is
 * the SparsityPattern class. As it turns out, however, this class has some
 * drawbacks when we try to fill it right away: its data structures are set
 * up in such a way that we need to have an estimate for the maximal number
 * of entries we may wish to have in each row. In two space dimensions,
 * reasonable values for this estimate are available through the
 * DoFHandler::max_couplings_between_dofs() function, but in three
 * dimensions the function almost always severely overestimates the true
 * number, leading to a lot of wasted memory, sometimes too much for the
 * machine used, even if the unused memory can be released immediately after
 * computing the sparsity pattern. In order to avoid this, we use an
 * intermediate object of type DynamicSparsityPattern that uses a
 * different %internal data structure and that we can later copy into the
 * SparsityPattern object without much overhead. (Some more information on
 * these data structures can be found in the @ref Sparsity module.) In order
 * to initialize this intermediate data structure, we have to give it the
 * size of the matrix, which in our case will be square with as many rows
 * and columns as there are degrees of freedom on the grid:
 * 
 * @code
 *   DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
 *                                                   dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * We then fill this object with the places where nonzero elements will be
 * located given the present numbering of degrees of freedom:
 * 
 * @code
 *   DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
 * 
 * @endcode
 * 
 * Now we are ready to create the actual sparsity pattern that we could
 * later use for our matrix. It will just contain the data already assembled
 * in the DynamicSparsityPattern.
 * 
 * @code
 *   SparsityPattern sparsity_pattern;
 *   sparsity_pattern.copy_from(dynamic_sparsity_pattern);
 * 
 * @endcode
 * 
 * With this, we can now write the results to a file:
 * 
 * @code
 *   std::ofstream out("sparsity_pattern1.svg");
 *   sparsity_pattern.print_svg(out);
 * @endcode
 * 
 * The result is stored in an <code>.svg</code> file, where each nonzero entry
 * in the matrix corresponds with a red square in the image. The output will
 * be shown below.
 *   

 * 
 * If you look at it, you will note that the sparsity pattern is
 * symmetric. This should not come as a surprise, since we have not given
 * the <code>DoFTools::make_sparsity_pattern</code> any information that
 * would indicate that our bilinear form may couple shape functions in a
 * non-symmetric way. You will also note that it has several distinct
 * region, which stem from the fact that the numbering starts from the
 * coarsest cells and moves on to the finer ones; since they are all
 * distributed symmetrically around the origin, this shows up again in the
 * sparsity pattern.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="RenumberingofDoFs"></a> 
 * <h3>Renumbering of DoFs</h3>
 * 

 * 
 * In the sparsity pattern produced above, the nonzero entries extended quite
 * far off from the diagonal. For some algorithms, for example for incomplete
 * LU decompositions or Gauss-Seidel preconditioners, this is unfavorable, and
 * we will show a simple way how to improve this situation.
 * 

 * 
 * Remember that for an entry $(i,j)$ in the matrix to be nonzero, the
 * supports of the shape functions i and j needed to intersect (otherwise in
 * the integral, the integrand would be zero everywhere since either the one
 * or the other shape function is zero at some point). However, the supports
 * of shape functions intersected only if they were adjacent to each other, so
 * in order to have the nonzero entries clustered around the diagonal (where
 * $i$ equals $j$), we would like to have adjacent shape functions to be
 * numbered with indices (DoF numbers) that differ not too much.
 * 

 * 
 * This can be accomplished by a simple front marching algorithm, where one
 * starts at a given vertex and gives it the index zero. Then, its neighbors
 * are numbered successively, making their indices close to the original
 * one. Then, their neighbors, if not yet numbered, are numbered, and so on.
 * 

 * 
 * One algorithm that adds a little bit of sophistication along these lines is
 * the one by Cuthill and McKee. We will use it in the following function to
 * renumber the degrees of freedom such that the resulting sparsity pattern is
 * more localized around the diagonal. The only interesting part of the
 * function is the first call to <code>DoFRenumbering::Cuthill_McKee</code>,
 * the rest is essentially as before:
 * 
 * @code
 * void renumber_dofs(DoFHandler<2> &dof_handler)
 * {
 *   DoFRenumbering::Cuthill_McKee(dof_handler);
 * 
 *   DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
 *                                                   dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
 * 
 *   SparsityPattern sparsity_pattern;
 *   sparsity_pattern.copy_from(dynamic_sparsity_pattern);
 * 
 *   std::ofstream out("sparsity_pattern2.svg");
 *   sparsity_pattern.print_svg(out);
 * }
 * 
 * @endcode
 * 
 * Again, the output is shown below. Note that the nonzero entries are
 * clustered far better around the diagonal than before. This effect is even
 * more distinguished for larger matrices (the present one has 1260 rows and
 * columns, but large matrices often have several 100,000s).
 * 

 * 
 * It is worth noting that the <code>DoFRenumbering</code> class offers a
 * number of other algorithms as well to renumber degrees of freedom. For
 * example, it would of course be ideal if all couplings were in the lower or
 * upper triangular part of a matrix, since then solving the linear system
 * would amount to only forward or backward substitution. This is of course
 * unachievable for symmetric sparsity patterns, but in some special
 * situations involving transport equations, this is possible by enumerating
 * degrees of freedom from the inflow boundary along streamlines to the
 * outflow boundary. Not surprisingly, <code>DoFRenumbering</code> also has
 * algorithms for this.
 * 

 * 
 * 

 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * Finally, this is the main program. The only thing it does is to allocate
 * and create the triangulation, then create a <code>DoFHandler</code> object
 * and associate it to the triangulation, and finally call above two functions
 * on it:
 * 
 * @code
 * int main()
 * {
 *   Triangulation<2> triangulation;
 *   make_grid(triangulation);
 * 
 *   DoFHandler<2> dof_handler(triangulation);
 * 
 *   distribute_dofs(dof_handler);
 *   renumber_dofs(dof_handler);
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


The program has, after having been run, produced two sparsity
patterns. We can visualize them by opening the <code>.svg</code> files in a web browser.

The results then look like this (every point denotes an entry which
might be nonzero; of course the fact whether the entry actually is
zero or not depends on the equation under consideration, but the
indicated positions in the matrix tell us which shape functions can
and which can't couple when discretizing a local, i.e. differential,
equation):
<table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-1.svg" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-2.svg" alt=""></td>
  </tr>
</table>

The different regions in the left picture, indicated by kinks in the lines and
single dots on the left and top, represent the degrees of
freedom on the different refinement levels of the triangulation.  As
can be seen in the right picture, the sparsity pattern is much better
clustered around the main diagonal of the matrix after
renumbering. Although this might not be apparent, the number of
nonzero entries is the same in both pictures, of course.



<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


Just as with step-1, you may want to play with the program a bit to
familiarize yourself with deal.II. For example, in the
<code>distribute_dofs</code> function, we use linear finite elements
(that's what the argument "1" to the FE_Q object is). Explore how the
sparsity pattern changes if you use higher order elements, for example
cubic or quintic ones (by using 3 and 5 as the respective arguments).

You could also explore how the sparsity pattern changes by refining
the mesh. You will see that not only the size of the matrix
changes, but also its bandwidth (the distance from the diagonal of
those nonzero elements of the matrix that are farthest away from the
diagonal), though the ratio of bandwidth to size typically shrinks,
i.e. the matrix clusters more around the diagonal.

Another idea of experiments would be to try other renumbering
strategies than Cuthill-McKee from the DoFRenumbering namespace and see how
they affect the sparsity pattern.

You can also visualize the output using <a
href="http://www.gnuplot.info/">GNUPLOT</a> (one of the simpler visualization
programs; maybe not the easiest to use since it is command line driven, but
also universally available on all Linux and other Unix-like systems) by changing from <code>print_svg()</code> to <code>print_gnuplot()</code> in <code>distribute_dofs()</code> and <code>renumber_dofs()</code>:
@code
examples/\step-2> gnuplot

        G N U P L O T
        Version 3.7 patchlevel 3
        last modified Thu Dec 12 13:00:00 GMT 2002
        System: Linux 2.6.11.4-21.10-default

        Copyright(C) 1986 - 1993, 1998 - 2002
        Thomas Williams, Colin Kelley and many others

        Type `help` to access the on-line reference manual
        The gnuplot FAQ is available from
        http://www.gnuplot.info/gnuplot-faq.html

        Send comments and requests for help to <info-gnuplot@dartmouth.edu>
        Send bugs, suggestions and mods to <bug-gnuplot@dartmouth.edu>


Terminal type set to 'x11'
gnuplot> set style data points
gnuplot> plot "sparsity_pattern.1"
@endcode

Another practice based on
<a href="http://www.gnuplot.info/">GNUPLOT</a> is trying to
print out the mesh with locations and numbering of the support
points. For that, you need to include header files for GridOut and MappingQ1.
The code for this is:
@code
  std::ofstream out("gnuplot.gpl");
  out << "plot '-' using 1:2 with lines, "
      << "'-' with labels point pt 2 offset 1,1"
      << std::endl;
  GridOut().write_gnuplot (triangulation, out);
  out << "e" << std::endl;
  const int dim = 2;
  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
                                        dof_handler,
                                        support_points);
  DoFTools::write_gnuplot_dof_support_point_info(out,
                                                 support_points);
  out << "e" << std::endl;
@endcode
After we run the code, we get a file called gnuplot.gpl. To view this
file, we can run the following code in the command line:
@code
gnuplot -p gnuplot.gpl
@endcode.
With that, you will get a picture similar to
@image html support_point_dofs1.png
depending on the mesh you are looking at. For more information, see DoFTools::write_gnuplot_dof_support_point_info.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-2.cc"
*/
