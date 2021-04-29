/**
@page step_3 The step-3 tutorial program
This tutorial depends on step-2.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thebasicsetupoffiniteelementmethods">The basic set up of finite element methods</a>
        <li><a href="#Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"> Should we multiply by a test function from the left or from the right? </a>
        <li><a href="#Computingthematrixandrighthandsidevector"> Computing the matrix and right hand side vector </a>
        <li><a href="#Abouttheimplementation">About the implementation</a>
        <li><a href="#Anoteontypes"> A note on types </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Manynewincludefiles">Many new include files</a>
        <li><a href="#ThecodeStep3codeclass">The <code>Step3</code> class</a>
      <ul>
        <li><a href="#Step3Step3">Step3::Step3</a>
        <li><a href="#Step3make_grid">Step3::make_grid</a>
        <li><a href="#Step3setup_system">Step3::setup_system</a>
        <li><a href="#Step3assemble_system">Step3::assemble_system</a>
        <li><a href="#Step3solve">Step3::solve</a>
        <li><a href="#Step3output_results">Step3::output_results</a>
        <li><a href="#Step3run">Step3::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
        <li><a href="#UsingHDF5tooutputthesolutionandadditionaldata">Using HDF5 to output the solution and additional data</a>
      <ul>
        <li><a href="#Changingtheoutputtoh5"> Changing the output to .h5</a>
        <li><a href="#Addingthepointvalueandthemeanseeextensionaboveintotheh5file"> Adding the point value and the mean (see extension above) into the .h5 file</a>
      </ul>
        <li><a href="#UsingRandggplot2togenerateplots"> Using R and ggplot2 to generate plots</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{10}

<a name="Thebasicsetupoffiniteelementmethods"></a><h3>The basic set up of finite element methods</h3>


This is the first example where we actually use finite elements to compute
something. We
will solve a simple version of Poisson's equation with zero boundary
values, but a nonzero right hand side:
@f{align*}
  -\Delta u &= f \qquad\qquad & \text{in}\ \Omega,
  \\
  u &= 0 \qquad\qquad & \text{on}\ \partial\Omega.
@f}
We will solve this equation on the square, $\Omega=[-1,1]^2$, for which
you've already learned how to generate a mesh in step-1 and step-2. In
this program, we will also only consider the particular case
$f(\mathbf x)=1$ and come back to how to implement the more general
case in the next tutorial program, step-4.

If you've learned about the basics of the finite element method, you will
remember the steps we need to take to approximate the solution $u$ by a finite
dimensional approximation. Specifically, we first need to derive the weak form
of the equation above, which we obtain by multiplying the equation by a test
function $\varphi$ <i>from the left</i> (we will come back to the reason for
multiplying from the left and not from the right below) and integrating over
the domain $\Omega$:
@f{align*}
  -\int_\Omega \varphi \Delta u = \int_\Omega \varphi f.
@f}
This can be integrated by parts:
@f{align*}
  \int_\Omega \nabla\varphi \cdot \nabla u
  -
  \int_{\partial\Omega} \varphi \mathbf{n}\cdot \nabla u
   = \int_\Omega \varphi f.
@f}
The test function $\varphi$ has to satisfy the same kind of boundary
conditions (in mathematical terms: it needs to come from the tangent space of
the set in which we seek the solution), so on the boundary $\varphi=0$ and
consequently the weak form we are looking for reads
@f{align*}
  (\nabla\varphi, \nabla u)
   = (\varphi, f),
@f}
where we have used the common notation $(a,b)=\int_\Omega a\; b$. The problem
then asks for a function $u$ for which this statement is true for all test
functions $\varphi$ from the appropriate space (which here is the space
$H^1$).

Of course we can't find such a function on a computer in the general case, and
instead we seek an approximation $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf
x)$, where the $U_j$ are unknown expansion coefficients we need to determine
(the "degrees of freedom" of this problem), and $\varphi_i(\mathbf x)$ are the
finite element shape functions we will use. To define these shape functions,
we need the following:

- A mesh on which to define shape functions. You have already seen how to
  generate and manipulate the objects that describe meshes in step-1 and
  step-2.
- A finite element that describes the shape functions we want to use on the
  reference cell (which in deal.II is always the unit interval $[0,1]$, the
  unit square $[0,1]^2$ or the unit cube $[0,1]^3$, depending on which space
  dimension you work in). In step-2, we had already used an object of type
  FE_Q<2>, which denotes the usual Lagrange elements that define shape
  functions by interpolation on support points. The simplest one is
  FE_Q<2>(1), which uses polynomial degree 1. In 2d, these are often referred
  to as <i>bilinear</i>, since they are linear in each of the two coordinates
  of the reference cell. (In 1d, they would be <i>linear</i> and in 3d
  <i>tri-linear</i>; however, in the deal.II documentation, we will frequently
  not make this distinction and simply always call these functions "linear".)
- A DoFHandler object that enumerates all the degrees of freedom on the mesh,
  taking the reference cell description the finite element object provides as
  the basis. You've also already seen how to do this in step-2.
- A mapping that tells how the shape functions on the real cell are obtained
  from the shape functions defined by the finite element class on the
  reference cell. By default, unless you explicitly say otherwise, deal.II
  will use a (bi-, tri-)linear mapping for this, so in most cases you don't
  have to worry about this step.

Through these steps, we now have a set of functions $\varphi_i$, and we can
define the weak form of the discrete problem: Find a function $u_h$, i.e., find
the expansion coefficients $U_j$ mentioned above, so that
@f{align*}
  (\nabla\varphi_i, \nabla u_h)
   = (\varphi_i, f),
   \qquad\qquad
   i=0\ldots N-1.
@f}
Note that we here follow the convention that everything is counted starting at
zero, as common in C and C++. This equation can be rewritten as a linear
system if you insert the representation $u_h(\mathbf x)=\sum_j U_j
\varphi_j(\mathbf x)$ and then observe that
@f{align*}{
  (\nabla\varphi_i, \nabla u_h)
  &= \left(\nabla\varphi_i, \nabla \Bigl[\sum_j U_j \varphi_j\Bigr]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \left[U_j \varphi_j\right]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \varphi_j \right) U_j.
@f}
With this, the problem reads: Find a vector $U$ so that
@f{align*}{
  A U = F,
@f}
where the matrix $A$ and the right hand side $F$ are defined as
@f{align*}
  A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j),
  \\
  F_i &= (\varphi_i, f).
@f}


<a name="Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"></a><h3> Should we multiply by a test function from the left or from the right? </h3>


Before we move on with describing how these quantities can be computed, note
that if we had multiplied the original equation from the <i>right</i> by a
test function rather than from the left, then we would have obtained a linear
system of the form
@f{align*}
  U^T A = F^T
@f}
with a row vector $F^T$. By transposing this system, this is of course
equivalent to solving
@f{align*}
  A^T U = F
@f}
which here is the same as above since $A=A^T$. But in general is not,
and in order to avoid
any sort of confusion, experience has shown that simply getting into the habit
of multiplying the equation from the left rather than from the right (as is
often done in the mathematical literature) avoids a common class of errors as
the matrix is automatically correct and does not need to be transposed when
comparing theory and implementation. See step-9 for the first example in this
tutorial where we have a non-symmetric bilinear form for which it makes a
difference whether we multiply from the right or from the left.


<a name="Computingthematrixandrighthandsidevector"></a><h3> Computing the matrix and right hand side vector </h3>


Now we know what we need (namely: objects that hold the matrix and
vectors, as well as ways to compute $A_{ij},F_i$), and we can look at what it
takes to make that happen:

- The object for $A$ is of type SparseMatrix while those for $U$ and $F$ are of
  type Vector. We will see in the program below what classes are used to solve
  linear systems.
- We need a way to form the integrals. In the finite element method, this is
  most commonly done using quadrature, i.e. the integrals are replaced by a
  weighted sum over a set of points on each cell. That is, we first split the
  integral over $\Omega$ into integrals over all cells,
  @f{align*}
    A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j)
    = \sum_{K \in {\mathbb T}} \int_K \nabla\varphi_i \cdot \nabla \varphi_j,
    \\
    F_i &= (\varphi_i, f)
    = \sum_{K \in {\mathbb T}} \int_K \varphi_i f,
  @f}
  and then approximate each cell's contribution by quadrature:
  @f{align*}
    A^K_{ij} &=
    \int_K \nabla\varphi_i \cdot \nabla \varphi_j
    \approx
    \sum_q \nabla\varphi_i(\mathbf x^K_q) \cdot \nabla
    \varphi_j(\mathbf x^K_q) w_q^K,
    \\
    F^K_i &=
    \int_K \varphi_i f
    \approx
    \sum_q \varphi_i(\mathbf x^K_q) f(\mathbf x^K_q) w^K_q,
  @f}
  where $\mathbf x^K_q$ is the $q$th quadrature point on cell $K$, and $w^K_q$
  the $q$th quadrature weight. There are different parts to what is needed in
  doing this, and we will discuss them in turn next.
- First, we need a way to describe the location $\mathbf x_q^K$ of quadrature
  points and their weights $w^K_q$. They are usually mapped from the reference
  cell in the same way as shape functions, i.e., implicitly using the
  MappingQ1 class or, if you explicitly say so, through one of the other
  classes derived from Mapping. The locations and weights on the reference
  cell are described by objects derived from the Quadrature base
  class. Typically, one chooses a quadrature formula (i.e. a set of points and
  weights) so that the quadrature exactly equals the integral in the matrix;
  this can be achieved because all factors in the integral are polynomial, and
  is done by Gaussian quadrature formulas, implemented in the QGauss class.
- We then need something that can help us evaluate $\varphi_i(\mathbf x^K_q)$
  on cell $K$. This is what the FEValues class does: it takes a finite element
  objects to describe $\varphi$ on the reference cell, a quadrature object to
  describe the quadrature points and weights, and a mapping object (or
  implicitly takes the MappingQ1 class) and provides values and derivatives of
  the shape functions on the real cell $K$ as well as all sorts of other
  information needed for integration, at the quadrature points located on $K$.

FEValues really is the central class in the assembly process. One way you can
view it is as follows: The FiniteElement and derived classes describe shape
<i>functions</i>, i.e., infinite dimensional objects: functions have values at
every point. We need this for theoretical reasons because we want to perform
our analysis with integrals over functions. However, for a computer, this is a
very difficult concept, since they can in general only deal with a finite
amount of information, and so we replace integrals by sums over quadrature
points that we obtain by mapping (the Mapping object) using  points defined on
a reference cell (the Quadrature object) onto points on the real cell. In
essence, we reduce the problem to one where we only need a finite amount of
information, namely shape function values and derivatives, quadrature weights,
normal vectors, etc, exclusively at a finite set of points. The FEValues class
is the one that brings the three components together and provides this finite
set of information on a particular cell $K$. You will see it in action when we
assemble the linear system below.

It is noteworthy that all of this could also be achieved if you simply created
these three objects yourself in an application program, and juggled the
information yourself. However, this would neither be simpler (the FEValues
class provides exactly the kind of information you actually need) nor faster:
the FEValues class is highly optimized to only compute on each cell the
particular information you need; if anything can be re-used from the previous
cell, then it will do so, and there is a lot of code in that class to make
sure things are cached wherever this is advantageous.

The final piece of this introduction is to mention that after a linear
system is obtained, it is solved using an iterative solver and then
postprocessed: we create an output file using the DataOut class that can then
be visualized using one of the common visualization programs.

@note The preceding overview of all the important steps of any finite element
implementation has its counterpart in deal.II: The library can naturally be
grouped into a number of "modules" that cover the basic concepts just
outlined. You can access these modules through the tab at the top of this
page. An overview of the most fundamental groups of concepts is also available
on the <a href="index.html">front page of the deal.II manual</a>.


<a name="Abouttheimplementation"></a><h3>About the implementation</h3>


Although this is the simplest possible equation you can solve using the finite
element method, this program shows the basic structure of most finite
element programs and also serves as the template that almost all of the
following programs will essentially follow. Specifically, the main class of
this program looks like this:
@code
class Step3
{
  public:
    Step3 ();
    void run ();

  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;

    Triangulation<2>     triangulation;
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};
@endcode

This follows the object oriented programming mantra of <a
href="http://en.wikipedia.org/wiki/Encapsulation_(object-oriented_programming)">data
encapsulation</a>, i.e. we do our best to hide almost all internal details of
this class in private members that are not accessible to the outside.

Let's start with the member variables: These follow the building blocks we
have outlined above in the bullet points, namely we need a Triangulation and a
DoFHandler object, and a finite element object that describes the kinds of
shape functions we want to use. The second group of objects relate to the
linear algebra: the system matrix and right hand side as well as the solution
vector, and an object that describes the sparsity pattern of the matrix. This
is all this class needs (and the essentials that any solver for a stationary
PDE requires) and that needs to survive throughout the entire program. In
contrast to this, the FEValues object we need for assembly is only required
throughout assembly, and so we create it as a local object in the function
that does that and destroy it again at its end.

Secondly, let's look at the member functions. These, as well, already form the
common structure that almost all following tutorial programs will use:
<ul>
  <li> <code>make_grid()</code>: This is what one could call a
       <i>preprocessing function</i>. As its name suggests, it sets up the
       object that stores the triangulation. In later examples, it could also
       deal with boundary conditions, geometries, etc.
  <li> <code>setup_system()</code>: This then is the function in which all the
       other data structures are set up that are needed to solve the
       problem. In particular, it will initialize the DoFHandler object and
       correctly size the various objects that have to do with the linear
       algebra. This function is often separated from the preprocessing
       function above because, in a time dependent program, it may be called
       at least every few time steps whenever the mesh
       is adaptively refined (something we will see how to do in step-6). On
       the other hand, setting up the mesh itself in the preprocessing
       function above is done only once at the beginning of the program and
       is, therefore, separated into its own function.
  <li> <code>assemble_system()</code>: This, then is where the contents of the
       matrix and right hand side are computed, as discussed at length in the
       introduction above. Since doing something with this linear system is
       conceptually very different from computing its entries, we separate it
       from the following function.
  <li> <code>solve()</code>: This then is the function in which we compute the
       solution $U$ of the linear system $AU=F$. In the current program, this
       is a simple task since the matrix is so simple, but it will become a
       significant part of a program's size whenever the problem is not so
       trivial any more (see, for example, step-20, step-22, or step-31 once
       you've learned a bit more about the library).
  <li> <code>output_results()</code>: Finally, when you have computed a
       solution, you probably want to do something with it. For example, you
       may want to output it in a format that can be visualized, or you may
       want to compute quantities you are interested in: say, heat fluxes in a
       heat exchanger, air friction coefficients of a wing, maximum bridge
       loads, or simply the value of the numerical solution at a point. This
       function is therefore the place for postprocessing your solution.
</ul>
All of this is held together by the single public function (other than the
constructor), namely the <code>run()</code> function. It is the one that is
called from the place where an object of this type is created, and it is the
one that calls all the other functions in their proper order. Encapsulating
this operation into the <code>run()</code> function, rather than calling all
the other functions from <code>main()</code> makes sure that you
can change how the separation of concerns within this class is
implemented. For example, if one of the functions becomes too big, you can
split it up into two, and the only places you have to be concerned about
changing as a consequence are within this very same class, and not anywhere
else.

As mentioned above, you will see this general structure &mdash; sometimes with
variants in spelling of the functions' names, but in essentially this order of
separation of functionality &mdash; again in many of the
following tutorial programs.


<a name="Anoteontypes"></a><h3> A note on types </h3>


deal.II defines a number of integral %types via alias in namespace dealii::types.
(In the previous sentence, the word "integral" is used as the <i>adjective</i>
that corresponds to the noun "integer". It shouldn't be confused with the
<i>noun</i> "integral" that represents the area or volume under a curve
or surface. The adjective "integral" is widely used in the C++ world in
contexts such as "integral type", "integral constant", etc.)
In particular, in this program you will see types::global_dof_index in a couple of
places: an integer type that is used to denote the <i>global</i> index of a
degree of freedom, i.e., the index of a particular degree of freedom within the
DoFHandler object that is defined on top of a triangulation (as opposed to the
index of a particular degree of freedom within a particular cell). For the
current program (as well as almost all of the tutorial programs), you will have
a few thousand to maybe a few million unknowns globally (and, for $Q_1$
elements, you will have 4 <i>locally on each cell</i> in 2d and 8 in 3d).
Consequently, a data type that allows to store sufficiently large numbers for
global DoF indices is <code>unsigned int</code> given that it allows to store
numbers between 0 and slightly more than 4 billion (on most systems, where
integers are 32-bit). In fact, this is what types::global_dof_index is.

So, why not just use <code>unsigned int</code> right away? deal.II used to do
this until version 7.3. However, deal.II supports very large computations (via
the framework discussed in step-40) that may have more than 4 billion unknowns
when spread across a few thousand processors. Consequently, there are
situations where <code>unsigned int</code> is not sufficiently large and we
need a 64-bit unsigned integral type. To make this possible, we introduced
types::global_dof_index which by default is defined as simply <code>unsigned
int</code> whereas it is possible to define it as <code>unsigned long long
int</code> if necessary, by passing a particular flag during configuration
(see the ReadMe file).

This covers the technical aspect. But there is also a documentation purpose:
everywhere in the library and codes that are built on it, if you see a place
using the data type types::global_dof_index, you immediately know that the
quantity that is being referenced is, in fact, a global dof index. No such
meaning would be apparent if we had just used <code>unsigned int</code> (which
may also be a local index, a boundary indicator, a material id,
etc.). Immediately knowing what a variable refers to also helps avoid errors:
it's quite clear that there must be a bug if you see an object of type
types::global_dof_index being assigned to variable of type
types::subdomain_id, even though they are both represented by unsigned
integers and the compiler will, consequently, not complain.

In more practical terms what the presence of this type means is that during
assembly, we create a $4\times 4$ matrix (in 2d, using a $Q_1$ element) of the
contributions of the cell we are currently sitting on, and then we need to add
the elements of this matrix to the appropriate elements of the global (system)
matrix. For this, we need to get at the global indices of the degrees of
freedom that are local to the current cell, for which we will always use the
following piece of the code:
@code
  cell->get_dof_indices (local_dof_indices);
@endcode
where <code>local_dof_indices</code> is declared as
@code
  std::vector<types::global_dof_index> local_dof_indices (fe.n_dofs_per_cell());
@endcode
The name of this variable might be a bit of a misnomer -- it stands for "the
global indices of those degrees of freedom locally defined on the current
cell" -- but variables that hold this information are universally named this
way throughout the library.

@note types::global_dof_index is not the only type defined in this namespace.
Rather, there is a whole family, including types::subdomain_id,
types::boundary_id, and types::material_id. All of these are alias for integer
data types but, as explained above, they are used throughout the library so that
(i) the intent of a variable becomes more easily discerned, and (ii) so that it
becomes possible to change the actual type to a larger one if necessary without
having to go through the entire library and figure out whether a particular use
of <code>unsigned int</code> corresponds to, say, a material indicator.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Manynewincludefiles"></a> 
 * <h3>Many new include files</h3>
 * 

 * 
 * These include files are already known to you. They declare the classes
 * which handle triangulations and enumeration of degrees of freedom:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/dofs/dof_handler.h>
 * @endcode
 * 
 * And this is the file in which the functions are declared that create grids:
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h>
 * 
 * @endcode
 * 
 * This file contains the description of the Lagrange interpolation finite
 * element:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * 
 * @endcode
 * 
 * And this file is needed for the creation of sparsity patterns of sparse
 * matrices, as shown in previous examples:
 * 
 * @code
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * The next two files are needed for assembling the matrix using quadrature on
 * each cell. The classes declared in them will be explained below:
 * 
 * @code
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * @endcode
 * 
 * The following three include files we need for the treatment of boundary
 * values:
 * 
 * @code
 * #include <deal.II/base/function.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * @endcode
 * 
 * We're now almost to the end. The second to last group of include files is
 * for the linear algebra which we employ to solve the system of equations
 * arising from the finite element discretization of the Laplace equation. We
 * will use vectors and full matrices for assembling the system of equations
 * locally on each cell, and transfer the results into a sparse matrix. We
 * will then use a Conjugate Gradient solver to solve the problem, for which
 * we need a preconditioner (in this program, we use the identity
 * preconditioner which does nothing, but we need to include the file anyway):
 * 
 * @code
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * @endcode
 * 
 * Finally, this is for output to a file and to the console:
 * 
 * @code
 * #include <deal.II/numerics/data_out.h>
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * ...and this is to import the deal.II namespace into the global scope:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep3codeclass"></a> 
 * <h3>The <code>Step3</code> class</h3>
 * 

 * 
 * Instead of the procedural programming of previous examples, we encapsulate
 * everything into a class for this program. The class consists of functions
 * which each perform certain aspects of a finite element program, a `main`
 * function which controls what is done first and what is done next, and a
 * list of member variables.
 * 

 * 
 * The public part of the class is rather short: it has a constructor and a
 * function `run` that is called from the outside and acts as something like
 * the `main` function: it coordinates which operations of this class shall be
 * run in which order. Everything else in the class, i.e. all the functions
 * that actually do anything, are in the private section of the class:
 * 
 * @code
 * class Step3
 * {
 * public:
 *   Step3();
 * 
 *   void run();
 * 
 * @endcode
 * 
 * Then there are the member functions that mostly do what their names
 * suggest and whose have been discussed in the introduction already. Since
 * they do not need to be called from outside, they are made private to this
 * class.
 * 

 * 
 * 
 * @code
 * private:
 *   void make_grid();
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void output_results() const;
 * 
 * @endcode
 * 
 * And finally we have some member variables. There are variables describing
 * the triangulation and the global numbering of the degrees of freedom (we
 * will specify the exact polynomial degree of the finite element in the
 * constructor of this class)...
 * 
 * @code
 *   Triangulation<2> triangulation;
 *   FE_Q<2>          fe;
 *   DoFHandler<2>    dof_handler;
 * 
 * @endcode
 * 
 * ...variables for the sparsity pattern and values of the system matrix
 * resulting from the discretization of the Laplace equation...
 * 
 * @code
 *   SparsityPattern      sparsity_pattern;
 *   SparseMatrix<double> system_matrix;
 * 
 * @endcode
 * 
 * ...and variables which will hold the right hand side and solution
 * vectors.
 * 
 * @code
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="Step3Step3"></a> 
 * <h4>Step3::Step3</h4>
 * 

 * 
 * Here comes the constructor. It does not much more than first to specify
 * that we want bi-linear elements (denoted by the parameter to the finite
 * element object, which indicates the polynomial degree), and to associate
 * the dof_handler variable to the triangulation we use. (Note that the
 * triangulation isn't set up with a mesh at all at the present time, but the
 * DoFHandler doesn't care: it only wants to know which triangulation it will
 * be associated with, and it only starts to care about an actual mesh once
 * you try to distribute degree of freedom on the mesh using the
 * distribute_dofs() function.) All the other member variables of the Step3
 * class have a default constructor which does all we want.
 * 
 * @code
 * Step3::Step3()
 *   : fe(1)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3make_grid"></a> 
 * <h4>Step3::make_grid</h4>
 * 

 * 
 * Now, the first thing we've got to do is to generate the triangulation on
 * which we would like to do our computation and number each vertex with a
 * degree of freedom. We have seen these two steps in step-1 and step-2
 * before, respectively.
 * 

 * 
 * This function does the first part, creating the mesh.  We create the grid
 * and refine all cells five times. Since the initial grid (which is the
 * square $[-1,1] \times [-1,1]$) consists of only one cell, the final grid
 * has 32 times 32 cells, for a total of 1024.
 * 

 * 
 * Unsure that 1024 is the correct number? We can check that by outputting the
 * number of cells using the <code>n_active_cells()</code> function on the
 * triangulation.
 * 
 * @code
 * void Step3::make_grid()
 * {
 *   GridGenerator::hyper_cube(triangulation, -1, 1);
 *   triangulation.refine_global(5);
 * 
 *   std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *             << std::endl;
 * }
 * 
 * @endcode
 * 
 * @note We call the Triangulation::n_active_cells() function, rather than
 * Triangulation::n_cells(). Here, <i>active</i> means the cells that aren't
 * refined any further. We stress the adjective "active" since there are more
 * cells, namely the parent cells of the finest cells, their parents, etc, up
 * to the one cell which made up the initial grid. Of course, on the next
 * coarser level, the number of cells is one quarter that of the cells on the
 * finest level, i.e. 256, then 64, 16, 4, and 1. If you called
 * <code>triangulation.n_cells()</code> instead in the code above, you would
 * consequently get a value of 1365 instead. On the other hand, the number of
 * cells (as opposed to the number of active cells) is not typically of much
 * interest, so there is no good reason to print it.
 * 

 * 
 * 

 * 
 * 
 * <a name="Step3setup_system"></a> 
 * <h4>Step3::setup_system</h4>
 * 

 * 
 * Next we enumerate all the degrees of freedom and set up matrix and vector
 * objects to hold the system data. Enumerating is done by using
 * DoFHandler::distribute_dofs(), as we have seen in the step-2 example. Since
 * we use the FE_Q class and have set the polynomial degree to 1 in the
 * constructor, i.e. bilinear elements, this associates one degree of freedom
 * with each vertex. While we're at generating output, let us also take a look
 * at how many degrees of freedom are generated:
 * 
 * @code
 * void Step3::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 *   std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 * @endcode
 * 
 * There should be one DoF for each vertex. Since we have a 32 times 32
 * grid, the number of DoFs should be 33 times 33, or 1089.
 * 

 * 
 * As we have seen in the previous example, we set up a sparsity pattern by
 * first creating a temporary structure, tagging those entries that might be
 * nonzero, and then copying the data over to the SparsityPattern object
 * that can then be used by the system matrix.
 * 
 * @code
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *   sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * Note that the SparsityPattern object does not hold the values of the
 * matrix, it only stores the places where entries are. The entries
 * themselves are stored in objects of type SparseMatrix, of which our
 * variable system_matrix is one.
 *   

 * 
 * The distinction between sparsity pattern and matrix was made to allow
 * several matrices to use the same sparsity pattern. This may not seem
 * relevant here, but when you consider the size which matrices can have,
 * and that it may take some time to build the sparsity pattern, this
 * becomes important in large-scale problems if you have to store several
 * matrices in your program.
 * 
 * @code
 *   system_matrix.reinit(sparsity_pattern);
 * 
 * @endcode
 * 
 * The last thing to do in this function is to set the sizes of the right
 * hand side vector and the solution vector to the right values:
 * 
 * @code
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Step3assemble_system"></a> 
 * <h4>Step3::assemble_system</h4>
 * 

 * 
 * 

 * 
 * The next step is to compute the entries of the matrix and right hand side
 * that form the linear system from which we compute the solution. This is the
 * central function of each finite element program and we have discussed the
 * primary steps in the introduction already.
 * 

 * 
 * The general approach to assemble matrices and vectors is to loop over all
 * cells, and on each cell compute the contribution of that cell to the global
 * matrix and right hand side by quadrature. The point to realize now is that
 * we need the values of the shape functions at the locations of quadrature
 * points on the real cell. However, both the finite element shape functions
 * as well as the quadrature points are only defined on the reference
 * cell. They are therefore of little help to us, and we will in fact hardly
 * ever query information about finite element shape functions or quadrature
 * points from these objects directly.
 * 

 * 
 * Rather, what is required is a way to map this data from the reference cell
 * to the real cell. Classes that can do that are derived from the Mapping
 * class, though one again often does not have to deal with them directly:
 * many functions in the library can take a mapping object as argument, but
 * when it is omitted they simply resort to the standard bilinear Q1
 * mapping. We will go this route, and not bother with it for the moment (we
 * come back to this in step-10, step-11, and step-12).
 * 

 * 
 * So what we now have is a collection of three classes to deal with: finite
 * element, quadrature, and mapping objects. That's too much, so there is one
 * type of class that orchestrates information exchange between these three:
 * the FEValues class. If given one instance of each three of these objects
 * (or two, and an implicit linear mapping), it will be able to provide you
 * with information about values and gradients of shape functions at
 * quadrature points on a real cell.
 * 

 * 
 * Using all this, we will assemble the linear system for this problem in the
 * following function:
 * 
 * @code
 * void Step3::assemble_system()
 * {
 * @endcode
 * 
 * Ok, let's start: we need a quadrature formula for the evaluation of the
 * integrals on each cell. Let's take a Gauss formula with two quadrature
 * points in each direction, i.e. a total of four points since we are in
 * 2D. This quadrature formula integrates polynomials of degrees up to three
 * exactly (in 1D). It is easy to check that this is sufficient for the
 * present problem:
 * 
 * @code
 *   QGauss<2> quadrature_formula(fe.degree + 1);
 * @endcode
 * 
 * And we initialize the object which we have briefly talked about above. It
 * needs to be told which finite element we want to use, and the quadrature
 * points and their weights (jointly described by a Quadrature object). As
 * mentioned, we use the implied Q1 mapping, rather than specifying one
 * ourselves explicitly. Finally, we have to tell it what we want it to
 * compute on each cell: we need the values of the shape functions at the
 * quadrature points (for the right hand side $(\varphi_i,f)$), their
 * gradients (for the matrix entries $(\nabla \varphi_i, \nabla
 * \varphi_j)$), and also the weights of the quadrature points and the
 * determinants of the Jacobian transformations from the reference cell to
 * the real cells.
 *   

 * 
 * This list of what kind of information we actually need is given as a
 * collection of flags as the third argument to the constructor of
 * FEValues. Since these values have to be recomputed, or updated, every
 * time we go to a new cell, all of these flags start with the prefix
 * <code>update_</code> and then indicate what it actually is that we want
 * updated. The flag to give if we want the values of the shape functions
 * computed is #update_values; for the gradients it is
 * #update_gradients. The determinants of the Jacobians and the quadrature
 * weights are always used together, so only the products (Jacobians times
 * weights, or short <code>JxW</code>) are computed; since we need them, we
 * have to list #update_JxW_values as well:
 * 
 * @code
 *   FEValues<2> fe_values(fe,
 *                         quadrature_formula,
 *                         update_values | update_gradients | update_JxW_values);
 * @endcode
 * 
 * The advantage of this approach is that we can specify what kind of
 * information we actually need on each cell. It is easily understandable
 * that this approach can significantly speed up finite element computations,
 * compared to approaches where everything, including second derivatives,
 * normal vectors to cells, etc are computed on each cell, regardless of
 * whether they are needed or not.
 *   

 * 
 * @note The syntax <code>update_values | update_gradients |
 * update_JxW_values</code> is not immediately obvious to anyone not
 * used to programming bit operations in C for years already. First,
 * <code>operator|</code> is the <i>bitwise or operator</i>, i.e.,
 * it takes two integer arguments that are interpreted as bit
 * patterns and returns an integer in which every bit is set for
 * which the corresponding bit is set in at least one of the two
 * arguments. For example, consider the operation
 * <code>9|10</code>. In binary, <code>9=0b1001</code> (where the
 * prefix <code>0b</code> indicates that the number is to be
 * interpreted as a binary number) and <code>10=0b1010</code>. Going
 * through each bit and seeing whether it is set in one of the
 * argument, we arrive at <code>0b1001|0b1010=0b1011</code> or, in
 * decimal notation, <code>9|10=11</code>. The second piece of
 * information you need to know is that the various
 * <code>update_*</code> flags are all integers that have <i>exactly
 * one bit set</i>. For example, assume that
 * <code>update_values=0b00001=1</code>,
 * <code>update_gradients=0b00010=2</code>,
 * <code>update_JxW_values=0b10000=16</code>. Then
 * <code>update_values | update_gradients | update_JxW_values =
 * 0b10011 = 19</code>. In other words, we obtain a number that
 * <i>encodes a binary mask representing all of the operations you
 * want to happen</i>, where each operation corresponds to exactly
 * one bit in the integer that, if equal to one, means that a
 * particular piece should be updated on each cell and, if it is
 * zero, means that we need not compute it. In other words, even
 * though <code>operator|</code> is the <i>bitwise OR operation</i>,
 * what it really represents is <i>I want this AND that AND the
 * other</i>. Such binary masks are quite common in C programming,
 * but maybe not so in higher level languages like C++, but serve
 * the current purpose quite well.
 * 

 * 
 * For use further down below, we define a shortcut for a value that will
 * be used very frequently. Namely, an abbreviation for the number of degrees
 * of freedom on each cell (since we are in 2D and degrees of freedom are
 * associated with vertices only, this number is four, but we rather want to
 * write the definition of this variable in a way that does not preclude us
 * from later choosing a different finite element that has a different
 * number of degrees of freedom per cell, or work in a different space
 * dimension).
 *   

 * 
 * In general, it is a good idea to use a symbolic name instead of
 * hard-coding these numbers even if you know them, since for example,
 * you may want to change the finite element at some time. Changing the
 * element would have to be done in a different function and it is easy
 * to forget to make a corresponding change in another part of the program.
 * It is better to not rely on your own calculations, but instead ask
 * the right object for the information: Here, we ask the finite element
 * to tell us about the number of degrees of freedom per cell and we
 * will get the correct number regardless of the space dimension or
 * polynomial degree we may have chosen elsewhere in the program.
 *   

 * 
 * The shortcut here, defined primarily to discuss the basic concept
 * and not because it saves a lot of typing, will then make the following
 * loops a bit more readable. You will see such shortcuts in many places in
 * larger programs, and `dofs_per_cell` is one that is more or less the
 * conventional name for this kind of object.
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 * @endcode
 * 
 * Now, we said that we wanted to assemble the global matrix and vector
 * cell-by-cell. We could write the results directly into the global matrix,
 * but this is not very efficient since access to the elements of a sparse
 * matrix is slow. Rather, we first compute the contribution of each cell in
 * a small matrix with the degrees of freedom on the present cell, and only
 * transfer them to the global matrix when the computations are finished for
 * this cell. We do the same for the right hand side vector. So let's first
 * allocate these objects (these being local objects, all degrees of freedom
 * are coupling with all others, and we should use a full matrix object
 * rather than a sparse one for the local operations; everything will be
 * transferred to a global sparse matrix later on):
 * 
 * @code
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 * @endcode
 * 
 * When assembling the contributions of each cell, we do this with the local
 * numbering of the degrees of freedom (i.e. the number running from zero
 * through dofs_per_cell-1). However, when we transfer the result into the
 * global matrix, we have to know the global numbers of the degrees of
 * freedom. When we query them, we need a scratch (temporary) array for
 * these numbers (see the discussion at the end of the introduction for
 * the type, types::global_dof_index, used here):
 * 
 * @code
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Now for the loop over all cells. We have seen before how this works for a
 * triangulation. A DoFHandler has cell iterators that are exactly analogous
 * to those of a Triangulation, but with extra information about the degrees
 * of freedom for the finite element you're using. Looping over the active
 * cells of a degree-of-freedom handler works the same as for a triangulation.
 *   

 * 
 * Note that we declare the type of the cell as `const auto &` instead of
 * `auto` this time around. In step 1, we were modifying the cells of the
 * triangulation by flagging them with refinement indicators. Here we're only
 * examining the cells without modifying them, so it's good practice to
 * declare `cell` as `const` in order to enforce this invariant.
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 * @endcode
 * 
 * We are now sitting on one cell, and we would like the values and
 * gradients of the shape functions be computed, as well as the
 * determinants of the Jacobian matrices of the mapping between
 * reference cell and true cell, at the quadrature points. Since all
 * these values depend on the geometry of the cell, we have to have the
 * FEValues object re-compute them on each cell:
 * 
 * @code
 *       fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Next, reset the local cell's contributions to global matrix and
 * global right hand side to zero, before we fill them:
 * 
 * @code
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 * @endcode
 * 
 * Now it is time to start integration over the cell, which we
 * do by looping over all quadrature points, which we will
 * number by q_index.
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         {
 * @endcode
 * 
 * First assemble the matrix: For the Laplace problem, the
 * matrix on each cell is the integral over the gradients of
 * shape function i and j. Since we do not integrate, but
 * rather use quadrature, this is the sum over all
 * quadrature points of the integrands times the determinant
 * of the Jacobian matrix at the quadrature point times the
 * weight of this quadrature point. You can get the gradient
 * of shape function $i$ at quadrature point with number q_index by
 * using <code>fe_values.shape_grad(i,q_index)</code>; this
 * gradient is a 2-dimensional vector (in fact it is of type
 * Tensor@<1,dim@>, with here dim=2) and the product of two
 * such vectors is the scalar product, i.e. the product of
 * the two shape_grad function calls is the dot
 * product. This is in turn multiplied by the Jacobian
 * determinant and the quadrature point weight (that one
 * gets together by the call to FEValues::JxW() ). Finally,
 * this is repeated for all shape functions $i$ and $j$:
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices())
 *             for (const unsigned int j : fe_values.dof_indices())
 *               cell_matrix(i, j) +=
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                  fe_values.JxW(q_index));           // dx
 * 
 * @endcode
 * 
 * We then do the same thing for the right hand side. Here,
 * the integral is over the shape function i times the right
 * hand side function, which we choose to be the function
 * with constant value one (more interesting examples will
 * be considered in the following programs).
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices())
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                             1. *                                // f(x_q)
 *                             fe_values.JxW(q_index));            // dx
 *         }
 * @endcode
 * 
 * Now that we have the contribution of this cell, we have to transfer
 * it to the global matrix and right hand side. To this end, we first
 * have to find out which global numbers the degrees of freedom on this
 * cell have. Let's simply ask the cell for that information:
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices);
 * 
 * @endcode
 * 
 * Then again loop over all shape functions i and j and transfer the
 * local elements to the global matrix. The global numbers can be
 * obtained using local_dof_indices[i]:
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices())
 *         for (const unsigned int j : fe_values.dof_indices())
 *           system_matrix.add(local_dof_indices[i],
 *                             local_dof_indices[j],
 *                             cell_matrix(i, j));
 * 
 * @endcode
 * 
 * And again, we do the same thing for the right hand side vector.
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices())
 *         system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *     }
 * 
 * 
 * @endcode
 * 
 * Now almost everything is set up for the solution of the discrete
 * system. However, we have not yet taken care of boundary values (in fact,
 * Laplace's equation without Dirichlet boundary values is not even uniquely
 * solvable, since you can add an arbitrary constant to the discrete
 * solution). We therefore have to do something about the situation.
 *   

 * 
 * For this, we first obtain a list of the degrees of freedom on the
 * boundary and the value the shape function shall have there. For
 * simplicity, we only interpolate the boundary value function, rather than
 * projecting it onto the boundary. There is a function in the library which
 * does exactly this: VectorTools::interpolate_boundary_values(). Its
 * parameters are (omitting parameters for which default values exist and
 * that we don't care about): the DoFHandler object to get the global
 * numbers of the degrees of freedom on the boundary; the component of the
 * boundary where the boundary values shall be interpolated; the boundary
 * value function itself; and the output object.
 *   

 * 
 * The component of the boundary is meant as follows: in many cases, you may
 * want to impose certain boundary values only on parts of the boundary. For
 * example, you may have inflow and outflow boundaries in fluid dynamics, or
 * clamped and free parts of bodies in deformation computations of
 * bodies. Then you will want to denote these different parts of the
 * boundary by indicators, and tell the interpolate_boundary_values
 * function to only compute the boundary values on a certain part of the
 * boundary (e.g. the clamped part, or the inflow boundary). By default,
 * all boundaries have a 0 boundary indicator, unless otherwise specified. If
 * sections of the boundary have different boundary conditions, you have to
 * number those parts with different boundary indicators. The function call
 * below will then only determine boundary values for those parts of the
 * boundary for which the boundary indicator is in fact the zero specified as
 * the second argument.
 *   

 * 
 * The function describing the boundary values is an object of type Function
 * or of a derived class. One of the derived classes is
 * Functions::ZeroFunction, which describes (not unexpectedly) a function
 * which is zero everywhere. We create such an object in-place and pass it to
 * the VectorTools::interpolate_boundary_values() function.
 *   

 * 
 * Finally, the output object is a list of pairs of global degree of freedom
 * numbers (i.e. the number of the degrees of freedom on the boundary) and
 * their boundary values (which are zero here for all entries). This mapping
 * of DoF numbers to boundary values is done by the <code>std::map</code>
 * class.
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values;
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<2>(),
 *                                            boundary_values);
 * @endcode
 * 
 * Now that we got the list of boundary DoFs and their respective boundary
 * values, let's use them to modify the system of equations
 * accordingly. This is done by the following function call:
 * 
 * @code
 *   MatrixTools::apply_boundary_values(boundary_values,
 *                                      system_matrix,
 *                                      solution,
 *                                      system_rhs);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3solve"></a> 
 * <h4>Step3::solve</h4>
 * 

 * 
 * The following function simply solves the discretized equation. As the
 * system is quite a large one for direct solvers such as Gauss elimination or
 * LU decomposition, we use a Conjugate Gradient algorithm. You should
 * remember that the number of variables here (only 1089) is a very small
 * number for finite element computations, where 100.000 is a more usual
 * number.  For this number of variables, direct methods are no longer usable
 * and you are forced to use methods like CG.
 * 
 * @code
 * void Step3::solve()
 * {
 * @endcode
 * 
 * First, we need to have an object that knows how to tell the CG algorithm
 * when to stop. This is done by using a SolverControl object, and as
 * stopping criterion we say: stop after a maximum of 1000 iterations (which
 * is far more than is needed for 1089 variables; see the results section to
 * find out how many were really used), and stop if the norm of the residual
 * is below $10^{-12}$. In practice, the latter criterion will be the one
 * which stops the iteration:
 * 
 * @code
 *   SolverControl solver_control(1000, 1e-12);
 * @endcode
 * 
 * Then we need the solver itself. The template parameter to the SolverCG
 * class is the type of the vectors, and leaving the empty angle brackets
 * would indicate that we are taking the default argument (which is
 * <code>Vector@<double@></code>). However, we explicitly mention the template
 * argument:
 * 
 * @code
 *   SolverCG<Vector<double>> solver(solver_control);
 * 
 * @endcode
 * 
 * Now solve the system of equations. The CG solver takes a preconditioner
 * as its fourth argument. We don't feel ready to delve into this yet, so we
 * tell it to use the identity operation as preconditioner:
 * 
 * @code
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 * @endcode
 * 
 * Now that the solver has done its job, the solution variable contains the
 * nodal values of the solution function.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3output_results"></a> 
 * <h4>Step3::output_results</h4>
 * 

 * 
 * The last part of a typical finite element program is to output the results
 * and maybe do some postprocessing (for example compute the maximal stress
 * values at the boundary, or the average flux across the outflow, etc). We
 * have no such postprocessing here, but we would like to write the solution
 * to a file.
 * 
 * @code
 * void Step3::output_results() const
 * {
 * @endcode
 * 
 * To write the output to a file, we need an object which knows about output
 * formats and the like. This is the DataOut class, and we need an object of
 * that type:
 * 
 * @code
 *   DataOut<2> data_out;
 * @endcode
 * 
 * Now we have to tell it where to take the values from which it shall
 * write. We tell it which DoFHandler object to use, and the solution vector
 * (and the name by which the solution variable shall appear in the output
 * file). If we had more than one vector which we would like to look at in
 * the output (for example right hand sides, errors per cell, etc) we would
 * add them as well:
 * 
 * @code
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 * @endcode
 * 
 * After the DataOut object knows which data it is to work on, we have to
 * tell it to process them into something the back ends can handle. The
 * reason is that we have separated the frontend (which knows about how to
 * treat DoFHandler objects and data vectors) from the back end (which knows
 * many different output formats) and use an intermediate data format to
 * transfer data from the front- to the backend. The data is transformed
 * into this intermediate format by the following function:
 * 
 * @code
 *   data_out.build_patches();
 * 
 * @endcode
 * 
 * Now we have everything in place for the actual output. Just open a file
 * and write the data into it, using VTK format (there are many other
 * functions in the DataOut class we are using here that can write the
 * data in postscript, AVS, GMV, Gnuplot, or some other file
 * formats):
 * 
 * @code
 *   std::ofstream output("solution.vtk");
 *   data_out.write_vtk(output);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3run"></a> 
 * <h4>Step3::run</h4>
 * 

 * 
 * Finally, the last function of this class is the main function which calls
 * all the other functions of the <code>Step3</code> class. The order in which
 * this is done resembles the order in which most finite element programs
 * work. Since the names are mostly self-explanatory, there is not much to
 * comment about:
 * 
 * @code
 * void Step3::run()
 * {
 *   make_grid();
 *   setup_system();
 *   assemble_system();
 *   solve();
 *   output_results();
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
 * This is the main function of the program. Since the concept of a
 * main function is mostly a remnant from the pre-object oriented era
 * before C++ programming, it often does not do much more than
 * creating an object of the top-level class and calling its principle
 * function.
 * 

 * 
 * Finally, the first line of the function is used to enable output of
 * some diagnostics that deal.II can generate.  The @p deallog
 * variable (which stands for deal-log, not de-allog) represents a
 * stream to which some parts of the library write output. For
 * example, iterative solvers will generate diagnostics (starting
 * residual, number of solver steps, final residual) as can be seen
 * when running this tutorial program.
 * 

 * 
 * The output of @p deallog can be written to the console, to a file,
 * or both. Both are disabled by default since over the years we have
 * learned that a program should only generate output when a user
 * explicitly asks for it. But this can be changed, and to explain how
 * this can be done, we need to explain how @p deallog works: When
 * individual parts of the library want to log output, they open a
 * "context" or "section" into which this output will be placed. At
 * the end of the part that wants to write output, one exits this
 * section again. Since a function may call another one from within
 * the scope where this output section is open, output may in fact be
 * nested hierarchically into these sections. The LogStream class of
 * which @p deallog is a variable calls each of these sections a
 * "prefix" because all output is printed with this prefix at the left
 * end of the line, with prefixes separated by colons. There is always
 * a default prefix called "DEAL" (a hint at deal.II's history as the
 * successor of a previous library called "DEAL" and from which the
 * LogStream class is one of the few pieces of code that were taken
 * into deal.II).
 * 

 * 
 * By default, @p logstream only outputs lines with zero prefixes --
 * i.e., all output is disabled because the default "DEAL" prefix is
 * always there. But one can set a different maximal number of
 * prefixes for lines that should be output to something larger, and
 * indeed here we set it to two by calling
 * LogStream::depth_console(). This means that for all screen output,
 * a context that has pushed one additional prefix beyond the default
 * "DEAL" is allowed to print its output to the screen ("console"),
 * whereas all further nested sections that would have three or more
 * prefixes active would write to @p deallog, but @p deallog does not
 * forward this output to the screen. Thus, running this example (or
 * looking at the "Results" section), you will see the solver
 * statistics prefixed with "DEAL:CG", which is two prefixes. This is
 * sufficient for the context of the current program, but you will see
 * examples later on (e.g., in step-22) where solvers are nested more
 * deeply and where you may get useful information by setting the
 * depth even higher.
 * 
 * @code
 * int main()
 * {
 *   deallog.depth_console(2);
 * 
 *   Step3 laplace_problem;
 *   laplace_problem.run();
 * 
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


The output of the program looks as follows:
@code
Number of active cells: 1024
Number of degrees of freedom: 1089
DEAL:cg::Starting value 0.121094
DEAL:cg::Convergence step 48 value 5.33692e-13
@endcode

The first two lines is what we wrote to <code>cout</code>. The last
two lines were generated without our intervention by the CG
solver. The first two lines state the residual at the start of the
iteration, while the last line tells us that the solver needed 47
iterations to bring the norm of the residual to 5.3e-13, i.e. below
the threshold 1e-12 which we have set in the `solve' function. We will
show in the next program how to suppress this output, which is
sometimes useful for debugging purposes, but often clutters up the
screen display.

Apart from the output shown above, the program generated the file
<code>solution.vtk</code>, which is in the VTK format that is widely
used by many visualization programs today -- including the two
heavy-weights <a href="https://www.llnl.gov/visit">VisIt</a> and
<a href="https://www.paraview.org">Paraview</a> that are the most
commonly used programs for this purpose today.

Using VisIt, it is not very difficult to generate a picture of the
solution like this:
<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-3.solution-3.png" alt="Visualization of the solution of step-3">
    </td>
  </tr>
</table>
It shows both the solution and the mesh, elevated above the $x$-$y$ plane
based on the value of the solution at each point. Of course the solution
here is not particularly exciting, but that is a result of both what the
Laplace equation represents and the right hand side $f(\mathbf x)=1$ we
have chosen for this program: The Laplace equation describes (among many
other uses) the vertical deformation of a membrane subject to an external
(also vertical) force. In the current example, the membrane's borders
are clamped to a square frame with no vertical variation; a constant
force density will therefore intuitively lead to a membrane that
simply bulges upward -- like the one shown above.

VisIt and Paraview both allow playing with various kinds of visualizations
of the solution. Several video lectures show how to use these programs.
@dealiiVideoLectureSeeAlso{11,32}



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


If you want to play around a little bit with this program, here are a few
suggestions:
</p>

<ul>
  <li>
  Change the geometry and mesh: In the program, we have generated a square
  domain and mesh by using the <code>GridGenerator::hyper_cube</code>
  function. However, the <code>GridGenerator</code> has a good number of other
  functions as well. Try an L-shaped domain, a ring, or other domains you find
  there.
  </li>

  <li>
  Change the boundary condition: The code uses the Functions::ZeroFunction
  function to generate zero boundary conditions. However, you may want to try
  non-zero constant boundary values using
  <code>ConstantFunction&lt;2&gt;(1)</code> instead of
  <code>ZeroFunction&lt;2&gt;()</code> to have unit Dirichlet boundary
  values. More exotic functions are described in the documentation of the
  Functions namespace, and you may pick one to describe your particular boundary
  values.
  </li>

  <li> Modify the type of boundary condition: Presently, what happens
  is that we use Dirichlet boundary values all around, since the
  default is that all boundary parts have boundary indicator zero, and
  then we tell the
  VectorTools::interpolate_boundary_values() function to
  interpolate boundary values to zero on all boundary components with
  indicator zero.  <p> We can change this behavior if we assign parts
  of the boundary different indicators. For example, try this
  immediately after calling GridGenerator::hyper_cube():
  @code
  triangulation.begin_active()->face(0)->set_boundary_id(1);
  @endcode

  What this does is it first asks the triangulation to
  return an iterator that points to the first active cell. Of course,
  this being the coarse mesh for the triangulation of a square, the
  triangulation has only a single cell at this moment, and it is
  active. Next, we ask the cell to return an iterator to its first
  face, and then we ask the face to reset the boundary indicator of
  that face to 1. What then follows is this: When the mesh is refined,
  faces of child cells inherit the boundary indicator of their
  parents, i.e. even on the finest mesh, the faces on one side of the
  square have boundary indicator 1. Later, when we get to
  interpolating boundary conditions, the
  VectorTools::interpolate_boundary_values() call will only produce boundary
  values for those faces that have zero boundary indicator, and leave
  those faces alone that have a different boundary indicator. What
  this then does is to impose Dirichlet boundary conditions on the
  former, and homogeneous Neumann conditions on the latter (i.e. zero
  normal derivative of the solution, unless one adds additional terms
  to the right hand side of the variational equality that deal with
  potentially non-zero Neumann conditions). You will see this if you
  run the program.

  An alternative way to change the boundary indicator is to label
  the boundaries based on the Cartesian coordinates of the face centers.
  For example, we can label all of the cells along the top and
  bottom boundaries with a boundary indicator 1 by checking to
  see if the cell centers' y-coordinates are within a tolerance
  (here 1e-12) of -1 and 1. Try this immediately after calling
  GridGenerator::hyper_cube(), as before:
  @code
  for (auto &face : triangulation.active_face_iterators())
    if (std::fabs(face->center()(1) - (-1.0)) < 1e-12 ||
        std::fabs(face->center()(1) - (1.0)) < 1e-12)
      face->set_boundary_id(1);
  @endcode
  Although this code is a bit longer than before, it is useful for
  complex geometries, as it does not require knowledge of face labels.

  <li>
  A slight variation of the last point would be to set different boundary
  values as above, but then use a different boundary value function for
  boundary indicator one. In practice, what you have to do is to add a second
  call to <code>interpolate_boundary_values</code> for boundary indicator one:
  @code
  VectorTools::interpolate_boundary_values(dof_handler,
					   1,
					   ConstantFunction<2>(1.),
					   boundary_values);
  @endcode
  If you have this call immediately after the first one to this function, then
  it will interpolate boundary values on faces with boundary indicator 1 to the
  unit value, and merge these interpolated values with those previously
  computed for boundary indicator 0. The result will be that we will get
  discontinuous boundary values, zero on three sides of the square, and one on
  the fourth.

  <li>
  Observe convergence: We will only discuss computing errors in norms in
  step-7, but it is easy to check that computations converge
  already here. For example, we could evaluate the value of the solution in a
  single point and compare the value for different %numbers of global
  refinement (the number of global refinement steps is set in
  <code>LaplaceProblem::make_grid</code> above). To evaluate the
  solution at a point, say at $(\frac 13, \frac 13)$, we could add the
  following code to the <code>LaplaceProblem::output_results</code> function:
  @code
    std::cout << "Solution at (1/3,1/3): "
              << VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3))
              << std::endl;
  @endcode
  For 1 through 9 global refinement steps, we then get the following sequence
  of point values:
  <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$u_h(\frac 13,\frac13)$</th> </tr>
    <tr> <td>1</td> <td>0.166667</td> </tr>
    <tr> <td>2</td> <td>0.227381</td> </tr>
    <tr> <td>3</td> <td>0.237375</td> </tr>
    <tr> <td>4</td> <td>0.240435</td> </tr>
    <tr> <td>5</td> <td>0.241140</td> </tr>
    <tr> <td>6</td> <td>0.241324</td> </tr>
    <tr> <td>7</td> <td>0.241369</td> </tr>
    <tr> <td>8</td> <td>0.241380</td> </tr>
    <tr> <td>9</td> <td>0.241383</td> </tr>
  </table>
  By noticing that the difference between each two consecutive values reduces
  by about a factor of 4, we can conjecture that the "correct" value may be
  $u(\frac 13, \frac 13)\approx 0.241384$. In fact, if we assumed this to be
  the correct value, we could show that the sequence above indeed shows ${\cal
  O}(h^2)$ convergence &mdash; theoretically, the convergence order should be
  ${\cal O}(h^2 |\log h|)$ but the symmetry of the domain and the mesh may lead
  to the better convergence order observed.

  A slight variant of this would be to repeat the test with quadratic
  elements. All you need to do is to set the polynomial degree of the finite
  element to two in the constructor
  <code>LaplaceProblem::LaplaceProblem</code>.

  <li>Convergence of the mean: A different way to see that the solution
  actually converges (to something &mdash; we can't tell whether it's really
  the correct value!) is to compute the mean of the solution. To this end, add
  the following code to <code>LaplaceProblem::output_results</code>:
  @code
    std::cout << "Mean value: "
              << VectorTools::compute_mean_value (dof_handler,
						  QGauss<2>(fe.degree + 1),
						  solution,
						  0)
              << std::endl;
  @endcode
  The documentation of the function explains what the second and fourth
  parameters mean, while the first and third should be obvious. Doing the same
  study again where we change the number of global refinement steps, we get
  the following result:
  <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$\int_\Omega u_h(x)\; dx$</th> </tr>
    <tr> <td>0</td> <td>0.09375000</td> </tr>
    <tr> <td>1</td> <td>0.12790179</td> </tr>
    <tr> <td>2</td> <td>0.13733440</td> </tr>
    <tr> <td>3</td> <td>0.13976069</td> </tr>
    <tr> <td>4</td> <td>0.14037251</td> </tr>
    <tr> <td>5</td> <td>0.14052586</td> </tr>
    <tr> <td>6</td> <td>0.14056422</td> </tr>
    <tr> <td>7</td> <td>0.14057382</td> </tr>
    <tr> <td>8</td> <td>0.14057622</td> </tr>
  </table>
  Again, the difference between two adjacent values goes down by about a
  factor of four, indicating convergence as ${\cal O}(h^2)$.
</ul>



<a name="UsingHDF5tooutputthesolutionandadditionaldata"></a><h3>Using %HDF5 to output the solution and additional data</h3>


%HDF5 is a commonly used format that can be read by many scripting
languages (e.g. R or Python). It is not difficult to get deal.II to
produce some %HDF5 files that can then be used in external scripts to
postprocess some of the data generated by this program. Here are some
ideas on what is possible.


<a name="Changingtheoutputtoh5"></a><h4> Changing the output to .h5</h4>


To fully make use of the automation we first need to introduce a private variable for the number of
global refinement steps <code>unsigned int n_refinement_steps </code>, which will be used for the output filename.
In <code>make_grid()</code> we then replace <code>triangulation.refine_global(5);</code> with
@code
n_refinement_steps = 5;
triangulation.refine_global(n_refinement_steps);
@endcode
The deal.II library has two different %HDF5 bindings, one in the HDF5
namespace (for interfacing to general-purpose data files)
and another one in DataOut (specifically for writing files for the
visualization of solutions).
Although the HDF5 deal.II binding supports both serial and MPI, the %HDF5 DataOut binding
only supports parallel output.
For this reason we need to initialize an MPI
communicator with only one processor. This is done by adding the following code.
@code
int main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  ...
}
@endcode
Next we change the `Step3::output_results()` output routine as
described in the DataOutBase namespace documentation:
@code
const std::string filename_h5 = "solution_" + std::to_string(n_refinement_steps) + ".h5";
DataOutBase::DataOutFilterFlags flags(true, true);
DataOutBase::DataOutFilter data_filter(flags);
data_out.write_filtered_data(data_filter);
data_out.write_hdf5_parallel(data_filter, filename_h5, MPI_COMM_WORLD);
@endcode
The resulting file can then be visualized just like the VTK file that
the original version of the tutorial produces; but, since %HDF5 is a
more general file format, it can also easily be processed in scripting
languages for other purposes.


<a name="Addingthepointvalueandthemeanseeextensionaboveintotheh5file"></a><h4> Adding the point value and the mean (see extension above) into the .h5 file</h4>


After outputting the solution, the file can be opened again to include
more datasets.  This allows us to keep all the necessary information
of our experiment in a single result file, which can then be read and
processed by some postprocessing script.
(Have a look at HDF5::Group::write_dataset() for further
information on the possible output options.)

To make this happen, we first include the necessary header into our file:
@code
#include <deal.II/base/hdf5.h>
@endcode
Adding the following lines to the end
of our output routine adds the information about the value of the
solution at a particular point, as well as the mean value of the
solution, to our %HDF5 file:
@code
HDF5::File data_file(filename_h5, HDF5::File::FileAccessMode::open, MPI_COMM_WORLD);
Vector<double> point_value(1);
point_value[0] = VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3));
data_file.write_dataset("point_value", point_value);
Vector<double> mean_value(1);
mean_value[0] = VectorTools::compute_mean_value(dof_handler,
                                                QGauss<2>(fe.degree + 1),
                                                solution, 0);
data_file.write_dataset("mean_value",mean_value);
@endcode



<a name="UsingRandggplot2togenerateplots"></a><h3> Using R and ggplot2 to generate plots</h3>


The data put into %HDF5 files above can then be used from scripting
languages for further postprocessing. In the following, let us show
how this can, in particular, be done with the
<a href="https://en.wikipedia.org/wiki/R_(programming_language)">R
programming language</a>, a widely used language in statistical data
analysis. (Similar things can also be done in Python, for example.)
If you are unfamiliar with R and ggplot2 you could check out the data carpentry course on R
<a href="https://datacarpentry.org/R-ecology-lesson/index.html">here</a>.
Furthermore, since most search engines struggle with searches of the form "R + topic",
we recommend using the specializes service <a
href="http://rseek.org">RSeek </a> instead.

The most prominent difference between R and other languages is that
the assignment operator (`a = 5`) is typically written as
`a <- 5`. As the latter is considered standard we will use it in our examples as well.
To open the `.h5` file in R you have to install the <a href="https://bioconductor.org/packages/release/bioc/html/rhdf5.html">rhdf5</a> package, which is a part of the Bioconductor package.

First we will include all necessary packages and have a look at how the data is structured in our file.
@code{.r}
library(rhdf5)     # library for handling HDF5 files
library(ggplot2)   # main plotting library
library(grDevices) # needed for output to PDF
library(viridis)   # contains good colormaps for sequential data

refinement <- 5
h5f <- H5Fopen(paste("solution_",refinement,".h5",sep=""))
print(h5f)
@endcode
This gives the following output
@code{.unparsed}
HDF5 FILE
   name /
filename

    name       otype  dclass     dim
0 cells       H5I_DATASET INTEGER  x 1024
1 mean_value  H5I_DATASET FLOAT   1
2 nodes       H5I_DATASET FLOAT    x 1089
3 point_value H5I_DATASET FLOAT   1
4 solution    H5I_DATASET FLOAT    x 1089
@endcode
The datasets can be accessed by <code>h5f\$name</code>. The function
<code>dim(h5f\$cells)</code> gives us the dimensions of the matrix
that is used to store our cells.
We can see the following three matrices, as well as the two
additional data points we added.
<ul>
<li> <code>cells</code>: a 4x1024 matrix that stores the  (C++) vertex indices for each cell
<li> <code>nodes</code>: a 2x1089 matrix storing the position values (x,y) for our cell vertices
<li> <code>solution</code>: a 1x1089 matrix storing the values of our solution at each vertex
</ul>
Now we can use this data to generate various plots. Plotting with ggplot2 usually splits into two steps.
At first the data needs to be manipulated and added to a <code>data.frame</code>.
After that, a <code>ggplot</code> object is constructed and manipulated by adding plot elements to it.

<code>nodes</code> and <code>cells</code> contain all the information we need to plot our grid.
The following code wraps all the data into one dataframe for plotting our grid:
@code{.r}
# Counting in R starts at 1 instead of 0, so we need to increment all
# vertex indices by one:
cell_ids <- h5f$cells+1

# Store the x and y positions of each vertex in one big vector in a
# cell by cell fashion (every 4 entries belong to one cell):
cells_x <- h5f$nodes[1,][cell_ids]
cells_y <- h5f$nodes[2,][cell_ids]

# Construct a vector that stores the matching cell by cell grouping
# (1,1,1,1,2,2,2,2,...):
groups <- rep(1:ncol(cell_ids),each=4)

# Finally put everything into one dataframe:
meshdata <- data.frame(x = cells_x, y = cells_y, id = groups)
@endcode

With the finished dataframe we have everything we need to plot our grid:
@code{.r}
pdf (paste("grid_",refinement,".pdf",sep=""),width = 5,height = 5) # Open new PDF file
plt <- ggplot(meshdata,aes(x=x,y=y,group=id))                      # Construction of our plot
                                                                   # object, at first only data

plt <- plt + geom_polygon(fill="white",colour="black")             # Actual plotting of the grid as polygons
plt <- plt + ggtitle(paste("grid at refinement level #",refinement))

print(plt)                                                         # Show the current state of the plot/add it to the pdf
dev.off()                                                          # Close PDF file
@endcode

The contents of this file then look as follows (not very exciting, but
you get the idea):
<table width="60%" align="center">
  <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.grid_5.png" alt="Grid after 5 refinement steps of step-3">
   </td>
  </tr>
</table>

We can also visualize the solution itself, and this is going to look
more interesting.
To make a 2D pseudocolor plot of our solution we will use <code>geom_raster</code>.
This function needs a structured grid, i.e. uniform in x and y directions.
Luckily our data at this point is structured in the right way.
The following code plots a pseudocolor representation of our surface into a new PDF:
@code{.r}
pdf (paste("pseudocolor_",refinement,".pdf",sep=""),width = 5,height = 4.2) # Open new PDF file
colordata <- data.frame(x = h5f$nodes[1,],y = h5f$nodes[2,] , solution = h5f$solution[1,])
plt <- ggplot(colordata,aes(x=x,y=y,fill=solution))
plt <- plt + geom_raster(interpolate=TRUE)
plt <- plt + scale_fill_viridis()
plt <- plt + ggtitle(paste("solution at refinement level #",refinement))

print(plt)
dev.off()
H5Fclose(h5f) # Close the HDF5 file
@endcode
This is now going to look as follows:
<table width="60%" align="center">
 <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.pseudocolor_5.png" alt="Solution after 5 refinement steps of step-3">
   </td>
 </tr>
</table>

For plotting the converge curves we need to re-run the C++ code multiple times with different values for <code>n_refinement_steps</code>
starting from 1.
Since every file only contains a single data point we need to loop over them and concatenate the results into a single vector.
@code{.r}
n_ref <- 8   # Maximum refinement level for which results are existing

# First we initiate all vectors with the results of the first level
h5f   <- H5Fopen("solution_1.h5")
dofs  <- dim(h5f$solution)[2]
mean  <- h5f$mean_value
point <- h5f$point_value
H5Fclose(h5f)

for (reflevel in 2:n_ref)
{
   h5f   <- H5Fopen(paste("solution_",reflevel,".h5",sep=""))
   dofs  <- c(dofs,dim(h5f\$solution)[2])
   mean  <- c(mean,h5f\$mean_value)
   point <- c(point,h5f\$point_value)
   H5Fclose(h5f)
}
@endcode
As we are not interested in the values themselves but rather in the error compared to a "exact" solution we will
assume our highest refinement level to be that solution and omit it from the data.
@code{.r}
# Calculate the error w.r.t. our maximum refinement step
mean_error  <- abs(mean[1:n_ref-1]-mean[n_ref])
point_error <- abs(point[1:n_ref-1]-point[n_ref])

# Remove the highest value from our DoF data
dofs     <- dofs[1:n_ref-1]
convdata <- data.frame(dofs = dofs, mean_value= mean_error, point_value = point_error)
@endcode
Now we have all the data available to generate our plots.
It is often useful to plot errors on a log-log scale, which is
accomplished in the following code:
@code
pdf (paste("convergence.pdf",sep=""),width = 5,height = 4.2)
plt <- ggplot(convdata,mapping=aes(x = dofs, y = mean_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "mean value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)

plt <- ggplot(convdata,mapping=aes(x = dofs, y = point_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "point value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)

dev.off()
@endcode
This results in the following plot that shows how the errors in the
mean value and the solution value at the chosen point nicely converge
to zero:
<table style="width:50%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_mean.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_point.png" alt=""></td>
  </tr>
</table>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-3.cc"
*/
