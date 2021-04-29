/**
@page step_4 The step-4 tutorial program
This tutorial depends on step-3.

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
        <li><a href="#ThecodeStep4codeclasstemplate">The <code>Step4</code> class template</a>
        <li><a href="#Righthandsideandboundaryvalues">Right hand side and boundary values</a>
        <li><a href="#ImplementationofthecodeStep4codeclass">Implementation of the <code>Step4</code> class</a>
      <ul>
        <li><a href="#Step4Step4">Step4::Step4</a>
        <li><a href="#Step4make_grid">Step4::make_grid</a>
        <li><a href="#Step4setup_system">Step4::setup_system</a>
        <li><a href="#Step4assemble_system">Step4::assemble_system</a>
        <li><a href="#Step4solve">Step4::solve</a>
        <li><a href="#Step4output_results">Step4::output_results</a>
        <li><a href="#Step4run">Step4::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


@dealiiVideoLecture{12,13}

deal.II has a unique feature which we call
``dimension independent programming''. You may have noticed in the
previous examples that many classes had a number in angle brackets
suffixed to them. This is to indicate that for example the
triangulation in two and three space dimensions are different, but
related data %types. We could as well have called them
<code>Triangulation2d</code> and <code>Triangulation3d</code> instead
of <code>Triangulation@<2@></code> and
<code>Triangulation@<3@></code> to name the two classes, but this
has an important drawback: assume you have a function which does
exactly the same functionality, but on 2d or 3d triangulations,
depending on which dimension we would like to solve the equation in
presently (if you don't believe that it is the common case that a
function does something that is the same in all dimensions, just take
a look at the code below - there are almost no distinctions between 2d
and 3d!). We would have to write the same function twice, once
working on <code>Triangulation2d</code> and once working with a
<code>Triangulation3d</code>. This is an unnecessary obstacle in
programming and leads to a nuisance to keep the two function in sync
(at best) or difficult to find errors if the two versions get out of
sync (at worst; this would probably the more common case).




Such obstacles can be circumvented by using some template magic as
provided by the C++ language: templatized classes and functions are
not really classes or functions but only a pattern depending on an
as-yet undefined data type parameter or on a numerical value which is
also unknown at the point of definition. However, the compiler can
build proper classes or functions from these templates if you provide
it with the information that is needed for that. Of course, parts of
the template can depend on the template parameters, and they will be
resolved at the time of compilation for a specific template
parameter. For example, consider the following piece of code:
@code
  template <int dim>
  void make_grid (Triangulation<dim> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode



At the point where the compiler sees this function, it does not know
anything about the actual value of <code>dim</code>. The only thing the compiler has is
a template, i.e. a blueprint, to generate
functions <code>make_grid</code> if given a particular value of
<code>dim</code>. Since <code>dim</code> has an unknown value, there is no
code the compiler can generate for the moment.



However, if later down the compiler would encounter code that looks, for
example, like this,
@code
  Triangulation<2> triangulation;
  make_grid (triangulation);
@endcode
then the compiler will deduce that the function <code>make_grid</code> for
<code>dim==2</code> was
requested and will compile the template above into a function with dim replaced
by 2 everywhere, i.e. it will compile the function as if it were defined
as
@code
  void make_grid (Triangulation<2> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode



However, it is worth to note that the function
<code>GridGenerator::hyper_cube</code> depends on the dimension as
well, so in this case, the compiler will call the function
<code>GridGenerator::hyper_cube@<2@></code> while if dim were 3,
it would call <code>GridGenerator::hyper_cube@<3@></code> which
might be (and actually is) a totally unrelated  function.



The same can be done with member variables. Consider the following
function, which might in turn call the above one:
@code
  template <int dim>
  void make_grid_and_dofs (Triangulation<dim> &triangulation)
  {
    make_grid (triangulation);

    DoFHandler<dim> dof_handler(triangulation);
    ...
  };
@endcode
This function has a member variable of type
<code>DoFHandler@<dim@></code>. Again, the compiler can't
compile this function until it knows for which dimension. If you call
this function for a specific dimension as above, the compiler will
take the template, replace all occurrences of dim by the dimension for
which it was called, and compile it. If you call the function several
times for different dimensions, it will compile it several times, each
time calling the right <code>make_grid</code> function and reserving the right
amount of memory for the member variable; note that the size of a
<code>DoFHandler</code> might, and indeed does, depend on the space dimension.



The deal.II library is built around this concept
of dimension-independent programming, and therefore allows you to program in
a way that will not need to
distinguish between the space dimensions. It should be noted that in
only a very few places is it necessary to actually compare the
dimension using <code>if</code>s or <code>switch</code>es. However, since the compiler
has to compile each function for each dimension separately, even there
it knows the value of <code>dim</code> at the time of compilation and will
therefore be able to optimize away the <code>if</code> statement along with the
unused branch.



In this example program, we will show how to program dimension
independently (which in fact is even simpler than if you had to take
care about the dimension) and we will extend the Laplace problem of
the last example to a program that runs in two and three space
dimensions at the same time. Other extensions are the use of a
non-constant right hand side function and of non-zero boundary values.


@note When using templates, C++ imposes all sorts of syntax constraints that
make it sometimes a bit difficult to understand why exactly something has to
be written this way. A typical example is the need to use the keyword
<code>typename</code> in so many places. If you are not entirely familiar with
this already, then several of these difficulties are explained in the deal.II
Frequently Asked Questions (FAQ) linked to from the <a
href="http://www.dealii.org/">deal.II homepage</a>.

<!--We need a blank line to end the above block properly.-->
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
 * The first few (many?) include files have already been used in the previous
 * example, so we will not explain their meaning here again.
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * This is new, however: in the previous example we got some unwanted output
 * from the linear solvers. If we want to suppress it, we have to include this
 * file and add a single line somewhere to the program (see the main()
 * function below for that):
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * 
 * @endcode
 * 
 * The final step, as in previous programs, is to import all the deal.II class
 * and function names into the global namespace:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep4codeclasstemplate"></a> 
 * <h3>The <code>Step4</code> class template</h3>
 * 

 * 
 * This is again the same <code>Step4</code> class as in the previous
 * example. The only difference is that we have now declared it as a class
 * with a template parameter, and the template parameter is of course the
 * spatial dimension in which we would like to solve the Laplace equation. Of
 * course, several of the member variables depend on this dimension as well,
 * in particular the Triangulation class, which has to represent
 * quadrilaterals or hexahedra, respectively. Apart from this, everything is
 * as before.
 * 
 * @code
 * template <int dim>
 * class Step4
 * {
 * public:
 *   Step4();
 *   void run();
 * 
 * private:
 *   void make_grid();
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void output_results() const;
 * 
 *   Triangulation<dim> triangulation;
 *   FE_Q<dim>          fe;
 *   DoFHandler<dim>    dof_handler;
 * 
 *   SparsityPattern      sparsity_pattern;
 *   SparseMatrix<double> system_matrix;
 * 
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideandboundaryvalues"></a> 
 * <h3>Right hand side and boundary values</h3>
 * 

 * 
 * In the following, we declare two more classes denoting the right hand side
 * and the non-homogeneous Dirichlet boundary values. Both are functions of a
 * dim-dimensional space variable, so we declare them as templates as well.
 * 

 * 
 * Each of these classes is derived from a common, abstract base class
 * Function, which declares the common interface which all functions have to
 * follow. In particular, concrete classes have to overload the
 * <code>value</code> function, which takes a point in dim-dimensional space
 * as parameters and returns the value at that point as a
 * <code>double</code> variable.
 * 

 * 
 * The <code>value</code> function takes a second argument, which we have here
 * named <code>component</code>: This is only meant for vector-valued
 * functions, where you may want to access a certain component of the vector
 * at the point <code>p</code>. However, our functions are scalar, so we need
 * not worry about this parameter and we will not use it in the implementation
 * of the functions. Inside the library's header files, the Function base
 * class's declaration of the <code>value</code> function has a default value
 * of zero for the component, so we will access the <code>value</code>
 * function of the right hand side with only one parameter, namely the point
 * where we want to evaluate the function. A value for the component can then
 * simply be omitted for scalar functions.
 * 

 * 
 * Function objects are used in lots of places in the library (for example, in
 * step-3 we used a Functions::ZeroFunction instance as an argument to
 * VectorTools::interpolate_boundary_values) and this is the first tutorial
 * where we define a new class that inherits from Function. Since we only ever
 * call Function::value(), we could get away with just a plain function (and
 * this is what is done in step-5), but since this is a tutorial we inherit from
 * Function for the sake of example.
 * 
 * @code
 * template <int dim>
 * class RightHandSide : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 * };
 * 
 * 
 * 
 * template <int dim>
 * class BoundaryValues : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 * };
 * 
 * @endcode
 * 
 * If you are not familiar with what the keywords `virtual` and `override` in
 * the function declarations above mean, you will probably want to take a look
 * at your favorite C++ book or an online tutorial such as
 * http://www.cplusplus.com/doc/tutorial/polymorphism/ . In essence, what is
 * happening here is that Function<dim> is an "abstract" base class that
 * declares a certain "interface" -- a set of functions one can call on
 * objects of this kind. But it does not actually *implement* these functions:
 * it just says "this is how Function objects look like", but what kind of
 * function it actually is, is left to derived classes that implement
 * the `value()` function.
 * 

 * 
 * Deriving one class from another is often called an "is-a" relationship
 * function. Here, the `RightHandSide` class "is a" Function class
 * because it implements the interface described by the Function base class.
 * (The actual implementation of the `value()` function is in the code block
 * below.) The `virtual` keyword then means "Yes, the
 * function here is one that can be overridden by derived classes",
 * and the `override` keyword means "Yes, this is in fact a function we know
 * has been declared as part of the base class". The `override` keyword is not
 * strictly necessary, but is an insurance against typos: If we get the name
 * of the function or the type of one argument wrong, the compiler will warn
 * us by stating "You say that this function overrides one in a base class,
 * but I don't actually know any such function with this name and these
 * arguments."
 * 

 * 
 * But back to the concrete case here:
 * For this tutorial, we choose as right hand side the function
 * $4(x^4+y^4)$ in 2D, or $4(x^4+y^4+z^4)$ in 3D. We could write this
 * distinction using an if-statement on the space dimension, but here is a
 * simple way that also allows us to use the same function in 1D (or in 4D, if
 * you should desire to do so), by using a short loop.  Fortunately, the
 * compiler knows the size of the loop at compile time (remember that at the
 * time when you define the template, the compiler doesn't know the value of
 * <code>dim</code>, but when it later encounters a statement or declaration
 * <code>RightHandSide@<2@></code>, it will take the template, replace all
 * occurrences of dim by 2 and compile the resulting function).  In other
 * words, at the time of compiling this function, the number of times the body
 * will be executed is known, and the compiler can minimize the overhead
 * needed for the loop; the result will be as fast as if we had used the
 * formulas above right away.
 * 

 * 
 * The last thing to note is that a <code>Point@<dim@></code> denotes a point
 * in dim-dimensional space, and its individual components (i.e. $x$, $y$,
 * ... coordinates) can be accessed using the () operator (in fact, the []
 * operator will work just as well) with indices starting at zero as usual in
 * C and C++.
 * 
 * @code
 * template <int dim>
 * double RightHandSide<dim>::value(const Point<dim> &p,
 *                                  const unsigned int /*component*/) const
 * {
 *   double return_value = 0.0;
 *   for (unsigned int i = 0; i < dim; ++i)
 *     return_value += 4.0 * std::pow(p(i), 4.0);
 * 
 *   return return_value;
 * }
 * 
 * 
 * @endcode
 * 
 * As boundary values, we choose $x^2+y^2$ in 2D, and $x^2+y^2+z^2$ in 3D. This
 * happens to be equal to the square of the vector from the origin to the
 * point at which we would like to evaluate the function, irrespective of the
 * dimension. So that is what we return:
 * 
 * @code
 * template <int dim>
 * double BoundaryValues<dim>::value(const Point<dim> &p,
 *                                   const unsigned int /*component*/) const
 * {
 *   return p.square();
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeStep4codeclass"></a> 
 * <h3>Implementation of the <code>Step4</code> class</h3>
 * 

 * 
 * Next for the implementation of the class template that makes use of the
 * functions above. As before, we will write everything as templates that have
 * a formal parameter <code>dim</code> that we assume unknown at the time we
 * define the template functions. Only later, the compiler will find a
 * declaration of <code>Step4@<2@></code> (in the <code>main</code> function,
 * actually) and compile the entire class with <code>dim</code> replaced by 2,
 * a process referred to as `instantiation of a template'. When doing so, it
 * will also replace instances of <code>RightHandSide@<dim@></code> by
 * <code>RightHandSide@<2@></code> and instantiate the latter class from the
 * class template.
 * 

 * 
 * In fact, the compiler will also find a declaration <code>Step4@<3@></code>
 * in <code>main()</code>. This will cause it to again go back to the general
 * <code>Step4@<dim@></code> template, replace all occurrences of
 * <code>dim</code>, this time by 3, and compile the class a second time. Note
 * that the two instantiations <code>Step4@<2@></code> and
 * <code>Step4@<3@></code> are completely independent classes; their only
 * common feature is that they are both instantiated from the same general
 * template, but they are not convertible into each other, for example, and
 * share no code (both instantiations are compiled completely independently).
 * 

 * 
 * 

 * 
 * 
 * <a name="Step4Step4"></a> 
 * <h4>Step4::Step4</h4>
 * 

 * 
 * After this introduction, here is the constructor of the <code>Step4</code>
 * class. It specifies the desired polynomial degree of the finite elements
 * and associates the DoFHandler to the triangulation just as in the previous
 * example program, step-3:
 * 
 * @code
 * template <int dim>
 * Step4<dim>::Step4()
 *   : fe(1)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4make_grid"></a> 
 * <h4>Step4::make_grid</h4>
 * 

 * 
 * Grid creation is something inherently dimension dependent. However, as long
 * as the domains are sufficiently similar in 2D or 3D, the library can
 * abstract for you. In our case, we would like to again solve on the square
 * $[-1,1]\times [-1,1]$ in 2D, or on the cube $[-1,1] \times [-1,1] \times
 * [-1,1]$ in 3D; both can be termed GridGenerator::hyper_cube(), so we may
 * use the same function in whatever dimension we are. Of course, the
 * functions that create a hypercube in two and three dimensions are very much
 * different, but that is something you need not care about. Let the library
 * handle the difficult things.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::make_grid()
 * {
 *   GridGenerator::hyper_cube(triangulation, -1, 1);
 *   triangulation.refine_global(4);
 * 
 *   std::cout << "   Number of active cells: " << triangulation.n_active_cells()
 *             << std::endl
 *             << "   Total number of cells: " << triangulation.n_cells()
 *             << std::endl;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Step4setup_system"></a> 
 * <h4>Step4::setup_system</h4>
 * 

 * 
 * This function looks exactly like in the previous example, although it
 * performs actions that in their details are quite different if
 * <code>dim</code> happens to be 3. The only significant difference from a
 * user's perspective is the number of cells resulting, which is much higher
 * in three than in two space dimensions!
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 * 
 *   std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 * 
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *   sparsity_pattern.copy_from(dsp);
 * 
 *   system_matrix.reinit(sparsity_pattern);
 * 
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4assemble_system"></a> 
 * <h4>Step4::assemble_system</h4>
 * 

 * 
 * Unlike in the previous example, we would now like to use a non-constant
 * right hand side function and non-zero boundary values. Both are tasks that
 * are readily achieved with only a few new lines of code in the assemblage of
 * the matrix and right hand side.
 * 

 * 
 * More interesting, though, is the way we assemble matrix and right hand side
 * vector dimension independently: there is simply no difference to the
 * two-dimensional case. Since the important objects used in this function
 * (quadrature formula, FEValues) depend on the dimension by way of a template
 * parameter as well, they can take care of setting up properly everything for
 * the dimension for which this function is compiled. By declaring all classes
 * which might depend on the dimension using a template parameter, the library
 * can make nearly all work for you and you don't have to care about most
 * things.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::assemble_system()
 * {
 *   QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 * @endcode
 * 
 * We wanted to have a non-constant right hand side, so we use an object of
 * the class declared above to generate the necessary data. Since this right
 * hand side object is only used locally in the present function, we declare
 * it here as a local variable:
 * 
 * @code
 *   RightHandSide<dim> right_hand_side;
 * 
 * @endcode
 * 
 * Compared to the previous example, in order to evaluate the non-constant
 * right hand side function we now also need the quadrature points on the
 * cell we are presently on (previously, we only required values and
 * gradients of the shape function from the FEValues object, as well as the
 * quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
 * do for us by also giving it the #update_quadrature_points flag:
 * 
 * @code
 *   FEValues<dim> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients |
 *                             update_quadrature_points | update_JxW_values);
 * 
 * @endcode
 * 
 * We then again define the same abbreviation as in the previous program.
 * The value of this variable of course depends on the dimension which we
 * are presently using, but the FiniteElement class does all the necessary
 * work for you and you don't have to care about the dimension dependent
 * parts:
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Next, we again have to loop over all cells and assemble local
 * contributions.  Note, that a cell is a quadrilateral in two space
 * dimensions, but a hexahedron in 3D. In fact, the
 * <code>active_cell_iterator</code> data type is something different,
 * depending on the dimension we are in, but to the outside world they look
 * alike and you will probably never see a difference. In any case, the real
 * type is hidden by using `auto`:
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       fe_values.reinit(cell);
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 * @endcode
 * 
 * Now we have to assemble the local matrix and right hand side. This is
 * done exactly like in the previous example, but now we revert the
 * order of the loops (which we can safely do since they are independent
 * of each other) and merge the loops for the local matrix and the local
 * vector as far as possible to make things a bit faster.
 *       

 * 
 * Assembling the right hand side presents the only significant
 * difference to how we did things in step-3: Instead of using a
 * constant right hand side with value 1, we use the object representing
 * the right hand side and evaluate it at the quadrature points:
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             for (const unsigned int j : fe_values.dof_indices())
 *               cell_matrix(i, j) +=
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                  fe_values.JxW(q_index));           // dx
 * 
 *             const auto x_q = fe_values.quadrature_point(q_index);
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                             right_hand_side.value(x_q) *        // f(x_q)
 *                             fe_values.JxW(q_index));            // dx
 *           }
 * @endcode
 * 
 * As a final remark to these loops: when we assemble the local
 * contributions into <code>cell_matrix(i,j)</code>, we have to multiply
 * the gradients of shape functions $i$ and $j$ at point number
 * q_index and
 * multiply it with the scalar weights JxW. This is what actually
 * happens: <code>fe_values.shape_grad(i,q_index)</code> returns a
 * <code>dim</code> dimensional vector, represented by a
 * <code>Tensor@<1,dim@></code> object, and the operator* that
 * multiplies it with the result of
 * <code>fe_values.shape_grad(j,q_index)</code> makes sure that the
 * <code>dim</code> components of the two vectors are properly
 * contracted, and the result is a scalar floating point number that
 * then is multiplied with the weights. Internally, this operator* makes
 * sure that this happens correctly for all <code>dim</code> components
 * of the vectors, whether <code>dim</code> be 2, 3, or any other space
 * dimension; from a user's perspective, this is not something worth
 * bothering with, however, making things a lot simpler if one wants to
 * write code dimension independently.
 * 

 * 
 * With the local systems assembled, the transfer into the global matrix
 * and right hand side is done exactly as before, but here we have again
 * merged some loops for efficiency:
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices);
 *       for (const unsigned int i : fe_values.dof_indices())
 *         {
 *           for (const unsigned int j : fe_values.dof_indices())
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               cell_matrix(i, j));
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     }
 * 
 * @endcode
 * 
 * As the final step in this function, we wanted to have non-homogeneous
 * boundary values in this example, unlike the one before. This is a simple
 * task, we only have to replace the Functions::ZeroFunction used there by an
 * object of the class which describes the boundary values we would like to
 * use (i.e. the <code>BoundaryValues</code> class declared above):
 *   

 * 
 * The function VectorTools::interpolate_boundary_values() will only work
 * on faces that have been marked with boundary indicator 0 (because that's
 * what we say the function should work on with the second argument below).
 * If there are faces with boundary id other than 0, then the function
 * interpolate_boundary_values will do nothing on these faces. For
 * the Laplace equation doing nothing is equivalent to assuming that
 * on those parts of the boundary a zero Neumann boundary condition holds.
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values;
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            BoundaryValues<dim>(),
 *                                            boundary_values);
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
 * <a name="Step4solve"></a> 
 * <h4>Step4::solve</h4>
 * 

 * 
 * Solving the linear system of equations is something that looks almost
 * identical in most programs. In particular, it is dimension independent, so
 * this function is copied verbatim from the previous example.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::solve()
 * {
 *   SolverControl            solver_control(1000, 1e-12);
 *   SolverCG<Vector<double>> solver(solver_control);
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 * 
 * @endcode
 * 
 * We have made one addition, though: since we suppress output from the
 * linear solvers, we have to print the number of iterations by hand.
 * 
 * @code
 *   std::cout << "   " << solver_control.last_step()
 *             << " CG iterations needed to obtain convergence." << std::endl;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4output_results"></a> 
 * <h4>Step4::output_results</h4>
 * 

 * 
 * This function also does what the respective one did in step-3. No changes
 * here for dimension independence either.
 * 

 * 
 * Since the program will run both 2d and 3d versions of the Laplace solver,
 * we use the dimension in the filename to generate distinct filenames for
 * each run (in a better program, one would check whether <code>dim</code> can
 * have other values than 2 or 3, but we neglect this here for the sake of
 * brevity).
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::output_results() const
 * {
 *   DataOut<dim> data_out;
 * 
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 * 
 *   data_out.build_patches();
 * 
 *   std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
 *   data_out.write_vtk(output);
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4run"></a> 
 * <h4>Step4::run</h4>
 * 

 * 
 * This is the function which has the top-level control over everything. Apart
 * from one line of additional output, it is the same as for the previous
 * example.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::run()
 * {
 *   std::cout << "Solving problem in " << dim << " space dimensions."
 *             << std::endl;
 * 
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
 * And this is the main function. It also looks mostly like in step-3, but if
 * you look at the code below, note how we first create a variable of type
 * <code>Step4@<2@></code> (forcing the compiler to compile the class template
 * with <code>dim</code> replaced by <code>2</code>) and run a 2d simulation,
 * and then we do the whole thing over in 3d.
 * 

 * 
 * In practice, this is probably not what you would do very frequently (you
 * probably either want to solve a 2d problem, or one in 3d, but not both at
 * the same time). However, it demonstrates the mechanism by which we can
 * simply change which dimension we want in a single place, and thereby force
 * the compiler to recompile the dimension independent class templates for the
 * dimension we request. The emphasis here lies on the fact that we only need
 * to change a single place. This makes it rather trivial to debug the program
 * in 2d where computations are fast, and then switch a single place to a 3 to
 * run the much more computing intensive program in 3d for `real'
 * computations.
 * 

 * 
 * Each of the two blocks is enclosed in braces to make sure that the
 * <code>laplace_problem_2d</code> variable goes out of scope (and releases
 * the memory it holds) before we move on to allocate memory for the 3d
 * case. Without the additional braces, the <code>laplace_problem_2d</code>
 * variable would only be destroyed at the end of the function, i.e. after
 * running the 3d problem, and would needlessly hog memory while the 3d run
 * could actually use it.
 * 
 * @code
 * int main()
 * {
 *   {
 *     Step4<2> laplace_problem_2d;
 *     laplace_problem_2d.run();
 *   }
 * 
 *   {
 *     Step4<3> laplace_problem_3d;
 *     laplace_problem_3d.run();
 *   }
 * 
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>



The output of the program looks as follows (the number of iterations
may vary by one or two, depending on your computer, since this is
often dependent on the round-off accuracy of floating point
operations, which differs between processors):
@code
Solving problem in 2 space dimensions.
   Number of active cells: 256
   Total number of cells: 341
   Number of degrees of freedom: 289
   26 CG iterations needed to obtain convergence.
Solving problem in 3 space dimensions.
   Number of active cells: 4096
   Total number of cells: 4681
   Number of degrees of freedom: 4913
   30 CG iterations needed to obtain convergence.
@endcode
It is obvious that in three spatial dimensions the number of cells and
therefore also the number of degrees of freedom is
much higher. What cannot be seen here, is that besides this higher
number of rows and columns in the matrix, there are also significantly
more entries per row of the matrix in three space
dimensions. Together, this leads to a much higher numerical effort for
solving the system of equation, which you can feel in the run time of the two
solution steps when you actually run the program.



The program produces two files: <code>solution-2d.vtk</code> and
<code>solution-3d.vtk</code>, which can be viewed using the programs
VisIt or Paraview (in case you do not have these programs, you can easily
change the
output format in the program to something which you can view more
easily). Visualizing solutions is a bit of an art, but it can also be fun, so
you should play around with your favorite visualization tool to get familiar
with its functionality. Here's what I have come up with for the 2d solution:

<p align="center">
  <img src="https://www.dealii.org/images/steps/developer/step-4.solution-2d.png" alt="">
</p>

(@dealiiVideoLectureSeeAlso{11,32})
The picture shows the solution of the problem under consideration as
a 3D plot. As can be seen, the solution is almost flat in the interior
of the domain and has a higher curvature near the boundary. This, of
course, is due to the fact that for Laplace's equation the curvature
of the solution is equal to the right hand side and that was chosen as
a quartic polynomial which is nearly zero in the interior and is only
rising sharply when approaching the boundaries of the domain; the
maximal values of the right hand side function are at the corners of
the domain, where also the solution is moving most rapidly.
It is also nice to see that the solution follows the desired quadratic
boundary values along the boundaries of the domain.
It can also be useful to verify a computed solution against an analytical
solution. For an explanation of this technique, see step-7.

On the other hand, even though the picture does not show the mesh lines
explicitly, you can see them as little kinks in the solution. This clearly
indicates that the solution hasn't been computed to very high accuracy and
that to get a better solution, we may have to compute on a finer mesh.

In three spatial dimensions, visualization is a bit more difficult. The left
picture shows the solution and the mesh it was computed on on the surface of
the domain. This is nice, but it has the drawback that it completely hides
what is happening on the inside. The picture on the right is an attempt at
visualizing the interior as well, by showing surfaces where the solution has
constant values (as indicated by the legend at the top left). Isosurface
pictures look best if one makes the individual surfaces slightly transparent
so that it is possible to see through them and see what's behind.

<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.solution-3d.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.contours-3d.png" alt="">
    </td>
  </tr>
</table>

@note
A final remark on visualization: the idea of visualization is to give insight,
which is not the same as displaying information. In particular, it is easy to
overload a picture with information, but while it shows more information it
makes it also more difficult to glean insight. As an example, the program I
used to generate these pictures, VisIt, by default puts tick marks on every
axis, puts a big fat label "X Axis" on the $x$ axis and similar for the other
axes, shows the file name from which the data was taken in the top left and
the name of the user doing so and the time and date on the bottom right. None
of this is important
here: the axes are equally easy to make out because the tripod at the bottom
left is still visible, and we know from the program that the domain is
$[-1,1]^3$, so there is no need for tick marks. As a consequence, I have
switched off all the extraneous stuff in the picture: the art of visualization
is to reduce the picture to those parts that are important to see what one
wants to see, but no more.



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>



Essentially the possibilities for playing around with the program are the same
as for the previous one, except that they will now also apply to the 3d
case. For inspiration read up on <a href="step_3.html#extensions"
target="body">possible extensions in the documentation of step 3</a>.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-4.cc"
*/
