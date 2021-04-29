/**
@page step_10 The step-10 tutorial program
This tutorial depends on step-7.

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



This is a rather short example which only shows some aspects of using
higher order mappings. By <em>mapping</em> we mean the transformation
between the unit cell (i.e. the unit line, square, or cube) to the
cells in real space. In all the previous examples, we have implicitly
used linear or d-linear mappings; you will not have noticed this at
all, since this is what happens if you do not do anything
special. However, if your domain has curved boundaries, there are
cases where the piecewise linear approximation of the boundary
(i.e. by straight line segments) is not sufficient, and you want that
your computational domain is an approximation to the real domain using
curved boundaries as well. If the boundary approximation uses
piecewise quadratic parabolas to approximate the true boundary, then
we say that this is a quadratic or $Q_2$ approximation. If we
use piecewise graphs of cubic polynomials, then this is a $Q_3$
approximation, and so on.



For some differential equations, it is known that piecewise linear
approximations of the boundary, i.e. $Q_1$ mappings, are not
sufficient if the boundary of the exact domain is curved. Examples are the
biharmonic equation using $C^1$ elements, or the Euler
equations of gas dynamics on domains with curved reflective boundaries. In these cases,
it is necessary to compute the integrals using a higher order
mapping. If we do not use such a higher
order mapping, the order of approximation of the boundary dominates
the order of convergence of the entire numerical scheme, irrespective
of the order of convergence of the discretization in the interior of
the domain.



Rather than demonstrating the use of higher order mappings with one of
these more complicated examples, we do only a brief computation:
calculating the value of $\pi=3.141592653589793238462643\ldots$ by two
different methods.



The first method uses a triangulated approximation of the circle with unit
radius and integrates a unit magnitude constant function ($f = 1$) over it. Of
course, if the domain were the exact unit circle, then the area would be $\pi$,
but since we only use an approximation by piecewise polynomial segments, the
value of the area we integrate over is not exactly $\pi$. However, it is known
that as we refine the triangulation, a $Q_p$ mapping approximates the boundary
with an order $h^{p+1}$, where $h$ is the mesh size. We will check the values
of the computed area of the circle and their convergence towards $\pi$ under
mesh refinement for different mappings. We will also find a convergence
behavior that is surprising at first, but has a good explanation.



The second method works similarly, but this time does not use the area
of the triangulated unit circle, but rather its perimeter. $\pi$ is then
approximated by half of the perimeter, as we choose the radius equal to one.


@note This tutorial shows in essence how to choose a particular
mapping for integrals, by attaching a particular geometry to the
triangulation (as had already been done in step-1, for example) and
then passing a mapping argument to the FEValues class that is used for
all integrals in deal.II. The geometry we choose is a circle, for
which deal.II already has a class (SphericalManifold) that can be
used. If you want to define your own geometry, for example because it
is complicated and cannot be described by the classes already
available in deal.II, you will want to read through step-53.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The first of the following include files are probably well-known by now and
 * need no further explanation.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * @endcode
 * 
 * This include file is new. Even if we are not solving a PDE in this tutorial,
 * we want to use a dummy finite element with zero degrees of freedoms provided
 * by the FE_Nothing class.
 * 
 * @code
 * #include <deal.II/fe/fe_nothing.h>
 * 
 * @endcode
 * 
 * The following header file is also new: in it, we declare the MappingQ class
 * which we will use for polynomial mappings of arbitrary order:
 * 
 * @code
 * #include <deal.II/fe/mapping_q.h>
 * 
 * @endcode
 * 
 * And this again is C++:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * #include <cmath>
 * 
 * @endcode
 * 
 * The last step is as in previous programs:
 * 
 * @code
 * namespace Step10
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Now, as we want to compute the value of $\pi$, we have to compare to
 * something. These are the first few digits of $\pi$, which we define
 * beforehand for later use. Since we would like to compute the difference
 * between two numbers which are quite accurate, with the accuracy of the
 * computed approximation to $\pi$ being in the range of the number of
 * digits which a double variable can hold, we rather declare the reference
 * value as a <code>long double</code> and give it a number of extra digits:
 * 
 * @code
 *   const long double pi = 3.141592653589793238462643L;
 * 
 * 
 * 
 * @endcode
 * 
 * Then, the first task will be to generate some output. Since this program
 * is so small, we do not employ object oriented techniques in it and do not
 * declare classes (although, of course, we use the object oriented features
 * of the library). Rather, we just pack the functionality into separate
 * functions. We make these functions templates on the number of space
 * dimensions to conform to usual practice when using deal.II, although we
 * will only use them for two space dimensions and throw an exception when
 * attempted to use for any other spatial dimension.
 *   

 * 
 * The first of these functions just generates a triangulation of a circle
 * (hyperball) and outputs the $Q_p$ mapping of its cells for different values
 * of <code>p</code>. Then, we refine the grid once and do so again.
 * 
 * @code
 *   template <int dim>
 *   void gnuplot_output()
 *   {
 *     std::cout << "Output of grids into gnuplot files:" << std::endl
 *               << "===================================" << std::endl;
 * 
 * @endcode
 * 
 * So first generate a coarse triangulation of the circle and associate a
 * suitable boundary description to it. By default,
 * GridGenerator::hyper_ball attaches a SphericalManifold to the boundary
 * (and uses FlatManifold for the interior) so we simply call that
 * function and move on:
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     GridGenerator::hyper_ball(triangulation);
 * 
 * @endcode
 * 
 * Then alternate between generating output on the current mesh
 * for $Q_1$, $Q_2$, and $Q_3$ mappings, and (at the end of the
 * loop body) refining the mesh once globally.
 * 
 * @code
 *     for (unsigned int refinement = 0; refinement < 2; ++refinement)
 *       {
 *         std::cout << "Refinement level: " << refinement << std::endl;
 * 
 *         std::string filename_base = "ball_" + std::to_string(refinement);
 * 
 *         for (unsigned int degree = 1; degree < 4; ++degree)
 *           {
 *             std::cout << "Degree = " << degree << std::endl;
 * 
 * @endcode
 * 
 * For this, first set up an object describing the mapping. This
 * is done using the MappingQ class, which takes as
 * argument to the constructor the polynomial degree which it
 * shall use.
 * 
 * @code
 *             const MappingQ<dim> mapping(degree);
 * @endcode
 * 
 * As a side note, for a piecewise linear mapping, you
 * could give a value of <code>1</code> to the constructor
 * of MappingQ, but there is also a class MappingQ1 that
 * achieves the same effect. Historically, it did a lot of
 * things in a simpler way than MappingQ but is today just
 * a wrapper around the latter. It is, however, still the
 * class that is used implicitly in many places of the
 * library if you do not specify another mapping
 * explicitly.
 * 

 * 
 * 

 * 
 * In order to actually write out the present grid with this
 * mapping, we set up an object which we will use for output. We
 * will generate Gnuplot output, which consists of a set of lines
 * describing the mapped triangulation. By default, only one line
 * is drawn for each face of the triangulation, but since we want
 * to explicitly see the effect of the mapping, we want to have
 * the faces in more detail. This can be done by passing the
 * output object a structure which contains some flags. In the
 * present case, since Gnuplot can only draw straight lines, we
 * output a number of additional points on the faces so that each
 * face is drawn by 30 small lines instead of only one. This is
 * sufficient to give us the impression of seeing a curved line,
 * rather than a set of straight lines.
 * 
 * @code
 *             GridOut               grid_out;
 *             GridOutFlags::Gnuplot gnuplot_flags(false, 60);
 *             grid_out.set_flags(gnuplot_flags);
 * 
 * @endcode
 * 
 * Finally, generate a filename and a file for output:
 * 
 * @code
 *             std::string filename =
 *               filename_base + "_mapping_q_" + std::to_string(degree) + ".dat";
 *             std::ofstream gnuplot_file(filename);
 * 
 * @endcode
 * 
 * Then write out the triangulation to this file. The last
 * argument of the function is a pointer to a mapping object. This
 * argument has a default value, and if no value is given a simple
 * MappingQ1 object is taken, which we briefly
 * described above. This would then result in a piecewise linear
 * approximation of the true boundary in the output.
 * 
 * @code
 *             grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping);
 *           }
 *         std::cout << std::endl;
 * 
 * @endcode
 * 
 * At the end of the loop, refine the mesh globally.
 * 
 * @code
 *         triangulation.refine_global();
 *       }
 *   }
 * 
 * @endcode
 * 
 * Now we proceed with the main part of the code, the approximation of
 * $\pi$. The area of a circle is of course given by $\pi r^2$, so having a
 * circle of radius 1, the area represents just the number that is searched
 * for. The numerical computation of the area is performed by integrating
 * the constant function of value 1 over the whole computational domain,
 * i.e. by computing the areas $\int_K 1 dx=\int_{\hat K} 1
 * \ \textrm{det}\ J(\hat x) d\hat x \approx \sum_i \textrm{det}
 * \ J(\hat x_i)w(\hat x_i)$,
 * where the sum extends over all quadrature points on all active cells in
 * the triangulation, with $w(x_i)$ being the weight of quadrature point
 * $x_i$. The integrals on each cell are approximated by numerical
 * quadrature, hence the only additional ingredient we need is to set up a
 * FEValues object that provides the corresponding `JxW` values of each
 * cell. (Note that `JxW` is meant to abbreviate <i>Jacobian determinant
 * times weight</i>; since in numerical quadrature the two factors always
 * occur at the same places, we only offer the combined quantity, rather
 * than two separate ones.) We note that here we won't use the FEValues
 * object in its original purpose, i.e. for the computation of values of
 * basis functions of a specific finite element at certain quadrature
 * points. Rather, we use it only to gain the `JxW` at the quadrature
 * points, irrespective of the (dummy) finite element we will give to the
 * constructor of the FEValues object. The actual finite element given to
 * the FEValues object is not used at all, so we could give any.
 * 
 * @code
 *   template <int dim>
 *   void compute_pi_by_area()
 *   {
 *     std::cout << "Computation of Pi by the area:" << std::endl
 *               << "==============================" << std::endl;
 * 
 * @endcode
 * 
 * For the numerical quadrature on all cells we employ a quadrature rule
 * of sufficiently high degree. We choose QGauss that is of order 8 (4
 * points), to be sure that the errors due to numerical quadrature are of
 * higher order than the order (maximal 6) that will occur due to the
 * order of the approximation of the boundary, i.e. the order of the
 * mappings employed. Note that the integrand, the Jacobian determinant,
 * is not a polynomial function (rather, it is a rational one), so we do
 * not use Gauss quadrature in order to get the exact value of the
 * integral as done often in finite element computations, but could as
 * well have used any quadrature formula of like order instead.
 * 
 * @code
 *     const QGauss<dim> quadrature(4);
 * 
 * @endcode
 * 
 * Now start by looping over polynomial mapping degrees=1..4:
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree)
 *       {
 *         std::cout << "Degree = " << degree << std::endl;
 * 
 * @endcode
 * 
 * First generate the triangulation, the boundary and the mapping
 * object as already seen.
 * 
 * @code
 *         Triangulation<dim> triangulation;
 *         GridGenerator::hyper_ball(triangulation);
 * 
 *         const MappingQ<dim> mapping(degree);
 * 
 * @endcode
 * 
 * We now create a finite element. Unlike the rest of the example
 * programs, we do not actually need to do any computations with shape
 * functions; we only need the `JxW` values from an FEValues
 * object. Hence we use the special finite element class FE_Nothing
 * which has exactly zero degrees of freedom per cell (as the name
 * implies, the local basis on each cell is the empty set). A more
 * typical usage of FE_Nothing is shown in step-46.
 * 
 * @code
 *         const FE_Nothing<dim> fe;
 * 
 * @endcode
 * 
 * Likewise, we need to create a DoFHandler object. We do not actually
 * use it, but it will provide us with `active_cell_iterators` that
 * are needed to reinitialize the FEValues object on each cell of the
 * triangulation.
 * 
 * @code
 *         DoFHandler<dim> dof_handler(triangulation);
 * 
 * @endcode
 * 
 * Now we set up the FEValues object, giving the Mapping, the dummy
 * finite element and the quadrature object to the constructor,
 * together with the update flags asking for the `JxW` values at the
 * quadrature points only. This tells the FEValues object that it
 * needs not compute other quantities upon calling the
 * <code>reinit</code> function, thus saving computation time.
 *         

 * 
 * The most important difference in the construction of the FEValues
 * object compared to previous example programs is that we pass a
 * mapping object as first argument, which is to be used in the
 * computation of the mapping from unit to real cell. In previous
 * examples, this argument was omitted, resulting in the implicit use
 * of an object of type MappingQ1.
 * 
 * @code
 *         FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values);
 * 
 * @endcode
 * 
 * We employ an object of the ConvergenceTable class to store all
 * important data like the approximated values for $\pi$ and the error
 * with respect to the true value of $\pi$. We will also use functions
 * provided by the ConvergenceTable class to compute convergence rates
 * of the approximations to $\pi$.
 * 
 * @code
 *         ConvergenceTable table;
 * 
 * @endcode
 * 
 * Now we loop over several refinement steps of the triangulation.
 * 
 * @code
 *         for (unsigned int refinement = 0; refinement < 6;
 *              ++refinement, triangulation.refine_global(1))
 *           {
 * @endcode
 * 
 * In this loop we first add the number of active cells of the
 * current triangulation to the table. This function automatically
 * creates a table column with superscription `cells`, in case
 * this column was not created before.
 * 
 * @code
 *             table.add_value("cells", triangulation.n_active_cells());
 * 
 * @endcode
 * 
 * Then we distribute the degrees of freedom for the dummy finite
 * element. Strictly speaking we do not need this function call in
 * our special case but we call it to make the DoFHandler happy --
 * otherwise it would throw an assertion in the FEValues::reinit
 * function below.
 * 
 * @code
 *             dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * We define the variable area as `long double` like we did for
 * the `pi` variable before.
 * 
 * @code
 *             long double area = 0;
 * 
 * @endcode
 * 
 * Now we loop over all cells, reinitialize the FEValues object
 * for each cell, and add up all the `JxW` values for this cell to
 * `area`...
 * 
 * @code
 *             for (const auto &cell : dof_handler.active_cell_iterators())
 *               {
 *                 fe_values.reinit(cell);
 *                 for (unsigned int i = 0; i < fe_values.n_quadrature_points; ++i)
 *                   area += static_cast<long double>(fe_values.JxW(i));
 *               }
 * 
 * @endcode
 * 
 * ...and store the resulting area values and the errors in the
 * table. We need a static cast to double as there is no
 * add_value(string, long double) function implemented. Note that
 * this also concerns the second call as the <code>fabs</code>
 * function in the <code>std</code> namespace is overloaded on its
 * argument types, so there exists a version taking and returning
 * a <code>long double</code>, in contrast to the global namespace
 * where only one such function is declared (which takes and
 * returns a double).
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(area));
 *             table.add_value("error", static_cast<double>(std::fabs(area - pi)));
 *           }
 * 
 * @endcode
 * 
 * We want to compute the convergence rates of the `error`
 * column. Therefore we need to omit the other columns from the
 * convergence rate evaluation before calling
 * `evaluate_all_convergence_rates`
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells");
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi");
 *         table.evaluate_all_convergence_rates(
 *           ConvergenceTable::reduction_rate_log2);
 * 
 * @endcode
 * 
 * Finally we set the precision and scientific mode for output of some
 * of the quantities...
 * 
 * @code
 *         table.set_precision("eval.pi", 16);
 *         table.set_scientific("error", true);
 * 
 * @endcode
 * 
 * ...and write the whole table to std::cout.
 * 
 * @code
 *         table.write_text(std::cout);
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * The following, second function also computes an approximation of $\pi$
 * but this time via the perimeter $2\pi r$ of the domain instead of the
 * area. This function is only a variation of the previous function. So we
 * will mainly give documentation for the differences.
 * 
 * @code
 *   template <int dim>
 *   void compute_pi_by_perimeter()
 *   {
 *     std::cout << "Computation of Pi by the perimeter:" << std::endl
 *               << "===================================" << std::endl;
 * 
 * @endcode
 * 
 * We take the same order of quadrature but this time a `dim-1`
 * dimensional quadrature as we will integrate over (boundary) lines
 * rather than over cells.
 * 
 * @code
 *     const QGauss<dim - 1> quadrature(4);
 * 
 * @endcode
 * 
 * We loop over all degrees, create the triangulation, the boundary, the
 * mapping, the dummy finite element and the DoFHandler object as seen
 * before.
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree)
 *       {
 *         std::cout << "Degree = " << degree << std::endl;
 *         Triangulation<dim> triangulation;
 *         GridGenerator::hyper_ball(triangulation);
 * 
 *         const MappingQ<dim>   mapping(degree);
 *         const FE_Nothing<dim> fe;
 * 
 *         DoFHandler<dim> dof_handler(triangulation);
 * 
 * @endcode
 * 
 * Then we create a FEFaceValues object instead of a FEValues object
 * as in the previous function. Again, we pass a mapping as first
 * argument.
 * 
 * @code
 *         FEFaceValues<dim> fe_face_values(mapping,
 *                                          fe,
 *                                          quadrature,
 *                                          update_JxW_values);
 *         ConvergenceTable  table;
 * 
 *         for (unsigned int refinement = 0; refinement < 6;
 *              ++refinement, triangulation.refine_global(1))
 *           {
 *             table.add_value("cells", triangulation.n_active_cells());
 * 
 *             dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * Now we run over all cells and over all faces of each cell. Only
 * the contributions of the `JxW` values on boundary faces are
 * added to the long double variable `perimeter`.
 * 
 * @code
 *             long double perimeter = 0;
 *             for (const auto &cell : dof_handler.active_cell_iterators())
 *               for (const auto &face : cell->face_iterators())
 *                 if (face->at_boundary())
 *                   {
 * @endcode
 * 
 * We reinit the FEFaceValues object with the cell
 * iterator and the number of the face.
 * 
 * @code
 *                     fe_face_values.reinit(cell, face);
 *                     for (unsigned int i = 0;
 *                          i < fe_face_values.n_quadrature_points;
 *                          ++i)
 *                       perimeter +=
 *                         static_cast<long double>(fe_face_values.JxW(i));
 *                   }
 * @endcode
 * 
 * Then store the evaluated values in the table...
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(perimeter / 2.0L));
 *             table.add_value(
 *               "error", static_cast<double>(std::fabs(perimeter / 2.0L - pi)));
 *           }
 * 
 * @endcode
 * 
 * ...and end this function as we did in the previous one:
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells");
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi");
 *         table.evaluate_all_convergence_rates(
 *           ConvergenceTable::reduction_rate_log2);
 * 
 *         table.set_precision("eval.pi", 16);
 *         table.set_scientific("error", true);
 * 
 *         table.write_text(std::cout);
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * } // namespace Step10
 * 
 * 
 * @endcode
 * 
 * The following main function just calls the above functions in the order of
 * their appearance. Apart from this, it looks just like the main functions of
 * previous tutorial programs.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       std::cout.precision(16);
 * 
 *       const unsigned int dim = 2;
 * 
 *       Step10::gnuplot_output<dim>();
 * 
 *       Step10::compute_pi_by_area<dim>();
 *       Step10::compute_pi_by_perimeter<dim>();
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



The program performs two tasks, the first being to generate a
visualization of the mapped domain, the second to compute pi by the
two methods described. Let us first take a look at the generated
graphics. They are generated in Gnuplot format, and can be viewed with
the commands
@code
set style data lines
set size ratio -1
unset key
unset xtics
unset ytics
plot [-1:1][-1:1] "ball_0_mapping_q_1.dat" lw 4 lt rgb "black"
@endcode
or using one of the other filenames. The second line makes sure that
the aspect ratio of the generated output is actually 1:1, i.e. a
circle is drawn as a circle on your screen, rather than as an
ellipse. The third line switches off the key in the graphic, as that
will only print information (the filename) which is not that important
right now. Similarly, the fourth and fifth disable tick marks. The plot
is then generated with a specific line width ("lw", here set to 4)
and line type ("lt", here chosen by saying that the line should be
drawn using the RGB color "black").

The following table shows the triangulated computational domain for $Q_1$,
$Q_2$, and $Q_3$ mappings, for the original coarse grid (left), and a once
uniformly refined grid (right).

<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q1.svg"
         alt="Five-cell discretization of the disk."
         width="400" height="400">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q1.svg"
         alt="20-cell discretization of the disk (i.e., five cells
              refined once)."
         width="400" height="400">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q2.svg"
         alt="Five-cell discretization of the disk with quadratic edges. The
              boundary is nearly indistinguishable from the actual circle."
         width="400" height="400"
         >
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q2.svg"
         alt="20-cell discretization with quadratic edges."
         width="400" height="400">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q3.svg"
         alt="Five-cell discretization of the disk with cubic edges. The
              boundary is nearly indistinguishable from the actual circle."
         width="400" height="400">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q3.svg"
         alt="20-cell discretization with cubic edges."
         width="400" height="400">
  </div>
</div>

These pictures show the obvious advantage of higher order mappings: they
approximate the true boundary quite well also on rather coarse meshes. To
demonstrate this a little further, here is part of the upper right quarter
circle of the coarse meshes with $Q_2$ and $Q_3$ mappings, where the dashed
red line marks the actual circle:

<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q2.svg"
         alt="Close-up of quadratic discretization. The distance between the
         quadratic interpolant and the actual circle is small."
         width="400" height="400">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q3.svg"
         alt="Close-up of cubic discretization. The distance between the
         cubic interpolant and the actual circle is very small."
         width="400" height="400">
  </div>
</div>

Obviously the quadratic mapping approximates the boundary quite well,
while for the cubic mapping the difference between approximated domain
and true one is hardly visible already for the coarse grid. You can
also see that the mapping only changes something at the outer
boundaries of the triangulation. In the interior, all lines are still
represented by linear functions, resulting in additional computations
only on cells at the boundary. Higher order mappings are therefore
usually not noticeably slower than lower order ones, because the
additional computations are only performed on a small subset of all
cells.



The second purpose of the program was to compute the value of pi to
good accuracy. This is the output of this part of the program:
@code
Output of grids into gnuplot files:
===================================
Refinement level: 0
Degree = 1
Degree = 2
Degree = 3

Refinement level: 1
Degree = 1
Degree = 2
Degree = 3

Computation of Pi by the area:
==============================
Degree = 1
cells      eval.pi            error
    5 1.9999999999999993 1.1416e+00    -
   20 2.8284271247461890 3.1317e-01 1.87
   80 3.0614674589207174 8.0125e-02 1.97
  320 3.1214451522580511 2.0148e-02 1.99
 1280 3.1365484905459380 5.0442e-03 2.00
 5120 3.1403311569547516 1.2615e-03 2.00

Degree = 2
cells      eval.pi            error
    5 3.1045694996615860 3.7023e-02    -
   20 3.1391475703122267 2.4451e-03 3.92
   80 3.1414377167038290 1.5494e-04 3.98
  320 3.1415829366419006 9.7169e-06 4.00
 1280 3.1415920457576898 6.0783e-07 4.00
 5120 3.1415926155921117 3.7998e-08 4.00

Degree = 3
cells      eval.pi            error
    5 3.1410033851499288 5.8927e-04    -
   20 3.1415830393583839 9.6142e-06 5.94
   80 3.1415925017363797 1.5185e-07 5.98
  320 3.1415926512106696 2.3791e-09 6.00
 1280 3.1415926535525927 3.7200e-11 6.00
 5120 3.1415926535892100 5.8302e-13 6.00

Degree = 4
cells      eval.pi            error
    5 3.1415871927401131 5.4608e-06    -
   20 3.1415926314742428 2.2116e-08 7.95
   80 3.1415926535026202 8.7173e-11 7.99
  320 3.1415926535894498 3.4350e-13 7.99
 1280 3.1415926535897896 3.4671e-15 6.63
 5120 3.1415926535897909 2.4009e-15 0.53

Computation of Pi by the perimeter:
===================================
Degree = 1
cells      eval.pi            error
    5 2.8284271247461898 3.1317e-01    -
   20 3.0614674589207178 8.0125e-02 1.97
   80 3.1214451522580520 2.0148e-02 1.99
  320 3.1365484905459389 5.0442e-03 2.00
 1280 3.1403311569547525 1.2615e-03 2.00
 5120 3.1412772509327724 3.1540e-04 2.00

Degree = 2
cells      eval.pi            error
    5 3.1248930668550594 1.6700e-02    -
   20 3.1404050605605449 1.1876e-03 3.81
   80 3.1415157631807009 7.6890e-05 3.95
  320 3.1415878042798613 4.8493e-06 3.99
 1280 3.1415923498174534 3.0377e-07 4.00
 5120 3.1415926345931995 1.8997e-08 4.00

Degree = 3
cells      eval.pi            error
    5 3.1414940401456048 9.8613e-05    -
   20 3.1415913432549156 1.3103e-06 6.23
   80 3.1415926341726910 1.9417e-08 6.08
  320 3.1415926532906897 2.9910e-10 6.02
 1280 3.1415926535851355 4.6578e-12 6.00
 5120 3.1415926535897190 7.4216e-14 5.97

Degree = 4
cells      eval.pi            error
    5 3.1415921029432572 5.5065e-07     -
   20 3.1415926513737595 2.2160e-09  7.96
   80 3.1415926535810712 8.7222e-12  7.99
  320 3.1415926535897576 3.5525e-14  7.94
 1280 3.1415926535897936 4.6729e-16  6.25
 5120 3.1415926535897918 1.4929e-15 -1.68
@endcode



One of the immediate observations from the output is that in all cases the
values converge quickly to the true value of
$\pi=3.141592653589793238462643$. Note that for the $Q_4$ mapping, we are
already in the regime of roundoff errors and the convergence rate levels off,
which is already quite a lot. However, also note that for the $Q_1$ mapping,
even on the finest grid the accuracy is significantly worse than on the coarse
grid for a $Q_3$ mapping!



The last column of the output shows the convergence order, in powers of the
mesh width $h$. In the introduction, we had stated that the convergence order
for a $Q_p$ mapping should be $h^{p+1}$. However, in the example shown, the
order is rather $h^{2p}$! This at first surprising fact is explained by the
properties of the $Q_p$ mapping. At order <i>p</i>, it uses support points
that are based on the <i>p</i>+1 point Gauss-Lobatto quadrature rule that
selects the support points in such a way that the quadrature rule converges at
order 2<i>p</i>. Even though these points are here only used for interpolation
of a <i>p</i>th order polynomial, we get a superconvergence effect when
numerically evaluating the integral, resulting in the observed high order of
convergence. (This effect is also discussed in detail in the following
publication: A. Bonito, A. Demlow, and J. Owen: "A priori error
estimates for finite element approximations to eigenvalues and
eigenfunctions of the Laplace-Beltrami operator", submitted, 2018.)
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-10.cc"
*/
