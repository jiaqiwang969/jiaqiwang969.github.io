/**
@page step_63 The step-63 tutorial program
This tutorial depends on step-16.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Equation">Equation</a>
      <ul>
        <li><a href="#Streamlinediffusion">Streamline diffusion</a>
      </ul>
        <li><a href="#Smoothers">Smoothers</a>
        <li><a href="#Testproblem">Test problem</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#MeshWorkerdata">MeshWorker data</a>
        <li><a href="#Problemparameters">Problem parameters</a>
        <li><a href="#Cellpermutations">Cell permutations</a>
        <li><a href="#Righthandsideandboundaryvalues">Right-hand side and boundary values</a>
        <li><a href="#Streamlinediffusionimplementation">Streamline diffusion implementation</a>
        <li><a href="#codeAdvectionProlemcodeclass"><code>AdvectionProlem</code> class</a>
      <ul>
        <li><a href="#codeAdvectionProblemsetup_systemcode"><code>AdvectionProblem::setup_system()</code></a>
        <li><a href="#codeAdvectionProblemassemble_cellcode"><code>AdvectionProblem::assemble_cell()</code></a>
        <li><a href="#codeAdvectionProblemassemble_system_and_multigridcode"><code>AdvectionProblem::assemble_system_and_multigrid()</code></a>
        <li><a href="#codeAdvectionProblemsetup_smoothercode"><code>AdvectionProblem::setup_smoother()</code></a>
        <li><a href="#codeAdvectionProblemsolvecode"><code>AdvectionProblem::solve()</code></a>
        <li><a href="#codeAdvectionProblemoutput_resultscode"><code>AdvectionProblem::output_results()</code></a>
        <li><a href="#codeAdvectionProblemruncode"><code>AdvectionProblem::run()</code></a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#GMRESIterationNumbers"> GMRES Iteration Numbers </a>
      <ul>
        <li><a href="#DoFCellRenumbering"> DoF/Cell Renumbering </a>
        <li><a href="#Pointvsblocksmoothers"> Point vs. block smoothers </a>
      </ul>
        <li><a href="#Cost"> Cost </a>
        <li><a href="#Additionalpoints"> Additional points </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#ConstantiterationsforQsub5sub"> Constant iterations for Q<sub>5</sub> </a>
        <li><a href="#Effectivenessofrenumberingforchangingepsilon"> Effectiveness of renumbering for changing epsilon </a>
        <li><a href="#Meshadaptivity"> Mesh adaptivity </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Thomas C. Clevenger and Timo Heister.

The creation of this tutorial was partially supported by NSF Award
DMS-1522191, DMS-1901529, OAC-1835452, by the Computational
Infrastructure in Geodynamics initiative (CIG), through the NSF under
Award EAR-0949446 and EAR-1550901 and The University of California -
Davis.
</i>

@dealiiTutorialDOI{10.5281/zenodo.3382899,https://zenodo.org/badge/DOI/10.5281/zenodo.3382899.svg}

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This program solves an advection-diffusion problem using a geometric multigrid
(GMG) preconditioner. The basics of this preconditioner are discussed in step-16;
here we discuss the necessary changes needed for a non-symmetric
PDE. Additionally, we introduce the idea of block smoothing (as compared to
point smoothing in step-16), and examine the effects of DoF renumbering for
additive and multiplicative smoothers.

<a name="Equation"></a><h3>Equation</h3>

The advection-diffusion equation is given by
@f{align*}{
-\varepsilon \Delta u + \boldsymbol{\beta}\cdot \nabla u & = f &
\text{ in } \Omega\\
u &= g & \text{ on } \partial\Omega
@f}
where $\varepsilon>0$, $\boldsymbol{\beta}$ is the <i>advection
direction</i>, and $f$ is a source. A few notes:

1. If $\boldsymbol{\beta}=\boldsymbol{0}$, this is the Laplace equation solved in step-16
(and many other places).

2. If $\varepsilon=0$ then this is the stationary advection equation solved in
step-9.

3. One can define a dimensionless number for this problem, called the
<i>Peclet number</i>: $\mathcal{P} \dealcoloneq \frac{\|\boldsymbol{\beta}\|
L}{\varepsilon}$, where $L$ is the length scale of the domain. It
characterizes the kind of equation we are
considering: If $\mathcal{P}>1$, we say the problem is
<i>advection-dominated</i>, else if $\mathcal{P}<1$ we will say the problem is
<i>diffusion-dominated</i>.

For the discussion in this tutorial we will be concerned with
advection-dominated flow. This is the complicated case: We know that
for diffusion-dominated problems, the standard Galerkin method works
just fine, and we also know that simple multigrid methods such as
those defined in step-16 are very efficient. On the other hand, for
advection-dominated problems, the standard Galerkin approach leads to
oscillatory and unstable discretizations, and simple solvers are often
not very efficient. This tutorial program is therefore intended to
address both of these issues.


<a name="Streamlinediffusion"></a><h4>Streamline diffusion</h4>


Using the standard Galerkin finite element method, for suitable test
functions $v_h$, a discrete weak form of the PDE would read
@f{align*}{
a(u_h,v_h) = F(v_h)
@f}
where
@f{align*}{
a(u_h,v_h) &= (\varepsilon \nabla v_h,\, \nabla u_h) +
(v_h,\,\boldsymbol{\beta}\cdot \nabla u_h),\\
F(v_h) &= (v_h,\,f).
@f}

Unfortunately, one typically gets oscillatory solutions with this
approach. Indeed, the following error estimate can be shown for this
formulation:
@f{align*}{
\|\nabla (u-u_h)\| \leq (1+\mathcal{P}) \inf_{v_h} \|\nabla (u-v_h)\|.
@f}
The infimum on the right can be estimated as follows if the exact
solution is sufficiently smooth:
@f{align*}{
  \inf_{v_h} \|\nabla (u-v_h)\|.
  \le
  \|\nabla (u-I_h u)\|
  \le
  h^k
  C
  \|\nabla^k u)\|
@f}
where $k$ is the polynomial degree of the finite elements used. As a
consequence, we obtain the estimate
@f{align*}{
\|\nabla (u-u_h)\|
\leq (1+\mathcal{P}) C h^k
  \|\nabla^k u)\|.
@f}
In other words, the numerical solution will converge. On the other hand,
given the definition of $\mathcal{P}$ above, we have to expect poor
numerical solutions with a large error when $\varepsilon \ll
\|\boldsymbol{\beta}\| L$, i.e., if the problem has only a small
amount of diffusion.

To combat this, we will consider the new weak form
@f{align*}{
a(u_h,\,v_h) + \sum_K (-\varepsilon \Delta u_h +
\boldsymbol{\beta}\cdot \nabla u_h-f,\,\delta_K
\boldsymbol{\beta}\cdot \nabla v_h)_K = F(v_h)
@f}
where the sum is done over all cells $K$ with the inner product taken
for each cell, and $\delta_K$ is a cell-wise constant
stabilization parameter defined in
@cite john2006discontinuity.

Essentially, adding in the
discrete strong form residual enhances the coercivity of the bilinear
form $a(\cdot,\cdot)$ which increases the stability of the discrete
solution. This method is commonly referred to as <i>streamline
diffusion</i> or <i>SUPG</i> (streamline upwind/Petrov-Galerkin).


<a name="Smoothers"></a><h3>Smoothers</h3>


One of the goals of this tutorial is to expand from using a simple
(point-wise) Gauss-Seidel (SOR) smoother that is used in step-16
(class PreconditionSOR) on each level of the multigrid hierarchy.
The term "point-wise" is traditionally used in solvers to indicate that one
solves at one "grid point" at a time; for scalar problems, this means
to use a solver that updates one unknown of the linear
system at a time, keeping all of the others fixed; one would then
iterate over all unknowns in the problem and, once done, start over again
from the first unknown until these "sweeps" converge. Jacobi,
Gauss-Seidel, and SOR iterations can all be interpreted in this way.
In the context of multigrid, one does not think of these methods as
"solvers", but as "smoothers". As such, one is not interested in
actually solving the linear system. It is enough to remove the high-frequency
part of the residual for the multigrid method to work, because that allows
restricting the solution to a coarser mesh.  Therefore, one only does a few,
fixed number of "sweeps" over all unknowns. In the code in this
tutorial this is controlled by the "Smoothing steps" parameter.

But these methods are known to converge rather slowly when used as
solvers. While as multigrid smoothers, they are surprisingly good,
they can also be improved upon. In particular, we consider
"cell-based" smoothers here as well. These methods solve for all
unknowns on a cell at once, keeping all other unknowns fixed; they
then move on to the next cell, and so on and so forth. One can think
of them as "block" versions of Jacobi, Gauss-Seidel, or SOR, but
because degrees of freedom are shared among multiple cells, these
blocks overlap and the methods are in fact
best be explained within the framework of additive and multiplicative
Schwarz methods.

In contrast to step-16, our test problem contains an advective
term. Especially with a small diffusion constant $\varepsilon$, information is
transported along streamlines in the given advection direction. This means
that smoothers are likely to be more effective if they allow information to
travel in downstream direction within a single smoother
application. If we want to solve one unknown (or block of unknowns) at
a time in the order in which these unknowns (or blocks) are
enumerated, then this information propagation property
requires reordering degrees of freedom or cells (for the cell-based smoothers)
accordingly so that the ones further upstream are treated earlier
(have lower indices) and those further downstream are treated later
(have larger indices). The influence of the ordering will be visible
in the results section.

Let us now briefly define the smoothers used in this tutorial.
For a more detailed introduction, we refer to
@cite KanschatNotesIterative and the books @cite smith2004domain and @cite toselli2006domain.
A Schwarz
preconditioner requires a decomposition
@f{align*}{
V = \sum_{j=1}^J V_j
@f}
of our finite element space $V$. Each subproblem $V_j$ also has a Ritz
projection $P_j: V \rightarrow V_j$ based on the bilinear form
$a(\cdot,\cdot)$. This projection induces a local operator $A_j$ for each
subproblem $V_j$. If $\Pi_j:V\rightarrow V_j$ is the orthogonal projector onto
$V_j$, one can show $A_jP_j=\Pi_j^TA$.

With this we can define an <i>additive Schwarz preconditioner</i> for the
operator $A$ as
@f{align*}{
 B^{-1} = \sum_{j=1}^J P_j A^{-1} = \sum_{j=1}^J A_j^{-1} \Pi_j^T.
@f}
In other words, we project our solution into each subproblem, apply the
inverse of the subproblem $A_j$, and sum the contributions up over all $j$.

Note that one can interpret the point-wise (one unknown at a time)
Jacobi method as an additive
Schwarz method by defining a subproblem $V_j$ for each degree of
freedom. Then, $A_j^{-1}$ becomes a multiplication with the inverse of a
diagonal entry of $A$.

For the "Block Jacobi" method used in this tutorial, we define a subproblem
$V_j$ for each cell of the mesh on the current level. Note that we use a
continuous finite element, so these blocks are overlapping, as degrees of
freedom on an interface between two cells belong to both subproblems. The
logic for the Schwarz operator operating on the subproblems (in deal.II they
are called "blocks") is implemented in the class RelaxationBlock. The "Block
Jacobi" method is implemented in the class RelaxationBlockJacobi. Many
aspects of the class (for example how the blocks are defined and how to invert
the local subproblems $A_j$) can be configured in the smoother data, see
RelaxationBlock::AdditionalData and DoFTools::make_cell_patches() for details.

So far, we discussed additive smoothers where the updates can be applied
independently and there is no information flowing within a single smoother
application. A <i>multiplicative Schwarz preconditioner</i> addresses this
and is defined by
@f{align*}{
 B^{-1} = \left( I- \prod_{j=1}^J \left(I-P_j\right) \right) A^{-1}.
@f}
In contrast to above, the updates on the subproblems $V_j$ are applied
sequentially. This means that the update obtained when inverting the
subproblem $A_j$ is immediately used in $A_{j+1}$. This becomes
visible when writing out the project:
@f{align*}{
 B^{-1}
 =
 \left(
   I
   -
   \left(I-P_1\right)\left(I-P_2\right)\cdots\left(I-P_J\right)
 \right)
 A^{-1}
 =
   A^{-1}
   -
   \left[ \left(I-P_1\right)
   \left[ \left(I-P_2\right)\cdots
     \left[\left(I-P_J\right) A^{-1}\right] \cdots \right] \right]
@f}

When defining the sub-spaces $V_j$ as whole blocks of degrees of
freedom, this method is implemented in the class RelaxationBlockSOR and used when you
select "Block SOR" in this tutorial. The class RelaxationBlockSOR is also
derived from RelaxationBlock. As such, both additive and multiplicative
Schwarz methods are implemented in a unified framework.

Finally, let us note that the standard Gauss-Seidel (or SOR) method can be
seen as a multiplicative Schwarz method with a subproblem for each DoF.


<a name="Testproblem"></a><h3>Test problem</h3>


We will be considering the following test problem: $\Omega =
[-1,\,1]\times[-1,\,1]\backslash B_{0.3}(0)$, i.e., a square
with a circle of radius 0.3 centered at the
origin removed. In addition, we use $\varepsilon=0.005$, $\boldsymbol{\beta} =
[-\sin(\pi/6),\,\cos(\pi/6)]$, $f=0$, and Dirichlet boundary values
@f{align*}{
g = \left\{\begin{array}{ll} 1 & \text{if } x=-1 \text{ or } y=-1,\,x\geq 0.5 \\
0 & \text{otherwise} \end{array}\right.
@f}

The following figures depict the solutions with (left) and without
(right) streamline diffusion. Without streamline diffusion we see large
oscillations around the boundary layer, demonstrating the instability
of the standard Galerkin finite element method for this problem.

<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-63-solution.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-63-solution-no-sd.png" alt="">
    </td>
  </tr>
</table>
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
 * Typical files needed for standard deal.II:
 * 
 * @code
 * #include <deal.II/base/tensor_function.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/relaxation_block.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_out.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * Include all relevant multilevel files:
 * 
 * @code
 * #include <deal.II/multigrid/mg_constrained_dofs.h>
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * @endcode
 * 
 * C++:
 * 
 * @code
 * #include <algorithm>
 * #include <fstream>
 * #include <iostream>
 * #include <random>
 * 
 * @endcode
 * 
 * We will be using MeshWorker::mesh_loop functionality for assembling matrices:
 * 
 * @code
 * #include <deal.II/meshworker/mesh_loop.h>
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MeshWorkerdata"></a> 
 * <h3>MeshWorker data</h3>
 * 

 * 
 * As always, we will be putting everything related to this program
 * into a namespace of its own.
 * 

 * 
 * Since we will be using the MeshWorker framework, the first step is
 * to define the following structures needed by the assemble_cell()
 * function used by MeshWorker::mesh_loop(): `ScratchData`
 * contains an FEValues object which is needed for assembling
 * a cell's local contribution, while `CopyData` contains the
 * output from a cell's local contribution and necessary information
 * to copy that to the global system. (Their purpose is also explained
 * in the documentation of the WorkStream class.)
 * 
 * @code
 * namespace Step63
 * {
 *   using namespace dealii;
 * 
 *   template <int dim>
 *   struct ScratchData
 *   {
 *     ScratchData(const FiniteElement<dim> &fe,
 *                 const unsigned int        quadrature_degree)
 *       : fe_values(fe,
 *                   QGauss<dim>(quadrature_degree),
 *                   update_values | update_gradients | update_hessians |
 *                     update_quadrature_points | update_JxW_values)
 *     {}
 * 
 *     ScratchData(const ScratchData<dim> &scratch_data)
 *       : fe_values(scratch_data.fe_values.get_fe(),
 *                   scratch_data.fe_values.get_quadrature(),
 *                   update_values | update_gradients | update_hessians |
 *                     update_quadrature_points | update_JxW_values)
 *     {}
 * 
 *     FEValues<dim> fe_values;
 *   };
 * 
 * 
 * 
 *   struct CopyData
 *   {
 *     CopyData() = default;
 * 
 *     unsigned int level;
 *     unsigned int dofs_per_cell;
 * 
 *     FullMatrix<double>                   cell_matrix;
 *     Vector<double>                       cell_rhs;
 *     std::vector<types::global_dof_index> local_dof_indices;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Problemparameters"></a> 
 * <h3>Problem parameters</h3>
 * 

 * 
 * The second step is to define the classes that deal with run-time
 * parameters to be read from an input file.
 *   

 * 
 * We will use ParameterHandler to pass in parameters at runtime. The
 * structure `Settings` parses and stores the parameters to be queried
 * throughout the program.
 * 
 * @code
 *   struct Settings
 *   {
 *     enum DoFRenumberingStrategy
 *     {
 *       none,
 *       downstream,
 *       upstream,
 *       random
 *     };
 * 
 *     void get_parameters(const std::string &prm_filename);
 * 
 *     double                 epsilon;
 *     unsigned int           fe_degree;
 *     std::string            smoother_type;
 *     unsigned int           smoothing_steps;
 *     DoFRenumberingStrategy dof_renumbering;
 *     bool                   with_streamline_diffusion;
 *     bool                   output;
 *   };
 * 
 * 
 * 
 *   void Settings::get_parameters(const std::string &prm_filename)
 *   {
 *     /* First declare the parameters... */
 *     ParameterHandler prm;
 * 
 *     prm.declare_entry("Epsilon",
 *                       "0.005",
 *                       Patterns::Double(0),
 *                       "Diffusion parameter");
 * 
 *     prm.declare_entry("Fe degree",
 *                       "1",
 *                       Patterns::Integer(1),
 *                       "Finite Element degree");
 *     prm.declare_entry("Smoother type",
 *                       "block SOR",
 *                       Patterns::Selection("SOR|Jacobi|block SOR|block Jacobi"),
 *                       "Select smoother: SOR|Jacobi|block SOR|block Jacobi");
 *     prm.declare_entry("Smoothing steps",
 *                       "2",
 *                       Patterns::Integer(1),
 *                       "Number of smoothing steps");
 *     prm.declare_entry(
 *       "DoF renumbering",
 *       "downstream",
 *       Patterns::Selection("none|downstream|upstream|random"),
 *       "Select DoF renumbering: none|downstream|upstream|random");
 *     prm.declare_entry("With streamline diffusion",
 *                       "true",
 *                       Patterns::Bool(),
 *                       "Enable streamline diffusion stabilization: true|false");
 *     prm.declare_entry("Output",
 *                       "true",
 *                       Patterns::Bool(),
 *                       "Generate graphical output: true|false");
 * 
 *     /* ...and then try to read their values from the input file: */
 *     if (prm_filename.empty())
 *       {
 *         prm.print_parameters(std::cout, ParameterHandler::Text);
 *         AssertThrow(
 *           false, ExcMessage("Please pass a .prm file as the first argument!"));
 *       }
 * 
 *     prm.parse_input(prm_filename);
 * 
 *     epsilon         = prm.get_double("Epsilon");
 *     fe_degree       = prm.get_integer("Fe degree");
 *     smoother_type   = prm.get("Smoother type");
 *     smoothing_steps = prm.get_integer("Smoothing steps");
 * 
 *     const std::string renumbering = prm.get("DoF renumbering");
 *     if (renumbering == "none")
 *       dof_renumbering = DoFRenumberingStrategy::none;
 *     else if (renumbering == "downstream")
 *       dof_renumbering = DoFRenumberingStrategy::downstream;
 *     else if (renumbering == "upstream")
 *       dof_renumbering = DoFRenumberingStrategy::upstream;
 *     else if (renumbering == "random")
 *       dof_renumbering = DoFRenumberingStrategy::random;
 *     else
 *       AssertThrow(false,
 *                   ExcMessage("The <DoF renumbering> parameter has "
 *                              "an invalid value."));
 * 
 *     with_streamline_diffusion = prm.get_bool("With streamline diffusion");
 *     output                    = prm.get_bool("Output");
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Cellpermutations"></a> 
 * <h3>Cell permutations</h3>
 *   

 * 
 * The ordering in which cells and degrees of freedom are traversed
 * will play a role in the speed of convergence for multiplicative
 * methods. Here we define functions which return a specific ordering
 * of cells to be used by the block smoothers.
 *   

 * 
 * For each type of cell ordering, we define a function for the
 * active mesh and one for a level mesh (i.e., for the cells at one
 * level of a multigrid hierarchy). While the only reordering
 * necessary for solving the system will be on the level meshes, we
 * include the active reordering for visualization purposes in
 * output_results().
 *   

 * 
 * For the two downstream ordering functions, we first create an
 * array with all of the relevant cells that we then sort in
 * downstream direction using a "comparator" object. The output of
 * the functions is then simply an array of the indices of the cells
 * in the just computed order.
 * 
 * @code
 *   template <int dim>
 *   std::vector<unsigned int>
 *   create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                                   const Tensor<1, dim>   direction,
 *                                   const unsigned int     level)
 *   {
 *     std::vector<typename DoFHandler<dim>::level_cell_iterator> ordered_cells;
 *     ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level));
 *     for (const auto &cell : dof_handler.cell_iterators_on_level(level))
 *       ordered_cells.push_back(cell);
 * 
 *     const DoFRenumbering::
 *       CompareDownstream<typename DoFHandler<dim>::level_cell_iterator, dim>
 *         comparator(direction);
 *     std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);
 * 
 *     std::vector<unsigned> ordered_indices;
 *     ordered_indices.reserve(dof_handler.get_triangulation().n_cells(level));
 * 
 *     for (const auto &cell : ordered_cells)
 *       ordered_indices.push_back(cell->index());
 * 
 *     return ordered_indices;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   std::vector<unsigned int>
 *   create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                                   const Tensor<1, dim>   direction)
 *   {
 *     std::vector<typename DoFHandler<dim>::active_cell_iterator> ordered_cells;
 *     ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells());
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       ordered_cells.push_back(cell);
 * 
 *     const DoFRenumbering::
 *       CompareDownstream<typename DoFHandler<dim>::active_cell_iterator, dim>
 *         comparator(direction);
 *     std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);
 * 
 *     std::vector<unsigned int> ordered_indices;
 *     ordered_indices.reserve(dof_handler.get_triangulation().n_active_cells());
 * 
 *     for (const auto &cell : ordered_cells)
 *       ordered_indices.push_back(cell->index());
 * 
 *     return ordered_indices;
 *   }
 * 
 * 
 * @endcode
 * 
 * The functions that produce a random ordering are similar in
 * spirit in that they first put information about all cells into an
 * array. But then, instead of sorting them, they shuffle the
 * elements randomly using the facilities C++ offers to generate
 * random numbers. The way this is done is by iterating over all
 * elements of the array, drawing a random number for another
 * element before that, and then exchanging these elements. The
 * result is a random shuffle of the elements of the array.
 * 
 * @code
 *   template <int dim>
 *   std::vector<unsigned int>
 *   create_random_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                               const unsigned int     level)
 *   {
 *     std::vector<unsigned int> ordered_cells;
 *     ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level));
 *     for (const auto &cell : dof_handler.cell_iterators_on_level(level))
 *       ordered_cells.push_back(cell->index());
 * 
 *     std::mt19937 random_number_generator;
 *     std::shuffle(ordered_cells.begin(),
 *                  ordered_cells.end(),
 *                  random_number_generator);
 * 
 *     return ordered_cells;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   std::vector<unsigned int>
 *   create_random_cell_ordering(const DoFHandler<dim> &dof_handler)
 *   {
 *     std::vector<unsigned int> ordered_cells;
 *     ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells());
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       ordered_cells.push_back(cell->index());
 * 
 *     std::mt19937 random_number_generator;
 *     std::shuffle(ordered_cells.begin(),
 *                  ordered_cells.end(),
 *                  random_number_generator);
 * 
 *     return ordered_cells;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideandboundaryvalues"></a> 
 * <h3>Right-hand side and boundary values</h3>
 * 

 * 
 * The problem solved in this tutorial is an adaptation of Ex. 3.1.3 found
 * on pg. 118 of <a
 * href="https://global.oup.com/academic/product/finite-elements-and-fast-iterative-solvers-9780199678808">
 * Finite Elements and Fast Iterative Solvers: with Applications in
 * Incompressible Fluid Dynamics by Elman, Silvester, and Wathen</a>. The
 * main difference being that we add a hole in the center of our domain with
 * zero Dirichlet boundary conditions.
 *   

 * 
 * For a complete description, we need classes that implement the
 * zero right-hand side first (we could of course have just used
 * Functions::ZeroFunction):
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual void value_list(const std::vector<Point<dim>> &points,
 *                             std::vector<double> &          values,
 *                             const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> &,
 *                                    const unsigned int component) const
 *   {
 *     Assert(component == 0, ExcIndexRange(component, 0, 1));
 *     (void)component;
 * 
 *     return 0.0;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                       std::vector<double> &          values,
 *                                       const unsigned int component) const
 *   {
 *     Assert(values.size() == points.size(),
 *            ExcDimensionMismatch(values.size(), points.size()));
 * 
 *     for (unsigned int i = 0; i < points.size(); ++i)
 *       values[i] = RightHandSide<dim>::value(points[i], component);
 *   }
 * 
 * 
 * @endcode
 * 
 * We also have Dirichlet boundary conditions. On a connected portion of the
 * outer, square boundary we set the value to 1, and we set the value to 0
 * everywhere else (including the inner, circular boundary):
 * 
 * @code
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual void value_list(const std::vector<Point<dim>> &points,
 *                             std::vector<double> &          values,
 *                             const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> & p,
 *                                     const unsigned int component) const
 *   {
 *     Assert(component == 0, ExcIndexRange(component, 0, 1));
 *     (void)component;
 * 
 * @endcode
 * 
 * Set boundary to 1 if $x=1$, or if $x>0.5$ and $y=-1$.
 * 
 * @code
 *     if (std::fabs(p[0] - 1) < 1e-8 ||
 *         (std::fabs(p[1] + 1) < 1e-8 && p[0] >= 0.5))
 *       {
 *         return 1.0;
 *       }
 *     else
 *       {
 *         return 0.0;
 *       }
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                        std::vector<double> &          values,
 *                                        const unsigned int component) const
 *   {
 *     Assert(values.size() == points.size(),
 *            ExcDimensionMismatch(values.size(), points.size()));
 * 
 *     for (unsigned int i = 0; i < points.size(); ++i)
 *       values[i] = BoundaryValues<dim>::value(points[i], component);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Streamlinediffusionimplementation"></a> 
 * <h3>Streamline diffusion implementation</h3>
 * 

 * 
 * The streamline diffusion method has a stabilization constant that
 * we need to be able to compute. The choice of how this parameter
 * is computed is taken from <a
 * href="https://link.springer.com/chapter/10.1007/978-3-540-34288-5_27">On
 * Discontinuity-Capturing Methods for Convection-Diffusion
 * Equations by Volker John and Petr Knobloch</a>.
 * 
 * @code
 *   template <int dim>
 *   double compute_stabilization_delta(const double         hk,
 *                                      const double         eps,
 *                                      const Tensor<1, dim> dir,
 *                                      const double         pk)
 *   {
 *     const double Peclet = dir.norm() * hk / (2.0 * eps * pk);
 *     const double coth =
 *       (1.0 + std::exp(-2.0 * Peclet)) / (1.0 - std::exp(-2.0 * Peclet));
 * 
 *     return hk / (2.0 * dir.norm() * pk) * (coth - 1.0 / Peclet);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProlemcodeclass"></a> 
 * <h3><code>AdvectionProlem</code> class</h3>
 * 

 * 
 * This is the main class of the program, and should look very similar to
 * step-16. The major difference is that, since we are defining our multigrid
 * smoother at runtime, we choose to define a function `create_smoother()` and
 * a class object `mg_smoother` which is a `std::unique_ptr` to a smoother
 * that is derived from MGSmoother. Note that for smoothers derived from
 * RelaxationBlock, we must include a `smoother_data` object for each level.
 * This will contain information about the cell ordering and the method of
 * inverting cell matrices.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class AdvectionProblem
 *   {
 *   public:
 *     AdvectionProblem(const Settings &settings);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 * 
 *     template <class IteratorType>
 *     void assemble_cell(const IteratorType &cell,
 *                        ScratchData<dim> &  scratch_data,
 *                        CopyData &          copy_data);
 *     void assemble_system_and_multigrid();
 * 
 *     void setup_smoother();
 * 
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     const FE_Q<dim>     fe;
 *     const MappingQ<dim> mapping;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     MGLevelObject<SparsityPattern> mg_sparsity_patterns;
 *     MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns;
 * 
 *     MGLevelObject<SparseMatrix<double>> mg_matrices;
 *     MGLevelObject<SparseMatrix<double>> mg_interface_in;
 *     MGLevelObject<SparseMatrix<double>> mg_interface_out;
 * 
 *     mg::Matrix<Vector<double>> mg_matrix;
 *     mg::Matrix<Vector<double>> mg_interface_matrix_in;
 *     mg::Matrix<Vector<double>> mg_interface_matrix_out;
 * 
 *     std::unique_ptr<MGSmoother<Vector<double>>> mg_smoother;
 * 
 *     using SmootherType =
 *       RelaxationBlock<SparseMatrix<double>, double, Vector<double>>;
 *     using SmootherAdditionalDataType = SmootherType::AdditionalData;
 *     MGLevelObject<SmootherAdditionalDataType> smoother_data;
 * 
 *     MGConstrainedDoFs mg_constrained_dofs;
 * 
 *     Tensor<1, dim> advection_direction;
 * 
 *     const Settings settings;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   AdvectionProblem<dim>::AdvectionProblem(const Settings &settings)
 *     : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     , dof_handler(triangulation)
 *     , fe(settings.fe_degree)
 *     , mapping(settings.fe_degree)
 *     , settings(settings)
 *   {
 *     advection_direction[0] = -std::sin(numbers::PI / 6.0);
 *     if (dim >= 2)
 *       advection_direction[1] = std::cos(numbers::PI / 6.0);
 *     if (dim >= 3)
 *       AssertThrow(false, ExcNotImplemented());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemsetup_systemcode"></a> 
 * <h4><code>AdvectionProblem::setup_system()</code></h4>
 * 

 * 
 * Here we first set up the DoFHandler, AffineConstraints, and
 * SparsityPattern objects for both active and multigrid level meshes.
 *   

 * 
 * We could renumber the active DoFs with the DoFRenumbering class,
 * but the smoothers only act on multigrid levels and as such, this
 * would not matter for the computations. Instead, we will renumber the
 * DoFs on each multigrid level below.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::setup_system()
 *   {
 *     const unsigned int n_levels = triangulation.n_levels();
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 0, BoundaryValues<dim>(), constraints);
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 1, BoundaryValues<dim>(), constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ false);
 * 
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     dof_handler.distribute_mg_dofs();
 * 
 * @endcode
 * 
 * Having enumerated the global degrees of freedom as well as (in
 * the last line above) the level degrees of freedom, let us
 * renumber the level degrees of freedom to get a better smoother
 * as explained in the introduction.  The first block below
 * renumbers DoFs on each level in downstream or upstream
 * direction if needed. This is only necessary for point smoothers
 * (SOR and Jacobi) as the block smoothers operate on cells (see
 * `create_smoother()`). The blocks below then also implement
 * random numbering.
 * 
 * @code
 *     if (settings.smoother_type == "SOR" || settings.smoother_type == "Jacobi")
 *       {
 *         if (settings.dof_renumbering ==
 *               Settings::DoFRenumberingStrategy::downstream ||
 *             settings.dof_renumbering ==
 *               Settings::DoFRenumberingStrategy::upstream)
 *           {
 *             const Tensor<1, dim> direction =
 *               (settings.dof_renumbering ==
 *                    Settings::DoFRenumberingStrategy::upstream ?
 *                  -1.0 :
 *                  1.0) *
 *               advection_direction;
 * 
 *             for (unsigned int level = 0; level < n_levels; ++level)
 *               DoFRenumbering::downstream(dof_handler,
 *                                          level,
 *                                          direction,
 *                                          /*dof_wise_renumbering = */ true);
 *           }
 *         else if (settings.dof_renumbering ==
 *                  Settings::DoFRenumberingStrategy::random)
 *           {
 *             for (unsigned int level = 0; level < n_levels; ++level)
 *               DoFRenumbering::random(dof_handler, level);
 *           }
 *         else
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 * @endcode
 * 
 * The rest of the function just sets up data structures. The last
 * lines of the code below is unlike the other GMG tutorials, as
 * it sets up both the interface in and out matrices. We need this
 * since our problem is non-symmetric.
 * 
 * @code
 *     mg_constrained_dofs.clear();
 *     mg_constrained_dofs.initialize(dof_handler);
 * 
 *     mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, {0, 1});
 * 
 *     mg_matrices.resize(0, n_levels - 1);
 *     mg_matrices.clear_elements();
 *     mg_interface_in.resize(0, n_levels - 1);
 *     mg_interface_in.clear_elements();
 *     mg_interface_out.resize(0, n_levels - 1);
 *     mg_interface_out.clear_elements();
 *     mg_sparsity_patterns.resize(0, n_levels - 1);
 *     mg_interface_sparsity_patterns.resize(0, n_levels - 1);
 * 
 *     for (unsigned int level = 0; level < n_levels; ++level)
 *       {
 *         {
 *           DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
 *                                      dof_handler.n_dofs(level));
 *           MGTools::make_sparsity_pattern(dof_handler, dsp, level);
 *           mg_sparsity_patterns[level].copy_from(dsp);
 *           mg_matrices[level].reinit(mg_sparsity_patterns[level]);
 *         }
 *         {
 *           DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
 *                                      dof_handler.n_dofs(level));
 *           MGTools::make_interface_sparsity_pattern(dof_handler,
 *                                                    mg_constrained_dofs,
 *                                                    dsp,
 *                                                    level);
 *           mg_interface_sparsity_patterns[level].copy_from(dsp);
 * 
 *           mg_interface_in[level].reinit(mg_interface_sparsity_patterns[level]);
 *           mg_interface_out[level].reinit(mg_interface_sparsity_patterns[level]);
 *         }
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemassemble_cellcode"></a> 
 * <h4><code>AdvectionProblem::assemble_cell()</code></h4>
 * 

 * 
 * Here we define the assembly of the linear system on each cell to
 * be used by the mesh_loop() function below. This one function
 * assembles the cell matrix for either an active or a level cell
 * (whatever it is passed as its first argument), and only assembles
 * a right-hand side if called with an active cell.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   template <class IteratorType>
 *   void AdvectionProblem<dim>::assemble_cell(const IteratorType &cell,
 *                                             ScratchData<dim> &  scratch_data,
 *                                             CopyData &          copy_data)
 *   {
 *     copy_data.level = cell->level();
 * 
 *     const unsigned int dofs_per_cell =
 *       scratch_data.fe_values.get_fe().n_dofs_per_cell();
 *     copy_data.dofs_per_cell = dofs_per_cell;
 *     copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 * 
 *     const unsigned int n_q_points =
 *       scratch_data.fe_values.get_quadrature().size();
 * 
 *     if (cell->is_level_cell() == false)
 *       copy_data.cell_rhs.reinit(dofs_per_cell);
 * 
 *     copy_data.local_dof_indices.resize(dofs_per_cell);
 *     cell->get_active_or_mg_dof_indices(copy_data.local_dof_indices);
 * 
 *     scratch_data.fe_values.reinit(cell);
 * 
 *     RightHandSide<dim>  right_hand_side;
 *     std::vector<double> rhs_values(n_q_points);
 * 
 *     right_hand_side.value_list(scratch_data.fe_values.get_quadrature_points(),
 *                                rhs_values);
 * 
 * @endcode
 * 
 * If we are using streamline diffusion we must add its contribution
 * to both the cell matrix and the cell right-hand side. If we are not
 * using streamline diffusion, setting $\delta=0$ negates this contribution
 * below and we are left with the standard, Galerkin finite element
 * assembly.
 * 
 * @code
 *     const double delta = (settings.with_streamline_diffusion ?
 *                             compute_stabilization_delta(cell->diameter(),
 *                                                         settings.epsilon,
 *                                                         advection_direction,
 *                                                         settings.fe_degree) :
 *                             0.0);
 * 
 *     for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             {
 * @endcode
 * 
 * The assembly of the local matrix has two parts. First
 * the Galerkin contribution:
 * 
 * @code
 *               copy_data.cell_matrix(i, j) +=
 *                 (settings.epsilon *
 *                  scratch_data.fe_values.shape_grad(i, q_point) *
 *                  scratch_data.fe_values.shape_grad(j, q_point) *
 *                  scratch_data.fe_values.JxW(q_point)) +
 *                 (scratch_data.fe_values.shape_value(i, q_point) *
 *                  (advection_direction *
 *                   scratch_data.fe_values.shape_grad(j, q_point)) *
 *                  scratch_data.fe_values.JxW(q_point))
 * @endcode
 * 
 * and then the streamline diffusion contribution:
 * 
 * @code
 *                 + delta *
 *                     (advection_direction *
 *                      scratch_data.fe_values.shape_grad(j, q_point)) *
 *                     (advection_direction *
 *                      scratch_data.fe_values.shape_grad(i, q_point)) *
 *                     scratch_data.fe_values.JxW(q_point) -
 *                 delta * settings.epsilon *
 *                   trace(scratch_data.fe_values.shape_hessian(j, q_point)) *
 *                   (advection_direction *
 *                    scratch_data.fe_values.shape_grad(i, q_point)) *
 *                   scratch_data.fe_values.JxW(q_point);
 *             }
 *           if (cell->is_level_cell() == false)
 *             {
 * @endcode
 * 
 * The same applies to the right hand side. First the
 * Galerkin contribution:
 * 
 * @code
 *               copy_data.cell_rhs(i) +=
 *                 scratch_data.fe_values.shape_value(i, q_point) *
 *                   rhs_values[q_point] * scratch_data.fe_values.JxW(q_point)
 * @endcode
 * 
 * and then the streamline diffusion contribution:
 * 
 * @code
 *                 + delta * rhs_values[q_point] * advection_direction *
 *                     scratch_data.fe_values.shape_grad(i, q_point) *
 *                     scratch_data.fe_values.JxW(q_point);
 *             }
 *         }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemassemble_system_and_multigridcode"></a> 
 * <h4><code>AdvectionProblem::assemble_system_and_multigrid()</code></h4>
 * 

 * 
 * Here we employ MeshWorker::mesh_loop() to go over cells and assemble the
 * system_matrix, system_rhs, and all mg_matrices for us.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::assemble_system_and_multigrid()
 *   {
 *     const auto cell_worker_active =
 *       [&](const decltype(dof_handler.begin_active()) &cell,
 *           ScratchData<dim> &                          scratch_data,
 *           CopyData &                                  copy_data) {
 *         this->assemble_cell(cell, scratch_data, copy_data);
 *       };
 * 
 *     const auto copier_active = [&](const CopyData &copy_data) {
 *       constraints.distribute_local_to_global(copy_data.cell_matrix,
 *                                              copy_data.cell_rhs,
 *                                              copy_data.local_dof_indices,
 *                                              system_matrix,
 *                                              system_rhs);
 *     };
 * 
 * 
 *     MeshWorker::mesh_loop(dof_handler.begin_active(),
 *                           dof_handler.end(),
 *                           cell_worker_active,
 *                           copier_active,
 *                           ScratchData<dim>(fe, fe.degree + 1),
 *                           CopyData(),
 *                           MeshWorker::assemble_own_cells);
 * 
 * @endcode
 * 
 * Unlike the constraints for the active level, we choose to create
 * constraint objects for each multigrid level local to this function
 * since they are never needed elsewhere in the program.
 * 
 * @code
 *     std::vector<AffineConstraints<double>> boundary_constraints(
 *       triangulation.n_global_levels());
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       {
 *         IndexSet locally_owned_level_dof_indices;
 *         DoFTools::extract_locally_relevant_level_dofs(
 *           dof_handler, level, locally_owned_level_dof_indices);
 *         boundary_constraints[level].reinit(locally_owned_level_dof_indices);
 *         boundary_constraints[level].add_lines(
 *           mg_constrained_dofs.get_refinement_edge_indices(level));
 *         boundary_constraints[level].add_lines(
 *           mg_constrained_dofs.get_boundary_indices(level));
 *         boundary_constraints[level].close();
 *       }
 * 
 *     const auto cell_worker_mg =
 *       [&](const decltype(dof_handler.begin_mg()) &cell,
 *           ScratchData<dim> &                      scratch_data,
 *           CopyData &                              copy_data) {
 *         this->assemble_cell(cell, scratch_data, copy_data);
 *       };
 * 
 *     const auto copier_mg = [&](const CopyData &copy_data) {
 *       boundary_constraints[copy_data.level].distribute_local_to_global(
 *         copy_data.cell_matrix,
 *         copy_data.local_dof_indices,
 *         mg_matrices[copy_data.level]);
 * 
 * @endcode
 * 
 * If $(i,j)$ is an `interface_out` dof pair, then $(j,i)$ is an
 * `interface_in` dof pair. Note: For `interface_in`, we load
 * the transpose of the interface entries, i.e., the entry for
 * dof pair $(j,i)$ is stored in `interface_in(i,j)`. This is an
 * optimization for the symmetric case which allows only one
 * matrix to be used when setting the edge_matrices in
 * solve(). Here, however, since our problem is non-symmetric,
 * we must store both `interface_in` and `interface_out`
 * matrices.
 * 
 * @code
 *       for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
 *         for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j)
 *           if (mg_constrained_dofs.is_interface_matrix_entry(
 *                 copy_data.level,
 *                 copy_data.local_dof_indices[i],
 *                 copy_data.local_dof_indices[j]))
 *             {
 *               mg_interface_out[copy_data.level].add(
 *                 copy_data.local_dof_indices[i],
 *                 copy_data.local_dof_indices[j],
 *                 copy_data.cell_matrix(i, j));
 *               mg_interface_in[copy_data.level].add(
 *                 copy_data.local_dof_indices[i],
 *                 copy_data.local_dof_indices[j],
 *                 copy_data.cell_matrix(j, i));
 *             }
 *     };
 * 
 *     MeshWorker::mesh_loop(dof_handler.begin_mg(),
 *                           dof_handler.end_mg(),
 *                           cell_worker_mg,
 *                           copier_mg,
 *                           ScratchData<dim>(fe, fe.degree + 1),
 *                           CopyData(),
 *                           MeshWorker::assemble_own_cells);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemsetup_smoothercode"></a> 
 * <h4><code>AdvectionProblem::setup_smoother()</code></h4>
 * 

 * 
 * Next, we set up the smoother based on the settings in the `.prm` file. The
 * two options that are of significance is the number of pre- and
 * post-smoothing steps on each level of the multigrid v-cycle and the
 * relaxation parameter.
 * 

 * 
 * Since multiplicative methods tend to be more powerful than additive method,
 * fewer smoothing steps are required to see convergence independent of mesh
 * size. The same holds for block smoothers over point smoothers. This is
 * reflected in the choice for the number of smoothing steps for each type of
 * smoother below.
 * 

 * 
 * The relaxation parameter for point smoothers is chosen based on trial and
 * error, and reflects values necessary to keep the iteration counts in
 * the GMRES solve constant (or as close as possible) as we refine the mesh.
 * The two values given for both "Jacobi" and "SOR" in the `.prm` files are
 * for degree 1 and degree 3 finite elements. If the user wants to change to
 * another degree, they may need to adjust these numbers. For block smoothers,
 * this parameter has a more straightforward interpretation, namely that for
 * additive methods in 2D, a DoF can have a repeated contribution from up to 4
 * cells, therefore we must relax these methods by 0.25 to compensate. This is
 * not an issue for multiplicative methods as each cell's inverse application
 * carries new information to all its DoFs.
 * 

 * 
 * Finally, as mentioned above, the point smoothers only operate on DoFs, and
 * the block smoothers on cells, so only the block smoothers need to be given
 * information regarding cell orderings. DoF ordering for point smoothers has
 * already been taken care of in `setup_system()`.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::setup_smoother()
 *   {
 *     if (settings.smoother_type == "SOR")
 *       {
 *         using Smoother = PreconditionSOR<SparseMatrix<double>>;
 * 
 *         auto smoother =
 *           std::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
 *                                                   Smoother,
 *                                                   Vector<double>>>();
 *         smoother->initialize(mg_matrices,
 *                              Smoother::AdditionalData(fe.degree == 1 ? 1.0 :
 *                                                                        0.62));
 *         smoother->set_steps(settings.smoothing_steps);
 *         mg_smoother = std::move(smoother);
 *       }
 *     else if (settings.smoother_type == "Jacobi")
 *       {
 *         using Smoother = PreconditionJacobi<SparseMatrix<double>>;
 *         auto smoother =
 *           std::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
 *                                                   Smoother,
 *                                                   Vector<double>>>();
 *         smoother->initialize(mg_matrices,
 *                              Smoother::AdditionalData(fe.degree == 1 ? 0.6667 :
 *                                                                        0.47));
 *         smoother->set_steps(settings.smoothing_steps);
 *         mg_smoother = std::move(smoother);
 *       }
 *     else if (settings.smoother_type == "block SOR" ||
 *              settings.smoother_type == "block Jacobi")
 *       {
 *         smoother_data.resize(0, triangulation.n_levels() - 1);
 * 
 *         for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
 *           {
 *             DoFTools::make_cell_patches(smoother_data[level].block_list,
 *                                         dof_handler,
 *                                         level);
 * 
 *             smoother_data[level].relaxation =
 *               (settings.smoother_type == "block SOR" ? 1.0 : 0.25);
 *             smoother_data[level].inversion = PreconditionBlockBase<double>::svd;
 * 
 *             std::vector<unsigned int> ordered_indices;
 *             switch (settings.dof_renumbering)
 *               {
 *                 case Settings::DoFRenumberingStrategy::downstream:
 *                   ordered_indices =
 *                     create_downstream_cell_ordering(dof_handler,
 *                                                     advection_direction,
 *                                                     level);
 *                   break;
 * 
 *                 case Settings::DoFRenumberingStrategy::upstream:
 *                   ordered_indices =
 *                     create_downstream_cell_ordering(dof_handler,
 *                                                     -1.0 * advection_direction,
 *                                                     level);
 *                   break;
 * 
 *                 case Settings::DoFRenumberingStrategy::random:
 *                   ordered_indices =
 *                     create_random_cell_ordering(dof_handler, level);
 *                   break;
 * 
 *                 case Settings::DoFRenumberingStrategy::none:
 *                   break;
 * 
 *                 default:
 *                   AssertThrow(false, ExcNotImplemented());
 *                   break;
 *               }
 * 
 *             smoother_data[level].order =
 *               std::vector<std::vector<unsigned int>>(1, ordered_indices);
 *           }
 * 
 *         if (settings.smoother_type == "block SOR")
 *           {
 *             auto smoother = std::make_unique<MGSmootherPrecondition<
 *               SparseMatrix<double>,
 *               RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>,
 *               Vector<double>>>();
 *             smoother->initialize(mg_matrices, smoother_data);
 *             smoother->set_steps(settings.smoothing_steps);
 *             mg_smoother = std::move(smoother);
 *           }
 *         else if (settings.smoother_type == "block Jacobi")
 *           {
 *             auto smoother = std::make_unique<
 *               MGSmootherPrecondition<SparseMatrix<double>,
 *                                      RelaxationBlockJacobi<SparseMatrix<double>,
 *                                                            double,
 *                                                            Vector<double>>,
 *                                      Vector<double>>>();
 *             smoother->initialize(mg_matrices, smoother_data);
 *             smoother->set_steps(settings.smoothing_steps);
 *             mg_smoother = std::move(smoother);
 *           }
 *       }
 *     else
 *       AssertThrow(false, ExcNotImplemented());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemsolvecode"></a> 
 * <h4><code>AdvectionProblem::solve()</code></h4>
 * 

 * 
 * Before we can solve the system, we must first set up the multigrid
 * preconditioner. This requires the setup of the transfer between levels,
 * the coarse matrix solver, and the smoother. This setup follows almost
 * identically to Step-16, the main difference being the various smoothers
 * defined above and the fact that we need different interface edge matrices
 * for in and out since our problem is non-symmetric. (In reality, for this
 * tutorial these interface matrices are empty since we are only using global
 * refinement, and thus have no refinement edges. However, we have still
 * included both here since if one made the simple switch to an adaptively
 * refined method, the program would still run correctly.)
 * 

 * 
 * The last thing to note is that since our problem is non-symmetric, we must
 * use an appropriate Krylov subspace method. We choose here to
 * use GMRES since it offers the guarantee of residual reduction in each
 * iteration. The major disavantage of GMRES is that, for each iteration,
 * the number of stored temporary vectors increases by one, and one also needs
 * to compute a scalar product with all previously stored vectors. This is
 * rather expensive. This requirement is relaxed by using the restarted GMRES
 * method which puts a cap on the number of vectors we are required to store
 * at any one time (here we restart after 50 temporary vectors, or 48
 * iterations). This then has the disadvantage that we lose information we
 * have gathered throughout the iteration and therefore we could see slower
 * convergence. As a consequence, where to restart is a question of balancing
 * memory consumption, CPU effort, and convergence speed.
 * However, the goal of this tutorial is to have very low
 * iteration counts by using a powerful GMG preconditioner, so we have picked
 * the restart length such that all of the results shown below converge prior
 * to restart happening, and thus we have a standard GMRES method. If the user
 * is interested, another suitable method offered in deal.II would be
 * BiCGStab.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::solve()
 *   {
 *     const unsigned int max_iters       = 200;
 *     const double       solve_tolerance = 1e-8 * system_rhs.l2_norm();
 *     SolverControl      solver_control(max_iters, solve_tolerance, true, true);
 *     solver_control.enable_history_data();
 * 
 *     using Transfer = MGTransferPrebuilt<Vector<double>>;
 *     Transfer mg_transfer(mg_constrained_dofs);
 *     mg_transfer.build(dof_handler);
 * 
 *     FullMatrix<double> coarse_matrix;
 *     coarse_matrix.copy_from(mg_matrices[0]);
 *     MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
 *     coarse_grid_solver.initialize(coarse_matrix);
 * 
 *     setup_smoother();
 * 
 *     mg_matrix.initialize(mg_matrices);
 *     mg_interface_matrix_in.initialize(mg_interface_in);
 *     mg_interface_matrix_out.initialize(mg_interface_out);
 * 
 *     Multigrid<Vector<double>> mg(
 *       mg_matrix, coarse_grid_solver, mg_transfer, *mg_smoother, *mg_smoother);
 *     mg.set_edge_matrices(mg_interface_matrix_out, mg_interface_matrix_in);
 * 
 *     PreconditionMG<dim, Vector<double>, Transfer> preconditioner(dof_handler,
 *                                                                  mg,
 *                                                                  mg_transfer);
 * 
 *     std::cout << "     Solving with GMRES to tol " << solve_tolerance << "..."
 *               << std::endl;
 *     SolverGMRES<Vector<double>> solver(
 *       solver_control, SolverGMRES<Vector<double>>::AdditionalData(50, true));
 * 
 *     Timer time;
 *     time.start();
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *     time.stop();
 * 
 *     std::cout << "          converged in " << solver_control.last_step()
 *               << " iterations"
 *               << " in " << time.last_wall_time() << " seconds " << std::endl;
 * 
 *     constraints.distribute(solution);
 * 
 *     mg_smoother.release();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemoutput_resultscode"></a> 
 * <h4><code>AdvectionProblem::output_results()</code></h4>
 * 

 * 
 * The final function of interest generates graphical output.
 * Here we output the solution and cell ordering in a .vtu format.
 * 

 * 
 * At the top of the function, we generate an index for each cell to
 * visualize the ordering used by the smoothers. Note that we do
 * this only for the active cells instead of the levels, where the
 * smoothers are actually used. For the point smoothers we renumber
 * DoFs instead of cells, so this is only an approximation of what
 * happens in reality. Finally, the random ordering is not the
 * random ordering we actually use (see `create_smoother()` for that).
 *   

 * 
 * The (integer) ordering of cells is then copied into a (floating
 * point) vector for graphical output.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     const unsigned int n_active_cells = triangulation.n_active_cells();
 *     Vector<double>     cell_indices(n_active_cells);
 *     {
 *       std::vector<unsigned int> ordered_indices;
 *       switch (settings.dof_renumbering)
 *         {
 *           case Settings::DoFRenumberingStrategy::downstream:
 *             ordered_indices =
 *               create_downstream_cell_ordering(dof_handler, advection_direction);
 *             break;
 * 
 *           case Settings::DoFRenumberingStrategy::upstream:
 *             ordered_indices =
 *               create_downstream_cell_ordering(dof_handler,
 *                                               -1.0 * advection_direction);
 *             break;
 * 
 *           case Settings::DoFRenumberingStrategy::random:
 *             ordered_indices = create_random_cell_ordering(dof_handler);
 *             break;
 * 
 *           case Settings::DoFRenumberingStrategy::none:
 *             ordered_indices.resize(n_active_cells);
 *             for (unsigned int i = 0; i < n_active_cells; ++i)
 *               ordered_indices[i] = i;
 *             break;
 * 
 *           default:
 *             AssertThrow(false, ExcNotImplemented());
 *             break;
 *         }
 * 
 *       for (unsigned int i = 0; i < n_active_cells; ++i)
 *         cell_indices(ordered_indices[i]) = static_cast<double>(i);
 *     }
 * 
 * @endcode
 * 
 * The remainder of the function is then straightforward, given
 * previous tutorial programs:
 * 
 * @code
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.add_data_vector(cell_indices, "cell_index");
 *     data_out.build_patches();
 * 
 *     const std::string filename =
 *       "solution-" + Utilities::int_to_string(cycle) + ".vtu";
 *     std::ofstream output(filename.c_str());
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeAdvectionProblemruncode"></a> 
 * <h4><code>AdvectionProblem::run()</code></h4>
 * 

 * 
 * As in most tutorials, this function creates/refines the mesh and calls
 * the various functions defined above to set up, assemble, solve, and output
 * the results.
 * 

 * 
 * In cycle zero, we generate the mesh for the on the square
 * <code>[-1,1]^dim</code> with a hole of radius 3/10 units centered
 * at the origin. For objects with `manifold_id` equal to one
 * (namely, the faces adjacent to the hole), we assign a spherical
 * manifold.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < (settings.fe_degree == 1 ? 7 : 5);
 *          ++cycle)
 *       {
 *         std::cout << "  Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,
 *                                                             0.3,
 *                                                             1.0);
 * 
 *             const SphericalManifold<dim> manifold_description(Point<dim>(0, 0));
 *             triangulation.set_manifold(1, manifold_description);
 *           }
 * 
 *         triangulation.refine_global();
 * 
 *         setup_system();
 * 
 *         std::cout << "     Number of active cells:       "
 *                   << triangulation.n_active_cells() << " ("
 *                   << triangulation.n_levels() << " levels)" << std::endl;
 *         std::cout << "     Number of degrees of freedom: "
 *                   << dof_handler.n_dofs() << std::endl;
 * 
 *         assemble_system_and_multigrid();
 * 
 *         solve();
 * 
 *         if (settings.output)
 *           output_results(cycle);
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * } // namespace Step63
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * Finally, the main function is like most tutorials. The only
 * interesting bit is that we require the user to pass a `.prm` file
 * as a sole command line argument. If no parameter file is given, the
 * program will output the contents of a sample parameter file with
 * all default values to the screen that the user can then copy and
 * paste into their own `.prm` file.
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       Step63::Settings settings;
 *       settings.get_parameters((argc > 1) ? (argv[1]) : "");
 * 
 *       Step63::AdvectionProblem<2> advection_problem_2d(settings);
 *       advection_problem_2d.run();
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


<a name="GMRESIterationNumbers"></a><h3> GMRES Iteration Numbers </h3>


The major advantage for GMG is that it is an $\mathcal{O}(n)$ method,
that is, the complexity of the problem increases linearly with the
problem size. To show then that the linear solver presented in this
tutorial is in fact $\mathcal{O}(n)$, all one needs to do is show that
the iteration counts for the GMRES solve stay roughly constant as we
refine the mesh.

Each of the following tables gives the GMRES iteration counts to reduce the
initial residual by a factor of $10^8$. We selected a sufficient number of smoothing steps
(based on the method) to get iteration numbers independent of mesh size. As
can be seen from the tables below, the method is indeed $\mathcal{O}(n)$.

<a name="DoFCellRenumbering"></a><h4> DoF/Cell Renumbering </h4>


The point-wise smoothers ("Jacobi" and "SOR") get applied in the order the
DoFs are numbered on each level. We can influence this using the
DoFRenumbering namespace. The block smoothers are applied based on the
ordering we set in `setup_smoother()`. We can visualize this numbering. The
following pictures show the cell numbering of the active cells in downstream,
random, and upstream numbering (left to right):

<img src="https://www.dealii.org/images/steps/developer/step-63-cell-order.png" alt="">

Let us start with the additive smoothers. The following table shows
the number of iterations necessary to obtain convergence from GMRES:

<table align="center" class="doxtable">
<tr>
  <th></th>
  <th></th>
  <th colspan="1">$Q_1$</th>
  <th colspan="7">Smoother (smoothing steps)</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th></th>
  <th colspan="3">Jacobi (6)</th>
  <th></th>
  <th colspan="3">Block Jacobi (3)</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th></th>
  <th colspan="3">Renumbering Strategy</th>
  <th></th>
  <th colspan="3">Renumbering Strategy</th>
</tr>
<tr>
  <th>Cells</th>
  <th></th>
  <th>DoFs</th>
  <th>Downstream</th>
  <th>Random</th>
  <th>Upstream</th>
  <th></th>
  <th>Downstream</th>
  <th>Random</th>
  <th>Upstream</th>
</tr>
<tr>
  <th>32</th>
  <th></th>
  <th>48</th>
  <td>3</th>
  <td>3</th>
  <td>3</th>
  <th></th>
  <td>3</th>
  <td>3</th>
  <td>3</th>
</tr>
<tr>
  <th>128</th>
  <th></th>
  <th>160</th>
  <td>6</th>
  <td>6</th>
  <td>6</th>
  <th></th>
  <td>6</th>
  <td>6</th>
  <td>6</th>
</tr>
<tr>
  <th>512</th>
  <th></th>
  <th>576</th>
  <td>11</th>
  <td>11</th>
  <td>11</th>
  <th></th>
  <td>9</th>
  <td>9</th>
  <td>9</th>
</tr>
<tr>
  <th>2048</th>
  <th></th>
  <th>2176</th>
  <td>15</th>
  <td>15</th>
  <td>15</th>
  <th></th>
  <td>13</th>
  <td>13</th>
  <td>13</th>
</tr>
<tr>
  <th>8192</th>
  <th></th>
  <th>8448</th>
  <td>18</th>
  <td>18</th>
  <td>18</th>
  <th></th>
  <td>15</th>
  <td>15</th>
  <td>15</th>
</tr>
<tr>
  <th>32768</th>
  <th></th>
  <th>33280</th>
  <td>20</th>
  <td>20</th>
  <td>20</th>
  <th></th>
  <td>16</th>
  <td>16</th>
  <td>16</th>
</tr>
<tr>
  <th>131072</th>
  <th></th>
  <th>132096</th>
  <td>20</th>
  <td>20</th>
  <td>20</th>
  <th></th>
  <td>16</th>
  <td>16</th>
  <td>16</th>
</tr>
</table>

We see that renumbering the
DoFs/cells has no effect on convergence speed. This is because these
smoothers compute operations on each DoF (point-smoother) or cell
(block-smoother) independently and add up the results. Since we can
define these smoothers as an application of a sum of matrices, and
matrix addition is commutative, the order at which we sum the
different components will not affect the end result.

On the other hand, the situation is different for multiplicative smoothers:

<table align="center" class="doxtable">
<tr>
  <th></th>
  <th></th>
  <th colspan="1">$Q_1$</th>
  <th colspan="7">Smoother (smoothing steps)</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th></th>
  <th colspan="3">SOR (3)</th>
  <th></th>
  <th colspan="3">Block SOR (1)</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th></th>
  <th colspan="3">Renumbering Strategy</th>
  <th></th>
  <th colspan="3">Renumbering Strategy</th>
</tr>
<tr>
  <th>Cells</th>
  <th></th>
  <th>DoFs</th>
  <th>Downstream</th>
  <th>Random</th>
  <th>Upstream</th>
  <th></th>
  <th>Downstream</th>
  <th>Random</th>
  <th>Upstream</th>
</tr>
<tr>
  <th>32</th>
  <th></th>
  <th>48</th>
  <td>2</th>
  <td>2</th>
  <td>3</th>
  <th></th>
  <td>2</th>
  <td>2</th>
  <td>3</th>
</tr>
<tr>
  <th>128</th>
  <th></th>
  <th>160</th>
  <td>5</th>
  <td>5</th>
  <td>7</th>
  <th></th>
  <td>5</th>
  <td>5</th>
  <td>7</th>
</tr>
<tr>
  <th>512</th>
  <th></th>
  <th>576</th>
  <td>7</th>
  <td>9</th>
  <td>11</th>
  <th></th>
  <td>7</th>
  <td>7</th>
  <td>12</th>
</tr>
<tr>
  <th>2048</th>
  <th></th>
  <th>2176</th>
  <td>10</th>
  <td>12</th>
  <td>15</th>
  <th></th>
  <td>8</th>
  <td>10</th>
  <td>17</th>
</tr>
<tr>
  <th>8192</th>
  <th></th>
  <th>8448</th>
  <td>11</th>
  <td>15</th>
  <td>19</th>
  <th></th>
  <td>10</th>
  <td>11</th>
  <td>20</th>
</tr>
<tr>
  <th>32768</th>
  <th></th>
  <th>33280</th>
  <td>12</th>
  <td>16</th>
  <td>20</th>
  <th></th>
  <td>10</th>
  <td>12</th>
  <td>21</th>
</tr>
<tr>
  <th>131072</th>
  <th></th>
  <th>132096</th>
  <td>12</th>
  <td>16</th>
  <td>19</th>
  <th></th>
  <td>11</th>
  <td>12</th>
  <td>21</th>
</tr>
</table>

Here, we can speed up
convergence by renumbering the DoFs/cells in the advection direction,
and similarly, we can slow down convergence if we do the renumbering
in the opposite direction. This is because advection-dominated
problems have a directional flow of information (in the advection
direction) which, given the right renumbering of DoFs/cells,
multiplicative methods are able to capture.

This feature of multiplicative methods is, however, dependent on the
value of $\varepsilon$. As we increase $\varepsilon$ and the problem
becomes more diffusion-dominated, we have a more uniform propagation
of information over the mesh and there is a diminished advantage for
renumbering in the advection direction. On the opposite end, in the
extreme case of $\varepsilon=0$ (advection-only), we have a 1st-order
PDE and multiplicative methods with the right renumbering become
effective solvers: A correct downstream numbering may lead to methods
that require only a single iteration because information can be
propagated from the inflow boundary downstream, with no information
transport in the opposite direction. (Note, however, that in the case
of $\varepsilon=0$, special care must be taken for the boundary
conditions in this case).


<a name="Pointvsblocksmoothers"></a><h4> %Point vs. block smoothers </h4>


We will limit the results to runs using the downstream
renumbering. Here is a cross comparison of all four smoothers for both
$Q_1$ and $Q_3$ elements:

<table align="center" class="doxtable">
<tr>
  <th></th>
  <td></th>
  <th colspan="1">$Q_1$</th>
  <th colspan="4">Smoother (smoothing steps)</th>
  <th></th>
  <th colspan="1">$Q_3$</th>
  <th colspan="4">Smoother (smoothing steps)</th>
</tr>
<tr>
  <th colspan="1">Cells</th>
  <td></th>
  <th colspan="1">DoFs</th>
  <th colspan="1">Jacobi (6)</th>
  <th colspan="1">Block Jacobi (3)</th>
  <th colspan="1">SOR (3)</th>
  <th colspan="1">Block SOR (1)</th>
  <th></th>
  <th colspan="1">DoFs</th>
  <th colspan="1">Jacobi (6)</th>
  <th colspan="1">Block Jacobi (3)</th>
  <th colspan="1">SOR (3)</th>
  <th colspan="1">Block SOR (1)</th>
</tr>
<tr>
  <th>32</th>
  <td></th>
  <th>48</th>
  <td>3</th>
  <td>3</th>
  <td>2</th>
  <td>2</th>
  <td></th>
  <th>336</th>
  <td>15</th>
  <td>14</th>
  <td>15</th>
  <td>6</th>
</tr>
<tr>
  <th>128</th>
  <td></th>
  <th>160</th>
  <td>6</th>
  <td>6</th>
  <td>5</th>
  <td>5</th>
  <td></th>
  <th>1248</th>
  <td>23</th>
  <td>18</th>
  <td>21</th>
  <td>9</th>
</tr>
<tr>
  <th>512</th>
  <td></th>
  <th>576</th>
  <td>11</th>
  <td>9</th>
  <td>7</th>
  <td>7</th>
  <td></th>
  <th>4800</th>
  <td>29</th>
  <td>21</th>
  <td>28</th>
  <td>9</th>
</tr>
<tr>
  <th>2048</th>
  <td></th>
  <th>2176</th>
  <td>15</th>
  <td>13</th>
  <td>10</th>
  <td>8</th>
  <td></th>
  <th>18816</th>
  <td>33</th>
  <td>22</th>
  <td>32</th>
  <td>9</th>
</tr>
<tr>
  <th>8192</th>
  <td></th>
  <th>8448</th>
  <td>18</th>
  <td>15</th>
  <td>11</th>
  <td>10</th>
  <td></th>
  <th>74496</th>
  <td>35</th>
  <td>22</th>
  <td>34</th>
  <td>10</th>
</tr>
<tr>
  <th>32768</th>
  <td></th>
  <th>33280</th>
  <td>20</th>
  <td>16</th>
  <td>12</th>
  <td>10</th>
  <td></th>
</tr>
<tr>
  <th>131072</th>
  <td></th>
  <th>132096</th>
  <td>20</th>
  <td>16</th>
  <td>12</th>
  <td>11</th>
  <td></th>
</tr>
</table>

We see that for $Q_1$, both multiplicative smoothers require a smaller
combination of smoothing steps and iteration counts than either
additive smoother. However, when we increase the degree to a $Q_3$
element, there is a clear advantage for the block smoothers in terms
of the number of smoothing steps and iterations required to
solve. Specifically, the block SOR smoother gives constant iteration
counts over the degree, and the block Jacobi smoother only sees about
a 38% increase in iterations compared to 75% and 183% for Jacobi and
SOR respectively.

<a name="Cost"></a><h3> Cost </h3>


Iteration counts do not tell the full story in the optimality of a one
smoother over another. Obviously we must examine the cost of an
iteration. Block smoothers here are at a disadvantage as they are
having to construct and invert a cell matrix for each cell. Here is a
comparison of solve times for a $Q_3$ element with 74,496 DoFs:

<table align="center" class="doxtable">
<tr>
  <th colspan="1">$Q_3$</th>
  <th colspan="4">Smoother (smoothing steps)</th>
</tr>
<tr>
  <th colspan="1">DoFs</th>
  <th colspan="1">Jacobi (6)</th>
  <th colspan="1">Block Jacobi (3)</th>
  <th colspan="1">SOR (3)</th>
  <th colspan="1">Block SOR (1)</th>
</tr>
<tr>
  <th>74496</th>
  <td>0.68s</th>
  <td>5.82s</th>
  <td>1.18s</th>
  <td>1.02s</th>
</tr>
</table>

The smoother that requires the most iterations (Jacobi) actually takes
the shortest time (roughly 2/3 the time of the next fastest
method). This is because all that is required to apply a Jacobi
smoothing step is multiplication by a diagonal matrix which is very
cheap. On the other hand, while SOR requires over 3x more iterations
(each with 3x more smoothing steps) than block SOR, the times are
roughly equivalent, implying that a smoothing step of block SOR is
roughly 9x slower than a smoothing step of SOR. Lastly, block Jacobi
is almost 6x more expensive than block SOR, which intuitively makes
sense from the fact that 1 step of each method has the same cost
(inverting the cell matrices and either adding or multiply them
together), and block Jacobi has 3 times the number of smoothing steps per
iteration with 2 times the iterations.


<a name="Additionalpoints"></a><h3> Additional points </h3>


There are a few more important points to mention:

<ol>
<li> For a mesh distributed in parallel, multiplicative methods cannot
be executed over the entire domain. This is because they operate one
cell at a time, and downstream cells can only be handled once upstream
cells have already been done. This is fine on a single processor: The
processor just goes through the list of cells one after the
other. However, in parallel, it would imply that some processors are
idle because upstream processors have not finished doing the work on
cells upstream from the ones owned by the current processor. Once the
upstream processors are done, the downstream ones can start, but by
that time the upstream processors have no work left. In other words,
most of the time during these smoother steps, most processors are in
fact idle. This is not how one obtains good parallel scalability!

One can use a hybrid method where
a multiplicative smoother is applied on each subdomain, but as you
increase the number of subdomains, the method approaches the behavior
of an additive method. This is a major disadvantage to these methods.
</li>

<li> Current research into block smoothers suggest that soon we will be
able to compute the inverse of the cell matrices much cheaper than
what is currently being done inside deal.II. This research is based on
the fast diagonalization method (dating back to the 1960s) and has
been used in the spectral community for around 20 years (see, e.g., <a
href="https://doi.org/10.1007/s10915-004-4787-3"> Hybrid
Multigrid/Schwarz Algorithms for the Spectral Element Method by Lottes
and Fischer</a>). There are currently efforts to generalize these
methods to DG and make them more robust. Also, it seems that one
should be able to take advantage of matrix-free implementations and
the fact that, in the interior of the domain, cell matrices tend to
look very similar, allowing fewer matrix inverse computations.
</li>
</ol>

Combining 1. and 2. gives a good reason for expecting that a method
like block Jacobi could become very powerful in the future, even
though currently for these examples it is quite slow.


<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="ConstantiterationsforQsub5sub"></a><h4> Constant iterations for Q<sub>5</sub> </h4>


Change the number of smoothing steps and the smoother relaxation
parameter (set in <code>Smoother::AdditionalData()</code> inside
<code>create_smoother()</code>, only necessary for point smoothers) so
that we maintain a constant number of iterations for a $Q_5$ element.

<a name="Effectivenessofrenumberingforchangingepsilon"></a><h4> Effectiveness of renumbering for changing epsilon </h4>


Increase/decrease the parameter "Epsilon" in the `.prm` files of the
multiplicative methods and observe for which values renumbering no
longer influences convergence speed.

<a name="Meshadaptivity"></a><h4> Mesh adaptivity </h4>


The code is set up to work correctly with an adaptively refined mesh (the
interface matrices are created and set). Devise a suitable refinement
criterium or try the KellyErrorEstimator class.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-63.cc"
*/
