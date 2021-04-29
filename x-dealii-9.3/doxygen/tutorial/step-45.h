/**
@page step_45 The step-45 tutorial program
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
        <li><a href="#Settingupperiodicityconstraintsondistributedtriangulations">Setting up periodicity constraints on distributed triangulations</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Daniel Arndt and Matthias Maier.</i>
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


In this example we present how to use periodic boundary conditions in
deal.II. Periodic boundary conditions are algebraic constraints that
typically occur in computations on representative regions of a larger
domain that repeat in one or more directions.

An example is the simulation of the electronic structure of photonic
crystals, because they have a lattice-like structure and, thus, it often
suffices to do the actual computation on only one box of the lattice. To
be able to proceed this way one has to assume that the model can be
periodically extended to the other boxes; this requires the solution to
have a periodic structure.

<a name="Procedure"></a>
<a name="Procedure"></a><h1>Procedure</h1>


deal.II provides a number of high level entry points to impose periodic
boundary conditions.
The general approach to apply periodic boundary conditions consists of
three steps (see also the
@ref GlossPeriodicConstraints "Glossary entry on periodic boundary conditions"):
-# Create a mesh
-# Identify those pairs of faces on different parts of the boundary across which
   the solution should be symmetric, using GridTools::collect_periodic_faces()
-# Add the periodicity information to the mesh
   using parallel::distributed::Triangulation::add_periodicity()
-# Add periodicity constraints using DoFTools::make_periodicity_constraints()

The second and third step are necessary for parallel meshes using the
parallel::distributed::Triangulation class
to ensure that cells on opposite sides of the domain but connected by periodic
faces are part of the ghost layer if one of them is stored on the local processor.
If the Triangulation is not a parallel::distributed::Triangulation,
these steps are not necessary.

The first step consists of collecting matching periodic faces and storing them in
a <code>std::vector</code> of GridTools::PeriodicFacePair. This is done with the
function GridTools::collect_periodic_faces() that can be invoked for example
like this:
@code
GridTools::collect_periodic_faces(dof_handler,
                                  b_id1,
                                  b_id2,
                                  direction,
                                  matched_pairs,
                                  offset = <default value>,
                                  matrix = <default value>,
                                  first_vector_components = <default value>);
@endcode

This call loops over all faces of the container dof_handler on the periodic
boundaries with boundary indicator @p b_id1 and @p b_id2,
respectively. (You can assign these boundary indicators by hand after
creating the coarse mesh, see
@ref GlossBoundaryIndicator "Boundary indicator". Alternatively, you
can also let many of the functions in namespace GridGenerator do this
for if you specify the "colorize" flag; in that case, these functions
will assign different boundary indicators to different parts of the
boundary, with the details typically spelled out in the documentation
of these functions.)

Concretely, if $\text{vertices}_{1/2}$ are the vertices of two faces
$\text{face}_{1/2}$, then the function call above will match pairs of
faces (and dofs) such that the difference between $\text{vertices}_2$
and $matrix\cdot \text{vertices}_1+\text{offset}$ vanishes in every
component apart from direction and stores the resulting pairs with
associated data in @p matched_pairs. (See
GridTools::orthogonal_equality() for detailed information about the
matching process.)

Consider, for example, the colored unit square $\Omega=[0,1]^2$ with boundary
indicator 0 on the left, 1 on the right, 2 on the bottom and 3 on the top
faces. (See the documentation of GridGenerator::hyper_cube() for this
convention on how boundary indicators are assigned.) Then,
@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs);
@endcode
would yield periodicity constraints such that $u(0,y)=u(1,y)$ for all
$y\in[0,1]$.

If we instead consider the parallelogram given by the convex hull of
$(0,0)$, $(1,1)$, $(1,2)$, $(0,1)$ we can achieve the constraints
$u(0,y)=u(1,y+1)$ by specifying an @p offset:
@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(0.,1.));
@endcode
or
@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*arbitrary direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(1.,1.));
@endcode
Here, again, the assignment of boundary indicators 0 and 1 stems from
what GridGenerator::parallelogram() documents.

The resulting @p matched_pairs can be used in
DoFTools::make_periodicity_constraints for populating an AffineConstraints
object with periodicity constraints:
@code
DoFTools::make_periodicity_constraints(matched_pairs, constraints);
@endcode

Apart from this high level interface there are also variants of
DoFTools::make_periodicity_constraints available that combine those two
steps (see the variants of DofTools::make_periodicity_constraints).

There is also a low level interface to
DoFTools::make_periodicity_constraints if more flexibility is needed. The
low level variant allows to directly specify two faces that shall be
constrained:
@code
using namespace DoFTools;
make_periodicity_constraints(face_1,
                             face_2,
                             affine_constraints,
                             component_mask = <default value>;
                             face_orientation = <default value>,
                             face_flip = <default value>,
                             face_rotation = <default value>,
                             matrix = <default value>);
@endcode
Here, we need to specify the orientation of the two faces using
@p face_orientation, @p face_flip and @p face_orientation. For a closer description
have a look at the documentation of DoFTools::make_periodicity_constraints.
The remaining parameters are the same as for the high level interface apart
from the self-explaining @p component_mask and @p affine_constraints.


<a name="problem"></a>
<a name="Apracticalexample"></a><h1>A practical example</h1>


In the following, we show how to use the above functions in a more involved
example. The task is to enforce rotated periodicity constraints for the
velocity component of a Stokes flow.

On a quarter-circle defined by $\Omega=\{{\bf x}\in(0,1)^2:\|{\bf x}\|\in (0.5,1)\}$ we are
going to solve the Stokes problem
@f{eqnarray*}
  -\Delta \; \textbf{u} + \nabla p &=& (\exp(-100\|{\bf x}-(.75,0.1)^T\|^2),0)^T, \\
  -\textrm{div}\;  \textbf{u}&=&0,\\
  \textbf{u}|_{\Gamma_1}&=&{\bf 0},
@f}
where the boundary $\Gamma_1$ is defined as $\Gamma_1 \dealcoloneq \{x\in \partial\Omega: \|x\|\in\{0.5,1\}\}$.
For the remaining parts of the boundary we are going to use periodic boundary conditions, i.e.
@f{align*}
  u_x(0,\nu)&=-u_y(\nu,0)&\nu&\in[0,1]\\
  u_y(0,\nu)&=u_x(\nu,0)&\nu&\in[0,1].
@f}

The mesh will be generated by GridGenerator::quarter_hyper_shell(),
which also documents how it assigns boundary indicators to its various
boundaries if its `colorize` argument is set to `true`.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * This example program is a slight modification of step-22 running in parallel
 * using Trilinos to demonstrate the usage of periodic boundary conditions in
 * deal.II. We thus omit to discuss the majority of the source code and only
 * comment on the parts that deal with periodicity constraints. For the rest
 * have a look at step-22 and the full source code at the bottom.
 * 

 * 
 * In order to implement periodic boundary conditions only two functions
 * have to be modified:
 * - <code>StokesProblem<dim>::setup_dofs()</code>:
 * To populate an AffineConstraints object with periodicity constraints
 * - <code>StokesProblem<dim>::create_mesh()</code>:
 * To supply a distributed triangulation with periodicity information.
 * 

 * 
 * The rest of the program is identical to step-22, so let us skip this part
 * and only show these two functions in the following. (The full program can be
 * found in the "Plain program" section below, though.)
 * 

 * 
 * 

 * 
 *   

 * 
 * 
 * <a name="Settingupperiodicityconstraintsondistributedtriangulations"></a> 
 * <h3>Setting up periodicity constraints on distributed triangulations</h3>
 * 
 * @code
 *   template <int dim>
 *   void StokesProblem<dim>::create_mesh()
 *   {
 *     Point<dim>   center;
 *     const double inner_radius = .5;
 *     const double outer_radius = 1.;
 * 
 *     GridGenerator::quarter_hyper_shell(
 *       triangulation, center, inner_radius, outer_radius, 0, true);
 * 
 * @endcode
 * 
 * Before we can prescribe periodicity constraints, we need to ensure that
 * cells on opposite sides of the domain but connected by periodic faces are
 * part of the ghost layer if one of them is stored on the local processor.
 * At this point we need to think about how we want to prescribe
 * periodicity. The vertices $\text{vertices}_2$ of a face on the left
 * boundary should be matched to the vertices $\text{vertices}_1$ of a face
 * on the lower boundary given by $\text{vertices}_2=R\cdot
 * \text{vertices}_1+b$ where the rotation matrix $R$ and the offset $b$ are
 * given by
 * @f{align*}
 * R=\begin{pmatrix}
 * 0&1\\-1&0
 * \end{pmatrix},
 * \quad
 * b=\begin{pmatrix}0&0\end{pmatrix}.
 * @f}
 * The data structure we are saving the resulting information into is here
 * based on the Triangulation.
 * 
 * @code
 *     std::vector<GridTools::PeriodicFacePair<
 *       typename parallel::distributed::Triangulation<dim>::cell_iterator>>
 *       periodicity_vector;
 * 
 *     FullMatrix<double> rotation_matrix(dim);
 *     rotation_matrix[0][1] = 1.;
 *     rotation_matrix[1][0] = -1.;
 * 
 *     GridTools::collect_periodic_faces(triangulation,
 *                                       2,
 *                                       3,
 *                                       1,
 *                                       periodicity_vector,
 *                                       Tensor<1, dim>(),
 *                                       rotation_matrix);
 * 
 * @endcode
 * 
 * Now telling the triangulation about the desired periodicity is
 * particularly easy by just calling
 * parallel::distributed::Triangulation::add_periodicity.
 * 
 * @code
 *     triangulation.add_periodicity(periodicity_vector);
 * 
 *     triangulation.refine_global(4 - dim);
 *   }
 * 
 * 
 *   template <int dim>
 *   void StokesProblem<dim>::setup_dofs()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::vector<unsigned int> block_component(dim + 1, 0);
 *     block_component[dim] = 1;
 *     DoFRenumbering::component_wise(dof_handler, block_component);
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
 *     const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];
 * 
 *     {
 *       owned_partitioning.clear();
 *       IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
 *       owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u));
 *       owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u + n_p));
 * 
 *       relevant_partitioning.clear();
 *       IndexSet locally_relevant_dofs;
 *       DoFTools::extract_locally_relevant_dofs(dof_handler,
 *                                               locally_relevant_dofs);
 *       relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u));
 *       relevant_partitioning.push_back(
 *         locally_relevant_dofs.get_view(n_u, n_u + n_p));
 * 
 *       constraints.clear();
 *       constraints.reinit(locally_relevant_dofs);
 * 
 *       FEValuesExtractors::Vector velocities(0);
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *       VectorTools::interpolate_boundary_values(mapping,
 *                                                dof_handler,
 *                                                0,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 *       VectorTools::interpolate_boundary_values(mapping,
 *                                                dof_handler,
 *                                                1,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 * 
 * @endcode
 * 
 * After we provided the mesh with the necessary information for the
 * periodicity constraints, we are now able to actual create them. For
 * describing the matching we are using the same approach as before, i.e.,
 * the $\text{vertices}_2$ of a face on the left boundary should be
 * matched to the vertices
 * $\text{vertices}_1$ of a face on the lower boundary given by
 * $\text{vertices}_2=R\cdot \text{vertices}_1+b$ where the rotation
 * matrix $R$ and the offset $b$ are given by
 * @f{align*}
 * R=\begin{pmatrix}
 * 0&1\\-1&0
 * \end{pmatrix},
 * \quad
 * b=\begin{pmatrix}0&0\end{pmatrix}.
 * @f}
 * These two objects not only describe how faces should be matched but
 * also in which sense the solution should be transformed from
 * $\text{face}_2$ to
 * $\text{face}_1$.
 * 
 * @code
 *       FullMatrix<double> rotation_matrix(dim);
 *       rotation_matrix[0][1] = 1.;
 *       rotation_matrix[1][0] = -1.;
 * 
 *       Tensor<1, dim> offset;
 * 
 * @endcode
 * 
 * For setting up the constraints, we first store the periodicity
 * information in an auxiliary object of type
 * <code>std::vector@<GridTools::PeriodicFacePair<typename
 * DoFHandler@<dim@>::%cell_iterator@> </code>. The periodic boundaries
 * have the boundary indicators 2 (x=0) and 3 (y=0). All the other
 * parameters we have set up before. In this case the direction does not
 * matter. Due to $\text{vertices}_2=R\cdot \text{vertices}_1+b$ this is
 * exactly what we want.
 * 
 * @code
 *       std::vector<
 *         GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
 *         periodicity_vector;
 * 
 *       const unsigned int direction = 1;
 * 
 *       GridTools::collect_periodic_faces(dof_handler,
 *                                         2,
 *                                         3,
 *                                         direction,
 *                                         periodicity_vector,
 *                                         offset,
 *                                         rotation_matrix);
 * 
 * @endcode
 * 
 * Next, we need to provide information on which vector valued components
 * of the solution should be rotated. Since we choose here to just
 * constraint the velocity and this starts at the first component of the
 * solution vector, we simply insert a 0:
 * 
 * @code
 *       std::vector<unsigned int> first_vector_components;
 *       first_vector_components.push_back(0);
 * 
 * @endcode
 * 
 * After setting up all the information in periodicity_vector all we have
 * to do is to tell make_periodicity_constraints to create the desired
 * constraints.
 * 
 * @code
 *       DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector,
 *                                                        constraints,
 *                                                        fe.component_mask(
 *                                                          velocities),
 *                                                        first_vector_components);
 * 
 *       VectorTools::interpolate_boundary_values(mapping,
 *                                                dof_handler,
 *                                                0,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 *       VectorTools::interpolate_boundary_values(mapping,
 *                                                dof_handler,
 *                                                1,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                fe.component_mask(velocities));
 *     }
 * 
 *     constraints.close();
 * 
 *     {
 *       TrilinosWrappers::BlockSparsityPattern bsp(owned_partitioning,
 *                                                  owned_partitioning,
 *                                                  relevant_partitioning,
 *                                                  mpi_communicator);
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (!((c == dim) && (d == dim)))
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler,
 *                                       coupling,
 *                                       bsp,
 *                                       constraints,
 *                                       false,
 *                                       Utilities::MPI::this_mpi_process(
 *                                         mpi_communicator));
 * 
 *       bsp.compress();
 * 
 *       system_matrix.reinit(bsp);
 *     }
 * 
 *     {
 *       TrilinosWrappers::BlockSparsityPattern preconditioner_bsp(
 *         owned_partitioning,
 *         owned_partitioning,
 *         relevant_partitioning,
 *         mpi_communicator);
 * 
 *       Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if ((c == dim) && (d == dim))
 *             preconditioner_coupling[c][d] = DoFTools::always;
 *           else
 *             preconditioner_coupling[c][d] = DoFTools::none;
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler,
 *                                       preconditioner_coupling,
 *                                       preconditioner_bsp,
 *                                       constraints,
 *                                       false,
 *                                       Utilities::MPI::this_mpi_process(
 *                                         mpi_communicator));
 * 
 *       preconditioner_bsp.compress();
 * 
 *       preconditioner_matrix.reinit(preconditioner_bsp);
 *     }
 * 
 *     system_rhs.reinit(owned_partitioning, mpi_communicator);
 *     solution.reinit(owned_partitioning,
 *                     relevant_partitioning,
 *                     mpi_communicator);
 *   }
 * 
 * @endcode
 * 
 * The rest of the program is then again identical to step-22. We will omit
 * it here now, but as before, you can find these parts in the "Plain program"
 * section below.
 * 

 * 
<a name="Results"></a><h1>Results</h1>


The created output is not very surprising. We simply see that the solution is
periodic with respect to the left and lower boundary:

<img src="https://www.dealii.org/images/steps/developer/step-45.periodic.png" alt="">

Without the periodicity constraints we would have ended up with the following solution:

<img src="https://www.dealii.org/images/steps/developer/step-45.non_periodic.png" alt="">
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-45.cc"
*/
