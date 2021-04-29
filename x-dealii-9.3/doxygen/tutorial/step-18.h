/**
@page step_18 The step-18 tutorial program
This tutorial depends on step-17.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Quasistaticelasticdeformation">Quasistatic elastic deformation</a>
      <ul>
        <li><a href="#Motivationofthemodel">Motivation of the model</a>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Updatingthestressvariable">Updating the stress variable</a>
      </ul>
        <li><a href="#Parallelgraphicaloutput">Parallel graphical output</a>
        <li><a href="#Atriangulationwithautomaticpartitioning">A triangulation with automatic partitioning</a>
        <li><a href="#Overallstructureoftheprogram">Overall structure of the program</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ThecodePointHistorycodeclass">The <code>PointHistory</code> class</a>
        <li><a href="#Thestressstraintensor">The stress-strain tensor</a>
        <li><a href="#Auxiliaryfunctions">Auxiliary functions</a>
        <li><a href="#ThecodeTopLevelcodeclass">The <code>TopLevel</code> class</a>
        <li><a href="#ThecodeBodyForcecodeclass">The <code>BodyForce</code> class</a>
        <li><a href="#ThecodeIncrementalBoundaryValuecodeclass">The <code>IncrementalBoundaryValue</code> class</a>
        <li><a href="#ImplementationofthecodeTopLevelcodeclass">Implementation of the <code>TopLevel</code> class</a>
      <ul>
        <li><a href="#Thepublicinterface">The public interface</a>
        <li><a href="#TopLevelcreate_coarse_grid">TopLevel::create_coarse_grid</a>
        <li><a href="#TopLevelsetup_system">TopLevel::setup_system</a>
        <li><a href="#TopLevelassemble_system">TopLevel::assemble_system</a>
        <li><a href="#TopLevelsolve_timestep">TopLevel::solve_timestep</a>
        <li><a href="#TopLevelsolve_linear_problem">TopLevel::solve_linear_problem</a>
        <li><a href="#TopLeveloutput_results">TopLevel::output_results</a>
        <li><a href="#TopLeveldo_initial_timestep">TopLevel::do_initial_timestep</a>
        <li><a href="#TopLeveldo_timestep">TopLevel::do_timestep</a>
        <li><a href="#TopLevelrefine_initial_grid">TopLevel::refine_initial_grid</a>
        <li><a href="#TopLevelmove_mesh">TopLevel::move_mesh</a>
        <li><a href="#TopLevelsetup_quadrature_point_history">TopLevel::setup_quadrature_point_history</a>
        <li><a href="#TopLevelupdate_quadrature_point_history">TopLevel::update_quadrature_point_history</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
      <ul>
        <li><a href="#Plasticitymodels">Plasticity models</a>
        <li><a href="#Stabilizationissues">Stabilization issues</a>
        <li><a href="#Refinementduringtimesteps">Refinement during timesteps</a>
        <li><a href="#Ensuringmeshregularity">Ensuring mesh regularity</a>
    </ul>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>



This tutorial program is another one in the series on the elasticity problem
that we have already started with step-8 and step-17. It extends it into two
different directions: first, it solves the quasistatic but time dependent
elasticity problem for large deformations with a Lagrangian mesh movement
approach. Secondly, it shows some more techniques for solving such problems
using %parallel processing with PETSc's linear algebra. In addition to this,
we show how to work around one of the two major bottlenecks of step-17, namely
that we generated graphical output from only one process, and that this scaled
very badly with larger numbers of processes and on large problems. (The other
bottleneck, namely that every processor has to hold the entire mesh and
DoFHandler, is addressed in step-40.) Finally, a
good number of assorted improvements and techniques are demonstrated that have
not been shown yet in previous programs.

As before in step-17, the program runs just as fine on a single sequential
machine as long as you have PETSc installed. Information on how to tell
deal.II about a PETSc installation on your system can be found in the deal.II
README file, which is linked to from the <a href="../../index.html">main
documentation page</a>
in your installation of deal.II, or on <a href="http://www.dealii.org/">the
deal.II webpage</a>.


<a name="Quasistaticelasticdeformation"></a><h3>Quasistatic elastic deformation</h3>


<a name="Motivationofthemodel"></a><h4>Motivation of the model</h4>


In general, time-dependent small elastic deformations are described by the
elastic wave equation
@f[
  \rho \frac{\partial^2 \mathbf{u}}{\partial t^2}
  + c \frac{\partial \mathbf{u}}{\partial t}
  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) = \mathbf{f}
  \qquad
  \textrm{in}\ \Omega,
@f]
where $\mathbf{u}=\mathbf{u} (\mathbf{x},t)$ is the deformation of the body, $\rho$
and $c$ the density and attenuation coefficient, and $\mathbf{f}$ external forces.
In addition, initial conditions
@f[
  \mathbf{u}(\cdot, 0) = \mathbf{u}_0(\cdot)
  \qquad
  \textrm{on}\ \Omega,
@f]
and Dirichlet (displacement) or Neumann (traction) boundary conditions need
to be specified for a unique solution:
@f{eqnarray*}
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega\backslash\Gamma_D.
@f}
In above formulation, $\varepsilon(\mathbf{u})= \frac 12 (\nabla \mathbf{u} + \nabla
\mathbf{u}^T)$ is the symmetric gradient of the displacement, also called the
<em>strain</em>. $C$ is a tensor of rank 4, called the <em>stress-strain
  tensor</em> (the inverse of the <a
  href="https://en.wikipedia.org/wiki/Hooke%27s_law#Hooke's_law_for_continuous_media"><em>compliance
  tensor</em></a>)
that contains knowledge of the elastic strength of the material; its
symmetry properties make sure that it maps symmetric tensors of rank 2
(&ldquo;matrices&rdquo; of dimension $d$, where $d$ is the spatial dimensionality) onto
symmetric tensors of the same rank. We will comment on the roles of the strain
and stress tensors more below. For the moment it suffices to say that we
interpret the term $\textrm{div}\  ( C \varepsilon(\mathbf{u}))$ as the vector with
components $\frac \partial{\partial x_j} C_{ijkl} \varepsilon(\mathbf{u})_{kl}$,
where summation over indices $j,k,l$ is implied.

The quasistatic limit of this equation is motivated as follows: each small
perturbation of the body, for example by changes in boundary condition or the
forcing function, will result in a corresponding change in the configuration
of the body. In general, this will be in the form of waves radiating away from
the location of the disturbance. Due to the presence of the damping term,
these waves will be attenuated on a time scale of, say, $\tau$. Now, assume
that all changes in external forcing happen on times scales that are
much larger than $\tau$. In that case, the dynamic nature of the change is
unimportant: we can consider the body to always be in static equilibrium,
i.e. we can assume that at all times the body satisfies
@f{eqnarray*}
  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega,
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N.
@f}
Note that the differential equation does not contain any time derivatives any
more -- all time dependence is introduced through boundary conditions and a
possibly time-varying force function $\mathbf{f}(\mathbf{x},t)$. The changes in
configuration can therefore be considered as being stationary
instantaneously. An alternative view of this is that $t$ is not really a time
variable, but only a time-like parameter that governs the evolution of the
problem.

While these equations are sufficient to describe small deformations, computing
large deformations is a little more complicated and, in general, leads
to nonlinear equations such as those treated in step-44. In the
following, let us consider some of the tools one would employ when
simulating problems in which the deformation becomes <i>large</i>.

@note The model we will consider below is not founded on anything that
would be mathematically sound: we will consider a model in which we
produce a small deformation, deform the physical coordinates of the
body by this deformation, and then consider the next loading step
again as a linear problem. This isn't consistent, since the assumption
of linearity implies that deformations are infinitesimal and so moving
around the vertices of our mesh by a finite amount before solving the
next linear problem is an inconsistent approach. We should therefore
note that it is not surprising that the equations discussed below
can't be found in the literature: <b>The model considered here has
little to do with reality!</b> On the other hand, the implementation
techniques we consider are very much what one would need to use when
implementing a <i>real</i> model, as we will see in step-44.


To come back to defining our "artificial" model, let us first
introduce a tensorial stress variable $\sigma$, and write the differential
equations in terms of the stress:
@f{eqnarray*}
  - \textrm{div}\  \sigma &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega(t),
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega(t),
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega(t)\backslash\Gamma_D.
@f}
Note that these equations are posed on a domain $\Omega(t)$ that
changes with time, with the boundary moving according to the
displacements $\mathbf{u}(\mathbf{x},t)$ of the points on the boundary. To
complete this system, we have to specify the incremental relationship between
the stress and the strain, as follows:
<a name="step_18.stress-strain"></a>
@f[
  \dot\sigma = C \varepsilon (\dot{\mathbf{u}}),
  \qquad
  \qquad
  \textrm{[stress-strain]}
@f]
where a dot indicates a time derivative. Both the stress $\sigma$ and the
strain $\varepsilon(\mathbf{u})$ are symmetric tensors of rank 2.


<a name="Timediscretization"></a><h4>Time discretization</h4>


Numerically, this system is solved as follows: first, we discretize
the time component using a backward Euler scheme. This leads to a
discrete equilibrium of force at time step $n$:
@f[
  -\textrm{div}\  \sigma^n = f^n,
@f]
where
@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n),
@f]
and $\Delta \mathbf{u}^n$ the incremental displacement for time step
$n$. In addition, we have to specify initial data $\mathbf{u}(\cdot,0)=\mathbf{u}_0$.
This way, if we want to solve for the displacement increment, we
have to solve the following system:
@f{align*}
  - \textrm{div}\   C \varepsilon(\Delta\mathbf{u}^n) &= \mathbf{f} + \textrm{div}\  \sigma^{n-1}
  \qquad
  &&\textrm{in}\ \Omega(t_{n-1}),
  \\
  \Delta \mathbf{u}^n(\mathbf{x},t) &= \mathbf{d}(\mathbf{x},t_n) - \mathbf{d}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_D\subset\partial\Omega(t_{n-1}),
  \\
  \mathbf{n} \ C \varepsilon(\Delta \mathbf{u}^n(\mathbf{x},t)) &= \mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_N=\partial\Omega(t_{n-1})\backslash\Gamma_D.
@f}
The weak form of this set of equations, which as usual is the basis for the
finite element formulation, reads as follows: find $\Delta \mathbf{u}^n \in
\{v\in H^1(\Omega(t_{n-1}))^d: v|_{\Gamma_D}=\mathbf{d}(\cdot,t_n) - \mathbf{d}(\cdot,t_{n-1})\}$
such that
<a name="step_18.linear-system"></a>
@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}
  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  \\
  &\qquad
  +(\mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1}), \varphi)_{\Gamma_N}
  +(\sigma^{n-1} \mathbf{n}, \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.
@f}
Using that $\sigma^{n-1} \mathbf{n}
            = [C \varepsilon(\mathbf{u}^{n-1})] \mathbf{n}
            = \mathbf{b}(\mathbf x, t_{n-1})$,
these equations can be simplified to
@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}
  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  +(\mathbf{b}(\mathbf{x},t_n),t_{n-1}), \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.
  \qquad
  \qquad
  \textrm{[linear-system]}
@f}

We note that, for simplicity, in the program we will always assume that there
are no boundary forces, i.e. $\mathbf{b} = 0$, and that the deformation of the
body is driven by body forces $\mathbf{f}$ and prescribed boundary displacements
$\mathbf{d}$ alone. It is also worth noting that when integrating by parts, we
would get terms of the form $(C \varepsilon(\Delta\mathbf{u}^n), \nabla \varphi
)_{\Omega(t_{n-1})}$, but that we replace them with the term involving the
symmetric gradient $\varepsilon(\varphi)$ instead of $\nabla\varphi$. Due to
the symmetry of $C$, the two terms are mathematically equivalent, but
the symmetric version avoids the potential for round-off errors making
the resulting matrix slightly non-symmetric.

The system at time step $n$, to be solved on the old domain
$\Omega(t_{n-1})$, has exactly the form of a stationary elastic
problem, and is therefore similar to what we have already implemented
in previous example programs. We will therefore not comment on the
space discretization beyond saying that we again use lowest order
continuous finite elements.

There are differences, however:
<ol>
  <li> We have to move (update) the mesh after each time step, in order to be
  able to solve the next time step on a new domain;

  <li> We need to know $\sigma^{n-1}$ to compute the next incremental
  displacement, i.e. we need to compute it at the end of the time step
  to make sure it is available for the next time step. Essentially,
  the stress variable is our window to the history of deformation of
  the body.
</ol>
These two operations are done in the functions <code>move_mesh</code> and
<code>update_quadrature_point_history</code> in the program. While moving
the mesh is only a technicality, updating the stress is a little more
complicated and will be discussed in the next section.


<a name="Updatingthestressvariable"></a><h4>Updating the stress variable</h4>


As indicated above, we need to have the stress variable $\sigma^n$ available
when computing time step $n+1$, and we can compute it using
<a name="step_18.stress-update"></a>
@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update]}
@f]
There are, despite the apparent simplicity of this equation, two questions
that we need to discuss. The first concerns the way we store $\sigma^n$: even
if we compute the incremental updates $\Delta\mathbf{u}^n$ using lowest-order
finite elements, then its symmetric gradient $\varepsilon(\Delta\mathbf{u}^n)$ is
in general still a function that is not easy to describe. In particular, it is
not a piecewise constant function, and on general meshes (with cells that are
not rectangles %parallel to the coordinate axes) or with non-constant
stress-strain tensors $C$ it is not even a bi- or trilinear function. Thus, it
is a priori not clear how to store $\sigma^n$ in a computer program.

To decide this, we have to see where it is used. The only place where we
require the stress is in the term
$(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}$. In practice, we of
course replace this term by numerical quadrature:
@f[
  (\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  =
  \sum_{K\subset {T}}
  (\sigma^{n-1},\varepsilon(\varphi))_K
  \approx
  \sum_{K\subset {T}}
  \sum_q
  w_q \ \sigma^{n-1}(\mathbf{x}_q) : \varepsilon(\varphi(\mathbf{x}_q),
@f]
where $w_q$ are the quadrature weights and $\mathbf{x}_q$ the quadrature points on
cell $K$. This should make clear that what we really need is not the stress
$\sigma^{n-1}$ in itself, but only the values of the stress in the quadrature
points on all cells. This, however, is a simpler task: we only have to provide
a data structure that is able to hold one symmetric tensor of rank 2 for each
quadrature point on all cells (or, since we compute in parallel, all
quadrature points of all cells that the present MPI process &ldquo;owns&rdquo;). At the
end of each time step we then only have to evaluate $\varepsilon(\Delta \mathbf{u}^n(\mathbf{x}_q))$, multiply it by the stress-strain tensor $C$, and use the
result to update the stress $\sigma^n(\mathbf{x}_q)$ at quadrature point $q$.

The second complication is not visible in our notation as chosen above. It is
due to the fact that we compute $\Delta u^n$ on the domain $\Omega(t_{n-1})$,
and then use this displacement increment to both update the stress as well as
move the mesh nodes around to get to $\Omega(t_n)$ on which the next increment
is computed. What we have to make sure, in this context, is that moving the
mesh does not only involve moving around the nodes, but also making
corresponding changes to the stress variable: the updated stress is a variable
that is defined with respect to the coordinate system of the material in the
old domain, and has to be transferred to the new domain. The reason for this
can be understood as follows: locally, the incremental deformation $\Delta\mathbf{u}$ can be decomposed into three parts, a linear translation (the constant part
of the displacement increment field in the neighborhood of a point), a
dilational
component (that part of the gradient of the displacement field that has a
nonzero divergence), and a rotation. A linear translation of the material does
not affect the stresses that are frozen into it -- the stress values are
simply translated along. The dilational or compressional change produces a
corresponding stress update. However, the rotational component does not
necessarily induce a nonzero stress update (think, in 2d, for example of the
situation where $\Delta\mathbf{u}=(y, -x)^T$, with which $\varepsilon(\Delta
\mathbf{u})=0$). Nevertheless, if the material was prestressed in a certain
direction, then this direction will be rotated along with the material.  To
this end, we have to define a rotation matrix $R(\Delta \mathbf{u}^n)$ that
describes, in each point the rotation due to the displacement increments. It
is not hard to see that the actual dependence of $R$ on $\Delta \mathbf{u}^n$ can
only be through the curl of the displacement, rather than the displacement
itself or its full gradient (as mentioned above, the constant components of
the increment describe translations, its divergence the dilational modes, and
the curl the rotational modes). Since the exact form of $R$ is cumbersome, we
only state it in the program code, and note that the correct updating formula
for the stress variable is then
<a name="step_18.stress-update+rot"></a>
@f[
  \sigma^n
  =
  R(\Delta \mathbf{u}^n)^T
  [\sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n)]
  R(\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update+rot]}
@f]

Both stress update and rotation are implemented in the function
<code>update_quadrature_point_history</code> of the example program.


<a name="Parallelgraphicaloutput"></a><h3>Parallel graphical output</h3>


In step-17, the main bottleneck for %parallel computations as far as run time
is concerned
was that only the first processor generated output for the entire domain.
Since generating graphical output is expensive, this did not scale well when
larger numbers of processors were involved. We will address this here. (For a
definition of what it means for a program to "scale", see
@ref GlossParallelScaling "this glossary entry".)

Basically, what we need to do is let every process
generate graphical output for that subset of cells that it owns, write them
into separate files and have a way to display all files for a certain timestep
at the same time. This way the code produces one <code>.vtu</code> file per process per
time step. The two common VTK file viewers ParaView and VisIt both support
opening more than one <code>.vtu</code> file at once. To simplify the process of picking
the correct files and allow moving around in time, both support record files
that reference all files for a given timestep. Sadly, the record files have a
different format between VisIt and Paraview, so we write out both formats.

The code will generate the files <code>solution-TTTT.NNN.vtu</code>,
where <code>TTTT</code> is the timestep number (starting from 1) and
<code>NNN</code> is the process rank (starting from
0). These files contain the locally owned cells for the timestep and
processor. The files <code>solution-TTTT.visit</code> is the visit record
for timestep <code>TTTT</code>, while <code>solution-TTTT.pvtu</code> is
the same for ParaView. (More recent versions of VisIt can actually read
<code>.pvtu</code> files as well, but it doesn't hurt to output both
kinds of record files.) Finally, the file
<code>solution.pvd</code> is a special record only supported by ParaView that references
all time steps. So in ParaView, only solution.pvd needs to be opened, while
one needs to select the group of all .visit files in VisIt for the same
effect.


<a name="Atriangulationwithautomaticpartitioning"></a><h3>A triangulation with automatic partitioning</h3>


In step-17, we used a regular triangulation that was simply replicated on
every processor, and a corresponding DoFHandler. Both had no idea that they
were used in a %parallel context -- they just existed in their entirety
on every processor, and we argued that this was eventually going to be a
major memory bottleneck.

We do not address this issue here (we will do so in step-40) but make
the situation slightly more automated. In step-17, we created the triangulation
and then manually "partitioned" it, i.e., we assigned
@ref GlossSubdomainId "subdomain ids" to every cell that indicated which
@ref GlossMPIProcess "MPI process" "owned" the cell. Here, we use a class
parallel::shared::Triangulation that at least does this part automatically:
whenever you create or refine such a triangulation, it automatically
partitions itself among all involved processes (which it knows about because
you have to tell it about the @ref GlossMPICommunicator "MPI communicator"
that connects these processes upon construction of the triangulation).
Otherwise, the parallel::shared::Triangulation looks, for all practical
purposes, like a regular Triangulation object.

The convenience of using this class does not only result from being able
to avoid the manual call to GridTools::partition(). Rather, the DoFHandler
class now also knows that you want to use it in a parallel context, and
by default automatically enumerates degrees of freedom in such a way
that all DoFs owned by process zero come before all DoFs owned by process 1,
etc. In other words, you can also avoid the call to
DoFRenumbering::subdomain_wise().

There are other benefits. For example, because the triangulation knows that
it lives in a %parallel universe, it also knows that it "owns" certain
cells (namely, those whose subdomain id equals its MPI rank; previously,
the triangulation only stored these subdomain ids, but had no way to
make sense of them). Consequently, in the assembly function, you can
test whether a cell is "locally owned" (i.e., owned by the current
process, see @ref GlossLocallyOwnedCell) when you loop over all cells
using the syntax
@code
  if (cell->is_locally_owned())
@endcode
This knowledge extends to the DoFHandler object built on such triangulations,
which can then identify which degrees of freedom are locally owned
(see @ref GlossLocallyOwnedDof) via calls such as
DoFHandler::compute_n_locally_owned_dofs_per_processor() and
DoFTools::extract_locally_relevant_dofs(). Finally, the DataOut class
also knows how to deal with such triangulations and will simply skip
generating graphical output on cells not locally owned.

Of course, as has been noted numerous times in the discussion in step-17,
keeping the entire triangulation on every process will not scale: large
problems may simply not fit into each process's memory any more, even if
we have sufficiently many processes around to solve them in a reasonable
time. In such cases, the parallel::shared::Triangulation is no longer
a reasonable basis for computations and we will show in step-40 how the
parallel::distributed::Triangulation class can be used to work around
this, namely by letting each process store only a <i>part</i> of the
triangulation.


<a name="Overallstructureoftheprogram"></a><h3>Overall structure of the program</h3>


The overall structure of the program can be inferred from the <code>run()</code>
function that first calls <code>do_initial_timestep()</code> for the first time
step, and then <code>do_timestep()</code> on all subsequent time steps. The
difference between these functions is only that in the first time step we
start on a coarse mesh, solve on it, refine the mesh adaptively, and then
start again with a clean state on that new mesh. This procedure gives us a
better starting mesh, although we should of course keep adapting the mesh as
iterations proceed -- this isn't done in this program, but commented on below.

The common part of the two functions treating time steps is the following
sequence of operations on the present mesh:
<ul>
<li> <code>assemble_system ()</code> [via <code>solve_timestep ()</code>]:
  This first function is also the most interesting one. It assembles the
  linear system corresponding to the discretized version of equation
  <a href="#step_18.linear-system">[linear-system]</a>. This leads to a system matrix $A_{ij} = \sum_K
  A^K_{ij}$ built up of local contributions on each cell $K$ with entries
  @f[
    A^K_{ij} = (C \varepsilon(\varphi_j), \varepsilon(\varphi_i))_K;
  @f]
  In practice, $A^K$ is computed using numerical quadrature according to the
  formula
  @f[
    A^K_{ij} = \sum_q w_q [\varepsilon(\varphi_i(\mathbf{x}_q)) : C :
                           \varepsilon(\varphi_j(\mathbf{x}_q))],
  @f]
  with quadrature points $\mathbf{x}_q$ and weights $w_q$. We have built these
  contributions before, in step-8 and step-17, but in both of these cases we
  have done so rather clumsily by using knowledge of how the rank-4 tensor $C$
  is composed, and considering individual elements of the strain tensors
  $\varepsilon(\varphi_i),\varepsilon(\varphi_j)$. This is not really
  convenient, in particular if we want to consider more complicated elasticity
  models than the isotropic case for which $C$ had the convenient form
  $C_{ijkl}  = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik} \delta_{jl}
  + \delta_{il} \delta_{jk})$. While we in fact do not use a more complicated
  form than this in the present program, we nevertheless want to write it in a
  way that would easily allow for this. It is then natural to introduce
  classes that represent symmetric tensors of rank 2 (for the strains and
  stresses) and 4 (for the stress-strain tensor $C$). Fortunately, deal.II
  provides these: the <code>SymmetricTensor<rank,dim></code> class template
  provides a full-fledged implementation of such tensors of rank <code>rank</code>
  (which needs to be an even number) and dimension <code>dim</code>.

  What we then need is two things: a way to create the stress-strain rank-4
  tensor $C$ as well as to create a symmetric tensor of rank 2 (the strain
  tensor) from the gradients of a shape function $\varphi_i$ at a quadrature
  point $\mathbf{x}_q$ on a given cell. At the top of the implementation of this
  example program, you will find such functions. The first one,
  <code>get_stress_strain_tensor</code>, takes two arguments corresponding to
  the Lam&eacute; constants $\lambda$ and $\mu$ and returns the stress-strain tensor
  for the isotropic case corresponding to these constants (in the program, we
  will choose constants corresponding to steel); it would be simple to replace
  this function by one that computes this tensor for the anisotropic case, or
  taking into account crystal symmetries, for example. The second one,
  <code>get_strain</code> takes an object of type <code>FEValues</code> and indices
  $i$ and $q$ and returns the symmetric gradient, i.e. the strain,
  corresponding to shape function $\varphi_i(\mathbf{x}_q)$, evaluated on the cell
  on which the <code>FEValues</code> object was last reinitialized.

  Given this, the innermost loop of <code>assemble_system</code> computes the
  local contributions to the matrix in the following elegant way (the variable
  <code>stress_strain_tensor</code>, corresponding to the tensor $C$, has
  previously been initialized with the result of the first function above):
  @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  for (unsigned int j=0; j<dofs_per_cell; ++j)
    for (unsigned int q_point=0; q_point<n_q_points;
         ++q_point)
      {
        const SymmetricTensor<2,dim>
          eps_phi_i = get_strain (fe_values, i, q_point),
          eps_phi_j = get_strain (fe_values, j, q_point);

        cell_matrix(i,j)
          += (eps_phi_i * stress_strain_tensor * eps_phi_j *
              fe_values.JxW (q_point));
      }
  @endcode
  It is worth noting the expressive power of this piece of code, and to
  compare it with the complications we had to go through in previous examples
  for the elasticity problem. (To be fair, the SymmetricTensor class
  template did not exist when these previous examples were written.) For
  simplicity, <code>operator*</code> provides for the (double summation) product
  between symmetric tensors of even rank here.

  Assembling the local contributions
  @f{eqnarray*}
      f^K_i &=&
      (\mathbf{f}, \varphi_i)_K -(\sigma^{n-1},\varepsilon(\varphi_i))_K
      \\
      &\approx&
      \sum_q
      w_q \left\{
        \mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q) -
        \sigma^{n-1}_q : \varepsilon(\varphi_i(\mathbf{x}_q))
      \right\}
  @f}
  to the right hand side of <a href="#step_18.linear-system">[linear-system]</a> is equally
  straightforward (note that we do not consider any boundary tractions $\mathbf{b}$ here). Remember that we only had to store the old stress in the
  quadrature points of cells. In the program, we will provide a variable
  <code>local_quadrature_points_data</code> that allows to access the stress
  $\sigma^{n-1}_q$ in each quadrature point. With this the code for the right
  hand side looks as this, again rather elegant:
  @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  {
    const unsigned int
      component_i = fe.system_to_component_index(i).first;

    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        const SymmetricTensor<2,dim> &old_stress
          = local_quadrature_points_data[q_point].old_stress;

        cell_rhs(i) += (body_force_values[q_point](component_i) *
                        fe_values.shape_value (i,q_point)
                        -
                        old_stress *
                        get_strain (fe_values,i,q_point)) *
                       fe_values.JxW (q_point);
      }
  }
  @endcode
  Note that in the multiplication $\mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q)$, we have made use of the fact that for the chosen finite element, only
  one vector component (namely <code>component_i</code>) of $\varphi_i$ is
  nonzero, and that we therefore also have to consider only one component of
  $\mathbf{f}(\mathbf{x}_q)$.

  This essentially concludes the new material we present in this function. It
  later has to deal with boundary conditions as well as hanging node
  constraints, but this parallels what we had to do previously in other
  programs already.

<li> <code>solve_linear_problem ()</code> [via <code>solve_timestep ()</code>]:
  Unlike the previous one, this function is not really interesting, since it
  does what similar functions have done in all previous tutorial programs --
  solving the linear system using the CG method, using an incomplete LU
  decomposition as a preconditioner (in the %parallel case, it uses an ILU of
  each processor's block separately). It is virtually unchanged
  from step-17.

<li> <code>update_quadrature_point_history ()</code> [via
  <code>solve_timestep ()</code>]: Based on the displacement field $\Delta \mathbf{u}^n$ computed before, we update the stress values in all quadrature points
  according to <a href="#step_18.stress-update">[stress-update]</a> and <a href="#step_18.stress-update+rot">[stress-update+rot]</a>,
  including the rotation of the coordinate system.

<li> <code>move_mesh ()</code>: Given the solution computed before, in this
  function we deform the mesh by moving each vertex by the displacement vector
  field evaluated at this particular vertex.

<li> <code>output_results ()</code>: This function simply outputs the solution
  based on what we have said above, i.e. every processor computes output only
  for its own portion of the domain. In addition to the solution, we also compute the norm of
  the stress averaged over all the quadrature points on each cell.
</ul>

With this general structure of the code, we only have to define what case we
want to solve. For the present program, we have chosen to simulate the
quasistatic deformation of a vertical cylinder for which the bottom boundary
is fixed and the top boundary is pushed down at a prescribed vertical
velocity. However, the horizontal velocity of the top boundary is left
unspecified -- one can imagine this situation as a well-greased plate pushing
from the top onto the cylinder, the points on the top boundary of the cylinder
being allowed to slide horizontally along the surface of the plate, but forced
to move downward by the plate. The inner and outer boundaries of the cylinder
are free and not subject to any prescribed deflection or traction. In
addition, gravity acts on the body.

The program text will reveal more about how to implement this situation, and
the results section will show what displacement pattern comes out of this
simulation.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * First the usual list of header files that have already been used in
 * previous example programs:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/multithread_info.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/petsc_vector.h>
 * #include <deal.II/lac/petsc_sparse_matrix.h>
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/petsc_precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparsity_tools.h>
 * #include <deal.II/distributed/shared_tria.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * And here the only three new things among the header files: an include file in
 * which symmetric tensors of rank 2 and 4 are implemented, as introduced in
 * the introduction:
 * 
 * @code
 * #include <deal.II/base/symmetric_tensor.h>
 * 
 * @endcode
 * 
 * And lastly a header that contains some functions that will help us compute
 * rotaton matrices of the local coordinate systems at specific points in the
 * domain.
 * 
 * @code
 * #include <deal.II/physics/transformations.h>
 * 
 * @endcode
 * 
 * This is then simply C++ again:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * #include <iomanip>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step18
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePointHistorycodeclass"></a> 
 * <h3>The <code>PointHistory</code> class</h3>
 * 

 * 
 * As was mentioned in the introduction, we have to store the old stress in
 * quadrature point so that we can compute the residual forces at this point
 * during the next time step. This alone would not warrant a structure with
 * only one member, but in more complicated applications, we would have to
 * store more information in quadrature points as well, such as the history
 * variables of plasticity, etc. In essence, we have to store everything
 * that affects the present state of the material here, which in plasticity
 * is determined by the deformation history variables.
 *   

 * 
 * We will not give this class any meaningful functionality beyond being
 * able to store data, i.e. there are no constructors, destructors, or other
 * member functions. In such cases of `dumb' classes, we usually opt to
 * declare them as <code>struct</code> rather than <code>class</code>, to
 * indicate that they are closer to C-style structures than C++-style
 * classes.
 * 
 * @code
 *   template <int dim>
 *   struct PointHistory
 *   {
 *     SymmetricTensor<2, dim> old_stress;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thestressstraintensor"></a> 
 * <h3>The stress-strain tensor</h3>
 * 

 * 
 * Next, we define the linear relationship between the stress and the strain
 * in elasticity. It is given by a tensor of rank 4 that is usually written
 * in the form $C_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il}
 * \delta_{jk}) + \lambda \delta_{ij} \delta_{kl}$. This tensor maps
 * symmetric tensor of rank 2 to symmetric tensors of rank 2. A function
 * implementing its creation for given values of the Lam&eacute; constants
 * $\lambda$ and $\mu$ is straightforward:
 * 
 * @code
 *   template <int dim>
 *   SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
 *                                                    const double mu)
 *   {
 *     SymmetricTensor<4, dim> tmp;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = 0; j < dim; ++j)
 *         for (unsigned int k = 0; k < dim; ++k)
 *           for (unsigned int l = 0; l < dim; ++l)
 *             tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
 *                                ((i == l) && (j == k) ? mu : 0.0) +
 *                                ((i == j) && (k == l) ? lambda : 0.0));
 *     return tmp;
 *   }
 * 
 * @endcode
 * 
 * With this function, we will define a static member variable of the main
 * class below that will be used throughout the program as the stress-strain
 * tensor. Note that in more elaborate programs, this will probably be a
 * member variable of some class instead, or a function that returns the
 * stress-strain relationship depending on other input. For example in
 * damage theory models, the Lam&eacute; constants are considered a function
 * of the prior stress/strain history of a point. Conversely, in plasticity
 * the form of the stress-strain tensor is modified if the material has
 * reached the yield stress in a certain point, and possibly also depending on
 * its prior history.
 *   

 * 
 * In the present program, however, we assume that the material is
 * completely elastic and linear, and a constant stress-strain tensor is
 * sufficient for our present purposes.
 * 

 * 
 * 

 * 
 * 

 * 
 * 
 * <a name="Auxiliaryfunctions"></a> 
 * <h3>Auxiliary functions</h3>
 * 

 * 
 * Before the rest of the program, here are a few functions that we need as
 * tools. These are small functions that are called in inner loops, so we
 * mark them as <code>inline</code>.
 *   

 * 
 * The first one computes the symmetric strain tensor for shape function
 * <code>shape_func</code> at quadrature point <code>q_point</code> by
 * forming the symmetric gradient of this shape function. We need that when
 * we want to form the matrix, for example.
 *   

 * 
 * We should note that in previous examples where we have treated
 * vector-valued problems, we have always asked the finite element object in
 * which of the vector component the shape function is actually non-zero,
 * and thereby avoided to compute any terms that we could prove were zero
 * anyway. For this, we used the <code>fe.system_to_component_index</code>
 * function that returns in which component a shape function was zero, and
 * also that the <code>fe_values.shape_value</code> and
 * <code>fe_values.shape_grad</code> functions only returned the value and
 * gradient of the single non-zero component of a shape function if this is
 * a vector-valued element.
 *   

 * 
 * This was an optimization, and if it isn't terribly time critical, we can
 * get away with a simpler technique: just ask the <code>fe_values</code>
 * for the value or gradient of a given component of a given shape function
 * at a given quadrature point. This is what the
 * <code>fe_values.shape_grad_component(shape_func,q_point,i)</code> call
 * does: return the full gradient of the <code>i</code>th component of shape
 * function <code>shape_func</code> at quadrature point
 * <code>q_point</code>. If a certain component of a certain shape function
 * is always zero, then this will simply always return zero.
 *   

 * 
 * As mentioned, using <code>fe_values.shape_grad_component</code> instead
 * of the combination of <code>fe.system_to_component_index</code> and
 * <code>fe_values.shape_grad</code> may be less efficient, but its
 * implementation is optimized for such cases and shouldn't be a big
 * slowdown. We demonstrate the technique here since it is so much simpler
 * and straightforward.
 * 
 * @code
 *   template <int dim>
 *   inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values,
 *                                             const unsigned int   shape_func,
 *                                             const unsigned int   q_point)
 *   {
 * @endcode
 * 
 * Declare a temporary that will hold the return value:
 * 
 * @code
 *     SymmetricTensor<2, dim> tmp;
 * 
 * @endcode
 * 
 * First, fill diagonal terms which are simply the derivatives in
 * direction <code>i</code> of the <code>i</code> component of the
 * vector-valued shape function:
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i)
 *       tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
 * 
 * @endcode
 * 
 * Then fill the rest of the strain tensor. Note that since the tensor is
 * symmetric, we only have to compute one half (here: the upper right
 * corner) of the off-diagonal elements, and the implementation of the
 * <code>SymmetricTensor</code> class makes sure that at least to the
 * outside the symmetric entries are also filled (in practice, the class
 * of course stores only one copy). Here, we have picked the upper right
 * half of the tensor, but the lower left one would have been just as
 * good:
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = i + 1; j < dim; ++j)
 *         tmp[i][j] =
 *           (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
 *            fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
 *           2;
 * 
 *     return tmp;
 *   }
 * 
 * 
 * @endcode
 * 
 * The second function does something very similar (and therefore is given
 * the same name): compute the symmetric strain tensor from the gradient of
 * a vector-valued field. If you already have a solution field, the
 * <code>fe_values.get_function_gradients</code> function allows you to
 * extract the gradients of each component of your solution field at a
 * quadrature point. It returns this as a vector of rank-1 tensors: one rank-1
 * tensor (gradient) per vector component of the solution. From this we have
 * to reconstruct the (symmetric) strain tensor by transforming the data
 * storage format and symmetrization. We do this in the same way as above,
 * i.e. we avoid a few computations by filling first the diagonal and then
 * only one half of the symmetric tensor (the <code>SymmetricTensor</code>
 * class makes sure that it is sufficient to write only one of the two
 * symmetric components).
 *   

 * 
 * Before we do this, though, we make sure that the input has the kind of
 * structure we expect: that is that there are <code>dim</code> vector
 * components, i.e. one displacement component for each coordinate
 * direction. We test this with the <code>Assert</code> macro that will
 * simply abort our program if the condition is not met.
 * 
 * @code
 *   template <int dim>
 *   inline SymmetricTensor<2, dim>
 *   get_strain(const std::vector<Tensor<1, dim>> &grad)
 *   {
 *     Assert(grad.size() == dim, ExcInternalError());
 * 
 *     SymmetricTensor<2, dim> strain;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       strain[i][i] = grad[i][i];
 * 
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = i + 1; j < dim; ++j)
 *         strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
 * 
 *     return strain;
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally, below we will need a function that computes the rotation matrix
 * induced by a displacement at a given point. In fact, of course, the
 * displacement at a single point only has a direction and a magnitude, it
 * is the change in direction and magnitude that induces rotations. In
 * effect, the rotation matrix can be computed from the gradients of a
 * displacement, or, more specifically, from the curl.
 *   

 * 
 * The formulas by which the rotation matrices are determined are a little
 * awkward, especially in 3d. For 2d, there is a simpler way, so we
 * implement this function twice, once for 2d and once for 3d, so that we
 * can compile and use the program in both space dimensions if so desired --
 * after all, deal.II is all about dimension independent programming and
 * reuse of algorithm thoroughly tested with cheap computations in 2d, for
 * the more expensive computations in 3d. Here is one case, where we have to
 * implement different algorithms for 2d and 3d, but then can write the rest
 * of the program in a way that is independent of the space dimension.
 *   

 * 
 * So, without further ado to the 2d implementation:
 * 
 * @code
 *   Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>> &grad_u)
 *   {
 * @endcode
 * 
 * First, compute the curl of the velocity field from the gradients. Note
 * that we are in 2d, so the rotation is a scalar:
 * 
 * @code
 *     const double curl = (grad_u[1][0] - grad_u[0][1]);
 * 
 * @endcode
 * 
 * From this, compute the angle of rotation:
 * 
 * @code
 *     const double angle = std::atan(curl);
 * 
 * @endcode
 * 
 * And from this, build the antisymmetric rotation matrix. We want this
 * rotation matrix to represent the rotation of the local coordinate system
 * with respect to the global Cartesian basis, to we construct it with a
 * negative angle. The rotation matrix therefore represents the rotation
 * required to move from the local to the global coordinate system.
 * 
 * @code
 *     return Physics::Transformations::Rotations::rotation_matrix_2d(-angle);
 *   }
 * 
 * 
 * @endcode
 * 
 * The 3d case is a little more contrived:
 * 
 * @code
 *   Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u)
 *   {
 * @endcode
 * 
 * Again first compute the curl of the velocity field. This time, it is a
 * real vector:
 * 
 * @code
 *     const Point<3> curl(grad_u[2][1] - grad_u[1][2],
 *                         grad_u[0][2] - grad_u[2][0],
 *                         grad_u[1][0] - grad_u[0][1]);
 * 
 * @endcode
 * 
 * From this vector, using its magnitude, compute the tangent of the angle
 * of rotation, and from it the actual angle of rotation with respect to
 * the Cartesian basis:
 * 
 * @code
 *     const double tan_angle = std::sqrt(curl * curl);
 *     const double angle     = std::atan(tan_angle);
 * 
 * @endcode
 * 
 * Now, here's one problem: if the angle of rotation is too small, that
 * means that there is no rotation going on (for example a translational
 * motion). In that case, the rotation matrix is the identity matrix.
 *     

 * 
 * The reason why we stress that is that in this case we have that
 * <code>tan_angle==0</code>. Further down, we need to divide by that
 * number in the computation of the axis of rotation, and we would get
 * into trouble when dividing doing so. Therefore, let's shortcut this and
 * simply return the identity matrix if the angle of rotation is really
 * small:
 * 
 * @code
 *     if (std::abs(angle) < 1e-9)
 *       {
 *         static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
 *         static const Tensor<2, 3> rot(rotation);
 *         return rot;
 *       }
 * 
 * @endcode
 * 
 * Otherwise compute the real rotation matrix. For this, again we rely on
 * a predefined function to compute the rotation matrix of the local
 * coordinate system.
 * 
 * @code
 *     const Point<3> axis = curl / tan_angle;
 *     return Physics::Transformations::Rotations::rotation_matrix_3d(axis,
 *                                                                    -angle);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeTopLevelcodeclass"></a> 
 * <h3>The <code>TopLevel</code> class</h3>
 * 

 * 
 * This is the main class of the program. Since the namespace already
 * indicates what problem we are solving, let's call it by what it does: it
 * directs the flow of the program, i.e. it is the toplevel driver.
 *   

 * 
 * The member variables of this class are essentially as before, i.e. it has
 * to have a triangulation, a DoF handler and associated objects such as
 * constraints, variables that describe the linear system, etc. There are a
 * good number of more member functions now, which we will explain below.
 *   

 * 
 * The external interface of the class, however, is unchanged: it has a
 * public constructor and destructor, and it has a <code>run</code>
 * function that initiated all the work.
 * 
 * @code
 *   template <int dim>
 *   class TopLevel
 *   {
 *   public:
 *     TopLevel();
 *     ~TopLevel();
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * The private interface is more extensive than in step-17. First, we
 * obviously need functions that create the initial mesh, set up the
 * variables that describe the linear system on the present mesh
 * (i.e. matrices and vectors), and then functions that actually assemble
 * the system, direct what has to be solved in each time step, a function
 * that solves the linear system that arises in each timestep (and returns
 * the number of iterations it took), and finally output the solution
 * vector on the correct mesh:
 * 
 * @code
 *     void create_coarse_grid();
 * 
 *     void setup_system();
 * 
 *     void assemble_system();
 * 
 *     void solve_timestep();
 * 
 *     unsigned int solve_linear_problem();
 * 
 *     void output_results() const;
 * 
 * @endcode
 * 
 * All, except for the first two, of these functions are called in each
 * timestep. Since the first time step is a little special, we have
 * separate functions that describe what has to happen in a timestep: one
 * for the first, and one for all following timesteps:
 * 
 * @code
 *     void do_initial_timestep();
 * 
 *     void do_timestep();
 * 
 * @endcode
 * 
 * Then we need a whole bunch of functions that do various things. The
 * first one refines the initial grid: we start on the coarse grid with a
 * pristine state, solve the problem, then look at it and refine the mesh
 * accordingly, and start the same process over again, again with a
 * pristine state. Thus, refining the initial mesh is somewhat simpler
 * than refining a grid between two successive time steps, since it does
 * not involve transferring data from the old to the new triangulation, in
 * particular the history data that is stored in each quadrature point.
 * 
 * @code
 *     void refine_initial_grid();
 * 
 * @endcode
 * 
 * At the end of each time step, we want to move the mesh vertices around
 * according to the incremental displacement computed in this time
 * step. This is the function in which this is done:
 * 
 * @code
 *     void move_mesh();
 * 
 * @endcode
 * 
 * Next are two functions that handle the history variables stored in each
 * quadrature point. The first one is called before the first timestep to
 * set up a pristine state for the history variables. It only works on
 * those quadrature points on cells that belong to the present processor:
 * 
 * @code
 *     void setup_quadrature_point_history();
 * 
 * @endcode
 * 
 * The second one updates the history variables at the end of each
 * timestep:
 * 
 * @code
 *     void update_quadrature_point_history();
 * 
 * @endcode
 * 
 * This is the new shared Triangulation:
 * 
 * @code
 *     parallel::shared::Triangulation<dim> triangulation;
 * 
 *     FESystem<dim> fe;
 * 
 *     DoFHandler<dim> dof_handler;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 * @endcode
 * 
 * One difference of this program is that we declare the quadrature
 * formula in the class declaration. The reason is that in all the other
 * programs, it didn't do much harm if we had used different quadrature
 * formulas when computing the matrix and the right hand side, for
 * example. However, in the present case it does: we store information in
 * the quadrature points, so we have to make sure all parts of the program
 * agree on where they are and how many there are on each cell. Thus, let
 * us first declare the quadrature formula that will be used throughout...
 * 
 * @code
 *     const QGauss<dim> quadrature_formula;
 * 
 * @endcode
 * 
 * ... and then also have a vector of history objects, one per quadrature
 * point on those cells for which we are responsible (i.e. we don't store
 * history data for quadrature points on cells that are owned by other
 * processors).
 * Note that, instead of storing and managing this data ourself, we
 * could use the CellDataStorage class like is done in step-44. However,
 * for the purpose of demonstration, in this case we manage the storage
 * manually.
 * 
 * @code
 *     std::vector<PointHistory<dim>> quadrature_point_history;
 * 
 * @endcode
 * 
 * The way this object is accessed is through a <code>user pointer</code>
 * that each cell, face, or edge holds: it is a <code>void*</code> pointer
 * that can be used by application programs to associate arbitrary data to
 * cells, faces, or edges. What the program actually does with this data
 * is within its own responsibility, the library just allocates some space
 * for these pointers, and application programs can set and read the
 * pointers for each of these objects.
 * 

 * 
 * 

 * 
 * Further: we need the objects of linear systems to be solved,
 * i.e. matrix, right hand side vector, and the solution vector. Since we
 * anticipate solving big problems, we use the same types as in step-17,
 * i.e. distributed %parallel matrices and vectors built on top of the
 * PETSc library. Conveniently, they can also be used when running on only
 * a single machine, in which case this machine happens to be the only one
 * in our %parallel universe.
 *     

 * 
 * However, as a difference to step-17, we do not store the solution
 * vector -- which here is the incremental displacements computed in each
 * time step -- in a distributed fashion. I.e., of course it must be a
 * distributed vector when computing it, but immediately after that we
 * make sure each processor has a complete copy. The reason is that we had
 * already seen in step-17 that many functions needed a complete
 * copy. While it is not hard to get it, this requires communication on
 * the network, and is thus slow. In addition, these were repeatedly the
 * same operations, which is certainly undesirable unless the gains of not
 * always having to store the entire vector outweighs it. When writing
 * this program, it turned out that we need a complete copy of the
 * solution in so many places that it did not seem worthwhile to only get
 * it when necessary. Instead, we opted to obtain the complete copy once
 * and for all, and instead get rid of the distributed copy
 * immediately. Thus, note that the declaration of
 * <code>incremental_displacement</code> does not denote a distribute
 * vector as would be indicated by the middle namespace <code>MPI</code>:
 * 
 * @code
 *     PETScWrappers::MPI::SparseMatrix system_matrix;
 * 
 *     PETScWrappers::MPI::Vector system_rhs;
 * 
 *     Vector<double> incremental_displacement;
 * 
 * @endcode
 * 
 * The next block of variables is then related to the time dependent
 * nature of the problem: they denote the length of the time interval
 * which we want to simulate, the present time and number of time step,
 * and length of present timestep:
 * 
 * @code
 *     double       present_time;
 *     double       present_timestep;
 *     double       end_time;
 *     unsigned int timestep_no;
 * 
 * @endcode
 * 
 * Then a few variables that have to do with %parallel processing: first,
 * a variable denoting the MPI communicator we use, and then two numbers
 * telling us how many participating processors there are, and where in
 * this world we are. Finally, a stream object that makes sure only one
 * processor is actually generating output to the console. This is all the
 * same as in step-17:
 * 
 * @code
 *     MPI_Comm mpi_communicator;
 * 
 *     const unsigned int n_mpi_processes;
 * 
 *     const unsigned int this_mpi_process;
 * 
 *     ConditionalOStream pcout;
 * 
 * @endcode
 * 
 * We are storing the locally owned and the locally relevant indices:
 * 
 * @code
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 * @endcode
 * 
 * Finally, we have a static variable that denotes the linear relationship
 * between the stress and strain. Since it is a constant object that does
 * not depend on any input (at least not in this program), we make it a
 * static variable and will initialize it in the same place where we
 * define the constructor of this class:
 * 
 * @code
 *     static const SymmetricTensor<4, dim> stress_strain_tensor;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBodyForcecodeclass"></a> 
 * <h3>The <code>BodyForce</code> class</h3>
 * 

 * 
 * Before we go on to the main functionality of this program, we have to
 * define what forces will act on the body whose deformation we want to
 * study. These may either be body forces or boundary forces. Body forces
 * are generally mediated by one of the four basic physical types of forces:
 * gravity, strong and weak interaction, and electromagnetism. Unless one
 * wants to consider subatomic objects (for which quasistatic deformation is
 * irrelevant and an inappropriate description anyway), only gravity and
 * electromagnetic forces need to be considered. Let us, for simplicity
 * assume that our body has a certain mass density, but is either
 * non-magnetic and not electrically conducting or that there are no
 * significant electromagnetic fields around. In that case, the body forces
 * are simply <code>rho g</code>, where <code>rho</code> is the material
 * density and <code>g</code> is a vector in negative z-direction with
 * magnitude 9.81 m/s^2.  Both the density and <code>g</code> are defined in
 * the function, and we take as the density 7700 kg/m^3, a value commonly
 * assumed for steel.
 *   

 * 
 * To be a little more general and to be able to do computations in 2d as
 * well, we realize that the body force is always a function returning a
 * <code>dim</code> dimensional vector. We assume that gravity acts along
 * the negative direction of the last, i.e. <code>dim-1</code>th
 * coordinate. The rest of the implementation of this function should be
 * mostly self-explanatory given similar definitions in previous example
 * programs. Note that the body force is independent of the location; to
 * avoid compiler warnings about unused function arguments, we therefore
 * comment out the name of the first argument of the
 * <code>vector_value</code> function:
 * 
 * @code
 *   template <int dim>
 *   class BodyForce : public Function<dim>
 *   {
 *   public:
 *     BodyForce();
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override;
 * 
 *     virtual void
 *     vector_value_list(const std::vector<Point<dim>> &points,
 *                       std::vector<Vector<double>> &  value_list) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   BodyForce<dim>::BodyForce()
 *     : Function<dim>(dim)
 *   {}
 * 
 * 
 *   template <int dim>
 *   inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/,
 *                                            Vector<double> &values) const
 *   {
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
 * 
 *     const double g   = 9.81;
 *     const double rho = 7700;
 * 
 *     values          = 0;
 *     values(dim - 1) = -rho * g;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BodyForce<dim>::vector_value_list(
 *     const std::vector<Point<dim>> &points,
 *     std::vector<Vector<double>> &  value_list) const
 *   {
 *     const unsigned int n_points = points.size();
 * 
 *     Assert(value_list.size() == n_points,
 *            ExcDimensionMismatch(value_list.size(), n_points));
 * 
 *     for (unsigned int p = 0; p < n_points; ++p)
 *       BodyForce<dim>::vector_value(points[p], value_list[p]);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeIncrementalBoundaryValuecodeclass"></a> 
 * <h3>The <code>IncrementalBoundaryValue</code> class</h3>
 * 

 * 
 * In addition to body forces, movement can be induced by boundary forces
 * and forced boundary displacement. The latter case is equivalent to forces
 * being chosen in such a way that they induce certain displacement.
 *   

 * 
 * For quasistatic displacement, typical boundary forces would be pressure
 * on a body, or tangential friction against another body. We chose a
 * somewhat simpler case here: we prescribe a certain movement of (parts of)
 * the boundary, or at least of certain components of the displacement
 * vector. We describe this by another vector-valued function that, for a
 * given point on the boundary, returns the prescribed displacement.
 *   

 * 
 * Since we have a time-dependent problem, the displacement increment of the
 * boundary equals the displacement accumulated during the length of the
 * timestep. The class therefore has to know both the present time and the
 * length of the present time step, and can then approximate the incremental
 * displacement as the present velocity times the present timestep.
 *   

 * 
 * For the purposes of this program, we choose a simple form of boundary
 * displacement: we displace the top boundary with constant velocity
 * downwards. The rest of the boundary is either going to be fixed (and is
 * then described using an object of type
 * <code>Functions::ZeroFunction</code>) or free (Neumann-type, in which case
 * nothing special has to be done).  The implementation of the class
 * describing the constant downward motion should then be obvious using the
 * knowledge we gained through all the previous example programs:
 * 
 * @code
 *   template <int dim>
 *   class IncrementalBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     IncrementalBoundaryValues(const double present_time,
 *                               const double present_timestep);
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override;
 * 
 *     virtual void
 *     vector_value_list(const std::vector<Point<dim>> &points,
 *                       std::vector<Vector<double>> &  value_list) const override;
 * 
 *   private:
 *     const double velocity;
 *     const double present_time;
 *     const double present_timestep;
 *   };
 * 
 * 
 *   template <int dim>
 *   IncrementalBoundaryValues<dim>::IncrementalBoundaryValues(
 *     const double present_time,
 *     const double present_timestep)
 *     : Function<dim>(dim)
 *     , velocity(.08)
 *     , present_time(present_time)
 *     , present_timestep(present_timestep)
 *   {}
 * 
 * 
 *   template <int dim>
 *   void
 *   IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/,
 *                                                Vector<double> &values) const
 *   {
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
 * 
 *     values    = 0;
 *     values(2) = -present_timestep * velocity;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void IncrementalBoundaryValues<dim>::vector_value_list(
 *     const std::vector<Point<dim>> &points,
 *     std::vector<Vector<double>> &  value_list) const
 *   {
 *     const unsigned int n_points = points.size();
 * 
 *     Assert(value_list.size() == n_points,
 *            ExcDimensionMismatch(value_list.size(), n_points));
 * 
 *     for (unsigned int p = 0; p < n_points; ++p)
 *       IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeTopLevelcodeclass"></a> 
 * <h3>Implementation of the <code>TopLevel</code> class</h3>
 * 

 * 
 * Now for the implementation of the main class. First, we initialize the
 * stress-strain tensor, which we have declared as a static const
 * variable. We chose Lam&eacute; constants that are appropriate for steel:
 * 
 * @code
 *   template <int dim>
 *   const SymmetricTensor<4, dim> TopLevel<dim>::stress_strain_tensor =
 *     get_stress_strain_tensor<dim>(/*lambda = */ 9.695e10,
 *                                   /*mu     = */ 7.617e10);
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thepublicinterface"></a> 
 * <h4>The public interface</h4>
 * 

 * 
 * The next step is the definition of constructors and destructors. There
 * are no surprises here: we choose linear and continuous finite elements
 * for each of the <code>dim</code> vector components of the solution, and a
 * Gaussian quadrature formula with 2 points in each coordinate
 * direction. The destructor should be obvious:
 * 
 * @code
 *   template <int dim>
 *   TopLevel<dim>::TopLevel()
 *     : triangulation(MPI_COMM_WORLD)
 *     , fe(FE_Q<dim>(1), dim)
 *     , dof_handler(triangulation)
 *     , quadrature_formula(fe.degree + 1)
 *     , present_time(0.0)
 *     , present_timestep(1.0)
 *     , end_time(10.0)
 *     , timestep_no(0)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
 *     , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
 *     , pcout(std::cout, this_mpi_process == 0)
 *   {}
 * 
 * 
 * 
 *   template <int dim>
 *   TopLevel<dim>::~TopLevel()
 *   {
 *     dof_handler.clear();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The last of the public functions is the one that directs all the work,
 * <code>run()</code>. It initializes the variables that describe where in
 * time we presently are, then runs the first time step, then loops over all
 * the other time steps. Note that for simplicity we use a fixed time step,
 * whereas a more sophisticated program would of course have to choose it in
 * some more reasonable way adaptively:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::run()
 *   {
 *     do_initial_timestep();
 * 
 *     while (present_time < end_time)
 *       do_timestep();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelcreate_coarse_grid"></a> 
 * <h4>TopLevel::create_coarse_grid</h4>
 * 

 * 
 * The next function in the order in which they were declared above is the
 * one that creates the coarse grid from which we start. For this example
 * program, we want to compute the deformation of a cylinder under axial
 * compression. The first step therefore is to generate a mesh for a
 * cylinder of length 3 and with inner and outer radii of 0.8 and 1,
 * respectively. Fortunately, there is a library function for such a mesh.
 *   

 * 
 * In a second step, we have to associated boundary conditions with the
 * upper and lower faces of the cylinder. We choose a boundary indicator of
 * 0 for the boundary faces that are characterized by their midpoints having
 * z-coordinates of either 0 (bottom face), an indicator of 1 for z=3 (top
 * face); finally, we use boundary indicator 2 for all faces on the inside
 * of the cylinder shell, and 3 for the outside.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::create_coarse_grid()
 *   {
 *     const double inner_radius = 0.8, outer_radius = 1;
 *     GridGenerator::cylinder_shell(triangulation, 3, inner_radius, outer_radius);
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             const Point<dim> face_center = face->center();
 * 
 *             if (face_center[2] == 0)
 *               face->set_boundary_id(0);
 *             else if (face_center[2] == 3)
 *               face->set_boundary_id(1);
 *             else if (std::sqrt(face_center[0] * face_center[0] +
 *                                face_center[1] * face_center[1]) <
 *                      (inner_radius + outer_radius) / 2)
 *               face->set_boundary_id(2);
 *             else
 *               face->set_boundary_id(3);
 *           }
 * 
 * @endcode
 * 
 * Once all this is done, we can refine the mesh once globally:
 * 
 * @code
 *     triangulation.refine_global(1);
 * 
 * @endcode
 * 
 * As the final step, we need to set up a clean state of the data that we
 * store in the quadrature points on all cells that are treated on the
 * present processor.
 * 
 * @code
 *     setup_quadrature_point_history();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_system"></a> 
 * <h4>TopLevel::setup_system</h4>
 * 

 * 
 * The next function is the one that sets up the data structures for a given
 * mesh. This is done in most the same way as in step-17: distribute the
 * degrees of freedom, then sort these degrees of freedom in such a way that
 * each processor gets a contiguous chunk of them. Note that subdivisions into
 * chunks for each processor is handled in the functions that create or
 * refine grids, unlike in the previous example program (the point where
 * this happens is mostly a matter of taste; here, we chose to do it when
 * grids are created since in the <code>do_initial_timestep</code> and
 * <code>do_timestep</code> functions we want to output the number of cells
 * on each processor at a point where we haven't called the present function
 * yet).
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 * @endcode
 * 
 * The next step is to set up constraints due to hanging nodes. This has
 * been handled many times before:
 * 
 * @code
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 * @endcode
 * 
 * And then we have to set up the matrix. Here we deviate from step-17, in
 * which we simply used PETSc's ability to just know about the size of the
 * matrix and later allocate those nonzero elements that are being written
 * to. While this works just fine from a correctness viewpoint, it is not
 * at all efficient: if we don't give PETSc a clue as to which elements
 * are written to, it is (at least at the time of this writing) unbearably
 * slow when we set the elements in the matrix for the first time (i.e. in
 * the first timestep). Later on, when the elements have been allocated,
 * everything is much faster. In experiments we made, the first timestep
 * can be accelerated by almost two orders of magnitude if we instruct
 * PETSc which elements will be used and which are not.
 *     

 * 
 * To do so, we first generate the sparsity pattern of the matrix we are
 * going to work with, and make sure that the condensation of hanging node
 * constraints add the necessary additional entries in the sparsity
 * pattern:
 * 
 * @code
 *     DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs);
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     sparsity_pattern,
 *                                     hanging_node_constraints,
 *                                     /*keep constrained dofs*/ false);
 *     SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
 *                                                locally_owned_dofs,
 *                                                mpi_communicator,
 *                                                locally_relevant_dofs);
 * @endcode
 * 
 * Note that we have used the <code>DynamicSparsityPattern</code> class
 * here that was already introduced in step-11, rather than the
 * <code>SparsityPattern</code> class that we have used in all other
 * cases. The reason for this is that for the latter class to work we have
 * to give an initial upper bound for the number of entries in each row, a
 * task that is traditionally done by
 * <code>DoFHandler::max_couplings_between_dofs()</code>. However, this
 * function suffers from a serious problem: it has to compute an upper
 * bound to the number of nonzero entries in each row, and this is a
 * rather complicated task, in particular in 3d. In effect, while it is
 * quite accurate in 2d, it often comes up with much too large a number in
 * 3d, and in that case the <code>SparsityPattern</code> allocates much
 * too much memory at first, often several 100 MBs. This is later
 * corrected when <code>DoFTools::make_sparsity_pattern</code> is called
 * and we realize that we don't need all that much memory, but at time it
 * is already too late: for large problems, the temporary allocation of
 * too much memory can lead to out-of-memory situations.
 *     

 * 
 * In order to avoid this, we resort to the
 * <code>DynamicSparsityPattern</code> class that is slower but does
 * not require any up-front estimate on the number of nonzero entries per
 * row. It therefore only ever allocates as much memory as it needs at any
 * given time, and we can build it even for large 3d problems.
 *     

 * 
 * It is also worth noting that due to the specifics of
 * parallel::shared::Triangulation, the sparsity pattern we construct is
 * global, i.e. comprises all degrees of freedom whether they will be
 * owned by the processor we are on or another one (in case this program
 * is run in %parallel via MPI). This of course is not optimal -- it
 * limits the size of the problems we can solve, since storing the entire
 * sparsity pattern (even if only for a short time) on each processor does
 * not scale well. However, there are several more places in the program
 * in which we do this, for example we always keep the global
 * triangulation and DoF handler objects around, even if we only work on
 * part of them. At present, deal.II does not have the necessary
 * facilities to completely distribute these objects (a task that, indeed,
 * is very hard to achieve with adaptive meshes, since well-balanced
 * subdivisions of a domain tend to become unbalanced as the mesh is
 * adaptively refined).
 *     

 * 
 * With this data structure, we can then go to the PETSc sparse matrix and
 * tell it to preallocate all the entries we will later want to write to:
 * 
 * @code
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          sparsity_pattern,
 *                          mpi_communicator);
 * @endcode
 * 
 * After this point, no further explicit knowledge of the sparsity pattern
 * is required any more and we can let the <code>sparsity_pattern</code>
 * variable go out of scope without any problem.
 * 

 * 
 * The last task in this function is then only to reset the right hand
 * side vector as well as the solution vector to its correct size;
 * remember that the solution vector is a local one, unlike the right hand
 * side that is a distributed %parallel one and therefore needs to know
 * the MPI communicator over which it is supposed to transmit messages:
 * 
 * @code
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *     incremental_displacement.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelassemble_system"></a> 
 * <h4>TopLevel::assemble_system</h4>
 * 

 * 
 * Again, assembling the system matrix and right hand side follows the same
 * structure as in many example programs before. In particular, it is mostly
 * equivalent to step-17, except for the different right hand side that now
 * only has to take into account internal stresses. In addition, assembling
 * the matrix is made significantly more transparent by using the
 * <code>SymmetricTensor</code> class: note the elegance of forming the
 * scalar products of symmetric tensors of rank 2 and 4. The implementation
 * is also more general since it is independent of the fact that we may or
 * may not be using an isotropic elasticity tensor.
 *   

 * 
 * The first part of the assembly routine is as always:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::assemble_system()
 *   {
 *     system_rhs    = 0;
 *     system_matrix = 0;
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     BodyForce<dim>              body_force;
 *     std::vector<Vector<double>> body_force_values(n_q_points,
 *                                                   Vector<double>(dim));
 * 
 * @endcode
 * 
 * As in step-17, we only need to loop over all cells that belong to the
 * present processor:
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 *           fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Then loop over all indices i,j and quadrature points and assemble
 * the system matrix contributions from this cell.  Note how we
 * extract the symmetric gradients (strains) of the shape functions
 * at a given quadrature point from the <code>FEValues</code>
 * object, and the elegance with which we form the triple
 * contraction <code>eps_phi_i : C : eps_phi_j</code>; the latter
 * needs to be compared to the clumsy computations needed in
 * step-17, both in the introduction as well as in the respective
 * place in the program:
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *                 {
 *                   const SymmetricTensor<2, dim>
 *                     eps_phi_i = get_strain(fe_values, i, q_point),
 *                     eps_phi_j = get_strain(fe_values, j, q_point);
 * 
 *                   cell_matrix(i, j) += (eps_phi_i *            
 *                                         stress_strain_tensor * 
 *                                         eps_phi_j              
 *                                         ) *                    
 *                                        fe_values.JxW(q_point); 
 *                 }
 * 
 * 
 * @endcode
 * 
 * Then also assemble the local right hand side contributions. For
 * this, we need to access the prior stress value in this quadrature
 * point. To get it, we use the user pointer of this cell that
 * points into the global array to the quadrature point data
 * corresponding to the first quadrature point of the present cell,
 * and then add an offset corresponding to the index of the
 * quadrature point we presently consider:
 * 
 * @code
 *           const PointHistory<dim> *local_quadrature_points_data =
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
 * @endcode
 * 
 * In addition, we need the values of the external body forces at
 * the quadrature points on this cell:
 * 
 * @code
 *           body_force.vector_value_list(fe_values.get_quadrature_points(),
 *                                        body_force_values);
 * @endcode
 * 
 * Then we can loop over all degrees of freedom on this cell and
 * compute local contributions to the right hand side:
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const unsigned int component_i =
 *                 fe.system_to_component_index(i).first;
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *                 {
 *                   const SymmetricTensor<2, dim> &old_stress =
 *                     local_quadrature_points_data[q_point].old_stress;
 * 
 *                   cell_rhs(i) +=
 *                     (body_force_values[q_point](component_i) *
 *                        fe_values.shape_value(i, q_point) -
 *                      old_stress * get_strain(fe_values, i, q_point)) *
 *                     fe_values.JxW(q_point);
 *                 }
 *             }
 * 
 * @endcode
 * 
 * Now that we have the local contributions to the linear system, we
 * need to transfer it into the global objects. This is done exactly
 * as in step-17:
 * 
 * @code
 *           cell->get_dof_indices(local_dof_indices);
 * 
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix,
 *                                                               cell_rhs,
 *                                                               local_dof_indices,
 *                                                               system_matrix,
 *                                                               system_rhs);
 *         }
 * 
 * @endcode
 * 
 * Now compress the vector and the system matrix:
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 * 
 * 
 * @endcode
 * 
 * The last step is to again fix up boundary values, just as we already
 * did in previous programs. A slight complication is that the
 * <code>apply_boundary_values</code> function wants to have a solution
 * vector compatible with the matrix and right hand side (i.e. here a
 * distributed %parallel vector, rather than the sequential vector we use
 * in this program) in order to preset the entries of the solution vector
 * with the correct boundary values. We provide such a compatible vector
 * in the form of a temporary vector which we then copy into the
 * sequential one.
 * 

 * 
 * We make up for this complication by showing how boundary values can be
 * used flexibly: following the way we create the triangulation, there are
 * three distinct boundary indicators used to describe the domain,
 * corresponding to the bottom and top faces, as well as the inner/outer
 * surfaces. We would like to impose boundary conditions of the following
 * type: The inner and outer cylinder surfaces are free of external
 * forces, a fact that corresponds to natural (Neumann-type) boundary
 * conditions for which we don't have to do anything. At the bottom, we
 * want no movement at all, corresponding to the cylinder being clamped or
 * cemented in at this part of the boundary. At the top, however, we want
 * a prescribed vertical downward motion compressing the cylinder; in
 * addition, we only want to restrict the vertical movement, but not the
 * horizontal ones -- one can think of this situation as a well-greased
 * plate sitting on top of the cylinder pushing it downwards: the atoms of
 * the cylinder are forced to move downward, but they are free to slide
 * horizontally along the plate.
 * 

 * 
 * The way to describe this is as follows: for boundary indicator zero
 * (bottom face) we use a dim-dimensional zero function representing no
 * motion in any coordinate direction. For the boundary with indicator 1
 * (top surface), we use the <code>IncrementalBoundaryValues</code> class,
 * but we specify an additional argument to the
 * <code>VectorTools::interpolate_boundary_values</code> function denoting
 * which vector components it should apply to; this is a vector of bools
 * for each vector component and because we only want to restrict vertical
 * motion, it has only its last component set:
 * 
 * @code
 *     FEValuesExtractors::Scalar                z_component(dim - 1);
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(dim),
 *                                              boundary_values);
 *     VectorTools::interpolate_boundary_values(
 *       dof_handler,
 *       1,
 *       IncrementalBoundaryValues<dim>(present_time, present_timestep),
 *       boundary_values,
 *       fe.component_mask(z_component));
 * 
 *     PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
 *     MatrixTools::apply_boundary_values(
 *       boundary_values, system_matrix, tmp, system_rhs, false);
 *     incremental_displacement = tmp;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_timestep"></a> 
 * <h4>TopLevel::solve_timestep</h4>
 * 

 * 
 * The next function is the one that controls what all has to happen within
 * a timestep. The order of things should be relatively self-explanatory
 * from the function names:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::solve_timestep()
 *   {
 *     pcout << "    Assembling system..." << std::flush;
 *     assemble_system();
 *     pcout << " norm of rhs is " << system_rhs.l2_norm() << std::endl;
 * 
 *     const unsigned int n_iterations = solve_linear_problem();
 * 
 *     pcout << "    Solver converged in " << n_iterations << " iterations."
 *           << std::endl;
 * 
 *     pcout << "    Updating quadrature point data..." << std::flush;
 *     update_quadrature_point_history();
 *     pcout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_linear_problem"></a> 
 * <h4>TopLevel::solve_linear_problem</h4>
 * 

 * 
 * Solving the linear system again works mostly as before. The only
 * difference is that we want to only keep a complete local copy of the
 * solution vector instead of the distributed one that we get as output from
 * PETSc's solver routines. To this end, we declare a local temporary
 * variable for the distributed vector and initialize it with the contents
 * of the local variable (remember that the
 * <code>apply_boundary_values</code> function called in
 * <code>assemble_system</code> preset the values of boundary nodes in this
 * vector), solve with it, and at the end of the function copy it again into
 * the complete local vector that we declared as a member variable. Hanging
 * node constraints are then distributed only on the local copy,
 * i.e. independently of each other on each of the processors:
 * 
 * @code
 *   template <int dim>
 *   unsigned int TopLevel<dim>::solve_linear_problem()
 *   {
 *     PETScWrappers::MPI::Vector distributed_incremental_displacement(
 *       locally_owned_dofs, mpi_communicator);
 *     distributed_incremental_displacement = incremental_displacement;
 * 
 *     SolverControl solver_control(dof_handler.n_dofs(),
 *                                  1e-16 * system_rhs.l2_norm());
 * 
 *     PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
 * 
 *     PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
 * 
 *     cg.solve(system_matrix,
 *              distributed_incremental_displacement,
 *              system_rhs,
 *              preconditioner);
 * 
 *     incremental_displacement = distributed_incremental_displacement;
 * 
 *     hanging_node_constraints.distribute(incremental_displacement);
 * 
 *     return solver_control.last_step();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveloutput_results"></a> 
 * <h4>TopLevel::output_results</h4>
 * 

 * 
 * This function generates the graphical output in .vtu format as explained
 * in the introduction. Each process will only work on the cells it owns,
 * and then write the result into a file of its own. Additionally, processor
 * 0 will write the record files the reference all the .vtu files.
 *   

 * 
 * The crucial part of this function is to give the <code>DataOut</code>
 * class a way to only work on the cells that the present process owns.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 * 
 * @endcode
 * 
 * Then, just as in step-17, define the names of solution variables (which
 * here are the displacement increments) and queue the solution vector for
 * output. Note in the following switch how we make sure that if the space
 * dimension should be unhandled that we throw an exception saying that we
 * haven't implemented this case yet (another case of defensive
 * programming):
 * 
 * @code
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 1:
 *           solution_names.emplace_back("delta_x");
 *           break;
 *         case 2:
 *           solution_names.emplace_back("delta_x");
 *           solution_names.emplace_back("delta_y");
 *           break;
 *         case 3:
 *           solution_names.emplace_back("delta_x");
 *           solution_names.emplace_back("delta_y");
 *           solution_names.emplace_back("delta_z");
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     data_out.add_data_vector(incremental_displacement, solution_names);
 * 
 * 
 * @endcode
 * 
 * The next thing is that we wanted to output something like the average
 * norm of the stresses that we have stored in each cell. This may seem
 * complicated, since on the present processor we only store the stresses
 * in quadrature points on those cells that actually belong to the present
 * process. In other words, it seems as if we can't compute the average
 * stresses for all cells. However, remember that our class derived from
 * <code>DataOut</code> only iterates over those cells that actually do
 * belong to the present processor, i.e. we don't have to compute anything
 * for all the other cells as this information would not be touched. The
 * following little loop does this. We enclose the entire block into a
 * pair of braces to make sure that the iterator variables do not remain
 * accidentally visible beyond the end of the block in which they are
 * used:
 * 
 * @code
 *     Vector<double> norm_of_stress(triangulation.n_active_cells());
 *     {
 * @endcode
 * 
 * Loop over all the cells...
 * 
 * @code
 *       for (auto &cell : triangulation.active_cell_iterators())
 *         if (cell->is_locally_owned())
 *           {
 * @endcode
 * 
 * On these cells, add up the stresses over all quadrature
 * points...
 * 
 * @code
 *             SymmetricTensor<2, dim> accumulated_stress;
 *             for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
 *               accumulated_stress +=
 *                 reinterpret_cast<PointHistory<dim> *>(cell->user_pointer())[q]
 *                   .old_stress;
 * 
 * @endcode
 * 
 * ...then write the norm of the average to their destination:
 * 
 * @code
 *             norm_of_stress(cell->active_cell_index()) =
 *               (accumulated_stress / quadrature_formula.size()).norm();
 *           }
 * @endcode
 * 
 * And on the cells that we are not interested in, set the respective
 * value in the vector to a bogus value (norms must be positive, and a
 * large negative value should catch your eye) in order to make sure
 * that if we were somehow wrong about our assumption that these
 * elements would not appear in the output file, that we would find out
 * by looking at the graphical output:
 * 
 * @code
 *         else
 *           norm_of_stress(cell->active_cell_index()) = -1e+20;
 *     }
 * @endcode
 * 
 * Finally attach this vector as well to be treated for output:
 * 
 * @code
 *     data_out.add_data_vector(norm_of_stress, "norm_of_stress");
 * 
 * @endcode
 * 
 * As a last piece of data, let us also add the partitioning of the domain
 * into subdomains associated with the processors if this is a parallel
 * job. This works in the exact same way as in the step-17 program:
 * 
 * @code
 *     std::vector<types::subdomain_id> partition_int(
 *       triangulation.n_active_cells());
 *     GridTools::get_subdomain_association(triangulation, partition_int);
 *     const Vector<double> partitioning(partition_int.begin(),
 *                                       partition_int.end());
 *     data_out.add_data_vector(partitioning, "partitioning");
 * 
 * @endcode
 * 
 * Finally, with all this data, we can instruct deal.II to munge the
 * information and produce some intermediate data structures that contain
 * all these solution and other data vectors:
 * 
 * @code
 *     data_out.build_patches();
 * 
 * @endcode
 * 
 * Let us call a function that opens the necessary output files and writes
 * the data we have generated into them. The function automatically
 * constructs the file names from the given directory name (the first
 * argument) and file name base (second argument). It augments the resulting
 * string by pieces that result from the time step number and a "piece
 * number" that corresponds to a part of the overall domain that can consist
 * of one or more subdomains.
 *     

 * 
 * The function also writes a record files (with suffix `.pvd`) for Paraview
 * that describes how all of these output files combine into the data for
 * this single time step:
 * 
 * @code
 *     const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", timestep_no, mpi_communicator, 4);
 * 
 * @endcode
 * 
 * The record files must be written only once and not by each processor,
 * so we do this on processor 0:
 * 
 * @code
 *     if (this_mpi_process == 0)
 *       {
 * @endcode
 * 
 * Finally, we write the paraview record, that references all .pvtu
 * files and their respective time. Note that the variable
 * times_and_names is declared static, so it will retain the entries
 * from the previous timesteps.
 * 
 * @code
 *         static std::vector<std::pair<double, std::string>> times_and_names;
 *         times_and_names.push_back(
 *           std::pair<double, std::string>(present_time, pvtu_filename));
 *         std::ofstream pvd_output("solution.pvd");
 *         DataOutBase::write_pvd_record(pvd_output, times_and_names);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_initial_timestep"></a> 
 * <h4>TopLevel::do_initial_timestep</h4>
 * 

 * 
 * This and the next function handle the overall structure of the first and
 * following timesteps, respectively. The first timestep is slightly more
 * involved because we want to compute it multiple times on successively
 * refined meshes, each time starting from a clean state. At the end of
 * these computations, in which we compute the incremental displacements
 * each time, we use the last results obtained for the incremental
 * displacements to compute the resulting stress updates and move the mesh
 * accordingly. On this new mesh, we then output the solution and any
 * additional data we consider important.
 *   

 * 
 * All this is interspersed by generating output to the console to update
 * the person watching the screen on what is going on. As in step-17, the
 * use of <code>pcout</code> instead of <code>std::cout</code> makes sure
 * that only one of the parallel processes is actually writing to the
 * console, without having to explicitly code an if-statement in each place
 * where we generate output:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::do_initial_timestep()
 *   {
 *     present_time += present_timestep;
 *     ++timestep_no;
 *     pcout << "Timestep " << timestep_no << " at time " << present_time
 *           << std::endl;
 * 
 *     for (unsigned int cycle = 0; cycle < 2; ++cycle)
 *       {
 *         pcout << "  Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           create_coarse_grid();
 *         else
 *           refine_initial_grid();
 * 
 *         pcout << "    Number of active cells:       "
 *               << triangulation.n_active_cells() << " (by partition:";
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p)
 *           pcout << (p == 0 ? ' ' : '+')
 *                 << (GridTools::count_cells_with_subdomain_association(
 *                      triangulation, p));
 *         pcout << ")" << std::endl;
 * 
 *         setup_system();
 * 
 *         pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (by partition:";
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p)
 *           pcout << (p == 0 ? ' ' : '+')
 *                 << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
 *                                                                     p));
 *         pcout << ")" << std::endl;
 * 
 *         solve_timestep();
 *       }
 * 
 *     move_mesh();
 *     output_results();
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_timestep"></a> 
 * <h4>TopLevel::do_timestep</h4>
 * 

 * 
 * Subsequent timesteps are simpler, and probably do not require any more
 * documentation given the explanations for the previous function above:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::do_timestep()
 *   {
 *     present_time += present_timestep;
 *     ++timestep_no;
 *     pcout << "Timestep " << timestep_no << " at time " << present_time
 *           << std::endl;
 *     if (present_time > end_time)
 *       {
 *         present_timestep -= (present_time - end_time);
 *         present_time = end_time;
 *       }
 * 
 * 
 *     solve_timestep();
 * 
 *     move_mesh();
 *     output_results();
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelrefine_initial_grid"></a> 
 * <h4>TopLevel::refine_initial_grid</h4>
 * 

 * 
 * The following function is called when solving the first time step on
 * successively refined meshes. After each iteration, it computes a
 * refinement criterion, refines the mesh, and sets up the history variables
 * in each quadrature point again to a clean state.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::refine_initial_grid()
 *   {
 * @endcode
 * 
 * First, let each process compute error indicators for the cells it owns:
 * 
 * @code
 *     Vector<float> error_per_cell(triangulation.n_active_cells());
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       incremental_displacement,
 *       error_per_cell,
 *       ComponentMask(),
 *       nullptr,
 *       MultithreadInfo::n_threads(),
 *       this_mpi_process);
 * 
 * @endcode
 * 
 * Then set up a global vector into which we merge the local indicators
 * from each of the %parallel processes:
 * 
 * @code
 *     const unsigned int n_local_cells =
 *       triangulation.n_locally_owned_active_cells();
 * 
 *     PETScWrappers::MPI::Vector distributed_error_per_cell(
 *       mpi_communicator, triangulation.n_active_cells(), n_local_cells);
 * 
 *     for (unsigned int i = 0; i < error_per_cell.size(); ++i)
 *       if (error_per_cell(i) != 0)
 *         distributed_error_per_cell(i) = error_per_cell(i);
 *     distributed_error_per_cell.compress(VectorOperation::insert);
 * 
 * @endcode
 * 
 * Once we have that, copy it back into local copies on all processors and
 * refine the mesh accordingly:
 * 
 * @code
 *     error_per_cell = distributed_error_per_cell;
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     error_per_cell,
 *                                                     0.35,
 *                                                     0.03);
 *     triangulation.execute_coarsening_and_refinement();
 * 
 * @endcode
 * 
 * Finally, set up quadrature point data again on the new mesh, and only
 * on those cells that we have determined to be ours:
 * 
 * @code
 *     setup_quadrature_point_history();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelmove_mesh"></a> 
 * <h4>TopLevel::move_mesh</h4>
 * 

 * 
 * At the end of each time step, we move the nodes of the mesh according to
 * the incremental displacements computed in this time step. To do this, we
 * keep a vector of flags that indicate for each vertex whether we have
 * already moved it around, and then loop over all cells and move those
 * vertices of the cell that have not been moved yet. It is worth noting
 * that it does not matter from which of the cells adjacent to a vertex we
 * move this vertex: since we compute the displacement using a continuous
 * finite element, the displacement field is continuous as well and we can
 * compute the displacement of a given vertex from each of the adjacent
 * cells. We only have to make sure that we move each node exactly once,
 * which is why we keep the vector of flags.
 *   

 * 
 * There are two noteworthy things in this function. First, how we get the
 * displacement field at a given vertex using the
 * <code>cell-@>vertex_dof_index(v,d)</code> function that returns the index
 * of the <code>d</code>th degree of freedom at vertex <code>v</code> of the
 * given cell. In the present case, displacement in the k-th coordinate
 * direction corresponds to the k-th component of the finite element. Using a
 * function like this bears a certain risk, because it uses knowledge of the
 * order of elements that we have taken together for this program in the
 * <code>FESystem</code> element. If we decided to add an additional
 * variable, for example a pressure variable for stabilization, and happened
 * to insert it as the first variable of the element, then the computation
 * below will start to produce nonsensical results. In addition, this
 * computation rests on other assumptions: first, that the element we use
 * has, indeed, degrees of freedom that are associated with vertices. This
 * is indeed the case for the present Q1 element, as would be for all Qp
 * elements of polynomial order <code>p</code>. However, it would not hold
 * for discontinuous elements, or elements for mixed formulations. Secondly,
 * it also rests on the assumption that the displacement at a vertex is
 * determined solely by the value of the degree of freedom associated with
 * this vertex; in other words, all shape functions corresponding to other
 * degrees of freedom are zero at this particular vertex. Again, this is the
 * case for the present element, but is not so for all elements that are
 * presently available in deal.II. Despite its risks, we choose to use this
 * way in order to present a way to query individual degrees of freedom
 * associated with vertices.
 *   

 * 
 * In this context, it is instructive to point out what a more general way
 * would be. For general finite elements, the way to go would be to take a
 * quadrature formula with the quadrature points in the vertices of a
 * cell. The <code>QTrapezoid</code> formula for the trapezoidal rule does
 * exactly this. With this quadrature formula, we would then initialize an
 * <code>FEValues</code> object in each cell, and use the
 * <code>FEValues::get_function_values</code> function to obtain the values
 * of the solution function in the quadrature points, i.e. the vertices of
 * the cell. These are the only values that we really need, i.e. we are not
 * at all interested in the weights (or the <code>JxW</code> values)
 * associated with this particular quadrature formula, and this can be
 * specified as the last argument in the constructor to
 * <code>FEValues</code>. The only point of minor inconvenience in this
 * scheme is that we have to figure out which quadrature point corresponds
 * to the vertex we consider at present, as they may or may not be ordered
 * in the same order.
 *   

 * 
 * This inconvenience could be avoided if finite elements have support
 * points on vertices (which the one here has; for the concept of support
 * points, see @ref GlossSupport "support points"). For such a case, one
 * could construct a custom quadrature rule using
 * FiniteElement::get_unit_support_points(). The first
 * <code>cell-&gt;n_vertices()*fe.dofs_per_vertex</code>
 * quadrature points will then correspond to the vertices of the cell and
 * are ordered consistent with <code>cell-@>vertex(i)</code>, taking into
 * account that support points for vector elements will be duplicated
 * <code>fe.dofs_per_vertex</code> times.
 *   

 * 
 * Another point worth explaining about this short function is the way in
 * which the triangulation class exports information about its vertices:
 * through the <code>Triangulation::n_vertices</code> function, it
 * advertises how many vertices there are in the triangulation. Not all of
 * them are actually in use all the time -- some are left-overs from cells
 * that have been coarsened previously and remain in existence since deal.II
 * never changes the number of a vertex once it has come into existence,
 * even if vertices with lower number go away. Secondly, the location
 * returned by <code>cell-@>vertex(v)</code> is not only a read-only object
 * of type <code>Point@<dim@></code>, but in fact a reference that can be
 * written to. This allows to move around the nodes of a mesh with relative
 * ease, but it is worth pointing out that it is the responsibility of an
 * application program using this feature to make sure that the resulting
 * cells are still useful, i.e. are not distorted so much that the cell is
 * degenerated (indicated, for example, by negative Jacobians). Note that we
 * do not have any provisions in this function to actually ensure this, we
 * just have faith.
 *   

 * 
 * After this lengthy introduction, here are the full 20 or so lines of
 * code:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::move_mesh()
 *   {
 *     pcout << "    Moving mesh..." << std::endl;
 * 
 *     std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
 *     for (auto &cell : dof_handler.active_cell_iterators())
 *       for (const auto v : cell->vertex_indices())
 *         if (vertex_touched[cell->vertex_index(v)] == false)
 *           {
 *             vertex_touched[cell->vertex_index(v)] = true;
 * 
 *             Point<dim> vertex_displacement;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               vertex_displacement[d] =
 *                 incremental_displacement(cell->vertex_dof_index(v, d));
 * 
 *             cell->vertex(v) += vertex_displacement;
 *           }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_quadrature_point_history"></a> 
 * <h4>TopLevel::setup_quadrature_point_history</h4>
 * 

 * 
 * At the beginning of our computations, we needed to set up initial values
 * of the history variables, such as the existing stresses in the material,
 * that we store in each quadrature point. As mentioned above, we use the
 * <code>user_pointer</code> for this that is available in each cell.
 *   

 * 
 * To put this into larger perspective, we note that if we had previously
 * available stresses in our model (which we assume do not exist for the
 * purpose of this program), then we would need to interpolate the field of
 * preexisting stresses to the quadrature points. Likewise, if we were to
 * simulate elasto-plastic materials with hardening/softening, then we would
 * have to store additional history variables like the present yield stress
 * of the accumulated plastic strains in each quadrature
 * points. Pre-existing hardening or weakening would then be implemented by
 * interpolating these variables in the present function as well.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::setup_quadrature_point_history()
 *   {
 * @endcode
 * 
 * For good measure, we set all user pointers of all cells, whether
 * ours of not, to the null pointer. This way, if we ever access the user
 * pointer of a cell which we should not have accessed, a segmentation
 * fault will let us know that this should not have happened:
 * 

 * 
 * 
 * @code
 *     triangulation.clear_user_data();
 * 
 * @endcode
 * 
 * Next, allocate the quadrature objects that are within the responsibility
 * of this processor. This, of course, equals the number of cells that
 * belong to this processor times the number of quadrature points our
 * quadrature formula has on each cell. Since the `resize()` function does
 * not actually shrink the amount of allocated memory if the requested new
 * size is smaller than the old size, we resort to a trick to first free all
 * memory, and then reallocate it: we declare an empty vector as a temporary
 * variable and then swap the contents of the old vector and this temporary
 * variable. This makes sure that the `quadrature_point_history` is now
 * really empty, and we can let the temporary variable that now holds the
 * previous contents of the vector go out of scope and be destroyed. In the
 * next step we can then re-allocate as many elements as we need, with the
 * vector default-initializing the `PointHistory` objects, which includes
 * setting the stress variables to zero.
 * 
 * @code
 *     {
 *       std::vector<PointHistory<dim>> tmp;
 *       quadrature_point_history.swap(tmp);
 *     }
 *     quadrature_point_history.resize(
 *       triangulation.n_locally_owned_active_cells() * quadrature_formula.size());
 * 
 * @endcode
 * 
 * Finally loop over all cells again and set the user pointers from the
 * cells that belong to the present processor to point to the first
 * quadrature point objects corresponding to this cell in the vector of
 * such objects:
 * 
 * @code
 *     unsigned int history_index = 0;
 *     for (auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell->set_user_pointer(&quadrature_point_history[history_index]);
 *           history_index += quadrature_formula.size();
 *         }
 * 
 * @endcode
 * 
 * At the end, for good measure make sure that our count of elements was
 * correct and that we have both used up all objects we allocated
 * previously, and not point to any objects beyond the end of the
 * vector. Such defensive programming strategies are always good checks to
 * avoid accidental errors and to guard against future changes to this
 * function that forget to update all uses of a variable at the same
 * time. Recall that constructs using the <code>Assert</code> macro are
 * optimized away in optimized mode, so do not affect the run time of
 * optimized runs:
 * 
 * @code
 *     Assert(history_index == quadrature_point_history.size(),
 *            ExcInternalError());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelupdate_quadrature_point_history"></a> 
 * <h4>TopLevel::update_quadrature_point_history</h4>
 * 

 * 
 * At the end of each time step, we should have computed an incremental
 * displacement update so that the material in its new configuration
 * accommodates for the difference between the external body and boundary
 * forces applied during this time step minus the forces exerted through
 * preexisting internal stresses. In order to have the preexisting
 * stresses available at the next time step, we therefore have to update the
 * preexisting stresses with the stresses due to the incremental
 * displacement computed during the present time step. Ideally, the
 * resulting sum of internal stresses would exactly counter all external
 * forces. Indeed, a simple experiment can make sure that this is so: if we
 * choose boundary conditions and body forces to be time independent, then
 * the forcing terms (the sum of external forces and internal stresses)
 * should be exactly zero. If you make this experiment, you will realize
 * from the output of the norm of the right hand side in each time step that
 * this is almost the case: it is not exactly zero, since in the first time
 * step the incremental displacement and stress updates were computed
 * relative to the undeformed mesh, which was then deformed. In the second
 * time step, we again compute displacement and stress updates, but this
 * time in the deformed mesh -- there, the resulting updates are very small
 * but not quite zero. This can be iterated, and in each such iteration the
 * residual, i.e. the norm of the right hand side vector, is reduced; if one
 * makes this little experiment, one realizes that the norm of this residual
 * decays exponentially with the number of iterations, and after an initial
 * very rapid decline is reduced by roughly a factor of about 3.5 in each
 * iteration (for one testcase I looked at, other testcases, and other
 * numbers of unknowns change the factor, but not the exponential decay).
 * 

 * 
 * In a sense, this can then be considered as a quasi-timestepping scheme to
 * resolve the nonlinear problem of solving large-deformation elasticity on
 * a mesh that is moved along in a Lagrangian manner.
 *   

 * 
 * Another complication is that the existing (old) stresses are defined on
 * the old mesh, which we will move around after updating the stresses. If
 * this mesh update involves rotations of the cell, then we need to also
 * rotate the updated stress, since it was computed relative to the
 * coordinate system of the old cell.
 *   

 * 
 * Thus, what we need is the following: on each cell which the present
 * processor owns, we need to extract the old stress from the data stored
 * with each quadrature point, compute the stress update, add the two
 * together, and then rotate the result together with the incremental
 * rotation computed from the incremental displacement at the present
 * quadrature point. We will detail these steps below:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::update_quadrature_point_history()
 *   {
 * @endcode
 * 
 * First, set up an <code>FEValues</code> object by which we will evaluate
 * the incremental displacements and the gradients thereof at the
 * quadrature points, together with a vector that will hold this
 * information:
 * 
 * @code
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients);
 * 
 *     std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads(
 *       quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim));
 * 
 * @endcode
 * 
 * Then loop over all cells and do the job in the cells that belong to our
 * subdomain:
 * 
 * @code
 *     for (auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 * @endcode
 * 
 * Next, get a pointer to the quadrature point history data local to
 * the present cell, and, as a defensive measure, make sure that
 * this pointer is within the bounds of the global array:
 * 
 * @code
 *           PointHistory<dim> *local_quadrature_points_history =
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
 *           Assert(local_quadrature_points_history >=
 *                    &quadrature_point_history.front(),
 *                  ExcInternalError());
 *           Assert(local_quadrature_points_history <=
 *                    &quadrature_point_history.back(),
 *                  ExcInternalError());
 * 
 * @endcode
 * 
 * Then initialize the <code>FEValues</code> object on the present
 * cell, and extract the gradients of the displacement at the
 * quadrature points for later computation of the strains
 * 
 * @code
 *           fe_values.reinit(cell);
 *           fe_values.get_function_gradients(incremental_displacement,
 *                                            displacement_increment_grads);
 * 
 * @endcode
 * 
 * Then loop over the quadrature points of this cell:
 * 
 * @code
 *           for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
 *             {
 * @endcode
 * 
 * On each quadrature point, compute the strain increment from
 * the gradients, and multiply it by the stress-strain tensor to
 * get the stress update. Then add this update to the already
 * existing strain at this point:
 * 
 * @code
 *               const SymmetricTensor<2, dim> new_stress =
 *                 (local_quadrature_points_history[q].old_stress +
 *                  (stress_strain_tensor *
 *                   get_strain(displacement_increment_grads[q])));
 * 
 * @endcode
 * 
 * Finally, we have to rotate the result. For this, we first
 * have to compute a rotation matrix at the present quadrature
 * point from the incremental displacements. In fact, it can be
 * computed from the gradients, and we already have a function
 * for that purpose:
 * 
 * @code
 *               const Tensor<2, dim> rotation =
 *                 get_rotation_matrix(displacement_increment_grads[q]);
 * @endcode
 * 
 * Note that the result, a rotation matrix, is in general an
 * antisymmetric tensor of rank 2, so we must store it as a full
 * tensor.
 * 

 * 
 * With this rotation matrix, we can compute the rotated tensor
 * by contraction from the left and right, after we expand the
 * symmetric tensor <code>new_stress</code> into a full tensor:
 * 
 * @code
 *               const SymmetricTensor<2, dim> rotated_new_stress =
 *                 symmetrize(transpose(rotation) *
 *                            static_cast<Tensor<2, dim>>(new_stress) * rotation);
 * @endcode
 * 
 * Note that while the result of the multiplication of these
 * three matrices should be symmetric, it is not due to floating
 * point round off: we get an asymmetry on the order of 1e-16 of
 * the off-diagonal elements of the result. When assigning the
 * result to a <code>SymmetricTensor</code>, the constructor of
 * that class checks the symmetry and realizes that it isn't
 * exactly symmetric; it will then raise an exception. To avoid
 * that, we explicitly symmetrize the result to make it exactly
 * symmetric.
 * 

 * 
 * The result of all these operations is then written back into
 * the original place:
 * 
 * @code
 *               local_quadrature_points_history[q].old_stress =
 *                 rotated_new_stress;
 *             }
 *         }
 *   }
 * 
 * @endcode
 * 
 * This ends the project specific namespace <code>Step18</code>. The rest is
 * as usual and as already shown in step-17: A <code>main()</code> function
 * that initializes and terminates PETSc, calls the classes that do the
 * actual work, and makes sure that we catch all exceptions that propagate
 * up to this point:
 * 
 * @code
 * } // namespace Step18
 * 
 * 
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step18;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       TopLevel<3> elastic_problem;
 *       elastic_problem.run();
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



Running the program takes a good while if one uses debug mode; it takes about
eleven minutes on my i7 desktop. Fortunately, the version compiled with
optimizations is much faster; the program only takes about a minute and a half
after recompiling with the command <tt>make release</tt> on the same machine, a
much more reasonable time.


If run, the program prints the following output, explaining what it is
doing during all that time:
@verbatim
\$ time make run
[ 66%] Built target \step-18
[100%] Run \step-18 with Release configuration
Timestep 1 at time 1
  Cycle 0:
    Number of active cells:       3712 (by partition: 3712)
    Number of degrees of freedom: 17226 (by partition: 17226)
    Assembling system... norm of rhs is 1.88062e+10
    Solver converged in 103 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       12812 (by partition: 12812)
    Number of degrees of freedom: 51738 (by partition: 51738)
    Assembling system... norm of rhs is 1.86145e+10
    Solver converged in 121 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 2 at time 2
    Assembling system... norm of rhs is 1.84169e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 3 at time 3
    Assembling system... norm of rhs is 1.82355e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 4 at time 4
    Assembling system... norm of rhs is 1.80728e+10
    Solver converged in 117 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 5 at time 5
    Assembling system... norm of rhs is 1.79318e+10
    Solver converged in 116 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 6 at time 6
    Assembling system... norm of rhs is 1.78171e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 7 at time 7
    Assembling system... norm of rhs is 1.7737e+10
    Solver converged in 112 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 8 at time 8
    Assembling system... norm of rhs is 1.77127e+10
    Solver converged in 111 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 9 at time 9
    Assembling system... norm of rhs is 1.78207e+10
    Solver converged in 113 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 10 at time 10
    Assembling system... norm of rhs is 1.83544e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...

[100%] Built target run
make run  176.82s user 0.15s system 198% cpu 1:28.94 total
@endverbatim
In other words, it is computing on 12,000 cells and with some 52,000
unknowns. Not a whole lot, but enough for a coupled three-dimensional
problem to keep a computer busy for a while. At the end of the day,
this is what we have for output:
@verbatim
\$ ls -l *vtu *visit
-rw-r--r-- 1 drwells users 1706059 Feb 13 19:36 solution-0010.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0010.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0010.visit
-rw-r--r-- 1 drwells users 1707907 Feb 13 19:36 solution-0009.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0009.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0009.visit
-rw-r--r-- 1 drwells users 1703771 Feb 13 19:35 solution-0008.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0008.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0008.visit
-rw-r--r-- 1 drwells users 1693671 Feb 13 19:35 solution-0007.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0007.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0007.visit
-rw-r--r-- 1 drwells users 1681847 Feb 13 19:35 solution-0006.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0006.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0006.visit
-rw-r--r-- 1 drwells users 1670115 Feb 13 19:35 solution-0005.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0005.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0005.visit
-rw-r--r-- 1 drwells users 1658559 Feb 13 19:35 solution-0004.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0004.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0004.visit
-rw-r--r-- 1 drwells users 1639983 Feb 13 19:35 solution-0003.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0003.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0003.visit
-rw-r--r-- 1 drwells users 1625851 Feb 13 19:35 solution-0002.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0002.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0002.visit
-rw-r--r-- 1 drwells users 1616035 Feb 13 19:34 solution-0001.000.vtu
-rw-r--r-- 1 drwells users     761 Feb 13 19:34 solution-0001.pvtu
-rw-r--r-- 1 drwells users      33 Feb 13 19:34 solution-0001.visit
@endverbatim


If we visualize these files with VisIt or Paraview, we get to see the full picture
of the disaster our forced compression wreaks on the cylinder (colors in the
images encode the norm of the stress in the material):


<div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0002.0000.png"
           alt="Time = 2"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 2
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0005.0000.png"
           alt="Time = 5"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 5
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0007.0000.png"
           alt="Time = 7"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 7
    </div>
  </div>
</div>


<div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0008.0000.png"
           alt="Time = 8"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 8
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0009.0000.png"
           alt="Time = 9"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 9
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0010.0000.png"
           alt="Time = 10"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 10
    </div>
  </div>
</div>


As is clearly visible, as we keep compressing the cylinder, it starts
to bow out near the fully constrained bottom surface and, after about eight
time units, buckle in an azimuthally symmetric manner.


Although the result appears plausible for the symmetric geometry and loading,
it is yet to be established whether or not the computation is fully converged.
In order to see whether it is, we ran the program again with one more global
refinement at the beginning and with the time step halved. This would have
taken a very long time on a single machine, so we used a proper workstation and
ran it on 16 processors in parallel. The beginning of the output now looks like
this:
@verbatim
Timestep 1 at time 0.5
  Cycle 0:
    Number of active cells:       29696 (by partition: 1808+1802+1894+1881+1870+1840+1884+1810+1876+1818+1870+1884+1854+1903+1816+1886)
    Number of degrees of freedom: 113100 (by partition: 6936+6930+7305+7116+7326+6869+7331+6786+7193+6829+7093+7162+6920+7280+6843+7181)
    Assembling system... norm of rhs is 1.10765e+10
    Solver converged in 209 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       102034 (by partition: 6387+6202+6421+6341+6408+6201+6428+6428+6385+6294+6506+6244+6417+6527+6299+6546)
    Number of degrees of freedom: 359337 (by partition: 23255+21308+24774+24019+22304+21415+22430+22184+22298+21796+22396+21592+22325+22553+21977+22711)
    Assembling system... norm of rhs is 1.35759e+10
    Solver converged in 268 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 2 at time 1
    Assembling system... norm of rhs is 1.34674e+10
    Solver converged in 267 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 3 at time 1.5
    Assembling system... norm of rhs is 1.33607e+10
    Solver converged in 265 iterations.
    Updating quadrature point data...
    Moving mesh...

Timestep 4 at time 2
    Assembling system... norm of rhs is 1.32558e+10
    Solver converged in 263 iterations.
    Updating quadrature point data...
    Moving mesh...

[...]

Timestep 20 at time 10
    Assembling system... norm of rhs is 1.47755e+10
    Solver converged in 425 iterations.
    Updating quadrature point data...
    Moving mesh...
@endverbatim
That's quite a good number of unknowns, given that we are in 3d. The output of
this program are 16 files for each time step:
@verbatim
\$ ls -l solution-0001*
-rw-r--r-- 1 wellsd2 user 761065 Feb 13 21:09 solution-0001.000.vtu
-rw-r--r-- 1 wellsd2 user 759277 Feb 13 21:09 solution-0001.001.vtu
-rw-r--r-- 1 wellsd2 user 761217 Feb 13 21:09 solution-0001.002.vtu
-rw-r--r-- 1 wellsd2 user 761605 Feb 13 21:09 solution-0001.003.vtu
-rw-r--r-- 1 wellsd2 user 756917 Feb 13 21:09 solution-0001.004.vtu
-rw-r--r-- 1 wellsd2 user 752669 Feb 13 21:09 solution-0001.005.vtu
-rw-r--r-- 1 wellsd2 user 735217 Feb 13 21:09 solution-0001.006.vtu
-rw-r--r-- 1 wellsd2 user 750065 Feb 13 21:09 solution-0001.007.vtu
-rw-r--r-- 1 wellsd2 user 760273 Feb 13 21:09 solution-0001.008.vtu
-rw-r--r-- 1 wellsd2 user 777265 Feb 13 21:09 solution-0001.009.vtu
-rw-r--r-- 1 wellsd2 user 772469 Feb 13 21:09 solution-0001.010.vtu
-rw-r--r-- 1 wellsd2 user 760833 Feb 13 21:09 solution-0001.011.vtu
-rw-r--r-- 1 wellsd2 user 782241 Feb 13 21:09 solution-0001.012.vtu
-rw-r--r-- 1 wellsd2 user 748905 Feb 13 21:09 solution-0001.013.vtu
-rw-r--r-- 1 wellsd2 user 738413 Feb 13 21:09 solution-0001.014.vtu
-rw-r--r-- 1 wellsd2 user 762133 Feb 13 21:09 solution-0001.015.vtu
-rw-r--r-- 1 wellsd2 user   1421 Feb 13 21:09 solution-0001.pvtu
-rw-r--r-- 1 wellsd2 user    364 Feb 13 21:09 solution-0001.visit
@endverbatim


Here are first the mesh on which we compute as well as the partitioning
for the 16 processors:


<div class="twocolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-000mesh.png"
           alt="Discretization"
           width="400">
    </div>
    <div class="text" align="center">
      Discretization
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.p.png"
           alt="Parallel partitioning"
           width="400">
    </div>
    <div class="text" align="center">
      Parallel partitioning
    </div>
  </div>
</div>


Finally, here is the same output as we have shown before for the much smaller
sequential case:

<div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.s.png"
           alt="Time = 2"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 2
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0005.s.png"
           alt="Time = 5"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 5
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0007.s.png"
           alt="Time = 7"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 7
    </div>
  </div>
</div>


<div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0008.s.png"
           alt="Time = 8"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 8
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0009.s.png"
           alt="Time = 9"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 9
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0010.s.png"
           alt="Time = 10"
           width="400">
    </div>
    <div class="text" align="center">
      Time = 10
    </div>
  </div>
</div>


As before, we observe that at high axial compression the cylinder begins
to buckle, but this time ultimately collapses on itself. In contrast to our
first run, towards the end of the simulation the deflection pattern becomes
nonsymmetric (the central bulge deflects laterally). The model clearly does not
provide for this (all our forces and boundary deflections are symmetric) but the
effect is probably physically correct anyway: in reality, small inhomogeneities
in the body's material properties would lead it to buckle to one side
to evade the forcing; in numerical simulations, small perturbations
such as numerical round-off or an inexact solution of a linear system
by an iterative solver could have the same effect. Another typical source for
asymmetries in adaptive computations is that only a certain fraction of cells
is refined in each step, which may lead to asymmetric meshes even if the
original coarse mesh was symmetric.


If one compares this with the previous run, the results both qualitatively
and quantitatively different. The previous computation was
therefore certainly not converged, though we can't say for sure anything about
the present one. One would need an even finer computation to find out. However,
the point may be moot: looking at the last picture in detail, it is pretty
obvious that not only is the linear small deformation model we chose completely
inadequate, but for a realistic simulation we would also need to make sure that
the body does not intersect itself during deformation (if we continued
compressing the cylinder we would observe some self-intersection).
Without such a formulation we cannot expect anything to make physical sense,
even if it produces nice pictures!


<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


The program as is does not really solve an equation that has many applications
in practice: quasi-static material deformation based on a purely elastic law
is almost boring. However, the program may serve as the starting point for
more interesting experiments, and that indeed was the initial motivation for
writing it. Here are some suggestions of what the program is missing and in
what direction it may be extended:

<a name="Plasticitymodels"></a><h5>Plasticity models</h5>


 The most obvious extension is to use a more
realistic material model for large-scale quasistatic deformation. The natural
choice for this would be plasticity, in which a nonlinear relationship between
stress and strain replaces equation <a href="#step_18.stress-strain">[stress-strain]</a>. Plasticity
models are usually rather complicated to program since the stress-strain
dependence is generally non-smooth. The material can be thought of being able
to withstand only a maximal stress (the yield stress) after which it starts to
&ldquo;flow&rdquo;. A mathematical description to this can be given in the form of a
variational inequality, which alternatively can be treated as minimizing the
elastic energy
@f[
  E(\mathbf{u}) =
  (\varepsilon(\mathbf{u}), C\varepsilon(\mathbf{u}))_{\Omega}
  - (\mathbf{f}, \mathbf{u})_{\Omega} - (\mathbf{b}, \mathbf{u})_{\Gamma_N},
@f]
subject to the constraint
@f[
  f(\sigma(\mathbf{u})) \le 0
@f]
on the stress. This extension makes the problem to be solved in each time step
nonlinear, so we need another loop within each time step.

Without going into further details of this model, we refer to the excellent
book by Simo and Hughes on &ldquo;Computational Inelasticity&rdquo; for a
comprehensive overview of computational strategies for solving plastic
models. Alternatively, a brief but concise description of an algorithm for
plasticity is given in an article by S. Commend, A. Truty, and Th. Zimmermann;
@cite CTZ04.


<a name="Stabilizationissues"></a><h5>Stabilization issues</h5>


The formulation we have chosen, i.e. using
piecewise (bi-, tri-)linear elements for all components of the displacement
vector, and treating the stress as a variable dependent on the displacement is
appropriate for most materials. However, this so-called displacement-based
formulation becomes unstable and exhibits spurious modes for incompressible or
nearly-incompressible materials. While fluids are usually not elastic (in most
cases, the stress depends on velocity gradients, not displacement gradients,
although there are exceptions such as electro-rheologic fluids), there are a
few solids that are nearly incompressible, for example rubber. Another case is
that many plasticity models ultimately let the material become incompressible,
although this is outside the scope of the present program.

Incompressibility is characterized by Poisson's ratio
@f[
  \nu = \frac{\lambda}{2(\lambda+\mu)},
@f]
where $\lambda,\mu$ are the Lam&eacute; constants of the material.
Physical constraints indicate that $-1\le \nu\le \frac 12$ (the condition
also follows from mathematical stability considerations). If $\nu$
approaches $\frac 12$, then the material becomes incompressible. In that
case, pure displacement-based formulations are no longer appropriate for the
solution of such problems, and stabilization techniques have to be employed
for a stable and accurate solution. The book and paper cited above give
indications as to how to do this, but there is also a large volume of
literature on this subject; a good start to get an overview of the topic can
be found in the references of the paper by H.-Y. Duan and Q. Lin; @cite DL05.


<a name="Refinementduringtimesteps"></a><h5>Refinement during timesteps</h5>


In the present form, the program
only refines the initial mesh a number of times, but then never again. For any
kind of realistic simulation, one would want to extend this so that the mesh
is refined and coarsened every few time steps instead. This is not hard to do,
in fact, but has been left for future tutorial programs or as an exercise, if
you wish.

The main complication one has to overcome is that one has to
transfer the data that is stored in the quadrature points of the cells of the
old mesh to the new mesh, preferably by some sort of projection scheme. The
general approach to this would go like this:

- At the beginning, the data is only available in the quadrature points of
  individual cells, not as a finite element field that is defined everywhere.

- So let us find a finite element field that <i>is</i> defined everywhere so
  that we can later interpolate it to the quadrature points of the new
  mesh. In general, it will be difficult to find a continuous finite element
  field that matches the values in the quadrature points exactly because the
  number of degrees of freedom of these fields does not match the number of
  quadrature points there are, and the nodal values of this global field will
  either be over- or underdetermined. But it is usually not very difficult to
  find a discontinuous field that matches the values in the quadrature points;
  for example, if you have a QGauss(2) quadrature formula (i.e. 4 points per
  cell in 2d, 8 points in 3d), then one would use a finite element of kind
  FE_DGQ(1), i.e. bi-/tri-linear functions as these have 4 degrees of freedom
  per cell in 2d and 8 in 3d.

- There are functions that can make this conversion from individual points to
  a global field simpler. The following piece of pseudo-code should help if
  you use a QGauss(2) quadrature formula. Note that the multiplication by the
  projection matrix below takes a vector of scalar components, i.e., we can only
  convert one set of scalars at a time from the quadrature points to the degrees
  of freedom and vice versa. So we need to store each component of stress separately,
  which requires <code>dim*dim</code> vectors. We'll store this set of vectors in a 2D array to
  make it easier to read off components in the same way you would the stress tensor.
  Thus, we'll loop over the components of stress on each cell and store
  these values in the global history field. (The prefix <code>history_</code>
  indicates that we work with quantities related to the history variables defined
  in the quadrature points.)
  @code
    FE_DGQ<dim>     history_fe (1);
    DoFHandler<dim> history_dof_handler (triangulation);
    history_dof_handler.distribute_dofs (history_fe);

    std::vector< std::vector< Vector<double> > >
                 history_field (dim, std::vector< Vector<double> >(dim)),
                 local_history_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_history_fe_values (dim, std::vector< Vector<double> >(dim));

    for (unsigned int i=0; i<dim; i++)
      for (unsigned int j=0; j<dim; j++)
      {
        history_field[i][j].reinit(history_dof_handler.n_dofs());
        local_history_values_at_qpoints[i][j].reinit(quadrature.size());
        local_history_fe_values[i][j].reinit(history_fe.n_dofs_per_cell());
      }

    FullMatrix<double> qpoint_to_dof_matrix (history_fe.dofs_per_cell,
                                             quadrature.size());
    FETools::compute_projection_from_quadrature_points_matrix
              (history_fe,
               quadrature, quadrature,
               qpoint_to_dof_matrix);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++dg_cell)
      {

        PointHistory<dim> *local_quadrature_points_history
          = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

        Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
                ExcInternalError());
        Assert (local_quadrature_points_history < &quadrature_point_history.back(),
                ExcInternalError());

        for (unsigned int i=0; i<dim; i++)
          for (unsigned int j=0; j<dim; j++)
          {
            for (unsigned int q=0; q<quadrature.size(); ++q)
              local_history_values_at_qpoints[i][j](q)
                = local_quadrature_points_history[q].old_stress[i][j];

            qpoint_to_dof_matrix.vmult (local_history_fe_values[i][j],
                                        local_history_values_at_qpoints[i][j]);

            dg_cell->set_dof_values (local_history_fe_values[i][j],
                                     history_field[i][j]);
          }
      }
  @endcode

- Now that we have a global field, we can refine the mesh and transfer the
  history_field vector as usual using the SolutionTransfer class. This will
  interpolate everything from the old to the new mesh.

- In a final step, we have to get the data back from the now interpolated
  global field to the quadrature points on the new mesh. The following code
  will do that:
  @code
    FullMatrix<double> dof_to_qpoint_matrix (quadrature.size(),
                                             history_fe.n_dofs_per_cell());
    FETools::compute_interpolation_to_quadrature_points_matrix
              (history_fe,
               quadrature,
               dof_to_qpoint_matrix);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();

    for (; cell != endc; ++cell, ++dg_cell)
    {
      PointHistory<dim> *local_quadrature_points_history
       = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

      Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
              ExcInternalError());
      Assert (local_quadrature_points_history < &quadrature_point_history.back(),
              ExcInternalError());

      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
        {
          dg_cell->get_dof_values (history_field[i][j],
                                   local_history_fe_values[i][j]);

          dof_to_qpoint_matrix.vmult (local_history_values_at_qpoints[i][j],
                                      local_history_fe_values[i][j]);

          for (unsigned int q=0; q<quadrature.size(); ++q)
            local_quadrature_points_history[q].old_stress[i][j]
              = local_history_values_at_qpoints[i][j](q);
      }
  @endcode

It becomes a bit more complicated once we run the program in parallel, since
then each process only stores this data for the cells it owned on the old
mesh. That said, using a parallel vector for <code>history_field</code> will
do the trick if you put a call to <code>compress</code> after the transfer
from quadrature points into the global vector.


<a name="Ensuringmeshregularity"></a><h5>Ensuring mesh regularity</h5>


At present, the program makes no attempt
to make sure that a cell, after moving its vertices at the end of the time
step, still has a valid geometry (i.e. that its Jacobian determinant is
positive and bounded away from zero everywhere). It is, in fact, not very hard
to set boundary values and forcing terms in such a way that one gets distorted
and inverted cells rather quickly. Certainly, in some cases of large
deformation, this is unavoidable with a mesh of finite mesh size, but in some
other cases this should be preventable by appropriate mesh refinement and/or a
reduction of the time step size. The program does not do that, but a more
sophisticated version definitely should employ some sort of heuristic defining
what amount of deformation of cells is acceptable, and what isn't.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-18.cc"
*/
