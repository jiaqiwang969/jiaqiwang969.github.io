/**
@page step_19 The step-19 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Spatialdiscretization">Spatial discretization</a>
        <li><a href="#Dealingwithparticlesprogrammatically">Dealing with particles programmatically</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Globaldefinitions">Global definitions</a>
        <li><a href="#Themainclass">The main class</a>
        <li><a href="#ThecodeCathodeRaySimulatorcodeclassimplementation">The <code>CathodeRaySimulator</code> class implementation</a>
      <ul>
        <li><a href="#ThecodeCathodeRaySimulatorcodeconstructor">The <code>CathodeRaySimulator</code> constructor</a>
        <li><a href="#ThecodeCathodeRaySimulatormake_gridcodefunction">The <code>CathodeRaySimulator::make_grid</code> function</a>
        <li><a href="#ThecodeCathodeRaySimulatorsetup_systemcodefunction">The <code>CathodeRaySimulator::setup_system</code> function</a>
        <li><a href="#ThecodeCathodeRaySimulatorassemble_systemcodefunction">The <code>CathodeRaySimulator::assemble_system</code> function</a>
        <li><a href="#CathodeRaySimulatorsolve">CathodeRaySimulator::solve</a>
        <li><a href="#CathodeRaySimulatorrefine_grid">CathodeRaySimulator::refine_grid</a>
        <li><a href="#CathodeRaySimulatorcreate_particles">CathodeRaySimulator::create_particles</a>
        <li><a href="#CathodeRaySimulatormove_particles">CathodeRaySimulator::move_particles</a>
        <li><a href="#CathodeRaySimulatortrack_lost_particle">CathodeRaySimulator::track_lost_particle</a>
        <li><a href="#CathodeRaySimulatorupdate_timestep_size">CathodeRaySimulator::update_timestep_size</a>
        <li><a href="#ThecodeCathodeRaySimulatoroutput_resultscodefunction">The <code>CathodeRaySimulator::output_results()</code> function</a>
        <li><a href="#CathodeRaySimulatorrun">CathodeRaySimulator::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Avoidingaperformancebottleneckwithparticles"> Avoiding a performance bottleneck with particles </a>
        <li><a href="#Morestatisticsaboutelectrons"> More statistics about electrons </a>
        <li><a href="#Abettersynchronizedvisualization"> A better-synchronized visualization </a>
        <li><a href="#Abettertimeintegrator"> A better time integrator </a>
        <li><a href="#Parallelization"> Parallelization </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly

<br>

<i>
This program was contributed by Wolfgang Bangerth, Rene Gassmoeller, and Peter Munch.

Wolfgang Bangerth acknowledges support through NSF
awards DMS-1821210, EAR-1550901, and OAC-1835673.
</i>

@note Support for particles exists in deal.II primarily due to the initial
  efforts of Rene Gassmoeller. Please acknowledge this work by citing
  the publication @cite GLHPW2018 if you use particle functionality in your
  own work.

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


The finite element method in general, and deal.II in particular, were invented
to solve partial differential equations -- in other words, to solve
[continuum mechanics](https://en.wikipedia.org/wiki/Continuum_mechanics) problems.
On the other hand, sometimes one wants to solve problems in which it is useful
to track individual objects ("particles") and how their positions evolve. If
this simply leads to a set of ordinary differential equations, for example
if you want to track the positions of the planets in the solar system over
time, then deal.II is clearly not your right tool. On the other hand, if
this evolution is due to the interaction with the solution of partial differential
equation, or if having a mesh to determine which particles interact
with others (such as in the
[smoothed particle hydrodynamics (SPH)](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
method), then deal.II has support for you.

The case we will consider here is how electrically charged particles move through
an electric field. As motivation, we will consider
[cathode rays](https://en.wikipedia.org/wiki/Cathode_ray): Electrons emitted by a
heated piece of metal that is negatively charged (the "cathode"), and that are
then accelerated by an electric field towards the positively charged electrode
(the "anode"). The anode is typically ring-shaped so that the majority of
electrons can fly through the hole in the form of an electron beam. In the olden
times, they might then have illuminated the screen of a TV built from a
[cathode ray tube](https://en.wikipedia.org/wiki/Cathode-ray_tube).
Today, instead, electron beams are useful in
[X-ray machines](https://en.wikipedia.org/wiki/X-ray_tube),
[electron beam lithography](https://en.wikipedia.org/wiki/Electron-beam_lithography),
[electron beam welding](https://en.wikipedia.org/wiki/Electron-beam_welding), and
a number of other areas.

The equations we will then consider are as follows: First, we need to describe
the electric field. This is most easily accomplished by noting that the electric
potential $V$ satisfied the equation
@f[
  -\epsilon_0 \Delta V = \rho
@f]
where $\epsilon_0$ is the dielectric constant of vacuum, and $\rho$ is the charge
density. This is augmented by boundary conditions that we will choose as follows:
@f{align*}{
  V &= -V_0 && \text{on}\; \Gamma_\text{cathode}\subset\partial\Omega \\
  V &= +V_0 && \text{on}\; \Gamma_\text{anode}\subset\partial\Omega \\
  \epsilon\frac{\partial V}{\partial n} &= 0
   && \text{on}\; \partial\Omega\setminus\Gamma_\text{cathode}\setminus\Gamma_\text{anode}.
@f}
In other words, we prescribe voltages $+V_0$ and $-V_0$ at the two electrodes
and insulating (Neumann) boundary conditions elsewhere. Since the dynamics of the
particles are purely due to the electric field $\mathbf E=\nabla V$, we could
as well have prescribed $2V_0$ and $0$ at the two electrodes -- all that matters
is the voltage difference at the two electrodes.

Given this electric potential $V$ and the electric field $\mathbf E=\nabla V$,
we can describe the trajectory of the $i$th particle using the differential
equation
@f[
  m {\ddot {\mathbf x}}_i = e\mathbf E,
@f]
where $m,e$ are the mass and electric charge of each particle. In practice, it
is convenient to write this as a system of first-order differential equations
in the position $\mathbf x$ and velocity $\mathbf v$:
@f{align*}{
  {\dot {\mathbf v}}_i &= \frac{e\mathbf E}{m}, \\
  {\dot {\mathbf x}}_i &= {\mathbf v}_i.
@f}
The deal.II class we will use to deal with particles, Particles::ParticleHandler,
stores particles in a way so that the position $\mathbf x_i$ is part of the
Particles::ParticleHandler data structures. (It stores particles sorted
by cell they are in, and consequently needs to know where each particle is.)
The velocity $\mathbf v_i$, on the other hand, is of no concern to
Particles::ParticleHandler and consequently we will store it as a
"property" of each particle that we will update in each time step. Properties
can also be used to store any other quantity we might care about each particle:
its charge, or if they were larger than just an electron, its color, mass,
attitude in space, chemical composition, etc.

There remain two things to discuss to complete the model:
Where particles start and what the charge density $\rho$ is.

First, historically, cathode rays used very large electric fields to pull
electrons out of the metal. This produces only a relatively small current. One
can do better by heating the cathode: a statistical fraction of electrons in that
case have enough thermal energy to leave the metal; the electric field then just
has to be strong enough to pull them away from the attraction of their host
body. We will model this in the following way: We will create a new particle if
(i) the electric field points away from the electrode, i.e., if
$\mathbf E \cdot \mathbf n < 0$ where $\mathbf n$ is the normal vector at a
face pointing out of the domain (into the electrode), and (ii) the electric
field exceeds a threshold value $|\mathbf E|\ge E_\text{threshold}$. This is
surely not a sufficiently accurate model for what really happens, but is good
enough for our current tutorial program.

Second, in principle we would have to model the charge density via
@f[
  \rho(\mathbf x) = \sum_i e\delta(\mathbf x-\mathbf x_i).
@f]

@note
The issue now is that in reality, a cathode ray tube in an old television
yields a current of somewhere around a few milli-Amperes. In the much higher
energy beams of particle accelerators, the current may only be a few
nano-Ampere. But an Ampere is $6\times 10^{18}$ electrons flowing per
second. Now, as you will see in the results section, we really only simulate
a few microseconds ($10^{-5}$ seconds), but that still results in very very
large numbers of electrons -- far more than we can hope to simulate
with a program as small as the current one. As a consequence, let us
presume that each particle represents $N$ electrons. Then the particle
mass and charge are also $Nm$ and $Ne$ and the equations we have to
solve are
@f[
  (Nm) {\ddot {\mathbf x}}_i = (Ne)\mathbf E,
@f]
which is of course exactly the same as above. On the other hand, the charge
density for these "clumps" of electrons is given by
@f[
  \rho(\mathbf x) = \sum_i (Ne)\delta(\mathbf x-\mathbf x_i).
@f]
It is this form that we will implement in the program, where $N$ is chosen
rather large in the program to ensure that the particles actually affect
the electric field. (This may not be realistic in practice: In most cases,
there are just not enough electrons to actually affect the overall
electric field. But realism is not our goal here.)


@note One may wonder why the equation for the electric field (or, rather,
the electric potential) has no time derivative whereas the equations for
the electron positions do. In essence, this is a modeling assumption: We
assume that the particles move so slowly that at any given time the
electric field is in equilibrium. This is saying, in other words, that
the velocity of the electrons is much less than the speed of light. In
yet other words, we can rephrase this in terms of the electrode voltage
$V_0$: Since every volt of electric potential accelerates electrons by
approximately 600 km/s (neglecting relativistic effects), requiring
$|\mathbf v_i\|\ll c$ is equivalent to saying that $2V_0 \ll 500 \text{V}$.
Under this assumption (and the assumption that the total number
of electrons is small), one can also neglect the creation of
magnetic fields by the moving charges, which would otherwise also affect
the movement of the electrons.


<a name="Timediscretization"></a><h3>Time discretization</h3>


The equations outlined above form a set of coupled differential equations.
Let us bring them all together in one place again to make that clear:
@f{align*}{
  -\epsilon_0 \Delta V &= \sum_i e\delta(\mathbf x-\mathbf x_i)
  \\
  {\dot {\mathbf x}}_i &= {\mathbf v}_i,
  \\
  {\dot {\mathbf v}}_i &= \frac{e\mathbf E}{m} = \frac{e\mathbf \nabla V}{m}.
@f}
Because of the awkward dependence of the electric potential on the
particle locations, we don't want to solve this as a coupled system
but instead use a decoupled approach where we first solve for the
potential in each time step and then the particle locations. (One
could also do it the other way around, of course.) This is very
much in the same spirit as we do in step-21, step-31, and step-32,
to name just a few, and can all be understood in the context of
the operator splitting methods discussed in step-58.

So, if we denote by an upper index $n$ the time step, and if we
use a simple time discretization for the ODE, then this means
that we have to solve the following set of equations in each time
step:
@f{align*}{
  -\epsilon_0 \Delta V^{(n)} &= \sum_i e\delta(\mathbf x-\mathbf x_i^{(n-1)})
  \\
  \frac{{\mathbf v}_i^{(n)}-{\mathbf v}_i^{(n-1)}}{\Delta t} &= \frac{e\nabla V^{(n)}}{m}
  \\
  \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} &= {\mathbf v}_i^{(n)}.
@f}
There are of course many better ways to do a time discretization (for
example the simple [leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration))
but this isn't the point of the tutorial program, and so we will be content
with what we have here. (We will comment on a piece of this puzzle in the
<a href="#extensions">possibilities for extensions</a> section of this program,
however.)

There remains the question of how we should choose the time step size $\Delta t$.
The limitation here is that the Particles::ParticleHandler class needs to
keep track of which cell each particle is in. This is particularly an issue if
we are running computations in parallel (say, in step-70) because in that case
every process only stores those cells it owns plus one layer of "ghost cells".
That's not relevant here, but in general we should make sure that over the
course of each time step, a particle moves only from one cell to any
of its immediate neighbors (face, edge, or vertex neighbors). If we can ensure
that, then Particles::ParticleHandler is guaranteed to be able to figure out
which cell a particle ends up in. To do this, a useful rule of thumb
is that we should choose the time step so that for all particles the expected
distance the particle moves by is less than one cell diameter:
@f[
  \Delta t \le \frac{h_i}{\|\mathbf v_i\|} \qquad\qquad \forall i,
@f]
or equivalently
@f[
  \Delta t \le \min_i \frac{h_i}{\|\mathbf v_i\|}.
@f]
Here, $h_i$ is the length of the shortest edge of the cell on which particle
$i$ is located -- in essence, a measure of the size of a cell.

On the other hand, a particle might already be at the boundary of one cell
and the neighboring cell might be once further refined. So then the time to
cross that *neighboring* cell would actually be half the amount above,
suggesting
@f[
  \Delta t \le \min_i \frac{\tfrac 12 h_i}{\|\mathbf v_i\|}.
@f]

But even that is not good enough: The formula above updates the particle
positions in each time using the formula
@f[
\frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} = {\mathbf v}_i^{(n)},
@f]
that is, using the *current* velocity ${\mathbf v}_i^{n}$. But we don't have
the current velocity yet at the time when we need to choose $\Delta t$ -- which
is after we have updated the potential $V^{(n)}$ but before we update the
velocity from ${\mathbf v}_i^{(n-1)}$ to ${\mathbf v}_i^{(n)}$. All we have is
${\mathbf v}_i^{(n-1)}$. So we need an additional safety factor for our final
choice:
@f[
  \Delta t^{(n)} =
  c_\text{safety} \min_i \frac{\tfrac 12 h_i}{\|\mathbf v_i^{(n-1)}\|}.
@f]
How large should $c_\text{safety}$ be? That depends on how much of underestimate
$\|\mathbf v_i^{(n-1)}\|$ might be compared to $\|\mathbf v_i^{(n)}\|$, and that
is actually quite easy to assess: A particle created in one time step with
zero velocity will roughly pick up equal velocity increments in each successive
time step if the electric field it encounters along the way were roughly
constant. So the maximal difference between $\|\mathbf v_i^{(n-1)}\|$ and
$\|\mathbf v_i^{(n)}\|$ would be a factor of two. As a consequence,
we will choose $c_\text{saftey}=0.5$.

There is only one other case we ought to consider: What happens in
the very first time step? There, any particles to be moved along have just
been created, but they have a zero velocity. So we don't know what
velocity we should choose for them. Of course, in all other time steps
there are also particles that have just been created, but in general,
the particles with the highest velocity limit the time step size and so the
newly created particles with their zero velocity don't matter. But if we *only*
have such particles?

In that case, we can use the following approximation: If a particle
starts at $\mathbf v^{(0)}=0$, then the update formula tells us that
@f[
  {\mathbf v}_i^{(1)} = \frac{e\nabla V^{(1)}}{m} \Delta t,
@f]
and consequently
@f[
    \frac{{\mathbf x}_i^{(1)}-{\mathbf x}_i^{(0)}}{\Delta t} = {\mathbf v}_i^{(1)},
@f]
which we can write as
@f[
    {\mathbf x}_i^{(1)} - {\mathbf x}_i^{(0)} = \frac{e\nabla V^{(1)}}{m} \Delta t^2.
@f]
Not wanting to move a particle by more than $\frac 12 h_i$ then implies that we should
choose the time step as
@f[
  \Delta t
  \le
  \min_i
  \sqrt{ \frac{h_i m}{e \|\nabla V^{(1)}\| }}.
@f]
Using the same argument about neighboring cells possibly being smaller by
a factor of two then leads to the final formula for time step zero:
@f[
  \Delta t
  =
  \min_i
  \sqrt{ \frac{\frac 12 h_i m}{e \|\nabla V^{(1)}\| } }.
@f]

Strictly speaking, we would have to evaluate the electric potential $V^{(1)}$ at
the location of each particle, but a good enough approximation is to use the
maximum of the values at the vertices of the respective cell. (Why the vertices
and not the midpoint? Because the gradient of the solution of the Laplace equation,
i.e., the electric field, is largest in corner singularities which are located
at the vertices of cells.) This has the advantage that we can make good use of the
FEValues functionality which can recycle pre-computed material as long as the
quadrature points are the same from one cell to the next.

We could always run this kind of scheme to estimate the difference between
$\mathbf v_i^{(n-1)}$ and $\mathbf v_i^{(n)}$, but it relies on evaluating the
electric field $\mathbf E$ on each cell, and that is expensive. As a
consequence, we will limit this approach to the very first time step.


<a name="Spatialdiscretization"></a><h3>Spatial discretization</h3>


Having discussed the time discretization, the discussion of the spatial
discretization is going to be short: We use quadratic finite elements,
i.e., the space $Q_2$, to approximate the electric potential $V$. The
mesh is adapted a couple of times during the initial time step. All
of this is entirely standard if you have read step-6, and the implementation
does not provide for any kind of surprise.



<a name="Dealingwithparticlesprogrammatically"></a><h3>Dealing with particles programmatically</h3>


Adding and moving particles is, in practice, not very difficult in deal.II.
To add one, the `create_particles()` function of this program simply
uses a code snippet of the following form:
@code
  Particles::Particle<dim> new_particle;
  new_particle.set_location(location);
  new_particle.set_reference_location
      (mapping.transform_real_to_unit_cell(cell, location));
  new_particle.set_id(n_current_particles);

  particle_handler.insert_particle(new_particle, cell);
@endcode
In other words, it is not all that different from inserting an object
into a `std::set` or `std::map`: Create the object, set its properties
(here, the current location, its reference cell location, and its id)
and call `insert_particle`. The only thing that may be surprising is
the reference location: In order to evaluate things such as
$\nabla V(\mathbf x_i)$, it is necessary to evaluate finite element
fields at locations $\mathbf x_i$. But this requires evaluating the
finite element shape functions at points on the refence cell
$\hat{\mathbf x}_i$. To make this efficient, every particle doesn't
just store its location and the cell it is on, but also what location
that point corresponds to in the cell's reference coordinate system.

Updating a particle's position is then no more difficult: One just has
to call
@code
  particle->set_location(new_location);
@endcode
We do this in the `move_particles()` function. The only difference
is that we then have to tell the Particles::ParticleHandler class
to also find what cell that position corresponds to (and, when computing
in parallel, which process owns this cell). For efficiency reason,
this is most easily done after updating all particles' locations,
and is achieved via the
Particles::ParticleHandler::sort_particles_into_subdomains_and_cells()
function.

There are, of course, times where a particle may leave the domain in
question. In that case,
Particles::ParticleHandler::sort_particles_into_subdomains_and_cells()
can not find a surrounding cell and simply deletes the particle. But, it
is often useful to track the number of particles that have been lost
this way, and for this the Particles::ParticleHandler class offers a
"signal" that one can attach to. We show how to do this in the
constructor of the main class to count how many particles were lost
in each time step. Specifically, the way this works is that
the Particles::ParticleHandler class has a "signal" to which one
can attach a function that will be executed whenever the signal
is triggered. Here, this looks as follows:
@code
    particle_handler.signals.particle_lost.connect(
      [this](const typename Particles::ParticleIterator<dim> &        particle,
             const typename Triangulation<dim>::active_cell_iterator &cell)
      {
        this->track_lost_particle(particle, cell);
      });
@endcode
That's a bit of a mouthful, but what's happening is this: We declare
a lambda function that "captures" the `this` pointer (so that we can access
member functions of the surrounding object inside the lambda function), and
that takes two arguments:
- A reference to the particle that has been "lost".
- A reference to the cell it was on last.
The lambda function then simply calls the `CathodeRaySimulator::track_lost_particle`
function with these arguments. When we attach this lambda function to the
signal, the Particles::ParticleHandler::sort_particles_into_subdomains_and_cells()
function will trigger the signal for every particle for which it can't
find a new home. This gives us the chance to record where the particle
is, and to record statistics on it.


@note In this tutorial program, we insert particles by hand and at
  locations we specifically choose based on conditions that include
  the solution of the electrostatic problem. But there are other cases
  where one primarily wants to use particles as passive objects, for
  example to trace and visualize the flow field of a fluid flow
  problem. In those cases, there are numerous functions in the
  Particles::Generators namespace that can generate particles
  automatically. One of the functions of this namespace is also used
  in the step-70 tutorial program, for example.


<a name="Thetestcase"></a><h3>The test case</h3>


The test case here is not meant to be a realistic depiction of a cathode
ray tube, but it has the right general characteristics and the point is,
in any case, only to demonstrate how one would implement deal.II codes
that use particles.

The following picture shows the geometry that we're going to use:

<p align="center">
  <img src="https://www.dealii.org/images/steps/developer/step-19.geometry.png"
       alt="The geometry used in this program"
       width="600">
</p>

In this picture, the parts of the boundary marked in red and blue are the
cathode, held at an electric potential $V=-V_0$. The part of the cathode shown
in red is the part that is heated, leading to electrons leaving the metal
and then being accelerated by the electric field (a few electric
field lines are also shown). The green part of the boundary is the anode,
held at $V=+V_0$. The rest of the boundary satisfies a Neumann boundary
condition.

This setup mimicks real devices. The re-entrant corner results in an
electric potential $V$ whose derivative (the electric field $\mathbf E$)
has a singularity -- in other words, it becomes very large in the vicinity
of the corner, allowing it to rip electrons away from the metal. These
electrons are then accelerated towards the (green) anode which has a
hole in the middle through which the electrons can escape the device and
fly on to hit the screen, where they excite the "phosphor" to then emit
the light that we see from these old-style TV screens. The non-heated
part of the cathode is not subject
to the emission of electrons -- in the code, we will mark this as the
"focussing element" of the tube, because its negative electric voltage
repels the electrons and makes sure that they do not just fly
away from the heated part of the cathode perpendicular to the boundary,
but in fact bend their paths towards the anode on the right.

The electric field lines also shown in the picture illustrate
that the electric field connects the negative and positive
electrodes, respectively. The accelerating force the electrons
experience is along these field lines. Finally, the picture shows the
mesh used in the computation, illustrating that there are
singularities at the tip of the re-rentrant corner as well
as at all places where the boundary conditions change; these
singularities are visible because the mesh is refined in these
locations.

Of practical interest is to figure out which fraction of the
electrons emitted from the cathode actually make it through the
hole in the anode -- electrons that just bounce into the anode
itself are not actually doing anything useful other than converting
eletricity into heat. As a consequence, in the `track_lost_particle()`
function (which is called for each particle that leaves the domain,
see above), we will estimate where it might have left the domain
and report this in the output.


@note It is worth repeating that neither the geometry used here,
nor in fact any other aspect of this program is intended to represent
anything even half-way realistic. Tutorial programs are our tools to
teach how deal.II works, and we often use situations for which we
have some kind of intuition since this helps us interpret the output
of a program, but that's about the extent to which we intend the
program to do anything of use besides being a teaching tool.
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
 * The majority of the include files used in this program are
 * well known from step-6 and similar programs:
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/fe/fe_point_evaluation.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * 
 * @endcode
 * 
 * The ones that are new are only the following three: The first declares the
 * DiscreteTime class that helps us keep track of time in a time-dependent
 * simulation. The latter two provide all of the particle functionality,
 * namely a way to keep track of particles located on a mesh (the
 * Particles::ParticleHandler class) and the ability to output these
 * particles' locations and their properties for the purposes of
 * visualization (the Particles::DataOut class).
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/particles/particle_handler.h>
 * #include <deal.II/particles/data_out.h>
 * 
 * #include <fstream>
 * 
 * using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Globaldefinitions"></a> 
 * <h3>Global definitions</h3>
 * 

 * 
 * As is customary, we put everything that corresponds to the details of the
 * program into a namespace of its own. At the top, we define a few constants
 * for which we would rather use symbolic names than hard-coded numbers.
 * 

 * 
 * Specifically, we define numbers for
 * @ref GlossBoundaryIndicator "boundary indicators"
 * for the various parts of the geometry, as well as the physical properties
 * of electrons and other specifics of the setup we use here.
 * 

 * 
 * For the boundary indicators, let us start enumerating at some
 * random value 101. The principle here is to use numbers that are
 * *uncommon*. If there are pre-defined boundary indicators previously
 * set by the `GridGenerator` functions, they will likely be small
 * integers starting from zero, but not in this rather randomly chosen
 * range. Using numbers such as those below avoids the possibility for
 * conflicts, and also reduces the temptation to just spell these
 * numbers out in the program (because you will probably never
 * remember which is which, whereas you might have been tempted if
 * they had started at 0).
 * 
 * @code
 * namespace Step19
 * {
 *   namespace BoundaryIds
 *   {
 *     constexpr types::boundary_id open          = 101;
 *     constexpr types::boundary_id cathode       = 102;
 *     constexpr types::boundary_id focus_element = 103;
 *     constexpr types::boundary_id anode         = 104;
 *   } // namespace BoundaryIds
 * 
 *   namespace Constants
 *   {
 *     constexpr double electron_mass   = 9.1093837015e-31;
 *     constexpr double electron_charge = 1.602176634e-19;
 * 
 *     constexpr double V0 = 1;
 * 
 *     constexpr double E_threshold = 0.05;
 * 
 *     constexpr double electrons_per_particle = 3e15;
 *   } // namespace Constants
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * The following is then the main class of this program. It has,
 * fundamentally, the same structure as step-6 and many other
 * tutorial programs. This includes the majority of the member
 * functions (with the purpose of the rest probably self-explanatory
 * from their names) as well as only a small number of member
 * variables beyond those of step-6, all of which are related to
 * dealing with particles.
 * 
 * @code
 *   template <int dim>
 *   class CathodeRaySimulator
 *   {
 *   public:
 *     CathodeRaySimulator();
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve_field();
 *     void refine_grid();
 * 
 *     void create_particles();
 *     void move_particles();
 *     void track_lost_particle(
 *       const typename Particles::ParticleIterator<dim> &        particle,
 *       const typename Triangulation<dim>::active_cell_iterator &cell);
 * 
 * 
 *     void update_timestep_size();
 *     void output_results() const;
 * 
 *     Triangulation<dim>        triangulation;
 *     MappingQ<dim>             mapping;
 *     FE_Q<dim>                 fe;
 *     DoFHandler<dim>           dof_handler;
 *     AffineConstraints<double> constraints;
 * 
 *     SparseMatrix<double> system_matrix;
 *     SparsityPattern      sparsity_pattern;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     Particles::ParticleHandler<dim> particle_handler;
 *     types::particle_index           next_unused_particle_id;
 *     types::particle_index           n_recently_lost_particles;
 *     types::particle_index           n_total_lost_particles;
 *     types::particle_index           n_particles_lost_through_anode;
 * 
 *     DiscreteTime time;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorcodeclassimplementation"></a> 
 * <h3>The <code>CathodeRaySimulator</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorcodeconstructor"></a> 
 * <h4>The <code>CathodeRaySimulator</code> constructor</h4>
 * 

 * 
 * So then let us get started on the implementation. What the constructor
 * does is really only a straight-forward initialization of all of the member
 * variables at the top. The only two worth mentioning are the
 * `particle_handler`, which is handed a reference to the triangulation
 * on which the particles will live (currently of course still empty,
 * but the particle handler stores the reference and will use it once
 * particles are added -- which happens after the triangulation is built).
 * The other piece of information it gets is how many "properties"
 * each particle needs to store. Here, all we need each particle to
 * remember is its current velocity, i.e., a vector with `dim`
 * components. There are, however, other intrinsic properties that
 * each particle has and that the Particles::ParticleHandler class
 * automatically and always makes sure are available; in particular,
 * these are the current location of a particle, the cell it is on,
 * it's reference location within that cell, and the particle's ID.
 *   

 * 
 * The only other variable of interest is `time`, an object of type
 * DiscreteTime. It keeps track of the current time we are in a
 * time-dependent simulation, and is initialized with the start time
 * (zero) and end time ($10^{-4}$). We will later set the time step
 * size in `update_timestep_size()`.
 *   

 * 
 * The body of the constructor consists of a piece of code we have
 * already discussed in the introduction. Namely, we make sure that the
 * `track_lost_particle()` function is called by the `particle_handler`
 * object every time a particle leaves the domain.
 * 
 * @code
 *   template <int dim>
 *   CathodeRaySimulator<dim>::CathodeRaySimulator()
 *     : mapping(1)
 *     , fe(2)
 *     , dof_handler(triangulation)
 *     , particle_handler(triangulation, mapping, /*n_properties=*/dim)
 *     , next_unused_particle_id(0)
 *     , n_recently_lost_particles(0)
 *     , n_total_lost_particles(0)
 *     , n_particles_lost_through_anode(0)
 *     , time(0, 1e-4)
 *   {
 *     particle_handler.signals.particle_lost.connect(
 *       [this](const typename Particles::ParticleIterator<dim> &        particle,
 *              const typename Triangulation<dim>::active_cell_iterator &cell) {
 *         this->track_lost_particle(particle, cell);
 *       });
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatormake_gridcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::make_grid</code> function</h4>
 * 

 * 
 * The next function is then responsible for generating the mesh on which
 * we want to solve. Recall how the domain looks like:
 * <p align="center">
 * <img
 * src="https://www.dealii.org/images/steps/developer/step-19.geometry.png"
 * alt="The geometry used in this program"
 * width="600">
 * </p>
 * We subdivide this geometry into a mesh of $4\times 2$ cells that looks
 * like this:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   *---*---*---*---*
 *   \   |   |   |   |
 *    *--*---*---*---*
 *   /   |   |   |   |
 *   *---*---*---*---*
 * @endcode
 * </div>
 * The way this is done is by first defining where the $15=5\times 3$
 * vertices are located -- here, we say that they are on integer points
 * with the middle one on the left side moved to the right by a value of
 * `delta=0.5`.
 *   

 * 
 * In the following, we then have to say which vertices together form
 * the 8 cells. The following code is then entirely equivalent to what
 * we also do in step-14:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::make_grid()
 *   {
 *     static_assert(dim == 2,
 *                   "This function is currently only implemented for 2d.");
 * 
 *     const double       delta = 0.5;
 *     const unsigned int nx    = 5;
 *     const unsigned int ny    = 3;
 * 
 *     const std::vector<Point<dim>> vertices 
 *       = {{0, 0},
 *          {1, 0},
 *          {2, 0},
 *          {3, 0},
 *          {4, 0},
 *          {delta, 1},
 *          {1, 1},
 *          {2, 1},
 *          {3, 1},
 *          {4, 1},
 *          {0, 2},
 *          {1, 2},
 *          {2, 2},
 *          {3, 2},
 *          {4, 2}};
 *     AssertDimension(vertices.size(), nx * ny);
 * 
 *     const std::vector<unsigned int> cell_vertices[(nx - 1) * (ny - 1)] = {
 *       {0, 1, nx + 0, nx + 1},
 *       {1, 2, nx + 1, nx + 2},
 *       {2, 3, nx + 2, nx + 3},
 *       {3, 4, nx + 3, nx + 4},
 * 
 *       {5, nx + 1, 2 * nx + 0, 2 * nx + 1},
 *       {nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2},
 *       {nx + 2, nx + 3, 2 * nx + 2, 2 * nx + 3},
 *       {nx + 3, nx + 4, 2 * nx + 3, 2 * nx + 4}};
 * 
 * @endcode
 * 
 * With these arrays out of the way, we can move to slightly higher
 * higher-level data structures. We create a vector of CellData
 * objects that store for each cell to be created the vertices in
 * question as well as the @ref GlossMaterialId "material id" (which
 * we will here simply set to zero since we don't use it in the program).
 *     

 * 
 * This information is then handed to the
 * Triangulation::create_triangulation() function, and the mesh is twice
 * globally refined.
 * 
 * @code
 *     std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>());
 *     for (unsigned int i = 0; i < cells.size(); ++i)
 *       {
 *         cells[i].vertices    = cell_vertices[i];
 *         cells[i].material_id = 0;
 *       }
 * 
 *     triangulation.create_triangulation(
 *       vertices,
 *       cells,
 *       SubCellData()); // No boundary information
 * 
 *     triangulation.refine_global(2);
 * 
 * @endcode
 * 
 * The remaining part of the function loops over all cells and their faces,
 * and if a face is at the boundary determines which boundary indicator
 * should be applied to it. The various conditions should make sense if
 * you compare the code with the picture of the geometry above.
 *     

 * 
 * Once done with this step, we refine the mesh once more globally.
 * 
 * @code
 *     for (auto &cell : triangulation.active_cell_iterators())
 *       for (auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             if ((face->center()[0] > 0) && (face->center()[0] < 0.5) &&
 *                 (face->center()[1] > 0) && (face->center()[1] < 2))
 *               face->set_boundary_id(BoundaryIds::cathode);
 *             else if ((face->center()[0] > 0) && (face->center()[0] < 2))
 *               face->set_boundary_id(BoundaryIds::focus_element);
 *             else if ((face->center()[0] > 4 - 1e-12) &&
 *                      ((face->center()[1] > 1.5) || (face->center()[1] < 0.5)))
 *               face->set_boundary_id(BoundaryIds::anode);
 *             else
 *               face->set_boundary_id(BoundaryIds::open);
 *           }
 * 
 *     triangulation.refine_global(1);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorsetup_systemcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::setup_system</code> function</h4>
 * 

 * 
 * The next function in this program deals with setting up the various
 * objects related to solving the partial differential equations. It is
 * in essence a copy of the corresponding function in step-6 and requires
 * no further discussion.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::cathode,
 *                                              Functions::ConstantFunction<dim>(
 *                                                -Constants::V0),
 *                                              constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::focus_element,
 *                                              Functions::ConstantFunction<dim>(
 *                                                -Constants::V0),
 *                                              constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::anode,
 *                                              Functions::ConstantFunction<dim>(
 *                                                +Constants::V0),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ false);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorassemble_systemcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::assemble_system</code> function</h4>
 * 

 * 
 * The function that computes
 * the matrix entries is again in essence a copy of the
 * corresponding function in step-6:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.dofs_per_cell;
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               for (const unsigned int j : fe_values.dof_indices())
 *                 cell_matrix(i, j) +=
 *                   (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                    fe_values.JxW(q_index));           // dx
 *             }
 * 
 * @endcode
 * 
 * The only interesting part of this function is how it forms the right
 * hand side of the linear system. Recall that the right hand side
 * of the PDE is
 * @f[
 * \sum_p (N e)\delta(\mathbf x-\mathbf x_p),
 * @f]
 * where we have used $p$ to index the particles here to avoid
 * confusion with the shape function $\varphi_i$; $\mathbf x_p$
 * is the position of the $p$th particle.
 *         

 * 
 * When multiplied by a test function $\varphi_i$ and integrated over
 * the domain results in a right hand side vector
 * @f{align*}{
 * F_i &= \int_\Omega \varphi_i (\mathbf x)\left[
 * \sum_p (N e)\delta(\mathbf x-\mathbf x_p) \right] dx
 * \\  &=  \sum_p (N e) \varphi_i(\mathbf x_p).
 * @f}
 * Note that the final line no longer contains an integral, and
 * consequently also no occurrence of $dx$ which would require the
 * appearance of the `JxW` symbol in our code.
 *         

 * 
 * For a given cell $K$, this cell's contribution to the right hand
 * side is then
 * @f{align*}{
 * F_i^K &= \sum_{p, \mathbf x_p\in K} (N e) \varphi_i(\mathbf x_p),
 * @f}
 * i.e., we only have to worry about those particles that are actually
 * located on the current cell $K$.
 *         

 * 
 * In practice, what we do here is the following: If there are any
 * particles on the current cell, then we first obtain an iterator range
 * pointing to the first particle of that cell as well as the particle
 * past the last one on this cell (or the end iterator) -- i.e., a
 * half-open range as is common for C++ functions. Knowing now the list
 * of particles, we query their reference locations (with respect to
 * the reference cell), evaluate the shape functions in these reference
 * locations, and compute the force according to the formula above
 * (without any FEValues::JxW).
 *         

 * 
 * @note It is worth pointing out that calling the
 * Particles::ParticleHandler::particles_in_cell() and
 * Particles::ParticleHandler::n_particles_in_cell() functions is not
 * very efficient on problems with a large number of particles. But it
 * illustrates the easiest way to write this algorithm, and so we are
 * willing to incur this cost for the moment for expository purposes.
 * We discuss the issue in more detail in the
 * <a href="#extensions">"possibilities for extensions" section</a>
 * below, and use a better approach in step-70, for example.
 * 
 * @code
 *         if (particle_handler.n_particles_in_cell(cell) > 0)
 *           for (const auto &particle : particle_handler.particles_in_cell(cell))
 *             {
 *               const Point<dim> reference_location =
 *                 particle.get_reference_location();
 *               for (const unsigned int i : fe_values.dof_indices())
 *                 cell_rhs(i) +=
 *                   (fe.shape_value(i, reference_location) * // phi_i(x_p)
 *                    (-Constants::electrons_per_particle *   // N
 *                     Constants::electron_charge));          // e
 *             }
 * 
 * @endcode
 * 
 * Finally, we can copy the contributions of this cell into
 * the global matrix and right hand side vector:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorsolve"></a> 
 * <h4>CathodeRaySimulator::solve</h4>
 * 

 * 
 * The function that solves the linear system is then again exactly as in
 * step-6:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::solve_field()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorrefine_grid"></a> 
 * <h4>CathodeRaySimulator::refine_grid</h4>
 * 

 * 
 * The final field-related function is the one that refines the grid. We will
 * call it a number of times in the first time step to obtain a mesh that
 * is well-adapted to the structure of the solution and, in particular,
 * resolves the various singularities in the solution that are due to
 * re-entrant corners and places where the boundary condition type
 * changes. You might want to refer to step-6 again for more details:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::refine_grid()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                        QGauss<dim - 1>(fe.degree + 1),
 *                                        {},
 *                                        solution,
 *                                        estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.1,
 *                                                     0.03);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorcreate_particles"></a> 
 * <h4>CathodeRaySimulator::create_particles</h4>
 * 

 * 
 * Let us now turn to the functions that deal with particles. The first one
 * is about the creation of particles. As mentioned in the introduction,
 * we want to create a particle at points of the cathode if the the electric
 * field $\mathbf E=\nabla V$ exceeds a certain threshold, i.e., if
 * $|\mathbf E| \ge E_\text{threshold}$, and if furthermore the electric field
 * points into the domain (i.e., if $\mathbf E \cdot \mathbf n < 0$). As is
 * common in the finite element method, we evaluate fields (and their
 * derivatives) at specific evaluation points; typically, these are
 * "quadrature points", and so we create a "quadrature formula" that we will
 * use to designate the points at which we want to evaluate the solution.
 * Here, we will simply take QMidpoint implying that we will only check the
 * threshold condition at the midpoints of faces. We then use this to
 * initialize an object of type FEFaceValues to evaluate the solution at these
 * points.
 *   

 * 
 * All of this will then be used in a loop over all cells, their faces, and
 * specifically those faces that are at the boundary and, moreover, the
 * cathode part of the boundary.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::create_particles()
 *   {
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      QMidpoint<dim - 1>(),
 *                                      update_quadrature_points |
 *                                        update_gradients |
 *                                        update_normal_vectors);
 * 
 *     std::vector<Tensor<1, dim>> solution_gradients(
 *       fe_face_values.n_quadrature_points);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary() &&
 *             (face->boundary_id() == BoundaryIds::cathode))
 *           {
 *             fe_face_values.reinit(cell, face);
 * 
 * @endcode
 * 
 * So we have found a face on the cathode. Next, we let the
 * FEFaceValues object compute the gradient of the solution at each
 * "quadrature" point, and extract the electric field vector from
 * the gradient in the form of a Tensor variable through the methods
 * discussed in the
 * @ref vector_valued "vector-valued problems" documentation module.
 * 
 * @code
 *             const FEValuesExtractors::Scalar electric_potential(0);
 *             fe_face_values[electric_potential].get_function_gradients(
 *               solution, solution_gradients);
 *             for (const unsigned int q_point :
 *                  fe_face_values.quadrature_point_indices())
 *               {
 *                 const Tensor<1, dim> E = solution_gradients[q_point];
 * 
 * @endcode
 * 
 * Electrons can only escape the cathode if the electric field
 * strength exceeds a threshold and,
 * crucially, if the electric field points *into* the domain.
 * Once we have that checked, we create a new
 * Particles::Particle object at this location and insert it
 * into the Particles::ParticleHandler object with a unique ID.
 *                 

 * 
 * The only thing that may be not obvious here is that we also
 * associate with this particle the location in the reference
 * coordinates of the cell we are currently on. This is done
 * because we will in downstream functions compute quantities
 * such as the electric field at the location of the particle
 * (e.g., to compute the forces that act on it when updating its
 * position in each time step). Evaluating a finite element
 * field at arbitrary coordinates is quite an expensive
 * operation because shape functions are really only defined on
 * the reference cell, and so when asking for the electric field
 * at an arbitrary point requires us first to determine what the
 * reference coordinates of that point are. To avoid having to
 * do this over and over, we determine these coordinates once
 * and for all and then store these reference coordinates
 * directly with the particle.
 * 
 * @code
 *                 if ((E * fe_face_values.normal_vector(q_point) < 0) &&
 *                     (E.norm() > Constants::E_threshold))
 *                   {
 *                     const Point<dim> location =
 *                       fe_face_values.quadrature_point(q_point);
 * 
 *                     Particles::Particle<dim> new_particle;
 *                     new_particle.set_location(location);
 *                     new_particle.set_reference_location(
 *                       mapping.transform_real_to_unit_cell(cell, location));
 *                     new_particle.set_id(next_unused_particle_id);
 *                     particle_handler.insert_particle(new_particle, cell);
 * 
 *                     ++next_unused_particle_id;
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * At the end of all of these insertions, we let the `particle_handler`
 * update some internal statistics about the particles it stores.
 * 
 * @code
 *     particle_handler.update_cached_numbers();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatormove_particles"></a> 
 * <h4>CathodeRaySimulator::move_particles</h4>
 * 

 * 
 * The second particle-related function is the one that moves the particles
 * in each time step. To do this, we have to loop over all cells, the
 * particles in each cell, and evaluate the electric field at each of the
 * particles' positions.
 *   

 * 
 * The approach used here is conceptually the same used in the
 * `assemble_system()` function: We loop over all cells, find the particles
 * located there (with the same caveat about the inefficiency of the algorithm
 * used here to find these particles), and use FEPointEvaluation object to
 * evaluate the gradient at these positions:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::move_particles()
 *   {
 *     const double dt = time.get_next_step_size();
 * 
 *     Vector<double>            solution_values(fe.n_dofs_per_cell());
 *     FEPointEvaluation<1, dim> evaluator(mapping, fe);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (particle_handler.n_particles_in_cell(cell) > 0)
 *         {
 *           const typename Particles::ParticleHandler<
 *             dim>::particle_iterator_range particles_in_cell =
 *             particle_handler.particles_in_cell(cell);
 * 
 *           std::vector<Point<dim>> particle_positions;
 *           for (const auto &particle : particles_in_cell)
 *             particle_positions.push_back(particle.get_reference_location());
 * 
 *           cell->get_dof_values(solution, solution_values);
 * 
 * @endcode
 * 
 * Then we can ask the FEPointEvaluation object for the gradients of
 * the solution (i.e., the electric field $\mathbf E$) at these
 * locations and loop over the individual particles:
 * 
 * @code
 *           evaluator.evaluate(cell,
 *                              particle_positions,
 *                              make_array_view(solution_values),
 *                              EvaluationFlags::gradients);
 * 
 *           {
 *             typename Particles::ParticleHandler<dim>::particle_iterator
 *               particle = particles_in_cell.begin();
 *             for (unsigned int particle_index = 0;
 *                  particle != particles_in_cell.end();
 *                  ++particle, ++particle_index)
 *               {
 *                 const Tensor<1, dim> E = evaluator.get_gradient(particle_index);
 * 
 * @endcode
 * 
 * Having now obtained the electric field at the location of one
 * of the particles, we use this to update first the velocity
 * and then the position. To do so, let us first get the old
 * velocity out of the properties stored with the particle,
 * compute the acceleration, update the velocity, and store this
 * new velocity again in the properties of the particle. Recall
 * that this corresponds to the first of the following set of
 * update equations discussed in the introduction:
 * @f{align*}{
 * \frac{{\mathbf v}_i^{(n)}
 * -{\mathbf v}_i^{(n-1)}}{\Delta t}
 * &= \frac{e\nabla V^{(n)}}{m}
 * \\ \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}
 * {\Delta t} &= {\mathbf v}_i^{(n)}.
 * @f}
 * 
 * @code
 *                 const Tensor<1, dim> old_velocity(particle->get_properties());
 * 
 *                 const Tensor<1, dim> acceleration =
 *                   Constants::electron_charge / Constants::electron_mass * E;
 * 
 *                 const Tensor<1, dim> new_velocity =
 *                   old_velocity + acceleration * dt;
 * 
 *                 particle->set_properties(make_array_view(new_velocity));
 * 
 * @endcode
 * 
 * With the new velocity, we can then also update the location
 * of the particle and tell the particle about it.
 * 
 * @code
 *                 const Point<dim> new_location =
 *                   particle->get_location() + dt * new_velocity;
 *                 particle->set_location(new_location);
 *               }
 *           }
 *         }
 * 
 * @endcode
 * 
 * Having updated the locations and properties (i.e., velocities) of all
 * particles, we need to make sure that the `particle_handler` again knows
 * which cells they are in, and what their locations in the coordinate
 * system of the reference cell are. The following function does that. (It
 * also makes sure that, in parallel computations, particles are moved from
 * one processor to another processor if a particle moves from the subdomain
 * owned by the former to the subdomain owned by the latter.)
 * 
 * @code
 *     particle_handler.sort_particles_into_subdomains_and_cells();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatortrack_lost_particle"></a> 
 * <h4>CathodeRaySimulator::track_lost_particle</h4>
 * 

 * 
 * The final particle-related function is the one that is called whenever a
 * particle is lost from the simulation. This typically happens if it leaves
 * the domain. If that happens, this function is called both the cell (which
 * we can ask for its new location) and the cell it was previously on. The
 * function then keeps track of updating the number of particles lost in this
 * time step, the total number of lost particles, and then estimates whether
 * the particle left through the hole in the middle of the anode. We do so by
 * first checking whether the cell it was in last had an $x$ coordinate to the
 * left of the right boundary (located at $x=4$) and the particle now has a
 * position to the right of the right boundary. If that is so, we compute a
 * direction vector of its motion that is normalized so that the $x$ component
 * of the direction vector is equal to $1$. With this direction vector, we can
 * compute where it would have intersected the line $x=4$. If this intersect
 * is between $0.5$ and $1.5$, then we claim that the particle left through
 * the hole and increment a counter.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::track_lost_particle(
 *     const typename Particles::ParticleIterator<dim> &        particle,
 *     const typename Triangulation<dim>::active_cell_iterator &cell)
 *   {
 *     ++n_recently_lost_particles;
 *     ++n_total_lost_particles;
 * 
 *     const Point<dim> current_location              = particle->get_location();
 *     const Point<dim> approximate_previous_location = cell->center();
 * 
 *     if ((approximate_previous_location[0] < 4) && (current_location[0] > 4))
 *       {
 *         const Tensor<1, dim> direction =
 *           (current_location - approximate_previous_location) /
 *           (current_location[0] - approximate_previous_location[0]);
 * 
 *         const double right_boundary_intercept =
 *           approximate_previous_location[1] +
 *           (4 - approximate_previous_location[0]) * direction[1];
 *         if ((right_boundary_intercept > 0.5) &&
 *             (right_boundary_intercept < 1.5))
 *           ++n_particles_lost_through_anode;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorupdate_timestep_size"></a> 
 * <h4>CathodeRaySimulator::update_timestep_size</h4>
 * 

 * 
 * As discussed at length in the introduction, we need to respect a time step
 * condition whereby particles can not move further than one cell in one time
 * step. To ensure that this is the case, we again first compute the maximal
 * speed of all particles on each cell, and divide the cell size by that
 * speed. We then compute the next time step size as the minimum of this
 * quantity over all cells, using the safety factor discussed in the
 * introduction, and set this as the desired time step size using the
 * DiscreteTime::set_desired_time_step_size() function.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::update_timestep_size()
 *   {
 *     if (time.get_step_number() > 0)
 *       {
 *         double min_cell_size_over_velocity = std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : dof_handler.active_cell_iterators())
 *           if (particle_handler.n_particles_in_cell(cell) > 0)
 *             {
 *               const double cell_size = cell->minimum_vertex_distance();
 * 
 *               double max_particle_velocity(0.0);
 * 
 *               for (const auto &particle :
 *                    particle_handler.particles_in_cell(cell))
 *                 {
 *                   const Tensor<1, dim> velocity(particle.get_properties());
 *                   max_particle_velocity =
 *                     std::max(max_particle_velocity, velocity.norm());
 *                 }
 * 
 *               if (max_particle_velocity > 0)
 *                 min_cell_size_over_velocity =
 *                   std::min(min_cell_size_over_velocity,
 *                            cell_size / max_particle_velocity);
 *             }
 * 
 *         constexpr double c_safety = 0.5;
 *         time.set_desired_next_step_size(c_safety * 0.5 *
 *                                         min_cell_size_over_velocity);
 *       }
 * @endcode
 * 
 * As mentioned in the introduction, we have to treat the very first
 * time step differently since there, particles are not available yet or
 * do not yet have the information associated that we need for the
 * computation of a reasonable step length. The formulas below follow the
 * discussion in the introduction.
 * 
 * @code
 *     else
 *       {
 *         const QTrapezoid<dim> vertex_quadrature;
 *         FEValues<dim> fe_values(fe, vertex_quadrature, update_gradients);
 * 
 *         std::vector<Tensor<1, dim>> field_gradients(vertex_quadrature.size());
 * 
 *         double min_timestep = std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : dof_handler.active_cell_iterators())
 *           if (particle_handler.n_particles_in_cell(cell) > 0)
 *             {
 *               const double cell_size = cell->minimum_vertex_distance();
 * 
 *               fe_values.reinit(cell);
 *               fe_values.get_function_gradients(solution, field_gradients);
 * 
 *               double max_E = 0;
 *               for (const auto q_point : fe_values.quadrature_point_indices())
 *                 max_E = std::max(max_E, field_gradients[q_point].norm());
 * 
 *               if (max_E > 0)
 *                 min_timestep =
 *                   std::min(min_timestep,
 *                            std::sqrt(0.5 * cell_size *
 *                                      Constants::electron_mass /
 *                                      Constants::electron_charge / max_E));
 *             }
 * 
 *         time.set_desired_next_step_size(min_timestep);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatoroutput_resultscodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::output_results()</code> function</h4>
 * 

 * 
 * The final function implementing pieces of the overall algorithm is the one
 * that generates graphical output. In the current context, we want to output
 * both the electric potential field as well as the particle locations and
 * velocities. But we also want to output the electric field, i.e., the
 * gradient of the solution.
 *   

 * 
 * deal.II has a general way how one can compute derived quantities from
 * the solution and output those as well. Here, this is the electric
 * field, but it could also be some other quantity -- say, the norm of the
 * electric field, or in fact anything else one could want to compute from
 * the solution $V_h(\mathbf x)$ or its derivatives. This general solution
 * uses the DataPostprocessor class and, in cases like the one here where we
 * want to output a quantity that represents a vector field, the
 * DataPostprocessorVector class.
 *   

 * 
 * Rather than try and explain how this class works, let us simply refer to
 * the documentation of the DataPostprocessorVector class that has essentially
 * this case as a well-documented example.
 * 
 * @code
 *   template <int dim>
 *   class ElectricFieldPostprocessor : public DataPostprocessorVector<dim>
 *   {
 *   public:
 *     ElectricFieldPostprocessor()
 *       : DataPostprocessorVector<dim>("electric_field", update_gradients)
 *     {}
 * 
 *     virtual void evaluate_scalar_field(
 *       const DataPostprocessorInputs::Scalar<dim> &input_data,
 *       std::vector<Vector<double>> &computed_quantities) const override
 *     {
 *       AssertDimension(input_data.solution_gradients.size(),
 *                       computed_quantities.size());
 * 
 *       for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
 *         {
 *           AssertDimension(computed_quantities[p].size(), dim);
 *           for (unsigned int d = 0; d < dim; ++d)
 *             computed_quantities[p][d] = input_data.solution_gradients[p][d];
 *         }
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * With this, the `output_results()` function becomes relatively
 * straightforward: We use the DataOut class as we have in almost every one of
 * the previous tutorial programs to output the solution (the "electric
 * potential") and we use the postprocessor defined above to also output its
 * gradient (the "electric field"). This all is then written into a file in
 * VTU format after also associating the current time and time step number
 * with this file.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::output_results() const
 *   {
 *     {
 *       ElectricFieldPostprocessor<dim> electric_field;
 *       DataOut<dim>                    data_out;
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "electric_potential");
 *       data_out.add_data_vector(solution, electric_field);
 *       data_out.build_patches();
 * 
 *       data_out.set_flags(
 *         DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));
 * 
 *       std::ofstream output("solution-" +
 *                            Utilities::int_to_string(time.get_step_number(), 4) +
 *                            ".vtu");
 *       data_out.write_vtu(output);
 *     }
 * 
 * @endcode
 * 
 * Output the particle positions and properties is not more complicated. The
 * Particles::DataOut class plays the role of the DataOut class for
 * particles, and all we have to do is tell that class where to take
 * particles from and how to interpret the `dim` components of the
 * properties -- namely, as a single vector indicating the velocity, rather
 * than as `dim` scalar properties. The rest is then the same as above:
 * 
 * @code
 *     {
 *       Particles::DataOut<dim, dim> particle_out;
 *       particle_out.build_patches(
 *         particle_handler,
 *         std::vector<std::string>(dim, "velocity"),
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>(
 *           dim, DataComponentInterpretation::component_is_part_of_vector));
 * 
 *       particle_out.set_flags(
 *         DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));
 * 
 *       std::ofstream output("particles-" +
 *                            Utilities::int_to_string(time.get_step_number(), 4) +
 *                            ".vtu");
 *       particle_out.write_vtu(output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorrun"></a> 
 * <h4>CathodeRaySimulator::run</h4>
 * 

 * 
 * The last member function of the principal class of this program is then the
 * driver. At the top, it refines the mesh a number of times by solving the
 * problem (with not particles yet created) on a sequence of finer and finer
 * meshes.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::run()
 *   {
 *     make_grid();
 * 
 * @endcode
 * 
 * do a few refinement cycles up front
 * 
 * @code
 *     const unsigned int n_pre_refinement_cycles = 3;
 *     for (unsigned int refinement_cycle = 0;
 *          refinement_cycle < n_pre_refinement_cycles;
 *          ++refinement_cycle)
 *       {
 *         setup_system();
 *         assemble_system();
 *         solve_field();
 *         refine_grid();
 *       }
 * 
 * 
 * @endcode
 * 
 * Now do the loop over time. The sequence of steps follows closely the
 * outline of the algorithm discussed in the introduction. As discussed in
 * great detail in the documentation of the DiscreteTime class, while we
 * move the field and particle information forward by one time step, the
 * time stored in the `time` variable is not consistent with where (some of)
 * these quantities are (in the diction of DiscreteTime, this is the "update
 * stage"). The call to `time.advance_time()` makes everything consistent
 * again by setting the `time` variable to the time at which the field and
 * particles already are, and once we are in this "consistent stage", we can
 * generate graphical output and write information about the current state
 * of the simulation to screen.
 * 
 * @code
 *     setup_system();
 *     do
 *       {
 *         std::cout << "Timestep " << time.get_step_number() + 1 << std::endl;
 *         std::cout << "  Field degrees of freedom:                 "
 *                   << dof_handler.n_dofs() << std::endl;
 * 
 *         assemble_system();
 *         solve_field();
 * 
 *         create_particles();
 *         std::cout << "  Total number of particles in simulation:  "
 *                   << particle_handler.n_global_particles() << std::endl;
 * 
 *         n_recently_lost_particles = 0;
 *         update_timestep_size();
 *         move_particles();
 * 
 *         time.advance_time();
 * 
 *         output_results();
 * 
 *         std::cout << "  Number of particles lost this time step:  "
 *                   << n_recently_lost_particles << std::endl;
 *         if (n_total_lost_particles > 0)
 *           std::cout << "  Fraction of particles lost through anode: "
 *                     << 1. * n_particles_lost_through_anode /
 *                          n_total_lost_particles
 *                     << std::endl;
 * 
 *         std::cout << std::endl
 *                   << "  Now at t=" << time.get_current_time()
 *                   << ", dt=" << time.get_previous_step_size() << '.'
 *                   << std::endl
 *                   << std::endl;
 *       }
 *     while (time.is_at_end() == false);
 *   }
 * } // namespace Step19
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The final function of the program is then again the `main()` function. It is
 * unchanged in all tutorial programs since step-6 and so there is nothing new
 * to discuss:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step19::CathodeRaySimulator<2> cathode_ray_simulator_2d;
 *       cathode_ray_simulator_2d.run();
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
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


When this program is run, it produces output that looks as follows:
```
Timestep 1
  Field degrees of freedom:                                 4989
  Total number of particles in simulation:  20
  Number of particles lost this time step:  0

  Now at t=2.12647e-07, dt=2.12647e-07.

Timestep 2
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  24
  Number of particles lost this time step:  0

  Now at t=4.14362e-07, dt=2.01715e-07.

Timestep 3
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  28
  Number of particles lost this time step:  0

  Now at t=5.96019e-07, dt=1.81657e-07.

Timestep 4
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  32
  Number of particles lost this time step:  0

  Now at t=7.42634e-07, dt=1.46614e-07.


...


  Timestep 1000
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  44
  Number of particles lost this time step:  6
  Fraction of particles lost through anode: 0.0601266

  Now at t=4.93276e-05, dt=4.87463e-08.

Timestep 1001
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  44
  Number of particles lost this time step:  0
  Fraction of particles lost through anode: 0.0601266

  Now at t=4.93759e-05, dt=4.82873e-08.


...


Timestep 2091
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  44
  Number of particles lost this time step:  0
  Fraction of particles lost through anode: 0.0503338

  Now at t=9.99237e-05, dt=4.26254e-08.

Timestep 2092
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  44
  Number of particles lost this time step:  0
  Fraction of particles lost through anode: 0.0503338

  Now at t=9.99661e-05, dt=4.24442e-08.

Timestep 2093
  Field degrees of freedom:                 4989
  Total number of particles in simulation:  44
  Number of particles lost this time step:  2
  Fraction of particles lost through anode: 0.050308

  Now at t=0.0001, dt=3.38577e-08.
```

Picking a random few time steps, we can visualize the solution in the
form of streamlines for the electric field and dots for the electrons:
<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-19.solution.0000.png"
         alt="The solution at time step 0 (t=0 seconds)."
         width="500">
    <br>
    Solution at time step 0 (t=0 seconds).
    <br>
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-19.solution.1400.png"
         alt="The solution at time step 1400 (t=0.000068 seconds)."
         width="500">
    <br>
    Solution at time step 1400 (t=0.000068 seconds).
    <br>
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-19.solution.0700.png"
         alt="The solution at time step 700 (t=0.000035 seconds)."
         width="500">
    <br>
    Solution at time step 700 (t=0.000035 seconds).
    <br>
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-19.solution.2092.png"
         alt="The solution at time step 2092 (t=0.0001 seconds)."
         width="500">
    <br>
    Solution at time step 2092 (t=0.0001 seconds).
    <br>
  </div>
</div>

That said, a more appropriate way to visualize the results of this
program are by creating a video that shows how these electrons move, and how
the electric field changes in response to their motion:

@htmlonly
<p align="center">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/HwUtE7xuteE"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly

What you can see here is how the "focus element" of the boundary with its negative
voltage repels the electrons and makes sure that they do not just fly away
perpendicular from the cathode (as they do in the initial part of their
trajectories). It also shows how the electric field lines
move around over time, in response to the charges flying by -- in other words,
the feedback the particles have on the electric field that itself drives the
motion of the electrons.

The movie suggests that electrons move in "bunches" or "bursts". One element of
this appearance is an artifact of how the movie was created: Every frame of the
movie corresponds to one time step, but the time step length varies. More specifically,
the fastest particle moving through the smallest cell determines the length of the
time step (see the discussion in the introduction), and consequently time steps
are small whenever a (fast) particle moves through the small cells at the right
edge of the domain; time steps are longer again once the particle has left
the domain. This slowing-accelerating effect can easily be visualized by plotting
the time step length shown in the screen output.

The second part of this is real, however: The simulation creates a large group
of particles in the beginning, and fewer after about the 300th time step. This
is probably because of the negative charge of the particles in the simulation:
They reduce the magnitude of the electric field at the (also negatively charged
electrode) and consequently reduce the number of points on the cathode at which
the magnitude exceeds the threshold necessary to draw an electron out of the
electrode.


<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Avoidingaperformancebottleneckwithparticles"></a><h4> Avoiding a performance bottleneck with particles </h4>


The `assemble_system()`, `move_particles()`, and `update_timestep_size()`
functions all call Particles::ParticleHandler::particles_in_cell() and
Particles::ParticleHandler::n_particles_in_cell() that query information
about the particles located on the current cell. While this is convenient,
it's also inefficient. To understand why this is so, one needs to know
how particles are stored in Particles::ParticleHandler: namely, in a
data structure in which particles are ordered in some kind of linear
fashion sorted by the cell they are on. Consequently, in order to find
the particles associated with a given cell, these functions need to
search for the first (and possibly last) particle on a given cell --
an effort that costs ${\cal O}(\log N)$ operations where $N$ is the
number of particles. But this is repeated on every cell; assuming that
for large computations, the number of cells and particles are roughly
proportional, the accumulated cost of these function calls is then
${\cal O}(N \log N)$ and consequently larger than the ${\cal O}(N)$
cost that we should shoot for with all parts of a program.

We can make this cheaper, though. First, instead of calling
Particles::ParticleHandler::n_particles_in_cell(), we might first call
Particles::ParticleHandler::particles_in_cell() and then compute the
number of particles on a cell by just computing the distance of the last
to the first particle on the current cell:
@code
  const typename Particles::ParticleHandler<dim, spacedim>::particle_iterator_range
    particles_in_cell = particle_handler.particles_in_cell(cell);
  const unsigned int
    n_particles_in_cell = std::distance (particles_in_cell.begin(),
                                         particles_in_cell.end());
@endcode
The first of these calls is of course still ${\cal O}(\log N)$,
but at least the second call only takes a compute time proportional to
the number of particles on the current cell and so, when accumulated
over all cells, has a cost of ${\cal O}(N)$.

But we can even get rid of the first of these calls with some proper algorithm
design. That's because particles are ordered in the same way as cells, and so
we can just walk them as we move along on the cells. The following outline
of an algorithm does this:
@code
  auto begin_particle_on_cell = particle_handler.begin();
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      unsigned int n_particles_on_cell = 0;
      auto end_particle_on_cell = begin_particle_on_cell;
      while (end_particle_on_cell->get_surrounding_cell(triangulation)
             == cell)
        {
          ++n_particles_on_cell;
          ++end_particle_on_cell;
        }

      ...now operate on the range of particles from begin_particle_on_cell
         to end_particle_on_cell, all of which are known to be on the current
         cell...;

      // Move the begin iterator forward so that it points to the first
      // particle on the next cell
      begin_particle_on_cell = end_particle_on_cell;
    }
@endcode

In this code, we touch every cell exactly once and we never have to search
the big data structure for the first or last particle on each cell. As a
consequence, the algorithm costs a total of ${\cal O}(N)$ for a complete
sweep of all particles and all cells.

It would not be very difficult to implement this scheme for all three of the
functions in this program that have this issue.


<a name="Morestatisticsaboutelectrons"></a><h4> More statistics about electrons </h4>


The program already computes the fraction of the electrons that leave the
domain through the hole in the anode. But there are other quantities one might be
interested in. For example, the average velocity of these particles. It would
not be very difficult to obtain each particle's velocity from its properties,
in the same way as we do in the `move_particles()` function, and compute
statistics from it.


<a name="Abettersynchronizedvisualization"></a><h4> A better-synchronized visualization </h4>


As discussed above, there is a varying time difference between different frames
of the video because we create output for every time step. A better way to
create movies would be to generate a new output file in fixed time intervals,
regardless of how many time steps lie between each such point.


<a name="Abettertimeintegrator"></a><h4> A better time integrator </h4>


The problem we are considering in this program is a coupled, multiphysics
problem. But the way we solve it is by first computing the (electric) potential
field and then update the particle locations. This is what is called an
"operator-splitting method", a concept we will investigate in more detail
in step-58.

While it is awkward to think of a way to solve this problem that does not involve
splitting the problem into a PDE piece and a particles piece, one
*can* (and probably should!) think of a better way to update the particle
locations. Specifically, the equations we use to update the particle location
are
@f{align*}{
  \frac{{\mathbf v}_i^{(n)}-{\mathbf v}_i^{(n-1)}}{\Delta t} &= \frac{e\nabla V^{(n)}}{m}
  \\
  \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} &= {\mathbf v}_i^{(n)}.
@f}
This corresponds to a simple forward-Euler time discretization -- a method of
first order accuracy in the time step size $\Delta t$ that we know we should
avoid because we can do better. Rather, one might want to consider a scheme such
as the
[leapfrog scheme](https://en.wikipedia.org/wiki/Leapfrog_integration)
or more generally
[symplectic integrators](https://en.wikipedia.org/wiki/Symplectic_integrator)
such as the
[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_integration).


<a name="Parallelization"></a><h4> Parallelization </h4>


In release mode, the program runs in about 3.5 minutes on one of the author's
laptops at the time of writing this. That's acceptable. But what if we wanted
to make the simulation three-dimensional? If we wanted to not use a maximum
of around 100 particles at any given time (as happens with the parameters
used here) but 100,000? If we needed a substantially finer mesh?

In those cases, one would want to run the program not just on a single processor,
but in fact on as many as one has available. This requires parallelization
both the PDE solution as well as over particles. In practice, while there
are substantial challenges to making this efficient and scale well, these
challenges are all addressed in deal.II itself. For example, step-40 shows
how to parallelize the finite element part, and step-70 shows how one can
then also parallelize the particles part.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-19.cc"
*/
