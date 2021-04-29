/**
@page step_68 The step-68 tutorial program
This tutorial depends on step-19.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Simulationofthemotionofmasslesstracerparticlesinavorticalflow">Simulation of the motion of massless tracer particles in a vortical flow</a>
        <li><a href="#ParticlesindealII">Particles in deal.II</a>
        <li><a href="#Challengesrelatedtodistributedparticlesimulations">Challenges related to distributed particle simulations</a>
      <ul>
        <li><a href="#Parallelparticlegeneration">Parallel particle generation</a>
        <li><a href="#Particleexchange">Particle exchange</a>
        <li><a href="#Balancingmeshandparticleload">Balancing mesh and particle load</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Runtimeparameterhandling">Run-time parameter handling</a>
        <li><a href="#Velocityprofile">Velocity profile</a>
        <li><a href="#ThecodeParticleTrackingcodeclassdeclaration">The <code>ParticleTracking</code> class declaration</a>
        <li><a href="#ThecodePatricleTrackingcodeclassimplementation">The <code>PatricleTracking</code> class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#Cellweight">Cell weight</a>
        <li><a href="#Particlesgeneration">Particles generation</a>
        <li><a href="#BackgroundDOFsandinterpolation">Background DOFs and interpolation</a>
        <li><a href="#Timeintegrationofthetrajectories">Time integration of the trajectories</a>
        <li><a href="#Dataoutput">Data output</a>
        <li><a href="#Runningthesimulation">Running the simulation</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Motionoftheparticles"> Motion of the particles </a>
        <li><a href="#Dynamicloadbalancing"> Dynamic load balancing </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program was contributed by
Bruno Blais (Polytechnique Montréal),
Toni El Geitani Nehme (Polytechnique Montréal),
Rene Gassmöller (University of California Davis),
and Peter Munch (Technical University of Munich and Helmholtz-Zentrum Geesthacht).
Bruno Blais was supported by NSERC Discovery grant
RGPIN-2020-04510, by Compute Canada and Calcul Québec.
</i>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Simulationofthemotionofmasslesstracerparticlesinavorticalflow"></a><h3>Simulation of the motion of massless tracer particles in a vortical flow</h3>


Particles play an important part in numerical models for a large
 number of applications. Particles are routinely used
 as massless tracers to visualize the dynamic of a transient flow. They
 can also play an intrinsic role as part of a more complex finite element
 model, as is the case for the Particle-In-Cell (PIC) method @cite GLHPW2018
 or they can even be used to simulate the motion of granular matter, as in
 the Discrete Element Method (DEM) @cite Blais2019. In the case
 of DEM, the resulting model is not related to the finite element method anymore,
 but just leads to a system of ordinary differential equation which describes
 the motion of the particles and the dynamic of their collisions. All of
 these models can be built using deal.II's particle handling capabilities.

In the present step, we use particles as massless tracers to illustrate
the dynamic of a vortical flow. Since the particles are massless tracers,
the position of each particle $i$ is described by the
following ordinary differential equation (ODE):
@f[
\frac{d \textbf{x}_i}{dt} =\textbf{u}(\textbf{x}_i)
@f]

where $\textbf{x}_i$ is the position of particle $i$ and $\textbf{u}(\textbf{x}_i)$ the flow velocity at its position.
In the present step, this ODE is solved using the explicit Euler method. The resulting scheme is:
@f[
\textbf{x}_{i}^{n+1} = \textbf{x}_{i}^{n} + \Delta t \; \textbf{u}(\textbf{x}_{i}^{n})
@f]

where $\textbf{x}_{i}^{n+1}$ and $\textbf{x}_{i}^{n}$ are the position
of particle $i$ at time $t+\Delta t$ and $t$, respectively and where $\Delta t$
is the time step. In the present step, the velocity at the location of particles
is obtained in two different fashions:
- By evaluating the velocity function at the location of the particles;
- By evaluating the velocity function on a background triangulation and, using
a  finite element support, interpolating at the position of the particle.

The first approach is not practical, since the velocity profile
is generally not known analytically. The second approach, based on interpolating a solution
at the position of the particles, mimics exactly what would be done in a
realistic computational fluid dynamic simulation, and this follows the way we have also evaluated
the finite element solution at particle locations in step-19. In this step, we illustrate both strategies.

We note that much greater accuracy could be obtained by using a fourth
order Runge-Kutta method or another appropriate scheme for the time integration
of the motion of the particles.  Implementing a more advanced time-integration scheme
would be a straightforward extension of this step.

<a name="ParticlesindealII"></a><h3>Particles in deal.II</h3>


In deal.II, Particles::Particle are very simple and flexible entities that can be used
to build PIC, DEM or any type of particle-based models. Particles have a location
in real space, a location in the reference space of the element in which they
are located and a unique ID. In the majority of cases, simulations that include
particles require a significant number of them. Thus, it becomes interesting
to handle all particles through an entity which agglomerates all particles.
In deal.II, this is achieved through the use of the Particles::ParticleHandler class.

By default, particles do not have a diameter,
a mass or any other physical properties which we would generally expect of physical particles. However, through
a ParticleHandler, particles have access to a Particles::PropertyPool. This PropertyPool is
an array which can be used to store an arbitrary number of properties
associated with the particles. Consequently, users can build their own
particle solver and attribute the desired properties to the particles (e.g., mass, charge,
diameter, temperature, etc.). In the present tutorial, this is used to
store the value of the fluid velocity and the process id to which the particles
belong.

<a name="Challengesrelatedtodistributedparticlesimulations"></a><h3>Challenges related to distributed particle simulations</h3>


Although the present step is not computationally intensive, simulations that
include many particles can be computationally demanding and require parallelization.
The present step showcases the distributed parallel capabilities of deal.II for particles.
In general, there are three main challenges
that specifically arise in parallel distributed simulations that include particles:
- Generating the particles on the distributed triangulation;
- Exchanging the particles that leave local domains between the processors;
- Load balancing the simulation so that every processor has a similar computational load.
These challenges and their solution in deal.II have been discussed in more detail in
@cite GLHPW2018, but we will summarize them below.

There are of course also questions on simply setting up a code that uses particles. These have largely already been
addressed in step-19. Some more advanced techniques will also be discussed in step-70.

<a name="Parallelparticlegeneration"></a><h4>Parallel particle generation</h4>


Generating distributed particles in a scalable way is not straightforward since
the processor to which they belong must first be identified before the cell in
which they are located is found.  deal.II provides numerous capabilities to
generate particles through the Particles::Generator namespace.  Some of these
particle generators create particles only on the locally owned subdomain. For
example, Particles::Generators::regular_reference_locations() creates particles
at the same reference locations within each cell of the local subdomain and
Particles::Generators::probabilistic_locations() uses a globally defined probability
density function to determine how many and where to generate particles locally.

In other situations, such as the present step, particles must be generated at
specific locations on cells that may be owned only by a subset of the processors.
In  most of these situations, the insertion of the particles is done for a very
limited number of time-steps and, consequently, does not constitute a large
portion of the computational cost. For these occasions, deal.II provides
convenient Particles::Generators that can globally insert the particles even if
the particle is not located in a cell owned by the parallel process on which the call to create the particle is initiated. The
generators first locate on which subdomain the particles are situated, identify
in which cell they are located and exchange the necessary information among
the processors to ensure that the particle is generated with the right
properties. Consequently, this type of particle generation can be communication
intensive. The Particles::Generators::dof_support_points and the
Particles::Generators::quadrature_points generate particles using a
triangulation and the points of an associated DoFHandler or quadrature
respectively. The triangulation that is used to generate the particles can be
the same triangulation that is used for the background mesh, in which case these
functions are very similar to the
Particles::Generators::regular_reference_locations() function described in the
previous paragraph. However, the triangulation used to generate particles can
also be different (non-matching) from the triangulation of the background grid,
which is useful to generate particles in particular shapes (as in this
example), or to transfer information between two different computational grids
(as in step-70).  Furthermore, the Particles::ParticleHandler class provides the
Particles::ParticleHandler::insert_global_particles() function which enables the
global insertion of particles from a vector of arbitrary points and a global
vector of bounding boxes. In the present step, we use the
Particles::Generators::quadrature_points() function on a non-matching triangulation to
insert particles located at positions in the shape of a disk.

<a name="Particleexchange"></a><h4>Particle exchange</h4>


As particles move around in parallel distributed computations they may leave
the locally owned subdomain and need to be transferred to their new owner
processes. This situation can arise in two very different ways: First, if the
previous owning process knows the new owner of the particles that were lost
(for example because the particles moved from the locally owned cell of one processor
into an adjacent ghost cells of a distributed
triangulation) then the transfer can be handled efficiently as a point-to-point
communication between each process and the new owners. This transfer happens
automatically whenever particles are sorted into their new cells. Secondly,
the previous owner may not know to which process the particle has moved. In
this case the particle is discarded by default, as a global search for the
owner can be expensive. step-19 shows how such a discarded particle can still
be collected, interpreted, and potentially reinserted by the user. In the
present example we prevent the second case by imposing a CFL criterion on the
timestep to ensure particles will at most move into the ghost layer of the
local process and can therefore be send to neighboring processes automatically.

<a name="Balancingmeshandparticleload"></a><h4>Balancing mesh and particle load</h4>


The last challenge that arises in parallel distributed computations using
particles is to balance the computational load between work that is done on the
grid, for example solving the finite-element problem, and the work that is done
on the particles, for example advecting the particles or computing the forces
between particles or between particles and grid. By default, for example in
step-40, deal.II distributes the background mesh as evenly as possible between
the available processes, that is it balances the number of cells on each
process. However, if some cells own many more particles than other cells, or if
the particles of one cell are much more computationally expensive than the
particles in other cells, then this problem no longer scales efficiently (for a
discussion of what we consider "scalable" programs, see
@ref GlossParallelScaling "this glossary entry"). Thus, we have to apply a form of
"load balancing", which means we estimate the computational load that is
associated with each cell and its particles. Repartitioning the mesh then
accounts for this combined computational load instead of the simplified
assumption of the number of cells @cite GLHPW2018.

In this section we only discussed the particle-specific challenges in distributed
computation. Parallel challenges that particles share with
finite-element solutions (parallel output, data transfer during mesh
refinement) can be addressed with the solutions found for
finite-element problems already discussed in other examples.

<a name="Thetestcase"></a><h3>The testcase</h3>


In the present step, we use particles as massless tracers to illustrate
the dynamics of a particular vortical flow: the Rayleigh--Kothe vortex. This flow pattern
is generally used as a complex test case for interface tracking methods
(e.g., volume-of-fluid and level set approaches) since
it leads to strong rotation and elongation of the fluid @cite Blais2013.

The stream function $\Psi$ of this Rayleigh-Kothe vortex is defined as:

@f[
\Psi = \frac{1}{\pi} \sin^2 (\pi x) \sin^2 (\pi y) \cos \left( \pi \frac{t}{T} \right)
@f]
where $T$ is half the period of the flow. The velocity profile in 2D ($\textbf{u}=[u,v]^T$) is :
@f{eqnarray*}
   u &=&  - \frac{\partial\Psi}{\partial y} = -2 \sin^2 (\pi x) \sin (\pi y) \cos (\pi y)  \cos \left( \pi \frac{t}{T} \right)\\
   v &=&  \frac{\partial\Psi}{\partial x} = 2 \cos(\pi x) \sin(\pi x) \sin^2 (\pi y) \cos \left( \pi \frac{t}{T} \right)
@f}

The velocity profile is illustrated in the following animation:

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/m6hQm7etji8"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly

It can be seen that this velocity reverses periodically due to the term
$\cos \left( \pi \frac{t}{T} \right)$ and that material will end up at its
starting position after every period of length $t=2T$. We will run this tutorial
program for exactly one period and compare the final particle location to the
initial location to illustrate this flow property. This example uses the testcase
to produce two models that handle the particles
slightly differently. The first model prescribes the exact analytical velocity
solution as the velocity for each particle. Therefore in this model there is no
error in the assigned velocity to the particles, and any deviation of particle
positions from the analytical position at a given time results from the error
in solving the equation of motion for the particle inexactly, using a time stepping method. In the second model the
analytical velocity field is first interpolated to a finite-element vector
space (to simulate the case that the velocity was obtained from solving a
finite-element problem, in the same way as the ODE for each particle in step-19 depends on a finite element
solution). This finite-element "solution" is then evaluated at
the locations of the particles to solve their equation of motion. The
difference between the two cases allows to assess whether the chosen
finite-element space is sufficiently accurate to advect the particles with the
optimal convergence rate of the chosen particle advection scheme, a question
that is important in practice to determine the accuracy of the combined
algorithm (see e.g. @cite Gassmoller2019).
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
 * 
 * @code
 * #include <deal.II/base/bounding_box.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/base/mpi.h>
 * #include <deal.II/base/parameter_acceptor.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/distributed/cell_weights.h>
 * #include <deal.II/distributed/solution_transfer.h>
 * #include <deal.II/distributed/tria.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * From the following include file we import the ParticleHandler class
 * that allows you to manage
 * a collection of particles (objects of type Particles::Particle), representing
 * a collection of points with some attached properties (e.g., an id) floating
 * on a parallel::distributed::Triangulation. The methods and classes in the
 * namespace Particles allows one to easily implement Particle-In-Cell methods
 * and particle tracing on distributed triangulations:
 * 
 * @code
 * #include <deal.II/particles/particle_handler.h>
 * 
 * @endcode
 * 
 * We import the particles generator
 * which allow us to insert the particles. In the present step, the particle
 * are globally inserted using a non-matching hyper-shell triangulation:
 * 
 * @code
 * #include <deal.II/particles/generators.h>
 * 
 * @endcode
 * 
 * Since the particles do not form a triangulation, they have their
 * own specific DataOut class which will enable us to write them
 * to commonly used parallel vtu format (or any number of other file formats):
 * 
 * @code
 * #include <deal.II/particles/data_out.h>
 * 
 * #include <cmath>
 * #include <iostream>
 * 
 * 
 * 
 * namespace Step68
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameterhandling"></a> 
 * <h3>Run-time parameter handling</h3>
 * 

 * 
 * Similarly to what is done in step-60, we set up a class that holds
 * all the parameters of our problem and derive it from the ParameterAcceptor
 * class to simplify the management and creation of parameter files.
 *   

 * 
 * The ParameterAcceptor paradigm requires all parameters to be writable by
 * the ParameterAcceptor methods. In order to avoid bugs that would be very
 * difficult to track down (such as writing things like `if (time = 0)`
 * instead of `if(time == 0)`), we declare all the parameters in an external
 * class, which is initialized before the actual `ParticleTracking` class, and
 * pass it to the main class as a `const` reference.
 *   

 * 
 * The constructor of the class is responsible for the connection between the
 * members of this class and the corresponding entries in the
 * ParameterHandler. Thanks to the use of the
 * ParameterHandler::add_parameter() method, this connection is trivial, but
 * requires all members of this class to be writable.
 * 
 * @code
 *   class ParticleTrackingParameters : public ParameterAcceptor
 *   {
 *   public:
 *     ParticleTrackingParameters();
 * 
 * @endcode
 * 
 * This class consists largely of member variables that
 * describe the details of the particle tracking simulation and its
 * discretization. The following parameters are about where output should
 * written to, the spatial discretization of the velocity (the default is
 * $Q_1$), the time step and the output frequency (how many time steps
 * should elapse before we generate graphical output again):
 * 
 * @code
 *     std::string output_directory = "./";
 * 
 *     unsigned int velocity_degree       = 1;
 *     double       time_step             = 0.002;
 *     double       final_time            = 4.0;
 *     unsigned int output_frequency      = 10;
 *     unsigned int repartition_frequency = 5;
 * 
 * @endcode
 * 
 * We allow every grid to be refined independently. In this tutorial, no
 * physics is resolved on the fluid grid, and its velocity is calculated
 * analytically.
 * 
 * @code
 *     unsigned int fluid_refinement              = 4;
 *     unsigned int particle_insertion_refinement = 3;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * There remains the task of declaring what run-time parameters we can accept
 * in input files. Since we have a very limited number of parameters, all
 * parameters are declared in the same section.
 * 
 * @code
 *   ParticleTrackingParameters::ParticleTrackingParameters()
 *     : ParameterAcceptor("Particle Tracking Problem/")
 *   {
 *     add_parameter(
 *       "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1));
 * 
 *     add_parameter("Output frequency",
 *                   output_frequency,
 *                   "Iteration frequency at which output results are written",
 *                   prm,
 *                   Patterns::Integer(1));
 * 
 *     add_parameter("Repartition frequency",
 *                   repartition_frequency,
 *                   "Iteration frequency at which the mesh is load balanced",
 *                   prm,
 *                   Patterns::Integer(1));
 * 
 *     add_parameter("Output directory", output_directory);
 * 
 *     add_parameter("Time step", time_step, "", prm, Patterns::Double());
 * 
 *     add_parameter("Final time",
 *                   final_time,
 *                   "End time of the simulation",
 *                   prm,
 *                   Patterns::Double());
 * 
 *     add_parameter("Fluid refinement",
 *                   fluid_refinement,
 *                   "Refinement level of the fluid domain",
 *                   prm,
 *                   Patterns::Integer(0));
 * 
 *     add_parameter(
 *       "Particle insertion refinement",
 *       particle_insertion_refinement,
 *       "Refinement of the volumetric mesh used to insert the particles",
 *       prm,
 *       Patterns::Integer(0));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Velocityprofile"></a> 
 * <h3>Velocity profile</h3>
 * 

 * 
 * The velocity profile is provided as a Function object.
 * This function is hard-coded within
 * the example.
 * 
 * @code
 *   template <int dim>
 *   class Vortex : public Function<dim>
 *   {
 *   public:
 *     Vortex()
 *       : Function<dim>(dim)
 *     {}
 * 
 * 
 *     virtual void vector_value(const Point<dim> &point,
 *                               Vector<double> &  values) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * The velocity profile for the Rayleigh-Kothe vertex is time-dependent.
 * Consequently, the current time in the
 * simulation (t) must be gathered from the Function object.
 * 
 * @code
 *   template <int dim>
 *   void Vortex<dim>::vector_value(const Point<dim> &point,
 *                                  Vector<double> &  values) const
 *   {
 *     const double T = 4;
 *     const double t = this->get_time();
 * 
 *     const double px = numbers::PI * point(0);
 *     const double py = numbers::PI * point(1);
 *     const double pt = numbers::PI / T * t;
 * 
 *     values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py);
 *     values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px);
 *     if (dim == 3)
 *       {
 *         values[2] = 0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeParticleTrackingcodeclassdeclaration"></a> 
 * <h3>The <code>ParticleTracking</code> class declaration</h3>
 * 

 * 
 * We are now ready to introduce the main class of our tutorial program.
 * 
 * @code
 *   template <int dim>
 *   class ParticleTracking
 *   {
 *   public:
 *     ParticleTracking(const ParticleTrackingParameters &par,
 *                      const bool                        interpolated_velocity);
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * This function is responsible for the initial
 * generation of the particles on top of the background grid.
 * 
 * @code
 *     void generate_particles();
 * 
 * @endcode
 * 
 * When the velocity profile is interpolated to the position of the
 * particles, it must first be stored using degrees of freedom.
 * Consequently, as is the case for other parallel case (e.g. step-40) we
 * initialize the degrees of freedom on the background grid.
 * 
 * @code
 *     void setup_background_dofs();
 * 
 * @endcode
 * 
 * In one of the test cases, the function is mapped to the background grid
 * and a finite element interpolation is used to calculate the velocity
 * at the particle location. This function calculates the value of the
 * function at the support point of the triangulation.
 * 
 * @code
 *     void interpolate_function_to_field();
 * 
 * @endcode
 * 
 * The next two functions are responsible for carrying out step of explicit
 * Euler time integration for the cases where the velocity field is
 * interpolated at the positions of the particles or calculated
 * analytically, respectively.
 * 
 * @code
 *     void euler_step_interpolated(const double dt);
 *     void euler_step_analytical(const double dt);
 * 
 * @endcode
 * 
 * The `cell_weight()` function indicates to the triangulation how much
 * computational work is expected to happen on this cell, and consequently
 * how the domain needs to be partitioned so that every MPI rank receives a
 * roughly equal amount of work (potentially not an equal number of cells).
 * While the function is called from the outside, it is connected to the
 * corresponding signal from inside this class, therefore it can be
 * `private`.
 * 
 * @code
 *     unsigned int cell_weight(
 *       const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *         &cell,
 *       const typename parallel::distributed::Triangulation<dim>::CellStatus
 *         status) const;
 * 
 * @endcode
 * 
 * The following two functions are responsible for outputting the simulation
 * results for the particles and for the velocity profile on the background
 * mesh, respectively.
 * 
 * @code
 *     void output_particles(const unsigned int it);
 *     void output_background(const unsigned int it);
 * 
 * @endcode
 * 
 * The private members of this class are similar to other parallel deal.II
 * examples. The parameters are stored as a `const` member. It is important
 * to note that we keep the `Vortex` class as a member since its time
 * must be modified as the simulation proceeds.
 * 

 * 
 * 
 * @code
 *     const ParticleTrackingParameters &par;
 * 
 *     MPI_Comm                                  mpi_communicator;
 *     parallel::distributed::Triangulation<dim> background_triangulation;
 *     Particles::ParticleHandler<dim>           particle_handler;
 * 
 *     DoFHandler<dim>                            fluid_dh;
 *     FESystem<dim>                              fluid_fe;
 *     MappingQ1<dim>                             mapping;
 *     LinearAlgebra::distributed::Vector<double> velocity_field;
 * 
 *     Vortex<dim> velocity;
 * 
 *     ConditionalOStream pcout;
 * 
 *     bool interpolated_velocity;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePatricleTrackingcodeclassimplementation"></a> 
 * <h3>The <code>PatricleTracking</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * The constructors and destructors are rather trivial. They are very similar
 * to what is done in step-40. We set the processors we want to work on
 * to all machines available (`MPI_COMM_WORLD`) and
 * initialize the <code>pcout</code> variable to only allow processor zero
 * to output anything to the standard output.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
 *                                           const bool interpolated_velocity)
 *     : par(par)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , background_triangulation(mpi_communicator)
 *     , fluid_dh(background_triangulation)
 *     , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *     , interpolated_velocity(interpolated_velocity)
 * 
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Cellweight"></a> 
 * <h4>Cell weight</h4>
 * 

 * 
 * This function is the key component that allow us to dynamically balance the
 * computational load for this example. The function attributes a weight to
 * every cell that represents the computational work on this cell. Here the
 * majority of work is expected to happen on the particles, therefore the
 * return value of this function (representing "work for this cell") is
 * calculated based on the number of particles in the current cell.
 * The function is
 * connected to the cell_weight() signal inside the triangulation, and will be
 * called once per cell, whenever the triangulation repartitions the domain
 * between ranks (the connection is created inside the
 * generate_particles() function of this class).
 * 
 * @code
 *   template <int dim>
 *   unsigned int ParticleTracking<dim>::cell_weight(
 *     const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *       &                                                                  cell,
 *     const typename parallel::distributed::Triangulation<dim>::CellStatus status)
 *     const
 *   {
 * @endcode
 * 
 * We do not assign any weight to cells we do not own (i.e., artificial
 * or ghost cells)
 * 
 * @code
 *     if (!cell->is_locally_owned())
 *       return 0;
 * 
 * @endcode
 * 
 * This determines how important particle work is compared to cell
 * work (by default every cell has a weight of 1000).
 * We set the weight per particle much higher to indicate that
 * the particle load is the only one that is important to distribute the
 * cells in this example. The optimal value of this number depends on the
 * application and can range from 0 (cheap particle operations,
 * expensive cell operations) to much larger than 1000 (expensive
 * particle operations, cheap cell operations, like presumed in this
 * example).
 * 
 * @code
 *     const unsigned int particle_weight = 10000;
 * 
 * @endcode
 * 
 * This example does not use adaptive refinement, therefore every cell
 * should have the status `CELL_PERSIST`. However this function can also
 * be used to distribute load during refinement, therefore we consider
 * refined or coarsened cells as well.
 * 
 * @code
 *     if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
 *         status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
 *       {
 *         const unsigned int n_particles_in_cell =
 *           particle_handler.n_particles_in_cell(cell);
 *         return n_particles_in_cell * particle_weight;
 *       }
 *     else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
 *       {
 *         unsigned int n_particles_in_cell = 0;
 * 
 *         for (unsigned int child_index = 0; child_index < cell->n_children();
 *              ++child_index)
 *           n_particles_in_cell +=
 *             particle_handler.n_particles_in_cell(cell->child(child_index));
 * 
 *         return n_particles_in_cell * particle_weight;
 *       }
 * 
 *     Assert(false, ExcInternalError());
 *     return 0;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Particlesgeneration"></a> 
 * <h4>Particles generation</h4>
 * 

 * 
 * This function generates the tracer particles and the background
 * triangulation on which these particles evolve.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::generate_particles()
 *   {
 * @endcode
 * 
 * We create a hyper cube triangulation which we globally refine. This
 * triangulation covers the full trajectory of the particles.
 * 
 * @code
 *     GridGenerator::hyper_cube(background_triangulation, 0, 1);
 *     background_triangulation.refine_global(par.fluid_refinement);
 * 
 * @endcode
 * 
 * In order to consider the particles when repartitioning the triangulation
 * the algorithm needs to know three things:
 *     

 * 
 * 1. How much weight to assign to each cell (how many particles are in
 * there);
 * 2. How to pack the particles before shipping data around;
 * 3. How to unpack the particles after repartitioning.
 *     

 * 
 * We attach the correct functions to the signals inside
 * parallel::distributed::Triangulation. These signal will be called every
 * time the repartition() function is called. These connections only need to
 * be created once, so we might as well have set them up in the constructor
 * of this class, but for the purpose of this example we want to group the
 * particle related instructions.
 * 
 * @code
 *     background_triangulation.signals.cell_weight.connect(
 *       [&](
 *         const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *           &cell,
 *         const typename parallel::distributed::Triangulation<dim>::CellStatus
 *           status) -> unsigned int { return this->cell_weight(cell, status); });
 * 
 *     background_triangulation.signals.pre_distributed_repartition.connect(
 *       [this]() { this->particle_handler.register_store_callback_function(); });
 * 
 *     background_triangulation.signals.post_distributed_repartition.connect(
 *       [&]() { this->particle_handler.register_load_callback_function(false); });
 * 
 * @endcode
 * 
 * This initializes the background triangulation where the particles are
 * living and the number of properties of the particles.
 * 
 * @code
 *     particle_handler.initialize(background_triangulation, mapping, 1 + dim);
 * 
 * @endcode
 * 
 * We create a particle triangulation which is solely used to generate
 * the points which will be used to insert the particles. This
 * triangulation is a hyper shell which is offset from the
 * center of the simulation domain. This will be used to generate a
 * disk filled with particles which will allow an easy monitoring
 * of the motion due to the vortex.
 * 
 * @code
 *     Point<dim> center;
 *     center[0] = 0.5;
 *     center[1] = 0.75;
 *     if (dim == 3)
 *       center[2] = 0.5;
 * 
 *     const double outer_radius = 0.15;
 *     const double inner_radius = 0.01;
 * 
 *     parallel::distributed::Triangulation<dim> particle_triangulation(
 *       MPI_COMM_WORLD);
 * 
 *     GridGenerator::hyper_shell(
 *       particle_triangulation, center, inner_radius, outer_radius, 6);
 *     particle_triangulation.refine_global(par.particle_insertion_refinement);
 * 
 * @endcode
 * 
 * We generate the necessary bounding boxes for the particles generator.
 * These bounding boxes are required to quickly identify in which
 * process's subdomain the inserted particle lies, and which cell owns it.
 * 
 * @code
 *     const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
 *       background_triangulation, IteratorFilters::LocallyOwnedCell());
 *     const auto global_bounding_boxes =
 *       Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);
 * 
 * @endcode
 * 
 * We generate an empty vector of properties. We will attribute the
 * properties to the particles once they are generated.
 * 
 * @code
 *     std::vector<std::vector<double>> properties(
 *       particle_triangulation.n_locally_owned_active_cells(),
 *       std::vector<double>(dim + 1, 0.));
 * 
 * @endcode
 * 
 * We generate the particles at the position of a single
 * point quadrature. Consequently, one particle will be generated
 * at the centroid of each cell.
 * 
 * @code
 *     Particles::Generators::quadrature_points(particle_triangulation,
 *                                              QMidpoint<dim>(),
 *                                              global_bounding_boxes,
 *                                              particle_handler,
 *                                              mapping,
 *                                              properties);
 * 
 *     pcout << "Number of particles inserted: "
 *           << particle_handler.n_global_particles() << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BackgroundDOFsandinterpolation"></a> 
 * <h4>Background DOFs and interpolation</h4>
 * 

 * 
 * This function sets up the background degrees of freedom used for the
 * velocity interpolation and allocates the field vector where the entire
 * solution of the velocity field is stored.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::setup_background_dofs()
 *   {
 *     fluid_dh.distribute_dofs(fluid_fe);
 *     const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
 *     IndexSet       locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);
 * 
 *     velocity_field.reinit(locally_owned_dofs,
 *                           locally_relevant_dofs,
 *                           mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This function takes care of interpolating the
 * vortex velocity field to the field vector. This is achieved rather easily
 * by using the VectorTools::interpolate() function.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::interpolate_function_to_field()
 *   {
 *     velocity_field.zero_out_ghost_values();
 *     VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
 *     velocity_field.update_ghost_values();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Timeintegrationofthetrajectories"></a> 
 * <h4>Time integration of the trajectories</h4>
 * 

 * 
 * We integrate the particle trajectories
 * using an analytically defined velocity field. This demonstrates a
 * relatively trivial usage of the particles.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::euler_step_analytical(const double dt)
 *   {
 *     const unsigned int this_mpi_rank =
 *       Utilities::MPI::this_mpi_process(mpi_communicator);
 *     Vector<double> particle_velocity(dim);
 * 
 * @endcode
 * 
 * Looping over all particles in the domain using a particle iterator
 * 
 * @code
 *     for (auto &particle : particle_handler)
 *       {
 * @endcode
 * 
 * We calculate the velocity of the particles using their current
 * location.
 * 
 * @code
 *         Point<dim> particle_location = particle.get_location();
 *         velocity.vector_value(particle_location, particle_velocity);
 * 
 * @endcode
 * 
 * This updates the position of the particles and sets the old position
 * equal to the new position of the particle.
 * 
 * @code
 *         for (int d = 0; d < dim; ++d)
 *           particle_location[d] += particle_velocity[d] * dt;
 * 
 *         particle.set_location(particle_location);
 * 
 * @endcode
 * 
 * We store the processor id (a scalar) and the particle velocity (a
 * vector) in the particle properties. In this example, this is done
 * purely for visualization purposes.
 * 
 * @code
 *         ArrayView<double> properties = particle.get_properties();
 *         for (int d = 0; d < dim; ++d)
 *           properties[d] = particle_velocity[d];
 *         properties[dim] = this_mpi_rank;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In contrast to the previous function in this function we
 * integrate the particle trajectories by interpolating the value of
 * the velocity field at the degrees of freedom to the position of
 * the particles.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::euler_step_interpolated(const double dt)
 *   {
 *     Vector<double> local_dof_values(fluid_fe.dofs_per_cell);
 * 
 * @endcode
 * 
 * We loop over all the local particles. Although this could be achieved
 * directly by looping over all the cells, this would force us
 * to loop over numerous cells which do not contain particles.
 * Rather, we loop over all the particles, but, we get the reference
 * of the cell in which the particle lies and then loop over all particles
 * within that cell. This enables us to gather the values of the velocity
 * out of the `velocity_field` vector once and use them for all particles
 * that lie within the cell.
 * 
 * @code
 *     auto particle = particle_handler.begin();
 *     while (particle != particle_handler.end())
 *       {
 *         const auto cell =
 *           particle->get_surrounding_cell(background_triangulation);
 *         const auto dh_cell =
 *           typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);
 * 
 *         dh_cell->get_dof_values(velocity_field, local_dof_values);
 * 
 * @endcode
 * 
 * Next, compute the velocity at the particle locations by evaluating
 * the finite element solution at the position of the particles.
 * This is essentially an optimized version of the particle advection
 * functionality in step 19, but instead of creating quadrature
 * objects and FEValues objects for each cell, we do the
 * evaluation by hand, which is somewhat more efficient and only
 * matters for this tutorial, because the particle work is the
 * dominant cost of the whole program.
 * 
 * @code
 *         const auto pic = particle_handler.particles_in_cell(cell);
 *         Assert(pic.begin() == particle, ExcInternalError());
 *         for (auto &p : pic)
 *           {
 *             const Point<dim> reference_location = p.get_reference_location();
 *             Tensor<1, dim>   particle_velocity;
 *             for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
 *               {
 *                 const auto comp_j = fluid_fe.system_to_component_index(j);
 * 
 *                 particle_velocity[comp_j.first] +=
 *                   fluid_fe.shape_value(j, reference_location) *
 *                   local_dof_values[j];
 *               }
 * 
 *             Point<dim> particle_location = particle->get_location();
 *             for (int d = 0; d < dim; ++d)
 *               particle_location[d] += particle_velocity[d] * dt;
 *             p.set_location(particle_location);
 * 
 * @endcode
 * 
 * Again, we store the particle velocity and the processor id in the
 * particle properties for visualization purposes.
 * 
 * @code
 *             ArrayView<double> properties = p.get_properties();
 *             for (int d = 0; d < dim; ++d)
 *               properties[d] = particle_velocity[d];
 * 
 *             properties[dim] =
 *               Utilities::MPI::this_mpi_process(mpi_communicator);
 * 
 *             ++particle;
 *           }
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Dataoutput"></a> 
 * <h4>Data output</h4>
 * 

 * 
 * The next two functions take care of writing both the particles
 * and the background mesh to vtu with a pvtu record. This ensures
 * that the simulation results can be visualized when the simulation is
 * launched in parallel.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::output_particles(const unsigned int it)
 *   {
 *     Particles::DataOut<dim, dim> particle_output;
 * 
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     solution_names.push_back("process_id");
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     particle_output.build_patches(particle_handler,
 *                                   solution_names,
 *                                   data_component_interpretation);
 *     const std::string output_folder(par.output_directory);
 *     const std::string file_name(interpolated_velocity ?
 *                                   "interpolated-particles" :
 *                                   "analytical-particles");
 * 
 *     pcout << "Writing particle output file: " << file_name << "-" << it
 *           << std::endl;
 * 
 *     particle_output.write_vtu_with_pvtu_record(
 *       output_folder, file_name, it, mpi_communicator, 6);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void ParticleTracking<dim>::output_background(const unsigned int it)
 *   {
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 * 
 *     DataOut<dim> data_out;
 * 
 * @endcode
 * 
 * Attach the solution data to data_out object
 * 
 * @code
 *     data_out.attach_dof_handler(fluid_dh);
 *     data_out.add_data_vector(velocity_field,
 *                              solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 *     Vector<float> subdomain(background_triangulation.n_active_cells());
 *     for (unsigned int i = 0; i < subdomain.size(); ++i)
 *       subdomain(i) = background_triangulation.locally_owned_subdomain();
 *     data_out.add_data_vector(subdomain, "subdomain");
 * 
 *     data_out.build_patches(mapping);
 * 
 *     const std::string output_folder(par.output_directory);
 *     const std::string file_name("background");
 * 
 *     pcout << "Writing background field file: " << file_name << "-" << it
 *           << std::endl;
 * 
 *     data_out.write_vtu_with_pvtu_record(
 *       output_folder, file_name, it, mpi_communicator, 6);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Runningthesimulation"></a> 
 * <h4>Running the simulation</h4>
 * This function orchestrates the entire simulation. It is very similar
 * to the other time dependent tutorial programs -- take step-21 or step-26 as
 * an example. Note that we use the DiscreteTime class to monitor the time,
 * the time-step and the step-number. This function is relatively
 * straightforward.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::run()
 *   {
 *     DiscreteTime discrete_time(0, par.final_time, par.time_step);
 * 
 *     generate_particles();
 * 
 *     pcout << "Repartitioning triangulation after particle generation"
 *           << std::endl;
 *     background_triangulation.repartition();
 * 
 * @endcode
 * 
 * We set the initial property of the particles by doing an
 * explicit Euler iteration with a time-step of 0 both in the case
 * of the analytical and the interpolated approach.
 * 
 * @code
 *     if (interpolated_velocity)
 *       {
 *         setup_background_dofs();
 *         interpolate_function_to_field();
 *         euler_step_interpolated(0.);
 *       }
 *     else
 *       euler_step_analytical(0.);
 * 
 *     output_particles(discrete_time.get_step_number());
 *     if (interpolated_velocity)
 *       output_background(discrete_time.get_step_number());
 * 
 * @endcode
 * 
 * The particles are advected by looping over time.
 * 
 * @code
 *     while (!discrete_time.is_at_end())
 *       {
 *         discrete_time.advance_time();
 *         velocity.set_time(discrete_time.get_previous_time());
 * 
 *         if ((discrete_time.get_step_number() % par.repartition_frequency) == 0)
 *           {
 *             background_triangulation.repartition();
 *             if (interpolated_velocity)
 *               setup_background_dofs();
 *           }
 * 
 *         if (interpolated_velocity)
 *           {
 *             interpolate_function_to_field();
 *             euler_step_interpolated(discrete_time.get_previous_step_size());
 *           }
 *         else
 *           euler_step_analytical(discrete_time.get_previous_step_size());
 * 
 * @endcode
 * 
 * After the particles have been moved, it is necessary to identify
 * in which cell they now reside. This is achieved by calling
 * <code>sort_particles_into_subdomains_and_cells</code>
 * 
 * @code
 *         particle_handler.sort_particles_into_subdomains_and_cells();
 * 
 *         if ((discrete_time.get_step_number() % par.output_frequency) == 0)
 *           {
 *             output_particles(discrete_time.get_step_number());
 *             if (interpolated_velocity)
 *               output_background(discrete_time.get_step_number());
 *           }
 *       }
 *   }
 * 
 * } // namespace Step68
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * The remainder of the code, the `main()` function, is standard.
 * We note that we run the particle tracking with the analytical velocity
 * and the interpolated velocity and produce both results
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   using namespace Step68;
 *   using namespace dealii;
 *   deallog.depth_console(1);
 * 
 *   try
 *     {
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       std::string prm_file;
 *       if (argc > 1)
 *         prm_file = argv[1];
 *       else
 *         prm_file = "parameters.prm";
 * 
 *       ParticleTrackingParameters par;
 *       ParameterAcceptor::initialize(prm_file);
 *       {
 *         Step68::ParticleTracking<2> particle_tracking(par, false);
 *         particle_tracking.run();
 *       }
 *       {
 *         Step68::ParticleTracking<2> particle_tracking(par, true);
 *         particle_tracking.run();
 *       }
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


The directory in which this program is run contains an example parameter file by default.
If you do not specify a parameter file as an argument on the command
line, the program will try to read the file "parameters.prm" by default, and
will execute the code.

On any number of cores, the simulation output will look like:

@code
bash$ mpirun -np 4 ./step-68 parameters.prm
Number of particles inserted: 606
Repartitioning triangulation after particle generation
Writing particle output file: analytical-particles-0
Writing particle output file: analytical-particles-10
Writing particle output file: analytical-particles-20
Writing particle output file: analytical-particles-30
...
Number of particles inserted: 606
Repartitioning triangulation after particle generation
Writing particle output file: interpolated-particles-0
Writing background field file: background-0
Writing particle output file: interpolated-particles-10
Writing background field file: background-10
Writing particle output file: interpolated-particles-20
Writing background field file: background-20
Writing particle output file: interpolated-particles-30
Writing background field file: background-30
...
Writing particle output file: interpolated-particles-1980
Writing background field file: background-1980
Writing particle output file: interpolated-particles-1990
Writing background field file: background-1990
Writing particle output file: interpolated-particles-2000
Writing background field file: background-2000
@endcode

We note that, by default, the simulation runs the particle tracking with
an analytical velocity for 2000 iterations, then restarts from the beginning and runs the particle tracking with
velocity interpolation for the same duration. The results are written every
10th iteration.

<a name="Motionoftheparticles"></a><h3> Motion of the particles </h3>


The following animation displays the trajectory of the particles as they
are advected by the flow field. We see that after the complete duration of the
flow, the particle go back to their initial configuration as is expected.

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/EbgS5Ch35Xs"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly

<a name="Dynamicloadbalancing"></a><h3> Dynamic load balancing </h3>


The following animation shows the impact of dynamic load balancing. We clearly
see that the subdomains adapt themselves to balance the number of particles per
subdomain. However, a perfect load balancing is not reached, in part due to
the coarseness of the background mesh.

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/ubUcsR4ECj4"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly


<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


This program highlights some of the main capabilities for handling particles in deal.II, notably their
capacity to be used in distributed parallel simulations. However, this step could
be extended in numerous manners:
- High-order time integration (for example using a Runge-Kutta 4 method) could be
used to increase the accuracy or allow for larger time-step sizes with the same accuracy.
- The full equation of motion (with inertia) could be solved for the particles. In
this case the particles would need to have additional properties such as their mass,
as in step-19, and if one wanted to also consider interactions with the fluid, their diameter.
- Coupling to a flow solver. This step could be straightforwardly coupled to any parallel
program in which the Stokes (step-32, step-70) or the Navier-Stokes equations are solved (e.g., step-57).
- Computing the difference in final particle positions between the two models
would allow to quantify the influence of the interpolation error on particle motion.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-68.cc"
*/
