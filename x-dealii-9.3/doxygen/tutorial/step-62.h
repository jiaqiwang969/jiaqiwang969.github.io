/**
@page step_62 The step-62 tutorial program
This tutorial depends on step-8, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Elasticequations">Elastic equations</a>
        <li><a href="#Simulationparameters">Simulation parameters</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Auxiliaryclassesandfunctions">Auxiliary classes and functions</a>
      <ul>
        <li><a href="#TheRightHandSideclass">The `RightHandSide` class</a>
        <li><a href="#ThePMLclass">The `PML` class</a>
        <li><a href="#TheRhoclass">The `Rho` class</a>
        <li><a href="#TheParametersclass">The `Parameters` class</a>
        <li><a href="#TheQuadratureCacheclass">The `QuadratureCache` class</a>
        <li><a href="#Theget_stiffness_tensorfunction">The `get_stiffness_tensor()` function</a>
      </ul>
        <li><a href="#TheElasticWaveclass">The `ElasticWave` class</a>
        <li><a href="#Implementationoftheauxiliaryclasses">Implementation of the auxiliary classes</a>
      <ul>
        <li><a href="#TheRightHandSideclassimplementation">The `RightHandSide` class implementation</a>
        <li><a href="#ThePMLclassimplementation">The `PML` class implementation</a>
        <li><a href="#TheRhoclassimplementation">The `Rho` class implementation</a>
        <li><a href="#TheParametersclassimplementation">The `Parameters` class implementation</a>
        <li><a href="#TheQuadratureCacheclassimplementation">The `QuadratureCache` class implementation</a>
      </ul>
        <li><a href="#ImplementationoftheElasticWaveclass">Implementation of the `ElasticWave` class</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#ElasticWavesetup_system">ElasticWave::setup_system</a>
        <li><a href="#ElasticWaveassemble_system">ElasticWave::assemble_system</a>
        <li><a href="#ElasticWavesolve">ElasticWave::solve</a>
        <li><a href="#ElasticWaveinitialize_position_vector">ElasticWave::initialize_position_vector</a>
        <li><a href="#ElasticWavestore_frequency_step_data">ElasticWave::store_frequency_step_data</a>
        <li><a href="#ElasticWaveoutput_results">ElasticWave::output_results</a>
        <li><a href="#ElasticWavesetup_quadrature_cache">ElasticWave::setup_quadrature_cache</a>
        <li><a href="#ElasticWavefrequency_sweep">ElasticWave::frequency_sweep</a>
        <li><a href="#ElasticWaverun">ElasticWave::run</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Resonancefrequencyandbandgap">Resonance frequency and bandgap</a>
        <li><a href="#Modeprofile">Mode profile</a>
        <li><a href="#Experimentalapplications">Experimental applications</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Daniel Garcia-Sanchez.</i>
<br>


@note As a prerequisite of this program, you need to have HDF5, complex PETSc,
and the p4est libraries installed. The installation of deal.II
together with these additional libraries is described in the <a
href="../../readme.html" target="body">README</a> file.

<a name="Introduction"></a><h1>Introduction</h1>

A phononic crystal is a periodic nanostructure that modifies the motion of
mechanical vibrations or [phonons](https://en.wikipedia.org/wiki/Phonon).
Phononic structures can be used to disperse, route and confine mechanical vibrations.
These structures have potential applications in
[quantum information](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.86.1391)
and have been used to study
[macroscopic quantum phenomena](https://science.sciencemag.org/content/358/6360/203).
Phononic crystals are usually fabricated in
[cleanrooms](https://en.wikipedia.org/wiki/Cleanroom).

In this tutorial we show how to a design a
[phononic superlattice cavity](https://doi.org/10.1103/PhysRevA.94.033813)
which is a particular type of phononic crystal that can be used to confine
mechanical vibrations. A phononic superlattice cavity is formed by two
[Distributed Bragg Reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector),
mirrors and a $\lambda/2$ cavity where $\lambda$ is the acoustic
wavelength. Acoustic DBRs are  periodic structures where a set of bilayer
stacks with contrasting physical properties (sound velocity index) is
repeated $N$ times.
Superlattice cavities are usually grown on a
[Gallium Arsenide](https://en.wikipedia.org/wiki/Gallium_arsenide)
wafer by
[Molecular Beam Epitaxy](https://en.wikipedia.org/wiki/Molecular-beam_epitaxy).
The bilayers correspond to GaAs/AlAs mirror pairs.
As shown below, the thickness of the mirror layers (brown and green) is
$\lambda/4$ and the thickness of the cavity (blue) is $\lambda/2$.

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.01.svg" height="200" />

In this tutorial we calculate the
[band gap](https://en.wikipedia.org/wiki/Band_gap) and the
mechanical resonance of a phononic superlattice cavity but the code presented here
can be easily used to design and calculate other types of
[phononic crystals](https://science.sciencemag.org/content/358/6360/203).

The device is a waveguide in which the wave goes from left to right.
The simulations of this tutorial are done in 2D, but the code is dimension
independent and can be easily used with 3D simulations.
The waveguide width is equal to the $y$ dimension of the domain and the
waveguide length is equal to the $x$ dimension of the domain.
There are two regimes that depend on the waveguide width:
- Single mode: In this case the width of the structure is much
  smaller than the wavelength.
  This case can be solved either with FEM (the approach that we take here) or with
  a simple semi-analytical
  [1D transfer matrix formalism](https://en.wikipedia.org/wiki/Transfer_matrix).
- Multimode: In this case the width of the structure is larger than the wavelength.
  This case can be solved using FEM
  or with a [scattering matrix formalism](https://doi.org/10.1103/PhysRevA.94.033813).
  Although we do not study this case in this tutorial, it is very easy to reach the multimode
  regime by increasing the parameter waveguide width (`dimension_y` in the jupyter notebook).

The simulations of this tutorial are performed in the frequency domain.
To calculate the transmission spectrum, we use a
[procedure](https://meep.readthedocs.io/en/latest/Python_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity/)
that is commonly used in time domain [FDTD](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method)
simulations. A pulse at a certain frequency is generated on the left side of the
structure and the transmitted energy is measured on the right side of the structure.
The simulation is run twice. First, we run the simulation with the phononic
structure and measure the transmitted energy:

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.02.svg" height="200" />

Then, we run the simulation without the phononic structure and measure the transmitted
energy. We use the simulation without the structure for the calibration:

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.03.svg" height="200" />

The transmission coefficient corresponds to the energy of the first simulation
divided by the calibration energy.
We repeat this procedure for each frequency step.


<a name="Elasticequations"></a><h3>Elastic equations</h3>

What we want to simulate here is the transmission of elastic
waves. Consequently, the right description of the problem uses the
elastic equations, which in the time domain are given by
@f[
\rho\partial_{tt} u_i - \partial_j (c_{ijkl} \varepsilon_{kl}) = f_i,
\qquad i=0,1,2
@f]
where the stiffness tensor $c_{ijkl}$ depends on the spatial coordinates and
the strain is the symmetrized gradient of the displacement, given by
@f[
\varepsilon_{kl} =\frac{1}{2}(\partial_k u_l + \partial_l u_k)
@f]

[A perfectly matched layer (PML)](https://en.wikipedia.org/wiki/Perfectly_matched_layer)
can be used to truncate the solution at the boundaries.
A PML is a transformation that results in a complex coordinate
stretching.

Instead of a time domain approach, this tutorial program converts the
equations above into the frequency domain by performing a Fourier
transform with regard to the time variable.
The elastic equations in the frequency domain then read as follows
@f{eqnarray*}
\nabla\cdot(\boldsymbol{\bar\sigma} \xi \boldsymbol{\Lambda})&=&-\omega^2\rho\xi\mathbf{\bar u}\\
\boldsymbol{\bar \sigma} &=&\mathbf{C}\boldsymbol{\bar\varepsilon}\\
\boldsymbol{\bar\varepsilon}&=&\frac{1}{2}[(\nabla\mathbf{\bar{u}}\boldsymbol{\Lambda}+\boldsymbol{\Lambda}^\mathrm{T}(\nabla\mathbf{\bar{u}})^\mathrm{T})]\\
\xi &=&\prod_i^\textrm{dim}s_i\\
\boldsymbol{\Lambda} &=& \operatorname{diag}(1/s_0,1/s_1,1/s_2)\qquad\textrm{for 3D}\\
\boldsymbol{\Lambda} &=& \operatorname{diag}(1/s_0,1/s_1)\qquad\textrm{for 2D}
@f}
where the coefficients $s_i = 1+is_i'(x,y,z)$ account for the absorption.
There are 3 $s_i$ coefficients in 3D and 2 in 2D.
The imaginary par of $s_i$ is equal to zero outside the PML.
The PMLs are reflectionless only for the exact wave equations.
When the set of equations is discretized the PML is no longer reflectionless.
The reflections can be made arbitrarily small as long as the
medium is slowly varying, see
[the adiabatic theorem](https://doi.org/10.1103/PhysRevE.66.066608).
In the code a quadratic turn-on of the PML has been used.
A linear and cubic turn-on is also
[known to work](https://doi.org/10.1364/OE.16.011376).
These equations can be expanded into
@f[
-\omega^2\rho \xi  u_m - \partial_n \left(\frac{\xi}{s_n}c_{mnkl}
\varepsilon_{kl}\right) = f_m
@f]
@f[
\varepsilon_{kl} =\frac{1}{2}\left(\frac{1}{s_k}\partial_k u_l
+ \frac{1}{s_l}\partial_l u_k\right)
@f]
where summation over repeated indices (here $n$, as well as $k$ and $l$) is as always implied.
Note that the strain is no longer symmetric after applying the complex coordinate
stretching of the PML.
This set of equations can be written as
@f[
-\omega^2\rho \xi  u_m - \partial_n \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k u_l
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l u_k\right) = f_m
@f]

The same as the strain, the stress tensor is not symmetric inside the PML ($s_j\neq 0$).
Indeed the fields inside the PML are not physical.
It is useful to introduce the tensors $\alpha_{mnkl}$ and $\beta_{mnkl}$.
@f[
-\omega^2\rho \xi  u_m - \partial_n \left(\alpha_{mnkl}\partial_k u_l
+  \beta_{mnkl}\partial_l u_k\right) = f_m
@f]

We can multiply by $\varphi_m$ and integrate over the domain $\Omega$ and integrate by parts.
@f{eqnarray*}
-\omega^2\int_\Omega\rho\xi\varphi_m u_m + \int_\Omega\partial_n\varphi_m \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k u_l
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l u_k\right) = \int_\Omega\varphi_m f_m
@f}
It is this set of equations we want to solve for a set of frequencies $\omega$ in order to compute the
transmission coefficient as function of frequency.
The linear system becomes
@f{eqnarray*}
AU&=&F\\
A_{ij} &=& -\omega^2\int_\Omega\rho \xi\varphi_m^i \varphi_m^j + \int_\Omega\partial_n\varphi_m^i \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k \varphi_l^j
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l \varphi_k^j\right)\\
F_i &=& \int_\Omega\varphi_m^i f_m
@f}

<a name="Simulationparameters"></a><h3>Simulation parameters</h3>

In this tutorial we use a python
[jupyter notebook](https://github.com/dealii/dealii/blob/phononic-cavity/examples/step-62/step-62.ipynb)
to set up the parameters and run the simulation.
First we create a HDF5 file where we store the parameters and the results of
the simulation.

Each of the simulations (displacement and calibration) is stored in a separate HDF5 group:
@code{.py}
import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import scipy.constants as constants
import scipy.optimize

# This considerably reduces the size of the svg data
plt.rcParams['svg.fonttype'] = 'none'

h5_file = h5py.File('results.h5', 'w')
data = h5_file.create_group('data')
displacement = data.create_group('displacement')
calibration = data.create_group('calibration')

# Set the parameters
for group in [displacement, calibration]:
    # Dimensions of the domain
    # The waveguide length is equal to dimension_x
    group.attrs['dimension_x'] = 2e-5
    # The waveguide width is equal to dimension_y
    group.attrs['dimension_y'] = 2e-8

    # Position of the probe that we use to measure the flux
    group.attrs['probe_pos_x']   = 8e-6
    group.attrs['probe_pos_y']   = 0
    group.attrs['probe_width_y'] = 2e-08

    # Number of points in the probe
    group.attrs['nb_probe_points'] = 5

    # Global refinement
    group.attrs['grid_level'] = 1

    # Cavity
    group.attrs['cavity_resonance_frequency'] = 20e9
    group.attrs['nb_mirror_pairs']            = 15

    # Material
    group.attrs['poissons_ratio'] = 0.27
    group.attrs['youngs_modulus'] = 270000000000.0
    group.attrs['material_a_rho'] = 3200
    if group == displacement:
        group.attrs['material_b_rho'] = 2000
    else:
        group.attrs['material_b_rho'] = 3200
    group.attrs['lambda'] = (group.attrs['youngs_modulus'] * group.attrs['poissons_ratio'] /
                           ((1 + group.attrs['poissons_ratio']) *
                           (1 - 2 * group.attrs['poissons_ratio'])))
    group.attrs['mu']= (group.attrs['youngs_modulus'] / (2 * (1 + group.attrs['poissons_ratio'])))

    # Force
    group.attrs['max_force_amplitude'] = 1e26
    group.attrs['force_sigma_x']       = 1e-7
    group.attrs['force_sigma_y']       = 1
    group.attrs['max_force_width_x']   = 3e-7
    group.attrs['max_force_width_y']   = 2e-8
    group.attrs['force_x_pos']         = -8e-6
    group.attrs['force_y_pos']         = 0

    # PML
    group.attrs['pml_x']            = True
    group.attrs['pml_y']            = False
    group.attrs['pml_width_x']      = 1.8e-6
    group.attrs['pml_width_y']      = 5e-7
    group.attrs['pml_coeff']        = 1.6
    group.attrs['pml_coeff_degree'] = 2

    # Frequency sweep
    group.attrs['center_frequency']    = 20e9
    group.attrs['frequency_range']     = 0.5e9
    group.attrs['start_frequency']     = group.attrs['center_frequency'] - group.attrs['frequency_range'] / 2
    group.attrs['stop_frequency']      = group.attrs['center_frequency'] + group.attrs['frequency_range'] / 2
    group.attrs['nb_frequency_points'] = 400

    # Other parameters
    if group == displacement:
        group.attrs['simulation_name'] = 'phononic_cavity_displacement'
    else:
        group.attrs['simulation_name'] = 'phononic_cavity_calibration'
    group.attrs['save_vtu_files'] = False

h5_file.close()
@endcode
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
 * Most of the include files we need for this program have already been
 * discussed in previous programs, in particular in step-40.
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/generic_linear_algebra.h>
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The following header provides the Tensor class that we use to represent the
 * material properties.
 * 
 * @code
 * #include <deal.II/base/tensor.h>
 * 
 * 
 * @endcode
 * 
 * The following header is necessary for the HDF5 interface of deal.II.
 * 
 * @code
 * #include <deal.II/base/hdf5.h>
 * 
 * @endcode
 * 
 * This header is required for the function VectorTools::point_value that we use
 * to evaluate the result of the simulation.
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * We need this header for the function GridTools::find_active_cell_around_point
 * that we use in the function `ElasticWave::store_frequency_step_data()`
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h>
 * 
 * namespace step62
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Auxiliaryclassesandfunctions"></a> 
 * <h3>Auxiliary classes and functions</h3>
 * The following classes are used to store the parameters of the simulation.
 * 

 * 
 * 
 * <a name="TheRightHandSideclass"></a> 
 * <h4>The `RightHandSide` class</h4>
 * This class is used to define the force pulse on the left side of the
 * structure:
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     RightHandSide(HDF5::Group &data);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 * 
 *   private:
 * @endcode
 * 
 * The variable `data` is the HDF5::Group in which all the simulation
 * results will be stored. Note that the variables `RightHandSide::data`,
 * `PML::data`, `Rho::data` and `Parameters::data` point to the same group
 * of the HDF5 file. When a HDF5::Group is copied, it will point to the same
 * group of the HDF5 file.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The simulation parameters are stored in `data` as HDF5 attributes. The
 * following attributes are defined in the jupyter notebook, stored in
 * `data` as HDF5 attributes and then read by the constructor.
 * 
 * @code
 *     const double     max_force_amplitude;
 *     const double     force_sigma_x;
 *     const double     force_sigma_y;
 *     const double     max_force_width_x;
 *     const double     max_force_width_y;
 *     const Point<dim> force_center;
 * 
 *   public:
 * @endcode
 * 
 * In this particular simulation the force has only a $x$ component,
 * $F_y=0$.
 * 
 * @code
 *     const unsigned int force_component = 0;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="ThePMLclass"></a> 
 * <h4>The `PML` class</h4>
 * This class is used to define the shape of the Perfectly Matches
 * Layer (PML) to absorb waves traveling towards the boundary:
 * 
 * @code
 *   template <int dim>
 *   class PML : public Function<dim, std::complex<double>>
 *   {
 *   public:
 *     PML(HDF5::Group &data);
 * 
 *     virtual std::complex<double>
 *     value(const Point<dim> &p, const unsigned int component) const override;
 * 
 *   private:
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const double pml_coeff;
 *     const int    pml_coeff_degree;
 *     const double dimension_x;
 *     const double dimension_y;
 *     const bool   pml_x;
 *     const bool   pml_y;
 *     const double pml_width_x;
 *     const double pml_width_y;
 *     const double a_coeff_x;
 *     const double a_coeff_y;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRhoclass"></a> 
 * <h4>The `Rho` class</h4>
 * This class is used to define the mass density.
 * 
 * @code
 *   template <int dim>
 *   class Rho : public Function<dim>
 *   {
 *   public:
 *     Rho(HDF5::Group &data);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const double       lambda;
 *     const double       mu;
 *     const double       material_a_rho;
 *     const double       material_b_rho;
 *     const double       cavity_resonance_frequency;
 *     const unsigned int nb_mirror_pairs;
 *     const double       dimension_y;
 *     const unsigned int grid_level;
 *     double             average_rho_width;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheParametersclass"></a> 
 * <h4>The `Parameters` class</h4>
 * This class contains all the parameters that will be used in the simulation.
 * 
 * @code
 *   template <int dim>
 *   class Parameters
 *   {
 *   public:
 *     Parameters(HDF5::Group &data);
 * 
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const std::string        simulation_name;
 *     const bool               save_vtu_files;
 *     const double             start_frequency;
 *     const double             stop_frequency;
 *     const unsigned int       nb_frequency_points;
 *     const double             lambda;
 *     const double             mu;
 *     const double             dimension_x;
 *     const double             dimension_y;
 *     const unsigned int       nb_probe_points;
 *     const unsigned int       grid_level;
 *     const Point<dim>         probe_start_point;
 *     const Point<dim>         probe_stop_point;
 *     const RightHandSide<dim> right_hand_side;
 *     const PML<dim>           pml;
 *     const Rho<dim>           rho;
 * 
 *   private:
 *     const double comparison_float_constant = 1e-12;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheQuadratureCacheclass"></a> 
 * <h4>The `QuadratureCache` class</h4>
 * The calculation of the mass and stiffness matrices is very expensive. These
 * matrices are the same for all the frequency steps. The right hand side
 * vector is also the same for all the frequency steps. We use this class to
 * store these objects and re-use them at each frequency step. Note that here
 * we don't store the assembled mass and stiffness matrices and right hand
 * sides, but instead the data for a single cell. `QuadratureCache` class is
 * very similar to the `PointHistory` class that has been used in step-18.
 * 
 * @code
 *   template <int dim>
 *   class QuadratureCache
 *   {
 *   public:
 *     QuadratureCache(const unsigned int dofs_per_cell);
 * 
 *   private:
 *     unsigned int dofs_per_cell;
 * 
 *   public:
 * @endcode
 * 
 * We store the mass and stiffness matrices in the variables
 * mass_coefficient and stiffness_coefficient. We store as well the
 * right_hand_side and JxW values which are going to be the same for all the
 * frequency steps.
 * 
 * @code
 *     FullMatrix<std::complex<double>>  mass_coefficient;
 *     FullMatrix<std::complex<double>>  stiffness_coefficient;
 *     std::vector<std::complex<double>> right_hand_side;
 *     double                            JxW;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theget_stiffness_tensorfunction"></a> 
 * <h4>The `get_stiffness_tensor()` function</h4>
 * 

 * 
 * This function returns the stiffness tensor of the material. For the sake of
 * simplicity we consider the stiffness to be isotropic and homogeneous; only
 * the density $\rho$ depends on the position. As we have previously shown in
 * step-8, if the stiffness is isotropic and homogeneous, the stiffness
 * coefficients $c_{ijkl}$ can be expressed as a function of the two
 * coefficients $\lambda$ and $\mu$. The coefficient tensor reduces to
 * @f[
 * c_{ijkl}
 * =
 * \lambda \delta_{ij} \delta_{kl} +
 * \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).
 * @f]
 * 
 * @code
 *   template <int dim>
 *   SymmetricTensor<4, dim> get_stiffness_tensor(const double lambda,
 *                                                const double mu)
 *   {
 *     SymmetricTensor<4, dim> stiffness_tensor;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = 0; j < dim; ++j)
 *         for (unsigned int k = 0; k < dim; ++k)
 *           for (unsigned int l = 0; l < dim; ++l)
 *             stiffness_tensor[i][j][k][l] =
 *               (((i == k) && (j == l) ? mu : 0.0) +
 *                ((i == l) && (j == k) ? mu : 0.0) +
 *                ((i == j) && (k == l) ? lambda : 0.0));
 *     return stiffness_tensor;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheElasticWaveclass"></a> 
 * <h3>The `ElasticWave` class</h3>
 * 

 * 
 * Next let's declare the main class of this program. Its structure is very
 * similar to the step-40 tutorial program. The main differences are:
 * - The sweep over the frequency values.
 * - We save the stiffness and mass matrices in `quadrature_cache` and
 * use them for each frequency step.
 * - We store the measured energy by the probe for each frequency step in the
 * HDF5 file.
 * 
 * @code
 *   template <int dim>
 *   class ElasticWave
 *   {
 *   public:
 *     ElasticWave(const Parameters<dim> &parameters);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system(const double omega,
 *                          const bool   calculate_quadrature_data);
 *     void solve();
 *     void initialize_probe_positions_vector();
 *     void store_frequency_step_data(const unsigned int frequency_idx);
 *     void output_results();
 * 
 * @endcode
 * 
 * This is called before every frequency step to set up a pristine state
 * for the cache variables.
 * 
 * @code
 *     void setup_quadrature_cache();
 * 
 * @endcode
 * 
 * This function loops over the frequency vector and runs the simulation for
 * each frequency step.
 * 
 * @code
 *     void frequency_sweep();
 * 
 * @endcode
 * 
 * The parameters are stored in this variable.
 * 
 * @code
 *     Parameters<dim> parameters;
 * 
 *     MPI_Comm mpi_communicator;
 * 
 *     parallel::distributed::Triangulation<dim> triangulation;
 * 
 *     QGauss<dim> quadrature_formula;
 * 
 * @endcode
 * 
 * We store the mass and stiffness matrices for each cell this vector.
 * 
 * @code
 *     std::vector<QuadratureCache<dim>> quadrature_cache;
 * 
 * 
 *     FESystem<dim>   fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 *     AffineConstraints<std::complex<double>> constraints;
 * 
 *     LinearAlgebraPETSc::MPI::SparseMatrix system_matrix;
 *     LinearAlgebraPETSc::MPI::Vector       locally_relevant_solution;
 *     LinearAlgebraPETSc::MPI::Vector       system_rhs;
 * 
 * 
 * @endcode
 * 
 * This vector contains the range of frequencies that we are going to
 * simulate.
 * 
 * @code
 *     std::vector<double> frequency;
 * 
 * @endcode
 * 
 * This vector contains the coordinates $(x,y)$ of the points of the
 * measurement probe.
 * 
 * @code
 *     FullMatrix<double> probe_positions;
 * 
 * @endcode
 * 
 * HDF5 datasets to store the frequency and `probe_positions` vectors.
 * 
 * @code
 *     HDF5::DataSet frequency_dataset;
 *     HDF5::DataSet probe_positions_dataset;
 * 
 * @endcode
 * 
 * HDF5 dataset that stores the values of the energy measured by the probe.
 * 
 * @code
 *     HDF5::DataSet displacement;
 * 
 * 
 *     ConditionalOStream pcout;
 *     TimerOutput        computing_timer;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Implementationoftheauxiliaryclasses"></a> 
 * <h3>Implementation of the auxiliary classes</h3>
 * 

 * 
 * 
 * <a name="TheRightHandSideclassimplementation"></a> 
 * <h4>The `RightHandSide` class implementation</h4>
 * 

 * 
 * The constructor reads all the parameters from the HDF5::Group `data` using
 * the HDF5::Group::get_attribute() function.
 * 
 * @code
 *   template <int dim>
 *   RightHandSide<dim>::RightHandSide(HDF5::Group &data)
 *     : Function<dim>(dim)
 *     , data(data)
 *     , max_force_amplitude(data.get_attribute<double>("max_force_amplitude"))
 *     , force_sigma_x(data.get_attribute<double>("force_sigma_x"))
 *     , force_sigma_y(data.get_attribute<double>("force_sigma_y"))
 *     , max_force_width_x(data.get_attribute<double>("max_force_width_x"))
 *     , max_force_width_y(data.get_attribute<double>("max_force_width_y"))
 *     , force_center(Point<dim>(data.get_attribute<double>("force_x_pos"),
 *                               data.get_attribute<double>("force_y_pos")))
 *   {}
 * 
 * @endcode
 * 
 * This function defines the spatial shape of the force vector pulse which
 * takes the form of a Gaussian function
 * @f{align*}
 * F_x &=
 * \left\{
 * \begin{array}{ll}
 * a \exp(- (\frac{(x-b_x)^2 }{ 2 \sigma_x^2}+\frac{(y-b_y)^2 }{ 2
 * \sigma_y^2}))
 * & \text{if}\, x_\textrm{min} <x<x_\textrm{max}\, \text{and}\,
 * y_\textrm{min} <y<y_\textrm{max}  \\ 0 & \text{otherwise},
 * \end{array}
 * \right.\\ F_y &= 0
 * @f}
 * where $a$ is the maximum amplitude that takes the force and $\sigma_x$ and
 * $\sigma_y$ are the standard deviations for the $x$ and $y$ components. Note
 * that the pulse has been cropped to $x_\textrm{min}<x<x_\textrm{max}$ and
 * $y_\textrm{min} <y<y_\textrm{max}$.
 * 
 * @code
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 *     if (component == force_component)
 *       {
 *         if (std::abs(p[0] - force_center[0]) < max_force_width_x / 2 &&
 *             std::abs(p[1] - force_center[1]) < max_force_width_y / 2)
 *           {
 *             return max_force_amplitude *
 *                    std::exp(-(std::pow(p[0] - force_center[0], 2) /
 *                                 (2 * std::pow(force_sigma_x, 2)) +
 *                               std::pow(p[1] - force_center[1], 2) /
 *                                 (2 * std::pow(force_sigma_y, 2))));
 *           }
 *         else
 *           {
 *             return 0;
 *           }
 *       }
 *     else
 *       {
 *         return 0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePMLclassimplementation"></a> 
 * <h4>The `PML` class implementation</h4>
 * 

 * 
 * As before, the constructor reads all the parameters from the HDF5::Group
 * `data` using the HDF5::Group::get_attribute() function. As we have
 * discussed, a quadratic turn-on of the PML has been defined in the jupyter
 * notebook. It is possible to use a linear, cubic or another power degree by
 * changing the parameter `pml_coeff_degree`. The parameters `pml_x` and
 * `pml_y` can be used to turn on and off the `x` and `y` PMLs.
 * 
 * @code
 *   template <int dim>
 *   PML<dim>::PML(HDF5::Group &data)
 *     : Function<dim, std::complex<double>>(dim)
 *     , data(data)
 *     , pml_coeff(data.get_attribute<double>("pml_coeff"))
 *     , pml_coeff_degree(data.get_attribute<int>("pml_coeff_degree"))
 *     , dimension_x(data.get_attribute<double>("dimension_x"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , pml_x(data.get_attribute<bool>("pml_x"))
 *     , pml_y(data.get_attribute<bool>("pml_y"))
 *     , pml_width_x(data.get_attribute<double>("pml_width_x"))
 *     , pml_width_y(data.get_attribute<double>("pml_width_y"))
 *     , a_coeff_x(pml_coeff / std::pow(pml_width_x, pml_coeff_degree))
 *     , a_coeff_y(pml_coeff / std::pow(pml_width_y, pml_coeff_degree))
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * The PML coefficient for the `x` component takes the form
 * $s'_x = a_x x^{\textrm{degree}}$
 * 
 * @code
 *   template <int dim>
 *   std::complex<double> PML<dim>::value(const Point<dim> & p,
 *                                        const unsigned int component) const
 *   {
 *     double calculated_pml_x_coeff = 0;
 *     double calculated_pml_y_coeff = 0;
 * 
 *     if ((component == 0) && pml_x)
 *       {
 *         const double pml_x_start_position = dimension_x / 2 - pml_width_x;
 *         if (std::abs(p[0]) > pml_x_start_position)
 *           {
 *             const double x_prime = std::abs(p[0]) - pml_x_start_position;
 *             calculated_pml_x_coeff =
 *               a_coeff_x * std::pow(x_prime, pml_coeff_degree);
 *           }
 *       }
 * 
 *     if ((component == 1) && pml_y)
 *       {
 *         const double pml_y_start_position = dimension_y / 2 - pml_width_y;
 *         if (std::abs(p[1]) > pml_y_start_position)
 *           {
 *             const double y_prime = std::abs(p[1]) - pml_y_start_position;
 *             calculated_pml_y_coeff =
 *               a_coeff_y * std::pow(y_prime, pml_coeff_degree);
 *           }
 *       }
 * 
 *     return 1. + std::max(calculated_pml_x_coeff, calculated_pml_y_coeff) *
 *                   std::complex<double>(0., 1.);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRhoclassimplementation"></a> 
 * <h4>The `Rho` class implementation</h4>
 * 

 * 
 * This class is used to define the mass density. As we have explaine before,
 * a phononic superlattice cavity is formed by two
 * [Distributed Reflector](https://en.wikipedia.org/wiki/Band_gap),
 * mirrors and a $\lambda/2$ cavity where $\lambda$ is the acoustic
 * wavelength. Acoustic DBRs are periodic structures where a set of bilayer
 * stacks with contrasting physical properties (sound velocity index) is
 * repeated $N$ times. The change of in the wave velocity is generated by
 * alternating layers with different density.
 * 
 * @code
 *   template <int dim>
 *   Rho<dim>::Rho(HDF5::Group &data)
 *     : Function<dim>(1)
 *     , data(data)
 *     , lambda(data.get_attribute<double>("lambda"))
 *     , mu(data.get_attribute<double>("mu"))
 *     , material_a_rho(data.get_attribute<double>("material_a_rho"))
 *     , material_b_rho(data.get_attribute<double>("material_b_rho"))
 *     , cavity_resonance_frequency(
 *         data.get_attribute<double>("cavity_resonance_frequency"))
 *     , nb_mirror_pairs(data.get_attribute<int>("nb_mirror_pairs"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , grid_level(data.get_attribute<int>("grid_level"))
 *   {
 * @endcode
 * 
 * In order to increase the precision we use
 * [subpixel
 * smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/).
 * 
 * @code
 *     average_rho_width = dimension_y / (std::pow(2.0, grid_level));
 *     data.set_attribute("average_rho_width", average_rho_width);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   double Rho<dim>::value(const Point<dim> &p,
 *                          const unsigned int /*component*/) const
 *   {
 * @endcode
 * 
 * The speed of sound is defined by
 * @f[
 * c = \frac{K_e}{\rho}
 * @f]
 * where $K_e$ is the effective elastic constant and $\rho$ the density.
 * Here we consider the case in which the waveguide width is much smaller
 * than the wavelength. In this case it can be shown that for the two
 * dimensional case
 * @f[
 * K_e = 4\mu\frac{\lambda +\mu}{\lambda+2\mu}
 * @f]
 * and for the three dimensional case $K_e$ is equal to the Young's modulus.
 * @f[
 * K_e = \mu\frac{3\lambda +2\mu}{\lambda+\mu}
 * @f]
 * 
 * @code
 *     double elastic_constant;
 *     if (dim == 2)
 *       {
 *         elastic_constant = 4 * mu * (lambda + mu) / (lambda + 2 * mu);
 *       }
 *     else if (dim == 3)
 *       {
 *         elastic_constant = mu * (3 * lambda + 2 * mu) / (lambda + mu);
 *       }
 *     else
 *       {
 *         Assert(false, ExcInternalError());
 *       }
 *     const double material_a_speed_of_sound =
 *       std::sqrt(elastic_constant / material_a_rho);
 *     const double material_a_wavelength =
 *       material_a_speed_of_sound / cavity_resonance_frequency;
 *     const double material_b_speed_of_sound =
 *       std::sqrt(elastic_constant / material_b_rho);
 *     const double material_b_wavelength =
 *       material_b_speed_of_sound / cavity_resonance_frequency;
 * 
 * @endcode
 * 
 * The density $\rho$ takes the following form
 * <img alt="Phononic superlattice cavity"
 * src="https://www.dealii.org/images/steps/developer/step-62.04.svg"
 * height="200" />
 * where the brown color represents material_a and the green color
 * represents material_b.
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_transition_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4);
 *         if (std::abs(p[0]) >=
 *               (layer_transition_center - average_rho_width / 2) &&
 *             std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2))
 *           {
 *             const double coefficient =
 *               (std::abs(p[0]) -
 *                (layer_transition_center - average_rho_width / 2)) /
 *               average_rho_width;
 *             return (1 - coefficient) * material_a_rho +
 *                    coefficient * material_b_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * Here we define the
 * [subpixel
 * smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/)
 * which improves the precision of the simulation.
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_transition_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 4;
 *         if (std::abs(p[0]) >=
 *               (layer_transition_center - average_rho_width / 2) &&
 *             std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2))
 *           {
 *             const double coefficient =
 *               (std::abs(p[0]) -
 *                (layer_transition_center - average_rho_width / 2)) /
 *               average_rho_width;
 *             return (1 - coefficient) * material_b_rho +
 *                    coefficient * material_a_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * then the cavity
 * 
 * @code
 *     if (std::abs(p[0]) <= material_a_wavelength / 2)
 *       {
 *         return material_a_rho;
 *       }
 * 
 * @endcode
 * 
 * the material_a layers
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 4 + material_a_wavelength / 8;
 *         const double layer_width = material_a_wavelength / 4;
 *         if (std::abs(p[0]) >= (layer_center - layer_width / 2) &&
 *             std::abs(p[0]) <= (layer_center + layer_width / 2))
 *           {
 *             return material_a_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * the material_b layers
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 8;
 *         const double layer_width = material_b_wavelength / 4;
 *         if (std::abs(p[0]) >= (layer_center - layer_width / 2) &&
 *             std::abs(p[0]) <= (layer_center + layer_width / 2))
 *           {
 *             return material_b_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * and finally the default is material_a.
 * 
 * @code
 *     return material_a_rho;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheParametersclassimplementation"></a> 
 * <h4>The `Parameters` class implementation</h4>
 * 

 * 
 * The constructor reads all the parameters from the HDF5::Group `data` using
 * the HDF5::Group::get_attribute() function.
 * 
 * @code
 *   template <int dim>
 *   Parameters<dim>::Parameters(HDF5::Group &data)
 *     : data(data)
 *     , simulation_name(data.get_attribute<std::string>("simulation_name"))
 *     , save_vtu_files(data.get_attribute<bool>("save_vtu_files"))
 *     , start_frequency(data.get_attribute<double>("start_frequency"))
 *     , stop_frequency(data.get_attribute<double>("stop_frequency"))
 *     , nb_frequency_points(data.get_attribute<int>("nb_frequency_points"))
 *     , lambda(data.get_attribute<double>("lambda"))
 *     , mu(data.get_attribute<double>("mu"))
 *     , dimension_x(data.get_attribute<double>("dimension_x"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , nb_probe_points(data.get_attribute<int>("nb_probe_points"))
 *     , grid_level(data.get_attribute<int>("grid_level"))
 *     , probe_start_point(data.get_attribute<double>("probe_pos_x"),
 *                         data.get_attribute<double>("probe_pos_y") -
 *                           data.get_attribute<double>("probe_width_y") / 2)
 *     , probe_stop_point(data.get_attribute<double>("probe_pos_x"),
 *                        data.get_attribute<double>("probe_pos_y") +
 *                          data.get_attribute<double>("probe_width_y") / 2)
 *     , right_hand_side(data)
 *     , pml(data)
 *     , rho(data)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheQuadratureCacheclassimplementation"></a> 
 * <h4>The `QuadratureCache` class implementation</h4>
 * 

 * 
 * We need to reserve enough space for the mass and stiffness matrices and the
 * right hand side vector.
 * 
 * @code
 *   template <int dim>
 *   QuadratureCache<dim>::QuadratureCache(const unsigned int dofs_per_cell)
 *     : dofs_per_cell(dofs_per_cell)
 *     , mass_coefficient(dofs_per_cell, dofs_per_cell)
 *     , stiffness_coefficient(dofs_per_cell, dofs_per_cell)
 *     , right_hand_side(dofs_per_cell)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationoftheElasticWaveclass"></a> 
 * <h3>Implementation of the `ElasticWave` class</h3>
 * 

 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * This is very similar to the constructor of step-40. In addition we create
 * the HDF5 datasets `frequency_dataset`, `position_dataset` and
 * `displacement`. Note the use of the `template` keyword for the creation of
 * the HDF5 datasets. It is a C++ requirement to use the `template` keyword in
 * order to treat `create_dataset` as a dependent template name.
 * 
 * @code
 *   template <int dim>
 *   ElasticWave<dim>::ElasticWave(const Parameters<dim> &parameters)
 *     : parameters(parameters)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , triangulation(mpi_communicator,
 *                     typename Triangulation<dim>::MeshSmoothing(
 *                       Triangulation<dim>::smoothing_on_refinement |
 *                       Triangulation<dim>::smoothing_on_coarsening))
 *     , quadrature_formula(2)
 *     , fe(FE_Q<dim>(1), dim)
 *     , dof_handler(triangulation)
 *     , frequency(parameters.nb_frequency_points)
 *     , probe_positions(parameters.nb_probe_points, dim)
 *     , frequency_dataset(parameters.data.template create_dataset<double>(
 *         "frequency",
 *         std::vector<hsize_t>{parameters.nb_frequency_points}))
 *     , probe_positions_dataset(parameters.data.template create_dataset<double>(
 *         "position",
 *         std::vector<hsize_t>{parameters.nb_probe_points, dim}))
 *     , displacement(
 *         parameters.data.template create_dataset<std::complex<double>>(
 *           "displacement",
 *           std::vector<hsize_t>{parameters.nb_probe_points,
 *                                parameters.nb_frequency_points}))
 *     , pcout(std::cout,
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *     , computing_timer(mpi_communicator,
 *                       pcout,
 *                       TimerOutput::summary,
 *                       TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesetup_system"></a> 
 * <h4>ElasticWave::setup_system</h4>
 * 

 * 
 * There is nothing new in this function, the only difference with step-40 is
 * that we don't have to apply boundary conditions because we use the PMLs to
 * truncate the domain.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::setup_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "setup");
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 *     locally_relevant_solution.reinit(locally_owned_dofs,
 *                                      locally_relevant_dofs,
 *                                      mpi_communicator);
 * 
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 * 
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(locally_relevant_dofs);
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 *     SparsityTools::distribute_sparsity_pattern(dsp,
 *                                                locally_owned_dofs,
 *                                                mpi_communicator,
 *                                                locally_relevant_dofs);
 * 
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          dsp,
 *                          mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveassemble_system"></a> 
 * <h4>ElasticWave::assemble_system</h4>
 * 

 * 
 * This function is also very similar to step-40, though there are notable
 * differences. We assemble the system for each frequency/omega step. In the
 * first step we set `calculate_quadrature_data = True` and we calculate the
 * mass and stiffness matrices and the right hand side vector. In the
 * subsequent steps we will use that data to accelerate the calculation.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::assemble_system(const double omega,
 *                                          const bool   calculate_quadrature_data)
 *   {
 *     TimerOutput::Scope t(computing_timer, "assembly");
 * 
 *     FEValues<dim>      fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<std::complex<double>>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Here we store the value of the right hand side, rho and the PML.
 * 
 * @code
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));
 *     std::vector<double>         rho_values(n_q_points);
 *     std::vector<Vector<std::complex<double>>> pml_values(
 *       n_q_points, Vector<std::complex<double>>(dim));
 * 
 * @endcode
 * 
 * We calculate the stiffness tensor for the $\lambda$ and $\mu$ that have
 * been defined in the jupyter notebook. Note that contrary to $\rho$ the
 * stiffness is constant among for the whole domain.
 * 
 * @code
 *     const SymmetricTensor<4, dim> stiffness_tensor =
 *       get_stiffness_tensor<dim>(parameters.lambda, parameters.mu);
 * 
 * @endcode
 * 
 * We use the same method of step-20 for vector-valued problems.
 * 
 * @code
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 * @endcode
 * 
 * We have to calculate the values of the right hand side, rho and
 * the PML only if we are going to calculate the mass and the
 * stiffness matrices. Otherwise we can skip this calculation which
 * considerably reduces the total calculation time.
 * 
 * @code
 *           if (calculate_quadrature_data)
 *             {
 *               fe_values.reinit(cell);
 * 
 *               parameters.right_hand_side.vector_value_list(
 *                 fe_values.get_quadrature_points(), rhs_values);
 *               parameters.rho.value_list(fe_values.get_quadrature_points(),
 *                                         rho_values);
 *               parameters.pml.vector_value_list(
 *                 fe_values.get_quadrature_points(), pml_values);
 *             }
 * 
 * @endcode
 * 
 * We have done this in step-18. Get a pointer to the quadrature
 * cache data local to the present cell, and, as a defensive
 * measure, make sure that this pointer is within the bounds of the
 * global array:
 * 
 * @code
 *           QuadratureCache<dim> *local_quadrature_points_data =
 *             reinterpret_cast<QuadratureCache<dim> *>(cell->user_pointer());
 *           Assert(local_quadrature_points_data >= &quadrature_cache.front(),
 *                  ExcInternalError());
 *           Assert(local_quadrature_points_data <= &quadrature_cache.back(),
 *                  ExcInternalError());
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             {
 * @endcode
 * 
 * The quadrature_data variable is used to store the mass and
 * stiffness matrices, the right hand side vector and the value
 * of `JxW`.
 * 
 * @code
 *               QuadratureCache<dim> &quadrature_data =
 *                 local_quadrature_points_data[q];
 * 
 * @endcode
 * 
 * Below we declare the force vector and the parameters of the
 * PML $s$ and $\xi$.
 * 
 * @code
 *               Tensor<1, dim>                       force;
 *               Tensor<1, dim, std::complex<double>> s;
 *               std::complex<double>                 xi(1, 0);
 * 
 * @endcode
 * 
 * The following block is calculated only in the first frequency
 * step.
 * 
 * @code
 *               if (calculate_quadrature_data)
 *                 {
 * @endcode
 * 
 * Store the value of `JxW`.
 * 
 * @code
 *                   quadrature_data.JxW = fe_values.JxW(q);
 * 
 *                   for (unsigned int component = 0; component < dim; ++component)
 *                     {
 * @endcode
 * 
 * Convert vectors to tensors and calculate xi
 * 
 * @code
 *                       force[component] = rhs_values[q][component];
 *                       s[component]     = pml_values[q][component];
 *                       xi *= s[component];
 *                     }
 * 
 * @endcode
 * 
 * Here we calculate the $\alpha_{mnkl}$ and $\beta_{mnkl}$
 * tensors.
 * 
 * @code
 *                   Tensor<4, dim, std::complex<double>> alpha;
 *                   Tensor<4, dim, std::complex<double>> beta;
 *                   for (unsigned int m = 0; m < dim; ++m)
 *                     for (unsigned int n = 0; n < dim; ++n)
 *                       for (unsigned int k = 0; k < dim; ++k)
 *                         for (unsigned int l = 0; l < dim; ++l)
 *                           {
 *                             alpha[m][n][k][l] = xi *
 *                                                 stiffness_tensor[m][n][k][l] /
 *                                                 (2.0 * s[n] * s[k]);
 *                             beta[m][n][k][l] = xi *
 *                                                stiffness_tensor[m][n][k][l] /
 *                                                (2.0 * s[n] * s[l]);
 *                           }
 * 
 *                   for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                     {
 *                       const Tensor<1, dim> phi_i =
 *                         fe_values[displacement].value(i, q);
 *                       const Tensor<2, dim> grad_phi_i =
 *                         fe_values[displacement].gradient(i, q);
 * 
 *                       for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                         {
 *                           const Tensor<1, dim> phi_j =
 *                             fe_values[displacement].value(j, q);
 *                           const Tensor<2, dim> grad_phi_j =
 *                             fe_values[displacement].gradient(j, q);
 * 
 * @endcode
 * 
 * calculate the values of the mass matrix.
 * 
 * @code
 *                           quadrature_data.mass_coefficient[i][j] =
 *                             rho_values[q] * xi * phi_i * phi_j;
 * 
 * @endcode
 * 
 * Loop over the $mnkl$ indices of the stiffness
 * tensor.
 * 
 * @code
 *                           std::complex<double> stiffness_coefficient = 0;
 *                           for (unsigned int m = 0; m < dim; ++m)
 *                             for (unsigned int n = 0; n < dim; ++n)
 *                               for (unsigned int k = 0; k < dim; ++k)
 *                                 for (unsigned int l = 0; l < dim; ++l)
 *                                   {
 * @endcode
 * 
 * Here we calculate the stiffness matrix.
 * Note that the stiffness matrix is not
 * symmetric because of the PMLs. We use the
 * gradient function (see the
 * [documentation](https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html))
 * which is a <code>Tensor@<2,dim@></code>.
 * The matrix $G_{ij}$ consists of entries
 * @f[
 * G_{ij}=
 * \frac{\partial\phi_i}{\partial x_j}
 * =\partial_j \phi_i
 * @f]
 * Note the position of the indices $i$ and
 * $j$ and the notation that we use in this
 * tutorial: $\partial_j\phi_i$. As the
 * stiffness tensor is not symmetric, it is
 * very easy to make a mistake.
 * 
 * @code
 *                                     stiffness_coefficient +=
 *                                       grad_phi_i[m][n] *
 *                                       (alpha[m][n][k][l] * grad_phi_j[l][k] +
 *                                        beta[m][n][k][l] * grad_phi_j[k][l]);
 *                                   }
 * 
 * @endcode
 * 
 * We save the value of the stiffness matrix in
 * quadrature_data
 * 
 * @code
 *                           quadrature_data.stiffness_coefficient[i][j] =
 *                             stiffness_coefficient;
 *                         }
 * 
 * @endcode
 * 
 * and the value of the right hand side in
 * quadrature_data.
 * 
 * @code
 *                       quadrature_data.right_hand_side[i] =
 *                         phi_i * force * fe_values.JxW(q);
 *                     }
 *                 }
 * 
 * @endcode
 * 
 * We loop again over the degrees of freedom of the cells to
 * calculate the system matrix. These loops are really quick
 * because we have already calculated the stiffness and mass
 * matrices, only the value of $\omega$ changes.
 * 
 * @code
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     {
 *                       std::complex<double> matrix_sum = 0;
 *                       matrix_sum += -std::pow(omega, 2) *
 *                                     quadrature_data.mass_coefficient[i][j];
 *                       matrix_sum += quadrature_data.stiffness_coefficient[i][j];
 *                       cell_matrix(i, j) += matrix_sum * quadrature_data.JxW;
 *                     }
 *                   cell_rhs(i) += quadrature_data.right_hand_side[i];
 *                 }
 *             }
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints.distribute_local_to_global(cell_matrix,
 *                                                  cell_rhs,
 *                                                  local_dof_indices,
 *                                                  system_matrix,
 *                                                  system_rhs);
 *         }
 * 
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesolve"></a> 
 * <h4>ElasticWave::solve</h4>
 * 

 * 
 * This is even more simple than in step-40. We use the parallel direct solver
 * MUMPS which requires less options than an iterative solver. The drawback is
 * that it does not scale very well. It is not straightforward to solve the
 * Helmholtz equation with an iterative solver. The shifted Laplacian
 * multigrid method is a well known approach to precondition this system, but
 * this is beyond the scope of this tutorial.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::solve()
 *   {
 *     TimerOutput::Scope              t(computing_timer, "solve");
 *     LinearAlgebraPETSc::MPI::Vector completely_distributed_solution(
 *       locally_owned_dofs, mpi_communicator);
 * 
 *     SolverControl                    solver_control;
 *     PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
 *     solver.solve(system_matrix, completely_distributed_solution, system_rhs);
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *           << std::endl;
 *     constraints.distribute(completely_distributed_solution);
 *     locally_relevant_solution = completely_distributed_solution;
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveinitialize_position_vector"></a> 
 * <h4>ElasticWave::initialize_position_vector</h4>
 * 

 * 
 * We use this function to calculate the values of the position vector.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::initialize_probe_positions_vector()
 *   {
 *     for (unsigned int position_idx = 0;
 *          position_idx < parameters.nb_probe_points;
 *          ++position_idx)
 *       {
 * @endcode
 * 
 * Because of the way the operator + and - are overloaded to subtract
 * two points, the following has to be done:
 * `Point_b<dim> + (-Point_a<dim>)`
 * 
 * @code
 *         const Point<dim> p =
 *           (position_idx / ((double)(parameters.nb_probe_points - 1))) *
 *             (parameters.probe_stop_point + (-parameters.probe_start_point)) +
 *           parameters.probe_start_point;
 *         probe_positions[position_idx][0] = p[0];
 *         probe_positions[position_idx][1] = p[1];
 *         if (dim == 3)
 *           {
 *             probe_positions[position_idx][2] = p[2];
 *           }
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavestore_frequency_step_data"></a> 
 * <h4>ElasticWave::store_frequency_step_data</h4>
 * 

 * 
 * This function stores in the HDF5 file the measured energy by the probe.
 * 
 * @code
 *   template <int dim>
 *   void
 *   ElasticWave<dim>::store_frequency_step_data(const unsigned int frequency_idx)
 *   {
 *     TimerOutput::Scope t(computing_timer, "store_frequency_step_data");
 * 
 * @endcode
 * 
 * We store the displacement in the $x$ direction; the displacement in the
 * $y$ direction is negligible.
 * 
 * @code
 *     const unsigned int probe_displacement_component = 0;
 * 
 * @endcode
 * 
 * The vector coordinates contains the coordinates in the HDF5 file of the
 * points of the probe that are located in locally owned cells. The vector
 * displacement_data contains the value of the displacement at these points.
 * 
 * @code
 *     std::vector<hsize_t>              coordinates;
 *     std::vector<std::complex<double>> displacement_data;
 *     for (unsigned int position_idx = 0;
 *          position_idx < parameters.nb_probe_points;
 *          ++position_idx)
 *       {
 *         Point<dim> point;
 *         for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *           {
 *             point[dim_idx] = probe_positions[position_idx][dim_idx];
 *           }
 *         bool point_in_locally_owned_cell;
 *         {
 * @endcode
 * 
 * First we have to find out if the point is in a locally owned cell.
 * 
 * @code
 *           auto mapping = StaticMappingQ1<dim>::mapping;
 *           const std::pair<typename DoFHandler<dim>::active_cell_iterator,
 *                           Point<dim>>
 *             cell_point = GridTools::find_active_cell_around_point(mapping,
 *                                                                   dof_handler,
 *                                                                   point);
 * 
 *           point_in_locally_owned_cell = cell_point.first->is_locally_owned();
 *         }
 *         if (point_in_locally_owned_cell)
 *           {
 * @endcode
 * 
 * Then we can store the values of the displacement in the points of
 * the probe in `displacement_data`.
 * 
 * @code
 *             Vector<std::complex<double>> tmp_vector(dim);
 *             VectorTools::point_value(dof_handler,
 *                                      locally_relevant_solution,
 *                                      point,
 *                                      tmp_vector);
 *             coordinates.emplace_back(position_idx);
 *             coordinates.emplace_back(frequency_idx);
 *             displacement_data.emplace_back(
 *               tmp_vector(probe_displacement_component));
 *           }
 *       }
 * 
 * @endcode
 * 
 * We write the displacement data in the HDF5 file. The call
 * HDF5::DataSet::write_selection() is MPI collective which means that all
 * the processes have to participate.
 * 
 * @code
 *     if (coordinates.size() > 0)
 *       {
 *         displacement.write_selection(displacement_data, coordinates);
 *       }
 * @endcode
 * 
 * Therefore even if the process has no data to write it has to participate
 * in the collective call. For this we can use HDF5::DataSet::write_none().
 * Note that we have to specify the data type, in this case
 * `std::complex<double>`.
 * 
 * @code
 *     else
 *       {
 *         displacement.write_none<std::complex<double>>();
 *       }
 * 
 * @endcode
 * 
 * If the variable `save_vtu_files` in the input file equals `True` then all
 * the data will be saved as vtu. The procedure to write `vtu` files has
 * been described in step-40.
 * 
 * @code
 *     if (parameters.save_vtu_files)
 *       {
 *         std::vector<std::string> solution_names(dim, "displacement");
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *           interpretation(
 *             dim, DataComponentInterpretation::component_is_part_of_vector);
 * 
 *         DataOut<dim> data_out;
 *         data_out.add_data_vector(dof_handler,
 *                                  locally_relevant_solution,
 *                                  solution_names,
 *                                  interpretation);
 *         Vector<float> subdomain(triangulation.n_active_cells());
 *         for (unsigned int i = 0; i < subdomain.size(); ++i)
 *           subdomain(i) = triangulation.locally_owned_subdomain();
 *         data_out.add_data_vector(subdomain, "subdomain");
 * 
 *         std::vector<Vector<double>> force(
 *           dim, Vector<double>(triangulation.n_active_cells()));
 *         std::vector<Vector<double>> pml(
 *           dim, Vector<double>(triangulation.n_active_cells()));
 *         Vector<double> rho(triangulation.n_active_cells());
 * 
 *         for (auto &cell : triangulation.active_cell_iterators())
 *           {
 *             if (cell->is_locally_owned())
 *               {
 *                 for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *                   {
 *                     force[dim_idx](cell->active_cell_index()) =
 *                       parameters.right_hand_side.value(cell->center(), dim_idx);
 *                     pml[dim_idx](cell->active_cell_index()) =
 *                       parameters.pml.value(cell->center(), dim_idx).imag();
 *                   }
 *                 rho(cell->active_cell_index()) =
 *                   parameters.rho.value(cell->center());
 *               }
 * @endcode
 * 
 * And on the cells that we are not interested in, set the
 * respective value to a bogus value in order to make sure that if
 * we were somehow wrong about our assumption we would find out by
 * looking at the graphical output:
 * 
 * @code
 *             else
 *               {
 *                 for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *                   {
 *                     force[dim_idx](cell->active_cell_index()) = -1e+20;
 *                     pml[dim_idx](cell->active_cell_index())   = -1e+20;
 *                   }
 *                 rho(cell->active_cell_index()) = -1e+20;
 *               }
 *           }
 * 
 *         for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *           {
 *             data_out.add_data_vector(force[dim_idx],
 *                                      "force_" + std::to_string(dim_idx));
 *             data_out.add_data_vector(pml[dim_idx],
 *                                      "pml_" + std::to_string(dim_idx));
 *           }
 *         data_out.add_data_vector(rho, "rho");
 * 
 *         data_out.build_patches();
 * 
 *         std::stringstream  frequency_idx_stream;
 *         const unsigned int nb_number_positions =
 *           ((unsigned int)std::log10(parameters.nb_frequency_points)) + 1;
 *         frequency_idx_stream << std::setw(nb_number_positions)
 *                              << std::setfill('0') << frequency_idx;
 *         std::string filename = (parameters.simulation_name + "_" +
 *                                 frequency_idx_stream.str() + ".vtu");
 *         data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveoutput_results"></a> 
 * <h4>ElasticWave::output_results</h4>
 * 

 * 
 * This function writes the datasets that have not already been written.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::output_results()
 *   {
 * @endcode
 * 
 * The vectors `frequency` and `position` are the same for all the
 * processes. Therefore any of the processes can write the corresponding
 * `datasets`. Because the call HDF5::DataSet::write is MPI collective, the
 * rest of the processes will have to call HDF5::DataSet::write_none.
 * 
 * @code
 *     if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *       {
 *         frequency_dataset.write(frequency);
 *         probe_positions_dataset.write(probe_positions);
 *       }
 *     else
 *       {
 *         frequency_dataset.write_none<double>();
 *         probe_positions_dataset.write_none<double>();
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesetup_quadrature_cache"></a> 
 * <h4>ElasticWave::setup_quadrature_cache</h4>
 * 

 * 
 * We use this function at the beginning of our computations to set up initial
 * values of the cache variables. This function has been described in step-18.
 * There are no differences with the function of step-18.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::setup_quadrature_cache()
 *   {
 *     triangulation.clear_user_data();
 * 
 *     {
 *       std::vector<QuadratureCache<dim>> tmp;
 *       quadrature_cache.swap(tmp);
 *     }
 * 
 *     quadrature_cache.resize(triangulation.n_locally_owned_active_cells() *
 *                               quadrature_formula.size(),
 *                             QuadratureCache<dim>(fe.n_dofs_per_cell()));
 *     unsigned int cache_index = 0;
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell->set_user_pointer(&quadrature_cache[cache_index]);
 *           cache_index += quadrature_formula.size();
 *         }
 *     Assert(cache_index == quadrature_cache.size(), ExcInternalError());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavefrequency_sweep"></a> 
 * <h4>ElasticWave::frequency_sweep</h4>
 * 

 * 
 * For clarity we divide the function `run` of step-40 into the functions
 * `run` and `frequency_sweep`. In the function `frequency_sweep` we place the
 * iteration over the frequency vector.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::frequency_sweep()
 *   {
 *     for (unsigned int frequency_idx = 0;
 *          frequency_idx < parameters.nb_frequency_points;
 *          ++frequency_idx)
 *       {
 *         pcout << parameters.simulation_name + " frequency idx: "
 *               << frequency_idx << '/' << parameters.nb_frequency_points - 1
 *               << std::endl;
 * 
 * 
 * 
 *         setup_system();
 *         if (frequency_idx == 0)
 *           {
 *             pcout << "   Number of active cells :       "
 *                   << triangulation.n_active_cells() << std::endl;
 *             pcout << "   Number of degrees of freedom : "
 *                   << dof_handler.n_dofs() << std::endl;
 *           }
 * 
 *         if (frequency_idx == 0)
 *           {
 * @endcode
 * 
 * Write the simulation parameters only once
 * 
 * @code
 *             parameters.data.set_attribute("active_cells",
 *                                           triangulation.n_active_cells());
 *             parameters.data.set_attribute("degrees_of_freedom",
 *                                           dof_handler.n_dofs());
 *           }
 * 
 * @endcode
 * 
 * We calculate the frequency and omega values for this particular step.
 * 
 * @code
 *         const double current_loop_frequency =
 *           (parameters.start_frequency +
 *            frequency_idx *
 *              (parameters.stop_frequency - parameters.start_frequency) /
 *              (parameters.nb_frequency_points - 1));
 *         const double current_loop_omega =
 *           2 * numbers::PI * current_loop_frequency;
 * 
 * @endcode
 * 
 * In the first frequency step we calculate the mass and stiffness
 * matrices and the right hand side. In the subsequent frequency steps
 * we will use those values. This improves considerably the calculation
 * time.
 * 
 * @code
 *         assemble_system(current_loop_omega,
 *                         (frequency_idx == 0) ? true : false);
 *         solve();
 * 
 *         frequency[frequency_idx] = current_loop_frequency;
 *         store_frequency_step_data(frequency_idx);
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 *         pcout << std::endl;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaverun"></a> 
 * <h4>ElasticWave::run</h4>
 * 

 * 
 * This function is very similar to the one in step-40.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::run()
 *   {
 * #ifdef DEBUG
 *     pcout << "Debug mode" << std::endl;
 * #else
 *     pcout << "Release mode" << std::endl;
 * #endif
 * 
 *     {
 *       Point<dim> p1;
 *       p1(0) = -parameters.dimension_x / 2;
 *       p1(1) = -parameters.dimension_y / 2;
 *       if (dim == 3)
 *         {
 *           p1(2) = -parameters.dimension_y / 2;
 *         }
 *       Point<dim> p2;
 *       p2(0) = parameters.dimension_x / 2;
 *       p2(1) = parameters.dimension_y / 2;
 *       if (dim == 3)
 *         {
 *           p2(2) = parameters.dimension_y / 2;
 *         }
 *       std::vector<unsigned int> divisions(dim);
 *       divisions[0] = int(parameters.dimension_x / parameters.dimension_y);
 *       divisions[1] = 1;
 *       if (dim == 3)
 *         {
 *           divisions[2] = 1;
 *         }
 *       GridGenerator::subdivided_hyper_rectangle(triangulation,
 *                                                 divisions,
 *                                                 p1,
 *                                                 p2);
 *     }
 * 
 *     triangulation.refine_global(parameters.grid_level);
 * 
 *     setup_quadrature_cache();
 * 
 *     initialize_probe_positions_vector();
 * 
 *     frequency_sweep();
 * 
 *     output_results();
 *   }
 * } // namespace step62
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h4>The main function</h4>
 * 

 * 
 * The main function is very similar to the one in step-40.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       const unsigned int dim = 2;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       HDF5::File data_file("results.h5",
 *                            HDF5::File::FileAccessMode::create,
 *                            MPI_COMM_WORLD);
 *       auto       data = data_file.create_group("data");
 * 
 * @endcode
 * 
 * Each of the simulations (displacement and calibration) is stored in a
 * separate HDF5 group:
 * 
 * @code
 *       const std::vector<std::string> group_names = {"displacement",
 *                                                     "calibration"};
 *       for (auto group_name : group_names)
 *         {
 * @endcode
 * 
 * For each of these two group names, we now create the group and put
 * attributes into these groups.
 * Specifically, these are:
 * - The dimensions of the waveguide (in $x$ and $y$ directions)
 * - The position of the probe (in $x$ and $y$ directions)
 * - The number of points in the probe
 * - The global refinement level
 * - The cavity resonance frequency
 * - The number of mirror pairs
 * - The material properties
 * - The force parameters
 * - The PML parameters
 * - The frequency parameters
 * 

 * 
 * 
 * @code
 *           auto group = data.create_group(group_name);
 * 
 *           group.set_attribute<double>("dimension_x", 2e-5);
 *           group.set_attribute<double>("dimension_y", 2e-8);
 *           group.set_attribute<double>("probe_pos_x", 8e-6);
 *           group.set_attribute<double>("probe_pos_y", 0);
 *           group.set_attribute<double>("probe_width_y", 2e-08);
 *           group.set_attribute<unsigned int>("nb_probe_points", 5);
 *           group.set_attribute<unsigned int>("grid_level", 1);
 *           group.set_attribute<double>("cavity_resonance_frequency", 20e9);
 *           group.set_attribute<unsigned int>("nb_mirror_pairs", 15);
 * 
 *           group.set_attribute<double>("poissons_ratio", 0.27);
 *           group.set_attribute<double>("youngs_modulus", 270000000000.0);
 *           group.set_attribute<double>("material_a_rho", 3200);
 * 
 *           if (group_name == std::string("displacement"))
 *             group.set_attribute<double>("material_b_rho", 2000);
 *           else
 *             group.set_attribute<double>("material_b_rho", 3200);
 * 
 *           group.set_attribute(
 *             "lambda",
 *             group.get_attribute<double>("youngs_modulus") *
 *               group.get_attribute<double>("poissons_ratio") /
 *               ((1 + group.get_attribute<double>("poissons_ratio")) *
 *                (1 - 2 * group.get_attribute<double>("poissons_ratio"))));
 *           group.set_attribute("mu",
 *                               group.get_attribute<double>("youngs_modulus") /
 *                                 (2 * (1 + group.get_attribute<double>(
 *                                             "poissons_ratio"))));
 * 
 *           group.set_attribute<double>("max_force_amplitude", 1e26);
 *           group.set_attribute<double>("force_sigma_x", 1e-7);
 *           group.set_attribute<double>("force_sigma_y", 1);
 *           group.set_attribute<double>("max_force_width_x", 3e-7);
 *           group.set_attribute<double>("max_force_width_y", 2e-8);
 *           group.set_attribute<double>("force_x_pos", -8e-6);
 *           group.set_attribute<double>("force_y_pos", 0);
 * 
 *           group.set_attribute<bool>("pml_x", true);
 *           group.set_attribute<bool>("pml_y", false);
 *           group.set_attribute<double>("pml_width_x", 1.8e-6);
 *           group.set_attribute<double>("pml_width_y", 5e-7);
 *           group.set_attribute<double>("pml_coeff", 1.6);
 *           group.set_attribute<unsigned int>("pml_coeff_degree", 2);
 * 
 *           group.set_attribute<double>("center_frequency", 20e9);
 *           group.set_attribute<double>("frequency_range", 0.5e9);
 *           group.set_attribute<double>(
 *             "start_frequency",
 *             group.get_attribute<double>("center_frequency") -
 *               group.get_attribute<double>("frequency_range") / 2);
 *           group.set_attribute<double>(
 *             "stop_frequency",
 *             group.get_attribute<double>("center_frequency") +
 *               group.get_attribute<double>("frequency_range") / 2);
 *           group.set_attribute<unsigned int>("nb_frequency_points", 400);
 * 
 *           if (group_name == std::string("displacement"))
 *             group.set_attribute<std::string>(
 *               "simulation_name", std::string("phononic_cavity_displacement"));
 *           else
 *             group.set_attribute<std::string>(
 *               "simulation_name", std::string("phononic_cavity_calibration"));
 * 
 *           group.set_attribute<bool>("save_vtu_files", false);
 *         }
 * 
 *       {
 * @endcode
 * 
 * Displacement simulation. The parameters are read from the
 * displacement HDF5 group and the results are saved in the same HDF5
 * group.
 * 
 * @code
 *         auto                    displacement = data.open_group("displacement");
 *         step62::Parameters<dim> parameters(displacement);
 * 
 *         step62::ElasticWave<dim> elastic_problem(parameters);
 *         elastic_problem.run();
 *       }
 * 
 *       {
 * @endcode
 * 
 * Calibration simulation. The parameters are read from the calibration
 * HDF5 group and the results are saved in the same HDF5 group.
 * 
 * @code
 *         auto                    calibration = data.open_group("calibration");
 *         step62::Parameters<dim> parameters(calibration);
 * 
 *         step62::ElasticWave<dim> elastic_problem(parameters);
 *         elastic_problem.run();
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


<a name="Resonancefrequencyandbandgap"></a><h3>Resonance frequency and bandgap</h3>


The results are analyzed in the
[jupyter notebook](https://github.com/dealii/dealii/blob/phononic-cavity/examples/step-62/step-62.ipynb)
with the following code
@code{.py}
h5_file = h5py.File('results.h5', 'r')
data = h5_file['data']

# Gaussian function that we use to fit the resonance
def resonance_f(freq, freq_m, quality_factor, max_amplitude):
    omega = 2 * constants.pi * freq
    omega_m = 2 * constants.pi * freq_m
    gamma = omega_m / quality_factor
    return max_amplitude * omega_m**2 * gamma**2 / (((omega_m**2 - omega**2)**2 + gamma**2 * omega**2))

frequency = data['displacement']['frequency'][...]
# Average the probe points
displacement = np.mean(data['displacement']['displacement'], axis=0)
calibration_displacement = np.mean(data['calibration']['displacement'], axis=0)
reflection_coefficient = displacement / calibration_displacement
reflectivity = (np.abs(np.mean(data['displacement']['displacement'][...]**2, axis=0))/
                np.abs(np.mean(data['calibration']['displacement'][...]**2, axis=0)))

try:
    x_data = frequency
    y_data = reflectivity
    quality_factor_guess = 1e3
    freq_guess = x_data[np.argmax(y_data)]
    amplitude_guess = np.max(y_data)
    fit_result, covariance = scipy.optimize.curve_fit(resonance_f, x_data, y_data,
                                                      [freq_guess, quality_factor_guess, amplitude_guess])
    freq_m = fit_result[0]
    quality_factor = np.abs(fit_result[1])
    max_amplitude = fit_result[2]
    y_data_fit = resonance_f(x_data, freq_m, quality_factor, max_amplitude)

    fig = plt.figure()
    plt.plot(frequency / 1e9, reflectivity, frequency / 1e9, y_data_fit)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('amplitude^2 (a.u.)')
    plt.title('Transmission\n' + 'freq = ' + "%.7g" % (freq_guess / 1e9) + 'GHz Q = ' + "%.6g" % quality_factor)
except:
    fig = plt.figure()
    plt.plot(frequency / 1e9, reflectivity)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('amplitude^2 (a.u.)')
    plt.title('Transmission')

fig = plt.figure()
plt.plot(frequency / 1e9, np.angle(reflection_coefficient))
plt.xlabel('frequency (GHz)')
plt.ylabel('phase (rad)')
plt.title('Phase (transmission coefficient)\n')

plt.show()
h5_file.close()
@endcode

A phononic cavity is characterized by the
[resonance frequency](https://en.wikipedia.org/wiki/Resonance) and the
[the quality factor](https://en.wikipedia.org/wiki/Q_factor).
The quality factor is equal to the ratio between the stored energy in the resonator and the energy
dissipated energy per cycle, which is approximately equivalent to the ratio between the
resonance frequency and the
[full width at half maximum (FWHM)](https://en.wikipedia.org/wiki/Full_width_at_half_maximum).
The FWHM is equal to the bandwidth over which the power of vibration is greater than half the
power at the resonant frequency.
@f[
Q = \frac{f_r}{\Delta f} = \frac{\omega_r}{\Delta \omega} =
2 \pi \times \frac{\text{energy stored}}{\text{energy dissipated per cycle}}
@f]

The square of the amplitude of the mechanical resonance $a^2$ as a function of the frequency
has a gaussian shape
@f[
a^2 = a_\textrm{max}^2\frac{\omega^2\Gamma^2}{(\omega_r^2-\omega^2)^2+\Gamma^2\omega^2}
@f]
where $f_r = \frac{\omega_r}{2\pi}$ is the resonance frequency and $\Gamma=\frac{\omega_r}{Q}$ is the dissipation rate.
We used the previous equation in the jupyter notebook to fit the mechanical resonance.

Given the values we have chosen for the parameters, one could estimate the resonance frequency
analytically. Indeed, this is then confirmed by what we get in this program:
the phononic superlattice cavity exhibits a mechanical resonance at 20GHz and a quality factor of 5046.
The following images show the transmission amplitude and phase as a function of frequency in the
vicinity of the resonance frequency:

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.05.png" height="400" />
<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.06.png" height="400" />

The images above suggest that the periodic structure has its intended effect: It really only lets waves of a very
specific frequency pass through, whereas all other waves are reflected. This is of course precisely what one builds
these sorts of devices for.
But it is not quite this easy. In practice, there is really only a "band gap", i.e., the device blocks waves other than
the desired one at 20GHz only within a certain frequency range. Indeed, to find out how large this "gap" is within
which waves are blocked, we can extend the frequency range to 16 GHz through the appropriate parameters in the
input file. We then obtain the following image:

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.07.png" height="400" />

What this image suggests is that in the range of around 18 to around 22 GHz, really only the waves with a frequency
of 20 GHz are allowed to pass through, but beyond this range, there are plenty of other frequencies that can pass
through the device.

<a name="Modeprofile"></a><h3>Mode profile</h3>


We can inspect the mode profile with Paraview or VisIt.
As we have discussed, at resonance all the mechanical
energy is transmitted and the amplitude of motion is amplified inside the cavity.
It can be observed that the PMLs are quite effective to truncate the solution.
The following image shows the mode profile at resonance:

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.08.png" height="400" />

On the other hand,  out of resonance all the mechanical energy is
reflected. The following image shows the profile at 19.75 GHz.
Note the interference between the force pulse and the reflected wave
at the position $x=-8\mu\textrm{m}$.

<img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.09.png" height="400" />

<a name="Experimentalapplications"></a><h3>Experimental applications</h3>


Phononic superlattice cavities find application in
[quantum optomechanics](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.86.1391).
Here we have presented the simulation of a 2D superlattice cavity,
but this code can be used as well to simulate "real world" 3D devices such as
[micropillar superlattice cavities](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.060101),
which are promising candidates to study macroscopic quantum phenomena.
The 20GHz mode of a micropillar superlattice cavity is essentially a mechanical harmonic oscillator that is very well isolated
from the environment. If the device is cooled down to 20mK in a dilution fridge, the mode would then become a
macroscopic quantum harmonic oscillator.


<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


Instead of setting the parameters in the C++ file we could set the parameters
using a python script and save them in the HDF5 file that we will use for
the simulations. Then the deal.II program will read the parameters from the
HDF5 file.

@code{.py}
import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import scipy.constants as constants
import scipy.optimize

# This considerably reduces the size of the svg data
plt.rcParams['svg.fonttype'] = 'none'

h5_file = h5py.File('results.h5', 'w')
data = h5_file.create_group('data')
displacement = data.create_group('displacement')
calibration = data.create_group('calibration')

# Set the parameters
for group in [displacement, calibration]:
    # Dimensions of the domain
    # The waveguide length is equal to dimension_x
    group.attrs['dimension_x'] = 2e-5
    # The waveguide width is equal to dimension_y
    group.attrs['dimension_y'] = 2e-8

    # Position of the probe that we use to measure the flux
    group.attrs['probe_pos_x']   = 8e-6
    group.attrs['probe_pos_y']   = 0
    group.attrs['probe_width_y'] = 2e-08

    # Number of points in the probe
    group.attrs['nb_probe_points'] = 5

    # Global refinement
    group.attrs['grid_level'] = 1

    # Cavity
    group.attrs['cavity_resonance_frequency'] = 20e9
    group.attrs['nb_mirror_pairs']            = 15

    # Material
    group.attrs['poissons_ratio'] = 0.27
    group.attrs['youngs_modulus'] = 270000000000.0
    group.attrs['material_a_rho'] = 3200
    if group == displacement:
        group.attrs['material_b_rho'] = 2000
    else:
        group.attrs['material_b_rho'] = 3200
    group.attrs['lambda'] = (group.attrs['youngs_modulus'] * group.attrs['poissons_ratio'] /
                           ((1 + group.attrs['poissons_ratio']) *
                           (1 - 2 * group.attrs['poissons_ratio'])))
    group.attrs['mu']= (group.attrs['youngs_modulus'] / (2 * (1 + group.attrs['poissons_ratio'])))

    # Force
    group.attrs['max_force_amplitude'] = 1e26
    group.attrs['force_sigma_x']       = 1e-7
    group.attrs['force_sigma_y']       = 1
    group.attrs['max_force_width_x']   = 3e-7
    group.attrs['max_force_width_y']   = 2e-8
    group.attrs['force_x_pos']         = -8e-6
    group.attrs['force_y_pos']         = 0

    # PML
    group.attrs['pml_x']            = True
    group.attrs['pml_y']            = False
    group.attrs['pml_width_x']      = 1.8e-6
    group.attrs['pml_width_y']      = 5e-7
    group.attrs['pml_coeff']        = 1.6
    group.attrs['pml_coeff_degree'] = 2

    # Frequency sweep
    group.attrs['center_frequency']    = 20e9
    group.attrs['frequency_range']     = 0.5e9
    group.attrs['start_frequency']     = group.attrs['center_frequency'] - group.attrs['frequency_range'] / 2
    group.attrs['stop_frequency']      = group.attrs['center_frequency'] + group.attrs['frequency_range'] / 2
    group.attrs['nb_frequency_points'] = 400

    # Other parameters
    if group == displacement:
        group.attrs['simulation_name'] = 'phononic_cavity_displacement'
    else:
        group.attrs['simulation_name'] = 'phononic_cavity_calibration'
    group.attrs['save_vtu_files'] = False

h5_file.close()
@endcode

In order to read the HDF5 parameters we have to use the
HDF5::File::FileAccessMode::open flag.
@code{.py}
      HDF5::File data_file("results.h5",
                           HDF5::File::FileAccessMode::open,
                           MPI_COMM_WORLD);
      auto       data = data_file.open_group("data");
@endcode
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-62.cc"
*/
