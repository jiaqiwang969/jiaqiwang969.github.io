/**
@page step_39 The step-39 tutorial program
This tutorial depends on step-12b.

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
        <li><a href="#Thelocalintegrators">The local integrators</a>
        <li><a href="#Themainclass">The main class</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Logfileoutput">Logfile output</a>
        <li><a href="#Postprocessingofthelogfile">Postprocessing of the logfile</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>

In this program, we use the interior penalty method and Nitsche's weak
boundary conditions to solve Poisson's equation. We use multigrid
methods on locally refined meshes, which are generated using a bulk
criterion and a standard error estimator based on cell and face
residuals. All operators are implemented using the MeshWorker interface.

Like in step-12, the discretization relies on finite element spaces,
which are polynomial inside the mesh cells $K\in \mathbb T_h$, but
have no continuity between cells. Since such functions have two values
on each interior face $F\in \mathbb F_h^i$, one from each side, we
define mean value and jump operators as follows: let
<i>K</i><sub>1</sub> and <i>K</i><sub>2</sub> be the two cells sharing
a face, and let the traces of functions <i>u<sub>i</sub></i> and the
outer normal vectors <b>n</b><i><sub>i</sub></i> be labeled
accordingly. Then, on the face, we let
@f[
	\average{ u } = \frac{u_1 + u_2}2
@f]

Note, that if such an expression contains a normal vector, the
averaging operator turns into a jump. The interior penalty method for the problem
@f[
  -\Delta u = f \text{ in }\Omega \qquad u = u^D \text{ on } \partial\Omega
@f]
becomes
@f{multline*}
  \sum_{K\in \mathbb T_h} (\nabla u, \nabla v)_K
  \\
  + \sum_{F \in F_h^i} \biggl\{4\sigma_F (\average{ u \mathbf n}, \average{ v \mathbf n })_F
  - 2 (\average{ \nabla u },\average{ v\mathbf n })_F
  - 2 (\average{ \nabla v },\average{ u\mathbf n })_F
  \biggr\}
  \\
  + \sum_{F \in F_h^b} \biggl\{2\sigma_F (u, v)_F
  - (\partial_n u,v)_F
  - (\partial_n v,u)_F
  \biggr\}
  \\
  = (f, v)_\Omega + \sum_{F \in F_h^b} \biggl\{
  2\sigma_F (u^D, v)_F - (\partial_n v,u^D)_F
  \biggr\}.
@f}

Here, $\sigma_F$ is the penalty parameter, which is chosen as follows:
for a face <i>F</i> of a cell <i>K</i>, compute the value
@f[
\sigma_{F,K} = p(p+1) \frac{|F|_{d-1}}{|K|_d},
@f]
where <i>p</i> is the polynomial degree of the finite element
functions and $|\cdot|_d$ and $|\cdot|_{d-1}$ denote the $d$ and $d-1$
dimensional Hausdorff measure of the corresponding
object. If the face is at the boundary, choose $\sigma_F = \sigma_{F,K}$.
For an interior face, we take the average of the two values at this face.

In our finite element program, we distinguish three different
integrals, corresponding to the sums over cells, interior faces and
boundary faces above. Since the MeshWorker::loop organizes the sums
for us, we only need to implement the integrals over each mesh
element. The class MatrixIntegrator below has these three functions
for the left hand side of the formula, the class RHSIntegrator for the
right.

As we will see below, even the error estimate is of the same
structure, since it can be written as
@f{align*}
  \eta^2 &= \eta_K^2 + \eta_F^2 + \eta_B^2
  \\
  \eta_K^2 &= \sum_{K\in \mathbb T_h} h^2 \|f + \Delta u_h\|^2
  \\
  \eta_F^2 &= \sum_{F \in F_h^i} \biggl\{
    4 \sigma_F \| \average{u_h\mathbf n} \|^2 + h \|\average{\partial_n u_h}\|^2 \biggr\}
  \\
  \eta_B^2 &= \sum_{F \in F_h^b} 2\sigma_F \| u_h-u^D \|^2.
@f}

Thus, the functions for assembling matrices, right hand side and error
estimates below exhibit that these loops are all generic and can be
programmed in the same way.

This program is related to step-12b, in that it uses MeshWorker and
discontinuous Galerkin methods. While there, we solved an advection
problem, here it is a diffusion problem. Here, we also use multigrid
preconditioning and a theoretically justified error estimator, see
Karakashian and Pascal (2003). The multilevel scheme was discussed in
detail in Kanschat (2004). The adaptive iteration and its convergence
have been discussed (for triangular meshes) in Hoppe, Kanschat, and
Warburton (2009).
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The include files for the linear algebra: A regular SparseMatrix, which in
 * turn will include the necessary files for SparsityPattern and Vector
 * classes.
 * 
 * @code
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/precondition_block.h>
 * #include <deal.II/lac/block_vector.h>
 * 
 * @endcode
 * 
 * Include files for setting up the mesh
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * @endcode
 * 
 * Include files for FiniteElement classes and DoFHandler.
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_dgp.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * The include files for using the MeshWorker framework
 * 
 * @code
 * #include <deal.II/meshworker/dof_info.h>
 * #include <deal.II/meshworker/integration_info.h>
 * #include <deal.II/meshworker/assembler.h>
 * #include <deal.II/meshworker/loop.h>
 * 
 * @endcode
 * 
 * The include file for local integrators associated with the Laplacian
 * 
 * @code
 * #include <deal.II/integrators/laplace.h>
 * 
 * @endcode
 * 
 * Support for multigrid methods
 * 
 * @code
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * #include <deal.II/multigrid/mg_transfer.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * 
 * @endcode
 * 
 * Finally, we take our exact solution from the library as well as quadrature
 * and additional tools.
 * 
 * @code
 * #include <deal.II/base/function_lib.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * @endcode
 * 
 * All classes of the deal.II library are in the namespace dealii. In order to
 * save typing, we tell the compiler to search names in there as well.
 * 
 * @code
 * namespace Step39
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * This is the function we use to set the boundary values and also the exact
 * solution we compare to.
 * 
 * @code
 *   Functions::SlitSingularityFunction<2> exact_solution;
 * 
 * @endcode
 * 
 * 
 * <a name="Thelocalintegrators"></a> 
 * <h3>The local integrators</h3>
 * 

 * 
 * MeshWorker separates local integration from the loops over cells and
 * faces. Thus, we have to write local integration classes for generating
 * matrices, the right hand side and the error estimator.
 * 

 * 
 * All these classes have the same three functions for integrating over
 * cells, boundary faces and interior faces, respectively. All the
 * information needed for the local integration is provided by
 * MeshWorker::IntegrationInfo<dim>. Note that the signature of the
 * functions cannot be changed, because it is expected by
 * MeshWorker::integration_loop().
 * 

 * 
 * The first class defining local integrators is responsible for computing
 * cell and face matrices. It is used to assemble the global matrix as well
 * as the level matrices.
 * 
 * @code
 *   template <int dim>
 *   class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
 *   {
 *   public:
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo,
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1,
 *               MeshWorker::DoFInfo<dim> &                 dinfo2,
 *               typename MeshWorker::IntegrationInfo<dim> &info1,
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * On each cell, we integrate the Dirichlet form. We use the library of
 * ready made integrals in LocalIntegrators to avoid writing these loops
 * ourselves. Similarly, we implement Nitsche boundary conditions and the
 * interior penalty fluxes between cells.
 *   

 * 
 * The boundary and flux terms need a penalty parameter, which should be
 * adjusted to the cell size and the polynomial degree. A safe choice of
 * this parameter for constant coefficients can be found in
 * LocalIntegrators::Laplace::compute_penalty() and we use this below.
 * 
 * @code
 *   template <int dim>
 *   void MatrixIntegrator<dim>::cell(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix,
 *                                            info.fe_values());
 *   }
 * 
 * 
 *   template <int dim>
 *   void MatrixIntegrator<dim>::boundary(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const unsigned int degree = info.fe_values(0).get_fe().tensor_degree();
 *     LocalIntegrators::Laplace::nitsche_matrix(
 *       dinfo.matrix(0, false).matrix,
 *       info.fe_values(0),
 *       LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, degree, degree));
 *   }
 * 
 * @endcode
 * 
 * Interior faces use the interior penalty method
 * 
 * @code
 *   template <int dim>
 *   void MatrixIntegrator<dim>::face(
 *     MeshWorker::DoFInfo<dim> &                 dinfo1,
 *     MeshWorker::DoFInfo<dim> &                 dinfo2,
 *     typename MeshWorker::IntegrationInfo<dim> &info1,
 *     typename MeshWorker::IntegrationInfo<dim> &info2) const
 *   {
 *     const unsigned int degree = info1.fe_values(0).get_fe().tensor_degree();
 *     LocalIntegrators::Laplace::ip_matrix(
 *       dinfo1.matrix(0, false).matrix,
 *       dinfo1.matrix(0, true).matrix,
 *       dinfo2.matrix(0, true).matrix,
 *       dinfo2.matrix(0, false).matrix,
 *       info1.fe_values(0),
 *       info2.fe_values(0),
 *       LocalIntegrators::Laplace::compute_penalty(
 *         dinfo1, dinfo2, degree, degree));
 *   }
 * 
 * @endcode
 * 
 * The second local integrator builds the right hand side. In our example,
 * the right hand side function is zero, such that only the boundary
 * condition is set here in weak form.
 * 
 * @code
 *   template <int dim>
 *   class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
 *   {
 *   public:
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo,
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1,
 *               MeshWorker::DoFInfo<dim> &                 dinfo2,
 *               typename MeshWorker::IntegrationInfo<dim> &info1,
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   void
 *   RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &,
 *                            typename MeshWorker::IntegrationInfo<dim> &) const
 *   {}
 * 
 * 
 *   template <int dim>
 *   void RHSIntegrator<dim>::boundary(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const FEValuesBase<dim> &fe           = info.fe_values();
 *     Vector<double> &         local_vector = dinfo.vector(0).block(0);
 * 
 *     std::vector<double> boundary_values(fe.n_quadrature_points);
 *     exact_solution.value_list(fe.get_quadrature_points(), boundary_values);
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree();
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() /
 *                            dinfo.cell->measure();
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
 *         local_vector(i) +=
 *           (-penalty * fe.shape_value(i, k)              // (-sigma * v_i(x_k)
 *            + fe.normal_vector(k) * fe.shape_grad(i, k)) // + n * grad v_i(x_k))
 *           * boundary_values[k] * fe.JxW(k);             // u^D(x_k) * dx
 *   }
 * 
 * 
 *   template <int dim>
 *   void
 *   RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
 *                            MeshWorker::DoFInfo<dim> &,
 *                            typename MeshWorker::IntegrationInfo<dim> &,
 *                            typename MeshWorker::IntegrationInfo<dim> &) const
 *   {}
 * 
 * 
 * @endcode
 * 
 * The third local integrator is responsible for the contributions to the
 * error estimate. This is the standard energy estimator due to Karakashian
 * and Pascal (2003).
 * 
 * @code
 *   template <int dim>
 *   class Estimator : public MeshWorker::LocalIntegrator<dim>
 *   {
 *   public:
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo,
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1,
 *               MeshWorker::DoFInfo<dim> &                 dinfo2,
 *               typename MeshWorker::IntegrationInfo<dim> &info1,
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * The cell contribution is the Laplacian of the discrete solution, since
 * the right hand side is zero.
 * 
 * @code
 *   template <int dim>
 *   void
 *   Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &                 dinfo,
 *                        typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const FEValuesBase<dim> &fe = info.fe_values();
 * 
 *     const std::vector<Tensor<2, dim>> &DDuh = info.hessians[0][0];
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         const double t = dinfo.cell->diameter() * trace(DDuh[k]);
 *         dinfo.value(0) += t * t * fe.JxW(k);
 *       }
 *     dinfo.value(0) = std::sqrt(dinfo.value(0));
 *   }
 * 
 * @endcode
 * 
 * At the boundary, we use simply a weighted form of the boundary residual,
 * namely the norm of the difference between the finite element solution and
 * the correct boundary condition.
 * 
 * @code
 *   template <int dim>
 *   void Estimator<dim>::boundary(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const FEValuesBase<dim> &fe = info.fe_values();
 * 
 *     std::vector<double> boundary_values(fe.n_quadrature_points);
 *     exact_solution.value_list(fe.get_quadrature_points(), boundary_values);
 * 
 *     const std::vector<double> &uh = info.values[0][0];
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree();
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() /
 *                            dinfo.cell->measure();
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         const double diff = boundary_values[k] - uh[k];
 *         dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
 *       }
 *     dinfo.value(0) = std::sqrt(dinfo.value(0));
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally, on interior faces, the estimator consists of the jumps of the
 * solution and its normal derivative, weighted appropriately.
 * 
 * @code
 *   template <int dim>
 *   void
 *   Estimator<dim>::face(MeshWorker::DoFInfo<dim> &                 dinfo1,
 *                        MeshWorker::DoFInfo<dim> &                 dinfo2,
 *                        typename MeshWorker::IntegrationInfo<dim> &info1,
 *                        typename MeshWorker::IntegrationInfo<dim> &info2) const
 *   {
 *     const FEValuesBase<dim> &          fe   = info1.fe_values();
 *     const std::vector<double> &        uh1  = info1.values[0][0];
 *     const std::vector<double> &        uh2  = info2.values[0][0];
 *     const std::vector<Tensor<1, dim>> &Duh1 = info1.gradients[0][0];
 *     const std::vector<Tensor<1, dim>> &Duh2 = info2.gradients[0][0];
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree();
 *     const double       penalty1 =
 *       degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
 *     const double penalty2 =
 *       degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
 *     const double penalty = penalty1 + penalty2;
 *     const double h       = dinfo1.face->measure();
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         const double diff1 = uh1[k] - uh2[k];
 *         const double diff2 =
 *           fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
 *         dinfo1.value(0) +=
 *           (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k);
 *       }
 *     dinfo1.value(0) = std::sqrt(dinfo1.value(0));
 *     dinfo2.value(0) = dinfo1.value(0);
 *   }
 * 
 * @endcode
 * 
 * Finally we have an integrator for the error. Since the energy norm for
 * discontinuous Galerkin problems not only involves the difference of the
 * gradient inside the cells, but also the jump terms across faces and at
 * the boundary, we cannot just use VectorTools::integrate_difference().
 * Instead, we use the MeshWorker interface to compute the error ourselves.
 * 

 * 
 * There are several different ways to define this energy norm, but all of
 * them are equivalent to each other uniformly with mesh size (some not
 * uniformly with polynomial degree). Here, we choose @f[ \|u\|_{1,h} =
 * \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
 * 4\sigma_F\|\average{ u \mathbf n}\|^2_F + \sum_{F \in F_h^b}
 * 2\sigma_F\|u\|^2_F @f]
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim>
 *   {
 *   public:
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo,
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override;
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1,
 *               MeshWorker::DoFInfo<dim> &                 dinfo2,
 *               typename MeshWorker::IntegrationInfo<dim> &info1,
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override;
 *   };
 * 
 * @endcode
 * 
 * Here we have the integration on cells. There is currently no good
 * interface in MeshWorker that would allow us to access values of regular
 * functions in the quadrature points. Thus, we have to create the vectors
 * for the exact function's values and gradients inside the cell
 * integrator. After that, everything is as before and we just add up the
 * squares of the differences.
 * 

 * 
 * Additionally to computing the error in the energy norm, we use the
 * capability of the mesh worker to compute two functionals at the same time
 * and compute the <i>L<sup>2</sup></i>-error in the same loop. Obviously,
 * this one does not have any jump terms and only appears in the integration
 * on cells.
 * 
 * @code
 *   template <int dim>
 *   void ErrorIntegrator<dim>::cell(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const FEValuesBase<dim> &   fe = info.fe_values();
 *     std::vector<Tensor<1, dim>> exact_gradients(fe.n_quadrature_points);
 *     std::vector<double>         exact_values(fe.n_quadrature_points);
 * 
 *     exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients);
 *     exact_solution.value_list(fe.get_quadrature_points(), exact_values);
 * 
 *     const std::vector<Tensor<1, dim>> &Duh = info.gradients[0][0];
 *     const std::vector<double> &        uh  = info.values[0][0];
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         double sum = 0;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           {
 *             const double diff = exact_gradients[k][d] - Duh[k][d];
 *             sum += diff * diff;
 *           }
 *         const double diff = exact_values[k] - uh[k];
 *         dinfo.value(0) += sum * fe.JxW(k);
 *         dinfo.value(1) += diff * diff * fe.JxW(k);
 *       }
 *     dinfo.value(0) = std::sqrt(dinfo.value(0));
 *     dinfo.value(1) = std::sqrt(dinfo.value(1));
 *   }
 * 
 * 
 *   template <int dim>
 *   void ErrorIntegrator<dim>::boundary(
 *     MeshWorker::DoFInfo<dim> &                 dinfo,
 *     typename MeshWorker::IntegrationInfo<dim> &info) const
 *   {
 *     const FEValuesBase<dim> &fe = info.fe_values();
 * 
 *     std::vector<double> exact_values(fe.n_quadrature_points);
 *     exact_solution.value_list(fe.get_quadrature_points(), exact_values);
 * 
 *     const std::vector<double> &uh = info.values[0][0];
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree();
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() /
 *                            dinfo.cell->measure();
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         const double diff = exact_values[k] - uh[k];
 *         dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
 *       }
 *     dinfo.value(0) = std::sqrt(dinfo.value(0));
 *   }
 * 
 * 
 *   template <int dim>
 *   void ErrorIntegrator<dim>::face(
 *     MeshWorker::DoFInfo<dim> &                 dinfo1,
 *     MeshWorker::DoFInfo<dim> &                 dinfo2,
 *     typename MeshWorker::IntegrationInfo<dim> &info1,
 *     typename MeshWorker::IntegrationInfo<dim> &info2) const
 *   {
 *     const FEValuesBase<dim> &  fe  = info1.fe_values();
 *     const std::vector<double> &uh1 = info1.values[0][0];
 *     const std::vector<double> &uh2 = info2.values[0][0];
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree();
 *     const double       penalty1 =
 *       degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
 *     const double penalty2 =
 *       degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
 *     const double penalty = penalty1 + penalty2;
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
 *       {
 *         const double diff = uh1[k] - uh2[k];
 *         dinfo1.value(0) += (penalty * diff * diff) * fe.JxW(k);
 *       }
 *     dinfo1.value(0) = std::sqrt(dinfo1.value(0));
 *     dinfo2.value(0) = dinfo1.value(0);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * This class does the main job, like in previous examples. For a
 * description of the functions declared here, please refer to the
 * implementation below.
 * 
 * @code
 *   template <int dim>
 *   class InteriorPenaltyProblem
 *   {
 *   public:
 *     using CellInfo = MeshWorker::IntegrationInfo<dim>;
 * 
 *     InteriorPenaltyProblem(const FiniteElement<dim> &fe);
 * 
 *     void run(unsigned int n_steps);
 * 
 *   private:
 *     void   setup_system();
 *     void   assemble_matrix();
 *     void   assemble_mg_matrix();
 *     void   assemble_right_hand_side();
 *     void   error();
 *     double estimate();
 *     void   solve();
 *     void   output_results(const unsigned int cycle) const;
 * 
 * @endcode
 * 
 * The member objects related to the discretization are here.
 * 
 * @code
 *     Triangulation<dim>        triangulation;
 *     const MappingQ1<dim>      mapping;
 *     const FiniteElement<dim> &fe;
 *     DoFHandler<dim>           dof_handler;
 * 
 * @endcode
 * 
 * Then, we have the matrices and vectors related to the global discrete
 * system.
 * 
 * @code
 *     SparsityPattern      sparsity;
 *     SparseMatrix<double> matrix;
 *     Vector<double>       solution;
 *     Vector<double>       right_hand_side;
 *     BlockVector<double>  estimates;
 * 
 * @endcode
 * 
 * Finally, we have a group of sparsity patterns and sparse matrices
 * related to the multilevel preconditioner.  First, we have a level
 * matrix and its sparsity pattern.
 * 
 * @code
 *     MGLevelObject<SparsityPattern>      mg_sparsity;
 *     MGLevelObject<SparseMatrix<double>> mg_matrix;
 * 
 * @endcode
 * 
 * When we perform multigrid with local smoothing on locally refined
 * meshes, additional matrices are required; see Kanschat (2004). Here is
 * the sparsity pattern for these edge matrices. We only need one, because
 * the pattern of the up matrix is the transpose of that of the down
 * matrix. Actually, we do not care too much about these details, since
 * the MeshWorker is filling these matrices.
 * 
 * @code
 *     MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
 * @endcode
 * 
 * The flux matrix at the refinement edge, coupling fine level degrees of
 * freedom to coarse level.
 * 
 * @code
 *     MGLevelObject<SparseMatrix<double>> mg_matrix_dg_down;
 * @endcode
 * 
 * The transpose of the flux matrix at the refinement edge, coupling
 * coarse level degrees of freedom to fine level.
 * 
 * @code
 *     MGLevelObject<SparseMatrix<double>> mg_matrix_dg_up;
 *   };
 * 
 * 
 * @endcode
 * 
 * The constructor simply sets up the coarse grid and the DoFHandler. The
 * FiniteElement is provided as a parameter to allow flexibility.
 * 
 * @code
 *   template <int dim>
 *   InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(
 *     const FiniteElement<dim> &fe)
 *     : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     , mapping()
 *     , fe(fe)
 *     , dof_handler(triangulation)
 *     , estimates(1)
 *   {
 *     GridGenerator::hyper_cube_slit(triangulation, -1, 1);
 *   }
 * 
 * 
 * @endcode
 * 
 * In this function, we set up the dimension of the linear system and the
 * sparsity patterns for the global matrix as well as the level matrices.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::setup_system()
 *   {
 * @endcode
 * 
 * First, we use the finite element to distribute degrees of freedom over
 * the mesh and number them.
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler.distribute_mg_dofs();
 *     unsigned int n_dofs = dof_handler.n_dofs();
 * @endcode
 * 
 * Then, we already know the size of the vectors representing finite
 * element functions.
 * 
 * @code
 *     solution.reinit(n_dofs);
 *     right_hand_side.reinit(n_dofs);
 * 
 * @endcode
 * 
 * Next, we set up the sparsity pattern for the global matrix. Since we do
 * not know the row sizes in advance, we first fill a temporary
 * DynamicSparsityPattern object and copy it to the regular
 * SparsityPattern once it is complete.
 * 
 * @code
 *     DynamicSparsityPattern dsp(n_dofs);
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
 *     sparsity.copy_from(dsp);
 *     matrix.reinit(sparsity);
 * 
 *     const unsigned int n_levels = triangulation.n_levels();
 * @endcode
 * 
 * The global system is set up, now we attend to the level matrices. We
 * resize all matrix objects to hold one matrix per level.
 * 
 * @code
 *     mg_matrix.resize(0, n_levels - 1);
 *     mg_matrix.clear_elements();
 *     mg_matrix_dg_up.resize(0, n_levels - 1);
 *     mg_matrix_dg_up.clear_elements();
 *     mg_matrix_dg_down.resize(0, n_levels - 1);
 *     mg_matrix_dg_down.clear_elements();
 * @endcode
 * 
 * It is important to update the sparsity patterns after <tt>clear()</tt>
 * was called for the level matrices, since the matrices lock the sparsity
 * pattern through the SmartPointer and Subscriptor mechanism.
 * 
 * @code
 *     mg_sparsity.resize(0, n_levels - 1);
 *     mg_sparsity_dg_interface.resize(0, n_levels - 1);
 * 
 * @endcode
 * 
 * Now all objects are prepared to hold one sparsity pattern or matrix per
 * level. What's left is setting up the sparsity patterns on each level.
 * 
 * @code
 *     for (unsigned int level = mg_sparsity.min_level();
 *          level <= mg_sparsity.max_level();
 *          ++level)
 *       {
 * @endcode
 * 
 * These are roughly the same lines as above for the global matrix,
 * now for each level.
 * 
 * @code
 *         DynamicSparsityPattern dsp(dof_handler.n_dofs(level));
 *         MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level);
 *         mg_sparsity[level].copy_from(dsp);
 *         mg_matrix[level].reinit(mg_sparsity[level]);
 * 
 * @endcode
 * 
 * Additionally, we need to initialize the transfer matrices at the
 * refinement edge between levels. They are stored at the index
 * referring to the finer of the two indices, thus there is no such
 * object on level 0.
 * 
 * @code
 *         if (level > 0)
 *           {
 *             DynamicSparsityPattern dsp;
 *             dsp.reinit(dof_handler.n_dofs(level - 1),
 *                        dof_handler.n_dofs(level));
 *             MGTools::make_flux_sparsity_pattern_edge(dof_handler, dsp, level);
 *             mg_sparsity_dg_interface[level].copy_from(dsp);
 *             mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
 *             mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
 *           }
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * In this function, we assemble the global system matrix, where by global
 * we indicate that this is the matrix of the discrete system we solve and
 * it is covering the whole mesh.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::assemble_matrix()
 *   {
 * @endcode
 * 
 * First, we need t set up the object providing the values we
 * integrate. This object contains all FEValues and FEFaceValues objects
 * needed and also maintains them automatically such that they always
 * point to the current cell. To this end, we need to tell it first, where
 * and what to compute. Since we are not doing anything fancy, we can rely
 * on their standard choice for quadrature rules.
 *     

 * 
 * Since their default update flags are minimal, we add what we need
 * additionally, namely the values and gradients of shape functions on all
 * objects (cells, boundary and interior faces). Afterwards, we are ready
 * to initialize the container, which will create all necessary
 * FEValuesBase objects for integration.
 * 
 * @code
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     UpdateFlags update_flags = update_values | update_gradients;
 *     info_box.add_update_flags_all(update_flags);
 *     info_box.initialize(fe, mapping);
 * 
 * @endcode
 * 
 * This is the object into which we integrate local data. It is filled by
 * the local integration routines in MatrixIntegrator and then used by the
 * assembler to distribute the information into the global matrix.
 * 
 * @code
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 * @endcode
 * 
 * Furthermore, we need an object that assembles the local matrix into the
 * global matrix. These assembler objects have all the knowledge
 * of the structures of the target object, in this case a
 * SparseMatrix, possible constraints and the mesh structure.
 * 
 * @code
 *     MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler;
 *     assembler.initialize(matrix);
 * 
 * @endcode
 * 
 * Now comes the part we coded ourselves, the local
 * integrator. This is the only part which is problem dependent.
 * 
 * @code
 *     MatrixIntegrator<dim> integrator;
 * @endcode
 * 
 * Now, we throw everything into a MeshWorker::loop(), which here
 * traverses all active cells of the mesh, computes cell and face matrices
 * and assembles them into the global matrix. We use the variable
 * <tt>dof_handler</tt> here in order to use the global numbering of
 * degrees of freedom.
 * 
 * @code
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
 *                                            dof_handler.end(),
 *                                            dof_info,
 *                                            info_box,
 *                                            integrator,
 *                                            assembler);
 *   }
 * 
 * 
 * @endcode
 * 
 * Now, we do the same for the level matrices. Not too surprisingly, this
 * function looks like a twin of the previous one. Indeed, there are only
 * two minor differences.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::assemble_mg_matrix()
 *   {
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     UpdateFlags update_flags = update_values | update_gradients;
 *     info_box.add_update_flags_all(update_flags);
 *     info_box.initialize(fe, mapping);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 * @endcode
 * 
 * Obviously, the assembler needs to be replaced by one filling level
 * matrices. Note that it automatically fills the edge matrices as well.
 * 
 * @code
 *     MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler;
 *     assembler.initialize(mg_matrix);
 *     assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down);
 * 
 *     MatrixIntegrator<dim> integrator;
 * @endcode
 * 
 * Here is the other difference to the previous function: we run
 * over all cells, not only the active ones. And we use functions
 * ending on <code>_mg</code> since we need the degrees of freedom
 * on each level, not the global numbering.
 * 
 * @code
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(),
 *                                            dof_handler.end_mg(),
 *                                            dof_info,
 *                                            info_box,
 *                                            integrator,
 *                                            assembler);
 *   }
 * 
 * 
 * @endcode
 * 
 * Here we have another clone of the assemble function. The difference to
 * assembling the system matrix consists in that we assemble a vector here.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::assemble_right_hand_side()
 *   {
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     UpdateFlags                         update_flags =
 *       update_quadrature_points | update_values | update_gradients;
 *     info_box.add_update_flags_all(update_flags);
 *     info_box.initialize(fe, mapping);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 * @endcode
 * 
 * Since this assembler allows us to fill several vectors, the interface is
 * a little more complicated as above. The pointers to the vectors have to
 * be stored in an AnyData object. While this seems to cause two extra
 * lines of code here, it actually comes handy in more complex
 * applications.
 * 
 * @code
 *     MeshWorker::Assembler::ResidualSimple<Vector<double>> assembler;
 *     AnyData                                               data;
 *     data.add<Vector<double> *>(&right_hand_side, "RHS");
 *     assembler.initialize(data);
 * 
 *     RHSIntegrator<dim> integrator;
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
 *                                            dof_handler.end(),
 *                                            dof_info,
 *                                            info_box,
 *                                            integrator,
 *                                            assembler);
 * 
 *     right_hand_side *= -1.;
 *   }
 * 
 * 
 * @endcode
 * 
 * Now that we have coded all functions building the discrete linear system,
 * it is about time that we actually solve it.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::solve()
 *   {
 * @endcode
 * 
 * The solver of choice is conjugate gradient.
 * 
 * @code
 *     SolverControl            control(1000, 1.e-12);
 *     SolverCG<Vector<double>> solver(control);
 * 
 * @endcode
 * 
 * Now we are setting up the components of the multilevel
 * preconditioner. First, we need transfer between grid levels. The object
 * we are using here generates sparse matrices for these transfers.
 * 
 * @code
 *     MGTransferPrebuilt<Vector<double>> mg_transfer;
 *     mg_transfer.build(dof_handler);
 * 
 * @endcode
 * 
 * Then, we need an exact solver for the matrix on the coarsest level.
 * 
 * @code
 *     FullMatrix<double> coarse_matrix;
 *     coarse_matrix.copy_from(mg_matrix[0]);
 *     MGCoarseGridHouseholder<double, Vector<double>> mg_coarse;
 *     mg_coarse.initialize(coarse_matrix);
 * 
 * @endcode
 * 
 * While transfer and coarse grid solver are pretty much generic, more
 * flexibility is offered for the smoother. First, we choose Gauss-Seidel
 * as our smoothing method.
 * 
 * @code
 *     GrowingVectorMemory<Vector<double>> mem;
 *     using RELAXATION = PreconditionSOR<SparseMatrix<double>>;
 *     mg::SmootherRelaxation<RELAXATION, Vector<double>> mg_smoother;
 *     RELAXATION::AdditionalData                         smoother_data(1.);
 *     mg_smoother.initialize(mg_matrix, smoother_data);
 * 
 * @endcode
 * 
 * Do two smoothing steps on each level.
 * 
 * @code
 *     mg_smoother.set_steps(2);
 * @endcode
 * 
 * Since the SOR method is not symmetric, but we use conjugate gradient
 * iteration below, here is a trick to make the multilevel preconditioner
 * a symmetric operator even for nonsymmetric smoothers.
 * 
 * @code
 *     mg_smoother.set_symmetric(true);
 * @endcode
 * 
 * The smoother class optionally implements the variable V-cycle, which we
 * do not want here.
 * 
 * @code
 *     mg_smoother.set_variable(false);
 * 
 * @endcode
 * 
 * Finally, we must wrap our matrices in an object having the required
 * multiplication functions.
 * 
 * @code
 *     mg::Matrix<Vector<double>> mgmatrix(mg_matrix);
 *     mg::Matrix<Vector<double>> mgdown(mg_matrix_dg_down);
 *     mg::Matrix<Vector<double>> mgup(mg_matrix_dg_up);
 * 
 * @endcode
 * 
 * Now, we are ready to set up the V-cycle operator and the multilevel
 * preconditioner.
 * 
 * @code
 *     Multigrid<Vector<double>> mg(
 *       mgmatrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
 * @endcode
 * 
 * Let us not forget the edge matrices needed because of the adaptive
 * refinement.
 * 
 * @code
 *     mg.set_edge_flux_matrices(mgdown, mgup);
 * 
 * @endcode
 * 
 * After all preparations, wrap the Multigrid object into another object,
 * which can be used as a regular preconditioner,
 * 
 * @code
 *     PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
 *       preconditioner(dof_handler, mg, mg_transfer);
 * @endcode
 * 
 * and use it to solve the system.
 * 
 * @code
 *     solver.solve(matrix, solution, right_hand_side, preconditioner);
 *   }
 * 
 * 
 * @endcode
 * 
 * Another clone of the assemble function. The big difference to the
 * previous ones is here that we also have an input vector.
 * 
 * @code
 *   template <int dim>
 *   double InteriorPenaltyProblem<dim>::estimate()
 *   {
 * @endcode
 * 
 * The results of the estimator are stored in a vector with one entry per
 * cell. Since cells in deal.II are not numbered, we have to create our
 * own numbering in order to use this vector. For the assembler used below
 * the information in which component of a vector the result is stored is
 * transmitted by the user_index variable for each cell. We need to set this
 * numbering up here.
 *     

 * 
 * On the other hand, somebody might have used the user indices
 * already. So, let's be good citizens and save them before tampering with
 * them.
 * 
 * @code
 *     std::vector<unsigned int> old_user_indices;
 *     triangulation.save_user_indices(old_user_indices);
 * 
 *     estimates.block(0).reinit(triangulation.n_active_cells());
 *     unsigned int i = 0;
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       cell->set_user_index(i++);
 * 
 * @endcode
 * 
 * This starts like before,
 * 
 * @code
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     const unsigned int                  n_gauss_points =
 *       dof_handler.get_fe().tensor_degree() + 1;
 *     info_box.initialize_gauss_quadrature(n_gauss_points,
 *                                          n_gauss_points + 1,
 *                                          n_gauss_points);
 * 
 * @endcode
 * 
 * but now we need to notify the info box of the finite element function we
 * want to evaluate in the quadrature points. First, we create an AnyData
 * object with this vector, which is the solution we just computed.
 * 
 * @code
 *     AnyData solution_data;
 *     solution_data.add<const Vector<double> *>(&solution, "solution");
 * 
 * @endcode
 * 
 * Then, we tell the Meshworker::VectorSelector for cells, that we need
 * the second derivatives of this solution (to compute the
 * Laplacian). Therefore, the Boolean arguments selecting function values
 * and first derivatives a false, only the last one selecting second
 * derivatives is true.
 * 
 * @code
 *     info_box.cell_selector.add("solution", false, false, true);
 * @endcode
 * 
 * On interior and boundary faces, we need the function values and the
 * first derivatives, but not second derivatives.
 * 
 * @code
 *     info_box.boundary_selector.add("solution", true, true, false);
 *     info_box.face_selector.add("solution", true, true, false);
 * 
 * @endcode
 * 
 * And we continue as before, with the exception that the default update
 * flags are already adjusted to the values and derivatives we requested
 * above.
 * 
 * @code
 *     info_box.add_update_flags_boundary(update_quadrature_points);
 *     info_box.initialize(fe, mapping, solution_data, solution);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 * @endcode
 * 
 * The assembler stores one number per cell, but else this is the same as
 * in the computation of the right hand side.
 * 
 * @code
 *     MeshWorker::Assembler::CellsAndFaces<double> assembler;
 *     AnyData                                      out_data;
 *     out_data.add<BlockVector<double> *>(&estimates, "cells");
 *     assembler.initialize(out_data, false);
 * 
 *     Estimator<dim> integrator;
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
 *                                            dof_handler.end(),
 *                                            dof_info,
 *                                            info_box,
 *                                            integrator,
 *                                            assembler);
 * 
 * @endcode
 * 
 * Right before we return the result of the error estimate, we restore the
 * old user indices.
 * 
 * @code
 *     triangulation.load_user_indices(old_user_indices);
 *     return estimates.block(0).l2_norm();
 *   }
 * 
 * @endcode
 * 
 * Here we compare our finite element solution with the (known) exact
 * solution and compute the mean quadratic error of the gradient and the
 * function itself. This function is a clone of the estimation function
 * right above.
 * 

 * 
 * Since we compute the error in the energy and the
 * <i>L<sup>2</sup></i>-norm, respectively, our block vector needs two
 * blocks here.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::error()
 *   {
 *     BlockVector<double> errors(2);
 *     errors.block(0).reinit(triangulation.n_active_cells());
 *     errors.block(1).reinit(triangulation.n_active_cells());
 * 
 *     std::vector<unsigned int> old_user_indices;
 *     triangulation.save_user_indices(old_user_indices);
 *     unsigned int i = 0;
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       cell->set_user_index(i++);
 * 
 *     MeshWorker::IntegrationInfoBox<dim> info_box;
 *     const unsigned int                  n_gauss_points =
 *       dof_handler.get_fe().tensor_degree() + 1;
 *     info_box.initialize_gauss_quadrature(n_gauss_points,
 *                                          n_gauss_points + 1,
 *                                          n_gauss_points);
 * 
 *     AnyData solution_data;
 *     solution_data.add<Vector<double> *>(&solution, "solution");
 * 
 *     info_box.cell_selector.add("solution", true, true, false);
 *     info_box.boundary_selector.add("solution", true, false, false);
 *     info_box.face_selector.add("solution", true, false, false);
 * 
 *     info_box.add_update_flags_cell(update_quadrature_points);
 *     info_box.add_update_flags_boundary(update_quadrature_points);
 *     info_box.initialize(fe, mapping, solution_data, solution);
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler);
 * 
 *     MeshWorker::Assembler::CellsAndFaces<double> assembler;
 *     AnyData                                      out_data;
 *     out_data.add<BlockVector<double> *>(&errors, "cells");
 *     assembler.initialize(out_data, false);
 * 
 *     ErrorIntegrator<dim> integrator;
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
 *                                            dof_handler.end(),
 *                                            dof_info,
 *                                            info_box,
 *                                            integrator,
 *                                            assembler);
 *     triangulation.load_user_indices(old_user_indices);
 * 
 *     deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl;
 *     deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * Create graphical output. We produce the filename by collating the
 * name from its various components, including the refinement cycle
 * that we output with two digits.
 * 
 * @code
 *   template <int dim>
 *   void
 *   InteriorPenaltyProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     const std::string filename =
 *       "sol-" + Utilities::int_to_string(cycle, 2) + ".gnuplot";
 * 
 *     deallog << "Writing solution to <" << filename << ">..." << std::endl
 *             << std::endl;
 *     std::ofstream gnuplot_output(filename);
 * 
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "u");
 *     data_out.add_data_vector(estimates.block(0), "est");
 * 
 *     data_out.build_patches();
 * 
 *     data_out.write_gnuplot(gnuplot_output);
 *   }
 * 
 * @endcode
 * 
 * And finally the adaptive loop, more or less like in previous examples.
 * 
 * @code
 *   template <int dim>
 *   void InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
 *   {
 *     deallog << "Element: " << fe.get_name() << std::endl;
 *     for (unsigned int s = 0; s < n_steps; ++s)
 *       {
 *         deallog << "Step " << s << std::endl;
 *         if (estimates.block(0).size() == 0)
 *           triangulation.refine_global(1);
 *         else
 *           {
 *             GridRefinement::refine_and_coarsen_fixed_fraction(
 *               triangulation, estimates.block(0), 0.5, 0.0);
 *             triangulation.execute_coarsening_and_refinement();
 *           }
 * 
 *         deallog << "Triangulation " << triangulation.n_active_cells()
 *                 << " cells, " << triangulation.n_levels() << " levels"
 *                 << std::endl;
 * 
 *         setup_system();
 *         deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
 *         for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
 *           deallog << ' ' << dof_handler.n_dofs(l);
 *         deallog << std::endl;
 * 
 *         deallog << "Assemble matrix" << std::endl;
 *         assemble_matrix();
 *         deallog << "Assemble multilevel matrix" << std::endl;
 *         assemble_mg_matrix();
 *         deallog << "Assemble right hand side" << std::endl;
 *         assemble_right_hand_side();
 *         deallog << "Solve" << std::endl;
 *         solve();
 *         error();
 *         deallog << "Estimate " << estimate() << std::endl;
 *         output_results(s);
 *       }
 *   }
 * } // namespace Step39
 * 
 * 
 * 
 * int main()
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step39;
 * 
 *       deallog.depth_console(2);
 *       std::ofstream logfile("deallog");
 *       deallog.attach(logfile);
 *       FE_DGQ<2>                 fe1(3);
 *       InteriorPenaltyProblem<2> test1(fe1);
 *       test1.run(12);
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


<a name="Logfileoutput"></a><h3>Logfile output</h3>

First, the program produces the usual logfile here stored in <tt>deallog</tt>. It reads (with omission of intermediate steps)

@code
DEAL::Element: FE_DGQ<2>(3)
DEAL::Step 0
DEAL::Triangulation 16 cells, 2 levels
DEAL::DoFHandler 256 dofs, level dofs 64 256
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 13 value 1.64974e-13
DEAL::energy-error: 0.297419
DEAL::L2-error:     0.00452447
DEAL::Estimate 0.990460
DEAL::Writing solution to <sol-00.gnuplot>...
DEAL::
DEAL::Step 1
DEAL::Triangulation 25 cells, 3 levels
DEAL::DoFHandler 400 dofs, level dofs 64 256 192
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 14 value 3.72262e-13
DEAL::energy-error: 0.258559
DEAL::L2-error:     0.00288510
DEAL::Estimate 0.738624
DEAL::Writing solution to <sol-01.gnuplot>...
DEAL::
DEAL::Step 2
DEAL::Triangulation 34 cells, 4 levels
DEAL::DoFHandler 544 dofs, level dofs 64 256 256 128
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 15 value 1.91610e-13
DEAL::energy-error: 0.189234
DEAL::L2-error:     0.00147954
DEAL::Estimate 0.657507
DEAL::Writing solution to <sol-02.gnuplot>...

...

DEAL::Step 10
DEAL::Triangulation 232 cells, 11 levels
DEAL::DoFHandler 3712 dofs, level dofs 64 256 896 768 768 640 512 256 256 256 256
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 51.1571
DEAL:cg::Convergence step 15 value 7.19599e-13
DEAL::energy-error: 0.0132475
DEAL::L2-error:     1.00423e-05
DEAL::Estimate 0.0470724
DEAL::Writing solution to <sol-10.gnuplot>...
DEAL::
DEAL::Step 11
DEAL::Triangulation 322 cells, 12 levels
DEAL::DoFHandler 5152 dofs, level dofs 64 256 1024 1024 896 768 768 640 448 320 320 320
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 52.2226
DEAL:cg::Convergence step 15 value 8.15195e-13
DEAL::energy-error: 0.00934891
DEAL::L2-error:     5.41095e-06
DEAL::Estimate 0.0329102
DEAL::Writing solution to <sol-11.gnuplot>...
DEAL::
@endcode

This log for instance shows that the number of conjugate gradient
iteration steps is constant at approximately 15.

<a name="Postprocessingofthelogfile"></a><h3>Postprocessing of the logfile</h3>


<img src="https://www.dealii.org/images/steps/developer/step-39-convergence.svg" alt="">
Using the perl script <tt>postprocess.pl</tt>, we extract relevant
data into <tt>output.dat</tt>, which can be used to plot graphs with
<tt>gnuplot</tt>. The graph above for instance was produced using the gnuplot
script <tt>plot_errors.gpl</tt> via

@code
perl postprocess.pl deallog &> output.dat
gnuplot plot_errors.gpl
@endcode

Reference data can be found in <tt>output.reference.dat</tt>.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-39.cc"
*/
