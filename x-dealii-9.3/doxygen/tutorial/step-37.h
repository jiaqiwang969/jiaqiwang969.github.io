/**
@page step_37 The step-37 tutorial program
This tutorial depends on step-16, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The test case</a>
        <li><a href="#Matrixvectorproductimplementation">Matrix-vector product implementation</a>
        <li><a href="#Combinationwithmultigrid">Combination with multigrid</a>
        <li><a href="#UsingCPUdependentinstructionsvectorization">Using CPU-dependent instructions (vectorization)</a>
        <li><a href="#Runningmultigridonlargescaleparallelcomputers">Running multigrid on large-scale parallel computers</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Matrixfreeimplementation">Matrix-free implementation</a>
      <ul>
        <li><a href="#Computationofcoefficient">Computation of coefficient</a>
        <li><a href="#LocalevaluationofLaplaceoperator">Local evaluation of Laplace operator</a>
      </ul>
        <li><a href="#LaplaceProblemclass">LaplaceProblem class</a>
      <ul>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_rhs">LaplaceProblem::assemble_rhs</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Comparisonwithasparsematrix">Comparison with a sparse matrix</a>
        <li><a href="#ResultsforlargescaleparallelcomputationsonSuperMUC"> Results for large-scale parallel computations on SuperMUC</a>
        <li><a href="#Adaptivity"> Adaptivity</a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions</a>
      <ul>
        <li><a href="#Kellyerrorestimator"> Kelly error estimator </a>
        <li><a href="#Sharedmemoryparallelization"> Shared-memory parallelization</a>
        <li><a href="#InhomogeneousDirichletboundaryconditions"> Inhomogeneous Dirichlet boundary conditions </a>
      <ul>
        <li><a href="#UseFEEvaluationread_dof_values_plaintoavoidresolvingconstraints"> Use FEEvaluation::read_dof_values_plain() to avoid resolving constraints </a>
        <li><a href="#UseLaplaceOperatorwithasecondAffineConstraintsobjectwithoutDirichletconditions"> Use LaplaceOperator with a second AffineConstraints object without Dirichlet conditions </a>
    </ul>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program was contributed by Katharina Kormann and Martin
Kronbichler.

The algorithm for the matrix-vector product is based on the article <a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic interface
for parallel cell-based finite element operator application</a> by Martin
Kronbichler and Katharina Kormann, Computers and Fluids 63:135&ndash;147,
2012, and the paper &quot;Parallel finite element operator application: Graph
partitioning and coloring&quot; by Katharina Kormann and Martin Kronbichler
in: Proceedings of the 7th IEEE International Conference on e-Science, 2011.

This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA). The
large-scale computations shown in the results section of this tutorial program
were supported by Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu)
by providing computing time on the GCS Supercomputer SuperMUC at Leibniz
Supercomputing Centre (LRZ, www.lrz.de) through project id pr83te. </i>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This example shows how to implement a matrix-free method, that is, a method
that does not explicitly store the matrix elements, for a second-order Poisson
equation with variable coefficients on a hypercube. The linear system will be
solved with a multigrid method and uses large-scale parallelism with MPI.

The major motivation for matrix-free methods is the fact that on today's
processors access to main memory (i.e., for objects that do not fit in the
caches) has become the bottleneck in many solvers for partial differential equations: To perform a
matrix-vector product based on matrices, modern CPUs spend far more time
waiting for data to arrive from memory than on actually doing the floating
point multiplications and additions. Thus, if we could substitute looking up
matrix elements in memory by re-computing them &mdash; or rather, the operator
represented by these entries &mdash;, we may win in terms of overall run-time
even if this requires a significant number of additional floating point
operations. That said, to realize this with a trivial implementation is not
enough and one needs to really look at the details to gain in
performance. This tutorial program and the papers referenced above show how
one can implement such a scheme and demonstrates the speedup that can be
obtained.


<a name="Thetestcase"></a><h3>The test case</h3>


In this example, we consider the Poisson problem @f{eqnarray*} -
\nabla \cdot a(\mathbf x) \nabla u &=& 1, \\ u &=& 0 \quad \text{on}\
\partial \Omega @f} where $a(\mathbf x)$ is a variable coefficient.
Below, we explain how to implement a matrix-vector product for this
problem without explicitly forming the matrix. The construction can,
of course, be done in a similar way for other equations as well.

We choose as domain $\Omega=[0,1]^3$ and $a(\mathbf x)=\frac{1}{0.05 +
2\|\mathbf x\|^2}$. Since the coefficient is symmetric around the
origin but the domain is not, we will end up with a non-symmetric
solution.


<a name="Matrixvectorproductimplementation"></a><h3>Matrix-vector product implementation</h3>


In order to find out how we can write a code that performs a matrix-vector
product, but does not need to store the matrix elements, let us start at
looking how a finite element matrix <i>A</i> is assembled:
@f{eqnarray*}
A = \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}}
P_{\mathrm{cell,{loc-glob}}}^T A_{\mathrm{cell}} P_{\mathrm{cell,{loc-glob}}}.
@f}
In this formula, the matrix <i>P</i><sub>cell,loc-glob</sub> is a rectangular
matrix that defines the index mapping from local degrees of freedom in the
current cell to the global degrees of freedom. The information from which this
operator can be built is usually encoded in the <code>local_dof_indices</code>
variable and is used in the assembly calls filling matrices in deal.II. Here,
<i>A</i><sub>cell</sub> denotes the cell matrix associated with <i>A</i>.

If we are to perform a matrix-vector product, we can hence use that
@f{eqnarray*}
y &=& A\cdot u = \left(\sum_{\text{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
A_\mathrm{cell} P_\mathrm{cell,{loc-glob}}\right) \cdot u
\\
&=& \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
A_\mathrm{cell} u_\mathrm{cell}
\\
&=& \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
v_\mathrm{cell},
@f}
where <i>u</i><sub>cell</sub> are the values of <i>u</i> at the degrees of freedom
of the respective cell, and
<i>v</i><sub>cell</sub>=<i>A</i><sub>cell</sub><i>u</i><sub>cell</sub>
correspondingly for the result.
A naive attempt to implement the local action of the Laplacian would hence be
to use the following code:
@code
Matrixfree<dim>::vmult (Vector<double>       &dst,
                        const Vector<double> &src) const
{
  dst = 0;

  QGauss<dim>  quadrature_formula(fe.degree+1);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_gradients | update_JxW_values|
                           update_quadrature_points);

  const unsigned int   dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_src (dofs_per_cell),
                       cell_dst (dofs_per_cell);
  const Coefficient<dim> coefficient;
  std::vector<double> coefficient_values(n_q_points);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  for (const auto & cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      fe_values.reinit (cell);
      coefficient.value_list(fe_values.get_quadrature_points(),
                             coefficient_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (fe_values.shape_grad(i,q) *
                                 fe_values.shape_grad(j,q) *
                                 fe_values.JxW(q)*
                                 coefficient_values[q]);

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        cell_src(i) = src(local_dof_indices(i));

      cell_matrix.vmult (cell_dst, cell_src);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        dst(local_dof_indices(i)) += cell_dst;
    }
}
@endcode

Here we neglected boundary conditions as well as any hanging nodes we may
have, though neither would be very difficult to include using the
AffineConstraints class. Note how we first generate the local matrix in the
usual way as a sum over all quadrature points for each local matrix entry.
To form the actual product as expressed in the above formula, we
extract the values of <code>src</code> of the cell-related degrees of freedom
(the action of <i>P</i><sub>cell,loc-glob</sub>), multiply by the local matrix
(the action of <i>A</i><sub>cell</sub>), and finally add the result to the
destination vector <code>dst</code> (the action of
<i>P</i><sub>cell,loc-glob</sub><sup>T</sup>, added over all the elements). It
is not more difficult than that, in principle.

While this code is completely correct, it is very slow. For every cell, we
generate a local matrix, which takes three nested loops with loop length equal
to the number of local degrees of freedom to compute. The
multiplication itself is then done by two nested loops, which means that it
is much cheaper.

One way to improve this is to realize that conceptually the local
matrix can be thought of as the product of three matrices,
@f{eqnarray*}
A_\mathrm{cell} = B_\mathrm{cell}^T D_\mathrm{cell} B_\mathrm{cell},
@f}
where for the example of the Laplace operator the (<i>q</i>*dim+<i>d,i</i>)-th
element of <i>B</i><sub>cell</sub> is given by
<code>fe_values.shape_grad(i,q)[d]</code>. This matrix consists of
<code>dim*n_q_points</code> rows and @p dofs_per_cell columns. The matrix
<i>D</i><sub>cell</sub> is diagonal and contains the values
<code>fe_values.JxW(q) * coefficient_values[q]</code> (or, rather, @p
dim copies of each of these values). This kind of representation of
finite element matrices can often be found in the engineering literature.

When the cell matrix is applied to a vector,
@f{eqnarray*}
A_\mathrm{cell}\cdot u_\mathrm{cell} = B_\mathrm{cell}^T
D_\mathrm{cell} B_\mathrm{cell} \cdot u_\mathrm{cell},
@f}
one would then not form the matrix-matrix products, but rather multiply one
matrix at a time with a vector from right to left so that only three
successive matrix-vector products are formed. This approach removes the three
nested loops in the calculation of the local matrix, which reduces the
complexity of the work on one cell from something like $\mathcal
{O}(\mathrm{dofs\_per\_cell}^3)$ to $\mathcal
{O}(\mathrm{dofs\_per\_cell}^2)$. An interpretation of this algorithm is that
we first transform the vector of values on the local DoFs to a vector of
gradients on the quadrature points. In the second loop, we multiply these
gradients by the integration weight and the coefficient. The third loop applies
the second gradient (in transposed form), so that we get back to a vector of
(Laplacian) values on the cell dofs.

The bottleneck in the above code is the operations done by the call to
FEValues::reinit for every <code>cell</code>, which take about as much time as
the other steps together (at least if the mesh is unstructured; deal.II can
recognize that the gradients are often unchanged on structured meshes). That
is certainly not ideal and we would like to do better than this. What the
reinit function does is to calculate the gradient in real space by
transforming the gradient on the reference cell using the Jacobian of the
transformation from real to reference cell. This is done for each basis
function on the cell, for each quadrature point. The Jacobian does not depend
on the basis function, but it is different on different quadrature points in
general. If you only build the matrix once as we've done in all previous
tutorial programs, there is nothing to be optimized since FEValues::reinit
needs to be called on every cell. In this process, the transformation is
applied while computing the local matrix elements.

In a matrix-free implementation, however, we will compute those integrals very
often because iterative solvers will apply the matrix many times during the
solution process. Therefore, we need to think about whether we may be able to
cache some data that gets reused in the operator applications, i.e., integral
computations. On the other hand, we realize that we must not cache too much
data since otherwise we get back to the situation where memory access becomes
the dominating factor. Therefore, we will not store the transformed gradients
in the matrix <i>B</i>, as they would in general be different for each basis
function and each quadrature point on every element for curved meshes.

The trick is to factor out the Jacobian transformation and first apply the
gradient on the reference cell only. This operation interpolates the vector of
values on the local dofs to a vector of (unit-coordinate) gradients on the
quadrature points. There, we first apply the Jacobian that we factored out
from the gradient, then apply the weights of the quadrature, and finally apply
the transposed Jacobian for preparing the third loop which tests by the
gradients on the unit cell and sums over quadrature points.

Let us again write this in terms of matrices. Let the matrix
<i>B</i><sub>cell</sub> denote the cell-related gradient matrix, with each row
containing the values on the quadrature points. It is constructed by a
matrix-matrix product as @f{eqnarray*} B_\mathrm{cell} =
J_\mathrm{cell}^{-\mathrm T} B_\mathrm{ref\_cell}, @f} where
<i>B</i><sub>ref_cell</sub> denotes the gradient on the reference cell and
<i>J</i><sup>-T</sup><sub>cell</sub> denotes the inverse transpose Jacobian of
the transformation from unit to real cell (in the language of transformations,
the operation represented by <i>J</i><sup>-T</sup><sub>cell</sub> represents a
covariant transformation). <i>J</i><sup>-T</sup><sub>cell</sub> is
block-diagonal, and the blocks size is equal to the dimension of the
problem. Each diagonal block is the Jacobian transformation that goes from the
reference cell to the real cell.

Putting things together, we find that
@f{eqnarray*}
A_\mathrm{cell} = B_\mathrm{cell}^T D B_\mathrm{cell}
                = B_\mathrm{ref\_cell}^T J_\mathrm{cell}^{-1}
                  D_\mathrm{cell}
                  J_\mathrm{cell}^{-\mathrm T} B_\mathrm{ref\_cell},
@f}
so we calculate the product (starting the local product from the right)
@f{eqnarray*}
v_\mathrm{cell} = B_\mathrm{ref\_cell}^T J_\mathrm{cell}^{-1} D J_\mathrm{cell}^{-\mathrm T}
B_\mathrm{ref\_cell} u_\mathrm{cell}, \quad
v = \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
v_\mathrm{cell}.
@f}
@code
  FEValues<dim> fe_values_reference (fe, quadrature_formula,
                                     update_gradients);
  Triangulation<dim> reference_cell;
  GridGenerator::hyper_cube(reference_cell, 0., 1.);
  fe_values_reference.reinit (reference_cell.begin());

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_inverse_jacobians | update_JxW_values |
                           update_quadrature_points);

  for (const auto & cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit (cell);
      coefficient.value_list(fe_values.get_quadrature_points(),
                             coefficient_values);

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        cell_src(i) = src(local_dof_indices(i));

      temp_vector = 0;
      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            temp_vector(q*dim+d) +=
              fe_values_reference.shape_grad(i,q)[d] * cell_src(i);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // apply the transposed inverse Jacobian of the mapping
          Tensor<1,dim> temp;
          for (unsigned int d=0; d<dim; ++d)
            temp[d] = temp_vector(q*dim+d);
          for (unsigned int d=0; d<dim; ++d)
            {
              double sum = 0;
              for (unsigned int e=0; e<dim; ++e)
                sum += fe_values.inverse_jacobian(q)[e][d] *
                               temp[e];
              temp_vector(q*dim+d) = sum;
            }

          // multiply by coefficient and integration weight
          for (unsigned int d=0; d<dim; ++d)
            temp_vector(q*dim+d) *= fe_values.JxW(q) * coefficient_values[q];

          // apply the inverse Jacobian of the mapping
          for (unsigned int d=0; d<dim; ++d)
            temp[d] = temp_vector(q*dim+d);
          for (unsigned int d=0; d<dim; ++d)
            {
              double sum = 0;
              for (unsigned int e=0; e<dim; ++e)
                sum += fe_values.inverse_jacobian(q)[d][e] *
                       temp[e];
              temp_vector(q*dim+d) = sum;
            }
        }

      cell_dst = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int d=0; d<dim; ++d)
            cell_dst(i) += fe_values_reference.shape_grad(i,q)[d] *
                                   temp_vector(q*dim+d);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        dst(local_dof_indices(i)) += cell_dst(i);
    }
}
@endcode

Note how we create an additional FEValues object for the reference cell
gradients and how we initialize it to the reference cell. The actual
derivative data is then applied by the inverse, transposed Jacobians (deal.II
calls the Jacobian matrix from real to unit cell inverse_jacobian, as the
forward transformation is from unit to real cell). The factor
$J_\mathrm{cell}^{-1} D_\mathrm{cell} J_\mathrm{cell}^{-\mathrm T}$ is
block-diagonal over quadrature. In this form, one realizes that variable
coefficients (possibly expressed through a tensor) and general grid topologies
with Jacobian transformations have a similar effect on the coefficient
transforming the unit-cell derivatives.

At this point, one might wonder why we store the matrix
$J_\mathrm{cell}^{-\mathrm T}$ and the coefficient separately, rather than
only the complete factor $J_\mathrm{cell}^{-1} D_\mathrm{cell}
J_\mathrm{cell}^{-\mathrm T}$. The latter would use less memory because the
tensor is symmetric with six independent values in 3D, whereas for the former
we would need nine entries for the inverse transposed Jacobian, one for the
quadrature weight and Jacobian determinant, and one for the coefficient,
totaling to 11 doubles. The reason is that the former approach allows for
implementing generic differential operators through a common framework of
cached data, whereas the latter specifically stores the coefficient for the
Laplacian. In case applications demand for it, this specialization could pay
off and would be worthwhile to consider. Note that the implementation in
deal.II is smart enough to detect Cartesian or affine geometries where the
Jacobian is constant throughout the cell and needs not be stored for every
cell (and indeed often is the same over different cells as well).

The final optimization that is most crucial from an operation count point of
view is to make use of the tensor product structure in the basis
functions. This is possible because we have factored out the gradient from the
reference cell operation described by <i>B</i><sub>ref_cell</sub>, i.e., an
interpolation operation over the completely regular data fields of the
reference cell. We illustrate the process of complexity reduction in two space
dimensions, but the same technique can be used in higher dimensions. On the
reference cell, the basis functions are of the tensor product form
$\phi(x,y,z) = \varphi_i(x) \varphi_j(y)$. The part of the matrix
<i>B</i><sub>ref_cell</sub> that computes the first component has the form
$B_\mathrm{sub\_cell}^x = B_\mathrm{grad,x} \otimes B_\mathrm{val,y}$, where
<i>B</i><sub>grad,x</sub> and <i>B</i><sub>val,y</sub> contain the evaluation
of all the 1D basis functions on all the 1D quadrature points. Forming a
matrix <i>U</i> with <i>U(j,i)</i> containing the coefficient belonging to
basis function $\varphi_i(x) \varphi_j(y)$, we get $(B_\mathrm{grad,x} \otimes
B_\mathrm{val,y})u_\mathrm{cell} = B_\mathrm{val,y} U B_\mathrm{grad,x}$. This
reduces the complexity for computing this product from $p^4$ to $2 p^3$, where
<i>p</i>-1 is the degree of the finite element (i.e., equivalently, <i>p</i>
is the number of shape functions in each coordinate direction), or $p^{2d}$ to
$d p^{d+1}$ in general. The reason why we look at the complexity in terms of
the polynomial degree is since we want to be able to go to high degrees and
possibly increase the polynomial degree <i>p</i> instead of the grid
resolution. Good algorithms for moderate degrees like the ones used here are
linear in the polynomial degree independent on the dimension, as opposed to
matrix-based schemes or naive evaluation through FEValues. The techniques used
in the implementations of deal.II have been established in the spectral
element community since the 1980s.

Implementing a matrix-free and cell-based finite element operator requires a
somewhat different program design as compared to the usual matrix assembly
codes shown in previous tutorial programs. The data structures for doing this
are the MatrixFree class that collects all data and issues a (parallel) loop
over all cells and the FEEvaluation class that evaluates finite element basis
functions by making use of the tensor product structure.

The implementation of the matrix-free matrix-vector product shown in this
tutorial is slower than a matrix-vector product using a sparse matrix for
linear elements, but faster for all higher order elements thanks to the
reduced complexity due to the tensor product structure and due to less memory
transfer during computations. The impact of reduced memory transfer is
particularly beneficial when working on a multicore processor where several
processing units share access to memory. In that case, an algorithm which is
computation bound will show almost perfect parallel speedup (apart from
possible changes of the processor's clock frequency through turbo modes
depending on how many cores are active), whereas an algorithm that is bound by
memory transfer might not achieve similar speedup (even when the work is
perfectly parallel and one could expect perfect scaling like in sparse
matrix-vector products). An additional gain with this implementation is that
we do not have to build the sparse matrix itself, which can also be quite
expensive depending on the underlying differential equation. Moreover, the
above framework is simple to generalize to nonlinear operations, as we
demonstrate in step-48.


<a name="Combinationwithmultigrid"></a><h3>Combination with multigrid</h3>


Above, we have gone to significant lengths to implement a matrix-vector
product that does not actually store the matrix elements. In many user codes,
however, one wants more than just doing a few matrix-vector products &mdash;
one wants to do as few of these operations as possible when solving linear
systems. In theory, we could use the CG method without preconditioning;
however, that would not be very efficient for the Laplacian. Rather,
preconditioners are used for increasing the speed of
convergence. Unfortunately, most of the more frequently used preconditioners
such as SSOR, ILU or algebraic multigrid (AMG) cannot be used here because
their implementation requires knowledge of the elements of the system matrix.

One solution is to use geometric multigrid methods as shown in step-16. They
are known to be very fast, and they are suitable for our purpose since all
ingredients, including the transfer between different grid levels, can be
expressed in terms of matrix-vector products related to a collection of
cells. All one needs to do is to find a smoother that is based on
matrix-vector products rather than all the matrix entries. One such candidate
would be a damped Jacobi iteration that requires access to the matrix
diagonal, but it is often not sufficiently good in damping all high-frequency
errors. The properties of the Jacobi method can be improved by iterating it a
few times with the so-called Chebyshev iteration. The Chebyshev iteration is
described by a polynomial expression of the matrix-vector product where the
coefficients can be chosen to achieve certain properties, in this case to
smooth the high-frequency components of the error which are associated to the
eigenvalues of the Jacobi-preconditioned matrix. At degree zero, the Jacobi
method with optimal damping parameter is retrieved, whereas higher order
corrections are used to improve the smoothing properties. The effectiveness of
Chebyshev smoothing in multigrid has been demonstrated, e.g., in the article
<a href="http://www.sciencedirect.com/science/article/pii/S0021999103001943">
<i>M. Adams, M. Brezina, J. Hu, R. Tuminaro. Parallel multigrid smoothers:
polynomial versus Gauss&ndash;Seidel, J. Comput. Phys. 188:593&ndash;610,
2003</i></a>. This publication also identifies one more advantage of
Chebyshev smoothers that we exploit here, namely that they are easy to
parallelize, whereas SOR/Gauss&ndash;Seidel smoothing relies on substitutions,
for which a naive parallelization works on diagonal sub-blocks of the matrix,
thereby decreases efficiency (for more detail see e.g. Y. Saad, Iterative
Methods for Sparse Linear Systems, SIAM, 2nd edition, 2003, chapters 11 & 12).

The implementation into the multigrid framework is then straightforward. The
multigrid implementation in this program is similar to step-16 and includes
adaptivity.


<a name="UsingCPUdependentinstructionsvectorization"></a><h3>Using CPU-dependent instructions (vectorization)</h3>


The computational kernels for evaluation in FEEvaluation are written in a way
to optimally use computational resources. To achieve this, they do not operate
on double data types, but something we call VectorizedArray (check e.g. the
return type of FEEvaluationBase::get_value, which is VectorizedArray for a
scalar element and a Tensor of VectorizedArray for a vector finite
element). VectorizedArray is a short array of doubles or float whose length
depends on the particular computer system in use. For example, systems based
on x86-64 support the streaming SIMD extensions (SSE), where the processor's
vector units can process two doubles (or four single-precision floats) by one
CPU instruction. Newer processors (from about 2012 and onwards) support the
so-called advanced vector extensions (AVX) with 256 bit operands, which can
use four doubles and eight floats, respectively. Vectorization is a
single-instruction/multiple-data (SIMD) concept, that is, one CPU instruction
is used to process multiple data values at once. Often, finite element
programs do not use vectorization explicitly as the benefits of this concept
are only in arithmetic intensive operations. The bulk of typical finite
element workloads are memory bandwidth limited (operations on sparse matrices
and vectors) where the additional computational power is useless.

Behind the scenes, optimized BLAS packages might heavily rely on
vectorization, though. Also, optimizing compilers might automatically
transform loops involving standard code into more efficient vectorized form
(deal.II uses OpenMP SIMD pragmas inside the regular loops of vector
updates). However, the data flow must be very regular in order for compilers
to produce efficient code. For example, already the automatic vectorization of
the prototype operation that benefits from vectorization, matrix-matrix
products, fails on most compilers (as of writing this tutorial in early 2012
and updating in late 2016, neither gcc nor the Intel compiler manage to
produce useful vectorized code for the FullMatrix::mmult function, and not
even on the simpler case where the matrix bounds are compile-time constants
instead of run-time constants as in FullMatrix::mmult). The main reason for
this is that the information to be processed at the innermost loop (that is
where vectorization is applied) is not necessarily a multiple of the vector
length, leaving parts of the resources unused. Moreover, the data that can
potentially be processed together might not be laid out in a contiguous way in
memory or not with the necessary alignment to address boundaries that are
needed by the processor. Or the compiler might not be able to prove that data
arrays do not overlap when loading several elements at once.

In the matrix-free implementation in deal.II, we have therefore chosen to
apply vectorization at the level which is most appropriate for finite element
computations: The cell-wise computations are typically exactly the same for
all cells (except for indices in the indirect addressing used while reading
from and writing to vectors), and hence SIMD can be used to process several
cells at once. In all what follows, you can think of a VectorizedArray to hold
data from several cells. Remember that it is not related to the spatial
dimension and the number of elements e.g. in a Tensor or Point.

Note that vectorization depends on the CPU the code is running on and for
which the code is compiled. In order to generate the fastest kernels of
FEEvaluation for your computer, you should compile deal.II with the so-called
<i>native</i> processor variant. When using the gcc compiler, it can be
enabled by setting the variable <tt>CMAKE_CXX_FLAGS</tt> to
<tt>"-march=native"</tt> in the cmake build settings (on the command line,
specify <tt>-DCMAKE_CXX_FLAGS="-march=native"</tt>, see the deal.II README for
more information). Similar options exist for other compilers. We output
the current vectorization length in the run() function of this example.


<a name="Runningmultigridonlargescaleparallelcomputers"></a><h3>Running multigrid on large-scale parallel computers</h3>


As mentioned above, all components in the matrix-free framework can easily be
parallelized with MPI using domain decomposition. Thanks to the easy access to
large-scale parallel meshes through p4est (see step-40 for details) in
deal.II, and the fact that cell-based loops with matrix-free evaluation
<i>only</i> need a decomposition of the mesh into chunks of roughly equal size
on each processor, there is relatively little to do to write a parallel
program working with distributed memory. While other tutorial programs using
MPI have relied on either PETSc or Trilinos, this program uses deal.II's own
parallel vector facilities.

The deal.II parallel vector class, LinearAlgebra::distributed::Vector, holds
the processor-local part of the solution as well as data fields for ghosted
DoFs, i.e. DoFs that are owned by a remote processor but accessed by cells
that are owned by the present processor. In the @ref GlossLocallyActiveDof
"glossary" these degrees of freedom are referred to as locally active degrees
of freedom. The function MatrixFree::initialize_dof_vector() provides a method
that sets this design. Note that hanging nodes can relate to additional
ghosted degrees of freedom that must be included in the distributed vector but
are not part of the locally active DoFs in the sense of the @ref
GlossLocallyActiveDof "glossary". Moreover, the distributed vector holds the
MPI metadata for DoFs that are owned locally but needed by other
processors. A benefit of the design of this vector class is the way ghosted
entries are accessed. In the storage scheme of the vector, the data array
extends beyond the processor-local part of the solution with further vector
entries available for the ghosted degrees of freedom. This gives a contiguous
index range for all locally active degrees of freedom. (Note that the index
range depends on the exact configuration of the mesh.) Since matrix-free
operations can be thought of doing linear algebra that is performance
critical, and performance-critical code cannot waste time on doing MPI-global
to MPI-local index translations, the availability of an index spaces local to
one MPI rank is fundamental. The way things are accessed here is a direct
array access. This is provided through
LinearAlgebra::distributed::Vector::local_element(), but it is actually rarely
needed because all of this happens internally in FEEvaluation.

The design of LinearAlgebra::distributed::Vector is similar to the
PETScWrappers::MPI::Vector and TrilinosWrappers::MPI::Vector data types we
have used in step-40 and step-32 before, but since we do not need any other
parallel functionality of these libraries, we use the
LinearAlgebra::distributed::Vector class of deal.II instead of linking in another
large library in this tutorial program. Also note that the PETSc and Trilinos
vectors do not provide the fine-grained control over ghost entries with direct
array access because they abstract away the necessary implementation details.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * First include the necessary files from the deal.II library.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer_matrix_free.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * This includes the data structures for the efficient implementation of
 * matrix-free methods or more generic finite element operators with the class
 * MatrixFree.
 * 
 * @code
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/operators.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * 
 * namespace Step37
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * To be efficient, the operations performed in the matrix-free
 * implementation require knowledge of loop lengths at compile time, which
 * are given by the degree of the finite element. Hence, we collect the
 * values of the two template parameters that can be changed at one place in
 * the code. Of course, one could make the degree of the finite element a
 * run-time parameter by compiling the computational kernels for all degrees
 * that are likely (say, between 1 and 6) and selecting the appropriate
 * kernel at run time. Here, we simply choose second order $Q_2$ elements
 * and choose dimension 3 as standard.
 * 
 * @code
 *   const unsigned int degree_finite_element = 2;
 *   const unsigned int dimension             = 3;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * We define a variable coefficient function for the Poisson problem. It is
 * similar to the function in step-5 but we use the form $a(\mathbf
 * x)=\frac{1}{0.05 + 2\|\bf x\|^2}$ instead of a discontinuous one. It is
 * merely to demonstrate the possibilities of this implementation, rather
 * than making much sense physically. We define the coefficient in the same
 * way as functions in earlier tutorial programs. There is one new function,
 * namely a @p value method with template argument @p number.
 * 
 * @code
 *   template <int dim>
 *   class Coefficient : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     template <typename number>
 *     number value(const Point<dim, number> &p,
 *                  const unsigned int        component = 0) const;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * This is the new function mentioned above: Evaluate the coefficient for
 * abstract type @p number. It might be just a usual double, but it can also
 * be a somewhat more complicated type that we call VectorizedArray. This
 * data type is essentially a short array of doubles as discussed in the
 * introduction that holds data from several cells. For example, we evaluate
 * the coefficient shown here not on a simple point as usually done, but we
 * hand it a Point<dim,VectorizedArray<double> > point, which is actually a
 * collection of four points in the case of AVX. Do not confuse the entries
 * in VectorizedArray with the different coordinates of the point. Indeed,
 * the data is laid out such that <code>p[0]</code> returns a
 * VectorizedArray, which in turn contains the x-coordinate for the first
 * point and the second point. You may access the coordinates individually
 * using e.g. <code>p[0][j]</code>, j=0,1,2,3, but it is recommended to
 * define operations on a VectorizedArray as much as possible in order to
 * make use of vectorized operations.
 *   

 * 
 * In the function implementation, we assume that the number type overloads
 * basic arithmetic operations, so we just write the code as usual. The base
 * class function @p value is then computed from the templated function with
 * double type, in order to avoid duplicating code.
 * 
 * @code
 *   template <int dim>
 *   template <typename number>
 *   number Coefficient<dim>::value(const Point<dim, number> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     return 1. / (0.05 + 2. * p.square());
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   double Coefficient<dim>::value(const Point<dim> & p,
 *                                  const unsigned int component) const
 *   {
 *     return value<double>(p, component);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Matrixfreeimplementation"></a> 
 * <h3>Matrix-free implementation</h3>
 * 

 * 
 * The following class, called <code>LaplaceOperator</code>, implements the
 * differential operator. For all practical purposes, it is a matrix, i.e.,
 * you can ask it for its size (member functions <code>m(), n()</code>) and
 * you can apply it to a vector (the <code>vmult()</code> function). The
 * difference to a real matrix of course lies in the fact that this class
 * does not actually store the <i>elements</i> of the matrix, but only knows
 * how to compute the action of the operator when applied to a vector.
 *   

 * 
 * The infrastructure describing the matrix size, the initialization from a
 * MatrixFree object, and the various interfaces to matrix-vector products
 * through vmult() and Tvmult() methods, is provided by the class
 * MatrixFreeOperator::Base from which this class derives. The
 * LaplaceOperator class defined here only has to provide a few interfaces,
 * namely the actual action of the operator through the apply_add() method
 * that gets used in the vmult() functions, and a method to compute the
 * diagonal entries of the underlying matrix. We need the diagonal for the
 * definition of the multigrid smoother. Since we consider a problem with
 * variable coefficient, we further implement a method that can fill the
 * coefficient values.
 *   

 * 
 * Note that the file <code>include/deal.II/matrix_free/operators.h</code>
 * already contains an implementation of the Laplacian through the class
 * MatrixFreeOperators::LaplaceOperator. For educational purposes, the
 * operator is re-implemented in this tutorial program, explaining the
 * ingredients and concepts used there.
 *   

 * 
 * This program makes use of the data cache for finite element operator
 * application that is integrated in deal.II. This data cache class is
 * called MatrixFree. It contains mapping information (Jacobians) and index
 * relations between local and global degrees of freedom. It also contains
 * constraints like the ones from hanging nodes or Dirichlet boundary
 * conditions. Moreover, it can issue a loop over all cells in %parallel,
 * making sure that only cells are worked on that do not share any degree of
 * freedom (this makes the loop thread-safe when writing into destination
 * vectors). This is a more advanced strategy compared to the WorkStream
 * class described in the @ref threads module. Of course, to not destroy
 * thread-safety, we have to be careful when writing into class-global
 * structures.
 *   

 * 
 * The class implementing the Laplace operator has three template arguments,
 * one for the dimension (as many deal.II classes carry), one for the degree
 * of the finite element (which we need to enable efficient computations
 * through the FEEvaluation class), and one for the underlying scalar
 * type. We want to use <code>double</code> numbers (i.e., double precision,
 * 64-bit floating point) for the final matrix, but floats (single
 * precision, 32-bit floating point numbers) for the multigrid level
 * matrices (as that is only a preconditioner, and floats can be processed
 * twice as fast). The class FEEvaluation also takes a template argument for
 * the number of quadrature points in one dimension. In the code below, we
 * hard-code it to <code>fe_degree+1</code>. If we wanted to change it
 * independently of the polynomial degree, we would need to add a template
 * parameter as is done in the MatrixFreeOperators::LaplaceOperator class.
 *   

 * 
 * As a sidenote, if we implemented several different operations on the same
 * grid and degrees of freedom (like a mass matrix and a Laplace matrix), we
 * would define two classes like the current one for each of the operators
 * (derived from the MatrixFreeOperators::Base class), and let both of them
 * refer to the same MatrixFree data cache from the general problem
 * class. The interface through MatrixFreeOperators::Base requires us to
 * only provide a minimal set of functions. This concept allows for writing
 * complex application codes with many matrix-free operations.
 *   

 * 
 * @note Storing values of type <code>VectorizedArray<number></code>
 * requires care: Here, we use the deal.II table class which is prepared to
 * hold the data with correct alignment. However, storing e.g. an
 * <code>std::vector<VectorizedArray<number> ></code> is not possible with
 * vectorization: A certain alignment of the data with the memory address
 * boundaries is required (essentially, a VectorizedArray that is 32 bytes
 * long in case of AVX needs to start at a memory address that is divisible
 * by 32). The table class (as well as the AlignedVector class it is based
 * on) makes sure that this alignment is respected, whereas std::vector does
 * not in general, which may lead to segmentation faults at strange places
 * for some systems or suboptimal performance for other systems.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   class LaplaceOperator
 *     : public MatrixFreeOperators::
 *         Base<dim, LinearAlgebra::distributed::Vector<number>>
 *   {
 *   public:
 *     using value_type = number;
 * 
 *     LaplaceOperator();
 * 
 *     void clear() override;
 * 
 *     void evaluate_coefficient(const Coefficient<dim> &coefficient_function);
 * 
 *     virtual void compute_diagonal() override;
 * 
 *   private:
 *     virtual void apply_add(
 *       LinearAlgebra::distributed::Vector<number> &      dst,
 *       const LinearAlgebra::distributed::Vector<number> &src) const override;
 * 
 *     void
 *     local_apply(const MatrixFree<dim, number> &                   data,
 *                 LinearAlgebra::distributed::Vector<number> &      dst,
 *                 const LinearAlgebra::distributed::Vector<number> &src,
 *                 const std::pair<unsigned int, unsigned int> &cell_range) const;
 * 
 *     void local_compute_diagonal(
 *       const MatrixFree<dim, number> &              data,
 *       LinearAlgebra::distributed::Vector<number> & dst,
 *       const unsigned int &                         dummy,
 *       const std::pair<unsigned int, unsigned int> &cell_range) const;
 * 
 *     Table<2, VectorizedArray<number>> coefficient;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * This is the constructor of the @p LaplaceOperator class. All it does is
 * to call the default constructor of the base class
 * MatrixFreeOperators::Base, which in turn is based on the Subscriptor
 * class that asserts that this class is not accessed after going out of scope
 * e.g. in a preconditioner.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
 *     : MatrixFreeOperators::Base<dim,
 *                                 LinearAlgebra::distributed::Vector<number>>()
 *   {}
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::clear()
 *   {
 *     coefficient.reinit(0, 0);
 *     MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
 *       clear();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computationofcoefficient"></a> 
 * <h4>Computation of coefficient</h4>
 * 

 * 
 * To initialize the coefficient, we directly give it the Coefficient class
 * defined above and then select the method
 * <code>coefficient_function.value</code> with vectorized number (which the
 * compiler can deduce from the point data type). The use of the
 * FEEvaluation class (and its template arguments) will be explained below.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::evaluate_coefficient(
 *     const Coefficient<dim> &coefficient_function)
 *   {
 *     const unsigned int n_cells = this->data->n_cell_batches();
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*this->data);
 * 
 *     coefficient.reinit(n_cells, phi.n_q_points);
 *     for (unsigned int cell = 0; cell < n_cells; ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           coefficient(cell, q) =
 *             coefficient_function.value(phi.quadrature_point(q));
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LocalevaluationofLaplaceoperator"></a> 
 * <h4>Local evaluation of Laplace operator</h4>
 * 

 * 
 * Here comes the main function of this class, the evaluation of the
 * matrix-vector product (or, in general, a finite element operator
 * evaluation). This is done in a function that takes exactly four
 * arguments, the MatrixFree object, the destination and source vectors, and
 * a range of cells that are to be worked on. The method
 * <code>cell_loop</code> in the MatrixFree class will internally call this
 * function with some range of cells that is obtained by checking which
 * cells are possible to work on simultaneously so that write operations do
 * not cause any race condition. Note that the cell range used in the loop
 * is not directly the number of (active) cells in the current mesh, but
 * rather a collection of batches of cells.  In other word, "cell" may be
 * the wrong term to begin with, since FEEvaluation groups data from several
 * cells together. This means that in the loop over quadrature points we are
 * actually seeing a group of quadrature points of several cells as one
 * block. This is done to enable a higher degree of vectorization.  The
 * number of such "cells" or "cell batches" is stored in MatrixFree and can
 * be queried through MatrixFree::n_cell_batches(). Compared to the deal.II
 * cell iterators, in this class all cells are laid out in a plain array
 * with no direct knowledge of level or neighborship relations, which makes
 * it possible to index the cells by unsigned integers.
 *   

 * 
 * The implementation of the Laplace operator is quite simple: First, we
 * need to create an object FEEvaluation that contains the computational
 * kernels and has data fields to store temporary results (e.g. gradients
 * evaluated on all quadrature points on a collection of a few cells). Note
 * that temporary results do not use a lot of memory, and since we specify
 * template arguments with the element order, the data is stored on the
 * stack (without expensive memory allocation). Usually, one only needs to
 * set two template arguments, the dimension as a first argument and the
 * degree of the finite element as the second argument (this is equal to the
 * number of degrees of freedom per dimension minus one for FE_Q
 * elements). However, here we also want to be able to use float numbers for
 * the multigrid preconditioner, which is the last (fifth) template
 * argument. Therefore, we cannot rely on the default template arguments and
 * must also fill the third and fourth field, consequently. The third
 * argument specifies the number of quadrature points per direction and has
 * a default value equal to the degree of the element plus one. The fourth
 * argument sets the number of components (one can also evaluate
 * vector-valued functions in systems of PDEs, but the default is a scalar
 * element), and finally the last argument sets the number type.
 *   

 * 
 * Next, we loop over the given cell range and then we continue with the
 * actual implementation: <ol> <li>Tell the FEEvaluation object the (macro)
 * cell we want to work on.  <li>Read in the values of the source vectors
 * (@p read_dof_values), including the resolution of constraints. This
 * stores $u_\mathrm{cell}$ as described in the introduction.  <li>Compute
 * the unit-cell gradient (the evaluation of finite element
 * functions). Since FEEvaluation can combine value computations with
 * gradient computations, it uses a unified interface to all kinds of
 * derivatives of order between zero and two. We only want gradients, no
 * values and no second derivatives, so we set the function arguments to
 * true in the gradient slot (second slot), and to false in the values slot
 * (first slot). There is also a third slot for the Hessian which is
 * false by default, so it needs not be given. Note that the FEEvaluation
 * class internally evaluates shape functions in an efficient way where one
 * dimension is worked on at a time (using the tensor product form of shape
 * functions and quadrature points as mentioned in the introduction). This
 * gives complexity equal to $\mathcal O(d^2 (p+1)^{d+1})$ for polynomial
 * degree $p$ in $d$ dimensions, compared to the naive approach with loops
 * over all local degrees of freedom and quadrature points that is used in
 * FEValues and costs $\mathcal O(d (p+1)^{2d})$.  <li>Next comes the
 * application of the Jacobian transformation, the multiplication by the
 * variable coefficient and the quadrature weight. FEEvaluation has an
 * access function @p get_gradient that applies the Jacobian and returns the
 * gradient in real space. Then, we just need to multiply by the (scalar)
 * coefficient, and let the function @p submit_gradient apply the second
 * Jacobian (for the test function) and the quadrature weight and Jacobian
 * determinant (JxW). Note that the submitted gradient is stored in the same
 * data field as where it is read from in @p get_gradient. Therefore, you
 * need to make sure to not read from the same quadrature point again after
 * having called @p submit_gradient on that particular quadrature point. In
 * general, it is a good idea to copy the result of @p get_gradient when it
 * is used more often than once.  <li>Next follows the summation over
 * quadrature points for all test functions that corresponds to the actual
 * integration step. For the Laplace operator, we just multiply by the
 * gradient, so we call the integrate function with the respective argument
 * set. If you have an equation where you test by both the values of the
 * test functions and the gradients, both template arguments need to be set
 * to true. Calling first the integrate function for values and then
 * gradients in a separate call leads to wrong results, since the second
 * call will internally overwrite the results from the first call. Note that
 * there is no function argument for the second derivative for integrate
 * step.  <li>Eventually, the local contributions in the vector
 * $v_\mathrm{cell}$ as mentioned in the introduction need to be added into
 * the result vector (and constraints are applied). This is done with a call
 * to @p distribute_local_to_global, the same name as the corresponding
 * function in the AffineConstraints (only that we now store the local vector
 * in the FEEvaluation object, as are the indices between local and global
 * degrees of freedom).  </ol>
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::local_apply(
 *     const MatrixFree<dim, number> &                   data,
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         AssertDimension(coefficient.size(0), data.n_cell_batches());
 *         AssertDimension(coefficient.size(1), phi.n_q_points);
 * 
 *         phi.reinit(cell);
 *         phi.read_dof_values(src);
 *         phi.evaluate(EvaluationFlags::gradients);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q);
 *         phi.integrate(EvaluationFlags::gradients);
 *         phi.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This function implements the loop over all cells for the
 * Base::apply_add() interface. This is done with the @p cell_loop of the
 * MatrixFree class, which takes the operator() of this class with arguments
 * MatrixFree, OutVector, InVector, cell_range. When working with MPI
 * parallelization (but no threading) as is done in this tutorial program,
 * the cell loop corresponds to the following three lines of code:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * src.update_ghost_values();
 * local_apply(*this->data, dst, src, std::make_pair(0U,
 *                                                   data.n_cell_batches()));
 * dst.compress(VectorOperation::add);
 * @endcode
 * </div>
 *   

 * 
 * Here, the two calls update_ghost_values() and compress() perform the data
 * exchange on processor boundaries for MPI, once for the source vector
 * where we need to read from entries owned by remote processors, and once
 * for the destination vector where we have accumulated parts of the
 * residuals that need to be added to the respective entry of the owner
 * processor. However, MatrixFree::cell_loop does not only abstract away
 * those two calls, but also performs some additional optimizations. On the
 * one hand, it will split the update_ghost_values() and compress() calls in
 * a way to allow for overlapping communication and computation. The
 * local_apply function is then called with three cell ranges representing
 * partitions of the cell range from 0 to MatrixFree::n_cell_batches(). On
 * the other hand, cell_loop also supports thread parallelism in which case
 * the cell ranges are split into smaller chunks and scheduled in an
 * advanced way that avoids access to the same vector entry by several
 * threads. That feature is explained in step-48.
 *   

 * 
 * Note that after the cell loop, the constrained degrees of freedom need to
 * be touched once more for sensible vmult() operators: Since the assembly
 * loop automatically resolves constraints (just as the
 * AffineConstraints::distribute_local_to_global() call does), it does not
 * compute any contribution for constrained degrees of freedom, leaving the
 * respective entries zero. This would represent a matrix that had empty
 * rows and columns for constrained degrees of freedom. However, iterative
 * solvers like CG only work for non-singular matrices. The easiest way to
 * do that is to set the sub-block of the matrix that corresponds to
 * constrained DoFs to an identity matrix, in which case application of the
 * matrix would simply copy the elements of the right hand side vector into
 * the left hand side. Fortunately, the vmult() implementations
 * MatrixFreeOperators::Base do this automatically for us outside the
 * apply_add() function, so we do not need to take further action here.
 *   

 * 
 * When using the combination of MatrixFree and FEEvaluation in parallel
 * with MPI, there is one aspect to be careful about &mdash; the indexing
 * used for accessing the vector. For performance reasons, MatrixFree and
 * FEEvaluation are designed to access vectors in MPI-local index space also
 * when working with multiple processors. Working in local index space means
 * that no index translation needs to be performed at the place the vector
 * access happens, apart from the unavoidable indirect addressing. However,
 * local index spaces are ambiguous: While it is standard convention to
 * access the locally owned range of a vector with indices between 0 and the
 * local size, the numbering is not so clear for the ghosted entries and
 * somewhat arbitrary. For the matrix-vector product, only the indices
 * appearing on locally owned cells (plus those referenced via hanging node
 * constraints) are necessary. However, in deal.II we often set all the
 * degrees of freedom on ghosted elements as ghosted vector entries, called
 * the @ref GlossLocallyRelevantDof "locally relevant DoFs described in the
 * glossary". In that case, the MPI-local index of a ghosted vector entry
 * can in general be different in the two possible ghost sets, despite
 * referring to the same global index. To avoid problems, FEEvaluation
 * checks that the partitioning of the vector used for the matrix-vector
 * product does indeed match with the partitioning of the indices in
 * MatrixFree by a check called
 * LinearAlgebra::distributed::Vector::partitioners_are_compatible. To
 * facilitate things, the MatrixFreeOperators::Base class includes a
 * mechanism to fit the ghost set to the correct layout. This happens in the
 * ghost region of the vector, so keep in mind that the ghost region might
 * be modified in both the destination and source vector after a call to a
 * vmult() method. This is legitimate because the ghost region of a
 * distributed deal.II vector is a mutable section and filled on
 * demand. Vectors used in matrix-vector products must not be ghosted upon
 * entry of vmult() functions, so no information gets lost.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::apply_add(
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src) const
 *   {
 *     this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The following function implements the computation of the diagonal of the
 * operator. Computing matrix entries of a matrix-free operator evaluation
 * turns out to be more complicated than evaluating the
 * operator. Fundamentally, we could obtain a matrix representation of the
 * operator by applying the operator on <i>all</i> unit vectors. Of course,
 * that would be very inefficient since we would need to perform <i>n</i>
 * operator evaluations to retrieve the whole matrix. Furthermore, this
 * approach would completely ignore the matrix sparsity. On an individual
 * cell, however, this is the way to go and actually not that inefficient as
 * there usually is a coupling between all degrees of freedom inside the
 * cell.
 *   

 * 
 * We first initialize the diagonal vector to the correct parallel
 * layout. This vector is encapsulated in a member called
 * inverse_diagonal_entries of type DiagonalMatrix in the base class
 * MatrixFreeOperators::Base. This member is a shared pointer that we first
 * need to initialize and then get the vector representing the diagonal
 * entries in the matrix. As to the actual diagonal computation, we again
 * use the cell_loop infrastructure of MatrixFree to invoke a local worker
 * routine called local_compute_diagonal(). Since we will only write into a
 * vector but not have any source vector, we put a dummy argument of type
 * <tt>unsigned int</tt> in place of the source vector to confirm with the
 * cell_loop interface. After the loop, we need to set the vector entries
 * subject to Dirichlet boundary conditions to one (either those on the
 * boundary described by the AffineConstraints object inside MatrixFree or
 * the indices at the interface between different grid levels in adaptive
 * multigrid). This is done through the function
 * MatrixFreeOperators::Base::set_constrained_entries_to_one() and matches
 * with the setting in the matrix-vector product provided by the Base
 * operator. Finally, we need to invert the diagonal entries which is the
 * form required by the Chebyshev smoother based on the Jacobi iteration. In
 * the loop, we assert that all entries are non-zero, because they should
 * either have obtained a positive contribution from integrals or be
 * constrained and treated by @p set_constrained_entries_to_one() following
 * cell_loop.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
 *   {
 *     this->inverse_diagonal_entries.reset(
 *       new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
 *     LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
 *       this->inverse_diagonal_entries->get_vector();
 *     this->data->initialize_dof_vector(inverse_diagonal);
 *     unsigned int dummy = 0;
 *     this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
 *                           this,
 *                           inverse_diagonal,
 *                           dummy);
 * 
 *     this->set_constrained_entries_to_one(inverse_diagonal);
 * 
 *     for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
 *       {
 *         Assert(inverse_diagonal.local_element(i) > 0.,
 *                ExcMessage("No diagonal entry in a positive definite operator "
 *                           "should be zero"));
 *         inverse_diagonal.local_element(i) =
 *           1. / inverse_diagonal.local_element(i);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In the local compute loop, we compute the diagonal by a loop over all
 * columns in the local matrix and putting the entry 1 in the <i>i</i>th
 * slot and a zero entry in all other slots, i.e., we apply the cell-wise
 * differential operator on one unit vector at a time. The inner part
 * invoking FEEvaluation::evaluate, the loop over quadrature points, and
 * FEEvalution::integrate, is exactly the same as in the local_apply
 * function. Afterwards, we pick out the <i>i</i>th entry of the local
 * result and put it to a temporary storage (as we overwrite all entries in
 * the array behind FEEvaluation::get_dof_value() with the next loop
 * iteration). Finally, the temporary storage is written to the destination
 * vector. Note how we use FEEvaluation::get_dof_value() and
 * FEEvaluation::submit_dof_value() to read and write to the data field that
 * FEEvaluation uses for the integration on the one hand and writes into the
 * global vector on the other hand.
 *   

 * 
 * Given that we are only interested in the matrix diagonal, we simply throw
 * away all other entries of the local matrix that have been computed along
 * the way. While it might seem wasteful to compute the complete cell matrix
 * and then throw away everything but the diagonal, the integration are so
 * efficient that the computation does not take too much time. Note that the
 * complexity of operator evaluation per element is $\mathcal
 * O((p+1)^{d+1})$ for polynomial degree $k$, so computing the whole matrix
 * costs us $\mathcal O((p+1)^{2d+1})$ operations, not too far away from
 * $\mathcal O((p+1)^{2d})$ complexity for computing the diagonal with
 * FEValues. Since FEEvaluation is also considerably faster due to
 * vectorization and other optimizations, the diagonal computation with this
 * function is actually the fastest (simple) variant. (It would be possible
 * to compute the diagonal with sum factorization techniques in $\mathcal
 * O((p+1)^{d+1})$ operations involving specifically adapted
 * kernels&mdash;but since such kernels are only useful in that particular
 * context and the diagonal computation is typically not on the critical
 * path, they have not been implemented in deal.II.)
 *   

 * 
 * Note that the code that calls distribute_local_to_global on the vector to
 * accumulate the diagonal entries into the global matrix has some
 * limitations. For operators with hanging node constraints that distribute
 * an integral contribution of a constrained DoF to several other entries
 * inside the distribute_local_to_global call, the vector interface used
 * here does not exactly compute the diagonal entries, but lumps some
 * contributions located on the diagonal of the local matrix that would end
 * up in a off-diagonal position of the global matrix to the diagonal. The
 * result is correct up to discretization accuracy as explained in <a
 * href="http://dx.doi.org/10.4208/cicp.101214.021015a">Kormann (2016),
 * section 5.3</a>, but not mathematically equal. In this tutorial program,
 * no harm can happen because the diagonal is only used for the multigrid
 * level matrices where no hanging node constraints appear.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
 *     const MatrixFree<dim, number> &             data,
 *     LinearAlgebra::distributed::Vector<number> &dst,
 *     const unsigned int &,
 *     const std::pair<unsigned int, unsigned int> &cell_range) const
 *   {
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);
 * 
 *     AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         AssertDimension(coefficient.size(0), data.n_cell_batches());
 *         AssertDimension(coefficient.size(1), phi.n_q_points);
 * 
 *         phi.reinit(cell);
 *         for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
 *               phi.submit_dof_value(VectorizedArray<number>(), j);
 *             phi.submit_dof_value(make_vectorized_array<number>(1.), i);
 * 
 *             phi.evaluate(EvaluationFlags::gradients);
 *             for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *               phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q),
 *                                   q);
 *             phi.integrate(EvaluationFlags::gradients);
 *             diagonal[i] = phi.get_dof_value(i);
 *           }
 *         for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
 *           phi.submit_dof_value(diagonal[i], i);
 *         phi.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemclass"></a> 
 * <h3>LaplaceProblem class</h3>
 * 

 * 
 * This class is based on the one in step-16. However, we replaced the
 * SparseMatrix<double> class by our matrix-free implementation, which means
 * that we can also skip the sparsity patterns. Notice that we define the
 * LaplaceOperator class with the degree of finite element as template
 * argument (the value is defined at the top of the file), and that we use
 * float numbers for the multigrid level matrices.
 *   

 * 
 * The class also has a member variable to keep track of all the detailed
 * timings for setting up the entire chain of data before we actually go
 * about solving the problem. In addition, there is an output stream (that
 * is disabled by default) that can be used to output details for the
 * individual setup operations instead of the summary only that is printed
 * out by default.
 *   

 * 
 * Since this program is designed to be used with MPI, we also provide the
 * usual @p pcout output stream that only prints the information of the
 * processor with MPI rank 0. The grid used for this programs can either be
 * a distributed triangulation based on p4est (in case deal.II is configured
 * to use p4est), otherwise it is a serial grid that only runs without MPI.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_rhs();
 *     void solve();
 *     void output_results(const unsigned int cycle) const;
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 * 
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     MappingQ1<dim> mapping;
 * 
 *     AffineConstraints<double> constraints;
 *     using SystemMatrixType =
 *       LaplaceOperator<dim, degree_finite_element, double>;
 *     SystemMatrixType system_matrix;
 * 
 *     MGConstrainedDoFs mg_constrained_dofs;
 *     using LevelMatrixType = LaplaceOperator<dim, degree_finite_element, float>;
 *     MGLevelObject<LevelMatrixType> mg_matrices;
 * 
 *     LinearAlgebra::distributed::Vector<double> solution;
 *     LinearAlgebra::distributed::Vector<double> system_rhs;
 * 
 *     double             setup_time;
 *     ConditionalOStream pcout;
 *     ConditionalOStream time_details;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * When we initialize the finite element, we of course have to use the
 * degree specified at the top of the file as well (otherwise, an exception
 * will be thrown at some point, since the computational kernel defined in
 * the templated LaplaceOperator class and the information from the finite
 * element read out by MatrixFree will not match). The constructor of the
 * triangulation needs to set an additional flag that tells the grid to
 * conform to the 2:1 cell balance over vertices, which is needed for the
 * convergence of the geometric multigrid routines. For the distributed
 * grid, we also need to specifically enable the multigrid hierarchy.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem()
 *     :
 * #ifdef DEAL_II_WITH_P4EST
 *     triangulation(
 *       MPI_COMM_WORLD,
 *       Triangulation<dim>::limit_level_difference_at_vertices,
 *       parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
 *     ,
 * #else
 *     triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     ,
 * #endif
 *     fe(degree_finite_element)
 *     , dof_handler(triangulation)
 *     , setup_time(0.)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *     ,
 * @endcode
 * 
 * The LaplaceProblem class holds an additional output stream that
 * collects detailed timings about the setup phase. This stream, called
 * time_details, is disabled by default through the @p false argument
 * specified here. For detailed timings, removing the @p false argument
 * prints all the details.
 * 
 * @code
 *     time_details(std::cout,
 *                  false && Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * The setup stage is in analogy to step-16 with relevant changes due to the
 * LaplaceOperator class. The first thing to do is to set up the DoFHandler,
 * including the degrees of freedom for the multigrid levels, and to
 * initialize constraints from hanging nodes and homogeneous Dirichlet
 * conditions. Since we intend to use this programs in %parallel with MPI,
 * we need to make sure that the constraints get to know the locally
 * relevant degrees of freedom, otherwise the storage would explode when
 * using more than a few hundred millions of degrees of freedom, see
 * step-40.
 * 

 * 
 * Once we have created the multigrid dof_handler and the constraints, we
 * can call the reinit function for the global matrix operator as well as
 * each level of the multigrid scheme. The main action is to set up the
 * <code> MatrixFree </code> instance for the problem. The base class of the
 * <code>LaplaceOperator</code> class, MatrixFreeOperators::Base, is
 * initialized with a shared pointer to MatrixFree object. This way, we can
 * simply create it here and then pass it on to the system matrix and level
 * matrices, respectively. For setting up MatrixFree, we need to activate
 * the update flag in the AdditionalData field of MatrixFree that enables
 * the storage of quadrature point coordinates in real space (by default, it
 * only caches data for gradients (inverse transposed Jacobians) and JxW
 * values). Note that if we call the reinit function without specifying the
 * level (i.e., giving <code>level = numbers::invalid_unsigned_int</code>),
 * MatrixFree constructs a loop over the active cells. In this tutorial, we
 * do not use threads in addition to MPI, which is why we explicitly disable
 * it by setting the MatrixFree::AdditionalData::tasks_parallel_scheme to
 * MatrixFree::AdditionalData::none. Finally, the coefficient is evaluated
 * and vectors are initialized as explained above.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 *     Timer time;
 *     setup_time = 0;
 * 
 *     system_matrix.clear();
 *     mg_matrices.clear_elements();
 * 
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler.distribute_mg_dofs();
 * 
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl;
 * 
 *     IndexSet locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
 *     constraints.close();
 *     setup_time += time.wall_time();
 *     time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *     time.restart();
 * 
 *     {
 *       typename MatrixFree<dim, double>::AdditionalData additional_data;
 *       additional_data.tasks_parallel_scheme =
 *         MatrixFree<dim, double>::AdditionalData::none;
 *       additional_data.mapping_update_flags =
 *         (update_gradients | update_JxW_values | update_quadrature_points);
 *       std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
 *         new MatrixFree<dim, double>());
 *       system_mf_storage->reinit(mapping,
 *                                 dof_handler,
 *                                 constraints,
 *                                 QGauss<1>(fe.degree + 1),
 *                                 additional_data);
 *       system_matrix.initialize(system_mf_storage);
 *     }
 * 
 *     system_matrix.evaluate_coefficient(Coefficient<dim>());
 * 
 *     system_matrix.initialize_dof_vector(solution);
 *     system_matrix.initialize_dof_vector(system_rhs);
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *     time.restart();
 * 
 * @endcode
 * 
 * Next, initialize the matrices for the multigrid method on all the
 * levels. The data structure MGConstrainedDoFs keeps information about
 * the indices subject to boundary conditions as well as the indices on
 * edges between different refinement levels as described in the step-16
 * tutorial program. We then go through the levels of the mesh and
 * construct the constraints and matrices on each level. These follow
 * closely the construction of the system matrix on the original mesh,
 * except the slight difference in naming when accessing information on
 * the levels rather than the active cells.
 * 
 * @code
 *     const unsigned int nlevels = triangulation.n_global_levels();
 *     mg_matrices.resize(0, nlevels - 1);
 * 
 *     std::set<types::boundary_id> dirichlet_boundary;
 *     dirichlet_boundary.insert(0);
 *     mg_constrained_dofs.initialize(dof_handler);
 *     mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
 *                                                        dirichlet_boundary);
 * 
 *     for (unsigned int level = 0; level < nlevels; ++level)
 *       {
 *         IndexSet relevant_dofs;
 *         DoFTools::extract_locally_relevant_level_dofs(dof_handler,
 *                                                       level,
 *                                                       relevant_dofs);
 *         AffineConstraints<double> level_constraints;
 *         level_constraints.reinit(relevant_dofs);
 *         level_constraints.add_lines(
 *           mg_constrained_dofs.get_boundary_indices(level));
 *         level_constraints.close();
 * 
 *         typename MatrixFree<dim, float>::AdditionalData additional_data;
 *         additional_data.tasks_parallel_scheme =
 *           MatrixFree<dim, float>::AdditionalData::none;
 *         additional_data.mapping_update_flags =
 *           (update_gradients | update_JxW_values | update_quadrature_points);
 *         additional_data.mg_level = level;
 *         std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level(
 *           new MatrixFree<dim, float>());
 *         mg_mf_storage_level->reinit(mapping,
 *                                     dof_handler,
 *                                     level_constraints,
 *                                     QGauss<1>(fe.degree + 1),
 *                                     additional_data);
 * 
 *         mg_matrices[level].initialize(mg_mf_storage_level,
 *                                       mg_constrained_dofs,
 *                                       level);
 *         mg_matrices[level].evaluate_coefficient(Coefficient<dim>());
 *       }
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_rhs"></a> 
 * <h4>LaplaceProblem::assemble_rhs</h4>
 * 

 * 
 * The assemble function is very simple since all we have to do is to
 * assemble the right hand side. Thanks to FEEvaluation and all the data
 * cached in the MatrixFree class, which we query from
 * MatrixFreeOperators::Base, this can be done in a few lines. Since this
 * call is not wrapped into a MatrixFree::cell_loop (which would be an
 * alternative), we must not forget to call compress() at the end of the
 * assembly to send all the contributions of the right hand side to the
 * owner of the respective degree of freedom.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_rhs()
 *   {
 *     Timer time;
 * 
 *     system_rhs = 0;
 *     FEEvaluation<dim, degree_finite_element> phi(
 *       *system_matrix.get_matrix_free());
 *     for (unsigned int cell = 0;
 *          cell < system_matrix.get_matrix_free()->n_cell_batches();
 *          ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_value(make_vectorized_array<double>(1.0), q);
 *         phi.integrate(EvaluationFlags::values);
 *         phi.distribute_local_to_global(system_rhs);
 *       }
 *     system_rhs.compress(VectorOperation::add);
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * The solution process is similar as in step-16. We start with the setup of
 * the transfer. For LinearAlgebra::distributed::Vector, there is a very
 * fast transfer class called MGTransferMatrixFree that does the
 * interpolation between the grid levels with the same fast sum
 * factorization kernels that get also used in FEEvaluation.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     Timer                            time;
 *     MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs);
 *     mg_transfer.build(dof_handler);
 *     setup_time += time.wall_time();
 *     time_details << "MG build transfer time     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *     time.restart();
 * 
 * @endcode
 * 
 * As a smoother, this tutorial program uses a Chebyshev iteration instead
 * of SOR in step-16. (SOR would be very difficult to implement because we
 * do not have the matrix elements available explicitly, and it is
 * difficult to make it work efficiently in %parallel.)  The smoother is
 * initialized with our level matrices and the mandatory additional data
 * for the Chebyshev smoother. We use a relatively high degree here (5),
 * since matrix-vector products are comparably cheap. We choose to smooth
 * out a range of $[1.2 \hat{\lambda}_{\max}/15,1.2 \hat{\lambda}_{\max}]$
 * in the smoother where $\hat{\lambda}_{\max}$ is an estimate of the
 * largest eigenvalue (the factor 1.2 is applied inside
 * PreconditionChebyshev). In order to compute that eigenvalue, the
 * Chebyshev initialization performs a few steps of a CG algorithm
 * without preconditioner. Since the highest eigenvalue is usually the
 * easiest one to find and a rough estimate is enough, we choose 10
 * iterations. Finally, we also set the inner preconditioner type in the
 * Chebyshev method which is a Jacobi iteration. This is represented by
 * the DiagonalMatrix class that gets the inverse diagonal entry provided
 * by our LaplaceOperator class.
 *     

 * 
 * On level zero, we initialize the smoother differently because we want
 * to use the Chebyshev iteration as a solver. PreconditionChebyshev
 * allows the user to switch to solver mode where the number of iterations
 * is internally chosen to the correct value. In the additional data
 * object, this setting is activated by choosing the polynomial degree to
 * @p numbers::invalid_unsigned_int. The algorithm will then attack all
 * eigenvalues between the smallest and largest one in the coarse level
 * matrix. The number of steps in the Chebyshev smoother are chosen such
 * that the Chebyshev convergence estimates guarantee to reduce the
 * residual by the number specified in the variable @p
 * smoothing_range. Note that for solving, @p smoothing_range is a
 * relative tolerance and chosen smaller than one, in this case, we select
 * three orders of magnitude, whereas it is a number larger than 1 when
 * only selected eigenvalues are smoothed.
 *     

 * 
 * From a computational point of view, the Chebyshev iteration is a very
 * attractive coarse grid solver as long as the coarse size is
 * moderate. This is because the Chebyshev method performs only
 * matrix-vector products and vector updates, which typically parallelize
 * better to the largest cluster size with more than a few tens of
 * thousands of cores than inner product involved in other iterative
 * methods. The former involves only local communication between neighbors
 * in the (coarse) mesh, whereas the latter requires global communication
 * over all processors.
 * 
 * @code
 *     using SmootherType =
 *       PreconditionChebyshev<LevelMatrixType,
 *                             LinearAlgebra::distributed::Vector<float>>;
 *     mg::SmootherRelaxation<SmootherType,
 *                            LinearAlgebra::distributed::Vector<float>>
 *                                                          mg_smoother;
 *     MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
 *     smoother_data.resize(0, triangulation.n_global_levels() - 1);
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       {
 *         if (level > 0)
 *           {
 *             smoother_data[level].smoothing_range     = 15.;
 *             smoother_data[level].degree              = 5;
 *             smoother_data[level].eig_cg_n_iterations = 10;
 *           }
 *         else
 *           {
 *             smoother_data[0].smoothing_range = 1e-3;
 *             smoother_data[0].degree          = numbers::invalid_unsigned_int;
 *             smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
 *           }
 *         mg_matrices[level].compute_diagonal();
 *         smoother_data[level].preconditioner =
 *           mg_matrices[level].get_matrix_diagonal_inverse();
 *       }
 *     mg_smoother.initialize(mg_matrices, smoother_data);
 * 
 *     MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
 *       mg_coarse;
 *     mg_coarse.initialize(mg_smoother);
 * 
 * @endcode
 * 
 * The next step is to set up the interface matrices that are needed for the
 * case with hanging nodes. The adaptive multigrid realization in deal.II
 * implements an approach called local smoothing. This means that the
 * smoothing on the finest level only covers the local part of the mesh
 * defined by the fixed (finest) grid level and ignores parts of the
 * computational domain where the terminal cells are coarser than this
 * level. As the method progresses to coarser levels, more and more of the
 * global mesh will be covered. At some coarser level, the whole mesh will
 * be covered. Since all level matrices in the multigrid method cover a
 * single level in the mesh, no hanging nodes appear on the level matrices.
 * At the interface between multigrid levels, homogeneous Dirichlet boundary
 * conditions are set while smoothing. When the residual is transferred to
 * the next coarser level, however, the coupling over the multigrid
 * interface needs to be taken into account. This is done by the so-called
 * interface (or edge) matrices that compute the part of the residual that
 * is missed by the level matrix with
 * homogeneous Dirichlet conditions. We refer to the @ref mg_paper
 * "Multigrid paper by Janssen and Kanschat" for more details.
 *     

 * 
 * For the implementation of those interface matrices, there is already a
 * pre-defined class MatrixFreeOperators::MGInterfaceOperator that wraps
 * the routines MatrixFreeOperators::Base::vmult_interface_down() and
 * MatrixFreeOperators::Base::vmult_interface_up() in a new class with @p
 * vmult() and @p Tvmult() operations (that were originally written for
 * matrices, hence expecting those names). Note that vmult_interface_down
 * is used during the restriction phase of the multigrid V-cycle, whereas
 * vmult_interface_up is used during the prolongation phase.
 *     

 * 
 * Once the interface matrix is created, we set up the remaining Multigrid
 * preconditioner infrastructure in complete analogy to step-16 to obtain
 * a @p preconditioner object that can be applied to a matrix.
 * 
 * @code
 *     mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
 *       mg_matrices);
 * 
 *     MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
 *       mg_interface_matrices;
 *     mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       mg_interface_matrices[level].initialize(mg_matrices[level]);
 *     mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
 *       mg_interface_matrices);
 * 
 *     Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
 *       mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
 *     mg.set_edge_matrices(mg_interface, mg_interface);
 * 
 *     PreconditionMG<dim,
 *                    LinearAlgebra::distributed::Vector<float>,
 *                    MGTransferMatrixFree<dim, float>>
 *       preconditioner(dof_handler, mg, mg_transfer);
 * 
 * @endcode
 * 
 * The setup of the multigrid routines is quite easy and one cannot see
 * any difference in the solve process as compared to step-16. All the
 * magic is hidden behind the implementation of the LaplaceOperator::vmult
 * operation. Note that we print out the solve time and the accumulated
 * setup time through standard out, i.e., in any case, whereas detailed
 * times for the setup operations are only printed in case the flag for
 * detail_times in the constructor is changed.
 * 

 * 
 * 
 * @code
 *     SolverControl solver_control(100, 1e-12 * system_rhs.l2_norm());
 *     SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);
 *     setup_time += time.wall_time();
 *     time_details << "MG build smoother time     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *     pcout << "Total setup time               (wall) " << setup_time << "s\n";
 * 
 *     time.reset();
 *     time.start();
 *     constraints.set_zero(solution);
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 * 
 *     pcout << "Time solve (" << solver_control.last_step() << " iterations)"
 *           << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) "
 *           << time.cpu_time() << "s/" << time.wall_time() << "s\n";
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemoutput_results"></a> 
 * <h4>LaplaceProblem::output_results</h4>
 * 

 * 
 * Here is the data output, which is a simplified version of step-5. We use
 * the standard VTU (= compressed VTK) output for each grid produced in the
 * refinement process. In addition, we use a compression algorithm that is
 * optimized for speed rather than disk usage. The default setting (which
 * optimizes for disk usage) makes saving the output take about 4 times as
 * long as running the linear solver, while setting
 * DataOutBase::VtkFlags::compression_level to
 * DataOutBase::VtkFlags::best_speed lowers this to only one fourth the time
 * of the linear solve.
 *   

 * 
 * We disable the output when the mesh gets too large. A variant of this
 * program has been run on hundreds of thousands MPI ranks with as many as
 * 100 billion grid cells, which is not directly accessible to classical
 * visualization tools.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     Timer time;
 *     if (triangulation.n_global_active_cells() > 1000000)
 *       return;
 * 
 *     DataOut<dim> data_out;
 * 
 *     solution.update_ghost_values();
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches(mapping);
 * 
 *     DataOutBase::VtkFlags flags;
 *     flags.compression_level = DataOutBase::VtkFlags::best_speed;
 *     data_out.set_flags(flags);
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", cycle, MPI_COMM_WORLD, 3);
 * 
 *     time_details << "Time write output          (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * The function that runs the program is very similar to the one in
 * step-16. We do few refinement steps in 3D compared to 2D, but that's
 * it.
 *   

 * 
 * Before we run the program, we output some information about the detected
 * vectorization level as discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     {
 *       const unsigned int n_vect_doubles = VectorizedArray<double>::size();
 *       const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;
 * 
 *       pcout << "Vectorization over " << n_vect_doubles
 *             << " doubles = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl;
 *     }
 * 
 *     for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, 0., 1.);
 *             triangulation.refine_global(3 - dim);
 *           }
 *         triangulation.refine_global(1);
 *         setup_system();
 *         assemble_rhs();
 *         solve();
 *         output_results(cycle);
 *         pcout << std::endl;
 *       };
 *   }
 * } // namespace Step37
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
 * Apart from the fact that we set up the MPI framework according to step-40,
 * there are no surprises in the main function.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step37;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
 * 
 *       LaplaceProblem<dimension> laplace_problem;
 *       laplace_problem.run();
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


<a name="Programoutput"></a><h3>Program output</h3>


Since this example solves the same problem as step-5 (except for
a different coefficient), there is little to say about the
solution. We show a picture anyway, illustrating the size of the
solution through both isocontours and volume rendering:

<img src="https://www.dealii.org/images/steps/developer/step-37.solution.png" alt="">

Of more interest is to evaluate some aspects of the multigrid solver.
When we run this program in 2D for quadratic ($Q_2$) elements, we get the
following output (when run on one core in release mode):
@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 81
Total setup time               (wall) 0.00159788s
Time solve (6 iterations)  (CPU/wall) 0.000951s/0.000951052s

Cycle 1
Number of degrees of freedom: 289
Total setup time               (wall) 0.00114608s
Time solve (6 iterations)  (CPU/wall) 0.000935s/0.000934839s

Cycle 2
Number of degrees of freedom: 1089
Total setup time               (wall) 0.00244665s
Time solve (6 iterations)  (CPU/wall) 0.00207s/0.002069s

Cycle 3
Number of degrees of freedom: 4225
Total setup time               (wall) 0.00678205s
Time solve (6 iterations)  (CPU/wall) 0.005616s/0.00561595s

Cycle 4
Number of degrees of freedom: 16641
Total setup time               (wall) 0.0241671s
Time solve (6 iterations)  (CPU/wall) 0.019543s/0.0195441s

Cycle 5
Number of degrees of freedom: 66049
Total setup time               (wall) 0.0967851s
Time solve (6 iterations)  (CPU/wall) 0.07457s/0.0745709s

Cycle 6
Number of degrees of freedom: 263169
Total setup time               (wall) 0.346374s
Time solve (6 iterations)  (CPU/wall) 0.260042s/0.265033s
@endcode

As in step-16, we see that the number of CG iterations remains constant with
increasing number of degrees of freedom. A constant number of iterations
(together with optimal computational properties) means that the computing time
approximately quadruples as the problem size quadruples from one cycle to the
next. The code is also very efficient in terms of storage. Around 2-4 million
degrees of freedom fit into 1 GB of memory, see also the MPI results below. An
interesting fact is that solving one linear system is cheaper than the setup,
despite not building a matrix (approximately half of which is spent in the
DoFHandler::distribute_dofs() and DoFHandler::distribute_mg_dofs()
calls). This shows the high efficiency of this approach, but also that the
deal.II data structures are quite expensive to set up and the setup cost must
be amortized over several system solves.

Not much changes if we run the program in three spatial dimensions. Since we
use uniform mesh refinement, we get eight times as many elements and
approximately eight times as many degrees of freedom with each cycle:

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 125
Total setup time               (wall) 0.00231099s
Time solve (6 iterations)  (CPU/wall) 0.000692s/0.000922918s

Cycle 1
Number of degrees of freedom: 729
Total setup time               (wall) 0.00289083s
Time solve (6 iterations)  (CPU/wall) 0.001534s/0.0024128s

Cycle 2
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0143182s
Time solve (6 iterations)  (CPU/wall) 0.010785s/0.0107841s

Cycle 3
Number of degrees of freedom: 35937
Total setup time               (wall) 0.087064s
Time solve (6 iterations)  (CPU/wall) 0.063522s/0.06545s

Cycle 4
Number of degrees of freedom: 274625
Total setup time               (wall) 0.596306s
Time solve (6 iterations)  (CPU/wall) 0.427757s/0.431765s

Cycle 5
Number of degrees of freedom: 2146689
Total setup time               (wall) 4.96491s
Time solve (6 iterations)  (CPU/wall) 3.53126s/3.56142s
@endcode

Since it is so easy, we look at what happens if we increase the polynomial
degree. When selecting the degree as four in 3D, i.e., on $\mathcal Q_4$
elements, by changing the line <code>const unsigned int
degree_finite_element=4;</code> at the top of the program, we get the
following program output:

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 729
Total setup time               (wall) 0.00633097s
Time solve (6 iterations)  (CPU/wall) 0.002829s/0.00379395s

Cycle 1
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0174279s
Time solve (6 iterations)  (CPU/wall) 0.012255s/0.012254s

Cycle 2
Number of degrees of freedom: 35937
Total setup time               (wall) 0.082655s
Time solve (6 iterations)  (CPU/wall) 0.052362s/0.0523629s

Cycle 3
Number of degrees of freedom: 274625
Total setup time               (wall) 0.507943s
Time solve (6 iterations)  (CPU/wall) 0.341811s/0.345788s

Cycle 4
Number of degrees of freedom: 2146689
Total setup time               (wall) 3.46251s
Time solve (7 iterations)  (CPU/wall) 3.29638s/3.3265s

Cycle 5
Number of degrees of freedom: 16974593
Total setup time               (wall) 27.8989s
Time solve (7 iterations)  (CPU/wall) 26.3705s/27.1077s
@endcode

Since $\mathcal Q_4$ elements on a certain mesh correspond to $\mathcal Q_2$
elements on half the mesh size, we can compare the run time at cycle 4 with
fourth degree polynomials with cycle 5 using quadratic polynomials, both at
2.1 million degrees of freedom. The surprising effect is that the solver for
$\mathcal Q_4$ element is actually slightly faster than for the quadratic
case, despite using one more linear iteration. The effect that higher-degree
polynomials are similarly fast or even faster than lower degree ones is one of
the main strengths of matrix-free operator evaluation through sum
factorization, see the <a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">matrix-free
paper</a>. This is fundamentally different to matrix-based methods that get
more expensive per unknown as the polynomial degree increases and the coupling
gets denser.

In addition, also the setup gets a bit cheaper for higher order, which is
because fewer elements need to be set up.

Finally, let us look at the timings with degree 8, which corresponds to
another round of mesh refinement in the lower order methods:

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0842004s
Time solve (8 iterations)  (CPU/wall) 0.019296s/0.0192959s

Cycle 1
Number of degrees of freedom: 35937
Total setup time               (wall) 0.327048s
Time solve (8 iterations)  (CPU/wall) 0.07517s/0.075999s

Cycle 2
Number of degrees of freedom: 274625
Total setup time               (wall) 2.12335s
Time solve (8 iterations)  (CPU/wall) 0.448739s/0.453698s

Cycle 3
Number of degrees of freedom: 2146689
Total setup time               (wall) 16.1743s
Time solve (8 iterations)  (CPU/wall) 3.95003s/3.97717s

Cycle 4
Number of degrees of freedom: 16974593
Total setup time               (wall) 130.8s
Time solve (8 iterations)  (CPU/wall) 31.0316s/31.767s
@endcode

Here, the initialization seems considerably slower than before, which is
mainly due to the computation of the diagonal of the matrix, which actually
computes a 729 x 729 matrix on each cell and throws away everything but the
diagonal. The solver times, however, are again very close to the quartic case,
showing that the linear increase with the polynomial degree that is
theoretically expected is almost completely offset by better computational
characteristics and the fact that higher order methods have a smaller share of
degrees of freedom living on several cells that add to the evaluation
complexity.

<a name="Comparisonwithasparsematrix"></a><h3>Comparison with a sparse matrix</h3>


In order to understand the capabilities of the matrix-free implementation, we
compare the performance of the 3d example above with a sparse matrix
implementation based on TrilinosWrappers::SparseMatrix by measuring both the
computation times for the initialization of the problem (distribute DoFs,
setup and assemble matrices, setup multigrid structures) and the actual
solution for the matrix-free variant and the variant based on sparse
matrices. We base the preconditioner on float numbers and the actual matrix
and vectors on double numbers, as shown above. Tests are run on an Intel Core
i7-5500U notebook processor (two cores and <a
href="http://en.wikipedia.org/wiki/Advanced_Vector_Extensions">AVX</a>
support, i.e., four operations on doubles can be done with one CPU
instruction, which is heavily used in FEEvaluation), optimized mode, and two
MPI ranks.

<table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="2">Sparse matrix</th>
    <th colspan="2">Matrix-free implementation</th>
  </tr>
  <tr>
    <th>n_dofs</th>
    <th>Setup + assemble</th>
    <th>&nbsp;Solve&nbsp;</th>
    <th>Setup + assemble</th>
    <th>&nbsp;Solve&nbsp;</th>
  </tr>
  <tr>
    <td align="right">125</td>
    <td align="center">0.0042s</td>
    <td align="center">0.0012s</td>
    <td align="center">0.0022s</td>
    <td align="center">0.00095s</td>
  </tr>
  <tr>
    <td align="right">729</td>
    <td align="center">0.012s</td>
    <td align="center">0.0040s</td>
    <td align="center">0.0027s</td>
    <td align="center">0.0021s</td>
  </tr>
  <tr>
    <td align="right">4,913</td>
    <td align="center">0.082s</td>
    <td align="center">0.012s</td>
    <td align="center">0.011s</td>
    <td align="center">0.0057s</td>
  </tr>
  <tr>
    <td align="right">35,937</td>
    <td align="center">0.73s</td>
    <td align="center">0.13s</td>
    <td align="center">0.048s</td>
    <td align="center">0.040s</td>
  </tr>
  <tr>
    <td align="right">274,625</td>
    <td align="center">5.43s</td>
    <td align="center">1.01s</td>
    <td align="center">0.33s</td>
    <td align="center">0.25s</td>
  </tr>
  <tr>
    <td align="right">2,146,689</td>
    <td align="center">43.8s</td>
    <td align="center">8.24s</td>
    <td align="center">2.42s</td>
    <td align="center">2.06s</td>
  </tr>
</table>

The table clearly shows that the matrix-free implementation is more than twice
as fast for the solver, and more than six times as fast when it comes to
initialization costs. As the problem size is made a factor 8 larger, we note
that the times usually go up by a factor eight, too (as the solver iterations
are constant at six). The main deviation is in the sparse matrix between 5k
and 36k degrees of freedom, where the time increases by a factor 12. This is
the threshold where the (L3) cache in the processor can no longer hold all
data necessary for the matrix-vector products and all matrix elements must be
fetched from main memory.

Of course, this picture does not necessarily translate to all cases, as there
are problems where knowledge of matrix entries enables much better solvers (as
happens when the coefficient is varying more strongly than in the above
example). Moreover, it also depends on the computer system. The present system
has good memory performance, so sparse matrices perform comparably
well. Nonetheless, the matrix-free implementation gives a nice speedup already
for the <i>Q</i><sub>2</sub> elements used in this example. This becomes
particularly apparent for time-dependent or nonlinear problems where sparse
matrices would need to be reassembled over and over again, which becomes much
easier with this class. And of course, thanks to the better complexity of the
products, the method gains increasingly larger advantages when the order of the
elements increases (the matrix-free implementation has costs
4<i>d</i><sup>2</sup><i>p</i> per degree of freedom, compared to
2<i>p<sup>d</sup></i> for the sparse matrix, so it will win anyway for order 4
and higher in 3d).

<a name="ResultsforlargescaleparallelcomputationsonSuperMUC"></a><h3> Results for large-scale parallel computations on SuperMUC</h3>


As explained in the introduction and the in-code comments, this program can be
run in parallel with MPI. It turns out that geometric multigrid schemes work
really well and can scale to very large machines. To the authors' knowledge,
the geometric multigrid results shown here are the largest computations done
with deal.II as of late 2016, run on up to 147,456 cores of the <a
href="https://www.lrz.de/services/compute/supermuc/systemdescription/">complete
SuperMUC Phase 1</a>. The ingredients for scalability beyond 1000 cores are
that no data structure that depends on the global problem size is held in its
entirety on a single processor and that the communication is not too frequent
in order not to run into latency issues of the network.  For PDEs solved with
iterative solvers, the communication latency is often the limiting factor,
rather than the throughput of the network. For the example of the SuperMUC
system, the point-to-point latency between two processors is between 1e-6 and
1e-5 seconds, depending on the proximity in the MPI network. The matrix-vector
products with @p LaplaceOperator from this class involves several
point-to-point communication steps, interleaved with computations on each
core. The resulting latency of a matrix-vector product is around 1e-4
seconds. Global communication, for example an @p MPI_Allreduce operation that
accumulates the sum of a single number per rank over all ranks in the MPI
network, has a latency of 1e-4 seconds. The multigrid V-cycle used in this
program is also a form of global communication. Think about the coarse grid
solve that happens on a single processor: It accumulates the contributions
from all processors before it starts. When completed, the coarse grid solution
is transferred to finer levels, where more and more processors help in
smoothing until the fine grid. Essentially, this is a tree-like pattern over
the processors in the network and controlled by the mesh. As opposed to the
@p MPI_Allreduce operations where the tree in the reduction is optimized to the
actual links in the MPI network, the multigrid V-cycle does this according to
the partitioning of the mesh. Thus, we cannot expect the same
optimality. Furthermore, the multigrid cycle is not simply a walk up and down
the refinement tree, but also communication on each level when doing the
smoothing. In other words, the global communication in multigrid is more
challenging and related to the mesh that provides less optimization
opportunities. The measured latency of the V-cycle is between 6e-3 and 2e-2
seconds, i.e., the same as 60 to 200 MPI_Allreduce operations.

The following figure shows a scaling experiments on $\mathcal Q_3$
elements. Along the lines, the problem size is held constant as the number of
cores is increasing. When doubling the number of cores, one expects a halving
of the computational time, indicated by the dotted gray lines. The results
show that the implementation shows almost ideal behavior until an absolute
time of around 0.1 seconds is reached. The solver tolerances have been set
such that the solver performs five iterations. This way of plotting data is
the <b>strong scaling</b> of the algorithm. As we go to very large core
counts, the curves flatten out a bit earlier, which is because of the
communication network in SuperMUC where communication between processors
farther away is slightly slower.

<img src="https://www.dealii.org/images/steps/developer/step-37.scaling_strong.png" alt="">

In addition, the plot also contains results for <b>weak scaling</b> that lists
how the algorithm behaves as both the number of processor cores and elements
is increased at the same pace. In this situation, we expect that the compute
time remains constant. Algorithmically, the number of CG iterations is
constant at 5, so we are good from that end. The lines in the plot are
arranged such that the top left point in each data series represents the same
size per processor, namely 131,072 elements (or approximately 3.5 million
degrees of freedom per core). The gray lines indicating ideal strong scaling
are by the same factor of 8 apart. The results show again that the scaling is
almost ideal. The parallel efficiency when going from 288 cores to 147,456
cores is at around 75% for a local problem size of 750,000 degrees of freedom
per core which takes 1.0s on 288 cores, 1.03s on 2304 cores, 1.19s on 18k
cores, and 1.35s on 147k cores. The algorithms also reach a very high
utilization of the processor. The largest computation on 147k cores reaches
around 1.7 PFLOPs/s on SuperMUC out of an arithmetic peak of 3.2 PFLOPs/s. For
an iterative PDE solver, this is a very high number and significantly more is
often only reached for dense linear algebra. Sparse linear algebra is limited
to a tenth of this value.

As mentioned in the introduction, the matrix-free method reduces the memory
consumption of the data structures. Besides the higher performance due to less
memory transfer, the algorithms also allow for very large problems to fit into
memory. The figure below shows the computational time as we increase the
problem size until an upper limit where the computation exhausts memory. We do
this for 1k cores, 8k cores, and 65k cores and see that the problem size can
be varied over almost two orders of magnitude with ideal scaling. The largest
computation shown in this picture involves 292 billion ($2.92 \cdot 10^{11}$)
degrees of freedom. On a DG computation of 147k cores, the above algorithms
were also run involving up to 549 billion (2^39) DoFs.

<img src="https://www.dealii.org/images/steps/developer/step-37.scaling_size.png" alt="">

Finally, we note that while performing the tests on the large-scale system
shown above, improvements of the multigrid algorithms in deal.II have been
developed. The original version contained the sub-optimal code based on
MGSmootherPrecondition where some MPI_Allreduce commands (checking whether
all vector entries are zero) were done on each smoothing
operation on each level, which only became apparent on 65k cores and
more. However, the following picture shows that the improvement already pay
off on a smaller scale, here shown on computations on up to 14,336 cores for
$\mathcal Q_5$ elements:

<img src="https://www.dealii.org/images/steps/developer/step-37.scaling_oldnew.png" alt="">


<a name="Adaptivity"></a><h3> Adaptivity</h3>


As explained in the code, the algorithm presented here is prepared to run on
adaptively refined meshes. If only part of the mesh is refined, the multigrid
cycle will run with local smoothing and impose Dirichlet conditions along the
interfaces which differ in refinement level for smoothing through the
MatrixFreeOperators::Base class. Due to the way the degrees of freedom are
distributed over levels, relating the owner of the level cells to the owner of
the first descendant active cell, there can be an imbalance between different
processors in MPI, which limits scalability by a factor of around two to five.

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions</h3>


<a name="Kellyerrorestimator"></a><h4> Kelly error estimator </h4>


As mentioned above the code is ready for locally adaptive h-refinement.
For the Poisson equation one can employ the Kelly error indicator,
implemented in the KellyErrorEstimator class. However one needs to be careful
with the ghost indices of parallel vectors.
In order to evaluate the jump terms in the error indicator, each MPI process
needs to know locally relevant DoFs.
However MatrixFree::initialize_dof_vector() function initializes the vector only with
some locally relevant DoFs.
The ghost indices made available in the vector are a tight set of only those indices
that are touched in the cell integrals (including constraint resolution).
This choice has performance reasons, because sending all locally relevant degrees
of freedom would be too expensive compared to the matrix-vector product.
Consequently the solution vector as-is is
not suitable for the KellyErrorEstimator class.
The trick is to change the ghost part of the partition, for example using a
temporary vector and LinearAlgebra::distributed::Vector::copy_locally_owned_data_from()
as shown below.

@code
IndexSet locally_relevant_dofs;
DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
LinearAlgebra::distributed::Vector<double> copy_vec(solution);
solution.reinit(dof_handler.locally_owned_dofs(),
                locally_relevant_dofs,
                triangulation.get_communicator());
solution.copy_locally_owned_data_from(copy_vec);
constraints.distribute(solution);
solution.update_ghost_values();
@endcode

<a name="Sharedmemoryparallelization"></a><h4> Shared-memory parallelization</h4>


This program is parallelized with MPI only. As an alternative, the MatrixFree
loop can also be issued in hybrid mode, for example by using MPI parallelizing
over the nodes of a cluster and with threads through Intel TBB within the
shared memory region of one node. To use this, one would need to both set the
number of threads in the MPI_InitFinalize data structure in the main function,
and set the MatrixFree::AdditionalData::tasks_parallel_scheme to
partition_color to actually do the loop in parallel. This use case is
discussed in step-48.

<a name="InhomogeneousDirichletboundaryconditions"></a><h4> Inhomogeneous Dirichlet boundary conditions </h4>


The presented program assumes homogeneous Dirichlet boundary conditions. When
going to non-homogeneous conditions, the situation is a bit more intricate. To
understand how to implement such a setting, let us first recall how these
arise in the mathematical formulation and how they are implemented in a
matrix-based variant. In essence, an inhomogeneous Dirichlet condition sets
some of the nodal values in the solution to given values rather than
determining them through the variational principles,
@f{eqnarray*}
u_h(\mathbf{x}) = \sum_{i\in \mathcal N} \varphi_i(\mathbf{x}) u_i =
\sum_{i\in \mathcal N \setminus \mathcal N_D} \varphi_i(\mathbf{x}) u_i +
\sum_{i\in \mathcal N_D} \varphi_i(\mathbf{x}) g_i,
@f}
where $u_i$ denotes the nodal values of the solution and $\mathcal N$ denotes
the set of all nodes. The set $\mathcal N_D\subset \mathcal N$ is the subset
of the nodes that are subject to Dirichlet boundary conditions where the
solution is forced to equal $u_i = g_i = g(\mathbf{x}_i)$ as the interpolation
of boundary values on the Dirichlet-constrained node points $i\in \mathcal
N_D$. We then insert this solution
representation into the weak form, e.g. the Laplacian shown above, and move
the known quantities to the right hand side:
@f{eqnarray*}
(\nabla \varphi_i, \nabla u_h)_\Omega &=& (\varphi_i, f)_\Omega \quad \Rightarrow \\
\sum_{j\in \mathcal N \setminus \mathcal N_D}(\nabla \varphi_i,\nabla \varphi_j)_\Omega \, u_j &=&
(\varphi_i, f)_\Omega
-\sum_{j\in \mathcal N_D} (\nabla \varphi_i,\nabla\varphi_j)_\Omega\, g_j.
@f}
In this formula, the equations are tested for all basis functions $\varphi_i$
with $i\in N \setminus \mathcal N_D$ that are not related to the nodes
constrained by Dirichlet conditions.

In the implementation in deal.II, the integrals $(\nabla \varphi_i,\nabla \varphi_j)_\Omega$
on the right hand side are already contained in the local matrix contributions
we assemble on each cell. When using
AffineConstraints::distribute_local_to_global() as first described in the
step-6 and step-7 tutorial programs, we can account for the contribution of
inhomogeneous constraints <i>j</i> by multiplying the columns <i>j</i> and
rows <i>i</i> of the local matrix according to the integrals $(\varphi_i,
\varphi_j)_\Omega$ by the inhomogeneities and subtracting the resulting from
the position <i>i</i> in the global right-hand-side vector, see also the @ref
constraints module. In essence, we use some of the integrals that get
eliminated from the left hand side of the equation to finalize the right hand
side contribution. Similar mathematics are also involved when first writing
all entries into a left hand side matrix and then eliminating matrix rows and
columns by MatrixTools::apply_boundary_values().

In principle, the components that belong to the constrained degrees of freedom
could be eliminated from the linear system because they do not carry any
information. In practice, in deal.II we always keep the size of the linear
system the same to avoid handling two different numbering systems and avoid
confusion about the two different index sets. In order to ensure that the
linear system does not get singular when not adding anything to constrained
rows, we then add dummy entries to the matrix diagonal that are otherwise
unrelated to the real entries.

In a matrix-free method, we need to take a different approach, since the @p
LaplaceOperator class represents the matrix-vector product of a
<b>homogeneous</b> operator (the left-hand side of the last formula).  It does
not matter whether the AffineConstraints object passed to the
MatrixFree::reinit() contains inhomogeneous constraints or not, the
MatrixFree::cell_loop() call will only resolve the homogeneous part of the
constraints as long as it represents a <b>linear</b> operator.

In our matrix-free code, the right hand side computation where the
contribution of inhomogeneous conditions ends up is completely decoupled from
the matrix operator and handled by a different function above. Thus, we need
to explicitly generate the data that enters the right hand side rather than
using a byproduct of the matrix assembly. Since we already know how to apply
the operator on a vector, we could try to use those facilities for a vector
where we only set the Dirichlet values:
@code
  // interpolate boundary values on vector solution
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(mapping,
                                           dof_handler,
                                           0,
                                           BoundaryValueFunction<dim>(),
                                           boundary_values);
  for (const std::pair<const types::global_dof_index, double> &pair : boundary_values)
    if (solution.locally_owned_elements().is_element(pair.first))
      solution(pair.first) = pair.second;
@endcode
or, equivalently, if we already had filled the inhomogeneous constraints into
an AffineConstraints object,
@code
  solution = 0;
  constraints.distribute(solution);
@endcode

We could then pass the vector @p solution to the @p
LaplaceOperator::vmult_add() function and add the new contribution to the @p
system_rhs vector that gets filled in the @p LaplaceProblem::assemble_rhs()
function. However, this idea does not work because the
FEEvaluation::read_dof_values() call used inside the vmult() functions assumes
homogeneous values on all constraints (otherwise the operator would not be a
linear operator but an affine one). To also retrieve the values of the
inhomogeneities, we could select one of two following strategies.

<a name="UseFEEvaluationread_dof_values_plaintoavoidresolvingconstraints"></a><h5> Use FEEvaluation::read_dof_values_plain() to avoid resolving constraints </h5>


The class FEEvaluation has a facility that addresses precisely this
requirement: For non-homogeneous Dirichlet values, we do want to skip the
implicit imposition of homogeneous (Dirichlet) constraints upon reading the
data from the vector @p solution. For example, we could extend the @p
LaplaceProblem::assemble_rhs() function to deal with inhomogeneous Dirichlet
values as follows, assuming the Dirichlet values have been interpolated into
the object @p constraints:
@code
template <int dim>
void LaplaceProblem<dim>::assemble_rhs()
{
  solution = 0;
  constraints.distribute(solution);
  solution.update_ghost_values();
  system_rhs = 0;

  const Table<2, VectorizedArray<double>> &coefficient = system_matrix.get_coefficient();
  FEEvaluation<dim, degree_finite_element> phi(*system_matrix.get_matrix_free());
  for (unsigned int cell = 0;
       cell < system_matrix.get_matrix_free()->n_cell_batches();
       ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(solution);
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          phi.submit_gradient(-coefficient(cell, q) * phi.get_gradient(q), q);
          phi.submit_value(make_vectorized_array<double>(1.0), q);
        }
      phi.integrate(EvaluationFlags::values|EvaluationFlags::gradients);
      phi.distribute_local_to_global(system_rhs);
    }
  system_rhs.compress(VectorOperation::add);
}
@endcode

In this code, we replaced the FEEvaluation::read_dof_values() function for the
tentative solution vector by FEEvaluation::read_dof_values_plain() that
ignores all constraints. Due to this setup, we must make sure that other
constraints, e.g. by hanging nodes, are correctly distributed to the input
vector already as they are not resolved as in
FEEvaluation::read_dof_values_plain(). Inside the loop, we then evaluate the
Laplacian and repeat the second derivative call with
FEEvaluation::submit_gradient() from the @p LaplaceOperator class, but with the
sign switched since we wanted to subtract the contribution of Dirichlet
conditions on the right hand side vector according to the formula above. When
we invoke the FEEvaluation::integrate() call, we then set both arguments
regarding the value slot and first derivative slot to true to account for both
terms added in the loop over quadrature points. Once the right hand side is
assembled, we then go on to solving the linear system for the homogeneous
problem, say involving a variable @p solution_update. After solving, we can
add @p solution_update to the @p solution vector that contains the final
(inhomogeneous) solution.

Note that the negative sign for the Laplacian alongside with a positive sign
for the forcing that we needed to build the right hand side is a more general
concept: We have implemented nothing else than Newton's method for nonlinear
equations, but applied to a linear system. We have used an initial guess for
the variable @p solution in terms of the Dirichlet boundary conditions and
computed a residual $r = f - Au_0$. The linear system was then solved as
$\Delta u = A^{-1} (f-Au)$ and we finally computed $u = u_0 + \Delta u$. For a
linear system, we obviously reach the exact solution after a single
iteration. If we wanted to extend the code to a nonlinear problem, we would
rename the @p assemble_rhs() function into a more descriptive name like @p
assemble_residual() that computes the (weak) form of the residual, whereas the
@p LaplaceOperator::apply_add() function would get the linearization of the
residual with respect to the solution variable.

<a name="UseLaplaceOperatorwithasecondAffineConstraintsobjectwithoutDirichletconditions"></a><h5> Use LaplaceOperator with a second AffineConstraints object without Dirichlet conditions </h5>


A second alternative to get the right hand side that re-uses the @p
LaplaceOperator::apply_add() function is to instead add a second LaplaceOperator
that skips Dirichlet constraints. To do this, we initialize a second MatrixFree
object which does not have any boundary value constraints. This @p matrix_free
object is then passed to a @p LaplaceOperator class instance @p
inhomogeneous_operator that is only used to create the right hand side:
@code
template <int dim>
void LaplaceProblem<dim>::assemble_rhs()
{
  system_rhs = 0;
  AffineConstraints<double> no_constraints;
  no_constraints.close();
  LaplaceOperator<dim, degree_finite_element, double> inhomogeneous_operator;

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_gradients | update_JxW_values | update_quadrature_points);
  std::shared_ptr<MatrixFree<dim, double>> matrix_free(
    new MatrixFree<dim, double>());
  matrix_free->reinit(dof_handler,
                      no_constraints,
                      QGauss<1>(fe.degree + 1),
                      additional_data);
  inhomogeneous_operator.initialize(matrix_free);

  solution = 0.0;
  constraints.distribute(solution);
  inhomogeneous_operator.evaluate_coefficient(Coefficient<dim>());
  inhomogeneous_operator.vmult(system_rhs, solution);
  system_rhs *= -1.0;

  FEEvaluation<dim, degree_finite_element> phi(
    *inhomogeneous_operator.get_matrix_free());
  for (unsigned int cell = 0;
       cell < inhomogeneous_operator.get_matrix_free()->n_cell_batches();
       ++cell)
    {
      phi.reinit(cell);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_value(make_vectorized_array<double>(1.0), q);
      phi.integrate(EvaluationFlags::values);
      phi.distribute_local_to_global(system_rhs);
    }
  system_rhs.compress(VectorOperation::add);
}
@endcode

A more sophisticated implementation of this technique could reuse the original
MatrixFree object. This can be done by initializing the MatrixFree object with
multiple blocks, where each block corresponds to a different AffineConstraints
object. Doing this would require making substantial modifications to the
LaplaceOperator class, but the MatrixFreeOperators::LaplaceOperator class that
comes with the library can do this. See the discussion on blocks in
MatrixFreeOperators::Base for more information on how to set up blocks.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-37.cc"
*/
