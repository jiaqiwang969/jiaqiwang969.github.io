/**
@page step_15 The step-15 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Foreword">Foreword</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Weakformulationoftheproblem">Weak formulation of the problem</a>
        <li><a href="#Questionsabouttheappropriatesolver"> Questions about the appropriate solver </a>
        <li><a href="#Choiceofsteplengthandglobalization"> Choice of step length and globalization </a>
        <li><a href="#Summaryofthealgorithmandtestcase"> Summary of the algorithm and testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a>
        <li><a href="#Boundarycondition">Boundary condition</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#MinimalSurfaceProblemMinimalSurfaceProblem">MinimalSurfaceProblem::MinimalSurfaceProblem</a>
        <li><a href="#MinimalSurfaceProblemsetup_system">MinimalSurfaceProblem::setup_system</a>
        <li><a href="#MinimalSurfaceProblemassemble_system">MinimalSurfaceProblem::assemble_system</a>
        <li><a href="#MinimalSurfaceProblemsolve">MinimalSurfaceProblem::solve</a>
        <li><a href="#MinimalSurfaceProblemrefine_mesh">MinimalSurfaceProblem::refine_mesh</a>
        <li><a href="#MinimalSurfaceProblemset_boundary_values">MinimalSurfaceProblem::set_boundary_values</a>
        <li><a href="#MinimalSurfaceProblemcompute_residual">MinimalSurfaceProblem::compute_residual</a>
        <li><a href="#MinimalSurfaceProblemdetermine_step_length">MinimalSurfaceProblem::determine_step_length</a>
        <li><a href="#MinimalSurfaceProblemoutput_results">MinimalSurfaceProblem::output_results</a>
        <li><a href="#MinimalSurfaceProblemrun">MinimalSurfaceProblem::run</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Steplengthcontrol"> Step length control </a>
        <li><a href="#Integratingmeshrefinementandnonlinearandlinearsolvers"> Integrating mesh refinement and nonlinear and linear solvers </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program grew out of a student project by Sven Wetterauer at the
University of Heidelberg, Germany. Most of the work for this program
is by him.
</i>
<br>


<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


<a name="Foreword"></a><h3>Foreword</h3>


This program deals with an example of a non-linear elliptic partial
differential equation, the minimal
surface equation. You can imagine the solution of this equation to describe
the surface spanned by a soap film that is enclosed by a
closed wire loop. We imagine the wire to not just be a planar loop, but in
fact curved. The surface tension of the soap film will then reduce the surface
to have minimal surface. The solution of the minimal surface equation
describes this shape with the wire's vertical displacement as a boundary
condition. For simplicity, we will here assume that the surface can be written
as a graph $u=u(x,y)$ although it is clear that it is not very hard to
construct cases where the wire is bent in such a way that the surface can only
locally be constructed as a graph but not globally.

Because the equation is non-linear, we can't solve it directly. Rather, we
have to use Newton's method to compute the solution iteratively.

@dealiiVideoLecture{31.5,31.55,31.6}
(@dealiiVideoLectureSeeAlso{31.65,31.7})



<a name="Classicalformulation"></a><h3>Classical formulation</h3>


In a classical sense, the problem is given in the following form:


  @f{align*}
    -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.
  @f}

$\Omega$ is the domain we get by projecting the wire's positions into $x-y$
space. In this example, we choose $\Omega$ as the unit disk.

As described above, we solve this equation using Newton's method in which we
compute the $n$th approximate solution from the $n$th$-1$ one, and use
a damping parameter $\alpha^n$ to get better global convergence behavior:
  @f{align*}
    F'(u^{n},\delta u^{n})&=- F(u^{n})
    \\
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}
with
  @f[
    F(u) \dealcoloneq -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)
  @f]
and $F'(u,\delta u)$ the derivative of F in direction of $\delta u$:
@f[
  F'(u,\delta u)=\lim \limits_{\epsilon \rightarrow 0}{\frac{F(u+\epsilon \delta u)-
  F(u)}{\epsilon}}.
@f]

Going through the motions to find out what $F'(u,\delta u)$ is, we find that
we have to solve a linear elliptic PDE in every Newton step, with $\delta u^n$
as the solution of:

  @f[
  - \nabla \cdot \left( \frac{1}{(1+|\nabla u^{n}|^{2})^{\frac{1}{2}}}\nabla
  \delta u^{n} \right) +
  \nabla \cdot \left( \frac{\nabla u^{n} \cdot
  \nabla \delta u^{n}}{(1+|\nabla u^{n}|^{2})^{\frac{3}{2}}} \nabla u^{n}
  \right)  =
  -\left( - \nabla \cdot \left( \frac{1}{(1+|\nabla u^{n}|^{2})^{\frac{1}{2}}}
  \nabla u^{n} \right) \right)
  @f]

In order to solve the minimal surface equation, we have to solve this equation
repeatedly, once per Newton step. To solve this, we have to take a look at the
boundary condition of this problem. Assuming that $u^{n}$ already has the
right boundary values, the Newton update $\delta u^{n}$ should have zero
boundary conditions, in order to have the right boundary condition after
adding both.  In the first Newton step, we are starting with the solution
$u^{0}\equiv 0$, the Newton update still has to deliver the right boundary
condition to the solution $u^{1}$.


Summing up, we have to solve the PDE above with the boundary condition $\delta
u^{0}=g$ in the first step and with $\delta u^{n}=0$ in all the following steps.


<a name="Weakformulationoftheproblem"></a><h3>Weak formulation of the problem</h3>


Starting with the strong formulation above, we get the weak formulation by multiplying
both sides of the PDE with a test function $\varphi$ and integrating by parts on both sides:
  @f[
  \left( \nabla \varphi , \frac{1}{(1+|\nabla u^{n}|^{2})^{\frac{1}{2}}}\nabla
  \delta u^{n} \right)-\left(\nabla \varphi ,\frac{\nabla u^{n} \cdot \nabla
  \delta u^{n}}{(1+|\nabla u^{n}|^{2})^{\frac{3}{2}}}\nabla u^{n}  \right)
  = -\left(\nabla \varphi , \frac{1}{(1+|\nabla u^{n}|^{2})^{\frac{1}{2}}} \nabla u^{n}
   \right).
  @f]
Here the solution $\delta u^{n}$ is a function in $H^{1}(\Omega)$, subject to
the boundary conditions discussed above.
Reducing this space to a finite dimensional space with basis $\left\{
\varphi_{0},\dots , \varphi_{N-1}\right\}$, we can write the solution:

@f[
  \delta u^{n}=\sum_{j=0}^{N-1} \delta U_{j} \varphi_{j}.
@f]

Using the basis functions as test functions and defining $a_{n} \dealcoloneq \frac{1}
{\sqrt{1+|\nabla u^{n}|^{2}}}$, we can rewrite the weak formulation:

@f[
  \sum_{j=0}^{N-1}\left[ \left( \nabla \varphi_{i} , a_{n} \nabla \varphi_{j} \right) -
  \left(\nabla u^{n}\cdot \nabla \varphi_{i} , a_{n}^{3} \nabla u^{n} \cdot \nabla
  \varphi_{j} \right) \right] \cdot \delta U_{j}=-\left( \nabla \varphi_{i} , a_{n}
  \nabla u^{n}\right) \qquad \forall i=0,\dots ,N-1,
@f]

where the solution $\delta u^{n}$ is given by the coefficients $\delta U^{n}_{j}$.
This linear system of equations can be rewritten as:

@f[
  A^{n}\; \delta U^{n}=b^{n},
@f]

where the entries of the matrix $A^{n}$ are given by:

@f[
  A^{n}_{ij} \dealcoloneq \left( \nabla \varphi_{i} , a_{n} \nabla \varphi_{j} \right) -
  \left(\nabla u^{n}\cdot \nabla \varphi_{i} , a_{n}^{3} \nabla u^{n} \cdot \nabla
  \varphi_{j} \right),
@f]

and the right hand side $b^{n}$ is given by:

@f[
  b^{n}_{i} \dealcoloneq -\left( \nabla \varphi_{i} , a_{n} \nabla u^{n}\right).
@f]


<a name="Questionsabouttheappropriatesolver"></a><h3> Questions about the appropriate solver </h3>


The matrix that corresponds to the Newton step above can be reformulated to
show its structure a bit better. Rewriting it slightly, we get that it has the
form
@f[
  A_{ij}
  =
  \left(
    \nabla \varphi_i,
    B
    \nabla \varphi_j
  \right),
@f]
where the matrix $B$ (of size $d \times d$ in $d$ space dimensions) is given
by the following expression:
@f[
  B
  =
  a_n \left\{
   \mathbf I
   -
   a_n^2 [\nabla u_n] \otimes [\nabla u_n]
  \right\}
  =
  a_n \left\{
   \mathbf I
   -
  \frac{\nabla u_n}{\sqrt{1+|\nabla u^{n}|^{2}}} \otimes
  \frac{\nabla u_n}{\sqrt{1+|\nabla u^{n}|^{2}}}
  \right\}.
@f]
From this expression, it is obvious that
$B$ is symmetric, and so $A$ is symmetric as well.
On the other hand, $B$ is also positive definite, which confers the same
property onto $A$. This can be seen by noting that the vector $v_1 =
\frac{\nabla u^n}{|\nabla u^n|}$ is an eigenvector of $B$ with eigenvalue
$\lambda_1=a_n \left(1-\frac{|\nabla u^n|^2}{1+|\nabla u^n|^2}\right) > 0$ while all vectors $v_2\ldots v_d$
that are perpendicular to $v_1$ and each other are eigenvectors with
eigenvalue $a_n$. Since all eigenvalues are positive, $B$ is positive definite
and so is $A$. We can thus use the CG method for solving the Newton steps.
(The fact that the matrix $A$ is symmetric and positive definite should not come
as a surprise. It results from taking the derivative of an operator that
results from taking the derivative of an energy functional: the minimal
surface equation simply minimizes some non-quadratic energy. Consequently,
the Newton matrix, as the matrix of second derivatives of a scalar energy,
must be symmetric since the derivative with regard to the $i$th and $j$th
degree of freedom should clearly commute. Likewise, if the energy functional
is convex, then the matrix of second derivatives must be positive definite,
and the direct calculation above simply reaffirms this.)

It is worth noting, however, that the positive definiteness degenerates for
problems where $\nabla u$ becomes large. In other words, if we simply multiply
all boundary values by 2, then to first order $u$ and $\nabla u$ will also be
multiplied by two, but as a consequence the smallest eigenvalue of $B$ will
become smaller and the matrix will become more ill-conditioned. (More
specifically, for $|\nabla u^n|\rightarrow\infty$ we have that
$\lambda_1 \propto a_n \frac{1}{|\nabla u^n|^2}$ whereas
$\lambda_2\ldots \lambda_d=a_n$; thus, the condition number of $B$,
which is a multiplicative factor in the condition number of $A$ grows
like ${\cal O}(|\nabla u^n|^2)$.) It is simple
to verify with the current program that indeed multiplying the boundary values
used in the current program by larger and larger values results in a problem
that will ultimately no longer be solvable using the simple preconditioned CG
method we use here.


<a name="Choiceofsteplengthandglobalization"></a><h3> Choice of step length and globalization </h3>


As stated above, Newton's method works by computing a direction
$\delta u^n$ and then performing the update $u^{n+1} = u^{n}+\alpha^n
\delta u^{n}$ with a step length $0 < \alpha^n \le 1$. It is a common
observation that for strongly nonlinear models, Newton's method does
not converge if we always choose $\alpha^n=1$ unless one starts with
an initial guess $u^0$ that is sufficiently close to the solution $u$
of the nonlinear problem. In practice, we don't always have such an
initial guess, and consequently taking full Newton steps (i.e., using
$\alpha=1$) does frequently not work.

A common strategy therefore is to use a smaller step length for the
first few steps while the iterate $u^n$ is still far away from the
solution $u$ and as we get closer use larger values for $\alpha^n$
until we can finally start to use full steps $\alpha^n=1$ as we are
close enough to the solution. The question is of course how to choose
$\alpha^n$. There are basically two widely used approaches: line
search and trust region methods.

In this program, we simply always choose the step length equal to
0.1. This makes sure that for the testcase at hand we do get
convergence although it is clear that by not eventually reverting to
full step lengths we forego the rapid, quadratic convergence that
makes Newton's method so appealing. Obviously, this is a point one
eventually has to address if the program was made into one that is
meant to solve more realistic problems. We will comment on this issue
some more in the <a href="#Results">results section</a>.


<a name="Summaryofthealgorithmandtestcase"></a><h3> Summary of the algorithm and testcase </h3>


Overall, the program we have here is not unlike step-6 in many regards. The
layout of the main class is essentially the same. On the other hand, the
driving algorithm in the <code>run()</code> function is different and works as
follows:
<ol>
<li>
  Start with the function $u^{0}\equiv 0$ and modify it in such a way
  that the values of $u^0$ along the boundary equal the correct
  boundary values $g$ (this happens in
  <code>MinimalSurfaceProblem::set_boundary_values</code>). Set
  $n=0$.
</li>

<li>
  Compute the Newton update by solving the system $A^{n}\;\delta
  U^{n}=b^{n}$
  with boundary condition $\delta u^{n}=0$ on $\partial \Omega$.
</li>

<li>
  Compute a step length $\alpha^n$. In this program, we always set
  $\alpha^n=0.1$. To make things easier to extend later on, this
  happens in a function of its own, namely in
  <code>MinimalSurfaceProblem::determine_step_length</code>.
</li>

<li>
  The new approximation of the solution is given by
  $u^{n+1}=u^{n}+\alpha^n \delta u^{n}$.
</li>

<li>
  If $n$ is a multiple of 5 then refine the mesh, transfer the
  solution $u^{n+1}$ to the new mesh and set the values of $u^{n+1}$
  in such a way that along the boundary we have
  $u^{n+1}|_{\partial\Gamma}=g$ (again in
  <code>MinimalSurfaceProblem::set_boundary_values</code>). Note that
  this isn't automatically
  guaranteed even though by construction we had that before mesh
  refinement $u^{n+1}|_{\partial\Gamma}=g$ because mesh refinement
  adds new nodes to the mesh where we have to interpolate the old
  solution to the new nodes upon bringing the solution from the old to
  the new mesh. The values we choose by interpolation may be close to
  the exact boundary conditions but are, in general, nonetheless not
  the correct values.
</li>

<li>
  Set $n\leftarrow n+1$ and go to step 2.
</li>
</ol>

The testcase we solve is chosen as follows: We seek to find the solution of
minimal surface over the unit disk $\Omega=\{\mathbf x: \|\mathbf
x\|<1\}\subset {\mathbb R}^2$ where the surface attains the values
$u(x,y)|{\partial\Omega} = g(x,y) \dealcoloneq \sin(2 \pi (x+y))$ along the
boundary.
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
 * The first few files have already been covered in previous examples and will
 * thus not be further commented on.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * We will use adaptive mesh refinement between Newton iterations. To do so,
 * we need to be able to work with a solution on the new mesh, although it was
 * computed on the old one. The SolutionTransfer class transfers the solution
 * from the old to the new mesh:
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * @endcode
 * 
 * We then open a namespace for this program and import everything from the
 * dealii namespace into it, as in previous programs:
 * 
 * @code
 * namespace Step15
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class template</h3>
 * 

 * 
 * The class template is basically the same as in step-6.  Three additions
 * are made:
 * - There are two solution vectors, one for the Newton update
 * $\delta u^n$, and one for the current iterate $u^n$.
 * - The <code>setup_system</code> function takes an argument that denotes
 * whether this is the first time it is called or not. The difference is
 * that the first time around we need to distribute the degrees of freedom
 * and set the solution vector for $u^n$ to the correct size. The following
 * times, the function is called after we have already done these steps as
 * part of refining the mesh in <code>refine_mesh</code>.
 * - We then also need new functions: <code>set_boundary_values()</code>
 * takes care of setting the boundary values on the solution vector
 * correctly, as discussed at the end of the
 * introduction. <code>compute_residual()</code> is a function that computes
 * the norm of the nonlinear (discrete) residual. We use this function to
 * monitor convergence of the Newton iteration. The function takes a step
 * length $\alpha^n$ as argument to compute the residual of $u^n + \alpha^n
 * \; \delta u^n$. This is something one typically needs for step length
 * control, although we will not use this feature here. Finally,
 * <code>determine_step_length()</code> computes the step length $\alpha^n$
 * in each Newton iteration. As discussed in the introduction, we here use a
 * fixed step length and leave implementing a better strategy as an
 * exercise.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class MinimalSurfaceProblem
 *   {
 *   public:
 *     MinimalSurfaceProblem();
 *     void run();
 * 
 *   private:
 *     void   setup_system(const bool initial_step);
 *     void   assemble_system();
 *     void   solve();
 *     void   refine_mesh();
 *     void   set_boundary_values();
 *     double compute_residual(const double alpha) const;
 *     double determine_step_length() const;
 *     void   output_results(const unsigned int refinement_cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     DoFHandler<dim> dof_handler;
 *     FE_Q<dim>       fe;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> current_solution;
 *     Vector<double> newton_update;
 *     Vector<double> system_rhs;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="Boundarycondition"></a> 
 * <h3>Boundary condition</h3>
 * 

 * 
 * The boundary condition is implemented just like in step-4.  It is chosen
 * as $g(x,y)=\sin(2 \pi (x+y))$:
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> &p,
 *                                     const unsigned int /*component*/) const
 *   {
 *     return std::sin(2 * numbers::PI * (p[0] + p[1]));
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="MinimalSurfaceProblemMinimalSurfaceProblem"></a> 
 * <h4>MinimalSurfaceProblem::MinimalSurfaceProblem</h4>
 * 

 * 
 * The constructor and destructor of the class are the same as in the first
 * few tutorials.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
 *     : dof_handler(triangulation)
 *     , fe(2)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemsetup_system"></a> 
 * <h4>MinimalSurfaceProblem::setup_system</h4>
 * 

 * 
 * As always in the setup-system function, we setup the variables of the
 * finite element method. There are same differences to step-6, because
 * there we start solving the PDE from scratch in every refinement cycle
 * whereas here we need to take the solution from the previous mesh onto the
 * current mesh. Consequently, we can't just reset solution vectors. The
 * argument passed to this function thus indicates whether we can
 * distributed degrees of freedom (plus compute constraints) and set the
 * solution vector to zero or whether this has happened elsewhere already
 * (specifically, in <code>refine_mesh()</code>).
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
 *   {
 *     if (initial_step)
 *       {
 *         dof_handler.distribute_dofs(fe);
 *         current_solution.reinit(dof_handler.n_dofs());
 * 
 *         hanging_node_constraints.clear();
 *         DoFTools::make_hanging_node_constraints(dof_handler,
 *                                                 hanging_node_constraints);
 *         hanging_node_constraints.close();
 *       }
 * 
 * 
 * @endcode
 * 
 * The remaining parts of the function are the same as in step-6.
 * 

 * 
 * 
 * @code
 *     newton_update.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 *     hanging_node_constraints.condense(dsp);
 * 
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemassemble_system"></a> 
 * <h4>MinimalSurfaceProblem::assemble_system</h4>
 * 

 * 
 * This function does the same as in the previous tutorials except that now,
 * of course, the matrix and right hand side functions depend on the
 * previous iteration's solution. As discussed in the introduction, we need
 * to use zero boundary values for the Newton updates; we compute them at
 * the end of this function.
 *   

 * 
 * The top of the function contains the usual boilerplate code, setting up
 * the objects that allow us to evaluate shape functions at quadrature
 * points and temporary storage locations for the local matrices and
 * vectors, as well as for the gradients of the previous solution at the
 * quadrature points. We then start the loop over all cells:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::assemble_system()
 *   {
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_gradients | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);
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
 * @endcode
 * 
 * For the assembly of the linear system, we have to obtain the values
 * of the previous solution's gradients at the quadrature
 * points. There is a standard way of doing this: the
 * FEValues::get_function_gradients function takes a vector that
 * represents a finite element field defined on a DoFHandler, and
 * evaluates the gradients of this field at the quadrature points of the
 * cell with which the FEValues object has last been reinitialized.
 * The values of the gradients at all quadrature points are then written
 * into the second argument:
 * 
 * @code
 *         fe_values.get_function_gradients(current_solution,
 *                                          old_solution_gradients);
 * 
 * @endcode
 * 
 * With this, we can then do the integration loop over all quadrature
 * points and shape functions.  Having just computed the gradients of
 * the old solution in the quadrature points, we are able to compute
 * the coefficients $a_{n}$ in these points.  The assembly of the
 * system itself then looks similar to what we always do with the
 * exception of the nonlinear terms, as does copying the results from
 * the local objects into the global ones:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double coeff =
 *               1.0 / std::sqrt(1 + old_solution_gradients[q] *
 *                                     old_solution_gradients[q]);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                   cell_matrix(i, j) +=
 *                     (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
 *                        * coeff                         //   * a_n
 *                        * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
 *                       -                                //  -
 *                       (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
 *                        * coeff * coeff * coeff         //   * a_n^3
 *                        * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
 *                           * old_solution_gradients[q]) //      * \nabla u_n)
 *                        * old_solution_gradients[q]))   //   * \nabla u_n)))
 *                      * fe_values.JxW(q));              // * dx
 * 
 *                 cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
 *                                 * coeff                     // * a_n
 *                                 * old_solution_gradients[q] // * u_n
 *                                 * fe_values.JxW(q));        // * dx
 *               }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               system_matrix.add(local_dof_indices[i],
 *                                 local_dof_indices[j],
 *                                 cell_matrix(i, j));
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *           }
 *       }
 * 
 * @endcode
 * 
 * Finally, we remove hanging nodes from the system and apply zero
 * boundary values to the linear system that defines the Newton updates
 * $\delta u^n$:
 * 
 * @code
 *     hanging_node_constraints.condense(system_matrix);
 *     hanging_node_constraints.condense(system_rhs);
 * 
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(),
 *                                              boundary_values);
 *     MatrixTools::apply_boundary_values(boundary_values,
 *                                        system_matrix,
 *                                        newton_update,
 *                                        system_rhs);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemsolve"></a> 
 * <h4>MinimalSurfaceProblem::solve</h4>
 * 

 * 
 * The solve function is the same as always. At the end of the solution
 * process we update the current solution by setting
 * $u^{n+1}=u^n+\alpha^n\;\delta u^n$.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(system_rhs.size(),
 *                                  system_rhs.l2_norm() * 1e-6);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
 * 
 *     hanging_node_constraints.distribute(newton_update);
 * 
 *     const double alpha = determine_step_length();
 *     current_solution.add(alpha, newton_update);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemrefine_mesh"></a> 
 * <h4>MinimalSurfaceProblem::refine_mesh</h4>
 * 

 * 
 * The first part of this function is the same as in step-6... However,
 * after refining the mesh we have to transfer the old solution to the new
 * one which we do with the help of the SolutionTransfer class. The process
 * is slightly convoluted, so let us describe it in detail:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::refine_mesh()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       current_solution,
 *       estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 * 
 * @endcode
 * 
 * Then we need an additional step: if, for example, you flag a cell that
 * is once more refined than its neighbor, and that neighbor is not
 * flagged for refinement, we would end up with a jump of two refinement
 * levels across a cell interface.  To avoid these situations, the library
 * will silently also have to refine the neighbor cell once. It does so by
 * calling the Triangulation::prepare_coarsening_and_refinement function
 * before actually doing the refinement and coarsening.  This function
 * flags a set of additional cells for refinement or coarsening, to
 * enforce rules like the one-hanging-node rule.  The cells that are
 * flagged for refinement and coarsening after calling this function are
 * exactly the ones that will actually be refined or coarsened. Usually,
 * you don't have to do this by hand
 * (Triangulation::execute_coarsening_and_refinement does this for
 * you). However, we need to initialize the SolutionTransfer class and it
 * needs to know the final set of cells that will be coarsened or refined
 * in order to store the data from the old mesh and transfer to the new
 * one. Thus, we call the function by hand:
 * 
 * @code
 *     triangulation.prepare_coarsening_and_refinement();
 * 
 * @endcode
 * 
 * With this out of the way, we initialize a SolutionTransfer object with
 * the present DoFHandler and attach the solution vector to it, followed
 * by doing the actual refinement and distribution of degrees of freedom
 * on the new mesh
 * 
 * @code
 *     SolutionTransfer<dim> solution_transfer(dof_handler);
 *     solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * Finally, we retrieve the old solution interpolated to the new
 * mesh. Since the SolutionTransfer function does not actually store the
 * values of the old solution, but rather indices, we need to preserve the
 * old solution vector until we have gotten the new interpolated
 * values. Thus, we have the new values written into a temporary vector,
 * and only afterwards write them into the solution vector object. Once we
 * have this solution we have to make sure that the $u^n$ we now have
 * actually has the correct boundary values. As explained at the end of
 * the introduction, this is not automatically the case even if the
 * solution before refinement had the correct boundary values, and so we
 * have to explicitly make sure that it now has:
 * 
 * @code
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     solution_transfer.interpolate(current_solution, tmp);
 *     current_solution = tmp;
 * 
 *     set_boundary_values();
 * 
 * @endcode
 * 
 * On the new mesh, there are different hanging nodes, which we have to
 * compute again. To ensure there are no hanging nodes of the old mesh in
 * the object, it's first cleared.  To be on the safe side, we then also
 * make sure that the current solution's vector entries satisfy the
 * hanging node constraints (see the discussion in the documentation of
 * the SolutionTransfer class for why this is necessary):
 * 
 * @code
 *     hanging_node_constraints.clear();
 * 
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 *     hanging_node_constraints.distribute(current_solution);
 * 
 * @endcode
 * 
 * We end the function by updating all the remaining data structures,
 * indicating to <code>setup_dofs()</code> that this is not the first
 * go-around and that it needs to preserve the content of the solution
 * vector:
 * 
 * @code
 *     setup_system(false);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemset_boundary_values"></a> 
 * <h4>MinimalSurfaceProblem::set_boundary_values</h4>
 * 

 * 
 * The next function ensures that the solution vector's entries respect the
 * boundary values for our problem.  Having refined the mesh (or just
 * started computations), there might be new nodal points on the
 * boundary. These have values that are simply interpolated from the
 * previous mesh (or are just zero), instead of the correct boundary
 * values. This is fixed up by setting all boundary nodes explicit to the
 * right value:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::set_boundary_values()
 *   {
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              BoundaryValues<dim>(),
 *                                              boundary_values);
 *     for (auto &boundary_value : boundary_values)
 *       current_solution(boundary_value.first) = boundary_value.second;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemcompute_residual"></a> 
 * <h4>MinimalSurfaceProblem::compute_residual</h4>
 * 

 * 
 * In order to monitor convergence, we need a way to compute the norm of the
 * (discrete) residual, i.e., the norm of the vector
 * $\left<F(u^n),\varphi_i\right>$ with $F(u)=-\nabla \cdot \left(
 * \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)$ as discussed in the
 * introduction. It turns out that (although we don't use this feature in
 * the current version of the program) one needs to compute the residual
 * $\left<F(u^n+\alpha^n\;\delta u^n),\varphi_i\right>$ when determining
 * optimal step lengths, and so this is what we implement here: the function
 * takes the step length $\alpha^n$ as an argument. The original
 * functionality is of course obtained by passing a zero as argument.
 *   

 * 
 * In the function below, we first set up a vector for the residual, and
 * then a vector for the evaluation point $u^n+\alpha^n\;\delta u^n$. This
 * is followed by the same boilerplate code we use for all integration
 * operations:
 * 
 * @code
 *   template <int dim>
 *   double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
 *   {
 *     Vector<double> residual(dof_handler.n_dofs());
 * 
 *     Vector<double> evaluation_point(dof_handler.n_dofs());
 *     evaluation_point = current_solution;
 *     evaluation_point.add(alpha, newton_update);
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_gradients | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double>              cell_residual(dofs_per_cell);
 *     std::vector<Tensor<1, dim>> gradients(n_q_points);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_residual = 0;
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * The actual computation is much as in
 * <code>assemble_system()</code>. We first evaluate the gradients of
 * $u^n+\alpha^n\,\delta u^n$ at the quadrature points, then compute
 * the coefficient $a_n$, and then plug it all into the formula for
 * the residual:
 * 
 * @code
 *         fe_values.get_function_gradients(evaluation_point, gradients);
 * 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double coeff =
 *               1. / std::sqrt(1 + gradients[q] * gradients[q]);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
 *                                    * coeff                    // * a_n
 *                                    * gradients[q]             // * u_n
 *                                    * fe_values.JxW(q));       // * dx
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           residual(local_dof_indices[i]) += cell_residual(i);
 *       }
 * 
 * @endcode
 * 
 * At the end of this function we also have to deal with the hanging node
 * constraints and with the issue of boundary values. With regard to the
 * latter, we have to set to zero the elements of the residual vector for
 * all entries that correspond to degrees of freedom that sit at the
 * boundary. The reason is that because the value of the solution there is
 * fixed, they are of course no "real" degrees of freedom and so, strictly
 * speaking, we shouldn't have assembled entries in the residual vector
 * for them. However, as we always do, we want to do exactly the same
 * thing on every cell and so we didn't not want to deal with the question
 * of whether a particular degree of freedom sits at the boundary in the
 * integration above. Rather, we will simply set to zero these entries
 * after the fact. To this end, we need to determine which degrees
 * of freedom do in fact belong to the boundary and then loop over all of
 * those and set the residual entry to zero. This happens in the following
 * lines which we have already seen used in step-11, using the appropriate
 * function from namespace DoFTools:
 * 
 * @code
 *     hanging_node_constraints.condense(residual);
 * 
 *     for (types::global_dof_index i :
 *          DoFTools::extract_boundary_dofs(dof_handler))
 *       residual(i) = 0;
 * 
 * @endcode
 * 
 * At the end of the function, we return the norm of the residual:
 * 
 * @code
 *     return residual.l2_norm();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemdetermine_step_length"></a> 
 * <h4>MinimalSurfaceProblem::determine_step_length</h4>
 * 

 * 
 * As discussed in the introduction, Newton's method frequently does not
 * converge if we always take full steps, i.e., compute $u^{n+1}=u^n+\delta
 * u^n$. Rather, one needs a damping parameter (step length) $\alpha^n$ and
 * set $u^{n+1}=u^n+\alpha^n\delta u^n$. This function is the one called
 * to compute $\alpha^n$.
 *   

 * 
 * Here, we simply always return 0.1. This is of course a sub-optimal
 * choice: ideally, what one wants is that the step size goes to one as we
 * get closer to the solution, so that we get to enjoy the rapid quadratic
 * convergence of Newton's method. We will discuss better strategies below
 * in the results section.
 * 
 * @code
 *   template <int dim>
 *   double MinimalSurfaceProblem<dim>::determine_step_length() const
 *   {
 *     return 0.1;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemoutput_results"></a> 
 * <h4>MinimalSurfaceProblem::output_results</h4>
 * 

 * 
 * This last function to be called from `run()` outputs the current solution
 * (and the Newton update) in graphical form as a VTU file. It is entirely the
 * same as what has been used in previous tutorials.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::output_results(
 *     const unsigned int refinement_cycle) const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(current_solution, "solution");
 *     data_out.add_data_vector(newton_update, "update");
 *     data_out.build_patches();
 * 
 *     const std::string filename =
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemrun"></a> 
 * <h4>MinimalSurfaceProblem::run</h4>
 * 

 * 
 * In the run function, we build the first grid and then have the top-level
 * logic for the Newton iteration.
 *   

 * 
 * As described in the introduction, the domain is the unit disk around
 * the origin, created in the same way as shown in step-6. The mesh is
 * globally refined twice followed later on by several adaptive cycles.
 *   

 * 
 * Before starting the Newton loop, we also need to do a bit of
 * setup work: We need to create the basic data structures and
 * ensure that the first Newton iterate already has the correct
 * boundary values, as discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_ball(triangulation);
 *     triangulation.refine_global(2);
 * 
 *     setup_system(/*first time=*/true);
 *     set_boundary_values();
 * 
 * @endcode
 * 
 * The Newton iteration starts next. We iterate until the (norm of the)
 * residual computed at the end of the previous iteration is less than
 * $10^{-3}$, as checked at the end of the `do { ... } while` loop that
 * starts here. Because we don't have a reasonable value to initialize
 * the variable, we just use the largest value that can be represented
 * as a `double`.
 * 
 * @code
 *     double       last_residual_norm = std::numeric_limits<double>::max();
 *     unsigned int refinement_cycle   = 0;
 *     do
 *       {
 *         std::cout << "Mesh refinement step " << refinement_cycle << std::endl;
 * 
 *         if (refinement_cycle != 0)
 *           refine_mesh();
 * 
 * @endcode
 * 
 * On every mesh we do exactly five Newton steps. We print the initial
 * residual here and then start the iterations on this mesh.
 *         

 * 
 * In every Newton step the system matrix and the right hand side have
 * to be computed first, after which we store the norm of the right
 * hand side as the residual to check against when deciding whether to
 * stop the iterations. We then solve the linear system (the function
 * also updates $u^{n+1}=u^n+\alpha^n\;\delta u^n$) and output the
 * norm of the residual at the end of this Newton step.
 *         

 * 
 * After the end of this loop, we then also output the solution on the
 * current mesh in graphical form and increment the counter for the
 * mesh refinement cycle.
 * 
 * @code
 *         std::cout << "  Initial residual: " << compute_residual(0) << std::endl;
 * 
 *         for (unsigned int inner_iteration = 0; inner_iteration < 5;
 *              ++inner_iteration)
 *           {
 *             assemble_system();
 *             last_residual_norm = system_rhs.l2_norm();
 * 
 *             solve();
 * 
 *             std::cout << "  Residual: " << compute_residual(0) << std::endl;
 *           }
 * 
 *         output_results(refinement_cycle);
 * 
 *         ++refinement_cycle;
 *         std::cout << std::endl;
 *       }
 *     while (last_residual_norm > 1e-3);
 *   }
 * } // namespace Step15
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h4>The main function</h4>
 * 

 * 
 * Finally the main function. This follows the scheme of all other main
 * functions:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step15;
 * 
 *       MinimalSurfaceProblem<2> laplace_problem_2d;
 *       laplace_problem_2d.run();
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



The output of the program looks as follows:
@code
Mesh refinement step
  Initial residual: 1.53143
  Residual: 1.08746
  Residual: 0.966748
  Residual: 0.859602
  Residual: 0.766462
  Residual: 0.685475

Mesh refinement step
  Initial residual: 0.868959
  Residual: 0.762125
  Residual: 0.677792
  Residual: 0.605762
  Residual: 0.542748
  Residual: 0.48704

Mesh refinement step
  Initial residual: 0.426445
  Residual: 0.382731
  Residual: 0.343865
  Residual: 0.30918
  Residual: 0.278147
  Residual: 0.250327

Mesh refinement step
  Initial residual: 0.282026
  Residual: 0.253146
  Residual: 0.227414
  Residual: 0.20441
  Residual: 0.183803
  Residual: 0.165319

Mesh refinement step
  Initial residual: 0.154404
  Residual: 0.138723
  Residual: 0.124694
  Residual: 0.112124
  Residual: 0.100847
  Residual: 0.0907222

....
@endcode

Obviously, the scheme converges, if not very fast. We will come back to
strategies for accelerating the method below.

One can visualize the solution after each set of five Newton
iterations, i.e., on each of the meshes on which we approximate the
solution. This yields the following set of images:

<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_15_solution_1.png"
         alt="Solution after zero cycles with contour lines." width="230" height="273">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_15_solution_2.png"
         alt="Solution after one cycle with contour lines." width="230" height="273">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_15_solution_3.png"
         alt="Solution after two cycles with contour lines." width="230" height="273">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_15_solution_4.png"
         alt="Solution after three cycles with contour lines." width="230" height="273">
  </div>
</div>

It is clearly visible, that the solution minimizes the surface
after each refinement. The solution converges to a picture one
would imagine a soap bubble to be that is located inside a wire loop
that is bent like
the boundary. Also it is visible, how the boundary
is smoothed out after each refinement. On the coarse mesh,
the boundary doesn't look like a sine, whereas it does the
finer the mesh gets.

The mesh is mostly refined near the boundary, where the solution
increases or decreases strongly, whereas it is coarsened on
the inside of the domain, where nothing interesting happens,
because there isn't much change in the solution. The ninth
solution and mesh are shown here:

<div class="onecolumn" style="width: 60%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_15_solution_9.png"
         alt="Grid and solution of the ninth cycle with contour lines." width="507" height="507">
  </div>
</div>



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


The program shows the basic structure of a solver for a nonlinear, stationary
problem. However, it does not converge particularly fast, for good reasons:

- The program always takes a step size of 0.1. This precludes the rapid,
  quadratic convergence for which Newton's method is typically chosen.
- It does not connect the nonlinear iteration with the mesh refinement
  iteration.

Obviously, a better program would have to address these two points.
We will discuss them in the following.


<a name="Steplengthcontrol"></a><h4> Step length control </h4>


Newton's method has two well known properties:
- It may not converge from arbitrarily chosen starting points. Rather, a
  starting point has to be close enough to the solution to guarantee
  convergence. However, we can enlarge the area from which Newton's method
  converges by damping the iteration using a <i>step length</i> 0<$\alpha^n\le
  1$.
- It exhibits rapid convergence of quadratic order if (i) the step length is
  chosen as $\alpha^n=1$, and (ii) it does in fact converge with this choice
  of step length.

A consequence of these two observations is that a successful strategy is to
choose $\alpha^n<1$ for the initial iterations until the iterate has come
close enough to allow for convergence with full step length, at which point we
want to switch to $\alpha^n=1$. The question is how to choose $\alpha^n$ in an
automatic fashion that satisfies these criteria.

We do not want to review the literature on this topic here, but only briefly
mention that there are two fundamental approaches to the problem: backtracking
line search and trust region methods. The former is more widely used for
partial differential equations and essentially does the following:
- Compute a search direction
- See if the resulting residual of $u^n + \alpha^n\;\delta u^n$ with
  $\alpha^n=1$ is "substantially smaller" than that of $u^n$ alone.
- If so, then take $\alpha^n=1$.
- If not, try whether the residual is "substantially smaller" with
  $\alpha^n=2/3$.
- If so, then take $\alpha^n=2/3$.
- If not, try whether the residual is "substantially smaller" with
  $\alpha^n=(2/3)^2$.
- Etc.
One can of course choose other factors $r, r^2, \ldots$ than the $2/3,
(2/3)^2, \ldots$ chosen above, for $0<r<1$. It is obvious where the term
"backtracking" comes from: we try a long step, but if that doesn't work we try
a shorter step, and ever shorter step, etc. The function
<code>determine_step_length()</code> is written the way it is to support
exactly this kind of use case.

Whether we accept a particular step length $\alpha^n$ depends on how we define
"substantially smaller". There are a number of ways to do so, but without
going into detail let us just mention that the most common ones are to use the
Wolfe and Armijo-Goldstein conditions. For these, one can show the following:
- There is always a step length $\alpha^n$ for which the conditions are
  satisfied, i.e., the iteration never gets stuck as long as the problem is
  convex.
- If we are close enough to the solution, then the conditions allow for
  $\alpha^n=1$, thereby enabling quadratic convergence.

We will not dwell on this here any further but leave the implementation of
such algorithms as an exercise. We note, however, that when implemented
correctly then it is a common observation that most reasonably nonlinear
problems can be solved in anywhere between 5 and 15 Newton iterations to
engineering accuracy &mdash; substantially fewer than we need with the current
version of the program.

More details on globalization methods including backtracking can be found,
for example, in @cite GNS08 and @cite NW99.

A separate point, very much worthwhile making, however, is that in practice
the implementation of efficient nonlinear solvers is about as complicated as
the implementation of efficient finite element methods. One should not
attempt to reinvent the wheel by implementing all of the necessary steps
oneself. Rather, just like building finite element solvers on libraries
such as deal.II, one should be building nonlinear solvers on libraries such
as [SUNDIALS](https://computing.llnl.gov/projects/sundials). In fact,
deal.II has interfaces to SUNDIALS and in particular to its nonlinear solver
sub-package KINSOL through the SUNDIALS::KINSOL class. It would not be
very difficult to base the current problem on that interface.



<a name="Integratingmeshrefinementandnonlinearandlinearsolvers"></a><h4> Integrating mesh refinement and nonlinear and linear solvers </h4>


We currently do exactly 5 iterations on each mesh. But is this optimal? One
could ask the following questions:
- Maybe it is worthwhile doing more iterations on the initial meshes since
  there, computations are cheap.
- On the other hand, we do not want to do too many iterations on every mesh:
  yes, we could drive the residual to zero on every mesh, but that would only
  mean that the nonlinear iteration error is far smaller than the
  discretization error.
- Should we use solve the linear systems in each Newton step with higher or
  lower accuracy?

Ultimately, what this boils down to is that we somehow need to couple the
discretization error on the current mesh with the nonlinear residual we want
to achieve with the Newton iterations on a given mesh, and to the linear
iteration we want to achieve with the CG method within each Newton
iterations.

How to do this is, again, not entirely trivial, and we again leave it as a
future exercise.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-15.cc"
*/
