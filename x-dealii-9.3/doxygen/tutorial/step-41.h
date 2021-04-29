/**
@page step_41 The step-41 tutorial program
This tutorial depends on step-15.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Introduction">Introduction</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Derivationofthevariationalinequality">Derivation of the variational inequality</a>
        <li><a href="#Formulationasasaddlepointproblem">Formulation as a saddle point problem</a>
        <li><a href="#ActiveSetmethodstosolvethesaddlepointproblem">Active Set methods to solve the saddle point problem</a>
        <li><a href="#Theprimaldualactivesetalgorithm">The primal-dual active set algorithm</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeObstacleProblemcodeclasstemplate">The <code>ObstacleProblem</code> class template</a>
        <li><a href="#Righthandsideboundaryvaluesandtheobstacle">Right hand side, boundary values, and the obstacle</a>
        <li><a href="#ImplementationofthecodeObstacleProblemcodeclass">Implementation of the <code>ObstacleProblem</code> class</a>
      <ul>
        <li><a href="#ObstacleProblemObstacleProblem">ObstacleProblem::ObstacleProblem</a>
        <li><a href="#ObstacleProblemmake_grid">ObstacleProblem::make_grid</a>
        <li><a href="#ObstacleProblemsetup_system">ObstacleProblem::setup_system</a>
        <li><a href="#ObstacleProblemassemble_system">ObstacleProblem::assemble_system</a>
        <li><a href="#ObstacleProblemassemble_mass_matrix_diagonal">ObstacleProblem::assemble_mass_matrix_diagonal</a>
        <li><a href="#ObstacleProblemupdate_solution_and_constraints">ObstacleProblem::update_solution_and_constraints</a>
        <li><a href="#ObstacleProblemsolve">ObstacleProblem::solve</a>
        <li><a href="#ObstacleProblemoutput_results">ObstacleProblem::output_results</a>
        <li><a href="#ObstacleProblemrun">ObstacleProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>This program was contributed by Jörg Frohne (University of Siegen,
Germany) while on a long-term visit to Texas A&amp;M University.
<br>
This material is based upon work partly supported by ThyssenKrupp Steel Europe.
</i>


<a name="Intro"></a>
<a name="Introduction"></a><h3>Introduction</h3>


This example is based on the Laplace equation in 2d and deals with the
question what happens if a membrane is deflected by some external force but is
also constrained by an obstacle. In other words, think of a elastic membrane
clamped at the boundary to a rectangular frame (we choose $\Omega =
\left[-1,1\right]^2$) and that sags through due to gravity acting on it. What
happens now if there is an obstacle under the membrane that prevents it from
reaching its equilibrium position if gravity was the only existing force? In
the current example program, we will consider that under the membrane is a
stair step obstacle against which gravity pushes the membrane.

This problem is typically called the "obstacle problem" (see also <a
href="http://en.wikipedia.org/wiki/Obstacle_problem">this Wikipedia article</a>), and it results in a
variational inequality, rather than a variational equation when put into the
weak form. We will below derive it from the classical formulation, but before we
go on to discuss the mathematics let us show how the solution of the problem we
will consider in this tutorial program looks to gain some intuition of what
we should expect:

<table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.png" alt="">
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.png" alt="">
    </td>
  </tr>
</table>

Here, at the left, we see the displacement of the membrane. The shape
of the obstacle underneath is clearly visible. On the right, we overlay which
parts of the membrane are in contact with the obstacle. We will later call
this set of points the "active set" to indicate that an inequality constraint
is active there.


<a name="Classicalformulation"></a><h3>Classical formulation</h3>


The classical formulation of the problem possesses the following form:
@f{align*}
 -\textrm{div}\ \sigma &\geq f & &\quad\text{in } \Omega,\\
 \sigma &= \nabla u & &\quad\text{in } \Omega,\\
 u(\mathbf x) &= 0 & &\quad\text{on }\partial\Omega,\\
(-\Delta u - f)(u - g) &= 0 & &\quad\text{in } \Omega,\\
 u(\mathbf x) &\geq g(\mathbf x) & &\quad\text{in } \Omega
@f}
with $u\in H^2(\Omega)$.  $u$ is a scalar valued function that denotes the
vertical displacement of the membrane. The first equation is called equilibrium
condition with a force of areal density $f$. Here, we will consider this force
to be gravity. The second one is known as Hooke's Law that says that the stresses
$\sigma$ are proportional to the gradient of the displacements $u$ (the
proportionality constant, often denoted by $E$, has been set to one here,
without loss of generality; if it is constant, it can be put into the right
hand side function). At the boundary we have zero Dirichlet
conditions. Obviously, the first two equations can be combined to yield
$-\Delta u \ge f$.

Intuitively, gravity acts downward and so $f(\mathbf x)$ is a negative
function (we choose $f=-10$ in this program). The first condition then means
that the total force acting on the membrane is gravity plus something
positive: namely the upward force that the obstacle exerts on the membrane at
those places where the two of them are in contact. How big is this additional
force? We don't know yet (and neither do we know "where" it actually acts) but
it must be so that the membrane doesn't penetrate the obstacle.

The fourth equality above together with the last inequality forms the obstacle
condition which has to hold at every point of the whole domain. The latter of
these two means that the membrane must be above the obstacle $g(\mathbf x)$
everywhere. The second to last equation, often called the "complementarity
condition" says that where the membrane is not in contact with the obstacle
(i.e., those $\mathbf x$ where $u(\mathbf x) - g(\mathbf x) \neq 0$), then
$-\Delta u=f$ at these locations; in other words, no additional forces act
there, as expected. On the other hand, where $u=g$ we can have $-\Delta u-f
\neq 0$, i.e., there can be additional forces (though there don't have to be:
it is possible for the membrane to just touch, not press against, the
obstacle).


<a name="Derivationofthevariationalinequality"></a><h3>Derivation of the variational inequality</h3>


An obvious way to obtain the variational formulation of the obstacle problem is to consider the total potential energy:
@f{equation*}
 E(u) \dealcoloneq \dfrac{1}{2}\int\limits_{\Omega} \nabla u \cdot \nabla u - \int\limits_{\Omega} fu.
@f}
We have to find a solution $u\in G$ of the following minimization problem:
@f{equation*}
 E(u)\leq E(v)\quad \forall v\in G,
@f}
with the convex set of admissible displacements:
@f{equation*}
 G \dealcoloneq \lbrace v\in V: v\geq g \text{ a.e. in } \Omega\rbrace,\quad V\dealcoloneq H^1_0(\Omega).
@f}
This set takes care of the third and fifth conditions above (the boundary
values and the complementarity condition).

Consider now the minimizer $u\in G$ of $E$ and any other function $v\in
G$. Then the function
@f{equation*}
 F(\varepsilon) \dealcoloneq E(u+\varepsilon(v-u)),\quad\varepsilon\in\left[0,1\right],
@f}
takes its minimum at $\varepsilon = 0$ (because $u$ is a minimizer of the
energy functional $E(\cdot)$), so that $F'(0)\geq 0$ for any choice
of $v$. Note that
$u+\varepsilon(v-u) = (1-\varepsilon)u+\varepsilon v\in G$ because of the
convexity of $G$. If we compute $F'(\varepsilon)\vert_{\varepsilon=0}$ it
yields the variational formulation we are searching for:

<i>Find a function $u\in G$ with</i>
@f{equation*}
 \left(\nabla u, \nabla(v-u)\right) \geq \left(f,v-u\right) \quad \forall v\in G.
@f}

This is the typical form of variational inequalities, where not just $v$
appears in the bilinear form but in fact $v-u$. The reason is this: if $u$ is
not constrained, then we can find test functions $v$ in $G$ so that $v-u$ can have
any sign. By choosing test functions $v_1,v_2$ so that $v_1-u = -(v_2-u)$ it
follows that the inequality can only hold for both $v_1$ and $v_2$ if the two
sides are in fact equal, i.e., we obtain a variational equality.

On the other hand, if $u=g$ then $G$ only allows test functions $v$ so that in fact
$v-u\ge 0$. This means that we can't test the equation with both $v-u$ and
$-(v-u)$ as above, and so we can no longer conclude that the two sides are in
fact equal. Thus, this mimics the way we have discussed the complementarity
condition above.



<a name="Formulationasasaddlepointproblem"></a><h3>Formulation as a saddle point problem</h3>


The variational inequality above is awkward to work with. We would therefore
like to reformulate it as an equivalent saddle point problem. We introduce a
Lagrange multiplier $\lambda$ and the convex cone $K\subset V'$, $V'$
dual space of $V$, $K \dealcoloneq \{\mu\in V': \langle\mu,v\rangle\geq 0,\quad \forall
v\in V, v \le 0 \}$ of
Lagrange multipliers, where $\langle\cdot,\cdot\rangle$ denotes the duality
pairing between $V'$ and $V$. Intuitively, $K$ is the cone of all "non-positive
functions", except that $K\subset (H_0^1)'$ and so contains other objects
besides regular functions as well.
This yields:

<i>Find $u\in V$ and $\lambda\in K$ such that</i>
@f{align*}
 a(u,v) + b(v,\lambda) &= f(v),\quad &&v\in V\\
 b(u,\mu - \lambda) &\leq \langle g,\mu - \lambda\rangle,\quad&&\mu\in K,
@f}
<i>with</i>
@f{align*}
 a(u,v) &\dealcoloneq \left(\nabla u, \nabla v\right),\quad &&u,v\in V\\
 b(u,\mu) &\dealcoloneq \langle u,\mu\rangle,\quad &&u\in V,\quad\mu\in V'.
@f}
In other words, we can consider $\lambda$ as the negative of the additional, positive force that the
obstacle exerts on the membrane. The inequality in the second line of the
statement above only appears to have the wrong sign because we have
$\mu-\lambda<0$ at points where $\lambda=0$, given the definition of $K$.

The existence and uniqueness of $(u,\lambda)\in V\times K$ of this saddle
point problem has been stated in Glowinski, Lions and Tr&eacute;moli&egrave;res: Numerical Analysis of Variational
Inequalities, North-Holland, 1981.



<a name="ActiveSetmethodstosolvethesaddlepointproblem"></a><h3>Active Set methods to solve the saddle point problem</h3>


There are different methods to solve the variational inequality. As one
possibility you can understand the saddle point problem as a convex quadratic program (QP) with
inequality constraints.

To get there, let us assume that we discretize both $u$ and $\lambda$ with the
same finite element space, for example the usual $Q_k$ spaces. We would then
get the equations
@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &[BU-G]_i \geq 0, \quad \Lambda_i \leq 0,\quad \Lambda_i[BU-G]_i = 0
\qquad \forall i.&
@f}
where $B$ is the mass matrix on the chosen finite element space and the
indices $i$ above are for all degrees of freedom in the set $\cal S$ of degrees of
freedom located in the interior of the domain
(we have Dirichlet conditions on the perimeter). However, we
can make our life simpler if we use a particular quadrature rule when
assembling all terms that yield this mass matrix, namely a quadrature formula
where quadrature points are only located at the interpolation points at
which shape functions are defined; since all but one shape function are zero
at these locations, we get a diagonal mass matrix with
@f{align*}
  B_{ii} = \int_\Omega \varphi_i(\mathbf x)^2\ \textrm{d}x,
  \qquad
  B_{ij}=0 \ \text{for } i\neq j.
@f}
To define $G$ we use the same technique as for $B$. In other words, we
define
@f{align*}
  G_{i} = \int_\Omega g_h(x) \varphi_i(\mathbf x)\ \textrm{d}x,
@f}
where $g_h$ is a suitable approximation of $g$. The integral in the definition
of $B_{ii}$ and $G_i$ are then approximated by the trapezoidal rule.
With this, the equations above can be restated as
@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &U_i-B_{ii}^{-1}G_i \ge 0, \quad \Lambda_i \leq 0,\quad \Lambda_i[U_i-B_{ii}^{-1}G_i] = 0
\qquad \forall i\in{\cal S}.&
@f}

Now we define for each degree of freedom $i$ the function
@f{equation*}
 C([BU]_i,\Lambda_i) \dealcoloneq -\Lambda_i + \min\lbrace 0, \Lambda_i + c([BU]_i - G_i) \rbrace,
@f}
with some $c>0$. (In this program we choose $c = 100$. It is a kind of a
penalty parameter which depends on the problem itself and needs to be chosen
large enough; for example there is no convergence for $c = 1$ using the
current program if we use 7 global refinements.)

After some head-scratching one can then convince oneself that the inequalities
above can equivalently be rewritten as
@f{equation*}
 C([BU]_i,\Lambda_i) = 0, \qquad \forall i\in{\cal S}.
@f}
The primal-dual active set strategy we will use here is an iterative scheme which is based on
this condition to predict the next active and inactive sets $\mathcal{A}_k$ and
$\mathcal{F}_k$ (that is, those complementary sets of indices $i$ for which
$U_i$ is either equal to or not equal to the value of the obstacle
$B^{-1}G$). For a more in depth treatment of this approach, see Hintermueller, Ito, Kunisch: The primal-dual active set
strategy as a semismooth newton method, SIAM J. OPTIM., 2003, Vol. 13, No. 3,
pp. 865-888.

<a name="Theprimaldualactivesetalgorithm"></a><h3>The primal-dual active set algorithm</h3>


The algorithm for the primal-dual active set method works as follows (NOTE: $B = B^T$):

1. Initialize $\mathcal{A}_k$ and $\mathcal{F}_k$, such that
 $\mathcal{S}=\mathcal{A}_k\cup\mathcal{F}_k$ and
 $\mathcal{A}_k\cap\mathcal{F}_k=\emptyset$ and set $k=1$.
2. Find the primal-dual pair $(U^k,\Lambda^k)$ that satisfies
 @f{align*}
  AU^k + B\Lambda^k &= F,\\
  [BU^k]_i &= G_i\quad&&\forall i\in\mathcal{A}_k,\\
  \Lambda_i^k &= 0\quad&&\forall i\in\mathcal{F}_k.
 @f}
 Note that the second and third conditions imply that exactly $|S|$ unknowns
 are fixed, with the first condition yielding the remaining $|S|$ equations
 necessary to determine both $U$ and $\Lambda$.
3. Define the new active and inactive sets by
 @f{equation*}
 \begin{split}
  \mathcal{A}_{k+1} \dealcoloneq \lbrace i\in\mathcal{S}:\Lambda^k_i + c([BU^k]_i - G_i)< 0\rbrace,\\
  \mathcal{F}_{k+1} \dealcoloneq \lbrace i\in\mathcal{S}:\Lambda^k_i + c([BU^k]_i - G_i)\geq 0\rbrace.
 \end{split}
 @f}
4. If $\mathcal{A}_{k+1}=\mathcal{A}_k$ (and then, obviously, also
 $\mathcal{F}_{k+1}=\mathcal{F}_k$) then stop, else set $k=k+1$ and go to step
 (2).

The method is called "primal-dual" because it uses both primal (the
displacement $U$) as well as dual variables (the Lagrange multiplier
$\Lambda$) to determine the next active set.

At the end of this section, let us add two observations. First,
for any primal-dual pair $(U^k,\Lambda^k)$ that satisfies these
condition, we can distinguish the following cases:

1. $\Lambda^k_i + c([BU^k]_i - G_i) < 0$ (i active):
  <br>
  Then either $[BU^k]_i<G_i$ and $\Lambda^k_i=0$ (penetration) or $\Lambda^k_i<0$ and $[BU^k]_i=G_i$ (pressing load).
2. $\Lambda^k_i + c([BU^k]_i - G_i)\geq 0$ (i inactive):
  <br>
  Then either $[BU^k]_i\geq G_i$ and $\Lambda^k_i=0$ (no contact) or $\Lambda^k_i\geq0$ and $[BU^k]_i=G_i$ (unpressing load).

Second, the method above appears intuitively correct and useful but a bit ad
hoc. However, it can be derived in a concisely in the following way. To this
end, note that we'd like to solve the nonlinear system
@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &C([BU-G]_i, \Lambda_i) = 0,
\qquad \forall i.&
@f}
We can iteratively solve this by always linearizing around the previous
iterate (i.e., applying a Newton method), but for this we need to linearize
the function $C(\cdot,\cdot)$ that is not differentiable. That said, it is
slantly differentiable, and in fact we have
@f{equation*}
 \dfrac{\partial}{\partial U^k_i}C([BU^k]_i,\Lambda^k_i) = \begin{cases}
                                   cB_{ii},& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)< 0\\
                                   0,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)\geq 0.
                                  \end{cases}
@f}
@f{equation*}
 \dfrac{\partial}{\partial\Lambda^k_i}C([BU^k]_i,\Lambda^k_i) = \begin{cases}
                                   0,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)< 0\\
                                   -1,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)\geq 0.
                                  \end{cases}
@f}
This suggest a semismooth Newton step of the form
@f{equation*}
 \begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & B_{\mathcal{F}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & 0 & B_{\mathcal{A}_k}\\
 0 & 0 & -Id_{\mathcal{F}_k} & 0\\
 0 & cB_{\mathcal{A}_k} & 0 & 0
\end{pmatrix}
\begin{pmatrix}
 \delta U^k_{\mathcal{F}_k}\\ \delta U^k_{\mathcal{A}_k}\\ \delta \Lambda^k_{\mathcal{F}_k}\\ \delta \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=
-\begin{pmatrix}
 (AU^k + \Lambda^k - F)_{\mathcal{F}_k}\\ (AU^k + \Lambda^k - F)_{\mathcal{A}_k}\\ -\Lambda^k_{\mathcal{F}_k}\\ c(B_{\mathcal{A}_k} U^k - G)_{\mathcal{A}_k}
\end{pmatrix},
@f}
where we have split matrices $A,B$ as well as vectors in the natural way into
rows and columns whose indices belong to either the active set
${\mathcal{A}_k}$ or the inactive set ${\mathcal{F}_k}$.

Rather than solving for updates $\delta U, \delta \Lambda$, we can also solve
for the variables we are interested in right away by setting $\delta U^k \dealcoloneq
U^{k+1} - U^k$ and $\delta \Lambda^k \dealcoloneq \Lambda^{k+1} - \Lambda^k$ and
bringing all known terms to the right hand side. This yields
@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & B_{\mathcal{F}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & 0 & B_{\mathcal{A}_k}\\
 0 & 0 & Id_{\mathcal{F}_k} & 0\\
 0 & B_{\mathcal{A}_k} & 0 & 0
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}\\ \Lambda^k_{\mathcal{F}_k}\\ \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}\\ 0\\ G_{\mathcal{A}_k}
\end{pmatrix}.
@f}
These are the equations outlined above in the description of the basic algorithm.

We could even drive this a bit further.
It's easy to see that we can eliminate the third row and the third column
because it implies $\Lambda_{\mathcal{F}_k} = 0$:
@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & B_{\mathcal{A}_k}\\
 0 & B_{\mathcal{A}_k} & 0
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}\\ \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}\\ G_{\mathcal{A}_k}
\end{pmatrix}.
@f}
This shows that one in fact only needs to solve for the Lagrange multipliers
located on the active set. By considering the second row one would then recover
the full Lagrange multiplier vector through
@f{equation*}
 \Lambda^k_S = B^{-1}\left(f_{\mathcal{S}} - A_{\mathcal{S}}U^k_{\mathcal{S}}\right).
@f}
Because of the third row and the fact that $B_{\mathcal{A}_k}$ is a diagonal matrix we are able
to calculate $U^k_{\mathcal{A}_k}=B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}$ directly. We can therefore also write the
linear system as follows:
@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & 0\\
 0 & Id_{\mathcal{A}_k} \\
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k} - A_{\mathcal{F}_k\mathcal{A}_k}B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}
 \\
 B_{\mathcal{A}_k}^{-1}G_{\mathcal{A}_k}
\end{pmatrix}.
@f}
Fortunately, this form is easy to arrive at: we simply build the usual Laplace
linear system
@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} \\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k}
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}
\end{pmatrix},
@f}
and then let the AffineConstraints class eliminate all constrained degrees of
freedom, namely $U^k_{\mathcal{A}_k}=B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}$,
in the same way as if the dofs in $\mathcal{A}_k$ were Dirichlet data. The
result linear system (the second to last one above) is symmetric and positive
definite and we solve it with a CG-method
and the AMG preconditioner from Trilinos.


<a name="Implementation"></a><h3>Implementation</h3>


This tutorial is quite similar to step-4. The general structure of the program
follows step-4 with minor differences:
- We need two new methods, <code>assemble_mass_matrix_diagonal</code> and
  <code>update_solution_and_constraints</code>.
- We need new member variables that denote the constraints we have here.
- We change the preconditioner for the solver.


You may want to read up on step-4 if you want to understand the
current program.
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
 * As usual, at the beginning we include all the header files we need in
 * here. With the exception of the various files that provide interfaces to
 * the Trilinos library, there are no surprises:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/index_set.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_vector.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * namespace Step41
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeObstacleProblemcodeclasstemplate"></a> 
 * <h3>The <code>ObstacleProblem</code> class template</h3>
 * 

 * 
 * This class supplies all function and variables needed to describe the
 * obstacle problem. It is close to what we had to do in step-4, and so
 * relatively simple. The only real new components are the
 * update_solution_and_constraints function that computes the active set and
 * a number of variables that are necessary to describe the original
 * (unconstrained) form of the linear system
 * (<code>complete_system_matrix</code> and
 * <code>complete_system_rhs</code>) as well as the active set itself and
 * the diagonal of the mass matrix $B$ used in scaling Lagrange multipliers
 * in the active set formulation. The rest is as in step-4:
 * 
 * @code
 *   template <int dim>
 *   class ObstacleProblem
 *   {
 *   public:
 *     ObstacleProblem();
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void
 *          assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);
 *     void update_solution_and_constraints();
 *     void solve();
 *     void output_results(const unsigned int iteration) const;
 * 
 *     Triangulation<dim>        triangulation;
 *     FE_Q<dim>                 fe;
 *     DoFHandler<dim>           dof_handler;
 *     AffineConstraints<double> constraints;
 *     IndexSet                  active_set;
 * 
 *     TrilinosWrappers::SparseMatrix system_matrix;
 *     TrilinosWrappers::SparseMatrix complete_system_matrix;
 * 
 *     TrilinosWrappers::MPI::Vector solution;
 *     TrilinosWrappers::MPI::Vector system_rhs;
 *     TrilinosWrappers::MPI::Vector complete_system_rhs;
 *     TrilinosWrappers::MPI::Vector diagonal_of_mass_matrix;
 *     TrilinosWrappers::MPI::Vector contact_force;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandtheobstacle"></a> 
 * <h3>Right hand side, boundary values, and the obstacle</h3>
 * 

 * 
 * In the following, we define classes that describe the right hand side
 * function, the Dirichlet boundary values, and the height of the obstacle
 * as a function of $\mathbf x$. In all three cases, we derive these classes
 * from Function@<dim@>, although in the case of <code>RightHandSide</code>
 * and <code>Obstacle</code> this is more out of convention than necessity
 * since we never pass such objects to the library. In any case, the
 * definition of the right hand side and boundary values classes is obvious
 * given our choice of $f=-10$, $u|_{\partial\Omega}=0$:
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       AssertIndexRange(component, 1);
 * 
 *       return -10;
 *     }
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       AssertIndexRange(component, 1);
 * 
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * We describe the obstacle function by a cascaded barrier (think: stair
 * steps):
 * 
 * @code
 *   template <int dim>
 *   class Obstacle : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override
 *     {
 *       (void)component;
 *       Assert(component == 0, ExcIndexRange(component, 0, 1));
 * 
 *       if (p(0) < -0.5)
 *         return -0.2;
 *       else if (p(0) >= -0.5 && p(0) < 0.0)
 *         return -0.4;
 *       else if (p(0) >= 0.0 && p(0) < 0.5)
 *         return -0.6;
 *       else
 *         return -0.8;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeObstacleProblemcodeclass"></a> 
 * <h3>Implementation of the <code>ObstacleProblem</code> class</h3>
 * 

 * 
 * 

 * 
 * 
 * <a name="ObstacleProblemObstacleProblem"></a> 
 * <h4>ObstacleProblem::ObstacleProblem</h4>
 * 

 * 
 * To everyone who has taken a look at the first few tutorial programs, the
 * constructor is completely obvious:
 * 
 * @code
 *   template <int dim>
 *   ObstacleProblem<dim>::ObstacleProblem()
 *     : fe(1)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemmake_grid"></a> 
 * <h4>ObstacleProblem::make_grid</h4>
 * 

 * 
 * We solve our obstacle problem on the square $[-1,1]\times [-1,1]$ in
 * 2D. This function therefore just sets up one of the simplest possible
 * meshes.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -1, 1);
 *     triangulation.refine_global(7);
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Total number of cells: " << triangulation.n_cells()
 *               << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemsetup_system"></a> 
 * <h4>ObstacleProblem::setup_system</h4>
 * 

 * 
 * In this first function of note, we set up the degrees of freedom handler,
 * resize vectors and matrices, and deal with the constraints. Initially,
 * the constraints are, of course, only given by boundary values, so we
 * interpolate them towards the top of the function.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     active_set.set_size(dof_handler.n_dofs());
 * 
 *     std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl
 *               << std::endl;
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              BoundaryValues<dim>(),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 * 
 *     system_matrix.reinit(dsp);
 *     complete_system_matrix.reinit(dsp);
 * 
 *     IndexSet solution_index_set = dof_handler.locally_owned_dofs();
 *     solution.reinit(solution_index_set, MPI_COMM_WORLD);
 *     system_rhs.reinit(solution_index_set, MPI_COMM_WORLD);
 *     complete_system_rhs.reinit(solution_index_set, MPI_COMM_WORLD);
 *     contact_force.reinit(solution_index_set, MPI_COMM_WORLD);
 * 
 * @endcode
 * 
 * The only other thing to do here is to compute the factors in the $B$
 * matrix which is used to scale the residual. As discussed in the
 * introduction, we'll use a little trick to make this mass matrix
 * diagonal, and in the following then first compute all of this as a
 * matrix and then extract the diagonal elements for later use:
 * 
 * @code
 *     TrilinosWrappers::SparseMatrix mass_matrix;
 *     mass_matrix.reinit(dsp);
 *     assemble_mass_matrix_diagonal(mass_matrix);
 *     diagonal_of_mass_matrix.reinit(solution_index_set);
 *     for (unsigned int j = 0; j < solution.size(); j++)
 *       diagonal_of_mass_matrix(j) = mass_matrix.diag_element(j);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemassemble_system"></a> 
 * <h4>ObstacleProblem::assemble_system</h4>
 * 

 * 
 * This function at once assembles the system matrix and right-hand-side and
 * applied the constraints (both due to the active set as well as from
 * boundary values) to our system. Otherwise, it is functionally equivalent
 * to the corresponding function in, for example, step-4.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::assemble_system()
 *   {
 *     std::cout << "   Assembling system..." << std::endl;
 * 
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const QGauss<dim>  quadrature_formula(fe.degree + 1);
 *     RightHandSide<dim> right_hand_side;
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
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 cell_matrix(i, j) +=
 *                   (fe_values.shape_grad(i, q_point) *
 *                    fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));
 * 
 *               cell_rhs(i) +=
 *                 (fe_values.shape_value(i, q_point) *
 *                  right_hand_side.value(fe_values.quadrature_point(q_point)) *
 *                  fe_values.JxW(q_point));
 *             }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraints.distribute_local_to_global(cell_matrix,
 *                                                cell_rhs,
 *                                                local_dof_indices,
 *                                                system_matrix,
 *                                                system_rhs,
 *                                                true);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemassemble_mass_matrix_diagonal"></a> 
 * <h4>ObstacleProblem::assemble_mass_matrix_diagonal</h4>
 * 

 * 
 * The next function is used in the computation of the diagonal mass matrix
 * $B$ used to scale variables in the active set method. As discussed in the
 * introduction, we get the mass matrix to be diagonal by choosing the
 * trapezoidal rule for quadrature. Doing so we don't really need the triple
 * loop over quadrature points, indices $i$ and indices $j$ any more and
 * can, instead, just use a double loop. The rest of the function is obvious
 * given what we have discussed in many of the previous tutorial programs.
 *   

 * 
 * Note that at the time this function is called, the constraints object
 * only contains boundary value constraints; we therefore do not have to pay
 * attention in the last copy-local-to-global step to preserve the values of
 * matrix entries that may later on be constrained by the active set.
 *   

 * 
 * Note also that the trick with the trapezoidal rule only works if we have
 * in fact $Q_1$ elements. For higher order elements, one would need to use
 * a quadrature formula that has quadrature points at all the support points
 * of the finite element. Constructing such a quadrature formula isn't
 * really difficult, but not the point here, and so we simply assert at the
 * top of the function that our implicit assumption about the finite element
 * is in fact satisfied.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::assemble_mass_matrix_diagonal(
 *     TrilinosWrappers::SparseMatrix &mass_matrix)
 *   {
 *     Assert(fe.degree == 1, ExcNotImplemented());
 * 
 *     const QTrapezoid<dim> quadrature_formula;
 *     FEValues<dim>         fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         cell_matrix = 0;
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             cell_matrix(i, i) +=
 *               (fe_values.shape_value(i, q_point) *
 *                fe_values.shape_value(i, q_point) * fe_values.JxW(q_point));
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraints.distribute_local_to_global(cell_matrix,
 *                                                local_dof_indices,
 *                                                mass_matrix);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemupdate_solution_and_constraints"></a> 
 * <h4>ObstacleProblem::update_solution_and_constraints</h4>
 * 

 * 
 * In a sense, this is the central function of this program.  It updates the
 * active set of constrained degrees of freedom as discussed in the
 * introduction and computes an AffineConstraints object from it that can then
 * be used to eliminate constrained degrees of freedom from the solution of
 * the next iteration. At the same time we set the constrained degrees of
 * freedom of the solution to the correct value, namely the height of the
 * obstacle.
 *   

 * 
 * Fundamentally, the function is rather simple: We have to loop over all
 * degrees of freedom and check the sign of the function $\Lambda^k_i +
 * c([BU^k]_i - G_i) = \Lambda^k_i + cB_i(U^k_i - [g_h]_i)$ because in our
 * case $G_i = B_i[g_h]_i$. To this end, we use the formula given in the
 * introduction by which we can compute the Lagrange multiplier as the
 * residual of the original linear system (given via the variables
 * <code>complete_system_matrix</code> and <code>complete_system_rhs</code>.
 * At the top of this function, we compute this residual using a function
 * that is part of the matrix classes.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::update_solution_and_constraints()
 *   {
 *     std::cout << "   Updating active set..." << std::endl;
 * 
 *     const double penalty_parameter = 100.0;
 * 
 *     TrilinosWrappers::MPI::Vector lambda(
 *       complete_index_set(dof_handler.n_dofs()));
 *     complete_system_matrix.residual(lambda, solution, complete_system_rhs);
 * 
 * @endcode
 * 
 * compute contact_force[i] = - lambda[i] * diagonal_of_mass_matrix[i]
 * 
 * @code
 *     contact_force = lambda;
 *     contact_force.scale(diagonal_of_mass_matrix);
 *     contact_force *= -1;
 * 
 * @endcode
 * 
 * The next step is to reset the active set and constraints objects and to
 * start the loop over all degrees of freedom. This is made slightly more
 * complicated by the fact that we can't just loop over all elements of
 * the solution vector since there is no way for us then to find out what
 * location a DoF is associated with; however, we need this location to
 * test whether the displacement of a DoF is larger or smaller than the
 * height of the obstacle at this location.
 *     

 * 
 * We work around this by looping over all cells and DoFs defined on each
 * of these cells. We use here that the displacement is described using a
 * $Q_1$ function for which degrees of freedom are always located on the
 * vertices of the cell; thus, we can get the index of each degree of
 * freedom and its location by asking the vertex for this information. On
 * the other hand, this clearly wouldn't work for higher order elements,
 * and so we add an assertion that makes sure that we only deal with
 * elements for which all degrees of freedom are located in vertices to
 * avoid tripping ourselves with non-functional code in case someone wants
 * to play with increasing the polynomial degree of the solution.
 *     

 * 
 * The price to pay for having to loop over cells rather than DoFs is that
 * we may encounter some degrees of freedom more than once, namely each
 * time we visit one of the cells adjacent to a given vertex. We will
 * therefore have to keep track which vertices we have already touched and
 * which we haven't so far. We do so by using an array of flags
 * <code>dof_touched</code>:
 * 
 * @code
 *     constraints.clear();
 *     active_set.clear();
 * 
 *     const Obstacle<dim> obstacle;
 *     std::vector<bool>   dof_touched(dof_handler.n_dofs(), false);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       for (const auto v : cell->vertex_indices())
 *         {
 *           Assert(dof_handler.get_fe().n_dofs_per_cell() == cell->n_vertices(),
 *                  ExcNotImplemented());
 * 
 *           const unsigned int dof_index = cell->vertex_dof_index(v, 0);
 * 
 *           if (dof_touched[dof_index] == false)
 *             dof_touched[dof_index] = true;
 *           else
 *             continue;
 * 
 * @endcode
 * 
 * Now that we know that we haven't touched this DoF yet, let's get
 * the value of the displacement function there as well as the value
 * of the obstacle function and use this to decide whether the
 * current DoF belongs to the active set. For that we use the
 * function given above and in the introduction.
 *           

 * 
 * If we decide that the DoF should be part of the active set, we
 * add its index to the active set, introduce an inhomogeneous
 * equality constraint in the AffineConstraints object, and reset the
 * solution value to the height of the obstacle. Finally, the
 * residual of the non-contact part of the system serves as an
 * additional control (the residual equals the remaining,
 * unaccounted forces, and should be zero outside the contact zone),
 * so we zero out the components of the residual vector (i.e., the
 * Lagrange multiplier lambda) that correspond to the area where the
 * body is in contact; at the end of the loop over all cells, the
 * residual will therefore only consist of the residual in the
 * non-contact zone. We output the norm of this residual along with
 * the size of the active set after the loop.
 * 
 * @code
 *           const double obstacle_value = obstacle.value(cell->vertex(v));
 *           const double solution_value = solution(dof_index);
 * 
 *           if (lambda(dof_index) + penalty_parameter *
 *                                     diagonal_of_mass_matrix(dof_index) *
 *                                     (solution_value - obstacle_value) <
 *               0)
 *             {
 *               active_set.add_index(dof_index);
 *               constraints.add_line(dof_index);
 *               constraints.set_inhomogeneity(dof_index, obstacle_value);
 * 
 *               solution(dof_index) = obstacle_value;
 * 
 *               lambda(dof_index) = 0;
 *             }
 *         }
 *     std::cout << "      Size of active set: " << active_set.n_elements()
 *               << std::endl;
 * 
 *     std::cout << "   Residual of the non-contact part of the system: "
 *               << lambda.l2_norm() << std::endl;
 * 
 * @endcode
 * 
 * In a final step, we add to the set of constraints on DoFs we have so
 * far from the active set those that result from Dirichlet boundary
 * values, and close the constraints object:
 * 
 * @code
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              BoundaryValues<dim>(),
 *                                              constraints);
 *     constraints.close();
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemsolve"></a> 
 * <h4>ObstacleProblem::solve</h4>
 * 

 * 
 * There is nothing to say really about the solve function. In the context
 * of a Newton method, we are not typically interested in very high accuracy
 * (why ask for a highly accurate solution of a linear problem that we know
 * only gives us an approximation of the solution of the nonlinear problem),
 * and so we use the ReductionControl class that stops iterations when
 * either an absolute tolerance is reached (for which we choose $10^{-12}$)
 * or when the residual is reduced by a certain factor (here, $10^{-3}$).
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::solve()
 *   {
 *     std::cout << "   Solving system..." << std::endl;
 * 
 *     ReductionControl                        reduction_control(100, 1e-12, 1e-3);
 *     SolverCG<TrilinosWrappers::MPI::Vector> solver(reduction_control);
 *     TrilinosWrappers::PreconditionAMG       precondition;
 *     precondition.initialize(system_matrix);
 * 
 *     solver.solve(system_matrix, solution, system_rhs, precondition);
 *     constraints.distribute(solution);
 * 
 *     std::cout << "      Error: " << reduction_control.initial_value() << " -> "
 *               << reduction_control.last_value() << " in "
 *               << reduction_control.last_step() << " CG iterations."
 *               << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemoutput_results"></a> 
 * <h4>ObstacleProblem::output_results</h4>
 * 

 * 
 * We use the vtk-format for the output.  The file contains the displacement
 * and a numerical representation of the active set.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::output_results(const unsigned int iteration) const
 *   {
 *     std::cout << "   Writing graphical output..." << std::endl;
 * 
 *     TrilinosWrappers::MPI::Vector active_set_vector(
 *       dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
 *     for (const auto index : active_set)
 *       active_set_vector[index] = 1.;
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "displacement");
 *     data_out.add_data_vector(active_set_vector, "active_set");
 *     data_out.add_data_vector(contact_force, "lambda");
 * 
 *     data_out.build_patches();
 * 
 *     std::ofstream output_vtk("output_" +
 *                              Utilities::int_to_string(iteration, 3) + ".vtk");
 *     data_out.write_vtk(output_vtk);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemrun"></a> 
 * <h4>ObstacleProblem::run</h4>
 * 

 * 
 * This is the function which has the top-level control over everything.  It
 * is not very long, and in fact rather straightforward: in every iteration
 * of the active set method, we assemble the linear system, solve it, update
 * the active set and project the solution back to the feasible set, and
 * then output the results. The iteration is terminated whenever the active
 * set has not changed in the previous iteration.
 *   

 * 
 * The only trickier part is that we have to save the linear system (i.e.,
 * the matrix and right hand side) after assembling it in the first
 * iteration. The reason is that this is the only step where we can access
 * the linear system as built without any of the contact constraints
 * active. We need this to compute the residual of the solution at other
 * iterations, but in other iterations that linear system we form has the
 * rows and columns that correspond to constrained degrees of freedom
 * eliminated, and so we can no longer access the full residual of the
 * original equation.
 * 
 * @code
 *   template <int dim>
 *   void ObstacleProblem<dim>::run()
 *   {
 *     make_grid();
 *     setup_system();
 * 
 *     IndexSet active_set_old(active_set);
 *     for (unsigned int iteration = 0; iteration <= solution.size(); ++iteration)
 *       {
 *         std::cout << "Newton iteration " << iteration << std::endl;
 * 
 *         assemble_system();
 * 
 *         if (iteration == 0)
 *           {
 *             complete_system_matrix.copy_from(system_matrix);
 *             complete_system_rhs = system_rhs;
 *           }
 * 
 *         solve();
 *         update_solution_and_constraints();
 *         output_results(iteration);
 * 
 *         if (active_set == active_set_old)
 *           break;
 * 
 *         active_set_old = active_set;
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * } // namespace Step41
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * And this is the main function. It follows the pattern of all other main
 * functions. The call to initialize MPI exists because the Trilinos library
 * upon which we build our linear solvers in this program requires it.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step41;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *         argc, argv, numbers::invalid_unsigned_int);
 * 
 * @endcode
 * 
 * This program can only be run in serial. Otherwise, throw an exception.
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
 *                   ExcMessage(
 *                     "This program can only be run in serial, use ./step-41"));
 * 
 *       ObstacleProblem<2> obstacle_problem;
 *       obstacle_problem.run();
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


Running the program produces output like this:
@code
Number of active cells: 16384
Total number of cells: 21845
Number of degrees of freedom: 16641

Newton iteration 0
   Assembling system...
   Solving system...
      Error: 0.310059 -> 5.16619e-05 in 5 CG iterations.
   Updating active set...
      Size of active set: 13164
   Residual of the non-contact part of the system: 1.61863e-05
   Writing graphical output...

Newton iteration 1
   Assembling system...
   Solving system...
      Error: 1.11987 -> 0.00109377 in 6 CG iterations.
   Updating active set...
      Size of active set: 12363
   Residual of the non-contact part of the system: 3.9373
   Writing graphical output...

...

Newton iteration 17
   Assembling system...
   Solving system...
      Error: 0.00713308 -> 2.29249e-06 in 4 CG iterations.
   Updating active set...
      Size of active set: 5399
   Residual of the non-contact part of the system: 0.000957525
   Writing graphical output...

Newton iteration 18
   Assembling system...
   Solving system...
      Error: 0.000957525 -> 2.8033e-07 in 4 CG iterations.
   Updating active set...
      Size of active set: 5399
   Residual of the non-contact part of the system: 2.8033e-07
   Writing graphical output...
@endcode

The iterations end once the active set doesn't change any more (it has
5,399 constrained degrees of freedom at that point). The algebraic
precondition is apparently working nicely since we only need 4-6 CG
iterations to solve the linear system (although this also has a lot to
do with the fact that we are not asking for very high accuracy of the
linear solver).

More revealing is to look at a sequence of graphical output files
(every third step is shown, with the number of the iteration in the
leftmost column):

<table align="center">
  <tr>
    <td valign="top">
      0 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.00.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.00.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.00.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      3 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.03.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.03.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      6 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.06.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.06.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.06.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      9 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.09.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      12 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.12.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.12.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.12.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      15 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.15.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.15.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.15.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      18 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.18.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.18.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.18.png" alt="">
    </td>
  </tr>
</table>

The pictures show that in the first step, the solution (which has been
computed without any of the constraints active) bends through so much
that pretty much every interior point has to be bounced back to the
stairstep function, producing a discontinuous solution. Over the
course of the active set iterations, this unphysical membrane shape is
smoothed out, the contact with the lower-most stair step disappears,
and the solution stabilizes.

In addition to this, the program also outputs the values of the
Lagrange multipliers. Remember that these are the contact forces and
so should only be positive on the contact set, and zero outside. If,
on the other hand, a Lagrange multiplier is negative in the active
set, then this degree of freedom must be removed from the active
set. The following pictures show the multipliers in iterations 1, 9
and 18, where we use red and browns to indicate positive values, and
blue for negative values.

<table align="center">
  <tr>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.01.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.18.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      Iteration 1
    </td>
    <td align="center">
      Iteration 9
    </td>
    <td align="center">
      Iteration 18
    </td>
  </tr>
</table>

It is easy to see that the positive values converge nicely to moderate
values in the interior of the contact set and large upward forces at
the edges of the steps, as one would expect (to support the large
curvature of the membrane there); at the fringes of the active set,
multipliers are initially negative, causing the set to shrink until,
in iteration 18, there are no more negative multipliers and the
algorithm has converged.



<a name="extensions"></a>
<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


As with any of the programs of this tutorial, there are a number of
obvious possibilities for extensions and experiments. The first one is
clear: introduce adaptivity. Contact problems are prime candidates for
adaptive meshes because the solution has lines along which it is less
regular (the places where contact is established between membrane and
obstacle) and other areas where the solution is very smooth (or, in
the present context, constant wherever it is in contact with the
obstacle). Adding this to the current program should not pose too many
difficulties, but it is not trivial to find a good error estimator for
that purpose.

A more challenging task would be an extension to 3d. The problem here
is not so much to simply make everything run in 3d. Rather, it is that
when a 3d body is deformed and gets into contact with an obstacle,
then the obstacle does not act as a constraining body force within the
domain as is the case here. Rather, the contact force only acts on the
boundary of the object. The inequality then is not in the differential
equation but in fact in the (Neumann-type) boundary conditions, though
this leads to a similar kind of variational
inequality. Mathematically, this means that the Lagrange multiplier
only lives on the surface, though it can of course be extended by zero
into the domain if that is convenient. As in the current program, one
does not need to form and store this Lagrange multiplier explicitly.

A further interesting problem for the 3d case is to consider contact problems
with friction. In almost every mechanical process friction has a big influence.
For the modelling we have to take into account tangential stresses at the contact
surface. Also we have to observe that friction adds another nonlinearity to
our problem.

Another nontrivial modification is to implement a more complex constitutive
law like nonlinear elasticity or elasto-plastic  material behavior.
The difficulty here is to handle the additional nonlinearity arising
through the nonlinear constitutive law.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-41.cc"
*/
