/**
@page step_47 The step-47 tutorial program
This tutorial depends on step-12.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Whatstheissue"> What's the issue? </a>
        <li><a href="#Whattodoinstead"> What to do instead? </a>
        <li><a href="#DerivationoftheC0IPmethod"> Derivation of the C0IP method </a>
      <ul>
        <li><a href="#ConvergenceRates">Convergence Rates </a>
      </ul>
        <li><a href="#OtherBoundaryConditions">Other Boundary Conditions</a>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Themainclass">The main class</a>
      <ul>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
        <li><a href="#Solvingthelinearsystemandpostprocessing">Solving the linear system and postprocessing</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestresultsoniQsub2subiwithigammapp1i">Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub3subiwithigammapp1i">Test results on <i>Q<sub>3</sub></i> with <i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub4subiwithigammapp1i">Test results on <i>Q<sub>4</sub></i> with <i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub2subiwithigamma1i">Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = 1</i> </a>
        <li><a href="#TestresultsoniQsub2subiwithigamma2i">Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = 2</i> </a>
        <li><a href="#Conclusionsforthechoiceofthepenaltyparameter"> Conclusions for the choice of the penalty parameter </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Derivationforthesimplysupportedplates"> Derivation for the simply supported plates </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<br>

<i>
This program was contributed by Natasha Sharma, Guido Kanschat, Timo
Heister, Wolfgang Bangerth, and Zhuoran Wang.

The first author would like to acknowledge the support of NSF Grant
No. DMS-1520862.
Timo Heister and Wolfgang Bangerth acknowledge support through NSF
awards DMS-1821210, EAR-1550901, and OAC-1835673.
</i>

<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


This program deals with the <a
href="https://en.wikipedia.org/wiki/Biharmonic_equation">biharmonic
equation</a>,
@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega.
@f}
This equation appears in the modeling of thin structures such as roofs
of stadiums. These objects are of course in reality
three-dimensional with a large aspect ratio of lateral extent to
perpendicular thickness, but one can often very accurately model these
structures as two dimensional by making assumptions about how internal
forces vary in the perpendicular direction. These assumptions lead to the
equation above.

The model typically comes in two different kinds, depending on what
kinds of boundary conditions are imposed. The first case,
@f{align*}{
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega,
@f}
corresponds to the edges of the thin structure attached to the top of
a wall of height $g(\mathbf x)$ in such a way that the bending forces
that act on the structure are $h(\mathbf x)$; in most physical
situations, one will have $h=0$, corresponding to the structure simply
sitting atop the wall.

In the second possible case of boundary values, one would have
@f{align*}{
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \frac{\partial u(\mathbf x)}{\partial \mathbf n} &= j(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.
@f}
This corresponds to a "clamped" structure for which a nonzero
$j(\mathbf x)$ implies a certain angle against the horizontal.

As with Dirichlet and Neumann boundary conditions for the Laplace
equation, it is of course possible to have one kind of boundary
conditions on one part of the boundary, and the other on the
remainder.


<a name="Whatstheissue"></a><h3> What's the issue? </h3>


The fundamental issue with the equation is that it takes four
derivatives of the solution. In the case of the Laplace equation
we treated in step-3, step-4, and several other tutorial programs,
one multiplies by a test function, integrates, integrates by parts,
and ends up with only one derivative on both the test function and
trial function -- something one can do with functions that are
continuous globally, but may have kinks at the interfaces between
cells: The derivative may not be defined at the interfaces, but
that is on a lower-dimensional manifold (and so doesn't show up
in the integrated value).

But for the biharmonic equation, if one followed the same procedure
using integrals over the entire domain (i.e., the union of all cells),
one would end up with two derivatives on the test functions and trial
functions each. If one were to use the usual piecewise polynomial
functions with their kinks on cell interfaces, the first derivative
would yield a discontinuous gradient, and the second derivative with
delta functions on the interfaces -- but because both the second
derivatives of the test functions and of the trial functions yield a
delta function, we would try to integrate the product of two delta
functions. For example, in 1d, where $\varphi_i$ are the usual
piecewise linear "hat functions", we would get integrals of the sort
@f{align*}{
  \int_0^L (\Delta \varphi_i) (\Delta \varphi_j)
  =
  \int_0^L
  \frac 1h \left[\delta(x-x_{i-1}) - 2\delta(x-x_i) + \delta(x-x_{i+1})\right]
  \frac 1h \left[\delta(x-x_{j-1}) - 2\delta(x-x_j) + \delta(x-x_{j+1})\right]
@f}
where $x_i$ is the node location at which the shape function
$\varphi_i$ is defined, and $h$ is the mesh size (assumed
uniform). The problem is that delta functions in integrals are defined
using the relationship
@f{align*}{
  \int_0^L \delta(x-\hat x) f(x) \; dx
  =
  f(\hat x).
@f}
But that only works if (i) $f(\cdot)$ is actually well defined at
$\hat x$, and (ii) if it is finite. On the other hand, an integral of
the form
@f{align*}{
\int_0^L \delta(x-x_i) \delta (x-x_i)
@f}
does not make sense. Similar reasoning can be applied for 2d and 3d
situations.

In other words: This approach of trying to integrate over the entire
domain and then integrating by parts can't work.

Historically, numerical analysts have tried to address this by
inventing finite elements that are "C<sup>1</sup> continuous", i.e., that use
shape functions that are not just continuous but also have continuous
first derivatives. This is the realm of elements such as the Argyris
element, the Clough-Tocher element and others, all developed in the
late 1960s. From a twenty-first century perspective, they can only be
described as bizarre in their construction. They are also exceedingly
cumbersome to implement if one wants to use general meshes. As a
consequence, they have largely fallen out of favor and deal.II currently
does not contain implementations of these shape functions.


<a name="Whattodoinstead"></a><h3> What to do instead? </h3>


So how does one approach solving such problems then? That depends a
bit on the boundary conditions. If one has the first set of boundary
conditions, i.e., if the equation is
@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega,
@f}
then the following trick works (at least if the domain is convex, see
below): In the same way as we obtained the
mixed Laplace equation of step-20 from the regular Laplace equation by
introducing a second variable, we can here introduce a variable
$v=\Delta u$ and can then replace the equations above by the
following, "mixed" system:
@f{align*}{
  -\Delta u(\mathbf x) +v(\mathbf x) &= 0
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  -\Delta v(\mathbf x) &= -f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  v(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.
@f}
In other words, we end up with what is in essence a system of two
coupled Laplace equations for $u,v$, each with Dirichlet-type boundary
conditions. We know how to solve such problems, and it should not be
very difficult to construct good solvers and preconditioners for this
system either using the techniques of step-20 or step-22. So this
case is pretty simple to deal with.

@note It is worth pointing out that this only works for domains whose
  boundary has corners if the domain is also convex -- in other words,
  if there are no re-entrant corners.
  This sounds like a rather random condition, but it makes
  sense in view of the following two facts: The solution of the
  original biharmonic equation must satisfy $u\in H^2(\Omega)$. On the
  other hand, the mixed system reformulation above suggests that both
  $u$ and $v$ satisfy $u,v\in H^1(\Omega)$ because both variables only
  solve a Poisson equation. In other words, if we want to ensure that
  the solution $u$ of the mixed problem is also a solution of the
  original biharmonic equation, then we need to be able to somehow
  guarantee that the solution of $-\Delta u=v$ is in fact more smooth
  than just $H^1(\Omega)$. This can be argued as follows: For convex
  domains,
  <a href="https://en.wikipedia.org/wiki/Elliptic_operator#Elliptic_regularity_theorem">"elliptic
  regularity"</a> implies that if the right hand side $v\in H^s$, then
  $u\in H^{s+2}$ if the domain is convex and the boundary is smooth
  enough. (This could also be guaranteed if the domain boundary is
  sufficiently smooth -- but domains whose boundaries have no corners
  are not very practical in real life.)
  We know that $v\in H^1$ because it solves the equation
  $-\Delta v=f$, but we are still left with the condition on convexity
  of the boundary; one can show that polygonal, convex domains are
  good enough to guarantee that $u\in H^2$ in this case (smoothly
  bounded, convex domains would result in $u\in H^3$, but we don't
  need this much regularity). On the other hand, if the domain is not
  convex, we can not guarantee that the solution of the mixed system
  is in $H^2$, and consequently may obtain a solution that can't be
  equal to the solution of the original biharmonic equation.

The more complicated situation is if we have the "clamped" boundary
conditions, i.e., if the equation looks like this:
@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \frac{\partial u(\mathbf x)}{\partial \mathbf n} &= j(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.
@f}
The same trick with the mixed system does not work here, because we
would end up with <i>both</i> Dirichlet and Neumann boundary conditions for
$u$, but none for $v$.


The solution to this conundrum arrived with the Discontinuous Galerkin
method wave in the 1990s and early 2000s: In much the same way as one
can use <i>discontinuous</i> shape functions for the Laplace equation
by penalizing the size of the discontinuity to obtain a scheme for an
equation that has one derivative on each shape function, we can use a
scheme that uses <i>continuous</i> (but not $C^1$ continuous) shape
functions and penalize the jump in the derivative to obtain a scheme
for an equation that has two derivatives on each shape function. In
analogy to the Interior Penalty (IP) method for the Laplace equation,
this scheme for the biharmonic equation is typically called the $C^0$ IP
(or C0IP) method, since it uses $C^0$ (continuous but not continuously
differentiable) shape functions with an interior penalty formulation.


<a name="DerivationoftheC0IPmethod"></a><h3> Derivation of the C0IP method </h3>


We base this program on the $C^0$ IP method presented by Susanne
Brenner and Li-Yeng Sung in the paper "C$^0$ Interior Penalty Method
for Linear Fourth Order Boundary Value Problems on polygonal
domains'' @cite Brenner2005 , where the method is
derived for the biharmonic equation with "clamped" boundary
conditions.

As mentioned, this method relies on the use of $C^0$ Lagrange finite
elements where the $C^1$ continuity requirement is relaxed and has
been replaced with interior penalty techniques. To derive this method,
we consider a $C^0$ shape function $v_h$ which vanishes on
$\partial\Omega$. We introduce notation $ \mathbb{F} $ as the set of
all faces of $\mathbb{T}$, $ \mathbb{F}^b $ as the set of boundary faces,
and $ \mathbb{F}^i $ as the set of interior faces for use further down below.
Since the higher order derivatives of $v_h$ have two
values on each interface $e\in \mathbb{F}$ (shared by the two cells
$K_{+},K_{-} \in \mathbb{T}$), we cope with this discontinuity by
defining the following single-valued functions on $e$:
@f{align*}{
  \jump{\frac{\partial^k v_h}{\partial \mathbf n^k}}
  &=
  \frac{\partial^k v_h|_{K_+}}{\partial \mathbf n^k} \bigg |_e
  - \frac{\partial^k v_h|_{K_-}}{\partial \mathbf n^k} \bigg |_e,
  \\
  \average{\frac{\partial^k v_h}{\partial \mathbf n^k}}
  &=
  \frac{1}{2}
  \bigg( \frac{\partial^k v_h|_{K_+}}{\partial \mathbf n^k} \bigg |_e
  + \frac{\partial^k v_h|_{K_-}}{\partial \mathbf n^k} \bigg |_e \bigg )
@f}
for $k =1,2$ (i.e., for the gradient and the matrix of second
derivatives), and where $\mathbf n$ denotes a unit vector normal to
$e$ pointing from $K_+$ to $K_-$. In the
literature, these functions are referred to as the "jump" and
"average" operations, respectively.

To obtain the $C^0$ IP approximation $u_h$, we left multiply the
biharmonic equation by $v_h$, and then integrate over $\Omega$. As
explained above, we can't do the integration by parts on all of
$\Omega$ with these shape functions, but we can do it on each cell
individually since the shape functions are just polynomials on each
cell. Consequently, we start by using the following
integration-by-parts formula on each mesh cell $K \in {\mathbb{T}}$:
@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K v_h (\nabla\cdot\nabla) (\Delta w_h)
  \\
  &= -\int_K \nabla v_h \cdot (\nabla \Delta w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n).
@f}
At this point, we have two options: We can integrate the domain term's
$\nabla\Delta w_h$ one more time to obtain
@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K (\Delta v_h) (\Delta w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)
     -\int_{\partial K} (\nabla v_h \cdot \mathbf n) \Delta w_h.
@f}
For a variety of reasons, this turns out to be a variation that is not
useful for our purposes.

Instead, what we do is recognize that
$\nabla\Delta w_h = \text{grad}\,(\text{div}\,\text{grad}\, w_h)$, and we
can re-sort these operations as
$\nabla\Delta w_h = \text{div}\,(\text{grad}\,\text{grad}\, w_h)$ where we
typically write $\text{grad}\,\text{grad}\, w_h = D^2 w_h$ to indicate
that this is the "Hessian" matrix of second derivatives. With this
re-ordering, we can now integrate the divergence, rather than the
gradient operator, and we get the following instead:
@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K (\nabla \nabla v_h) : (\nabla \nabla w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)
     -\int_{\partial K} (\nabla v_h \otimes \mathbf n) : (\nabla\nabla w_h)
  \\
  &= \int_K (D^2 v_h) : (D^2 w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)
     -\int_{\partial K} (\nabla v_h) \cdot (D^2 w_h \mathbf n).
@f}
Here, the colon indicates a double-contraction over the indices of the
matrices to its left and right, i.e., the scalar product between two
tensors. The outer product of two vectors $a \otimes b$ yields the
matrix $(a \otimes b)_{ij} = a_i b_j$.

Then, we sum over all cells $K \in  \mathbb{T}$, and take into account
that this means that every interior face appears twice in the
sum. If we therefore split everything into a sum of integrals over
cell interiors and a separate sum over cell interfaces, we can use
the jump and average operators defined above. There are two steps
left: First, because our shape functions are continuous, the gradients
of the shape functions may be discontinuous, but the continuity
guarantees that really only the normal component of the gradient is
discontinuous across faces whereas the tangential component(s) are
continuous. Second, the discrete formulation that results is not
stable as the mesh size goes to zero, and to obtain a stable
formulation that converges to the correct solution, we need to add
the following terms:
@f{align*}{
-\sum_{e \in \mathbb{F}} \int_{e}
  \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
  \jump{\frac{\partial u_h}{\partial \mathbf n}}
+ \sum_{e \in \mathbb{F}}
  \frac{\gamma}{h_e}\int_e
  \jump{\frac{\partial v_h}{\partial \mathbf n}}
  \jump{\frac{\partial u_h}{\partial \mathbf n}}.
@f}
Then, after making cancellations that arise, we arrive at the following
C0IP formulation of the biharmonic equation: find $u_h$ such that $u_h =
g$ on $\partial \Omega$ and
@f{align*}{
\mathcal{A}(v_h,u_h)&=\mathcal{F}(v_h) \quad \text{holds for all test functions } v_h,
@f}
where
@f{align*}{
\mathcal{A}(v_h,u_h):=&\sum_{K \in \mathbb{T}}\int_K D^2v_h:D^2u_h \ dx
\\
&
 -\sum_{e \in \mathbb{F}} \int_{e}
  \jump{\frac{\partial v_h}{\partial \mathbf n}}
  \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
 -\sum_{e \in \mathbb{F}} \int_{e}
 \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
\\
&+ \sum_{e \in \mathbb{F}}
 \frac{\gamma}{h_e}
 \int_e
 \jump{\frac{\partial v_h}{\partial \mathbf n}}
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds,
@f}
and
@f{align*}{
\mathcal{F}(v_h)&:=\sum_{K \in \mathbb{T}}\int_{K} v_h f \ dx
-
\sum_{e \in \mathbb{F}, e\subset\partial\Omega}
\int_e \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}} j \ ds
+
\sum_{e \in \mathbb{F}, e\subset\partial\Omega}
\frac{\gamma}{h_e}
\int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds.
@f}
Here, $\gamma$ is the penalty parameter which both weakly enforces the
boundary condition
@f{align*}{
\frac{\partial u(\mathbf x)}{\partial \mathbf n} = j(\mathbf x)
@f}
on the boundary interfaces $e \in \mathbb{F}^b$, and also ensures that
in the limit $h\rightarrow 0$, $u_h$ converges to a $C^1$ continuous
function. $\gamma$ is chosen to be large enough to guarantee the
stability of the method. We will discuss our choice in the program below.


<a name="ConvergenceRates"></a><h4>Convergence Rates </h4>

On polygonal domains, the weak solution $u$ to the biharmonic equation
lives in $H^{2 +\alpha}(\Omega)$ where $\alpha \in(1/2, 2]$ is
determined by the interior angles at the corners of $\Omega$. For
instance, whenever $\Omega$ is convex, $\alpha=1$; $\alpha$ may be less
than one if the domain has re-entrant corners but
$\alpha$ is close to $1$ if one of all interior angles is close to
$\pi$.

Now suppose that the $C^0$ IP solution $u_h$ is approximated by $C^0$
shape functions with polynomial degree $p \ge 2$. Then the
discretization outlined above yields the convergence rates as
discussed below.


<b>Convergence in the $C^0$ IP-norm</b>

Ideally, we would like to measure convergence in the "energy norm"
$\|D^2(u-u_h)\|$. However, this does not work because, again, the
discrete solution $u_h$ does not have two (weak) derivatives. Instead,
one can define a discrete ($C^0$ IP) seminorm that is "equivalent" to the
energy norm, as follows:
@f{align*}{
 |u_h|_{h}^2 :=
 \sum\limits_{K \in \mathbb{T}} \big|u_h\big|_{H^2(K)}^2
 +
 \sum\limits_{e \in \mathbb{F} }
 \frac{\gamma }{h_e} \left\|
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \right\|_{L^2(e)}^2.
@f}

In this seminorm, the theory in the paper mentioned above yields that we
can expect
@f{align*}{
 |u-u_h|_{h}^2 = {\cal O}(h^{p-1}),
@f}
much as one would expect given the convergence rates we know are true
for the usual discretizations of the Laplace equation.

Of course, this is true only if the exact solution is sufficiently
smooth. Indeed, if $f \in H^m(\Omega)$ with $m \ge 0$,
$u \in H^{2+\alpha}(\Omega)$ where $ 2 < 2+\alpha  \le m+4$,
then the convergence rate of the $C^0$ IP method is
$\mathcal{O}(h^{\min\{p-1, \alpha\}})$. In other words, the optimal
convergence rate can only be expected if the solution is so smooth
that $\alpha\ge p-1$; this can
only happen if (i) the domain is convex with a sufficiently smooth
boundary, and (ii) $m\ge p-3$. In practice, of course, the solution is
what it is (independent of the polynomial degree we choose), and the
last condition can then equivalently be read as saying that there is
definitely no point in choosing $p$ large if $m$ is not also
large. In other words, the only reasonably choices for $p$ are $p\le
m+3$ because larger polynomial degrees do not result in higher
convergence orders.

For the purposes of this program, we're a bit too lazy to actually
implement this equivalent seminorm -- though it's not very difficult and
would make for a good exercise. Instead, we'll simply check in the
program what the "broken" $H^2$ seminorm
@f{align*}{
 \left(|u_h|^\circ_{H^2}\right)^2
 :=
 \sum\limits_{K \in \mathbb{T}} \big|u_h\big|_{H^2(K)}^2
 =
 \sum\limits_{K \in \mathbb{T}} \big|D^2 u_h\big|_{L_2}^2
@f}
yields. The convergence rate in this norm can, from a theoretical
perspective, of course not be <i>worse</i> than the one for
$|\cdot|_h$ because it contains only a subset of the necessary terms,
but it could at least conceivably be better. It could also be the case that
we get the optimal convergence rate even though there is a bug in the
program, and that that bug would only show up in sub-optimal rates for
the additional terms present in $|\cdot|_h$. But, one might hope
that if we get the optimal rate in the broken norm and the norms
discussed below, then the program is indeed correct. The results
section will demonstrate that we obtain optimal rates in all norms
shown.


<b>Convergence in the $L_2$-norm</b>

The optimal convergence rate in the $L_2$-norm is $\mathcal{O}(h^{p+1})$
provided $p \ge 3$. More details can be found in Theorem 4.6 of
@cite Engel2002 .

The default in the program below is to choose $p=2$. In that case, the
theorem does not apply, and indeed one only gets $\mathcal{O}(h^2)$
instead of $\mathcal{O}(h^3)$ as we will show in the results section.


<b>Convergence in the $H^1$-seminorm</b>

Given that we expect
$\mathcal{O}(h^{p-1})$ in the best of cases for a norm equivalent to
the $H^2$ seminorm, and $\mathcal{O}(h^{p+1})$ for the $L_2$ norm, one
may ask about what happens in the $H^1$ seminorm that is intermediate
to the two others. A reasonable guess is that one should expect
$\mathcal{O}(h^{p})$. There is probably a paper somewhere that proves
this, but we also verify that this conjecture is experimentally true
below.



<a name="OtherBoundaryConditions"></a><h3>Other Boundary Conditions</h3>


We remark that the derivation of the $C^0$ IP method for the
biharmonic equation with other boundary conditions -- for instance,
for the first set of boundary conditions namely $u(\mathbf x) =
g(\mathbf x)$ and $\Delta u(\mathbf x)= h(\mathbf x)$ on
$\partial\Omega$ -- can be obtained with suitable modifications to
$\mathcal{A}(\cdot,\cdot)$ and $\mathcal{F}(\cdot)$ described in
the book chapter @cite Brenner2011 .


<a name="Thetestcase"></a><h3>The testcase</h3>


The last step that remains to describe is what this program solves
for. As always, a trigonometric function is both a good and a bad
choice because it does not lie in any polynomial space in which we may
seek the solution while at the same time being smoother than real
solutions typically are (here, it is in $C^\infty$ while real
solutions are typically only in $H^3$ or so on convex polygonal
domains, or somewhere between $H^2$ and $H^3$ if the domain is not
convex). But, since we don't have the means to describe solutions of
realistic problems in terms of relatively simple formulas, we just go
with the following, on the unit square for the domain $\Omega$:
@f{align*}{
  u = \sin(\pi x) \sin(\pi y).
@f}
As a consequence, we then need choose as boundary conditions the following:
@f{align*}{
  g &= u|_{\partial\Omega} = \sin(\pi x) \sin(\pi y)|_{\partial\Omega},
  \\
  j &= \frac{\partial u}{\partial\mathbf n}|_{\partial\Omega}
  \\
    &= \left.\begin{pmatrix}
                \pi\cos(\pi x) \sin(\pi y) \\
                \pi\sin(\pi x) \cos(\pi y)
             \end{pmatrix}\right|_{\partial\Omega} \cdot \mathbf n.
@f}
The right hand side is easily computed as
@f{align*}{
  f = \Delta^2 u = 4 \pi^4 \sin(\pi x) \sin(\pi y).
@f}
The program has classes `ExactSolution::Solution` and
`ExactSolution::RightHandSide` that encode this information.
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
 * The first few include files have already been used in the previous
 * example, so we will not explain their meaning here again. The principal
 * structure of the program is very similar to that of, for example, step-4
 * and so we include many of the same header files.
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/sparse_direct.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * The two most interesting header files will be these two:
 * 
 * @code
 * #include <deal.II/fe/fe_interface_values.h>
 * #include <deal.II/meshworker/mesh_loop.h>
 * @endcode
 * 
 * The first of these is responsible for providing the class FEInterfaceValues
 * that can be used to evaluate quantities such as the jump or average
 * of shape functions (or their gradients) across interfaces between cells.
 * This class will be quite useful in evaluating the penalty terms that appear
 * in the C0IP formulation.
 * 

 * 
 * 

 * 
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * #include <cmath>
 * 
 * 
 * namespace Step47
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * In the following namespace, let us define the exact solution against
 * which we will compare the numerically computed one. It has the form
 * $u(x,y) = \sin(\pi x) \sin(\pi y)$ (only the 2d case is implemented),
 * and the namespace also contains a class that corresponds to the right
 * hand side that produces this solution.
 * 
 * @code
 *   namespace ExactSolution
 *   {
 *     using numbers::PI;
 * 
 *     template <int dim>
 *     class Solution : public Function<dim>
 *     {
 *     public:
 *       static_assert(dim == 2, "Only dim==2 is implemented.");
 * 
 *       virtual double value(const Point<dim> &p,
 *                            const unsigned int /*component*/ = 0) const override
 *       {
 *         return std::sin(PI * p[0]) * std::sin(PI * p[1]);
 *       }
 * 
 *       virtual Tensor<1, dim>
 *       gradient(const Point<dim> &p,
 *                const unsigned int /*component*/ = 0) const override
 *       {
 *         Tensor<1, dim> r;
 *         r[0] = PI * std::cos(PI * p[0]) * std::sin(PI * p[1]);
 *         r[1] = PI * std::cos(PI * p[1]) * std::sin(PI * p[0]);
 *         return r;
 *       }
 * 
 *       virtual void
 *       hessian_list(const std::vector<Point<dim>> &       points,
 *                    std::vector<SymmetricTensor<2, dim>> &hessians,
 *                    const unsigned int /*component*/ = 0) const override
 *       {
 *         for (unsigned i = 0; i < points.size(); ++i)
 *           {
 *             const double x = points[i][0];
 *             const double y = points[i][1];
 * 
 *             hessians[i][0][0] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
 *             hessians[i][0][1] = PI * PI * std::cos(PI * x) * std::cos(PI * y);
 *             hessians[i][1][1] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
 *           }
 *       }
 *     };
 * 
 * 
 *     template <int dim>
 *     class RightHandSide : public Function<dim>
 *     {
 *     public:
 *       static_assert(dim == 2, "Only dim==2 is implemented");
 * 
 *       virtual double value(const Point<dim> &p,
 *                            const unsigned int /*component*/ = 0) const override
 * 
 *       {
 *         return 4 * std::pow(PI, 4.0) * std::sin(PI * p[0]) *
 *                std::sin(PI * p[1]);
 *       }
 *     };
 *   } // namespace ExactSolution
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
 * The following is the principal class of this tutorial program. It has
 * the structure of many of the other tutorial programs and there should
 * really be nothing particularly surprising about its contents or
 * the constructor that follows it.
 * 
 * @code
 *   template <int dim>
 *   class BiharmonicProblem
 *   {
 *   public:
 *     BiharmonicProblem(const unsigned int fe_degree);
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void compute_errors();
 *     void output_results(const unsigned int iteration) const;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     MappingQ<dim> mapping;
 * 
 *     FE_Q<dim>                 fe;
 *     DoFHandler<dim>           dof_handler;
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   BiharmonicProblem<dim>::BiharmonicProblem(const unsigned int fe_degree)
 *     : mapping(1)
 *     , fe(fe_degree)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * Next up are the functions that create the initial mesh (a once refined
 * unit square) and set up the constraints, vectors, and matrices on
 * each mesh. Again, both of these are essentially unchanged from many
 * previous tutorial programs.
 * 
 * @code
 *   template <int dim>
 *   void BiharmonicProblem<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0., 1.);
 *     triangulation.refine_global(1);
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Total number of cells: " << triangulation.n_cells()
 *               << std::endl;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BiharmonicProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              ExactSolution::Solution<dim>(),
 *                                              constraints);
 *     constraints.close();
 * 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, true);
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Assemblingthelinearsystem"></a> 
 * <h4>Assembling the linear system</h4>
 *   

 * 
 * The following pieces of code are more interesting. They all relate to the
 * assembly of the linear system. While assembling the cell-interior terms
 * is not of great difficulty -- that works in essence like the assembly
 * of the corresponding terms of the Laplace equation, and you have seen
 * how this works in step-4 or step-6, for example -- the difficulty
 * is with the penalty terms in the formulation. These require the evaluation
 * of gradients of shape functions at interfaces of cells. At the least,
 * one would therefore need to use two FEFaceValues objects, but if one of the
 * two sides is adaptively refined, then one actually needs an FEFaceValues
 * and one FESubfaceValues objects; one also needs to keep track which
 * shape functions live where, and finally we need to ensure that every
 * face is visited only once. All of this is a substantial overhead to the
 * logic we really want to implement (namely the penalty terms in the
 * bilinear form). As a consequence, we will make use of the
 * FEInterfaceValues class -- a helper class in deal.II that allows us
 * to abstract away the two FEFaceValues or FESubfaceValues objects and
 * directly access what we really care about: jumps, averages, etc.
 *   

 * 
 * But this doesn't yet solve our problem of having to keep track of
 * which faces we have already visited when we loop over all cells and
 * all of their faces. To make this process simpler, we use the
 * MeshWorker::mesh_loop() function that provides a simple interface
 * for this task: Based on the ideas outlined in the WorkStream
 * namespace documentation, MeshWorker::mesh_loop() requires three
 * functions that do work on cells, interior faces, and boundary
 * faces. These functions work on scratch objects for intermediate
 * results, and then copy the result of their computations into
 * copy data objects from where a copier function copies them into
 * the global matrix and right hand side objects.
 *   

 * 
 * The following structures then provide the scratch and copy objects
 * that are necessary for this approach. You may look up the WorkStream
 * namespace as well as the
 * @ref threads "Parallel computing with multiple processors"
 * module for more information on how they typically work.
 * 
 * @code
 *   template <int dim>
 *   struct ScratchData
 *   {
 *     ScratchData(const Mapping<dim> &      mapping,
 *                 const FiniteElement<dim> &fe,
 *                 const unsigned int        quadrature_degree,
 *                 const UpdateFlags         update_flags,
 *                 const UpdateFlags         interface_update_flags)
 *       : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
 *       , fe_interface_values(mapping,
 *                             fe,
 *                             QGauss<dim - 1>(quadrature_degree),
 *                             interface_update_flags)
 *     {}
 * 
 * 
 *     ScratchData(const ScratchData<dim> &scratch_data)
 *       : fe_values(scratch_data.fe_values.get_mapping(),
 *                   scratch_data.fe_values.get_fe(),
 *                   scratch_data.fe_values.get_quadrature(),
 *                   scratch_data.fe_values.get_update_flags())
 *       , fe_interface_values(scratch_data.fe_values.get_mapping(),
 *                             scratch_data.fe_values.get_fe(),
 *                             scratch_data.fe_interface_values.get_quadrature(),
 *                             scratch_data.fe_interface_values.get_update_flags())
 *     {}
 * 
 *     FEValues<dim>          fe_values;
 *     FEInterfaceValues<dim> fe_interface_values;
 *   };
 * 
 * 
 * 
 *   struct CopyData
 *   {
 *     CopyData(const unsigned int dofs_per_cell)
 *       : cell_matrix(dofs_per_cell, dofs_per_cell)
 *       , cell_rhs(dofs_per_cell)
 *       , local_dof_indices(dofs_per_cell)
 *     {}
 * 
 * 
 *     CopyData(const CopyData &) = default;
 * 
 * 
 *     struct FaceData
 *     {
 *       FullMatrix<double>                   cell_matrix;
 *       std::vector<types::global_dof_index> joint_dof_indices;
 *     };
 * 
 *     FullMatrix<double>                   cell_matrix;
 *     Vector<double>                       cell_rhs;
 *     std::vector<types::global_dof_index> local_dof_indices;
 *     std::vector<FaceData>                face_data;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * The more interesting part is where we actually assemble the linear system.
 * Fundamentally, this function has five parts:
 * - The definition of the `cell_worker` lambda function, a small
 * function that is defined within the `assemble_system()`
 * function and that will be responsible for computing the local
 * integrals on an individual cell. It will work on a copy of the
 * `ScratchData` class and put its results into the corresponding
 * `CopyData` object.
 * - The definition of the `face_worker` lambda function that does
 * the integration of all terms that live on the interfaces between
 * cells.
 * - The definition of the `boundary_worker` function that does the
 * same but for cell faces located on the boundary of the domain.
 * - The definition of the `copier` function that is responsible
 * for copying all of the data the previous three functions have
 * put into copy objects for a single cell, into the global matrix
 * and right hand side.
 *   

 * 
 * The fifth part is the one where we bring all of this together.
 *   

 * 
 * Let us go through each of these pieces necessary for the assembly
 * in turns.
 * 
 * @code
 *   template <int dim>
 *   void BiharmonicProblem<dim>::assemble_system()
 *   {
 *     using Iterator = typename DoFHandler<dim>::active_cell_iterator;
 * 
 * @endcode
 * 
 * The first piece is the `cell_worker` that does the assembly
 * on the cell interiors. It is a (lambda) function that takes
 * a cell (input), a scratch object, and a copy object (output)
 * as arguments. It looks like the assembly functions of many
 * other of the tutorial programs, or at least the body of the
 * loop over all cells.
 *     

 * 
 * The terms we integrate here are the cell contribution
 * @f{align*}{
 * A^K_{ij} = \int_K \nabla^2\varphi_i(x) : \nabla^2\varphi_j(x) dx
 * @f}
 * to the global matrix, and
 * @f{align*}{
 * f^K_i = \int_K \varphi_i(x) f(x) dx
 * @f}
 * to the right hand side vector.
 *     

 * 
 * We use the same technique as used in the assembly of step-22
 * to accelerate the function: Instead of calling
 * `fe_values.shape_hessian(i, qpoint)` in the innermost loop,
 * we create a variable `hessian_i` that evaluates this
 * value once in the loop over `i` and re-use the so-evaluated
 * value in the loop over `j`. For symmetry, we do the same with a
 * variable `hessian_j`, although it is indeed only used once and
 * we could have left the call to `fe_values.shape_hessian(j,qpoint)`
 * in the instruction that computes the scalar product between
 * the two terms.
 * 
 * @code
 *     auto cell_worker = [&](const Iterator &  cell,
 *                            ScratchData<dim> &scratch_data,
 *                            CopyData &        copy_data) {
 *       copy_data.cell_matrix = 0;
 *       copy_data.cell_rhs    = 0;
 * 
 *       FEValues<dim> &fe_values = scratch_data.fe_values;
 *       fe_values.reinit(cell);
 * 
 *       cell->get_dof_indices(copy_data.local_dof_indices);
 * 
 *       const ExactSolution::RightHandSide<dim> right_hand_side;
 * 
 *       const unsigned int dofs_per_cell =
 *         scratch_data.fe_values.get_fe().n_dofs_per_cell();
 * 
 *       for (unsigned int qpoint = 0; qpoint < fe_values.n_quadrature_points;
 *            ++qpoint)
 *         {
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const Tensor<2, dim> hessian_i =
 *                 fe_values.shape_hessian(i, qpoint);
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const Tensor<2, dim> hessian_j =
 *                     fe_values.shape_hessian(j, qpoint);
 * 
 *                   copy_data.cell_matrix(i, j) +=
 *                     scalar_product(hessian_i,   // nabla^2 phi_i(x)
 *                                    hessian_j) * // nabla^2 phi_j(x)
 *                     fe_values.JxW(qpoint);      // dx
 *                 }
 * 
 *               copy_data.cell_rhs(i) +=
 *                 fe_values.shape_value(i, qpoint) * // phi_i(x)
 *                 right_hand_side.value(
 *                   fe_values.quadrature_point(qpoint)) * // f(x)
 *                 fe_values.JxW(qpoint);                  // dx
 *             }
 *         }
 *     };
 * 
 * 
 * @endcode
 * 
 * The next building block is the one that assembles penalty terms on each
 * of the interior faces of the mesh. As described in the documentation of
 * MeshWorker::mesh_loop(), this function receives arguments that denote
 * a cell and its neighboring cell, as well as (for each of the two
 * cells) the face (and potentially sub-face) we have to integrate
 * over. Again, we also get a scratch object, and a copy object
 * for putting the results in.
 *     

 * 
 * The function has three parts itself. At the top, we initialize
 * the FEInterfaceValues object and create a new `CopyData::FaceData`
 * object to store our input in. This gets pushed to the end of the
 * `copy_data.face_data` variable. We need to do this because
 * the number of faces (or subfaces) over which we integrate for a
 * given cell differs from cell to cell, and the sizes of these
 * matrices also differ, depending on what degrees of freedom
 * are adjacent to the face or subface. As discussed in the documentation
 * of MeshWorker::mesh_loop(), the copy object is reset every time a new
 * cell is visited, so that what we push to the end of
 * `copy_data.face_data()` is really all that the later `copier` function
 * gets to see when it copies the contributions of each cell to the global
 * matrix and right hand side objects.
 * 
 * @code
 *     auto face_worker = [&](const Iterator &    cell,
 *                            const unsigned int &f,
 *                            const unsigned int &sf,
 *                            const Iterator &    ncell,
 *                            const unsigned int &nf,
 *                            const unsigned int &nsf,
 *                            ScratchData<dim> &  scratch_data,
 *                            CopyData &          copy_data) {
 *       FEInterfaceValues<dim> &fe_interface_values =
 *         scratch_data.fe_interface_values;
 *       fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);
 * 
 *       copy_data.face_data.emplace_back();
 *       CopyData::FaceData &copy_data_face = copy_data.face_data.back();
 * 
 *       copy_data_face.joint_dof_indices =
 *         fe_interface_values.get_interface_dof_indices();
 * 
 *       const unsigned int n_interface_dofs =
 *         fe_interface_values.n_current_interface_dofs();
 *       copy_data_face.cell_matrix.reinit(n_interface_dofs, n_interface_dofs);
 * 
 * @endcode
 * 
 * The second part deals with determining what the penalty
 * parameter should be. By looking at the units of the various
 * terms in the bilinear form, it is clear that the penalty has
 * to have the form $\frac{\gamma}{h_K}$ (i.e., one over length
 * scale), but it is not a priori obvious how one should choose
 * the dimension-less number $\gamma$. From the discontinuous
 * Galerkin theory for the Laplace equation, one might
 * conjecture that the right choice is $\gamma=p(p+1)$ is the
 * right choice, where $p$ is the polynomial degree of the
 * finite element used. We will discuss this choice in a bit
 * more detail in the results section of this program.
 *       

 * 
 * In the formula above, $h_K$ is the size of cell $K$. But this
 * is not quite so straightforward either: If one uses highly
 * stretched cells, then a more involved theory says that $h$
 * should be replaced by the diameter of cell $K$ normal to the
 * direction of the edge in question.  It turns out that there
 * is a function in deal.II for that. Secondly, $h_K$ may be
 * different when viewed from the two different sides of a face.
 *       

 * 
 * To stay on the safe side, we take the maximum of the two values.
 * We will note that it is possible that this computation has to be
 * further adjusted if one were to use hanging nodes resulting from
 * adaptive mesh refinement.
 * 
 * @code
 *       const unsigned int p = fe.degree;
 *       const double       gamma_over_h =
 *         std::max((1.0 * p * (p + 1) /
 *                   cell->extent_in_direction(
 *                     GeometryInfo<dim>::unit_normal_direction[f])),
 *                  (1.0 * p * (p + 1) /
 *                   ncell->extent_in_direction(
 *                     GeometryInfo<dim>::unit_normal_direction[nf])));
 * 
 * @endcode
 * 
 * Finally, and as usual, we loop over the quadrature points and
 * indices `i` and `j` to add up the contributions of this face
 * or sub-face. These are then stored in the
 * `copy_data.face_data` object created above. As for the cell
 * worker, we pull the evaluation of averages and jumps out of
 * the loops if possible, introducing local variables that store
 * these results. The assembly then only needs to use these
 * local variables in the innermost loop. Regarding the concrete
 * formula this code implements, recall that the interface terms
 * of the bilinear form were as follows:
 * @f{align*}{
 * -\sum_{e \in \mathbb{F}} \int_{e}
 * \jump{ \frac{\partial v_h}{\partial \mathbf n}}
 * \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
 * -\sum_{e \in \mathbb{F}} \int_{e}
 * \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
 * \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
 * + \sum_{e \in \mathbb{F}}
 * \frac{\gamma}{h_e}
 * \int_e
 * \jump{\frac{\partial v_h}{\partial \mathbf n}}
 * \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds.
 * @f}
 * 
 * @code
 *       for (unsigned int qpoint = 0;
 *            qpoint < fe_interface_values.n_quadrature_points;
 *            ++qpoint)
 *         {
 *           const auto &n = fe_interface_values.normal(qpoint);
 * 
 *           for (unsigned int i = 0; i < n_interface_dofs; ++i)
 *             {
 *               const double av_hessian_i_dot_n_dot_n =
 *                 (fe_interface_values.average_hessian(i, qpoint) * n * n);
 *               const double jump_grad_i_dot_n =
 *                 (fe_interface_values.jump_gradient(i, qpoint) * n);
 * 
 *               for (unsigned int j = 0; j < n_interface_dofs; ++j)
 *                 {
 *                   const double av_hessian_j_dot_n_dot_n =
 *                     (fe_interface_values.average_hessian(j, qpoint) * n * n);
 *                   const double jump_grad_j_dot_n =
 *                     (fe_interface_values.jump_gradient(j, qpoint) * n);
 * 
 *                   copy_data_face.cell_matrix(i, j) +=
 *                     (-av_hessian_i_dot_n_dot_n       // - {grad^2 v n n }
 *                        * jump_grad_j_dot_n           // [grad u n]
 *                      - av_hessian_j_dot_n_dot_n      // - {grad^2 u n n }
 *                          * jump_grad_i_dot_n         // [grad v n]
 *                      +                               // +
 *                      gamma_over_h *                  // gamma/h
 *                        jump_grad_i_dot_n *           // [grad v n]
 *                        jump_grad_j_dot_n) *          // [grad u n]
 *                     fe_interface_values.JxW(qpoint); // dx
 *                 }
 *             }
 *         }
 *     };
 * 
 * 
 * @endcode
 * 
 * The third piece is to do the same kind of assembly for faces that
 * are at the boundary. The idea is the same as above, of course,
 * with only the difference that there are now penalty terms that
 * also go into the right hand side.
 *     

 * 
 * As before, the first part of the function simply sets up some
 * helper objects:
 * 
 * @code
 *     auto boundary_worker = [&](const Iterator &    cell,
 *                                const unsigned int &face_no,
 *                                ScratchData<dim> &  scratch_data,
 *                                CopyData &          copy_data) {
 *       FEInterfaceValues<dim> &fe_interface_values =
 *         scratch_data.fe_interface_values;
 *       fe_interface_values.reinit(cell, face_no);
 *       const auto &q_points = fe_interface_values.get_quadrature_points();
 * 
 *       copy_data.face_data.emplace_back();
 *       CopyData::FaceData &copy_data_face = copy_data.face_data.back();
 * 
 *       const unsigned int n_dofs =
 *         fe_interface_values.n_current_interface_dofs();
 *       copy_data_face.joint_dof_indices =
 *         fe_interface_values.get_interface_dof_indices();
 * 
 *       copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);
 * 
 *       const std::vector<double> &JxW = fe_interface_values.get_JxW_values();
 *       const std::vector<Tensor<1, dim>> &normals =
 *         fe_interface_values.get_normal_vectors();
 * 
 * 
 *       const ExactSolution::Solution<dim> exact_solution;
 *       std::vector<Tensor<1, dim>>        exact_gradients(q_points.size());
 *       exact_solution.gradient_list(q_points, exact_gradients);
 * 
 * 
 * @endcode
 * 
 * Positively, because we now only deal with one cell adjacent to the
 * face (as we are on the boundary), the computation of the penalty
 * factor $\gamma$ is substantially simpler:
 * 
 * @code
 *       const unsigned int p = fe.degree;
 *       const double       gamma_over_h =
 *         (1.0 * p * (p + 1) /
 *          cell->extent_in_direction(
 *            GeometryInfo<dim>::unit_normal_direction[face_no]));
 * 
 * @endcode
 * 
 * The third piece is the assembly of terms. This is now
 * slightly more involved since these contains both terms for
 * the matrix and for the right hand side. The former is exactly
 * the same as for the interior faces stated above if one just
 * defines the jump and average appropriately (which is what the
 * FEInterfaceValues class does). The latter requires us to
 * evaluate the boundary conditions $j(\mathbf x)$, which in the
 * current case (where we know the exact solution) we compute
 * from $j(\mathbf x) = \frac{\partial u(\mathbf x)}{\partial
 * {\mathbf n}}$. The term to be added to the right hand side
 * vector is then
 * $\frac{\gamma}{h_e}\int_e
 * \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds$.
 * 
 * @code
 *       for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
 *         {
 *           const auto &n = normals[qpoint];
 * 
 *           for (unsigned int i = 0; i < n_dofs; ++i)
 *             {
 *               const double av_hessian_i_dot_n_dot_n =
 *                 (fe_interface_values.average_hessian(i, qpoint) * n * n);
 *               const double jump_grad_i_dot_n =
 *                 (fe_interface_values.jump_gradient(i, qpoint) * n);
 * 
 *               for (unsigned int j = 0; j < n_dofs; ++j)
 *                 {
 *                   const double av_hessian_j_dot_n_dot_n =
 *                     (fe_interface_values.average_hessian(j, qpoint) * n * n);
 *                   const double jump_grad_j_dot_n =
 *                     (fe_interface_values.jump_gradient(j, qpoint) * n);
 * 
 *                   copy_data_face.cell_matrix(i, j) +=
 *                     (-av_hessian_i_dot_n_dot_n  // - {grad^2 v n n}
 *                        * jump_grad_j_dot_n      //   [grad u n]
 *                                                 
 *                      - av_hessian_j_dot_n_dot_n // - {grad^2 u n n}
 *                          * jump_grad_i_dot_n    //   [grad v n]
 *                                                 
 *                      + gamma_over_h             //  gamma/h
 *                          * jump_grad_i_dot_n    // [grad v n]
 *                          * jump_grad_j_dot_n    // [grad u n]
 *                      ) *
 *                     JxW[qpoint]; // dx
 *                 }
 * 
 *               copy_data.cell_rhs(i) +=
 *                 (-av_hessian_i_dot_n_dot_n *       // - {grad^2 v n n }
 *                    (exact_gradients[qpoint] * n)   //   (grad u_exact . n)
 *                  +                                 // +
 *                  gamma_over_h                      //  gamma/h
 *                    * jump_grad_i_dot_n             // [grad v n]
 *                    * (exact_gradients[qpoint] * n) // (grad u_exact . n)
 *                  ) *
 *                 JxW[qpoint]; // dx
 *             }
 *         }
 *     };
 * 
 * @endcode
 * 
 * Part 4 is a small function that copies the data produced by the
 * cell, interior, and boundary face assemblers above into the
 * global matrix and right hand side vector. There really is not
 * very much to do here: We distribute the cell matrix and right
 * hand side contributions as we have done in almost all of the
 * other tutorial programs using the constraints objects. We then
 * also have to do the same for the face matrix contributions
 * that have gained content for the faces (interior and boundary)
 * and that the `face_worker` and `boundary_worker` have added
 * to the `copy_data.face_data` array.
 * 
 * @code
 *     auto copier = [&](const CopyData &copy_data) {
 *       constraints.distribute_local_to_global(copy_data.cell_matrix,
 *                                              copy_data.cell_rhs,
 *                                              copy_data.local_dof_indices,
 *                                              system_matrix,
 *                                              system_rhs);
 * 
 *       for (auto &cdf : copy_data.face_data)
 *         {
 *           constraints.distribute_local_to_global(cdf.cell_matrix,
 *                                                  cdf.joint_dof_indices,
 *                                                  system_matrix);
 *         }
 *     };
 * 
 * 
 * @endcode
 * 
 * Having set all of this up, what remains is to just create a scratch
 * and copy data object and call the MeshWorker::mesh_loop() function
 * that then goes over all cells and faces, calls the respective workers
 * on them, and then the copier function that puts things into the
 * global matrix and right hand side. As an additional benefit,
 * MeshWorker::mesh_loop() does all of this in parallel, using
 * as many processor cores as your machine happens to have.
 * 
 * @code
 *     const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
 *     ScratchData<dim>   scratch_data(mapping,
 *                                   fe,
 *                                   n_gauss_points,
 *                                   update_values | update_gradients |
 *                                     update_hessians | update_quadrature_points |
 *                                     update_JxW_values,
 *                                   update_values | update_gradients |
 *                                     update_hessians | update_quadrature_points |
 *                                     update_JxW_values | update_normal_vectors);
 *     CopyData           copy_data(dof_handler.get_fe().n_dofs_per_cell());
 *     MeshWorker::mesh_loop(dof_handler.begin_active(),
 *                           dof_handler.end(),
 *                           cell_worker,
 *                           copier,
 *                           scratch_data,
 *                           copy_data,
 *                           MeshWorker::assemble_own_cells |
 *                             MeshWorker::assemble_boundary_faces |
 *                             MeshWorker::assemble_own_interior_faces_once,
 *                           boundary_worker,
 *                           face_worker);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solvingthelinearsystemandpostprocessing"></a> 
 * <h4>Solving the linear system and postprocessing</h4>
 *   

 * 
 * The show is essentially over at this point: The remaining functions are
 * not overly interesting or novel. The first one simply uses a direct
 * solver to solve the linear system (see also step-29):
 * 
 * @code
 *   template <int dim>
 *   void BiharmonicProblem<dim>::solve()
 *   {
 *     std::cout << "   Solving system..." << std::endl;
 * 
 *     SparseDirectUMFPACK A_direct;
 *     A_direct.initialize(system_matrix);
 *     A_direct.vmult(solution, system_rhs);
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function evaluates the error between the computed solution
 * and the exact solution (which is known here because we have chosen
 * the right hand side and boundary values in a way so that we know
 * the corresponding solution). In the first two code blocks below,
 * we compute the error in the $L_2$ norm and the $H^1$ semi-norm.
 * 
 * @code
 *   template <int dim>
 *   void BiharmonicProblem<dim>::compute_errors()
 *   {
 *     {
 *       Vector<float> norm_per_cell(triangulation.n_active_cells());
 *       VectorTools::integrate_difference(mapping,
 *                                         dof_handler,
 *                                         solution,
 *                                         ExactSolution::Solution<dim>(),
 *                                         norm_per_cell,
 *                                         QGauss<dim>(fe.degree + 2),
 *                                         VectorTools::L2_norm);
 *       const double error_norm =
 *         VectorTools::compute_global_error(triangulation,
 *                                           norm_per_cell,
 *                                           VectorTools::L2_norm);
 *       std::cout << "   Error in the L2 norm           :     " << error_norm
 *                 << std::endl;
 *     }
 * 
 *     {
 *       Vector<float> norm_per_cell(triangulation.n_active_cells());
 *       VectorTools::integrate_difference(mapping,
 *                                         dof_handler,
 *                                         solution,
 *                                         ExactSolution::Solution<dim>(),
 *                                         norm_per_cell,
 *                                         QGauss<dim>(fe.degree + 2),
 *                                         VectorTools::H1_seminorm);
 *       const double error_norm =
 *         VectorTools::compute_global_error(triangulation,
 *                                           norm_per_cell,
 *                                           VectorTools::H1_seminorm);
 *       std::cout << "   Error in the H1 seminorm       : " << error_norm
 *                 << std::endl;
 *     }
 * 
 * @endcode
 * 
 * Now also compute an approximation to the $H^2$ seminorm error. The actual
 * $H^2$ seminorm would require us to integrate second derivatives of the
 * solution $u_h$, but given the Lagrange shape functions we use, $u_h$ of
 * course has kinks at the interfaces between cells, and consequently second
 * derivatives are singular at interfaces. As a consequence, we really only
 * integrate over the interior of cells and ignore the interface
 * contributions. This is *not* an equivalent norm to the energy norm for
 * the problem, but still gives us an idea of how fast the error converges.
 *     

 * 
 * We note that one could address this issue by defining a norm that
 * is equivalent to the energy norm. This would involve adding up not
 * only the integrals over cell interiors as we do below, but also adding
 * penalty terms for the jump of the derivative of $u_h$ across interfaces,
 * with an appropriate scaling of the two kinds of terms. We will leave
 * this for later work.
 * 
 * @code
 *     {
 *       const QGauss<dim>            quadrature_formula(fe.degree + 2);
 *       ExactSolution::Solution<dim> exact_solution;
 *       Vector<double> error_per_cell(triangulation.n_active_cells());
 * 
 *       FEValues<dim> fe_values(mapping,
 *                               fe,
 *                               quadrature_formula,
 *                               update_values | update_hessians |
 *                                 update_quadrature_points | update_JxW_values);
 * 
 *       FEValuesExtractors::Scalar scalar(0);
 *       const unsigned int         n_q_points = quadrature_formula.size();
 * 
 *       std::vector<SymmetricTensor<2, dim>> exact_hessians(n_q_points);
 *       std::vector<Tensor<2, dim>>          hessians(n_q_points);
 *       for (auto &cell : dof_handler.active_cell_iterators())
 *         {
 *           fe_values.reinit(cell);
 *           fe_values[scalar].get_function_hessians(solution, hessians);
 *           exact_solution.hessian_list(fe_values.get_quadrature_points(),
 *                                       exact_hessians);
 * 
 *           double local_error = 0;
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             {
 *               local_error +=
 *                 ((exact_hessians[q_point] - hessians[q_point]).norm_square() *
 *                  fe_values.JxW(q_point));
 *             }
 *           error_per_cell[cell->active_cell_index()] = std::sqrt(local_error);
 *         }
 * 
 *       const double error_norm = error_per_cell.l2_norm();
 *       std::cout << "   Error in the broken H2 seminorm: " << error_norm
 *                 << std::endl;
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Equally uninteresting is the function that generates graphical output.
 * It looks exactly like the one in step-6, for example.
 * 
 * @code
 *   template <int dim>
 *   void
 *   BiharmonicProblem<dim>::output_results(const unsigned int iteration) const
 *   {
 *     std::cout << "   Writing graphical output..." << std::endl;
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches();
 * 
 *     const std::string filename =
 *       ("output_" + Utilities::int_to_string(iteration, 6) + ".vtu");
 *     std::ofstream output_vtu(filename);
 *     data_out.write_vtu(output_vtu);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The same is true for the `run()` function: Just like in previous
 * programs.
 * 
 * @code
 *   template <int dim>
 *   void BiharmonicProblem<dim>::run()
 *   {
 *     make_grid();
 * 
 *     const unsigned int n_cycles = 4;
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         std::cout << "Cycle " << cycle << " of " << n_cycles << std::endl;
 * 
 *         triangulation.refine_global(1);
 *         setup_system();
 * 
 *         assemble_system();
 *         solve();
 * 
 *         output_results(cycle);
 * 
 *         compute_errors();
 *         std::cout << std::endl;
 *       }
 *   }
 * } // namespace Step47
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
 * Finally for the `main()` function. There is, again, not very much to see
 * here: It looks like the ones in previous tutorial programs. There
 * is a variable that allows selecting the polynomial degree of the element
 * we want to use for solving the equation. Because the C0IP formulation
 * we use requires the element degree to be at least two, we check with
 * an assertion that whatever one sets for the polynomial degree actually
 * makes sense.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step47;
 * 
 *       const unsigned int fe_degree = 2;
 *       Assert(fe_degree >= 2,
 *              ExcMessage("The C0IP formulation for the biharmonic problem "
 *                         "only works if one uses elements of polynomial "
 *                         "degree at least 2."));
 * 
 *       BiharmonicProblem<2> biharmonic_problem(fe_degree);
 *       biharmonic_problem.run();
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


We run the program with right hand side and boundary values as
discussed in the introduction. These will produce the
solution $u = \sin(\pi x) \sin(\pi y)$ on the domain $\Omega = (0,1)^2$.
We test this setup using $Q_2$, $Q_3$, and $Q_4$ elements, which one can
change via the `fe_degree` variable in the `main()` function. With mesh
refinement, the $L_2$ convergence rates, $H^1$-seminorm rate,
and $H^2$-seminorm convergence of $u$
should then be around 2, 2, 1 for $Q_2$ (with the $L_2$ norm
sub-optimal as discussed in the introduction); 4, 3, 2 for
$Q_3$; and 5, 4, 3 for $Q_4$, respectively.

From the literature, it is not immediately clear what
the penalty parameter $\gamma$ should be. For example,
@cite Brenner2009 state that it needs to be larger than one, and
choose $\gamma=5$. The FEniCS/Dolphin tutorial chooses it as
$\gamma=8$, see
https://fenicsproject.org/docs/dolfin/1.6.0/python/demo/documented/biharmonic/python/documentation.html
. @cite Wells2007 uses a value for $\gamma$ larger than the
number of edges belonging to an element for Kirchhoff plates (see
their Section 4.2). This suggests that maybe
$\gamma = 1$, $2$, are too small; on the other hand, a value
$p(p+1)$ would be reasonable,
where $p$ is the degree of polynomials. The last of these choices is
the one one would expect to work by comparing
to the discontinuous Galerkin formulations for the Laplace equation
(see, for example, the discussions in step-39 and step-74),
and it will turn out to also work here.
But we should check what value of $\gamma$ is right, and we will do so
below; changing $\gamma$ is easy in the two `face_worker` and
`boundary_worker` functions defined in `assemble_system()`.


<a name="TestresultsoniQsub2subiwithigammapp1i"></a><h3>Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = p(p+1)</i> </h3>


We run the code with differently refined meshes
and get the following convergence rates.

<table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   8.780e-03 </td><td>       </td><td>  7.095e-02   </td><td>           </td><td>  1.645 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   3.515e-03   </td><td>  1.32 </td><td> 2.174e-02  </td><td>     1.70     </td><td> 8.121e-01  </td><td>  1.018  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   1.103e-03   </td><td>  1.67   </td><td> 6.106e-03    </td><td>  1.83        </td><td>   4.015e-01 </td><td> 1.016  </td>
  </tr>
  <tr>
   <td>   5                  </td><td>  3.084e-04  </td><td>  1.83   </td><td>  1.622e-03   </td><td>    1.91        </td><td> 1.993e-01 </td><td>  1.010   </td>
  </tr>
</table>
We can see that the $L_2$ convergence rates are around 2,
$H^1$-seminorm convergence rates are around 2,
and $H^2$-seminorm convergence rates are around 1. The latter two
match the theoretically expected rates; for the former, we have no
theorem but are not surprised that it is sub-optimal given the remark
in the introduction.


<a name="TestresultsoniQsub3subiwithigammapp1i"></a><h3>Test results on <i>Q<sub>3</sub></i> with <i>&gamma; = p(p+1)</i> </h3>



<table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    2.045e-04 </td><td>       </td><td>   4.402e-03   </td><td>           </td><td> 1.641e-01 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   1.312e-05   </td><td> 3.96  </td><td>  5.537e-04  </td><td>   2.99     </td><td> 4.096e-02 </td><td>  2.00  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   8.239e-07 </td><td>  3.99  </td><td> 6.904e-05   </td><td> 3.00     </td><td> 1.023e-02 </td><td> 2.00 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   5.158e-08  </td><td>  3.99 </td><td> 8.621e-06 </td><td>  3.00      </td><td> 2.558e-03  </td><td>  2.00  </td>
  </tr>
</table>
We can see that the $L_2$ convergence rates are around 4,
$H^1$-seminorm convergence rates are around 3,
and $H^2$-seminorm convergence rates are around 2.
This, of course, matches our theoretical expectations.


<a name="TestresultsoniQsub4subiwithigammapp1i"></a><h3>Test results on <i>Q<sub>4</sub></i> with <i>&gamma; = p(p+1)</i> </h3>


<table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    6.510e-06 </td><td>       </td><td> 2.215e-04   </td><td>           </td><td>  1.275e-02 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   2.679e-07  </td><td>  4.60  </td><td> 1.569e-05  </td><td>   3.81    </td><td> 1.496e-03 </td><td>  3.09  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   9.404e-09  </td><td> 4.83   </td><td> 1.040e-06    </td><td> 3.91       </td><td> 1.774e-04 </td><td> 3.07 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   7.943e-10 </td><td>  3.56  </td><td>   6.693e-08 </td><td> 3.95     </td><td> 2.150e-05  </td><td> 3.04    </td>
  </tr>
</table>
We can see that the $L_2$ norm convergence rates are around 5,
$H^1$-seminorm convergence rates are around 4,
and $H^2$-seminorm convergence rates are around 3.
On the finest mesh, the $L_2$ norm convergence rate
is much smaller than our theoretical expectations
because the linear solver becomes the limiting factor due
to round-off. Of course the $L_2$ error is also very small already in
that case.


<a name="TestresultsoniQsub2subiwithigamma1i"></a><h3>Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = 1</i> </h3>


For comparison with the results above, let us now also consider the
case where we simply choose $\gamma=1$:

<table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   7.350e-02 </td><td>       </td><td>   7.323e-01   </td><td>           </td><td> 10.343 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   6.798e-03   </td><td> 3.43  </td><td> 1.716e-01   </td><td>   2.09    </td><td>4.836 </td><td>  1.09 </td>
  </tr>
  <tr>
   <td>   4                  </td><td>  9.669e-04   </td><td> 2.81   </td><td> 6.436e-02    </td><td> 1.41      </td><td>  3.590 </td><td> 0.430 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   1.755e-04 </td><td> 2.46 </td><td>  2.831e-02  </td><td>    1.18      </td><td>3.144  </td><td>  0.19  </td>
  </tr>
</table>
Although $L_2$ norm convergence rates of $u$ more or less
follows the theoretical expectations,
the $H^1$-seminorm and $H^2$-seminorm do not seem to converge as expected.
Comparing results from $\gamma = 1$ and $\gamma = p(p+1)$, it is clear that
$\gamma = p(p+1)$ is a better penalty.
Given that $\gamma=1$ is already too small for $Q_2$ elements, it may not be surprising that if one repeated the
experiment with a $Q_3$ element, the results are even more disappointing: One again only obtains convergence
rates of 2, 1, zero -- i.e., no better than for the $Q_2$ element (although the errors are smaller in magnitude).
Maybe surprisingly, however, one obtains more or less the expected convergence orders when using $Q_4$
elements. Regardless, this uncertainty suggests that $\gamma=1$ is at best a risky choice, and at worst an
unreliable one and that we should choose $\gamma$ larger.


<a name="TestresultsoniQsub2subiwithigamma2i"></a><h3>Test results on <i>Q<sub>2</sub></i> with <i>&gamma; = 2</i> </h3>


Since $\gamma=1$ is clearly too small, one might conjecture that
$\gamma=2$ might actually work better. Here is what one obtains in
that case:

<table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   4.133e-02 </td><td>       </td><td>  2.517e-01   </td><td>           </td><td> 3.056 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>  6.500e-03   </td><td>2.66  </td><td> 5.916e-02  </td><td>  2.08    </td><td>1.444 </td><td>  1.08 </td>
  </tr>
  <tr>
   <td>   4                  </td><td> 6.780e-04   </td><td> 3.26  </td><td> 1.203e-02    </td><td> 2.296      </td><td> 6.151e-01 </td><td> 1.231 </td>
  </tr>
  <tr>
   <td>   5                  </td><td> 1.622e-04 </td><td> 2.06 </td><td>  2.448e-03  </td><td>   2.297     </td><td> 2.618e-01  </td><td> 1.232  </td>
  </tr>
</table>
In this case, the convergence rates more or less follow the
theoretical expectations, but, compared to the results from $\gamma =
p(p+1)$, are more variable.
Again, we could repeat this kind of experiment for $Q_3$ and $Q_4$ elements. In both cases, we will find that we
obtain roughly the expected convergence rates. Of more interest may then be to compare the absolute
size of the errors. While in the table above, for the $Q_2$ case, the errors on the finest grid are comparable between
the $\gamma=p(p+1)$ and $\gamma=2$ case, for $Q_3$ the errors are substantially larger for $\gamma=2$ than for
$\gamma=p(p+1)$. The same is true for the $Q_4$ case.


<a name="Conclusionsforthechoiceofthepenaltyparameter"></a><h3> Conclusions for the choice of the penalty parameter </h3>


The conclusions for which of the "reasonable" choices one should use for the penalty parameter
is that $\gamma=p(p+1)$ yields the expected results. It is, consequently, what the code
uses as currently written.


<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


There are a number of obvious extensions to this program that would
make sense:

- The program uses a square domain and a uniform mesh. Real problems
  don't come this way, and one should verify convergence also on
  domains with other shapes and, in particular, curved boundaries. One
  may also be interested in resolving areas of less regularity by
  using adaptive mesh refinement.

- From a more theoretical perspective, the convergence results above
  only used the "broken" $H^2$ seminorm $|\cdot|^\circ_{H^2}$ instead
  of the "equivalent" norm $|\cdot|_h$. This is good enough to
  convince ourselves that the program isn't fundamentally
  broken. However, it might be interesting to measure the error in the
  actual norm for which we have theoretical results. Implementing this
  addition should not be overly difficult using, for example, the
  FEInterfaceValues class combined with MeshWorker::mesh_loop() in the
  same spirit as we used for the assembly of the linear system.


<a name="Derivationforthesimplysupportedplates"></a>  <h4> Derivation for the simply supported plates </h4>


  Similar to the "clamped" boundary condition addressed in the implementation,
  we will derive the $C^0$ IP finite element scheme for simply supported plates:
  @f{align*}{
    \Delta^2 u(\mathbf x) &= f(\mathbf x)
    \qquad \qquad &&\forall \mathbf x \in \Omega,
    u(\mathbf x) &= g(\mathbf x) \qquad \qquad
    &&\forall \mathbf x \in \partial\Omega, \\
    \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
    &&\forall \mathbf x \in \partial\Omega.
  @f}
  We multiply the biharmonic equation by the test function $v_h$ and integrate over $ K $ and get:
  @f{align*}{
    \int_K v_h (\Delta^2 u_h)
     &= \int_K (D^2 v_h) : (D^2 u_h)
       + \int_{\partial K} v_h \frac{\partial (\Delta u_h)}{\partial \mathbf{n}}
       -\int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}).
  @f}

  Summing up over all cells $K \in  \mathbb{T}$,since normal directions of $\Delta u_h$ are pointing at
  opposite directions on each interior edge shared by two cells and $v_h = 0$ on $\partial \Omega$,
  @f{align*}{
  \sum_{K \in \mathbb{T}} \int_{\partial K} v_h \frac{\partial (\Delta u_h)}{\partial \mathbf{n}} = 0,
  @f}
  and by the definition of jump over cell interfaces,
  @f{align*}{
  -\sum_{K \in \mathbb{T}} \int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}) = -\sum_{e \in \mathbb{F}} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}}).
  @f}
  We separate interior faces and boundary faces of the domain,
  @f{align*}{
  -\sum_{K \in \mathbb{T}} \int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}) = -\sum_{e \in \mathbb{F}^i} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}})
  - \sum_{e \in \partial \Omega} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} h,
  @f}
  where $\mathbb{F}^i$ is the set of interior faces.
  This leads us to
  @f{align*}{
  \sum_{K \in \mathbb{T}} \int_K (D^2 v_h) : (D^2 u_h) \ dx - \sum_{e \in \mathbb{F}^i} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}}) \ ds
  = \sum_{K \in \mathbb{T}}\int_{K} v_h f  \ dx + \sum_{e\subset\partial\Omega} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} h \ ds.
  @f}

  In order to symmetrize and stabilize the discrete problem,
  we add symmetrization and stabilization term.
  We finally get the $C^0$ IP finite element scheme for the biharmonic equation:
  find $u_h$ such that $u_h =g$ on $\partial \Omega$ and
  @f{align*}{
  \mathcal{A}(v_h,u_h)&=\mathcal{F}(v_h) \quad \text{holds for all test functions } v_h,
  @f}
  where
  @f{align*}{
  \mathcal{A}(v_h,u_h):=&\sum_{K \in \mathbb{T}}\int_K D^2v_h:D^2u_h \ dx
  \\
  &
   -\sum_{e \in \mathbb{F}^i} \int_{e}
    \jump{\frac{\partial v_h}{\partial \mathbf n}}
    \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
   -\sum_{e \in \mathbb{F}^i} \int_{e}
   \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
   \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
  \\
  &+ \sum_{e \in \mathbb{F}^i}
   \frac{\gamma}{h_e}
   \int_e
   \jump{\frac{\partial v_h}{\partial \mathbf n}}
   \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds,
  @f}
  and
  @f{align*}{
  \mathcal{F}(v_h)&:=\sum_{K \in \mathbb{T}}\int_{K} v_h f \ dx
  +
  \sum_{e\subset\partial\Omega}
  \int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} h \ ds.
  @f}
  The implementation of this boundary case is similar to the "clamped" version
  except that `boundary_worker` is no longer needed for system assembling
  and the right hand side is changed according to the formulation.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-47.cc"
*/
