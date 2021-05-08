  /**  @page step_27 The step-27 tutorial program  。 

本教程取决于  step-6  。

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Finiteelementcollections">Finite element collections</a><a href="#Finiteelementcollections">Finite element collections</a>
        <li><a href="#TheDoFHandlerclassinihpimodeassociatingcellswithfiniteelementsandconstraints">The DoFHandler class in <i>hp</i><a href="#TheDoFHandlerclassinihpimodeassociatingcellswithfiniteelementsandconstraints">The DoFHandler class in <i>hp</i>-mode, associating cells with finite elements, and constraints</a><i>hp</i>-mode, associating cells with finite elements, and constraints</a>
        <li><a href="#Assemblingmatricesandvectorswithhpobjects">Assembling matrices and vectors with hp-objects</a><a href="#Assemblingmatricesandvectorswithhpobjects">Assembling matrices and vectors with hp-objects</a>
        <li><a href="#Asimpleindicatorforhprefinementandestimatingsmoothness">A simple indicator for hp-refinement and estimating smoothness</a><a href="#Asimpleindicatorforhprefinementandestimatingsmoothness">A simple indicator for hp-refinement and estimating smoothness</a>
      <ul>
        <li><a href="#Theidea">The idea</a><a href="#Theidea">The idea</a>
        <li><a href="#Whatwehavetodo">What we have to do</a><a href="#Whatwehavetodo">What we have to do</a>
        <li><a href="#Compensatingforanisotropy">Compensating for anisotropy</a><a href="#Compensatingforanisotropy">Compensating for anisotropy</a>
        <li><a href="#Questionsaboutcellsizes">Questions about cell sizes</a><a href="#Questionsaboutcellsizes">Questions about cell sizes</a>
      </ul>
        <li><a href="#Complicationswithlinearsystemsforhpdiscretizations">Complications with linear systems for hp-discretizations</a><a href="#Complicationswithlinearsystemsforhpdiscretizations">Complications with linear systems for hp-discretizations</a>
      <ul>
        <li><a href="#Creatingthesparsitypattern">Creating the sparsity pattern</a><a href="#Creatingthesparsitypattern">Creating the sparsity pattern</a>
        <li><a href="#Eliminatingconstraineddegreesoffreedom">Eliminating constrained degrees of freedom</a><a href="#Eliminatingconstraineddegreesoffreedom">Eliminating constrained degrees of freedom</a>
      </ul>
        <li><a href="#Thetestcase">The test case</a><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#Themainclass">The main class</a> ]<a href="#Themainclass">The main class</a>
        <li><a href="#Equationdata">Equation data</a><a href="#Equationdata">Equation data</a>
        <li><a href="#Implementationofthemainclass">Implementation of the main class</a><a href="#Implementationofthemainclass">Implementation of the main class</a>
      <ul>
        <li><a href="#LaplaceProblemLaplaceProblemconstructor">LaplaceProblem::LaplaceProblem constructor</a><a href="#LaplaceProblemLaplaceProblemconstructor">LaplaceProblem::LaplaceProblem constructor</a>
        <li><a href="#LaplaceProblemLaplaceProblemdestructor">LaplaceProblem::~LaplaceProblem destructor</a><a href="#LaplaceProblemLaplaceProblemdestructor">LaplaceProblem::~LaplaceProblem destructor</a>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblempostprocess">LaplaceProblem::postprocess</a><a href="#LaplaceProblempostprocess">LaplaceProblem::postprocess</a>
        <li><a href="#LaplaceProblemcreate_coarse_grid">LaplaceProblem::create_coarse_grid</a><a href="#LaplaceProblemcreate_coarse_grid">LaplaceProblem::create_coarse_grid</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main function</a><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Differenthpdecisionstrategies">Different hp-decision strategies</a><a href="#Differenthpdecisionstrategies">Different hp-decision strategies</a>
        <li><a href="#Parallelhpadaptivefiniteelements">Parallel hp-adaptive finite elements</a><a href="#Parallelhpadaptivefiniteelements">Parallel hp-adaptive finite elements</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly 

<a name="Intro"></a> <a name="Introduction"></a><h1>Introduction</h1>


本教程程序试图展示如何使用 $hp$ --deal.II的有限元方法。它解决的是拉普拉斯方程，因此只建立在前几个教程程序的基础上，特别是建立在 step-4 的尺寸独立编程和 step-6 的自适应网格细化上。

 $hp$ 有限元方法是由Babu&scaron;ka和Guo在20世纪80年代初提出的，作为(i)网格细化（即减少有限元计算中的网格参数 $h$ ）或(ii)增加用于形状函数的多项式程度 $p$ 的替代方法。它是基于这样的观察：如果解足够光滑，增加形状函数的多项式程度可以减少近似误差。另一方面，众所周知，即使对于一般良好的椭圆问题，在边界附近、角落或系数不连续的地方也不能保证较高的规则度；因此，在这些地方不能通过增加多项式度来改善近似  $p$  ，而只能通过细化网格，即减少网格大小  $h$  。这些减少误差的不同手段导致了 $hp$ 有限元的概念，即在解足够光滑的地方，近似的有限元空间被调整为具有高多项式度 $p$ ，而在解缺乏规则性的地方，网格宽度 $h$ 被减少。在关于这种方法的第一篇论文中已经意识到， $hp$ --无限元可以是一个强有力的工具，它可以保证误差不仅以自由度数的某个负数减少，而且实际上是以指数形式减少。

为了实现这种方法，我们需要一些高于一般有限元程序所需的东西，特别是高于我们在通往 step-6 的教程程序中介绍的东西。特别是，我们将不得不讨论以下几个方面。  <ul>   <li>  我们现在不是在所有单元上使用同一个有限元，而是希望有一个有限元对象的集合，并将每个单元与这个集合中的一个对象相关联。 </li>   

    <li>  然后必须根据与该特定单元相关的有限元在每个单元上分配自由度。约束条件的生成方式与悬挂节点的生成方式相同，但我们现在还必须处理两个相邻单元分配不同有限元的情况。 </li>   

    <li>  我们将需要能够汇集单元和面对全局矩阵和右手向量的贡献。 </li>   

    <li>  在解决了所产生的线性系统后，我们将想分析该解决方案。特别是，我们要计算误差指标，告诉我们是否应该细化一个给定的单元和/或是否应该增加对其使用的形状函数的多项式程度。 </li>   </ul>  。 

我们将在本介绍的以下几个小节中讨论所有这些方面。我们不会感到惊讶的是，这些任务中的大部分已经被deal.II提供的功能很好地支持了，我们只需要提供程序应该做什么的逻辑，而不是确切地说明这一切将如何发生。

在deal.II中， $hp$ 的功能主要被打包到hp-namespace中。这个命名空间提供了处理 $hp$ 分解、集合矩阵和向量以及其他任务的类。我们将在下文中进一步了解它们中的许多。此外，DoFTools和VectorTools命名空间中的大多数函数除了接受非 $hp$ 的对象外，还接受 $hp$ 的对象。许多 $hp$ 的实现也在 @ref hp 文档模块和那里的链接中讨论。

也许值得在这第一部分介绍的最后给出一个稍大的观点。  $hp$  -功能已经在许多不同的有限元软件包中实现（例如，见 @ref hp_paper  "hp-paper "中引用的参考文献列表）。然而，总的来说，这些软件包大多只在(i)2d情况下和/或(ii)不连续的Galerkin方法中实现了它。后者是一个重要的简化，因为根据定义，不连续的有限元不要求单元之间的面的连续性，因此，当不同的多项式程度的有限元在一个共同的面相遇时，不需要特殊处理。相比之下，deal.II实现了最普遍的情况，即它允许1d、2d和3d的连续和不连续元素，并自动处理由此产生的复杂性。特别是，它可以处理不同程度的元素在一个面或边缘相遇时的约束（类似于悬挂节点约束）的计算。许多必要的算法和数据结构技术在 @ref hp_paper "hp-paper "中为那些对这些细节感兴趣的人进行了描述。

我们希望提供这样一个通用的实现方法将有助于进一步探索 $hp$ -方法的潜力。




<a name="Finiteelementcollections"></a><h3>Finite element collections</h3>


现在再来谈谈如何在deal.II中使用 $hp$ 的功能的细节。我们必须处理的第一个方面是，现在我们不再只有一个用于所有单元的单一有限元，而是有许多不同的单元可以选择使用。为此，deal.II引入了<i>finite element
collection</i>的概念，在 hp::FECollection. 类中实现。实质上，这样的集合就像 <code>std::vector@<FiniteElement@></code> 类型的对象，但多了一些小功能和内存管理，更适合手头的工作。正如我们以后所看到的，我们还将使用类似的正交集合，以及&mdash; 尽管我们在这里没有使用它们&mdash; 还有映射集合的概念。所有这些类在 @ref hpcollection 概述中都有描述。

在这个教程程序中，我们将使用2到7阶（在2D中）或2到5阶（在3D中）的连续Lagrange元素。然后可以按如下方式创建所用元素的集合。

@code
  hp::FECollection<dim> fe_collection;
  for (unsigned int degree = 2; degree <= max_degree; ++degree)
    fe_collection.push_back(FE_Q<dim>(degree));
@endcode 






<a name="TheDoFHandlerclassinihpimodeassociatingcellswithfiniteelementsandconstraints"></a><h3>The DoFHandler class in <i>hp</i><i>hp</i>-mode, associating cells with finite elements, and constraints</h3> 


我们要考虑的下一个任务是如何处理我们要使用的有限元对象的列表。在以前的教程程序中，从 step-2 开始，我们已经看到DoFHandler类负责在网格（由三角形对象描述）和有限元之间建立联系，为网格的每个顶点、面、边和单元分配正确的自由度数。

这里的情况有点复杂，因为我们不只是有一个单一的有限元对象，而是可能想在不同的单元上使用不同的元素。因此，我们需要两样东西。(i) 一个可以处理这种情况的DoFHandler类的版本，以及(ii) 告诉DoFHandler在哪个单元上使用哪个元素的方法。

这两件事中的第一件是在DoFHandler类的<i>hp</i>模式下实现的：不是将其与一个三角形和一个有限元对象相关联，而是与一个三角形和一个有限元集合相关联。第二部分是通过在DoFHandler的所有单元上的循环来实现的，并为每个单元设置集合中的有限元的索引，该索引将被用于该单元。我们将集合中用于某个单元的有限元对象的索引称为该单元的<i>active FE
index</i>，以表明这是在该单元上活动的有限元，而集合中的所有其他元素在该单元上是不活动的。其大致内容是这样的。

@code
  DoFHandler<dim> dof_handler(triangulation);
  for (auto &cell: dof_handler.active_cell_iterators())
    cell->set_active_fe_index(...);
  dof_handler.distribute_dofs(fe_collection);
@endcode 



调用 <code>set_active_fe_index()</code> 中的点表示我们以后要有某种策略来决定在哪个单元上使用哪个元素；我们以后会再来讨论这个问题。这里的重点是，这个代码片段的第一行和最后一行与非 $hp$ 的情况几乎完全相同。

另一个复杂的情况是，这次我们不只是有来自局部网格细化的悬挂节点，我们还必须处理这样的情况：如果有两个具有不同活动有限元指数的单元在一个面上相遇（例如一个Q2和一个Q3单元），那么我们必须计算有限元场的额外约束，以确保它是连续的。这在概念上与我们计算悬挂节点约束的方式非常相似，事实上，代码看起来完全一样。

@code
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
@endcode 

换句话说， DoFTools::make_hanging_node_constraints 不仅要处理悬挂节点约束，还要同时处理 $hp$  -约束。




<a name="Assemblingmatricesandvectorswithhpobjects"></a><h3>Assembling matrices and vectors with hp-objects</h3>


在这之后，我们必须为正确大小的线性系统设置矩阵和向量，并将它们组合起来。设置它们的方法与非 $hp$ 情况下的设置方法完全相同。组装则需要更多的思考。

主要的想法当然是不变的：我们必须在所有的单元中循环，集合局部贡献，然后将它们复制到全局对象中。正如第一次在 step-3 中详细讨论的那样，deal.II有FEValues类，它将有限元描述、映射和正交公式拉到一起，帮助评估形状函数的值和梯度，以及映射到单元实际位置的每个正交点的其他信息。每当我们转到一个新的单元时，我们就会重新初始化这个FEValues对象，从而要求它重新计算从一个单元到另一个单元的那部分信息。然后，它可以被用来总结对双线性形式和右手边的局部贡献。

在 $hp$ 有限元方法的背景下，我们必须处理这样一个事实，即我们不在每个单元上使用相同的有限元对象。事实上，我们甚至不应该对所有单元使用相同的正交对象，而应该对使用高阶有限元的单元使用高阶正交公式。同样地，我们可能也想在这些单元中使用高阶映射。

为了方便这些考虑，deal.II有一个类 hp::FEValues ，可以做我们在当前背景下需要的事情。不同的是，它不是一个单一的有限元、正交公式和映射，而是这些对象的集合。它的使用非常类似于常规的FEValues类，也就是说，在所有单元中循环的有趣部分看起来像这样。

@code
  hp::FEValues<dim> hp_fe_values(mapping_collection,
                                 fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      hp_fe_values.reinit(cell);


      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();


      ...  // assemble local contributions and copy them into global object
    }
@endcode 



在这个教程程序中，我们将始终使用Q1映射，所以 hp::FEValues 结构的映射集合参数将被省略。在循环内部，我们首先为当前单元格初始化 hp::FEValues 对象。第二个、第三个和第四个参数表示我们希望在该单元上使用的正交、映射和有限元对象在各自集合中的索引。这些参数可以省略（在下面的程序中也是如此），在这种情况下， <code>cell-@>active_fe_index()</code> 被用于该索引。之所以这样选择这些参数的顺序，是因为有时人们可能想从各自的集合中挑选不同的正交或映射对象，但几乎不可能从这个单元上使用的有限元中挑选不同的有限元，也就是说，索引与 <code>cell-@>active_fe_index()</code> 不同的有限元。因此，有限元集合索引是最后一个默认参数，这样可以方便地省略它。

这个 <code>reinit</code> 调用的作用如下： hp::FEValues 类检查它是否已经为这个有限元、正交和映射对象的组合分配了一个非 $hp$  -FEValues对象。如果没有，它将分配一个。然后，它为当前单元重新初始化这个对象，之后，现在有一个FEValues对象可用于当前单元上选定的有限元、正交和映射。然后通过调用 <code>hp_fe_values.get_present_fe_values()</code> 获得对该对象的引用，并将以通常的方式用于集合局部贡献。




<a name="Asimpleindicatorforhprefinementandestimatingsmoothness"></a><h3>A simple indicator for hp-refinement and estimating smoothness</h3>


自适应有限元方法的核心之一是，我们用一个指标来检查计算出的解（后验），告诉我们哪些是误差最大的单元，然后对它们进行细化。在其他许多教程程序中，我们使用KellyErrorEstimator类来获得单元格上的误差大小指示，尽管我们也在一些程序中讨论了更复杂的策略，最重要的是在  step-14  。

在任何情况下，只要决定只是 "细化这个单元 "或 "不细化这个单元"，实际的细化步骤就不是特别具有挑战性。然而，在这里，我们有一个能够进行hp细化的代码，也就是说，每当我们检测到某个单元上的误差太大，我们突然有两个选择：我们可以通过把它分割成几个小单元来细化这个单元，或者我们可以增加在它上面使用的形状函数的多项式程度。我们如何知道哪个是更有希望的策略？回答这个问题是本文写作时 $hp$ -无限元研究的核心问题。

简而言之，这个问题在目前的文献中似乎还没有解决。有许多或多或少复杂的方案来解决这个问题，但没有什么像KellyErrorEstimator那样被普遍接受为一个好的、即使不是最优的误差指标。大多数建议采用这样的事实：只要解是局部光滑的，增加多项式的度数是有益的，而只要网格是粗糙的，就应该细化。然而，如何确定解的局部光滑度以及决定一个解何时足够光滑以允许增加 $p$ 的问题无疑是大而重要的问题。

在下文中，我们提出了一个简单的解的局部光滑度的估计方法。正如我们将在结果部分看到的，这个估计器有缺陷，特别是就有局部悬空节点的单元而言。因此，我们不打算把下面的想法作为问题的完整解决方案。相反，它是作为一个值得进一步研究和调查的思路来处理的。换句话说，我们不打算在关于一般问题的答案的争论中进入一个复杂的建议。然而，为了证明我们对 $hp$ -无限元素的方法，我们需要一个简单的指标，它确实产生一些有用的信息，能够驱动本教程程序将进行的简单计算。


<a name="Theidea"></a><h4>The idea</h4>


我们的方法很简单：对于一个函数 $u({\bf x})$ 来说，它在单元格 $K$ 的Sobolev空间 $H^s(K)$ 中，必须满足条件@f[
   \int_K |\nabla^s u({\bf x})|^2 \; d{\bf x} < \infty.
@f] 假设单元格 $K$ 不是退化的，即。从单元格到单元格 $K$ 的映射足够规则，上述条件当然等同于@f[
   \int_{\hat K} |\nabla^s \hat u(\hat{\bf x})|^2 \; d\hat{\bf x} < \infty\,,
@f]，其中 $\hat u(\hat{\bf x})$ 是映射回单元格 $\hat K$ 的函数 $u({\bf x})$  。从这里，我们可以做以下工作：首先，让我们定义 $\hat u$ 的傅里叶级数为@f[
   \hat u(\hat{\bf x})
   = \sum_{\bf k} \hat U_{\bf k}\,e^{-i {\bf k}\cdot \hat{\bf x}},
@f]，傅里叶向量 ${\bf k}=(k_x,k_y)$ 在2d， ${\bf k}=(k_x,k_y,k_z)$ 在3d，等等，和 $k_x,k_y,k_z=0,2\pi,4\pi,\ldots$  。膨胀系数 $\hat U_{\bf k}$ 可以利用指数基@f[
\int_{\hat K} e^{-i {\bf m}\cdot \hat{\bf x}} e^{i {\bf n}\cdot \hat{\bf x}} d\hat{\bf x} = \delta_{\bf m \bf n},
@f]的 $L^2$ 正交性得到，这导致了以下表达式@f[
   \hat U_{\bf k}
   = \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat u(\hat{\bf x}) d\hat{\bf x} \,.
@f] 很明显，我们可以把 $H^s$ 的 $\hat u$  ] 的规范写成@f[
  \int_{\hat K} |\nabla^s \hat u(\hat{\bf x})|^2 \; d\hat{\bf x}
  =
  \int_{\hat K}
  \left|
    \sum_{\bf k} |{\bf k}|^s e^{-i{\bf k}\cdot \hat{\bf x}} \hat U_{\bf k}
  \right|^2 \; d\hat{\bf x}
  =
  \sum_{\bf k}
    |{\bf k}|^{2s}
    |\hat U_{\bf k}|^2.
@f] 换句话说，如果这个规范是有限的（即，为了使 $\hat u(\hat{\bf x})$ 在 $H^s(\hat K)$ 中），我们需要@f[
   |\hat U_{\bf k}| = {\cal O}\left(|{\bf k}|^{-\left(s+1/2+\frac{d-1}{2}+\epsilon\right)}\right).
@f] 换句话说：我们想要的规则性 $s$ 越高，傅里叶系数归零的速度越快。如果你想知道额外的指数 $\frac{d-1}2$ 从何而来：我们想利用这样一个事实： $\sum_l a_l < \infty$ 如果序列 $a_l =
{\cal O}(l^{-1-\epsilon})$ 为任何 $\epsilon>0$  。问题是，我们在这里不仅有一个单一变量的求和，而且有 $2\pi$ 的所有整数倍，这些整数倍位于 $d$ 维球内，因为我们有矢量成分 $k_x, k_y,
\ldots$  。就像我们通过用整个直线上的积分代替总和来证明上面的序列 $a_l$ 收敛一样，我们可以用 $d$ -维空间的积分代替我们的 $d$ -维总和。现在我们必须注意，在距离 $|{\bf k}|$ 和 $|{\bf k}|+d|{\bf k}|$ 之间，在一个常数以内，有 $|{\bf k}|^{d-1}$ 个模式，这与我们可以将体积元素 $dx\;dy$ 转化为 $2\pi r\; dr$ 的方式相同。因此，不再是 $|{\bf k}|^{2s}|\hat
U_{\bf k}|^2$ 必须衰减为 ${\cal O}(|{\bf k}|^{-1-\epsilon})$ ，而实际上是 $|{\bf k}|^{2s}|\hat U_{\bf k}|^2 |{\bf k}|^{d-1}$  。指数的比较产生了这个结果。

我们可以把这个问题转过来。假设我们得到了一个未知平滑度的函数 $\hat u$ 。让我们计算一下它的傅里叶系数 $\hat U_{\bf k}$ ，看看它们衰减的速度。如果它们以@f[
   |\hat U_{\bf k}| = {\cal O}(|{\bf k}|^{-\mu-\epsilon}),
@f]的速度衰减，那么我们这里的函数就在 $H^{\mu-d/2}$ 中。


<a name="Whatwehavetodo"></a><h4>What we have to do</h4>


那么，我们要做什么来估计 $u({\bf x})$ 在 $K$ 单元上的局部平稳性？显然，第一步是计算我们解决方案的傅里叶系数。傅里叶级数是无限级数，我们通过只计算级数的前几项来简化我们的任务，比如 $|{\bf k}|\le 2\pi N$ 有一个截止点 $N$  。让我们顺便说一下，我们希望选择 $N$ 足够大，这样我们至少可以捕获那些变化最大的形状函数的变化。另一方面，我们不应该把 $N$ 选得太大：显然，一个有限元函数，作为一个多项式，在任何给定的单元上都在 $C^\infty$ 中，所以系数将不得不在某一点上指数衰减；由于我们想估计这个多项式所近似的函数的平稳性，而不是多项式本身，我们需要为 $N$ 选择一个合理的截止点。无论怎样，计算这个数列并不特别难：从定义@f[
   \hat U_{\bf k}
   = \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat u(\hat{\bf x}) d\hat{\bf x}
@f]中我们可以看到，我们可以将系数 $\hat U_{\bf k}$ 计算为@f[
   \hat U_{\bf k}
   =
   \sum_{i=0}^{\textrm{dofs per cell}}
   \left[\int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat \varphi_i(\hat{\bf x})
   d\hat{\bf x} \right] u_i,
@f]，其中 $u_i$ 是这个单元上的 $i$ 个自由度的值。换句话说，我们可以把它写成@f[
   \hat U_{\bf k}
   = {\cal F}_{{\bf k},j} u_j,
@f]与矩阵@f[
   {\cal F}_{{\bf k},j}
   =
   \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat \varphi_j(\hat{\bf x}) d\hat{\bf x}.
@f]的矩阵-向量乘积。对于给定数量的形状函数 $\varphi_j$ 和傅里叶模式 $N$ ，这个矩阵很容易计算。因此，寻找系数 $\hat U_{\bf k}$ 是一个相当琐碎的工作。为了进一步简化我们的生活，我们将使用 FESeries::Fourier 类，它正是这样做的。

下一个任务是，我们必须估计这些系数的衰减速度  $|{\bf k}|$  。问题是，当然，我们首先只有有限的这些系数。换句话说，我们能做的最好的事情就是将一个函数 $\alpha |{\bf k}|^{-\mu}$ 拟合到我们的数据点 $\hat U_{\bf k}$ ，例如通过最小二乘法程序确定 $\alpha,\mu$ 。@f[
   \min_{\alpha,\mu}
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( |\hat U_{\bf k}| - \alpha |{\bf k}|^{-\mu}\right)^2
@f] 然而，这样做的问题是，它导致了一个非线性问题，这是我们想避免的事实。另一方面，如果我们试图将我们的系数的对数与 $\alpha |{\bf k}|^{-\mu}$ 的对数进行拟合，我们可以将这个问题转化为一个更简单的问题，就像这样。@f[
   \min_{\alpha,\mu}
   Q(\alpha,\mu) =
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( \ln |\hat U_{\bf k}| - \ln (\alpha |{\bf k}|^{-\mu})\right)^2.
@f] 利用通常的对数事实，我们可以看到，这将产生问题 @f[
   \min_{\beta,\mu}
   Q(\beta,\mu) =
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( \ln |\hat U_{\bf k}| - \beta + \mu \ln |{\bf k}|\right)^2,
@f] 其中  $\beta=\ln \alpha$  。现在这是一个问题，对于这个问题，优化条件  $\frac{\partial Q}{\partial\beta}=0,
\frac{\partial Q}{\partial\mu}=0$  ，在  $\beta,\mu$  中是线性的。我们可以把这些条件写成如下。@f[
   \left(\begin{array}{cc}
   \sum_{{\bf k}, |{\bf k}|\le N} 1 &
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|
   \\
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}| &
   \sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2
   \end{array}\right)
   \left(\begin{array}{c}
   \beta \\ -\mu
   \end{array}\right)
   =
   \left(\begin{array}{c}
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|
   \\
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}|
   \end{array}\right)
@f] 这个线性系统很容易倒置，得到@f[
   \beta =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }
@f] 和 @f[
   \mu =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }.
@f] 。

这无非是线性回归拟合，为了做到这一点，我们将使用  FESeries::linear_regression().  虽然我们对  $\beta$  的实际值不是特别感兴趣，但上面的公式给了我们一个计算指数  $\mu$  的平均值，然后我们可以用它来确定  $\hat u(\hat{\bf x})$  在  $H^s(\hat K)$  与  $s=\mu-\frac d2$  之间。

上面概述的这些步骤适用于许多不同的情况，这促使我们在deal.II中引入了一个通用函数 SmoothnessEstimator::Fourier::coefficient_decay() ，它将本节中描述的所有任务结合在一个简单的函数调用中。我们将在这个程序的实现中使用它。


<a name="Compensatingforanisotropy"></a><h4>Compensating for anisotropy</h4>


在上面的公式中，我们已经得出了傅里叶系数  $\hat U_{\bf
k}$  。因为 ${\bf k}$ 是一个矢量，对于相同的绝对值 $|{\bf k}|$ ，我们将得到许多傅里叶系数 $\hat U_{{\bf k}}$ ，对应于不同方向的傅里叶变换。如果我们现在考虑一个像 $|x|y^2$ 这样的函数，那么我们会发现在 $x$ 方向有很多大的傅里叶系数，因为这个方向的函数是不平滑的，但在 $y$ 方向有快速衰减的傅里叶系数，因为那里的函数是平滑的。由此产生的问题是：如果我们简单地将我们的多项式衰减 $\alpha |{\bf k}|^\mu$ 与<i>all</i>的傅里叶系数拟合，我们将把它拟合成一个光滑度<i>averaged in all spatial directions</i>。这就是我们想要的吗？还是只考虑所有 ${\bf k}$ 中幅度最大的系数 $\hat U_{{\bf k}}$ ，本质上是试图确定解在那个空间方向上的平滑度，在这个方向上解显得最粗糙？

人们也许可以为这两种情况争辩。如果deal.II有能力使用各向异性的有限元，即在不同的空间方向使用不同的多项式度数的有限元，那么这个问题将更有意义，因为它们能够更好地利用方向可变的光滑度。唉，在写这个教程程序的时候，这种能力并不存在。

无论如何，由于我们只有同位素的有限元类，我们采取的观点是，我们应该将多项式的程度调整到最低的规则性，以保持低的数值努力。因此，与其使用公式@f[
   \mu =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }.
@f] 计算如上所示的 $\mu$ ，我们必须稍微修改所有的和：与其在所有的傅里叶模式上求和，我们只求那些傅里叶系数是所有 $\hat U_{{\bf k}}$ 中幅度最大的 $|{\bf k}|$ ，即上述所有和必须由以下和代替。@f[
  \sum_{{\bf k}, |{\bf k}|\le N}
  \longrightarrow
  \sum_{\begin{matrix}{{\bf k}, |{\bf k}|\le N} \\ {|\hat U_{{\bf k}}| \ge |\hat U_{{\bf k}'}|
  \ \textrm{for all}\ {\bf k}'\ \textrm{with}\ |{\bf k}'|=|{\bf k}|}\end{matrix}}.
@f] 这就是我们将在程序中实现的形式。


<a name="Questionsaboutcellsizes"></a><h4>Questions about cell sizes</h4>


有人可能会问，我们只计算解的<i>reference cell</i>（而不是实数单元）上的傅里叶变换，这是否有问题？毕竟，在变换过程中，我们将解决方案拉伸了一个系数 $\frac 1h$ ，从而使傅里叶频率移动了一个系数 $h$ 。这是一个特别值得关注的问题，因为我们可能有相邻的单元，其网格大小 $h$ 相差2倍，如果其中一个单元比另一个更精细。正如我们在下面的结果部分所看到的那样，解决方案的估计平滑度应该是一个或多或少的连续函数，但在网格大小跳跃的地方却表现出跳跃性。因此，我们似乎很自然地要问，我们是否必须对变换进行补偿。

简短的回答是 "不"。在上述过程中，我们试图找到系数 $\beta,\mu$ ，使项@f[
   \ln |\hat U_{{\bf k}}| - \beta + \mu \ln |{\bf k}|.
@f]的平方之和最小化，以补偿变换意味着不试图适应衰减 $|{\bf k}|^\mu$  ]相对于傅里叶频率 ${\bf k}$ <i>on the unit
cell</i>，而是拟合在参考单元<i>to the Fourier frequencies on the real cell $|\bf
k|h$</i>上计算的系数 $\hat U_{{\bf k}}$ ，其中 $h$ 是变换算子的规范（即类似单元的直径）。换句话说，我们将不得不最小化@f[
   \ln |\hat U_{{\bf k}}| - \beta + \mu \ln (|{\bf k}|h).
@f]项的平方之和。然而，利用对数的基本属性，这只是相当于最小化 @f[
   \ln |\hat U_{{\bf k}}| - (\beta - \mu \ln h) + \mu \ln (|{\bf k}|).
@f] 换句话说，这和原始的最小二乘法问题将产生相同的最佳拟合指数  $\mu$  ，尽管在一种情况下偏移量是  $\beta$  而在另一种情况下是  $\beta-\mu \ln h$  。然而，由于我们对偏移量根本不感兴趣，而只对指数感兴趣，所以无论我们是否为了考虑网格大小的影响而对傅里叶频率进行缩放，都没有关系，在两种情况下估计的平滑度指数都是一样的。




<a name="Complicationswithlinearsystemsforhpdiscretizations"></a><h3>Complications with linear systems for hp-discretizations</h3>


<a name="Creatingthesparsitypattern"></a><h4>Creating the sparsity pattern</h4>


 $hp$ -方法的问题之一是，形状函数的高多项式程度与大量受限自由度一起导致矩阵在某些行有大量非零条目。同时，因为有些地方我们使用的是低多项式度，因此矩阵行的非零项相对较少。因此，为这些矩阵分配稀疏性模式是一个挑战：我们不能简单地从带宽的估计开始组装一个SparsityPattern，而不使用大量的额外内存。

我们为底层线性系统创建SparsityPattern的方式与我们用来执行约束的策略是紧密相连的。  <ol>   <li>  在不考虑约束条件的情况下组装矩阵，然后用 AffineConstraints::condense, 或 </li>   <li>  在我们用 AffineConstraints::distribute_local_to_global.</li>   </ol> 组装系统时应用约束条件，大多数基于deal.II的程序使用 DoFTools::make_sparsity_pattern  函数来分配一个考虑到约束的动态SparsityPattern。然后系统矩阵使用从DynamicSparsityPattern复制过来的SparsityPattern。这种方法在 step-2 中作了解释，并在大多数教程程序中使用。

早期的教程程序使用一阶或二阶有限元，因此删除疏散模式中对应于受限自由度的条目对矩阵明确存储的零的总体数量影响不大。然而，由于多达三分之一的自由度在hp分解中可能受到约束（而且，对于更高程度的元素，这些约束可以将一个自由度与多达10或20个其他自由度相联系），值得考虑这些约束，因为所产生的矩阵将更加稀疏（因此，矩阵-向量乘积或因子化也将大大加快）。


<a name="Eliminatingconstraineddegreesoffreedom"></a><h4>Eliminating constrained degrees of freedom</h4> 


 $hp$ 方法的第二个问题是，我们有很多受约束的自由度：通常有三分之一的自由度（在三维中）是受约束的，因为它们要么属于有悬空节点的单元，要么位于与具有更高或更低多项式度的单元相邻的单元上。事实上，这并不比非 $hp$ 模式中受约束自由度的比例多多少，但不同的是，每个受约束的悬空节点不仅受制于相邻的两个自由度，而且还受制于更多的自由度。

事实证明，在用 step-6 计算元素矩阵和向量时，首先提出的消除约束的策略，对于这种情况也是最有效的方法。另一种策略是首先建立没有约束的矩阵，然后 "浓缩 "掉有约束的自由度，这种策略要昂贵得多。事实证明，用这种低效率的算法建立稀疏模式至少需要 ${\cal O}(N \log N)$ 个未知数，而理想的有限元程序当然只有与未知数成线性的算法。对稀疏模式的创建以及矩阵的装配进行计时显示， step-6 中提出的算法（并在下面的代码中使用）确实更快。

在我们的程序中，我们还将把边界条件作为（可能是不均匀的）约束条件，并将矩阵的行和列也消除。为此我们所要做的就是在设置阶段调用插值Dirichlet边界条件的函数，以便告诉AffineConstraints对象关于它们的情况，然后同时在矩阵和向量上做从局部到全局的数据转移。这正是我们在  step-6  中所展示的内容。




<a name="Thetestcase"></a><h3>The test case</h3>


我们将用这个程序解决的测试案例是我们在 step-14 中已经看过的案例的重考：我们用 $f=(x+1)(y+1)$ 解决2d中的拉普拉斯方程@f[


   -\Delta u = f
@f]，并且用 $u$ 的零迪里希特边界值。我们在域 $[-1,1]^2\backslash[-\frac 12,\frac 12]^2$ 上这样做，即一个中间有一个正方形孔的正方形。

当然，与 step-14 不同的是，我们使用 $hp$ 的无限元来求解。这个测试案例是有意义的，因为它在洞的四角有重入角，在这些地方解有奇异性。因此，我们期望解在域的内部是平滑的，而在奇点附近是粗糙的。希望我们的细化和光滑度指标能够看到这种行为，并在远离奇点的地方细化网格，同时提高多项式的度数。正如我们将在结果部分看到的，情况确实如此。<a name="CommProg"></a> <h1> The commented program</h1>


<a name="Includefiles"></a><h3>Include files</h3>




前面几个文件已经在前面的例子中讲过了，因此不再进一步评论。

@code
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


@endcode 



这些是我们需要的新文件。第一个和第二个提供了FECollection和<i>hp</i>版本的FEValues类，如本程序介绍中所述。下一个文件提供了自动 $hp$ 适应的功能，为此我们将使用基于衰变系列扩展系数的估计算法，这是最后两个文件的一部分。

@code
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/numerics/smoothness_estimator.h>


@endcode 



最后一组include文件是标准的C++头文件。

@code
#include <fstream>
#include <iostream>



@endcode 



最后，这和以前的程序一样。

@code
namespace Step27
{
  using namespace dealii;



@endcode 




<a name="Themainclass"></a> <h3>The main class</h3> 




这个程序的主类看起来非常像前几个教程程序中已经使用过的主类，例如  step-6  中的那个。主要的区别是我们将refine_grid和output_results函数合并为一个，因为我们还想输出一些用于决定如何细化网格的量（特别是估计的解的光滑度）。   


就成员变量而言，我们使用与 step-6 中已经使用的相同的结构，但我们需要集合，而不是单个的有限元、正交和面状正交对象。我们将在类的构造函数中填充这些集合。最后一个变量， <code>max_degree</code>  ，表示使用的形状函数的最大多项式程度。

@code
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    ~LaplaceProblem();


    void run();


  private:
    void setup_system();
    void assemble_system();
    void solve();
    void create_coarse_grid();
    void postprocess(const unsigned int cycle);


    Triangulation<dim> triangulation;


    DoFHandler<dim>          dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;


    AffineConstraints<double> constraints;


    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;


    Vector<double> solution;
    Vector<double> system_rhs;


    const unsigned int max_degree;
  };





@endcode 




<a name="Equationdata"></a> <h3>Equation data</h3>   


接下来，让我们为这个问题定义右手边的函数。它在1d中是 $x+1$ ，在2d中是 $(x+1)(y+1)$ ，以此类推。

@code
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component) const override;
  };



  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    double product = 1;
    for (unsigned int d = 0; d < dim; ++d)
      product *= (p[d] + 1);
    return product;
  }





@endcode 




<a name="Implementationofthemainclass"></a> <h3>Implementation of the main class</h3>





<a name="LaplaceProblemLaplaceProblemconstructor"></a> <h4>LaplaceProblem::LaplaceProblem constructor</h4>




这个类的构造函数是相当直接的。它将DoFHandler对象与三角形关联起来，然后将最大多项式度数设置为7（在1d和2d中）或5（在3d及以上）。我们这样做是因为使用高阶多项式度数会变得过于昂贵，特别是在更高的空间维度上。   


在这之后，我们填补了有限元的集合，以及单元和面的四分法对象。我们从二次元开始，每个正交公式的选择是为了适合 hp::FECollection 对象中的匹配有限元。

@code
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    : dof_handler(triangulation)
    , max_degree(dim <= 2 ? 7 : 5)
  {
    for (unsigned int degree = 2; degree <= max_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }
  }



@endcode 




<a name="LaplaceProblemLaplaceProblemdestructor"></a> <h4>LaplaceProblem::~LaplaceProblem destructor</h4>




解构器与我们在  step-6  中已经做的没有变化。

@code
  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    dof_handler.clear();
  }



@endcode 




<a name="LaplaceProblemsetup_system"></a> <h4>LaplaceProblem::setup_system</h4>   


这个函数再次逐字复制了我们在  step-6  中已经做过的事情。尽管函数调用的名称和参数完全相同，但内部使用的算法在某些方面是不同的，因为这里的dof_handler变量是在  $hp$  -mode。

@code
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe_collection);


    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());


    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();


    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);


    system_matrix.reinit(sparsity_pattern);
  }





@endcode 




<a name="LaplaceProblemassemble_system"></a> <h4>LaplaceProblem::assemble_system</h4>




这是一个从每个单元的局部贡献组装全局矩阵和右手向量的函数。它的主要工作与之前的许多教程程序中描述的一样。重要的差异是<i>hp</i>有限元方法所需要的。特别是，我们需要使用FEValues对象的集合（通过 hp::FEValues 类实现），并且在将局部贡献复制到全局对象时，我们必须消除受限自由度。这两点在本程序的介绍中都有详细解释。   


还有一个轻微的复杂情况是，由于我们在不同的单元上使用了不同的多项式度数，所以持有局部贡献的矩阵和向量在所有单元上的大小不一样。因此，在所有单元的循环开始时，我们每次都必须将它们的大小调整到正确的大小（由 <code>dofs_per_cell</code> 给出）。因为这些类的实现方式是减少矩阵或向量的大小不会释放当前分配的内存（除非新的大小为零），所以在循环开始时调整大小的过程只需要在最初几次迭代中重新分配内存。一旦我们在一个单元中找到了最大的有限元度，就不会再发生重新分配，因为所有后续的 <code>reinit</code> 调用只会将大小设置为适合当前分配的内存。这一点很重要，因为分配内存是很昂贵的，而且每次我们访问一个新的单元时都这样做会花费大量的计算时间。

@code
  template <int dim>
  void LaplaceProblem<dim>::assemble_system()
  {
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);


    RightHandSide<dim> rhs_function;


    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;


    std::vector<types::global_dof_index> local_dof_indices;


    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();


        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;


        cell_rhs.reinit(dofs_per_cell);
        cell_rhs = 0;


        hp_fe_values.reinit(cell);


        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();


        std::vector<double> rhs_values(fe_values.n_quadrature_points);
        rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values);


        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
             ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_point) * // grad phi_j(x_q)
                   fe_values.JxW(q_point));           // dx


              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f(x_q)
                              fe_values.JxW(q_point));            // dx
            }


        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);


        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }





@endcode 




<a name="LaplaceProblemsolve"></a><h4>LaplaceProblem::solve</h4>




解决线性系统的函数与之前的例子完全没有变化。我们只是试图将初始残差（相当于右手边的 $l_2$ 准则）减少一定的系数。

@code
  template <int dim>
  void LaplaceProblem<dim>::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 1e-12 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);


    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);


    cg.solve(system_matrix, solution, system_rhs, preconditioner);


    constraints.distribute(solution);
  }





@endcode 




<a name="LaplaceProblempostprocess"></a> <h4>LaplaceProblem::postprocess</h4>




在解决了线性系统之后，我们要对解决方案进行后处理。在这里，我们所做的就是估计误差，估计解的局部平滑度，如介绍中所述，然后写出图形输出，最后根据之前计算的指标，在 $h$ 和 $p$ 中细化网格。我们在同一个函数中完成这一切，因为我们希望估计的误差和平滑度指标不仅用于细化，而且还包括在图形输出中。

@code
  template <int dim>
  void LaplaceProblem<dim>::postprocess(const unsigned int cycle)
  {
@endcode 



让我们从计算估计误差和平滑度指标开始，这两个指标对于我们三角测量的每个活动单元都是一个数字。对于误差指标，我们一如既往地使用KellyErrorEstimator类。

@code
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);


@endcode 



估算平滑度是用介绍中所述的衰减膨胀系数的方法进行的。我们首先需要创建一个对象，能够将每一个单元上的有限元解转化为一串傅里叶级数系数。SmoothnessEstimator命名空间为这样一个 FESeries::Fourier 对象提供了一个工厂函数，它为估计平滑度的过程进行了优化。然后在最后一个函数中实际确定每个单独单元上的傅里叶系数的衰减情况。

@code
    Vector<float> smoothness_indicators(triangulation.n_active_cells());
    FESeries::Fourier<dim> fourier =
      SmoothnessEstimator::Fourier::default_fe_series(fe_collection);
    SmoothnessEstimator::Fourier::coefficient_decay(fourier,
                                                    dof_handler,
                                                    solution,
                                                    smoothness_indicators);


@endcode 



接下来我们要生成图形输出。除了上面得出的两个估计量之外，我们还想输出网格上每个元素所使用的有限元的多项式程度。     


要做到这一点，我们需要在所有单元上循环，并使用  <code>cell-@>active_fe_index()</code>  轮询它们的活动有限元索引。然后我们使用这个操作的结果，在有限元集合中查询具有该索引的有限元，最后确定该元素的多项式程度。我们将结果放入一个矢量，每个单元有一个元素。DataOut类要求这是一个 <code>float</code> or <code>double</code> 的向量，尽管我们的值都是整数，所以我们就用这个。

@code
    {
      Vector<float> fe_degrees(triangulation.n_active_cells());
      for (const auto &cell : dof_handler.active_cell_iterators())
        fe_degrees(cell->active_cell_index()) =
          fe_collection[cell->active_fe_index()].degree;


@endcode 



现在有了所有的数据向量--解决方案、估计误差和平滑度指标以及有限元度--我们创建一个用于图形输出的DataOut对象并附加所有数据。

@code
      DataOut<dim> data_out;


      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.add_data_vector(estimated_error_per_cell, "error");
      data_out.add_data_vector(smoothness_indicators, "smoothness");
      data_out.add_data_vector(fe_degrees, "fe_degree");
      data_out.build_patches();


@endcode 



生成输出的最后一步是确定一个文件名，打开文件，并将数据写入其中（这里，我们使用VTK格式）。

@code
      const std::string filename =
        "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk";
      std::ofstream output(filename);
      data_out.write_vtk(output);
    }


@endcode 



在这之后，我们想实际细化网格，在  $h$  和  $p$  。我们要做的是：首先，我们使用估计的误差来标记那些误差最大的单元，以便进行细化。这就是我们一直以来的做法。

@code
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      0.3,
                                                      0.03);


@endcode 



接下来，我们要弄清楚哪些被标记为细化的单元格实际上应该增加 $p$ 而不是减少 $h$ 。我们选择的策略是，我们看一下那些被标记为细化的单元格的平滑度指标，对那些平滑度大于某个相对阈值的单元格增加 $p$ 。换句话说，对于每一个(i)细化标志被设置，(ii)平滑度指标大于阈值，以及(iii)我们在有限元集合中仍有一个多项式程度高于当前的有限元的单元，我们将分配一个未来的FE指数，对应于一个比当前高一个程度的多项式。下面的函数正是能够做到这一点。在没有更好的策略的情况下，我们将通过在标记为细化的单元上的最小和最大平滑度指标之间进行插值来设置阈值。由于角部奇点具有很强的局部性，我们将在数量上支持 $p$ -而不是 $h$ -细化。我们通过设置0.2的小插值因子，以低阈值实现这一点。以同样的方式，我们处理要被粗化的单元，当它们的平滑度指标低于在要被粗化的单元上确定的相应阈值时，减少它们的多项式程度。

@code
      hp::Refinement::p_adaptivity_from_relative_threshold(
        dof_handler, smoothness_indicators, 0.2, 0.2);


@endcode 



上述函数只决定了多项式程度是否会通过未来的FE指数发生变化，但并不操作 $h$  -细化标志。所以对于被标记为两种细化类别的单元格，我们更倾向于 $p$  -而不是 $h$  -细化。下面的函数调用确保只有 $p$  -或 $h$  -细化中的一个被施加，而不是同时施加。

@code
      hp::Refinement::choose_p_over_h(dof_handler);


@endcode 



对于网格自适应细化，我们通过调用 Triangulation::prepare_coarsening_and_refinement(). 将相邻单元的细化水平差异限制为1，以确保2:1的网格平衡。 我们希望对相邻单元的p水平实现类似的效果：未来有限元的水平差异不允许超过指定的差异。通过其默认参数，调用 hp::Refinement::limit_p_level_difference() 可以确保它们的级差被限制在1以内。这不一定会减少域中的悬挂节点的数量，但可以确保高阶多项式不会被限制在面的低得多的多项式上，例如，五阶多项式对二阶多项式。

@code
      triangulation.prepare_coarsening_and_refinement();
      hp::Refinement::limit_p_level_difference(dof_handler);


@endcode 



在这个过程的最后，我们再细化网格。在这个过程中，正在进行二分的单元的子单元会继承其母单元的有限元指数。此外，未来的有限元指数将变成活动的，因此新的有限元将在下一次调用 DoFHandler::distribute_dofs(). 后被分配给单元。 

@code
      triangulation.execute_coarsening_and_refinement();
    }
  }



@endcode 




<a name="LaplaceProblemcreate_coarse_grid"></a> <h4>LaplaceProblem::create_coarse_grid</h4>




在创建初始网格时，会使用以下函数。我们想要创建的网格实际上类似于 step-14 中的网格，即中间有方孔的方形域。它可以由完全相同的函数生成。然而，由于它的实现只是2d情况下的一种特殊化，我们将介绍一种不同的创建这个域的方法，它是独立于维度的。   


我们首先创建一个有足够单元的超立方体三角形，这样它就已经包含了我们想要的域 $[-1,1]^d$ ，并细分为 $4^d$ 单元。然后，我们通过测试每个单元上顶点的坐标值来移除域中心的那些单元。最后，我们像往常一样对如此创建的网格进行全局细化。

@code
  template <int dim>
  void LaplaceProblem<dim>::create_coarse_grid()
  {
    Triangulation<dim> cube;
    GridGenerator::subdivided_hyper_cube(cube, 4, -1., 1.);


    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
    for (const auto &cell : cube.active_cell_iterators())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        if (cell->vertex(v).square() < .1)
          cells_to_remove.insert(cell);


    GridGenerator::create_triangulation_with_removed_cells(cube,
                                                           cells_to_remove,
                                                           triangulation);


    triangulation.refine_global(3);
  }





@endcode 




<a name="LaplaceProblemrun"></a> <h4>LaplaceProblem::run</h4> 




这个函数实现了程序的逻辑，就像以前大多数程序中的相应函数一样，例如见  step-6  。   


基本上，它包含了自适应循环：在第一次迭代中创建一个粗大的网格，然后建立线性系统，组装它，求解，并对解决方案进行后处理，包括网格细化。然后再重新开始。同时，也输出一些信息给那些盯着屏幕试图搞清楚程序是干什么的人。

@code
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;


        if (cycle == 0)
          create_coarse_grid();


        setup_system();


        std::cout << "   Number of active cells      : "
                  << triangulation.n_active_cells() << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl
                  << "   Number of constraints       : "
                  << constraints.n_constraints() << std::endl;


        assemble_system();
        solve();
        postprocess(cycle);
      }
  }
} // namespace Step27



@endcode 




<a name="Themainfunction"></a> <h3>The main function</h3> 




主函数又是逐字逐句地执行我们之前的功能：把创建和运行一个主类的对象包装成一个 <code>try</code> 块，并捕捉任何抛出的异常，从而在出现问题时产生有意义的输出。

@code
int main()
{
  try
    {
      using namespace Step27;


      LaplaceProblem<2> laplace_problem;
      laplace_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;


      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }


  return 0;
}
@endcode 

<a name="Results"></a><h1>Results</h1>


在本节中，我们讨论了运行当前教程程序所产生的一些结果。更多的结果，特别是对三维计算的扩展和确定程序的各个部分需要多少计算时间，在 @ref hp_paper "hp-paper "中给出。

当运行时，这是该程序产生的结果。

@code
> make run
[ 66%] Built target step-27
[100%] Run step-27 with Release configuration
Cycle 0:
   Number of active cells      : 768
   Number of degrees of freedom: 3264
   Number of constraints       : 384
Cycle 1:
   Number of active cells      : 807
   Number of degrees of freedom: 4764
   Number of constraints       : 756
Cycle 2:
   Number of active cells      : 927
   Number of degrees of freedom: 8226
   Number of constraints       : 1856
Cycle 3:
   Number of active cells      : 978
   Number of degrees of freedom: 12146
   Number of constraints       : 2944
Cycle 4:
   Number of active cells      : 1104
   Number of degrees of freedom: 16892
   Number of constraints       : 3998
Cycle 5:
   Number of active cells      : 1149
   Number of degrees of freedom: 22078
   Number of constraints       : 5230
@endcode 



我们从中了解到的第一件事是，受限自由度的数量是总自由度的20-25%，至少在后来的网格上，当我们有相对高阶的元素时（在3D中，受限自由度的比例可以达到30%）。事实上，这与非 $hp$ 分化的数量级相同。例如，在 step-6 程序的最后一步，我们有18353个自由度，其中4432个是受约束的。不同的是，在后一个程序中，每个受约束的悬挂节点只对相邻的两个自由度进行约束，而在 $hp$ -案例中，受约束的节点对更多的自由度进行约束。还要注意的是，目前的程序在约束列表中还包括受迪里希特边界条件约束的节点。在第0周期，所有的约束条件实际上都是因为边界条件。

也许更有趣的是看一下图形输出。首先，这是该问题的解决方案。

<img src="https://www.dealii.org/images/steps/developer/  step-27  -solution.png" alt="解决方案的高程图，显示内部（重入）角附近缺乏规则性。" width="200" height="200">。

其次，让我们看一下生成的网格序列。

  <div class="threecolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-00.svg"
         alt="Triangulation containing reentrant corners without adaptive refinement."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-01.svg"
         alt="Triangulation containing reentrant corners with one level of
         refinement. New cells are placed near the corners."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-02.svg"
         alt="Triangulation containing reentrant corners with two levels of
         refinement. New cells are placed near the corners."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-03.svg"
         alt="Triangulation containing reentrant corners with three levels of
         refinement. New cells are placed near the corners."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-04.svg"
         alt="Triangulation containing reentrant corners with four levels of
         refinement. New cells are placed near the corners."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-05.svg"
         alt="Triangulation containing reentrant corners with five levels of
         refinement. New cells are placed near the corners."
         width="200" height="200">
  </div>
</div>   

可以清楚地看到网格在角部奇点附近是如何被细化的，正如人们所期望的那样。更有趣的是，我们应该好奇地看看这些网格单元的有限元多项式度数的分布，其中最浅的颜色对应的是2度，最深的对应的是7度。

  <div class="threecolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-00.svg"
         alt="Initial grid where all cells contain just biquadratic functions."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-01.svg"
         alt="Depiction of local approximation degrees after one refinement."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-02.svg"
         alt="Depiction of local approximation degrees after two refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-03.svg"
         alt="Depiction of local approximation degrees after three refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-04.svg"
         alt="Depiction of local approximation degrees after four refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-05.svg"
         alt="Depiction of local approximation degrees after five refinements."
         width="200" height="200">
  </div>
</div>   

虽然这肯定不是一个完美的安排，但它确实有一定的意义：我们在靠近边界和角落的地方使用低阶元素，那里的规则性很低。另一方面，在以下情况下使用高阶元素：(i) 误差一度相当大，即主要在角部奇点周围的一般区域和右上角解大的地方；(ii) 解是平滑的，即远离边界的地方。

这种多项式程度的安排当然是由我们的光滑性估计器得出的。这里是对解的光滑度的估计，深色表示最不光滑，浅色表示最光滑的区域。

  <div class="threecolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-00.svg"
         alt="Estimated regularity per cell on the initial grid."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-01.svg"
         alt="Depiction of the estimated regularity per cell after one refinement."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-02.svg"
         alt="Depiction of the estimated regularity per cell after two refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-03.svg"
         alt="Depiction of the estimated regularity per cell after three refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-04.svg"
         alt="Depiction of the estimated regularity per cell after four refinements."
         width="200" height="200">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-05.svg"
         alt="Depiction of the estimated regularity per cell after five refinements."
         width="200" height="200">
  </div>
</div>   

从这里可以得出的主要结论是，内角的规则性损失是一个高度局部的现象；它似乎只影响到与角本身相邻的单元，所以当我们细化网格时，黑色的颜色就不再明显。除了角部，这一连串的图示意味着平滑度估计在某种程度上与网格细化无关，特别是当我们远离边界时。同样明显的是，平滑度估计值与解的实际尺寸无关（见上面的解的图片），这也是应该的。然而，更值得关注的一点是，人们在仔细观察后发现，我们的估计器高估了有悬空节点的单元格的解决方案的平滑度。这反过来又导致了这些区域的多项式度数更高，使有限元在单元上的分配出现偏差。

目前我们对这种影响没有很好的解释。一种理论是，有悬空节点的单元上的数值解当然是受限制的，因此不能完全自由地探索函数空间以接近精确解。这种自由度的缺乏可能表现为在这些单元上产生具有抑制振荡的数值解，意味着更高的平滑度。估计器会捕捉到这个信号，估计的平滑度会高估实际值。然而，这个程序的作者目前还没有找到发生了什么的明确答案。

当然，更大的问题是如何避免这个问题。可能的方法包括不是在单个细胞上估计平滑度，而是在每个细胞周围的细胞集合体或斑块上估计。也可以根据每个细胞所具有的受限自由度的数量，为其找到简单的修正系数。无论哪种情况，都有大量的机会来进一步研究寻找好的 $hp$ -精化标准。另一方面，目前方案的主要内容是在deal.II中演示使用 $hp$ -技术，这不受我们使用可能的次优细化标准的影响。




<a name="extensions"></a> <a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Differenthpdecisionstrategies"></a><h4>Different hp-decision strategies</h4>


本教程只演示了决定 $h$ -和 $p$ -适应的一种特殊策略。事实上，还有许多自动决定适应类型的方法，其中一些已经在deal.II中实现： <ul>   <li>  <i>Fourier coefficient decay:</i> 这是本教程中目前实现的策略。关于这个策略的更多信息，请参见 SmoothnessEstimator::Fourier 命名空间的一般文档。 </li>   

    <li>  <i>Legendre coefficient decay:</i> 这个策略与目前的策略很相似，但使用了Legendre级数展开，而不是傅里叶级数：这个策略使用Legendre多项式而不是正弦波作为基础函数。当然，由于我们在每个单元上使用有限维度的多项式来近似解，因此解在Legendre多项式中的扩展也是有限的，因此，当我们谈论这个扩展的 "衰减 "时，我们只能考虑这个扩展的有限多个非零系数，而不是用渐进的方式来思考。  但是，如果我们有足够多的这些系数，我们当然可以把这些系数的衰减看作是精确解的系数衰减的特征（一般来说，精确解不是多项式的，所以会有一个无限的Legendre扩展），考虑我们拥有的系数应该可以揭示一些关于精确解的特性。

  从傅里叶策略到Legendre策略的过渡非常简单。  你只需要改变数列扩展类和相应的平滑度估计函数，使其成为适当命名空间 FESeries::Legendre 和 SmoothnessEstimator::Legendre. 的一部分。 关于这一策略的理论背景，请查阅 SmoothnessEstimator::Legendre 命名空间的一般文档，以及 @cite mavriplis1994hp 、 @cite eibner2007hp 和 @cite davydov2017hp 。 </li>   

    <li>  <i>Refinement history:</i> 最后一个策略与其他两个有很大不同。在理论上，我们知道在改变函数空间离散化后，误差会如何收敛。通过 $h$ -细化，解决方案以代数方式收敛，正如在 step-7 中已经指出的。不过，如果解足够平滑，我们预计解将随着有限元的多项式程度的增加而指数级收敛。我们可以在下面的步骤中把对误差的适当预测与实际误差进行比较，以了解我们对适应类型的选择是否合理。

  过渡到这种策略的过程比较复杂。为此，我们需要一个纯 $h$ -或 $p$ -细化的初始化步骤，我们需要将预测的误差转移到适应的网格上。 hp::Refinement::predict_error() 函数的大量文档不仅描述了这种方法的理论细节，而且还介绍了如何在你的代码中实现这种策略的蓝图。欲了解更多信息，请参见  @cite melenk2001hp  。

  请注意，在时间依赖性问题上，用这个特殊的函数你无法预测下一个时间步骤的误差。因此，这种策略不能不加思索地应用于这种类型的问题。另外，也可以使用下面的方法，它也适用于所有其他的策略：从每个时间步长的粗网格开始，不断细化，直到对结果满意，然后才进入下一个时间步长。 </li>  </ul>  。 

试着在本教程中实施这些策略之一，观察结果的细微变化。你会注意到，所有的策略都能够识别重心角附近的奇点，并将在这些区域执行 $h$ -精简，而在体域中更倾向于 $p$ -精简。这些策略的详细比较见 @cite fehling2020  。


<a name="Parallelhpadaptivefiniteelements"></a><h4>Parallel hp-adaptive finite elements</h4>


本教程中介绍的所有功能都已经适用于顺序和并行应用。不费吹灰之力就可以改成  parallel::shared::Triangulation  或  parallel::distributed::Triangulation  类。如果你觉得急于尝试，我们建议先阅读前者的 step-18 和后者的 step-40 ，以了解有关该主题的进一步背景信息，然后再回来看本教程，试试你新获得的技能。<a name="PlainProg"></a> <h1> The plain program</h1>  @include "step-27.cc" 。 

  */  
