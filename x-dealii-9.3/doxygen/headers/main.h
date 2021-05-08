// ---------------------------------------------------------------------
//ØØ
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//∏∏
// ---------------------------------------------------------------------


/**
 *@mainpage
 *Deal.ii pre 9.3 中文翻译版，如有任何建议，请联系 jiaqiwang969@gmail.com
 *
 *这是deal.II类和函数的主要起始页文件。其他方面的文档，例如构建系统，可以在别处找到。此外，还有关于使用库的<a href="Tutorial.html">教程</a>。
 *
 *deal.II库中基本涵盖几大模块（请参阅 <a href="modules.html">模块页面</a>或中的相应条目此页顶部的菜单）。这些模块最终拼装为完整的有限元代码。deal.II interact中的类数量由以下图表给出。可以了解deal.II的主类的组织架构，还可以单击进入类目（灰色框表示可选外部库，灰色椭圆可选外部库的子集deal.II可以与之交互的应用程序）：
 * 
 
 * @dot
 digraph G
{
  graph[rankdir="TB",bgcolor="transparent"];

  node [fontname="FreeSans",fontsize=15,
        shape=record,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];
  edge [color="black", weight=10];

  tria       [label="Triangulation",    URL="\ref grid"];
  fe         [label="Finite elements",    URL="\ref feall"];
  mapping    [label="Mapping",          URL="\ref mapping"];
  quadrature [label="Quadrature",       URL="\ref Quadrature"];
  dh         [label="DoFHandler",       URL="\ref dofs"];
  fevalues   [label="FEValues",         URL="\ref feaccess"];
  systems    [label="Linear systems",   URL="\ref LAC"];
  solvers    [label="Linear solvers",   URL="\ref Solvers"];
  output     [label="Graphical output", URL="\ref output"];
  manifold   [label="Manifold",         URL="\ref manifold"];

  tria -> dh              [color="black",style="solid"];
  fe -> dh                [color="black",style="solid"];
  fe -> fevalues          [color="black",style="solid"];
  mapping -> fevalues     [color="black",style="solid"];
  quadrature -> fevalues  [color="black",style="solid"];
  dh -> systems           [color="black",style="solid"];
  fevalues -> systems     [color="black",style="solid"];
  systems -> solvers      [color="black",style="solid"];
  solvers -> output       [color="black",style="solid"];
  manifold -> tria        [color="black",style="solid"];
  manifold -> mapping     [color="black",style="solid"];



  node [fontname="FreeSans",fontsize=12,
        shape=record,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];

  opencascade [label="OpenCASCADE"];
  subgraph linalglibs {
    rank="same";
    petsc       [label="PETSc",    URL="\ref PETScWrappers"];
    trilinos    [label="Trilinos", URL="\ref TrilinosWrappers"];
    cuda        [label="CUDA",     URL="\ref CUDAWrappers"];
  }
  umfpack     [label="UMFPACK"];

  petsc -> systems        [dir="none"];
  petsc -> solvers        [dir="none"];
  trilinos -> systems     [dir="none"];
  trilinos -> solvers     [dir="none"];
  cuda -> systems         [dir="none"];
  cuda -> solvers         [dir="none"];
  umfpack -> solvers      [dir="none"];
  opencascade -> manifold [dir="none"];


  node [fontname="FreeSans",fontsize=12,
        shape=ellipse,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];

  gmsh        [label="gmsh", URL="\ref Gmsh"];
  visit       [label="VisIt"]
  paraview    [label="ParaView"]

  gmsh -> tria       [dir="none"];
  output -> visit    [dir="none"];
  output -> paraview [dir="none"];
}
 * @enddot
 *
 *这些模块都包括在程序教程中，可以参考setp-3，其中给出了它们如何组合在一起的概述。以下是一个指南，包括分类的基本介绍，以及链接与他们每个人有关的文件:
 *
 *<ol>
 *<li> <b>%Triangulation</b>: 三角剖分是单元格cells及其低维边界对象的集合。单元是参考超立方体[0,1]<sup>dim</sup> 的映像，在参考单元和实单元之间的映射模块中有适当
 *的 @ref mapping 映射。
 * 
 *三角剖分存储网格的几何和拓扑特性：单元如何连接以及它们的顶点在哪里。三角剖分不知道您可能要在该网格上使用的有限元的任何信息，三角剖分甚至不知道其单元的形状：在2d中，它只知道单元有4个面（线）和4个顶点（在3d中，它有6个面（四边形）、12条线和8个顶点），但其他一切都是由映射类定义的。
 * 
 *三角剖分的属性和数据几乎总是通过所有单元格上的循环进行查询，也可能查询每个单元格的所有面。因此，有关网格的大部分知识都隐藏在迭代器 @em iterators 后面，即指针式结构，可以从一个单元格迭代到下一个单元格，并且可以询问有关它当前指向的单元格的信息。
 * 
 *描述三角剖分和单元的类位于"栅格和三角剖分"模块中并记录在 @ref grid。迭代器在"类网格容器上的迭代器"模块中 @ref Iterators 进行了描述。
 * 
 *<li> <b>%Manifold</b>: 流形描述细胞的形状，更一般地说，描述求解方程的域的几何结构。他们使用微分几何的语言。更多信息可以在三角剖分的流形描述中找到 @ref manifold 。
 * 
 *<li> <b>Finite Element</b>:有限元类描述单元上定义的有限元空间的属性。例如，这包括顶点、直线或单元内部的自由度。除此之外，有限元类当然必须在单元上的点上提供单个形状函数的值和梯度。
 * 
 *有限元模块中描述了有限元类 @ref feall。
 * 
 *<li> <b>%Quadrature</b>: 与有限元一样，求积对象定义在单元上。它们只描述了正交点在单元上的位置，以及正交点在单元上的权重。
 * 
 *描述特定求积公式的类的文档可以在求积公式模块中找到 @ref Quadrature。
 * 
 *<li> <b>%DoFHandler</b>: DoFHandler 对象是三角剖分和有限元的汇合点：有限元类描述每个顶点、线或单元需要多少自由度，DoFHandler 类分配这个空间，以便三角剖分的每个顶点、线或单元都有正确的自由度数。它还提供了一个全局编号。
 * 
 *另一种观点是：虽然网格和有限元描述了有限维空间 $V_h$ 的抽象性质，我们在其中寻找离散解，但 DoFHandler 类列举了这个空间的具体基础，以便我们可以用一组有序的系数 $U_j$ 来表示离散解 $u_h(\mathbf x)= \sum_j U_j \varphi_i(\mathbf x)$。
 * 
 *就像三角剖分对象一样，DoFHandler 上的大多数操作都是通过在所有单元格上循环并对每个单元格或其中的一个子集执行操作来完成的。因此，这两个类的接口非常相似：它们允许获取第一个和最后一个单元格（或面或行等）的迭代器，并通过这些迭代器提供信息。可以从这些迭代器获得的信息是已经可以从三角剖分迭代器（它们实际上是派生类）获得的几何和拓扑信息，以及诸如当前单元上的自由度的全局数之类的信息。还可以要求迭代器从存储与三角剖分相关联的所有自由度值的数据向量中提取与当前单元格上的自由度对应的值。
 * 
 *值得注意的是，就像三角剖分一样，DoFHandler 类不知道从单位单元到单个单元的映射。它也不知道与它管理的自由度相对应的形状函数：它只知道，例如，每个顶点有2个自由度，每个单元内部有4个自由度。它们的细节与DoFHandler类无关，但它们存在的事实除外。
 * 
 *DoFHandler 类及其关联在自由度模块 @ref dofs 中进行了描述。此外，还有专门的版本可以处理多级和hp离散化 。
 *这些在"多级支持" @ref mg 和 "hp有限元支持" @ref hp 模块中进行了描述。有限元方法通常意味着对自由度的约束，例如悬挂节点或边界条件适用的节点；自由度约束模块 @ref constraints 中描述了处理此类约束的方法。
 * 
 *<li> <b>%Mapping</b>: 有限元程序的下一步是使用有限元的形状函数和由求积规则定义的求积点来计算矩阵和右手边项或三角形剖分的每个单元上的其他量。为此，有必要将形状函数、正交点和正交权重从单位单元映射到三角剖分的每个单元。这不是由映射和派生类直接完成的，而是由映射 Mapping 和派生类促进的：它们描述了如何将点从单位映射到实空间并返回，还提供了这种导数和雅可比行列式的梯度。 
 * 
 *这些类都在 @ref mapping 映射中描述。
 * 
 *<li> <b>%FEValues</b>: 下一步是实际采取一个有限元，并评估其形状函数和梯度点定义的求积公式时，映射到真正的细胞。这是 FEValues 类和兄弟类的领域：在某种意义上，它们提供了有限元函数空间的逐点视图。
 * 
 *这似乎是有限制的：在数学分析中，我们总是用单元上的积分或单元的面来写公式，涉及有限元形状函数。因此，有人认为有必要将有限元空间描述为连续空间。然而，在实践中，这并不是必须的：所有的积分在实际计算中都被使用求积公式的近似值所代替，因此真正需要的只是在一个域内有限个给定位置处计算形状函数的能力。FEValues 类提供的正是这些信息：给定有限元、求积和映射对象，它们计算连续函数空间（与离散相对，而不是与不连续相对）对离散点数的限制。有许多对象可以执行此操作：FEValues 用于对单元格求值， FEFaceValues 用于对单元格面求值，FESubfaceValues 用于对单元格面部分求值。所有这些类都在 @ref feaccess 模块中描述。
 * 
 *<li> <b>Linear Systems</b>: 线性系统：如果一个人知道如何用 FEValues 和friends计算单个单元上的形状函数的值和梯度，并且知道如何用 DoFHandler 迭代器得到单元上自由度的全局数，然后，下一步是使用问题的双线性形式来组装线性系统的系统矩阵（和右侧）。然后我们将从这个线性系统中确定问题的解。
 * 
 *为此，我们需要有存储和管理矩阵和向量项的类。deal.II为此提供了一整套类，以及与提供类似功能的其他软件包的接口。这方面的文档可以在线性代数类模块 @ref LAC 中找到。
 * 
 * 
 *<li> <b>Linear Solvers</b>:  为了确定有限维线性方程组的解，需要线性解算器。在有限元应用中，它们经常是迭代的，但有时也可能需要使用直接或稀疏的直接求解器。deal.II有很多这样的东西。它们记录在"线性解算器类" @ref Solvers 模块中。
 * 
 *<li> <b>Output</b>:  最后，若已经获得了一个有限元问题的解决方案，对给定的三角剖分，人们往往会希望它的后处理使用可视化程序。这个库本身并不进行可视化，而是生成各种图形格式的输出文件，这些图形格式被广泛使用的可视化工具所理解。
 * 
 *图形输出模块 @ref output 中给出了这样做的类的描述。
 *</ol>
 * *此外，deal.II还有许多类组，它们超出了这里列出的类组。它们涉及到上述层次结构的更精细的概念，或者涉及到诸如输入和输出的处理等不一定特定于有限元程序的相关方面，但也出现在那里。这些类都列在本页顶部菜单栏上的"类和名称空间"视图中，并且也被分组到各自的模块中（请参阅本页顶部的"模块"链接 <a href="modules.html">Modules link</a> ）。
 * 
 *我们为那些希望将应用程序文档直接链接到deal.II在线文档的用户提供Doxygen标记文件。标记文件位于 <a href="../deal.tag"><code>deal.tag</code></a> 。对于deal.II的每个版本，它都位于Doxygen参考文档正上方的目录中。为了使用标记文件，您必须将其下载到Doxygen可以找到它的地方。之后，在Doxygen选项文件中找到关键标记文件 <code>TAGFILES</code> ，并编写如下代码
 *<pre>
 *TAGFILES = deal.tag=http://www.dealii.org/X.Y.Z/doxygen/deal.II
 *</pre>
 *其中 <code>X.Y.Z</code> 指向要链接到的版本。确保使用匹配的标记文件。理论上，您也可以针对deal.II的修订版本进行链接，但是您必须担心，如果deal.II结构发生变化，您的链接可能会失效。
 */
