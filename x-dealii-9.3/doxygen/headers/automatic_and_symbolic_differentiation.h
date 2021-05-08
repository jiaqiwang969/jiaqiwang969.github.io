// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/**
 * @defgroup auto_symb_diff Automatic and symbolic differentiation
 *
 * @brief 专用于实现与自动和符号微分有关的函数和类的模块。
 *
 *下面我们将简要介绍什么是自动和符号微分，这些计算/数值格式存在哪些变化，以及它们如何集成到deal.II的框架中。所有这些方案的目的都是在不想手工计算的情况下，自动计算函数的导数或近似值。常见的例子是在有限元环境中，想要解决一个非线性问题，这个问题一般转化为求残差$F(u,\nabla u)=0$ 的形式，需要通过牛顿法求复杂函数的微分进行优化迭代求解；或者给定一些参数下的优化问题，${\cal A}(q,u,\nabla u)=f$，并希望获得对参数 $q$ 的导数，例如，优化关于q的输出泛函，或用于关于q的灵敏度分析。我们应把 $q$ 看作是设计参数：比如说，机翼的宽度或形状，用来建造物体的材料的刚度系数，传送到设备的功率，传送到燃烧器的气体的化学成分。在所有这些情况下，我们应该想到 $F$ 和 $\cal A$ 区分起来既复杂又麻烦——至少在手工计算时是这样。一个相对简单的非线性问题已经凸显了手工计算导数的繁琐，如 step-15 所示。然而，在现实中，人们可能会思考一些问题，如化学反应流，其中流体方程具有密度和粘度等系数，这些系数强烈而非线性地依赖于流体在每个点的化学成分、温度和压力；化学物质之间的反应是基于反应系数的，反应系数也是非线性的，并且以复杂的方式依赖于化学成分、温度和压力。在许多情况下，所有这些系数的精确公式可以用几行来写出，可以包括几个非线性项的指数和（调和或几何）平均值，亦或可以包含数据点之间的查表和插值。仅仅把这些条件弄对就够难了；计算这些项的导数在大多数应用中是不切实际的，实际上，不可能得到正确的结果。如果没有计算机的帮助，更高的导数就更不可能实现。自动或符号微分是一种解决方法：只需实现一个函数，该函数只需根据输入计算这些系数一次，就可以得到（正确的！）无需进一步编码工作的导数（尽管在运行时、编译时或两者都有不可忽略的计算成本）。
 *
 *
 * @section auto_diff_1 Automatic differentiation
 *
 *自动微分 <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">Automatic differentiation </a> (通常也称为算法微分）是一种数值方法，可用于“自动”计算函数相对于一个或多个输入变量的第一阶（可能是更高阶）导数。尽管这需要一定的计算成本，但使用这种工具的好处可能是巨大的。当正确使用时，通常复杂函数的导数可以计算到非常高的精度。尽管这些框架所能达到的精确精度在很大程度上取决于其底层的数学公式，但有些实现的计算精度取决于机器精度的数量级。注意，这不同于经典的数值微分法（例如，通过在不同的点对函数求值来使用有限差分近似），后者的精度取决于微扰大小以及所选的有限差分格式；这些方法的误差比公式化的自动微分方法要大得多。
 *
 *接下来，我们将举在有限元环境中使用自动微分的三个实际例子：
 * - 一个新的非线性公式的快速原型，无需手动计算线性化本身,
 * - 在复杂的多物理框架内，有限元残差的自动线性化，
 * - 验证基于单元的计算（如残差）和基于连续点的计算（如非线性本构关系的切线）的线性化的用户实现.
 *
 *自动可微数的实现有很多种。它们主要分为两大类，即源代码转换 <em>source code transformation</em> 和运算符重载 <em>operator overloading</em> 。第一个方法基于某个输入函数生成新的可编译代码，该函数在执行时返回输入函数的导数。第二个是利用C++操作符定义的能力来为自定义类类型重载。因此，表示这样一个自可微的类，在对其执行每一个数学运算之后，原则上可以计算并跟踪它的值及其方向导数的值。由于专门实现源代码转换方法的库共同描述了用作函数预处理器的高度专业化的工具，因此它们在deal.II本身中没有直接的支持。但是，后者表示通过在适当的上下文中使用模板元编程可以支持的专用数字类型。给出上面的例子，这意味着 FEValues 类（和friends）以及张量 Tensor 和 SymmetricTensor 类应该支持使用这些特殊数字进行的计算(理论上，一个完整的程序是可以微分的。例如，这在解决方案对输入参数的敏感性分析中是有用的。但是，到目前为止，这还没有经过测试。）
 *
 *基于操作符重载 <em>operator overloading</em> 的专用框架的实现通常分为三类。在每一个函数中，都有一些自定义的数据类，表示被求值函数的浮点值及其导数
 *
 * -# 利用双重 <em>dual</em> /复杂阶跃 <em>complex-step</em> / 超双重公式（有时称为无带方法 <em>tapeless</em> ），
 * -# 利用 有带方法 <em>taping</em> , 
 * -# 通过表达式模板 <em>expression templates</em> 使用编译时优化的策略。
 *
 *为了初步了解这些不同的实现在实践中可能是什么样子，我们提供以下这些方法的通用摘要：
 * -# 上面列出的前两种无带方法<em>tapeless</em> (双数法和复步长法）使用截断泰勒级数的一些变化，以及摄动参数定义的特定选择，使用基于有限差分的方法计算函数导数。“对偶”数构成在计算函数值时同时计算的累积方向导数；在复步长方法中，虚值有效地服务于此目的。摄动参数的选择决定了格式的数值性质，如泰勒格式截断的影响；对偶数的一阶导数中不包含任何高阶项，而复步长法忽略了已有的高阶项。结果表明，这两种方法都不受相减对消误差的影响，而且在它们的有限差分格式中，它们对为数值扰动选择的内部步长不敏感。因此，对偶数方法产生精确的一阶导数，而复阶近似则不产生。然而，对偶数的标准实现不能产生二阶导数的精确值。超对偶数对这一观点有不同的看法，数字以类似于四元数的形式表示（即携带额外的非实数分量），导数是从泰勒级数的所有四个分量的高阶截断计算出来的。结果是，通过适当的实现，一阶导数和二阶导数都可以精确计算。
 * 
 * -# 在有带方法 <em>taping</em> 中，选择指定的代码子区域作为一个子区域，对该子区域使用活动（标记）输入变量执行的所有操作进行跟踪，并将其记录在称为磁带的数据结构中。在磁带区域的末尾，可以通过使用不同的输入变量集“重放”磁带来重新评估记录的函数，而不是直接重新计算函数。假设磁带区域代表一个光滑函数，则可以通过参考磁带上跟踪和存储的代码路径来计算函数的任意高阶导数(例如，可以通过对感兴趣点附近的函数进行求值来实现这一点。）存在一些策略，用于处理在求值点处录制的函数不平滑或不具有分析性的情况。此外，可能需要考虑分支函数的情况，其中磁带不再是连续的，而是由于原始记录的输入而在不同的求值路径上分叉。
* 
* -# 基于表达式模板 <a href="https://en.wikipedia.org/wiki/Expression_templates"> 的方法利用由抽象语法树（AST）构造的计算图（在本例中是一个有向无环图（DAG）<a href="https://en.wikipedia.org/wiki/Directed_acyclic_graph">directed acyclic graph (DAG)</a>)），它从函数的输入值解析函数输出。树上最外层的叶子代表自变量或常量，由一元运算符转换并由二元运算符连接（在最简单的情况下）。因此，在编译时对函数输入执行的操作是已知的，并且与之相关联的导数操作也可以同时使用计算操作的导数的众所周知的规则（例如加法和减法下的导数的关联性、乘积规则、，以及链式法则）。此运算符返回的编译输出类型不需要是泛型的，而是可以基于给定给DAG顶点上的特定运算符的特定输入（可能包含微分历史）进行专门化。通过这种方式，可以为用于评估依赖函数的每个中间结果的非常专门的单个操作生成编译时优化的指令集。
*
*
*当然，每一种方法都有其优缺点，对于要解决的给定问题，一种方法可能比另一种方法更合适。由于上述实现细节（以及其他未讨论的细节）可能对用户隐藏，因此理解使用这些“黑盒”自动可微数中的任何一个的含义、运行时成本和潜在限制可能仍然很重要。
 *
 *除了提供的链接文章外，用于提供此处提供的详细信息的资源还包括：
 * @code{.bib}
 * @InProceedings{Fike2011a,
 *   author    = {Fike, Jeffrey A and Alonso, Juan J},
 *   title     = {The Development of Hyper-Dual Numbers for Exact Second-Derivative Calculations},
 *   booktitle = {49th {AIAA} Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition},
 *   year      = {2011},
 *   volume    = {886},
 *   pages     = {124},
 *   month     = {jan},
 *   publisher = {American Institute of Aeronautics and Astronautics},
 *   doi       = {10.2514/6.2011-886},
 * }
 * @endcode
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 *   title     = {Getting Started with ADOL-C},
 *   author    = {Walther, Andrea and Griewank, Andreas},
 *   year      = {2009},
 *   booktitle = {Combinatorial scientific computing},
 *   doi       = {10.1.1.210.4834},
 *   pages     = {181--202}
 * }
 * @endcode
 *
 * ### 链式规则的利用 Exploitation of the chain-rule
 *
 *在最实际的意义上，上述任何类别都利用链式规则来计算复合函数的总导数。为了执行这个操作，他们通常使用两种机制中的一种来计算导数，特别是 
 * - 正向模式（或正向累积）自动微分 <em>forward-mode</em> (or <em>forward accumulation</em>) auto-differentiation, 
 * - 反向模式（或反向累积）自动微分 <em>reverse-mode</em> (or <em>reverse accumulation</em>) auto-differentiation.
 *
 *作为一个关注点，最优雅可比积累 <em>optimal Jacobian accumulation</em>，它执行一个最小的计算集，位于这两个极限情况之间。一般复合函数的计算在图论中仍是一个开放的问题。
 *
 *借助下面的图表（它和维基百科文章中列出的一些细节 <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">Wikipedia article</a> ），让我们考虑一下函数 $f (\mathbf{x}) = \sin (x_{1}) + x_{1} x_{2}$ 以及它带导数的计算方法:
 *
 * <div class="twocolumn" style="width: 80%">
 *   <div class="parent">
 *     <div class="img" align="center">
 *       <img src="https://upload.wikimedia.org/wikipedia/commons/a/a4/ForwardAccumulationAutomaticDifferentiation.png"
 *            alt="Forward mode automatic differentiation"
 *            width="400">
 *     </div>
 *     <div class="text" align="center">
 *       Forward mode automatic differentiation
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       <img src="https://upload.wikimedia.org/wikipedia/commons/a/a0/ReverseaccumulationAD.png"
 *            alt="Reverse mode automatic differentiation"
 *            width="400">
 *     </div>
 *     <div class="text" align="center">
 *       Reverse mode automatic differentiation
 *     </div>
 *   </div>
 * </div>
 *
 *具体来说，我们将简要描述什么是正向和反向自动微分。注意，在图中，沿文本中图形的边缘是函数 $w$ 相对于第 $i$ 个变量的方向导数，用符号 $\dot{w} = \dfrac{d w}{d x_{i}}$ 表示。在源代码文章中 <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">source article</a> 列出了用于呈现函数值及其方向导数的具体计算。对于第二个示例，感兴趣的读者参阅此链接 <a href="http://www.columbia.edu/~ahd2125/post/2015/12/5/">this article</a>。 
 *
 *首先考虑任何复合函数 $f(x)$，这里表示为有两个独立变量，可以被分解成它的基本函数的组合 
 * @f[
 *   f (\mathbf{x})
 *   = f_{0} \circ f_{1} \circ f_{2} \circ \ldots \circ f_{n} (\mathbf{x})
 *   \quad .
 * @f]
 * 如前所述，如果每一个基算子$f_{n}$都是光滑可微的，那么链式法则可以普遍用于计算$f$的总导数，即$\dfrac{d f(x)}{d \mathbf{x}}$。“正向”与“反向”模式的区别在于链式规则是如何计算的，但最终两者都计算总导数
 * @f[
 *   \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \dfrac{d f_{0}}{d f_{1}} \dfrac{d f_{1}}{d f_{2}} \dfrac{d f_{2}}{d f_{3}} \ldots \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}}
 *   \quad .
 * @f]
 *
 *在正向模式下，链式规则是从“内到外”自然计算出来的。因此自变量是固定的，每个子函数$f'_{i} \vert_{f'_{i+1}}$递归计算，其结果作为输入返回给父函数。使用括号封装和固定操作顺序，这意味着我们计算
 * @f[
 *   \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \dfrac{d f_{0}}{d f_{1}} \left( \dfrac{d f_{1}}{d f_{2}} \left(\dfrac{d f_{2}}{d f_{3}} \left(\ldots \left( \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)\right)\right)\right)
 *   \quad .
 * @f]
 *前向扫描的计算复杂度与输入函数的计算复杂度成正比。然而，对于要计算的每个方向导数，需要计算图的一次扫描。
 *
 *在反向模式下，链式规则是从“由外而内”的角度计算出来的，有些不自然。首先计算并固定因变量的值，然后计算前面的微分运算，并从左到右依次与前面的结果相乘。同样，如果我们使用括号封装并固定操作顺序，这意味着反向计算由
 * @f[
 * \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \left( \left( \left( \left( \left( \dfrac{d f_{0}}{d f_{1}} \right) \dfrac{d f_{1}}{d f_{2}} \right) \dfrac{d f_{2}}{d f_{3}} \right) \ldots \right) \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)
 *   \quad .
 * @f]
 * 中间值$\dfrac{d f_{i-1}}{d f_{i}}$ 被叫做 <em>adjoints</em>, 必须在遍历计算图时计算并存储。然而，对于每个相依的标量函数，计算图的一次扫描会同时呈现所有的方向导数。
 *
 *总的来说，每种模式的效率由独立（输入）变量和相依（输出）变量的数量决定。如果输出在数量上大大超过输入，那么正向模式可以比反向模式更有效。当输入变量的数量大大超过输出变量的数量时，情况正好相反。这一点可以用来帮助通知哪种数字类型最适合使用自动微分执行哪一组操作。例如，在许多要计算二阶导数的应用中，结合反向和正向模式是合适的。前者通常用于计算一阶导数，后者用于计算二阶导数。
 *
 *
 * @subsection auto_diff_1_1 已支持 automatic differentiation 现有库
 *
 *我们目前已经验证了以下一些类型和组合的实现：
 *
 *  - Taped ADOL-C (n-differentiable, in theory, but internal drivers for up to second-order
 *    derivatives will be implemented)
 *  - Tapeless ADOL-C (once differentiable)
 *  - Forward-mode Sacado with dynamic memory allocation using expression templates (once differentiable)
 *  - Nested forward-mode Sacado using expression templates (twice differentiable)
 *  - Reverse-mode Sacado (once differentiable)
 *  - Nested reverse and dynamically-allocated forward-mode Sacado (twice differentiable, but results memory leak described in Differentiation::AD::NumberTypes)
 *
 *注意，在上面，“动态内存分配”指的是在编译时不需要指定自变量的数量。 
 *
 * ADOL-C 用户手册 <a href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">ADOL-C user manual</a>
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 *   title     = {Getting Started with ADOL-C},
 *   author    = {Walther, Andrea and Griewank, Andreas},
 *   year      = {2009},
 *   booktitle = {Combinatorial scientific computing},
 *   doi       = {10.1.1.210.4834},
 *   pages     = {181--202},
 *   url       = {https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf}
 * }
 * @endcode
 *
 *提供了对其磁带化和无磁带化实现的基本见解，以及如何将ADOL-C合并到用户代码中。对于理解ADOL-C的实现以及如何在数字代码中使用ADOL-C的可能性，还有一些有用的资源，包括：
 *
 * @code{.bib}
 * @Article{Griewank1996a,
 *   author    = {Griewank, Andreas and Juedes, David and Utke, Jean},
 *   title     = {Algorithm 755: {ADOL-C}: a package for the automatic differentiation of algorithms written in {C/C++}},
 *   journal   = {ACM Transactions on Mathematical Software (TOMS)},
 *   year      = {1996},
 *   volume    = {22},
 *   number    = {2},
 *   pages     = {131--167},
 *   doi       = {10.1145/229473.229474},
 *   publisher = {ACM}
 * }
 * @endcode
 * @code{.bib}
 * @InCollection{Bischof2008a,
 *   author =    {Bischof, Christian and Guertler, Niels and Kowarz, Andreas and Walther, Andrea},
 *   title =     {Parallel reverse mode automatic differentiation for OpenMP programs with ADOL-C},
 *   booktitle = {Advances in Automatic Differentiation},
 *   publisher = {Springer},
 *   year =      {2008},
 *   pages =     {163--173}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Kulshreshtha2012a,
 *   chapter   = {Computing Derivatives in a Meshless Simulation Using Permutations in {ADOL}-C},
 *   pages     = {321--331},
 *   title     = {Recent Advances in Algorithmic Differentiation},
 *   publisher = {Springer Berlin Heidelberg},
 *   year      = {2012},
 *   author    = {Kshitij Kulshreshtha and Jan Marburger},
 *   editor    = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 *   series    = {Lecture Notes in Computational Science and Engineering},
 *   doi       = {10.1007/978-3-642-30023-3_29},
 * }
 * @endcode
 * @code{.bib}
 * @InProceedings{Kulshreshtha2013a,
 *   author    = {Kulshreshtha, Kshitij and Koniaeva, Alina},
 *   title     = {Vectorizing the forward mode of ADOL-C on a GPU using CUDA},
 *   booktitle = {13th European AD Workshop},
 *   year      = {2013},
 *   month     = jun
 * }
 * @endcode
 *
 *类似地，为理解Sacado数字类型的实现（特别是如何使用和利用表达式模板）而选择的有用资源包括：
 *
 * @code{.bib}
 * @InCollection{Bartlett2006a,
 *   author        = {Bartlett, R. A. and Gay, D. M. and Phipps, E. T.},
 *   title         = {Automatic Differentiation of C++ Codes for Large-Scale Scientific Computing},
 *   booktitle     = {International Conference on Computational Science {\textendash} {ICCS} 2006},
 *   publisher     = {Springer Berlin Heidelberg},
 *   year          = {2006},
 *   editor        = {Alexandrov, V.N. and van Albada, G.D. and Sloot, P.M.A. amd Dongarra, J.},
 *   pages         = {525--532},
 *   doi           = {10.1007/11758549_73},
 *   organization  = {Springer}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Gay2012a,
 *   chapter   = {Using expression graphs in optimization algorithms},
 *   pages     = {247--262},
 *   title     = {Mixed Integer Nonlinear Programming},
 *   publisher = {Springer New York},
 *   year      = {2012},
 *   author    = {Gay, D. M.},
 *   editor    = {Lee, J. and Leyffer, S.},
 *   isbn      = {978-1-4614-1927-3},
 *   doi       = {10.1007/978-1-4614-1927-3_8}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Phipps2012a,
 *   chapter     = {Efficient Expression Templates for Operator Overloading-based Automatic Differentiation},
 *   pages       = {309--319},
 *   title       = {Recent Advances in Algorithmic Differentiation},
 *   publisher   = {Springer},
 *   year        = {2012},
 *   author      = {Eric Phipps and Roger Pawlowski},
 *   editor      = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 *   series      = {Lecture Notes in Computational Science and Engineering},
 *   volume      = {73},
 *   date        = {2012-05-15},
 *   doi         = {10.1007/978-3-642-30023-3_28},
 *   eprint      = {1205.3506v1},
 *   eprintclass = {cs.MS},
 *   eprinttype  = {arXiv}
 * }
 * @endcode
 *
 *正向和反向模式Sacado数的实现相当复杂。从trilinos12.12开始，数学运算的实现涉及到大量的预处理器指令和宏编程。因此，代码可能很难理解，而且这些类没有有意义的配套文档。因此，可以在 <a href="https://trilinos.org/docs/dev/packages/sacado/doc/html/classSacado_1_1Fad_1_1SimpleFad.html">this link for the Sacado::Fad::SimpleFad class</a> 的链接中找到理解这些数字的原理实现的有用资源，该类概述了不使用表达式模板的前向模式自动可微数的引用（尽管据报道效率很低）实现(虽然没有明确说明，但似乎Sacado::Fad::SimpleFad类是以双数的精神实现的。）
 *
 * @subsection auto_diff_1_2 如何将自动差异化集成到deal.II中
 *
 *由于每个自动微分库的接口差别很大，不久的将来将为每个数字建立一个统一的内部接口。我们的目标是允许一些驱动程序类（提供核心功能，稍后将在下一节中介绍）具有与不同的自动差异库交互的一致机制。具体来说，他们需要能够正确初始化和最终确定数据，这些数据将被解释为公式的因变量和自变量。
 *
 *实现支持的自动可微数接口的文件摘要如下：
 *
 * - ad_drivers.h: 提供充当内部支持的自动差异库接口驱动程序的类。它们在内部用作我们提供的helper类的中介。
 * - ad_helpers.h: 提供一组类来帮助在许多不同的上下文中执行自动区分。 具体细节请查看 @ref auto_diff_1_3.
 * - ad_number_types.h: 为驱动程序类支持的自动可微数组合引入枚举（称为类型代码）。下面讨论使用这种有点限制的机制背后的理由。
 * - ad_number_traits.h: 声明一些内部类，这些类将针对每个自动区分库和/或数字类型进行专门化。它们随后通过NumberTraits和ADNumberTraits类为类提供统一的接口，这些类在整个驱动程序中广泛使用。我们还提供了一些机制来方便地查询这些数字的select属性，即一些类型特征。
 * - adolc_math.h: ADOL-C数学运算的扩展，允许在整个库中一致地使用这些数字。
 * - adolc_number_types.h: 内部类的实现，定义了如何使用ADOL-C数。
 * - adolc_product_types.h: 定义一些乘积和标量类型，允许将ADOL-C数与 Tensor 和 SymmetricTensor 类结合使用。
 * - sacado_math.h: Sacado数学操作的扩展，允许在整个库中一致地使用这些数字。
 * - sacado_number_types.h: 内部类的实现，定义如何使用支持的Sacado数。
 * - sacado_product_types.h: 定义一些乘积和标量类型，允许将支持的Sacado数与Tensor和SymmetricTensor类结合使用。
 *
 *通过对每个支持的数字类型使用类型代码，我们人为地限制了库中可以使用的自动可微数字的类型。这种设计选择是由于这样一个事实，即确保每个数字类型都正确初始化，并且嵌套（模板化）类型的所有组合对于库执行的所有操作都保持有效，这一点非常重要。此外，库中还有一些冗长的函数，它们被实例化为支持的数字类型，并且只有在使用库知道的自动可微数时才满足内部检查。这再次确保了所有计算的完整性。最后，使用一个简单的枚举作为类模板参数，最终可以很容易地在生产代码中使用的类型之间进行切换，而不需要对用户代码进行任何进一步的修改。
 *
 *
 * @subsubsection auto_diff_1_3 自动微分库的用户界面
 *
 *deal.II库为我们支持的自动差异化库提供了统一的接口。迄今为止，已经为以下上下文开发了 helper 类：
 *
 * - 设计在正交点水平（或任何一般连续点）运行的类:
 *   - Differentiation::AD::ScalarFunction: %标量值函数的微分。一个典型的应用是直接从应变能函数发展本构关系。
 *   - Differentiation::AD::VectorFunction: %向量值函数的微分。这可以用来线性化本构关系的运动变量，或协助求解局部内变量的演化方程。
 * - 设计用于单元级操作的类：:
 *   - Differentiation::AD::EnergyFunctional: %标量值能量泛函的微分，如由变分公式产生的微分。
 *   - Differentiation::AD::ResidualLinearization: %向量值有限元残差的微分，导致其一致线性化。
 *
 *当然，用户也可以自己管理初始化和派生计算。
 *
 *关于如何使用ADOL-C实现这一点的最新例子见下面几个链接
 * - 用户手册 <a href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">user manual</a>,
 * - github项目 <a href="https://github.com/coin-or/ADOL-C">development repository</a>, and
 * - deal.ii 测试例子 <a href="https://github.com/dealii/dealii/tree/master/tests/adolc">test-suite</a>,
 *
 *而对于Sacado，可以在
 * - <a href="https://github.com/trilinos/Trilinos/tree/master/packages/sacado/example">development repository</a>,
 * - deal.ii 展示例子 <a href="https://github.com/dealii/code-gallery/tree/master/Quasi_static_Finite_strain_Compressible_Elasticity">code-gallery example</a>, and
 * - deal.ii 测试例子 <a href="https://github.com/dealii/dealii/tree/master/tests/sacado">test-suite</a>.
 *
 *
 * @section symb_diff_1 基于符号表达的微分
 *
 * <a href="https://en.wikipedia.org/wiki/Symbolic_differentiation">Symbolic differentiation</a> 符号微分在设计和使用上与自动微分有很大的不同。任何符号库的底层都是一个计算机代数系统（CAS），它实现了一种语言和一组算法来操作符号（或“类似字符串”）表达式。从哲学的角度来看，这与代数运算是如何用手工进行的最为相似。
 
 
 *为了更好地区分符号微分法和数值方法，比如自动微分法，让我们考虑一个非常简单的例子。例如给定一个函数 $f(x,y) = [2x+1]^{y}$, 其中 $x$ and $y$ 相互独立。通过应用链式法则，这个函数的导数可以简单地表示出来 $\dfrac{d f(x,y)}{d x} = 2y[2x+1]^{y-1}$ 和 $\dfrac{d f(x,y)}{d y} = [2x+1]^{y} \ln(2x+1)$.这些正是在定义符号变量x和y之后，从CAS得到的结果，定义符号表达式 `f = pow(2x+1, y)` 然后计算导数 `diff(f, x)` 和 `diff(f, y)`. 在这一点上，没有假设x和y代表什么；它们后来可能被解释为普通（标量）数、复数或幂函数和自然对数函数定义良好的其他数。显然，这意味着也没有关于表达式或其导数的计算点的假设。我们可以很容易地得到$\dfrac{d f(x, y)}{d x}$的表达式，并在$x=1, y=2.5$时对其求值，然后在不重新计算导数表达式本身的情况下，在$x=3.25, y=-6$. 事实上，任何符号变量或表达式的解释，以及变量之间的相互依赖性，都可以在其操作过程中的任何时候定义或重新定义；这使得计算具有一定程度的灵活性，而这种灵活性是自动微分所不能比拟的。例如，可以执行$g(x) = \dfrac{d f(x, y)}{d x} \vert_{y=1}$，然后为 $x$ 的几个不同值重新计算$g(x)$。也可以用后后验来表达 `x` 和 `y` 之间的相互依赖关系,例如 $y \rightarrow y(x) := 2x$ 。对于这种情况，这意味着最初计算的导数 $\dfrac{d f(x, y)}{d x} \rightarrow \dfrac{\partial f(x, y(x))}{\partial x} = 2y(x) [2x+1]^{y(x)-1} = 4x[2x+1]^{2x-1}$ 和 $\dfrac{d f(x, y)}{d y} \rightarrow \dfrac{\partial f(x, y(x))}{\partial y} = [2x+1]^{y(x)} \ln(2x+1) = [2x+1]^{2x} \ln(2x+1)$ 真正代表偏导数而不是全导数。当然，如果在计算导数 $\dfrac{d f(x, y(x))}{d x}$ 和 $\dfrac{d f(x, y(x))}{d y}$ 之前明确定义了这种相互依赖关系，这才可能对应于总导数（这是本例中自动微分能够实现的唯一结果）。
*
*由于复杂的CAS形成符号运算的基础，操作的类型不一定仅限于微分，而是可以跨越与离散微分学有关的操作谱，纯数学中的主题，等等。SymPy库 <a href="https://www.sympy.org/en/index.html">SymPy</a> 的文档提供了大量的例子，突出了一个成熟的CAS的能力。通过 Differentiation::SD::Expression 类, 以及 Differentiation::SD 命名空间的相关函数, 我们为高性能SymEngine <a href="https://github.com/symengine/symengine">SymEngine</a> 符号操作库提供了一个包装器，它丰富了运算符重载，并提供了一个一致的接口，使其易于使用。事实上，在许多情况下，这个类可以作为算术类型的“插入式”替换，将运算从数字性质转换为符号性质；当在基础数字类型上模板化类时，这一点变得特别容易。由于专注于偏微分方程的数值模拟，deal.II中公开的CAS功能侧重于符号表达的创建、操作和区分。
*
*SymEngine功能的便利包装器主要集中在仅涉及基于字典（即，类似于“基于字符串”）操作的操作上。尽管 Symeengine 以高效的方式执行这些操作，但是仍然知道它们的计算开销很大，特别是在对大型表达式执行这些操作时。因此，当在生产代码中使用此代码时，应该预期执行差分、符号替换等的代码部分的性能可能是一个限制因素。因此，deal.II提供了一个接口，可以通过 @p BatchOptimizer 类（通常利用 SymEngine 提供的功能）加速对冗长符号表达式的求值。具体而言, @p BatchOptimizer 使用公共子表达式消除（CSE）等方法同时优化符号表达式集合，以及通过使用自定义生成的 `std::function` 函数或使用LLVM JIT编译器编译表达式来生成高性能代码路径来计算这些表达式。
*
*作为最后一点，重要的是要认识到deal.II当前对支持的符号库接口的实现中仍然存在的主要缺陷。目前实现的功能级别有效地将符号代数的使用限制在传统的用例中（即标量和张量代数，这可能有助于定义本构关系或复杂函数以作为边界条件或源项）。在将来，我们还将实现类来帮助执行汇编操作，其思想与 Differentication::AD 命名空间中所做的相同。
*
 * 实现支持的符号可微数接口的文件总结如下：
 * - symengine_math.h: 数学运算的实现，允许实现符号表达式的类在整个库和用户代码中一致地使用。它为标准名称空间中的许多数学函数提供了对应的定义。
 * - symengine_number_traits.h: 提供一些机制来轻松查询符号数的select属性，即某些类型特征。
 * - symengine_number_types.h: Differentiation::SD::Expression类的实现，可用于表示标量符号变量、标量符号表达式等。这个表达式类被赋予了一整套重载的运算符，用于SymEngine库支持的所有数学和逻辑运算，并且在数值建模的上下文中被认为是有用的。
 * - symengine_optimizer.h: 执行 Differentiation::SD::BatchOptimizer 类，可用于使用各种技术加速（在某些情况下，显著地）符号表达式的计算。
 * - symengine_product_types.h: 定义一些乘积和标量类型，允许将符号表达式与 Tensor 和 SymmetricTensor 类结合使用。 
 * - symengine_scalar_operations.h: 定义可以对标量符号表达式或变量执行或使用标量符号表达式或变量执行的许多操作。这包括（但不限于）标量符号的创建、对标量执行微分以及标量表达式中的符号替换。
 * - symengine_tensor_operations.h: 定义可以对符号表达式或变量的张量执行或使用张量执行的许多操作。这包括（但不限于）符号张量的创建、针对符号张量执行微分、符号张量的微分以及张量表达式内的符号替换。
 * - symengine_types.h: 为符号计算上下文中常用的某些类型提供别名。
 * - symengine_utilities.h: 提供一些在符号计算上下文中有用的实用函数。
 */
