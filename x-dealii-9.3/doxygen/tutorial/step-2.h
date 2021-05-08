  /**   @page step_2 The step-2 tutorial program   

本教程取决于  step-1  。

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Sparsity"> Sparsity </a><a href="#Sparsity"> Sparsity </a>
        <li><a href="#Howdegreesoffreedomareenumerated"> How degrees of freedom are enumerated </a><a href="#Howdegreesoffreedomareenumerated"> How degrees of freedom are enumerated </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Meshgeneration">Mesh generation</a><a href="#Meshgeneration">Mesh generation</a>
        <li><a href="#CreationofaDoFHandler">Creation of a DoFHandler</a><a href="#CreationofaDoFHandler">Creation of a DoFHandler</a>
        <li><a href="#RenumberingofDoFs">Renumbering of DoFs</a><a href="#RenumberingofDoFs">Renumbering of DoFs</a>
        <li><a href="#Themainfunction">The main function</a><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly 

<a name="Intro"></a> <a name="Introduction"></a><h1>Introduction</h1>


  @dealiiVideoLecture{9}   

在前面的例子中，我们已经创建了一个网格，现在我们展示如何在这个网格上定义自由度。在这个例子中，我们将使用最低阶的( $Q_1$ )有限元，自由度与网格的顶点相关。以后的例子将展示更高阶的元素，自由度不一定与顶点相关，但可以与边、面或单元相关。

术语 "自由度 "在有限元界通常用来表示两个略有不同但相关的东西。首先，我们希望将有限元解表示为形状函数的线性组合，形式为 $u_h(\mathbf x) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf
x)$  。这里， $U_j$ 是一个膨胀系数的向量。因为我们还不知道它们的值（我们将计算它们作为线性或非线性系统的解），它们被称为 "未知数 "或 "自由度"。该术语的第二个含义可以解释如下。对有限元问题的数学描述通常是说，我们正在寻找一个满足某些方程组的有限维函数 $u_h \in V_h$ （例如， $a(u_h,\varphi_h)=(f,\varphi_h)$ 为所有测试函数 $\varphi_h\in
V_h$ ）。换句话说，我们在这里说的是，解决方案需要位于某个空间  $V_h$  中。然而，为了在计算机上实际解决这个问题，我们需要选择这个空间的一个基；这就是我们在上面用系数 $U_j$ 展开 $u_h(\mathbf x)$ 时使用的形状函数 $\varphi_j(\mathbf x)$ 的集合。当然，空间 $V_h$ 的基数有很多，但我们将特别选择由传统上在网格单元上局部定义的有限元函数描述的基数。在这种情况下描述 "自由度 "需要我们简单地 <i>enumerate</i> 空间的基函数  $V_h$  。对于 $Q_1$ 元素，这意味着简单地以某种方式列举网格的顶点，但对于高阶元素，还必须列举与网格的边、面或单元内部相关的形状函数。换句话说，自由度的枚举是一个完全独立于我们用于顶点的索引的东西。提供这种列举 $V_h$ 的基础函数的类被称为DoFHandler。

在一个网格上定义自由度（简称 "DoF"）是一个相当简单的任务，因为这个库为你做了所有的工作。基本上，你所要做的就是创建一个有限元对象（从deal.II已有的众多有限元类中选取，例如参见 @ref fe 文档），并通过 DoFHandler::distribute_dofs 函数将其交给DoFHandler对象（"分配DoF "是我们用来描述<i>enumerating</i>基础函数过程的术语，如上所述）。DoFHandler是一个知道哪些自由度住在哪里的类，也就是说，它可以回答 "全局有多少个自由度 "和 "在这个单元上，给我住在这里的形状函数的全局索引 "这样的问题。当你决定你的系统矩阵应该有多大时，以及当把单个单元的贡献复制到全局矩阵时，你需要这种信息。

<a name="Sparsity"></a><h3> Sparsity </h3>


然后，下一步将是使用这个有限元和网格来计算对应于特定微分方程的矩阵和右手。我们将为 step-3 程序保留这一步骤，而是谈论有限元程序的一个实际问题，即有限元矩阵总是非常稀疏的：这些矩阵中几乎所有条目都是零。

更准确地说，如果矩阵中的非零项<i>per row</i>的数量与整个自由度的数量无关，我们就说该矩阵是稀疏的。例如，拉普拉斯方程的有限差分近似的简单5点模版导致了一个稀疏矩阵，因为每行的非零条目数是5，因此与矩阵的总大小无关。对于更复杂的问题--例如， step-22 的斯托克斯问题--特别是在三维中，每行的条目数可能是几百个。但重要的一点是，这个数字与问题的总体大小无关：如果你细化网格，每行的最大未知数保持不变。

稀疏性是有限元方法的显著特征之一，与之相比，例如使用泰勒扩展和匹配系数来逼近偏微分方程的解，或使用傅里叶基础。

在实践中，正是由于矩阵的稀疏性，使我们能够解决有数百万或数十亿未知数的问题。为了理解这一点，请注意，一个有 $N$ 行的矩阵，每个非零项的数量都有一个固定的上限，需要 ${\cal O}(N)$ 个内存位置来存储，而矩阵-向量乘法也只需要 ${\cal O}(N)$ 次操作。因此，如果我们有一个线性求解器，只需要固定数量的矩阵-向量乘法就可以得出这个矩阵的线性系统的解，那么我们就会有一个求解器，可以以最佳的复杂度找到所有 $N$ 个未知数的值，即总共需要 ${\cal O}(N)$ 次运算。很明显，如果矩阵不是稀疏的，这是不可能的（因为那样的话，矩阵中的条目数必须是 ${\cal O}(N^s)$ 与一些 $s>1$ ，做固定数量的矩阵-向量乘积将需要 ${\cal O}(N^s)$ 次操作），但这也需要非常专业的求解器，如多网格方法，以满足求解只需要固定数量的矩阵-向量乘法的要求。在本教程的其余程序中，我们将经常研究使用何种求解器的问题。

稀疏性是由以下事实产生的：有限元形状函数是在单个单元上定义的<i>locally</i>，而不是全局的，而且双线性形式中的局部微分算子只对支持度重叠的形状函数进行耦合。一个函数的 "支持 "是指它的非零区域。对于有限元方法，形状函数的支持通常是指与它所定义的顶点、边或面相邻的单元。) 换句话说，不在同一单元上定义的自由度 $i$ 和 $j$ 不会重叠，因此，矩阵条目 $A_{ij}$ 将为零。 (在某些情况下，如非连续加尔金法，形状函数也可以通过面积分连接到相邻的单元。但有限元方法一般不会将形状函数与定义了该函数的单元的近邻相联系。) 


<a name="Howdegreesoffreedomareenumerated"></a><h3> How degrees of freedom are enumerated </h3>


默认情况下，DoFHandler类以一种相当随机的方式列举网格上的自由度；因此，稀疏模式也没有为任何特定目的进行优化。为了说明这一点，下面的代码将演示一个简单的方法来输出对应于DoFHandler的 "稀疏模式"，即一个对象代表了在网格上离散偏微分方程时可能建立的矩阵的所有潜在的非零元素及其DoFHandler。这种缺乏结构的稀疏模式将从我们下面展示的图片中显现出来。

对于大多数应用和算法来说，自由度的确切编号方式并不重要。例如，我们用来解决线性系统的共轭梯度方法并不关心。另一方面，有些算法确实关心：特别是一些预处理程序，如SSOR，如果它们能以特定的顺序走过自由度，就能更好地工作，如果我们能以这样的方式排序，使SSOR能以这样的顺序从零到 $N$ 迭代它们，那就太好了。其他的例子包括计算不完整的LU或Cholesky因式分解，或者如果我们关心矩阵的块结构（见 step-20 的例子）。因此，deal.II在命名空间DoFRenumbering中有可以以特定方式重新列举自由度的算法。重新编号可以被认为是选择了一个不同的、排列过的有限元空间的基础。因此，这种重新编号所产生的稀疏模式和矩阵与我们没有明确的重新编号所得到的行和列相比，也只是一种排列组合。

在下面的程序中，我们将使用Cuthill和McKee的算法来做这件事。我们将在下面的<a href="#Results">results section</a>中展示原始自由度列举和重新编号版本的稀疏模式。<a name="CommProg"></a> <h1> The commented program</h1>

前面的几个包括就像前面的程序一样，所以不需要额外的注释。

@code
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>


@endcode 



然而，下一个文件是新的。我们需要这个包含文件来将自由度（"DoF"）与顶点、线和单元联系起来。

@code
#include <deal.II/dofs/dof_handler.h>


@endcode 



下面的包含文件包含了对双线性有限元的描述，包括它在三角形的每个顶点上都有一个自由度，但在面和单元的内部没有自由度。




事实上，该文件包含了对Lagrange元素的一般描述，即还有二次、三次等版本，而且不仅是2d，还有1d和3d。

@code
#include <deal.II/fe/fe_q.h>
@endcode 



在以下文件中，可以找到几个操纵自由度的工具。

@code
#include <deal.II/dofs/dof_tools.h>
@endcode 



我们将使用一个稀疏矩阵来可视化自由度在网格上的分布所产生的非零项的模式。这个类目可以在这里找到。

@code
#include <deal.II/lac/sparse_matrix.h>
@endcode 



我们还将需要使用一个中间的稀疏模式结构，在这个文件中可以找到。

@code
#include <deal.II/lac/dynamic_sparsity_pattern.h>


@endcode 



我们要使用一种特殊的算法对自由度进行重新编号。它在这里被声明。

@code
#include <deal.II/dofs/dof_renumbering.h>


@endcode 



而这也是C++输出所需要的。

@code
#include <fstream>


@endcode 



最后，和 step-1 一样，我们将deal.II命名空间导入到全局范围。

@code
using namespace dealii;


@endcode 




<a name="Meshgeneration"></a> <h3>Mesh generation</h3>




这是之前 step-1 示例程序中产生圆形网格的函数，细化步骤较少。唯一的区别是它通过其参数返回它所产生的网格。

@code
void make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);


  for (unsigned int step = 0; step < 3; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        for (const auto v : cell->vertex_indices())
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));


            if (std::fabs(distance_from_center - inner_radius) <=
                1e-6 * inner_radius)
              {
                cell->set_refine_flag();
                break;
              }
          }


      triangulation.execute_coarsening_and_refinement();
    }
}


@endcode 




<a name="CreationofaDoFHandler"></a> <h3>Creation of a DoFHandler</h3>




到目前为止，我们只有一个网格，即一些几何信息（顶点的位置）和一些拓扑信息（顶点如何连接到线，线如何连接到单元，以及哪些单元与其他单元相邻）。要使用数值算法，还需要一些逻辑信息：我们希望将自由度数字与每个顶点（或线，或单元，如果我们使用高阶元素的话）联系起来，以便随后生成描述三角形上的有限元场的矩阵和向量。




这个函数显示了如何做到这一点。要考虑的对象是 <code>DoFHandler</code> 类模板。 然而，在我们这样做之前，我们首先需要一些东西来描述有多少自由度要与这些对象中的每一个相关联。由于这是有限元空间定义的一个方面，有限元基类存储了这个信息。因此，在目前情况下，我们创建了一个描述拉格朗日元素的派生类 <code>FE_Q</code> 的对象。它的构造函数需要一个参数，说明元素的多项式程度，这里是1（表示双线性元素）；这就对应于每个顶点的一个自由度，而在线和四边形内部没有自由度。如果给构造函数的值是3，我们就会得到一个双立方体元素，每个顶点有一个自由度，每条线有两个自由度，单元内有四个自由度。一般来说， <code>FE_Q</code> 表示具有完整多项式（即张量积多项式）的连续元素家族，直到指定的阶数。




我们首先需要创建一个该类的对象，然后将其传递给 <code>DoFHandler</code> 对象，为自由度分配存储空间（用交易二的行话说：我们<i>distribute degrees of
freedom</i>）。

@code
void distribute_dofs(DoFHandler<2> &dof_handler)
{
  const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);


@endcode 



现在我们已经将自由度与每个顶点的全局数字联系起来，我们想知道如何将其可视化？ 没有简单的方法可以直接将与每个顶点相关的自由度数字可视化。然而，这样的信息几乎不会真正重要，因为编号本身或多或少是任意的。还有更重要的因素，我们将在下文中展示其中一个。   


与三角形的每个顶点相关的是一个形状函数。假设我们想解决类似拉普拉斯方程的问题，那么不同的矩阵条目将是每一对此类形状函数的梯度的积分。显然，由于形状函数只在与它们相关的顶点相邻的单元格上是非零的，所以只有当与该列和行%号相关的形状函数的支持相交时，矩阵条目才是非零的。这只是相邻形状函数的情况，因此也只是相邻顶点的情况。现在，由于顶点被上述函数 (DoFHandler::distribute_dofs), 或多或少地随机编号，矩阵中非零项的模式将有些参差不齐，我们现在就来看一下。   


首先，我们必须创建一个结构，用来存储非零元素的位置。然后，这可以被一个或多个稀疏矩阵对象使用，这些对象在这个稀疏模式所存储的位置上存储条目的值。存储这些位置的类是SparsityPattern类。然而，事实证明，当我们试图立即填充这个类时，它有一些缺点：它的数据结构的设置方式是，我们需要对我们可能希望在每一行的最大条目数有一个估计。在两个空间维度上，通过 DoFHandler::max_couplings_between_dofs() 函数可以得到合理的估计值，但是在三个维度上，该函数几乎总是严重高估真实的数字，导致大量的内存浪费，有时对于使用的机器来说太多，即使未使用的内存可以在计算稀疏模式后立即释放。为了避免这种情况，我们使用了一个中间对象DynamicSparsityPattern，该对象使用了一个不同的%内部数据结构，我们可以随后将其复制到SparsityPattern对象中，而不需要太多的开销。关于这些数据结构的一些更多信息可以在 @ref Sparsity 模块中找到）。为了初始化这个中间数据结构，我们必须给它提供矩阵的大小，在我们的例子中，矩阵是正方形的，行和列的数量与网格上的自由度相同。

@code
  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());


@endcode 



然后我们用非零元素的位置来填充这个对象，鉴于目前自由度的编号，非零元素将位于这个对象中。

@code
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);


@endcode 



现在我们已经准备好创建实际的稀疏模式，以后可以用于我们的矩阵。它将包含已经在DynamicSparsityPattern中集合的数据。

@code
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);


@endcode 



有了这个，我们现在可以把结果写到一个文件里。

@code
  std::ofstream out("sparsity_pattern1.svg");
  sparsity_pattern.print_svg(out);
@endcode 



结果被存储在一个 <code>.svg</code> 文件中，矩阵中的每个非零条目都对应于图像中的一个红色方块。输出结果将显示在下面。   


如果你看一下，你会注意到稀疏性模式是对称的。这不应该是一个惊喜，因为我们没有给 <code>DoFTools::make_sparsity_pattern</code> 任何信息，表明我们的双线性形式可能以非对称的方式耦合形状函数。你还会注意到，它有几个明显的区域，这源于编号从最粗的单元开始，然后到较细的单元；由于它们都是围绕原点对称分布的，这在稀疏模式中又显示出来。

@code
}



@endcode 




<a name="RenumberingofDoFs"></a><h3>Renumbering of DoFs</h3>




在上面产生的稀疏性模式中，非零条目从对角线上延伸得相当远。对于某些算法，例如不完全LU分解或Gauss-Seidel预处理，这是不利的，我们将展示一个简单的方法来改善这种情况。




请记住，为了使矩阵中的一个条目 $(i,j)$ 不为零，形状函数i和j的支持需要相交（否则在积分中，积分项将到处为零，因为在某个点上，一个或另一个形状函数为零）。然而，形状函数的支持只有在它们彼此相邻的情况下才会相交，所以为了让非零条目聚集在对角线周围（其中 $i$ 等于 $j$ ），我们希望相邻的形状函数的索引（DoF编号）相差不太大。




这可以通过一个简单的前行算法来实现，即从一个给定的顶点开始，给它的索引为0。然后，它的邻居被连续编号，使它们的指数接近于原始指数。然后，他们的邻居，如果还没有被编号，也被编号，以此类推。




Cuthill和McKee的算法是一种沿着这些思路增加了一点复杂性的算法。我们将在下面的函数中使用它来对自由度进行重新编号，从而使产生的稀疏模式在对角线周围更加局部化。该函数唯一有趣的部分是对 <code>DoFRenumbering::Cuthill_McKee</code> 的第一次调用，其余部分基本与之前一样。

@code
void renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);


  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);


  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);


  std::ofstream out("sparsity_pattern2.svg");
  sparsity_pattern.print_svg(out);
}


@endcode 



同样，输出结果显示如下。请注意，非零项在对角线附近的聚集情况比以前好得多。这种效果对于较大的矩阵来说甚至更加突出（目前的矩阵有1260行和列，但大的矩阵往往有几十万行）。




值得注意的是， <code>DoFRenumbering</code> 类也提供了一些其他的算法来重新编号自由度。例如，如果所有的耦合都在矩阵的下三角或上三角部分，那当然是最理想的，因为那样的话，解决线性系统就只需要向前或向后替换。当然，这对于对称稀疏模式来说是无法实现的，但在一些涉及传输方程的特殊情况下，通过列举从流入边界沿流线到流出边界的自由度，这是可能的。毫不奇怪， <code>DoFRenumbering</code> 也有这方面的算法。










<a name="Themainfunction"></a> <h3>The main function</h3>




最后，这是主程序。它所做的唯一事情是分配和创建三角形，然后创建一个 <code>DoFHandler</code> 对象并将其与三角形相关联，最后对其调用上述两个函数。

@code
int main()
{
  Triangulation<2> triangulation;
  make_grid(triangulation);


  DoFHandler<2> dof_handler(triangulation);


  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
@endcode 

<a name="Results"></a><h1>Results</h1>


该程序运行后，产生了两个疏散模式。我们可以通过在网络浏览器中打开 <code>.svg</code> 文件来可视化它们。

结果是这样的（每一个点都表示一个可能为非零的条目；当然，这个条目是否真的为零取决于所考虑的方程，但矩阵中的指示位置告诉我们，在离散化局部，即微分方程时，哪些形状函数可以，哪些不可以耦合）。  <table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-1.svg" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-2.svg" alt=""></td>
  </tr>
</table>   

左图中的不同区域，由线条中的结点和左边和上面的单点表示，代表了三角化不同细化层次上的自由度。 从右图中可以看出，重新编号后，稀疏模式在矩阵的主对角线附近的聚类情况要好得多。虽然这可能不明显，但两张图片中的非零项数量当然是一样的。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3> 


就像 step-1 一样，你可能想玩一下程序，熟悉一下deal.II。例如，在 <code>distribute_dofs</code> 函数中，我们使用线性有限元（这就是FE_Q对象的参数 "1"）。如果你使用高阶元素，例如立方或五元元素（使用3和5作为各自的参数），探索稀疏模式如何变化。

你也可以通过细化网格来探索稀疏性模式的变化。你会发现不仅矩阵的大小会发生变化，而且其带宽也会发生变化（矩阵中离对角线最远的非零元素与对角线的距离），不过带宽与大小的比例通常会缩小，也就是说，矩阵在对角线周围聚集得更多。

另一个实验想法是尝试DoFRenumbering命名空间中除Cuthill-McKee之外的其他重编号策略，看看它们如何影响稀疏性模式。

你也可以用<a
href="http://www.gnuplot.info/">GNUPLOT</a>（较简单的可视化程序之一；也许不是最容易使用的，因为它是命令行驱动的，但在所有的Linux和其他类似Unix的系统上也是普遍可用的）通过改变 <code>print_svg()</code> to <code>print_gnuplot()</code> in <code>distribute_dofs()</code> and <code>renumber_dofs()</code>  来可视化输出。

@code
examples/\step-2> gnuplot


        G N U P L O T
        Version 3.7 patchlevel 3
        last modified Thu Dec 12 13:00:00 GMT 2002
        System: Linux 2.6.11.4-21.10-default


        Copyright(C) 1986 - 1993, 1998 - 2002
        Thomas Williams, Colin Kelley and many others


        Type `help` to access the on-line reference manual
        The gnuplot FAQ is available from
        http://www.gnuplot.info/gnuplot-faq.html


        Send comments and requests for help to <info-gnuplot@dartmouth.edu>
        Send bugs, suggestions and mods to <bug-gnuplot@dartmouth.edu>



Terminal type set to 'x11'
gnuplot> set style data points
gnuplot> plot "sparsity_pattern.1"
@endcode 



另一个基于<a href="http://www.gnuplot.info/">GNUPLOT</a>的做法是尝试打印出带有支撑点位置和编号的网格。为此，你需要包含GridOut和MappingQ1的头文件。这方面的代码是。

@code
  std::ofstream out("gnuplot.gpl");
  out << "plot '-' using 1:2 with lines, "
      << "'-' with labels point pt 2 offset 1,1"
      << std::endl;
  GridOut().write_gnuplot (triangulation, out);
  out << "e" << std::endl;
  const int dim = 2;
  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
                                        dof_handler,
                                        support_points);
  DoFTools::write_gnuplot_dof_support_point_info(out,
                                                 support_points);
  out << "e" << std::endl;
@endcode 

在我们运行代码后，我们得到一个叫做gnuplot.gpl的文件。要查看这个文件，我们可以在命令行中运行以下代码。

@code
gnuplot -p gnuplot.gpl
@endcode. 有了它，你会得到一个类似于 @image html support_point_dofs1.png 的图片，这取决于你正在看的网格。更多信息，请参见  DoFTools::write_gnuplot_dof_support_point_info.  <a name="PlainProg"></a> <h1> The plain program</h1>  @include "step-2.cc"  。 

  */  
