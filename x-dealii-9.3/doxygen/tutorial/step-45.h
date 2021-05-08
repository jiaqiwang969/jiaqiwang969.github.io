  /**   @page step_45 The step-45 tutorial program 。 

本教程取决于  step-6  。

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Settingupperiodicityconstraintsondistributedtriangulations">Setting up periodicity constraints on distributed triangulations</a><a href="#Settingupperiodicityconstraintsondistributedtriangulations">Setting up periodicity constraints on distributed triangulations</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly 

  <br>   

<i>This program was contributed by Daniel Arndt and Matthias Maier.</i> <a name="Intro"></a> <a name="Introduction"></a><h1>Introduction</h1>


在这个例子中，我们介绍了如何在deal.II中使用周期性边界条件。周期性边界条件是代数约束，通常发生在一个大域的代表性区域的计算中，在一个或多个方向上重复。

一个例子是模拟光子晶体的电子结构，因为它们有一个类似格子的结构，因此，往往只需要在格子的一个盒子上进行实际计算。为了能够以这种方式进行，我们必须假设该模型可以周期性地扩展到其他盒子；这就要求解具有周期性结构。

<a name="Procedure"></a> <a name="Procedure"></a><h1>Procedure</h1>


deal.II提供了一些高水平的入口来施加周期性边界条件。应用周期性边界条件的一般方法包括三个步骤（也可参见 @ref GlossPeriodicConstraints "关于周期性边界条件的词汇 "条目）。

-# 创建一个网格 

-# 使用 GridTools::collect_periodic_faces() 确定边界不同部分的一对面，这些面的解应该是对称的。 

-# 使用 parallel::distributed::Triangulation::add_periodicity() 将周期性信息添加到网格中。 

-# 使用 DoFTools::make_periodicity_constraints() 添加周期性约束。 

第二和第三步对于使用 parallel::distributed::Triangulation 类的并行网格是必要的，以确保位于域的对面但由周期性面连接的单元是幽灵层的一部分，如果其中一个单元被存储在本地处理器上。如果Triangulation不是 parallel::distributed::Triangulation, 类，这些步骤就没有必要。

第一步包括收集匹配的周期性面，并将它们存储在 GridTools::PeriodicFacePair. 的 <code>std::vector</code> 中，这是通过函数 GridTools::collect_periodic_faces() 完成的，例如可以这样调用。

@code
GridTools::collect_periodic_faces(dof_handler,
                                  b_id1,
                                  b_id2,
                                  direction,
                                  matched_pairs,
                                  offset = <default value>,
                                  matrix = <default value>,
                                  first_vector_components = <default value>);
@endcode 



这个调用在周期性边界上的容器dof_handler的所有面进行循环，边界指标分别为 @p b_id1 和 @p b_id2, 。(你可以在创建粗略的网格后手动分配这些边界指标，见 @ref GlossBoundaryIndicator  "边界指标"。另外，如果你指定了 "着色 "标志，你也可以让命名空间GridGenerator中的许多函数来做这件事；在这种情况下，这些函数会给边界的不同部分分配不同的边界指标，细节通常在这些函数的文档中详细说明）。) 

具体来说，如果 $\text{vertices}_{1/2}$ 是两个面 $\text{face}_{1/2}$ 的顶点，那么上面的函数调用将匹配成对的面（和道夫），使得 $\text{vertices}_2$ 和 $matrix\cdot \text{vertices}_1+\text{offset}$ 之间的差异在除方向之外的每个分量中都消失，并将产生的对与相关数据存储在 @p matched_pairs. 中（关于匹配过程的详细信息，参见 GridTools::orthogonal_equality() ）。

例如，考虑彩色单位方块 $\Omega=[0,1]^2$ ，其边界指标0在左边，1在右边，2在下面，3在上面的面。见 GridGenerator::hyper_cube() 的文件中关于如何分配边界指标的这个约定）。然后。

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs);
@endcode 

将产生周期性约束，即 $u(0,y)=u(1,y)$ 对所有 $y\in[0,1]$ 。

如果我们考虑由 $(0,0)$ ,  $(1,1)$ ,  $(1,2)$ ,  $(0,1)$ 的凸壳给出的平行四边形，我们可以通过指定一个 @p offset: 来实现约束 $u(0,y)=u(1,y+1)$  。 

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(0.,1.));
@endcode 

或 

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*arbitrary direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(1.,1.));
@endcode 

这里，同样，边界指标0和1的分配源于 GridGenerator::parallelogram() 的文件。

由此产生的 @p matched_pairs 可以用在 DoFTools::make_periodicity_constraints 中，用于用周期性约束填充AffineConstraints对象。

@code
DoFTools::make_periodicity_constraints(matched_pairs, constraints);
@endcode 



除了这个高级接口外，还有一些 DoFTools::make_periodicity_constraints 的变体可以结合这两个步骤（见 DofTools::make_periodicity_constraints). 的变体）。 

如果需要更多的灵活性，也有一个 DoFTools::make_periodicity_constraints 的低级接口。该低级变体允许直接指定两个应被约束的面。

@code
using namespace DoFTools;
make_periodicity_constraints(face_1,
                             face_2,
                             affine_constraints,
                             component_mask = <default value>;
                             face_orientation = <default value>,
                             face_flip = <default value>,
                             face_rotation = <default value>,
                             matrix = <default value>);
@endcode 

这里，我们需要使用 @p face_orientation,  @p face_flip 和 @p face_orientation. 来指定两个面的方向。要想了解更详细的描述，请看 DoFTools::make_periodicity_constraints. 的文档。除了自我解释的 @p component_mask 和 @p affine_constraints. ，其余参数与高层接口相同。 


<a name="problem"></a> <a name="Apracticalexample"></a><h1>A practical example</h1>


在下文中，我们将展示如何在一个更复杂的例子中使用上述函数。任务是对斯托克斯流的速度分量实施旋转的周期性约束。

在由 $\Omega=\{{\bf x}\in(0,1)^2:\|{\bf x}\|\in (0.5,1)\}$ 定义的四分之一圆上，我们要解决斯托克斯问题@f{eqnarray*}


  -\Delta \; \textbf{u} + \nabla p &=& (\exp(-100\|{\bf x}-(.75,0.1)^T\|^2),0)^T, \\


  -\textrm{div}\;  \textbf{u}&=&0,\\
  \textbf{u}|_{\Gamma_1}&=&{\bf 0},
@f} 

其中边界 $\Gamma_1$ 被定义为 $\Gamma_1 \dealcoloneq \{x\in \partial\Omega: \|x\|\in\{0.5,1\}\}$  。对于边界的其余部分，我们将使用周期性边界条件，即 

@f{align*}
  u_x(0,\nu)&=-u_y(\nu,0)&\nu&\in[0,1]\\
  u_y(0,\nu)&=u_x(\nu,0)&\nu&\in[0,1].
@f} 



网格将由 GridGenerator::quarter_hyper_shell(), 生成，它也记录了如果它的`colorize'参数设置为`true'，它是如何给它的各个边界分配边界指标的。<a name="CommProg"></a> <h1> The commented program</h1>

这个例子程序是对 step-22 的轻微修改，使用Trilinos并行运行，以演示交易.II中周期性边界条件的用法。因此我们不讨论大部分的源代码，只对处理周期性约束的部分进行评论。其余的请看 step-22 和底部的完整源代码。




为了实现周期性边界条件，只有两个函数需要修改。

-  <code>StokesProblem<dim>::setup_dofs()</code>  : 用周期性约束来填充AffineConstraints对象 

-  <code>StokesProblem<dim>::create_mesh()</code>  : 为分布式三角形提供周期性信息。




程序的其余部分与 step-22 相同，所以让我们跳过这一部分，只在下面展示这两个函数。不过完整的程序可以在下面的 "普通程序 "部分找到）。








   




<a name="Settingupperiodicityconstraintsondistributedtriangulations"></a> <h3>Setting up periodicity constraints on distributed triangulations</h3>

@code
  template <int dim>
  void StokesProblem<dim>::create_mesh()
  {
    Point<dim>   center;
    const double inner_radius = .5;
    const double outer_radius = 1.;


    GridGenerator::quarter_hyper_shell(
      triangulation, center, inner_radius, outer_radius, 0, true);


@endcode 



在我们可以规定周期性约束之前，我们需要确保位于域的对面但由周期性面连接的单元是幽灵层的一部分，如果其中一个单元存储在本地处理器上。在这一点上，我们需要考虑我们要如何规定周期性。左边边界上的面的顶点 $\text{vertices}_2$ 应该与下面边界上的面的顶点 $\text{vertices}_1$ 相匹配，由 $\text{vertices}_2=R\cdot
\text{vertices}_1+b$ 给出，其中旋转矩阵 $R$ 和偏移量 $b$ 给出 

@f{align*}
R=\begin{pmatrix}
0&1\\-1&0
\end{pmatrix},
\quad
b=\begin{pmatrix}0&0\end{pmatrix}.
@f} 

我们将所得信息保存到这里的数据结构是基于三角法的。

@code
    std::vector<GridTools::PeriodicFacePair<
      typename parallel::distributed::Triangulation<dim>::cell_iterator>>
      periodicity_vector;


    FullMatrix<double> rotation_matrix(dim);
    rotation_matrix[0][1] = 1.;
    rotation_matrix[1][0] = -1.;


    GridTools::collect_periodic_faces(triangulation,
                                      2,
                                      3,
                                      1,
                                      periodicity_vector,
                                      Tensor<1, dim>(),
                                      rotation_matrix);


@endcode 



现在，只要调用 parallel::distributed::Triangulation::add_periodicity. ，告诉三角法关于所需的周期性就特别容易。 

@code
    triangulation.add_periodicity(periodicity_vector);


    triangulation.refine_global(4 - dim);
  }



  template <int dim>
  void StokesProblem<dim>::setup_dofs()
  {
    dof_handler.distribute_dofs(fe);


    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);


    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];


    {
      owned_partitioning.clear();
      IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
      owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u));
      owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u + n_p));


      relevant_partitioning.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
      relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u));
      relevant_partitioning.push_back(
        locally_relevant_dofs.get_view(n_u, n_u + n_p));


      constraints.clear();
      constraints.reinit(locally_relevant_dofs);


      FEValuesExtractors::Vector velocities(0);


      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));


@endcode 



在我们为网格提供了周期性约束的必要信息后，我们现在就可以实际创建它们了。对于描述匹配，我们使用与之前相同的方法，也就是说，左边边界上的一个面的 $\text{vertices}_2$ 应该与下面边界上的一个面的顶点 $\text{vertices}_1$ 相匹配，这些顶点的旋转矩阵 $R$ 和偏移量 $b$ 是通过以下方式给出的 

@f{align*}
R=\begin{pmatrix}
0&1\\-1&0
\end{pmatrix},
\quad
b=\begin{pmatrix}0&0\end{pmatrix}.
@f} 

这两个对象不仅描述了面的匹配方式，而且还描述了解决方案应该从 $\text{face}_2$ 转换到 $\text{face}_1$ 的哪个意义上。

@code
      FullMatrix<double> rotation_matrix(dim);
      rotation_matrix[0][1] = 1.;
      rotation_matrix[1][0] = -1.;


      Tensor<1, dim> offset;


@endcode 



为了设置约束，我们首先将周期性信息存储在一个类型为 <code>std::vector@<GridTools::PeriodicFacePair<typename   DoFHandler@<dim@>::%cell_iterator@>  </code>的辅助对象中。周期性边界的边界指标为2（x=0）和3（y=0）。所有其他的参数我们之前已经设置好了。在这种情况下，方向并不重要。由于 $\text{vertices}_2=R\cdot \text{vertices}_1+b$ 这正是我们想要的。

@code
      std::vector<
        GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
        periodicity_vector;


      const unsigned int direction = 1;


      GridTools::collect_periodic_faces(dof_handler,
                                        2,
                                        3,
                                        direction,
                                        periodicity_vector,
                                        offset,
                                        rotation_matrix);


@endcode 



接下来，我们需要提供关于解决方案中哪些矢量值的分量应该被旋转的信息。由于我们在这里选择只约束速度，并且从解向量的第一个分量开始，我们简单地插入一个0。

@code
      std::vector<unsigned int> first_vector_components;
      first_vector_components.push_back(0);


@endcode 



在周期性_vector中设置了所有的信息后，我们要做的就是告诉make_periodicity_constraints来创建所需的约束。

@code
      DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector,
                                                       constraints,
                                                       fe.component_mask(
                                                         velocities),
                                                       first_vector_components);


      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
    }


    constraints.close();


    {
      TrilinosWrappers::BlockSparsityPattern bsp(owned_partitioning,
                                                 owned_partitioning,
                                                 relevant_partitioning,
                                                 mpi_communicator);


      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if (!((c == dim) && (d == dim)))
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;


      DoFTools::make_sparsity_pattern(dof_handler,
                                      coupling,
                                      bsp,
                                      constraints,
                                      false,
                                      Utilities::MPI::this_mpi_process(
                                        mpi_communicator));


      bsp.compress();


      system_matrix.reinit(bsp);
    }


    {
      TrilinosWrappers::BlockSparsityPattern preconditioner_bsp(
        owned_partitioning,
        owned_partitioning,
        relevant_partitioning,
        mpi_communicator);


      Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if ((c == dim) && (d == dim))
            preconditioner_coupling[c][d] = DoFTools::always;
          else
            preconditioner_coupling[c][d] = DoFTools::none;


      DoFTools::make_sparsity_pattern(dof_handler,
                                      preconditioner_coupling,
                                      preconditioner_bsp,
                                      constraints,
                                      false,
                                      Utilities::MPI::this_mpi_process(
                                        mpi_communicator));


      preconditioner_bsp.compress();


      preconditioner_matrix.reinit(preconditioner_bsp);
    }


    system_rhs.reinit(owned_partitioning, mpi_communicator);
    solution.reinit(owned_partitioning,
                    relevant_partitioning,
                    mpi_communicator);
  }


@endcode 



然后，程序的其余部分又与  step-22  相同。我们现在省略它，但和以前一样，你可以在下面的 "普通程序 "部分找到这些部分。




<a name="Results"></a><h1>Results</h1>


创建的输出并不十分令人惊讶。我们只是看到，相对于左边界和下边界而言，解是周期性的。

  <img src="https://www.dealii.org/images/steps/developer/step-45.periodic.png" alt="">   

如果没有周期性约束，我们最终会得到以下的解决方案。

  <img src="https://www.dealii.org/images/steps/developer/step-45.non_periodic.png" alt="">  <a name="PlainProg"></a> <h1> The plain program</h1>  @include "step-45.cc" 。 

  */  
