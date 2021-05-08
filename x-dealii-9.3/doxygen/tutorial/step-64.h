  /**   @page step_64 The step-64 tutorial program   

本教程取决于  step-7  ,  step-37  。

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The test case</a><a href="#Thetestcase">The test case</a>
        <li><a href="#Movingdatatoandfromthedevice">Moving data to and from the device</a><a href="#Movingdatatoandfromthedevice">Moving data to and from the device</a>
        <li><a href="#Matrixvectorproductimplementation">Matrix-vector product implementation</a><a href="#Matrixvectorproductimplementation">Matrix-vector product implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ClasscodeVaryingCoefficientFunctorcode">Class <code>VaryingCoefficientFunctor</code></a><a href="#ClasscodeVaryingCoefficientFunctorcode">Class <code>VaryingCoefficientFunctor</code></a>
        <li><a href="#ClasscodeHelmholtzOperatorQuadcode">Class <code>HelmholtzOperatorQuad</code></a><a href="#ClasscodeHelmholtzOperatorQuadcode">Class <code>HelmholtzOperatorQuad</code></a>
        <li><a href="#ClasscodeLocalHelmholtzOperatorcode">Class <code>LocalHelmholtzOperator</code></a><a href="#ClasscodeLocalHelmholtzOperatorcode">Class <code>LocalHelmholtzOperator</code></a>
        <li><a href="#ClasscodeHelmholtzOperatorcode">Class <code>HelmholtzOperator</code></a><a href="#ClasscodeHelmholtzOperatorcode">Class <code>HelmholtzOperator</code></a>
        <li><a href="#ClasscodeHelmholtzProblemcode">Class <code>HelmholtzProblem</code></a><a href="#ClasscodeHelmholtzProblemcode">Class <code>HelmholtzProblem</code></a>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly 

  <br>   

<i>
This program was contributed by Bruno Turcksin and Daniel Arndt, Oak Ridge National Laboratory.
</i> 


<a name="Introduction"></a><h1>Introduction</h1>


这个例子展示了如何使用CUDA在GPU上实现一个无矩阵的方法，用于超立方体上系数可变的亥姆霍兹方程。该线性系统将使用共轭梯度法进行求解，并使用MPI进行并行化。

在过去的几年里，异构计算，特别是GPU，已经获得了很大的普及。这是因为在给定的功率预算下，GPU比CPU提供更好的计算能力和内存带宽。在2019年初的架构中，对于PDE相关的任务，GPU的功率效率约为服务器CPU的2-3倍，宽<a
href="https://en.wikipedia.org/wiki/SIMD">SIMD</a>。GPU也是机器学习中最受欢迎的架构。另一方面，GPU并不容易编程。这个程序探索了deal.II的能力，看看这样的程序可以如何有效地实现。

虽然我们试图让CPU和GPU的无矩阵类的界面尽可能接近，但还是有一些区别。当在GPU上使用无矩阵框架时，人们必须编写一些CUDA代码。然而，数量相当少，而且对CUDA的使用仅限于几个关键词。


<a name="Thetestcase"></a><h3>The test case</h3>


在这个例子中，我们考虑Helmholtz问题@f{eqnarray*} - \nabla \cdot
\nabla u + a(\mathbf x) u &=&1,\\ u &=& 0 \quad \text{on } \partial \Omega @f}

其中 $a(\mathbf x)$ 是一个可变系数。

我们选择  $\Omega=[0,1]^3$  和  $a(\mathbf x)=\frac{10}{0.05 +
2\|\mathbf x\|^2}$  作为域。由于系数是围绕原点对称的，但域不是，我们最终会得到一个非对称的解决方案。

如果你在本教程中读到这里，你就会知道这个问题的弱表述是什么样子的，以及原则上如何为它组装线性系统。当然，在这个程序中，我们实际上不会形成矩阵，而只是表示它与之相乘时的作用。


<a name="Movingdatatoandfromthedevice"></a><h3>Moving data to and from the device</h3> 


GPU（我们从现在开始用 "设备 "一词来指代GPU）有自己的内存，与CPU（我们从现在开始用 "主机 "一词）可访问的内存分开。设备上的正常计算可以分为三个独立的步骤。

 -# 数据从主机移到设备上。

 -# 计算在设备上完成。

 -# 结果从设备上移回主机。

数据移动可以由用户代码明确完成，也可以使用UVM（统一虚拟内存）自动完成。在deal.II中，只支持第一种方法。虽然这意味着用户有额外的负担，但这可以更好地控制数据移动，更重要的是可以避免在主机上而不是设备上错误地运行重要的内核。

deal.II中的数据移动是通过 LinearAlgebra::ReadWriteVector. 完成的，这些向量可以被看作是主机上的缓冲区，用来存储从设备上收到的数据或发送数据到设备上。有两种类型的向量可以在设备上使用。

-  LinearAlgebra::CUDAWrappers::Vector,  它类似于更常见的Vector<Number>，和 

-  LinearAlgebra::distributed::Vector<Number,   MemorySpace::CUDA>,  这是一个普通的 LinearAlgebra::distributed::Vector ，我们指定了要使用哪个内存空间。

如果没有指定内存空间，默认是 MemorySpace::Host. 。 

接下来，我们展示如何使用 LinearAlgebra::CUDAWrappers::Vector: 将数据移入/移出设备。 

@code
  unsigned int size = 10;
  LinearAlgebra::ReadWriteVector<double> rw_vector(size);


  ...do something with the rw_vector...


  // Move the data to the device:
  LinearAlgebra::CUDAWrappers::Vector<double> vector_dev(size);
  vector_dev.import(rw_vector, VectorOperations::insert);


  ...do some computations on the device...


  // Move the data back to the host:
  rw_vector.import(vector_dev, VectorOperations::insert);
@endcode 

这里使用的两个向量类都只在一台机器上工作，也就是说，一个内存空间在主机上，一个在设备上。

但在有些情况下，人们希望在若干台机器上的多个MPI进程之间运行并行计算，而每台机器都配备有GPU。在这种情况下，人们希望使用 `LinearAlgebra::distributed::Vector<Number,MemorySpace::CUDA>`, ，它是类似的，但`import()`阶段可能涉及MPI通信。

@code
  IndexSet locally_owned_dofs, locally_relevant_dofs;
  ...fill the two IndexSet objects...


  // Create the ReadWriteVector using an IndexSet instead of the size
  LinearAlgebra::ReadWriteVector<double> owned_rw_vector(locally_owned_dofs);


  ...do something with the rw_vector...


  // Move the data to the device:
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
    distributed_vector_dev(locally_owned_dofs, MPI_COMM_WORLD);
  distributed_vector_dev.import(owned_rw_vector, VectorOperations::insert);


  ...do something with the dev_vector...


  // Create a ReadWriteVector with a different IndexSet:
  LinearAlgebra::ReadWriteVector<double>
    relevant_rw_vector(locally_relevant_dofs);


  // Move the data to the host, possibly using MPI communication:
  relevant_rw_vector.import(distributed_vector_dev, VectorOperations::insert);
@endcode 

`relevant_rw_vector`是一个存储向量所有元素的子集的对象。通常，这些是 @ref GlossLocallyRelevantDof "本地相关的DoF"，这意味着它们在不同的MPI进程之间是重叠的。因此，一台机器上存储在该向量中的元素可能与该机器上的GPU存储的元素不一致，需要MPI通信来导入它们。

在所有这些情况下，在导入一个向量时，可以插入数值（使用 VectorOperation::insert) 或添加到向量的先前内容中（使用 VectorOperation::add).   


<a name="Matrixvectorproductimplementation"></a><h3>Matrix-vector product implementation</h3>


在设备上评估无矩阵算子所需的代码与主机上的代码非常相似。然而，也有一些区别，主要是 Step-37 中的`local_apply()`函数和正交点的循环都需要封装在自己的函数中。<a name="CommProg"></a> <h1> The commented program</h1>

首先包括从以前的教程中知道的deal.II库中的必要文件。

@code
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>


#include <deal.II/dofs/dof_tools.h>


#include <deal.II/fe/fe_q.h>


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>


#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>


#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>


@endcode 



下面的包括在GPU上实现无矩阵方法的数据结构。

@code
#include <deal.II/base/cuda.h>


#include <deal.II/matrix_free/cuda_fe_evaluation.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>
#include <deal.II/matrix_free/operators.h>


#include <fstream>



@endcode 



像往常一样，我们把所有的东西都包围在一个自己的命名空间中。

@code
namespace Step64
{
  using namespace dealii;



@endcode 




<a name="ClasscodeVaryingCoefficientFunctorcode"></a> <h3>Class <code>VaryingCoefficientFunctor</code></h3>




接下来，我们定义一个类，实现我们想在亥姆霍兹算子中使用的变化系数。后来，我们想把这个类型的对象传递给一个 CUDAWrappers::MatrixFree 对象，该对象希望该类有一个`运算器()`，为给定的单元填充构造器中提供的值。这个操作符需要在设备上运行，所以它需要为编译器标记为`__device__`。

@code
  template <int dim, int fe_degree>
  class VaryingCoefficientFunctor
  {
  public:
    VaryingCoefficientFunctor(double *coefficient)
      : coef(coefficient)
    {}


    __device__ void operator()(
      const unsigned int                                          cell,
      const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data);


@endcode 



由于 CUDAWrappers::MatrixFree::Data 不知道其数组的大小，我们需要在这个类中存储正交点的数量和自由度的数量，以进行必要的索引转换。

@code
    static const unsigned int n_dofs_1d = fe_degree + 1;
    static const unsigned int n_local_dofs =
      dealii::Utilities::pow(n_dofs_1d, dim);
    static const unsigned int n_q_points =
      dealii::Utilities::pow(n_dofs_1d, dim);


  private:
    double *coef;
  };





@endcode 



下面的函数实现了这个系数。记得在介绍中，我们曾将其定义为 $a(\mathbf
x)=\frac{10}{0.05 + 2\|\mathbf x\|^2}$   

@code
  template <int dim, int fe_degree>
  __device__ void VaryingCoefficientFunctor<dim, fe_degree>::operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data)
  {
    const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
      cell, gpu_data, n_dofs_1d, n_q_points);
    const Point<dim> q_point =
      CUDAWrappers::get_quadrature_point<dim, double>(cell,
                                                      gpu_data,
                                                      n_dofs_1d);


    double p_square = 0.;
    for (unsigned int i = 0; i < dim; ++i)
      {
        const double coord = q_point[i];
        p_square += coord * coord;
      }
    coef[pos] = 10. / (0.05 + 2. * p_square);
  }



@endcode 




<a name="ClasscodeHelmholtzOperatorQuadcode"></a> <h3>Class <code>HelmholtzOperatorQuad</code></h3>




类`HelmholtzOperatorQuad`实现了Helmholtz算子在每个正交点的评估。它使用了类似于  step-37  中介绍的 MatrixFree 框架的机制。与那里不同的是，实际的正交点索引是通过转换当前线程索引隐式处理的。和以前一样，这个类的函数需要在设备上运行，所以需要为编译器标记为`__device__`。

@code
  template <int dim, int fe_degree>
  class HelmholtzOperatorQuad
  {
  public:
    __device__ HelmholtzOperatorQuad(double coef)
      : coef(coef)
    {}


    __device__ void
    operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const;


  private:
    double coef;
  };



@endcode 



我们在这里要解决的Helmholtz问题以弱的形式写成如下。@f{eqnarray*}
(\nabla v, \nabla u)+ (v, a(\mathbf x) u) &=&(v,1) \quad \forall v.
@f} 

如果你看过  step-37  ，那么很明显，左边的两个项对应着这里的两个函数调用。

@code
  template <int dim, int fe_degree>
  __device__ void HelmholtzOperatorQuad<dim, fe_degree>::
                  operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const
  {
    fe_eval->submit_value(coef * fe_eval->get_value());
    fe_eval->submit_gradient(fe_eval->get_gradient());
  }



@endcode 




<a name="ClasscodeLocalHelmholtzOperatorcode"></a> <h3>Class <code>LocalHelmholtzOperator</code></h3>




最后，我们需要定义一个类，实现整个运算符的评估，在基于矩阵的方法中对应于矩阵-向量积。

@code
  template <int dim, int fe_degree>
  class LocalHelmholtzOperator
  {
  public:
    LocalHelmholtzOperator(double *coefficient)
      : coef(coefficient)
    {}


    __device__ void operator()(
      const unsigned int                                          cell,
      const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
      CUDAWrappers::SharedData<dim, double> *                     shared_data,
      const double *                                              src,
      double *                                                    dst) const;


@endcode 



同样， CUDAWrappers::MatrixFree 对象不知道自由度的数量和正交点的数量，所以我们需要在调用算子中存储这些用于索引计算。

@code
    static const unsigned int n_dofs_1d    = fe_degree + 1;
    static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
    static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);


  private:
    double *coef;
  };



@endcode 



这是一个调用算子，在给定的单元上执行亥姆霍兹算子的评估，类似于CPU上的MatrixFree框架。特别是，我们需要访问源向量的值和梯度，我们将值和梯度信息写入目标向量。

@code
  template <int dim, int fe_degree>
  __device__ void LocalHelmholtzOperator<dim, fe_degree>::operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
    CUDAWrappers::SharedData<dim, double> *                     shared_data,
    const double *                                              src,
    double *                                                    dst) const
  {
    const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
      cell, gpu_data, n_dofs_1d, n_q_points);


    CUDAWrappers::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double>
      fe_eval(cell, gpu_data, shared_data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(true, true);
    fe_eval.apply_for_each_quad_point(
      HelmholtzOperatorQuad<dim, fe_degree>(coef[pos]));
    fe_eval.integrate(true, true);
    fe_eval.distribute_local_to_global(dst);
  }



@endcode 




<a name="ClasscodeHelmholtzOperatorcode"></a> <h3>Class <code>HelmholtzOperator</code></h3> 




HelmholtzOperator "类作为 "LocalHelmholtzOperator "的封装器，定义了一个可以与SolverCG等线性求解器一起使用的接口。特别是，像每一个实现线性算子接口的类一样，它需要有一个`vmult()'函数来执行线性算子对源向量的操作。

@code
  template <int dim, int fe_degree>
  class HelmholtzOperator
  {
  public:
    HelmholtzOperator(const DoFHandler<dim> &          dof_handler,
                      const AffineConstraints<double> &constraints);


    void
    vmult(LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &dst,
          const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
            &src) const;


    void initialize_dof_vector(
      LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const;


  private:
    CUDAWrappers::MatrixFree<dim, double>       mf_data;
    LinearAlgebra::CUDAWrappers::Vector<double> coef;
  };





@endcode 



下面是这个类的构造函数的实现。在第一部分，我们初始化`mf_data`成员变量，该变量将在评估算子时为我们提供必要的信息。   


在第二部分中，我们需要在每个活动的、本地拥有的单元中存储每个正交点的系数值。我们可以向平行三角法询问活动的、本地拥有的单元的数量，但手头只有一个DoFHandler对象。由于 DoFHandler::get_triangulation() 返回的是Triangulation对象，而不是 parallel::TriangulationBase 对象，我们必须对返回值进行下移。在这里这样做是安全的，因为我们知道三角形实际上是一个 parallel:distributed::Triangulation 对象。

@code
  template <int dim, int fe_degree>
  HelmholtzOperator<dim, fe_degree>::HelmholtzOperator(
    const DoFHandler<dim> &          dof_handler,
    const AffineConstraints<double> &constraints)
  {
    MappingQGeneric<dim> mapping(fe_degree);
    typename CUDAWrappers::MatrixFree<dim, double>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags = update_values | update_gradients |
                                           update_JxW_values |
                                           update_quadrature_points;
    const QGauss<1> quad(fe_degree + 1);
    mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);



    const unsigned int n_owned_cells =
      dynamic_cast<const parallel::TriangulationBase<dim> *>(
        &dof_handler.get_triangulation())


        ->n_locally_owned_active_cells();
    coef.reinit(Utilities::pow(fe_degree + 1, dim) * n_owned_cells);


    const VaryingCoefficientFunctor<dim, fe_degree> functor(coef.get_values());
    mf_data.evaluate_coefficients(functor);
  }



@endcode 



然后，关键的一步是使用前面所有的类来循环所有的单元格来执行矩阵-向量乘积。我们在下一个函数中实现这一点。   


在应用亥姆霍兹算子时，我们必须注意正确处理边界条件。因为本地算子不知道约束条件，所以我们必须在事后将正确的值从源头复制到目的向量上。

@code
  template <int dim, int fe_degree>
  void HelmholtzOperator<dim, fe_degree>::vmult(
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &      dst,
    const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &src)
    const
  {
    dst = 0.;
    LocalHelmholtzOperator<dim, fe_degree> helmholtz_operator(
      coef.get_values());
    mf_data.cell_loop(helmholtz_operator, src, dst);
    mf_data.copy_constrained_values(src, dst);
  }





  template <int dim, int fe_degree>
  void HelmholtzOperator<dim, fe_degree>::initialize_dof_vector(
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const
  {
    mf_data.initialize_dof_vector(vec);
  }



@endcode 




<a name="ClasscodeHelmholtzProblemcode"></a><h3>Class <code>HelmholtzProblem</code></h3>




这是本程序的主类。它定义了我们用于教程程序的通常框架。唯一值得评论的一点是`solve()`函数和向量类型的选择。

@code
  template <int dim, int fe_degree>
  class HelmholtzProblem
  {
  public:
    HelmholtzProblem();


    void run();


  private:
    void setup_system();


    void assemble_rhs();


    void solve();


    void output_results(const unsigned int cycle) const;


    MPI_Comm mpi_communicator;


    parallel::distributed::Triangulation<dim> triangulation;


    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;


    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;


    AffineConstraints<double>                          constraints;
    std::unique_ptr<HelmholtzOperator<dim, fe_degree>> system_matrix_dev;


@endcode 



由于`solve()`函数中的所有操作都是在显卡上执行的，因此所使用的向量也有必要在GPU上存储其值。  LinearAlgebra::distributed::Vector 可以被告知要使用哪个内存空间。还有 LinearAlgebra::CUDAWrappers::Vector ，总是使用GPU内存存储，但不与MPI一起工作。值得注意的是，如果MPI实现是CUDA感知的，并且配置标志`DEAL_II_MPI_WITH_CUDA_SUPPORT`被启用，不同MPI进程之间的通信可以得到改善。(这个标志的值需要在安装deal.II时调用`cmake`时设置)。     


此外，我们还保留了一个带有CPU存储的解决方案向量，这样我们就可以像往常一样查看和显示解决方案。

@code
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
                                                                  ghost_solution_host;
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> solution_dev;
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
      system_rhs_dev;


    ConditionalOStream pcout;
  };



@endcode 



除了 `Helmholtzproblem::solve()` 之外，这个类的所有其余函数的实现并不包含任何新的内容，我们不会对整体方法进一步发表过多的评论。

@code
  template <int dim, int fe_degree>
  HelmholtzProblem<dim, fe_degree>::HelmholtzProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator)
    , fe(fe_degree)
    , dof_handler(triangulation)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {}





  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::setup_system()
  {
    dof_handler.distribute_dofs(fe);


    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    system_rhs_dev.reinit(locally_owned_dofs, mpi_communicator);


    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();


    system_matrix_dev.reset(
      new HelmholtzOperator<dim, fe_degree>(dof_handler, constraints));


    ghost_solution_host.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               mpi_communicator);
    system_matrix_dev->initialize_dof_vector(solution_dev);
    system_rhs_dev.reinit(solution_dev);
  }





@endcode 



与 step-4 或 step-6 等程序不同，我们不必组装整个线性系统，只需组装右手边的向量。这在本质上与我们在 step-4 中所做的一样，例如，但我们必须注意在将局部贡献复制到全局矢量时使用正确的约束对象。特别是，我们需要确保对应于边界节点的条目被正确地清零了。这对于CG的收敛是必要的。 另一个解决方案是修改上面的`vmult()`函数，我们假装源向量的条目为零，在矩阵-向量乘积中不考虑它们。但这里使用的方法更简单）。)    


在函数的最后，我们不能直接将数值从主机复制到设备上，而是需要使用一个类型为 LinearAlgebra::ReadWriteVector 的中间对象来构建必要的正确通信模式。

@code
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::assemble_rhs()
  {
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
                      system_rhs_host(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
    const QGauss<dim> quadrature_formula(fe_degree + 1);


    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);


    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();


    Vector<double> cell_rhs(dofs_per_cell);


    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_rhs = 0;


          fe_values.reinit(cell);


          for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1.0 *
                                fe_values.JxW(q_index));
            }


          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 system_rhs_host);
        }
    system_rhs_host.compress(VectorOperation::add);


    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import(system_rhs_host, VectorOperation::insert);
    system_rhs_dev.import(rw_vector, VectorOperation::insert);
  }





@endcode 



这个solve()函数最后包含了对之前讨论的新类的调用。这里我们不使用任何预处理程序，即通过身份矩阵进行预处理，以只关注 CUDAWrappers::MatrixFree 框架的特殊性。当然，在实际应用中，选择一个合适的预处理程序是至关重要的，但我们至少有与 step-37 中相同的限制，因为矩阵条目是即时计算而不是存储的。   


在函数的第一部分解出线性系统后，我们将解决方案从设备上复制到主机上，以便能够查看其值并在`output_results()`中显示。这种转移的工作方式与前一个函数的结尾相同。

@code
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::solve()
  {
    PreconditionIdentity preconditioner;


    SolverControl solver_control(system_rhs_dev.size(),
                                 1e-12 * system_rhs_dev.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>> cg(
      solver_control);
    cg.solve(*system_matrix_dev, solution_dev, system_rhs_dev, preconditioner);


    pcout << "  Solved in " << solver_control.last_step() << " iterations."
          << std::endl;


    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import(solution_dev, VectorOperation::insert);
    ghost_solution_host.import(rw_vector, VectorOperation::insert);


    constraints.distribute(ghost_solution_host);


    ghost_solution_host.update_ghost_values();
  }


@endcode 



输出结果函数和平时一样，因为我们已经将数值从GPU复制回CPU。   


既然我们已经在用这个函数做事情了，我们不妨计算一下解决方案的 $L_2$ 规范。我们通过调用 VectorTools::integrate_difference(). 来做到这一点，该函数旨在通过评估数值解（由自由度值的向量给出）和代表精确解的对象之间的差异来计算误差。但是我们可以通过传入一个零函数来轻松地计算解决方案的 $L_2$ 准则。也就是说，我们不是在评估误差 $\|u_h-u\|_{L_2(\Omega)}$ ，而是在评估 $\|u_h-0\|_{L_2(\Omega)}=\|u_h\|_{L_2(\Omega)}$ 。

@code
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::output_results(
    const unsigned int cycle) const
  {
    DataOut<dim> data_out;


    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(ghost_solution_host, "solution");
    data_out.build_patches();


    DataOutBase::VtkFlags flags;
    flags.compression_level = DataOutBase::VtkFlags::best_speed;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2);


    Vector<float> cellwise_norm(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      ghost_solution_host,
                                      Functions::ZeroFunction<dim>(),
                                      cellwise_norm,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm);
    const double global_norm =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_norm,
                                        VectorTools::L2_norm);
    pcout << "  solution norm: " << global_norm << std::endl;
  }



@endcode 



`run()`函数中也没有什么令人惊讶的地方。我们只是在一系列（全局）细化的网格上计算解决方案。

@code
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::run()
  {
    for (unsigned int cycle = 0; cycle < 7 - dim; ++cycle)
      {
        pcout << "Cycle " << cycle << std::endl;


        if (cycle == 0)
          GridGenerator::hyper_cube(triangulation, 0., 1.);
        triangulation.refine_global(1);


        setup_system();


        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;


        assemble_rhs();
        solve();
        output_results(cycle);
        pcout << std::endl;
      }
  }
} // namespace Step64



@endcode 




<a name="Thecodemaincodefunction"></a><h3>The <code>main()</code> function</h3>




最后为`main()`函数。 默认情况下，所有的MPI等级将尝试访问编号为0的设备，我们假设它是与某一MPI等级运行的CPU相关的GPU设备。这是可行的，但是如果我们在运行MPI支持时，可能会有多个MPI进程在同一台机器上运行（例如，每个CPU核心一个），然后它们都想访问该机器上的同一个GPU。如果机器上只有一个GPU，我们对此无能为力。该机器上的所有MPI行列都需要共享它。但是如果有不止一个GPU，那么最好为不同的进程解决不同的显卡。下面的选择是基于MPI进程的ID，通过将GPU轮流分配给GPU行列。(为了正确地工作，这个方案假定一台机器上的MPI等级是连续的。如果不是这样的话，那么等级与GPU的关联可能就不是最佳的了）。) 为了使其正常工作，在使用这个函数之前，需要对MPI进行初始化。

@code
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step64;


      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);


      int         n_devices       = 0;
      cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
      AssertCuda(cuda_error_code);
      const unsigned int my_mpi_id =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      const int device_id = my_mpi_id % n_devices;
      cuda_error_code     = cudaSetDevice(device_id);
      AssertCuda(cuda_error_code);


      HelmholtzProblem<3, 3> helmholtz_problem;
      helmholtz_problem.run();
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


由于本教程的主要目的是演示如何使用 CUDAWrappers::MatrixFree 接口，而不是计算任何有用的东西本身，我们只是在这里显示预期的输出。

@code
Cycle 0
   Number of active cells:       8
   Number of degrees of freedom: 343
  Solved in 27 iterations.
  solution norm: 0.0205439


Cycle 1
   Number of active cells:       64
   Number of degrees of freedom: 2197
  Solved in 60 iterations.
  solution norm: 0.0205269


Cycle 2
   Number of active cells:       512
   Number of degrees of freedom: 15625
  Solved in 114 iterations.
  solution norm: 0.0205261


Cycle 3
   Number of active cells:       4096
   Number of degrees of freedom: 117649
  Solved in 227 iterations.
  solution norm: 0.0205261
@endcode 



在这里，人们可以提出两个看法。首先，数值解的规范收敛了，大概是收敛到精确（但未知）解的规范。其次，每次细化网格时，迭代次数大约增加一倍。这与CG迭代次数随矩阵条件数的平方根增长的预期一致；而且我们知道二阶微分运算的矩阵条件数的增长方式为 ${\cal O}(h^{-2})$  。这当然是相当低效的，因为一个最佳解算器的迭代次数与问题的大小无关。但是要有这样一个求解器，就需要使用比我们在这里使用的身份矩阵更好的预处理程序。


<a name="extensions"></a> <a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


目前，这个程序完全没有使用预处理程序。这主要是因为构造一个高效的无矩阵预处理程序是不容易的。 然而，只需要相应矩阵的对角线的简单选择是很好的选择，这些也可以用无矩阵的方式计算。另外，也许更好的是，我们可以扩展教程，使用类似于  step-37  的切比雪夫平滑器的多网格。<a name="PlainProg"></a> <h1> The plain program</h1>  @include "step-64.cu"  。 

  */  
