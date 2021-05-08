  /**   @page step_54 The step-54 tutorial program 。 

本教程取决于  step-53  。

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#CADsurfaces"> CAD surfaces </a><a href="#CADsurfaces"> CAD surfaces </a>
        <li><a href="#TheCADboundaryprojectorclasses"> The CAD boundary projector classes </a><a href="#TheCADboundaryprojectorclasses"> The CAD boundary projector classes </a>
        <li><a href="#Thetestcase"> The testcase </a><a href="#Thetestcase"> The testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#TheTriangulationOnCADclass">The TriangulationOnCAD class</a><a href="#TheTriangulationOnCADclass">The TriangulationOnCAD class</a>
      <ul>
        <li><a href="#TriangulationOnCADTriangulationOnCAD">TriangulationOnCAD::TriangulationOnCAD</a><a href="#TriangulationOnCADTriangulationOnCAD">TriangulationOnCAD::TriangulationOnCAD</a>
        <li><a href="#TriangulationOnCADread_domain">TriangulationOnCAD::read_domain</a><a href="#TriangulationOnCADread_domain">TriangulationOnCAD::read_domain</a>
        <li><a href="#TriangulationOnCADrefine_mesh">TriangulationOnCAD::refine_mesh</a><a href="#TriangulationOnCADrefine_mesh">TriangulationOnCAD::refine_mesh</a>
        <li><a href="#TriangulationOnCADoutput_results">TriangulationOnCAD::output_results</a><a href="#TriangulationOnCADoutput_results">TriangulationOnCAD::output_results</a>
        <li><a href="#TriangulationOnCADrun">TriangulationOnCAD::run</a><a href="#TriangulationOnCADrun">TriangulationOnCAD::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly 

  <br>   

<i>This program was contributed by Andrea Mola and Luca Heltai.</i> 

  @note  这个程序阐述了工业几何的概念，使用与OpenCASCADE库（http://www.opencascade.org）接口的工具，允许指定任意的IGES文件来描述你的几何图形的边界。

  @dealiiTutorialDOI{10.5281/zenodo.546220,https://zenodo.org/badge/DOI/10.5281/zenodo.546220.svg}   

<a name="Intro"></a> <a name="Introduction"></a><h1>Introduction</h1> 




在之前的一些教程中（ step-1 ,  step-3 ,  step-5 ,  step-6 和 step-49 等），我们已经学会了如何使用 deal.II 中提供的网格细化方法。这些教程展示了如何使用这些工具为一次模拟产生一个精细的网格，如  step-3  ；或者从一个粗大的网格开始，在自适应细化的网格上进行一系列的模拟，如  step-6  的情况。无论采取哪种方法，网格细化都需要对计算域边界进行适当的几何描述，以便在每次细化时将新的网格节点置于边界表面。例如， step-5 展示了如何创建一个圆形网格，将一个圆形流形对象自动附加到计算域上，从而使位于边界上的面被细化到圆上。  step-53 展示了如何用一个由实验获得的数据定义的流形做这件事。但是，至少就基本边界形状而言，deal.II实际上只提供了圆、球、盒和其他基本组合。在本教程中，我们将展示如何使用一套开发的类来导入任意的CAD几何图形，将它们分配到计算域的所需边界，并在这种复杂的形状上细化计算网格。


<a name="CADsurfaces"></a><h3> CAD surfaces </h3> 


在最常见的工业实践中，任意形状的物体的几何模型是通过计算机辅助设计（CAD）工具实现的。在过去的几十年里，CAD建模器的使用得到了推广，因为它们可以为每个设计对象生成一个完整的虚拟模型，通过计算机可以在实物制作之前对其最精细的细节进行可视化、检查和分析。 从数学的角度来看，CAD建模人员的引擎是由分析几何学来代表的，特别是由参数化的曲线和曲面，如B-splines和NURBS，它们足够丰富，可以代表大多数实际利益的表面。 一旦一个虚拟模型准备好了，所需物体的所有几何特征都被存储在文件中，这些文件实质上包含了构成该物体的参数化曲面和曲线的系数。根据用于定义几何模型的具体CAD工具，当然有几种不同的文件格式，可以组织CAD模型的信息。为了提供一个跨CAD工具交换数据的共同基础，美国国家标准局在1980年发布了初始图形交换表示法（IGES）中性文件格式，在本例中使用。

<a name="TheCADboundaryprojectorclasses"></a><h3> The CAD boundary projector classes </h3>


为了导入和查询CAD模型，deal.II库为CAD建模的OpenCASCADE开源库实现了一系列封装函数。这些函数允许将IGES文件导入OpenCASCADE本地对象，并将其包裹在一系列Manifold类中。

一旦从IGES文件导入，模型就被存储在 <code>TopoDS_Shape</code> 中，这是OpenCASCADE框架中定义的通用拓扑实体。从 <code>TopoDS_Shape</code> 中，就可以访问构成它的所有子形状（如顶点、边和面）以及它们的几何描述。在deal.II框架中，组成一个形状的拓扑实体被用来创建一个相应的Mifold表示。在 step-6 中，我们看到如何使用 GridGenerator::hyper_sphere() 来创建一个超球体，它自动地将一个SphericalManifold附加到所有的边界面。这保证了边界面在网格细化过程中保持在球体或圆上。CAD建模界面的功能被设计为保留相同的结构，允许用户使用导入的CAD形状建立一个投影仪对象，保持我们在其他教程程序中使用的相同程序，即把这种投影仪对象分配给粗略网格的单元、面或边。在每个细化周期，新的网格节点将通过将现有对象的中点投影到指定的几何体上而自动生成。

与球形或圆形边界不同的是，具有复杂几何形状的边界带来的问题是，将细化后产生的新节点放在规定形状上最合适。例如，PolarManifold将周围的点转换为极坐标，计算该坐标系中的平均值（对每个坐标单独计算），最后将点转换回直角坐标。

不过，在一个任意的复杂形状的情况下，一个适当的选择来放置一个新的节点是不可能那么容易确定的。deal.II中的OpenCASCADE封装器提供了几个采用不同投影策略的投影仪类。第一个投影仪，在 OpenCASCADE::ArclengthProjectionLineManifold 类中实现，只用于边缘细化。它的建立是给它分配一个维度为1的拓扑形状，或者是一个 <code>TopoDS_Edge</code> or a <code>TopoDS_Wire</code> （这是一个复合形状，由几个相连的 <code>TopoDS_Edge</code> 组成），并细化一个网格边缘，找到新的顶点作为点，将位于原始边缘顶点之间的CAD曲线部分的曲线长度分成两个均匀部分。

  <img src="https://www.dealii.org/images/steps/developer/step-54.CurveSplit.png" alt="" width="500">   


在 OpenCASCADE::NormalProjectionBoundary 类中实现了一个不同的投影策略。在构造时分配的 <code>TopoDS_Shape</code> 可以是任意的（图形、面、边的集合或单个面或边都可以）。新的单元格节点首先通过对周围的点进行平均计算，方法与FlatManifold相同。在第二步中，所有的新节点将沿着形状的法线方向被投射到 <code>TopoDS_Shape</code> 。如果没有法线投影，则选择最接近形状的点--通常位于形状的边界上--。 如果形状是由几个子形状组成的，则投影到每个子形状上，并选择最近的投影点。

  <img src="https://www.dealii.org/images/steps/developer/step-54.NormalProjectionEdge.png" alt="" width="500">   <img src="https://www.dealii.org/images/steps/developer/step-54.NormalProjection.png" alt="" width="500"> 。 

正如我们即将体验到的，对于某些形状，将投影方向设置为CAD表面的法线，将不会导致合适质量的表面网格元素。这是因为CAD表面的法线方向原则上与网格需要新节点的方向没有关系。在这种情况下， OpenCASCADE::DirectionalProjectionBoundary 类可以提供帮助。这个类的构造是指定一个 <code>TopoDS_Shape</code> （至少包含一个面）和一个方向，所有的投影将沿着这个方向进行。新的点将通过首先对周围的点进行平均计算（就像FlatManifold的情况一样），然后沿着构建时使用的方向，在拓扑形状和通过所产生的点的线之间取最接近的交点。 通过这种方式，用户将对要执行的投影方向有更高的控制，以确保良好的网格质量。

  <img src="https://www.dealii.org/images/steps/developer/step-54.DirectionalProjection.png" alt="" width="500">   


当然，后一种方法只有在曲面的方向比较统一的情况下才有效，这样就可以确定一个单一的投影方向。在曲面方向接近于投影方向的情况下，甚至有可能找不到方向性的投影。为了克服这些问题， OpenCASCADE::NormalToMeshProjectionBoundary 类实现了第三个投影算法。 OpenCASCADE::NormalToMeshProjectionBoundary 类是将一个 <code>TopoDS_Shape</code> （至少包含一个面）分配给构造函数建立的，其工作原理与 OpenCASCADE::DirectionalProjection. 完全一样。但是，正如该类的名字所暗示的， OpenCASCADE::NormalToMeshProjectionBoundary 试图想出一个合适的对要精化的网格元素的法线方向的估计，并将其用于新节点在CAD曲面上的投影。如果我们考虑二维空间中的网格边缘，其轴线方向是一个方向，沿着这个方向分割，以产生两个相同长度的新单元。我们在此将这一概念扩展到三维空间，并将所有新的节点投射到一个近似于单元格法线的方向。

在下图中，受本教程中考虑的几何形状的启发，我们试图比较所考虑的三种投影仪的行为。从左边可以看出，给定原始单元（蓝色），用法线投影找到的新点的位置不允许生成均匀的新元素（红色）。这种情况在进一步的细化步骤中会变得更糟。 由于我们考虑的几何体在某种程度上垂直于水平方向，以水平方向为投影方向的方向性投影（中心图像）在获得新的网格点方面做得相当好。然而，由于图片底部的表面几乎是水平的，我们可以预期在这些区域进行进一步细化步骤时，会出现问题。最后，右边的图片显示，位于单元轴上的节点将导致两个新单元具有相同的长度。当然，三维的情况会比这个简单的二维案例中描述的情况更复杂一些。然而，这个测试的结果证实，当考虑到任意形状的表面时，除非你有一个几何形状，已知有一个更具体的方法是合适的，否则法线方向是测试的三种方法中最好的。


  <img src="https://www.dealii.org/images/steps/developer/step-54.ProjectionComparisons.png" alt="" width="700">   


<a name="Thetestcase"></a><h3> The testcase </h3>


在这个程序中，我们将考虑为一个描述船头的真实几何体创建一个曲面网格（这个几何体经常被用于CAD和网格生成的比较中，并且可以免费获得）。我们得到的曲面网格可以用来求解边界元素方程，以模拟水在船舶周围的流动（类似于 step-34 的方式），但我们在这里不会尝试这样做。为了让你了解我们所考虑的几何形状，这里有一张图片。

  <img src="https://www.dealii.org/images/steps/developer/step-54.bare.png" alt="" width="500">   

在这个程序中，我们从文件中读取几何体和粗略的网格，然后使用上面讨论的几个选项来放置新的顶点，进行一系列的网格细化步骤。<a name="CommProg"></a> <h1> The commented program</h1>


<a name="Includefiles"></a> <h3>Include files</h3> 




我们首先包括一堆文件，我们将在程序的各个部分使用这些文件。它们中的大多数已经在以前的教程中讨论过了。

@code
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>


@endcode 



这些是opencascade支持类和函数的头文件。注意，只有当你在编译deal.II库时支持OpenCASCADE，即在deal.II配置过程中调用 <code>-DDEAL_II_WITH_OPENCASCADE=ON</code> 和 <code>-DOPENCASCADE_DIR=/path/to/your/opencascade/installation</code> 时，这些才会包含合理的数据。

@code
#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>



@endcode 



最后，几个C++标准头文件 

@code
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


@endcode 



我们将程序的其余部分隔离在它自己的命名空间中 

@code
namespace Step54
{
  using namespace dealii;





@endcode 




<a name="TheTriangulationOnCADclass"></a> <h3>The TriangulationOnCAD class</h3> 




这是主类。它真正做的是存储输入和输出文件的名称，以及一个三角图。然后，它提供了一个函数，可以从粗略的网格中生成这样一个三角形，使用介绍中讨论的策略之一，并列在类的顶部的枚举类型中。   


这个类的成员函数类似于你可以在其他大多数教程程序中找到的模拟用网格的设置阶段的内容。







@code
  class TriangulationOnCAD
  {
  public:
    enum ProjectionType
    {
      NormalProjection       = 0,
      DirectionalProjection  = 1,
      NormalToMeshProjection = 2
    };



    TriangulationOnCAD(
      const std::string &  initial_mesh_filename,
      const std::string &  cad_file_name,
      const std::string &  output_filename,
      const ProjectionType surface_projection_kind = NormalProjection);


    void run();


  private:
    void read_domain();


    void refine_mesh();


    void output_results(const unsigned int cycle);


    Triangulation<2, 3> tria;


    const std::string initial_mesh_filename;
    const std::string cad_file_name;
    const std::string output_filename;


    const ProjectionType surface_projection_kind;
  };



@endcode 




<a name="TriangulationOnCADTriangulationOnCAD"></a> <h4>TriangulationOnCAD::TriangulationOnCAD</h4> 




TriangulationOnCAD类的构造函数非常简单。输入参数是输入和输出文件名的字符串，以及决定在网格细化循环中使用哪种曲面投影仪的枚举类型（详见下文）。







@code
  TriangulationOnCAD::TriangulationOnCAD(
    const std::string &  initial_mesh_filename,
    const std::string &  cad_file_name,
    const std::string &  output_filename,
    const ProjectionType surface_projection_kind)
    : initial_mesh_filename(initial_mesh_filename)
    , cad_file_name(cad_file_name)
    , output_filename(output_filename)
    , surface_projection_kind(surface_projection_kind)
  {}



@endcode 




<a name="TriangulationOnCADread_domain"></a> <h4>TriangulationOnCAD::read_domain</h4> 








下面的函数代表了这个程序的核心。 在这个函数中，我们导入CAD形状，在此基础上生成并细化我们的三角剖面。我们假设CAD表面包含在 @p cad_file_name 文件中（我们在输入目录中提供了一个名为 "input/DTMB-5415_bulbous_bow.iges "的IGES文件的例子，它代表了一艘船的球状船头）。由于存在几个凸和凹的高曲率区域，我们提供的几何体是一个特别有意义的例子。   


在导入船首表面后，我们提取了组成它的一些曲线和曲面，并使用它们来生成一组投影仪。这些投影仪定义了三角法在细化单元时必须遵循的规则，以定位每个新节点。   


为了初始化Triangulation，就像以前的教程程序一样，我们导入一个以VTK格式保存的已有网格。在这里我们假设用户已经在外部生成了一个粗略的网格，与IGES的几何图形相匹配。在编写本教程的时候，deal.II库并不自动支持生成这样的网格，但是有一些工具可以从CAD文件开始为你提供合理的初始网格。在我们的例子中，导入的网格是由一个四边形单元组成的，其顶点被放置在CAD形状上。   


在导入IGES几何体和初始网格后，我们将之前讨论过的投影仪分配给每个需要在CAD表面进行细化的边和单元。   


在本教程中，我们将测试介绍中所描述的三种不同的CAD表面投影器，并将分析每一种投影器所得到的结果。 如前所述，这些投影策略中的每一种都是在不同的类中实现的，这些类型的对象可以用 Triangulation::set_manifold 的方法分配给一个三角形。   


然后，下面的函数首先导入给定的CAD文件。函数的参数是一个包含所需文件名的字符串，以及一个比例因子。在这个例子中，比例因子被设置为1e-3，因为原始几何体是以毫米为单位的（这是大多数IGES文件的典型计量单位），而我们更喜欢以米为单位工作。 该函数的输出是一个OpenCASCADE通用拓扑形状类的对象，即一个 @p TopoDS_Shape. 。 

@code
  void TriangulationOnCAD::read_domain()
  {
    TopoDS_Shape bow_surface = OpenCASCADE::read_IGES(cad_file_name, 1e-3);


@endcode 



每个CAD几何对象都被定义了一个公差，这表明其位置可能不准确。例如，顶点的公差 @p tol 表示它可以位于以标称位置为中心、半径为 @p tol. 的球体中的任何一点。 当把一个点投影到一个表面（这反过来也有其公差）时，我们必须记住，投影的精度将受到建立表面时的公差的限制。




下面的方法是提取给定形状的公差，并使其稍大一些，以避免麻烦。

@code
    const double tolerance = OpenCASCADE::get_shape_tolerance(bow_surface) * 5;


@endcode 



我们现在要从通用形状中提取一组复合子形状。特别是，CAD文件的每个面是由类型为 @p TopoDS_Wire, 的修剪曲线组成的，它是构成曲面边界的 @p TopoDS_Edges 的集合，以及对曲面本身的NURBS描述。我们将使用一个线条投影仪来将我们的三角形的边界与划定曲面的线联系起来。 为了提取所有的复合子形状，如线、壳或实体，我们求助于OpenCASCADE命名空间的一种方法。  OpenCASCADE::extract_compound_shapes 的输入是一个形状和一组空的 std::vectors 子形状，它将被填入在给定拓扑形状中发现的所有复合形状。

@code
    std::vector<TopoDS_Compound>  compounds;
    std::vector<TopoDS_CompSolid> compsolids;
    std::vector<TopoDS_Solid>     solids;
    std::vector<TopoDS_Shell>     shells;
    std::vector<TopoDS_Wire>      wires;


    OpenCASCADE::extract_compound_shapes(
      bow_surface, compounds, compsolids, solids, shells, wires);


@endcode 



接下来的几个步骤比较熟悉，允许我们从外部VTK文件导入一个现有的网格，并将其转换为交易三角。

@code
    std::ifstream in;


    in.open(initial_mesh_filename);


    GridIn<2, 3> gi;
    gi.attach_triangulation(tria);
    gi.read_vtk(in);


@endcode 



我们输出这个初始网格，将其保存为细化步骤0。

@code
    output_results(0);


@endcode 



导入的网格有一个位于三维空间的单一二维单元。我们现在要确保它是根据上面导入的CAD几何图形进行细化的。为此，我们得到一个单元的迭代器，并给它分配manifold_id 1（见 @ref GlossManifoldIndicator  "这个词汇表条目"）。我们还得到了一个指向其四个面的迭代器，并为每个面分配了manifold_id 2。

@code
    Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
    cell->set_manifold_id(1);


    for (const auto &face : cell->face_iterators())
      face->set_manifold_id(2);


@endcode 



一旦CAD几何体和初始网格都被导入和消化，我们就用CAD的面和曲线来定义投影仪，并将它们分配给刚才指定的流形ID。




第一个投影器是用我们的CAD文件中的单线定义的。 ArclengthProjectionLineManifold将确保位于导线上的每条网格边缘都被细化为一个位于导线上的点，并将其分割为两个位于边缘顶点之间的相等弧线。我们首先检查导线矢量是否至少包含一个元素，然后为其创建一个Manifold对象。     


一旦创建了投影仪，我们就把它分配给三角形的所有部分，manifold_id=2。

@code
    Assert(
      wires.size() > 0,
      ExcMessage(
        "I could not find any wire in the CAD file you gave me. Bailing out."));


    OpenCASCADE::ArclengthProjectionLineManifold<2, 3> line_projector(
      wires[0], tolerance);


    tria.set_manifold(2, line_projector);


@endcode 



根据构造函数的 @p surface_projection_kind 选项所指定的内容来创建曲面投影仪。特别是，如果surface_projection_kind值等于 @p NormalProjection, ，我们选择 OpenCASCADE::NormalProjectionManifold. ，然后新的网格点将在所考虑的单元格/边的arycenter处最初生成，然后沿其法线方向投影到CAD表面。 NormalProjectionManifold构造器只需要一个形状和一个公差，然后我们把它分配给三角结构，用于所有具有id 1的流形的零件。

@code
    switch (surface_projection_kind)
      {
        case NormalProjection:
          {
            OpenCASCADE::NormalProjectionManifold<2, 3> normal_projector(
              bow_surface, tolerance);
            tria.set_manifold(1, normal_projector);


            break;
          }


@endcode 



  @p If  surface_projection_kind的值是 @p DirectionalProjection,  我们选择 OpenCASCADE::DirectionalProjectionManifold 类。新的网格点将在所考虑的单元格/边的arycenter处初始生成，然后沿着指定给 OpenCASCADE::DirectionalProjectionManifold 构造函数的方向在CAD表面上投影。在这种情况下，投影是沿着Y轴进行的。

@code
        case DirectionalProjection:
          {
            OpenCASCADE::DirectionalProjectionManifold<2, 3>
              directional_projector(bow_surface,
                                    Point<3>(0.0, 1.0, 0.0),
                                    tolerance);
            tria.set_manifold(1, directional_projector);


            break;
          }


@endcode 



作为第三个选项，如果 @p surface_projection_kind 的值是 @p NormalToMeshProjection, ，我们选择 OpenCASCADE::NormalToMeshProjectionManifold. ，新的网格点将再次在所考虑的单元格/边的arycenter处初始生成，然后沿着一个估计的网格法线方向投影到CAD表面。 OpenCASCADE::NormalToMeshProjectionManifold  构造函数只需要一个形状（至少包含一个面）和一个公差。

@code
        case NormalToMeshProjection:
          {
            OpenCASCADE::NormalToMeshProjectionManifold<2, 3>
              normal_to_mesh_projector(bow_surface, tolerance);
            tria.set_manifold(1, normal_to_mesh_projector);


            break;
          }


@endcode 



最后，我们使用良好的软件清洁性，确保这真的涵盖了 @p case 语句的所有可能选项。如果我们得到任何其他的值，我们就直接中止程序。

@code
        default:
          AssertThrow(false, ExcInternalError());
      }
  }



@endcode 




<a name="TriangulationOnCADrefine_mesh"></a> <h4>TriangulationOnCAD::refine_mesh</h4>. 




这个函数是全局细化网格的。在其他教程中，它通常也会分配自由度，并调整矩阵和向量的大小。这些任务在这里没有进行，因为我们没有在生成的三角网格上运行任何模拟。   


虽然这个函数看起来很简单，但这是我们对这个教程程序感兴趣的大部分工作实际发生的地方。特别是在细化定义船体表面的四边形和直线时，Triangulation类会询问我们分配给处理单个流形ID的各种对象，以确定新顶点应该位于何处。

@code
  void TriangulationOnCAD::refine_mesh()
  {
    tria.refine_global(1);
  }





@endcode 




<a name="TriangulationOnCADoutput_results"></a> <h4>TriangulationOnCAD::output_results</h4> 




输出我们的计算结果是一项相当机械的任务。这个函数的所有组成部分在之前已经讨论过了。

@code
  void TriangulationOnCAD::output_results(const unsigned int cycle)
  {
    const std::string filename =
      (output_filename + "_" + Utilities::int_to_string(cycle) + ".vtk");
    std::ofstream logfile(filename);
    GridOut       grid_out;
    grid_out.write_vtk(tria, logfile);
  }



@endcode 




<a name="TriangulationOnCADrun"></a> <h4>TriangulationOnCAD::run</h4> 




这是主函数。它的简短性应该是不言自明的。

@code
  void TriangulationOnCAD::run()
  {
    read_domain();


    const unsigned int n_cycles = 5;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        refine_mesh();
        output_results(cycle + 1);
      }
  }
} // namespace Step54



@endcode 




<a name="Themainfunction"></a> <h3>The main() function</h3> 




这是本程序的主要功能。它的基本结构与以前所有的教程程序一样，但通过新顶点放置的三种可能性来运行主类。

@code
int main()
{
  try
    {
      using namespace Step54;


      const std::string in_mesh_filename = "input/initial_mesh_3d.vtk";
      const std::string cad_file_name    = "input/DTMB-5415_bulbous_bow.iges";


      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in direction normal to CAD surface"
                << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::string        out_mesh_filename = ("3d_mesh_normal_projection");
      TriangulationOnCAD tria_on_cad_norm(in_mesh_filename,
                                          cad_file_name,
                                          out_mesh_filename,
                                          TriangulationOnCAD::NormalProjection);
      tria_on_cad_norm.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;


      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in y-axis direction" << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      out_mesh_filename = ("3d_mesh_directional_projection");
      TriangulationOnCAD tria_on_cad_dir(
        in_mesh_filename,
        cad_file_name,
        out_mesh_filename,
        TriangulationOnCAD::DirectionalProjection);
      tria_on_cad_dir.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;


      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in direction normal to mesh elements"
                << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      out_mesh_filename = ("3d_mesh_normal_to_mesh_projection");
      TriangulationOnCAD tria_on_cad_norm_to_mesh(
        in_mesh_filename,
        cad_file_name,
        out_mesh_filename,
        TriangulationOnCAD::NormalToMeshProjection);
      tria_on_cad_norm_to_mesh.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
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


程序的执行会产生一系列的网格文件 <code>3d_mesh_*.vtk</code> ，我们可以用任何可以读取VTK文件格式的常用可视化程序进行可视化。

下表说明了采用正常投影策略得到的结果。表中的前两行显示的是逐步细化的网格的侧视图，覆盖在精确几何图形的非常精细的渲染上。深红色和浅红色的区域只是表示当前的网格或精细的几何体更接近观察者；这种区别没有任何特别深刻的意义。最后一排图片描述了第二排所示的相同网格的正面视图（镜像到几何体的两侧）。


  <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_5.png" alt="" width="400"></td>
  </tr>
</table>   

从图片中可以看出---也正如我们所预料的---正常的细化策略在应用于具有明显曲率变化的表面时，无法产生形状良好的元素。这一点在船体的球体上尤其明显，所有的新点都被放置在球体的上部，而下部则完全没有解决。

下表与上表的安排一样，说明了采用方向性投影方法得到的结果，其中选择的投影方向是Y轴（在每张图片的左下方用一个小的黄色箭头表示）。


  <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_5.png" alt="" width="400"></td>
  </tr>
</table>   

这些图像证实，用定向投影法得到的网格质量明显高于沿表面法线投影法得到的网格。然而，在球体底部观察到一些在Y方向上拉长的元素，那里的表面几乎与选择的投影方向平行。

最后的测试显示了使用面的法线投影的结果。

  <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_5.png" alt="" width="400"></td>
  </tr>
</table>   

这些图片证实了法线对网格投影的方法导致网格在整个细化步骤中保持均匀的间隔。同时，这些网格很好地表现了原始的几何形状，甚至在球体的底部区域也是如此，而使用定向投影仪或法线投影仪则不能很好地恢复。<a name="PlainProg"></a> <h1> The plain program</h1>  @include "step-54.cc"  。 

  */  
