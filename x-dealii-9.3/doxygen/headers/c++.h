// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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
 * @defgroup CPP11 deal.II and the C++11 standard
 *
 *从版本9起，Del.II需要至少支持 <a href="http://en.wikipedia.org/wiki/C%2B%2B11">C++11</a> 的编译器。作为这一部分的一部分，Del.II内部实现的许多地方现在使用的是仅在C++ 11中引入的特性。也就是说，Del.II还具有使C++与C++ 11更容易使用的功能和类。
 *
 *一个例子是支持基于C++的11范围的循环<a href="http://en.wikipedia.org/wiki/C++11#Range-based_for_loop">range-based for loops</a>。基于deal.II的代码通常有许多这样的循环
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   typename Triangulation<dim>::active_cell_iterator
 *     cell = triangulation.begin_active(),
 *     endc = triangulation.end();
 *   for (; cell!=endc; ++cell)
 *     cell->set_refine_flag();
 * @endcode
 *使用C++ 11的基于循环的范围，现在可以将其写入如下：
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   for (auto &cell : triangulation.active_cell_iterators())
 *     cell->set_refine_flag();
 * @endcode
 *这依赖于如下函数 Triangulation::active_cell_iterators(),类似的还有DoFHandler类中的
 DoFHandler::active_cell_iterators(), hp::DoFHandler::active_cell_iterators().这些函数的一些变体为所有单元格（不仅仅是活动单元格）和各个级别的单元格提供迭代器范围。
 *
 *库中还有许多其他函数允许习惯性地使用基于范围的for循环。例如，GeometryInfo::face_indices(), GeometryInfo::vertex_indices(), FEValuesBase::quadrature_point_indices(), 等等.
 *
 *C++ 11还介绍了[constexpr](https://en.cppreference.com/w/cpp/language/constexpr)和函数的概念。定义为constexpr的变量是在程序编译期间计算的常量值，因此与初始化相关的运行时成本为零。此外，constexpr常量还正确定义了生存期，从而完全防止了所谓的“静态初始化顺序失败”。函数可以标记为constexpr，表示如果它们的输入参数是常量表达式，则可以生成编译时常量返回值。此外，至少有一个constexpr构造函数的类可以初始化为`constexpr`。
 * 
 *例如，由于构造函数 Tensor::Tensor(const array_type &) 是`constexpr`，我们可以在编译时用数组初始化一个张量，如下所示：
* @code
 * constexpr double[2][2] entries = {{1., 0.}, {0., 1.}};
 * constexpr Tensor<2, 2> A(entries);
 * @endcode
 *这里，A的内容不存储在堆栈上。相反，它们在编译时初始化并插入到可执行程序的.data部分。程序可以在运行时使用这些值，而无需花费时间进行初始化。初始化张量可以简化为一行。
* @code
 * constexpr Tensor<2, 2> A({{1., 0.}, {0., 1.}});
 * @endcode
 *某些函数，如 determinant()被指定为CONTXPR，但它们需要具有C++ 14能力的编译器。因此，该函数在内部声明为：
* @code
 * template <int dim, typename Number>
 * DEAL_II_CONSTEXPR Number determinant(const Tensor<2, dim, Number> &t);
 * @endcode
 *如果编译器支持C++14，那么宏macro @ref DEAL_II_CONSTEXPR 将被简化为 `constexpr` 。不然, 会被忽略。因此，使用较新的编译器，用户可以编写
 * @code
 * constexpr double det_A = determinant(A);
 * @endcode
 *假设`A`用 `constexpr`说明符声明了。这个例子展示了使用constexpr的性能提升，因为在这里我们在编译时执行了一个复杂度为 $O(\text{dim}^3)$ 的操作，避免了任何运行时开销。
 */



/**
 *Del.II目前只需要一个C++ 11兼容编译器，但是有许多来自C++ 14标准的函数和类，在编译器只支持C++ 11的情况下，也很容易提供。这些是在当前命名空间中收集的。
 *
 *最为著名的例子是 <a
 *href="https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique">`std::make_unique`</a>
 *函数，它可以被忽略为没有包含在C++11中(因为在C++11中存在 <a
 * href="https://en.cppreference.com/w/cpp/memory/shared_ptr/make_shared">`std::make_shared`</a>).
 *
 *在这个命名空间中还有其他的小的添加，允许我们在这一点上使用C++ 14的特性，即使我们不需要C++ 14兼容编译器。
 *
 * @note 如果使用中的编译器实际上支持C++ 14，那么这个命名空间的内容就是从命名空间“STD”中简单导入的类和函数。也就是说，我们退回到编译器提供的内容，而不是我们自己的实现。
 */
namespace std_cxx14
{}



/**
 *Del.II目前只需要一个C++ 11兼容编译器，但是有许多来自C++ 17标准的函数和类，在编译器只支持C++ 11的情况下，也很容易提供。这些是在当前命名空间中收集的。
 *
 *最显著的例子是<a
 * href="https://en.cppreference.com/w/cpp/utility/optional">`std::optional`</a>从C++ 17标准开始的C++类的可选类。
 *
 *在这个命名空间中还有其他的小的添加，允许我们在这一点上使用C++ 17的特性，即使我们不需要C++ 17兼容编译器。
 *
 * @note 如果使用中的编译器实际上支持C++ 17，那么这个命名空间的内容只是从命名空间STD导入的类和函数，也就是说，我们回到编译器所提供的，而不是我们自己的实现。
 */
namespace std_cxx17
{}



/**
 *Del.II目前只需要一个C++ 11兼容编译器，但是有许多来自C++ 20标准的函数和类，在编译器只支持C++ 11的情况下，也很容易提供。这些是在当前命名空间中收集的。
 *
 *其中的一个例子是 <a
 * href="https://en.cppreference.com/w/cpp/ranges/iota_view">`std::ranges::iota_view`</a>
 *这个类从C++ 20标准开始，将其引入到C++中。它用作 GeometryInfo::face_indices(), GeometryInfo::vertex_indices(), 和 FEValuesBase::quadrature_point_indices() functions, 支持基于范围的for循环。(参考 @ref CPP11 的介绍）
 *
 *在这个命名空间中还有其他的小的添加，允许我们在这一点上使用C++ 20的特性，即使我们不需要C++ 20兼容编译器。
 *
 * @note 如果使用中的编译器实际上支持C++ 20，那么这个命名空间的内容只是从命名空间STD导入的类和函数，也就是说，我们回到编译器所提供的，而不是我们自己的实现。
 */
namespace std_cxx20
{}
