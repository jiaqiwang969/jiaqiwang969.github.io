// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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
 * @page CodingConventions Coding conventions used throughout deal.II

在deal.II中，我们努力保持编程风格和提供的接口类型尽可能一致。为此，我们采用了一套编码约定，尽可能遵循这些约定。它们有两个部分：风格问题和我们称之为“防御编程”的东西，后者是让我们的代码帮助我们发现bug的一种尝试。在阅读时，重要的是要记住，风格不是上帝赋予的，也不是比其他任何一套惯例更好；他们的目的仅仅是尽可能保持deal.II的一致性。一致性减少了我们产生的错误的数量，因为我们可以，例如，总是假设输入参数在函数调用的输出参数之前。它们还简化了代码的阅读，因为通过查看一段代码的编写风格，一些事情已经变得很清楚了，而不必查看某个东西的确切定义。

<h3>Notes on deal.II 缩进 indentation</h3>

<p>deal.II 采用 <code>clang-format</code> 6.0 进行规范化缩进. 样式文件位于
@code
  \${SOURCE_DIR}/.clang-format
@endcode

<p>在提交之前，您应该运行
@code
  clang-format -i <file>
@endcode
在你的每个档案上。这将确保缩进符合本页中概述的样式准则。

这很麻烦。因此，更容易的是，你可以直接跑
@code
  make indent
@endcode

在设置要编译的库的任何目录中，缩进最近更改的所有源文件。如果要确保缩进对于所有提交都是正确的，那么可能需要设置一个预提交。一种方法是复制 <code>\${SOURCE_DIR}/contrib/git-hooks/pre-commit</code> 到
<code>\${SOURCE_DIR}/.git/hooks/pre-commit</code> 并确保它是可执行的。

如果您正在使用的系统安装了多个版本的 <code>clang-format</code> （或者如果它不在路径中），您应该将上面的 <code>make indent</code> 命令替换为
@code
  make DEAL_II_CLANG_FORMAT=/path/to/clang-6.0/clang-format indent
@endcode
指向正确的可执行文件。
</p>

<h3>风格问题</h3>

<ol>
<li> %返回某个数（单元数、自由度等）的函数应以<code>n_*</code>开头。 例如:
  SparsityPatternBase::n_nonzero_elements().</li>

<li> %设置位或标志的函数应以<code>set_*</code>开始；清除位或标志的函数应命名为<code>clear_*</code>。例如: CellAccessor::set_refine_flag().</li>

<li> 应该使用传统的逻辑运算符，而不是它们的英文等效运算符(i.e., 用 <code>&&</code>, <code>||</code>, 或 <code>!</code> 而非 <code>and</code>, <code>or</code>, and <code>not</code>).

<li> 在实现文件中，在每个函数之后，期望有三行空行，以实现更好的可读性。一个空行出现在对代码块进行分组的函数中，因为两个空行不足以明显区分代码是否属于两个不同的函数。</li>

<li> 每当整数变量只能采用非负值时，它就被标记为无符号 unsigned。这同样适用于只能返回正值或零值的函数。 例如: Triangulation::n_active_cells().</li>

<li> 每当函数的参数不改变时，就应该将其标记为const，即使是按值传递。通常，我们将输入参数标记为const。这有助于作为一个额外的文档工具来澄清参数（输入、输出或两者）的意图，并允许编译器在此类参数更改时发出警告，这通常是非自愿的或风格不佳的。</li>

<li> 每当函数不更改嵌入类/对象的任何成员变量时，都应将其标记为const。</li>

<li> %函数名和变量名不能只包含一个或两个字母，除非变量是纯计数索引。</li>

<li> 类型别名（<code>using</code> 声明）优先于 typedef 声明。</li>

<li> 使用GeometryInfo中的几何信息来获得每个单元的面数、每个单元的子单元数、与面3相邻的子单元的子索引等，而不是将它们直接作为 <code>2*dim</code>、<code>(1@<@<dim)</code> 和<code>{0,3}</code> 写入代码。这降低了出错的可能性并增强了代码的可读性。</li>

<li> 类声明的布局如下：首先是公共函数块，从构造函数开始，然后是析构函数。如果存在公共成员变量，则这些变量必须出现在构造函数之前。公共变量只能在常量（特别是静态和常量）或不可避免的情况下使用。
<br>
在公共成员之后，是受保护的成员，最后是私有成员。顺序如下：首先是变量，然后是函数。
<br>
Exceptions 特殊情况应在公共部分结束时在非公共部分开始前宣布。

对于既不是<code>static const</code>也不是<code>static constexpr</code>，我们 C++ 11风格的类成员初始化；而是这样代替
@code
  class Foo
  {
    int a = 42;
    int *b = nullptr;
  };
@endcode
  写成
@code
  class Foo
  {
    Foo();

    int a;
    int *b;
  };



  inline Foo::Foo()
  : a(42)
  , b(nullptr)
  {}
@endcode
  </li>

<li> 如果一个函数既有输入参数又有输出参数，通常输入参数应在输出参数之前，除非有充分的理由改变这个顺序(最常见的原因是使用默认值跟踪输入参数。）</li>

<li> Exceptions 异常用于内部参数检查和通过Assert宏进行一致性检查。如C++语言所做的异常处理（尝试/抛出/捕捉 try/throw/catch，并使用AssertThrow Macro）用于处理运行时错误（如I/O故障），这在任何情况下都必须进行，而不仅仅是在调试模式下。</li>

<li> 有时，通过使用几个非成员函数来实现一个类是有意义的，这些函数不是公共接口的一部分，并且只在当前源文件中调用。此类自由函数应放在内部命名空间中，其结构如下：
  @code
  namespace internal
  {
    namespace ClassNameImplementation
    {
      // free functions go here
    }
  }
  @endcode
  其中 <code>ClassName</code> 是调用类的名称.
</li>

<li> 类、名称空间 和 类型 通常使用大写字母来表示单词的开头（例如triiterator）&mdash; 有时称为<a href="http://en.wikipedia.org/wiki/Camel_case"><i>camel case</i></a> &mdash 函数和变量使用小写字母和下划线分隔单词。唯一的例外是三角剖分中的迭代器别名 Triangulation 和 DoFHandler（命名为cell_iterator, active_line_iterator等），以明确与标准库容器类的连接。</li>

<li> 对于具有多个模板参数的类，维度通常放在数据类型说明符之前, 即书写为 Point<dim,number> 而非 Point<number,dim>. </li>

<li> 在deal.II中有几个地方我们在头文件中使用前向声明。这样做的原因是，当我们只需要将某个类型标记为函数的参数时，我们可以通过不使用头来提高编译速度。deal.II中使用的约定是，如果我们只需要一个类型名，那么可以在我们需要的头中向前声明该类型；如果一个函数（或成员函数）可以返回一个值，那么该值类型的声明应该是可用的（通过包含必要的头）。例如，<code>deal.II/dofs/dof_handler.h</code> 包含<code>deal.II/dofs/dof_accessor.h</code>，这样就可以编写<code>dof_handler.begin_active()->is_active()</code>之类的内容，而不必显式包含声明<code>begin_active()</code>返回的对象类型的头。</li>

<li> 每个类必须有至少200页的文档；-） ;-)</li>

</ol>


<h3> 实例化 Instantiation of templated functions/classes</h3>

<p>deal.II中的大多数类和函数都是模板化的。这就带来了这样一个问题：如果有的话，这些对象是如何实例化的，在哪里实例化的。在整个deal.II中，我们采用以下公约：</p>

<ol>

<li> 如果我们可以枚举所有可能的模板参数（例如，维度只能是1、2或3），那么函数模板进入.cc文件，我们显式实例化所有可能的参数。用户将不需要看到这些函数模板，因为他们无论如何也不想为任何其他模板参数实例化这些函数。</li>

<li> 如果我们不能列举所有可能的模板参数（例如，向量类型-因为用户可能想定义自己的向量类型），但至少知道一些常见的使用情况，那么该函数就被放入.templates.h文件中。我们将其 \#include 到.cc文件中，并实例化所有公共参数的函数。对于几乎所有的用户来说，这都很好--他们只使用我们已经实例化的( vector，matrix, ...) 类型，对他们来说.templates.h文件不会有任何意义。它也不会减慢他们的编译速度，因为他们看到的任何东西都不会 \#include .templates.h文件。但是定义自己类型(vector, matrix, ...) 的用户可以通过包含.templates.h文件用自己的用户定义类型实例化模板函数。</li>

<li> 最后，如果我们不能预先假定模板参数将采用哪些值（例如，从Subscriptor派生的任何类都可以用作参数），那么函数的定义将在头文件的底部提供声明。定义应该用<code>\#ifndef DOXYGEN ... \#endif</code> 以防止 Doxygen 误判。</li>

</ol>

<p> 对于前两种情况，实例化指令在 <code>.inst.in</code> 文件中定义。它们由一个名为 expand_instantiations 的二进制文件（从<code>cmake/scripts/expand_instantiations.cc</code> 构建）进行处理，并根据您的配置通过cmake动态定义参数（请参阅构建目录中的<code>cmake/config/template-arguments.in</code> ）。正是这些<code>.inst</code>文件最终包含在相应的<code>.cc</code>文件中。</p>


<h3>Defensive programming</h3>

<p> 防御性编程是一个术语，我们在谈论编写代码时经常使用这个术语，而我们的思维定势是错误会发生。在这里，错误有两种表现：第一，我自己在写函数的时候会犯错误；第二，其他人在调用我的函数时可能会出错。在任何一种情况下，我都希望我的代码能够（I）尽可能地避免错误，（ii）编译器已经可以找到一些错误，以及（iii）剩余的错误相对容易找到，例如因为程序中止。防御性编程是一套使这些目标更有可能实现的策略。</p>

<p>
随着时间的推移，我们已经学会了一些技巧，其中一些我们在这里列出：
<ol>
<li> <i>断言参数的前提条件:</i> 人们总是用错误或荒谬的参数调用函数。作为原型示例，考虑向量加法的一个简单实现:
  @code
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
  @endcode
  虽然正确，但如果两个向量的大小不相同，则此函数将遇到问题。你认为用不同大小的向量来调用这个函数是愚蠢的吗？是的，当然是。但这种情况经常发生：人们忘记重新初始化一个向量，或者在不同的函数中重置它，等等。因此，如果你在这样一个不幸的情况下，可能需要很长时间才能弄清楚发生了什么，因为你可能只是读取未初始化的内存，或者你正在写入lhs向量实际上并不拥有的内存。两者都不会导致程序的立即终止，但您可能会在以后的某个时间出现随机错误。如果程序马上停在这里就容易多了。下面的实现正好可以做到这一点：
  @code
    Vector &
    operator+=(Vector       &lhs,
               const Vector &rhs)
    {
      Assert (lhs.size() == rhs.size(),
              ExcDimensionMismatch(lhs.size(), rhs.size());
      for (unsigned int i=0; i<lhs.size(); ++i)
        lhs(i) += rhs(i);
      return lhs;
    }
  @endcode
  <code>Assert</code> 宏用来确保条件在运行时为true，否则将打印包含由第二个参数编码的信息的字符串并中止程序。这样，当您编写一个新的程序来调用这个函数时，您将立即了解到您的错误，并且有机会在不必认真调试任何东西的情况下修复它。
  <p>
  一般来说，无论何时实现一个新函数，都要考虑参数的前提条件<i>preconditions</i>，即函数希望它们中的每一个或它们的组合都是正确的。然后为所有这些前提条件编写断言。在某些情况下，这可能是六个断言，但请记住，每个断言都是通过琐碎的方法已经发现的潜在bug。
  <p>
  在最后一点中，让我们注意断言当然是昂贵的：当您将程序链接到库的调试版本时，断言可能会使程序慢3到5倍。但是，如果考虑到您的总体开发时间，快速发现bug的能力可能远远超过您等待程序完成所花费的时间。此外，对Assert宏的调用将在优化模式下从程序中删除（假设只有在知道调试模式下一切正常运行时才使用该模式）。优化后的库比调试库快3-5倍，但代价是发现bug要困难得多。
</li>

<li> <i> 断言后置条件 :</i> 如果一个函数计算一些非平凡的东西，那么代码中可能会有一个bug。要找到这些，请使用后置条件：就像您对输入参数的有用值有一定的了解一样，您也知道您期望的可能返回值是什么。例如，计算向量范数的函数期望范数为正。你可以这样写：
  @code
    double norm(const Vector &v)
    {
      double s = 0;
      for (unsigned int i=0; i<v.size(); ++i)
        s += v(i) * v(i);

      Assert (s >= 0, ExcInternalError());
      return std::sqrt(s);
    }
  @endcode
  这个函数太简单了，无法真正证明这个断言是正确的，但是想象一下计算要长一些，您可以看到断言如何帮助您确保（或避免）自己出错。请注意，有人可能会认为，一旦我们运行了多次程序并发现该条件从未触发，就应该删除断言。但最好还是把它放在原处：它为将来（和读者）对函数的知识编码；如果有人出现并用更有效的算法替换函数的实现，断言可以帮助确保函数继续执行它应该执行的操作。
</li>

<li> <i>断言内部状态:</i> 同样的道理，如果你有一个复杂的算法，使用断言来确保你的心理模型与事实相符。例如，假设您正在编写一个函数，以确保网格大小不会在局部发生太大变化。您可能会得到以下类型的代码：
  @code
    for (const auto &cell = triangulation.active_cell_iterators())
      for (unsigned int face=0; ...)
        {
          if (something)
            { ... }
          else
            {
              // we have a cell whose neighbor must
              // be at the boundary if we got here
            }
        }
  @endcode
  使我们进入else分支的条件可能很复杂，虽然我们认为这里唯一的可能性是邻居在边界上，但我们的实现中可能有一个bug。我们的想法中可能也有一个bug，或者有人在相同的函数中更改了上面的代码，忘记了这里的问题，或者在库中完全不同的位置进行了更改，使得这个假设站不住脚。在所有这些情况下，我们断言的明确声明确保这些问题很容易被发现。
  </li>

<li> <i>如果变量位于堆栈上，则在其声明点初始化变量:</i>
  传统的C语言要求在函数的开头声明变量，即使它们只在下面的更进一步使用。这就产生了我们可以想象的一维代码：:
  @code
    template <int dim>
    void foo ()
    {
      Point<dim> cell_center;
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  问题是，如果声明和初始化之间的代码既长又复杂，您就无法在一个页面上查找变量的类型和值。事实上，甚至可能不太清楚这个变量是用来初始化的，或者它是否被意外地忽略了初始化。
  <p>
  更好的方法如下：
  @code
    template <int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          Point<dim> cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  这样就更清楚了变量的类型是什么，而且它实际上只在初始化时使用。此外，如果有人想阅读代码来查看变量实际上在做什么，那么在最内部可能的作用域中声明和初始化它会使这项任务变得更容易：我们不必在声明之外向上查找它，也不必在当前作用域的末尾向下查找，因为这是变量消亡的地方。
  <p>
  最后一点，很明显，您只能对完全位于堆栈上的变量执行这类操作，而无需在堆上分配内存。在deal.II中，这仅适用于<code>int, double, char</code>等内置类型，以及 Point 和 Tensor 类。其他所有东西都有类似<code>std::vector</code>的成员变量，这需要内存分配 &mdash;您不希望在循环内声明这些，至少在循环频繁遍历的情况下是这样。
  </li>

<li> <i>使变量成为 const:</i> 为了学习上面的例子，请注意，在大多数情况下，我们永远不会再更改这样初始化的变量。换言之，如果是这样的话，我们不妨这样写：
  @code
    template <int dim>
    void foo ()
    {
      ... // something lengthy and complicated
      for (const auto &cell = dof_handler.active_cell_iterators())
        {
          const Point<dim> cell_center = (cell->vertex(0) +
                                          cell->vertex(1)) / 2;
          ...
        }
      ...
    }
  @endcode
  通过将变量标记为常量，我们可以确保不会意外更改它。例如，编译器可以捕获如下代码：
  @code
        if (cell_center[0] = 0)
          ...
  @endcode
  这很可能是一个 <code>==</code> 而不是一个赋值。通过将变量标记为const，编译器就会告诉我们这个bug。也许同样重要的是，代码的人类读者不需要进一步了解变量的值是否真的在声明和使用之间发生了更改—如果它被标记为const，就不可能发生更改。
</li>

<li> <i> 使函数的输入参数为常量:</i> 对于函数参数也是如此：如果您不想更改变量（通常是输入参数的情况），那么将其标记为常量。例如，以下函数应将其参数作为常量值：
  @code
     template <int dim>
     typename Triangulation<dim>::cell_iterator
     CellAccessor<dim>::child(const unsigned int child_no)
     {
       ...
       return something;
     }
  @endcode
  例如，在这里，用户调用<code>cell-@>child(3)</code>。函数没有理由改变子参数的值，所以把它标记为常量：这可以帮助代码的读者理解这是函数的一个输入参数，我们不需要在下面搜索它是否被改变，如果我们不小心更改了值，它可以帮助编译器发现错误。
</ol>

 */
