#ifndef DUNE_PDELAB_BASIS_UTILITY_HH
#define DUNE_PDELAB_BASIS_UTILITY_HH

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>

#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

namespace Dune::PDELab
{

  namespace Impl
  {

    std::string formatMultiIndex(const auto &mi)
    {
      std::stringstream s;
      s << "\"(";
      for (std::size_t i = 0; i < mi.size(); ++i)
      {
        if (i > 0)
          s << ",";
        s << mi[i];
      }
      s << ")\"";
      return s.str();
    };

    template <class Basis>
    void printIndexTree(std::ofstream &out,
                        const Basis &basis,
                        typename Basis::SizePrefix prefix,
                        const std::set<typename Basis::SizePrefix> &indices)
    {

      using SizePrefix = typename Basis::SizePrefix;
      SizePrefix child_prefix = prefix;

      auto size = basis.size(prefix);
      if (size == 0)
      {
        if (not indices.contains(prefix))
          out << "  " << formatMultiIndex(prefix) << " [shape=box, style=dashed] " << std::endl;
        return;
      }
      out << "  " << formatMultiIndex(prefix) << " [style=filled, fillcolor=deepskyblue, shape=box] " << std::endl;
      child_prefix.push_back(0);
      for (std::size_t i = 0; i < size; ++i)
      {
        child_prefix.back() = i;
        std::string arrow_sty = "";
        if (basis.size(child_prefix) == 0 and not indices.contains(child_prefix))
          arrow_sty = " [style=dashed]";
        out << "  " << formatMultiIndex(prefix) << " -> " << formatMultiIndex(child_prefix) << arrow_sty << std::endl;
        printIndexTree(out, basis, child_prefix, indices);
      }
    }

  } // namespace Impl

  /**
   * @brief Write the index tree of a given basis coefficients in Graphviz DOT format.
   *
   * This function traverses the index tree of the provided basis and outputs its structure
   * to the given output stream in DOT language, suitable for visualization with Graphviz.
   * The nodes represent multi-indices in the basis, and edges represent parent-child relationships
   * in the index tree. Nodes are styled and colored based on their type (leaf or non-leaf).
   *
   * Inner nodes are filled with blue color, leaf nodes with green color, and ghost or empty
   * nodes are represented with dashed borders. Ghost or empty nodes are those that do not
   * correspond to any actual basis function but are included due to unremoved nodes in
   * the index tree after topological interleavings.
   *
   * @tparam Basis The type of the basis, which must model Dune::Functions::Concept::GlobalBasis.
   * @param filename The name of the file to write the DOT graph to.
   * @param basis The basis whose index tree will be printed.
   * @param prefix The current multi-index prefix (used for recursion; defaults to empty).
   * @param indices Set of multi-indices representing all leaf nodes (used for styling; defaults to empty).
   * @return The number of children for the current prefix node.
   *
   * @note The output DOT graph can be rendered using Graphviz tools (e.g., `cat output.txt | dot -Tsvg > output.svg`).
   * @note Use this function for debugging and visualizing the structure of a small basis.
   */
  template <class Basis>
    requires(bool(Dune::models<Dune::Functions::Concept::GlobalBasis<typename Basis::GridView>, Basis>()))
  void writeIndexTree(const std::string &filename, const Basis &basis)
  {
    std::set<typename Basis::SizePrefix> indices;
    std::ofstream out(filename);

    out << "digraph {" << std::endl;
    indices.clear();
    auto localView = basis.localView();
    for (const auto &e : elements(basis.gridView()))
    {
      localView.bind(e);
      for (typename Basis::size_type i = 0; i < localView.size(); ++i)
      {
        auto mi = localView.index(i);
        typename Basis::SizePrefix si;
        for (auto x : mi)
          si.push_back(x);
        indices.insert(si);
      }
    }
    for (const auto &mi : indices)
      out << "  " << Impl::formatMultiIndex(mi) << " [style=filled, fillcolor=darkolivegreen2, shape=box] " << std::endl;

    Impl::printIndexTree(out, basis, {}, indices);

    out << "}" << std::endl;
    out.close();
  }

  namespace ContainerDescriptors
  {

    using namespace Dune::Functions::ContainerDescriptors;

    namespace Impl
    {

      template <class, class = void>
      struct BlockType
      {
        using type = Unknown;
      };

      template <class ChildDescriptor>
      struct BlockType<Vector<ChildDescriptor>> : std::type_identity<ChildDescriptor>
      {
      };

      template <class ChildDescriptor>
      struct BlockType<UniformVector<ChildDescriptor>> : std::type_identity<ChildDescriptor>
      {
      };

      template <class ChildDescriptor, std::size_t i>
      struct BlockType<Array<ChildDescriptor, i>> : std::type_identity<ChildDescriptor>
      {
      };

      template <class ChildDescriptor, std::size_t i>
      struct BlockType<UniformArray<ChildDescriptor, i>> : std::type_identity<ChildDescriptor>
      {
      };

      template <class... ChildDescriptors>
      struct BlockType<Tuple<ChildDescriptors...>, std::void_t<std::common_type_t<ChildDescriptors...>>>
          : std::type_identity<std::common_type_t<ChildDescriptors...>>
      {
      };

      template <class Descriptor>
      struct IsCompileTimeUniform : std::false_type
      {
      };

      template <>
      struct IsCompileTimeUniform<Value> : std::true_type
      {
      };

      template <class Child, std::size_t n>
      struct IsCompileTimeUniform<UniformArray<Child, n>> : std::true_type
      {
      };
    }

    //! The block type of a descriptor is the type of a single block of the descriptor, if any.
    template <class Descriptor>
    using BlockType = typename Impl::BlockType<Descriptor>::type;

    //! Whether the descriptor has compile-time uniform size.
    template <class Descriptor>
    constexpr bool IsCompileTimeUniform = Impl::IsCompileTimeUniform<Descriptor>::value;

  } // namespace ContainerDescriptors

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_UTILITY_HH
