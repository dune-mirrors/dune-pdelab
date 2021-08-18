// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_TREE_TO_GRAPHVIZ_HH
#define DUNE_PDELAB_TREE_TO_GRAPHVIZ_HH

#include <type_traits>
#include <iostream>

#include <dune/common/typetraits.hh>
#include <dune/common/classname.hh>

#include <dune/typetree/visitor.hh>
#include <dune/typetree/traversal.hh>

namespace Dune::PDELab {

  namespace Impl {

    template<class F>
    struct TreeToGraphvizVisitor : public TypeTree::TreeVisitor
                                 , public TypeTree::DynamicTraversal
    {
      template<class F_>
      TreeToGraphvizVisitor(std::ostream& out, F_&& f) : _out{out}, _f{std::forward<F_>(f)}{
        _out << "digraph {" << std::endl;
      }

      ~TreeToGraphvizVisitor() {
        _out << "}" << std::endl;
      }

      template<typename T, typename TreePath>
      void pre(const T& node, TreePath treePath) const {
        _out << "\"" << treePath <<"\" [label=\"" << _f(node, treePath) << "\" shape=box]"  << std::endl;
      }

      template<typename T, typename Child, typename TreePath, typename ChildIndex>
      void beforeChild(T&& t, Child&& child, const TreePath& treePath, ChildIndex childIndex) {
        _out << "\"" << treePath << "\" -> \"" << push_back(treePath,childIndex) << "\"" << std::endl;
      }

      template<typename T, typename TreePath>
      void leaf(const T& node, const TreePath& treePath) {
        _out << "\"" << treePath <<"\" [label=\"" << _f(node, treePath) << "\" shape=box]"  << std::endl;
      }

      std::ostream& _out;
      F _f;
    };

      template<class Node>
      using HasNameConcept = decltype((std::string{std::declval<Node>().name()},true));

  }

  template<class Tree, class F>
  void writeTreeToGraphviz(std::ostream& out, const Tree& tree, F&& f){
    TypeTree::applyToTree(tree, Impl::TreeToGraphvizVisitor<std::decay_t<F>>{out, std::forward<F>(f)});
  }

  template<class Tree>
  void writeTreeToGraphviz(std::ostream& out, const Tree& tree){
    writeTreeToGraphviz(out, tree, [](const auto& node, auto& path) -> std::string {
      using Node = std::decay_t<decltype(node)>;
      std::stringstream ss;
      ss << path;
      if constexpr (Std::is_detected<Impl::HasNameConcept,Node>{}){
        std::string name = node.name();
        return name.empty() ? ss.str() : name;
      } else
        return ss.str();
    });
  }

} // end namespace Dune::PDELab

#endif // DUNE_PDELAB_TREE_TO_GRAPHVIZ_HH
