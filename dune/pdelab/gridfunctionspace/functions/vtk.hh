// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_VTK_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_VTK_HH

#include <cstddef>
#include <memory>
#include <sstream>
#include <string>

#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/typetree/childaccess.hh>
#include <dune/common/typetree/nodeconcepts.hh>
#include <dune/common/typetree/treepath.hh>

#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridfunctionspace/functions/discretegridfunction.hh>
// primary templates; this header is included at its end
#include <dune/pdelab/gridfunctionspace/vtk.hh>

/** \file
 * \brief VTK output of coefficient vectors of dune-functions bases.
 *
 * Overload of Dune::PDELab::addSolutionToVTKWriter() for dune-functions
 * bases, with the same output as for grid function spaces: every leaf of the
 * basis tree is written as a scalar function, except for the leaves below a
 * power node with a static number of leaf children, which are written as one
 * vector-valued function (as the children of a VectorGridFunctionSpace).
 *
 * \warning This is experimental.
 */

namespace Dune::PDELab {

  namespace Impl::Functions {

    // A power node with compile-time degree and leaf children, written as a vector
    template<class Node>
    constexpr bool isVTKVectorNode()
    {
      if constexpr (Dune::TypeTree::Concept::StaticDegreeInnerTreeNode<Node>
                    and Dune::TypeTree::Concept::UniformInnerTreeNode<Node>)
        return Dune::TypeTree::Concept::LeafTreeNode<Dune::TypeTree::Child<Node,0>>;
      else
        return false;
    }

    // Component names of bases are built from the tree path: solution, solution_0, solution_0_1, ...
    template<class TreePath>
    std::string vtkComponentName(const TreePath& tp)
    {
      std::stringstream name;
      name << "solution";
      for (std::size_t i = 0; i < tp.size(); ++i)
        name << "_" << tp[i];
      return name.str();
    }

    template<class VTKWriter, class Basis, class X, class NameGenerator, class Predicate, class Node, class TreePath>
    void addBasisSolutionToVTKWriter(VTKWriter& vtkWriter, const Basis& basis, const std::shared_ptr<const X>& x,
                                     const NameGenerator& nameGenerator, const Predicate& predicate,
                                     const Node& node, const TreePath& tp)
    {
      auto addVertexData = [&](auto dgf) {
        using DGF = typename decltype(dgf)::element_type;
        const std::string name = nameGenerator(vtkComponentName(tp), tp);
        vtkWriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF>>(dgf, name.c_str()));
      };
      if constexpr (Dune::TypeTree::Concept::LeafTreeNode<Node>)
      {
        if (predicate(node, tp))
        {
          auto subBasis = std::make_shared<const decltype(Dune::Functions::subspaceBasis(basis, tp))>(Dune::Functions::subspaceBasis(basis, tp));
          using SubBasis = typename decltype(subBasis)::element_type;
          addVertexData(std::make_shared<DiscreteGridFunction<SubBasis,X>>(subBasis, x));
        }
      }
      else if constexpr (isVTKVectorNode<Node>())
      {
        if (predicate(node, tp))
        {
          auto subBasis = std::make_shared<const decltype(Dune::Functions::subspaceBasis(basis, tp))>(Dune::Functions::subspaceBasis(basis, tp));
          using SubBasis = typename decltype(subBasis)::element_type;
          addVertexData(std::make_shared<VectorDiscreteGridFunction<SubBasis,X>>(subBasis, x));
        }
      }
      else
        Dune::Hybrid::forEach(Dune::range(node.degree()), [&](auto i) {
          addBasisSolutionToVTKWriter(vtkWriter, basis, x, nameGenerator, predicate,
                                      node.child(i), Dune::TypeTree::push_back(tp, i));
        });
    }

  } // namespace Impl::Functions

  /** \brief Add the functions of a coefficient vector of a dune-functions basis to a VTKWriter
   *
   * \param vtkWriter     VTKWriter or SubsamplingVTKWriter
   * \param basis         dune-functions global basis
   * \param x             coefficient vector, e.g. Backend::Vector<Basis,...>; referenced until the file is written
   * \param nameGenerator name scheme, e.g. vtk::defaultNameScheme(); as bases have no names,
   *                      it is called with the component name "solution_<tree path>"
   * \param predicate     predicate(node, treePath) selects the functions to write
   */
  template<typename VTKWriter, typename Basis, typename X,
           typename NameGenerator = vtk::DefaultFunctionNameGenerator,
           typename Predicate = vtk::DefaultPredicate>
    requires Impl::DuneFunctionsBasis<Basis>
  void addSolutionToVTKWriter(VTKWriter& vtkWriter,
                              const Basis& basis,
                              const X& x,
                              const NameGenerator& nameGenerator = vtk::defaultNameScheme(),
                              const Predicate& predicate = Predicate())
  {
    auto localView = basis.localView();
    Impl::Functions::addBasisSolutionToVTKWriter(vtkWriter, basis, Dune::stackobject_to_shared_ptr(x),
                                                 nameGenerator, predicate,
                                                 localView.tree(), Dune::TypeTree::treePath());
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_VTK_HH
