// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ADAPTIVITY_ADAPTIVITY_HH
#define DUNE_PDELAB_ADAPTIVITY_ADAPTIVITY_HH

#include<dune/common/exceptions.hh>

#include<array>
#include<limits>
#include<vector>
#include<map>
#include<unordered_map>
#include<dune/common/dynmatrix.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include<dune/pdelab/backend/common/spacetraits.hh>
#include<dune/pdelab/gridoperator/functions/localfunctionspace.hh>

#include<dune/common/typetree/childaccess.hh>
#include<dune/common/typetree/traversal.hh>
#include<dune/functions/backends/istlvectorbackend.hh>
#include<dune/functions/functionspacebases/nodes.hh>

#include<dune/pdelab/common/function.hh>
// for InterpolateBackendStandard
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
// for intersectionoperator
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

namespace Dune {
  namespace PDELab {


    template<typename GFS>
    struct LeafOffsetCache
    {

      typedef typename GFS::Traits::GridView::template Codim<0>::Entity Cell;
      typedef LocalFunctionSpace<GFS> LFS;

        // we need an additional entry because we store offsets and we also want the
        // offset after the last leaf for size calculations
      typedef std::array<std::size_t,TypeTree::TreeInfo<GFS>::leafCount + 1> LeafOffsets;

      const LeafOffsets& operator[](GeometryType gt) const
      {
        const LeafOffsets& leaf_offsets = _leaf_offset_cache[GlobalGeometryTypeIndex::index(gt)];
        // make sure we have data for this geometry type
        assert(leaf_offsets.back() > 0);
        return leaf_offsets;
      }

      void update(const Cell& e)
      {
        LeafOffsets& leaf_offsets = _leaf_offset_cache[GlobalGeometryTypeIndex::index(e.type())];
        if (leaf_offsets.back() == 0)
          {
            _lfs.bind(e);
            extract_lfs_leaf_sizes(_lfs,leaf_offsets.begin()+1);
            // convert to offsets
            std::partial_sum(leaf_offsets.begin(),leaf_offsets.end(),leaf_offsets.begin());
            // sanity check
            assert(leaf_offsets.back() == _lfs.size());
          }
      }

      explicit LeafOffsetCache(const GFS& gfs)
        : _lfs(gfs)
        , _leaf_offset_cache(GlobalGeometryTypeIndex::size(Cell::dimension))
      {}

      LFS _lfs;
      std::vector<LeafOffsets> _leaf_offset_cache;

    };


#ifndef DOXYGEN
    namespace Impl::Adapt {

      /* Access to the space-specific parts of the solution transfer. The
       * transfer algorithm below works on grid function spaces and on
       * dune-functions bases; it only needs
       * - the grid and the element range of the space,
       * - the local coefficients of a leaf element (read and add),
       * - a tree of local finite elements for any element of the hierarchy
       *   (leaves with finiteElement(), size() and localIndex(i)),
       * - the update of the space after grid.adapt().
       */

      template<class Space>
      auto elementRange(const Space& space)
      {
        if constexpr (Impl::DuneFunctionsBasis<Space>)
          return elements(space.gridView(), Partitions::interior);
        else
          return elements(space.entitySet(), Partitions::interior);
      }

      template<class Space>
      const auto& grid(const Space& space)
      {
        return space.gridView().grid();
      }

      //! update the space after the grid has been adapted
      template<class Space, class Grid>
      void update(Space& space, Grid& grid)
      {
        if constexpr (Impl::DuneFunctionsBasis<Space>)
          space.update(grid.leafGridView());
        else
          space.update(true);
      }

      //! local coefficients of a leaf element: grid function space
      template<class Space>
      class LocalCoefficients
      {
        using LFS = LocalFunctionSpace<Space>;
        using Cache = LFSIndexCache<LFS>;

      public:
        explicit LocalCoefficients(const Space& space)
          : lfs_(space), cache_(lfs_)
        {}

        template<class Element>
        void bind(const Element& e)
        {
          lfs_.bind(e);
          cache_.update();
        }

        const LFS& tree() const { return lfs_; }
        std::size_t size() const { return lfs_.size(); }

        template<class X, class C>
        void read(const X& x, C& local) const
        {
          typename X::template ConstLocalView<Cache> view(x);
          view.bind(cache_);
          local.resize(size());
          view.read(local);
          view.unbind();
        }

        template<class X, class C>
        void add(X& x, const C& local) const
        {
          typename X::template LocalView<Cache> view(x);
          view.bind(cache_);
          view.add(local);
          view.commit();
          view.unbind();
        }

      private:
        LFS lfs_;
        Cache cache_;
      };

      //! local coefficients of a leaf element: dune-functions basis
      template<class Space>
        requires Impl::DuneFunctionsBasis<Space>
      class LocalCoefficients<Space>
      {
        using LFS = Impl::Functions::LocalFunctionSpace<Space>;

      public:
        explicit LocalCoefficients(const Space& space)
          : lfs_(space)
        {}

        template<class Element>
        void bind(const Element& e) { lfs_.bind(e); }

        const auto& tree() const { return lfs_.tree(); }
        std::size_t size() const { return lfs_.size(); }

        template<class X, class C>
        void read(const X& x, C& local) const
        {
          auto xb = Dune::Functions::istlVectorBackend(Backend::native(x));
          local.resize(size());
          for (std::size_t i = 0; i < size(); ++i)
            local[i] = xb[lfs_.index(i)];
        }

        template<class X, class C>
        void add(X& x, const C& local) const
        {
          auto xb = Dune::Functions::istlVectorBackend(Backend::native(x));
          for (std::size_t i = 0; i < size(); ++i)
            xb[lfs_.index(i)] += local[i];
        }

      private:
        LFS lfs_;
      };

      //! local finite elements of any element of the grid hierarchy: grid function space
      template<class Space>
      class ElementTree
      {
        using Element = typename Space::Traits::EntitySet::Element;

      public:
        explicit ElementTree(const Space& space) : lfs_(space) {}

        // A local function space only provides finite elements on elements of
        // its entity set. The layout (leaf sizes and local indices) is taken from
        // a leaf descendant of the same geometry type, the finite elements from
        // the finite element maps (as the original implementation did).
        void bind(const Element& e)
        {
          element_ = e;
          const int maxLevel = lfs_.gridFunctionSpace().gridView().grid().maxLevel();
          for (const auto& d : descendantElements(e, maxLevel))
            if (d.isLeaf() and d.type() == e.type())
            {
              lfs_.bind(d);
              return;
            }
          DUNE_THROW(Exception, "No leaf descendant of the same geometry type found");
        }

        const LocalFunctionSpace<Space>& tree() const { return lfs_; }

        template<class Leaf>
        decltype(auto) finiteElement(const Leaf& leaf) const
        {
          return leaf.gridFunctionSpace().finiteElementMap().find(element_);
        }

      private:
        LocalFunctionSpace<Space> lfs_;
        Element element_;
      };

      //! local finite elements of any element of the grid hierarchy: dune-functions basis
      template<class Space>
        requires Impl::DuneFunctionsBasis<Space>
      class ElementTree<Space>
      {
        using Tree = std::decay_t<decltype(std::declval<const Space&>().preBasis().makeNode())>;

      public:
        explicit ElementTree(const Space& space)
          : tree_(space.preBasis().makeNode())
        {
          Dune::Functions::initializeTree(tree_);
        }

        // only the finite elements are set up, no global indices are computed,
        // so this also works for elements that are not in the leaf grid view
        template<class Element>
        void bind(const Element& e) { Dune::Functions::bindTree(tree_, e); }

        const Tree& tree() const { return tree_; }

        template<class Leaf>
        const auto& finiteElement(const Leaf& leaf) const { return leaf.finiteElement(); }

      private:
        Tree tree_;
      };

    } // namespace Impl::Adapt
#endif // DOXYGEN


    /*! @class L2Projection
     *
     * @brief Cached inverse local mass matrices for the local L2 projection
     *        onto the elements of a coarser grid level
     *
     * @tparam GFS Type of ansatz space (grid function space or dune-functions basis)
     * @tparam U   Container class for the solution
     */
    template<class GFS, class U>
    class L2Projection
    {
      typedef typename U::ElementType DF;

    public:

      typedef DynamicMatrix<typename U::ElementType> MassMatrix;

      /*! @brief The constructor.
       *
       * @param gfs      The ansatz space
       * @param intorder Integration order of the mass matrices and the projection
       */
      explicit L2Projection(const GFS& gfs, int intorder = 2)
        : _intorder(intorder)
      {}

      //! Integration order of the mass matrices and the projection
      int integrationOrder() const { return _intorder; }

      /*! @brief Inverse local mass matrix of a leaf of the space tree on an element
       *
       * @param e         The element; the matrix is cached per geometry type
       * @param leafIndex Index of the leaf in the space tree
       * @param fe        Local finite element of this leaf on e
       */
      template<class Element, class FiniteElement>
      const MassMatrix& inverseMassMatrix(const Element& e, std::size_t leafIndex, const FiniteElement& fe)
      {
        auto& matrices = _inverse_mass_matrices[GlobalGeometryTypeIndex::index(e.type())];
        if (matrices.size() <= leafIndex)
          matrices.resize(leafIndex+1);
        MassMatrix& mass_matrix = matrices[leafIndex];
        // if the matrix isn't empty, it has already been cached
        if (mass_matrix.N() > 0)
          return mass_matrix;

        constexpr int dim = Element::dimension;
        using Range = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
        const std::size_t local_size = fe.localBasis().size();
        mass_matrix.resize(local_size,local_size);
        std::vector<Range> phi(local_size);
        for (const auto& ip : QuadratureRules<DF,dim>::rule(e.type(),_intorder))
          {
            fe.localBasis().evaluateFunction(ip.position(),phi);
            const DF factor = ip.weight();
            for (std::size_t i = 0; i < local_size; ++i)
              for (std::size_t j = 0; j < local_size; ++j)
                mass_matrix[i][j] += phi[i] * phi[j] * factor;
          }
        mass_matrix.invert();
        return mass_matrix;
      }

    private:

      int _intorder;
      std::unordered_map<std::size_t,std::vector<MassMatrix>> _inverse_mass_matrices;
    };


#ifndef DOXYGEN
    namespace Impl::Adapt {

      /* Store the local coefficients of x for every interior leaf element in
       * transfer_map (keyed by the local id of the element), and for elements
       * that might vanish during adaptation the local L2 projection of x onto
       * all ancestors that might vanish.
       */
      template<typename Space, typename X, typename TransferMap, typename Projection>
      void backup(const Space& space, Projection& projection, const X& x, TransferMap& transfer_map)
      {
        using RF = typename X::ElementType;
        using Dune::TypeTree::child;

        const auto& g = grid(space);
        const auto& id_set = g.localIdSet();
        const int max_level = g.maxLevel();
        const int int_order = projection.integrationOrder();

        LocalCoefficients<Space> local(space);
        ElementTree<Space> coarse(space);
        std::vector<RF> u_fine;

        for (const auto& element : elementRange(space))
          {
            local.bind(element);
            local.read(x,transfer_map[id_set.id(element)]);

            auto ancestor = element;
            while (ancestor.mightVanish())
              {
                // work around UG bug!
                if (!ancestor.hasFather())
                  break;

                ancestor = ancestor.father();

                auto& u_coarse = transfer_map[id_set.id(ancestor)];
                // don't project more than once
                if (u_coarse.size() > 0)
                  continue;

                coarse.bind(ancestor);
                u_coarse.assign(coarse.tree().size(),RF(0));
                const auto coarse_geometry = ancestor.geometry();
                constexpr int dim = std::decay_t<decltype(element)>::dimension;
                using DF = typename std::decay_t<decltype(element)>::Geometry::ctype;

                for (const auto& descendant : descendantElements(ancestor,max_level))
                  {
                    // only evaluate on entities with data
                    if (!descendant.isLeaf())
                      continue;
                    local.bind(descendant);
                    local.read(x,u_fine);
                    const auto fine_geometry = descendant.geometry();

                    // project leaf by leaf of the space tree
                    std::size_t leaf_index = 0;
                    Dune::TypeTree::forEachLeafNode(coarse.tree(), [&](const auto& coarse_leaf, auto tree_path) {
                      const auto& fine_leaf = child(local.tree(),tree_path);
                      const auto& fine_fe = fine_leaf.finiteElement();
                      const auto& coarse_fe = coarse.finiteElement(coarse_leaf);
                      using Range = typename std::decay_t<decltype(coarse_fe)>::Traits::LocalBasisType::Traits::RangeType;
                      const auto& inverse_mass_matrix = projection.inverseMassMatrix(ancestor,leaf_index,coarse_fe);

                      std::vector<Range> fine_phi, coarse_phi;
                      for (const auto& ip : QuadratureRules<DF,dim>::rule(descendant.type(),int_order))
                        {
                          const auto coarse_local = coarse_geometry.local(fine_geometry.global(ip.position()));
                          fine_fe.localBasis().evaluateFunction(ip.position(),fine_phi);
                          coarse_fe.localBasis().evaluateFunction(coarse_local,coarse_phi);
                          const DF factor = ip.weight()
                            * fine_geometry.integrationElement(ip.position())
                            / coarse_geometry.integrationElement(coarse_local);

                          auto val = Range{0.0};
                          for (std::size_t i = 0; i < fine_phi.size(); ++i)
                            val.axpy(u_fine[fine_leaf.localIndex(i)],fine_phi[i]);

                          assert(inverse_mass_matrix.M()==coarse_phi.size());
                          for (std::size_t i = 0; i < coarse_phi.size(); ++i)
                            {
                              auto y = Range{0.0};
                              for (std::size_t j = 0; j < inverse_mass_matrix.M(); ++j)
                                y.axpy(inverse_mass_matrix[i][j],coarse_phi[j]);
                              u_coarse[coarse_leaf.localIndex(i)] += factor * (y * val);
                            }
                        }
                      ++leaf_index;
                    });
                  }
              }
          }
      }

      //! The function of a coarse element, evaluated in local coordinates of a fine element
      template<typename FiniteElement, typename Geometry, typename LocalDOFVector, typename LocalIndex>
      struct CoarseFunction
      {
        using Range = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
        using Domain = typename FiniteElement::Traits::LocalBasisType::Traits::DomainType;

        //! Traits class containing decayed types
        struct Traits
        {
          using RangeType = std::decay_t<Range>;
        };

        Range operator()(const Domain& x) const
        {
          _phi.resize(_finite_element.localBasis().size());
          _finite_element.localBasis().evaluateFunction(_coarse_geometry.local(_fine_geometry.global(x)),_phi);
          Range y = 0;
          for (std::size_t i = 0; i < _phi.size(); ++i)
            y.axpy(_dofs[_local_index(i)],_phi[i]);
          return y;
        }

        const FiniteElement& _finite_element;
        Geometry _coarse_geometry;
        Geometry _fine_geometry;
        const LocalDOFVector& _dofs;
        LocalIndex _local_index;
        mutable std::vector<Range> _phi;
      };

      /* Rebuild x after adaptation from transfer_map: for every interior leaf
       * element the stored coefficients are copied (element unchanged) or the
       * function of the nearest stored ancestor is interpolated. DOFs shared
       * by several elements are averaged.
       */
      template<typename Space, typename X, typename TransferMap>
      void replay(const Space& space, X& x, const TransferMap& transfer_map)
      {
        using RF = typename X::ElementType;
        using Dune::TypeTree::child;
        using CountVector = Backend::Vector<Space,int>;
        CountVector uc(space,0);

        const auto& id_set = grid(space).localIdSet();
        LocalCoefficients<Space> local(space);
        ElementTree<Space> coarse(space);
        std::vector<RF> u_fine, u_tmp;
        std::vector<int> counts;

        for (const auto& element : elementRange(space))
          {
            auto ancestor = element;
            typename TransferMap::const_iterator map_it;
            while ((map_it = transfer_map.find(id_set.id(ancestor))) == transfer_map.end())
              {
                if (!ancestor.hasFather())
                  DUNE_THROW(Exception,
                             "transferMap of GridAdaptor didn't contain ancestor of element with id " << id_set.id(ancestor));
                ancestor = ancestor.father();
              }
            const auto& u_coarse = map_it->second;

            local.bind(element);
            // test identity using ids
            if (id_set.id(element) == id_set.id(ancestor))
              {
                // no interpolation necessary, just copy the saved data
                local.add(x,u_coarse);
              }
            else
              {
                coarse.bind(ancestor);
                u_fine.assign(local.size(),RF(0));
                Dune::TypeTree::forEachLeafNode(coarse.tree(), [&](const auto& coarse_leaf, auto tree_path) {
                  const auto& fine_leaf = child(local.tree(),tree_path);
                  const auto& coarse_fe = coarse.finiteElement(coarse_leaf);
                  auto local_index = [&](std::size_t i) { return coarse_leaf.localIndex(i); };
                  CoarseFunction<std::decay_t<decltype(coarse_fe)>,typename std::decay_t<decltype(element)>::Geometry,
                                 std::vector<RF>,decltype(local_index)>
                    f{coarse_fe,ancestor.geometry(),element.geometry(),u_coarse,local_index,{}};
                  const auto& fine_fe = fine_leaf.finiteElement();
                  u_tmp.assign(fine_fe.localBasis().size(),RF(0));
                  fine_fe.localInterpolation().interpolate(f,u_tmp);
                  for (std::size_t i = 0; i < u_tmp.size(); ++i)
                    u_fine[fine_leaf.localIndex(i)] = u_tmp[i];
                });
                local.add(x,u_fine);
              }

            counts.assign(local.size(),1);
            local.add(uc,counts);
          }

        if constexpr (not Impl::DuneFunctionsBasis<Space>)
          {
            typedef Dune::PDELab::AddDataHandle<Space,X> DOFHandle;
            DOFHandle addHandle1(space,x);
            space.entitySet().gridView().communicate(addHandle1,
                                                     Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
            typedef Dune::PDELab::AddDataHandle<Space,CountVector> CountHandle;
            CountHandle addHandle2(space,uc);
            space.entitySet().gridView().communicate(addHandle2,
                                                     Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
          }
        else if (space.gridView().comm().size() > 1)
          DUNE_THROW(NotImplemented,"Grid adaptation with dune-functions bases is sequential only");

        // normalize multiple-interpolated DOFs by taking the arithmetic average
        typename CountVector::iterator ucit = uc.begin();
        for (typename X::iterator uit = x.begin(), uend = x.end(); uit != uend; ++uit, ++ucit)
          (*uit) /= ((*ucit) > 0 ? (*ucit) : 1.0);
      }

    } // namespace Impl::Adapt
#endif // DOXYGEN


    /*! @class GridAdaptor
     *
     * @brief Class for automatic adaptation of the grid.
     *
     *        The GridAdaptor capsules the act of deciding which Elems to refine and coarsen,
     *        adapting the grid, and transfering the solution from the old grid to the new one.
     *        Works for grid function spaces and dune-functions bases (sequential only for
     *        bases); trees of spaces are transferred leaf by leaf.
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFSU       Type of ansatz space, we need to update it after adaptation
     * @tparam U          Container class of the solution
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFSU, class U, class Projection>
    class GridAdaptor
    {
      typedef typename Grid::LocalIdSet IDSet;
      typedef typename IDSet::IdType ID;

    public:
      typedef std::unordered_map<ID,std::vector<typename U::ElementType> > MapType;


      /*! @brief The constructor.
       *
       * @param gfs         The ansatz space, we need to update it
       */
      explicit GridAdaptor(const GFSU& gfs)
      {}

      /* @brief Save the solution before adaptation
       *
       * @param[in]  u           The solution that will be saved
       * @param[out] transferMap The map containing the solution during adaptation
       */
      void backupData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, MapType& transfer_map)
      {
        Impl::Adapt::backup(gfsu,projection,u,transfer_map);
      }

      /* @brief Restore the solution after adaptation
       *
       * @param[out] u           The solution after adaptation
       * @param[in]  transferMap The map that contains the information for the rebuild of u
       */
      void replayData(Grid& grid, GFSU& gfsu, Projection& projection, U& u, const MapType& transfer_map)
      {
        Impl::Adapt::replay(gfsu,u,transfer_map);
      }

    };

    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space (grid function space or dune-functions basis),
     *                    we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     *
     * @param grid        The grid we want to adapt
     * @param gfs         The ansatz space; a dune-functions basis is updated with
     *                    the new leaf grid view (sequential only)
     * @param x1          The DOF container
     * @param int_order   The integration order used in inversion of the mass matrix.
     *                    For scalar problems 2*k with k being the polynomial degree
     *                    of ansatz and test space is a suitable choice. For systems,
     *                    the same order is currently used for all components.
     */
    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, int int_order)
    {
      typedef L2Projection<GFS,X> Projection;
      Projection projection(gfs,int_order);

      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor(gfs);

      // prepare the grid for refinement
      grid.preAdapt();

      // save u
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);

      // adapt the grid
      grid.adapt();

      // update the function spaces
      Impl::Adapt::update(gfs,grid);

      // reset u
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);

      // clean up
      grid.postAdapt();
    }

    /*! grid adaptation as a function
     *
     * @brief adapt a grid, corresponding function space and solution vectors
     *
     * Assumes that the grid's elements have been marked for refinement and coarsening appropriately before
     *
     * @tparam Grid       Type of the grid we want to adapt
     * @tparam GFS        Type of ansatz space (grid function space or dune-functions basis),
     *                    we need to update it after adaptation
     * @tparam X          Container class for DOF vectors
     * @tparam Projection Projection used when Elems vanish
     */
    template<class Grid, class GFS, class X>
    void adapt_grid (Grid& grid, GFS& gfs, X& x1, X& x2, int int_order)
    {
      typedef L2Projection<GFS,X> Projection;
      Projection projection(gfs,int_order);

      GridAdaptor<Grid,GFS,X,Projection> grid_adaptor(gfs);

      // prepare the grid for refinement
      grid.preAdapt();

      // save solution
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap1;
      grid_adaptor.backupData(grid,gfs,projection,x1,transferMap1);
      typename GridAdaptor<Grid,GFS,X,Projection>::MapType transferMap2;
      grid_adaptor.backupData(grid,gfs,projection,x2,transferMap2);

      // adapt the grid
      grid.adapt();

      // update the function spaces
      Impl::Adapt::update(gfs,grid);

      // interpolate solution
      x1 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x1,transferMap1);
      x2 = X(gfs,0.0);
      grid_adaptor.replayData(grid,gfs,projection,x2,transferMap2);

      // clean up
      grid.postAdapt();
    }

#ifndef DOXYGEN
    namespace impl{

      // Struct for storing a GridFunctionSpace, corrosponding vectors and integration order
      template <typename G, typename... X>
      struct GFSWithVectors
      {
        // Export types
        using GFS = G;
        using Tuple = std::tuple<X&...>;

        GFSWithVectors (GFS& gfs, int integrationOrder, X&... x) :
          _gfs(gfs),
          _integrationOrder(integrationOrder),
          _tuple(x...)
        {}

        GFS& _gfs;
        int _integrationOrder;
        Tuple _tuple;
      };

      // Forward declarations needed for the recursion
      template <typename Grid>
      void iteratePacks(Grid& grid);
      template <typename Grid, typename X, typename... XS>
      void iteratePacks(Grid& grid, X& x, XS&... xs);

      // This function is called after the last vector of the tuple.  Here
      // the next pack is called.  On the way back we update the current
      // function space.
      template<std::size_t I = 0, typename Grid, typename X, typename... XS>
      inline typename std::enable_if<I == std::tuple_size<typename X::Tuple>::value, void>::type
      iterateTuple(Grid& grid, X& x, XS&... xs)
      {
        // Iterate next pack
        iteratePacks(grid,xs...);

        // On our way back we need to update the current function space
        Impl::Adapt::update(x._gfs,grid);
      }

      /* In this function we store the data of the current vector (indicated
       * by template parameter I) of the current pack. After recursively
       * iterating through the other packs and vectors we replay the data.
       *
       * @tparam I      std:size_t used for tmp
       * @tparam Grid   Grid type
       * @tparam X      Current  pack
       * @tparam ...XS  Remaining packs
       */
      template<std::size_t I = 0, typename Grid, typename X, typename... XS>
      inline typename std::enable_if<I < std::tuple_size<typename X::Tuple>::value, void>::type
      iterateTuple(Grid& grid, X& x, XS&... xs)
      {
        // Get some basic types
        using GFS = typename X::GFS;
        using Tuple = typename X::Tuple;
        using V  = typename std::decay<typename std::tuple_element<I,Tuple>::type>::type;
        // // alternative:
        // auto v = std::get<I>(x._tuple);
        // using V = decltype(v);

        // Setup classes for data restoring
        typedef Dune::PDELab::L2Projection <GFS,V> Projection;
        Projection projection(x._gfs,x._integrationOrder);
        GridAdaptor<Grid,GFS,V,Projection> gridAdaptor(x._gfs);

        // Store vector data
        typename GridAdaptor<Grid,GFS,V,Projection>::MapType transferMap;
        gridAdaptor.backupData(grid,x._gfs,projection,std::get<I>(x._tuple),transferMap);

        // Recursively iterate through remaining vectors (and packs). Grid
        // adaption will be done at the end of recursion.
        iterateTuple<I + 1, Grid, X, XS...>(grid,x,xs...);

        // Play back data. Note: At this point the function space was
        // already updatet.
        std::get<I>(x._tuple) = V(x._gfs,0.0);
        gridAdaptor.replayData(grid,x._gfs,projection,std::get<I>(x._tuple),transferMap);
      }

      // This gets called after the last pack.  After this function call we
      // have visited every vector of every pack and we will go back through
      // the recursive function calls.
      template <typename Grid>
      void iteratePacks(Grid& grid)
      {
        // Adapt the grid
        grid.adapt();
      }

      /* Use template meta programming to iterate over packs at compile time
       *
       * In order to adapt our grid and all vectors of all packs we need to
       * do the following:
       * - Iterate over all vectors of all packs.
       * - Store the data from the vectors where things could change.
       * - Adapt our grid.
       * - Update function spaces and restore data.
       *
       * The key point is that we need the object that stores the data to
       * replay it.  Because of that we can not just iterate over the packs
       * and within each pack iterate over the vectors but we have to make
       * one big recursion.  Therefore we iterate over the vectors of the
       * current pack.
       */
       template <typename Grid, typename X, typename... XS>
       void iteratePacks(Grid& grid, X& x, XS&... xs)
       {
         iterateTuple(grid,x,xs...);
       }

    } // namespace impl
#endif // DOXYGEN

    /*! \brief Pack function space and vectors for grid adaption
     *
     * This function packs a GridFunctionSpace an integration order and an
     * arbitrary number of vectors in a single struct.
     *
     * Important: You have to make sure that all vectors belong to the
     * same function space.
     *
     * @tparam GFS   Grid function space
     * @tparam ...X  Arbitrary number of corresponding vectors
     */
    template <typename GFS, typename... X>
    impl::GFSWithVectors<GFS,X...> transferSolutions(GFS& gfs, int integrationOrder, X&... x)
    {
      impl::GFSWithVectors<GFS,X...> gfsWithVectors(gfs, integrationOrder, x...);
      return gfsWithVectors;
    }

    /*! \brief Adapt grid and multiple function spaces with corresponding vectors
     *
     * Assumes that the grid's elements have been marked for refinement and
     * coarsement appropriately befor
     *
     * @tparam Grid   Type of the Grid
     * @tparam X      Packed GFS with vectors that should be adapted
     *
     * Note: A pack can be created using the transferSolution function.
     */
    template <typename Grid, typename... X>
    void adaptGrid(Grid& grid, X&... x)
    {
      // Prepare the grid for refinement
      grid.preAdapt();

      // Iterate over packs
      impl::iteratePacks(grid,x...);

      // Clean up
      grid.postAdapt();
    }


    template<typename T>
    void error_fraction(const T& x, typename T::ElementType alpha, typename T::ElementType beta,
                        typename T::ElementType& eta_alpha, typename T::ElementType& eta_beta, int verbose=0)
    {
      if (verbose>0)
        std::cout << "+++ error fraction: alpha=" << alpha << " beta=" << beta << std::endl;
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error = x.one_norm();
      NumberType max_error = x.infinity_norm();
      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;
          for (const auto& error : x)
            {
              if (error >=eta_alpha) { sum_alpha += error; alpha_count++;}
              if (error < eta_beta) { sum_beta += error; beta_count++;}
            }
          if (verbose>1)
            {
              std::cout << "+++ " << j << " eta_alpha=" << eta_alpha << " alpha_fraction=" << sum_alpha/total_error
                        << " elements: " << alpha_count << " of " << x.N() << std::endl;
              std::cout << "+++ " << j << " eta_beta=" << eta_beta << " beta_fraction=" << sum_beta/total_error
                        << " elements: " << beta_count << " of " << x.N() << std::endl;
            }
          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          if (sum_alpha>alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          if (sum_beta>beta*total_error)
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;
        }
      if (verbose>0)
        {
          std::cout << "+++ refine_threshold=" << eta_alpha
                    << " coarsen_threshold=" << eta_beta << std::endl;
        }
    }


    template<typename T>
    void element_fraction(const T& x, typename T::ElementType alpha, typename T::ElementType beta,
                          typename T::ElementType& eta_alpha, typename T::ElementType& eta_beta, int verbose=0)
    {
      const int steps=20; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error =x.N();
      NumberType max_error = x.infinity_norm();
      NumberType eta_alpha_left = 0.0;
      NumberType eta_alpha_right = max_error;
      NumberType eta_beta_left = 0.0;
      NumberType eta_beta_right = max_error;
      for (int j=1; j<=steps; j++)
        {
          eta_alpha = 0.5*(eta_alpha_left+eta_alpha_right);
          eta_beta = 0.5*(eta_beta_left+eta_beta_right);
          NumberType sum_alpha=0.0;
          NumberType sum_beta=0.0;
          unsigned int alpha_count = 0;
          unsigned int beta_count = 0;

          for (const auto& error : x)
            {
              if (error>=eta_alpha) { sum_alpha += 1.0; alpha_count++;}
              if (error< eta_beta) { sum_beta +=1.0; beta_count++;}
            }
          if (verbose>1)
            {
              std::cout << j << " eta_alpha=" << eta_alpha << " alpha_fraction=" << sum_alpha/total_error
                        << " elements: " << alpha_count << " of " << x.N() << std::endl;
              std::cout << j << " eta_beta=" << eta_beta << " beta_fraction=" << sum_beta/total_error
                        << " elements: " << beta_count << " of " << x.N() << std::endl;
            }
          if (std::abs(alpha-sum_alpha/total_error) <= 0.01 && std::abs(beta-sum_beta/total_error) <= 0.01) break;
          if (sum_alpha>alpha*total_error)
            eta_alpha_left = eta_alpha;
          else
            eta_alpha_right = eta_alpha;
          if (sum_beta>beta*total_error)
            eta_beta_right = eta_beta;
          else
            eta_beta_left = eta_beta;
        }
      if (verbose>0)
        {
          std::cout << "+++ refine_threshold=" << eta_alpha
                    << " coarsen_threshold=" << eta_beta << std::endl;
        }
    }

    /** Compute error distribution
     */
    template<typename T>
    void error_distribution(const T& x, unsigned int bins)
    {
      const int steps=30; // max number of bisection steps
      typedef typename T::ElementType NumberType;
      NumberType total_error = x.one_norm();
      NumberType total_elements = x.N();
      NumberType max_error = x.infinity_norm();
      std::vector<NumberType> left(bins,0.0);
      std::vector<NumberType> right(bins,max_error*(1.0+1e-8));
      std::vector<NumberType> eta(bins);
      std::vector<NumberType> target(bins);
      for (unsigned int k=0; k<bins; k++)
        target[k]= (k+1)/((NumberType)bins);
      for (int j=1; j<=steps; j++)
        {
          for (unsigned int k=0; k<bins; k++)
            eta[k]= 0.5*(left[k]+right[k]);
          std::vector<NumberType> sum(bins,0.0);
          std::vector<int> count(bins,0);

          for (typename T::const_iterator it = x.begin(),
                 end = x.end();
               it != end;
               ++it)
            {
              for (unsigned int k=0; k<bins; k++)
                if (*it<=eta[k])
                  {
                    sum[k] += *it;
                    count[k] += 1;
                  }
            }
          // std::cout << std::endl;
          // std::cout << "// step " << j << std::endl;
          // for (unsigned int k=0; k<bins; k++)
          //    std::cout << k+1 << " " << count[k] << " " << eta[k] << " " << right[k]-left[k]
          //          << " " << sum[k]/total_error << " " << target[k] << std::endl;
          for (unsigned int k=0; k<bins; k++)
            if (sum[k]<=target[k]*total_error)
              left[k] = eta[k];
            else
              right[k] = eta[k];
        }
      std::vector<NumberType> sum(bins,0.0);
      std::vector<int> count(bins,0);
      for (unsigned int i=0; i<x.N(); i++)
        for (unsigned int k=0; k<bins; k++)
          if (x[i]<=eta[k])
            {
              sum[k] += x[i];
              count[k] += 1;
            }
      std::cout << "+++ error distribution" << std::endl;
      std::cout << "+++ number of elements: " << x.N() << std::endl;
      std::cout << "+++ max element error:  " << max_error << std::endl;
      std::cout << "+++ total error:        " << total_error << std::endl;
      std::cout << "+++ bin #elements eta sum/total " << std::endl;
      for (unsigned int k=0; k<bins; k++)
        std::cout << "+++ " << k+1 << " " << count[k] << " " << eta[k] << " " << sum[k]/total_error << std::endl;
    }

    template<typename Grid, typename X>
    void mark_grid (Grid &grid, const X& x, typename X::ElementType refine_threshold,
                    typename X::ElementType coarsen_threshold, int min_level = 0, int max_level = std::numeric_limits<int>::max(), int verbose=0)
    {
      typedef typename Grid::LeafGridView GV;

      GV gv = grid.leafGridView();

      unsigned int refine_cnt=0;
      unsigned int coarsen_cnt=0;

      typedef typename X::GridFunctionSpace GFS;
      Impl::Adapt::LocalCoefficients<GFS> local(x.gridFunctionSpace());
      std::vector<typename X::ElementType> x_view;

      for(const auto& cell : elements(gv))
        {
          local.bind(cell);
          local.read(x,x_view);

          if (x_view[0]>=refine_threshold && cell.level() < max_level)
            {
              grid.mark(1,cell);
              refine_cnt++;
            }
          if (x_view[0]<=coarsen_threshold && cell.level() > min_level)
            {
              grid.mark(-1,cell);
              coarsen_cnt++;
            }
        }
      if (verbose>0)
        std::cout << "+++ mark_grid: " << refine_cnt << " marked for refinement, "
                  << coarsen_cnt << " marked for coarsening" << std::endl;
    }


    template<typename Grid, typename X>
    void mark_grid_for_coarsening (Grid &grid, const X& x, typename X::ElementType refine_threshold,
                                   typename X::ElementType coarsen_threshold, int verbose=0)
    {
      typedef typename Grid::LeafGridView GV;

      GV gv = grid.leafGridView();

      unsigned int coarsen_cnt=0;

      typedef typename X::GridFunctionSpace GFS;
      Impl::Adapt::LocalCoefficients<GFS> local(x.gridFunctionSpace());
      std::vector<typename X::ElementType> x_view;

      for(const auto& cell : elements(gv))
        {
          local.bind(cell);
          local.read(x,x_view);

          if (x_view[0]>=refine_threshold)
            {
              grid.mark(-1,cell);
              coarsen_cnt++;
            }
          if (x_view[0]<=coarsen_threshold)
            {
              grid.mark(-1,cell);
              coarsen_cnt++;
            }
        }
      if (verbose>0)
        std::cout << "+++ mark_grid_for_coarsening: "
                  << coarsen_cnt << " marked for coarsening" << std::endl;
    }


    class TimeAdaptationStrategy
    {
      // strategy parameters
      double scaling;
      double optimistic_factor;
      double coarsen_limit;
      double balance_limit;
      double tol;
      double T;
      int verbose;
      bool no_adapt;
      double refine_fraction_while_refinement;
      double coarsen_fraction_while_refinement;
      double coarsen_fraction_while_coarsening;
      double timestep_decrease_factor;
      double timestep_increase_factor;

      // results to be reported to the user after evaluating the error
      bool accept;
      bool adapt_dt;
      bool adapt_grid;
      double newdt;
      double q_s, q_t;

      // state variables
      bool have_decreased_time_step;
      bool have_refined_grid;

      // the only state variable: accumulated error
      double accumulated_estimated_error_squared;
      double minenergy_rate;

    public:
      TimeAdaptationStrategy (double tol_, double T_, int verbose_=0)
        : scaling(16.0), optimistic_factor(1.0), coarsen_limit(0.5), balance_limit(0.33333),
          tol(tol_), T(T_), verbose(verbose_), no_adapt(false),
          refine_fraction_while_refinement(0.7),
          coarsen_fraction_while_refinement(0.2),
          coarsen_fraction_while_coarsening(0.2),
          timestep_decrease_factor(0.5), timestep_increase_factor(1.5),
          accept(false), adapt_dt(false), adapt_grid(false), newdt(1.0),
          have_decreased_time_step(false), have_refined_grid(false),
          accumulated_estimated_error_squared(0.0),
          minenergy_rate(0.0)
      {
      }

      void setTimeStepDecreaseFactor (double s)
      {
        timestep_decrease_factor=s;
      }

      void setTimeStepIncreaseFactor (double s)
      {
        timestep_increase_factor=s;
      }

      void setRefineFractionWhileRefinement (double s)
      {
        refine_fraction_while_refinement=s;
      }

      void setCoarsenFractionWhileRefinement (double s)
      {
        coarsen_fraction_while_refinement=s;
      }

      void setCoarsenFractionWhileCoarsening (double s)
      {
        coarsen_fraction_while_coarsening=s;
      }

      void setMinEnergyRate (double s)
      {
        minenergy_rate=s;
      }

      void setCoarsenLimit (double s)
      {
        coarsen_limit=s;
      }

      void setBalanceLimit (double s)
      {
        balance_limit=s;
      }

      void setTemporalScaling (double s)
      {
        scaling=s;
      }

      void setOptimisticFactor (double s)
      {
        optimistic_factor=s;
      }

      void setAdaptationOn ()
      {
        no_adapt = false;
      }

      void setAdaptationOff ()
      {
        no_adapt = true;
      }

      bool acceptTimeStep () const
      {
        return accept;
      }

      bool adaptDT () const
      {
        return adapt_dt;
      }

      bool adaptGrid () const
      {
        return adapt_grid;
      }

      double newDT () const
      {
        return newdt;
      }

      double qs () const
      {
        return q_s;
      }

      double qt () const
      {
        return q_t;
      }

      double endT() const
      {
        return T;
      }

      double accumulatedErrorSquared () const
      {
        return accumulated_estimated_error_squared;
      }

      // to be called when new time step is done
      void startTimeStep ()
      {
        have_decreased_time_step = false;
        have_refined_grid = false;
      }

      template<typename GM, typename X>
      void evaluate_estimators (GM& grid, double time, double dt, const X& eta_space,  const X& eta_time, double energy_timeslab)
      {
        using std::sqrt;
        accept=false;
        adapt_dt=false;
        adapt_grid=false;
        newdt=dt;

        double spatial_error = eta_space.one_norm();
        double temporal_error = scaling*eta_time.one_norm();
        double sum = spatial_error + temporal_error;
        //double allowed = optimistic_factor*(tol*tol-accumulated_estimated_error_squared)*dt/(T-time);
        double allowed = tol*tol*(energy_timeslab+minenergy_rate*dt);
        q_s = spatial_error/sum;
        q_t = temporal_error/sum;

        // print some statistics
        if (verbose>0)
          std::cout << "+++"
                    << " q_s=" << q_s
                    << " q_t=" << q_t
                    << " sum/allowed=" << sum/allowed
                    // << " allowed=" << allowed
                    << " estimated error=" << sqrt(accumulated_estimated_error_squared+sum)
                    << " energy_rate=" << energy_timeslab/dt
                    << std::endl;

        // for simplicity: a mode that does no adaptation at all
        if (no_adapt)
          {
            accept = true;
            accumulated_estimated_error_squared += sum;
            if (verbose>1) std::cout << "+++ no adapt mode" << std::endl;
            return;
          }

        // the adaptation strategy
        if (sum<=allowed)
          {
            // we will accept this time step
            accept = true;
            if (verbose>1) std::cout << "+++ accepting time step" << std::endl;
            accumulated_estimated_error_squared += sum;

            // check if grid size or time step needs to be adapted
            if (sum<coarsen_limit*allowed)
              {
                // the error is too small, i.e. the computation is inefficient
                if (q_t<balance_limit)
                  {
                    // spatial error is dominating => increase time step
                    newdt = timestep_increase_factor*dt;
                    adapt_dt = true;
                    if (verbose>1) std::cout << "+++ spatial error dominates: increase time step" << std::endl;
                  }
                else
                  {
                    if (q_s>balance_limit)
                      {
                        // step sizes balanced: coarsen in time
                        newdt = timestep_increase_factor*dt;
                        adapt_dt = true;
                        if (verbose>1) std::cout << "+++ increasing time step" << std::endl;
                      }
                    // coarsen grid in space
                    double eta_refine, eta_coarsen;
                    if (verbose>1) std::cout << "+++ mark grid for coarsening" << std::endl;
                    //error_distribution(eta_space,20);
                    Dune::PDELab::error_fraction(eta_space,coarsen_fraction_while_coarsening,
                                                 coarsen_fraction_while_coarsening,eta_refine,eta_coarsen);
                    Dune::PDELab::mark_grid_for_coarsening(grid,eta_space,eta_refine,eta_coarsen,verbose);
                    adapt_grid = true;
                  }
              }
            else
              {
                // modify at least the time step
                if (q_t<balance_limit)
                  {
                    // spatial error is dominating => increase time step
                    newdt = timestep_increase_factor*dt;
                    adapt_dt = true;
                    if (verbose>1) std::cout << "+++ spatial error dominates: increase time step" << std::endl;
                  }
              }
          }
        else
          {
            // error is too large, we need to do something
            if (verbose>1) std::cout << "+++ will redo time step" << std::endl;
            if (q_t>1-balance_limit)
              {
                // temporal error is dominating => decrease time step only
                newdt = timestep_decrease_factor*dt;
                adapt_dt = true;
                have_decreased_time_step = true;
                if (verbose>1) std::cout << "+++ decreasing time step only" << std::endl;
              }
            else
              {
                if (q_t<balance_limit)
                  {
                    if (!have_decreased_time_step)
                      {
                        // time steps size not balanced (too small)
                        newdt = timestep_increase_factor*dt;
                        adapt_dt = true;
                        if (verbose>1) std::cout << "+++ increasing time step" << std::endl;
                      }
                  }
                else
                  {
                    // step sizes balanced: refine in time as well
                    newdt = timestep_decrease_factor*dt;
                    adapt_dt = true;
                    have_decreased_time_step = true;
                    if (verbose>1) std::cout << "+++ decreasing time step" << std::endl;
                  }
                // refine grid in space
                double eta_refine, eta_coarsen;
                if (verbose>1) std::cout << "+++ BINGO mark grid for refinement and coarsening" << std::endl;
                //error_distribution(eta_space,20);
                Dune::PDELab::error_fraction(eta_space,refine_fraction_while_refinement,
                                             coarsen_fraction_while_refinement,eta_refine,eta_coarsen,0);
                Dune::PDELab::mark_grid(grid,eta_space,eta_refine,eta_coarsen,verbose);
                adapt_grid = true;
              }
          }
      }
    };



  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ADAPTIVITY_ADAPTIVITY_HH
