# PDELab

@tableofcontents

[PDELab][pdelab] is a discretization module for the [DUNE framework][dune]
@cite BastianEtAl2008a @cite BastianEtAl2008b @cite BastianEtAl2021 @cite BlattEtAl2025
for solving partial differential equations (PDEs) @cite BastianHeimannMarnach2010.
It provides high-level abstractions that allow for quick development of
simulation codes, while keeping the flexibility to implement new numerical
methods and to optimize for specific simulation scenarios.
All components are generic: the same code works for different grids,
dimensions, polynomial degrees, and number types, and in sequential as well as
in parallel runs.

@button{@ref pdelab-getting-started "Getting started"}
@button{@ref pdelab_recipes "Recipes"}
@button{[Tutorials][tutorials]}

## What PDELab offers {#pdelab-features}

A PDELab simulation is typically composed of the following building blocks.

### Function spaces
A @ref Dune::PDELab::GridFunctionSpace "GridFunctionSpace" combines a grid view
with a [finite element map](@ref FiniteElementMap), which provides the local
finite element on each grid element. Finite element maps are available, e.g., for
conforming Lagrange elements on simplices and cubes (Pk, Qk),
discontinuous Galerkin elements (Lagrange, Legendre, Lobatto, monomial, and
orthonormal bases), finite volumes (P0), Raviart-Thomas and
Brezzi-Douglas-Marini elements, lowest-order Nédélec (Whitney) edge elements
on simplices, and Rannacher-Turek elements.
Systems of PDEs are described by combining spaces into
@ref Dune::PDELab::PowerGridFunctionSpace "PowerGridFunctionSpace" and
@ref Dune::PDELab::CompositeGridFunctionSpace "CompositeGridFunctionSpace",
e.g., for Taylor-Hood elements.
The [ordering](@ref Ordering) of the degrees of freedom and the blocking
of the resulting vectors and matrices can be chosen freely.
See the [Grid Function Space](@ref GridFunctionSpace) topic for details.

### Local operators
The weak formulation of a PDE is implemented as a *local operator*, which
computes the residual and its Jacobian on single grid elements and faces
(methods `alpha_volume`, `alpha_skeleton`, `alpha_boundary`, `lambda_volume`, ...).
Local operators only need to implement the residual: Jacobians can be
computed by numerical differentiation.
PDELab comes with a library of local operators, e.g., for
* convection-diffusion-reaction equations with conforming finite elements, discontinuous Galerkin methods, and cell-centered finite volumes,
* linear elasticity,
* Stokes and Navier-Stokes equations with Taylor-Hood and discontinuous Galerkin discretizations,
* Maxwell's equations and linear acoustics with discontinuous Galerkin methods,
* two-phase flow in porous media with finite volumes.

See the [Local Operator](@ref LocalOperator) topic for details.

### Assembly
A @ref Dune::PDELab::GridOperator "GridOperator" assembles the global residual
and Jacobian from a local operator, a trial and a test space.
It also provides matrix-free application of the Jacobian, and a fast assembler
for discontinuous Galerkin methods.
See the @ref GridOperatorDocPage "Grid Operator" page for an overview.

### Constraints
Constraints restrict the function spaces, e.g., to incorporate
essential boundary conditions (Dirichlet constraints), hanging nodes on
nonconforming grids, and the degrees of freedom of other processes in parallel runs.

### Solvers
PDELab provides solver backends based on [dune-istl][istl]: iterative solvers
with various preconditioners including algebraic multigrid (AMG), and direct
solvers, for sequential runs as well as for overlapping and nonoverlapping
domain decompositions in parallel runs.
On top of these, PDELab provides solvers for linear stationary problems
(@ref Dune::PDELab::StationaryLinearProblemSolver "StationaryLinearProblemSolver")
and a Newton method with line search for nonlinear problems
(@ref Dune::PDELab::NewtonMethod "NewtonMethod").

### Time stepping
Instationary problems are solved with the method of lines.
Implicit and explicit [one-step methods](@ref OneStepMethod) are available,
including implicit and explicit Euler, the one-step θ-scheme (e.g., Crank-Nicolson),
Alexander's diagonally implicit Runge-Kutta schemes, the fractional-step θ-scheme, Heun's method, a third-order
strong stability preserving method by Shu and Osher, and the classical
fourth-order Runge-Kutta method.

### Adaptivity, parallelism, and output
Grids can be adapted based on error estimators, with the solution transferred
to the new grid (@ref Dune::PDELab::adapt_grid "adapt_grid").
Parallel computations use the domain decomposition of the grid, with
communication handled transparently (see the [Entity Set Partitions](@ref Partition) topic).
Solutions can be written for visualization, e.g., in ParaView, with
@ref Dune::PDELab::addSolutionToVTKWriter "addSolutionToVTKWriter".

## A first look {#pdelab-first-look}

The following code solves the Laplace equation with Q2 finite elements on
a structured grid of the unit square and writes the solution to a VTK file.
The weak formulation is provided by a predefined local operator for
convection-diffusion-reaction problems.

```cpp
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/pdelab.hh>

using namespace Dune::PDELab;

// A structured grid of the unit square with 32x32 cells
Dune::YaspGrid<2> grid({1.0,1.0}, {32,32});
auto gv = grid.leafGridView();
using GV = decltype(gv);

// The PDE: a convection-diffusion-reaction model problem
// (here the Laplace equation with Dirichlet data g(x) = |x|)
using Problem = ConvectionDiffusionModelProblem<GV,double>;
Problem problem;
ConvectionDiffusionBoundaryConditionAdapter bctype(problem);

// Function space: Q2 finite elements with Dirichlet constraints
using FEM = QkLocalFiniteElementMap<GV,double,double,2>;
FEM fem(gv);
using GFS = GridFunctionSpace<GV,FEM,ConformingDirichletConstraints,
                              ISTL::VectorBackend<>>;
GFS gfs(gv,fem);
gfs.name("u");
using CC = GFS::ConstraintsContainer<double>::Type;
CC cc;
constraints(bctype,gfs,cc);

// Local operator (weak form on one element) and
// grid operator (assembly over the whole grid)
using LOP = ConvectionDiffusionFEM<Problem,FEM>;
LOP lop(problem);
using MBE = ISTL::BCRSMatrixBackend<>;
using GO = GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC>;
GO go(gfs,cc,gfs,cc,lop,MBE(25));

// Coefficient vector, initialized with the Dirichlet values
using X = Backend::Vector<GFS,double>;
X x(gfs,0.0);
ConvectionDiffusionDirichletExtensionAdapter g(gv,problem);
interpolate(g,gfs,x);

// Assemble and solve with AMG-preconditioned CG from dune-istl
using LS = ISTLBackend_SEQ_CG_AMG_SSOR<GO>;
LS ls(100,1);
StationaryLinearProblemSolver<GO,LS,X> slp(go,ls,x,1e-10);
slp.apply();

// Write the solution for visualization, e.g. in ParaView
Dune::VTKWriter vtkwriter(gv);
addSolutionToVTKWriter(vtkwriter,gfs,x);
vtkwriter.write("poisson");
```

The [recipe on linear system assembly](@ref recipe-linear-system-assembly)
and the [tutorials][tutorials] explain these steps in detail,
including how to implement your own local operator.

## Getting started {#pdelab-getting-started}

### Tutorials
The [dune-pdelab-tutorials][tutorials] module contains complete, documented
examples together with the mathematical background:

| Tutorial | Topic |
|----------|-------|
| [tutorial00][t00] | Poisson equation, P1 conforming finite elements on simplicial meshes |
| [tutorial01][t01] | Nonlinear Poisson equation with conforming finite elements of various orders on simplicial and cube meshes |
| [tutorial02][t02] | Nonlinear Poisson equation with the cell-centered finite volume method |
| [tutorial03][t03] | Instationary nonlinear Poisson equation |
| [tutorial04][t04] | Wave equation, solved as a system of two equations |
| [tutorial05][t05] | Nonlinear Poisson equation with adaptive finite elements |
| [tutorial06][t06] | Nonlinear Poisson equation in parallel |
| [tutorial07][t07] | Discontinuous Galerkin methods for hyperbolic conservation laws |
| [tutorial08][t08] | Unsteady incompressible Navier-Stokes equations |
| [tutorial09][t09] | Generating local operators with code generation |

The tutorials module also contains introductions to C++, to the simulation
workflow (build system, meshes, and visualization), and to the DUNE grid interface.

### Recipes
The @ref pdelab_recipes "recipes" are short, compiled examples showing how to
achieve common tasks, e.g., operations on grid functions, setting up blocked
data structures, assembling and solving linear systems, and communication in
parallel programs.

### Further resources
* The [beginners resources for PDELab][beginners] on the DUNE website.
* The [DUNE book][dunebook] (O. Sander, Springer, 2020) gives an in-depth
  introduction to the DUNE core modules, dune-typetree, dune-functions, and PDELab.
* The [installation instructions][installation] of the DUNE website.
  The directory `doc/installation` of this module contains scripts to set up
  DUNE and PDELab on Ubuntu and macOS.

## Documentation map {#pdelab-documentation}

* @reftopics: the API documentation grouped by topic, e.g.,
  [Grid Function Space](@ref GridFunctionSpace),
  [Local Operator](@ref LocalOperator),
  [Backends](@ref Backend), and
  [Entity Set Partitions](@ref Partition).
* The @ref GridOperatorDocPage "Grid Operator" page explains the assembly
  and the composition of grid operators for instationary problems.
* The [Ordering](@ref Ordering) topic explains how degrees of freedom are
  mapped to vector and matrix entries.
* The @ref pdelab_recipes "recipes" page lists all recipes.

## Dependencies and building {#pdelab-building}

PDELab depends on the DUNE core modules
([dune-common][common], [dune-geometry][geometry], [dune-grid][grid],
[dune-localfunctions][localfunctions], [dune-istl][istl])
and on the staging modules [dune-typetree][typetree] and [dune-functions][functions].
For the versions required, see the file `dune.module`.

PDELab is built like any other DUNE module with `dunecontrol`, see the
[installation instructions][installation].
This documentation is built with `make doc` and can then be found in
`[build-directory]/doc/doxygen/html/`.

## Getting help and contributing {#pdelab-help}

* Questions can be asked on the [PDELab mailing list][mailinglist].
* Bugs can be reported in the [issue tracker][issues] (requires an account on the DUNE GitLab).
* The source code is developed in the [PDELab GitLab repository][gitlab].
  Contributions in the form of merge requests are welcome.

## Citing PDELab {#pdelab-citing}

If you use PDELab in your scientific work, please cite
@cite BastianHeimannMarnach2010 and the DUNE framework
@cite BastianEtAl2021 @cite BlattEtAl2025.

## License {#pdelab-license}

PDELab is free software. It can be used, modified, and redistributed under the
terms of either the GNU Lesser General Public License, version 3 or later,
or the GNU General Public License, version 2, with a special runtime exception.
See the file `LICENSE.md` of this module for details.

[pdelab]: https://www.dune-project.org/modules/dune-pdelab/
[dune]: https://www.dune-project.org
[tutorials]: https://www.dune-project.org/modules/dune-pdelab-tutorials/
[t00]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial00
[t01]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial01
[t02]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial02
[t03]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial03
[t04]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial04
[t05]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial05
[t06]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial06
[t07]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial07
[t08]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial08
[t09]: https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials/-/tree/master/tutorial09
[beginners]: https://www.dune-project.org/doc/gettingstarted/beginners-resources-pdelab/
[installation]: https://www.dune-project.org/installation/installation-buildsrc/
[dunebook]: https://doi.org/10.1007/978-3-030-59702-3
[mailinglist]: https://lists.dune-project.org/mailman/listinfo/dune-pdelab
[issues]: https://gitlab.dune-project.org/pdelab/dune-pdelab/-/issues
[gitlab]: https://gitlab.dune-project.org/pdelab/dune-pdelab
[common]: https://gitlab.dune-project.org/core/dune-common
[geometry]: https://gitlab.dune-project.org/core/dune-geometry
[grid]: https://gitlab.dune-project.org/core/dune-grid
[localfunctions]: https://gitlab.dune-project.org/core/dune-localfunctions
[istl]: https://gitlab.dune-project.org/core/dune-istl
[typetree]: https://gitlab.dune-project.org/staging/dune-typetree
[functions]: https://gitlab.dune-project.org/staging/dune-functions
