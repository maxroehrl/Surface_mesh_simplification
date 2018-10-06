# Surface_mesh_simplification
A CGAL Polyhedron demo plugin for surface mesh simplification with color based constraints

## Inspiration
This project is based on these CGAL Polyhedron demo files:
* [Mesh_simplification_plugin.cpp](https://github.com/CGAL/cgal/blob/master/Polyhedron/demo/Polyhedron/Plugins/Surface_mesh/Mesh_simplification_plugin.cpp)
* [Mesh_simplification_dialog.ui](https://github.com/CGAL/cgal/blob/master/Polyhedron/demo/Polyhedron/Plugins/Surface_mesh/Mesh_simplification_dialog.ui)

## Build requirements
* [Cmake](https://cmake.org/download/)
* [Boost](https://www.boost.org/users/download/)
* [CGAL](https://github.com/CGAL/cgal/releases)
* [QT5](https://www.qt.io/download)
* [Eigen3](https://github.com/eigenteam/eigen-git-mirror)
* [Intel TBB](https://github.com/01org/tbb/releases)

## How to build
In the CGAL directory run:
```
cd demo\Polyhedron\Plugins
git clone https://github.com/maxroehrl/Surface_mesh_simplification.git
```
Make sure all environment variables for the build requirements are set like shown [here](https://www.cgal.org/download/windows.html) for Windows.

Run cmake on the Polyhedron demo and open the generated solution in Visual Studio.

Build the ALL_BUILD project to build the whole demo.

## How to use
Start the Polyhedron demo and open an OFF file with color information with the "surface_mesh_io_plugin".

Right click on the mesh and under "Operations" select "Triangulated Surface Mesh Simplification" and "Color Constrained Simplification".

In the dialog select your color source and start the simplification process.
