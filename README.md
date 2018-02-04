## INTRODUCTION
GEL is a small library with a good support for 2,3,4 D vectors, simple linear system solving (via Lapack), Mesh representation of surfaces, and the volume representation. There are also several other things some of which pertain to simple user interfaces. For instance, GEL contains a console library.

## DOCUMENTATION
Some installation instructions below. But for more documentation please see the doc directory. There is a doxygen script for creating a reference manual and a latex file intro.tex which explains the basics. Please doxygen or pdflatex your documentation. A license is also found in the intro document.

### Building on XCode/OSX
An XCode project has been created and is found in the GEL_MAC directory. The XCode project produces a framework for the GEL library.

### Building on Windows
An Visual Studio 2013 project has been created and is found in the GEL_WIN directory. The project produces a framework for the GEL library.
Look in the doc folder for installation of glut and glew

### Building on Linux
TODO

### Building with CMake
There is a basic CMakeLists.txt - it compiles the GEL library and the MeshEditGlut demo.
```
#$> mkdir build; cd build; cmake ..
#$> make -j8
```
