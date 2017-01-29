## General

## How to use it
Prerequisites are <a href = "https://cmake.org/"> Cmake</a> and <a href = "https://github.com/libigl/libigl"> libigl </a>. The test unit test program uses <a href = "https://github.com/philsquared/Catch"> catch</a>, which is included with the project.

This project uses the data representation and input output routines of libigl.

The project can be used as a library, here are some routines that can be useful.

* `interpolate.h` contains routines that find the distance from one surface to another at each vertex, and routines to modify surfaces based on these distances.
* `laplace_smooth.h` contains a routine for smoothing the normal vectors of a surfaces.
* `skirt.h` contains routines for finding and organizing border vertices, and routines for creating a skirt around an open mesh.
* `thicken.h` contains routines for turning an open mesh into a solid.

To build the test and the sample applications, use Cmake. You may have to set `libigl_location`, which is predefined as a sibling directory to the meat and bone source code directory.

The sample application can be invoked as follows:

```console
libmeat_and_bone_bin <work surface> <reference surface> <output file>
```
