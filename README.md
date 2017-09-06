## General
This library contains routines for adding a skirt to smooth the edges of surfaces, turning surfaces into 3d shapes, and interpolating surfaces. It was made for modifying orthopedic devices to be more comfortable pressing against soft tissue and bone--- or in plain language, to reduce chafing. <a href = "https://eklarlund.github.io/meat-and-bone/">Click here for more information. </a>

## How to use it
Prerequisites are <a href = "https://cmake.org/"> Cmake</a> and <a href = "https://github.com/libigl/libigl"> libigl </a>. The test unit test program uses <a href = "https://github.com/philsquared/Catch"> catch</a>, which is included with the project.

This project uses the data representation and input output routines of libigl.

The project can be used as a library, here are some routines that can be useful.

* `interpolate.h`  displaces a work surface towards or away from a reference surface along a normal field -- just two intuitive intervals, four reals, define the transformation
* `laplace_smooth.h` displaces vertices for smoother appearance
* `make_skirt.h` fillets a border of a surface, that is extends the edge in a nice rounded manner
* `make_bitangents.h` smooths the curve vectors of a border for even smoother skirts
* `make_solid.h` thickens a surface into a solid

To build the test and the sample applications, use Cmake. You may have to set `libigl_location`, which is predefined as a sibling directory to the meat and bone source code directory.

The sample application can be invoked from the build directory as follows:

```console
./libmeat_and_bone_bin testdata/bone_example/work.obj testdata/bone_example/ref.obj result.obj
```

Inspiration and Guidance: Nils Klarlund
