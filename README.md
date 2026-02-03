# Classic marker detection in SingleR

![Unit tests](https://github.com/SingleR-inc/singler_classic_markers/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/SingleR-inc/singler_classic_markers/actions/workflows/doxygenate.yaml/badge.svg)
[![codecov](https://codecov.io/gh/SingleR-inc/singler_classic_markers/branch/master/graph/badge.svg?token=q3AJjQtKtj)](https://codecov.io/gh/SingleR-inc/singler_classic_markers)

## Overview

This repository implements the classic marker detection algorithm from [**SingleR**](https://bioconductor.org/packages/SingleR). 
For each pairwise comparison between labels, genes are ranked based on the difference in the medians of log-expression values.
The top genes are then used as markers in the **SingleR** algorithm.
We also provide methods for combining information across blocks in the presence of, e.g., batch effects.

## Quick start

**singler_classic_markers** is a header-only library, so it can be easily used by just `#include`ing the relevant source files.

```cpp
#include "singler_classic_markers/singler_classic_markers.hpp"

// Prepare the reference matrix as a tatami::NumericMatrix.
ref_mat;

// Prepare a vector of labels, one per column of ref_mat.
ref_labels;

// Computing the top markers, which can be directly used in singlepp::train_single().
singler_classic_markers::ChooseOptions opts;
auto markers = singler_classic_markers::choose_index(ref_mat, ref_labels.data(), opts); 
```

See the [reference documentation](https://singler-inc.github.io/singler_classic_markers) for more details.

## Alternatives

As its name suggests, the classic method was part of the original implementation of the **SingleR** method.
It is simple and works well for the built-in reference datasets, but is less appropriate for single-cell references where the median is often zero within each label.
In such cases, we suggest using the `score_markers_best()` function from the [**scran_markers**](https://github.com/libscran/scran_markers) library,
which ranks genes on effect sizes that are more robust to a large number of zeros.

For bulk RNA-seq references, we could also consider using tools like [**edgeR**](https://bioconductor.org/packages/edgeR) or [**DESeq2**](https://bioconductor.org/packages/DESeq2).
Briefly, we could perform differential expression comparisons between all labels (possibly using TREAT to enforce a minimum log-fold change) and then take the top genes as markers.
This would likely yield better markers than **SingleR**'s classic approach,
as **edgeR** and friends consider the variance within each label and will favor genes that have more consistent behavior.

## Building projects 

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  singler_classic_markers
  GIT_REPOSITORY https://github.com/SingleR-inc/singler_classic_markers
  GIT_TAG master # replace with a pinned release
)

FetchContent_MakeAvailable(singler_classic_markers)
```

Then you can link to **singler_classic_markers** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe singler_classic_markers)

# For libaries
target_link_libraries(mylib INTERFACE singler_classic_markers)
```

By default, this will use `FetchContent` to fetch all external dependencies.
Applications are advised to pin the versions of all dependencies - see [`extern/CMakeLists.txt`](extern/CMakeLists.txt) for suggested versions.
If you want to install them manually, use `-DSINGLER_CLASSIC_MARKERS_FETCH_EXTERN=OFF`.

### CMake with `find_package()`

```cmake
find_package(singler_singler_classic_markers CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE singler::singler_classic_markers)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DSINGLER_CLASSIC_MARKERS_TESTS=OFF
cmake --build . --target install
```

Again, this will use `FetchContent` to retrieve dependencies, see comments above.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This assumes that the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) are available during compilation.

## References

Aran D et al. (2019). 
Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
_Nat. Immunol._ 20, 163-172
