#ifndef SINGLER_CLASSIC_MARKERS_UTILS_HPP
#define SINGLER_CLASSIC_MARKERS_UTILS_HPP

#include <type_traits>

namespace singler_classic_markers {

template<typename Input_>
using I = std::remove_cv_t<std::remove_reference_t<Input_> >;

}

#endif
