#include <gtest/gtest.h>

#include "singler_classic_markers/number.hpp"

TEST(DefaultNumber, Classic) {
    // Just getting some test coverage here.
    EXPECT_EQ(singler_classic_markers::default_number(0), 0);
    EXPECT_GT(singler_classic_markers::default_number(100), 0);
    EXPECT_GT(singler_classic_markers::default_number(1000), 0);
    EXPECT_EQ(singler_classic_markers::get_num_keep<int>(100, {}), singler_classic_markers::default_number(100));
    EXPECT_EQ(singler_classic_markers::get_num_keep<int>(100, 21), 21);
}
