#include <gtest/gtest.h>

#include <cmath>
#include <vector>
#include <cstddef>

#include "utils.h"
#include "spawn_matrix.h"

#include "singler_classic_markers/blocked.hpp"
#include "singler_classic_markers/choose.hpp"

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "topicks/topicks.hpp"

class BlockedTest : public ::testing::TestWithParam<int> {};

TEST_P(BlockedTest, Single) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    auto mat = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 0.3);
    size_t nlabels = 5;
    auto labels = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = requested;
    mopt.keep_ties = false;
    auto simple = singler_classic_markers::choose(*mat, labels.data(), mopt);

    singler_classic_markers::ChooseBlockedOptions bopt;
    bopt.number = requested;
    bopt.keep_ties = false;
    std::vector<int> blocks(nsamples);
    auto blocked = singler_classic_markers::choose_blocked(*mat, labels.data(), blocks.data(), bopt);
    EXPECT_EQ(simple, blocked);

    // Same result with the minimum.
    auto min_bopt = bopt;
    min_bopt.use_minimum = true;
    auto min_blocked = singler_classic_markers::choose_blocked(*mat, labels.data(), blocks.data(), min_bopt);
    EXPECT_EQ(simple, min_blocked);

    // Works with indices only.
    auto iblocked = singler_classic_markers::choose_blocked_index(*mat, labels.data(), blocks.data(), bopt);
    EXPECT_EQ(iblocked, strip_to_indices(blocked));

    // Same result when parallelized.
    bopt.num_threads = 3;
    auto blockedp = singler_classic_markers::choose_blocked(*mat, labels.data(), blocks.data(), bopt);
    EXPECT_EQ(simple, blockedp);
}

TEST_P(BlockedTest, Multiple) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    // We need higher densities to get enough genes for a meaningful test with
    // use_minimum=true, otherwise nothing is selected because the minimum is
    // zero for most genes and negative for half of the test.
    auto mat1 = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 1);
    auto mat2 = spawn_matrix(ngenes, nsamples, /* seed = */ 9876 * requested, /* density = */ 1);

    size_t nlabels = 5;
    auto labels1 = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);
    auto labels2 = spawn_labels(nsamples, nlabels, /* seed = */ 4321 * requested);

    auto compute_reference = [&](bool use_minimum) -> auto {
        auto medians1 = tatami_stats::grouped_medians::by_row(*mat1, labels1.data(), {});
        auto medians2 = tatami_stats::grouped_medians::by_row(*mat2, labels2.data(), {});
        std::vector<std::vector<std::vector<std::pair<int, double> > > > output(nlabels);

        std::vector<double> buffer(ngenes);
        for (std::size_t l = 0; l < nlabels; ++l) {
            output[l].resize(nlabels);

            for (std::size_t l2 = 0; l2 < nlabels; ++l2) {
                if (use_minimum) {
                    for (size_t r = 0; r < ngenes; ++r) {
                        buffer[r] = std::min(medians1[l][r] - medians1[l2][r], medians2[l][r] - medians2[l2][r]);
                    }
                } else {
                    for (size_t r = 0; r < ngenes; ++r) {
                        buffer[r] = ((medians1[l][r] - medians1[l2][r]) + (medians2[l][r] - medians2[l2][r]))/2.0;
                    }
                }

                topicks::PickTopGenesOptions<double> opt;
                opt.keep_ties = false; 
                opt.bound = 0;

                auto keep = topicks::pick_top_genes_index<int>(ngenes, buffer.data(), requested, true, opt);
                auto& result = output[l][l2];
                for (auto k : keep) {
                    output[l][l2].emplace_back(k, buffer[k]);
                }

                std::sort(result.begin(), result.end(), [](const std::pair<int, double>& left, const std::pair<int, double>& right) -> bool {
                    if (left.second == right.second) {
                        return left.first < right.first;
                    } else {
                        return left.second > right.second;
                    }
                });
            }
        }

        return output;
    };

    tatami::DelayedBind<double, int> combined(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{ mat1, mat2 }, false);
    std::vector<int> blocks(nsamples);
    blocks.insert(blocks.end(), nsamples, 1);
    std::vector<int> labels = labels1;
    labels.insert(labels.end(), labels2.begin(), labels2.end());

    singler_classic_markers::ChooseBlockedOptions bopt;
    bopt.number = requested;
    bopt.keep_ties = false;
    auto mean_blocked = singler_classic_markers::choose_blocked(combined, labels.data(), blocks.data(), bopt);
    auto mean_ref = compute_reference(false);
    EXPECT_EQ(mean_blocked, mean_ref);

    // Same result with the minimum.
    auto min_bopt = bopt;
    min_bopt.use_minimum = true;
    auto min_blocked = singler_classic_markers::choose_blocked(combined, labels.data(), blocks.data(), min_bopt);
    auto min_ref = compute_reference(true);
    EXPECT_EQ(min_blocked, min_ref);

    // Same result when parallelized.
    bopt.num_threads = 3;
    auto blockedp = singler_classic_markers::choose_blocked(combined, labels.data(), blocks.data(), bopt);
    EXPECT_EQ(mean_blocked, blockedp);
}

TEST_P(BlockedTest, Overlap) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    // We need higher densities to get enough genes for a meaningful test with
    // use_minimum=true, otherwise nothing is selected because the minimum is
    // zero for most genes and negative for half of the test.
    auto mat1 = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 1);
    auto mat2 = spawn_matrix(ngenes, nsamples, /* seed = */ 9876 * requested, /* density = */ 1);

    size_t nlabels = 5;
    auto labels1 = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);
    auto labels2 = spawn_labels(nsamples, nlabels, /* seed = */ 4321 * requested);

    std::shared_ptr<tatami::Matrix<double, int> > combined(
        new tatami::DelayedBind<double, int>(std::vector<std::shared_ptr<tatami::Matrix<double, int> > >{ mat1, mat2 }, false)
    );
    std::vector<int> blocks(nsamples);
    blocks.insert(blocks.end(), nsamples, 1);
    std::vector<int> labels = labels1;
    for (auto ll : labels2) {
        labels.push_back(ll + 1);
    }

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = requested;
    auto basic1 = singler_classic_markers::choose(*mat1, labels1.data(), mopt);
    auto basic2 = singler_classic_markers::choose(*mat2, labels2.data(), mopt);

    std::vector<int> subset, sublabels, subblocks;
    for (std::size_t i = 0, end = labels.size(); i < end; ++i) {
        if (labels[i] >= 1 && labels[i] < 5) {
            subset.push_back(i);
            sublabels.push_back(labels[i] - 1);
            subblocks.push_back(blocks[i]);
        }
    }
    auto subcombined = tatami::make_DelayedSubset<double, int>(combined, std::move(subset), false);

    auto compare_results = [&](const auto& blocked, const auto& subblocked) -> void {
        for (std::size_t l = 0; l <= nlabels; ++l) {
            for (std::size_t l2 = 0; l2 <= nlabels; ++l2) {
                if (l == l2) {
                    EXPECT_TRUE(blocked[l][l2].empty());
                    continue;
                }

                if (l == 0) {
                    if (l2 == 5) {
                        EXPECT_TRUE(blocked[l][l2].empty());
                    } else {
                        EXPECT_EQ(blocked[l][l2], basic1[l][l2]);
                    }
                } else if (l == 5) {
                    if (l2 == 0) {
                        EXPECT_TRUE(blocked[l][l2].empty());
                    } else {
                        EXPECT_EQ(blocked[l][l2], basic2[l - 1][l2 - 1]);
                    }
                } else {
                    if (l2 == 0) {
                        EXPECT_EQ(blocked[l][l2], basic1[l][l2]);
                    } else if (l2 == 5) {
                        EXPECT_EQ(blocked[l][l2], basic2[l - 1][l2 - 1]);
                    } else {
                        EXPECT_EQ(blocked[l][l2], subblocked[l - 1][l2 - 1]);
                    }
                }
            }
        }
    };

    // For the mean.
    {
        singler_classic_markers::ChooseBlockedOptions bopt;
        bopt.number = requested;
        bopt.keep_ties = false;
        auto mean_blocked = singler_classic_markers::choose_blocked(*combined, labels.data(), blocks.data(), bopt);
        auto subblocked = singler_classic_markers::choose_blocked(*subcombined, sublabels.data(), subblocks.data(), bopt);
        compare_results(mean_blocked, subblocked);
    }

    // For the minimum.
    {
        singler_classic_markers::ChooseBlockedOptions bopt;
        bopt.number = requested;
        bopt.keep_ties = false;
        bopt.use_minimum = true;
        auto mean_blocked = singler_classic_markers::choose_blocked(*combined, labels.data(), blocks.data(), bopt);
        auto subblocked = singler_classic_markers::choose_blocked(*subcombined, sublabels.data(), subblocks.data(), bopt);
        compare_results(mean_blocked, subblocked);
    }
}

TEST_P(BlockedTest, Impossible) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    auto mat = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 0.3);
    size_t nlabels = 5;
    auto labels = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);

    auto blocks = labels;
    for (auto& b : blocks) {
        b = (b == 3);
    }

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = requested;
    mopt.keep_ties = false;

    std::vector<int> subset, sublabels;
    for (std::size_t i = 0, end = labels.size(); i < end; ++i) {
        if (labels[i] != 3) {
            subset.push_back(i);
            sublabels.push_back(labels[i]);
        }
    }
    auto submat = tatami::make_DelayedSubset<double, int>(mat, std::move(subset), false);

    auto blocked = singler_classic_markers::choose_blocked(*mat, labels.data(), blocks.data(), {});
    auto subblocked = singler_classic_markers::choose(*submat, sublabels.data(), {});
    for (std::size_t l = 0; l < nlabels; ++l) {
        for (std::size_t l2 = 0; l2 < nlabels; ++l2) {
            if (l == l2 || l == 3 || l2 == 3) {
                EXPECT_TRUE(blocked[l][l2].empty());
                EXPECT_TRUE(subblocked[l][l2].empty());
                continue;
            }

            EXPECT_EQ(blocked[l][l2], subblocked[l][l2]);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    Blocked,
    BlockedTest,
    ::testing::Values(1, 20, 50, 1000) // number of top genes.
);
