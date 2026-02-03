#include <gtest/gtest.h>

#include "utils.h"
#include "spawn_matrix.h"

#include "singler_classic_markers/choose.hpp"

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "topicks/topicks.hpp"

class ChooseTest : public ::testing::TestWithParam<int> {
protected:
    static std::vector<std::vector<std::vector<std::pair<int, double> > > > reference(
        const tatami::Matrix<double, int>& matrix,
        const int* label,
        int num_keep
    ) {
        auto medians = tatami_stats::grouped_medians::by_row(matrix, label, {});
        std::size_t nlabels = medians.size();

        const auto NR = matrix.nrow();
        std::vector<double> buffer(NR);
        std::vector<std::vector<std::vector<std::pair<int, double> > > > output(nlabels);

        for (std::size_t l = 0; l < nlabels; ++l) {
            output[l].resize(nlabels);
            for (std::size_t l2 = 0; l2 < nlabels; ++l2) {
                for (int r = 0; r < NR; ++r) {
                    buffer[r] = medians[l][r] - medians[l2][r];
                }

                topicks::PickTopGenesOptions<double> opt;
                opt.keep_ties = false; 
                opt.bound = 0;

                auto keep = topicks::pick_top_genes_index(NR, buffer.data(), num_keep, true, opt);
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
    }
};

TEST_P(ChooseTest, Basic) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    auto mat = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 0.3);
    size_t nlabels = 5;
    auto labels = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = requested;
    mopt.keep_ties = false;

    auto output = singler_classic_markers::choose(*mat, labels.data(), mopt);
    auto ref = reference(*mat, labels.data(), requested);
    EXPECT_EQ(output, ref);

    // Same result with sparse.
    auto smat = tatami::convert_to_compressed_sparse<double, int>(*mat, true, {});
    auto soutput = singler_classic_markers::choose(*mat, labels.data(), mopt);
    EXPECT_EQ(output, soutput);

    // Works with indices only.
    auto ioutput = singler_classic_markers::choose_index(*mat, labels.data(), mopt);
    EXPECT_EQ(ioutput, strip_to_indices(output));

    // Same result when parallelized.
    mopt.num_threads = 3;
    auto outputp = singler_classic_markers::choose(*mat, labels.data(), mopt);
    EXPECT_EQ(output, outputp);
}

TEST_P(ChooseTest, Missing) { 
    size_t ngenes = 500;
    size_t nsamples = 50;
    int requested = GetParam();

    auto mat = spawn_matrix(ngenes, nsamples, /* seed = */ 1234 * requested, /* density = */ 0.3);
    size_t nlabels = 5;
    auto labels = spawn_labels(nsamples, nlabels, /* seed = */ 6789 * requested);

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = requested;
    auto output = singler_classic_markers::choose(*mat, labels.data(), mopt);

    // Check that we are robust to missing labels.
    for (auto& l : labels) {
        ++l;
    }
    auto output2 = singler_classic_markers::choose(*mat, labels.data(), mopt);

    for (std::size_t l = 0; l <= nlabels; ++l) {
        for (std::size_t l2 = 0; l2 <= nlabels; ++l2) {
            if (l == l2 || l2 == 0 || l == 0) {
                EXPECT_TRUE(output2[l][l2].empty());
                continue;
            }

            EXPECT_EQ(output[l-1][l2-1], output2[l][l2]);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    Choose,
    ChooseTest,
    ::testing::Values(1, 20, 50, 1000) // number of top genes.
);

TEST(Choose, Ties) { 
    size_t ngenes = 100;
    size_t nlabels = 2;
    std::vector<double> expression(nlabels * ngenes);
    std::fill(expression.begin(), expression.begin() + 50, 1);
    std::fill(expression.begin() + 150, expression.end(), 1);
    tatami::DenseColumnMatrix<double, int> mat(ngenes, nlabels, std::move(expression));

    singler_classic_markers::ChooseOptions mopt;
    mopt.number = 10;
    std::vector<int> grouping { 0, 1 };
    auto tied = singler_classic_markers::choose(mat, grouping.data(), mopt);

    // Ties are broken in a stable way.
    std::vector<std::pair<int, double> > expected1 { {0,1.}, {1,1.0}, {2,1.0}, {3,1.0}, {4,1.0}, {5,1.0}, {6,1.0}, {7,1.0}, {8,1.0}, {9,1.0} };
    EXPECT_EQ(tied[0][1], expected1);
    std::vector<std::pair<int, double> > expected2 { {50,1.}, {51,1.0}, {52,1.0}, {53,1.0}, {54,1.0}, {55,1.0}, {56,1.0}, {57,1.0}, {58,1.0}, {59,1.0} };
    EXPECT_EQ(tied[1][0], expected2);

    // Unless all ties are requested.
    mopt.keep_ties = true;
    tied = singler_classic_markers::choose(mat, grouping.data(), mopt);
    EXPECT_EQ(tied[0][1].size(), 50);
    EXPECT_EQ(tied[1][0].size(), 50);
}
