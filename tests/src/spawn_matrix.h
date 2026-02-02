#ifndef SPAWN_MATRIX_H
#define SPAWN_MATRIX_H

#include <vector>
#include <random>
#include <memory>

#include "tatami/tatami.hpp"

inline std::shared_ptr<tatami::Matrix<double, int> > spawn_matrix(size_t nr, size_t nc, int seed, double density) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<> dist;
    std::uniform_real_distribution<> udist;

    std::vector<double> contents(nr * nc);
    for (auto& c : contents) {
        c = (udist(rng) <= density ? dist(rng) : 0.0);
    }

    return std::shared_ptr<tatami::Matrix<double, int> >(new tatami::DenseColumnMatrix<double, int>(nr, nc, std::move(contents)));
}

inline std::vector<int> spawn_labels(size_t nc, size_t nlabels, int seed) {
    std::vector<int> labels(nc);
    std::iota(labels.begin(), labels.begin() + nlabels, 0); // at least one entry per label.
    std::mt19937_64 rng(seed);
    for (size_t i = nlabels; i < nc; ++i) {
        labels[i] = rng() % nlabels;
    }
    std::shuffle(labels.begin(), labels.end(), rng);
    return labels;
}

#endif
