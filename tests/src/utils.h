#ifndef UTILS_H
#define UTILS_H

#include <vector>

template<typename Index_, typename Stat_>
std::vector<std::vector<std::vector<Index_> > > strip_to_indices(const std::vector<std::vector<std::vector<std::pair<Index_, Stat_> > > >& input) {
    std::vector<std::vector<std::vector<Index_> > > output;
    output.reserve(input.size());

    for (const auto& x : input) {
        std::vector<std::vector<Index_> > out2;
        out2.reserve(x.size());

        for (const auto& y : x) {
            std::vector<Index_> out3;
            out3.reserve(y.size());

            for (const auto& z : y){
                out3.push_back(z.first);
            }

            out2.push_back(out3);
        }

        output.push_back(out2);
    }

    return output;
}

#endif
