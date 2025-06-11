#pragma once

#include <memory>
#include <string>
#include <vector>

#include "spoa-sys/spoa/include/spoa/spoa.hpp"

namespace spoa {

std::unique_ptr<spoa::AlignmentEngine>
create_alignment_engine(spoa::AlignmentType type, std::int8_t m, std::int8_t n,
                        std::int8_t g, std::int8_t e, std::int8_t q,
                        std::int8_t c);

std::unique_ptr<spoa::Alignment> align(spoa::AlignmentEngine &engine,
                                       const char *sequence,
                                       std::uint32_t sequence_len,
                                       const Graph &graph);

std::unique_ptr<spoa::Graph> create_graph();

void add_alignment_with_qual(spoa::Graph &graph,
                             const spoa::Alignment &alignment,
                             const char *sequence, std::uint32_t sequence_len,
                             const char *quality, std::uint32_t quality_len);

void add_alignment(spoa::Graph &graph, const spoa::Alignment &alignment,
                   const char *sequence, std::uint32_t sequence_len,
                   std::uint32_t weight);

void graph_clear(spoa::Graph &graph);

std::unique_ptr<std::string> generate_consensus(spoa::Graph &graph);

std::unique_ptr<std::string>
generate_consensus_with_min_coverage(spoa::Graph &graph,
                                     std::int32_t min_coverage);

std::unique_ptr<std::vector<std::string>>
generate_multiple_sequence_alignment(spoa::Graph &graph,
                                     bool include_consensus);

} // namespace spoa
