#pragma once
#include "branching_rules_primitives.hpp"
#include "reduction_rules.hpp"
#include "minimal_solutions_computation.hpp"
#include "dominance_free_solutions_computation.hpp"
#include "adjusted_solutions_computation.hpp"

BranchingRule generate_branching_rule(const _Graph & g) {
    vector<int> red_vertices;
    for(int i = 0; i < g.n; ++i) {
        if(g.red_vertices_mask & (1<<i)) {
            red_vertices.push_back(i);
        }
    }

    if(is_g_rv_handled_by_red_component_reduction(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_rcr");
    }

    if(is_g_rv_handled_by_red_star_reduction(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_rsr");
    }

    const auto & minimal_solutions = compute_all_minimal_solutions(g);
    double bf_minimal_solutions = compute_branching_factor_for_solutions(minimal_solutions);

    if(bf_minimal_solutions <= DPVC_BF) {
        return BranchingRule(g, {}, minimal_solutions, {}, minimal_solutions, bf_minimal_solutions, "minimal");
    }

    const auto & dominance_free_solutions = compute_all_dominance_free_solutions(minimal_solutions, g, red_vertices);

    const auto & adjusted_solutions = compute_adjusted_solutions_bit_masks(g, dominance_free_solutions);
    double bf_adjusted_solutions = compute_branching_factor_for_solutions(adjusted_solutions);

    return BranchingRule(g, red_vertices, minimal_solutions, dominance_free_solutions, adjusted_solutions, bf_adjusted_solutions, "adjusted_dominance");
}
