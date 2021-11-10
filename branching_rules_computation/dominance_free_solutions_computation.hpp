#pragma once
#include "branching_rules_primitives.hpp"

vector<bool>  _is_solution_i_dominated_by_solution_j_domination_masks(int s_i, const vector<int> & solutions_bit_masks,
    const vector<int> & g_adj_lists_bit_masks, const vector<int> & safe_red_vertices_subsets_bit_masks) {
    vector<bool> solution_j_domination_masks(solutions_bit_masks.size(), false);

    int solution_i = solutions_bit_masks[s_i];

    // Go over all safe red_vertices subsets.
    for(int safe_rv_subset : safe_red_vertices_subsets_bit_masks) {
        auto solution_i_rv_subset_intersection = set_intersection_bit_masks(solution_i, safe_rv_subset);

        // If the solution does not intersect the rv_subset, it will not help.
        if(size_bit_masks(solution_i_rv_subset_intersection)==0) {
            continue;
        }
        // Assume that everything from the solution outside the rv_subset is deleted.
        auto deleted_vertices = set_minus_bit_masks(solution_i, solution_i_rv_subset_intersection);

        // Compute the "blue" neighborhood of the rv_subset after deleting the vertices.
        int solution_i_safe_rv_subset_neighborhood = 0;
        for(int v : bit_mask_to_vector(safe_rv_subset)) solution_i_safe_rv_subset_neighborhood = set_union_bit_masks(
            solution_i_safe_rv_subset_neighborhood,
            set_minus_bit_masks(g_adj_lists_bit_masks[v], deleted_vertices));
        solution_i_safe_rv_subset_neighborhood = set_minus_bit_masks(solution_i_safe_rv_subset_neighborhood, safe_rv_subset);

        // If the "blue" neighborhood is larger than the intersection, no bottleneck occurs, it will not help.
        if(size_bit_masks(solution_i_rv_subset_intersection) < size_bit_masks(solution_i_safe_rv_subset_neighborhood)) {
            continue;
        }
        if(size_bit_masks(solution_i_safe_rv_subset_neighborhood) == 0) {
            continue;
        }

        // We go over all solutions other than solution_i.
        for(int s_j = 0; s_j < solutions_bit_masks.size(); ++s_j) {
            if(s_j==s_i) continue;

            // If the rest of the solution_j is a subset of the bottleneck, then the solution_j dominates solution_i.
            // That is because one can create a better solution solution_i' from solution_i and solution_j is a subset of solution_i'.
            if(is_subset_bit_masks(set_minus_bit_masks(solutions_bit_masks[s_j], deleted_vertices), solution_i_safe_rv_subset_neighborhood)) {
                solution_j_domination_masks[s_j] = true;
            }
        }
    }

    return solution_j_domination_masks;
}

vector<vector<int>> construct_dominance_graph(const vector<vector<int>> & solutions, const _Graph & g, const vector<int> & red_vertices) {
    int dg_n = solutions.size();
    vector<vector<int>> dg_adj(dg_n);
    map<vector<int>, int> solution_to_v;
    for(int i = 0; i < solutions.size(); ++i) {
        solution_to_v[solutions[i]] = i;
    }

    if(red_vertices.size()==0) {
        return dg_adj;
    }

    vector<vector<int>> g_adj_lists = g.get_adjacent_lists();
    vector<vector<int>> safe_red_vertices_subsets;
    for(auto rv_subset : ordered_powerset_nonempty(red_vertices)) {
        if(is_g_solved_by_solution(g, set_minus(range(g.n), rv_subset))) {
            safe_red_vertices_subsets.push_back(rv_subset);
        }
    }

    vector<int> g_adj_lists_bit_masks = convert_vectors_to_bit_masks(g_adj_lists);
    vector<int> safe_red_vertices_subsets_bit_masks = convert_vectors_to_bit_masks(safe_red_vertices_subsets);
    vector<int> solutions_bit_masks = convert_vectors_to_bit_masks(solutions);

    for(int s_i = 0; s_i < solutions.size(); ++s_i) {
        vector<bool> solution_j_domination_masks = _is_solution_i_dominated_by_solution_j_domination_masks(
            s_i, solutions_bit_masks, g_adj_lists_bit_masks, safe_red_vertices_subsets_bit_masks);
        for(int s_j = 0; s_j < solutions.size(); ++ s_j) {
            if(s_i == s_j) continue;
            if(solution_j_domination_masks[s_j]) {
                dg_adj[solution_to_v[solutions[s_i]]].push_back(solution_to_v[solutions[s_j]]);
            }
        }
    }

    return dg_adj;
}

void _strongly_connected_components_dfs1(int v, const vector<vector<int>> & dg_adj, vector<int> & dg_dfs_closed, int & dg_dfs_closed_num) {
    dg_dfs_closed[v] = 0;
    for(int nv : dg_adj[v]) {
        if(dg_dfs_closed[nv]!=-1) continue;
        _strongly_connected_components_dfs1(nv, dg_adj, dg_dfs_closed, dg_dfs_closed_num);
    }
    dg_dfs_closed[v] = dg_dfs_closed_num++;
}

void _strongly_connected_components_dfs2(int v, const vector<vector<int>> & dg_adj, vector<int> & dg_dfs_components, int dg_dfs_components_num) {
    dg_dfs_components[v] = dg_dfs_components_num;
    for(int nv : dg_adj[v]) {
        if(dg_dfs_components[nv]!=-1) continue;
        _strongly_connected_components_dfs2(nv, dg_adj, dg_dfs_components, dg_dfs_components_num);
    }
}

vector<vector<int>> compute_all_dominance_free_solutions(const vector<vector<int>> & solutions, const _Graph & g, const vector<int> & red_vertices) {
    int dg_n = solutions.size();
    vector<vector<int>> dg_adj(dg_n);
    map<vector<int>, int> solution_to_v;
    for(int i = 0; i < solutions.size(); ++i) {
        solution_to_v[solutions[i]] = i;
    }

    dg_adj = construct_dominance_graph(solutions, g, red_vertices);

    vector<int> dg_dfs_closed(dg_n, -1);
    int dg_dfs_closed_num = 0;
    for(int i = 0; i < dg_n; ++i) {
        if(dg_dfs_closed[i]==-1) {
            _strongly_connected_components_dfs1(i, dg_adj, dg_dfs_closed, dg_dfs_closed_num);
        }
    }
    vector<vector<int>> dg_adj_reversed(dg_n);
    for(int i = 0; i < dg_n; ++i) {
        for(int nv : dg_adj[i]) {
            dg_adj_reversed[nv].push_back(i);
        }
    }
    vector<int> dg_dfs_closed_order(dg_n);
    for(int i = 0; i < dg_n; ++i) {
        dg_dfs_closed_order[dg_n-1-dg_dfs_closed[i]] = i;
    }
    vector<int> dg_dfs_components(dg_n, -1);
    int dg_dfs_components_num = 0;
    for(int v : dg_dfs_closed_order) {
        if(dg_dfs_components[v]==-1) {
            _strongly_connected_components_dfs2(v, dg_adj_reversed, dg_dfs_components, dg_dfs_components_num++);
        }
    }
    int cdg_n = dg_dfs_components_num;
    vector<vector<int>> cdg_adj(cdg_n);
    std::set<pair<int,int>> contraction_edges;
    for(int i = 0; i < dg_n; ++i) {
        for(int nv : dg_adj[i]) {
            int ci = dg_dfs_components[i];
            int cj = dg_dfs_components[nv];
            if(ci==cj)continue;
            if(!contraction_edges.count({ci, cj})) {
                contraction_edges.insert({ci, cj});
                cdg_adj[ci].push_back(cj);
            }
        }
    }

    vector<vector<int>> csink_solutions;

    for(int i = 0; i < cdg_n; ++i) {
        if(cdg_adj[i].size()==0) {
            for(int j = 0; j < dg_n; ++j) {
                if(dg_dfs_components[j]==i) {
                    csink_solutions.push_back(solutions[j]);
                    break;
                }
            }
        }
    }

    return csink_solutions;
}
