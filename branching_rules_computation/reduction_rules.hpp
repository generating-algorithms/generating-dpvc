#pragma once
#include "branching_rules_primitives.hpp"

bool is_g_rv_handled_by_red_component_reduction(const _Graph & g, const vector<int> & red_vertices) {
    vector<bool> red_vertices_mask(g.n, false);
    for(int v : red_vertices) red_vertices_mask[v] = true;
    auto adj_lists = g.get_adjacent_lists();

    vector<bool> closed(g.n, false);
    vector<vector<int>> red_components;
    for(int rv : red_vertices) {
        if(closed[rv] == false) {
            vector<int> rc;
            queue<int> q;
            rc.push_back(rv);
            q.push(rv);
            closed[rv] = true;
            while(q.size()) {
                int v = q.front();
                q.pop();
                for(int u : adj_lists[v]) {
                    if(red_vertices_mask[u] == true && closed[u] == false) {
                        rc.push_back(u);
                        q.push(u);
                        closed[u] = true;
                    }
                }
            }

            red_components.push_back(std::move(rc));
        }
    }

    vector<bool> star_mask(g.n, false);

    for(const auto & red_component : red_components) {
        // To avoid expensive check that the red_component is d-path free, we just check if it smaller than d.
        if(red_component.size() >= DPVC_PATH_LEN) continue;

        std::set<int> red_component_blue_neighborhood;
        for(int v : red_component) {
            for(int u : adj_lists[v]) {
                if(red_vertices_mask[u] == false) {
                    red_component_blue_neighborhood.insert(u);
                }
            }
        }

        if(red_component_blue_neighborhood.size() == 0) {
            return false;
        }
        // If there are two red components with the same blue_neighborhood which is just one vertex, it is handled.
        if(red_component_blue_neighborhood.size() == 1) {
            int star_center = *red_component_blue_neighborhood.begin();
            if(star_mask[star_center]) return true;
            star_mask[star_center] = true;
        }

    }

    return false;
}


bool is_g_rv_handled_by_red_star_reduction(const _Graph & g, const vector<int> & red_vertices) {
    if(DPVC_PATH_LEN < 4) {
        return false;
    }

    std::map<int, int> kcenter_red_rays_count;
    for(int v : red_vertices) {
        int kcenter = g.get_adjacent_vertices_mask(v);
        kcenter_red_rays_count[kcenter]++;
    }

    for(const pair<int, int> & krrc : kcenter_red_rays_count) {
        int kcenter_size = __builtin_popcount(krrc.first);
        int red_rays_count = krrc.second;

        if(DPVC_PATH_LEN >= 4) {
            int kcenter_size_threshold = DPVC_PATH_LEN / 2 - 1;
            if(kcenter_size <= kcenter_size_threshold) {
                if(red_rays_count >= 2*kcenter_size) return true;
            }
        }
    }

    return false;
}

