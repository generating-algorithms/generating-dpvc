#include <chrono>
#include "../common_primitives/primitive_operations.hpp"

struct _ForbiddenInducedSubgraphs {
    std::set<_NautyCanonicalRep> canonical_reps;
    vector<int> sizes;
    int min_size;
    int max_size;
    void add(_NautyCanonicalRep canonical_rep, int s) {
        canonical_reps.insert(canonical_rep);
        if(sizes.size()==0) {
            min_size=s;
            max_size=s;
        }
        min_size = std::min(min_size, s);
        max_size = std::max(max_size, s);
        for(int ss : sizes){if(ss==s){s=-1;}}
        if(s!=-1)sizes.push_back(s);
    }
    void add(const _Graph & g) {
        assert(g.canonical_rep.first);
        add(g.canonical_rep, g.n);
    }
    void clear() {
        canonical_reps.clear();
        sizes.clear();
    }
};

bool is_g_containing_forbidden_induced_subgraph(const _Graph & g, const _ForbiddenInducedSubgraphs & fis, _Nauty & nt) {
    for(int s_idx = fis.sizes.size()-1; s_idx >= 0; --s_idx) {
        int s = fis.sizes[s_idx];
        if(s>g.n) continue;

        if(g.n == s) {
            if(fis.canonical_reps.count(g.canonical_rep)) {
                return true;
            }
        }
        else {
            for(int vs_mask : ordered_combinations_of_bits(g.n, s)) {
                auto edges = g.get_induced_edges_list(vs_mask);

                _NautyCanonicalRep g_prime_canonical = nt.get_canonical_rep(s, edges);
                if(fis.canonical_reps.count(g_prime_canonical)) {
                    return true;
                }
            }
        }
    }
    return false;
}

// OPTIMIZATION:
// The red_vertices are propagated between generations. No expansion will be connected to some red_vertex, so we may skip them.
vector<_Graph> expand_graph(const _Graph & g, _Nauty & nt) {
    vector<int> non_red_vertices_edges;
    for(int v = 0; v < g.n; ++v) {
        if(!(g.red_vertices_mask&(1<<v))){
            non_red_vertices_edges.push_back(v);
        }
    }
    vector<_Graph> egs;
    for(const auto & rvs : ordered_powerset_nonempty(non_red_vertices_edges)) {
        _Graph eg(g.n+1, g.edges, make_empty_nauty_canonical_rep(), g.red_vertices_mask);
        for(int v : rvs) {
            eg.set_edge(v, g.n);
        }
        eg.canonical_rep = nt.get_canonical_rep(g.n+1, eg.get_edges_list());
        egs.push_back(eg);
    }

    return egs;
}

vector<_Graph> filter_graphs_containing_forbidden_induced_subgraph(const vector<_Graph> & gs, const _ForbiddenInducedSubgraphs & fis, _Nauty & nt) {
    vector<_Graph> filtered_gs;
    std::map<_NautyCanonicalRep, bool> decided_graphs;
    for(const _Graph & g : gs) {
        if(decided_graphs.count(g.canonical_rep)) {
            if(decided_graphs[g.canonical_rep]) {
                filtered_gs.push_back(g);
            }
        }
        else {
            if(!is_g_containing_forbidden_induced_subgraph(g, fis, nt)) {
                filtered_gs.push_back(g);
                decided_graphs[g.canonical_rep] = true;
            }
            else {
                decided_graphs[g.canonical_rep] = false;
            }
        }
    }
    return filtered_gs;
}

void parallel_filter_graphs_containing_forbidden_induced_subgraph(vector<vector<_Graph>> & gss_to_filter, const _ForbiddenInducedSubgraphs & fis, vector<_Nauty> & nts) {

    #pragma omp parallel default(none) shared(gss_to_filter, fis, nts)
    {
        _Nauty & nt = nts[omp_get_thread_num()];
        #pragma omp for schedule(dynamic)
        for(int gss_idx = 0; gss_idx < gss_to_filter.size(); ++gss_idx) {
            gss_to_filter[gss_idx] = filter_graphs_containing_forbidden_induced_subgraph(gss_to_filter[gss_idx], fis, nt);
        }
    }
}

// Implementation of https://link.springer.com/chapter/10.1007/978-3-319-21840-3_49 602  - enumerating all induced connected subgraphs.
bool is_g_containing_forbidden_induced_subgraph_containing_v(int v, const _Graph & g, const _ForbiddenInducedSubgraphs & fis, _Nauty & nt, std::map<std::pair<int, _GraphEdges>, bool> & cache) {
    vector<uint32_t> adj_masks(g.n, 0);
    for(int i = 0; i < g.n; i++) {
        for(int j = i + 1 ; j < g.n; j++) {
            if(g.has_edge(i, j)){
                adj_masks[i] |= 1<<j;
                adj_masks[j] |= 1<<i;
            }
        }
    }
    assert(v<g.n);

    vector<std::tuple<uint32_t, uint32_t, uint32_t>> _stack;

    // cur_ics, cur_closed, cur_neighborhood,
    _stack.push_back(std::make_tuple(((uint32_t)1)<<v, ((uint32_t)1)<<v, adj_masks[v]));

    while(_stack.size()) {
        auto _frame = _stack.back();
        _stack.pop_back();

        const uint32_t cur_ics = std::get<0>(_frame);
        const uint32_t cur_closed = std::get<1>(_frame);
        const uint32_t cur_neighborhood = std::get<2>(_frame);
        int s = __builtin_popcount(cur_ics);

        if(cur_neighborhood == 0) {
            if(fis.min_size<=s && s <= fis.max_size) {
                // OPTIMIZATION:
                // As the computation of the canonical rep is costly, do fast identity check first to see if we have already looked at this one.
                std::pair<int, _GraphEdges> cache_key(s, g.get_induced_edges_mask(cur_ics));
                auto cache_entry = cache.find(cache_key);
                if(cache_entry != cache.end()) {
                    if(cache_entry->second) return true;
                }
                else {
                    const auto & g_prime_canonical = nt.get_canonical_rep(s, g.get_induced_edges_list(cur_ics));
                    if(fis.canonical_reps.count(g_prime_canonical)) {
                        cache[cache_key] = true;
                        return true;
                    }
                    cache[cache_key] = false;
                }
            }
            continue;
        }
        if(s>fis.max_size) continue;


        const int nv = __builtin_ctz(cur_neighborhood);

        const uint32_t new_cur_ics_1 = cur_ics | (1<<nv);
        const uint32_t new_cur_closed_1 = cur_closed | (1<<nv);
        const uint32_t new_cur_neighborhood_1 = (cur_neighborhood | adj_masks[nv]) & ~new_cur_closed_1;

        _stack.push_back(std::make_tuple(new_cur_ics_1, new_cur_closed_1, new_cur_neighborhood_1));

        const uint32_t new_cur_ics_2 = cur_ics;
        const uint32_t new_cur_closed_2 = cur_closed | (1<<nv);
        const uint32_t new_cur_neighborhood_2 = cur_neighborhood & ~new_cur_closed_1;

        _stack.push_back(std::make_tuple(new_cur_ics_2, new_cur_closed_2, new_cur_neighborhood_2));
    }

    return false;
}

vector<_Graph> expand_graph_and_filter_forbidden_subgraphs(const _Graph & g,  const _ForbiddenInducedSubgraphs & fis, _Nauty & nt) {
    auto egs = expand_graph(g, nt);

    std::map<std::pair<int, _GraphEdges>, bool> cache;

    vector<_Graph> filtered_gs;
    std::map<_NautyCanonicalRep, bool> decided_graphs;
    for(const _Graph & eg : egs) {
        if(decided_graphs.count(eg.canonical_rep)) {
            if(decided_graphs[eg.canonical_rep]) {
                filtered_gs.push_back(eg);
            }
        }
        else {
            // OPTIMIZATION:
            // We know that g does not contain any forbidden induced subgraph (because it would be filtered earlier on).
            // Therefore we only need to test all connected induced subgraphs of the expansion that contain the new expanding vertex.
            if(!is_g_containing_forbidden_induced_subgraph_containing_v(g.n, eg, fis, nt, cache)) {
                filtered_gs.push_back(eg);
                decided_graphs[eg.canonical_rep] = true;
            }
            else {
                decided_graphs[eg.canonical_rep] = false;
            }
        }
    }
    return filtered_gs;
}

vector<vector<_Graph>> parallel_expand_graph_without_forbidden_subgraph(const vector<_Graph> & gs, const _ForbiddenInducedSubgraphs & fis, vector<_Nauty> & nts) {
    vector<vector<_Graph>> egs(gs.size());
    int _progress = 0;
    int _progress_s = (gs.size()+100-1) / 100;
    auto _progress_clock = std::chrono::steady_clock::now();
    int64_t _progress_duration = 0;
    cerr<<"-- expansion------------------------------------->"<<gs.size()<<endl;
    #pragma omp parallel default(none) shared(gs, fis, nts, egs, _progress, _progress_s, _progress_clock, _progress_duration, cerr)
    {
        _Nauty & nt = nts[omp_get_thread_num()];
        #pragma omp for schedule(dynamic)
        for(int gs_idx = 0; gs_idx < gs.size(); ++gs_idx) {
            egs[gs_idx] = expand_graph_and_filter_forbidden_subgraphs(gs[gs_idx], fis, nt);
            #pragma omp atomic update
            _progress ++;
            if((_progress % _progress_s) == 0){
                #pragma omp critical
                {
                _progress_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - _progress_clock).count();
                int64_t eta_nanos = (gs.size() - _progress) * (_progress_duration / _progress);
                cerr<<eta_nanos/1000/1000/1000<<"s ";
                _progress_clock = std::chrono::steady_clock::now();
                }
            }
        }
    }
    cerr<<endl;
    return egs;
}

vector<_Graph> gather_nonisomorphic_graphs(const vector<vector<_Graph>> & egs) {
    vector<_Graph> all_graphs;
    for(const vector<_Graph> & eg : egs) {
        all_graphs.insert(all_graphs.end(), eg.begin(), eg.end());
    }
    if(!all_graphs.size()) {
        return {};
    }
    sort(all_graphs.begin(), all_graphs.end());
    vector<_Graph> nonisomorphic_graphs;
    nonisomorphic_graphs.push_back(all_graphs[0]);
    for(int i = 1; i < all_graphs.size(); ++i) {
        if(all_graphs[i-1].canonical_rep != all_graphs[i].canonical_rep) {
            nonisomorphic_graphs.push_back(all_graphs[i]);
        }
    }
    return nonisomorphic_graphs;
}

vector<_Graph> get_connected_graphs_of_size(int n, vector<_Nauty> & nts) {
    vector<_Graph> expanded_graphs;
    _ForbiddenInducedSubgraphs fis;
    expanded_graphs.push_back(_Graph(1, _GraphEdges(), make_empty_nauty_canonical_rep(), 0));
    for(int i = 0; i < n - 1; ++i) {
        expanded_graphs = gather_nonisomorphic_graphs(parallel_expand_graph_without_forbidden_subgraph(expanded_graphs, fis, nts));
    }
    return expanded_graphs;
}
