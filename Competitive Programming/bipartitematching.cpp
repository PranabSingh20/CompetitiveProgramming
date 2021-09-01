//https://judge.yosupo.jp/submission/52092
// consider shuffling adjacency lists if this gives TLE
template <bool ToShuffle = false>
struct bipartite_matching {
    int n_left, n_right;
    std::vector<std::vector<int>> g;
    std::vector<int> match_from_left, match_from_right;
    std::vector<int> vis;
    int iteration;

    bipartite_matching(int _n_left, int _n_right)
        : n_left(_n_left),
          n_right(_n_right),
          g(_n_left),
          vis(_n_left, 0),
          match_from_left(_n_left, -1),
          match_from_right(_n_right, -1),
          iteration(0) {}

    // u on left, v on right
    void add(int u, int v) { g[u].push_back(v); }

    bool dfs(int u) {
        vis[u] = iteration;
        for (auto v : g[u]) {
            if (match_from_right[v] == -1) {
                match_from_right[v] = u;
                match_from_left[u] = v;
                return true;
            }
        }
        for (auto v : g[u]) {
            if (vis[match_from_right[v]] != iteration &&
                    dfs(match_from_right[v])) {
                match_from_right[v] = u;
                match_from_left[u] = v;
                return true;
            }
        }
        return false;
    }

    int get_max_matching() {
        if constexpr (ToShuffle) {
            for (int i = 0; i < n_left; ++i)
                std::random_shuffle(std::begin(g[i]), std::end(g[i]));
        }
        int new_matchings = 0;
        int matchings = 0;
        do {
            iteration++;
            new_matchings = 0;
            for (int u = 0; u < n_left; ++u)
                if (match_from_left[u] == -1 && dfs(u)) new_matchings++;
            matchings += new_matchings;
        } while (new_matchings > 0);
        return matchings;
    }

    std::vector<std::pair<int, int>> get_edges() {
        std::vector<std::pair<int, int>> ans;
        for (int i = 0; i < n_left; ++i)
            if (match_from_left[i] != -1)
                ans.emplace_back(i, match_from_left[i]);
        return ans;
    }
};

void solve() {
    int n_left, n_right, m;
    cin >> n_left >> n_right >> m;
    bipartite_matching<true> matching(n_left, n_right);
    while (m--) {
        int u, v;
        cin >> u >> v;
        matching.add(u, v); //0-based
    }
    cout << matching.get_max_matching() << endl;
    auto edges = matching.get_edges();
    for (auto e : edges) cout << e.first << ' ' << e.second << '\n';
}