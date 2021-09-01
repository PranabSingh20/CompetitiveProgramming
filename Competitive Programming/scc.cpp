//https://judge.yosupo.jp/submission/52102
auto FindSCC(const std::vector<std::vector<int>> &G, int const Base = 0) {
    std::vector<std::vector<int>> SCC;
    std::vector<unsigned> S, P, depth(G.size());
    const auto dfs = [&](const auto & self, auto u) -> void {
        auto d = S.size();
        S.push_back(u);
        P.push_back(d + 1);
        depth[u] = d + 1;
        for (auto v : G[u]) {
            if (!depth[v])
                self(self, v);
            else
                while (P.back() > depth[v]) P.pop_back();
        }
        if (P.back() > d) {
            SCC.emplace_back(S.begin() + d, S.end());
            for (auto v : SCC.back()) depth[v] = -1;
            S.erase(S.begin() + d, S.end());
            P.pop_back();
        }
    };
    for (auto u = Base; u < G.size(); ++u)
        if (!depth[u]) dfs(dfs, u);
    reverse(SCC.begin(), SCC.end());
    return SCC;
}

void solve() {
    int n, m;
    cin >> n >> m;
    vector<vector<int>> g(n);
    while (m--) {
        int u, v;
        cin >> u >> v;
        g[u].push_back(v);
    }
    auto scc = FindSCC(g);
    cout << scc.size() << '\n';
    for (const auto& component : scc) {
        cout << component.size() << ' ';
        for (auto x : component) cout << x << ' ';
        cout << '\n';
    }
}