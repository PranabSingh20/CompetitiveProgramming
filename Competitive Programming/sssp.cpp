//https://judge.yosupo.jp/submission/52076
template <typename T>
using min_heap = std::priority_queue<T, std::vector<T>, std::greater<T>>;
template <class T, class G>
auto dijkstra(const G& g, int s) {
    std::vector d(g.size(), std::numeric_limits<T>::max());
    std::vector prv(g.size(), -1);
    min_heap<std::pair<T, int>> rh;
    rh.emplace(d[s] = 0, s);
    while (!rh.empty()) {
        auto [dv, v] = rh.top();
        rh.pop();
        if (dv != d[v]) continue;
        for (auto && [to, w] : g[v]) {
            if (d[to] > dv + w) {
                d[to] = dv + w;
                rh.emplace(d[to], to);
                prv[to] = v;
            }
        }
    }
    return std::make_pair(d, prv);
}

void solve() {
    int n, m, s, t;
    cin >> n >> m >> s >> t;
    vector<vector<pair<int, int>>> g(n);
    while (m--) {
        int u, v, w;
        cin >> u >> v >> w; //0-based
        g[u].emplace_back(v, w);
    }
    auto [d, prv] = dijkstra<int64_t>(g, s);
    if (d[t] == numeric_limits<int64_t>::max()) {
        cout << -1 << '\n';
    } else {
        vector<pair<int, int>> res;
        for (int cur = t; cur != s; cur = prv[cur])
            res.emplace_back(prv[cur], cur);
        reverse(begin(res), end(res));
        cout << d[t] << ' ' << res.size() << '\n';
        for (auto [a, b] : res) cout << a << ' ' << b << '\n';
    }
}