//https://judge.yosupo.jp/submission/52004
struct GraphEdgePointers {
    struct edge {
        int to, nxt;
        edge(int to, int nxt) : to(to), nxt(nxt) {}
    };
    std::vector<int> head;
    std::vector<edge> edges;
    int cur_edges;
    GraphEdgePointers(int n, int m) : head(n + 1, 0), cur_edges(1) {
        edges.reserve(m + 1);
        edges.emplace_back(0, 0);
    }
    void add_edge(int u, int v) {
        edges.emplace_back(v, head[u]);
        head[u] = cur_edges++;
    }
};

struct GeneralMatching {
    int n, m;
    GraphEdgePointers g;
    int n_matches, q_n, book_mark;
    std::vector<int> mate, q, book, type, fa, bel;

    GeneralMatching(int n, int m)
        : n(n),
          m(m),
          g(n, 2 * m),
          n_matches(0),
          q_n(0),
          book_mark(0),
          mate(n + 1),
          q(n + 1),
          book(n + 1),
          type(n + 1),
          fa(n + 1),
          bel(n + 1) {}

    void add_edge(int u, int v) {
        g.add_edge(u, v);
        g.add_edge(v, u);
    }

    int maximum_matching() {
        n_matches = 0;
        for (int u = 1; u <= n; ++u)
            if (!mate[u] && match(u)) ++n_matches;
        return n_matches;
    }

    std::vector<std::pair<int, int>> get_edges() {
        std::vector<std::pair<int, int>> ans;
        for (int u = 1; u <= n; ++u)
            if (mate[u] > u) ans.emplace_back(u, mate[u]);
        return ans;
    }

private:
    void augment(int u) {
        while (u) {
            int nu = mate[fa[u]];
            mate[mate[u] = fa[u]] = u;
            u = nu;
        }
    }
    int get_lca(int u, int v) {
        ++book_mark;
        while (true) {
            if (u) {
                if (book[u] == book_mark) return u;
                book[u] = book_mark;
                u = bel[fa[mate[u]]];
            }
            std::swap(u, v);
        }
    }
    void go_up(int u, int v, const int& mv) {
        while (bel[u] != mv) {
            fa[u] = v;
            v = mate[u];
            if (type[v] == 1) type[q[++q_n] = v] = 0;
            if (bel[u] == u) bel[u] = mv;
            if (bel[v] == v) bel[v] = mv;
            u = fa[v];
        }
    }
    void after_go_up() {
        for (int u = 1; u <= n; ++u) bel[u] = bel[bel[u]];
    }
    bool match(const int& sv) {
        for (int u = 1; u <= n; ++u) bel[u] = u, type[u] = -1;
        type[q[q_n = 1] = sv] = 0;
        for (int i = 1; i <= q_n; ++i) {
            int u = q[i];
            for (int e = g.head[u]; e; e = g.edges[e].nxt) {
                int v = g.edges[e].to;
                if (!~type[v]) {
                    fa[v] = u, type[v] = 1;
                    int nu = mate[v];
                    if (!nu) {
                        augment(v);
                        return true;
                    }
                    type[q[++q_n] = nu] = 0;
                } else if (!type[v] && bel[u] != bel[v]) {
                    int lca = get_lca(u, v);
                    go_up(u, v, lca);
                    go_up(v, u, lca);
                    after_go_up();
                }
            }
        }
        return false;
    }
};

void solve() {
    int n, m;
    cin >> n >> m;
    GeneralMatching matching(n, m);
    while (m--) {
        int u, v;
        cin >> u >> v;
        ++u, ++v;
        matching.add_edge(u, v);
    }
    cout << matching.maximum_matching() << endl;
    auto edges = matching.get_edges();
    for (auto e : edges) cout << e.first << ' ' << e.second << '\n';
}