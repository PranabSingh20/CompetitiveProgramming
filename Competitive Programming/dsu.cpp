struct DSU {
    //nodes are 1 base
    int SetCnt, FirstSetCnt;
    vector<int> par, sz;

    //constructor
    DSU (const int &x) {
        SetCnt = x;
        FirstSetCnt = x;
        par.resize(x + 1);
        sz.resize(x + 1);
        fill(par.begin(), par.end(), -1);
        fill(sz.begin(), sz.end(), 1);
    }
    //finding root of a vertex
    int find_root(int x) {
        return (par[x] == -1 ? x : par[x] = find_root(par[x]));
    }

    //merging two component
    bool setunion(int x, int y) {
        int p1 = find_root(x), p2 = find_root(y);
        SetCnt -= (p1 != p2);
        if (sz[p1] > sz[p2])
            swap(p1, p2);
        if (p1 != p2)
            sz[p2] += sz[p1];
        return (p1 == p2 ? false : (par[p1] = p2) | 1);
    }

    //meging two DSU (adding every edge in second DSU to first DSU) ---> O(n.log*(n))
    bool merge(DSU &x) {
        if (x.size() != FirstSetCnt)
            return false;
        for (int i = 1; i <= FirstSetCnt; i++) {
            int p1 = find_root(i), p2 = x.find_root(i);
            if (p1 != p2)
                setunion(p1, p2);
        }
        return true;
    }

    //clear DSU
    void clear() {
        fill(par.begin(), par.end(), -1);
        fill(sz.begin(), sz.end(), 1);
        SetCnt = FirstSetCnt;
    }

    //number of nodes in component of node x
    int size_cmp(const int &x) {
        return sz[find_root(x)];
    }

    //number of components
    int num_cmp() {
        return SetCnt;
    }

    //number of nodes
    int size() {
        return FirstSetCnt;
    }
};