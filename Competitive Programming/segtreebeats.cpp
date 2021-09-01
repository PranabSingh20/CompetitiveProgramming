//https://judge.yosupo.jp/submission/2700
//usage: https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
//Desc : set (min(ai, x), max(ai, x)), increment ai += x, in [l, r) 0-based
template<class T> struct segment_tree_chminmaxadd_rsq {
    const T ide = 0;
    static const T INF = std::numeric_limits<T>::max();
    static const T NINF = std::numeric_limits<T>::min();
    struct node {
        T sum;
        T gst1, gst2, gcnt;
        T lst1, lst2, lcnt;
        T len;
        T add;
        node(): gst2(NINF), gcnt(1), lst2(INF), lcnt(1), len(1), add(0) {}
    };

    std::vector<node> v;

    int n;
    int h;

    void fix(int k) {
        node& p = v[k];
        node& l = v[k * 2 + 0];
        node& r = v[k * 2 + 1];
        p.sum = l.sum + r.sum;

        if (r.gst1 < l.gst1) {
            p.gst1 = l.gst1;
            p.gcnt = l.gcnt;
            p.gst2 = std::max(l.gst2, r.gst1);
        }
        else if (l.gst1 < r.gst1) {
            p.gst1 = r.gst1;
            p.gcnt = r.gcnt;
            p.gst2 = std::max(l.gst1, r.gst2);
        }
        else {
            p.gst1 = l.gst1;
            p.gcnt = l.gcnt + r.gcnt;
            p.gst2 = std::max(l.gst2, r.gst2);
        }

        if (r.lst1 > l.lst1) {
            p.lst1 = l.lst1;
            p.lcnt = l.lcnt;
            p.lst2 = std::min(l.lst2, r.lst1);
        }
        else if (l.lst1 > r.lst1) {
            p.lst1 = r.lst1;
            p.lcnt = r.lcnt;
            p.lst2 = std::min(l.lst1, r.lst2);
        }
        else {
            p.lst1 = l.lst1;
            p.lcnt = l.lcnt + r.lcnt;
            p.lst2 = std::min(l.lst2, r.lst2);
        }

    }

    segment_tree_chminmaxadd_rsq() {}
    segment_tree_chminmaxadd_rsq(const std::vector<T>& vec) {
        n = 1;
        h = 1;
        while (n < vec.size()) n <<= 1, h++;
        v.resize(2 * n);
        for (int i = 0; i < vec.size(); i++) {
            v[i + n].sum = vec[i];
            v[i + n].gst1 = vec[i];
            v[i + n].lst1 = vec[i];
        }
        for (int i = n; i -- > 1;) {
            fix(i);
            v[i].len = v[i * 2 + 0].len + v[i * 2 + 1].len;
        }
    }

    void eff_add(int k, T x) {
        auto& p = v[k];
        p.sum  += x * p.len;
        p.gst1 += x;
        p.lst1 += x;
        p.add += x;
        if (p.gst2 != NINF) p.gst2 += x;
        if (p.lst2 !=  INF) p.lst2 += x;
    }

    void eff_chmin(int k, T x) {
        auto& p = v[k];
        p.sum += (x - p.gst1) * p.gcnt;
        if (p.gst1 == p.lst1) {
            p.gst1 = p.lst1 = x;
        }
        else if (p.gst1 == p.lst2) {
            p.gst1 = p.lst2 = x;
        }
        else {
            p.gst1 = x;
        }
    }

    void eff_chmax(int k, T x) {
        auto& p = v[k];
        p.sum += (x - p.lst1) * p.lcnt;
        if (p.lst1 == p.gst1) {
            p.lst1 = p.gst1 = x;
        }
        else if (p.lst1 == p.gst2) {
            p.lst1 = p.gst2 = x;
        }
        else {
            p.lst1 = x;
        }
    }

    void push(int k) {
        auto& p = v[k];
        if (p.add != 0) {
            eff_add(k * 2 + 0, p.add);
            eff_add(k * 2 + 1, p.add);
            p.add = 0;
        }
        if (p.gst1 < v[k * 2 + 0].gst1) eff_chmin(k * 2 + 0, p.gst1);
        if (p.gst1 < v[k * 2 + 1].gst1) eff_chmin(k * 2 + 1, p.gst1);

        if (p.lst1 > v[k * 2 + 0].lst1) eff_chmax(k * 2 + 0, p.lst1);
        if (p.lst1 > v[k * 2 + 1].lst1) eff_chmax(k * 2 + 1, p.lst1);
    }

    void infuse(int k) {
        k = k >> __builtin_ctz(k);
        while (k >>= 1) fix(k);
    }

    void infiltrate(int k) {
        if (k == n << 1) return;
        for (int i = h; i -- > 1;) push(k >> i);
    }

    void subtree_chmin(int k, T x) {
        if (v[k].gst1 <= x) return;
        if (v[k].gst2 < x) {
            eff_chmin(k, x);
            return;
        }
        push(k);
        subtree_chmin(k * 2 + 0, x);
        subtree_chmin(k * 2 + 1, x);
        fix(k);
    }

    void subtree_chmax(int k, T x) {
        if (x <= v[k].lst1) return;
        if (x < v[k].lst2) {
            eff_chmax(k, x);
            return;
        }
        push(k);
        subtree_chmax(k * 2 + 0, x);
        subtree_chmax(k * 2 + 1, x);
        fix(k);
    }

    void range_chmin(int a, int b, T x) {
        infiltrate(a + n);
        infiltrate(b + n);
        int l = a + n;
        int r = b + n;
        while (l < r) {
            if (l & 1) subtree_chmin(l++, x);
            if (r & 1) subtree_chmin(--r, x);
            l >>= 1;
            r >>= 1;
        }
        infuse(a + n);
        infuse(b + n);
    }

    void range_chmax(int a, int b, T x) {
        infiltrate(a + n);
        infiltrate(b + n);
        int l = a + n;
        int r = b + n;
        while (l < r) {
            if (l & 1) subtree_chmax(l++, x);
            if (r & 1) subtree_chmax(--r, x);
            l >>= 1;
            r >>= 1;
        }
        infuse(a + n);
        infuse(b + n);
    }

    void range_add(int a, int b, T x) {
        infiltrate(a + n);
        infiltrate(b + n);
        int l = a + n;
        int r = b + n;
        while (l < r) {
            if (l & 1) eff_add(l++, x);
            if (r & 1) eff_add(--r, x);
            l >>= 1;
            r >>= 1;
        }
        infuse(a + n);
        infuse(b + n);
    }

    T range_sum(int l, int r) {
        l += n;
        r += n;
        infiltrate(l);
        infiltrate(r);
        T lx = ide;
        T rx = ide;
        while (l < r) {
            if (l & 1) lx = lx + v[l++].sum;
            if (r & 1) rx = v[--r].sum + rx;
            l >>= 1;
            r >>= 1;
        }
        return lx + rx;
    }
    T range_max(int l, int r) {
        l += n;
        r += n;
        infiltrate(l);
        infiltrate(r);
        T lx = NINF;
        T rx = NINF;
        while (l < r) {
            if (l & 1) lx = max(lx, v[l++].gst1);
            if (r & 1) rx = max(v[--r].gst1, rx);
            l >>= 1;
            r >>= 1;
        }
        return max(lx, rx);
    }
    T range_min(int l, int r) {
        l += n;
        r += n;
        infiltrate(l);
        infiltrate(r);
        T lx = INF;
        T rx = INF;
        while (l < r) {
            if (l & 1) lx = min(lx, v[l++].lst1);
            if (r & 1) rx = min(v[--r].lst1, rx);
            l >>= 1;
            r >>= 1;
        }
        return min(lx, rx);
    }
};