template<class T> struct LazySegmentTree {
public : const T ID {}; vector<T> t, d; int n, h;
    T operation(T a, T b) {return a + b;}
    LazySegmentTree(const vector<T> &v) {
        n = (int)v.size();
        h = sizeof(int) * 8 - __builtin_clz(n);
        t.assign(2 * n, ID), d.assign(n, ID);
        for (int i = 0; i < n; ++i)t[n + i] = v[i];
        for (int i = n - 1; i > 0; --i) t[i] = operation(t[i << 1], t[i << 1 | 1]);
    }
    void apply(int p, T value) {
        t[p] += value;
        if (p < n) d[p] += value;
    }
    void build(int p) {
        while (p > 1) p >>= 1, t[p] = operation(t[p << 1], t[p << 1 | 1]) + d[p];
    }
    void push(int p) {
        for (int s = h; s > 0; --s) {
            int i = p >> s;
            if (d[i] != 0) {
                apply(i << 1, d[i]);
                apply(i << 1 | 1, d[i]);
                d[i] = 0;
            }
        }
    }
    void modify(int l, int r, T value) {
        l += n, r += n + 1;
        int l0 = l, r0 = r;
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) apply(l++, value);
            if (r & 1) apply(--r, value);
        }
        build(l0), build(r0 - 1);
    }
    T query(int l, int r) {
        push(l + n), push(r + n);
        if (l == r)return t[l + n];
        T ra = t[n + l++], rb = t[n + r--];
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) ra = operation(ra, t[l++]);
            if (r & 1) rb = operation(t[--r], rb);
        }
        return operation(ra, rb);
    }
};