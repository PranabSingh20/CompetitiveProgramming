template<class T> struct SegmentTree {
public : const T ID {}; vector<T> t; int n;
    T operation(T a, T b) {return a + b;}
    SegmentTree(const vector<T> &v) {
        n = (int)v.size();
        t.assign(2 * n, ID);
        for (int i = 0; i < n; ++i)t[n + i] = v[i];
        for (int i = n - 1; i > 0; --i) t[i] = operation(t[i << 1], t[i << 1 | 1]);
    }
    void modify(int p, T val) {
        for (t[p += n] = val; p > 1; p >>= 1) t[p >> 1] = operation(t[p], t[p ^ 1]);
    }
    T query(int l, int r) {
        if (l == r)return t[n + l];
        T ra = t[n + l], rb = t[n + r];
        l++, r--;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) ra = operation(ra, t[l++]);
            if (r & 1) rb = operation(t[--r], rb);
        }
        return operation(ra, rb);
    }
};