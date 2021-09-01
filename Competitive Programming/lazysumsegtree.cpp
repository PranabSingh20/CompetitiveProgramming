//https://github.com/NyaanNyaan/library/blob/master/data-structure/range-sum-range-add-bit.hpp
template <typename T>
struct BinaryIndexedTree {
    int N;
    vector<T> data;

    BinaryIndexedTree() = default;

    BinaryIndexedTree(int size) { init(size); }

    void init(int size) {
        N = size + 2;
        data.assign(N + 1, 0);
    }

    // get sum of [0,k]
    T sum(int k) const {
        if (k < 0) return 0;  // return 0 if k < 0
        T ret = 0;
        for (++k; k > 0; k -= k & -k) ret += data[k];
        return ret;
    }

    // getsum of [l,r]
    inline T sum(int l, int r) const { return sum(r) - sum(l - 1); }

    // get value of k
    inline T operator[](int k) const { return sum(k) - sum(k - 1); }

    // data[k] += x
    void add(int k, T x) {
        for (++k; k < N; k += k & -k) data[k] += x;
    }
};

template <typename T>
struct RangeAddRangeSumBIT {
    BinaryIndexedTree<T> a, b;
    RangeAddRangeSumBIT(const vector<T>&v) : a((int)v.size() + 1), b((int)v.size() + 1) {
        for (int i = 0; i < (int)v.size(); ++i) {
            a.add(i, v[i]);
            a.add(i + 1, -v[i]);
            b.add(i, v[i] * (1 - i));
            b.add(i + 1, v[i] * i);
        }
    }

    void add(int l, int r, T x) {
        a.add(l, x);
        a.add(r, -x);
        b.add(l, x * (1 - l));
        b.add(r, x * (r - 1));
    }

    // return sum of [l, r)
    T query(T l, T r) {
        --r, --l;
        return a.sum(r) * r + b.sum(r) - a.sum(l) * l - b.sum(l);
    }
};