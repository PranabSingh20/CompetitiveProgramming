#include "bits/stdc++.h"
using namespace std;
using namespace chrono;
typedef int64_t ll;
typedef long double ld;
typedef unsigned long long ULL;
#define endl "\n"
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define pb push_back
void read(vector<int> &a) {for (auto &x : a)cin >> x;}
void read(vector<ll> &a) {for (auto &x : a)cin >> x;}
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
#define sim template < class c
#define ris return * this
#define dor > debug & operator <<
#define eni(x) sim > typename \
  enable_if<sizeof dud<c>(0) x 1, debug&>::type operator<<(c i) {
sim > struct rge { c b, e; };
sim > rge<c> range(c i, c j) { return rge<c> {i, j}; }
sim > auto dud(c * x) -> decltype(cerr << *x, 0);
sim > char dud(...);
struct debug {
#ifndef ONLINE_JUDGE
    ~debug() { cerr << endl; }
    eni( != ) cerr << boolalpha << i; ris;
}
eni( == ) ris << range(begin(i), end(i));
}
sim, class b dor(pair < b, c > d) {
    ris << "(" << d.first << ", " << d.second << ")";
}
sim dor(rge<c> d) {
    *this << "[";
    for (auto it = d.b; it != d.e; ++it)
        *this << ", " + 2 * (it == d.b) << *it;
    ris << "]";
}
#else
    sim dor(const c&) { ris; }
#endif
};
#define imie(...) " [" << #__VA_ARGS__ ": " << (__VA_ARGS__) << "] "

const int MOD = 1e9 + 7;
const int INF = (int)2e9 + 7;
const ll LINF = (ll)1e18;
const ld PI = 3.1415926535897932384626433832795;
template< int mod >
struct NumberTheoreticTransform {

    vector< int > rev, rts;
    int base, max_base, root;

    NumberTheoreticTransform() : base(1), rev{0, 1}, rts{0, 1} {
        assert(mod >= 3 && mod % 2 == 1);
        auto tmp = mod - 1;
        max_base = 0;
        while (tmp % 2 == 0) tmp >>= 1, max_base++;
        root = 2;
        while (mod_pow(root, (mod - 1) >> 1) == 1) ++root;
        assert(mod_pow(root, mod - 1) == 1);
        root = mod_pow(root, (mod - 1) >> max_base);
    }

    inline int mod_pow(int x, int n) {
        int ret = 1;
        while (n > 0) {
            if (n & 1) ret = mul(ret, x);
            x = mul(x, x);
            n >>= 1;
        }
        return ret;
    }

    inline int inverse(int x) {
        return mod_pow(x, mod - 2);
    }

    inline unsigned add(unsigned x, unsigned y) {
        x += y;
        if (x >= mod) x -= mod;
        return x;
    }

    inline unsigned mul(unsigned a, unsigned b) {
        return 1ull * a * b % (unsigned long long) mod;
    }

    void ensure_base(int nbase) {
        if (nbase <= base) return;
        rev.resize(1 << nbase);
        rts.resize(1 << nbase);
        for (int i = 0; i < (1 << nbase); i++) {
            rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
        }
        assert(nbase <= max_base);
        while (base < nbase) {
            int z = mod_pow(root, 1 << (max_base - 1 - base));
            for (int i = 1 << (base - 1); i < (1 << base); i++) {
                rts[i << 1] = rts[i];
                rts[(i << 1) + 1] = mul(rts[i], z);
            }
            ++base;
        }
    }


    void ntt(vector< int > &a) {
        const int n = (int) a.size();
        assert((n & (n - 1)) == 0);
        int zeros = __builtin_ctz(n);
        ensure_base(zeros);
        int shift = base - zeros;
        for (int i = 0; i < n; i++) {
            if (i < (rev[i] >> shift)) {
                swap(a[i], a[rev[i] >> shift]);
            }
        }
        for (int k = 1; k < n; k <<= 1) {
            for (int i = 0; i < n; i += 2 * k) {
                for (int j = 0; j < k; j++) {
                    int z = mul(a[i + j + k], rts[j + k]);
                    a[i + j + k] = add(a[i + j], mod - z);
                    a[i + j] = add(a[i + j], z);
                }
            }
        }
    }


    vector< int > multiply(vector< int > a, vector< int > b) {
        int need = a.size() + b.size() - 1;
        int nbase = 1;
        while ((1 << nbase) < need) nbase++;
        ensure_base(nbase);
        int sz = 1 << nbase;
        a.resize(sz, 0);
        b.resize(sz, 0);
        ntt(a);
        ntt(b);
        int inv_sz = inverse(sz);
        for (int i = 0; i < sz; i++) {
            a[i] = mul(a[i], mul(b[i], inv_sz));
        }
        reverse(a.begin() + 1, a.end());
        ntt(a);
        a.resize(need);
        return a;
    }

};
void solve() {
    int n, m;
    const int N = 1e6;
    n = rng() % N, m = rng() % N;
    vector<int> a(n), b(m);
    for (int i = 0; i < n; i++) {
        a[i] = rng() % (int)1e9;
    }
    for (int i = 0; i < m; i++) {
        b[i] = rng() % (int)1e9;
    }
    auto t1 = high_resolution_clock::now();
    NumberTheoreticTransform<998244353> ntt;
    vector<int> c = ntt.multiply(a, b);
    auto t2 = high_resolution_clock::now();
    auto time = duration_cast<duration<double>>(t2 - t1);
    cout << "Time elapsed = " << time.count() << endl;
}

int32_t main() {
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int test_cases = 1;
    // cin >> test_cases;
    for (int tt = 1; tt <= test_cases; tt++)
    {
        //cout << "Case #" << tt << ": ";
        solve();
    }
}

