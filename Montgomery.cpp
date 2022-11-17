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
struct NTRU {

    int max_base, root;

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

    vector<int> Montgomery(vector<int> a) {
        int n = a.size();
        vector<int> b = a;
        for (int i = 1; i < n; i++) {
            b[i] = mul(b[i - 1], a[i]);
        }
        vector<int> t(n), c(n);
        t[n - 1] = inverse(b[n - 1]);
        for (int i = n - 1; i > 0; i--) {
            c[i] = mul(t[i], b[i - 1]);
            t[i - 1] = mul(t[i], a[i]);
        }
        c[0] = t[0];
        return c;
    }

};

void solve() {
    int n;
    cin >> n ;
    vector<int> a(n);
    read(a);
    NTRU<998244353> ntru;
    vector<int> inv = ntru.Montgomery(a);
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
    auto t1 = high_resolution_clock::now();
    for (int tt = 1; tt <= test_cases; tt++)
    {
        //cout << "Case #" << tt << ": ";
        solve();
    }
    auto t2 = high_resolution_clock::now();
    auto time = duration_cast<duration<double>>(t2 - t1);
    cerr << "Time elapsed = " << time.count() << endl;
}

