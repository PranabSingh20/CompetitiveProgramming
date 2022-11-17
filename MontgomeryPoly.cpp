/**
 *  date : 2021-05-19 02:14:34
 */

#define NDEBUG
#include <immintrin.h>

//
#include <algorithm>

#include <cassert>

#include <chrono>

#include <cstdint>

#include <cstdio>

#include <cstring>

#include <iostream>

#include <type_traits>

#include <utility>

#include <vector>


using namespace std;
using namespace chrono;

namespace fastio {
static constexpr int SZ = 1 << 17;
char inbuf[SZ], outbuf[SZ];
int in_left = 0, in_right = 0, out_right = 0;

struct Pre {
    char num[40000];
    constexpr Pre() : num() {
        for (int i = 0; i < 10000; i++) {
            int n = i;
            for (int j = 3; j >= 0; j--) {
                num[i * 4 + j] = n % 10 + '0';
                n /= 10;
            }
        }
    }
} constexpr pre;

inline void load() {
    int len = in_right - in_left;
    memcpy(inbuf, inbuf + in_left, len);
    in_right = len + fread(inbuf + len, 1, SZ - len, stdin);
    in_left = 0;
}

inline void flush() {
    fwrite(outbuf, 1, out_right, stdout);
    out_right = 0;
}

inline void skip_space() {
    if (in_left + 32 > in_right) load();
    while (inbuf[in_left] <= ' ') in_left++;
}

inline void rd(char& c) {
    if (in_left + 32 > in_right) load();
    c = inbuf[in_left++];
}
template <typename T>
inline void rd(T& x) {
    if (in_left + 32 > in_right) load();
    char c;
    do c = inbuf[in_left++];
    while (c < '-');
    [[maybe_unused]] bool minus = false;
    if constexpr (is_signed<T>::value == true) {
        if (c == '-') minus = true, c = inbuf[in_left++];
    }
    x = 0;
    while (c >= '0') {
        x = x * 10 + (c & 15);
        c = inbuf[in_left++];
    }
    if constexpr (is_signed<T>::value == true) {
        if (minus) x = -x;
    }
}
inline void rd() {}
template <typename Head, typename... Tail>
inline void rd(Head& head, Tail&... tail) {
    rd(head);
    rd(tail...);
}

inline void wt(char c) {
    if (out_right > SZ - 32) flush();
    outbuf[out_right++] = c;
}
inline void wt(bool b) {
    if (out_right > SZ - 32) flush();
    outbuf[out_right++] = b ? '1' : '0';
}
template <typename T>
inline void wt(T x) {
    if (out_right > SZ - 32) flush();
    if (!x) {
        outbuf[out_right++] = '0';
        return;
    }
    if constexpr (is_signed<T>::value == true) {
        if (x < 0) outbuf[out_right++] = '-', x = -x;
    }
    int i = 12;
    char buf[16];
    while (x >= 10000) {
        memcpy(buf + i, pre.num + (x % 10000) * 4, 4);
        x /= 10000;
        i -= 4;
    }
    if (x < 100) {
        if (x < 10) {
            outbuf[out_right] = '0' + x;
            ++out_right;
        } else {
            uint32_t q = (uint32_t(x) * 205) >> 11;
            uint32_t r = uint32_t(x) - q * 10;
            outbuf[out_right] = '0' + q;
            outbuf[out_right + 1] = '0' + r;
            out_right += 2;
        }
    } else {
        if (x < 1000) {
            memcpy(outbuf + out_right, pre.num + (x << 2) + 1, 3);
            out_right += 3;
        } else {
            memcpy(outbuf + out_right, pre.num + (x << 2), 4);
            out_right += 4;
        }
    }
    memcpy(outbuf + out_right, buf + i + 4, 12 - i);
    out_right += 12 - i;
}
inline void wt() {}
template <typename Head, typename... Tail>
inline void wt(Head&& head, Tail&&... tail) {
    wt(head);
    wt(forward<Tail>(tail)...);
}
template <typename... Args>
inline void wtn(Args&&... x) {
    wt(forward<Args>(x)...);
    wt('\n');
}

struct Dummy {
    Dummy() { atexit(flush); }
} dummy;

}  // namespace fastio
using fastio::rd;
using fastio::skip_space;
using fastio::wt;
using fastio::wtn;

//

constexpr unsigned int constexpr_primitive_root(unsigned int mod) {
    using u32 = unsigned int;
    using u64 = unsigned long long;
    if (mod == 2) return 1;
    u64 m = mod - 1, ds[32] = {}, idx = 0;
    for (u64 i = 2; i * i <= m; ++i) {
        if (m % i == 0) {
            ds[idx++] = i;
            while (m % i == 0) m /= i;
        }
    }
    if (m != 1) ds[idx++] = m;
    for (u32 _pr = 2, flg = true;; _pr++, flg = true) {
        for (u32 i = 0; i < idx && flg; ++i) {
            u64 a = _pr, b = (mod - 1) / ds[i], r = 1;
            for (; b; a = a * a % mod, b >>= 1)
                if (b & 1) r = r * a % mod;
            if (r == 1) flg = false;
        }
        if (flg == true) return _pr;
    }
}





template <uint32_t mod>
struct LazyMontgomeryModInt {
    using mint = LazyMontgomeryModInt;
    using i32 = int32_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    static constexpr u32 get_r() {
        u32 ret = mod;
        for (i32 i = 0; i < 4; ++i) ret *= 2 - mod * ret;
        return ret;
    }

    static constexpr u32 r = get_r();
    static constexpr u32 n2 = -u64(mod) % mod;
    static_assert(r * mod == 1, "invalid, r * mod != 1");
    static_assert(mod < (1 << 30), "invalid, mod >= 2 ^ 30");
    static_assert((mod & 1) == 1, "invalid, mod % 2 == 0");

    u32 a;

    constexpr LazyMontgomeryModInt() : a(0) {}
    constexpr LazyMontgomeryModInt(const int64_t &b)
        : a(reduce(u64(b % mod + mod) * n2)) {};

    static constexpr u32 reduce(const u64 &b) {
        return (b + u64(u32(b) * u32(-r)) * mod) >> 32;
    }

    constexpr mint &operator+=(const mint &b) {
        if (i32(a += b.a - 2 * mod) < 0) a += 2 * mod;
        return *this;
    }

    constexpr mint &operator-=(const mint &b) {
        if (i32(a -= b.a) < 0) a += 2 * mod;
        return *this;
    }

    constexpr mint &operator*=(const mint &b) {
        a = reduce(u64(a) * b.a);
        return *this;
    }

    constexpr mint &operator/=(const mint &b) {
        *this *= b.inverse();
        return *this;
    }

    constexpr mint operator+(const mint &b) const { return mint(*this) += b; }
    constexpr mint operator-(const mint &b) const { return mint(*this) -= b; }
    constexpr mint operator*(const mint &b) const { return mint(*this) *= b; }
    constexpr mint operator/(const mint &b) const { return mint(*this) /= b; }
    constexpr bool operator==(const mint &b) const {
        return (a >= mod ? a - mod : a) == (b.a >= mod ? b.a - mod : b.a);
    }
    constexpr bool operator!=(const mint &b) const {
        return (a >= mod ? a - mod : a) != (b.a >= mod ? b.a - mod : b.a);
    }
    constexpr mint operator-() const { return mint() - mint(*this); }

    constexpr mint pow(u64 n) const {
        mint ret(1), mul(*this);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    constexpr mint inverse() const { return pow(mod - 2); }

    friend ostream &operator<<(ostream &os, const mint &b) {
        return os << b.get();
    }

    friend istream &operator>>(istream &is, mint &b) {
        int64_t t;
        is >> t;
        b = LazyMontgomeryModInt<mod>(t);
        return (is);
    }

    constexpr u32 get() const {
        u32 ret = reduce(a);
        return ret >= mod ? ret - mod : ret;
    }

    static constexpr u32 get_mod() { return mod; }
};
//

#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")

using m256 = __m256i;
struct alignas(32) mmint {
    m256 x;
    static mmint R, M0, M1, M2, N2;

    mmint() : x() {}
    inline mmint(const m256& _x) : x(_x) {}
    inline mmint(unsigned int a) : x(_mm256_set1_epi32(a)) {}
    inline mmint(unsigned int a0, unsigned int a1, unsigned int a2,
                 unsigned int a3, unsigned int a4, unsigned int a5,
                 unsigned int a6, unsigned int a7)
        : x(_mm256_set_epi32(a7, a6, a5, a4, a3, a2, a1, a0)) {}
    inline operator m256&() { return x; }
    inline operator const m256&() const { return x; }
    inline int& operator[](int i) { return *(reinterpret_cast<int*>(&x) + i); }
    inline const int& operator[](int i) const {
        return *(reinterpret_cast<const int*>(&x) + i);
    }

    friend ostream& operator<<(ostream& os, const mmint& m) {
        unsigned r = R[0], mod = M1[0];
        auto reduce1 = [&](const uint64_t& b) {
            unsigned res = (b + uint64_t(unsigned(b) * unsigned(-r)) * mod) >> 32;
            return res >= mod ? res - mod : res;
        };
        for (int i = 0; i < 8; i++) {
            os << reduce1(m[i]) << (i == 7 ? "" : " ");
        }
        return os;
    }

    template <typename mint>
    static void set_mod() {
        R = _mm256_set1_epi32(mint::r);
        M0 = _mm256_setzero_si256();
        M1 = _mm256_set1_epi32(mint::get_mod());
        M2 = _mm256_set1_epi32(mint::get_mod() * 2);
        N2 = _mm256_set1_epi32(mint::n2);
    }

    static inline mmint reduce(const mmint& prod02, const mmint& prod13) {
        m256 unpalo = _mm256_unpacklo_epi32(prod02, prod13);
        m256 unpahi = _mm256_unpackhi_epi32(prod02, prod13);
        m256 prodlo = _mm256_unpacklo_epi64(unpalo, unpahi);
        m256 prodhi = _mm256_unpackhi_epi64(unpalo, unpahi);
        m256 hiplm1 = _mm256_add_epi32(prodhi, M1);
        m256 prodlohi = _mm256_shuffle_epi32(prodlo, 0xF5);
        m256 lmlr02 = _mm256_mul_epu32(prodlo, R);
        m256 lmlr13 = _mm256_mul_epu32(prodlohi, R);
        m256 prod02_ = _mm256_mul_epu32(lmlr02, M1);
        m256 prod13_ = _mm256_mul_epu32(lmlr13, M1);
        m256 unpalo_ = _mm256_unpacklo_epi32(prod02_, prod13_);
        m256 unpahi_ = _mm256_unpackhi_epi32(prod02_, prod13_);
        m256 prod = _mm256_unpackhi_epi64(unpalo_, unpahi_);
        return _mm256_sub_epi32(hiplm1, prod);
    }

    static inline mmint itom(const mmint& A) { return A * N2; }

    static inline mmint mtoi(const mmint& A) {
        m256 A13 = _mm256_shuffle_epi32(A, 0xF5);
        m256 lmlr02 = _mm256_mul_epu32(A, R);
        m256 lmlr13 = _mm256_mul_epu32(A13, R);
        m256 prod02_ = _mm256_mul_epu32(lmlr02, M1);
        m256 prod13_ = _mm256_mul_epu32(lmlr13, M1);
        m256 unpalo_ = _mm256_unpacklo_epi32(prod02_, prod13_);
        m256 unpahi_ = _mm256_unpackhi_epi32(prod02_, prod13_);
        m256 prod = _mm256_unpackhi_epi64(unpalo_, unpahi_);
        m256 cmp = _mm256_cmpgt_epi32(prod, M0);
        m256 dif = _mm256_and_si256(cmp, M1);
        return _mm256_sub_epi32(dif, prod);
    }

    friend inline mmint operator+(const mmint& A, const mmint& B) {
        m256 apb = _mm256_add_epi32(A, B);
        m256 ret = _mm256_sub_epi32(apb, M2);
        m256 cmp = _mm256_cmpgt_epi32(M0, ret);
        m256 add = _mm256_and_si256(cmp, M2);
        return _mm256_add_epi32(add, ret);
    }

    friend inline mmint operator-(const mmint& A, const mmint& B) {
        m256 ret = _mm256_sub_epi32(A, B);
        m256 cmp = _mm256_cmpgt_epi32(M0, ret);
        m256 add = _mm256_and_si256(cmp, M2);
        return _mm256_add_epi32(add, ret);
    }

    friend inline mmint operator*(const mmint& A, const mmint& B) {
        m256 a13 = _mm256_shuffle_epi32(A, 0xF5);
        m256 b13 = _mm256_shuffle_epi32(B, 0xF5);
        m256 prod02 = _mm256_mul_epu32(A, B);
        m256 prod13 = _mm256_mul_epu32(a13, b13);
        return reduce(prod02, prod13);
    }

    inline mmint& operator+=(const mmint& A) { return (*this) = (*this) + A; }
    inline mmint& operator-=(const mmint& A) { return (*this) = (*this) - A; }
    inline mmint& operator*=(const mmint& A) { return (*this) = (*this) * A; }

    bool operator==(const mmint& A) {
        m256 sub = _mm256_sub_epi32(x, A.x);
        return _mm256_testz_si256(sub, sub) == 1;
    }
    bool operator!=(const mmint& A) { return !((*this) == A); }
};
__attribute__((aligned(32))) mmint mmint::R;
__attribute__((aligned(32))) mmint mmint::M0, mmint::M1, mmint::M2, mmint::N2;

/**
 * @brief vectorize modint
 */

//

using mint = LazyMontgomeryModInt<998244353>;

constexpr uint32_t mod = mint::get_mod();
constexpr uint32_t pr = constexpr_primitive_root(mod);
constexpr int level = __builtin_ctzll(mod - 1);
mint dw[level], dy[level];

void setwy(int k) {
    mint w[level], y[level];
    w[k - 1] = mint(pr).pow((mod - 1) / (1 << k));
    y[k - 1] = w[k - 1].inverse();
    for (int i = k - 2; i > 0; --i)
        w[i] = w[i + 1] * w[i + 1], y[i] = y[i + 1] * y[i + 1];
    dw[1] = w[1], dy[1] = y[1], dw[2] = w[2], dy[2] = y[2];
    for (int i = 3; i < k; ++i) {
        dw[i] = dw[i - 1] * y[i - 2] * w[i];
        dy[i] = dy[i - 1] * w[i - 2] * y[i];
    }
}

void ntt(mmint* a, int k) {
    assert(k > 4);
    if (k & 1) {
        int v = 1 << (k - 1 - 3);
        for (int j = 0; j < v; ++j) {
            mmint ajv = a[j + v];
            a[j + v] = a[j] - ajv;
            a[j] += ajv;
        }
    }
    int u = 1 << (2 + (k & 1));
    int v = 1 << (k - 2 - (k & 1));
    mint one = mint(1), imag = dw[1];
    mmint ONE{one.a}, IMAG{imag.a};

    while (v) {
        if (v == 2) exit(1);
        if (v == 1) {
            /**/
            mint W = one, X = one, Y = imag;
            mmint A = a[0];
            mmint B = {one.a, one.a, W.a,           W.a,
                       one.a, one.a, (W * dw[1]).a, (W * dw[1]).a
                      };
            int jh = 0;
            for (; jh < u - 8; jh += 8) {
                mmint T0 = A * B;
                mmint T0m = _mm256_sub_epi32(mmint::M2, T0);
                mmint T1 = m256(_mm256_shuffle_ps((__m256)T0m.x, (__m256)T0.x, 0x4E));
                mmint C = {one.a, Y.a,           one.a, X.a,
                           one.a, (Y * dw[2]).a, one.a, (X * dw[2]).a
                          };
                X = X * dw[2];
                mmint T2 = (T0 + T1) * C;
                X *= dw[__builtin_ctzll(jh + 8)];
                W = X * X, Y = X * imag;
                mmint T2m = _mm256_sub_epi32(mmint::M2, T2);
                mmint T3 = _mm256_shuffle_epi32(T2, 0x88);
                A = a[jh / 8 + 1];
                B = {one.a, one.a, W.a,           W.a,
                     one.a, one.a, (W * dw[1]).a, (W * dw[1]).a
                    };
                mmint T4 = m256(_mm256_shuffle_ps((__m256)T2.x, (__m256)T2m.x, 0xDD));
                mmint T5 = T3 + T4;
                a[jh / 8] = _mm256_shuffle_epi32(T5, 0x8D);
            }
            mmint T0 = A * B;
            mmint T0m = _mm256_sub_epi32(mmint::M2, T0);
            mmint T1 = m256(_mm256_shuffle_ps((__m256)T0m.x, (__m256)T0.x, 0x4E));
            mmint C = {one.a, Y.a,           one.a, X.a,
                       one.a, (Y * dw[2]).a, one.a, (X * dw[2]).a
                      };
            mmint T2 = (T0 + T1) * C;
            mmint T2m = _mm256_sub_epi32(mmint::M2, T2);
            mmint T3 = _mm256_shuffle_epi32(T2, 0x88);
            mmint T4 = m256(_mm256_shuffle_ps((__m256)T2.x, (__m256)T2m.x, 0xDD));
            mmint T5 = T3 + T4;
            a[jh / 8] = _mm256_shuffle_epi32(T5, 0x8D);
            //*/
        } else if (v == 4) {
            /**/
            mint W = one, X = one, WX = one;
            mmint IMAG1{one.a, one.a, one.a, one.a, imag.a, imag.a, imag.a, imag.a};
            mmint t01 = a[0], t23 = a[1], c01{one.a}, c23{one.a};
            int jh = 0;
            for (; jh < u - 4; jh += 4) {
                X *= dw[__builtin_ctzll(jh + 4)];
                mmint tp = t01 + t23, tm = (t01 - t23) * IMAG1;
                mmint tpm = _mm256_sub_epi32(mmint::M2, tp);
                mmint tmm = _mm256_sub_epi32(mmint::M2, tm);
                t01 = a[(jh >> 1) + 2];
                t23 = a[(jh >> 1) + 3];
                W = X * X, WX = W * X;
                c01 = {one.a, one.a, one.a, one.a, X.a, X.a, X.a, X.a};
                c23 = {W.a, W.a, W.a, W.a, WX.a, WX.a, WX.a, WX.a};
                mmint u0 = _mm256_permute2x128_si256(tp, tpm, 0x00);
                mmint u1 = _mm256_permute2x128_si256(tp, tpm, 0x31);
                mmint u2 = _mm256_permute2x128_si256(tm, tmm, 0x00);
                mmint u3 = _mm256_permute2x128_si256(tm, tmm, 0x31);
                t01 *= c01, t23 *= c23;
                a[(jh >> 1) | 0] = u0 + u1;
                a[(jh >> 1) | 1] = u2 + u3;
            }
            mmint tp = t01 + t23, tm = (t01 - t23) * IMAG1;
            mmint tpm = _mm256_sub_epi32(mmint::M2, tp);
            mmint tmm = _mm256_sub_epi32(mmint::M2, tm);
            mmint u0 = _mm256_permute2x128_si256(tp, tpm, 0x00);
            mmint u1 = _mm256_permute2x128_si256(tp, tpm, 0x31);
            mmint u2 = _mm256_permute2x128_si256(tm, tmm, 0x00);
            mmint u3 = _mm256_permute2x128_si256(tm, tmm, 0x31);
            a[(jh >> 1) | 0] = u0 + u1;
            a[(jh >> 1) | 1] = u2 + u3;
            //*/
        } else {
            int v8 = v / 8;
            {
                int j0 = 0, j1 = v8, j2 = j1 + v8, j3 = j2 + v8;
                int je = v8 - 1;
                mmint t0 = a[j0], t1 = a[j1];
                mmint t2 = a[j2], t3 = a[j3];
                for (; j0 < je; ++j0, ++j1, ++j2, ++j3) {
                    mmint t0p2 = t0 + t2, t1p3 = t1 + t3;
                    mmint t0m2 = t0 - t2, t1m3 = (t1 - t3) * IMAG;
                    t0 = a[j0 + 1], t1 = a[j1 + 1];
                    t2 = a[j2 + 1], t3 = a[j3 + 1];
                    a[j0] = t0p2 + t1p3, a[j1] = t0p2 - t1p3;
                    a[j2] = t0m2 + t1m3, a[j3] = t0m2 - t1m3;
                }
                mmint t0p2 = t0 + t2, t1p3 = t1 + t3;
                mmint t0m2 = t0 - t2, t1m3 = (t1 - t3) * IMAG;
                a[j0] = t0p2 + t1p3, a[j1] = t0p2 - t1p3;
                a[j2] = t0m2 + t1m3, a[j3] = t0m2 - t1m3;
            }
            mmint W{one.a}, X{dw[2].a}, Y{one.a};
            for (int jh = 4; jh < u;) {
                W = X * X, Y = X * IMAG;
                int j0 = jh * v8, j1 = j0 + v8, j2 = j1 + v8, j3 = j2 + v8;
                int je = j0 + v8 - 1;
                mmint t0, t1, t2, t3;
                t0 = a[j0], t1 = a[j1];
                t2 = a[j2] * W, t3 = a[j3] * W;
                for (; j0 < je; ++j0, ++j1, ++j2, ++j3) {
                    mmint t0p2 = t0 + t2, t1p3 = (t1 + t3) * X;
                    mmint t0m2 = t0 - t2, t1m3 = (t1 - t3) * Y;
                    t0 = a[j0 + 1], t1 = a[j1 + 1];
                    t2 = a[j2 + 1] * W, t3 = a[j3 + 1] * W;
                    a[j0] = t0p2 + t1p3, a[j1] = t0p2 - t1p3;
                    a[j2] = t0m2 + t1m3, a[j3] = t0m2 - t1m3;
                }
                mmint t0p2 = t0 + t2, t1p3 = (t1 + t3) * X;
                mmint t0m2 = t0 - t2, t1m3 = (t1 - t3) * Y;
                a[j0] = t0p2 + t1p3, a[j1] = t0p2 - t1p3;
                a[j2] = t0m2 + t1m3, a[j3] = t0m2 - t1m3;
                X *= mmint{dw[__builtin_ctzll((jh += 4))].a};
            }
        }
        u <<= 2;
        v >>= 2;
    }
}

void intt(mmint* a, int k) {
    assert(k > 4);
    int u = 1 << (k - 2);
    int v = 1;
    mint one = mint(1), imag = dy[1];
    mmint ONE{one.a}, IMAG{imag.a};
    while (u) {
        assert(v != 2);
        if (v == 1) {
            u <<= 2;
            /**/
            mint W = one, X = one, Y = one;
            for (int jh = 0; jh < u;) {
                W = X * X, Y = X * imag;
                mmint& A = a[jh / 8];
                mint t0, t1, t2, t3, t4, t5, t6, t7;
                t0.a = A[0], t1.a = A[1], t2.a = A[2], t3.a = A[3];
                t4.a = A[4], t5.a = A[5], t6.a = A[6], t7.a = A[7];
                mint t0p1 = t0 + t1, t2p3 = t2 + t3;
                mint t0m1 = (t0 - t1) * X, t2m3 = (t2 - t3) * Y;
                A[0] = (t0p1 + t2p3).a;
                A[1] = (t0m1 + t2m3).a;
                A[2] = ((t0p1 - t2p3) * W).a;
                A[3] = ((t0m1 - t2m3) * W).a;
                X *= dy[2], W *= dy[1], Y *= dy[2];
                mint t4p5 = t4 + t5, t6p7 = t6 + t7;
                mint t4m5 = (t4 - t5) * X, t6m7 = (t6 - t7) * Y;
                A[4] = (t4p5 + t6p7).a;
                A[5] = (t4m5 + t6m7).a;
                A[6] = ((t4p5 - t6p7) * W).a;
                A[7] = ((t4m5 - t6m7) * W).a;
                X *= dy[__builtin_ctzll(jh += 8)];
            }
            //*/
        } else if (v == 4) {
            u <<= 2;
            /**/
            mint W = one, X = one, Y = one;

            for (int jh = 0; jh < u;) {
                W = X * X, Y = X * imag;
                mmint c23{X.a, X.a, X.a, X.a, Y.a, Y.a, Y.a, Y.a};
                mmint Ws{W.a};
                mmint t01 = a[(jh >> 1) | 0];
                mmint t23 = a[(jh >> 1) | 1];
                mmint t02 = _mm256_permute2x128_si256(t01, t23, 0x20);
                mmint t13 = _mm256_permute2x128_si256(t01, t23, 0x31);
                mmint tp = t02 + t13, tm = (t02 - t13) * c23;
                mmint u0 = _mm256_permute2x128_si256(tp, tm, 0x20);
                mmint u1 = _mm256_permute2x128_si256(tp, tm, 0x31);
                a[(jh >> 1) | 0] = u0 + u1;
                a[(jh >> 1) | 1] = (u0 - u1) * Ws;
                X *= dy[__builtin_ctzll((jh += 4))];
            }
            //*/
        } else {
            int v8 = v / 8;
            u <<= 2;
            {
                int j0 = 0, j1 = v8, j2 = j1 + v8, j3 = j2 + v8;
                int je = v8 - 1;
                mmint t0 = a[j0], t1 = a[j1];
                mmint t2 = a[j2], t3 = a[j3];
                for (; j0 < je; ++j0, ++j1, ++j2, ++j3) {
                    mmint t0p1 = t0 + t1, t2p3 = t2 + t3;
                    mmint t0m1 = t0 - t1, t2m3 = (t2 - t3) * IMAG;
                    t0 = a[j0 + 1], t1 = a[j1 + 1];
                    t2 = a[j2 + 1], t3 = a[j3 + 1];
                    a[j0] = t0p1 + t2p3, a[j1] = t0m1 + t2m3;
                    a[j2] = t0p1 - t2p3, a[j3] = t0m1 - t2m3;
                }
                mmint t0p1 = t0 + t1, t2p3 = t2 + t3;
                mmint t0m1 = t0 - t1, t2m3 = (t2 - t3) * IMAG;
                a[j0] = t0p1 + t2p3, a[j1] = t0m1 + t2m3;
                a[j2] = t0p1 - t2p3, a[j3] = t0m1 - t2m3;
            }
            mmint W{one.a}, X{dy[2].a}, Y{one.a};
            for (int jh = 4; jh < u;) {
                W = X * X, Y = X * IMAG;
                int j0 = jh * v8, j1 = j0 + v8, j2 = j1 + v8, j3 = j2 + v8;
                int je = j0 + v8 - 1;
                mmint t0 = a[j0], t1 = a[j1];
                mmint t2 = a[j2], t3 = a[j3];
                for (; j0 < je; ++j0, ++j1, ++j2, ++j3) {
                    mmint t0p1 = t0 + t1, t2p3 = t2 + t3;
                    mmint t0m1 = (t0 - t1) * X, t2m3 = (t2 - t3) * Y;
                    t0 = a[j0 + 1], t1 = a[j1 + 1];
                    t2 = a[j2 + 1], t3 = a[j3 + 1];
                    a[j0] = t0p1 + t2p3, a[j1] = t0m1 + t2m3;
                    a[j2] = (t0p1 - t2p3) * W, a[j3] = (t0m1 - t2m3) * W;
                }
                mmint t0p1 = t0 + t1, t2p3 = t2 + t3;
                mmint t0m1 = (t0 - t1) * X, t2m3 = (t2 - t3) * Y;
                a[j0] = t0p1 + t2p3, a[j1] = t0m1 + t2m3;
                a[j2] = (t0p1 - t2p3) * W, a[j3] = (t0m1 - t2m3) * W;
                X *= mmint{dy[__builtin_ctzll(jh += 4)].a};
            }
        }
        u >>= 4;
        v <<= 2;
    }
    if (k & 1) {
        u = 1 << (k - 1 - 3);
        for (int j = 0; j < u; ++j) {
            mmint ajv = a[j] - a[j + u];
            a[j] += a[j + u];
            a[j + u] = ajv;
        }
    }
}

vector<mint> multiply(const vector<mint>& a, const vector<mint>& b) {
    int l = a.size() + b.size() - 1;
    if (min<int>(a.size(), b.size()) < 20) {
        vector<mint> s(l);
        for (int i = 0; i < (int)a.size(); ++i)
            for (int j = 0; j < (int)b.size(); ++j) s[i + j] += a[i] * b[j];
        return s;
    }
    int k = 2, M = 4;
    while (M < l) M <<= 1, ++k;
    setwy(k);
    mmint *s = new mmint[M / 8], *t = new mmint[M / 8];
    memset((__attribute__((aligned(32))) int*)s, 0, sizeof(int) * M);
    memset((__attribute__((aligned(32))) int*)t, 0, sizeof(int) * M);
    for (int i = 0; i < (int)a.size(); ++i) s[i / 8][i % 8] = a[i].a;
    for (int i = 0; i < (int)b.size(); ++i) t[i / 8][i % 8] = b[i].a;
    ntt(s, k), ntt(t, k);
    for (int i = 0; i < M / 8; ++i) s[i] *= t[i];
    intt(s, k);
    mmint invm{(mint(M).inverse()).a};
    for (int i = 0; i < M / 8; i++) s[i] *= invm;
    vector<mint> res(M);
    std::memcpy(res.data(), (mint*)s, sizeof(int) * M);
    return res;
}

mmint a[1 << 17], b[1 << 17];

void inline_multiply(int M, int k) {
    setwy(k);
    ntt(a, k), ntt(b, k);
    for (int i = 0; i < M / 8; i += 4) {
        a[i + 0] *= b[i + 0];
        a[i + 1] *= b[i + 1];
        a[i + 2] *= b[i + 2];
        a[i + 3] *= b[i + 3];
    }
    intt(a, k);
}

void LC_convolution() {
    int n, m;
    unsigned int x, *p;

    rd(n, m);
    // exampleでサーバーを温める
    if (n == 4 or n == 1) {
        for (int N : vector<int> {270000, 70000, 17000}) {
            vector<mint> s(N), t(N);
            auto u = multiply(s, t);
        }
    }

    p = (unsigned int*)a;
    for (int i = 0; i < n; ++i, ++p) rd(x), *p = x;
    p = (unsigned int*)b;
    for (int i = 0; i < m; ++i, ++p) rd(x), *p = x;

    if (n + m < 32) {
        long long c[32] = {};
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                c[i + j] += 1LL * a[i / 8][i % 8] * b[j / 8][j % 8];
                c[i + j] %= 998244353;
            }
        }
        for (int i = 0; i < n + m - 1; i++) wt(c[i], " \n"[i == n + m - 2]);
        return;
    }

    int k = 2, M = 4, l = n + m - 1;
    while (M < l) M <<= 1, ++k;
    for (int i = 0; i < (n + 7) / 8; ++i) {
        a[i + 0] = mmint::itom(a[i + 0]);
    }
    for (int i = 0; i < (m + 7) / 8; ++i) {
        b[i + 0] = mmint::itom(b[i + 0]);
    }
    inline_multiply(M, k);
    mmint im{(mint(M).inverse()).a};
    for (int i = 0; i < M / 8; ++i) {
        a[i + 0] *= im;
        a[i + 0] = mmint::mtoi(a[i + 0]);
    }

    p = (unsigned int*)a;
    wt(*p);
    p++;
    for (int i = 1; i < l; i++, p++) {
        wt(' ');
        wt(*p);
    }
    wt('\n');
    return;
}

#include <unordered_set>



namespace my_rand {

// [0, 2^64 - 1)
uint64_t rng() {
    static uint64_t x_ =
        uint64_t(chrono::duration_cast<chrono::nanoseconds>(
                     chrono::high_resolution_clock::now().time_since_epoch())
                 .count()) *
        10150724397891781847ULL;
    x_ ^= x_ << 7;
    return x_ ^= x_ >> 9;
}

// [l, r)
int64_t randint(int64_t l, int64_t r) {
    assert(l < r);
    return l + rng() % (r - l);
}

// choose n numbers from [l, r) without overlapping
vector<int64_t> randset(int64_t l, int64_t r, int64_t n) {
    assert(l <= r && n <= r - l);
    unordered_set<int64_t> s;
    for (int64_t i = n; i; --i) {
        int64_t m = randint(l, r + 1 - i);
        if (s.find(m) != s.end()) m = r - i;
        s.insert(m);
    }
    vector<int64_t> ret;
    for (auto& x : s) ret.push_back(x);
    return ret;
}

// [0.0, 1.0)
double rnd() {
    union raw_cast {
        double t;
        uint64_t u;
    };
    constexpr uint64_t p = uint64_t(1023 - 64) << 52;
    return rng() * ((raw_cast*)(&p))->t;
}

template <typename T>
void randshf(vector<T>& v) {
    int n = v.size();
    for (int loop = 0; loop < 2; loop++)
        for (int i = 0; i < n; i++) swap(v[i], v[randint(0, n)]);
}

}  // namespace my_rand

using my_rand::randint;
using my_rand::randset;
using my_rand::randshf;
using my_rand::rnd;
using my_rand::rng;


struct Timer {
    chrono::high_resolution_clock::time_point st;

    Timer() { reset(); }

    void reset() { st = chrono::high_resolution_clock::now(); }

    chrono::milliseconds::rep elapsed() {
        auto ed = chrono::high_resolution_clock::now();
        return chrono::duration_cast<chrono::milliseconds>(ed - st).count();
    }
};
void calc_time() {
    /**/
    {
        vector<mint> s(21), t(21);
        for (int i = 0; i < 20; i++) s[i] = t[i] = 1;
        auto u = multiply(s, t);
        for (auto& x : u) std::cerr << x << " ";
        std::cerr << endl;
        for (int i = 0; i <= 39; i++) assert(u[i] == min(1 + i, 39 - i));
        assert(u[39] == u[40] and u[40] == 0);
    }

    for (int N : vector<int> {64, 128, 256, 512, 1024, 114514}) {
        vector<mint> s(N), t(N);
        auto u = multiply(s, t);
        for (auto& x : u) assert(x == 0);
    }

    std::cerr << "verify OK" << endl;
    //*/

    int k = 20;
    for (int i = 0; i < 1 << k; i++) {
        a[i / 8][i & 8] = unsigned(i) * 123456789u;
    }

    {
        Timer timer;
        mmint sm = 0;
        for (int i = 0; i < 1000; i++) {
            intt(a, k);
            sm += a[rng() & 131071];
        }
        std::cerr << sm << endl;
        std::cerr << timer.elapsed() << endl;
    }

    {
        Timer timer;
        mmint sm = 0;
        for (int i = 0; i < 1000; i++) {
            ntt(a, k);
            sm += a[rng() & 131071];
        }
        std::cerr << sm << endl;
        std::cerr << timer.elapsed() << endl;
    }
}

int main() {
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif
    auto t1 = high_resolution_clock::now();
    mmint::set_mod<mint>();

#if defined(NyaanLocal) || defined(PROFILER)
    calc_time();
    return 0;
#endif
    LC_convolution();
    auto t2 = high_resolution_clock::now();
    auto time = duration_cast<duration<double>>(t2 - t1);
    cerr << "Time elapsed = " << time.count() << endl;
}