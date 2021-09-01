const int N = (int)1e5 + 5;
ll factorialNumInverse[N], naturalNumInverse[N], fact[N];

void InverseofNumber(ll p) {
    naturalNumInverse[0] = naturalNumInverse[1] = 1;
    for (int i = 2; i < N; ++i) {
        naturalNumInverse[i] = naturalNumInverse[p % i] * (p - p / i) % p;
    }
}
void InverseofFactorial(ll p) {
    factorialNumInverse[0] = factorialNumInverse[1] = 1;
    for (int i = 2; i < N; ++i) {
        factorialNumInverse[i] = (naturalNumInverse[i] * factorialNumInverse[i - 1]) % p;
    }
}
void factorial(ll p)
{
    fact[0] = 1;
    for (int i = 1; i < N; i++) {
        fact[i] = (fact[i - 1] * i) % p;
    }
}

ll ncr(ll N, ll R, ll p)
{
    ll ans = ((fact[N] * factorialNumInverse[R]) % p * factorialNumInverse[N - R]) % p;
    return ans;
}
ll p = 1000000007;
InverseofNumber(p);
InverseofFactorial(p);
factorial(p);