const int N = 1e6 + 5;
bool composite[N];
int prime[N], lp[N], sz;

void sieve() {
    for (int i = 2; i < N; i++) {
        if (!composite[i])prime[sz++] = i, lp[i] = i;
        for (int j = 0; j < sz && i * prime[j] < N; j++) {
            composite[i * prime[j]] = true;
            lp[i * prime[j]] = prime[j];
            if (i % prime[j] == 0) {
                break;
            }
        }
    }
}

vector<int> getFactorization(int x) {
    vector<int> ret;
    while (x != 1)
    {
        ret.push_back(lp[x]);
        x = x / lp[x];
    }
    return ret;
}