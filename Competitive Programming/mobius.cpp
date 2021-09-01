const int N = 1e6 + 5;
bool composite[N];
int prime[N], mf[N], sz;

void sieve() {
    mf[1] = 1;
    for (int i = 2; i < N; i++) {
        if (!composite[i])prime[sz++] = i, mf[i] = - 1;
        for (int j = 0; j < sz && i * prime[j] < N; j++) {
            composite[i * prime[j]] = true;
            if (i % prime[j] == 0) {
                mf[i * prime[j]] = 0;
                break;
            }
            else {
                mf[i * prime[j]] = -mf[i];
            }
        }
    }
}