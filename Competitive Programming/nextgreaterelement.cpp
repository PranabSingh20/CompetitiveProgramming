vector<pair<int, int>> next_greater_element(vector<int> v) {
    int n = (int)v.size();
    stack<pair<int, int>> s;
    vector<pair<int, int>> ret(n);
    for (int i = n - 1; i >= 0; i--) {
        while (!s.empty() && s.top().first <= v[i]) {
            s.pop();
        }
        if (s.empty()) {
            ret[i] = { -1, - 1};
        }
        else {
            ret[i] = s.top();
        }
        s.push({v[i], i});
    }
    return ret;
}