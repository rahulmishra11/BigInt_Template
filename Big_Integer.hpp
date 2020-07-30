/*Convert int main to int_32 main
 as int is declared here as long long int.

 Use the header  #include "Big_Integer.h"
 */

#define int long long int
const int base = 1e10;
const int base_digits = 9;
#define vector<int> vll
class BigInt
{
public:
    // A vector to carry 9 digits
    vector<int> z;
    // sign of the number
    // 1 -> Positive
    // -1 ->Negative
    int sign;
    BigInt()
    {
        sign = 1;
    }
    BigInt(int a) {
        *this = a;
    }

    BigInt(string &str) {
        read_str(str);
    }

    void operator=(const BigInt &a) {
        sign = a.sign;
        z = a.z;
    }

    void operator=(int x) {
        sign = 1;
        if (x < 0) {
            sign = -1, x = -x;
        }
        z.clear();
        while ( x > 0) {
            z.push_back(x % base);
            x = x / base;
        }
    }

    BigInt operator+(const BigInt &x) const {
        if (sign == x.sign) {
            BigInt res = x;
            for (int i = 0, carry = 0; i < (int) max(z.size(), x.z.size()) || carry; ++i) {
                if (i == (int) res.z.size()) {
                    res.z.push_back(0);
                }
                res.z[i] += carry + (i < (int) z.size() ? z[i] : 0);
                carry = res.z[i] >= base;
                if (carry) {
                    res.z[i] -= base;
                }
            }
            return res;
        }
        else {
            return *this - (-x);
        }
    }

    void operator*=(int x) {
        if (x < 0) {
            sign = -sign, x = -x;
        }
        for (int i = 0, carry = 0; i < (int) z.size() || carry; ++i) {
            if (i == (int) z.size()) {
                z.push_back(0);
            }
            int cur = z[i] * x + carry;
            carry = cur / base;
            z[i] = cur % base;
        }
        trim();
    }

    BigInt operator*(int x) const {
        BigInt res = *this;
        res *= x;
        return res;
    }

    BigInt operator-(const BigInt &x) const {
        if (sign == x.sign) {
            if (abs() >= x.abs()) {
                BigInt res = *this;
                for (int i = 0, carry = 0; i < (int) x.z.size() || carry; ++i) {
                    res.z[i] -= carry + (i < (int) x.z.size() ? x.z[i] : 0);
                    carry = res.z[i] < 0;
                    if (carry) {
                        res.z[i] += base;
                    }
                }
                res.trim();
                return res;
            }
            else {
                return -(x - *this);
            }
        }
        else {
            return *this + (-x);
        }
    }


    void trim() {
        while (!z.empty() && z.back() == 0) {
            z.pop_back();
        }
        if (z.empty()) {
            sign = 1;
        }
    }


    void read_str(const string &s) {
        sign = 1;
        z.clear();
        int pos = 0;
        while (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) {
            if (s[pos] == '-') {
                sign = -sign;
            }
            ++pos;
        }
        for (int i = s.size() - 1; i >= pos; i -= base_digits) {
            int x = 0;
            for (int j = max(pos, i - base_digits + 1); j <= i; j++) {
                x = x * 10 + s[j] - '0';
            }
            z.push_back(x);
        }
        trim();
    }



    BigInt operator/(const BigInt &x) const {
        return divmod(*this, x).first;
    }

    BigInt operator%(const BigInt &x) const {
        return divmod(*this, x).second;
    }

    void operator/=(int x) {
        if (x < 0) {
            sign = -sign, x = -x;
        }
        for (int i = z.size() - 1, rem = 0; i >= 0; --i) {
            int cur = z[i] + rem * base;
            z[i] = cur / x;
            rem = cur % x;
        }
        trim();
    }

    BigInt operator/(int x) const {
        BigInt res = *this;
        res /= x;
        return res;
    }

    int operator%(int x) const {
        if (x < 0) {
            x = -x;
        }
        int m = 0;
        for (int i = z.size() - 1; i >= 0; --i) {
            m = (m * base + z[i]) % x;
        }
        return m * sign;
    }

    void operator+=(const BigInt &x) {
        *this = *this + x;
    }
    void operator-=(const BigInt &x) {
        *this = *this - x;
    }
    void operator*=(const BigInt &x) {
        *this = *this * x;
    }
    void operator/=(const BigInt &x) {
        *this = *this / x;
    }
    void operator%=(const BigInt &x) {
        *this = *this % x;
    }
    void operator%=(int x) {
        *this = *this % x;
    }

    bool operator<(const BigInt &x) const {
        if (sign != x.sign) {
            return sign < x.sign;
        }
        if (z.size() != x.z.size()) {
            return z.size() * sign < x.z.size() * x.sign;
        }
        for (int i = z.size() - 1; i >= 0; i--) {
            if (z[i] != x.z[i]) {
                return z[i] * sign < x.z[i] * sign;
            }
        }
        return false;
    }

    bool operator>(const BigInt &x) const {
        return x < *this;
    }
    bool operator<=(const BigInt &x) const {
        return !(x < *this);
    }
    bool operator>=(const BigInt &x) const {
        return !(*this < x);
    }
    bool operator==(const BigInt &x) const {
        return !(*this < x) && !(x < *this);
    }
    bool operator!=(const BigInt &x) const {
        return *this < x || x < *this;
    }

    bool isZero() const {
        return z.empty() || ((int) z.size() == 1 && !z[0]);
    }

    BigInt operator-() const {
        BigInt res = *this;
        res.sign = -sign;
        return res;
    }

    BigInt abs() const {
        BigInt res = *this;
        res.sign *= res.sign;
        return res;
    }

    int longValue() const {
        int res = 0;
        for (int i = z.size() - 1; i >= 0; i--) {
            res = res * base + z[i];
        }
        return res * sign;
    }

    friend BigInt gcd(const BigInt &a, const BigInt &b) {
        return b.isZero() ? a : gcd(b, a % b);
    }
    friend BigInt lcm(const BigInt &a, const BigInt &b) {
        return a / gcd(a, b) * b;
    }

    friend istream &operator>>(istream &stream, BigInt &x) {
        string s;
        stream >> s;
        x.read_str(s);
        return stream;
    }

    friend ostream &operator<<(ostream &stream, const BigInt &x) {
        if (x.sign == -1) {
            stream << '-';
        }
        stream << (x.z.empty() ? 0 : x.z.back());
        for (int i = x.z.size() - 2; i >= 0; --i) {
            stream << setw(base_digits) << setfill('0') << x.z[i];
        }
        return stream;
    }

    static vector<int> convert_base(const vector<int> &a, int old_digits, int new_digits) {
        vector<int> p(max(old_digits, new_digits) + 1);
        p[0] = 1;
        for (int i = 1; i < (int) p.size(); i++) {
            p[i] = p[i - 1] * 10;
        }
        vector<int> res;
        int cur = 0;
        int cur_digits = 0;
        for (int i = 0; i < (int) a.size(); i++) {
            cur += a[i] * p[cur_digits];
            cur_digits += old_digits;
            while (cur_digits >= new_digits) {
                res.push_back(cur % p[new_digits]);
                cur /= p[new_digits];
                cur_digits -= new_digits;
            }
        }
        res.push_back(cur);
        while (!res.empty() && res.back() == 0) {
            res.pop_back();
        }
        return res;
    }



    static vll karatsubaMultiply(const vll &a, const vll &b) {
        int n = a.size();
        vll res(n + n);
        if (n <= 32) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }
        int k = n >> 1;
        vll a1(a.begin(), a.begin() + k);
        vll a2(a.begin() + k, a.end());
        vll b1(b.begin(), b.begin() + k);
        vll b2(b.begin() + k, b.end());
        vll a1b1 = karatsubaMultiply(a1, b1);
        vll a2b2 = karatsubaMultiply(a2, b2);
        for (int i = 0; i < k; i++) {
            a2[i] += a1[i];
        }
        for (int i = 0; i < k; i++) {
            b2[i] += b1[i];
        }
        vll r = karatsubaMultiply(a2, b2);
        for (int i = 0; i < (int) a1b1.size(); i++) {
            r[i] -= a1b1[i];
        }
        for (int i = 0; i < (int) a2b2.size(); i++) {
            r[i] -= a2b2[i];
        }
        for (int i = 0; i < (int) r.size(); i++) {
            res[i + k] += r[i];
        }
        for (int i = 0; i < (int) a1b1.size(); i++) {
            res[i] += a1b1[i];
        }
        for (int i = 0; i < (int) a2b2.size(); i++) {
            res[i + n] += a2b2[i];
        }
        return res;
    }

    BigInt operator*(const BigInt &x) const {
        vector<int> a6 = convert_base(this->z, base_digits, 6);
        vector<int> b6 = convert_base(x.z, base_digits, 6);
        vll a(a6.begin(), a6.end());
        vll b(b6.begin(), b6.end());
        while (a.size() < b.size()) {
            a.push_back(0);
        }
        while (b.size() < a.size()) {
            b.push_back(0);
        }
        while (a.size() & (a.size() - 1)) {
            a.push_back(0);
            b.push_back(0);
        }
        vll c = karatsubaMultiply(a, b);
        BigInt res;
        res.sign = sign * x.sign;
        for (int i = 0, carry = 0; i < (int) c.size(); i++) {
            int cur = c[i] + carry;
            res.z.push_back(cur % 1000000);
            carry = cur / 1000000;
        }
        res.z = convert_base(res.z, 6, base_digits);
        res.trim();
        return res;
    }


    friend pair<BigInt, BigInt> divmod(const BigInt &a1, const BigInt &b1) {
        int norm = base / (b1.z.back() + 1);
        BigInt a = a1.abs() * norm;
        BigInt b = b1.abs() * norm;
        BigInt q, r;
        q.z.resize(a.z.size());
        for (int i = a.z.size() - 1; i >= 0; i--) {
            r *= base;
            r += a.z[i];
            int s1 = b.z.size() < r.z.size() ? r.z[b.z.size()] : 0;
            int s2 = b.z.size() - 1 < r.z.size() ? r.z[b.z.size() - 1] : 0;
            int d = (s1 * base + s2) / b.z.back();
            r -= b * d;
            while (r < 0) {
                r += b, --d;
            }
            q.z[i] = d;
        }
        q.sign = a1.sign * b1.sign;
        r.sign = a1.sign;
        q.trim();
        r.trim();
        return make_pair(q, r / norm);
    }

    friend BigInt sqrt(const BigInt &a1) {
        BigInt a = a1;
        while (a.z.empty() || (int) a.z.size() % 2 == 1) {
            a.z.push_back(0);
        }
        int n = a.z.size();
        int firstDigit = sqrt(a.z[n - 1] * base + a.z[n - 2]);
        int norm = base / (firstDigit + 1);
        a *= norm;
        a *= norm;
        while (a.z.empty() || (int) a.z.size() % 2 == 1) {
            a.z.push_back(0);
        }
        BigInt r = a.z[n - 1] * base + a.z[n - 2];
        firstDigit = sqrt(a.z[n - 1] * base + a.z[n - 2]);
        int q = firstDigit;
        BigInt res;
        for (int j = n / 2 - 1; j >= 0; j--) {
            for (;; --q) {
                BigInt r1 = (r - (res * 2 * base + q) * q) * base * base +
                            (j > 0 ? a.z[2 * j - 1] * base + a.z[2 * j - 2] : 0);
                if (r1 >= 0) {
                    r = r1;
                    break;
                }
            }
            res *= base;
            res += q;
            if (j > 0) {
                int d1 = res.z.size() + 2 < r.z.size() ? r.z[res.z.size() + 2] : 0;
                int d2 = res.z.size() + 1 < r.z.size() ? r.z[res.z.size() + 1] : 0;
                int d3 = res.z.size() < r.z.size() ? r.z[res.z.size()] : 0;
                q = (d1 * base * base + d2 * base + d3) / (firstDigit * 2);
            }
        }
        res.trim();
        return res / norm;
    }


};