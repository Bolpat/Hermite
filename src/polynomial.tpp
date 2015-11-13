// Compile with clang++-3.0 -std=c++11

// Polynomial<T> //
namespace Modulus
{

// Constructors //

template <typename T, typename deg_type> inline
Polynomial<T, deg_type>::Polynomial(T const & t, deg_type deg)
{
    if (t != T()) coeffs[deg] = t;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type>
Polynomial<T, deg_type>::fromCoeffVector(std::vector<T> const & coeffs)
{
    Polynomial res;
    for (size_t i = 0; i < coeffs.size(); ++i)
        if (coeffs[i] != T())
            res.coeffs[i] = coeffs[i];
    return res;
}


// Methods //

template <typename T, typename deg_type> inline
deg_type deg(Polynomial<T, deg_type> const & p)
{
    if (p.coeffs.empty()) return 0;
    return p.coeffs.rbegin()->first;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type>
Polynomial<T, deg_type>::with_monic(deg_type dg) const
{
    Polynomial result = *this;
    const deg_type d = std::max(deg(result) + 1, dg);
    result.coeffs[d] = T(1);
    return result;
}


// Shift //

template <typename T, typename deg_type>
Polynomial<T, deg_type>
    operator <<(Polynomial<T, deg_type> const & p, deg_type n)
{
    Polynomial<T, deg_type> res;
    for (auto & dcp : p.coeffs) res.coeffs[dcp.first + n] = dcp.second;
    return res;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type>
    operator >>(Polynomial<T, deg_type> const & p, deg_type n)
{
    Polynomial<T, deg_type> res;
    for (auto & dcp : p.coeffs) if (dcp.first >= n) res.coeffs[dcp.first - n] = dcp.second;
    return res;
}


// Arithmetric //

template <typename T, typename deg_type>
Polynomial<T, deg_type>
Polynomial<T, deg_type>::operator - () const
{
    Polynomial<T, deg_type> res = *this;
    for (auto & dcp : res.coeffs) dcp.second = -dcp.second;
    return res;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type>
    operator + (Polynomial<T, deg_type> const & p, Polynomial<T, deg_type> const & q)
{
    if (deg(p) > deg(q)) return q + p;
    Polynomial<T, deg_type> res = p;
    res += q;
    return res;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type>
    operator * (Polynomial<T, deg_type> const & p, Polynomial<T, deg_type> const & q)
{
    Polynomial<T, deg_type> res;
    for (auto & dcp : q.coeffs) res += dcp.second * p << dcp.first;
    return res;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type> &
Polynomial<T, deg_type>::operator *=(T const & t) &
{
    for (auto & coeff : coeffs) coeff.second *= t;
    return *this;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type> &
Polynomial<T, deg_type>::operator +=(Polynomial<T, deg_type> const & q) &
{
    for (auto & dcp : q.coeffs) if ((coeffs[dcp.first] += dcp.second) == T()) coeffs.erase(dcp.first);
    return *this;
}


template <typename T, typename deg_type>
std::pair<Polynomial<T, deg_type>, Polynomial<T, deg_type>>
Polynomial<T, deg_type>::divmod(Polynomial<T, deg_type> a, Polynomial<T, deg_type> const & b)
{
    if (b.coeffs.empty()) return std::make_pair(Polynomial(), Polynomial());
    
    using namespace std;
    Polynomial<T, deg_type> q;
    while (deg(a) >= deg(b))
    {
        auto const d = deg(a) - deg(b);
        auto const c = a.leading_coeff() / b.leading_coeff();
        q.coeffs[d] = c;
        auto const subtr = b * c << d;
        a -= subtr;
    }
    return std::make_pair(std::move(q), std::move(a));
}

template <typename T, typename deg_type> inline
void Polynomial<T, deg_type>::divmod(Polynomial<T, deg_type> const & a,
                                     Polynomial<T, deg_type> const & b,
                                     Polynomial<T, deg_type>       & q,
                                     Polynomial<T, deg_type>       & r)
{
    auto qr_pair = divmod(a, b);
    q = std::move(qr_pair.first);
    r = std::move(qr_pair.second);
}

// Access //

template <typename T, typename deg_type>
T const & Polynomial<T, deg_type>::at(deg_type k) const &
{
    auto const it = coeffs.find(k);
    if (it == coeffs.end()) return T();
    return it->second;
}

template <typename T, typename deg_type>
T         Polynomial<T, deg_type>::at(deg_type k)      &&
{
    auto it = coeffs.find(k);
    if (it == coeffs.end()) return T();
    return move(it->second);
}

template <typename T, typename deg_type>
T squareMultiply(T const & v, deg_type deg)
{   
    if (deg == 0) return T(1);
    deg_type zero = 0;
    deg_type hiBit = ~((~zero) >> 1U);
    while ((deg & hiBit) == zero) hiBit >>= 1U;
    // hiBit masks first bit of deg.
    
    hiBit >>= 1U;
    T prod = v; // consume the bit.
    
    while (hiBit != 0)
    {
        prod *= prod;
        if (deg & hiBit) prod *= v;
        hiBit >>= 1U;
    }
    return prod;
}

template <typename T, typename deg_type>
T power(T const & value, unsigned k, std::map<deg_type, T> & power_memo)
{
    auto pow_it = power_memo.find(k);
    if (pow_it == power_memo.end())
    {
        auto p = squareMultiply(value, k);
        power_memo[k] = p;
        return p;
    }
    else return pow_it->second;
}

template <typename T, typename deg_type>
T       Polynomial<T, deg_type>::operator ()(T const & value) const
{
    if (coeffs.empty()) return T();
    // under here coeffs is nonempty.

    if (coeffs.size() == 1)
        return coeffs.begin()->second * squareMultiply(value, coeffs.begin()->first);
    // under here coeffs has at least 2 elements.
    std::map<deg_type, T> power_memo;
    power_memo[0] = T(1);
    power_memo[1] = value;

    std::vector<deg_type> power_diffs;
    for (auto it2 =  coeffs.begin(), it = it2++;
              it2 != coeffs.end();
            ++it2,                ++it)
    {
      power_diffs.push_back(it2->first - it->first);
    }

    T sum = T();

    auto it2 = coeffs.rbegin(),
         it  = it2++;
    auto itd = power_diffs.rbegin();
    for (; it2 != coeffs.rend(); ++it, ++it2, ++itd)
    {
        sum += it->second;
        sum *= power(value, *itd, power_memo);
    }
    sum += it->second;
    sum *= power(value, it->first, power_memo);

    return sum;
}

template <typename T, typename deg_type>
Polynomial<T, deg_type> Polynomial<T, deg_type>::deriv() const
{
    // checkd: the operation creates a valid polynmonial.
    Polynomial d;
    for (auto it = coeffs.begin(); ++it != coeffs.end();)
        d.coeffs[it->first - 1] = it->first * it->second;
    return d;
}

// Helper functions //

template<typename T> bool equal_minus_one(const T & value)
{
    if (std::numeric_limits<T>::is_signed)
        return value == T(-1);
    else
        return false;
}

template <typename T, typename deg_type>
std::string monom_string(std::pair<deg_type, T> const & m)
{
    std::ostringstream os;
    switch (m.first)
    {
    case  0:
        os << m.second;
        break;
        
    case  1:
        if      (m.second == 1)                 os << 'x';
        else if (equal_minus_one<T>(m.second))  os << "-x";
        else                                    os << m.second << 'x';
        break;
        
    default:
        if      (m.second ==  1)               os << "x^" << m.first;
        else if (equal_minus_one<T>(m.second)) os << "-x^" << m.first;
        else                                   os << m.second << "x^" << m.first;
    }
    return os.str();
}

// IO streams //

template <typename T, typename deg_type>
std::ostream &
    operator <<(std::ostream & o, Polynomial<T, deg_type> const & p)
{
    auto it = p.coeffs.rbegin();
    
    if (it == p.coeffs.rend()) return o << T();
    
    o << monom_string(*it);
    while (++it != p.coeffs.rend())
    {
        if (it->second < T()) o << " - " << monom_string(std::make_pair(it->first, -it->second));
        else                  o << " + " << monom_string(*it);
    }
    return o;
}

template <typename T, typename deg_type>
bool read_monom(std::istringstream & is, deg_type & deg, T & coeff)
{
    using std::ios_base;
    using std::istringstream;
    
    std::string inp;
    is >> inp;

    if (inp == "") is.setstate(ios_base::failbit);
    else if (inp[0] == 'x' || (inp.size() >= 2 && inp[0] == '+' && inp[1] == 'x')) // starts with "x" or "+x"
    {
        auto x_pos = inp.find('x');
        if      (inp.back() == 'x')     { deg = 1; coeff = T(1); }
        else if (inp[x_pos + 1] != '^') is.setstate(ios_base::failbit);
        else
        {
            if (inp[0] == 'x') inp[0] = inp[1]          = ' '; // overwrite "x^"
            else               inp[0] = inp[1] = inp[2] = ' '; // overwrite "+x^"

            istringstream iss_inp(inp);
            if (iss_inp >> deg) coeff = T(1);
            else                is.setstate(ios_base::failbit);
        }
    }
    else if (inp[0] == '-' && inp[1] == 'x') // starts with "-x"
    {
        if (inp.back() == 'x') { deg = 1; coeff = T(-1); }
        else
        {
            inp[0] = inp[1] = inp[2] = ' '; // overwrite "-x^"
            istringstream iss_inp(inp);
            if (iss_inp >> deg) coeff = T(-1);
            else                is.setstate(ios_base::failbit);
        }
    }
    else
    {
        auto x_pos = inp.find('x');
        if (x_pos == std::string::npos)
        {
            istringstream iss_inp(inp);
            if (iss_inp >> coeff) deg = 0;
            else                  is.setstate(ios_base::failbit);
        }
        else if (inp.back() == 'x')
        {
            inp[x_pos] = ' ';
            istringstream iss_inp(inp);
            if (iss_inp >> coeff) deg = 1;
            else                  is.setstate(ios_base::failbit);
        }
        else
        {
            if (inp[x_pos + 1] != '^') is.setstate(ios_base::failbit);
            else
            {
                inp[x_pos] = inp[x_pos + 1] = ' ';
                istringstream iss_inp(inp);
                if (!(iss_inp >> coeff >> deg)) is.setstate(ios_base::failbit);
            }
        }
    }

    return is;
}

void str_replace(std::string & subject, std::string const & search, std::string const & replace)
{
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos)
    {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
}

template <typename T, typename deg_type>
std::istream &
    operator >>(std::istream & i, Polynomial<T, deg_type>       & p)
{
    p = T();
    
    std::string inp;
    if (i >> inp)
    {
        str_replace(inp, "+", " +");
        str_replace(inp, "-", " -");
        
        std::istringstream is(inp);
        deg_type deg;
        T coeff;
        while (read_monom(is, deg, coeff))
        {
            if (coeff != T())
                p.coeffs[deg] = std::move(coeff);
        }
        if (is.bad()) i.setstate(std::ios_base::failbit);
    }
    return i;
}

// Hash //

template <typename T, typename deg_type>
size_t Polynomial<T, deg_type>::hash() const noexcept
{
    std::hash<deg_type> dhash;
    std::hash<T>        Thash;
    
    size_t res = 0;
    for (auto & dcp : coeffs)
    {
        auto fsthash = dhash(dcp.first);
        res ^= fsthash << (4 * sizeof(size_t));
        res ^= fsthash >> (4 * sizeof(size_t));
        res ^= Thash(dcp.second);
    }
    return res;
}

} // namespace Modulus
