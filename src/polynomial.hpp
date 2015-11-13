#pragma once

// Compile with clang++-3.0 -std=c++11

#include <iostream>
#include <sstream>

#include <limits>

#include <vector>
#include <map>

#include <algorithm>
#include <functional>

namespace Modulus
{
// Forward declaration of the types and non-inline template friend functions.

template <typename T, typename deg_type = size_t>   class Polynomial;

template <typename T, typename deg_type>    deg_type                deg         (Polynomial<T, deg_type> const &);
template <typename T, typename deg_type>    Polynomial<T, deg_type> operator << (Polynomial<T, deg_type> const &,  deg_type);
template <typename T, typename deg_type>    Polynomial<T, deg_type> operator >> (Polynomial<T, deg_type> const &,  deg_type);
template <typename T, typename deg_type>    Polynomial<T, deg_type> operator +  (Polynomial<T, deg_type> const &,  Polynomial<T, deg_type> const &);
template <typename T, typename deg_type>    Polynomial<T, deg_type> operator *  (Polynomial<T, deg_type> const &,  Polynomial<T, deg_type> const &);
template <typename T, typename deg_type>    std::ostream &          operator << (std::ostream &,                   Polynomial<T, deg_type> const &);
template <typename T, typename deg_type>    std::istream &          operator >> (std::istream &,                   Polynomial<T, deg_type>       &);

template <typename T, typename deg_type = size_t>
class Eks
{
    deg_type deg;
    constexpr Eks(deg_type deg) : deg(deg) {  }
    
public:
    constexpr Eks() : deg(1) { }
    
    operator Polynomial<T, deg_type>() const noexcept
    {
        return Polynomial<T, deg_type>(T(1), deg);
    }
    
    friend constexpr deg_type deg(Eks e) { return e.deg; }
    
    constexpr Eks operator ^ (deg_type d) const { return Eks(deg * d); }
    constexpr Eks operator * (Eks      b) const { return Eks(deg + b.deg); }

    Polynomial<T, deg_type> operator-() const { return -Polynomial<T, deg_type>(*this); }
    
    friend Polynomial<T, deg_type> operator+(Eks a, T const & b) { return static_cast<Polynomial<T, deg_type>>(a) + Polynomial<T, deg_type>(b); }
    friend Polynomial<T, deg_type> operator-(Eks a, T const & b) { return static_cast<Polynomial<T, deg_type>>(a) - Polynomial<T, deg_type>(b); }
    
    friend Polynomial<T, deg_type> operator+(Eks a, Eks b) { return static_cast<Polynomial<T, deg_type>>(a) + static_cast<Polynomial<T, deg_type>>(b); }
    friend Polynomial<T, deg_type> operator-(Eks a, Eks b) { return static_cast<Polynomial<T, deg_type>>(a) - static_cast<Polynomial<T, deg_type>>(b); }
};


template <typename T, typename deg_type> /* default deg_type = size_t*/
class Polynomial // non-literal type (non-trivial destructor)
{
    static_assert(std::numeric_limits<deg_type>::is_integer && !std::numeric_limits<deg_type>::is_signed,
            "deg_type must be an unsigned integer.");
    
private:
    std::map<deg_type, T> coeffs;
    
public:
    //static constexpr Eks<T, deg_type> X = Eks<T, deg_type>();

    Polynomial() { }
    Polynomial(T const &, deg_type deg = 0);
    
    static  Polynomial      fromCoeffVector(std::vector<T> const & coeffs);
    
    friend  deg_type        deg<>(Polynomial const &);
            Polynomial      with_monic(deg_type dg = 0) const;
    
            T               leading_coeff() const { return coeffs.rbegin()->second; }
    
            bool is_zero() const                                          { return coeffs.empty(); }
    friend  bool operator == (Polynomial const & p, Polynomial const & q) { return p.coeffs == q.coeffs; }
    friend  bool operator != (Polynomial const & p, Polynomial const & q) { return p.coeffs != q.coeffs; }
    friend  bool operator <  (Polynomial const & p, Polynomial const & q) { return p.coeffs <  q.coeffs; }
    
    friend  Polynomial      operator << <>  (Polynomial const & p,  deg_type d);
    friend  Polynomial      operator >> <>  (Polynomial const & p,  deg_type d);
            Polynomial &    operator<<=     (                       deg_type d) & { return *this = *this << d; }
            Polynomial &    operator>>=     (                       deg_type d) & { return *this = *this >> d; }
    
    static  std::pair<Polynomial, Polynomial> divmod(Polynomial         a, Polynomial const & b);
    static  void                              divmod(Polynomial const & a, Polynomial const & b, Polynomial & q, Polynomial & r);
    
    // Must be friend, see http://stackoverflow.com/q/29318021/3273130?sem=2
    friend  Polynomial      operator +      (Polynomial p) { return p; }
            Polynomial      operator -      () const;
    
    friend  Polynomial      operator +      (T          const & t,   Polynomial         q) { return q += Polynomial(t); }
    friend  Polynomial      operator +      (Polynomial         p,   T          const & t) { return p += Polynomial(t); }
    friend  Polynomial      operator -      (T          const & t,   Polynomial         q) { return q -= Polynomial(t); }
    friend  Polynomial      operator -      (Polynomial         p,   T          const & t) { return p -= Polynomial(t); }
    friend  Polynomial      operator *      (T          const & t,   Polynomial         q) { return q *= t; }
    friend  Polynomial      operator *      (Polynomial         p,   T          const & t) { return p *= t; }
    
    friend  Polynomial      operator + <>   (Polynomial const & p,   Polynomial const & q);
    friend  Polynomial      operator * <>   (Polynomial const & p,   Polynomial const & q);
    friend  Polynomial      operator -      (Polynomial const & p,   Polynomial const & q) { return p + (-q); }
    friend  Polynomial      operator /      (Polynomial const & p,   Polynomial const & q) { return divmod(p, q).first; }
    friend  Polynomial      operator %      (Polynomial const & p,   Polynomial const & q) { return divmod(p, q).second; }
    
            Polynomial &    operator *=     (T          const & t) &;
    
            Polynomial &    operator +=     (Polynomial const & q) &;
            Polynomial &    operator -=     (Polynomial const & q) & { return *this += (-q); }
            Polynomial &    operator *=     (Polynomial const & q) & { return *this = *this * q; }
            Polynomial &    operator /=     (Polynomial const & q) & { return *this = *this / q; }
            Polynomial &    operator %=     (Polynomial const & q) & { return *this = *this % q; }
    
            T const &       at              (deg_type) const &;
            T               at              (deg_type)      &&;
            
            T               operator ()     (T const &) const;
            
            Polynomial      deriv           () const;
    
    friend  std::ostream &  operator << <>  (std::ostream &, Polynomial const &);
    friend  std::istream &  operator >> <>  (std::istream &, Polynomial       &);
    
            size_t       hash() const noexcept;
};

} // namespace Modulus

namespace std
{

template<typename T, typename deg_type>
struct hash<Modulus::Polynomial<T, deg_type>>
{
    inline size_t operator()(Modulus::Polynomial<T, deg_type> const & p) const noexcept
    { return p.hash(); }
};

} // namespace std


// As polynomial is a temlate class, it is necessary to include the tpp-file as well.
#include "Polynomial.tpp"
