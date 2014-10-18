#ifndef __POLY_MATH__
#define __POLY_MATH__

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <functional>

namespace polymath
{

template <typename> struct check_for_number ;
template <> struct check_for_number<int> {} ;
template <> struct check_for_number<long int> {} ;
template <> struct check_for_number<long long int> {} ;
template <> struct check_for_number<float> {} ;
template <> struct check_for_number<double> {} ;
template <> struct check_for_number<unsigned int> {} ;
template <> struct check_for_number<unsigned long int> {} ;
template <> struct check_for_number<unsigned long long int> {} ;
template <> struct check_for_number<short int> {} ;

unsigned int factorial(unsigned int n)
{
	auto ret = 1u ;
	if(n == 1u || n == 0u)
		return ret ;
	for(auto i = 2u; i <= n; ++i)
		ret *= i ;
	return ret ;
}

unsigned int combination(unsigned int n, unsigned int r)
{
	assert(n >= r && "n must be greater than r");
	auto ret = 1u ;
	if(r == 0 || n == r)
		return ret ;
	for(auto i = n; i >= n - r + 1u; --i)
		ret *= i ;
	return ret / factorial(r);
}

template <typename, char> class Polynomial ;

template <typename Number, char var = 'x'> class Term
{
	check_for_number<Number> check;
	Number _coefficient ;
	unsigned int _power ;
	
public:
	Term(Number num, unsigned int exp) :
		_coefficient{num}, _power{exp} {}
		
	bool operator<(Term term) const { return _power < term._power; }
	
	bool similar(Term term) { return _power == term._power; }
	Number getCoefficient() const { return _coefficient; }
	unsigned int getPower() const { return _power; }
		
	friend class Polynomial<Number, var> ;
	friend std::ostream& operator<<(std::ostream& stream, const Term<Number, var>& term) {
		//if(term._coefficient != 0)
		//{
			char sign = (term._coefficient >= 0) ? '+' : '\0' ;
			if(term._power != 0)
				stream << sign << term._coefficient << var << "^" << term._power ;
			else
				stream << sign << term._coefficient ;
		//}
		return stream ;
	}
};

template <typename Number, char var = 'x'>
class Polynomial
{
	std::vector<Term<Number, var>> _terms;
	unsigned int _degree ;
	void adjust();
	bool binomialCapable() const ;
	static Number zero ;
	
public:
	Polynomial();
	Polynomial(const Polynomial<Number, var>&) = delete ;
	Polynomial(Polynomial<Number, var>&&) = default ;
	Polynomial& operator=(Polynomial<Number, var>&&) = default ;
	Polynomial(const std::vector<Term<Number, var>>&);
	Polynomial(const std::initializer_list<Term<Number, var>>);
		
	void addTerm(const Term<Number, var>&);
	
	unsigned int degree() const { return _degree; }
	unsigned int validTerms() const ;
	Number operator()(Number) const ;
	Term<Number, var> operator[](unsigned int) const ;
	Polynomial<Number, var> getDerivative(unsigned int) const ;
	Polynomial<Number, var> getIntegral(Number) const ;
	Polynomial<Number, var> clone() const ;
	std::vector<Number> getRealRoots(unsigned int, Number, Number, Number) const ;
	Polynomial<Number, var> transform(const std::function<Number(Number)>&) const;
	
	template <typename number, char v> friend Polynomial<number, v> operator+(const Polynomial<number, v>&, const Polynomial<number, v>&);
	template <typename number, char v> friend Polynomial<number, v> operator-(const Polynomial<number, v>&, const Polynomial<number, v>&);
	template <typename number, char v> friend Polynomial<number, v> operator*(const Polynomial<number, v>&, const Polynomial<number, v>&);
	template <typename number, char v> friend Polynomial<number, v> operator/(const Polynomial<number, v>&, const Polynomial<number, v>&);
	template <typename number, char v> friend Polynomial<number, v> operator%(const Polynomial<number, v>&, const Polynomial<number, v>&);
	template <typename number, char v> friend Polynomial<number, v> operator^(const Polynomial<number, v>&, unsigned int);
	
	friend std::ostream& operator<<(std::ostream& stream, const Polynomial<Number, var>& poly) {
		for(auto& term : poly._terms)
			stream << term << " " ;
		return stream ;
	}
};

Number Polynomial<Number, var>::zero = static_cast<Number>(0);

template <typename Number, char var>
Polynomial<Number, var>::Polynomial()
{
	_degree = 0 ;
	_terms.push_back(Term<Number, var>(0, 0));
}

template <typename Number, char var>
Polynomial<Number, var>::Polynomial(const std::vector<Term<Number, var>>& _terms_)
{
	std::vector<Term<Number, var>> temp(_terms_);
	std::sort(temp.begin(), temp.end());
	_degree = 0 ;
	_terms.push_back(Term<Number, var>(0, 0));
	for(auto i = 0; i < temp.size(); ++i)
		addTerm(temp[i]);
}

template <typename Number, char var>
Polynomial<Number, var>::Polynomial(std::initializer_list<Term<Number, var>> _terms_)
{
	std::vector<Term<Number, var>> temp(_terms_);
	std::sort(temp.begin(), temp.end());
	_degree = 0 ;
	_terms.push_back(Term<Number, var>(0, 0));
	for(auto i = 0; i < temp.size(); ++i)
		addTerm(temp[i]);
}

template <typename Number, char var>
void Polynomial<Number, var>::adjust()
{
	while(_degree > 0)
	{
		if(_terms[_degree]._coefficient == 0)
		{
			_terms.erase(_terms.begin() + _degree);
			--_degree ;
		}
		else
			break ;
	}
	if(_degree == 0)
		_terms.push_back(Term<Number, var>(0, 0));
}

template <typename Number, char var>
Polynomial<Number, var> Polynomial<Number, var>::clone() const 
{
	Polynomial<Number, var> poly ;
	poly._terms = _terms ;
	poly._degree = _degree ;
	return poly ;
}

template <typename Number, char var>
Polynomial<Number, var> Polynomial<Number, var>::getDerivative(unsigned int order) const
{
	Polynomial<Number, var> derivative ;
	for(auto i = _degree; i >= order; --i)
	{
		auto coefficient = _terms[i]._coefficient ;
		for(auto j = 0u; j < order; ++j)
			coefficient *= static_cast<Number>(i - j);
		derivative.addTerm(Term<Number, var>(coefficient, i-order));
	}
	return derivative ;
}

template <typename Number, char var>
Polynomial<Number, var> Polynomial<Number, var>::getIntegral(Number constant) const
{
	Polynomial<Number, var> integral ;
	for(auto i = 0; i < _degree; ++i)
	{
		integral.addTerm(Term<Number, var>(
			_terms[i]._coefficient / static_cast<Number>(_terms[i]._exponent + 1),
			_terms[i]._exponent + 1));
	}
	integral.addTerm(Term<Number, var>(constant, 0u));
	return integral;
}

template <typename Number, char var>
bool Polynomial<Number, var>::binomialCapable() const
{
	bool ret = _terms[0]._coefficient > 0 ;
	auto count = 0u ;
	for(auto i = 1u; i < _terms.size(); ++i)
		if(_terms[i]._coefficient > static_cast<Number>(0))
			++count ;
	return ret && (count == 1u);
}

template <typename Number, char var>
void Polynomial<Number, var>::addTerm(const Term<Number, var>& term)
{
	if(_degree >= term._power)
	{
		_terms[term._power]._coefficient += term._coefficient ;
	}
	else
	{
		auto diff = term._power - _degree ;
		for(auto i = 1; i < diff; ++i)
			_terms.push_back(Term<Number, var>(0, _degree + i));
		_terms.push_back(term);
		_degree = term._power ;
	}
}

template <typename Number, char var>
unsigned int Polynomial<Number, var>::validTerms() const
{
	auto count = 0u ;
	for(auto& term : _terms)
		count += (term._coefficient != static_cast<Number>(0)) ;
	return count ;
}

template <typename Number, char var>
Polynomial<Number, var> Polynomial<Number, var>::transform(const std::function<Number(Number)>& f) const
{
	Polynomial<Number, var> ret = clone();
	for(auto term : ret._terms)
		term._coefficient = f(term._coefficient);
	ret.adjust();
	return ret ;
}

template <typename Number, char var>
std::vector<Number> Polynomial<Number, var>::getRealRoots(unsigned int iter_count, Number min, Number start, Number end) const
{
	std::vector<Number> roots ;
	Polynomial<Number, var> temp = clone(), d1, d2;
	Number tempval, current, approx ; 
	if(_terms[0]._coefficient == static_cast<Number>(0))
	{
		roots.push_back(static_cast<Number>(0));
		for(auto term : roots._terms)
			term._power -= 1u ;
		roots.adjust();
	}
	tempval = temp(start);
	for(auto i = (++start); i < end; ++i) {
		current = temp(i);
		if(tempval < 0 && current > 0 ||
		   tempval > 0 && current < 0)
		{
			approx = tempval ;
			d1 = temp.getDerivative(1);
			d2 = temp.getDerivative(2);
			for(auto i = 0u; i < iter_count && a > min; ++i)
			{
				auto val1 = temp(approx) ;
				if(val1 < min)
					break ;
				auto val2 = d1(approx), val3 = d2(approx);
				auto G = val2/val1 ;
				auto H = (G*G) - val3/val1 ;
				auto n = static_cast<Number>(temp._degree) ;
				auto t = static_cast<Number>(std::sqrt(static_cast<double>((n - 
					static_cast<Number>(1)) * ((n*H) - (G*G)))));
				auto a = n/(G + (G >= 0)*t - (G < 0)*t);
				approx -= a ;
			}
			roots.push_back(approx);
		}
		else
			tempval = current;
	}
	return roots ;
}

template <typename Number, char var>
Number Polynomial<Number, var>::operator()(Number n) const
{
	Number result = _terms[0]._coefficient;
	Number inc = static_cast<Number>(1);
	for(auto i = 1u; i < _terms.size(); ++i)
	{
		inc *= n ;
		result += _terms[i]._coefficient * inc ;
	}
	return result ;
}

template <typename number, char v>
Polynomial<number, v> operator+(const Polynomial<number, v>& lhs, const Polynomial<number, v>& rhs)
{
	Polynomial<number, v> sum;
	auto min = std::min(lhs.degree(), rhs.degree()) ;
	auto i = 0 ;
	for(; i <= min; ++i)
		sum.addTerm(Term<number, v>(lhs._terms[i].getCoefficient() + rhs._terms[i].getCoefficient() , i)) ;
	while(i < lhs._terms.size())
	{
		sum.addTerm(Term<number, v>(lhs._terms[i].getCoefficient() , i));
		++i ;
	}
	while(i < rhs._terms.size())
	{
		sum.addTerm(Term<number, v>(rhs._terms[i].getCoefficient() , i));
		++i ;
	}
	sum._degree = i - 1 ;
	sum.adjust();
	return sum ;
}

template <typename number, char v>
Polynomial<number, v> operator-(const Polynomial<number, v>& lhs, const Polynomial<number, v>& rhs)
{
	Polynomial<number, v> diff;
	auto min = std::min(lhs.degree(), rhs.degree()) ;
	auto i = 0 ;
	for(; i <= min; ++i)
		diff.addTerm(Term<number, v>(lhs._terms[i].getCoefficient() - rhs._terms[i].getCoefficient() , i)) ;
	while(i < lhs._terms.size())
	{
		diff.addTerm(Term<number, v>(lhs._terms[i].getCoefficient() , i));
		++i ;
	}
	while(i < rhs._terms.size())
	{
		diff.addTerm(Term<number, v>(-rhs._terms[i].getCoefficient() , i));
		++i ;
	}
	diff._degree = i - 1 ;
	diff.adjust();
	return diff ;
}

template <typename number, char v>
Polynomial<number, v> operator*(const Polynomial<number, v>& lhs, const Polynomial<number, v>& rhs)
{
	Polynomial<number, v> product;
	for(auto& rterm : rhs._terms)
		for(auto& lterm : lhs._terms)
			product.addTerm(Term<number, v>(lterm.getCoefficient() * rterm.getCoefficient(),
			                                lterm.getPower() + rterm.getPower()));
	return product ;
}

template <typename number, char v>
Polynomial<number, v> operator/(const Polynomial<number, v>& lhs, const Polynomial<number, v>& rhs)
{
	Polynomial<number, v> quotient, dividend(lhs._terms);
	Term<number, v> last = rhs._terms[rhs._terms.size() - 1] ;
	if(rhs.degree() > lhs.degree())
		return quotient ;
	// monomial division is simplified.
	if(rhs.validTerms() == 1u)
	{
		auto i = dividend._degree ;
		auto power = rhs._terms[rhs._degree].getPower() ;
		auto coefficient = rhs._terms[rhs._degree].getCoefficient() ;
		for(; i > 0; --i)
			if(dividend._terms[i].getPower() >= power)
				quotient.addTerm(Term<number, v>(
					dividend._terms[i].getCoefficient() / coefficient,
					dividend._terms[i].getPower() - power
				));
		if(dividend._terms[i].getPower() >= rhs._terms[i].getPower())
			quotient.addTerm(Term<number, v>(
				dividend._terms[i].getCoefficient() / coefficient,
				dividend._terms[i].getPower() - power
			));
		return quotient ;
	}
	while(dividend.degree() >= rhs.degree())
	{
		auto term = dividend._terms[dividend._terms.size() - 1]; // term to cancel
		auto multiplier = term.getCoefficient()/last.getCoefficient() ;
		auto power = term.getPower() - last.getPower() ;
		quotient.addTerm(Term<number, v>(multiplier, power));
		dividend = dividend - (Polynomial<number, v>({ Term<number, v>(multiplier, power) }) * rhs) ;
	}
	return quotient ;
}

template <typename number, char v>
Polynomial<number, v> operator%(const Polynomial<number, v>& lhs, const Polynomial<number, v>& rhs)
{
	Polynomial<number, v> remainder, dividend(lhs._terms);
	Term<number, v> last = rhs._terms[rhs._terms.size() - 1] ;
	if(rhs.degree() >= lhs.degree())
		return remainder;
	// monomial remainder is simplified.
	if(rhs.validTerms() == 1u)
	{
		auto mpower = rhs._terms[rhs._degree].getPower() ;
		for(auto i = 0; i < dividend._degree; ++i)
			if(dividend._terms[i].getPower() < mpower)
				remainder.addTerm(dividend._terms[i]);
		return remainder ;
	}
	while(dividend.degree() >= rhs.degree())
	{
		auto term = dividend._terms[dividend._terms.size() - 1]; // term to cancel
		auto multiplier = term.getCoefficient()/last.getCoefficient() ;
		auto power = term.getPower() - last.getPower() ;
		dividend = dividend - (Polynomial<number, v>({ Term<number, v>(multiplier, power) }) * rhs) ;
	}
	return dividend ;
}


template <typename number, char v>
Polynomial<number, v> operator^(const Polynomial<number, v>& poly, unsigned int power)
{
	Polynomial<number, v> result;
	if(poly.validTerms() == 1u)
	{
		for(auto i = 0u; i < poly._degree; ++i)
			result.addTerm(Term<number, v>(
				static_cast<number>(std::pow(poly._terms[0].getCoefficient(), power)),
				poly._terms[0].getPower() + power
			));
	}
	else if(poly.binomialCapable())
	{
		for(auto i = 0; i <= power; ++i)
			result.addTerm(Term<number, v>(
				std::pow(poly._terms[0].getCoefficient(), power - i) *
				static_cast<number>(combination(power, i)),
				i * poly._terms[1].getPower()
			));
	}
	else
	{
		result = poly ;
		for(auto i = 1u; i < power; ++i)
			result = result * poly ;
	}
	return result ;
}

}

#endif
