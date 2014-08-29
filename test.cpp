#include <vector>
#include <iostream>
#include "src\poly.cpp"

using namespace polymath ;

int main() {
	std::vector<Term<float, 'y'>> terms1 = { Term<float, 'y'>(1, 2) };
	Polynomial<float, 'y'> p1(terms1), p2 = { Term<float, 'y'>(3, 3), Term<float, 'y'>(3, 2), Term<float, 'y'>(4, 1) };
	Polynomial<float, 'y'> p = p2 % p1 ;
	std::cout << p2.getDerivative(2) ;
	return 0 ;
}

