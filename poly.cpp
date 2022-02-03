#include <iostream>
#include <cmath>
#include <climits>
#include "poly.h"
#include "assert.h"


namespace coen79_lab4 {
    
    //HELPER MEMBER FUNCTION
    void polynomial::update_current_degree() {
        //POSTCONDITION: The current degree is set to the highest exponent
        //with a non-zero coeficient.
        current_degree = degree();
    }

    //CONSTRUCTOR
    polynomial::polynomial(double c, unsigned int exponent) {
        //PRECONDITION: exponent <= MAXIMUM_DEGREE
        //POSTCONDITION: This polynomial has been created with all zero
        //     coefficients, except for coefficient c for the specified exponent.
        //     When used as a default constructor (using default values for
        //     both arguments), the result is a polynomial with all zero
        //     coefficients.
        assert(exponent <= MAXIMUM_DEGREE);
        poly[MAXIMUM_DEGREE + 1];
        for(int i = 0; i <= MAXIMUM_DEGREE + 1; ++i) {
            poly[i] = 0;
        }
        poly[exponent] = c;
    }

    //MODIFICATION MEMBER FUNCTIONS
    void polynomial::assign_coef(double coefficient, unsigned int exponent) {
        //PRECONDITION: exponent <= MAXIMUM_DEGREE.
        //POSTCONDITION: Sets the coefficient for the specified exponent.
        assert(exponent <= MAXIMUM_DEGREE);
        poly[exponent] = coefficient;
    }

    void polynomial::add_to_coef(double amount, unsigned int exponent) {
        //PRECONDITION: exponent <= MAXIMUM_DEGREE.
        //POSTCONDITION: Adds the given amount to the coefficient of the
        //     specified exponent.
        assert(exponent <= MAXIMUM_DEGREE);
        poly[exponent] += amount;
    }

    void polynomial::clear() {
        //POSTCONDITION: All coefficients of this polynomial are set to zero.
        for(int i = 0; i < MAXIMUM_DEGREE + 1; ++i) {
            poly[i] = 0;
        }
        current_degree = 0;
    }

    //CONSTANT MEMBER FUNCTIONS
    polynomial polynomial::antiderivative() const {
        //PRECONDITION: degree( ) < MAXIMUM_DEGREE
        //POSTCONDITION: The return value is the antiderivative (indefinite integral)
        //     of this polynomial.
        //     NOTE: The return polynomial always has a constant term of zero.
        assert(degree() < MAXIMUM_DEGREE);
        polynomial antiderivative; 
        for(int i = (degree() + 1); i > 0; --i) {
            if(i != 0) {
                antiderivative.assign_coef((this->poly[i - 1])/i, i);
            }
        }
        return antiderivative;
    }

    double polynomial::coefficient(unsigned int exponent) const {
        //POSTCONDITION: Returns coefficient at specified exponent of this
        //     polynomial.
        //     NOTE: For exponents > MAXIMUM_DEGREE, the return value is always zero.
        if(exponent > MAXIMUM_DEGREE) {
            return 0;
        }
        return poly[exponent];
    }

    double polynomial::definite_integral(double x0, double x1) const {
        //POSTCONDITION: Returns the value of the definite integral computed from
        //     x0 to x1.  The answer is computed analytically.
        //first get the antiderivative of the polynomial calling the definite_integral function
        polynomial antiderivative = this->antiderivative(); 
        double sumx1 = 0;
        double sumx0 = 0;
        //evalutate the integral with the bounds given
        for(int i = 0; i < antiderivative.degree(); ++i) {
            sumx1 += (antiderivative.coefficient(i) * (pow(x1, i)));
            sumx0 += (antiderivative.coefficient(i) * (pow(x0, i)));
        }
        return (sumx1 - sumx0);
    }
    
    unsigned int polynomial::degree() const {
        //POSTCONDITION: The function returns the value of the largest exponent
        //     with a non-zero coefficient.
        //     If all coefficients are zero, then the function returns zero (even
        //     though, technically, this polynomial does not have a degree).
        int highest_degree = 0;
        for(int i = 0; i < MAXIMUM_DEGREE + 1; ++i) {
            if(poly[i] != 0 && i > highest_degree) {
                highest_degree = i;
            }
        }
        return highest_degree;
    }

    polynomial polynomial::derivative () const {
        //POSTCONDITION: The return value is the first derivative of this
        //     polynomial.
        polynomial derivative;
        for(int e = 0; e < MAXIMUM_DEGREE + 1; ++e) {
            derivative.assign_coef(poly[e + 1] * (e + 1), e);
            
        }
        return derivative;
    }

    double polynomial::eval(double x) const {
        //POSTCONDITION: The return value is the value of this polynomial with the
        //     given value for the variable x.
        double sum = 0;
        for(int e = 0; e <= this->degree(); ++e) {
            if(poly[e] != 0) {
                sum += (poly[e] * pow(x,e));
            }
        }
        return sum;
    }

    bool polynomial::is_zero() const {
        //POSTCONDITION: The return value is true if and only if the polynomial
        //     is the zero polynomial.
        return(degree() == 0);
    }

    unsigned int polynomial::next_term(unsigned int e) const {
        //POSTCONDITION: The return value is the next exponent n which is LARGER
        //     than e such that coefficient(n) != 0.
        //     If there is no such term, then the return value is zero.
        ++e;
        //following loop skips over 0 coefficients until it finds the next term in the array
        while(e < MAXIMUM_DEGREE + 1) {
            if(poly[e] != 0) {
                return e; //return the next exponent value (array location) if it exists
            }
            ++e;
        }
        return 0;
    }

    unsigned int polynomial::previous_term(unsigned int e) const {
        //POSTCONDITION: The return value is the next exponent n which is SMALLER
        //     than e such that coefficient(n) != 0.
        //     If there is no such term, then the return value is UINT_MAX
        //     from <climits>.
        assert(e <= MAXIMUM_DEGREE);
        //following loop skips over 0 coefficients until it finds the previous term in the array
        for(int i = e-1; i >= 0; --i) {
            if(poly[i] != 0) {
                return i; //return the previous exponent value (array location) if it exists
            }
        }
        return UINT_MAX;

    }

    double polynomial::operator ()(double x) const {
        //double operator( ) (double x) const
        //     Same as the eval member function.
        return this->eval(x);
    }

    polynomial operator +(const polynomial& p1, const polynomial& p2) {
        //POSTCONDITION: return-value is a polynomial with each coefficient
        //     equal to the sum of the coefficients of p1 & p2 for any given
        //     exponent.
        polynomial added = p1; //set the polynomial to be returned equal to the first polynomial passed 
        int e = p2.degree();
        //loop through the second polynomial adding like terms together
        while(e != UINT_MAX) {
            added.add_to_coef(p2.coefficient(e), e);
            e = p2.previous_term(e);
        }
        return added;
    }

    polynomial operator -(const polynomial& p1, const polynomial& p2) {
        //POSTCONDITION: return-value is a polynomial with each coefficient
        //     equal to the difference of the coefficients of p1 & p2 for any given
        //     exponent.
        //see comments for + operator
        polynomial sub = p1;
        int e = p2.degree();
        while(e != UINT_MAX) {
            sub.add_to_coef((p2.coefficient(e) * -1), e);
            e = p2.previous_term(e);
        }
        return sub;
    }

    polynomial operator *(const polynomial& p1, const polynomial& p2) {
        //PRECONDITION: p1.degree( ) + p2.degree( ) <= MAXIMUM_DEGREE.
        //POSTCONDITION: Each term of p1 has been multiplied by each term of p2,
        //     and the answer is the sum of all these term-by-term products.
        assert((p1.degree() + p2.degree()) <= polynomial::MAXIMUM_DEGREE);
        polynomial multiplied; //polynomial that is returned 

        //double loop. first will iterate through the first polynomial, 
        //second loop will multiply the current term of the first polynomial by each term of the second
        //      then will add like terms together and place them into polynomial multiplied
        for(int p1e = p1.degree(); p1e != UINT_MAX; p1e = p1.previous_term(p1e)) {
            for(int p2e = p2.degree(); p2e != UINT_MAX; p2e = p2.previous_term(p2e)) {
                //following if else statments will test if there is already a like term in the array or not
                if(multiplied.coefficient(p1e + p2e) == 0) {
                    multiplied.assign_coef((p1.coefficient(p1e) * p2.coefficient(p2e)), (p1e + p2e));
                }
                else {
                    multiplied.add_to_coef((p1.coefficient(p1e) * p2.coefficient(p2e)), (p1e + p2e));
                }
            }
        }
        return multiplied;
    }

    std::ostream& operator <<(std::ostream& out, const polynomial& p) {
        //POSTCONDITION: The polynomial has been printed to ostream out, which,
        //     in turn, has been returned to the calling function.
        int e = p.degree();
        if(e > 0) {
            out <<  p.coefficient(e) << "x^" << e << " ";
            e = p.previous_term(e);
            while(e != UINT_MAX) {
                //each element will be tested for the degree to output correctly 
                if(e > 1) { 
                    //if the degree n is greater than 1, the degree will be printed "^n"
                    //also test each coefficient to print either a - or + between each term
                    if(p.coefficient(e) < 0) {
                        out << " - " << (p.coefficient(e) * -1) << "x^" << e;
                    }
                    else {
                        out << " + " << p.coefficient(e) << "x^" << e;
                    }
                }
                //if n = 1. do not print the degree, just the variable x
                else if (e == 1) {
                    if(p.coefficient(e) < 0) {
                        out << " - " << (p.coefficient(e) * -1) << "x";
                    }
                    else {
                        out << " + " << p.coefficient(e) << "x";
                    }
                }
                e = p.previous_term(e); //iterator
            }
        }
        //for the final term, print out only the coefficient.
        if(p.coefficient(0) < 0) {
            out << " - " << (p.coefficient(0) * -1);
        }
        else if (p.coefficient(0) > 0) {
            out << " + " << p.coefficient(0);
        }
        return out;
    }
    

}
