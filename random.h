// FILE: rand_gen.h
// CLASS PROVIDED: rand_gen
//
//  COEN 79
//  --- Behnam Dezfouli, COEN, SCU ---
//
//   (a class to generate pseudorandom numbers)
//   This class is part of the namespace coen79_2C
//
// CONSTRUCTOR for the rand_gen class:
//   rand_gen( );
//     Postcondition: The object has been initialized, and is ready to generate pseudorandom numbers.
//
// PUBLIC MODIFICATION member functions for the rand_gen class:
//  void set_seed(int newSeed)
//    Postcondition: the number new_seed is the new seed
//  int next()
//    Postcondition: The return value is the next random integer, and the return value is set as the new seed
//
// PUBLIC CONSTANT member functions for the rand_gen class:
//  void print_info()
//    Postcondition: Prints the values of the private member variables



#ifndef RAND_NUM_GEN_H
#define RAND_NUM_GEN_H

#include <iostream>
#include <cassert>
#include <stdio.h>

using namespace std;

namespace coen79_lab2 {
    
	class rand_gen
    	{
	public:
		//sets the seed, multiplier, increment and modulus
		rand_gen(int temp_seed, int temp_multiplier, int temp_increment, int temp_modulus);//generates the next seed and requires all parameters to be filled
        int next();//sets new seed
		void set_seed(int temp_seed);//changes value of the current seed. won't work if not an integer
		int get_seed() const; //getters to access private info
		int get_multiplier() const;
		int get_increment() const;
		int get_modulus() const;
		
		void print_info() const; //prints out names of private variables
	private:
		int seed;
		int multiplier;
		int increment;
		int modulus;
    	};
}


#endif