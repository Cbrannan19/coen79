//random.cpp
//class that creates pseudorandom numbers
//each time a new random number is computed, the value of the seed is changed to the new number

#include <iostream>
#include "random.h"

using namespace std;
using namespace coen79_lab2;

namespace coen79_lab2{

	
	rand_gen::rand_gen(int temp_seed, int temp_multiplier, int temp_increment, int temp_modulus){
		this->seed = temp_seed;//temporary seed
		this->multiplier = temp_multiplier;//temporary multiplier
		this->increment = temp_increment;
		this->modulus = temp_modulus;
	}
	
	void rand_gen::set_seed(int temp_seed){
		seed = temp_seed;
		return;
	}

	int rand_gen::next(){
		int new_seed;
		new_seed = (((multiplier * seed) + increment) % modulus);
		seed = new_seed;
		return new_seed;
	}
	
	int rand_gen::get_seed() const{
		return seed;
	}

	int rand_gen::get_multiplier() const{
		return multiplier;
	}

	int rand_gen::get_increment() const{
		return increment;
	}

	int rand_gen::get_modulus() const{
		return modulus;
	}

	void rand_gen::print_info() const{
		cout << "Seed: " << get_seed() << endl;
		cout << "Multiplier: " << get_multiplier() << endl;
		cout << "Increment: " << get_increment() << endl;
		cout << "Modulus: " << get_modulus() << endl;
	}
}
