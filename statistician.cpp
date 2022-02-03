#include <cassert> 
#include <iostream>
#include <algorithm>
#include "statistician.h"

using namespace std;
using namespace coen79_lab2;

namespace coen79_lab2{                       //header definition

	void statistician::next(double r){       //add number to sequence
		total += r;                          //received is r, add to total number
		used++;                              //my counter
		
		if(length()==1){                     //if only one value
			min = r; 
			max = r;
		}else{
			if(r > max){
				max = r;                     //change max and min here so the functions can just return the values
			}else if(r < min){
				min = r;
			}
		}
	}

	double statistician::sum() const{
		return total;
	}
	
	double statistician::mean() const{
		assert(length()>0);
		return sum()/used;                   //returns sum
	}
	
	double statistician::minimum() const{
		assert(length()>0);
		return min;
	}

	double statistician::maximum() const{
		return max;
	}
	
	void statistician::reset(){
		total = 0.0; used = 0.0;             //resets back to 0.0
	}
	
	bool operator ==(const statistician& s1, const statistician& s2){
		return(s1.length()==s2.length() &&   //makes sure all are equal
			s1.mean()==s2.mean() &&
			s1.sum()==s2.sum() &&
			s1.maximum()==s2.maximum() &&
			s1.minimum()==s2.minimum());
	}
	
	statistician operator+(const statistician& s1, const statistician&s2){
		statistician c;
		c.total = s1.total + s2.total;
		c.used = s1.used + s2.used;
		c.min = (s1.min < s2.min) ?s1.min:s2.min;
		c.max = (s1.max > s2.max)?s1.max:s2.max;

		return c;
	}
	
	statistician operator *(double scale, const statistician& s){
		statistician c;
		c.total = (scale * s.total);
		c.used = s.used;
		c.min = (scale * s.max);
		c.max = (scale * s.min);
		return c;
	}

}
	
void print_statistician_info(statistician & s){//function to print statistician info
	cout << "Sum: " << s.sum() << endl;
	cout << "Mean: " << s.mean() << endl;
	cout << "Smallest: " << s.minimum() << endl;
	cout << "Largest: " << s.maximum() << endl;
	cout << "Length: " << s.length() << endl;
	cout << endl;
}
