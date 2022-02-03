#include "sequence1.h"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cmath>

namespace coen79_lab3{

	sequence::sequence() {
		current_index = 0;
		used = 0;
	}
	
	bool sequence::is_item() const { //checks if there is an item and either returns true or false
		if(current_index<used){
			return true;
		}
		else
			return false;
	}
	void sequence::start(){ //starts the sequence
		current_index=0;
	}
	    
	void sequence::end(){ //ends the sequence
		current_index=used-1;
	}
	
	void sequence::last(){ //last item in sequence
		current_index=CAPACITY-1;
	}

	void sequence::advance(){
		assert(is_item());
		if(current_index<used){
			current_index++;//advance index
		}
		
	}
	void sequence::retreat(){ 
		if(current_index>0)
		{
			current_index--;//decrement index
		}
	}

	void sequence::insert(const value_type & entry){//insert new entry in sequence
		assert(size()<CAPACITY);
		if(used==0){
			data[used]=entry;
			current_index=0;
		}
		else //insert a new value
		{
			for(int i = used; i >= current_index+1; --i){
				data[i] = data[i-1];
			}
			data[current_index]=entry;
		}
		used++; //increments size
	}

	void sequence::attach(const value_type & entry){ //attach new item to sequence at the end instead of start
		assert(size()<CAPACITY);

		if(!is_item()){ //if no current item
			data[used] = entry; //places in end of sequence
		}
		else 
		{ 
			for(int i = used; i < current_index+1; --i){//traverses the sequence
				data[i] = data[i-1];
			}
			data[current_index+1] = entry;
			current_index = current_index+1;
		}
		used++;
	}

	void sequence::remove_current(){//remove current variable
		assert(is_item()==true);
		for (int i=current_index+1; i<used; ++i){
			data[i-1]=data[i]; //shifts items to left
		}
		used--; //decrement used by 1
	}

	void sequence::insert_front(const value_type& entry) {//insert back into the front of sequence
		assert(size()<CAPACITY);
		current_index=0;
		insert(entry);
		used++; //increment size
	}

	void sequence::attach_back(const value_type& entry){
		assert(size()<CAPACITY);
		data[used]=entry; //attach at end
		used++; //increment size
	}

	void sequence::remove_front(){
		for (int i=0; i<used-1; ++i){
			data[i]=data[i+1]; //shift everything to left
		}
		used--;//decrement used
	}


	sequence::value_type sequence::current() const{
		assert(is_item());
		return data[current_index]; //return item
	}

	sequence::value_type sequence::operator[] (int index)const {
		assert(index<used);
		return data[index]; //return item in desired index
	}

	double sequence::mean() const{ //sum of items
		double sum=0.0;
		for (int i=0; i < used; ++i){
			sum+=data[i];
		}
		return sum/size();
	}

	double sequence::standardDeviation() const { 
		double stanDev=0.0;
		for (int i=0; i < used; ++i){
			stanDev+=pow(data[i]-mean(),2);
		}
		return sqrt(stanDev/used); 
	}

	sequence operator +=(const sequence& lhs, const sequence& rhs){
		sequence n(lhs); //copy of lhs
		int start=lhs.used;
		n.used=lhs.used+rhs.used;
		for (int i=0; i < rhs.used; ++i){ //add everything from rhs to lhs
			n.data[start]=rhs.data[i];
			start++;
		}
		return n;//return copy of lhs, executes
	}

	sequence operator +(const sequence& lhs, const sequence& rhs){
		sequence answer;
		assert(lhs.size() + rhs.size() <= sequence::CAPACITY);
		answer+=lhs;
		answer+=rhs;
		return answer;
	}

	sequence::value_type summation(const sequence &s){
		double total=0.0;
		for (int i=0; i< s.size(); ++i){
			total+=s.data[i];
		}
		return total;
	}
}








