#include "word.h"

/******************************************************************
* Constructors of Word class.                                     *
******************************************************************/

Word::Word(){
	this->len = 0;
	err = false;
}

/* Create a Word from an array of boolean values */
Word::Word(std::vector<int> data){
	this->len = data.size();
	this->data = data;

	for(int i = 0; i < len; i++){
		this->data[i] = data[i];
	}

	err = false;
}

/******************************************************************
* Copy constructors of Word class.                                *
******************************************************************/
Word::Word(const Word& rhs){
	len = rhs.len;
	data = rhs.data;            
	for(int i = 0; i < len; i++){
		data[i] = rhs.data[i];
	}
	err = rhs.err;
}

/******************************************************************
* Destructor of Word class.                                       *
******************************************************************/
Word::~Word() = default;

/******************************************************************
* Assignment operator of Word class.                              *
******************************************************************/
Word& Word::operator = (const Word rhs) {
	len = rhs.len;
	data = rhs.data;		
	for(int i = 0; i < len; i++){
		data[i] = rhs.data[i];
	}
	err = rhs.err;
	return *this;
}

/******************************************************************
* Member functions of Word class.                                 *
******************************************************************/
int Word::get(int i){
	if(i < len) return data[i];
	err = true;
	return false;	
}

void Word::set(bool v, int i){
	if(i < len) data[i] = v;
	else err = true;
}

int Word::length(void) const{
	return len;
}

bool Word::get_err(void) const{
	
return err;
}

void Word::set_err(void)
{
	err = true;
}

/******************************************************************
* Function to compute Hamming Distance between two words          *
******************************************************************/
int hamming_distance(Word& w1, Word& w2)
{
	int hd = 0;
	if(w1.length()!= w2.length()){
		return -1;
	}

	for(int i = 0; i < w1.length(); i++){
		if(w1.get(i) != w2.get(i)) hd++;
	}

	return hd;

}

/******************************************************************
* Stream Print operator of Word class.                            *
******************************************************************/
std::ostream& operator << (std::ostream& os, const Word& v) {

	for(int i =0; i < v.len; i++){
        if(v.data[i]) os<<"1";
        else os<<"0";
	}
	return os;
}