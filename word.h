#ifndef __WORD_H__
#define __WORD_H__

#include <iostream>
#include <string>
#include <vector>


class Word
{
	std::vector<int> data;
	int len;
    bool err;

	public:
	
	/* constructors */
	Word(std::vector<int> data);

	/* copy constructor */
	Word(const Word& rhs);

	/* destructor */
	~Word();

    /* assignment operator */
    Word& operator = (const Word rhs);
	
	/*members*/
	int get(int i);
	void set(bool v, int i);
	int length(void) const;
    bool get_err(void) const;
	void set_err(void);

	/* hamming distance */
	friend int hamming_distance(Word& w1, Word& w2);

	/* stream print operator */
	friend std::ostream& operator << (std::ostream& os, const Word& v);

};

#endif