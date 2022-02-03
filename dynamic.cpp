#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <bits/stdc++.h> 

//Declaring standard namespace std to make functions and declarations shorter.
using namespace std; 

//Size of all characters with one additional character for punctuation.
const int SIZE = 27;

//The Trie Node each hold an array of children (the children nodes of the tree) and frequency value of 
//each of the words. With one boolean being whether a given node is the end of a word.
struct TrieNode 
{ 
    struct TrieNode *children[SIZE]; 
    bool isEndOfWord; 
    int frequency;
}; 

//Function: This function is the constructor for the Trie tree data structure 
//that will be used in the program.
//
//Efficiency: O(SIZE) with SIZE being all the possible lowercase and uppercase
//letters along with any punctuation.
//
//Return type: A pointer to the newly created Trie tree data structure.
//
//Precondition: None, enough memory to allocate for the Trie Tree.
//
//Postcondition: A new Trie tree data structure with children set to NULL and
//the frequency set to 0.
struct TrieNode *getNode(void)
{ 
    struct TrieNode *pNode =  new TrieNode; 
    pNode->isEndOfWord = false;
    for (int i = 0; i < SIZE; i++)
	{
        pNode->children[i] = NULL;
	}
    return pNode; 
} 

//Function: This function checks whether a character is a punctutation mark or 
//not.
//
//Efficiency: O(1) by only checking one character.
//
//Return type: A boolean of true or false.
//
//Precondition: A valid character "c".
//
//Postcondition: True or false if a given character is punctuation or not.
bool punct(char c) 
{
	return ((c == '.')||(c == '!')||(c == '?')||(c == ',')||(c== ';'));
}

//Function: This function adds words into the Trie tree and marks whether a given node
//is the last letter of a given word in a Trie Tree and then sets the frequency of the 
//word in the Trie tree. 
//
//Efficiency: O(line length) with line length being how long the string being passed to 
//the function is.
//
//Return type: Boolean, True or False whether a word has been added to the Trie Tree. 
//
//Precondition: A valid root to a Trie tree, a valid string, and the given frequency value.
//
//Postcondition: True or false whether or not the word is an English word.
bool insert(struct TrieNode *root, string line, double frequency)
{
	struct TrieNode *it = root; 
	if(line.length() == 0)
		return false;
	for (int i = 0; i < line.length(); i++) 
	{ 
		if((line.length() == 1)) //only one letter words are not added to the Trie Tree (they are ignored) 
		{
			return false;
		}
		if(isalpha(line[i])||(line[i] == '\'')) //the if statement goes to the correct array location with its given ascii value
		{
			int index = line[i];
			if (isupper(line[i]))
				index = line[i] - 'A';
			else if (isalpha(line[i]))
				index = (line[i] - 'a');
			else if (line[i] == '\'')
				index = 26;
			if (!it->children[index])  //if there is no node in the child position, create a Trie Node
			    it->children[index] = getNode(); 
			it = it->children[index]; 
			
		}
		else //if it isn't an alphanumeric character, do not add it to the tree at all
		{
			return false;
		}
	}
	if(it->isEndOfWord) //if at the end of a word, then do not add anything
	{
		return true;
	}
	//once at the end of a given word, set the frequency value 
	it->isEndOfWord = true;
	it->children[27] = getNode(); 
	it = it->children[27];
	it->isEndOfWord = true; 
	it->frequency = frequency;
	return true;
} 

//Function: This function finds whether the given string is a word in the Trie Tree.
//
//Efficiency: O(line length) with line length being how long the string being passed to 
//the function is.
//
//Return type: Boolean, True or False whether a word is in the Trie Tree or not. Test
//
//Precondition: A valid root to a Trie tree, and valid string.
//
//Postcondition: True or false whether or not the word is in the Trie Tree.
bool search(struct TrieNode *root, string line) 
{ 
	struct TrieNode *it = root; 
	if((line == "I")||(line == "a")||(line == "A")) //if a single letter, then return true in these cases
	{
		return true;
	}
	for (int i = 0; i < line.length(); i++) //checks through each possible character 
	{
		if(i > 0)
		{	if(isupper(line[i]))
			{
				return false;
			}
		}
		if(isalpha(line[i])||(punct(line[i]))||(line[i] == '\'')) //checks through the same cases as insertion to find array location
		{
			int index = line[i];
			if (isupper(line[i]))
				index = line[i] - 'A';
			else if (islower(line[i]))
				index = line[i] - 'a';
			else if (line[i] == '\'')
				index = 26;
			else if (punct(line[i]))
				index = 27;
			if (!it->children[index]) //if the char of the string isn't created, then the word is not in Trie Tree, return false
			{
			    return false;
			}
			it = it->children[index]; 
		}
		else
		{
			return false;
		}
	} 
	return (it != NULL && it->isEndOfWord); 
} 

//Function: This function finds whether the given string is a word in the Trie Tree and then 
//returns its frequency.
//
//Efficiency: O(line length) with line length being how long the string being passed to 
//the function is.
//
//Return type: An integer, the given frequency of a word in the Trie Tree. 
//
//Precondition: A valid root to a Trie tree, and valid string.
//
//Postcondition: The frequency of the word in the Trie Tree.
int searchint(struct TrieNode *root, string line) 
{ 
    struct TrieNode *it = root; 
    assert(search(root, line));
    for (int i = 0; i < line.length(); i++) //checks through each possible character 
    { 
		if((line == "I")||(line == "a")||(line == "A")) //returns a frequency of "0" for these three characters (since they are veyr frequent
		{
			return 0;
		}			
     	if(isalpha(line[i])||(line[i] == '\'')||(punct(line[i]))) //checks through the same cases as insertion to find array location
		{
			int index = line[i];
			if (isupper(line[i]))
				index = line[i] - 'A';
			else if (islower(line[i]))
				index = line[i] - 'a';
			else if (line[i] == '\'')
				index = 26;
			else if (punct(line[i]))
				break;
			it = it->children[index]; 
		}
    }
	it = it->children[27]; 
    return it->frequency; 
}


//Function: This function breaks a given string into each of its given components and tries to 
//find all the variations that exist for that given string. (this is the dynamic programming function)
//
//Return type: Void, no return type.
//
//Precondition: A valid root to a Trie tree, and valid string, string length, a vector address for the rank and 
//the current rank of the word.
//
//Postcondition: Files created for each possible variation of the given corrupt text file.
void sentenceBreakUtil(struct TrieNode *root, string str, int n, string result, vector <int> &rank, int current_rank) 
{ 
	for (int i = 1; i <= n; i++) 
	{ 
		string prefix = str.substr(0, i); 
		if (search(root, prefix)) 
		{ 
			int new_rank = current_rank + searchint(root, prefix);
			if (i == n)
			{ 
				result += prefix; 
				rank.push_back(new_rank); //pushing the frequency value into a vector object 
				ostringstream ss;
				ss << "solution_" << new_rank << ".txt"; //giving the file name the frequency value
				string name = ss.str();
   				ofstream outfile(name);
				outfile << result << endl; //putting the possible output string into a file
				return;
			}
			sentenceBreakUtil(root, str.substr(i, n - i), n - i, result + prefix + " ", rank, new_rank); 
		} 
	}
}



//Function: The function called by main to start the process of dynamic programming and 
//create all the possible output files.
//
//Return type: Void, no return type.
//
//Precondition: A valid root to a Trie tree, and valid string, and the address of a vector object.
//
//Postcondition: Same as "sentenceBreakUtil" function, the files created for each 
//possible variation of the given corrupt text file.
void sentenceBreak(struct TrieNode *root, string str, vector <int> &rank) 
{
	return sentenceBreakUtil(root, str, str.size(), "", rank, 0);
}


//Function: This function checks to see whether a given string can be broken into different words.
//
//Return type: Boolean, True or False whether a given string can be broken into different words.
//
//Precondition: A valid root to a Trie tree, and valid string.
//
//Postcondition: Files created for each possible variation of the given corrupt text file.
bool sentenceCanBreak(struct TrieNode *root, string str) 
{ 
    int size = str.size(); 
    if (size == 0) //base case of recursion
	{
		return true;
	}
    if((str == "I")||(str == "A")||(str == "a")) //check for single characters
	{
		return true;
	}
    for (int i = 1; i <= size; i++) //continually search and then reduce size to see if string can have multiple possible words
    { 
        if (search(root, str.substr(0, i) ) && sentenceCanBreak(root, str.substr(i, size) ))			
		return true;
    } 
    return false; 
} 


//Main function
int main(int argc, const char * argv[]) 
{ 
	string word; 
	ifstream file1;
 	struct TrieNode *root = getNode(); //create Trie tree root
	file1.open(argv[1]); //open dictionary as first given parameter of the program
	if(file1 == NULL) //checking if dictionary is valid
	{
		cout << "Error Opening Dictionary!" << endl;
		return 0;
	}
	else //creating the Trie tree, each line is a word, and its frequency is the line number starting from 1 to 
		 //the end of the file (which should be around 40,000 given the number of words in dictionary3.txt)
	{
		int counter = 1;
		string word;
		cout << "Creating Trie Tree..." << endl;
		while(getline(file1, word))
		{
			if(word[1] != '\n')
			{
				word.erase(word.length()-1);
				bool temp = insert(root, word, counter);
				if(temp)
				{
					counter++;
				}
			}
		}
		cout << "Trie Tree created!" << endl;
	}
	string file_text;
	ifstream file2;	
	file2.open(argv[2]);
	if(file2 == NULL) //checking if corrupt file is valid
	{
		cout << endl << "Error Opening Input File" << endl;
		return 0;
	}
	else
	{
		cout << "Generating all possible files..." << endl; 
		getline(file2, file_text); //putting all characters into one string from corrupt file
		file2.ignore(); //stop putting characters into "file_text" once you reach End of File
		if(word.length() != 0)
			word.erase(word.length() - 1);
		if(sentenceCanBreak(root, file_text))
		{
			vector <int> rank; //create a vector object for the ranking of the words
			sentenceBreak(root, file_text, rank); //do the dynamic programming
			cout << endl << "The file extensions are the sum of the rankings of the words in the solution file." << endl;
			cout << endl << "Files Generated." << endl << "The rankings from most likely to less likely are: " << endl;
			sort(rank.begin(), rank.end()); //sort the rankings of the different files
			int i = 1; //Have the number ranking start at 1.
			for(std::vector<int>::iterator it = rank.begin(); it != rank.end(); ++it) //for loop through each ranking and print it out on its own line.
			{	
				cout << i << ": " << "\t";
				cout << "solution_" << *it << endl;
				++i;
			}
		}
		else
		{
			cout << "Invalid text!" << endl;		
		}
	}
    return 0; 
}