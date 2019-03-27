#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;
int main () {

	ofstream myfile;
	myfile.open ("row2.txt");

	long int row = 2000;
	long int col = 4000;


	double r = ((double) rand() / (RAND_MAX));

	int k=0;

	for(long int i=0;i<row;i++)
	{

		myfile<<i+1;
		myfile<<' ';

		k = rand();
		k %= 50;
		k *= 100;

		for(long int j=0;j<col;j++)
		{

			r = ((double) rand() / (RAND_MAX))*100;


			myfile<<r+k<<' ';
		}
		myfile<<'\n';

	}

	myfile.close();
	return 0;
}
