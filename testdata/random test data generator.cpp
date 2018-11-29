//random test data generator
#include<iostream>
#include <stdlib.h>
#include <time.h>
#include<fstream>
#include<cstring>
#include<cmath>
using namespace std;

int main(){
	int cnt = 1;
	fstream input;
		
	for(int i = 0; i < 999; i++){
		if(i % 10 == 0)
			srand((unsigned)time(NULL));
		char filename[25] = {"testdata_"};
		char* index = new char [5];
		itoa(cnt, index, 10);
		strcat(filename, index);
		strcat(filename, ".txt");
		input.open(filename, ios::out);
		int n = rand() % 1000 + 1;
		input << n << " ";
		int m = static_cast<int>(sqrt(n) * 0.9 + rand() % static_cast<int>(sqrt(n) * 0.2));
		input << m << " ";
		int w = static_cast<int>(sqrt(n) * 0.7 + rand() % static_cast<int>(sqrt(n) * 0.6));
		input << w << " ";
		int d = rand() % (9*n) + n;
		input << d << "\n";
		for(int j = 0; j < m; j++){
			input << rand() % n << " ";
		}
		input << "\n";
		for(int j = 0; j < m; j++){
			input << rand() % n << " ";
		}
		input << "\n";
		for(int j = 0; j < m; j++){
			input << rand() % n << " ";
		}
		input << "\n";
		for(int j = 0; j < m; j++){
			input << rand() % 1000 << " ";
		}
		input << "\n";
		input << 0 << " " << 0 << " " << n << " " << n;
 		input.close();
 		cout << "(" << cnt << "/999)\n";
 		cnt++;
	}
	return 0; 
}
