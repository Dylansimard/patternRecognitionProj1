#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

double ranf(double m);
float generateSamples(int mu, int sigma);
void printToFile(vector<float> &valX, vector<float> &valY, string outputfile);
void bayesClassifier(vector<float> &valX, vector<float> &valY, float muOne, float muTwo, float varianceOne, float varianceTwo);
int main()
{
	vector<float> x1, y1, x2, y2;

	for (int i = 0; i < 50000; ++i) {
		x1.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(4, 1));
		y2.push_back(generateSamples(4, 1));
        
	}
    printToFile(x1, y1, "Sample1.txt");
	printToFile(x2, y2, "Sample2.txt");

}





void printToFile(vector<float> &valX, vector<float> &valY, string outputfile){
    ofstream outfile(outputfile);
        outfile << " X VALUE    Y VALUE" << endl;
    for(int i = 0; i < 50000; ++i){
        outfile << "(" << valX[i] << ", " << valY[i] << ")" << endl;
    }
}

void bayesClassifier(vector<float> &valX, vector<float> &valY, float muOne, float muTwo, float varianceOne, float varianceTwo){
    
}