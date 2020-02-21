#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include "boxmuller.cpp"
using namespace std;


float generateSamples(int mu, int sigma);
void printToFile(vector<float> &valX, vector<float> &valY, string outputfile);
bool bayesClassifier(vector<float> &valX, vector<float> &valY, float muOne, float muTwo, float varianceOne, float varianceTwo);

float gaussianDescriminant(float valX, float valY, vector<float> mu, vector<vector<float>> sigma);
float calculateDenominator(vector<vector<float>> sigma);

float calculateExponent(float valX, float valY, vector<float> mu, vector<vector<float>> sigma);

int main()
{
	vector<float> x1, y1, x2, y2;

	for (int i = 0; i < 50000; ++i) {
		x1.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(4, 1));
		y2.push_back(generateSamples(4, 1));
        
	}
    //printToFile(x1, y1, "Sample1.txt");
	//printToFile(x2, y2, "Sample2.txt");

	
	vector<float> muOne = {1.0, 1.0};
	vector<vector<float>> sigmaOne = {{1.0, 0.0}, {0.0 , 1.0}};

	vector<float> muTwo = {4.0, 4.0};
	vector<vector<float>> sigmaTwo = {{1.0, 0.0}, {0.0, 1.0}};

	vector<pair<int, int>> truePositives;
	vector<pair<int, int>> falseNegatives;


	for(int i = 0; i < 50000; ++i){
		float class1 = gaussianDescriminant(x1[i], y1[i], muOne, sigmaOne);
		class1 *= .5;

		float class2 = gaussianDescriminant(x1[i], y1[i], muTwo, sigmaTwo);
		class2 *= .5;

		pair<int, int> temp = make_pair(x1[i], x2[i]);

		if(class1 >= class2){
			truePositives.push_back(temp);
		}
		else{
			falseNegatives.push_back(temp);
		}
	}

	cout << "True positives: " << truePositives.size() << endl;
	cout << "False negatives: " << falseNegatives.size() << endl;
}





void printToFile(vector<float> &valX, vector<float> &valY, string outputfile){
    ofstream outfile(outputfile);
        outfile << " X VALUE    Y VALUE" << endl;
    for(int i = 0; i < 50000; ++i){
        outfile << "(" << valX[i] << ", " << valY[i] << ")" << endl;
    }
}




bool bayesClassifier(vector<float> &valX, vector<float> &valY, float muOne, float muTwo, float varianceOne, float varianceTwo){
    
}

/**
 * Thinking run this function for each input vector x = [x_i, y_i]
 * With respective mu vector -- 
 * 
 * 								[mu_x, mu_y]
 * 
 * and respectove sigma matrix -- 
 * 
 * 								[sigma_x, 0]
 * 								[0, sigma_y]
 * 
 * Thinking splitting the descriminant into 2 parts --
 * 
 * 		The denominator -- 
 * 
 * 								((2 * pi) ^ d/2) * |sigma| ^ (1/2)
 * 
 * 		The exponent --
 * 
 * 								-1/2 * ((x - mu) ^ t) * (sigma ^ -1) * (x - mu)
 * 				
 **/

float gaussianDescriminant(float valX, float valY, vector<float> mu, vector<vector<float>> sigma){

	float denominator = calculateDenominator(sigma);

	denominator = 1/denominator;

	float exponent = calculateExponent(valX, valY, mu, sigma);

	float result = pow(denominator, exponent);

	return result;

}


// Function to calculate: 		  "left"           "right"
//							((2 * pi) ^ d/2) * |sigma| ^ (1/2)
//
//	|sigma| is the determinant of sigma
//
// Assuming d = dimensionality, which == 2
// 
// returns the denominator value -- not 1/denominator -- still have to do 1/denominator in prev func.

float calculateDenominator(vector<vector<float>> sigma){
	
	float left = 2 * M_PI;

	// sigma will always be 2x2 in this assignment
	float determinant = sigma[0][0] * sigma[1][1] - sigma[0][1] * sigma[1][0];

	float right = pow(determinant, .5);

	float result = left * right;

	return result;

}

// Function to calculate: 				"left"			   "middle"      "right"
//								-1/2 * ((x - mu) ^ t) * (sigma ^ -1) * (x - mu)
float calculateExponent(float valX, float valY, vector<float> mu, vector<vector<float>> sigma){
	
	// (x - mu) -- vector
	vector<float> left = {valX - mu[0], valY - mu[1]};

	// above * 1/2 as a scalar
	left[0] *= .5;
	left[1] *=.5;


	// following for "middle" sigma ^ -1 -- have to inverse matrix
	// done by the following:

	// assume matrix = 
	//					[a, b]
	//					[c, d]
	// the inverse is:
	//								
	//				(1 / determinent) * [d, -b]
	//									[-c, a]


	float determinant = sigma[0][0] * sigma[1][1] - sigma[0][1] * sigma[1][0];
	determinant = 1/determinant;
	vector<vector<float>> middle = {
		{(determinant * sigma[1][1]), (determinant * -sigma[0][1])},
		{(determinant * -sigma[1][0]), (determinant * sigma[0][0])}
	};


	vector<float> right = {(valX - mu[0]), (valY - mu[1])};

	// now left ^ t * middle
	// the left side is technically [x
	//								 y]
	// so it is already transposed as it is stored
	// makes the equation look like:
	//		
	//			 "left"				"middle"
	//		[left_0, left_1] * [middle_0, middle_2
	//						    middle_1, middle_3]

	// will result in the following:
	//
	//			[leftMiddle_0, leftMiddle_1]

	vector<float> leftMiddle = {
		(left[0] * middle[0][0] + left[1] * middle[1][0]), (left[0] * middle[0][1] + left[1] * middle[1][1])
	};


	// now left * middle is calculated -- calculate leftmiddle * right --> result
	// leftMiddle stored as:
	//
	//			[leftMiddle_0, leftMiddle_1]
	//
	// right is stored as:
	//
	//						[right_0
	//						 right_1]
	//
	// so resulting value should be a singular value --> which is the return value

	float result = ((leftMiddle[0] * right[0]) + (leftMiddle[1] * right[1]));

	return result;

}