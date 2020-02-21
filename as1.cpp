#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include "boxmuller.cpp"
using namespace std;

void printToFile(vector<float> &valX, vector<float> &valY, string outputfile);
float gaussianDescriminant(float valX, float valY, vector<float> mu, vector<vector<float>> sigma);
float calculateDenominator(vector<vector<float>> sigma);
float calculateExponent(float valX, float valY, vector<float> mu, vector<vector<float>> sigma);
float bhattacharyyaBound(vector<float> mu1, vector<float> mu2, vector<vector<float>> sigma1, vector<vector<float>> sigma2);
float calculateBhattacharyyaDenominator(vector<vector<float>> sigma1, vector<vector<float>> sigma2);


int main()
{
	vector<float> x1, y1, x2, y2;

	for (int i = 0; i < 50000; ++i) {
		x1.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(4, 1));
		y2.push_back(generateSamples(4, 1));
        
	}
	
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

	float bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo);
	cout <<  "BhattacharyyaBound: " << bhat << endl;
}





void printToFile(vector<float> &valX, vector<float> &valY, string outputfile){
    ofstream outfile(outputfile);
        outfile << " X VALUE    Y VALUE" << endl;
    for(int i = 0; i < 50000; ++i){
        outfile << "(" << valX[i] << ", " << valY[i] << ")" << endl;
    }
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



// Calculate denominator in last portion of Bhattacharyya EQ
//
//	1/2 * ln [(1-B)*Sigma1 + B*Sigma2]
//           ------------------------
//   ======> Sigma1^(1-B) * Sigma2^(B) <=======
//				 ^^^          ^^^
//               det          det
float calculateBhattacharyyaDenominator(vector<vector<float>> sigmaClass1, vector<vector<float>> sigmaClass2){
	float B = 0.5;
	float detClass1 = sigmaClass1[0][0] * sigmaClass1[1][1] - sigmaClass1[0][1] * sigmaClass1[1][0];

	float detClass2 = sigmaClass2[0][0] * sigmaClass2[1][1] - sigmaClass2[0][1] * sigmaClass2[1][0];

	return (pow(detClass1, 0.5) * pow(detClass2, 0.5));
}

//
// Calculate Bhattacharyya Bound EQ
//
// First: B(1 - B)
//        --------
//            2
//
// Second: (mu1 - mu2)^t = (mu1[0] - mu2[0])
//
// Third: [(1-B) * Sigma1 + (B * Sigma2)]^-1
//
// Fourth: (mu1 - mu2)
//
// Fifth: 1/2 * ln [(1 - B) * Sigma1 + (B * Sigma2)]
//                 ---------------------------------
//                     Sigma1^(1-B) * Simga2^B
//
// Combine All Parts: First * Second * Third * Fourth + Fifth
float bhattacharyyaBound(vector<float> mu1, vector<float> mu2, vector<vector<float>> sigma1, vector<vector<float>> sigma2){
	float first = (.5 * .5)/2;
	
	vector<float> second = {mu1[0] - mu2[0], mu1[1] - mu2[1]};

	float detClass1 = sigma1[0][0] * sigma1[1][1] - sigma1[0][1] * sigma1[1][0];
	float detClass2 = sigma2[0][0] * sigma2[1][1] - sigma2[0][1] * sigma2[1][0];
	//cout << detClass1 << " "  << detClass2 << endl;
	
	float third = (.5 * detClass1) + (1/(detClass2 * 0.5));
	float thirdInverse = 1/third;
	//cout << thirdInverse << endl;

	vector<float> fourth = {mu1[0] - mu2[0], mu1[1] - mu2[0]};
	//cout << fourth[0] << fourth[1] <<  endl;

	float numerator = (0.5 * detClass1) + (0.5 * detClass2);
	//cout << numerator << endl;
	float insideLN = numerator/calculateBhattacharyyaDenominator(sigma1, sigma2);
	//cout << insideLN << endl;
	float fifth = 0.5 * log(insideLN);
	
	float firstSecond = first * second[0] + first * second[1];

	float secondThird = thirdInverse * firstSecond + thirdInverse * firstSecond;

	float thirdFourth = {secondThird * fourth[0] + secondThird * fourth[1]};
	
	float fourthFifth = thirdFourth + fifth;

	float result = exp(-fourthFifth);
	
	return result;


}