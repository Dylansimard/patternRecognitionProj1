#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include "boxmuller.cpp"
using namespace std;

void printToFile(vector<pair<float, float>> &truePositives, vector<pair<float, float>> &falseNegatives, string outputfile);

float gaussianDescriminant(float valX, float valY, vector<float> mu, int sigma, float probability);
float gaussianDescriminant3(float valX, float valY, vector<float> mu, vector<vector<float>> sigma, float probability);

float minDistClassifier(float valX, float valY, vector<float> mu);

float calculateDecisionBound(float valX, float valY, vector<float> mu1, vector<float> mu2, float sigma, float probability1, float probability2);

float bhattacharyyaBound(vector<float> mu1, vector<float> mu2, vector<vector<float>> sigma1, vector<vector<float>> sigma2, float beta);
float calculateBhattacharyyaDenominator(vector<vector<float>> sigma1, vector<vector<float>> sigma2);





int main()
{



	vector<float> x1, y1, x2, y2;

	for (int i = 0; i < 100000; ++i) {
		x1.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(4, 1));
		y2.push_back(generateSamples(4, 1));
        
	}
	
	vector<float> muOne = {1.0, 1.0};
	vector<vector<float>> sigmaOne = {{1.0, 0}, {0 , 1.0}};

	vector<float> muTwo = {4.0, 4.0};
	vector<vector<float>> sigmaTwo = {{1.0, 0.0}, {0.0, 1.0}};

	vector<pair<float, float>> truePositives, truePositives2;
	vector<pair<float, float>> falseNegatives, falseNegatives2;


	// for part a -- class 1
	vector<pair<float, float>> decisionBoundx, decisionBoundy;

	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant(x1[i], y1[i], muOne, 1, .5);

		float class2 = gaussianDescriminant(x1[i], y1[i], muTwo, 1, .5);

		//decisionBoundx.push_back(calculateDecisionBound(x1[i], y1[i], muOne, muTwo, 1, .5, .5));

		pair<float, float> temp = make_pair(x1[i], y1[i]);

		if(class1 == class2){
			cout << "success?" << endl;
			decisionBoundx.push_back(temp);
		}
		else if(class1 > class2){
			truePositives.push_back(temp);
		}
		else{
			falseNegatives.push_back(temp);
		}




		

	}

	string fileName = "num1class1PartA.csv";
	//printToFile(truePositives, falseNegatives, fileName);

	

	cout << "---------------------------" << endl;
	cout << "Number 1 Class 1" << endl;
	cout << "Part A" << endl;
	cout << "True positives: " << truePositives.size() << endl;
	cout << "False negatives: " << falseNegatives.size() << endl;

	float bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;

	

	truePositives.clear();
	falseNegatives.clear();


	// for part a -- class 2
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant(x2[i], y2[i], muOne, 1, .5);

		float class2 = gaussianDescriminant(x2[i], y2[i], muTwo, 1, .5);

		//decisionBoundy.push_back(calculateDecisionBound(x2[i], y2[i], muOne, muTwo, 1, .5, .5));

		pair<float, float> temp = make_pair(x2[i], y2[i]);

		if(class1 <= class2){
			truePositives.push_back(temp);
		}
		else{
			falseNegatives.push_back(temp);
		}
	}

	//cout << "test" << endl <<  calculateDecisionBound(x1[0], y1[0], muOne, muTwo, 1, .5, .5) << endl;

	fileName = "num1class2PartA.csv";
	printToFile(truePositives, falseNegatives, fileName);

	cout << "---------------------------" << endl;
	cout << "Number 1 Class 2" << endl;
	cout << "Part A" << endl;
	cout << "True positives: " << truePositives.size() << endl;
	cout << "False negatives: " << falseNegatives.size() << endl;

	cout <<  "Bhattacharyya Bound: " << bhat << endl;

	cout << "---------------------------" << endl;

	ofstream outFile("test.csv");

	for(int i = 0; i < 100000; ++i){
		if(i < truePositives.size()){
			outFile << truePositives[i].first << "," << truePositives[i].second << ",";
		}
		else{
			outFile << "," << ",";
		}

		if(i < falseNegatives.size()){
			outFile << falseNegatives[i].first << "," << falseNegatives[i].second << ",";
		}
		else{
			outFile << "," << ",";
		}

		if(i < truePositives2.size()){
			outFile << truePositives2[i].first << "," << truePositives2[i].second << ",";
		}
		else{
			outFile << "," << ",";
		}

		if(i < falseNegatives2.size()){
			outFile << falseNegatives2[i].first << "," << falseNegatives2[i].second << "," ;
		}
		else{
			outFile << "," << ",";
		}
		
		//outFile << decisionBoundx[i] << "," << decisionBoundy[i] << "\n";

	}

	ofstream outFile2("test2.csv");
	for(int i = 0; i < 100000; ++i){
		//outFile2 << x1[i] << "," << y1[i] << "," << x2[i] << "," << y2[i] << "," << decisionBoundx[i] << "\n";

	}

	// for part b -- class 1	
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant(x1[i], y1[i], muOne, 1, .2);

		float class2 = gaussianDescriminant(x1[i], y1[i], muTwo, 1, .8);

		pair<float, float> temp = make_pair(x1[i], y1[i]);

		if(class1 >= class2){
			truePositives2.push_back(temp);
		}
		else{
			falseNegatives2.push_back(temp);
		}
	}

	fileName = "num1class1PartB.csv";
	printToFile(truePositives2, falseNegatives2, fileName);	

	cout << "Number 1 Class 1" << endl;
	cout << "Part B" << endl;
	cout << "True positives: " << truePositives2.size() << endl;
	cout << "False negatives: " << falseNegatives2.size() << endl;
	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;
	
	cout << "---------------------------" << endl;



	truePositives2.clear();
	falseNegatives2.clear();

	// for part b -- class 2	
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant(x2[i], y2[i], muOne, 1, .2);

		float class2 = gaussianDescriminant(x2[i], y2[i], muTwo, 1, .8);

		pair<float, float> temp = make_pair(x1[i], y2[i]);

		if(class1 <= class2){
			truePositives2.push_back(temp);
		}
		else{
			falseNegatives2.push_back(temp);
		}
	}

	fileName = "num1class2PartB.csv";
	printToFile(truePositives2, falseNegatives2, fileName);	

	cout << "Number 1 Class 2" << endl;
	cout << "Part B" << endl;
	cout << "True positives: " << truePositives2.size() << endl;
	cout << "False negatives: " << falseNegatives2.size() << endl;
	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;
	
	cout << "---------------------------" << endl;

/*
*
*	START PART TWO FOR as1.cpp
*
*	new mu1 = [1, 1]
*	new mu2 = [4, 4]
*	new sigma1 = [1, 0
*				  0, 1]
*	new sigma2 = [4, 0
*				  0, 8]
*
*	REPEAT ALL STEPS FROM PART 1
*
*
*/
	

	muOne = {1.0, 1.0};
	sigmaOne = {{1.0, 0.0}, {0.0 , 1.0}};

	muTwo = {4.0, 4.0};
	sigmaTwo = {{4.0, 0.0}, {0.0, 8.0}};

	
	//Clearing generated data vectors 
	//Clearing truePositive and falseNegative vectors
		x1.clear();
		x2.clear();
		y1.clear();
		y2.clear();
		truePositives.clear();
		falseNegatives.clear();
		truePositives2.clear();
		falseNegatives2.clear();
		
	// Generate new data set for Part 2
	for (int i = 0; i < 100000; ++i) {
		x1.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(4, 4));
		y2.push_back(generateSamples(4, 8));
        
	}

	// for number 2 part a class 1
	for(int i = 0; i < 100000; ++i){
		//cout << "class1" << endl;
		float class1 = gaussianDescriminant3(x1[i], y1[i], muOne, sigmaOne, .5);
		//cout << "total = " << class1 << endl;

		//cout << "class2" << endl;
		float class2 = gaussianDescriminant3(x1[i], y1[i], muTwo, sigmaTwo, .5);
		//cout << "total = " << class2 << endl;

		pair<float, float> temp = make_pair(x1[i], y1[i]);

		if(class1 >= class2){
			truePositives.push_back(temp);
		}
		else{
			falseNegatives.push_back(temp);
		}
	}
	
	cout << "Number 2 Class 1" << endl;
	cout << "Part A" << endl;
	cout << "True positives: " << truePositives.size() << endl;
	cout << "False negatives: " << falseNegatives.size() << endl;

	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;

	cout << "---------------------------" << endl;

	fileName = "num2class1PartA.csv";
	printToFile(truePositives, falseNegatives, fileName);	

	truePositives.clear();
	falseNegatives.clear();


	// for number 2 part a class 2
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant3(x2[i], y2[i], muOne, sigmaOne, .5);

		float class2 = gaussianDescriminant3(x2[i], y2[i], muTwo, sigmaTwo, .5);
		pair<float, float> temp = make_pair(x2[i], y2[i]);

		if(class1 <= class2){
			truePositives.push_back(temp);
		}
		else{
			falseNegatives.push_back(temp);
		}
	}
	
	cout << "Number 2 Class 2" << endl;
	cout << "Part A" << endl;
	cout << "True positives: " << truePositives.size() << endl;
	cout << "False negatives: " << falseNegatives.size() << endl;

	cout <<  "Bhattacharyya Bound: " << bhat << endl;

	cout << "---------------------------" << endl;

	fileName = "num2class2PartA.csv";
	printToFile(truePositives, falseNegatives, fileName);	



	// for number 2 part b class 1
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant3(x1[i], y1[i], muOne, sigmaOne, .2);

		float class2 = gaussianDescriminant3(x1[i], y1[i], muTwo, sigmaTwo, .8);

		pair<float, float> temp = make_pair(x1[i], y1[i]);

		if(class1 >= class2){
			truePositives2.push_back(temp);
		}
		else{
			falseNegatives2.push_back(temp);
		}
	}

	cout << "Number 2 Class 1" << endl;
	cout << "Part B" << endl;
	cout << "True positives: " << truePositives2.size() << endl;
	cout << "False negatives: " << falseNegatives2.size() << endl;
	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;
	
	cout << "---------------------------" << endl;

	fileName = "num2class1PartB.csv";
	printToFile(truePositives2, falseNegatives2, fileName);

	truePositives2.clear();
	falseNegatives2.clear();

	// for number 2 part b class 2
	for(int i = 0; i < 100000; ++i){
		float class1 = gaussianDescriminant3(x2[i], y2[i], muOne, sigmaOne, .2);

		float class2 = gaussianDescriminant3(x2[i], y2[i], muTwo, sigmaTwo, .8);

		pair<float, float> temp = make_pair(x2[i], y2[i]);

		

		if(class1 <= class2){
			truePositives2.push_back(temp);
		}
		else{
			falseNegatives2.push_back(temp);
		}
	}

	cout << "Number 2 Class 2" << endl;
	cout << "Part B" << endl;
	cout << "True positives: " << truePositives2.size() << endl;
	cout << "False negatives: " << falseNegatives2.size() << endl;
	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;
	
	cout << "---------------------------" << endl;

	fileName = "num2class2PartB.csv";
	printToFile(truePositives2, falseNegatives2, fileName);


	truePositives2.clear();
	falseNegatives2.clear();

	// for number 3
	for(int i = 0; i < 100000; ++i){
		float class1 = minDistClassifier(x2[i], y2[i], muOne);

		float class2 = minDistClassifier(x2[i], y2[i], muTwo);

		pair<float, float> temp = make_pair(x2[i], y2[i]);

		

		if(class1 <= class2){
			truePositives2.push_back(temp);
		}
		else{
			falseNegatives2.push_back(temp);
		}
	}

	cout << "Number 3" << endl;
	cout << "True positives: " << truePositives2.size() << endl;
	cout << "False negatives: " << falseNegatives2.size() << endl;
	bhat = bhattacharyyaBound(muOne, muTwo, sigmaOne, sigmaTwo, 0.5);
	cout <<  "Bhattacharyya Bound: " << bhat << endl;
	
	cout << "---------------------------" << endl;


	
}


/**
 * 	
 * 	output to csv file
 * 	, delimits rows, \n signifies a new line
 * 	truePositives are always bigger, so loops through that entire vector
 * 	only writes false negatives for how many there are, then only writes the truepositives.
 * 	will have to do multiple times -- once for each class -- I think 4 total
 * 
 * 
 **/


void printToFile(vector<pair<float, float>> &truePositives, vector<pair<float, float>> &falseNegatives, string outputfile){
    ofstream outfile(outputfile);
    outfile << "TP X VALUE, TP Y VALUE, FN X VALUE, FN Y VALUE\n" << endl;

	for(int i = 0; i < truePositives.size(); ++i){

		if(i <= falseNegatives.size()){
			outfile << truePositives[i].first << "," << truePositives[i].second << "," << falseNegatives[i].first << "," << falseNegatives[i].second << "\n";
		}
		else{
			outfile << truePositives[i].first << "," << truePositives[i].second << "\n";
		}

	}
    
}


/**
 * gi(x) = - ((x-mu)^2 / 2 * (little Sigma)^2) + ln p(wi)
 * 
 * 		
 * 
 * 
 * 				
 **/

float gaussianDescriminant(float valX, float valY, vector<float> mu, int sigma, float probability){

	// (x - mu) ^ 2
	vector<float> temp;
	temp.push_back(abs(valX - mu[0]));
	temp.push_back(abs(valY - mu[1]));

	float top = (temp[0] * temp[0]) + (temp[1] * temp[1]);



	// 2 * (little sigma) ^ 2
	float bottom = pow(sigma, 2);
	bottom *= 2;

	float left = top / bottom;

	left *= -1;

	float probLog = log(probability);

	return left + probLog;

}

/**				left			 middle			right
 * 	g_i(x) = (x^t * W_i * x) + (w_i^t * x) + w_i0
 * 
 * 	W_i = -.5 * sigma_i ^ -1
 * 
 * 	w_i^t = (sigma_i ^ -1) * mu_i
 * 
 * 	w_i0 = (-.5 * (mu_i ^ t) * (sigma_i ^ -1) * mu_i) - (-.5 * ln|sigma_i|) + ln(P(w_i))
 * 															
 * 
 * */
float gaussianDescriminant3(float valX, float valY, vector<float> mu, vector<vector<float>> sigma, float probability){


	// calculating first term x^t * W_i * x
	// W_i = -1/2 * sigma^-1
	// final value storeed in left
	float determinant = (sigma[0][0] * sigma[1][1] - sigma[1][0] * sigma[0][1]);


	vector<vector<float>> inverseSigma{
		{((1/determinant) * sigma[1][1]), (-1 * (1/determinant) * sigma[1][0])},
		{(-1 * (1/determinant) * sigma[0][1]), ((1/determinant) * sigma[0][0])}
	};

	vector<vector<float>> W_i = inverseSigma;
	W_i[0][0] *= -.5;
	W_i[0][1] *= -.5;
	W_i[1][0] *= -.5;
	W_i[1][1] *= -.5;

	vector<float> temp = {(valX * W_i[0][0]) + (valY * W_i[1][0]), (valX * W_i[0][1]) + (valY * W_i[1][1])};

	float left = temp[0] * valX + temp[1] * valY;

	//w_i^t * x
	//w_i = inverseSigma * mu_i

	vector<float> w_i {(inverseSigma[0][0] * mu[0]) + (inverseSigma[1][0] * mu[1]), (inverseSigma[0][1] * mu[0]) + (inverseSigma[1][1] * mu[1])};

	float middle = w_i[0] * valX + w_i[1] * valY;

	// calculate w_i0
	// w_i0 = -.5 * mu^t * inverseSigma * mu - 1/2 ln|sigma| + ln prob

	vector<float> next1 {((float)-.5 * mu[0]), ((float)-.5 * mu[1])};

	vector<float> next2 {(next1[0] * inverseSigma[0][0]) + (next1[1] * inverseSigma[1][0]), (next1[0] * inverseSigma[0][1]) + (next1[0] * inverseSigma[1][1])};

	float next3 = next2[0] * mu[0] + next2[1] * mu[1];

	float right = next3 - (.5 * log(determinant)) + log(probability);

	return left + middle + right;

}

// g_i(x) = - ||x - u_i|| ^ 2
// ||x - u_i|| ^ 2 = 
//						(x-u_i)^t * (x-U_i)
//

float minDistClassifier(float valX, float valY, vector<float> mu){

	vector<float> temp1 {(valX - mu[0]), (valY - mu[1])};
	vector<float> temp2 = temp1;

	return -((temp1[0] * temp2[0]) + (temp1[1] * temp2[1]));

}


/**
 * 	w^t * (x - x_0)
 * 
 * 	w = mu_i - mu_j
 * 
 * 	x_0 = .5 * (mu_i + mu_j) - ( (sigma^2) / || mu_i - mu_j|| ^ 2) * ln ( P(w_i) / P(w_j) ) * ( mu_i - mu_j )
 * 			 first						second								third					fourth
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 **/

float calculateDecisionBound(float valX, float valY, vector<float> mu1, vector<float> mu2, float sigma, float probability1, float probability2){

	// w^t
	vector<float> w;
	w.push_back(mu1[0] - mu2[0]);
	w.push_back(mu1[1] - mu2[1]);


	// x_0

	vector<float> first;
	first.push_back(.5 * (mu1[0] + mu2[0]));
	first.push_back(.5 * (mu1[1] + mu2[1]));



	float second = sigma * sigma;



	// (x - mu) ^ 2
	vector<float> temp;
	temp.push_back(abs(mu1[0] - mu2[0]));
	temp.push_back(abs(mu1[1] - mu2[1]));

	float next = (temp[0] * temp[0]) + (temp[1] * temp[1]);

	second = second / next;

	

	float third = log(probability1 / probability2);


	vector<float> fourth;
	fourth.push_back(mu1[0] - mu2[0]);
	fourth.push_back(mu1[1] - mu2[1]);

	// second * third
	float secondThird = second * third;
	
	// times fourth as scalar
	fourth[0] *= secondThird;
	fourth[1] *= secondThird;

	//first - the rest
	first[0] -= fourth[0];
	first[1] -= fourth[1];

	// (x - x_0)
	vector<float> right;
	right.push_back(valX - first[0]);
	right.push_back(valY - first[1]);


	// w * (x - x_0)

	float result = (w[0] * right[0]) + (w[1] * right[1]);

	return result;
	
}


//
// Calculate Bhattacharyya Bound EQ
// Value for B = 0.5
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
float bhattacharyyaBound(vector<float> mu1, vector<float> mu2, vector<vector<float>> sigma1, vector<vector<float>> sigma2, float beta){
	float first = (beta * (1 - beta))/2;
	
	float detClass1 = sigma1[0][0] * sigma1[1][1] - sigma1[0][1] * sigma1[1][0];
	float detClass2 = sigma2[0][0] * sigma2[1][1] - sigma2[0][1] * sigma2[1][0];
	
	vector<float> second = {mu1[0] - mu2[0], mu1[1] - mu2[1]};

	
	vector<vector<float>> third = {{(sigma1[0][0] + sigma2[0][0]) * (beta), 0}, {0, (sigma1[1][1] + sigma2[1][1]) * (beta)}};
	float detThird = third[0][0] * third[1][1] - third[0][1] * third[1][0];
	float inverseDetThird = 1/detThird;

	vector<vector<float>> inverseThird = {{inverseDetThird * third[1][1], 0}, {0, inverseDetThird * third[0][0]}};

	vector<float> fourth = {mu1[0] - mu2[0], mu1[1] - mu2[1]};

	float numerator = detThird;

	float denominator = (pow(detClass1, 1 - beta) * pow(detClass2, beta));

	float fifth = 0.5 * log(numerator / denominator);

	//
	//	first: float value first = B(1 - B)
	//     	                       --------
	//     	                          2
	//	vector second: Using class1 and class2 respectively mu1 and mu2
	//		 second = [mu1_x - mu2_x, mu1_y - mu2_y]
	//
	//	firstSecond = first * second = [first * (mu1_x - mu2_x), first * (mu1_y - mu2_y)]
	//
	vector<float> firstSecond = {first * second[0], first * second[1]};

	vector<float> secondThird = {inverseThird[0][0] * firstSecond[0], inverseThird[1][1] * firstSecond[1]};

	float thirdFourth = (secondThird[0] * fourth[0]) + (secondThird[1] * fourth[1]);
	
	float fourthFifth = thirdFourth + fifth;

	float result = exp(-fourthFifth) * .5;

	return result;


}