#include <iostream>
#include <vector>
using namespace std;


double ranf(double m) {
	return (m * rand()) / (double)RAND_MAX;
}

float generateSamples(int mu, int sigma) {

	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last) {
		y1 = y2;
		use_last = 0;
	}

	else {
		do {
			x1 = 2.0 * ranf(1) - 1.0;
			x2 = 2.0 * ranf(1) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return (mu + y1 * sigma);

}




int main()
{

	vector<float> x1, y1, x2, y2;

	for (int i = 0; i < 50000; ++i) {
		x1.push_back(generateSamples(1, 1));
		x2.push_back(generateSamples(1, 1));
		y1.push_back(generateSamples(4, 1));
		y2.push_back(generateSamples(4, 1));
	}

}

