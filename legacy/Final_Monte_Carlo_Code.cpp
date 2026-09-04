
#if __has_include(<cstdlib>)
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#if __has_include(<cmath>)
#include <cmath>
#else
#include <math.h>
#endif

#if __has_include(<cstdio>)
#include <cstdio>
#else
#include <stdio.h>
#endif

#if __has_include(<ctime>)
#include <ctime>
#else
#include <time.h>
#endif

#include <iostream>
#include <vector>

double myFunction(long double x);
std::vector<double> monteCarloEstimate(double lowBound, double upBound, int iterations);
std::vector<double> HitOrMiss(double lowBound, double upBound, double WMax, int Nevents);
void CreateEvents(std::vector<double> CosThetas, double ECOM);

using namespace std;

constexpr double PI = 3.14159265358979323846;
double ECM = 80; // This is the Giga electronvolts measured

int main()
{
	srand(static_cast<unsigned int>(time(nullptr)));

	double lowerBound, upperBound;
	int iterations;

	lowerBound = -1;
	upperBound = 1;
	iterations = 100;
	double Nevents = 100;
	std::vector<double> MCResult = monteCarloEstimate(lowerBound, upperBound, iterations);
	std::vector<double> Events = HitOrMiss(lowerBound, upperBound, MCResult[2], static_cast<int>(Nevents));
	// Loop Over Vector Events and Plot values should be distributed over 1+x^2
	printf("Estimate for %.1f -> %.1f is %.2f, (%i iterations)\n",
			lowerBound, upperBound, MCResult[0], iterations);
	cout << " Error = " << MCResult[1] << endl;
	cout << " Maximimum Value  = " << MCResult[2] << endl;
	cout << " Generated " << Events.size() << endl;
	CreateEvents(Events, ECM);
	return 0;
}


double myFunction(long double x)
//Function to integrate
{
    double alpha = 1./137. ;
    
	return pow(alpha/(2*ECM),2)*(1 + x*x);
}

std::vector<double> monteCarloEstimate(double lowBound, double upBound, int iterations)
// Function to execute Monte Carlo integration on predefined function
{
	std::vector<double> Result;
	double totalSum = 0;
	double randNum, functionVal;
	double SumSqr = 0;
	int iter = 0;
	double WMax = -1;
	while (iter < iterations)
	{
		// Select a random number within the limits of integration
		randNum = lowBound + (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (upBound - lowBound); // this is x

		// Sample the function's values
		functionVal = myFunction(randNum) * (upBound - lowBound);
		if (functionVal > WMax) WMax = functionVal;
		// Add the f(x) value to the running sum
		totalSum += functionVal;
		SumSqr += pow(functionVal, 2);
		iter++;
	}

	double estimate = totalSum / iterations;
	double er = sqrt(((SumSqr / iterations) - pow(estimate, 2)) / iterations);
	Result.push_back(estimate);
	Result.push_back(er);
	Result.push_back(WMax);
	return Result;
}

std::vector<double> HitOrMiss(double lowBound, double upBound, double WMax, int Nevents)
// Function to Execute Hit Or Miss
{
	std::vector<double> Events;
	double x1, functionVal;
	int iter = 0;

	while (iter < Nevents) {
		x1 = lowBound + (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (upBound - lowBound); // this is x

		// x1 is the cosine of the angle
		functionVal = myFunction(x1) * (upBound - lowBound);
		double Prob = functionVal / WMax;
		double RandNumb = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		if (RandNumb < Prob) {
			Events.push_back(x1);
			iter++;
		}
	}

	return Events;
}

void CreateEvents(std::vector<double> CosThetas, double ECOM)
{
	int thetaSize = static_cast<int>(CosThetas.size());
	int iter2 = 0;
	while (iter2 < thetaSize) {
		double cosTH = CosThetas[iter2];
		double sinTH = sqrt(1.0 - pow(cosTH, 2));
		double Phi = (2.0 * PI) * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
		cout << " Event " << iter2 << endl;
		cout << " P1 Vector = " << ECM / 2 << " " << 0 << " " << 0 << " " << ECM / 2 << "k" << endl;
		cout << " P2 Vector = " << ECM / 2 << " " << 0 << " " << 0 << " " << -ECM / 2 << "k" << endl;
		cout << " P3 Vector = " << ECM / 2 << " " << (ECM / 2) * sinTH * cos(Phi) << "i" << " " << (ECM / 2) * sinTH * sin(Phi) << "j" << " " << (ECM / 2) * cosTH << "k" << endl;
		cout << " P4 Vector = " << ECM / 2 << " " << -(ECM / 2) * sinTH * cos(Phi) << "i" << " " << -(ECM / 2) * sinTH * sin(Phi) << "j" << " " << -(ECM / 2) * cosTH << "k" << endl;
		cout << " Cosine(Theta)= " << cosTH << endl;
		cout << " Sine(Theta)= " << sinTH << endl;
		cout << " Angle Theta = " << acos(cosTH) * (180.0 / PI) << " " << "Degrees" << endl;
		cout << " The Transverse Momentum= " << (ECM / 2) * sinTH << "GeV" << endl;
		cout << "_________________________________________________________" << endl;
		iter2++;
	}
}


