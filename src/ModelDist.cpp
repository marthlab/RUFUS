/*By ANDREW FARRELL
 * ModelDist.cpp
 * ----------------------------------------------------------
 * Takes raw k-mer frequency histogram produced by Jellyfish,
 * and fits models of the underlying copy number distribution
 * to the raw k-mer histogram for each sample
 * ----------------------------------------------------------
 */

#include <bitset>
#include <fstream>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "Util.h"

using namespace std;

double pi = 3.14159L;
bool Diploid = true;

double norm(double x, double mu, double sigma, double skew, double p) {
	if (x < mu) {
		sigma = sigma + pow((mu - x) * skew, p);
	}
	return (1 / (sqrt(2 * pi * pow(sigma, 2)))) *
		exp(-(pow(x - mu, 2) / (2 * pow(sigma, 2))));	// normal distribution
}

double calcError(vector<long>& histo, vector<double>& ModelSums, double Reads,
								 double ReadLength, int HS, double SC,
								 vector<double>& ErrorModel, vector<double>& ErrorDist) {
	double eH = 10;
	double ErrorTotal = 0;
	double Pe = histo[1] / ((double)Reads * (double)ReadLength * (double)eH);
	double TotalMode = 0;

	for (int i = 2; i < histo.size(); i++) {
		TotalMode += histo[i] * i;
	}

	for (double i = 1; i < histo.size(); i++) {
		double sum = 0;

		for (double d = 2; d < histo.size(); d++) {
			sum += (pow(Pe, (double)i) * (((d - i) * (1 - Pe)) / (d - (i - 1)))) *
	(Reads * (((double)histo[i] * i) / TotalMode)) * eH * ReadLength;
		}

		if (sum > 1e-200) {
			ErrorModel.push_back(sum);
			ErrorTotal += sum;
		} else
			ErrorModel.push_back(0);
	}
	for (double i = 1; i < histo.size(); i++) {
		ErrorDist.push_back(ErrorModel[i] / ErrorTotal);
	}
}

//TODO: refactor variables: i, , many more
double testModelLog(double SC, double stdev, double factor, double skew,
										double power, vector<long>& histo, int Inflection,
										int MaxCopy, double Ybar, double& Rsq) {
	vector<vector<double> > dist;
	vector<vector<double> > prob;
	vector<double> placeHolder;
	dist.push_back(placeHolder);
	prob.push_back(placeHolder);
	int i;
	int j;

	for (i = 1; i < histo.size(); i++) {
		vector<double> d;
		d.push_back(0);
		if (Diploid) {
			d.push_back(
			norm(i, SC / 2,
					 stdev * (1 - ((1 - (stdev / (stdev + ((1) * factor)))) / 2)),
					 skew, power));
		}

		for (long j = 1; j < histo.size() / SC; j++) {
			d.push_back(norm(i, SC * j, stdev + ((j - 1) * factor), skew, power));
		}
		dist.push_back(d);
	}

	for (long j = 1; j < histo.size() / SC; j++) {
		double sum = 0.0;

		for (i = 1; i < histo.size(); i++) {
			sum += dist[i][j];
		}

		for (i = 1; i < histo.size(); i++) {
			dist[i][j] = dist[i][j] / sum;
		}
	}

	for (i = 1; i < histo.size(); i++) {
		vector<double> p;
		double total = 0;
		p.push_back(0);

		for (long j = 1; j < histo.size() / SC; j++) {
			total += dist[i][j];
		}

		for (long j = 1; j < dist[i].size(); j++) {
			p.push_back(dist[i][j] / total);
		}
		prob.push_back(p);
	}

	vector<double> RC;
	RC.push_back(0);
	if (Diploid) {
		double tSC = histo[SC] / dist[SC][2];
		double hetKmers =
			(histo[SC / 2] - (dist[SC / 2][2] * tSC)) / dist[SC / 2][1];

		if (hetKmers > 0) {
			RC.push_back(hetKmers);
		} else {
			RC.push_back(0);
		}
		RC.push_back(histo[SC] / dist[SC][2]);
		
	} else {
		RC.push_back(histo[SC] / dist[SC][1]);
	}

	int max = histo.size() / SC;
	if (Diploid) {
		max = max - 1;
	}
	
	for (long a = 2; a < histo.size() / SC; a++) {
		double Reads = 0;
		if (Diploid) {
			Reads =
	((double)histo[SC * a] / dist[SC * a][a + 1] * prob[SC * a][a + 1]);
		} else {
			Reads = ((double)histo[SC * a] / dist[SC * a][a] * prob[SC * a][a]);
		}
		RC.push_back(Reads);
	}

	vector<vector<double> > model;
	vector<double> ModelSums;
	ModelSums.push_back(0);

	for (long i = 1; i < histo.size(); i++) {
		vector<double> m;
		double sum = 0;
		m.push_back(0);

		for (long j = 1; j < histo.size() / SC; j++) {
			m.push_back(dist[i][j] * RC[j]);
			sum += m[j];
		}
		model.push_back(m);
		ModelSums.push_back(sum);
	}

	double SSresLog = 0;
	double SSres = 0;
	double SStot = 0;

	for (double i = Inflection; i < SC * MaxCopy; i++) {
		SSresLog += (pow(log((double)histo[i]) - log((double)ModelSums[i]), 2));
		SSres += (pow((double)histo[i] - ModelSums[i], 2));
		SStot += (pow((double)ModelSums[i] - Ybar, 2));
	}
	Rsq = SSres / SStot;
	return SSresLog;
}

double testModel(double SC, double stdev, double factor, double skew,
								 double power, vector<long>& histo, int Inflection, int MaxCopy,
								 double Ybar, double& Rsq) {
	vector<vector<double> > dist;
	vector<vector<double> > prob;
	vector<double> placeHolder;
	dist.push_back(placeHolder);
	prob.push_back(placeHolder);
	int i;
	int j;

	for (i = 1; i < histo.size(); i++) {
		vector<double> d;
		d.push_back(0);
		if (Diploid) {
			d.push_back(
			norm(i, SC / 2,
					 stdev * (1 - ((1 - (stdev / (stdev + ((1) * factor)))) / 2)),
					 skew, power));
		}

		for (long j = 1; j < histo.size() / SC; j++) {
			d.push_back(norm(i, SC * j, stdev + ((j - 1) * factor), skew, power));
		}
		dist.push_back(d);
	}

	for (long j = 1; j < histo.size() / SC; j++) {
		double sum = 0.0;

		for (i = 1; i < histo.size(); i++) {
			sum += dist[i][j];
		}

		for (i = 1; i < histo.size(); i++) {
			dist[i][j] = dist[i][j] / sum;
		}
	}

	for (i = 1; i < histo.size(); i++) {
		vector<double> p;
		double total = 0;
		p.push_back(0);
		for (long j = 1; j < histo.size() / SC; j++) {
			total += dist[i][j];
		}

		for (long j = 1; j < dist[i].size(); j++) {
			p.push_back(dist[i][j] / total);
		}
		prob.push_back(p);
	}

	vector<double> RC;
	RC.push_back(0);

	if (Diploid) {
		double tSC = histo[SC] / dist[SC][2];
		double hetKmers =
			(histo[SC / 2] - (dist[SC / 2][2] * tSC)) / dist[SC / 2][1];

		if (hetKmers > 0) {
			RC.push_back(hetKmers);
		} else {
			RC.push_back(0);
		}
		RC.push_back(histo[SC] / dist[SC][2]);

	} else {
		RC.push_back(histo[SC] / dist[SC][1]);
	}
	int max = histo.size() / SC;

	if (Diploid) {
		max = max - 1;
	}

	for (long a = 2; a < histo.size() / SC; a++) {
		double Reads = 0;
		if (Diploid) {
			Reads =
	((double)histo[SC * a] / dist[SC * a][a + 1] * prob[SC * a][a + 1]);
		} else {
			Reads = ((double)histo[SC * a] / dist[SC * a][a] * prob[SC * a][a]);
		}
		RC.push_back(Reads);
	}

	vector<vector<double> > model;
	vector<double> ModelSums;
	ModelSums.push_back(0);

	for (long i = 1; i < histo.size(); i++) {
		vector<double> m;
		double sum = 0;
		m.push_back(0);

		for (long j = 1; j < histo.size() / SC; j++) {
			m.push_back(dist[i][j] * RC[j]);
			sum += m[j];
		}

		model.push_back(m);
		ModelSums.push_back(sum);
	}
	double SSres = 0;
	double SStot = 0;

	for (double i = Inflection; i < SC * MaxCopy; i++) {
		SSres += (pow((double)histo[i] - ModelSums[i], 2));
		SStot += (pow((double)ModelSums[i] - Ybar, 2));
	}

	Rsq = SSres / SStot;
	return SSres;
}

void correctErrorModel(vector<double>& Error, double& Total, int max) {
	double bestP = 0;
	bool flag = false;

	for (double p = 7; p >= 1; p += -.01) {

		for (int i = 1; i < max; i++) {
			if ((Error[i] - ((1 / (pow(i, p))) * Error[1])) < 0) {
				flag = true;
				break;
			}
		}

		if (flag) {
			break;
		}
		bestP = p;
	}
	Total = 0;

	for (int i = 1; i < Error.size(); i++) {
		Error[i] = (1 / (pow(i, bestP))) * Error[1];
		Total += Error[i];
	}
	cout << "best error is 1/x^" << bestP << endl;
}

void FitErrorModel(vector<double>& Error, double& Total, int max) {
	double LastSSQ = 0;
	float bestP = 0;

	for (int i = 1; i < max; i++) {
		LastSSQ += pow(log(Error[i]) - log((1 / (pow(i, 100))) * Error[1]), 2);
	}

	for (float p = 7; p > .1; p += -.001) {
		double SSQ = 0;
		for (int i = 1; i < max; i++) {
			SSQ += pow(log(Error[i]) - log((1 / (pow(i, p))) * Error[1]), 2);
		}

		if (SSQ < LastSSQ) {
			LastSSQ = SSQ;
			bestP = p;
		}
	}
	Total = 0;

	for (int i = 1; i < Error.size(); i++) {
		Error[i] = (1 / (pow(i, bestP))) * Error[1];
		Total += Error[i];
	}
	cout << "best error is 1/x^" << bestP << endl;
}

void correctContSubtract(vector<double>& Cont, double& Total, int inflect) {
	double last = 0;
	for (int i = 1; i < Cont.size(); i++) {
		if (i > inflect && Cont[i] > last) {
			Cont[i] = last;
			Total += last;
		} else {
			Cont[i] = Cont[i];
			last = Cont[i];
			Total += Cont[i];
		}
	}
}

int main(int argc, char* argv[]) {
	 cout << "Call is histoFile HS ReadLength Threads" << endl;
	ifstream HistoFile;
	HistoFile.open(argv[1]);

	if (HistoFile.is_open()) {
		cout << "Parent File open - " << argv[1] << endl;
	} else {
		cout << "Error, HistoFile could not be opened";
		return 0;
	}

	string path;


	ofstream ModelFile;
	path = argv[1];
	path += ".7.7.model";
	ModelFile.open(path.c_str());

	if (ModelFile.is_open()) {
	} else {
		cout << "Error, Model file could not be opened";
		return 0;
	}

	ofstream DistFile;
	path = argv[1];
	path += ".7.7.dist";
	DistFile.open(path.c_str());
	if (DistFile.is_open()) {
	} else {
		cout << "Error, Model file could not be opened";
		return 0;
	}
 	
	double SC = 1;
	double SCvalue = -1;
	double stdi = -1;
	double stdvalue = -1;
	double value = -1;
	double HistoSum = 0;
	int HS = -1;
	string par = argv[2];
	HS = atoi(par.c_str());
	int Inflection = -1;
	string line;
	line = argv[3];
	int ReadLength = atoi(line.c_str());
	int NumberOfReads = -1;
	line = argv[4];
	int Threads = atoi(line.c_str());
	vector<string> temp;
	vector<long> histo;
	histo.push_back(0);
	getline(HistoFile, line);	// burn the 0 0 line
	temp = Util::Split(line, '\t');
	cout << "first line = " << temp[0] << " - " << temp[1] << endl;
	int count = 0;

	while (atoi(temp[1].c_str()) == 0 or atoi(temp[0].c_str()) == 0) {
		cout << "getting another " << endl;
		getline(HistoFile, line);
		temp = Util::Split(line, '\t');
		cout << "got " << temp[0] << " - " << temp[1] << endl;
		count++;
		if (count > 10) {
			cout << "ERROR there are no kmers in this file" << endl;
			return 1;
		}
	}

	cout << "going with " << temp[0] << " - " << temp[1] << endl;
	value = atol(temp[1].c_str());
	histo.push_back(value);
	long last = value;
	bool PastInflection = false;
	long i = 1;
	long total = 0;
	long TotalKmers = 0;

	while (getline(HistoFile, line)) {
		i++;
		temp = Util::Split(line, '\t');
		value = atol(temp[1].c_str());
		histo.push_back(value);
		total += value;
		TotalKmers += value * atoi(temp[0].c_str());
		HistoSum += value;

		if (value - last > 0 && PastInflection == false) {
			Inflection = i - 1;
			PastInflection = true;
		}

		if (PastInflection) {
			if (SCvalue < value) {
				SCvalue = value;
				SC = i;
			}
		}
		last = histo[i];
	}

	NumberOfReads = TotalKmers / (ReadLength - HS + 1);
	cout << "Number of reads = " << NumberOfReads << endl;
	double Ybar = (double)total / (double)i;

	for (i = 0; i < 10; i++) {
		cout << "I = " << i << " \t " << histo[i] << endl;
	}
	cout << "SC = " << SC << " vlaue = " << SCvalue << endl;
	int rawSC = SC; 
	stdvalue = SCvalue * exp(-.5);

	for (i = SC; i < histo.size(); i++) {
		if (histo[i] - stdvalue < 0) {
			stdi = i;
			break;
		}
	}

	double stdev = i - SC;
	cout << "stdi = " << i << " stdev = " << stdev << endl;
	vector<long> histo2;
	vector<double> ErrorModel;
	vector<double> ErrorDist;

	for (int i = 0; i < histo.size(); i++) {
		histo2.push_back(histo[i]);
		ErrorModel.push_back(histo[i]);
	}

	double burner;
	FitErrorModel(ErrorModel, burner, Inflection);

	for (int i = 0; i < histo.size(); i++) {
		ErrorDist.push_back(ErrorModel[i] / burner);
		if (histo[i] - (ErrorModel[i]) > 0)
			histo2[i] = histo[i] - (ErrorModel[i]);
		else
			histo2[i] = 0;
	}

	double factor = 1;
	double skew = 0;
	double Power = 1;
	double Rsq = -1;
	double bestS = stdev;
	double bestF = factor;
	double bestSC = SC;
	double bestSK = skew;
	double bestP = Power;

	for (i = 0; i <= 2; i++) {
		cout << "On " << i + 1 << " pass" << endl;
		double count;
		double Flow = 1;
		double Fhigh = 20;
		count = 0;

		while (Flow / Fhigh < .999 and Fhigh > 1e-10) {
			count++;
			double values[9];

#pragma omp parallel for num_threads(11)
			for (int x = 0; x <= 10; x++) {
				values[x] =
				testModelLog(bestSC, bestS, Flow + (((Fhigh - Flow) / 10) * x),
				 bestSK, bestP, histo2, Inflection, 5, Ybar, Rsq);
			}

			double lowest = values[0];
			int lowestX = 0;

			for (int x = 1; x <= 10; x++) {
				if (values[x] < lowest) {
					lowestX = x;
					lowest = values[x];
				}
			}

			if (Flow + ((Fhigh - Flow) / 10) * (lowestX - 1) >= 0) {
				Flow = Flow + ((Fhigh - Flow) / 10) * (lowestX - 1);
			} else {
				Flow = 0;
			}

			Fhigh = Flow + ((Fhigh - Flow) / 10) * (lowestX + 1);
			bestF = Flow + ((Fhigh - Flow) / 10) * (lowestX);
		}
		cout << "	 best Factor = " << bestF << " steps = " << count << endl;
		double SClow = SC * .9;
		double SChigh = SC * 1.1;
		count = 0;

		while (SClow / SChigh < .999 and SChigh > 1e-50) {
			count++;
			double values[9];

#pragma omp parallel for num_threads(11)
			for (int x = 0; x <= 10; x++) {
				values[x] = testModel(SClow + ((SChigh - SClow) / 10) * x, bestS, bestF,
															bestSK, bestP, histo2, Inflection, 5, Ybar, Rsq);
			}
			double lowest = values[0];
			int lowestX = 0;
			for (int x = 1; x <= 10; x++) {
				if (values[x] < lowest) {
					lowestX = x;
					lowest = values[x];
				}
			}

			if (SClow + ((SChigh - SClow) / 10) * (lowestX - 1) >= 0) {
				SClow = SClow + ((SChigh - SClow) / 10) * (lowestX - 1);
			} else {
				SClow = 0;
			}

			SChigh = SClow + ((SChigh - SClow) / 10) * (lowestX + 1);
			bestSC = SClow + ((SChigh - SClow) / 10) * (lowestX);
		}
		cout << "		bestSC = " << bestSC << " steps = " << count << endl;
		double STlow = stdev * .9;
		double SThigh = stdev * 1.1;
		count = 0;

		while (STlow / SThigh < .99 and SThigh > 1e-50) {
			count++;
			double values[9];

#pragma omp parallel for num_threads(11)
			for (int x = 0; x <= 10; x++) {
				values[x] =
		testModel(bestSC, STlow + ((SThigh - STlow) / 10) * x, bestF,
				bestSK, bestP, histo2, Inflection, 5, Ybar, Rsq);
			}

			double lowest = values[0];
			int lowestX = 0;

			for (int x = 1; x <= 10; x++) {
				if (values[x] < lowest) {
					lowestX = x;
					lowest = values[x];
				}
			}

			if (STlow + ((SThigh - STlow) / 10) * (lowestX - 1) >= 0) {
				STlow = STlow + ((SThigh - STlow) / 10) * (lowestX - 1);
			} else {
				STlow = 0;
			}

			SThigh = STlow + ((SThigh - STlow) / 10) * (lowestX + 1);
			bestS = STlow + ((SThigh - STlow) / 10) * (lowestX);
		}

		cout << "		best StdDev = " << bestS << " steps = " << count << endl;
		double SKlow = 0;
		double SKhigh = 2;
		count = 0;

		while (SKlow / SKhigh < .999 and SKhigh < 1e-50) {
			count++;
			double values[9];

#pragma omp parallel for num_threads(11)
			for (int x = 0; x <= 10; x++) {
				values[x] = testModelLog(bestSC, bestS, bestF,
																 SKlow + ((SKhigh - SKlow) / 10) * x, bestP,
																 histo2, Inflection, 5, Ybar, Rsq);
			}

			double lowest = values[0];
			int lowestX = 0;

			for (int x = 1; x <= 10; x++) {
				if (values[x] < lowest) {
					lowestX = x;
					lowest = values[x];
				}
			}

			if (SKlow + ((SKhigh - SKlow) / 10) * (lowestX - 1) >= 0) {
				SKlow = SKlow + ((SKhigh - SKlow) / 10) * (lowestX - 1);
			} else {
				SKlow = 0;
			}

			SKhigh = SKlow + ((SKhigh - SKlow) / 10) * (lowestX + 1);
			bestSK = SKlow + ((SKhigh - SKlow) / 10) * (lowestX);

			if (SKlow < 1e-20 and SKhigh < 1e-20) {
				break;
			}
		}

		cout << "		best skew factor = " << bestSK << " steps = " << count << endl;
		double Plow = 1;
		double Phigh = 2;
		count = 0;

		while (Plow / Phigh < .999 and Phigh > 1e-50) {
			count++;
			double values[9];

#pragma omp parallel for num_threads(11)
			for (int x = 0; x <= 10; x++) {
				values[x] = testModelLog(bestSC, bestS, bestF, bestSK,
																 Plow + ((Phigh - Plow) / 10) * x, histo2,
																 Inflection, 5, Ybar, Rsq);
			}

			double lowest = values[0];
			int lowestX = 0;

			for (int x = 1; x <= 10; x++) {
				if (values[x] < lowest) {
					lowestX = x;
					lowest = values[x];
				}
			}

			if (Plow + ((Phigh - Plow) / 10) * (lowestX - 1) >= 1) {
				Plow = Plow + ((Phigh - Plow) / 10) * (lowestX - 1);
			} else {
				Plow = 1;
			}

			Phigh = Plow + ((Phigh - Plow) / 10) * (lowestX + 1);
			bestP = Plow + ((Phigh - Plow) / 10) * (lowestX);
		}

		cout << "		best Power factor = " << bestP << " steps = " << count << endl;
		stdev = bestS;
		factor = bestF;
		SC = bestSC;
		factor = bestF;
		skew = bestSK;
		Power = bestP;
	}

	vector<vector<double> > model;
	vector<vector<double> > dist;
	vector<vector<double> > prob;
	vector<double> placeHolder;
	prob.push_back(placeHolder);
	stdev = bestS;
	factor = bestF;
	SC = bestSC;
	skew = bestSK;
	Power = bestP;
	cout << "Best Model is SC = " << SC << " StdDev = " << stdev
			 << " F = " << factor << " skew = " << skew << " bestP = " << Power
			 << endl;

	for (int i = 0; i < histo.size(); i++) {
		vector<double> d;
		d.push_back(0);
		if (Diploid) {
			d.push_back(
			norm(i, SC / 2,
					 stdev * (1 - ((1 - (stdev / (stdev + ((1) * factor)))) / 2)),
					 skew, Power));
		}

		for (long j = 1; j < histo.size() / SC; j++) {
			d.push_back(norm(i, SC * j, stdev + ((j - 1) * factor), skew, Power));
		}
		dist.push_back(d);
	}
	// data liklyhood correction for non symettical distribution
	for (long j = 1; j < histo.size() / SC; j++) {
		double sum = 0.0;
		for (int i = 0; i < histo.size(); i++) {
			sum += dist[i][j];
		}
		for (int i = 0; i < histo.size(); i++) {
			dist[i][j] = dist[i][j] / sum;
		}
	}
	// marginal probabilites
	for (int i = 0; i < histo.size(); i++) {
		vector<double> p;
		double total = 0;
		p.push_back(0);

		for (long j = 1; j < histo.size() / SC; j++) {
			total += dist[i][j];
		}

		for (long j = 1; j < histo.size() / SC; j++) {
			p.push_back(dist[i][j] / total);
		}
		prob.push_back(p);
	}

	vector<double> RC;
	RC.push_back(0);

	if (Diploid) {
		double tSC = histo[SC] / dist[SC][2];
		double hetKmers =
			(histo[SC / 2] - (dist[SC / 2][2] * tSC)) / dist[SC / 2][1];
		if (hetKmers > 0) {
			RC.push_back(hetKmers);
		} else {
			RC.push_back(0);
		}

		RC.push_back(histo[SC] / dist[SC][2]);

	} else {
		RC.push_back(histo[SC] / dist[SC][1]);
	}

	int max = histo.size() / SC;

	if (Diploid) {
		max = max - 1;
	}

	for (long a = 2; a < histo.size() / SC; a++) {
		double Reads = 0;
		if (Diploid) {
			Reads =
	((double)histo[SC * a] / dist[SC * a][a + 1] * prob[SC * a][a + 1]);
		} else {
			Reads = ((double)histo[SC * a] / dist[SC * a][a] * prob[SC * a][a]);
		}

		RC.push_back(Reads);
	}

	vector<double> ModelSums;

	for (long i = 0; i < histo.size(); i++) {
		vector<double> m;
		double sum = 0;
		m.push_back(0);

		for (long j = 1; j < histo.size() / SC; j++) {
			m.push_back(dist[i][j] * RC[j]);
			sum += m[j];
		}

		model.push_back(m);
		ModelSums.push_back(sum);
	}

	double GenomeSize = 0;
	for (int i = 1; i < RC.size(); i++) {
		GenomeSize += RC[i] * i;
	}

	cout << "GenomeSize = " << GenomeSize << endl;

//	for (long i = 1; i < SC; i++) {
//		prob[i][1] = model[i][1] / histo[i];
//	}

	int CutOff = -1;

	for (long KmerCount = 2; KmerCount < SC * 5; KmerCount++) {
		float ModelSum = 0;

		for (long CopyNumber = 1; CopyNumber < 10; CopyNumber++) {
			ModelSum += model[KmerCount][CopyNumber];
		}

		if (ModelSum / ((float)ErrorModel[KmerCount] + ModelSum) > 0.5) {
			CutOff = KmerCount;
		}
	}

	int kcutoff = -1;

	for (long k = 1; k < dist.size(); k++) {
		float num = 0;

		for (long c = 1; c < dist[1].size(); c++) {
			num += dist[k][c];
		}

		cout << "prob not error = " << num / (num + ErrorDist[k]) << endl;

		if (num / (num + ErrorDist[k]) > 0.5) {
			kcutoff = k ;//- 1;
			cout << "this one" << endl;
			break;
		}
	}
	
	ofstream ProbFile;
        path = argv[1];
        path += ".7.7.prob";
        ProbFile.open(path.c_str());

        if (ProbFile.is_open()) {
        } else {
                cout << "Error, Prob file could not be opened";
                return 0;
        }
	
	cout << "here1" << endl; 
	ModelFile << 3 << endl << kcutoff << endl;	// Inflection << endl;
	DistFile << 3 << endl << kcutoff << endl; 
	ProbFile <<  3 << endl << kcutoff << endl;
	ModelFile << HistoSum << endl;
	DistFile << HistoSum << endl;
	ProbFile << HistoSum << endl;
	ModelFile << rawSC << endl;  
	DistFile << rawSC << endl;
	ProbFile << rawSC << endl;
	cout << "here" << endl; 

	for (long CopyNumber = 1; CopyNumber < model[1].size(); CopyNumber++) {
		long LocalSum = 0;

		for (long KmerCount = 1; KmerCount < histo.size(); KmerCount++) {
			LocalSum += model[KmerCount][CopyNumber];
		}

		ModelFile << ((double)LocalSum) / ((double)HistoSum) << '\t';
	}

	DistFile << burner << '\t' << 0 << '\t';

	for (int i = 1; i < RC.size(); i++) {
		DistFile << RC[i] << '\t';
	}

	ModelFile << endl;
	DistFile << endl;
	ModelFile << "K\tRawCount\tErrorModel\tContSubtract\tModelSum\t1x\t2x\t3x\t4x\t5x\t6x\t7x\t8x\t9x"<< endl;
	ModelFile << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0;
	ModelFile << '\t' << 0;

	for (long CopyNumber = 1; CopyNumber < 10; CopyNumber++) {
		ModelFile << '\t' << model[0][CopyNumber];
	}

	ModelFile << endl;

	for (long KmerCount = 1; KmerCount < SC * 5; KmerCount++) {
		ModelFile << KmerCount << '\t' << histo[KmerCount] << '\t'
							<< ErrorModel[KmerCount] << '\t' << 0;
		ModelFile << '\t' << ModelSums[KmerCount];

		for (long CopyNumber = 1; CopyNumber < 10; CopyNumber++) {
			ModelFile << '\t' << model[KmerCount][CopyNumber];
		}
		ModelFile << endl;
	}
	ModelFile.close();
	DistFile << SC << endl;
	DistFile << 0 << '\t' << 0 << '\t' << 0;

	for (long c = 1; c < dist[1].size(); c++) {
		DistFile << '\t' << dist[0][c];
	}

	DistFile << endl;

	for (long k = 1; k < dist.size(); k++) {
		DistFile << k << '\t' << ErrorDist[k] << '\t' << 0;

		for (long c = 1; c < dist[1].size(); c++) {
			DistFile << '\t' << dist[k][c];
		}
		DistFile << endl;
	}

	for (long c = 1; c < dist[1].size(); c++) {
                DistFile << '\t' << dist[0][c];
        }
	 DistFile.close();
        
	
	 ProbFile << endl << endl; 
        for (long k = 1; k < prob.size(); k++) {
                ProbFile << k << '\t' << ErrorDist[k] << '\t' << 0;

                for (long c = 1; c < prob[1].size(); c++) {
                        ProbFile << '\t' <<prob[k][c];
                }
                ProbFile << endl;
        }

	ProbFile.close();
	cout << "GenomeSize = " << GenomeSize << endl;
	cout << "Inflection point = " << Inflection << endl;
	cout << "Recomended RUFUS cutoff = " << SC - (5 * stdev) << endl;
	cout << "-1std = " << SC - (1 * stdev) << "	-2std = " << SC - (2 * stdev)
			 << "	-3std = " << SC - (3 * stdev) << "	-4std = " << SC - (4 * stdev)
			 << endl;
}
