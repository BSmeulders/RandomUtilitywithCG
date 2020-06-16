#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <stack>
// #include <boost/math/distributions.hpp>
#include <cassert>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
using namespace std;

#include "Typedef.h"
#include "Structures.h"

vector<double> vector_median(TwoDoubleMatrix doublevector);
double vector_median(vector<double> vector); // In QPColumn2.cpp
ThreeDoubleMatrix pi_hat_read_dirty(const configuration config, const datastructure data);
TwoIntMatrix Generate_Pattern(const configuration &config, const ThreeIntMatrix &patches);
TwoIntMatrix Generate_Pattern_Pure_Random(const configuration &config, const ThreeIntMatrix &patches);
int Random_Pattern_test(const configuration &config, const ThreeIntMatrix &patches);
TwoIntMatrix strongly_connected_components(int Graph_Size, TwoIntMatrix RP);

vector<double> vector_median(TwoDoubleMatrix doublevector)
{
	vector<double> median_vector(doublevector.size());

	for (int i = 0; i < doublevector.size(); i++)
		median_vector[i] = vector_median(doublevector[i]);
	
	return median_vector;
}

double vector_median(vector<double> vector)
{
	double median = 0;

	sort(vector.begin(), vector.end());
	if (vector.size() % 2 == 1)
		median = vector[floor(0.5*vector.size())];
	else
		median = (vector[0.5*vector.size()] + vector[0.5*vector.size() - 1]) / 2;

	return median;
}

vector<vector<double>> variance_pi_hat(const configuration & config, const ThreeDoubleMatrix & bs_pi_hat)
{
	vector<vector<double>> pi_hat_var(config.size);
	for (int i = 0; i < config.size; i++)
		pi_hat_var[i].resize(bs_pi_hat[0][i].size());

	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < pi_hat_var[i].size(); j++)
		{
			double mean = 0;
			for (int k = 0; k < bs_pi_hat.size(); k++)
			{
				mean = mean + bs_pi_hat[k][i][j];
			}
			mean = mean / bs_pi_hat.size();
			double temp = 0;
			for (int k = 0; k < bs_pi_hat.size(); k++)
			{
				temp = temp + (bs_pi_hat[k][i][j]-mean) * (bs_pi_hat[k][i][j] - mean);
			}
			pi_hat_var[i][j] = temp / bs_pi_hat.size();
		}
	}
	
	return pi_hat_var;
}

ThreeDoubleMatrix pi_hat_read_dirty(const configuration config, const datastructure data)
{
	ThreeDoubleMatrix pi_hat(1);
	string filename = "pi_hat.csv";
	ifstream inputfile(filename);

	if (!inputfile) // Check for existence of file, if not goto end.
	{
		cout << endl << "Failed to open file " << filename << endl;;
		cin.get();
	}

	pi_hat[0].resize(config.size);
	for (int i = 0; i < config.size; i++)
	{
		pi_hat[0][i].resize(data.patches[i].size());
		for (int j = 0; j < data.patches[i].size(); j++)
		{
			inputfile >> pi_hat[0][i][j];
			cout << pi_hat[0][i][j] << "\t";
		}
		cout << endl;
	}
	return pi_hat;
}

ThreeIntMatrix Generate_Tau_Set(const configuration & config, const ThreeIntMatrix & patches)
{
	ThreeIntMatrix pattern_set;
	srand(config.startyear*config.endyear*config.nr_products);
	while (pattern_set.size() < config.tau_vars)
	{
		TwoIntMatrix new_pattern = Generate_Pattern(config, patches);
		/*for (int i = 0; i < pattern_set.size() + 1; i++)
		{
			if (i == pattern_set.size()) // If the loop has not been broken and the last pattern already in the set has been tested, add the new pattern to the set.
			{
				//cout << "Adding new pattern in position " << i << endl;
				pattern_set.push_back(new_pattern);
				break;
			}
			if (new_pattern == pattern_set[i]) // If the pattern is the same as the previous pattern, break out of the loop.
			{
				//cout << "New pattern is equal to pattern " << i << endl;
				break;
			}
		}*/
		pattern_set.push_back(new_pattern);
	}

	return pattern_set;
}

ThreeIntMatrix Pure_Generate_Tau_Set(const configuration & config, const ThreeIntMatrix & patches)
{
	ThreeIntMatrix pattern_set;
	srand(config.startyear*config.endyear*config.nr_products);
	while (pattern_set.size() < config.tau_vars)
	{
		TwoIntMatrix new_pattern = Generate_Pattern_Pure_Random(config, patches);
		pattern_set.push_back(new_pattern);
	}

	return pattern_set;
}

int Random_Pattern_test(const configuration & config, const ThreeIntMatrix & patches)
{
	int number_tested = 0;
	// The pattern that is for each budget the choice of one alternative.
	TwoIntMatrix pattern; pattern.resize(config.size); for (int i = 0; i < config.size; i++) { pattern[i].resize(patches[i].size()); }
	// The Revealed preference matrix is a zero/one matrix showing the (reversed) RP relations. RP[0][1] == 1 means the patch chosen on 
	// budget 1 is below budget 2.
	TwoIntMatrix rp_relations; rp_relations.resize(config.size); for (int i = 0; i < config.size; i++) { rp_relations[i].resize(config.size); }
	TwoIntMatrix SCC;
	// chosen_patch gives for each budget the number of the chosen patch.
	vector<int> chosen_patches(config.size);

	// Loop which fixes the type.
	int correct = 0;
	while (correct < 100)
	{
		number_tested++;
		// First, pick a random starting patch for each time period.
		for (int i = 0; i < config.size; i++)
		{
			int nr_patches = patches[i].size();
			int chosen_patch = rand() % nr_patches;
			chosen_patches[i] = chosen_patch;
			rp_relations[i] = patches[i][chosen_patch];
		}
		
		// Get the Strongly Connected Components.
		SCC = strongly_connected_components(config.size, rp_relations);
		// If there is no SCC with more than one element, break.
		if (SCC.size() == config.size)
		{
			correct++;
		}
		if (number_tested % 1000 == 0)
			cout << number_tested << "\t" << correct << endl;
	}
	return number_tested;
}

TwoIntMatrix Generate_Pattern_Pure_Random(const configuration & config, const ThreeIntMatrix & patches)
{
	TwoIntMatrix rp_relations; rp_relations.resize(config.size); for (int i = 0; i < config.size; i++) { rp_relations[i].resize(config.size); }
	vector<int> chosen_patches(config.size);
	TwoIntMatrix SCC;
	while (SCC.size() < config.size)
	{
		// First, pick a random starting patch for each time period.
		for (int i = 0; i < config.size; i++)
		{
			int nr_patches = patches[i].size();
			int chosen_patch = rand() % nr_patches;
			chosen_patches[i] = chosen_patch;
			rp_relations[i] = patches[i][chosen_patch];
		}
		// Test for strongly connected components
		SCC = strongly_connected_components(config.size, rp_relations);
	}
	TwoIntMatrix pattern(config.size);
	for (int i = 0; i < config.size; i++)
	{
		pattern[i].resize(patches[i].size(),0);
		pattern[i][chosen_patches[i]] = 1;
	}
	return pattern;
}