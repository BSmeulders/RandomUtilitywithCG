// RUM_Computation_Source.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// QPColumn2.cpp : Defines the entry point for the console application.
//
using namespace std;

#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <stack>
#include <chrono>
#include <cassert>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "Typedef.h"
#include "Structures.h"

configuration Readconfig(char configname[], int* error);
datastructure Readdata(const configuration& config, int* error);
ThreeIntMatrix Create_Patches(const configuration& config, TwoDoubleMatrix prices);
void Test_PatchReverse(ThreeIntMatrix& all_patches, const configuration& config, TwoDoubleMatrix p, vector<int> patch, int depth, int budget);
TwoDoubleMatrix Initial_RHS(datastructure& data, const configuration& config);


// EndogenousStastic: Code of Kitamura and Stoyes RUM_7X Series.
void Endogenous_Stastic(datastructure& data, const configuration& config); // In Endogenous.cpp

// Functions for random generation of patterns
TwoIntMatrix Generate_Pattern(const configuration& config, const ThreeIntMatrix& patches);
ThreeIntMatrix Generate_Pattern_Set(const configuration& config, const ThreeIntMatrix& patches);
void change_budget(vector<int> SCC, const ThreeIntMatrix& patches, const configuration& config, TwoIntMatrix& rp_relations, vector<int>& chosen_patches);
TwoIntMatrix strongly_connected_components(int Graph_Size, TwoIntMatrix RP); // In SCC.cpp

// Misc functions
double vector_median(vector<double> vector);
ThreeDoubleMatrix pi_hat_read_dirty(const configuration config, const datastructure data);


int main(int argc, char* argv[])
{
	int error = 0; // Flag for errors.

	configuration config = Readconfig(argv[1], &error);
	if (error == 1)
	{
		cout << "Error in reading config, ending program.";
		return 0;
	}
	
	datastructure data = Readdata(config, &error);
	if (error == 1)
	{
		cout << "Error in reading data, ending program";
		return 0;
	}
	Endogenous_Stastic(data, config);
		
	return 0;
}


TwoIntMatrix Generate_Pattern(const configuration& config, const ThreeIntMatrix& patches)
{
	// The pattern that is for each budget the choice of one alternative.
	TwoIntMatrix pattern; pattern.resize(config.size); for (int i = 0; i < config.size; i++) { pattern[i].resize(patches[i].size()); }
	// The Revealed preference matrix is a zero/one matrix showing the (reversed) RP relations. RP[0][1] == 1 means the patch chosen on 
	// budget 1 is below budget 2.
	TwoIntMatrix rp_relations; rp_relations.resize(config.size); for (int i = 0; i < config.size; i++) { rp_relations[i].resize(config.size); }
	TwoIntMatrix SCC;
	// chosen_patch gives for each budget the number of the chosen patch.
	vector<int> chosen_patches(config.size);

	// First, pick a random starting patch for each time period.
	for (int i = 0; i < config.size; i++)
	{
		int nr_patches = patches[i].size();
		int chosen_patch = rand() % nr_patches;
		chosen_patches[i] = chosen_patch;
		rp_relations[i] = patches[i][chosen_patch];
	}

	chrono::high_resolution_clock::time_point start_time, end_time;
	start_time = chrono::high_resolution_clock::now();
	double time_elapsed = 0;

	// Loop which fixes the type.
	for (bool fixed = 0; time_elapsed < 10 && fixed == 0;)
	{
		// Get the Strongly Connected Components.
		SCC = strongly_connected_components(config.size, rp_relations);
		// If there is no SCC with more than one element, break.
		if (SCC.size() == config.size) break;
		// If there is a SCC with more than one element, pick a random budget in it.
		// For that random budget, pick a patch with minimal changes, but less incoming arcs within the SCC.
		for (int i = 0; i < SCC.size(); i++)
		{
			if (SCC[i].size() > 1)
				change_budget(SCC[i], patches, config, rp_relations, chosen_patches);
		}
		end_time = chrono::high_resolution_clock::now();
		time_elapsed = chrono::duration_cast<chrono::duration<double>>(end_time - start_time).count();
	}
	if (time_elapsed < 10)
	{
		// Go from chosen patches to pattern.
		for (int i = 0; i < config.size; i++)
		{
			for (int j = 0; j < pattern[i].size(); j++)
			{
				if (chosen_patches[i] == j)
					pattern[i][j] = 1;
				else
					pattern[i][j] = 0;
				//cout << pattern[i][j] << " ";
			}
			//cout << endl;
		}
	}
	else // Else we are stuck in a loop, and we just try again to build a pattern.
	{
		pattern = Generate_Pattern(config, patches);
	}



	return pattern;
}
ThreeIntMatrix Generate_Pattern_Set(const configuration& config, const ThreeIntMatrix& patches)
{
	ThreeIntMatrix pattern_set;

	while (pattern_set.size() < config.start_vars)
	{
		TwoIntMatrix new_pattern = Generate_Pattern(config, patches);
		for (int i = 0; i < pattern_set.size() + 1; i++)
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

		}
	}

	return pattern_set;
}
void change_budget(vector<int> SCC, const ThreeIntMatrix& patches, const configuration& config, TwoIntMatrix& rp_relations, vector<int>& chosen_patches)
{
	// We choose a random budget in the SCC to change.
	int budget_to_change = SCC[rand() % SCC.size()];

	vector<int> patch_scores; patch_scores.resize(patches[budget_to_change].size());
	// For each patch, we compute a score. A patch gets 1 point per removed RP relation, 5 points per added RP relation 
	// and 999 points if all relations within the SCC are the same as the current patch.

	for (int i = 0; i < patches[budget_to_change].size(); i++)
	{
		for (int j = 0; j < config.size; j++)
		{
			if (rp_relations[budget_to_change][j] == 1 && patches[budget_to_change][i][j] == 0)
				patch_scores[i] += 1;
			if (rp_relations[budget_to_change][j] == 0 && patches[budget_to_change][i][j] == 1)
				patch_scores[i] += 5;
		}
		int equal_in_SCC = 0;
		for (int j = 0; j < SCC.size(); j++)
		{
			if (rp_relations[budget_to_change][SCC[j]] == patches[budget_to_change][i][SCC[j]])
				equal_in_SCC++;
		}
		if (equal_in_SCC == SCC.size()) patch_scores[i] += 999;
	}

	// Find the minimum patch score
	int min_score = 9999;
	int new_patch = 0;
	for (int i = 0; i < patches[budget_to_change].size(); i++)
	{
		if (patch_scores[i] < min_score)
		{
			min_score = patch_scores[i];
			new_patch = i;
		}
	}
	// Change the patches.
	chosen_patches[budget_to_change] = new_patch;
	rp_relations[budget_to_change] = patches[budget_to_change][new_patch];
	/*for (int i = 0; i < chosen_patches.size(); i++)
	{
		cout << chosen_patches[i] << "\t";
	}
	cout << endl;*/
}

configuration Readconfig(char configname[], int* error)
{
	configuration config;
	namespace po = boost::program_options;
	ifstream configfile(configname);
	if (!configfile)
	{
		cout << "No configuration file found." << endl;
		*error = 1;
	}
	po::options_description file_options(configname);
	file_options.add_options()
		("Quadratic Program", po::value<bool>(&config.QP)->default_value(1), "Quadratic Program Switch")
		("Linear Program", po::value<bool>(&config.LP)->default_value(0), "Linear Program Switch")
		("Heuristics", po::value<int>(&config.heuristics)->default_value(0), "Heuristics Choice")
		("Repetitions", po::value<int>(&config.repetitions)->default_value(100), "The number of repetitions in the bootstrap procedure")
		("Startyear", po::value<int>(&config.startyear)->required(), "The first year used for the analysis")
		("Endyear", po::value<int>(&config.endyear)->required(), "The last year used in the analysis")
		("Products", po::value<int>(&config.nr_products)->default_value(3), "The number of product categories in the analysis.")
		("Start Variables", po::value<int>(&config.start_vars)->default_value(1), "The initial size of the array of the variables in the Column Generation models")
		("Tau Variables", po::value<int>(&config.tau_vars)->default_value(0), "The number of rational types generated for the tightening.")
		("Polynomial Degree", po::value<int>(&config.poly_degree)->default_value(3), "Polynomial degree used for endogneous estimation.")
		("Recycle", po::value<int>(&config.recycle)->default_value(0), "Switch for variable recycling in the CG.")
		("Recycle Credits", po::value<int>(&config.recycle_creds)->default_value(50), "Number of times variable can be unused before it is a candidate for recycling.")
		("Upper Bound", po::value<int>(&config.UB_Threshold)->default_value(0), "1 if the computation of Jstat for bootstrap is terminated if the upper bound for this value is below Jstat.")
		("Lower Bound", po::value<int>(&config.LB_Threshold)->default_value(0), "Depreciated")
		("Time Limit", po::value<int>(&config.time_limit)->default_value(3600));
	po::variables_map vm;
	ifstream ifs(configname);
	if (!ifs)
	{
		cout << "Failed to open Configuration File" << endl;
		*error = 1;
		return config;
	}
	else
	{
		try {
			store(po::parse_config_file(ifs, file_options), vm);
			notify(vm);
		}
		catch (exception e)
		{
			cout << "Required values not included in config file" << endl; cout << "Press enter to continue" << endl; cin.get(); *error = 1;
		}
	}
	if (config.QP + config.LP != 1)
	{
		cout << "Choice QP or LP not clear." << endl; cout << "Press enter to continue" << endl; cin.get(); *error = 1;
	}
	if (config.startyear < 75 || config.endyear < 75 || config.startyear > 99 || config.endyear > 99 || config.startyear >= config.endyear)
	{
		cout << "Start and/or Endyear are out of range for the dataset, or they are inconsistent (startyear > endyear)" << endl << "Please try again." << endl;
		*error = 1;
	}
	config.size = config.endyear - config.startyear + 1;
	return config;
}
datastructure Readdata(const configuration& config, int* error)
{
	datastructure data;
	vector<int> dominatingbudgets;
	string filename;
	ostringstream convert;
	bool flag = 0;

	data.prices.resize(config.size);
	data.quant.resize(config.size);
	data.shares.resize(config.size);
	data.income.resize(config.size);
	data.instrument.resize(config.size);
	data.patches.resize(config.size); // We expand to make room for all budgets
	data.frequency.resize(config.size);
	data.HH_year.resize(config.size);
	dominatingbudgets.resize(config.size);

	for (int i = 0; i < config.size; i++) // Reading in of all the prices
	{
		data.prices[i].resize(config.nr_products);
		convert << "input/" << config.nr_products << "goods-p" << config.startyear + i << ".csv";
		filename = convert.str(); convert.str("");
		ifstream inputfile(filename);
		if (!inputfile) // Check for existence of file, if not goto end.
		{
			cout << endl << "Failed to open file " << filename << endl;;
			*error = 1;
			return data;
		}
		for (int j = 0; j < config.nr_products; j++) inputfile >> data.prices[i][j];

		// Reading of the HH specific data.
		convert << "input/" << config.nr_products << "goods-q" << config.startyear + i << ".csv";
		filename = convert.str(); convert.str("");
		ifstream HHinputfile(filename);
		if (!HHinputfile) // Check for existence of file, if not goto end.
		{
			cout << endl << "Failed to open file " << filename << endl;;
			*error = 1;
			return data;
		}
		HHinputfile >> data.HH_year[i];
		data.shares[i].resize(data.HH_year[i]); data.income[i].resize(data.HH_year[i]); data.instrument[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			data.shares[i][j].resize(config.nr_products);
			for (int k = 0; k < config.nr_products; k++)
			{
				HHinputfile >> data.shares[i][j][k]; HHinputfile.ignore(1);
			}
			HHinputfile >> data.income[i][j]; HHinputfile.ignore(1);
			HHinputfile >> data.instrument[i][j];
		}
	}

	data.patches = Create_Patches(config, data.prices);

	if (flag == 1)
	{
		*error = 1;
	}
	return data;
}

ThreeIntMatrix Create_Patches(const configuration& config, TwoDoubleMatrix prices)
{
	vector<int> current_patch; current_patch.resize(config.size);
	ThreeIntMatrix all_patches(config.size);


	for (int i = 0; i < config.size; i++) // For each budget we create using Reverse patches.
	{
		if (i != config.size - 1) // We will always put patch[i] = 0 for budget i.
		{
			current_patch[config.size - 1] = 1;
			Test_PatchReverse(all_patches, config, prices, current_patch, 1, i);
		}

		current_patch[config.size - 1] = 0;
		Test_PatchReverse(all_patches, config, prices, current_patch, 1, i);
	}

	ofstream outputpatches;
	outputpatches.open("zPatches.txt");
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < all_patches[i].size(); j++)
		{
			for (int k = 0; k < all_patches[i][j].size(); k++)
			{
				outputpatches << all_patches[i][j][k] << " ";
			}
			outputpatches << endl;
		}
		outputpatches << endl << endl << endl;
	}


	return all_patches;
}
void Test_PatchReverse(ThreeIntMatrix& all_patches, const configuration& config, TwoDoubleMatrix p, vector<int> patch, int depth, int budget)
{
	// Test_Patch works iteratively. We test whether a patch can exist with the given relations to other budgets, up to the given depth. (i.e. if the depth is 3, we will only look a the last 3 budgets).
	// If a solution exists, we will continue by setting the budget just beyond the considered depth to 0 or 1. This determines whether we test for the patch being above or below that budget.
	// We then increase the depth and run Test_Patch again. 

	// If the depth is equal to the number of budgets, a positive result to the test will result in the patch being added to the list.

	// We will always set a constraint that the p*q should be equal to 1 for the budget for which we are identifying the patches, no matter the depth.

	// Setting up the LP model
	IloEnv env;
	IloModel model(env);
	IloNumVar slack_var(env, 0, IloInfinity, ILOFLOAT); // Slack_var to allow for some small numerical instability.
	IloNumVarArray q(env, config.nr_products, 0, IloInfinity, ILOFLOAT);
	IloObjective obj = (IloMinimize(env, slack_var)); model.add(obj);
	IloExpr equal(env);
	for (int i = 0; i < config.nr_products; i++) equal += q[i] * p[budget][i];
	IloRange equalrange = IloRange(equal == 1); model.add(equalrange);
	IloRangeArray unequal(env, depth);
	for (int i = config.size - 1; i > config.size - 1 - depth; i--)
	{
		if (patch[i] == 0)
		{
			unequal[config.size - 1 - i] = IloRange(slack_var >= 1);
			for (int j = 0; j < config.nr_products; j++) unequal[config.size - 1 - i].setLinearCoef(q[j], p[i][j]);
		}
		else if (patch[i] == 1)
		{
			unequal[config.size - 1 - i] = IloRange(-slack_var <= 1);
			for (int j = 0; j < config.nr_products; j++) unequal[config.size - 1 - i].setLinearCoef(q[j], p[i][j]);
		}
		else
		{
			cout << "Error" << endl;
			cin.get();
		}
	}
	model.add(unequal);
	IloCplex cplex_model(model);
	cplex_model.setOut(env.getNullStream());

	cplex_model.solve(); // We solve the model

	double slack = cplex_model.getObjValue();
	cplex_model.clear(); q.end(); slack_var.end(); obj.end(); /*equal1.end(); equal2.end();*/ unequal.end();
	cplex_model.end();
	env.end();

	// First, the case where the depth is lower than the number of budgets we consider.
	if (depth < config.size && slack <= 0.00000001) // We leave some slack for cases with numerical instability.
	{
		if (config.size - 1 - depth != budget) // We will never set patch[depth] = 1 for budget "budget" = depth.
		{
			patch[config.size - 1 - depth] = 1;
			Test_PatchReverse(all_patches, config, p, patch, depth + 1, budget);
		}

		patch[config.size - 1 - depth] = 0;
		Test_PatchReverse(all_patches, config, p, patch, depth + 1, budget);
	}

	if (depth == config.size && slack <= 0.0000001)
	{
		all_patches[budget].push_back(patch);
	}

}

TwoDoubleMatrix Initial_RHS(datastructure& data, const configuration& config)
{
	TwoDoubleMatrix lp_rhs;
	lp_rhs.resize(config.size);

	// This loop creates the data that will be given to the solver. From the data frequencies, it converts to ratios.
	for (int i = 0; i < config.size; i++)
	{
		int sum = 0;
		for (int j = 0; j < data.frequency[i].size(); j++)
		{
			sum = sum + data.frequency[i][j];
		}
		lp_rhs[i].resize(data.frequency[i].size());
		for (int j = 0; j < data.frequency[i].size(); j++)
		{
			lp_rhs[i][j] = (double)data.frequency[i][j] / (double)sum;
		}
	}

	return lp_rhs;
}
