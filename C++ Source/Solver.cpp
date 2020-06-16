#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <stack>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cassert>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
using namespace std;

#include "Typedef.h"
#include "Structures.h"

// Depending on the pi_hat tested (untightened, tightened, bootstrap), different things need to be returned.
TwoDoubleMatrix QP_CG_Untightened(const configuration & config, const TwoDoubleMatrix & pi, ThreeIntMatrix & starting_set, const ThreeIntMatrix & patches, const vector<int> & HH_year, double & Jstat);
vector<double> QP_CG_Bootstrap(const configuration & config, const ThreeDoubleMatrix & pi, const ThreeIntMatrix & patches, const vector<int> & HH_year, double Jstat);
TwoIntMatrix Best_Insertion_Heur(const configuration & config, const ThreeIntMatrix & patches, const TwoDoubleMatrix & distance_vals, double & price_obj);

IloModel IP_Model2(const ThreeIntMatrix &patches, const configuration &config, IloEnv &env, NumVarMatrix patch_vars, NumVarMatrix preference_vars, IloObjective ip_obj);
IloModel IP_Model3(const ThreeIntMatrix &patches, const configuration &config, IloEnv &env, NumVarMatrix patch_vars, NumVarMatrix preference_vars, TwoDoubleMatrix distance_vals);

TwoDoubleMatrix QP_CG_Untightened(const configuration & config, const TwoDoubleMatrix & pi, ThreeIntMatrix & starting_set, const ThreeIntMatrix & patches, const vector<int> & HH_year, double & Jstat)
{
	vector<double> distances(pi.size());
	TwoDoubleMatrix nu_hat; nu_hat.resize(pi.size()); for (int i = 0; i < config.size; i++) { nu_hat[i].resize(pi[i].size()); }

	ofstream outputtypes;
	outputtypes.open("zPatterns.txt");

	// First we set up the Quadratic Program and the Pricing Problem. For the QP parameters, we use those of the first pi input.
	// The QP model.
	IloEnv env;
	NumVarMatrix slack_vars(env, config.size); for (int i = 0; i < config.size; i++) { slack_vars[i] = IloNumVarArray(env, patches[i].size(), -1, 1, ILOFLOAT); }
	ConstraintMatrix qp_cons(env, config.size);
	for (int i = 0; i < config.size; i++)
	{
		qp_cons[i] = IloRangeArray(env, patches[i].size());
		for (int j = 0; j < patches[i].size(); j++)
		{
			qp_cons[i][j] = IloRange(slack_vars[i][j] == pi[i][j]); // If, for a particular dimension, the coordinate of the solution point is higher than the coordinate of the data point, the slack_var is negative.
		}
	}
	IloExpr objective(env); for (int i = 0; i < config.size; i++) { for (int j = 0; j < patches[i].size(); j++) { objective += slack_vars[i][j] * slack_vars[i][j]; } }
	IloObjective obj = IloMinimize(env, objective);
	IloModel qp_model(env); qp_model.add(obj);
	for (int i = 0; i < config.size; i++) { qp_model.add(qp_cons[i]); }
	// Adding the initial variables
	IloNumVarArray sarp_vars(env, starting_set.size());
	for (int i = 0; i < starting_set.size(); i++)
	{
		sarp_vars[i] = IloNumVar(env, 0, 1, ILOFLOAT);
		for (int j = 0; j < config.size; j++)
		{
			for (int k = 0; k < patches[j].size(); k++)
			{
				if (starting_set[i][j][k] == 1)
				{
					qp_cons[j][k].setLinearCoef(sarp_vars[i], 1);
				}
			}
		}
	}

	IloCplex cplex_qp(qp_model);
	cplex_qp.setOut(env.getNullStream());

	// We set up the Pricing Problem. 
	NumVarMatrix patch_vars(env, config.size); 
	for (int i = 0; i < config.size; i++) { patch_vars[i] = IloNumVarArray(env, patches[i].size(), 0, 1, ILOINT); }
	NumVarMatrix preference_vars(env, config.size); 
	for (int i = 0; i < config.size; i++) { preference_vars[i] = IloNumVarArray(env, config.size, 0, 1, ILOINT); }
	
	// A matrix used later to store the output of the QP
	TwoDoubleMatrix distance_vals;
	distance_vals.resize(patches.size()); for (int i = 0; i < patches.size(); i++) { distance_vals[i].resize(patches[i].size()); }
	// Start solving!
	// Solve the Quadratic Model
	cplex_qp.solve();
	int patterns = starting_set.size();
	double distance = cplex_qp.getObjValue();
	while (distance > 0.00000001)
	{
		double target = 0;
		
		// If current distance > threshold, update the pricing problem such that the objective is the hyperplane tangent to the current solution.
		for (int i = 0; i < distance_vals.size(); i++)
		{
			for (int j = 0; j < distance_vals[i].size(); j++)
			{
				distance_vals[i][j] = cplex_qp.getValue(slack_vars[i][j]);
				target += distance_vals[i][j] * (pi[i][j] - distance_vals[i][j]); // pi - distance_vals is current solution.
			}
		}
		// Solve Pricing Problem.
		// First we use a heuristic.
		double price_sol_val;
		TwoIntMatrix new_pattern;
		if (config.heuristics == 1)
		{
			new_pattern = Best_Insertion_Heur(config, patches, distance_vals, price_sol_val);
		}
		else
		{
			price_sol_val = -99;
		}
		
		double ip_sol_val = 0;
		// If the heuristic does not return a good solution, we run the IP Model.
		if (price_sol_val < target + 0.000001)
		{
			chrono::high_resolution_clock::time_point mini_start_time = chrono::high_resolution_clock::now();
			IloModel ip_model = IP_Model3(patches, config, env, patch_vars, preference_vars, distance_vals);
			IloCplex cplex_ip(ip_model);
			cplex_ip.setOut(env.getNullStream());
			cplex_ip.solve();
			chrono::high_resolution_clock::time_point mini_end_time = chrono::high_resolution_clock::now();
			cout << "IP Solve = " << chrono::duration_cast<chrono::duration<double>>(mini_end_time - mini_start_time).count() << endl;
			ip_sol_val = cplex_ip.getObjValue();
			new_pattern.resize(config.size);
			for (int i = 0; i < config.size; i++)
			{
				new_pattern[i].resize(0);
				for (int j = 0; j < patches[i].size(); j++)
				{
					if (cplex_ip.getValue(patch_vars[i][j]) > 0.99)
					{
						new_pattern[i].push_back(1);
					}
					else
					{
						new_pattern[i].push_back(0);
					}
				}
			}
			cplex_ip.end();
			ip_model.end();
		}
		
		// Add new variable to the Quadratic Problem if a new column is found
		if (price_sol_val >= target + 0.000001 || ip_sol_val >= target + 0.000001)
		{
			sarp_vars.add(IloNumVar(env, 0, 1, ILOFLOAT));
			patterns++;
			starting_set.push_back(new_pattern);
			for (int i = 0; i < config.size; i++)
			{
				for (int j = 0; j < patches[i].size(); j++)
				{
					if (new_pattern[i][j] > 0.99)
					{
						qp_cons[i][j].setLinearCoef(sarp_vars[patterns - 1], 1);
						outputtypes << "1" << "\t";
					}
					else
						outputtypes << "0" << "\t";
				}
				outputtypes << endl;
			}
		}
		else
			break;
		
		
		// Resolve
		cplex_qp.solve();
		distance = cplex_qp.getObjValue();
	}
	int sum_HH = accumulate(HH_year.begin(), HH_year.end(), 0);
	cout << "Distance = " << distance << endl;
	Jstat = distance*(double)sum_HH;

	for (int i = 0; i < pi.size(); i++)
	{
		for (int j = 0; j < pi[i].size(); j++)
		{
			nu_hat[i][j] = pi[i][j] - distance_vals[i][j];
		}
	}

	return nu_hat;
}

vector<double> QP_CG_Bootstrap(const configuration & config, const ThreeDoubleMatrix & pi_vector, const ThreeIntMatrix & patches, const vector<int> & HH_year, double Jstat)
{
	// Time variables
	chrono::high_resolution_clock::time_point start_time, end_time;
	chrono::high_resolution_clock::time_point point_start_time, point_end_time;
	chrono::high_resolution_clock::time_point mini_start_time, mini_end_time;
	start_time = chrono::high_resolution_clock::now();
	int recycled_patterns = 0;

	vector<int> variable_credits(0);

	ostringstream convert;
	convert << "Results-" << config.nr_products << "-" << config.endyear - config.startyear << "-" << config.startyear << "-" << config.endyear;
	if (config.heuristics == 1)
		convert << "-BI";
	else
		convert << "EXACT";
	if (config.UB_Threshold == 1 && config.LB_Threshold == 0)
		convert << "-UBBounds";
	else if (config.UB_Threshold == 0 && config.LB_Threshold == 1)
		convert << "-LBBounds";
	else if (config.UB_Threshold == 1 && config.LB_Threshold == 1)
		convert << "-AllBounds";
	else
		convert << "-NoBounds";
	convert << ".txt";

	ofstream output;
	string filename = convert.str();
	output.open(filename);
	ostringstream convert2;
	convert2 << "TPP-" << config.nr_products << "-" << config.endyear - config.startyear << "-" << config.startyear << "-" << config.endyear;
	if (config.heuristics == 1)
		convert2 << "-BI";
	else
		convert2 << "EXACT";
	if (config.UB_Threshold == 1 && config.LB_Threshold == 0)
		convert2 << "-UBBounds";
	else if (config.UB_Threshold == 0 && config.LB_Threshold == 1)
		convert2 << "-LBBounds";
	else if (config.UB_Threshold == 1 && config.LB_Threshold == 1)
		convert2 << "-AllBounds";
	else
		convert2 << "-NoBounds";
	convert2 << ".txt";
	ofstream timeoutput;
	filename = convert2.str();
	timeoutput.open(filename);

	
	vector<double> JStats(config.repetitions);
	int sum_HH = accumulate(HH_year.begin(), HH_year.end(), 0);

	int CPLEX_calls = 0;
	int CPLEX_patterns = 0;
	cout << "Start bootstrap Jstat Calc" << endl;

	// First we set up the Quadratic Program and the Pricing Problem.
	// The QP model.
	IloEnv env;
	NumVarMatrix slack_vars(env, config.size); for (int i = 0; i < config.size; i++) { slack_vars[i] = IloNumVarArray(env, patches[i].size(), -1, 1, ILOFLOAT); }
	ConstraintMatrix qp_cons(env, config.size);
	for (int i = 0; i < config.size; i++)
	{
		qp_cons[i] = IloRangeArray(env, patches[i].size());
		for (int j = 0; j < patches[i].size(); j++)
		{
			qp_cons[i][j] = IloRange(slack_vars[i][j] == 0); // If the current point is above the data point, the slack_var is negative.
		}
	}
	IloExpr objective(env); for (int i = 0; i < config.size; i++) { for (int j = 0; j < patches[i].size(); j++) { objective += slack_vars[i][j] * slack_vars[i][j]; } }
	IloObjective obj = IloMinimize(env, objective);
	IloModel qp_model(env); qp_model.add(obj);
	for (int i = 0; i < config.size; i++) { qp_model.add(qp_cons[i]); }
	// Adding the initial variables
	IloNumVarArray sarp_vars(env);
	int patterns = 0;

	IloCplex cplex_qp(qp_model);
	cplex_qp.setOut(env.getNullStream());

	// We set up the Pricing Problem. 
	NumVarMatrix patch_vars(env, config.size); for (int i = 0; i < config.size; i++) { patch_vars[i] = IloNumVarArray(env, patches[i].size(), 0, 1, ILOINT); }
	NumVarMatrix preference_vars(env, config.size); for (int i = 0; i < config.size; i++) { preference_vars[i] = IloNumVarArray(env, config.size, 0, 1, ILOINT); }
	
	IloObjective ip_obj = IloMaximize(env);
	IloModel ip_model = IP_Model2(patches, config, env, patch_vars, preference_vars, ip_obj);
	IloCplex cplex_ip(ip_model);
	cplex_ip.setOut(env.getNullStream());

	// A matrix used later to store the output of the QP
	TwoDoubleMatrix distance_vals;
	distance_vals.resize(patches.size()); for (int i = 0; i < patches.size(); i++) { distance_vals[i].resize(patches[i].size()); }

	ofstream outputtypes;
	outputtypes.open("zPatterns.txt");

	double time_elapsed_double = 0;
	int finished_points = 0;

	// Set up the Lower Bound problem
	IloModel lb_model(env);
	NumVarMatrix lb_point(env, config.size); for (int i = 0; i < config.size; i++) { lb_point[i] = IloNumVarArray(env, patches[i].size(), 0, 1, ILOFLOAT); }
	NumVarMatrix lb_distance(env, config.size); for (int i = 0; i < config.size; i++) { lb_distance[i] = IloNumVarArray(env, patches[i].size(), -1, 1, ILOFLOAT); }
	IloExpr lb_expr_obj(env);
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < patches[i].size(); j++)
			lb_expr_obj += lb_distance[i][j] * lb_distance[i][j];
	}
	IloObjective lb_obj = IloMinimize(env, lb_expr_obj);
	lb_model.add(lb_obj);
	ConstraintMatrix lb_cons(env, config.size);
	for (int i = 0; i < config.size; i++)
	{
		lb_cons[i] = IloRangeArray(env, patches[i].size());
		for (int j = 0; j < patches[i].size(); j++)
		{
			lb_cons[i][j] = IloRange(lb_point[i][j] - lb_distance[i][j] == 0); // RHS must be set let on.
		}
		lb_model.add(lb_cons[i]);
	}
	IloRange lb_plane_cons(env, -99, 0); // Coefficients and RHS must also be set.
	lb_model.add(lb_plane_cons);
	IloCplex lb_problem(lb_model);
	lb_problem.setOut(env.getNullStream());
	int LB_not_checked = 0;
	int LB_breaks = 0;
	double lower_bound = -99;

	// Start solving!
	for (int k = 0; k < pi_vector.size() && time_elapsed_double < config.time_limit; k++)
	{
		cout << "Bootstrap " << k << " Jstat = ";
		point_start_time = chrono::high_resolution_clock::now();
		TwoIntMatrix previous_pattern; // A pattern to check whether there are no loops that keep adding the same pattern.
		// Set the pi (with Tau-adjustment)
		TwoDoubleMatrix pi = pi_vector[k];
		for (int i = 0; i < pi.size(); i++)
		{
			for (int j = 0; j < pi[i].size(); j++)
			{
				qp_cons[i][j].setBounds(pi[i][j], pi[i][j]); // Update the QP problem.
			}
		}
		// Set the RHS for the LB problem
		if (config.LB_Threshold == 1)
		{
			for (int i = 0; i < distance_vals.size(); i++)
			{
				for (int j = 0; j < distance_vals[i].size(); j++)
				{
					lb_cons[i][j].setBounds(pi[i][j], pi[i][j]);
				}
			}
		}
		// Solve the Quadratic Model
		mini_start_time = chrono::high_resolution_clock::now();
		cplex_qp.solve();
		mini_end_time = chrono::high_resolution_clock::now();
		cout << "QP Solve = " << chrono::duration_cast<chrono::duration<double>>(mini_end_time - mini_start_time).count() << endl;
		double distance = cplex_qp.getObjValue();

		double distance_threshold;
		if (config.UB_Threshold == 0)
			distance_threshold = 0;
		else
			distance_threshold = Jstat / (double)sum_HH - 0.0000002;

		while (distance > distance_threshold + 0.0000001)
		{
			double target = 0;

			// If current distance > threshold, update the pricing problem such that the objective is the hyperplane tangent to the current solution.
			for (int i = 0; i < distance_vals.size(); i++)
			{
				for (int j = 0; j < distance_vals[i].size(); j++)
				{
					distance_vals[i][j] = cplex_qp.getValue(slack_vars[i][j]);
					ip_obj.setLinearCoef(patch_vars[i][j], distance_vals[i][j]);
					target += distance_vals[i][j] * (pi[i][j] - distance_vals[i][j]); // pi - distance_vals is current solution.
				}
				//cout << endl;
			}
			// Update the LB problem
			if(config.LB_Threshold == 1)
			{
				for (int i = 0; i < distance_vals.size(); i++)
				{
					for (int j = 0; j < distance_vals[i].size(); j++)
					{
						lb_plane_cons.setLinearCoef(lb_point[i][j], distance_vals[i][j]);
					}
				}
			}
			bool replace_var_flag = 0;
			int replace_var;
			if (config.recycle == 1)
			{
				for (int var = 0; var < sarp_vars.getSize(); var++)
				{
					double variable_value = cplex_qp.getValue(sarp_vars[var]);
					if (cplex_qp.getValue(sarp_vars[var]) <= 0.00000001)
					{
						variable_credits[var] = variable_credits[var]-1;
					}
					else
					{
						variable_credits[var] = config.recycle_creds;
					}
					if (variable_credits[var] <= 0 && replace_var_flag == 0)
					{
						replace_var_flag = 1;
						replace_var = var;
					}
				}
			}
			// Solve Pricing Problem.
			// First we use a heuristic.
			double price_sol_val;
			TwoIntMatrix new_pattern;
			if (config.heuristics == 1)
			{
				cout << "Heuristics" << endl;
				mini_start_time = chrono::high_resolution_clock::now();
				new_pattern = Best_Insertion_Heur(config, patches, distance_vals, price_sol_val);
				mini_end_time = chrono::high_resolution_clock::now();
				cout << "Heur Solve = " << chrono::duration_cast<chrono::duration<double>>(mini_end_time - mini_start_time).count() << endl;
				/*if (config.LB_Threshold == 1)
				{
					lb_plane_cons.setBounds(price_sol_val, price_sol_val);
					lb_problem.solve();
					lower_bound = lb_problem.getObjValue();
					if (lower_bound*(double)sum_HH > Jstat) // If the "heuristic" lower bound is bigger than JSTAT, we note this. If this happens often enough, we will calculate the exact lower bound.
						LB_not_checked++;
				}*/
			}
			else
			{
				price_sol_val = -99;
			}
			double ip_sol_val = 0;
			// If the heuristic does not return a good solution, we run the IP Model. We also do this if we want to check the exact lower bound.
			if (price_sol_val < target + 0.000001 || (LB_not_checked > 0 && LB_not_checked % 50 == 0))
			{
				CPLEX_calls++;
				mini_start_time = chrono::high_resolution_clock::now();
				//cplex_ip.exportModel("IPmodel.lp");
				cplex_ip.solve();
				ip_sol_val = cplex_ip.getObjValue();
				mini_end_time = chrono::high_resolution_clock::now();
				cout << "\t IP Solve = " << chrono::duration_cast<chrono::duration<double>>(mini_end_time - mini_start_time).count() << endl;
				if (ip_sol_val >= target)
					CPLEX_patterns++;
				new_pattern.resize(config.size);
				outputtypes << "IP Solution" << endl << endl;
				for (int i = 0; i < config.size; i++)
				{
					new_pattern[i].resize(0);
					for (int j = 0; j < patches[i].size(); j++)
					{
						if (cplex_ip.getValue(patch_vars[i][j]) > 0.99)
						{
							new_pattern[i].push_back(1);
						}
						else
						{
							new_pattern[i].push_back(0);
						}
					}
				}
				if (config.LB_Threshold == 1 && ip_sol_val > target + 0.000001)
				{
					//cout << "Before LB problem" << endl;
					cout << ip_sol_val << endl;
					lb_plane_cons.setUB(ip_sol_val);
					//cout << "Bounds set" << endl;
					lb_problem.solve();
					//lb_problem.exportModel("lb.lp");
					//cout << "Problem Solved" << endl;
					lower_bound = lb_problem.getObjValue();
					if (lower_bound*(double)sum_HH > Jstat) // If the exact LB is still greater than Jstat, there is no use going on with the calculation, and it is terminated.
					{
						cout << "Lower Bound Break" << endl;
						LB_breaks++;
						break;
					}
					//cout << "After LB Problem" << endl;
				}
			}

			// Add new variable to the Quadratic Problem if a new column is found
			if (price_sol_val >= target + 0.000001 || ip_sol_val >= target + 0.000001)
			{
				if (replace_var_flag == 0)
				{
					cout << "New Pattern \t";
					sarp_vars.add(IloNumVar(env, 0, 1, ILOFLOAT));
					variable_credits.push_back(config.recycle_creds);
					patterns++;
					for (int i = 0; i < config.size; i++)
					{
						for (int j = 0; j < patches[i].size(); j++)
						{
							if (new_pattern[i][j] > 0.99)
							{
								qp_cons[i][j].setLinearCoef(sarp_vars[patterns - 1], 1);
							}
						}
					}
				}
				else
				{
					cout << "Recycle Pattern \t";
					recycled_patterns++;
					for (int i = 0; i < config.size; i++)
					{
						for (int j = 0; j < patches[i].size(); j++)
						{
							if (new_pattern[i][j] > 0.99)
							{
								qp_cons[i][j].setLinearCoef(sarp_vars[replace_var], 1);
							}
							else
								qp_cons[i][j].setLinearCoef(sarp_vars[replace_var], 0);
						}
					}
				}
				cout << "Added Var" << endl;
			}
			else
			{
				cout << "No new variables" << endl;
				break;
			}

			// Resolve
			mini_start_time = chrono::high_resolution_clock::now();
			cplex_qp.solve();
			distance = cplex_qp.getObjValue();
			mini_end_time = chrono::high_resolution_clock::now();
			cout << "QP Solve = " << chrono::duration_cast<chrono::duration<double>>(mini_end_time - mini_start_time).count() << endl;
		}
		point_end_time = chrono::high_resolution_clock::now();
		timeoutput << k << ": " << chrono::duration_cast<chrono::duration<double>>(point_end_time - point_start_time).count();
		JStats[k] = distance*(double)sum_HH;
		timeoutput << "\t " << JStats[k] << endl; 
		time_elapsed_double = chrono::duration_cast<chrono::duration<double>>(point_end_time - start_time).count();
		finished_points = k + 1;
		//cout << JStats[k] << endl;
	}
	int bs_count = 0;
	for (int i = 0; i < JStats.size(); i++)
	{
		if (JStats[i] < Jstat)
			bs_count++;
	}
	end_time = chrono::high_resolution_clock::now();

	output << "Years " << config.startyear << "-" << config.endyear << endl;
	output << "Jstat = " << Jstat << endl;
	output << "p-value = " << 1 - ((double)bs_count / config.repetitions) << endl;
	output << "Computation Time = " << chrono::duration_cast<chrono::duration<double>>(end_time - start_time).count() << endl;
	output << "Patterns = " << patterns << endl;
	output << "Cplex Calls = " << CPLEX_calls << endl;
	output << "Lower Bound Breaks = " << LB_breaks << endl;
	output << "Potential Lower Bound Breaks (Heur) = " << LB_not_checked << endl;
	output << "Recycled Variables = " << recycled_patterns << endl;
	output << "Finished bootstrap points = " << finished_points << endl;
	output << endl;
	output << "Heuristics = " << config.heuristics << endl;
	output << "Recycle = " << config.recycle << endl;
	output << "Recycle Credits = " << config.recycle_creds << endl;
	output << "Upper Bound = " << config.UB_Threshold << endl;
	output << "Lower Bound = " << config.LB_Threshold << endl;
	return JStats;
}

IloModel IP_Model2(const ThreeIntMatrix & patches, const configuration & config, IloEnv & env, NumVarMatrix patch_vars, NumVarMatrix preference_vars, IloObjective ip_obj)
{
	IloModel model(env);

	// Objective Function
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < patches[i].size(); j++)
		{
			ip_obj.setLinearCoef(patch_vars[i][j], 1);
		}
	}
	model.add(ip_obj);

	// Assignment constraint (one patch on every budget is chosen)
	for (int i = 0; i < config.size; i++)
	{
		IloExpr assign_constraint(env);
		for (int j = 0; j < patches[i].size(); j++)
		{
			assign_constraint += patch_vars[i][j];
		}
		model.add(assign_constraint == 1);
		assign_constraint.end();
	}

	// Link Constraints, from choice of patches to preferences.
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < config.size; j++)
		{
			IloExpr constraint(env);
			for (int k = 0; k < patches[i].size(); k++)
			{
				if (patches[i][k][j] == 1)
					constraint += patch_vars[i][k];
			}
			constraint += -preference_vars[j][i];
			model.add(constraint <= 0),
				constraint.end();
		}
	}

	// Triangle Constraints
	for (int i = 0; i<config.size; i++)
	{
		for (int j = i + 1; j<config.size; j++)
		{
			for (int k = j + 1; k<config.size; k++)
			{
				model.add(preference_vars[i][j] + preference_vars[j][k] + preference_vars[k][i] <= 2);
				model.add(preference_vars[j][i] + preference_vars[i][k] + preference_vars[k][j] <= 2);
			}
		}
	}

	// Enforce anti-symmetry
	for (int i = 0; i<config.size; i++)
	{
		for (int j = i; j<config.size; j++)
		{
			if (i == j)
				model.add(preference_vars[i][i] == 0);
			else
				model.add(preference_vars[i][j] + preference_vars[j][i] == 1);

		}
	}

	return model;
}

IloModel IP_Model3(const ThreeIntMatrix & patches, const configuration & config, IloEnv & env, NumVarMatrix patch_vars, NumVarMatrix preference_vars, TwoDoubleMatrix distance_vals)
{
	IloModel model(env);
	// Objective Function
	IloObjective ip_obj = IloMaximize(env);
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < patches[i].size(); j++)
		{
			ip_obj.setLinearCoef(patch_vars[i][j], distance_vals[i][j]);
		}
	}
	model.add(ip_obj);

	// Assignment constraint (one patch on every budget is chosen)
	for (int i = 0; i < config.size; i++)
	{
		IloExpr assign_constraint(env);
		for (int j = 0; j < patches[i].size(); j++)
		{
			assign_constraint += patch_vars[i][j];
		}
		model.add(assign_constraint == 1);
		assign_constraint.end();
	}

	// Link Constraints, from choice of patches to preferences.
	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < config.size; j++)
		{
			IloExpr constraint(env);
			for (int k = 0; k < patches[i].size(); k++)
			{
				if (patches[i][k][j] == 1)
					constraint += patch_vars[i][k];
			}
			constraint += -preference_vars[j][i];
			model.add(constraint <= 0),
				constraint.end();
		}
	}

	// Triangle Constraints
	for (int i = 0; i<config.size; i++)
	{
		for (int j = i + 1; j<config.size; j++)
		{
			for (int k = j + 1; k<config.size; k++)
			{
				model.add(preference_vars[i][j] + preference_vars[j][k] + preference_vars[k][i] <= 2);
				model.add(preference_vars[j][i] + preference_vars[i][k] + preference_vars[k][j] <= 2);
			}
		}
	}

	// Enforce anti-symmetry
	for (int i = 0; i<config.size; i++)
	{
		for (int j = i; j<config.size; j++)
		{
			if (i == j)
				model.add(preference_vars[i][i] == 0);
			else
				model.add(preference_vars[i][j] + preference_vars[j][i] == 1);

		}
	}
	
	//IloCplex cplex2(model);
	//cplex2.exportModel("IPModel2.txt");
	return model;
}

TwoIntMatrix Best_Insertion_Heur(const configuration & config, const ThreeIntMatrix & patches, const TwoDoubleMatrix & distance_vals, double & price_obj)
{
	/*for (int i = 0; i < distance_vals.size(); i++)
	{
		for (int j = 0; j < distance_vals[i].size(); j++)
		{
			cout << distance_vals[i][j] << "\t";
		}
		cout << endl;
	}*/
	
	
	TwoIntMatrix random_orders(10); // The random orders will determine in which order the budgets will be added in the best insertion.
	for (int order = 0; order < 10; order++) // Make the random orders
	{
		for (int i = 0; i < patches.size(); i++)
		{
			random_orders[order].push_back(i);
		}
		random_shuffle(random_orders[order].begin(), random_orders[order].end());
		/*cout << "Order " << order << endl;
		for (int i = 0; i < patches.size(); i++)
		{
			cout << random_orders[order][i] << "\t";
		}
		cout << endl;*/
	} 


	// The best solutions so far.
	TwoIntMatrix best_patches_solution(patches.size()); // Saves which patches are chosen in the best solution.
	for (int i = 0; i < patches.size(); i++) { best_patches_solution[i].resize(patches[i].size()); }
	double best_sol_val = -99;

	for (int order = 0; order < 10; order++)
	{
		/*cout << "Next order" << endl;
		for (int i = 0; i < patches.size(); i++)
		{
			cout << random_orders[order][i] << "\t";
		}*/
		vector<int> Incumbentorder(patches.size(),-1);
		double Incumbentordervalue = -99;
		vector<double> Incumbentorderindvalue(patches.size(), -99);
		vector<int> Incumbentpatches(patches.size(), -1);
		

		Incumbentorder[0] = random_orders[order][0];
		
		for (int patch = 0; patch < patches[Incumbentorder[0]].size(); patch++)
		{
			if (distance_vals[Incumbentorder[0]][patch] > Incumbentorderindvalue[0])
			{
				Incumbentorderindvalue[0] = distance_vals[Incumbentorder[0]][patch];
				Incumbentpatches[0] = patch;
			}
		}
		Incumbentordervalue = Incumbentorderindvalue[0];

		for (int budget = 1; budget < random_orders[order].size(); budget++) // Loop for adding each budget to the order.
		{
			// The NewIncumb is the incumbent at this stage in the insertion. Can be changed if a better insertion point is found.
			// The normal Incumb gives the best solution after the insertion of the previous budget.
			vector<int> NewIncumbentorder(patches.size(), -1);
			double NewIncumbentordervalue = -99;
			vector<double> NewIncumbentorderindvalue(patches.size(), -99);
			vector<int> NewIncumbentpatches(patches.size(), -1);
			
			// The loop for different insertion locations.
			for(int insertloc = 0; insertloc <= budget; insertloc++)
			{ 
				//cout << "Inserting " << random_orders[order][budget] << " at position " << insertloc << endl;

				// The Temporder is for the current insertion point.
				vector<int> Temporder(patches.size(), -1);
				double Tempordervalue = -99;
				vector<double> Temporderindvalue(patches.size(), -99);
				vector<int> Temppatches(patches.size(), -1);

				// First, we look at the positions behind the insertion location. There are no changes compared to the incumbent. (No risk of adding new violations).
				for (int position = insertloc+1; position <= budget; position++)
				{
					Temporder[position] = Incumbentorder[position - 1];
					Temporderindvalue[position] = Incumbentorderindvalue[position - 1];
					Temppatches[position] = Incumbentpatches[position -1];
					/*cout << Temporder[position] << " remains at position " << position << " with patch " << Temppatches[position] << " and value "
						<< Temporderindvalue[position] << endl;*/
				}
				// Next, we look at the insertion location.
				// We go through the patches, and pick the patch with the highest value that is not below a budget further back in the ordering.
				Temporder[insertloc] = random_orders[order][budget];
				for (int patch = 0; patch < patches[random_orders[order][budget]].size(); patch++)
				{
					// Only evaluate feasibility if better value that incumbent.
					if (distance_vals[random_orders[order][budget]][patch] > Temporderindvalue[insertloc])
					{
						//cout << "Patch " << patch << " is interesting with value " << distance_vals[random_orders[order][budget]][patch] << endl;
						bool feasible = 1;
						for (int position = insertloc + 1; position <= budget; position++)
						{
							if (patches[random_orders[order][budget]][patch][Temporder[position]] == 1)
							{
								feasible = 0;
								//cout << "It is infeasible due to a conflict with " << Temporder[position] << " at position " << position << endl;
								break;
							}
						}
						if (feasible == 1)
						{
							Temppatches[insertloc] = patch;
							Temporderindvalue[insertloc] = distance_vals[random_orders[order][budget]][patch];
							//cout << "It becomes the incumbent" << endl;
						}
					}
				}
				// Finally, we look a the patches in front of the insertion location. The current incumbent patch can become infeasible.
				// First put everything correctly in the temporder before we start comparing stuff.
				for (int position = 0; position < insertloc; position++)
					Temporder[position] = Incumbentorder[position];
				// Now start the checking of the incumbent patches.
				for (int position = 0; position < insertloc; position++)
				{
					//cout << "Do we need to change " << Incumbentorder[position] << " at position " << position << "?" << endl;
					Temporder[position] = Incumbentorder[position];
					// If the incumbent patch is still feasible, then this is still the best option.
					if (patches[Incumbentorder[position]][Incumbentpatches[position]][random_orders[order][budget]] == 0)
					{
						//cout << "Still feasible" << random_orders[order][budget] << endl;
						Temporderindvalue[position] = Incumbentorderindvalue[position];
						Temppatches[position] = Incumbentpatches[position];
					}
					// If not, we again look for the best patch.
					else
					{
						/*cout << "---------------------------------------------------------------------- ISSUE!" << endl;
						cout << "Non-feasible Budget = " << Incumbentorder[position] << "\t Patch = " << Incumbentpatches[position] << endl;
						cout << "Inserted budget = " << random_orders[order][budget] << endl;
						*/
						Temporderindvalue[position] = -99;
						Temppatches[position] = -1;
						for (int patch = 0; patch < patches[Incumbentorder[position]].size(); patch++)
						{
							// Only evaluate feasibility if better value that incumbent.
							if (distance_vals[Incumbentorder[position]][patch] > Temporderindvalue[position])
							{
								//cout << "Patch " << patch << " is interesting" << endl;
								bool feasible = 1;
								for (int position2 = position + 1; position2 <= budget; position2++)
								{
									//cout << "Comparing against " << Temporder[position2] << endl;
									if (patches[Temporder[position]][patch][Temporder[position2]] == 1)
									{
										//cout << "Violation" << endl;
										feasible = 0;
										break;
									}
								}
								if (feasible == 1)
								{
									//cout << "No conflict, accepting patch" << endl;
									Temppatches[position] = patch;
									Temporderindvalue[position] = distance_vals[Incumbentorder[position]][patch];
								}
							}
						}
					}	
				}
				// We now have information for the best patches given the insert location.
				/*for (int i = 0; i < Temporder.size(); i++) cout << Temporder[i] << "\t"; cout << endl;
				for (int i = 0; i < Temporder.size(); i++) cout << Temppatches[i] << "\t"; cout << endl;
				for (int i = 0; i < Temporder.size(); i++) cout << Temporderindvalue[i] << "\t"; cout << endl;*/
				Tempordervalue = accumulate(Temporderindvalue.begin(), Temporderindvalue.begin()+budget+1, 0.0);
				//cout << Tempordervalue << endl;
				if (Tempordervalue > NewIncumbentordervalue)
				{
					NewIncumbentorder = Temporder;
					NewIncumbentorderindvalue = Temporderindvalue;
					NewIncumbentordervalue = Tempordervalue;
					NewIncumbentpatches = Temppatches;
				}
			}
			// After checking all insertion locations, we are now sure what the best one is. We move on to inserting the next budget, and update
			// the incumbent.
			Incumbentorder = NewIncumbentorder;
			Incumbentorderindvalue = NewIncumbentorderindvalue;
			Incumbentordervalue = NewIncumbentordervalue;
			Incumbentpatches = NewIncumbentpatches;
			/*cout << endl << endl << endl << endl;
			cout << "New Incumbent Order = " << endl;
			for (int i = 0; i < Incumbentorder.size(); i++) cout << Incumbentorder[i] << "\t";
			cout << endl;
			for (int i = 0; i < Incumbentorder.size(); i++) cout << Incumbentpatches[i] << "\t";
			cout << endl;
			for (int i = 0; i < budget; i++)
			{
				for (int j = i + 1; j <= budget; j++)
				{
					if (patches[Incumbentorder[i]][Incumbentpatches[i]][Incumbentorder[j]] == 1)
					{
						cout << "---------ERROR-----------------ERROR-----------------ERROR-----------ERROR-------" << endl;
						cin.get();
					}
				}
			}
			//cin.get();*/
		}
		//cout << "Testing new solution" << endl;
		if (Incumbentordervalue > best_sol_val)
		{
			// First reset the best_patches_solution, so no patches are chosen.
			for (int i = 0; i < best_patches_solution.size(); i++)
			{
				for (int j = 0; j < best_patches_solution[i].size(); j++)
				{
					best_patches_solution[i][j] = 0;
				}
			}
			
			best_sol_val = Incumbentordervalue;
			price_obj = best_sol_val;
			for (int i = 0; i < Incumbentorder.size(); i++)
			{
				best_patches_solution[Incumbentorder[i]][Incumbentpatches[i]] = 1;
			}
		}
	}
	/*for (int i = 0; i < best_patches_solution.size(); i++)
	{
		for (int j = 0; j < best_patches_solution[i].size(); j++)
			cout << best_patches_solution[i][j] << "\t";
		cout << endl;
	}*/
	return best_patches_solution;
}
