#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>
//#include <ppl.h>
#include <Eigen/Dense>
// #include <boost/math/distributions.hpp>
#include <cassert>
#include <stack>
using namespace std;

#include "Typedef.h"
#include "Structures.h"

void Endogenous_Stastic(datastructure &data, const configuration &config);
TwoDoubleMatrix Endogenous_PiHat(const configuration & config, const ThreeDoubleMatrix & shares, const TwoDoubleMatrix & income, 
	const vector<double> &median_log_income, const TwoDoubleMatrix & instrument, const datastructure &data);
bbdatastructure Endogenous_Randomize(const configuration & config, const datastructure & data);

// Subfunctions of Endogenous_PiHat.
double compute_epsilon_tilde(const configuration &config, const TwoDoubleMatrix &r_year,const TwoDoubleMatrix &R_Hat_year, const vector<int> &income_comparison_year_HH, int HH_number);
ThreeDoubleMatrix compute_alpha(const configuration & config, const datastructure &data, const TwoIntMatrix &HH_patch, const ThreeDoubleMatrix &s, const ThreeDoubleMatrix &S_basis);
double compute_tau(const configuration & config, const datastructure &data, const vector<double> & median_log_income);

// Calling Solver
TwoDoubleMatrix QP_CG_Untightened(const configuration & config, const TwoDoubleMatrix & pi, ThreeIntMatrix & starting_set, const ThreeIntMatrix & patches, const vector<int> & HH_year, double & Jstat);
TwoDoubleMatrix Tighten(const configuration &config, const TwoDoubleMatrix &pi_hat, const ThreeIntMatrix &tau_set, double tau);
vector<double> QP_CG_Bootstrap(const configuration & config, const ThreeDoubleMatrix & pi_vector, const ThreeIntMatrix & patches, const vector<int> & HH_year, double JStat);


// Misc functions
vector<double> vector_median(TwoDoubleMatrix doublevector);
double vector_median(vector<double> vector); // In QPColumn2.cpp
vector<vector<double>> variance_pi_hat(const configuration & config, const ThreeDoubleMatrix & bs_pi_hat); // In Misc.cpp
ThreeIntMatrix Generate_Tau_Set(const configuration & config, const ThreeIntMatrix &patches);
TwoDoubleMatrix Tighten(const configuration& config, const TwoDoubleMatrix& pi_hat, const ThreeIntMatrix& tau_set, double tau);

// Output
void Endogenous_output(const configuration & config, double Jstat, double p);

void Endogenous_Stastic(datastructure & data, const configuration & config)
{
	ThreeIntMatrix pattern_set;
	TwoDoubleMatrix log_income(config.size);
	ThreeIntMatrix tau_set = Generate_Tau_Set(config, data.patches); // Pattern set will be added to later on, set of patterns for tau adjustment is fixed.
	vector<double> median_log_income(config.size);
	srand(10000 * (config.repetitions / 1000) + 1000 * config.poly_degree + 100 * (config.startyear * 10 - config.endyear) + 10 * config.nr_products);

	// Get the median_log_income.
	for (int i = 0; i < config.size; i++)
	{
		log_income[i].resize(data.HH_year[i]);
		for (int j = 0; j < log_income[i].size(); j++)
		{
			log_income[i][j] = log(data.income[i][j]);
		}
		median_log_income[i] = vector_median(log_income[i]);
	}
	// Debug
	/*cout << "Median Log Incomes" << endl;
	for (int i = 0; i < median_log_income.size(); i++)
		cout << median_log_income[i] << "\t";
	system("PAUSE");*/
	// Create PiHat.
	TwoDoubleMatrix pi_hat = Endogenous_PiHat(config, data.shares, data.income, median_log_income, data.instrument, data);

	// Time variables
	time_t start_time;
	time_t end_time;
	double time_elapsed;

	// Time Measurement
	time(&start_time);

	// Solve Untightened (
	double Jstat = 0;
	cout << "QP Untightened Solve" << endl;
	TwoDoubleMatrix nu_hat = QP_CG_Untightened(config, pi_hat, pattern_set, data.patches, data.HH_year, Jstat);
	cout << "Jstat = " << Jstat << endl;

	time(&end_time);
	time_elapsed = difftime(end_time, start_time);
	cout << "Time = " << time << endl;
	// Compute Tau
	double tau = compute_tau(config, data, median_log_income);
	cout << "Tau = " << tau << endl;
	TwoDoubleMatrix pi_hat_tight = Tighten(config, pi_hat, tau_set, tau);
	double Jstat2 = 0;
	cout << "QP Tightened Solve" << endl;
	TwoDoubleMatrix nu_tight = QP_CG_Untightened(config, pi_hat_tight, pattern_set, data.patches, data.HH_year, Jstat2);

	// Bootstrap
	ThreeDoubleMatrix bs_pi_hat(config.repetitions);
//#pragma omp parallel
	{
//#pragma omp for
		for (int i = 0; i < config.repetitions; i++)
		{
			bbdatastructure bbstruct = Endogenous_Randomize(config, data);
			bs_pi_hat[i] = Endogenous_PiHat(config, bbstruct.shares, bbstruct.income, median_log_income, bbstruct.instrument, data);
			for (int j = 0; j < bs_pi_hat[i].size(); j++)
			{
				for (int k = 0; k < bs_pi_hat[i][j].size(); k++)
				{
					bs_pi_hat[i][j][k] = bs_pi_hat[i][j][k] - pi_hat[j][k] + nu_tight[j][k];
				}
			}
		}
	}
	cout << "After Bootstrap" << endl;
	vector<double> bs_Jstat;

	bs_Jstat = QP_CG_Bootstrap(config, bs_pi_hat, data.patches, data.HH_year, Jstat);
		
	cout << "Jstat" << endl;
	cout << "bs_Jstat size = " << bs_Jstat.size() << endl;
	int bs_count = 0;
	for (int i = 0; i < bs_Jstat.size(); i++)
	{
		//cout << bs_Jstat[i] << endl;
		if (bs_Jstat[i] > Jstat)
			bs_count++;
	}
	double p = 1 - (bs_count / config.repetitions);
	cout << "Count" << bs_count << endl;
		
	Endogenous_output(config, Jstat, p);
}

TwoDoubleMatrix Endogenous_PiHat(const configuration & config, const ThreeDoubleMatrix & shares, const TwoDoubleMatrix & income,
	const vector<double> &median_log_income, const TwoDoubleMatrix & instrument, const datastructure &data)
{
	// Construct the basis function.
	ThreeDoubleMatrix r(config.size); // Basis function for instrument [i,j,k] with i = year, j = HH, k = poly degree.
	for (int i = 0; i < config.size; i++)
	{
		r[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			r[i][j].resize(config.poly_degree + 1);
			for (int k = 0; k < config.poly_degree + 1; k++)
			{
				if (k == 0)
					r[i][j][k] = 1;
				else
					r[i][j][k] = r[i][j][k - 1] * instrument[i][j];
			}
		}
	}
	
	// Variance of basis function.
	ThreeDoubleMatrix R_Hat(config.size);
	for (int i = 0; i < config.size; i++)
	{
		R_Hat[i].resize(config.poly_degree + 1);
		for (int j = 0; j < config.poly_degree + 1; j++)
		{
			R_Hat[i][j].resize(config.poly_degree + 1);
			for (int k = 0; k < config.poly_degree + 1; k++)
			{
				for (int l = 0; l < r[i].size(); l++)
				{
					R_Hat[i][j][k] += r[i][l][j] * r[i][l][k];
				}
			}
		}
	}
	
	ThreeIntMatrix income_comparison(config.size);
	for (int i = 0; i < config.size; i++)
	{
		income_comparison[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			income_comparison[i][j].resize(data.HH_year[i]);
			for (int k = 0; k < data.HH_year[i]; k++)
			{
				if (income[i][j] >= income[i][k])
					income_comparison[i][j][k] = 1;
				else
					income_comparison[i][j][k] = 0;
			}
		}
	}

	TwoDoubleMatrix epsilon_tilde(config.size);
	for (int i = 0; i < config.size; i++)
	{
		epsilon_tilde[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			epsilon_tilde[i][j] = compute_epsilon_tilde(config, r[i], R_Hat[i], income_comparison[i][j], j);
		}
	}

	// Truncation.
	vector<double> v_N(config.size);
	TwoDoubleMatrix iota_U(config.size);
	TwoDoubleMatrix iota_L(config.size);
	for (int i = 0; i < config.size; i++)
	{
		v_N[i] = 1 / (pow(data.HH_year[i], 1.0 / 3));
		iota_U[i].resize(data.HH_year[i]); iota_L[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			iota_U[i][j] = pow((1 - epsilon_tilde[i][j] + v_N[i]), 2) / (4 * v_N[i]);
			iota_L[i][j] = pow((epsilon_tilde[i][j] + v_N[i]), 2) / (4 * v_N[i]);
		}
	}
	

	for (int i = 0; i < config.size; i++)
	{
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			if (epsilon_tilde[i][j] > (1 + v_N[i]))
			{
				epsilon_tilde[i][j] = 1;
			}	
			else if (epsilon_tilde[i][j] > (1 - v_N[i]) && epsilon_tilde[i][j] <= (1 + v_N[i]))
			{
				epsilon_tilde[i][j] = 1 - iota_U[i][j];
			}
			else if (epsilon_tilde[i][j] > (-v_N[i]) && epsilon_tilde[i][j] <= (v_N[i]))
			{
				epsilon_tilde[i][j] = iota_L[i][j];
			}	
			else if (epsilon_tilde[i][j] < -v_N[i])
			{
				epsilon_tilde[i][j] = 0;
			}
		}
	}

	// Compute alpha
	// Compute implied quantities of goods
	ThreeDoubleMatrix quantity(config.size);
	for (int i = 0; i < config.size; i++)
	{
		quantity[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			quantity[i][j].resize(config.nr_products);
			for (int k = 0; k < config.nr_products; k++)
			{
				quantity[i][j][k] = shares[i][j][k] / data.prices[i][k];
			}
		}
	}

	// Check on which patch each person is given the quantities.
	TwoIntMatrix HH_patch(config.size); // First coordinate is the year, second the HH. Value is the number of the patch the HH is on (0,1,..)
	TwoIntMatrix nr_on_patch(config.size); // The number observations on the patch.
	for (int i = 0; i < config.size; i++) // Year
	{
		HH_patch[i].resize(data.HH_year[i]);
		nr_on_patch[i].resize(data.patches[i].size());
		for (int j = 0; j < data.HH_year[i]; j++)	// Household
		{
			vector<int> budget_relation(config.size);
			for (int k = 0; k < config.size; k++)		// Year to compare against.
			{
				if (i == k)
				{
					budget_relation[k] = 0;
				}
				else
				{
					double pq = 0;
					for (int l = 0; l < config.nr_products; l++) // Goods
					{
						pq += quantity[i][j][l] * data.prices[k][l]; // Compute the cost of bundle in period l.
					}
					if (pq <= 1) // Check to see whether bundle is below or above budget line l.
						budget_relation[k] = 1;
					else
						budget_relation[k] = 0;
				}
			}
			// Check on which patch the bundle is.
			for (int k2 = 0; k2 < data.patches[i].size(); k2++)
			{
				if (budget_relation == data.patches[i][k2])
				{
					HH_patch[i][j] = k2; // If patch is found, save it and break.
					nr_on_patch[i][k2]++;
					break;
				}
			}
			
		}
	}

	// Construct a second basis function.
	ThreeDoubleMatrix  s(config.size); // Basis function with i = year, j = HH, k = position.
	for (int i = 0; i < config.size; i++)
	{
		s[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			for (int k = 0; k < config.poly_degree + 1; k++)
			{
				for (int l = 0; l < config.poly_degree + 1; l++)
				{
					if ((k + l) <= config.poly_degree)
					{
						s[i][j].push_back(pow(log(income[i][j]), k)*pow(epsilon_tilde[i][j], l));
					}
				}
			}
		}
	}

	ThreeDoubleMatrix S_basis(config.size);
	for (int i = 0; i < config.size; i++)
	{
		S_basis[i].resize(s[i][0].size());
		for (int a = 0; a < S_basis[i].size(); a++)
		{
			S_basis[i][a].resize(S_basis[i].size());
			for (int b = 0; b < S_basis[i][a].size(); b++)
			{
				for (int k = 0; k < s[i].size(); k++)
				{
					S_basis[i][a][b] += s[i][k][a] * s[i][k][b];
				}
			}
		}
	}

	// Compute alpha
	ThreeDoubleMatrix alpha = compute_alpha(config, data, HH_patch, s, S_basis);

	// Construct basis function.
	TwoDoubleMatrix d(config.size);
	for (int i = 0; i < config.size; i++)
	{
		for (int k = 0; k < config.poly_degree + 1; k++)
		{
			for (int l = 0; l < config.poly_degree + 1; l++)
			{
				if ((k + l) <= config.poly_degree)
				{
					d[i].push_back(pow(median_log_income[i], k) / (l + 1));
				}
			}
		}
	}

	// Compute pi_hat.
	TwoDoubleMatrix pi_hat(config.size); // Year & Patch
	for (int i = 0; i < config.size; i++)
	{
		pi_hat[i].resize(data.patches[i].size());
		for (int j = 0; j < data.patches[i].size(); j++)
		{
			for (int k = 0; k < d[i].size(); k++)
			{
				pi_hat[i][j] += d[i][k] * alpha[i][j][k];
			}
			if (pi_hat[i][j] > 1)
				pi_hat[i][j] = 1;
			if (pi_hat[i][j] < 0)
				pi_hat[i][j] = 0;
		}
	}

	return pi_hat;
}

bbdatastructure Endogenous_Randomize(const configuration & config, const datastructure & data)
{
	bbdatastructure bootstrap;
	bootstrap.income.resize(config.size);
	bootstrap.instrument.resize(config.size);
	bootstrap.shares.resize(config.size);
	

	for (int i = 0; i < config.size; i++)
	{
		bootstrap.income[i].resize(data.HH_year[i]);
		bootstrap.instrument[i].resize(data.HH_year[i]);
		bootstrap.shares[i].resize(data.HH_year[i]);
		for (int j = 0; j < data.HH_year[i]; j++)
		{
			int household = rand() % data.HH_year[i];
			bootstrap.income[i][j] = data.income[i][household];
			bootstrap.instrument[i][j] = data.instrument[i][household];
			bootstrap.shares[i][j] = data.shares[i][household];
		}
	}

	return bootstrap;
}

double compute_epsilon_tilde(const configuration &config, const TwoDoubleMatrix &r_year, const TwoDoubleMatrix &R_Hat_year, const vector<int> &income_comparison_year_HH, int HH_number)
{
	// Debug
	bool debug = 0;
	bool allow_debug = 0;
	if (r_year.size() == 1033 && HH_number == 0 && allow_debug == 1)
		debug = 1;
	
	Eigen::MatrixXd r_matrix(r_year.size(), config.poly_degree + 1);
	for (int i = 0; i < r_year.size(); i++)
	{
		for (int j = 0; j < config.poly_degree + 1; j++)
		{
			r_matrix(i, j) = r_year[i][j];
		}
	}
	if (debug == 1)
	{
		cout << r_matrix << endl;
		system("PAUSE");
	}
		
	Eigen::MatrixXd R_Hat_matrix(config.poly_degree + 1, config.poly_degree + 1);
	for (int i = 0; i < config.poly_degree + 1; i++)
	{
		for (int j = 0; j < config.poly_degree + 1; j++)
		{
			R_Hat_matrix(i, j) = R_Hat_year[i][j];
		}
	}
	if (debug == 1)
	{
		cout << "R_Hat_Matrix" << endl;
		cout << R_Hat_matrix << endl;
		system("PAUSE");
	}
	Eigen::MatrixXd income_comparison_year_HH_matrix(income_comparison_year_HH.size(), 1);
	for (int i = 0; i < r_year.size(); i++)
	{
		income_comparison_year_HH_matrix(i) = income_comparison_year_HH[i];
	}
	
	Eigen::MatrixXd r_matrix_trans = r_matrix.transpose();
	Eigen::MatrixXd r_x_income_comp = r_matrix_trans*income_comparison_year_HH_matrix;
	if (debug == 1)
	{
		cout << r_x_income_comp << endl << endl;
		system("PAUSE");
	}
	
	Eigen::MatrixXd R_Hat_Inverse = R_Hat_matrix.inverse();
	Eigen::VectorXd x = R_Hat_Inverse*r_x_income_comp;

		
	Eigen::MatrixXd r_HH_row(1, config.poly_degree + 1);
	for (int i = 0; i < config.poly_degree + 1; i++) r_HH_row(i) = r_year[HH_number][i];
	Eigen::MatrixXd result = r_HH_row*x;
	double epsilon_tilde = result(0,0);
	
	return epsilon_tilde;
}
ThreeDoubleMatrix compute_alpha(const configuration & config, const datastructure &data, const TwoIntMatrix &HH_patch, const ThreeDoubleMatrix & s, const ThreeDoubleMatrix & S_basis)
{
	ThreeDoubleMatrix spatch(config.size);
	for (int i = 0; i < config.size; i++)
	{
		spatch[i].resize(data.patches[i].size());
		for (int j = 0; j < data.patches[i].size(); j++) // Resize all year-patch combinations.
		{
			spatch[i][j].resize(s[i][0].size(),0);
		}
		for (int j = 0; j < data.HH_year[i]; j++) // Go through all HH to sum.
		{
			for (int k = 0; k < s[i][j].size(); k++)
			{
				spatch[i][HH_patch[i][j]][k] += s[i][j][k];
			}
		}
	}	
	Eigen::MatrixXd S_Matrix(S_basis[0].size(), S_basis[0].size());
	Eigen::MatrixXd spatch_vector(S_basis[0].size(),1);
	Eigen::VectorXd alpha_vector(S_basis[0].size());
	ThreeDoubleMatrix alpha(config.size);

	for (int i = 0; i < config.size; i++)
	{
		alpha[i].resize(data.patches[i].size());
		// Set up the values for S_Matrix for this year.
		for (int j = 0; j < S_basis[i].size(); j++)
		{
			for (int k = 0; k < S_basis[i].size(); k++)
			{
				S_Matrix(j, k) = S_basis[i][j][k];
			}
		}
		// Go on to the computation for each patch in the year.
		for (int j = 0; j < data.patches[i].size(); j++)
		{
			alpha[i][j].resize(S_basis[0].size());
			if (spatch[i][j][0] != 0.0) // If equal to 0, each entry in alpha will also be 0, no computation necessary.
			{
				// Set up te spatch_vector for the patch.
				for(int k = 0; k < S_basis[0].size(); k++)
				{
					spatch_vector(k) = spatch[i][j][k];
				}
				// Solve for alpha
				alpha_vector = S_Matrix.colPivHouseholderQr().solve(spatch_vector);
				for (int k = 0; k < S_basis[0].size(); k++)
				{
					alpha[i][j][k] = alpha_vector(k);
				}
			}
		}
	}


	return alpha;
}

double compute_tau(const configuration & config, const datastructure & data, const vector<double>& median_log_income)
{
	ThreeDoubleMatrix bs_pi_hat(config.repetitions);
	vector<vector<double>> bs_pi_hat_var;
	
//#pragma omp parallel
	{
//#pragma omp for
		for (int i = 0; i < config.repetitions; i++)
		{
			bbdatastructure bbstruct = Endogenous_Randomize(config, data);
			bs_pi_hat[i] = Endogenous_PiHat(config, bbstruct.shares, bbstruct.income, median_log_income, bbstruct.instrument, data);
		}
	}
	bs_pi_hat_var = variance_pi_hat(config, bs_pi_hat);
	vector<double> n_K(config.size);
	double min_n_K;
	for (int i = 0; i < config.size; i++)
	{
		double sum_v_N = 0;
		for (int j = 0; j < bs_pi_hat_var[i].size(); j++)
		{
			sum_v_N += (bs_pi_hat_var[i][j] * data.HH_year[i]);
		}
		n_K[i] = data.HH_year[i] * data.patches[i].size() / sum_v_N;
		if (i == 0)
			min_n_K = n_K[i];
		else if (n_K[i] < min_n_K)
			min_n_K = n_K[i];
	}
	double tau = sqrt(log(min_n_K) / min_n_K);
	return tau;
}

void Endogenous_output(const configuration & config, double Jstat, double p)
{
	ostringstream convertp;

	convertp << "L:" << config.endyear - config.startyear << ":" << config.startyear << " - " << config.endyear << ".csv";
	string filename = convertp.str(); convertp.str("");
	ofstream output;
	output.open(filename);
	output << Jstat << "," << p;
}

TwoDoubleMatrix Tighten(const configuration& config, const TwoDoubleMatrix& pi_hat, const ThreeIntMatrix& tau_set, double tau)
{
	TwoDoubleMatrix tau_set_sum(config.size);
	for (int i = 0; i < config.size; i++)
	{
		tau_set_sum[i].resize(pi_hat[i].size(), 0);
		for (int j = 0; j < tau_set_sum[i].size(); j++)
		{
			for (int set_element = 0; set_element < tau_set.size(); set_element++)
			{
				tau_set_sum[i][j] += tau_set[set_element][i][j];
			}
			tau_set_sum[i][j] = tau_set_sum[i][j] / tau_set.size();
		}
	}
	TwoDoubleMatrix pi_hat_tight(config.size);
	for (int i = 0; i < pi_hat_tight.size(); i++)
	{
		pi_hat_tight[i].resize(pi_hat[i].size());
		for (int j = 0; j < pi_hat_tight[i].size(); j++)
		{
			pi_hat_tight[i][j] = pi_hat[i][j] - tau_set_sum[i][j] * tau;
		}
	}
	return pi_hat_tight;
}