using namespace std;

#ifndef Struct_H
#define Struct_H

struct datastructure {
	ThreeDoubleMatrix shares; // The share of consumption for each class for each person. [i,j,k] with i the Period, j the HH and k the class.
	ThreeDoubleMatrix quant;
	TwoDoubleMatrix prices; // The prices of all classes in all periods [i,j] with i the period, j the class of products.
	TwoDoubleMatrix income; // The amount of money spent on products in the considered classes of each HH in period [i,j] with i the period, j the HH.
	TwoDoubleMatrix instrument; // An instrument representing the total income (Z in the K&S code) of the HH.  
	ThreeIntMatrix patches; // A 3-dimensional matrix representing the patches. The first coordinate is the number of the budget, the second the number of the patch. patches[2][3] thus refers to the fourth patch on the third budget.
							// The third coordinate again refers to a budget. The values is 0/1, with 1 if the patch is BELOW the budget. I.E. if patches[2][3][0] == 1, then the patch is BELOW budget 1 (and any patch on that budget is preferred to
							// the patch 3/4.
	TwoIntMatrix frequency; // The number of times a bundle on the given patch is observed in the dataset. frequency[2][3] gives the frequency of the fourth patch on the second budget.
	vector<int> HH_year; // The number of Households in a given year.
};

struct bbdatastructure {
	ThreeDoubleMatrix shares; // The share of consumption for each class for each person. Index is the class.
	TwoDoubleMatrix income; // The amount of money spent on products in the considered classes of each HH in period [i,j] with i the period, j the HH.
	TwoDoubleMatrix instrument; // An instrument representing the total income (Z in the K&S code) of the HH.  
};

struct DKSQdata {
	ThreeDoubleMatrix shares; // The share of consumption for each class for each person [i,j,k] with i the Period, j the HH and k the class. 
	TwoDoubleMatrix prices; // The prices of all classes in all periods [i,j] with i the period, j the class of products.
	ThreeIntMatrix patches; // A 3-dimensional matrix representing the patches. The first coordinate is the number of the budget, the second the number of the patch. patches[2][3] thus refers to the fourth patch on the third budget.
							// The third coordinate again refers to a budget. The values is 0/1, with 1 if the patch is BELOW the budget. I.E. if patches[2][3][0] == 1, then the patch is BELOW budget 1 (and any patch on that budget is preferred to
							// the patch 3/4.
	vector<int> HH_year; // The number of Households in a given year.
};

struct configuration
{
	bool QP;
	bool LP;
	int heuristics;
	int repetitions;
	int startyear;// The first year of the considered subset, two final digits.
	int endyear;// The final year of the considered subset, two final digits.
	int size; // The number of years considered in the dataset.
	int nr_products;// The number of considered classes of goods.
	int start_vars;
	int tau_vars; // The number of types generated to generate the cones origin.
	int poly_degree;
	int recycle; // Whether or not recycling of variables is done.
	int recycle_creds; // Number of times a qp variable can be unused before it is a candidate for recycling.
	int UB_Threshold;
	int LB_Threshold;
	bool type_estimation;
	int estimation_number;
	int dataset; // 1 is KS, 2 is DKSQ
	int time_limit;
};

#endif 

