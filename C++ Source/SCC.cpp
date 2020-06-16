#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
// #include <boost/math/distributions.hpp>
#include <cassert>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <stack>
using namespace std;

typedef vector< vector<int>> TwoIntMatrix;

TwoIntMatrix strongly_connected_components(int Graph_Size, TwoIntMatrix RP);
void strongconnect(int v, int& index, int *v_index, int *lowlink, stack<int>& S, int& SCCindex, bool *stack_ell, int Graph_Size, TwoIntMatrix RP, vector<int> &SCC);


TwoIntMatrix strongly_connected_components(int Graph_Size, TwoIntMatrix RP)							// Using Tarjan's algorithm
{
	int vertex = 0;
	int index = 0;
	int SCCindex = 0;
	int nr_SCC = 0;
	int* lowlink = NULL;
	int* v_index = NULL;
	vector<int> SCC;
	bool* stack_ell = NULL;
	stack<int> S;

	SCC.resize(Graph_Size);
	lowlink = new int[Graph_Size];
	v_index = new int[Graph_Size];
	stack_ell = new bool[Graph_Size];


	for (index = 0; index < Graph_Size; index++)					// Initialize the values of the arrays of every vertex to 0
	{
		v_index[index] = 0;
		stack_ell[index] = 0;
		lowlink[index] = 0;
	}

	for (index = 1, vertex = 0; vertex < Graph_Size; vertex++)	// Check whether a vertex is already added to the stack (v_index defined)
	{
		if (v_index[vertex] == 0)								// 
		{
			//cout << "Strongconnect because unassigned for node " << vertex << endl;
			strongconnect(vertex, index, v_index, lowlink, S, SCCindex, stack_ell, Graph_Size, RP, SCC);
		}
	}
	delete[] lowlink;
	delete[] v_index;
	delete[] stack_ell;

	// SCCout has one row per SCC, and the elements of that row are the observations in that SCC.
	for (int i = 0; i < Graph_Size; i++)
	{
		if (SCC[i]+1 > nr_SCC)
			nr_SCC = SCC[i]+1;
	}
	TwoIntMatrix SCCout; SCCout.resize(nr_SCC);

	for (int i = 0; i < Graph_Size; i++)
	{
		SCCout[SCC[i]].push_back(i);
	}
	return SCCout;
}

void strongconnect(int v, int& index, int *v_index, int *lowlink, stack<int>& S, int& SCCindex, bool *stack_ell, int Graph_Size, TwoIntMatrix RP, vector<int> &SCC)
{
	int i;
	int pop_ver = -1;

	v_index[v] = index;
	lowlink[v] = index;
	index++;
	S.push(v);
	stack_ell[v] = 1;
	//cout << "In strongconnect for node " << v << " v_index = " << v_index[v] << " lowlink = " << lowlink[v] << endl;


	for (i = 0; i < Graph_Size; i++)
	{
		if (RP[v][i] == 1 && v != i)
		{
			if (v_index[i] == 0)
			{
				//cout << "In strongconnect for node " << v << " starting a new stronconnect for node " << i << endl;
				strongconnect(i, index, v_index, lowlink, S, SCCindex, stack_ell, Graph_Size, RP, SCC);
				lowlink[v] = min(lowlink[v], lowlink[i]);
				//cout << "After strongconnect for node " << i << " the lowlink for " << v << " has become " << lowlink[v] << endl;
			}
			else if (stack_ell[i] == 1) //Check whether i is in the stack
			{
				lowlink[v] = min(lowlink[v], v_index[i]);
				//cout << "Succesor " << i << " in stack for node " << v << " this changes lowlink to " << lowlink[v] << endl;
			}
		}
	}

	if (lowlink[v] == v_index[v])
	{
		//cout << "Node " << v << endl;
		for (; v != pop_ver;)
		{
			pop_ver = S.top();
			//cout << "Pop out " << pop_ver << endl;
			//cout << SCCindex << endl;
			SCC[pop_ver] = SCCindex;
			S.pop();
			stack_ell[pop_ver] = 0;
		}
		SCCindex++;
	}
}