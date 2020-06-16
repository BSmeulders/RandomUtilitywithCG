#pragma once
#include <vector>
#include <ilcplex/ilocplex.h>

using namespace std;

#ifndef Typedef_H
#define Typedef_H

typedef vector< vector<int>> TwoIntMatrix;
typedef vector< vector< vector<int>>> ThreeIntMatrix;
typedef vector< vector<double>> TwoDoubleMatrix;
typedef vector< vector < vector<double>>> ThreeDoubleMatrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloRangeArray> ConstraintMatrix;

#endif