README:

1. Source code.
Only source code is provided, which must be compiled before use. 

2. External dependencies.
The code relies on several external libraries.
- CPLEX Optmizer (https://www.ibm.com/analytics/cplex-optimizer).
Used for solving the master problem (Quadratic Program) and the pricing problem (IP) to optimality.
CPLEX requires a license. Free academic licenses are available. CPLEX versions 12.7 and 12.10 have been used with this code. Other versions should work as well.

- Boost libraries (https://www.boost.org/)
The boost library is used for reading the configuration file.

- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
Used for matrix calculations.

3. Data
The folder "input" contains all input files (prices and quantities) in the right format for the program to read. The "input" folder should be placed in the same folder as the compiled executable. 

4. Configuration file
The configuration file can be used to set the period to be analyzed, number of goods, solver options, etc... 
Any non-comment lines in the configuration file that do not correspond to a program options can make crash the program. Please make minimal changes inside this file.

5. Batch file
A batch file is provided. It calls the compiled executable and passes the configuration file as an argument. This file should be in the same folder as the executable, configuration file and input folder.
