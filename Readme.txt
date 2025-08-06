This repository contains the code for all algorithms discussed in the paper "Extending the reflect flow formulation to variable-sized one-dimensional cutting and skiving stock problems" by Maxence Delorme and John Martinovic. 

Our algorithms are coded in C++ and use the commercial solver Gurobi for the ILP models. 
This repository contains the 30 solution methods tested in Section 7. 
In folder CSP-SSP, 2 solution methods for the CSP, 7 for the SSP (with suffix "_SSP"), and the instance files.
In folder CSP-SSP, 5 solution methods for the VSCSP, 11 for the VSSSP (with suffix "_SSP"), and the instance files.
In folder MKP, 5 solution methods for the MKP and the instance files.

The different folders correspond to the following methods in our paper:
- ARCFLOW					| The arcflow formulation
- LARCFLOW					| The arcflow formulation with reversed loss arcs (only for the SSP and VSSSP)
- REFLECT					| The reflect formulation 
- REFLECTB 					| The reflect formulation with backward loss arcs (only for the SSP and VSSSP)
- REFLECTBSB				| The reflect formulation with backward loss arcs and the adapted reduction procedure (only for the SSP and VSSSP)
- REFLECTBTVAR				| The reflect formulation with backward loss arcs and conversion variables (only for the SSP and VSSSP)
- REFLECTBSBTVAR			| The reflect formulation with backward loss arcs, the adapted reduction procedure, and conversion variables (only for the SSP and VSSSP)
- REFLECTNOSMALLB			| The reflect formulation with small-roll arcs (only for the VSCSP, MKP, and VSSSP)
- REFLECTNOSMALLBTRANSA		| The reflect formulation with small-roll and extended reflected connection arcs (only for the VSCSP, MKP, and VSSSP)
- REFLECTNOSMALLBTRANSADA	| The reflect formulation with small-roll, extended reflected connection, and downgrade arcs (only for the VSCSP, MKP, and VSSSP)
- REFLECTBFULL				| The reflect formulation with backward loss arcs, the adapted reduction procedure, conversion variables, small-roll, extended reflected connection, and downgrade arcs (only for the VSSSP)

Each folder contains the same substructure. For example, 1_CYCLE contains the following files:
- Allocation.cpp	| Contains a number of secondary functions (this file is usually the same for each subfolder)
- Allocation.h		| The header file corresponding to Allocation.cpp (this file is usually the same for each subfolder)
- main.cpp			| The front-end code for using the method  
- main.h			| The header file corresponding to main.cpp 
- makefile			| Used for compiling under linux (it needs to be updated by the user)
- time.cpp			| A generic file used to measure the computation time 
- time.h			| The header file corresponding to time.cpp 

********

Once compiled, the following command can be used to run the algorithm:
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" "V" for the 30 approaches
where
- PROGRAM is the name of the compiled software 
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance
- V is the configuration used (0, 1, or 2) related to the use of dummy variables

********

## License
This project is licensed under the [CC BY-NC-ND 4.0](LICENSE.md) license.

Any questions, comments, or issues with the code can be reported to Maxence Delorme's academic email address.