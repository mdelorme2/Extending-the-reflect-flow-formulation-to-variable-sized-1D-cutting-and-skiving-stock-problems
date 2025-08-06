#ifndef ALLOCATION_H
	#define ALLOCATION_H
	
	using namespace std;
	#include <iostream> 
	#include <iomanip> 
	#include <fstream>
	#include <sstream>
	#include <vector>
	#include <string>
	#include <set>
	#include <iostream> 
	#include <math.h> 
	#include <cstdlib>
	#include <algorithm>
	
	class Info;
	class Allocation;

/*	*************************************************************************************
	************************************* INFO ******************************************
	************************************************************************************* */

	class Info{
	public:
		bool opt;
		vector<double> timeCPU;
		int ObjVal;
		int ObjBound;
		int nbCons;
		int nbVar;
		int nbNZ;
		int correct;
	};

/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */

	class Allocation{
	public:
		
		// Data read from the file
		string name;
		int nbItemTypes;
		vector<vector<int> > items;
		int cap;
		
		// Given by the ILP model
		Info infos;

		void load(const string& path, const string& filein);
		void printProb();
		void printInfo(const string& pathAndFileout);
	};
	
	
#endif 