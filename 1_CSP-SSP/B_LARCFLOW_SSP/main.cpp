#include "main.h"

/*	*************************************************************************************
	*************************************  MAIN *****************************************
	************************************************************************************* */

double initTimeModelCPU;

int main(int argc, char **argv){
	
	initTimeModelCPU = getCPUTime();
	
	// local variables
	Allocation allo ;
	string filein = argv[2];
	string path = argv[1];
	string pathAndFileout = argv[3];
	int version = atoi(argv[4]); // related to dedicated variables for the total number of bins: without this variable (0), with this variable (1), with this variable and priority (2) 
		
	// functions
	allo.load(path,filein);
	allo.printProb();
	allo.infos.timeCPU.push_back(0);

	ILP(allo, version);
	allo.printInfo(pathAndFileout);
}

int ILP(Allocation& allo, const int& version){

	// Create arcs
	vector<bool> isNodeActive (allo.cap,false);
	isNodeActive[0] = true;
	vector<vector<int> > arcs;

	// Item arcs		
	vector<int> isItemTypeShifted(allo.nbItemTypes,0); 
	for (int i = 0; i < allo.nbItemTypes; i++){
		vector<bool> done (allo.cap,false);
		for (int j = allo.cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){
				for (int k = 1; k <= allo.items[i][1]; k++){ 
					if(j + (k-1) * allo.items[i][0] >= allo.cap || done[j + (k-1) * allo.items[i][0]]) break;
					// if tail between cap - weight[i] and cap 
					if(j + (k-1) * allo.items[i][0] > allo.cap - allo.items[i][0]){
						isItemTypeShifted[i] = 1;
						if(!done[allo.cap - allo.items[i][0]]){
							arcs.push_back({allo.cap - allo.items[i][0], allo.cap, i});
							done[allo.cap - allo.items[i][0]] = true;
						}
					}
					else{						
						arcs.push_back({j + (k-1) * allo.items[i][0], j + k * allo.items[i][0],i}); 
						isNodeActive[j + k * allo.items[i][0]] = true;
					}
					done[j + (k-1) * allo.items[i][0]] = true;
				}
			}
		}
	}

	// Reversed Loss Arcs 
	int lit = allo.nbItemTypes - 1;
	for (int i = 0; i < allo.nbItemTypes; i++){
		if(isItemTypeShifted[i]){
			lit = min(lit,i);
			isNodeActive[allo.cap - allo.items[i][0]] = 1;
		}
	}
	
	vector<int> onlyActiveNodes;
	for (int i = allo.cap - allo.items[lit][0]; i < allo.cap; i++){
		if(isNodeActive[i])
			onlyActiveNodes.push_back(i);
	}	

	for (int i = 0; i < onlyActiveNodes.size() - 1; i++){
		arcs.push_back({onlyActiveNodes[i+1], onlyActiveNodes[i],-1});
	}
	
	// Print arcs
	cout << arcs.size() << " arcs " << " --- lit is " << lit << endl;
	/*for(int i=0;i<arcs.size();i++){
		for(int j=0;j<arcs[i].size();j++){
			cout << arcs[i][j] << " ";
		}
		cout << endl;
	}*/
	
	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
		
	GRBEnv env;
	env.set(GRB_DoubleParam_MemLimit, 60);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBVar nbBinUsed;
		if(version >= 1) nbBinUsed = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER); 
		vector<GRBVar>  isArcUsed(arcs.size());
		vector<GRBLinExpr> cIn(allo.cap + 1,0);
		vector<GRBLinExpr> cOut(allo.cap + 1,0);
		vector<GRBLinExpr> nbItemUsed(allo.nbItemTypes,0);
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			if (arcs[i][2] >= 0) 
				isArcUsed[i] = model.addVar(0,  allo.items[arcs[i][2]][1], 0, GRB_INTEGER);
			else 
				isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
		}
		
		model.update();
	
		// Perform values
		for (int i = 0; i < arcs.size(); i++){
			cIn[arcs[i][1]] += isArcUsed[i];
			cOut[arcs[i][0]] += isArcUsed[i];
			if (arcs[i][2] >= 0) 
				nbItemUsed[arcs[i][2]] += isArcUsed[i];
		}
		
		// Flow conservation
		for (int i = 1; i < allo.cap; i++){
			if (isNodeActive[i] == true)
				model.addConstr(cIn[i] == cOut[i]);
		}	
		model.addConstr(cOut[0] == cIn[allo.cap]);

		// Item availability constraints
		for (int k = 0; k < allo.nbItemTypes; k++){
			model.addConstr(nbItemUsed[k] <= allo.items[k][1]); 
		}
					
		// Objective function
		if(version >= 1){
			objFun += nbBinUsed;
			model.addConstr(cOut[0] >= nbBinUsed); 
		}			
		else{
			objFun += cOut[0];
		}
		model.setObjective(objFun, GRB_MAXIMIZE); 
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		if(version >= 2){
			nbBinUsed.set(GRB_IntAttr_BranchPriority, 5);
		}
		model.optimize();
		
		// Filling Info
		allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
		allo.infos.ObjBound = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		allo.infos.opt = false;

		// Get Info pre preprocessing
		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);

		// If no solution found
		if (model.get(GRB_IntAttr_SolCount) < 1){
			cout << "Failed to optimize ILP. " << endl;
			allo.infos.ObjBound  = -1;
			allo.infos.ObjVal  = 0;
			allo.infos.correct = 1;
			return -1;
		}

		// If solution found
		allo.infos.ObjVal = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);	
		if(allo.infos.ObjVal == allo.infos.ObjBound) allo.infos.opt = true;

		// Filling Solution 	
		vector<vector<int> > bins;
		vector<vector<vector<int> > > arcsUsed(allo.cap+1);
		for (int i = 0; i < arcs.size(); i++){
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
			}
		}
		while (arcsUsed[0].size() > 0){
			vector<int> bin;
			int currTail = 0;
			for (;;){
				int nextTail = arcsUsed[currTail].back()[1];
				if(arcsUsed[currTail].back()[2] >= 0)
					bin.push_back(arcsUsed[currTail].back()[2]);
				arcsUsed[currTail].pop_back();
				currTail = nextTail;
				if (currTail == allo.cap){
					break;
				}
			}
			bins.push_back(bin);
		}
		
		// Check the correctness and print solution
		allo.infos.correct = 1;
		vector<int> itemTaken (allo.nbItemTypes,0);
		vector<int> binWeight (bins.size(),0);
		for (int i = 0; i < bins.size(); i++){
			cout << "Bin " << i << ": ";
			for (int j = 0; j < bins[i].size(); j++){
				binWeight[i] += allo.items[bins[i][j]][0];
				itemTaken[bins[i][j]] += 1;
				cout << bins[i][j] << " ";
			}
			cout << endl;
		}
		// Number of bins
		if(bins.size() != allo.infos.ObjVal){
			allo.infos.correct = 0;
			cout << "Solution uses " << bins.size() << " bins and not " <<  allo.infos.ObjVal << endl;
		}
		// Number of items
		for (int k = 0; k < allo.nbItemTypes; k++){
			if(itemTaken[k] > allo.items[k][1]){
				allo.infos.correct = 0;
				cout << "Item " << k << " packed " << itemTaken[k] << "/" << allo.items[k][1] << endl;
			}
		}
		// Bin weight
		for (int i = 0; i < bins.size(); i++){
			if(binWeight[i] < allo.cap){
				allo.infos.correct = 0;
				cout << "Bin " << i << " has capacity " <<  binWeight[i] << "/" << allo.cap << endl;
			}
		}		
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	// End
	return 0;
}
