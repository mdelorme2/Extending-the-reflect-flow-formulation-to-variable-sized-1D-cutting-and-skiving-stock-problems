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
	int version = atoi(argv[4]); // related to dedicated variables for the total number of bins of a certain type: without these variables (0), with these variables (1), with these variables and priority (2) 
		
	// functions
	allo.load(path,filein);
	allo.printProb();
	allo.infos.timeCPU.push_back(0);

	ILP(allo, version);
	allo.printInfo(pathAndFileout);
}

int ILP(Allocation& allo, const int& version){

	allo.infos.nbVarS = 0;
	allo.infos.nbConsS = 0;
	
	// Create arcs
	int cap = allo.bins[0][0];
	vector<bool> isNodeActive (cap+1,false);
	isNodeActive[0] = true;
	vector<vector<int> > arcs;

	// Item arcs		
	for (int i = 0; i < allo.nbItemTypes; i++){
		vector<bool> done (cap,false);
		for (int j = cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){
				for (int k = 1; k <= allo.items[i][1]; k++){ 
					if(j + (k-1) * allo.items[i][0] >= cap || done[j + (k-1) * allo.items[i][0]]) break; 			
					arcs.push_back({j + (k-1) * allo.items[i][0], min(cap, j + k * allo.items[i][0]),i,-1}); 
					isNodeActive[min(cap, j + k * allo.items[i][0])] = true;
					done[j + (k-1) * allo.items[i][0]] = true;
				}
			}
		}
	}
	
	// Bin arcs	
	for (int i = 0; i < allo.nbBinTypes; i++){
		arcs.push_back({allo.bins[i][0], cap + 1,-1,i});
		isNodeActive[allo.bins[i][0]] = true;
	}

	// Loss arcs
	vector<int> onlyActiveNodes;
	for (int i = allo.bins.back()[0]; i <= cap; i++){
		if(isNodeActive[i])
			onlyActiveNodes.push_back(i);
	}	
	for (int i = 0; i < onlyActiveNodes.size() - 1; i++){
		arcs.push_back({onlyActiveNodes[i+1], onlyActiveNodes[i],-1,-1});
	}
	
	// Print arcs
	cout << arcs.size() << " arcs " << endl;
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
		vector<GRBVar> nbBinUsed(allo.nbBinTypes);
		vector<GRBLinExpr> nbBinUsedExpr(allo.nbBinTypes);
		vector<GRBVar>  isArcUsed(arcs.size());
		vector<GRBLinExpr> cIn(cap + 2,0);
		vector<GRBLinExpr> cOut(cap + 2,0);
		vector<GRBLinExpr> nbItemUsed(allo.nbItemTypes,0);
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			if (arcs[i][2] >= 0){
				isArcUsed[i] = model.addVar(0, allo.items[arcs[i][2]][1], 0, GRB_INTEGER);
			}
			else{
				if(arcs[i][3] >= 0){ 
					isArcUsed[i] = model.addVar(0, allo.bins[arcs[i][3]][1], 0, GRB_INTEGER);
				}
				else{
					isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
				}
			}
		}
		
		if(version >= 1){
			for (int i = 0; i < allo.nbBinTypes; i++){
				nbBinUsed[i] = model.addVar(0, allo.bins[i][1], 0, GRB_INTEGER);
			}	
		}
		
		model.update();
	
		// Perform values
		for (int i = 0; i < arcs.size(); i++){
			cIn[arcs[i][1]] += isArcUsed[i];
			cOut[arcs[i][0]] += isArcUsed[i];
			if(arcs[i][2] >= 0) 
				nbItemUsed[arcs[i][2]] += isArcUsed[i];
			if(arcs[i][3] >= 0) 
				nbBinUsedExpr[arcs[i][3]] += isArcUsed[i];
		}
		
		// Flow conservation
		for (int i = 1; i <= cap; i++){
			if (isNodeActive[i] == true)
				model.addConstr(cIn[i] == cOut[i]);
		}	
		model.addConstr(cOut[0] == cIn[cap+1]);

		// Item availability constraints
		for (int k = 0; k < allo.nbItemTypes; k++){
			model.addConstr(nbItemUsed[k] <= allo.items[k][1]); 
		}
					
		// Objective function
		if(version >= 1){
			for (int i = 0; i < allo.nbBinTypes; i++){
				objFun += nbBinUsed[i] * allo.bins[i][2];
				model.addConstr(nbBinUsedExpr[i] >= nbBinUsed[i]); 
			}	
		}			
		else{
			for (int i = 0; i < allo.nbBinTypes; i++){
				model.addConstr(nbBinUsedExpr[i] <= allo.bins[i][1]); ;
				objFun += nbBinUsedExpr[i] * allo.bins[i][2];
			}	
		}
		model.setObjective(objFun, GRB_MAXIMIZE); 
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		if(version >= 2){
			for (int i = 0; i < allo.nbBinTypes; i++){
				nbBinUsed[i].set(GRB_IntAttr_BranchPriority, 5);
			}
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
		if(model.get(GRB_IntAttr_Status) == 2) allo.infos.ObjBound = allo.infos.ObjVal;		
		if(allo.infos.ObjVal == allo.infos.ObjBound) allo.infos.opt = true;

		// Filling Solution
		vector<vector<vector<int> > > bins(allo.nbBinTypes);
		vector<vector<vector<int> > > arcsUsed(cap+1);
		for (int i = 0; i < arcs.size(); i++){
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
			}
		}
		while (arcsUsed[0].size() > 0){
			vector<int> bin;
			int currTail = 0;
			int type = -1;
			for (;;){
				int nextTail = arcsUsed[currTail].back()[1];
				if(arcsUsed[currTail].back()[2] >= 0)
					bin.push_back(arcsUsed[currTail].back()[2]);
				type = arcsUsed[currTail].back()[3];
				arcsUsed[currTail].pop_back();
				currTail = nextTail;
				if (type >= 0){
					break;
				}
			}
			bins[type].push_back(bin);
		}
		
		// Check the correctness and print solution
		allo.infos.correct = 1;
		int cost = 0;
		vector<int> itemTaken (allo.nbItemTypes,0);
		vector<int> binTaken (allo.nbBinTypes,0);
		vector<vector<int> > binWeight (bins.size()); 
		for(int i = 0; i < bins.size(); i++){
			binWeight[i].resize(bins[i].size(),0);
			for (int j = 0; j < bins[i].size(); j++){
				cout << "Bin type " << i << " nb " << j << ": ";
				cost += allo.bins[i][2];
				binTaken[i] += 1;
				for (int k = 0; k < bins[i][j].size(); k++){
					binWeight[i][j] += allo.items[bins[i][j][k]][0];
					itemTaken[bins[i][j][k]] += 1;
					cout << bins[i][j][k] << " ";
				}
				cout << endl;
			}
		}
		// Objective value
		if(cost != allo.infos.ObjVal){
			allo.infos.correct = 0;
			cout << "Solution has value " << cost << " and not " <<  allo.infos.ObjVal << endl;
		}
		// Number of bins
		for(int i = 0; i < allo.nbBinTypes; i++){
			if(binTaken[i] > allo.bins[i][1]){
				allo.infos.correct = 0;
				cout << "Bin type " << i << " used " << binTaken[i] << "/" << allo.bins[i][1] << endl;
			}
		}
		// Number of items
		for (int k = 0; k < allo.nbItemTypes; k++){
			if(itemTaken[k] > allo.items[k][1]){
				allo.infos.correct = 0;
				cout << "Item type " << k << " packed " << itemTaken[k] << "/" << allo.items[k][1] << endl;
			}
		}
		// Bin weight
		for(int i = 0; i < binWeight.size(); i++){
			for (int j = 0; j < binWeight[i].size(); j++){
				if(binWeight[i][j] < allo.bins[i][0]){
					allo.infos.correct = 0;
					cout << "Bin type " << i << " index " << j << " has capacity " <<  binWeight[i][j] << "/" << allo.bins[i][0] << endl;
				}
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
