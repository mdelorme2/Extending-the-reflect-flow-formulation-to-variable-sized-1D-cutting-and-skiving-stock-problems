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
	int version = atoi(argv[4]); // related to dedicated variables for the total number of items of a certain type: without these variables (0), with these variables (1), with these variables and priority (2) 
		
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
	int cap = allo.bins[0];
	vector<bool> isNodeActive (cap + 1,false); 
	vector<bool> isNodeActiveWR (cap + 1,false); 
	isNodeActive[0] = true;
	isNodeActiveWR[0] = true;
	vector<vector<int> > arcs;

	// Item arcs (goes from 0 to cap, all item lentghs are multiplied by 2)		
	for (int i = 0; i < allo.nbItems; i++){
		for (int j = cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){			
				if (j + allo.items[i][0] * 2 <= cap){
					arcs.push_back({j, j +allo.items[i][0] * 2,i,-1}); 
					isNodeActive[j + allo.items[i][0] * 2] = true;
					isNodeActiveWR[j + allo.items[i][0] * 2] = true;
				}
				for(int l = 0; l < allo.nbBins;l++){				
					if (2*allo.bins[l] > cap && j < allo.bins[l] &&  j + allo.items[i][0] * 2 > allo.bins[l] && j <= allo.bins[l] * 2 - (j + allo.items[i][0] * 2)){
						allo.infos.nbVarS += 1;
						arcs.push_back({j, allo.bins[l]* 2 - (j + allo.items[i][0] * 2),i,l}); 
						isNodeActiveWR[allo.bins[l] * 2 - (j + allo.items[i][0] * 2)] = true;
					}
				}
			}
		}
	}
	
	// Bin arcs	
	for (int i = 0; i < allo.nbBins; i++){
		if(2*allo.bins[i] > cap){
			allo.infos.nbVarS += 1;
			arcs.push_back({allo.bins[i], allo.bins[i],-1,i});
			isNodeActiveWR[allo.bins[i]] = true;
		}
		else{
			allo.infos.nbVarS += 1;
			arcs.push_back({0, 2*allo.bins[i],-2,i});
			isNodeActiveWR[2*allo.bins[i]] = true;
		}
	}

	// Loss arcs
	vector<int> onlyActiveNodesWR;
	for (int i = 0; i <= cap; i++){
		if(isNodeActiveWR[i])
			onlyActiveNodesWR.push_back(i);
	}	
	
	for (int i = 0; i < onlyActiveNodesWR.size() - 1; i++){
		arcs.push_back({onlyActiveNodesWR[i], onlyActiveNodesWR[i+1],-1,-1});
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
		vector<GRBLinExpr> nbBinUsed(allo.nbBins,0);
		vector<GRBVar>  isArcUsed(arcs.size());
		vector<GRBLinExpr> cIn(cap + 1,0);
		vector<GRBLinExpr> cOut(cap + 1,0);
		vector<GRBLinExpr> cRIn(cap + 1,0);
		vector<GRBVar> nbItemUsed(allo.nbItems);
		vector<GRBLinExpr> nbItemUsedExpr(allo.nbItems,0);
		GRBLinExpr	nbRef = 0;
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			if (arcs[i][2] >= 0){
				isArcUsed[i] = model.addVar(0, 1, 0, GRB_INTEGER);
			}
			else{
				if(arcs[i][3] >= 0){ 
					isArcUsed[i] = model.addVar(0, 1, 0, GRB_INTEGER);
				}
				else{
					isArcUsed[i] = model.addVar(0, allo.nbBins*2, 0, GRB_INTEGER);
				}
			}
		}
		
		if(version >= 1){
			for (int i = 0; i < allo.nbItems; i++){
				nbItemUsed[i] = model.addVar(0, 1, 0, GRB_INTEGER);
			}	
		}
		
		model.update();
	
		// Perform values
		for (int i = 0; i < arcs.size(); i++){
			if(arcs[i][2] >= 0) // The arc contains an item 
				nbItemUsedExpr[arcs[i][2]] += isArcUsed[i];
			cOut[arcs[i][0]] += isArcUsed[i];
			if(arcs[i][3] == -1) // The arc is a standard arc
				cIn[arcs[i][1]] += isArcUsed[i];
			else{ // The arc is a reflected arc 
				cRIn[arcs[i][1]] += isArcUsed[i];
				nbRef += isArcUsed[i];
				nbBinUsed[arcs[i][3]] += isArcUsed[i];
			}	
		}
		
		// Flow conservation
		for (int i = 1; i <= cap; i++){
			if (isNodeActiveWR[i] == true)
				model.addConstr(cOut[i] + cRIn[i] == cIn[i]);
		}	
		model.addConstr(cOut[0] + cRIn[0] == 2 * nbRef);

		// Bin availability constraints
		for (int k = 0; k < allo.nbBins; k++){
			model.addConstr(nbBinUsed[k] <= 1); 
		}
					
		// Objective function
		if(version >= 1){
			for (int i = 0; i < allo.nbItems; i++){
				objFun += nbItemUsed[i] * allo.items[i][1];
				model.addConstr(nbItemUsed[i] <= nbItemUsedExpr[i]); 
			}	
		}			
		else{
			for (int i = 0; i < allo.nbItems; i++){
				model.addConstr(nbItemUsedExpr[i] <= 1);
				objFun += nbItemUsedExpr[i] * allo.items[i][1];
			}	
		}
		model.setObjective(objFun, GRB_MAXIMIZE); 
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		if(version >= 2){
			for (int i = 0; i < allo.nbItems; i++){
				nbItemUsed[i].set(GRB_IntAttr_BranchPriority, 5);
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
			allo.infos.opt = true;
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
		vector<vector<vector<int> > > bins(allo.nbBins);
		vector<vector<vector<int> > > arcsUsed(cap+1);
		vector<int> itemUsed(allo.nbItems,0);
		for (int i = 0; i < arcs.size(); i++){
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
			}
		}
		
		// Step 1 ==> Create half knapsacks 
		vector<vector<int> > halfBins;

		while (arcsUsed[0].size() > 0){
			int currTail = 0;
			vector<int> halfBin(2);
			halfBin[0] = 0; // the weight in the knapsack
			halfBin[1] = -1; // the type of knapsack
			while (arcsUsed[currTail].size() != 0){			
				int nextTail = arcsUsed[currTail].back()[1];
				if (arcsUsed[currTail].back()[2] >= 0 && itemUsed[arcsUsed[currTail].back()[2]] == 0){
					halfBin.push_back(arcsUsed[currTail].back()[2]);
					halfBin[0] += allo.items[arcsUsed[currTail].back()[2]][0];
					itemUsed[arcsUsed[currTail].back()[2]] += 1;
				}
				halfBin[1] = arcsUsed[currTail].back()[3];
				arcsUsed[currTail].pop_back();
				currTail = nextTail;
				if(halfBin[1] >= 0){
					halfBin[0] += cap - allo.bins[halfBin[1]];
					break;
				}
			}
			if(currTail != 0) halfBins.push_back(halfBin);
			else{
				vector<int> bin;
				//cout << "Selfsufficient bin ";
				//for(int j = 2; j < halfBin.size(); j++) cout << halfBin[j] << " "; 
				for(int j = 2; j < halfBin.size(); j++) bin.push_back(halfBin[j]); 
				bins[halfBin[1]].push_back(bin);
			}
		}

		// Step 2 ==> Order the half bins and match the largest with the smallest
		sort(halfBins.begin(), halfBins.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
			if (a[1] == -1 && b[1] != -1) return false;
			if (b[1] == -1 && a[1] != -1) return true;
			return a[0] > b[0];
		});
		int firstHalfIndex = 0; int secondHalfIndex = halfBins.size()-1;
		while (firstHalfIndex <= secondHalfIndex){
			vector<int> bin;
			/*cout << firstHalfIndex << " ";
			for(int j = 0; j < halfBins[firstHalfIndex].size(); j++){
				cout << halfBins[firstHalfIndex][j] << " "; 
			}*/
			for(int j = 2; j < halfBins[firstHalfIndex].size(); j++) 
				bin.push_back(halfBins[firstHalfIndex][j]); 
			/*cout << " matched with " << secondHalfIndex << " ";
			for(int j = 0; j < halfBins[secondHalfIndex].size(); j++){
				cout <<  halfBins[secondHalfIndex][j]; 
			}
			cout << endl;*/
			for(int j = 2; j < halfBins[secondHalfIndex].size(); j++) 
				bin.push_back(halfBins[secondHalfIndex][j]); 
			secondHalfIndex--;
			bins[halfBins[firstHalfIndex][1]].push_back(bin);
			firstHalfIndex++;		
		}

		// Check the correctness and print solution
		allo.infos.correct = 1;
		int cost = 0;
		vector<int> itemTaken (allo.nbItems,0);
		vector<int> binTaken (allo.nbBins,0);
		vector<vector<int> > binWeight (bins.size()); 
		for(int i = 0; i < bins.size(); i++){
			binWeight[i].resize(bins[i].size(),0);
			for (int j = 0; j < bins[i].size(); j++){
				cout << "Knapsack " << i << " nb " << j << ": ";
				binTaken[i] += 1;
				for (int k = 0; k < bins[i][j].size(); k++){
					binWeight[i][j] += allo.items[bins[i][j][k]][0];
					cost += allo.items[bins[i][j][k]][1];
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
		// Number of knapsacks
		for(int i = 0; i < allo.nbBins; i++){
			if(binTaken[i] > 1){
				allo.infos.correct = 0;
				cout << "Knapsack " << i << " used " << binTaken[i] << "/1" << endl;
			}
		}
		// Number of items
		for (int k = 0; k < allo.nbItems; k++){
			if(itemTaken[k] > 1){
				allo.infos.correct = 0;
				cout << "Item " << k << " packed " << itemTaken[k] << "/1" << endl;
			}
		}
		// Knapsack weight
		for(int i = 0; i < binWeight.size(); i++){
			for (int j = 0; j < binWeight[i].size(); j++){
				if(binWeight[i][j] > allo.bins[i]){
					allo.infos.correct = 0;
					cout << "Knapsack " << i << " index " << j << " has capacity " <<  binWeight[i][j] << "/" << allo.bins[i] << endl;
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
