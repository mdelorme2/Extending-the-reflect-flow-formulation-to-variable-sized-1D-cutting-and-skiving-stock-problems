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
	vector<bool> isNodeActive (cap + 1,false); 
	vector<bool> isNodeActiveWR (cap + 1,false); 
	isNodeActive[0] = true;
	isNodeActiveWR[0] = true;
	vector<vector<int> > arcs;

	// Item arcs (goes from 0 to cap, all item lentghs are multiplied by 2)		
	vector<int> minVerR (allo.nbBinTypes); // leftmost node that is the head of a reflected arc for each bin
	for (int i = 0; i < allo.nbBinTypes; i++) minVerR[i] = allo.bins[i][0];
	for (int i = 0; i < allo.nbItemTypes; i++){
		vector<bool> done (cap + 1,false);
		for (int j = cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){
				for (int k = 1; k <= allo.items[i][1]; k++){ 
					if(j + (k-1) * allo.items[i][0] * 2 > cap - 1 || done[j + (k-1) * allo.items[i][0] * 2]) break; 
					done[j + (k-1) * allo.items[i][0] * 2] = true;					
					if (j + k * allo.items[i][0] * 2 <= cap){
						arcs.push_back({j + (k-1) * allo.items[i][0] * 2, j + k * allo.items[i][0] * 2,i,-1}); 
						isNodeActive[j + k * allo.items[i][0] * 2] = true;
						isNodeActiveWR[j + k * allo.items[i][0] * 2] = true;
					}
					for(int l = 0; l < allo.nbBinTypes;l++){				
						if (j + (k-1) * allo.items[i][0] * 2 < allo.bins[l][0] &&  j + k * allo.items[i][0] * 2 > allo.bins[l][0]){
							int dest = max(0,allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2));
							arcs.push_back({j + (k-1) * allo.items[i][0] * 2,dest,i,l}); 
							isNodeActiveWR[dest] = true;
							minVerR[l] = min(minVerR[l], dest);
							allo.infos.nbVarS += 1;
						}
					}
				}
			}
		}
	}

	// Loss arcs
	vector<vector<bool> > isNodeActiveWRperBin (allo.nbBinTypes,vector<bool>(cap + 1, false)); 
	for(int l = 0; l < allo.nbBinTypes;l++){
		isNodeActiveWR[allo.bins[l][0]] = true;		
		vector<int> onlyActiveNodesWR;	
		for (int i = max(1,minVerR[l]); i < allo.bins[l][0]; i++){
			if (isNodeActiveWR[i] == true){
				isNodeActiveWRperBin[l][i] = true;
				onlyActiveNodesWR.push_back(i);
			}
		}	
		onlyActiveNodesWR.push_back(allo.bins[l][0]); 
		isNodeActiveWRperBin[l][allo.bins[l][0]] = true;
		for (int i = 0; i < onlyActiveNodesWR.size() - 1; i++){
			arcs.push_back({onlyActiveNodesWR[i], onlyActiveNodesWR[i+1],-1,l});
			allo.infos.nbVarS += 1;
		}
	}

	// Bin arcs	
	for (int i = 0; i < allo.nbBinTypes; i++){
		arcs.push_back({allo.bins[i][0], allo.bins[i][0],-2,i});	
		allo.infos.nbVarS += 1;		
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
		vector<GRBLinExpr> cIn(cap + 1,0);
		vector<GRBLinExpr> cOut(cap + 1,0);
		vector<GRBLinExpr> cRIn(cap + 1,0); 
		vector<vector<GRBLinExpr> > cRInPB(allo.nbBinTypes,vector<GRBLinExpr>(cap + 1, 0));
		vector<GRBLinExpr> nbItemUsed(allo.nbItemTypes,0);
		GRBLinExpr	nbRef = 0;
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			if (arcs[i][2] >= 0){
				isArcUsed[i] = model.addVar(0, allo.items[arcs[i][2]][1], 0, GRB_INTEGER);
			}
			if (arcs[i][2] == -1){
				isArcUsed[i] = model.addVar(0, allo.bins[arcs[i][3]][1], 0, GRB_INTEGER);
			}
			if (arcs[i][2] == -2){
				isArcUsed[i] = model.addVar(-allo.bins[arcs[i][3]][1], allo.bins[arcs[i][3]][1], 0, GRB_INTEGER);
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
			if(arcs[i][2] >= 0) // The arc contains an item 
				nbItemUsed[arcs[i][2]] += isArcUsed[i];
			if(arcs[i][2] == -1){ // The arc is a forward loss arc
				cRInPB[arcs[i][3]][arcs[i][1]] += isArcUsed[i];
				cRIn[arcs[i][1]] += isArcUsed[i];
				cRInPB[arcs[i][3]][arcs[i][0]] -= isArcUsed[i];
				cRIn[arcs[i][0]] -= isArcUsed[i];
			}
			else{ // The arc is not a forward loss arc
				cOut[arcs[i][0]] += isArcUsed[i];
				if(arcs[i][3] == -1) // The arc is a standard arc
					cIn[arcs[i][1]] += isArcUsed[i];
				else{ // The arc is a reflected arc 
					cRInPB[arcs[i][3]][arcs[i][1]] += isArcUsed[i];
					cRIn[arcs[i][1]] += isArcUsed[i];
					nbRef += isArcUsed[i];
					nbBinUsedExpr[arcs[i][3]] += isArcUsed[i];
				}	
			}
		}
		
		// Flow conservation
		for (int i = 1; i <= cap; i++){
			if (isNodeActiveWR[i] == true)
				model.addConstr(cOut[i] + cRIn[i] == cIn[i]);
			for(int l = 0; l < allo.nbBinTypes;l++){
				if (isNodeActiveWRperBin[l][i] == true){
					model.addConstr(cRInPB[l][i] >= 0);
					allo.infos.nbConsS += 1;
				}
			}
		}	
		model.addConstr(cOut[0] + cRIn[0] == 2 * nbRef);
		
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
		vector<vector<vector<int> > > transfer(allo.nbBinTypes,vector<vector<int> >(cap+1));
		for (int i = 0; i < arcs.size(); i++){
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
				cout << arcs[i][0] << " " << arcs[i][1] << " " << arcs[i][2] << " " << arcs[i][3] << endl;
			}
			if(arcs[i][2] >= 0){ // The arc is an item arc 
				for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
					arcsUsed[arcs[i][0]].push_back(arcs[i]);
				}
			}
			if(arcs[i][2] == -1){ // The arc is a forward loss arc 
				for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
					transfer[arcs[i][3]][arcs[i][0]].push_back(arcs[i][1]);
				}
			}
			if(arcs[i][2] == -2){ // The arc is a connection arc
				for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
					arcsUsed[arcs[i][0]].push_back(arcs[i]);
				}	
				for (int j = ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j < 0; j++){
					transfer[arcs[i][3]][arcs[i][0]].push_back(arcs[i][1]);
				}				
			}
		}
		
		// Step 1 ==> Create half bins 
		vector<vector<int> > halfBins;

		while (arcsUsed[0].size() > 0){
			int currTail = 0;
			vector<int> halfBin(2);
			halfBin[0] = 0; // the weight in the bin
			halfBin[1] = -1; // the type of bin
			while (arcsUsed[currTail].size() != 0){
			//	cout << currTail << endl;
				int nextTail = arcsUsed[currTail].back()[1];
				if(arcsUsed[currTail].back()[2] >= 0){
					halfBin.push_back(arcsUsed[currTail].back()[2]);
					halfBin[0] += allo.items[arcsUsed[currTail].back()[2]][0];
				}
				halfBin[1] = arcsUsed[currTail].back()[3];
				arcsUsed[currTail].pop_back();
				currTail = nextTail;
				if(halfBin[1] >= 0){
					int addLoss = 0; 
					while(transfer[halfBin[1]][currTail].size() > 0){
						nextTail = transfer[halfBin[1]][currTail].back();
						addLoss += nextTail - currTail;
						transfer[halfBin[1]][currTail].pop_back();
						if(currTail == nextTail){
							halfBin[1] = -1;
							halfBin[0] -= 0.5*addLoss;
							break;
						}
						currTail = nextTail;
					}
					if(halfBin[1] >= 0){
						halfBin[0] += cap - allo.bins[halfBin[1]][0];
						break;
					}
				}
			}
		//	cout << "STOP" << endl;
			halfBins.push_back(halfBin);
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
			cout << firstHalfIndex << " ";
			for(int j = 0; j < halfBins[firstHalfIndex].size(); j++){
				cout << halfBins[firstHalfIndex][j] << " "; 
			}
			for(int j = 2; j < halfBins[firstHalfIndex].size(); j++) 
				bin.push_back(halfBins[firstHalfIndex][j]); 
			if(halfBins[firstHalfIndex][0] < cap){ // handle the case where an item has weight allo.cap
				cout << " matched with " << secondHalfIndex << " ";
				for(int j = 0; j < halfBins[secondHalfIndex].size(); j++){
					cout <<  halfBins[secondHalfIndex][j] << " "; 
				}
				cout << endl;
				for(int j = 2; j < halfBins[secondHalfIndex].size(); j++) 
					bin.push_back(halfBins[secondHalfIndex][j]); 
				secondHalfIndex--;
			}
			bins[halfBins[firstHalfIndex][1]].push_back(bin);
			firstHalfIndex++;		
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
