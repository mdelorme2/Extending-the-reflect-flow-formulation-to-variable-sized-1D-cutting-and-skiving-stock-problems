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
	vector<bool> isNodeActive (allo.cap + 1,false); 
	vector<bool> isNodeActiveWR (allo.cap + 1,false); 
	isNodeActive[0] = true;
	isNodeActiveWR[0] = true;
	vector<vector<int> > arcs;

	// Item arcs (goes from 0 to cap, all item lentghs are multiplied by 2)		
	int minVerR = allo.cap; // leftmost node that is the head of a reflected arc, or the tail of a transposed arc
	for (int i = 0; i < allo.nbItemTypes; i++){
		vector<bool> done (allo.cap,false);
		bool doneR = false;
		for (int j = allo.cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){
				for (int k = 1; k <= allo.items[i][1]; k++){ 
					if(j + (k-1) * allo.items[i][0] * 2 >= allo.cap || done[j + (k-1) * allo.items[i][0] * 2]) break;
					done[j + (k-1) * allo.items[i][0] * 2] = true;					
					if(j + k * allo.items[i][0] * 2 <= allo.cap){
						arcs.push_back({j + (k-1) * allo.items[i][0] * 2, j + k * allo.items[i][0] * 2,i,-1}); 
						isNodeActive[j + k * allo.items[i][0] * 2] = true;
						isNodeActiveWR[j + k * allo.items[i][0] * 2] = true;
					}
					else{
						arcs.push_back({j + (k-1) * allo.items[i][0] * 2, allo.cap * 2 - (j + k * allo.items[i][0] * 2),i,0}); 
						isNodeActiveWR[allo.cap * 2 - (j + k * allo.items[i][0] * 2)] = true;
						minVerR = min(minVerR, allo.cap * 2 - (j + k * allo.items[i][0] * 2));
						if(!doneR){
							int org = max(0,allo.cap-2*allo.items[i][0]);
							arcs.push_back({org,allo.cap,i,-1});
							isNodeActiveWR[org] = true;
							minVerR = min(minVerR, org);
							doneR = true;		
						}
					}
				}
			}
		}
	}

	cout << arcs.size() << " item arcs " << endl;
	
	// Loss arcs
	isNodeActiveWR[allo.cap] = true;
	vector<int> onlyActiveNodes;
	for (int i = max(1,minVerR); i <= allo.cap; i++){
		if (isNodeActiveWR[i] == true){
			onlyActiveNodes.push_back(i);
		}
	}

	for (int i = 0; i < onlyActiveNodes.size() - 1; i++){
		arcs.push_back({onlyActiveNodes[i+1], onlyActiveNodes[i],-1,-1});
	}
	
	// Connection arcs
	arcs.push_back({allo.cap, allo.cap, -1, 0});
	
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
		GRBVar nbBinUsed;
		if(version >= 1) nbBinUsed = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER); 
		vector<GRBVar>  isArcUsed(arcs.size());
		vector<GRBVar>  isItemDowngraded(allo.nbItemTypes - 1); // The last item cannot be downgraded
		vector<GRBLinExpr> cIn(allo.cap + 1,0);
		vector<GRBLinExpr> cOut(allo.cap + 1,0);
		vector<GRBLinExpr> cRIn(allo.cap + 1,0);
		vector<GRBLinExpr> nbItemUsed(allo.nbItemTypes,0);
		GRBLinExpr	nbRef = 0;
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			if (arcs[i][2] >= 0) 
				isArcUsed[i] = model.addVar(0, allo.items[arcs[i][2]][1], 0, GRB_INTEGER);
			else 
				isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
		}
		
		for (int i = 0; i < allo.nbItemTypes - 1; i++){
			isItemDowngraded[i] = model.addVar(0, allo.items[i][1], 0, GRB_INTEGER);
		}	
		
		model.update();
	
		// Perform values
		for (int i = 0; i < arcs.size(); i++){		
			if(arcs[i][2] >= 0) // The arc contains an item 
				nbItemUsed[arcs[i][2]] += isArcUsed[i];
			cOut[arcs[i][0]] += isArcUsed[i];
			if(arcs[i][3] == -1) // The arc is a standard arc
				cIn[arcs[i][1]] += isArcUsed[i];
			else{ // The arc is a reflected arc 
				cRIn[arcs[i][1]] += isArcUsed[i];
				nbRef += isArcUsed[i];
			}	
		}
		
		// Flow conservation
		for (int i = 1; i <= allo.cap; i++){
			if (isNodeActiveWR[i] == true)
				model.addConstr(cOut[i] + cRIn[i] ==  cIn[i]);
		}	
		model.addConstr(cOut[0] + cRIn[0] == 2 * nbRef);

		// Item availability constraints
		for (int k = 0; k < allo.nbItemTypes; k++){
			if(k == 0)
				model.addConstr(nbItemUsed[k] + isItemDowngraded[k] <= allo.items[k][1]); 
			else{
				if(k == allo.nbItemTypes - 1)
					model.addConstr(nbItemUsed[k] <= allo.items[k][1] + isItemDowngraded[k-1]);
				else 
					model.addConstr(nbItemUsed[k] + isItemDowngraded[k] <= allo.items[k][1] + isItemDowngraded[k-1]);
			}
		}
					
		// Objective function
		if(version >= 1){
			objFun += nbBinUsed;
			model.addConstr(nbRef >= nbBinUsed); 
		}			
		else{
			objFun += nbRef;
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
		if(model.get(GRB_IntAttr_Status) == 2) allo.infos.ObjBound = allo.infos.ObjVal;
		if(allo.infos.ObjVal == allo.infos.ObjBound) allo.infos.opt = true;

		// Filling Solution 	
		vector<vector<int> > bins;
		vector<int> downgraded (allo.nbItemTypes,0);
		vector<vector<vector<int> > > arcsUsed(allo.cap+1);
		for (int i = 0; i < arcs.size(); i++){ 
			for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
			}
		}
		for (int k = 0; k < allo.nbItemTypes - 1; k++){
			downgraded[k] = ceil(isItemDowngraded[k].get(GRB_DoubleAttr_X) - EPSILON);
		}
		
		// Step 1 ==> Create half bins 
		vector<vector<int> > halfBins;

		while (arcsUsed[0].size() > 0){
			int currTail = 0;
			vector<int> halfBin(2);
			halfBin[0] = 0; // the weight in the bin
			halfBin[1] = -1; // the type of bin
			while (arcsUsed[currTail].size() != 0){			
				int nextTail = arcsUsed[currTail].back()[1];
				if (arcsUsed[currTail].back()[2] >= 0){
					halfBin.push_back(arcsUsed[currTail].back()[2]);
					halfBin[0] += allo.items[arcsUsed[currTail].back()[2]][0];
				}
				halfBin[1] = arcsUsed[currTail].back()[3];
				arcsUsed[currTail].pop_back();
				currTail = nextTail;
				if(halfBin[1] == 0) break;
			}
			if(currTail != 0) halfBins.push_back(halfBin);
			else{
				vector<int> bin;
				for(int j = 2; j < halfBin.size(); j++) bin.push_back(halfBin[j]); 
				bins.push_back(bin);
			}
		}

		// Step 2 ==> Order the half bins and match the largest with the smallest
		sort(halfBins.begin(), halfBins.end(), [](const std::vector<int>& a, const std::vector<int>& b) {return a[0] > b[0];}); // Sort
		int firstHalfIndex = 0; int secondHalfIndex = halfBins.size()-1;
		while (firstHalfIndex <= secondHalfIndex){
			vector<int> bin;
			for(int j = 2; j < halfBins[firstHalfIndex].size(); j++){
				bin.push_back(halfBins[firstHalfIndex][j]); 
				while(bin.back() != 0 && downgraded[bin.back() - 1] > 0){
					downgraded[bin.back() - 1] -= 1;
					bin.back() -= 1;
				}
			}
			for(int j = 2; j < halfBins[secondHalfIndex].size(); j++){
				bin.push_back(halfBins[secondHalfIndex][j]); 
				while(bin.back() != 0 && downgraded[bin.back() - 1] > 0){
					downgraded[bin.back() - 1] -= 1;
					bin.back() -= 1;
				}
			}
			secondHalfIndex--;
			firstHalfIndex++;
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
