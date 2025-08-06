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
	vector<vector<bool> > isNodeActiveTran (allo.nbBinTypes,vector<bool>(cap + 1,false)); 
	vector<vector<int> > isNodeActiveDown (allo.nbBinTypes,vector<int>(cap + 1,cap + 1)); 
	vector<bool> isNodeActiveDownT (cap + 1,false);
	isNodeActive[0] = true;
	isNodeActiveWR[0] = true;
	vector<vector<int> > arcs;

	// Item arcs (goes from 0 to cap, all item lentghs are multiplied by 2)
	int minVerR = cap; // leftmost node that is the head of a reflected arc, or the tail of a transposed arc	
	for (int i = 0; i < allo.nbItemTypes; i++){
		vector<bool> done (cap + 1,false);
		bool doneR = false;
		for (int j = cap - 1; j >= 0; j--){
			if (isNodeActive[j] == true){
				for (int k = 1; k <= allo.items[i][1]; k++){ 
					if(j + (k-1) * allo.items[i][0] * 2 > cap - 1 || done[j + (k-1) * allo.items[i][0] * 2]) break; 
					done[j + (k-1) * allo.items[i][0] * 2] = true;	
					if (j + k * allo.items[i][0] * 2 <= cap){
					//	if(j + k * allo.items[i][0] * 2 != cap || !doneR) to uncomment to avoid rare duplicates, but left there to ensure fairness with other methods
							arcs.push_back({j + (k-1) * allo.items[i][0] * 2, j + k * allo.items[i][0] * 2,i,-1}); 
						isNodeActive[j + k * allo.items[i][0] * 2] = true;
						isNodeActiveWR[j + k * allo.items[i][0] * 2] = true;
						// if the reflected arc would have been created for a smaller bin
						for(int l = 1; l < allo.nbBinTypes;l++){
							if (2*allo.bins[l][0] > cap && j + (k-1) * allo.items[i][0] * 2 <= allo.bins[l][0] &&  j + k * allo.items[i][0] * 2 > allo.bins[l][0]){
								if(i == allo.nbItemTypes - 1 || 2*(allo.bins[l][0] - (j + (k-1) * allo.items[i][0] * 2)) > allo.items[i+1][0] * 2){
									// Compute the head of the reflected transfer arc, which is the maximum between d and e
									int head = max(j + (k-1) * allo.items[i][0] * 2,allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2));
									isNodeActiveTran[l][allo.bins[l][0] * 2 - head] = true;
								}		
							}
						}
					}
					else{ // only reflected arc for the largest bin
						// If arc (d,e,i,0) such that d <= e 
						if(j + (k-1) * allo.items[i][0] * 2 <= cap * 2 - (j + k * allo.items[i][0] * 2)){
							arcs.push_back({j + (k-1) * allo.items[i][0] * 2, cap * 2 - (j + k * allo.items[i][0] * 2),i,0}); 
							isNodeActiveWR[cap * 2 - (j + k * allo.items[i][0] * 2)] = true;
							minVerR = min(minVerR, cap * 2 - (j + k * allo.items[i][0] * 2));
							allo.infos.nbVarS += 1;
						}
						else{ // If arc (d,e,i,0) such that e < d, , add (d,d,i,0) instead, and only add the arc if i is the last item or if item i + 1 would not produce the same arc (d,d,i,0) 
							if(i == allo.nbItemTypes - 1 || 2*(cap - (j + (k-1) * allo.items[i][0] * 2)) > allo.items[i+1][0] * 2){
								arcs.push_back({j + (k-1) * allo.items[i][0] * 2, j + (k-1) * allo.items[i][0] * 2,i,0}); 
								isNodeActiveWR[j + (k-1) * allo.items[i][0] * 2] = true;
								minVerR = min(minVerR, j + (k-1) * allo.items[i][0] * 2);	
								allo.infos.nbVarS += 1;					
							}									
						}						
						if(!doneR && allo.items[i][0] < cap){
							int org = max(0,cap-2*allo.items[i][0]);
							// If org = 0, only add the arc if i is the last item or if item i - 1 would not produce the same arc (0,c,i,0) 
							if((org != 0) /*&& !done[org]) to uncomment to avoid rare duplicates, but left there to ensure fairness with other methods */ || i == allo.nbItemTypes - 1 || allo.items[i+1][0] * 2 < cap){
								arcs.push_back({org,cap,i,-1});
								isNodeActiveWR[org] = true;
								minVerR = min(minVerR, org);
								doneR = true;	
								allo.infos.nbVarS += 1;	
							}							
						}			
						// if the reflected arc would also have been created for a smaller bin						
						for(int l = 1; l < allo.nbBinTypes;l++){				
							if (2*allo.bins[l][0] > cap && j + (k-1) * allo.items[i][0] * 2 < allo.bins[l][0] &&  j + k * allo.items[i][0] * 2 > allo.bins[l][0]){
								if(i == allo.nbItemTypes - 1 || 2*(allo.bins[l][0] - (j + (k-1) * allo.items[i][0] * 2)) > allo.items[i+1][0] * 2){
									cout << "arc " << j + (k-1) * allo.items[i][0] * 2 << " " << allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2) << " " << i << " " << l << " should have been created " << endl;
									// Compute the head of the downgrade transfer arc, which is the maximum between d and e
									int head = max(j + (k-1) * allo.items[i][0] * 2,allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2));
									if(allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2) < 0) head = allo.bins[l][0] * 2 - (j + k * allo.items[i][0] * 2);
									isNodeActiveDown[l][cap * 2 - (j + k * allo.items[i][0] * 2)] = min(isNodeActiveDown[l][cap * 2 - (j + k * allo.items[i][0] * 2)],head);
									isNodeActiveDownT[cap * 2 - (j + k * allo.items[i][0] * 2)] = true;
								}
							}
						}
					}
				}
			}
		}
	}

	// Bin arcs	
	for (int i = 0; i < allo.nbBinTypes; i++){
		if(2*allo.bins[i][0] > cap){
			isNodeActiveTran[i][allo.bins[i][0]] = allo.bins[i][0];
			for(int j = 0; j <= cap; j++){
				if(isNodeActiveTran[i][j]){
					cout << "Reflected transfer arc added " << j << " " << 2*allo.bins[i][0]-j << " -1 " << i << endl;
					arcs.push_back({j,2*allo.bins[i][0]-j,-1,i});
					isNodeActiveWR[j] = true;
					isNodeActiveWR[2*allo.bins[i][0]-j] = true;
					minVerR = min(minVerR, 2*allo.bins[i][0]-j);
					allo.infos.nbVarS += 1;				
				}
				if(isNodeActiveDown[i][j] < cap + 1){
					int dest = max(0,isNodeActiveDown[i][j]);
					cout << "Downgrade arc added " << j << " " << dest << " -3 " << i << endl;
					arcs.push_back({j,dest,-3,i});
					isNodeActiveWR[dest] = true;
					minVerR = min(minVerR,dest);
					allo.infos.nbVarS += 1;	
				}
			}
		}
		else{
			arcs.push_back({0, 2*allo.bins[i][0],-2,i});
			isNodeActiveWR[2*allo.bins[i][0]] = true;
			minVerR = min(minVerR, 2*allo.bins[i][0]);
			allo.infos.nbVarS += 1;
		}			
	}

	// Loss arcs
	cout << "minVerR is " <<  minVerR << endl;
	vector<int> onlyActiveNodesWR;	
	for (int i = max(1,minVerR); i <= cap; i++){
		if (isNodeActiveWR[i] == true){
			onlyActiveNodesWR.push_back(i);
		}
	}	
	for (int i = 0; i < onlyActiveNodesWR.size() - 1; i++){
		arcs.push_back({onlyActiveNodesWR[i+1], onlyActiveNodesWR[i],-1,-1});
		allo.infos.nbVarS += 1;
	}
	
	// Print arcs
	cout << arcs.size() << " arcs " << endl;
	for(int i=0;i<arcs.size();i++){
		for(int j=0;j<arcs[i].size();j++){
			cout << arcs[i][j] << " ";
		}
		cout << endl;
	}
	
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
		vector<GRBVar>  isItemDowngraded(allo.nbItemTypes - 1); // The last item cannot be downgraded
		vector<GRBLinExpr> cIn(cap + 1,0);
		vector<GRBLinExpr> cOut(cap + 1,0);
		vector<GRBLinExpr> cRIn(cap + 1,0); 
		vector<GRBLinExpr> cDAl(cap + 1,0);
		vector<GRBLinExpr> nbItemUsed(allo.nbItemTypes,0);
		GRBLinExpr	nbRef = 0;
		
		// Initialization
		for (int i = 0; i < arcs.size(); i++){
			isArcUsed[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
		}
	
		for (int i = 0; i < allo.nbItemTypes - 1; i++){
			isItemDowngraded[i] = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER);
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
			if(arcs[i][2] == -3){ // The arc is a downgrade arc
				cDAl[arcs[i][0]] -= isArcUsed[i];
				cRIn[arcs[i][0]] -= isArcUsed[i];	
				cRIn[arcs[i][1]] += isArcUsed[i];	
				nbBinUsedExpr[arcs[i][3]] += isArcUsed[i];
				nbBinUsedExpr[0] -= isArcUsed[i];
			}
			else{ // The arc is not a downgrade arc
				cOut[arcs[i][0]] += isArcUsed[i];
				if(arcs[i][3] == -1) // The arc is a standard arc
					cIn[arcs[i][1]] += isArcUsed[i];
				else{ // The arc is a reflected arc 
					if(arcs[i][3] == 0) // The arc is a reflected arc for bin 0 
						cDAl[arcs[i][1]] += isArcUsed[i];
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
			if(isNodeActiveDownT[i])
				model.addConstr(cDAl[i] >= 0);
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
		vector<int> downgraded (allo.nbItemTypes,0); // for the items
		vector<vector<vector<int> > > arcsUsed(cap+1);
		vector<vector<int > > downgrade(cap+1); // for the bins
		for (int i = 0; i < arcs.size(); i++){
			if(arcs[i][2] != -3){
				for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
					cout << arcs[i][0] << " " << arcs[i][1] << " " << arcs[i][2] << " " << arcs[i][3] << endl;
					arcsUsed[arcs[i][0]].push_back(arcs[i]);
				}
			}
			else{
				for (int j = 0; j < ceil(isArcUsed[i].get(GRB_DoubleAttr_X) - EPSILON); j++){
					cout << arcs[i][0] << " " << arcs[i][1] << " " << arcs[i][2] << " " << arcs[i][3] << endl;
					downgrade[arcs[i][0]].push_back(arcs[i][3]);
				}
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
					halfBin[0] += cap - allo.bins[halfBin[1]][0];
					if(halfBin[1] == 0 && downgrade[currTail].size() >= 1){
						halfBin[1] = downgrade[currTail].back();
						halfBin[0] += cap - allo.bins[halfBin[1]][0];
						downgrade[currTail].pop_back();
						currTail -= 2*(cap - allo.bins[halfBin[1]][0]);
						currTail = max(0,currTail);
					}
					break;
				}
			}
		//	cout << "STOP" << endl;
			if(currTail != 0) halfBins.push_back(halfBin);
			else{
				vector<int> bin;
				cout << "Selfsufficient bin ";
				for(int j = 2; j < halfBin.size(); j++) cout << halfBin[j] << " "; 
				for(int j = 2; j < halfBin.size(); j++){
					bin.push_back(halfBin[j]); 
					while(bin.back() != 0 && downgraded[bin.back() - 1] > 0){
						downgraded[bin.back() - 1] -= 1;
						bin.back() -= 1;
					}
				}
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
			cout << firstHalfIndex << " ";
			for(int j = 0; j < halfBins[firstHalfIndex].size(); j++){
				cout << halfBins[firstHalfIndex][j] << " "; 
			}
			for(int j = 2; j < halfBins[firstHalfIndex].size(); j++){ 
				bin.push_back(halfBins[firstHalfIndex][j]); 
				while(bin.back() != 0 && downgraded[bin.back() - 1] > 0){
					downgraded[bin.back() - 1] -= 1;
					bin.back() -= 1;
				}
			}			
			cout << " matched with " << secondHalfIndex << " ";
			for(int j = 0; j < halfBins[secondHalfIndex].size(); j++){
				cout <<  halfBins[secondHalfIndex][j] << " "; 
			}
			cout << endl;
			for(int j = 2; j < halfBins[secondHalfIndex].size(); j++){ 
				bin.push_back(halfBins[secondHalfIndex][j]); 
				while(bin.back() != 0 && downgraded[bin.back() - 1] > 0){
					downgraded[bin.back() - 1] -= 1;
					bin.back() -= 1;
				}
			}
			secondHalfIndex--;
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
