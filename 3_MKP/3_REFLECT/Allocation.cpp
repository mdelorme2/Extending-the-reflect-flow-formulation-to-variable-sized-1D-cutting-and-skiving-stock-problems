#include "Allocation.h" 

/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */
	
void Allocation::load(const string& path, const string& filein){
	// Local variables 
	istringstream iss;
	istringstream tempIss;
	string parser;
	string garbage;
	string nameFile = path + filein;
	string tempString;
	char tempChar;
	
	// File opening
	ifstream file(nameFile.c_str(), ios::in);

	// File lecture
	if (file){
		// Name of the instance is filein
		name = filein;

		// Number of knapsacks
		getline(file, parser); iss.str(parser);
		iss >> nbBins;  
		iss.clear();

		// Number of items
		getline(file, parser); iss.str(parser);
		iss >> nbItems;  
		iss.clear();

		// For each knapsack
		for(int i = 0; i < nbBins; i++){
			getline(file, parser); iss.str(parser);
			int bin;
			iss >> bin; 
			iss.clear();
			bins.push_back(bin);
		}

		// For each item
		for(int i = 0; i < nbItems; i++){
			getline(file, parser); iss.str(parser);
			vector<int> item(2);
			iss >> item[0]; iss >> item[1]; 
			iss.clear();
			items.push_back(item);
		}
		
		// Knapsack and items are not sorted, so we need to sort them
		sort(bins.begin(), bins.end(), [](const int& a, const int& b) {return a > b;}); 
		sort(items.begin(), items.end(), [](const std::vector<int>& a, const std::vector<int>& b) {return a[0]*100000 + a[1] > b[0]*100000 + b[1];}); // Sort
		file.close();
			
	}
	else cout << "Could not open the file " << nameFile << endl;
}

void Allocation::printProb(){
	cout << "Instance " << name << endl;
	for(int i = 0; i < nbBins;i++){
		cout << "Knapsack " << i << " capacity " << bins[i] << endl;
	}
	for(int i = 0; i < nbItems;i++){
		cout << "Item " << i << " weight " << items[i][0] << " price " << items[i][1] << endl;
	}
}

void Allocation::printInfo(const string& pathAndFileout){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << name << "\t" << infos.opt << "\t";
	for (int i = 0; i < infos.timeCPU.size();i++)
		file << infos.timeCPU[i] << "\t";
	file << min(infos.ObjVal,infos.ObjBound) << "\t" << max(infos.ObjVal,infos.ObjBound) <<  "\t" << infos.nbVar << "\t" << infos.nbCons << "\t" << infos.nbNZ << "\t" << infos.nbVarS << "\t" << infos.nbConsS <<"\t" << infos.correct << endl;
	file.close();
}