#include "Allocation.h" 

/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */
	
void Allocation::load(const string& path, const string& filein){
	// Local variables 
	istringstream iss;
	string parser;
	string nameFile = path + filein;
	
	// File opening
	ifstream file(nameFile.c_str(), ios::in);

	// File lecture
	if (file){
		// Name of the instance is filein
		name = filein;

		// Number of item types
		getline(file, parser); iss.str(parser);
		iss >> nbItemTypes;
		iss.clear();	

		// Capacity
		getline(file, parser); iss.str(parser);
		iss >> cap;  
		iss.clear();

		// For each item type
		for(int i = 0; i < nbItemTypes; i++){
			getline(file, parser); iss.str(parser);
			vector<int> item(2);
			iss >> item[0]; iss >> item[1]; 
			iss.clear();
			items.push_back(item);
		}
		
		file.close();
			
	}
	else cout << "Could not open the file " << nameFile << endl;
}

void Allocation::printProb(){
	cout << "Instance " << name << endl;
	cout << "Bin capacity " << cap << endl;
	for(int i = 0; i < nbItemTypes;i++){
		cout << "Item " << i << " weight " << items[i][0] << " demand " << items[i][1] << endl;
	}
}

void Allocation::printInfo(const string& pathAndFileout){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << name << "\t" << infos.opt << "\t";
	for (int i = 0; i < infos.timeCPU.size();i++)
		file << infos.timeCPU[i] << "\t";
	file << min(infos.ObjVal,infos.ObjBound) << "\t" << max(infos.ObjVal,infos.ObjBound) <<  "\t" << infos.nbVar << "\t" << infos.nbCons << "\t" << infos.nbNZ << "\t" << infos.correct << endl;
	file.close();
}