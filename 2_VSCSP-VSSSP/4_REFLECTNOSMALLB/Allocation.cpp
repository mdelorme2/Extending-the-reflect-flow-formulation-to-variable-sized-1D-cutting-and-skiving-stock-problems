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

		// Number of bin types
		getline(file, parser); iss.str(parser);
		iss >> nbBinTypes;  
		iss.clear();

		// For each bin
		for(int i = 0; i < nbBinTypes; i++){
			getline(file, parser); iss.str(parser);
			vector<int> bin(3);
			iss >> bin[0]; iss >> bin[1]; iss >> bin[2]; 
			iss.clear();
			bins.push_back(bin);
		}
		
		// Number of item types
		getline(file, parser); iss.str(parser);
		iss >> nbItemTypes;
		iss.clear();	

		// For each item
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
	for(int i = 0; i < nbBinTypes;i++){
		cout << "Bin " << i << " capacity " << bins[i][0] << " availability " << bins[i][1] << " price " << bins[i][2] << endl;
	}
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