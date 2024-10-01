#pragma once
#include <string>
#include <vector>

class Zone {

public:
	Zone();

	std::vector<int> node; 
	std::vector<int> cell; 
	std::vector<int> face; 
	
	std::string type, type_init = "None";
	std::string name, name_init = "None";
};

