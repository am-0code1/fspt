#pragma once
#include <vector>

class Bin
{

public:
	Bin(){};
	~Bin() {delete[] vertex;}

	/*
	 * A spatial bin/box/window to store the id of nodes/cells that reside within the bin.
	 */

	std::vector<int> node; // container for node ids that belong to the bin.

	// cells can be added to bins to be used in the original ASG variant
	std::vector<int> cell; //container for cell ids that belong to the bin.

	double *vertex;
};
