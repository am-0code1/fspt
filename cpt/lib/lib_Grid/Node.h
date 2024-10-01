#pragma once
#include "Props.h"
#include <vector>

class Node
{

public:
	Node();

	Props props;

	/*
	 * Note:
	 *
	 * For a boundary node, index of one of its cells is zero.
	 */

	std::vector<int> zone; // container for zone ids that the node is part of.
	std::vector<int> cell; // container for cell ids that the node is part of.
	std::vector<int> face; // container for face ids that the node is part of.

	int type;
	int type_init = -1;
	int *bin; // bin that node resides in -- bin[0]: bin.x, bin[1]: bin.y
};