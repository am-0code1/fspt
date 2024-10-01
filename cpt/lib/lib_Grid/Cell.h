#pragma once
#include "Node.h"

class Cell
{

public:
	Cell();

	Props props;
	Props_aux props_aux;

	/*
	 * Note:
	 *
	 * Index of boundary nodes on a cell is zero.
	 */

	std::vector<int> zone; // container for zone ids that the cell is part of.
	std::vector<int> node; // container for node ids belonging to the cell.
	std::vector<int> face; // container for face ids that the cell is across one of the side of.

	double area;
	int type, type_init = -1;
	int *bin; // bin that cell belongs to -- bin[0]: bin.x, bin[1]: bin.y
};