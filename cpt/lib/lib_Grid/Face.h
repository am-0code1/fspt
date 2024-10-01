#pragma once
#include <vector>

class Face
{

public:
	Face();

	std::vector<int> zone; // container for zone ids that the face is part of.
	std::vector<int> node; // container for node ids that belong to the face.
	std::vector<int> cell; // container for cell ids that belong to the face.

	int bc_type, face_type,
		bc_type_init = -1,
		face_type_init = -1;
};
