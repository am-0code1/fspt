#include "Cell.h"

Cell::Cell()
{
	props.pos = new double[props.NDIM];
	props.vel = new double[props.NDIM];
	bin = new int[props.NDIM];

	for (int i = 0; i < props.NDIM; i++)
	{
		props.pos[i] = 0;
		props.vel[i] = 0;
		bin[i] = -1;
	}

	props.pressure = 0;
	props.density = 0;

	//-----------------------------
	//aux props
	//-----------------------------
	props_aux.vel_mag=0;
	
	//-----------------------------
	//Ansys fluent mesh properties
	//-----------------------------
	type = type_init;
	area = 0;

}
