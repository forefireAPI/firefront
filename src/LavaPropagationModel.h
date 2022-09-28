/*

Copyright (C) 2012 ForeFire Team, SPE, Universit� de Corse.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 US

*/

#ifndef LAVAPROPAGATIONMODEL_H_
#define LAVAPROPAGATIONMODEL_H_

#include "PropagationModel.h"
#include "FireDomain.h"

namespace libforefire {

class LavaPropagationModel : public PropagationModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t effectiveSlope;
	size_t viscosity;
	size_t flowSpeed;
	size_t rugosity;
	size_t curvature;
	size_t fastInSection;

	/*! coefficients needed by the model */
	double slowFactor;
	int lavaviscosity;
	double lavadensity;
	int lavayield;
	/*! local variables */

	/*! result of the model */
	double getSpeed(double*);

public:

	LavaPropagationModel(const int& = 0, DataBroker* db=0);
	virtual ~LavaPropagationModel();

	string getName();

};

PropagationModel* getLavaPropagationModel(const int& = 0, DataBroker* db=0);

}

#endif /* LAVAPROPAGATIONMODEL_H_ */
