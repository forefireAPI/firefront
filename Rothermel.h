/*

Copyright (C) 2012 ForeFire Team, SPE, Université de Corse.

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

#ifndef ROTHERMEL_H_
#define ROTHERMEL_H_

#include "PropagationModel.h"
#include "FireDomain.h"
#include <math.h>

namespace libforefire {

class Rothermel: public PropagationModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t effectiveSlope;
	size_t normalWind20ft;
	size_t sigma;
	size_t rhop;
	size_t delta;
	size_t betaop;
	size_t wl;
	size_t Mf;
	size_t Mx;
	size_t St;
	size_t Se;
	size_t h;
	size_t a;

	/*! coefficients needed by the model */

	/*! local variables */

	/*! result of the model */
	double getSpeed(double*);

public:
	Rothermel(const int& = 0, DataBroker* db=0);
	virtual ~Rothermel();

	string getName();

};

PropagationModel* getRothermelModel(const int& = 0, DataBroker* db=0);

} /* namespace libforefire */
#endif /* ROTHERMEL_H_ */
