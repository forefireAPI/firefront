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

#ifndef Balbi2015_H_
#define Balbi2015_H_

#include "PropagationModel.h"
#include "FireDomain.h"

using namespace std;

namespace libforefire {

class Balbi2015 : public PropagationModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t slope;
	size_t fdepth;
	size_t curvature;
	size_t normalWind;
	size_t moisture;
	size_t temperature;
	size_t Rhod;
	size_t Rhol;
	size_t Md;
	size_t Ml;
	size_t sd;
	size_t sl;
	size_t e;
	size_t Sigmad;
	size_t Sigmal;
	size_t stoch;
	size_t RhoA;
	size_t Ta;
	size_t Tau0;
	size_t Deltah;
	size_t DeltaH;
	size_t Cp;
	size_t Ti;
	size_t X0;
	size_t r00;
	size_t Blai;

	/*! coefficients needed by the model */
	double cooling;
	double Cpa;

	/*! local variables */

	/*! result of the model */
	double getSpeed(double*);

public:

	Balbi2015(const int& = 0, DataBroker* db=0);
	virtual ~Balbi2015();

	string getName();

};

PropagationModel* getBalbi2015Model(const int& = 0, DataBroker* db=0);

} /* namespace libforefire */
#endif /* Balbi2015_H_ */
