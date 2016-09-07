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

#ifndef BALBINOV2011_H_
#define BALBINOV2011_H_

#include "PropagationModel.h"
#include "FireDomain.h"

using namespace std;

namespace libforefire {

class BalbiNov2011 : public PropagationModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t slope;
	size_t normalWind;
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
	double Cpa;
	double cooling;
	double adjustementSlope;
	double adjustementWind;

	/*! local variables */

	/*! result of the model */
	double getSpeed(double*);

public:

	BalbiNov2011(const int& = 0, DataBroker* db=0);
	virtual ~BalbiNov2011();

	string getName();

};

PropagationModel* getBalbiNov2011Model(const int& = 0, DataBroker* db=0);

}

#endif /* BALBINOV2011_H_ */
