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

#ifndef FOREFIREV1VAPORFLUXMODEL_H_
#define FOREFIREV1VAPORFLUXMODEL_H_

#include "FluxModel.h"
#include "FireDomain.h"
#include <math.h>

using namespace std;

namespace libforefire {

class ForeFireV1VaporFluxModel: public FluxModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t sd;
	size_t md;
	size_t sigmad;
	size_t ml;
	size_t sigmal;
	size_t lai;
	size_t m_e;

	/*! coefficients needed by the model */
	double BurrowsCst;
	double liveTimeScaleRatio;
	double timeScalesRatio;
	vector<double> LAICoeffs;
	vector<double> mCoeffs;
	double chi_b;

	/*! local variables */

	/*! local functions */
	double texp(double&, double&);
	double tauexp(double&, double&);

	/*! result of the model */
	double getValue(double*, const double&
			, const double&, const double&);

public:
	ForeFireV1VaporFluxModel(const int& = 0, DataBroker* = 0);
	virtual ~ForeFireV1VaporFluxModel();

	string getName();
};

FluxModel* getForeFireV1VaporFluxModel(const int& = 0, DataBroker* = 0);

} /* namespace libforefire */
#endif /* FOREFIREV1VAPORFLUXMODEL_H_ */
