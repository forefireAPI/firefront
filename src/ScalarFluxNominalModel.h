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

#ifndef SCALARFLUXNOMINALMODEL_H_
#define SCALARFLUXNOMINALMODEL_H_

#include "FluxModel.h"
#include "FireDomain.h"

namespace libforefire {

class ScalarFluxNominalModel: public FluxModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t tau0;
	size_t sd;

	/*! coefficients needed by the model */
	double nominalScalarFlux;

	/*! local variables */

	/*! result of the model */
	double getValue(double*, const double&
			, const double&, const double&);

public:
	ScalarFluxNominalModel(const int& = 0, DataBroker* = 0);
	virtual ~ScalarFluxNominalModel();

	string getName();
};

FluxModel* getScalarFluxNominalModel(const int& = 0, DataBroker* = 0);

} /* namespace libforefire */
#endif /* SCALARFLUXNOMINALMODEL_H_ */
