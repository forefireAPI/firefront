/*
 * LavaCO2FluxModel.cpp
 *
 *  Created on: 30 nov. 2012
 *      Author: jdurand
 */

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

#include "LavaCO2FluxModel.h"

using namespace std;

namespace libforefire {

/* name of the model */
const string LavaCO2FluxModel::name = "LavaCO2Flux";

/* instantiation */
FluxModel* getLavaCO2FluxModel(const int& index, DataBroker* db) {
	return new LavaCO2FluxModel(index, db);
}

/* registration */
int LavaCO2FluxModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getLavaCO2FluxModel );

/* constructor */
LavaCO2FluxModel::LavaCO2FluxModel(
		const int & mindex, DataBroker* db)
	: FluxModel(mindex, db) {
/* defining the properties needed for the model */
/* allocating the vector for the values of these properties */
if ( numProperties > 0 ) properties =  new double[numProperties];
/* registering the model in the data broker */

dataBroker->registerFluxModel(this);

}
/* Definition of the coefficients */

/* destructor (shoudn't be modified) */
LavaCO2FluxModel::~LavaCO2FluxModel() {
if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string LavaCO2FluxModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double LavaCO2FluxModel::getValue(double* valueOf
			, const double& bt, const double& et, const double& at){

	return 1;

}
} /* namespace libforefire */

