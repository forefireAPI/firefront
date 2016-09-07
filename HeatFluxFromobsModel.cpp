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

#include "HeatFluxFromobsModel.h"

namespace libforefire {

/* name of the model */
const string HeatFluxFromobsModel::name = "heatFluxFromobs";

/* registration */
int HeatFluxFromobsModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getHeatFluxFromobsModel );

/* instantiation */
FluxModel* getHeatFluxFromobsModel(const int& index, DataBroker* db) {
	return new HeatFluxFromobsModel(index, db);
}

/* constructor */
HeatFluxFromobsModel::HeatFluxFromobsModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	residenceTime = registerProperty("rtime");

	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
	burningDuration = 30.;
	if ( params->isValued("burningDuration") )
		burningDuration = params->getDouble("burningDuration");
	nominalHeatFlux = 150000.;
	if ( params->isValued("nominalHeatFlux") )
		nominalHeatFlux = params->getDouble("nominalHeatFlux");
}

/* destructor (shoudn't be modified) */
HeatFluxFromobsModel::~HeatFluxFromobsModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string HeatFluxFromobsModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double HeatFluxFromobsModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	/* Mean heat flux released between the time interval [bt, et] */
	/* The heat flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of 'burningDuration', constant of the model */


	return valueOf[residenceTime];

}

} /* namespace libforefire */
