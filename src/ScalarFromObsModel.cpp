/*

Copyright (C) 2012 ForeFire Team, SPE, UniversitŽ de Corse.

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

#include "ScalarFromObsModel.h"

namespace libforefire {

/* name of the model */
const string ScalarFromObsModel::name = "SFObs";

/* registration */
int ScalarFromObsModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getScalarFromObsModel );

/* instantiation */
FluxModel* getScalarFromObsModel(const int& index, DataBroker* db) {
	return new ScalarFromObsModel(index, db);
}

/* constructor */
ScalarFromObsModel::ScalarFromObsModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	nominalHeatFlux_bmap = registerProperty("FromObs_NominalHeatFlux");
	EFScalar_bmap = registerProperty("FromObs_EFScalar");
	residenceTime_bmap = registerProperty("FromObs_residenceTime");
	radiation_fraction_bmap = registerProperty("FromObs_radiationFraction");
	conversion_factor_bmap = registerProperty("FromObs_conversionFactor");

    /* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
}

/* destructor (shoudn't be modified) */
ScalarFromObsModel::~ScalarFromObsModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string ScalarFromObsModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */
double ScalarFromObsModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	/* Mean heat flux released between the time interval [bt, et] */
	/* The heat flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of 'burningDuration', constant of the model */
    
    double burningDuration = valueOf[residenceTime_bmap];
    double nominalHeatFlux = valueOf[nominalHeatFlux_bmap];
    double EFScalar = valueOf[EFScalar_bmap];

    double HeatFlux = computeHeatFLuxFromBmap(burningDuration,nominalHeatFlux,bt,et,at);
    double radiation_fraction = valueOf[radiation_fraction_bmap];
    double conversion_factor = valueOf[conversion_factor_bmap];

    /*double ScalarFlux = EFScalar * (HeatFlux/(1.-radiation_fraction));      */
    double ScalarFlux = EFScalar * conversion_factor * radiation_fraction * HeatFlux / (1.-radiation_fraction);

	return ScalarFlux;


}

} /* namespace libforefire */
