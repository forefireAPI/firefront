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

#include "VaporFluxFromObsModel.h"

namespace libforefire {

/* name of the model */
const string VaporFluxFromObsModel::name = "vaporFluxFromObs";

/* instantiation */
FluxModel* getVaporFluxFromObsModel(const int& index, DataBroker* db) {
	return new VaporFluxFromObsModel(index, db);
}

/* registration */
int VaporFluxFromObsModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getVaporFluxFromObsModel );

/* constructor */
VaporFluxFromObsModel::VaporFluxFromObsModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	nominalHeatFlux_bmap = registerProperty("FromObs_NominalHeatFlux");
	Moisture_bmap = registerProperty("FromObs_Moisture");
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
VaporFluxFromObsModel::~VaporFluxFromObsModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string VaporFluxFromObsModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double VaporFluxFromObsModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	/* Mean vapor flux released between the time interval [bt, et] */
	/* The vapor flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of given by fuel properties tau0/sd */

    double burningDuration = valueOf[residenceTime_bmap];
    double nominalHeatFlux = valueOf[nominalHeatFlux_bmap];
    double Moisture = valueOf[Moisture_bmap];
    
    double HeatFlux = computeHeatFLuxFromBmap(burningDuration,nominalHeatFlux,bt,et,at);
    double radiation_fraction = valueOf[radiation_fraction_bmap];
    double conversion_factor = valueOf[conversion_factor_bmap];
	
	/*double VaporFlux = EFVapor * (HeatFlux/(1.-radiation_fraction));*/
    double VaporFlux = Moisture * conversion_factor * radiation_fraction * HeatFlux/(1.-radiation_fraction);
	
	return VaporFlux;

}

} /* namespace libforefire */
