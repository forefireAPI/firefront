/*

Copyright (C) 2012 ForeFire Team, SPE, Universite de Corse.

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

#include "HeatFluxFromObsModel.h"

namespace libforefire {

/* name of the model */
const string HeatFluxFromObsModel::name = "heatFluxFromObs";

/* registration */
int HeatFluxFromObsModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getHeatFluxFromObsModel );

/* instantiation */
FluxModel* getHeatFluxFromObsModel(const int& index, DataBroker* db) {
	return new HeatFluxFromObsModel(index, db);
}

/* constructor */
HeatFluxFromObsModel::HeatFluxFromObsModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	residenceTime_bmap = registerProperty("FromObs_residenceTime");
	nominalHeatFlux_bmap = registerProperty("FromObs_NominalHeatFlux");

    /* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
	/*nominalHeatFlux = 150000.;
	if ( params->isValued("nominalHeatFlux") )
		nominalHeatFlux = params->getDouble("nominalHeatFlux");*/

}

/* destructor (shoudn't be modified) */
HeatFluxFromObsModel::~HeatFluxFromObsModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string HeatFluxFromObsModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double computeHeatFLuxFromBmap(const double& burningDuration, const double&  nominalHeatFlux, 
		 const double& bt, const double& et, const double& at){
	
    /* Mean heat flux released between the time interval [bt, et] */
	/* The heat flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of 'burningDuration', constant of the model */

	/* Instantaneous flux */
	/* ------------------ */

	if ( bt == et ){
		if ( bt < at ) return 0;
		if ( bt < at + burningDuration ) return nominalHeatFlux;
		return 0;
	}

	/* Averaged flux */
	/* ------------- */

	/* looking outside burning interval */
	if ( et < at or bt > at + burningDuration ) return 0;
	/* begin time outside interval, end time inside */
	if ( bt < at and et <= at + burningDuration ) return nominalHeatFlux*(et-at)/(et-bt);
	/* begin time outside interval, end time outside */
	if ( bt < at and et > at + burningDuration ) return nominalHeatFlux*burningDuration/(et-bt);
	/* begin time inside interval, end time inside */
	if ( bt >= at and et <= at + burningDuration ) return nominalHeatFlux;
	/* begin time inside interval, end time outside */
	return nominalHeatFlux*(at+burningDuration-bt)/(et-bt);
}

double HeatFluxFromObsModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	
    /* Mean heat flux released between the time interval [bt, et] */
	/* The heat flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of 'burningDuration', constant of the model */


    double burningDuration = valueOf[residenceTime_bmap];
    double nominalHeatFlux = valueOf[nominalHeatFlux_bmap];
	
	double mm; 
	mm =  computeHeatFLuxFromBmap(burningDuration,nominalHeatFlux,bt,et,at);
	
	/*if (mm>0)
		cout << et << ' ' << bt << ' ' << at << ' ' << burningDuration << ' ' << nominalHeatFlux << ' ' << mm << '\n';
	*/

	return mm;
	/*return computeHeatFLuxFromBmap(burningDuration,nominalHeatFlux,bt,et,at);*/

}

} /* namespace libforefire */
