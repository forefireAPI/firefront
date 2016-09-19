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

#include "FRPModel.h"

namespace libforefire {

/* name of the model */
const string FRPModel::name = "FRP";

/* registration */
int FRPModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getFRPModel );

/* instantiation */
FluxModel* getFRPModel(const int& index, DataBroker* db) {
	return new FRPModel(index, db);
}

/* constructor */
FRPModel::FRPModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	heatFlux = registerProperty("heatFlux");

    /* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
    FRP_max = 240.e3 ;
    FRP_ratio = 0.1;
	burningDuration = 30.;
	if ( params->isValued("burningDuration") )
		burningDuration = params->getDouble("burningDuration");
	/*nominalHeatFlux = 150000.;
	if ( params->isValued("nominalHeatFlux") )
		nominalHeatFlux = params->getDouble("nominalHeatFlux");*/
}

/* destructor (shoudn't be modified) */
FRPModel::~FRPModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string FRPModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */
double FRPModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	/* Mean heat flux released between the time interval [bt, et] */
	/* The heat flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of 'burningDuration', constant of the model */
    
    double heatFlux_here = valueOf[heatFlux];

	/* Instantaneous flux */
	/* ------------------ */
	if ( bt == et ){
		if ( bt < at ) return 0;
		if ( bt < at + burningDuration ) return FRP_ratio * heatFlux_here / FRP_max;
		return 0;
	}

	/* Averaged flux */
	/* ------------- */

	/* looking outside burning interval */
	if ( et < at or bt > at + burningDuration ) return 0;
	/* begin time outside interval, end time inside */
	if ( bt < at and et <= at + burningDuration ) return FRP_ratio * heatFlux_here*(et-at)/(et-bt) / FRP_max;
	/* begin time outside interval, end time outside */
	if ( bt < at and et > at + burningDuration ) return FRP_ratio * heatFlux_here*burningDuration/(et-bt) / FRP_max;
	/* begin time inside interval, end time inside */
	if ( bt >= at and et <= at + burningDuration ) return FRP_ratio * heatFlux_here / FRP_max;
	/* begin time inside interval, end time outside */
	return FRP_ratio * heatFlux_here*(at+burningDuration-bt)/(et-bt) / FRP_max;

}

} /* namespace libforefire */
