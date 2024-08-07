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


#include "FluxModel.h"
#include "FireDomain.h"

namespace libforefire {

class VaporFluxFromObsModel: public FluxModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
    size_t evaporationTime_data;
    size_t nominalVaporFlux_data; 
	
    /*! coefficients needed by the model */
	//double FromObsVaporFlux;

	/*! local variables */

	/*! result of the model */
	double getValue(double*, const double&
			, const double&, const double&);

public:
	VaporFluxFromObsModel(const int& = 0, DataBroker* = 0);
	virtual ~VaporFluxFromObsModel();

	string getName();
};

FluxModel* getVaporFluxFromObsModel(const int& = 0, DataBroker* = 0);

double computeVaporFLuxFromBmap(const double&, const double&, 
                                 const double&, const double&, const double&);

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
	evaporationTime_data   = registerProperty("FromObs_evaporationTime");
    nominalVaporFlux_data  = registerProperty("FromObs_VaporFlux_evaporation");

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

double computeVaporFLuxFromBmap(const double& evaporationTime, const double&  nominalVaporFlux, 
                                 const double& bt, const double& et, const double& at){

    /* Instantaneous flux */
	/* ------------------ */
	if ( bt == et ){
		if ( bt < at ) return 0;
		if ( bt < at + evaporationTime ) return nominalVaporFlux;
		return 0;
	}

	/* Averaged flux */
	/* ------------- */

	/* looking outside burning interval */
	if ( et < at or bt > at + evaporationTime ) return 0;
	/* begin time outside interval, end time inside */
	if ( bt < at and et <= at + evaporationTime ) return nominalVaporFlux*(et-at)/(et-bt);
	/* begin time outside interval, end time outside */
	if ( bt < at and et > at + evaporationTime ) return nominalVaporFlux*evaporationTime/(et-bt);
	/* begin time inside interval, end time inside */
	if ( bt >= at and et <= at + evaporationTime ) return nominalVaporFlux;
	/* begin time inside interval, end time outside */
	if ( bt >= at and et > at + evaporationTime ) return nominalVaporFlux*(at+evaporationTime-bt)/(et-bt);

    return -999;
}

double VaporFluxFromObsModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	/* Mean vapor flux released between the time interval [bt, et] */
	/* The vapor flux is supposed to be constant from the arrival time (at)
	 * and for a period of time of given by fuel properties tau0/sd */

    double evaporationTime = valueOf[evaporationTime_data];
    double nominalVaporFlux = valueOf[nominalVaporFlux_data];
   

    double VaporFlux = computeVaporFLuxFromBmap(evaporationTime,nominalVaporFlux, bt, et, at);
	
	return VaporFlux;

}

} /* namespace libforefire */
