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

#ifndef HEATFLUXFROMOBSMODEL_H_
#define HEATFLUXFROMOBSMODEL_H_

#include "FluxModel.h"
#include "FireDomain.h"

using namespace std;

namespace libforefire {

class HeatFluxFromObsModel: public FluxModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
    size_t evaporationTime_data;
    size_t residenceTime_data;
    size_t burningTime_data;
    size_t nominalHeatFlux_f_data;
    size_t nominalHeatFlux_s_data;

	/*! coefficients needed by the model */
	//double burningTime;
	//double residenceTime;

	/*double nominalHeatFlux;*/

	/*! local variables */

	/*! result of the model */
	double getValue(double*, const double&
			, const double&, const double&);

public:
	HeatFluxFromObsModel(const int& = 0, DataBroker* = 0);
	virtual ~HeatFluxFromObsModel();

	string getName();
    
};

struct SensibleheatFlux{
    double flaming;
    double smoldering; 

    SensibleheatFlux(double f, double s) : flaming(f), smoldering(s) {}
};



FluxModel* getHeatFluxFromObsModel(const int& = 0, DataBroker* = 0);
    
/*! \compute heat flux from local input */
SensibleheatFlux computeHeatFLuxFromBmap(const double&, const double&, const double&, const double&, const double&, const double&, const double&);

} /* namespace libforefire */
#endif /* HEATFLUXFROMOBSMODEL_H_ */
