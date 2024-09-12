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
 


#include "PropagationModel.h"
#include "FireDomain.h"
#include <math.h>
using namespace std;
namespace libforefire {

class RothermelAndrews2018: public PropagationModel {
	
	// ROS model as defined in Table 5 and Table 6 in Andrews, 2018.
    // @book{Andrews_2018, 
    //     title={The Rothermel surface fire spread model and associated developments: A comprehensive explanation}, 
    //     url={http://dx.doi.org/10.2737/RMRS-GTR-371}, 
    //     DOI={10.2737/rmrs-gtr-371}, 
    //     institution={U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station}, 
    //     author={Andrews, Patricia L.}, 
    //     year={2018}
    //     }

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;
	double windReductionFactor;
	/*! properties needed by the model */
	size_t wv;
	size_t slope;
	size_t wo;
	size_t fd;
	size_t fpsa;
	size_t mf;
	size_t pp;
	size_t h;
	size_t st;
	size_t se;
	size_t me;

	/*! result of the model */
	double getSpeed(double*);

public:
	RothermelAndrews2018(const int& = 0, DataBroker* db=0);
	virtual ~RothermelAndrews2018();

	string getName();

};

PropagationModel* getRothermelAndrews2018Model(const int& = 0, DataBroker* db=0);


/* name of the model */
const string RothermelAndrews2018::name = "RothermelAndrews2018";

/* instantiation */
PropagationModel* getRothermelAndrews2018Model(const int & mindex, DataBroker* db) {
	return new RothermelAndrews2018(mindex, db);
}

/* registration */
int RothermelAndrews2018::isInitialized =
		FireDomain::registerPropagationModelInstantiator(name, getRothermelAndrews2018Model );

/* constructor */
RothermelAndrews2018::RothermelAndrews2018(const int & mindex, DataBroker* db)
: PropagationModel(mindex, db) {
	/* defining the properties needed for the model */
	windReductionFactor = params->getDouble("windReductionFactor");

	wv = registerProperty("normalWind");
	slope = registerProperty("slope");

	wo = registerProperty("fuel.fl1h_tac");
	fd = registerProperty("fuel.fd_ft");
	fpsa = registerProperty("fuel.SAVcar_ftinv");
	mf = registerProperty("fuel.mdOnDry1h_r");
	pp = registerProperty("fuel.fuelDens_lbft3");
	h = registerProperty("fuel.H_BTUlb");
	me = registerProperty("fuel.Dme_pc");
	st = registerProperty("fuel.totMineral_r");
	se = registerProperty("fuel.effectMineral_r");
	
	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerPropagationModel(this);

	/* Definition of the coefficients */
}

/* destructor (shoudn't be modified) */
RothermelAndrews2018::~RothermelAndrews2018() {
}

/* accessor to the name of the model */
string RothermelAndrews2018::getName(){
	return name;
}

/* *********************************************** */
/* Model for the propagation velovity of the front */
/* *********************************************** */

double RothermelAndrews2018::getSpeed(double* Z){

	// Constants
	double msToftmin = 196.85039;
	double ftminToms = 0.00508;
	double lbft2Tokgm2 = 4.88243;
	double tacTokgm2 = 0.224;
	double pcTor = 0.01;

	double wvC = msToftmin * Z[wv]; // convert m/s to feat/min
	double meC = Z[me] * pcTor; // convert percentage to ratio
	double woC = Z[wo] * tacTokgm2 / lbft2Tokgm2; // convert US short tons per acre to pounds per square foot

	if (wvC < 0) wvC = 0;

	if(woC > 0){
		double Beta_op = 3.348 * pow(Z[fpsa], -0.8189);  // Optimum packing ratio
		double ODBD = woC / Z[fd]; // Ovendry bulk density
        double Beta = ODBD / pp; // Packing ratio
		double WN = woC / (1 + Z[st]); // Net fuel loading
        double A =  133.0 / pow(fpsa, 0.7913); // updated A
        double T_max = pow(fpsa,1.5) * pow(495.0 + 0.0594 * pow(fpsa, 1.5),-1.0); // Maximum reaction velocity
		double T = T_max * pow((Beta / Beta_op), A) * exp(A * (1 - Beta / Beta_op));  // Optimum reaction velocity
        double NM = 1. - 2.59 * (Z[mf] / meC) + 5.11 * pow(mf / meC, 2.) - 3.52 * pow(mf / meC,3.);  // Moisture damping coeff.
        double NS = 0.174 * pow(Z[se], -0.19);  // Mineral damping coefficient
        double RI = T * WN * Z[h] * NM * NS;
        double PFR = pow(192.0 + 0.2595 * fpsa, -1) * exp((0.792 + 0.681 * pow(fpsa, 0.5)) * (Beta + 0.1));  // Propogating flux ratio
        // Wind Coefficien t
        double B = 0.02526 * pow(fpsa, 0.54);
        double C = 7.47 * exp(-0.1333 * pow(fpsa, 0.55));
        double E = 0.715 * exp(-3.59 * pow(10, -4) * fpsa);
        if (wvC > (0.9 * RI)) wvC = 0.9 * RI;
		double WC = (C * pow(wvC, B)) * pow((Beta / Beta_op), (-E));

		double SC = 0;
		if (Z[slope] >0){ // Rothermel only for upslope, assuming no slope if negative slope is encountered
			5.275*(pow(Beta, -0.3))*pow(Z[slope], 2);
		}
        // Heat sink
        double EHN = exp(-138. / fpsa);  // Effective Heating Number = f(surface are volume ratio)
        double QIG = 250. + 1116. * mf;  // Heat of preignition= f(moisture content)
        // rate of spread (ft per minute) // RI = BTU/ft^2
        double numerator = (RI * PFR * (1 + WC + SC));
        double denominator = (ODBD * EHN * QIG);
        double R = numerator / denominator; // WC and SC will be zero at slope = wind = 0

		if(R <= 0.0) {
			return 0;
		}else{
			return R * ftminToms;
		}
	}else{
		return 0;
	}
} /* namespace libforefire */
}