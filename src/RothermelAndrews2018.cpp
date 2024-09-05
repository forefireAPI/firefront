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

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;
	double windReductionFactor;
	/*! properties needed by the model */
	size_t wv_;
	size_t slope_;
	size_t wo_;
	size_t fd_;
	size_t fpsa_;
	size_t mf_;
	size_t pp_;
	size_t h_;
	size_t mois_ext;
	size_t st;
	size_t se;

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

	wv_ = registerProperty("normalWind_ms");
	slope_ = registerProperty("slope_deg");

	wo_ = registerProperty("fl1h_tac");
	fd_ = registerProperty("fd_ft");
	fpsa_ = registerProperty("SAVcar_ftinv");
	mf_ = registerProperty("mdOnDry1h_r");
	pp_ = registerProperty("fuelDens_lbft3");
	h_ = registerProperty("H_BTUlb");
	mois_ext = registerProperty("mois_ext_r");
	st = registerProperty("st_r");
	se = registerProperty("se_r");
	
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

double RothermelAndrews2018::getSpeed(double* valueOf){

	double mois_ext = 0.3;
	double st = 0.0555;
	double se = 0.01;
	double Pi = 3.14159265358;
	double ms2ftmin = 196.85039;
	double ftmin2ms = 0.00508;

	double slope = Pi * valueOf[slope_] / 180; // convert degrees to radians
	double wv = ms2ftmin * valueOf[wv_]; // convert m/s to feat/min
	double pp = valueOf[pp_];
	double mf = valueOf[mf_];
	double fpsa = valueOf[fpsa_];
	double fd = valueOf[fd_];
	double wo = valueOf[wo_];
	double h = valueOf[h_];

	if (wv < 0) wv = 0;

	if(wv < 0){
		double tan_slope = tan(slope);
		double Beta_op = 3.348 * pow(fpsa, -0.8189);  // Optimum packing ratio
		double ODBD = wo / fd; // Ovendry bulk density
        double Beta = ODBD / pp; // Packing ratio
		double WN = wo / (1 + st); // Net fuel loading
        double A =  133.0 / pow(fpsa, 0.7913); // updated A
        double T_max = pow(fpsa,1.5) * pow(495.0 + 0.0594 * pow(fpsa, 1.5),-1.0); // Maximum reaction velocity
		double T = T_max * pow((Beta / Beta_op), A) * exp(A * (1 - Beta / Beta_op));  // Optimum reaction velocity
        double NM = 1. - 2.59 * (mf / mois_ext) + 5.11 * pow(mf / mois_ext, 2.) - 3.52 * pow(mf / mois_ext,3.);  // Moisture damping coeff.
        double NS = 0.174 * pow(se, -0.19);  // Mineral damping coefficient
        double RI = T * WN * h * NM * NS;
        double PFR = pow(192.0 + 0.2595 * fpsa, -1) * exp((0.792 + 0.681 * pow(fpsa, 0.5)) * (Beta + 0.1));  // Propogating flux ratio
        // Wind Coefficien t
        double B = 0.02526 * pow(fpsa, 0.54);
        double C = 7.47 * exp(-0.1333 * pow(fpsa, 0.55));
        double E = 0.715 * exp(-3.59 * pow(10, -4) * fpsa);
        if (wv > (0.9 * RI)) wv = 0.9 * RI;
		double WC = (C * pow(wv, B)) * pow((Beta / Beta_op), (-E));
		double SC = 5.275*(pow(Beta, -0.3))*pow(tan_slope, 2);
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
			return R * ftmin2ms;
		}
	}else{
		return 0;
	}
} /* namespace libforefire */
}