/*

Copyright (C) 2012 ForeFire Team, SPE, Université de Corse.

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

#include "Rothermel.h"

using namespace std;

namespace libforefire {

/* name of the model */
const string Rothermel::name = "Rothermel";

/* instantiation */
PropagationModel* getRothermelModel(const int & mindex, DataBroker* db) {
	return new Rothermel(mindex, db);
}

/* registration */
int Rothermel::isInitialized =
		FireDomain::registerPropagationModelInstantiator(name, getRothermelModel );

/* constructor */
Rothermel::Rothermel(const int & mindex, DataBroker* db)
: PropagationModel(mindex, db) {
	/* defining the properties needed for the model */

	slope = registerProperty("slope");
	normalWind = registerProperty("normalWind");
	Rhod = registerProperty("fuel.Rhod");
	Rhol = registerProperty("fuel.Rhol");
	Md = registerProperty("fuel.Md");
	Ml = registerProperty("fuel.Ml");
	sd = registerProperty("fuel.sd");
	sl = registerProperty("fuel.sl");
	e = registerProperty("fuel.e");
	Sigmad = registerProperty("fuel.Sigmad");
	Sigmal = registerProperty("fuel.Sigmal");
	stoch = registerProperty("fuel.stoch");
	RhoA = registerProperty("fuel.RhoA");
	Ta = registerProperty("fuel.Ta");
	Tau0 = registerProperty("fuel.Tau0");
	Deltah = registerProperty("fuel.Deltah");
	DeltaH = registerProperty("fuel.DeltaH");
	Cp = registerProperty("fuel.Cp");
	Ti = registerProperty("fuel.Ti");
	X0 = registerProperty("fuel.X0");
	r00 = registerProperty("fuel.r00");
	Blai = registerProperty("fuel.Blai");




	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerPropagationModel(this);

	/* Definition of the coefficients */
}

/* destructor (shoudn't be modified) */
Rothermel::~Rothermel() {
}

/* accessor to the name of the model */
string Rothermel::getName(){
	return name;
}

/* *********************************************** */
/* Model for the propagation velovity of the front */
/* *********************************************** */

double Rothermel::getSpeed(double* valueOf){


	double lRhod = valueOf[Rhod]* 0.06;// conversion  Pound per cubic foot
	double lRhol = 120;
	double lMd  = valueOf[Md];
	double lsd  = valueOf[sd] / 3.2808399 ;
	double le   = valueOf[e];
	if (le==0) return 0;
	double lSigmad = valueOf[Sigmad] * 2088.54342;// conversion   livre par pieds carré
	double lDeltaH = 8000;// conversion  BTU
	double normal_wind  = valueOf[normalWind]* 196.850394 ; //conversion ft/min
	double localngle =  valueOf[slope];


	normal_wind *= 0.4; // factor in the data seen in 2013

	if (normal_wind < 0) normal_wind = 0;


	double tanangle = tan(localngle);
	if (tanangle<0) tanangle=0;


	double Mchi = 0.3; // Moisture of extinction
	double Etas = 1; // no mineral damping

	double Wn = lSigmad;

	double Etam = 1 - 2.59*(lMd/Mchi) + 5.11*pow((lMd/Mchi),2) - 2.59*pow((lMd/Mchi),3);

	double A = 1/(4.774*pow(lsd, 0.1) - 7.27);

	double Beta = (lRhod)/lRhol;

	double Betaop = 3.338*pow(lsd, -0.8189);

	double RprimeMax = pow(lsd, 1.5)*(1/(495+0.0594*pow(lsd, 1.5)));

	double Rprime = pow(RprimeMax*(Beta/Betaop),A) *  exp(A*(1-(Beta/Betaop))) ;

	double chi = pow(192 + 0.2595*lsd, -1) * exp((0.792 + 0.681*pow(lsd, 0.5))*(Beta+0.1));

	double epsilon  = exp(-138/lsd);

	double Qig = 250 + 1.116*lMd;

	double C = 7.47 * exp(-0.133*pow(lsd, 0.55));

	double B = 0.02526*pow(lsd, 0.54);

	double E = 0.715* exp(-3.59*(10E-4 * lsd));

	double Ir = Rprime*Wn*lDeltaH*Etam*Etas;

	/// wind limit 2013 //
	double Uf = 96.81*pow(Ir, 1/3);

	if (normal_wind>Uf) {
		normal_wind = Uf;
	}


	double phiV = C* pow((Beta/Betaop),E)*pow(normal_wind,B) ;

	double phiP = 5.275*pow(Beta, -0.3)*pow(tanangle,2);

	double R0 = (Ir*chi)/(lRhod*epsilon*Qig);


	double R = R0*(1+phiV+phiP);


	if(R < R0)  R = R0;

	if(R > 0.0) {

		return R*0.00508; //piedsmin en Ms
	}
	return 0;








}

} /* namespace libforefire */
