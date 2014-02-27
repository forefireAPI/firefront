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

#include "BalbiNov2011.h"

namespace libforefire {

/* name of the model */
const string BalbiNov2011::name = "BalbiNov2011";

/* instantiation */
PropagationModel* getBalbiNov2011Model(const int & mindex, DataBroker* db) {
	return new BalbiNov2011(mindex, db);
}

/* registration */
int BalbiNov2011::isInitialized =
        FireDomain::registerPropagationModelInstantiator(name, getBalbiNov2011Model );

/* constructor */
BalbiNov2011::BalbiNov2011(const int & mindex, DataBroker* db)
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
	cooling = 0.;
	if ( params->isValued("BalbiNov2011.cooling") )
		cooling = params->getDouble("BalbiNov2011.cooling");
	Cpa = 1004.;
	if ( params->isValued("BalbiNov2011.Cpa") )
		Cpa = params->getDouble("BalbiNov2011.Cpa");
}

/* destructor (shoudn't be modified) */
BalbiNov2011::~BalbiNov2011() {
}

/* accessor to the name of the model */
string BalbiNov2011::getName(){
	return name;
}

/* *********************************************** */
/* Model for the propagation velovity of the front */
/* *********************************************** */

double BalbiNov2011::getSpeed(double* valueOf){

	double lRhod = valueOf[Rhod];
	double lRhol = valueOf[Rhol];
	double lMd  = valueOf[Md];
	double lMl  = valueOf[Ml];
	double lsd  = valueOf[sd];
	double lsl  = valueOf[sl];
	double le   = valueOf[e];
	double lSigmad = valueOf[Sigmad];
	double lSigmal = valueOf[Sigmal];
	double lstoch = valueOf[stoch];
	double lRhoA  = valueOf[RhoA];
	double lTa  = valueOf[Ta];
	double lTau0  = valueOf[Tau0];
	double lDeltah   = valueOf[Deltah];
	double lDeltaH = valueOf[DeltaH];
	double lCp = valueOf[Cp];
	double lTi  = valueOf[Ti];
	double lX0  = valueOf[X0];
	double lr00  = valueOf[r00];
	double lai = valueOf[Blai];

	double cosCurv = 1;

	if(le <= 0 ) return 0;

	double Betad =   lSigmad /(le* lRhod);
	double Betal =   lSigmal /(le* lRhol);
	double Sd = lsd * le * Betad;
	double Sl = lsl * le * Betal;
	double nu = min((Sd) / lai, 1.);
	double normal_wind = valueOf[normalWind];
	double B = 5.6E-8;
	double a   = lDeltah/ ((lCp*(lTi-lTa)));
	double r0 = lsd * lr00;
	double A0 = (lX0*lDeltaH)/(4*lCp*(lTi-lTa));
	double xsi = ((lMl-lMd)*((lSigmal/lSigmad)*(lDeltah/lDeltaH)));
	double A  = cosCurv * (nu*A0 / (1 + a * lMd)) * (1-xsi);
	double T = lTa + ( lDeltaH*(1-lX0)*(1-xsi) )     / ((lstoch+1)*Cpa);
	double R00 = (B*T*T*T*T)   / (lCp*(lTi-lTa));
	double R = 0;
	double R0 =0;
	double u00 = (2*lai*(lstoch+1)*T*lRhod)/(lRhoA*lTa*lTau0);
	double u0 = nu * u00;

	double  tanGamma =  valueOf[slope] + (normal_wind/u0);
	double gamma = atan(tanGamma);

	R0 = (le / lSigmad)   * (R00) / (1 + a * lMd) * Sd/(Sd+Sl) * Sd/(Sd+Sl);
	if ( gamma > 0. ) {
		double geomFactor = r0*((1+sin(gamma)-cos(gamma))/(1.+cos(gamma)));
		double Rt = R0 + A*geomFactor;
		R = 0.5*( Rt + sqrt( Rt*Rt + 4.*r0*R0/cos(gamma) ) );
	} else {
		R = R0;
	}

	return R ;

}

}
