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
	sigma = registerProperty("sigma");
	rhop = registerProperty("rhop");
	delta = registerProperty("delta");
	betaop = registerProperty("betaop");
	wl = registerProperty("wl");
	Mf = registerProperty("Mf");
	Mx = registerProperty("Mx");
	St = registerProperty("St");
	Se = registerProperty("Se");
	h = registerProperty("h");
	a = registerProperty("a");
	effectiveSlope = registerProperty("effectiveSlope");
	normalWind20ft = registerProperty("normalWind20ft");

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
	/* Computing the spread rate without wind and slope R0 */
	double GammaMax = pow(valueOf[sigma],1.5)/(495.+0.594*pow(valueOf[sigma],1.5));
	double w0 = valueOf[wl]/(1.+valueOf[Mf]);
	double beta = w0/(valueOf[rhop]*valueOf[delta]);
	double A = 1./(4.77*pow(valueOf[sigma],0.1)-7.27);
	double Gamma = GammaMax*pow(beta/valueOf[betaop], A)*exp(A*(1-beta/valueOf[betaop]));
	double wn = w0/(1.+valueOf[St]);
	double ratio = valueOf[Mf]/valueOf[Mx];
	double etam = 1.-2.59*ratio + 5.11*ratio*ratio - 3.52*ratio*ratio*ratio;
	double etas = 0.174*pow(valueOf[Se],-0.19);
	double Ir = Gamma*wn*valueOf[h]*etam*etas;
	double xsi = exp((1.+beta)*(0.792+0.681*pow(valueOf[sigma],0.5)));
	double rhob = w0/valueOf[delta];
	double epsilon = exp(-138/valueOf[sigma]);
	double Qig = 250.*beta + 1116*valueOf[Mf];
	double R0 = Ir*xsi/(rhob*epsilon*Qig);
	/* computing the wind factor */
	double C = 7.47*exp(-0.133*pow(valueOf[sigma], 0.55));
	double Ua = (valueOf[normalWind20ft] < 30. )?
			valueOf[a]*valueOf[normalWind20ft] : valueOf[a]*30.;
	double E = 0.715*exp(-0.000359*valueOf[sigma]);
	double phiw = C*pow(Ua, beta)*pow(beta/betaop, -E);
	/* computing the slope factor */
	double tanphi = valueOf[effectiveSlope];
	double phis = (tanphi > 0.)? 5.275*pow(beta, -0.3)*tanphi*tanphi : 0;
	/* computing the rate of spread (Rothermel 1972) */
	return R0*(1. + phiw + phis);
}

} /* namespace libforefire */
