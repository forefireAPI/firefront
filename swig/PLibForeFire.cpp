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


#include "PLibForeFire.h"

using namespace std;
using namespace libforefire;
#include <iostream>


Command* pyxecutor;
Command::Session* session ;
SimulationParameters* params ;



PLibForeFire::PLibForeFire(){
	pyxecutor = new Command();
	session = &(pyxecutor->currentSession);
	params = session->params;
}



void PLibForeFire::createDomain( int id
		,  int year,  int month
		,  int day,  double t
		,  double lat,  double lon
		,  int mdimx,  double* meshx
		,  int mdimy,  double* meshy
		,  int mdimz,  double* zgrid
		,  double dt){


	/* Defining the Fire Domain */
		if (session->fd) delete session->fd;

		session->fd = new FireDomain(id, year, month, day, t, lat, lon
				, mdimx, meshx, mdimy, meshy, mdimz, dt);

		pyxecutor->domain = session->fd;

		// A FireDomain has been created, the level is increased
		pyxecutor->increaseLevel();
		session->ff = session->fd->getDomainFront();
		// Defining the timetable of the events to be be in the domain
		if (session->tt) delete session->tt;
		session->tt = new TimeTable();
		// Associating this timetable to the domain
		session->fd->setTimeTable(session->tt);
		// Defining the simulator
		if (session->sim) delete session->sim;
		session->sim = new Simulator(session->tt, session->fd->outputs);


		session->outStrRep = new StringRepresentation(pyxecutor->domain);
		if ( SimulationParameters::GetInstance()->getInt("outputsUpdate") != 0 ){
			session->tt->insert(new FFEvent(session->outStrRep));
		}

		double deltaT = session->fd->getSecondsFromReferenceTime(year, month, day, t);

		pyxecutor->setReferenceTime(deltaT);
		pyxecutor->setStartTime(deltaT);


}





void PLibForeFire::addLayer(char *type, char* layername, char* keyname){

	pyxecutor->domain->addLayer(string(type),string(layername),string(keyname));
}

void PLibForeFire::setInt(char* name, int val){
	string lname(name);
	params->setInt(lname,val);
}
int PLibForeFire::getInt(char* name ){
	string lname(name);
	return params->getInt(lname);
}
void PLibForeFire::setDouble(char* name, double val){
	string lname(name);
	params->setDouble(lname,val);
}
double PLibForeFire::getDouble(char* name ){
	string lname(name);
	return params->getDouble(lname);
}
void PLibForeFire::setString(char* name, char* val){
	string lname(name);
	string lval(val);
		params->setParameter(lname,val);
}

string PLibForeFire::getString(char *name)
{

	return params->getParameter(string(name));
}
string PLibForeFire::execute(char *command)
{
	ostringstream stringOut;
	pyxecutor->setOstringstream(&stringOut);
	string smsg(command);
	pyxecutor->ExecuteCommand(smsg);
	return stringOut.str();
}


void PLibForeFire::addScalarLayer(char *type, char *name, double x0 , double y0, double t0, double width , double height, double timespan, int nnx, int nny, int nnz,  double* values){

	FFPoint *p0 = new FFPoint(x0,y0,0);
	FFPoint *pe = new FFPoint(width,height,0);
	string lname(name);
	string ltype(type);

	size_t ni = nnx;
	size_t nj = nny;
	size_t nk = nnz;
	size_t nl = 1;

 	pyxecutor->domain->addScalarLayer(type, lname,  x0,y0,  t0,                    width ,            height,  timespan,       ni,   nj,  nk, nl,values);

}

void PLibForeFire::addIndexLayer(char *type,char *name, double x0 , double y0, double t0, double width , double height, double timespan, int nnx, int nny, int nnz,  int* values){
	FFPoint *p0 = new FFPoint(x0,y0,0);
	FFPoint *pe = new FFPoint(width,height,0);
	string lname(name);
	string ltype(type);

	size_t ni = nnx;
	size_t nj = nny;
	size_t nk = nnz;
	size_t nl = 1;

 	pyxecutor->domain->addIndexLayer(ltype, lname,  x0,y0,  t0,                    width ,            height,  timespan,       ni,   nj,  nk, nl,values);
}

void PLibForeFire::getDoubleArray(char* name, double** outA, int* outNI, int* outNJ, int* outNK){
	double lTime = pyxecutor->domain->getSimulationTime();


	PLibForeFire::getDoubleArray(name,lTime,outA,outNI,  outNJ,  outNK);
}

void PLibForeFire::getDoubleArray(char* name, double t, double** outA, int* outNI, int* outNJ, int* outNK){
	string lname(name);
	FluxLayer<double>* myFluxLayer = pyxecutor->domain->getFluxLayer(lname);

		if ( myFluxLayer ){
			FFArray<double>* srcD;
			myFluxLayer->getMatrix(&srcD, t);
			if (srcD != 0) {
				*outNI = srcD->getDim("x");
				*outNJ = srcD->getDim("y");
				*outNK = srcD->getDim("z");
				*outA =  srcD->getData();
				return;
			}
		}

	DataLayer<double>* myDataLayer = pyxecutor->domain->getDataLayer(lname);

		if ( myDataLayer ){
			FFArray<double>* srcD;
			myDataLayer->getMatrix(&srcD, t);
			if (srcD != 0) {
				*outNI = srcD->getDim("x");
				*outNJ = srcD->getDim("y");
				*outNK = srcD->getDim("z");
				*outA =  srcD->getData();
				return;
			}
		}

		*outNI = 0;
		*outNJ = 0;
		*outNK = 0;
		*outA =  NULL;
		return ;
}

