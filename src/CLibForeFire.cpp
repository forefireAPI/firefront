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

#include "CLibForeFire.h"
#include "SimulationParameters.h"

using namespace std;

namespace libforefire {

Command executor;
static Command::Session* session = &(executor.currentSession);



Command* getLauncher(){
	return &executor;
}

void MNHInit(const double t){
	executor.setReferenceTime(t);


	/* Reading all the information on the parameters of ForeFire */
	ostringstream paramsfile;
	paramsfile<<SimulationParameters::GetInstance()->getParameter("caseDirectory")<<'/'
			<<SimulationParameters::GetInstance()->getParameter("ForeFireDataDirectory")<<'/'
			<<SimulationParameters::GetInstance()->getParameter("paramsFile");
	ifstream inputParams(paramsfile.str().c_str());
	if ( inputParams ) {
		string line;
		size_t numLine = 0;
		while ( getline( inputParams, line ) ) {
			numLine++;
			// checking for comments or newline
			if((line[0] == '#')||(line[0] == '*')||(line[0] == '\n'))
				continue;
			// treating the command of the current line
			executor.ExecuteCommand(line);
		}
	} else {
		cout<<"ERROR: File for ForeFire initialization not found !!"<<endl;
		cout<<'\t'<<"looked in "<<paramsfile.str()<<endl;
	}

	executor.outputDirs = SimulationParameters::GetInstance()->getParameterArray("atmoOutputDirectories");

}

void MNHCreateDomain(const int id
		, const int year, const int month
		, const int day, const double t
		, const double lat, const double lon
		, const int mdimx, const double* meshx
		, const int mdimy, const double* meshy
		, const int mdimz, const double* zgrid
		, const double dt){

	/* Defining the Fire Domain */
	if (session->fd) delete session->fd;

	session->fd = new FireDomain(id, year, month, day, t, lat, lon
			, mdimx, meshx, mdimy, meshy, mdimz, dt);

	executor.domain = session->fd;

	// A FireDomain has been created, the level is increased
	executor.increaseLevel();
	session->ff = session->fd->getDomainFront();
	// Defining the timetable of the events to be be in the domain
	if (session->tt) delete session->tt;
	session->tt = new TimeTable();
	// Associating this timetable to the domain
	session->fd->setTimeTable(session->tt);
	// Defining the simulator
	if (session->sim) delete session->sim;
	session->sim = new Simulator(session->tt, session->fd->outputs);

	/* Managing the outputs */
	ostringstream ffOutputsPattern;
	ffOutputsPattern<<SimulationParameters::GetInstance()->getParameter("caseDirectory")<<'/'
			<<SimulationParameters::GetInstance()->getParameter("fireOutputDirectory")<<'/'
			<<SimulationParameters::GetInstance()->getParameter("outputFiles")
			<<"."<<session->fd->getDomainID();
	SimulationParameters::GetInstance()->setParameter("ffOutputsPattern", ffOutputsPattern.str());
	session->outStrRep = new StringRepresentation(executor.domain);
	if ( SimulationParameters::GetInstance()->getInt("outputsUpdate") != 0 ){
		session->tt->insert(new FFEvent(session->outStrRep));
	}
// TODO	createNcFile(params->getParameter("ForeFireDataDirectory"), mdimx, mdimy, mdimz, meshx, meshy, zgrid);

	// Reading all the information on the initialization of ForeFire
	ostringstream initfile;
	if ( SimulationParameters::GetInstance()->getInt("parallelInit") != 1 ) {
		// Mono-file case: one file for all the processors
		initfile<<SimulationParameters::GetInstance()->getParameter("caseDirectory")<<'/'
								<<SimulationParameters::GetInstance()->getParameter("ForeFireDataDirectory")<<'/'
								<<SimulationParameters::GetInstance()->getParameter("InitFile");
	} else {
		// Multi-proc initial conditions
		initfile<<SimulationParameters::GetInstance()->getParameter("caseDirectory")<<'/'
								<<SimulationParameters::GetInstance()->getParameter("ForeFireDataDirectory")<<'/'
								<<SimulationParameters::GetInstance()->getParameter("InitFiles")
								<<"."<<id<<"."<<SimulationParameters::GetInstance()->getParameter("InitTime");
	}
	ifstream inputInit(initfile.str().c_str());
	if ( inputInit ) {
		string line;
		size_t numLine = 0;
		// skip the firest "firedomain" line for // init with multiple files
		if ( SimulationParameters::GetInstance()->getInt("parallelInit") == 1 ) getline( inputInit, line );
		while ( getline( inputInit, line ) ) {
			numLine++;
			// checking for comments or newline
			if((line[0] == '#')||(line[0] == '*')||(line[0] == '\n'))
				continue;
			// treating the command of the current line
			executor.ExecuteCommand(line);
		}
	} else {
		cout<<"File for ForeFire initialization "<<initfile.str()
				<<" not found !! No markers will be advected" << endl;
	}

	// completing the last front
	executor.completeFront(executor.currentSession.ff);

	// Managing the communication matrices for parallel computation
	if ( session->fd->parallel ) session->fd->createFirenodesMatrices();

	// advancing the simulation to the beginning of the atmospheric simulation
	double deltaT = session->fd->getSecondsFromReferenceTime(year, month, day, t);

	executor.setReferenceTime(deltaT);
	executor.setStartTime(deltaT);



}

void CheckLayer(const char* lname){
	string tmpname(lname);
	// searching for concerned layer
	FluxLayer<double>* myLayer = session->fd->getFluxLayer(tmpname);

	if ( myLayer == 0 ){
		if ( !session->fd->addFluxLayer(tmpname) ){
			cout<<"WARNING: layer for "<<tmpname
					<<" could not be found within ForeFire framework, "
					<<"this should cause serious problems when coupling"<<endl;
		}
	}
}

void MNHStep(double dt){
	ostringstream cmd;
	cmd << "step[dt=" << dt <<"]";
	string scmd = cmd.str();
	executor.ExecuteCommand(scmd);
}

void MNHGoTo(double time){
	ostringstream cmd;
	cmd << "goTo[t=" << time <<"]";
	string scmd = cmd.str();
	executor.ExecuteCommand(scmd);
}

void executeMNHCommand(const char* cmd){
	string scmd(cmd);
	executor.ExecuteCommand(scmd);
}

void FFPutString(const char* mname, char* str){
	// NOT TO BE CALLED
}

void FFGetString(const char* mname, const char* str){
	string name(mname);
	string val(str);
	SimulationParameters::GetInstance()->setParameter(name, val);
}

void FFPutInt(const char* mname, int* n){
	string name(mname);
	*n = SimulationParameters::GetInstance()->getInt(name);
}

void FFGetInt(const char* mname, int* n){
	string name(mname);
	SimulationParameters::GetInstance()->setInt(name, *n);
}

void FFPutIntArray(const char* mname, int* x,
		size_t sizein, size_t sizeout){
	// TODO
}

void FFGetIntArray(const char* mname, double time
		, int* x, int sizein, int sizeout){
	// TODO
}

void FFPutDouble(const char* mname, double* x){
	string name(mname);
	*x = SimulationParameters::GetInstance()->getDouble(name);
}

void FFGetDouble(const char* mname, double* x){
	string name(mname);
	SimulationParameters::GetInstance()->setDouble(name, *x);
}

void FFPutDoubleArray(const char* mname, double* x,
		size_t sizein, size_t sizeout){
	string tmpname(mname);
	// searching for concerned layer
	DataLayer<double>* myLayer = session->fd->getDataLayer(tmpname);
	if ( myLayer ){
		FFArray<double>* myMatrix;
		// getting the pointer
		myLayer->getMatrix(&myMatrix, executor.getTime());

		myMatrix->copyDataToFortran(x);
	} else {
		cout<<"Error trying to put data from unknown layer "<<tmpname<<endl;
	}
}

void FFGetDoubleArray(const char* mname, double t
		, double* x, size_t sizein, size_t sizeout){
	string tmpname(mname);
	double ct = executor.refTime + t;
	// searching for the layer to put data
	DataLayer<double>* myLayer = session->fd->getDataLayer(tmpname);
	if ( myLayer ){
		myLayer->setMatrix(tmpname, x, sizein, sizeout, ct);
	} else {
		cout<<"Error trying to get data for unknown layer "<<tmpname<<endl;
	}
}

void FFDumpDoubleArray(size_t nmodel, size_t nip, const char* mname, double t
		, double* x, size_t sizein, size_t ni, size_t nj, size_t nk, size_t sizeout){

	string tmpname(mname);
	ostringstream outputfile;
	double ct = executor.refTime + t;
	outputfile<<SimulationParameters::GetInstance()->getParameter("caseDirectory")<<'/'
			<<executor.outputDirs[nmodel-1]<<'/'
			<<SimulationParameters::GetInstance()->getParameter("outputFiles")
			<<"."<<nip<<"."<<tmpname;

	size_t niC = (size_t)(ni+0);
	size_t njC = (size_t)(nj+0);
	size_t nkC = (size_t)(nk+0);

	ofstream FileOut(outputfile.str().c_str(), ios_base::binary | ios::app);
	FileOut.write(reinterpret_cast<const char*>(&niC), sizeof(size_t));
	FileOut.write(reinterpret_cast<const char*>(&njC), sizeof(size_t));
	FileOut.write(reinterpret_cast<const char*>(&nkC), sizeof(size_t));
	FileOut.write(reinterpret_cast<const char*>(&ct), sizeof(double));

	size_t indF = 0;
	try {
		for ( indF = 0; indF < sizein; indF++ ) {
                     FileOut.write(reinterpret_cast<const char*>(x+indF), sizeof(double));
		}
	} catch (...) {
		cout << "Problem in passing a Fortran array to C array "
				<<tmpname<<endl;
	}

	FileOut.close();

}

#ifdef NETCDF_NOT_LEGACY
void saveNcRecord(int rec){cout << "CLibforefire:: saveNcRecord " << " newCDF Not Implemented" << endl;}
void createNcFile(string filename
		, const int& consted_ni, const int& consted_nj, const int& consted_nk
		, const double* meshx, const double* meshy, const double* zgrid){
			cout << "CLibforefire:: createNcFile " << " newCDF Not Implemented" << endl;
		}
#else

void saveNcRecord(int rec){
		size_t i = 0;
		size_t j = 0;
		size_t k = 0;
        ostringstream oss;

        oss<<session->params->getParameter("caseDirectory")<<'/'
        	<<session->params->getParameter("fireOutputDirectory")<<'/'
        	<<"data.nc";

		size_t   nx 	= 20;  // Nb de lignes de cells
		size_t   ny 	= 30;   // Nb de colonnes de cell
		size_t   nz 	= 40;   // Nb de colonnes de cell

		cout << "reopening "<< oss.str()<< endl;

		NcFile dataFile(oss.str().c_str(), NcFile::Write);
		 cout << "file opened Select"<< endl;
	   if (!dataFile.is_valid())
	   {
	      cout << "Couldn't open file!"<< endl;
	   }

		NcVar *cell_active = dataFile.get_var("values");

		NcDim* ctDim = dataFile.get_dim("time");
		  double  value[nz][ny][nx];

		   cout << "***  dumping record " <<oss.str()<< endl;

	//   cell_active->put(&value[0][0][0], nz, ny,nx);
       for (i = 0; i < nx ; i++){	for (j = 0; j < ny ; j++){for (k = 0; k < nz ; k++){value[k][j][i] = k;}}}
       cell_active->put_rec(ctDim,&value[0][0][0],0);
	   for (i = 0; i < nx ; i++){	for (j = 0; j < ny ; j++){for (k = 0; k < nz ; k++){value[k][j][i] = i;}}}
	   cell_active->put_rec(ctDim,&value[0][0][0],1);
	   for (i = 0; i < nx ; i++){	for (j = 0; j < ny ; j++){for (k = 0; k < nz ; k++){value[k][j][i] = j*i;}}}
	   cell_active->put_rec(ctDim,&value[0][0][0],2);
	   for (i = 0; i < nx ; i++){	for (j = 0; j < ny ; j++){for (k = 0; k < nz ; k++){value[k][j][i] = j;}}}
	   cell_active->put_rec(ctDim,&value[0][0][0],3);
	   for (i = 0; i < nx ; i++){	for (j = 0; j < ny ; j++){for (k = 0; k < nz ; k++){value[k][j][i] = 5;}}}
	   cell_active->put_rec(ctDim,&value[0][0][0],4);


	   cout << "*** SUCCESS dumping record " <<oss.str()<< endl;
	   dataFile.close();
}


void createNcFile(string filename
		, const int& consted_ni, const int& consted_nj, const int& consted_nk
		, const double* meshx, const double* meshy, const double* zgrid){
	ostringstream oss;

	size_t   ndim 	= 3;   // Nb de dimension

	// TODO transpose zgrid

	size_t i = 0;
	size_t j = 0;
	size_t k = 0;

	NcFile dataFile(filename.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);

	if (!dataFile.is_valid())
	{
		cout << "Couldn't open file for creation"<< endl;
	}

	NcVar *cell_geom = NULL;

	size_t ni = (size_t)consted_ni;
	size_t nj = (size_t)consted_nj;
	size_t nk = (size_t)consted_nk;

	NcDim* cxDim = dataFile.add_dim("ni", ni);
	NcDim* cyDim = dataFile.add_dim("nj", nj);
	NcDim* czDim = dataFile.add_dim("nk", nk);
	NcDim* cxyzDim = dataFile.add_dim("xyz", ndim);

	NcDim* ctDim = dataFile.add_dim("time");
	cout << "*** dims added " <<oss.str()<< endl;

	double  rv[nk][nj][ni][ndim];
	cout << "*** in Heap " <<oss.str()<< endl;

	for (i = 0; i < ni ; i++){
		for (j = 0; j < nj ; j++){
			for (k = 0; k < nk ; k++){
				rv[k][j][i][0] = meshx[i];
				rv[k][j][i][1] = meshy[j];
				rv[k][j][i][2] = zgrid[i + j*ni + k*ni*nj];
			}
		}
	}

	cout << "*** Geom " <<oss.str()<< endl;

	cell_geom = dataFile.add_var("grid_geometry", ncDouble, czDim, cyDim,cxDim,cxyzDim);
	cell_geom->put(&rv[0][0][0][0], nk, nj, ni, ndim);

	NcVar *cell_active = NULL;
	if ( !SimulationParameters::GetInstance()->isValued("outputs.variables") )
		cout<<"ERROR: vector of parameters outputs.variables should be valued"<<endl;
	vector<string> variables = SimulationParameters::GetInstance()->getParameterArray("outputs.variables");
	vector<string>::iterator variable;
	for ( variable = variables.begin(); variable != variables.end(); ++variable ){
		cell_active = dataFile.add_var(variable->c_str(), ncDouble,  ctDim,czDim, cyDim,cxDim);
	}

	cout << "*** SUCCESS creating " <<oss.str()<< endl;

	dataFile.close();
}
 
#endif



}
