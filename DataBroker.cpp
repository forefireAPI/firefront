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

#include "DataBroker.h"
#include "MultiplicativeLayer.h"
#include "FireDomain.h"

namespace libforefire {

const DataBroker::propGetterMap DataBroker::propPropertiesGetters =
		DataBroker::makePGmap();
const DataBroker::fluxGetterMap DataBroker::fluxPropertiesGetters =
		DataBroker::makeFGmap();

DataLayer<double>* DataBroker::fuelLayer = 0;
DataLayer<double>* DataBroker::moistureLayer = 0;
DataLayer<double>* DataBroker::dummyLayer = 0;
DataLayer<double>* DataBroker::altitudeLayer = 0;
DataLayer<double>* DataBroker::slopeLayer = 0;
DataLayer<double>* DataBroker::windULayer = 0;
DataLayer<double>* DataBroker::windVLayer = 0;

XYZTDataLayer<double>* DataBroker::PwindULayer = 0;
XYZTDataLayer<double>* DataBroker::PwindVLayer = 0;
FluxLayer<double>* DataBroker::heatFluxLayer = 0;

SimulationParameters* DataBroker::params = SimulationParameters::GetInstance();
double frontScanDistance = 1000;

DataBroker::DataBroker(FireDomain* fd) :
		domain(fd) {
	commonInitialization();
}

DataBroker::~DataBroker() {
	// Deleting layers
	while (!layers.empty()) {
		delete layers.back();
		layers.pop_back();
	}
	// Deleting atmospheric and parallel data
	if (atmosphericData)
		delete atmosphericData;
	if (parallelData)
		delete parallelData;
	// Deleting optimized data brokers
	delete[] optimizedPropDataBroker;
	delete[] optimizedFluxDataBroker;
}

void DataBroker::commonInitialization() {
	/* atmospheric data */
	atmosphericData = new AtmosphericData();
	/* parallel data */
	parallelData = new ParallelData();
	/* data brokers for propagation models */
	propDataGetters.resize(FireDomain::NUM_MAX_PROPMODELS);
	numPropDataGetters.resize(FireDomain::NUM_MAX_PROPMODELS);
	optimizedPropDataBroker = new bool[FireDomain::NUM_MAX_PROPMODELS];
	for (size_t i = 0; i < FireDomain::NUM_MAX_PROPMODELS; i++) {
		numPropDataGetters[i] = 0;
		optimizedPropDataBroker[i] = true;
	}
	frontScanDistance = params->getDouble("frontScanDistance");
	/* data brokers for flux models */
	fluxDataGetters.resize(FireDomain::NUM_MAX_FLUXMODELS);
	numFluxDataGetters.resize(FireDomain::NUM_MAX_FLUXMODELS);
	optimizedFluxDataBroker = new bool[FireDomain::NUM_MAX_FLUXMODELS];
	for (size_t i = 0; i < FireDomain::NUM_MAX_FLUXMODELS; i++) {
		numFluxDataGetters[i] = 0;
		optimizedFluxDataBroker[i] = true;
	}

	/* Getting the fuel parameters' table */
	/*------------------------------------*/
	ostringstream infile;
    /*
	infile << params->getParameter("caseDirectory") << '/'
			<< params->getParameter("ForeFireDataDirectory") << '/'
			<< params->getParameter("fuelsTableFile");
    */
    infile << params->GetPath(params->getParameter("fuelsTableFile"));
	readTableFromAsciiFile(infile.str(), fuelPropertiesTable);
}

void DataBroker::updateFuelValues(PropagationModel* model, string key, double value){

	size_t curfuel = 0;
	map<string, double>::iterator iindex;
	map<string, double>::iterator iparam;
	// second dimension
	for (curfuel = 0; curfuel < fuelPropertiesTable.size(); curfuel++) {
		// Getting the index of the fuel
		iindex = fuelPropertiesTable[curfuel].find(key);
		if(iindex != fuelPropertiesTable[curfuel].end()){
			iindex->second = value;

		}

	}



	extractFuelProperties(fuelPropertiesTable, model);




}


void DataBroker::registerPropagationModel(PropagationModel* model) {
	/* Constructing the vector of property getters */

	propGetterMap::const_iterator pg;
	bool fuelAsked = false;


	for (size_t prop = 0; prop < model->numProperties; prop++) {
		try {
			if ((model->wantedProperties)[prop].substr(0, 4) == "fuel") {

				if (!fuelAsked) {
					propDataGetters[model->index].push_back(&getFuelProperties);
					numPropDataGetters[model->index]++;
					neededProperties.push_back("fuel");
					fuelAsked = true;
				}
			} else {
				pg = propPropertiesGetters.find(
						(model->wantedProperties)[prop]);
				if (pg == propPropertiesGetters.end()) {
					optimizedPropDataBroker[model->index] = false;

				} else {
					propDataGetters[model->index].push_back(pg->second);
					numPropDataGetters[model->index]++;
				}
				if (find(neededProperties.begin(), neededProperties.end(),
						(model->wantedProperties)[prop])
						== neededProperties.end())
					neededProperties.push_back((model->wantedProperties)[prop]);
			}
		} catch (...) {
			cout
					<< "Encountered an error while constructing the function to compute "
					<< (model->wantedProperties)[prop] << endl;
		}
	}

	/* registering the prop model in the fire domain */
	domain->registerPropagationModel(model->index, model);

	/* constructing the table of fuel parameters values */
	if (model->numFuelProperties > 0)
		extractFuelProperties(fuelPropertiesTable, model);

}


void DataBroker::registerFluxModel(FluxModel* model) {

	/* Constructing the vector of property getters */
	fluxGetterMap::const_iterator pg;
	bool fuelAsked = false;
	for (size_t prop = 0; prop < model->numProperties; prop++) {
		try {
			if ((model->wantedProperties)[prop].substr(0, 4) == "fuel") {
				if (!fuelAsked) {
					fluxDataGetters[model->index].push_back(&getFuelProperties);
					numFluxDataGetters[model->index]++;
					if (find(neededProperties.begin(), neededProperties.end(),
							"fuel") == neededProperties.end())
						neededProperties.push_back("fuel");
					fuelAsked = true;
				}
			} else {
				pg = fluxPropertiesGetters.find(
						(model->wantedProperties)[prop]);
				if (pg == fluxPropertiesGetters.end()) {
					cout << "could not find an optimized property getter for "
							<< (model->wantedProperties)[prop] << endl;
					cout << "databroker is switched to an un-optimized mode"
							<< endl;
					optimizedFluxDataBroker[model->index] = false;
				} else {
					fluxDataGetters[model->index].push_back(pg->second);
					numFluxDataGetters[model->index]++;
				}
				if (find(neededProperties.begin(), neededProperties.end(),
						(model->wantedProperties)[prop])
						== neededProperties.end())
					neededProperties.push_back((model->wantedProperties)[prop]);
			}
		} catch (...) {
			cout
					<< "Encountered an error while constructing the function to compute "
					<< (model->wantedProperties)[prop] << endl;
		}
	}

	/* registering the flux model in the fire domain */
	domain->registerFluxModel(model->index, model);

	/* constructing the table of fuel parameters values */
	if (model->numFuelProperties > 0)
		extractFuelProperties(fuelPropertiesTable, model);

}

void DataBroker::extractFuelProperties(vector<map<string, double> > propsTable,
		ForeFireModel* model) {

	try {
		delete(model->fuelPropertiesTable);
		model->fuelPropertiesTable = new FFArray<double>("fuelProperties", 0.,
				FuelDataLayer<double>::MAXNUMFUELS, model->numFuelProperties);
		size_t curfuel = 0;
		map<string, double>::iterator iindex;
		map<string, double>::iterator iparam;
		// second dimension
		for (curfuel = 0; curfuel < propsTable.size(); curfuel++) {
			// Getting the index of the fuel
			iindex = propsTable[curfuel].find("Index");
			size_t ind = (size_t) iindex->second;
			// Getting the values of the wanted parameters
			for (size_t param = 0; param < model->numFuelProperties; param++) {
				iparam = propsTable[curfuel].find(
						(model->fuelPropertiesNames)[param]);
				if (iparam != propsTable[curfuel].end()) {
					(*(model->fuelPropertiesTable))(ind, param) =
							iparam->second;
				} else {
					cout << "Parameter " << (model->fuelPropertiesNames)[param]
							<< " could not be found " << "for fuel " << ind
							<< " in the fuel table" << endl;
				}
			}
		}
	} catch (const bad_alloc &) {
		// deleting what has been allocated
		cout << "ERROR : while retrieving the table of wanted fuel parameters"
				<< " in FuelDataLayer<T>::extractParameters" << endl;
	}
}

void DataBroker::registerLayer(string name, DataLayer<double>* layer) {


	/* inserting the layer into the map of layers */
	//cout << "registering " << name << " !" << endl;
	ilayer = layersMap.find(name);
	if (ilayer != layersMap.end()) {

		if (params->getParameter("runmode") == "standalone") {
			cout << "Redefining layer for variable " << name << " !" << endl;
			DataLayer<double>* oldlayer = ilayer->second;
			layersMap.erase(name);
			layers.remove(oldlayer);
			delete oldlayer;
		}
	}
	layersMap.insert(make_pair(name, layer));
	layers.push_back(layer);

	/* looking for possible match with predefined layers */
	if (name.find("altitude") != string::npos) {
		altitudeLayer = layer;
		// along with the topography comes the slope
		slopeLayer = new GradientDataLayer<double>("slope", altitudeLayer,
				params->getDouble("spatialIncrement"));
		registerLayer("slope", slopeLayer);
	}
	if (name.find("moisture") != string::npos)
		moistureLayer = layer;
	if (name.find("windU") != string::npos)
		windULayer = layer;
	if (name.find("windV") != string::npos)
		windVLayer = layer;
	if (name.find("fuel") != string::npos){
		fuelLayer = layer;
	}

}

void DataBroker::registerFluxLayer(string name, FluxLayer<double>* layer) {

	flayer = fluxLayersMap.find(name);
	if (flayer != fluxLayersMap.end()) {
		FluxLayer<double>* oldlayer = flayer->second;
		fluxLayersMap.erase(name);
		fluxLayers.remove(oldlayer);
		delete oldlayer;
	}
	fluxLayersMap.insert(make_pair(name, layer));
	fluxLayers.push_back(layer);
	registerLayer(name, layer);

	/* looking for possible match with predefined layers */
	if (name.find("heatFlux") != string::npos){

		heatFluxLayer = layer;
	}
}

void DataBroker::setAtmosphericDomain(const FFPoint& SWCorner,
		const FFPoint& NECorner, const size_t& nx, const size_t& ny) {
	atmoSWCorner = SWCorner;
	atmoNECorner = NECorner;
	atmosphericNx = nx;
	atmosphericNy = ny;
}

void DataBroker::initializeAtmosphericLayers(const double& time,
		const size_t& bnx, const size_t& bny) {

	atmosphericData->setSize(atmosphericNx, atmosphericNy);

	double dx = (atmoNECorner.getX() - atmoSWCorner.getX()) / atmosphericNx;
	double dy = (atmoNECorner.getY() - atmoSWCorner.getY()) / atmosphericNy;

	// Loading the atmospheric layers
	FFPoint windUOrigin = FFPoint(atmoSWCorner.getX() - dx,
			atmoSWCorner.getY() - 0.5 * dy);
	TwoTimeArrayLayer<double>* wul = new TwoTimeArrayLayer<double>("windU",
			atmosphericData->windU, time, atmosphericData->oldWindU, time,
			windUOrigin, dx, dy);
	registerLayer("windU", wul);
	registerLayer("outerWindU", wul);
	FFPoint windVOrigin = FFPoint(atmoSWCorner.getX() - 0.5 * dx,
			atmoSWCorner.getY() - dy);
	TwoTimeArrayLayer<double>* wvl = new TwoTimeArrayLayer<double>("windV",
			atmosphericData->windV, time, atmosphericData->oldWindV, time,
			windVOrigin, dx, dy);
	registerLayer("windV", wvl);
	registerLayer("outerWindV", wvl);

	// Loading the topography
	Array2DdataLayer<double>* alt = new Array2DdataLayer<double>("altitude",
			atmosphericData->topography, atmoSWCorner, atmoNECorner);

	registerLayer("altitude", alt);
}

void DataBroker::initializePropagativeLayer(string filename) {

	NcFile* NcdataFile = new NcFile(filename.c_str(), NcFile::ReadOnly);
	if (NcdataFile->is_valid()) {
		const char* type = "type";
		string layerType;
		NcAtt* att;
		// Loading the propagative layer from a netCDF file
		for (int layer = 0; layer < NcdataFile->num_vars(); layer++) {
			string varName(NcdataFile->get_var(layer)->name());
			att = NcdataFile->get_var(layer)->get_att(type);
			att == 0 ? layerType.erase() : layerType.assign(att->as_string(0));
			if (layerType.find("propagative") != string::npos) {
				PropagativeLayer<double>* propLayer =
						constructPropagativeLayerFromFile(NcdataFile, varName);
				domain->setPropagativeLayer(propLayer);
			}
			delete att;
			att = 0;
		}
	} else {
		cout << "WARNING: could not find netCDF file " << filename << endl;
	}
	if (NcdataFile)
		delete NcdataFile;

	if (domain->getPropagativeLayer() == 0) {
		if (!domain->addPropagativeLayer(
				params->getParameter("propagationModel"))) {
			cout
					<< "PROBLEM: it was not possible to retrieve propagation model "
					<< params->getParameter("propagationModel") << endl;
		}
	}

	if (domain->getPropagativeLayer() != 0) registerLayer(domain->getPropagativeLayer()->getKey(), domain->getPropagativeLayer());

}

void DataBroker::initializeFluxLayers(string filename) {

	NcFile* NcdataFile = new NcFile(filename.c_str(), NcFile::ReadOnly);
	if (NcdataFile->is_valid()) {
		const char* type = "type";
		string layerType;
		NcAtt* att;
		// Loading the flux layers from a netCDF file
		for (int layer = 0; layer < NcdataFile->num_vars(); layer++) {
			string varName(NcdataFile->get_var(layer)->name());
			att = NcdataFile->get_var(layer)->get_att(type);
			att == 0 ? layerType.erase() : layerType.assign(att->as_string(0));
			if (layerType.find("flux") != string::npos) {
				FluxLayer<double>* newlayer = constructFluxLayerFromFile(
						NcdataFile, varName);
				registerFluxLayer(varName, newlayer);
			}
			delete att;
			att = 0;
		}
	} else {
		cout << "WARNING: could not find netCDF file " << filename << endl;
	}

	if (NcdataFile != 0)
		delete NcdataFile;

}

void DataBroker::initializeParallelProperties(const size_t& nx,
		const size_t& ny, const size_t& nz, const double& initval) {

	parallelData->setSize(nx, ny, nz, initval);
	// Loading the parallel layers
	Array3DdataLayer<double>* newlayer = new Array3DdataLayer<double>("posX",
			parallelData->FireNodesInCellPosX, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesPosX", newlayer);
	newlayer = new Array3DdataLayer<double>("posY",
			parallelData->FireNodesInCellPosY, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesPosY", newlayer);
	newlayer = new Array3DdataLayer<double>("velX",
			parallelData->FireNodesInCellVelX, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesVelX", newlayer);
	newlayer = new Array3DdataLayer<double>("velY",
			parallelData->FireNodesInCellVelY, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesVelY", newlayer);
	newlayer = new Array3DdataLayer<double>("time",
			parallelData->FireNodesInCellTime, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesTime", newlayer);
	newlayer = new Array3DdataLayer<double>("id",
			parallelData->FireNodesInCellId, atmoSWCorner, atmoNECorner);
	registerLayer("FireNodesId", newlayer);
}

void DataBroker::addConstantLayer(string name, const double& val) {
	/* Adding a simple constant layer for variable name */
	XYZTDataLayer<double>* newLayer = new XYZTDataLayer<double>(name, val);
	registerLayer(name, newLayer);
}



void DataBroker::loadFromNCFile(string filename) {

	/* Loading all the properties in the netCDF file */
	NcFile* NcdataFile = new NcFile(filename.c_str(), NcFile::ReadOnly);
	propGetterMap::const_iterator pg;
	if (NcdataFile->is_valid()) {
		const char* type = "type";
		/* loading the properties except for the fuel */
		for (int layer = 0; layer < NcdataFile->num_vars(); layer++) {
			string varName(NcdataFile->get_var(layer)->name());
			NcAtt* att;
            att = NcdataFile->get_var(layer)->get_att(type);
			string layerType(att->as_string(0));
			delete att;
			if (varName == "wind") {
								// wind is given by the NetCDF file

								XYZTDataLayer<double> * wul = constructXYZTLayerFromFile(
										NcdataFile, varName.c_str(),0);
								registerLayer("windU", wul);
								XYZTDataLayer<double> * wvl = constructXYZTLayerFromFile(
																		NcdataFile, varName.c_str(),1);
								registerLayer("windV", wvl);
								PwindULayer = wul;
								PwindVLayer = wvl;

							}
			pg = propPropertiesGetters.find(varName);

			if (pg != propPropertiesGetters.end()) {

				if (varName == "altitude") {
					// altitude is given by the NetCDF file
					XYZTDataLayer<double> * alt = constructXYZTLayerFromFile(
							NcdataFile, varName.c_str(),-1);

					registerLayer(varName, alt);
				/*	XYZTDataLayer<double> * valt = constructXYZTLayerFromFile(
							NcdataFile, varName.c_str());

					registerLayer("lavadepot", valt);*/
				} else if (varName == "windU") {
					// wind is given by the NetCDF file
					XYZTDataLayer<double> * wul = constructXYZTLayerFromFile(
							NcdataFile, varName.c_str(),-1);
					registerLayer(varName, wul);
				} else if (varName == "windV") {
					// wind is given by the NetCDF file
					XYZTDataLayer<double> * wvl = constructXYZTLayerFromFile(
							NcdataFile, varName.c_str(),-1);
					registerLayer(varName, wvl);
				} else if (varName == "fieldSpeed") {
					dummyLayer = constructXYZTLayerFromFile(NcdataFile,
							varName.c_str(),-1);
					registerLayer(varName, dummyLayer);
				}
			} else if (layerType.find("data") != string::npos) {

				if (find(neededProperties.begin(), neededProperties.end(),
						varName) != neededProperties.end()) {
					XYZTDataLayer<double>* newlayer =
							constructXYZTLayerFromFile(NcdataFile,
									varName.c_str(),-1);
					registerLayer(varName, newlayer);
				}
			} else if (layerType.find("parameter") != string::npos) {

				/* loading the parameters */
				NcVar* ncparams = NcdataFile->get_var(varName.c_str());
				int numParams = (size_t) ncparams->num_atts();
				NcAtt* ncparam;
				for (int i = 0; i < numParams; i++) {
					ncparam = ncparams->get_att(i);
					string pname(ncparam->name());
					string sval(ncparam->as_string(0));
					params->setParameter(pname, sval);
				}
				if (ncparam != 0)
					delete ncparam;
			}
		}

		/* loading the fuel layer */
		fuelLayer = constructFuelLayerFromFile(NcdataFile);
		registerLayer("fuel", fuelLayer);

	}

}

void DataBroker::dropProperty(string name) {
	vector<string>::iterator it = find(neededProperties.begin(),
			neededProperties.end(), name);
	if (it != neededProperties.end())
		neededProperties.erase(it);
}

void DataBroker::insureLayersExistence() {
	/* checking the layers needed for the propagation and flux models */
	while (!neededProperties.empty()) {
		if (getLayer(neededProperties.back()) == 0) {
			if (neededProperties.back().find("normalWind") != string::npos) {
				if (windULayer == 0 or windVLayer == 0)
					cout << "layer for normal wind doesn't rely on existing"
							<< " windU and windV layers, creating them if needed"
							<< endl;
				/* special treatment for normal wind */
				if (windULayer == 0) {
					double val = 0.;
					if (params->isValued("windU"))
						val = params->getDouble("windU");
					addConstantLayer("windU", val);
				}
				if (windVLayer == 0) {
					double val = 0.;
					if (params->isValued("windV"))
						val = params->getDouble("windV");
					addConstantLayer("windV", val);
				}
			} else if (neededProperties.back().find("uel") != string::npos) {
				cout
						<< "layer for fuel doesn't rely on existing fuel layer, creating it"
						<< endl;
				int fuelIndex = 0;
				if (params->isValued("fuel"))
					fuelIndex = params->getDouble("fuel");
				fuelLayer = new FuelDataLayer<double>("fuel", fuelIndex);
			} else if (neededProperties.back().find("slope") != string::npos
					and altitudeLayer == 0) {
				cout
						<< "layer for slope doesn't rely on existing altitude layer, creating it"
						<< endl;
				double val = 0.;
				if (params->isValued("altitude"))
					val = params->getDouble("altitude");
				addConstantLayer("altitude", val);
			} else if (neededProperties.back().find("epth") != string::npos) {
				int fdepth = 1;
				FireNode::setFrontDepthComputation(fdepth);
				/* checking that a heat flux layer is present for burn checks.
				 * If not, instantiating a HeatFluxBasicModel */
				if (heatFluxLayer == 0 and !domain->addFluxLayer("heatFlux")) {
					cout
							<< "WARNING: heat flux layer could not be found within ForeFire framework, this should"
							<< " cause serious problems for the simulation as front depth is required"
							<< endl;
				}

			} else if (neededProperties.back().find("urvature")
					!= string::npos) {
				int curv = 1;
				FireNode::setCurvatureComputation(curv);
			} else {
				/* else instantiating a constant layer for the property */
				double val = 0.;
				if (params->isValued(neededProperties.back()))
					val = params->getDouble(neededProperties.back());
				addConstantLayer(neededProperties.back(), val);
			}
		}
		/* erasing the property from the list of properties to be treated */
		dropProperty(neededProperties.back());
	}

	/* checking the layers that are always needed */
	if (altitudeLayer == 0) {
		double val = 0.;
		if (params->isValued("altitude"))
			val = params->getDouble("altitude");
		addConstantLayer("altitude", val);
	}

}

void DataBroker::initFluxLayers(const double& t) {
	list<FluxLayer<double>*>::iterator flay;
	for (flay = fluxLayers.begin(); flay != fluxLayers.end(); ++flay)
		(*flay)->setFirstCall(t);
}

bool DataBroker::isRelevantData(FFPoint& SW, FFPoint& ext) {
	if (SW.getX() > domain->NECornerX()
			or SW.getX() + ext.getX() < domain->SWCornerX()
			or SW.getY() > domain->NECornerY()
			or SW.getY() + ext.getY() < domain->SWCornerY())
		return false;
	return true;
}

double DataBroker::getNetCDFFileVersion(NcVar* var) {
	if(true) return 1;
	NcAtt* xsw = var->get_att("version");
	double version = xsw == 0 ? -1 : xsw->as_double(0);
	delete xsw;
	return version;
}

FFPoint DataBroker::getNetCDFSWCorner(NcVar* var) {
	const char* SWx = "SWx";
	const char* SWy = "SWy";
	const char* SWz = "SWz";
	NcAtt* xsw = var->get_att(SWx);
	double Xorigin = xsw == 0 ? 0 : xsw->as_double(0);
	NcAtt* ysw = var->get_att(SWy);
	double Yorigin = ysw == 0 ? 0 : ysw->as_double(0);
	NcAtt* zsw = var->get_att(SWz);
	double Zorigin = zsw == 0 ? 0 : zsw->as_double(0);
	delete xsw;
	delete ysw;
	delete zsw;
	FFPoint relSWC = FFPoint(Xorigin, Yorigin, Zorigin);
	FFPoint shiftPos = FFPoint(params->getDouble("SHIFT_ALL_DATA_ABSCISSA_BY"),
			params->getDouble("SHIFT_ALL_DATA_ORDINATES_BY"));
	return shiftPos + relSWC;
}

double DataBroker::getNetCDFTimeOrigin(NcVar* var) {
	const char* t0 = "t0";
	NcAtt* attt0 = var->get_att(t0);
	double timeOrigin = attt0 == 0 ? 0 : attt0->as_double(0);
	delete attt0;
	return timeOrigin;
}

FFPoint DataBroker::getNetCDFSpatialSpan(NcVar* var) {
	const char* llx = "Lx";
	const char* lly = "Ly";
	const char* llz = "Lz";
	NcAtt* attlx = var->get_att(llx);
	double Lx = attlx == 0 ? numeric_limits<double>::infinity() : attlx->as_double(0);
	NcAtt* attly = var->get_att(lly);
	double Ly = attly == 0 ? numeric_limits<double>::infinity() : attly->as_double(0);
	NcAtt* attlz = var->get_att(llz);
	double Lz = attlz == 0 ? numeric_limits<double>::infinity() : attlz->as_double(0);
	delete attlx;
	delete attly;
	delete attlz;
	return FFPoint(Lx, Ly, Lz);
}

double DataBroker::getNetCDFTimeSpan(NcVar* var) {
	const char* llt = "Lt";
	NcAtt* attlt = var->get_att(llt);
	double Lt = attlt == 0 ? numeric_limits<double>::infinity() : attlt->as_double(0);
	delete attlt;
	return Lt;
}

XYZTDataLayer<double>* DataBroker::constructXYZTLayerFromFile(
		NcFile* NcdataFile, string property, int dimTSelected = -1) {

	/* Getting the information on the mesh */
	/*-------------------------------------*/
	// Getting information on the extension of the domain
	const char* sdom = "domain";
	NcVar* domvar = NcdataFile->get_var(sdom);
	FFPoint SWCorner = getNetCDFSWCorner(domvar);
	FFPoint spatialExtent = getNetCDFSpatialSpan(domvar);
	double timeOrigin = getNetCDFTimeOrigin(domvar);
	double Lt = getNetCDFTimeSpan(domvar);
	double version = getNetCDFFileVersion(domvar);

	/* Getting the desired variable */
	/*------------------------------*/
	NcVar* values = NcdataFile->get_var(property.c_str());
	if (values->num_dims() > 4) {
		cout << "Variable " << property << " doesn't have the right dimension ("
				<< values->num_dims() << " instead of 4 maximum)" << endl;
	}

	size_t nx = 0;
	size_t ny = 0;
	size_t nz = 0;
	size_t nt = 0;

	if (values->num_dims() == 2) {
		ny = (size_t) values->get_dim(0)->size();
		nx = (size_t) values->get_dim(1)->size();
	}
	if (values->num_dims() == 3) {
		nz = (size_t) values->get_dim(0)->size();
		ny = (size_t) values->get_dim(1)->size();
		nx = (size_t) values->get_dim(2)->size();

	}
	if (values->num_dims() == 4) {
		nt = (size_t) values->get_dim(0)->size();
		nz = (size_t) values->get_dim(1)->size();
		ny = (size_t) values->get_dim(2)->size();
		nx = (size_t) values->get_dim(3)->size();
	}

	/* Getting the data */
	/*------------------*/

		 if (dimTSelected>-1){
			 	 nt = 1;

		    }
	double* data = readAndTransposeFortranProjectedField(values, nt, nz, ny, nx,
			version>0,dimTSelected);

	if (isRelevantData(SWCorner, spatialExtent)) {

		/* Instanciating the data layer */
		/*------------------------------*/

		XYZTDataLayer<double>* newlayer = new XYZTDataLayer<double>(property,
				SWCorner, timeOrigin, spatialExtent, Lt, nx, ny, nz, nt, data);
		delete[] data;
				return newlayer;



	} else {

		cout << "WARNING, spatial domain of validity for variable " << property
				<< ": " << SWCorner.print() << " -> ("
				<< SWCorner.getX() + spatialExtent.getX() << ","
				<< SWCorner.getY() + spatialExtent.getY() << ","
				<< SWCorner.getZ() + spatialExtent.getZ()
				<< ") has no intersection with the fire domain: ("
				<< domain->SWCornerX() << "," << domain->SWCornerY() << ") -> ("
				<< domain->NECornerX() << "," << domain->NECornerY() << ")"
				<< endl << "...resizing to the fire domain..." << endl;

		FFPoint SWC = FFPoint(domain->SWCornerX(), domain->SWCornerY(), 0.);
		FFPoint ext = FFPoint(domain->NECornerX() - domain->SWCornerX(),
				domain->NECornerY() - domain->SWCornerY(), 10000.);
		double t0 = 0;
		double Dt = 10000000.;
		XYZTDataLayer<double>* newlayer = new XYZTDataLayer<double>(property,
				SWC, t0, ext, Dt, nx, ny, nz, nt, data);
		delete[] data;
		return newlayer;

	}
}

FuelDataLayer<double>* DataBroker::constructFuelLayerFromFile(
		NcFile* NcdataFile) {

	/* Getting the information on the mesh */
	/*-------------------------------------*/
	// Getting information on the extension of the domain
	const char* sdom = "domain";
	NcVar* domVar = NcdataFile->get_var(sdom);
	FFPoint SWCorner = getNetCDFSWCorner(domVar);
	FFPoint spatialExtent = getNetCDFSpatialSpan(domVar);
	double timeOrigin = getNetCDFTimeOrigin(domVar);
	double Lt = getNetCDFTimeSpan(domVar);
    double version = getNetCDFFileVersion(domVar);
	/* Getting the fuel Map */
	/*----------------------*/
	string fmap = "fuel";
	delete domVar;

	NcVar* values = NcdataFile->get_var(fmap.c_str());

	if (values->num_dims() > 4) {
		cout << "Variable " << fmap << " doesn't have the right dimension ("
				<< values->num_dims() << " instead of 4 maximum)" << endl;
	}

	size_t nx = 0;
	size_t ny = 0;
	size_t nz = 0;
	size_t nt = 0;
	if (values->num_dims() == 2) {
		ny = (size_t) values->get_dim(0)->size();
		nx = (size_t) values->get_dim(1)->size();
	}
	if (values->num_dims() == 3) {
		nz = (size_t) values->get_dim(0)->size();
		ny = (size_t) values->get_dim(1)->size();
		nx = (size_t) values->get_dim(2)->size();

	}
	if (values->num_dims() == 4) {
		nt = (size_t) values->get_dim(0)->size();
		nz = (size_t) values->get_dim(1)->size();
		ny = (size_t) values->get_dim(2)->size();
		nx = (size_t) values->get_dim(3)->size();
	}

	/* Getting the data */
	/*------------------*/
	int* fuelMap = readAndTransposeIntFortranProjectedField(values, nt, nz, ny,
			nx, version>0,-1);

	if (isRelevantData(SWCorner, spatialExtent)) {

		FuelDataLayer<double>* newlayer = new FuelDataLayer<double>(fmap,
				SWCorner, timeOrigin, spatialExtent, Lt, nx, ny, nz, nt,
				fuelMap);
		delete[] fuelMap;
		return newlayer;

	} else {

		cout << "WARNING, spatial domain of validity for fuel table: "
				<< SWCorner.print() << " -> ("
				<< SWCorner.getX() + spatialExtent.getX() << ","
				<< SWCorner.getY() + spatialExtent.getY() << ","
				<< SWCorner.getZ() + spatialExtent.getZ()
				<< ") has no intersection with the fire domain: ("
				<< domain->SWCornerX() << "," << domain->SWCornerY() << ") -> ("
				<< domain->NECornerX() << "," << domain->NECornerY() << ")"
				<< endl << "...resizing to the fire domain..." << endl;

		FFPoint SWC = FFPoint(domain->SWCornerX(), domain->SWCornerY(), 0.);
		FFPoint ext = FFPoint(domain->NECornerX() - domain->SWCornerX(),
				domain->NECornerY() - domain->SWCornerY(), 10000.);
		double t0 = 0;
		double Dt = 10000000.;
		FuelDataLayer<double>* newlayer = new FuelDataLayer<double>(fmap, SWC,
				t0, ext, Dt, nx, ny, nz, nt, fuelMap);
		delete[] fuelMap;
		return newlayer;

	}

	return 0;
}

void DataBroker::readTableFromAsciiFile(string filename,
		vector<map<string, double> >& table) {


	vector<string> paramNames;

	/* Opening the file */
	ifstream file(filename.c_str());

	if (!file) {
		cout << "WARNING: could not load file for fuel properties " << filename
				<< endl;
		return;
	}

	string line;
	const string delimiter = ";";

	/* getting lines of parameter one after each other.
	 * First line corresponds to the keys of the parameters. */
	getline(file, line);
	tokenize(line, paramNames, delimiter);
    
	/* Retrieving the data in the other lines */
	size_t fuelNum = 0;
	vector<string> vals;
	double dval;
	while (getline(file, line)) {
		/* Affecting the map */
		table.push_back(map<string, double>());
		/* cutting the line according to delimiter */
		tokenize(line, vals, delimiter);
		if (vals.size() != paramNames.size()) {
			cout << "Number of parameters for fuel " << vals[0]
					<< " isn't in accordance with the number of parameters in the table"
					<< " (" << vals.size() << " vs " << paramNames.size() << ")"
					<< endl;
			;
		} else {
			for (size_t pos = 0; pos < vals.size(); pos++) {
				istringstream iss(vals[pos]);
				if (iss >> dval) {
					table[fuelNum].insert(make_pair(paramNames[pos], dval));
					//cout << "at "<<fuelNum<< "  param "<<paramNames[pos] <<" = "<<dval<<endl;
				} else {
					cout << "could not cast " << vals[pos]
							<< " into a suitable value " << "for parameter "
							<< paramNames[pos] << " of fuel " << vals[0]
							<< endl;
				}
			}
			fuelNum++;
		}
		vals.clear();
	}
	file.close();

}

FluxLayer<double>* DataBroker::constructFluxLayerFromFile(NcFile* NcdataFile,
		string property) {


	/* Getting the information on the mesh */
	/*-------------------------------------*/
	// Getting information on the extension of the domain
	const char* sdom = "domain";
	NcVar* dom = NcdataFile->get_var(sdom);


	FFPoint SWCorner = getNetCDFSWCorner(dom);
	FFPoint spatialExtent = getNetCDFSpatialSpan(dom);
	double timeOrigin = getNetCDFTimeOrigin(dom);
	double Lt = getNetCDFTimeSpan(dom);
    double version = getNetCDFFileVersion(dom);
	/* Getting the desired flux layer */
	/*--------------------------------*/
	NcVar* flux = NcdataFile->get_var(property.c_str());
	if (flux->num_dims() > 4) {
		cout << "Variable " << property << " doesn't have the right dimension ("
				<< flux->num_dims() << " instead of 4 maximum)" << endl;
	}

	/* Sending the information on the models to the domain */
	/*-----------------------------------------------------*/
	const char* mindex = "indices";
	NcAtt* indices = flux->get_att(mindex);
	NcAtt* tmpName;
	size_t numModels = (size_t) indices->num_vals();
	ostringstream mname;
	string modelName;
	for (size_t i = 0; i < numModels; i++) {
		mname.str("");
		int ind = indices->as_int(i);
		mname << "model" << ind << "name";
		tmpName = flux->get_att(mname.str().c_str());
		modelName = tmpName->as_string(0);
		delete tmpName;
		domain->fluxModelInstanciation(ind, modelName);
	}
	delete indices;

	size_t nx = 0;
	size_t ny = 0;
	size_t nz = 0;
	size_t nt = 0;
	if (flux->num_dims() == 2) {
		ny = (size_t) flux->get_dim(0)->size();
		nx = (size_t) flux->get_dim(1)->size();
	}
	if (flux->num_dims() == 3) {
		nz = (size_t) flux->get_dim(0)->size();
		ny = (size_t) flux->get_dim(1)->size();
		nx = (size_t) flux->get_dim(2)->size();

	}
	if (flux->num_dims() == 4) {
		nt = (size_t) flux->get_dim(0)->size();
		nz = (size_t) flux->get_dim(1)->size();
		ny = (size_t) flux->get_dim(2)->size();
		nx = (size_t) flux->get_dim(3)->size();
	}



	/* Getting the data */
	/*------------------*/
	int* data = readAndTransposeIntFortranProjectedField(flux, nt, nz, ny, nx,
			version>0,-1);

	if (isRelevantData(SWCorner, spatialExtent)) {

		/* Instanciating the data layer */
		/*------------------------------*/
		FluxLayer<double>* newlayer = new FluxLayer<double>(property,
				atmoSWCorner, atmoNECorner, atmosphericNx, atmosphericNy,
				domain->getCells(), data, SWCorner, timeOrigin, spatialExtent,
				Lt, nx, ny, nz, nt);
		delete[] data;
		return newlayer;

	} else {

		cout << "WARNING, spatial domain of validity for variable " << property
				<< ": " << SWCorner.print() << " -> ("
				<< SWCorner.getX() + spatialExtent.getX() << ","
				<< SWCorner.getY() + spatialExtent.getY() << ","
				<< SWCorner.getZ() + spatialExtent.getZ()
				<< ") has no intersection with the fire domain: ("
				<< domain->SWCornerX() << "," << domain->SWCornerY() << ") -> ("
				<< domain->NECornerX() << "," << domain->NECornerY() << ")"
				<< endl << "...resizing to the fire domain..." << endl;

		FFPoint SWC = FFPoint(domain->SWCornerX(), domain->SWCornerY(), 0.);
		FFPoint ext = FFPoint(domain->NECornerX() - domain->SWCornerX(),
				domain->NECornerY() - domain->SWCornerY(), 10000.);
		double t0 = 0;
		double Dt = 10000000.;
		FluxLayer<double>* newlayer = new FluxLayer<double>(property,
				atmoSWCorner, atmoNECorner, atmosphericNx, atmosphericNy,
				domain->getCells(), data, SWC, t0, ext, Dt, nx, ny, nz, nt);
		delete[] data;
		return newlayer;

	}

}

int* DataBroker::readAndTransposeIntFortranProjectedField(NcVar* val,
		const size_t &nt, const size_t &nz, const size_t &ny, const size_t &nx,
		bool transpose = true, int selectedT =-1) {

	int* tmp = new int[nt * nz * ny * nx];

	if (!val->get(tmp, nt, nz, ny, nx)) {
		cout << "error in getting the int data NC field" << endl;
		delete tmp;
		return NULL;
	}
	if (!transpose)
		return tmp;

	int* data = new int[nt * nz * ny * nx];

	size_t size = nt * nz * ny * nx;
	size_t indC, indF;
	size_t ii, jj, kk, rest;

	for (indF = 0; indF < size; indF++) {
		/* first compute the indices
		 in the Fortran representation of the array*/
		kk = indF / (nx * ny);
		rest = indF - kk * nx * ny;
		jj = rest / nx;
		ii = rest % nx;
		/* then compute the corresponding index
		 in C representation */
		indC = ii * ny * nz + jj * nz + kk;
		data[indC] = tmp[indF];
	}
	free(tmp);
	return data;
}
double* DataBroker::readAndTransposeFortranProjectedField(NcVar* val,
		const size_t &nt, const size_t &nz, const size_t &ny, const size_t &nx,
		bool transpose = true,  int selectedT = -1) {
	size_t nnt = nt;
	 if (selectedT>-1){
		 	 nnt = 1;
	    	val->set_cur(selectedT,0,0,0);
	    }


	double* tmp = new double[nnt * nz * ny * nx];


	if (!transpose){

		cout<< "not transposing field "<<endl;
		if (!val->get(tmp, nnt, nz, ny, nx)) {
				cout << "error in getting the NC field" << endl;
				delete tmp;
				return NULL;
			}
		return tmp;

	}

	if (!val->get(tmp, nnt, nz, ny, nx)) {
		cout << "error in getting the double data NC field" << endl;
		delete tmp;
		return NULL;
	}
	double* data = new double[nnt * nz * ny * nx];
	size_t size = nnt * nz * ny * nx;
	size_t indC, indF;
	size_t ii, jj, kk, rest;

	for (indF = 0; indF < size; indF++) {
		/* first compute the indices
		 in the Fortran representation of the array*/
		kk = indF / (nx * ny);
		rest = indF - kk * nx * ny;
		jj = rest / nx;
		ii = rest % nx;
		/* then compute the corresponding index
		 in C representation */
		indC = ii * ny * nz + jj * nz + kk;
		data[indC] = tmp[indF];
	}

	free(tmp);
	return data;
}

PropagativeLayer<double>* DataBroker::constructPropagativeLayerFromFile(
		NcFile* NcdataFile, string property) {

   /* Getting the information on the mesh */
	/*-------------------------------------*/
	// Getting information on the extension of the domain
	const char* sdom = "domain";
	NcVar* dom = NcdataFile->get_var(sdom);
	FFPoint SWCorner = getNetCDFSWCorner(dom);
	FFPoint spatialExtent = getNetCDFSpatialSpan(dom);
	double timeOrigin = getNetCDFTimeOrigin(dom);
	double Lt = getNetCDFTimeSpan(dom);
    double version = getNetCDFFileVersion(dom);
	/* Getting the desired flux layer */
	/*--------------------------------*/
	NcVar* prop = NcdataFile->get_var(property.c_str());
	if (prop->num_dims() > 4) {
		cout << "Variable " << property << " doesn't have the right dimension ("
				<< prop->num_dims() << " instead of 4 maximum)" << endl;
	}

	/* Sending the information on the models to the domain */
	/*-----------------------------------------------------*/
	const char* mindex = "indices";
	NcAtt* indices = prop->get_att(mindex);
	NcAtt* tmpName;
	size_t numModels = (size_t) indices->num_vals();
	ostringstream mname;
	string modelName;
	for (size_t i = 0; i < numModels; i++) {
		mname.str("");
		int ind = indices->as_int(i);
		mname << "model" << ind << "name";
		tmpName = prop->get_att(mname.str().c_str());
		modelName = tmpName->as_string(0);
		delete tmpName;
		domain->propModelInstanciation(ind, modelName);
	}
	delete indices;

	size_t nx = 0;
	size_t ny = 0;
	size_t nz = 0;
	size_t nt = 0;
	if (prop->num_dims() == 2) {
		ny = (size_t) prop->get_dim(0)->size();
		nx = (size_t) prop->get_dim(1)->size();
	}
	if (prop->num_dims() == 3) {
		nz = (size_t) prop->get_dim(0)->size();
		ny = (size_t) prop->get_dim(1)->size();
		nx = (size_t) prop->get_dim(2)->size();

	}
	if (prop->num_dims() == 4) {
		nt = (size_t) prop->get_dim(0)->size();
		nz = (size_t) prop->get_dim(1)->size();
		ny = (size_t) prop->get_dim(2)->size();
		nx = (size_t) prop->get_dim(3)->size();
	}

	/* Getting the data */
	/*------------------*/
	int* data = readAndTransposeIntFortranProjectedField(prop, nt, nz, ny, nx, version >0);

	delete prop;
	if (isRelevantData(SWCorner, spatialExtent)) {

		/* Instanciating the data layer */
		/*------------------------------*/
		PropagativeLayer<double>* newlayer = new PropagativeLayer<double>(
				property, data, SWCorner, timeOrigin, spatialExtent, Lt, nx, ny,
				nz, nt);
		delete[] data;
		return newlayer;

	} else {

		cout << "WARNING, spatial domain of validity for variable " << property
				<< ": " << SWCorner.print() << " -> ("
				<< SWCorner.getX() + spatialExtent.getX() << ","
				<< SWCorner.getY() + spatialExtent.getY() << ","
				<< SWCorner.getZ() + spatialExtent.getZ()
				<< ") has no intersection with the fire domain: ("
				<< domain->SWCornerX() << "," << domain->SWCornerY() << ") -> ("
				<< domain->NECornerX() << "," << domain->NECornerY() << ")"
				<< endl << "...resizing to the fire domain..." << endl;

		FFPoint SWC = FFPoint(domain->SWCornerX(), domain->SWCornerY(), 0.);
		FFPoint ext = FFPoint(domain->NECornerX() - domain->SWCornerX(),
				domain->NECornerY() - domain->SWCornerY(), 10000.);
		double t0 = 0;
		double Dt = 10000000.;
		PropagativeLayer<double>* newlayer = new PropagativeLayer<double>(
				property, data, SWC, t0, ext, Dt, nx, ny, nz, nt);
		delete[] data;
		return newlayer;

	}

}

void DataBroker::tokenize(const string& str, vector<string>& tokens,
		const string& delimiter = " ") {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiter, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter, lastPos);
	}
}


void DataBroker::getPropagationData(PropagationModel* model, FireNode* fn) {
	size_t nfilled = 0;
	if (optimizedPropDataBroker[model->index]) {
		for (size_t prop = 0; prop < numPropDataGetters[model->index]; prop++) {
			nfilled += (propDataGetters[model->index][prop])(fn, model,
					nfilled);
		}
	} else {
		for (size_t prop = 0; prop < model->numProperties; prop++) {
			nfilled += getLayer(model->wantedProperties[prop])->getValuesAt(fn,
					model, nfilled);
		}
	}
}

void DataBroker::getFluxData(FluxModel* model, FFPoint& loc, const double& t) {
	size_t nfilled = 0;
	if (optimizedFluxDataBroker[model->index]) {
		for (size_t prop = 0; prop < numFluxDataGetters[model->index]; prop++) {
			nfilled += (fluxDataGetters[model->index][prop])(loc, t, model,
					nfilled);
		}
	} else {
		for (size_t prop = 0; prop < model->numProperties; prop++) {
			nfilled += getLayer(model->wantedProperties[prop])->getValuesAt(loc,
					t, model, nfilled);
		}
	}
}

double DataBroker::getProperty(FireNode* fn, string property) {
	return getLayer(property)->getValueAt(fn);
}

DataLayer<double>* DataBroker::getLayer(const string& property) {
	// Scanning the scalar layers

	ilayer = layersMap.find(property);
	if (ilayer != layersMap.end())
		return ilayer->second;
	return 0;
}

FluxLayer<double>* DataBroker::getFluxLayer(const string& property) {
	// Scanning the scalar layers

	flayer = fluxLayersMap.find(property);
	if (flayer != fluxLayersMap.end())
		return flayer->second;

	return 0;
}
void DataBroker::computeActiveSurfacesFlux(const double &t) {
	// Scanning the scalar layers
	map<string, FluxLayer<double>*>::iterator iter = fluxLayersMap.begin();
	int numFluxModelsMax = 10;
	for (; iter != fluxLayersMap.end(); ++iter) {

		FluxLayer<double>* flayer = iter->second;
		int modelCount[numFluxModelsMax];
		for (int i = 0; i < numFluxModelsMax; i++)
			modelCount[i] = 0;
		flayer->computeActiveMatrix(t, modelCount);

	}

}

/* *************************************** */
/* Property getters for propagation models */
/* *************************************** */

int DataBroker::getDummy(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = dummyLayer->getValueAt(fn);
	return 1;
}

int DataBroker::getFuelProperties(FireNode* fn, PropagationModel* model,
		int start) {
	int numberOfValuesFilled = fuelLayer->getValuesAt(fn, model, start);
	return numberOfValuesFilled;
}

int DataBroker::getMoisture(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = moistureLayer->getValueAt(fn);
	return 1;
}

int DataBroker::getAltitude(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = altitudeLayer->getValueAt(fn);
	return 1;
}

int DataBroker::getSlope(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = slopeLayer->getValueAt(fn);
	return 1;
}

int DataBroker::getWindU(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = windULayer->getValueAt(fn);
	return 1;
}

int DataBroker::getWindV(FireNode* fn, PropagationModel* model, int keynum) {
	(model->properties)[keynum] = windVLayer->getValueAt(fn);
	return 1;
}

int DataBroker::getNormalWind(FireNode* fn, PropagationModel* model,
		int keynum) {
	double u = windULayer->getValueAt(fn);
	double v = windVLayer->getValueAt(fn);
	FFVector wind = FFVector(u, v);
	(model->properties)[keynum] = wind.scalarProduct(fn->getNormal());
	return 1;
}

int DataBroker::getFrontDepth(FireNode* fn, PropagationModel* model,
		int keynum) {
	(model->properties)[keynum] = fn->getFrontDepth();
	return 1;
}

int DataBroker::getFrontCurvature(FireNode* fn, PropagationModel* model,
		int keynum) {
	(model->properties)[keynum] = fn->getCurvature();
	return 1;
}
int DataBroker::getFrontFastestInSection(FireNode* fn, PropagationModel* model,
		int keynum) {
	(model->properties)[keynum] = fn->getLowestNearby(frontScanDistance);
	return 1;
}

/* ******************************** */
/* Property getters for flux models */
/* ******************************** */

int DataBroker::getFuelProperties(FFPoint loc, const double&t, FluxModel* model,
		int start) {
	int numberOfValuesFilled = fuelLayer->getValuesAt(loc, t, model, start);
	return numberOfValuesFilled;
}

int DataBroker::getMoisture(FFPoint loc, const double& t, FluxModel* model,
		int keynum) {
	(model->properties)[keynum] = moistureLayer->getValueAt(loc, t);
	return 1;
}

int DataBroker::getAltitude(FFPoint loc, const double& t, FluxModel* model,
		int keynum) {
	(model->properties)[keynum] = altitudeLayer->getValueAt(loc, t);
	return 1;
}

int DataBroker::getWindU(FFPoint loc, const double& t, FluxModel* model,
		int keynum) {
	(model->properties)[keynum] = windULayer->getValueAt(loc, t);
	return 1;
}

int DataBroker::getWindV(FFPoint loc, const double& t, FluxModel* model,
		int keynum) {
	(model->properties)[keynum] = windVLayer->getValueAt(loc, t);
	return 1;
}

void DataBroker::getMatrix(string matrixName, const FFPoint& SWCorner,
		const FFPoint& NECorner, const double& t, FFArray<double>** matrix) {
	ilayer = layersMap.find(matrixName);
	if (ilayer != layersMap.end()) {
		try {
			ilayer->second->getMatrix(matrix, t);
		} catch (...) {
			cout << "Encountered an error while retrieving matrix "
					<< matrixName << " with the data broker" << endl;
		}
	} else {

		 cout<<"Encountered an error while retrieving matrix "
		 <<matrixName<<" with the data broker: unknown matrix "<< matrixName<<endl;

	}
}

string DataBroker::printLayers() {
	ostringstream oss;
	oss << "Layers referenced in the data broker are:" << endl;
	for (ilayer = layersMap.begin(); ilayer != layersMap.end(); ++ilayer) {
		if (ilayer->second == 0) {
			oss << '\t' << "problem with a de-referenced layer" << endl;
			return oss.str();
		}
		oss << '\t' << "layer for property " << ilayer->second->getKey()
				<< " is at " << ilayer->second << endl;
	}
	return oss.str();
}

}
