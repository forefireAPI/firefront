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
#include "Command.h"
using namespace std;

namespace libforefire {

size_t Command::currentLevel = 0;

bool Command::init = true;
bool Command::currentFrontCompleted = false;

double Command::startTime = 0;
double Command::endTime = 0;

bool Command::firstCommand = true;
size_t Command::refTabs = 0;

FFPoint* Command::lastReadLoc = 0;
FireNode* Command::previousNode = 0;
FireNode* Command::leftLinkNode = 0;
FireNode* Command::rightLinkNode = 0;

double Command::bmapOutputUpdate = 0;
int Command::numBmapOutputs = 0;
double Command::refTime = 0;
int Command::numAtmoIterations = 0;

const string Command::stringError = "1234567890";
const FFPoint Command::pointError = FFPoint(1234567890., 1234567890.);
const FFVector Command::vectorError = FFVector(1234567890., 1234567890.);

vector<string> Command::outputDirs;


Command::Session Command::currentSession =
{
	/*struct Session{
		SimulationParameters* params;
		FireDomain* fd;
		FireDomain* fdp;
		FireFront* ff;
		StringRepresentation* outStrRep;
		StringRepresentation* outStrRepp;
		TimeTable* tt;
		Simulator* sim;
		ostream* outStream;
		int debugMode;
	};*/

		SimulationParameters::GetInstance(),
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		&cout,
		0,
};
 

// Construction of the dictionary
Command::CmdDictEntry Command::cmdDict[] =	{
		CmdDictEntry("FireDomain[...]","creates a fire domain")
		,CmdDictEntry("FireFront[...]","creates a fire front with the following fire nodes")
		,CmdDictEntry("FireNode[...]","adds a fire node to the given location")
		,CmdDictEntry("step[...]","advances the simulation for a user-defined amount of time")
		,CmdDictEntry("goTo[...]","advances the simulation till a user-defined time")
		,CmdDictEntry("setParameter[...]","sets a given parameter")
		,CmdDictEntry("setParameters[...]","sets a given list of parameters")
		,CmdDictEntry("getParameter[...]","return the parameter from a key")
		,CmdDictEntry("trigger[...]","set a value at a certain date to change model ")
		,CmdDictEntry("include[...]","includes commands from a file")
		,CmdDictEntry("print[...]","prints the state of the simulation")
		,CmdDictEntry("save[]","saves the simulation in hdf format")
		,CmdDictEntry("help","displays messages about the usage of commands")
		,CmdDictEntry("man[command]","displays the man page of the desired 'command'")
		,CmdDictEntry("loadData[...]","load a NC data file")
		,CmdDictEntry("systemExec[...]","runs a system command")
		,CmdDictEntry("clear[]","clears the simulation")
		,CmdDictEntry("quit","terminates the simulation")
};

const Command::commandMap Command::translator = Command::makeCmds();
const Command::commandMan Command::manpages = Command::makeMan();

 
// Defaults constructor and destructor for the 'Command' abstract class
Command::Command() {

}

Command::~Command(){
}

void Command::increaseLevel(){
	currentLevel++;
}


void Command::decreaseLevel(){
	currentLevel--;
}

int Command::createDomain(const string& arg, size_t& numTabs){

	if ( currentSession.fd != 0 ){
	 
	}

	size_t n = argCount(arg);
	if ( n >= 3 ) {
		FFPoint SW = getPoint("sw", arg);
		FFPoint NE = getPoint("ne", arg);
		double t = getFloat("t",arg);

		setStartTime(t);
		setReferenceTime(t);
		/* creating the domain */
		if ( currentSession.fd != 0 ){
 
			if(currentSession.fd->getDomainID() == 1){
				
				currentSession.params->setParameter("runmode", "masterMNH");
				currentSession.fdp = new FireDomain(t, SW, NE);
				currentSession.fdp->setTimeTable(currentSession.tt);
				//currentSession.ff = currentSession.fdp->getDomainFront();
				
				ostringstream ffOutputsFNAME;
				ffOutputsFNAME<<currentSession.params->getParameter("caseDirectory")<<'/'
						<<currentSession.params->getParameter("fireOutputDirectory")<<'/'
						<<currentSession.params->getParameter("outputFiles")
						<<"."<<currentSession.fdp->getDomainID();
				currentSession.params->setParameter("PffOutputsPattern", ffOutputsFNAME.str());
				currentSession.outStrRepp = new StringRepresentation(currentSession.fdp);
				currentSession.outStrRepp->setOutPattern(ffOutputsFNAME.str());
				if ( currentSession.params->getInt("outputsUpdate") != 0 ){
					currentSession.tt->insert(new FFEvent(currentSession.outStrRepp));
					
				}
				
				
			}
		}else{
			currentSession.fd = new FireDomain(t, SW, NE);
			
		
			/* setting up the pointer to the domain */
		

			/* creating the related timetable and simulator */
			currentSession.tt = new TimeTable(new FFEvent(getDomain()));
			currentSession.sim = new Simulator(currentSession.tt
					, currentSession.fd->outputs);
			getDomain()->setTimeTable(currentSession.tt);
			/* linking to the domain front */
			currentSession.ff = getDomain()->getDomainFront();
			/* managing the outputs */
			ostringstream ffOutputsPattern;
			ffOutputsPattern<<currentSession.params->getParameter("caseDirectory")<<'/'
					<<currentSession.params->getParameter("fireOutputDirectory")<<'/'
					<<currentSession.params->getParameter("outputFiles")
					<<"."<<getDomain()->getDomainID();
			currentSession.params->setParameter("ffOutputsPattern", ffOutputsPattern.str());
			currentSession.outStrRep = new StringRepresentation(currentSession.fd);
			if ( currentSession.params->getInt("outputsUpdate") != 0 ){
				currentSession.tt->insert(new FFEvent(currentSession.outStrRep));
			}
		/* increasing the level */
			increaseLevel();
		}
		/* local copy of the reference time */
		return normal;
	} else {
		throw MissingOption(3-n);
	}
}

int Command::startFire(const string& arg, size_t& numTabs){
	
	size_t n = argCount(arg);
    if (n < 2) {
        throw MissingTime();
    }

	if( getDomain()->getDomainID() > 0 ){
		 
			return normal;
	}
    
	FireDomain* refDomain = getDomain();
	//if (currentSession.fdp != 0) {
	//		cout<<"Starting fire reference from "<<refDomain->getDomainID()<<" to "<<currentSession.fdp->getDomainID()<<endl;
     //       refDomain = currentSession.fdp;
      //  } 

		
	  
    SimulationParameters *simParam = SimulationParameters::GetInstance();
    
    FFPoint pos = getPoint("loc", arg);
    if ((refDomain->striclyWithinDomain(pos)))
    {
        if (currentSession.ff != 0) {
            completeFront(currentSession.ff);
        }

		
        double t = getFloat("t", arg);
        
        if (t == FLOATERROR) 
        {
            string date = getString("date", arg);
            
            if (date == stringError) {
                return normal;
            }
            
            int year, yday;
            double secs;
            
             simParam->ISODateDecomposition(date, secs, year, yday) ;

            t = simParam->SecsBetween(simParam->getDouble("refTime"), simParam->getInt("refYear"), simParam->getInt("refDay"), secs, year, yday);
            
            cout << endl << endl << "REF DATE : " << simParam->FormatISODate(simParam->getDouble("refTime"), simParam->getInt("refYear"), simParam->getInt("refDay"))
                    << endl << "DATE OF FIRE : " << date << endl << "SECS BETWEEN BOTH : " << t << endl << endl;
        }
        
        double perimRes = refDomain->getPerimeterResolution() * 2;
        int fdom = getInt("domain", arg);
        if ( fdom == INTERROR ) fdom = 0;

        double fdepth = currentSession.params->getDouble("initialFrontDepth");
        double kappa = 0.;
        string state = "init";
        FireFront* contfront = currentSession.ff->getContFront();
        currentSession.ff = refDomain->addFireFront(t, contfront);

        FFVector vel1 = FFVector(0, 1);
        FFVector vel2 = FFVector(1, -1);
        FFVector vel3 = FFVector(-1, -1);
        FFVector diffP1 = perimRes*vel1;
        FFVector diffP2 = perimRes*vel2;
        FFVector diffP3 = perimRes*vel3;

        FFPoint pos1 = pos+diffP1.toPoint();
        FFPoint pos2 = pos+diffP2.toPoint();
        FFPoint pos3 = pos+diffP3.toPoint();
        
        vel1 *= 0.1;
        vel2 *= 0.1;
        vel3 *= 0.1; 

        FireNode* lastnode = refDomain->addFireNode(pos1, vel1, t, fdepth, kappa, currentSession.ff, 0);
        lastnode = refDomain->addFireNode(pos2, vel2, t, fdepth, kappa, currentSession.ff, lastnode);
        refDomain->addFireNode(pos3, vel3, t, fdepth, kappa, currentSession.ff, lastnode); 

		completeFront(currentSession.ff); 
    }
    
    //simParam->setParameter("ISOdate", simParam->FormatISODate(simParam->getInt("refTime") + refDomain->getSimulationTime(), refDomain->getReferenceYear(), refDomain->getReferenceDay()));
    
	return normal;
}
int Command::createFireFront(const string& arg, size_t& numTabs){

	size_t n = argCount(arg);
	if ( n >= 1 ) {
		double t = getFloat("t",arg);
		if ( numTabs == currentLevel ) {

			/* first completing the previous front if needed */
			if ( currentSession.ff != 0 ){
				completeFront(currentSession.ff);
				currentFrontCompleted = false;
			}
			/* creating the new front */
			FireFront* contfront = currentSession.ff->getContFront();
			currentSession.ff = getDomain()->addFireFront(t, contfront);

		} else if (  numTabs > currentLevel ) {

			/* first completing the previous front */
			completeFront(currentSession.ff);
			currentFrontCompleted = false;
			/* creation of an inner front to the current one */
			FireFront* contfront = currentSession.ff;
			currentSession.ff = getDomain()->addFireFront(t, contfront);
			currentLevel = numTabs;

		} else if ( numTabs < currentLevel ){

			/* first completing the previous front */
			completeFront(currentSession.ff);
			currentFrontCompleted = false;
			/* creation of a new front at a higher level then the current one */
			FireFront* contfront = currentSession.ff;
			for ( size_t k=0; k < currentLevel-numTabs; k++ ) {
				contfront = contfront->getContFront();
			}
			currentSession.ff = getDomain()->addFireFront(t, contfront);
			currentLevel = numTabs;

		}
	} else {
		throw MissingTime();
	}

	return normal;
}

int Command::addFireNode(const string& arg, size_t& numTabs){

	if ( lastReadLoc == 0 ) lastReadLoc = new FFPoint(-numeric_limits<double>::infinity(),-numeric_limits<double>::infinity());

	size_t n = argCount(arg);
	double perimRes = getDomain()->getPerimeterResolution();
	if ( n >= 3 ) {
		FFPoint pos = getPoint("loc", arg);
		FFVector vel = getVector("vel", arg);
		double t = getFloat("t",arg);
		int fdom = getInt("domain", arg);
		if ( fdom == INTERROR ) fdom = 0;
		int id = getInt("id", arg);
		if ( id == INTERROR ) id = 0;
		double fdepth = getFloat("fdepth", arg);
		if ( fdepth == FLOATERROR )
			fdepth = currentSession.params->getDouble("initialFrontDepth");
		double kappa = getFloat("kappa", arg);
		if ( kappa == FLOATERROR ) kappa = 0.;
		string state = getString("state", arg);
		if ( state == stringError or state == "moving" ) state = "init";
		if ( numTabs != currentLevel + 1 ) {
			cout<<getDomain()->getDomainID()<<": WARNING : asked for a FireNode "
					<<" with wrong indentation, treating it the current fire front"<<endl;
		}

		if ( state == "link" ){

			/* Creating a link node */
			previousNode = getDomain()->addFireNode(pos, vel, t
					, fdepth, kappa, currentSession.ff, previousNode
					, fdom, id , FireNode::link);
			if ( getDomain()->commandOutputs ) cout<<getDomain()->getDomainID()<<": INIT -> added "<<previousNode->toString()<<endl;

		} else if ( state == "final" ){

				/* Creating a link node */
				previousNode = getDomain()->addFireNode(pos, vel, t
						, fdepth, kappa, currentSession.ff, previousNode
						, fdom, id , FireNode::final);
				if ( getDomain()->commandOutputs ) cout<<getDomain()->getDomainID()<<": INIT -> added "<<previousNode->toString()<<endl;

		} else if ( getDomain()->striclyWithinDomain(pos)
				and getDomain()->striclyWithinDomain(*lastReadLoc) ){

			/* Both nodes are within the domain and represent real firenodes */
			/* checking the distance between the two nodes */
			double distanceBetweenNodes = lastReadLoc->distance2D(pos);
			int interNodes = (int) floor(distanceBetweenNodes/(2.*perimRes));
			FFPoint posinc = (1./(interNodes+1))*(pos-previousNode->getLoc());
			FFVector velinc = (1./(interNodes+1))*(vel-previousNode->getVel());
			double timeinc = (1./(interNodes+1))*(t-previousNode->getTime());
			FFPoint ipos;
			FFVector ivel;
			double itime;
			for ( int k = 0; k < interNodes; k++ ){
				ipos = previousNode->getLoc() + posinc;
				ivel = previousNode->getVel() + velinc;
				itime = previousNode->getTime() + timeinc;
				previousNode = getDomain()->addFireNode(ipos, ivel, itime
						, fdepth, kappa, currentSession.ff, previousNode);
				if ( getDomain()->commandOutputs )
					cout<<getDomain()->getDomainID()<<": INIT -> added inter-node "
					<<previousNode->toString()<<endl;
			}
			// creating the firenode
			previousNode = getDomain()->addFireNode(pos, vel, t
					, fdepth, kappa, currentSession.ff, previousNode, fdom, id);
			if ( getDomain()->commandOutputs ) cout<<getDomain()->getDomainID()
					<<": INIT -> added "<<previousNode->toString()<<endl;

		} else if ( !getDomain()->striclyWithinDomain(pos)
				and getDomain()->striclyWithinDomain(*lastReadLoc) ){

			FFPoint linkPoint = getDomain()->findIntersectionWithFrontiers(
					pos, *lastReadLoc);
			rightLinkNode = getDomain()->addLinkNode(linkPoint);
			/* checking the distance between the two nodes */
			double distanceBetweenNodes = linkPoint.distance2D(previousNode->getLoc());
			if ( distanceBetweenNodes > 2.*perimRes ){
				int interNodes = (int) floor(distanceBetweenNodes/(2.*perimRes));
				FFPoint posinc = (1./(interNodes+1))*(linkPoint-previousNode->getLoc());
				FFPoint ipos;
				for ( int k = 0; k < interNodes; k++ ){
					ipos = previousNode->getLoc() + posinc;
					previousNode = getDomain()->addFireNode(ipos, vel, t
							, fdepth, kappa, currentSession.ff, previousNode);
					if ( getDomain()->commandOutputs )
						cout<<getDomain()->getDomainID()<<": INIT -> added inter-node "
						<<previousNode->toString()<<endl;
				}
			}
			previousNode->insertAfter(rightLinkNode);
			if ( leftLinkNode != 0 ) getDomain()->relateLinkNodes(rightLinkNode, leftLinkNode);

		} else if ( getDomain()->striclyWithinDomain(pos)
				and !getDomain()->striclyWithinDomain(*lastReadLoc) ){

			if ( lastReadLoc->getX() == -numeric_limits<double>::infinity() ){

				/* First node to be created */
				previousNode = getDomain()->addFireNode(pos, vel, t, fdepth, kappa
						, currentSession.ff, previousNode, fdom, id);
				if ( getDomain()->commandOutputs )
					cout<<getDomain()->getDomainID()<<": INIT -> added "<<previousNode->toString()<<endl;

			} else {

				/* Creating a link node */
				FFPoint linkPoint =
						getDomain()->findIntersectionWithFrontiers(
								pos, *lastReadLoc);
				leftLinkNode = getDomain()->addLinkNode(linkPoint);
				if ( rightLinkNode != 0 ) {
					currentSession.ff->addFireNode(leftLinkNode, rightLinkNode);
				} else {
					leftLinkNode->setFront(currentSession.ff);
					currentSession.ff->addFireNode(leftLinkNode);
				}
				previousNode = leftLinkNode;
				/* checking the distance between the two nodes */
				double distanceBetweenNodes = leftLinkNode->getLoc().distance2D(pos);
				if ( distanceBetweenNodes > 2.*perimRes ){
					int interNodes = (int) floor(distanceBetweenNodes/(2.*perimRes));
					FFPoint posinc = (1./(interNodes+1))*(pos-leftLinkNode->getLoc());
					FFPoint ipos;
					for ( int k = 0; k < interNodes; k++ ){
						ipos = previousNode->getLoc() + posinc;
						previousNode = getDomain()->addFireNode(ipos, vel, t
								, fdepth, kappa, currentSession.ff, previousNode);
						if ( getDomain()->commandOutputs )
							cout<<getDomain()->getDomainID()<<": INIT -> added inter-node "
							<<previousNode->toString()<<endl;
					}
				}
				// creating the firenode
				previousNode = getDomain()->addFireNode(pos, vel, t
						, fdepth, kappa, currentSession.ff, previousNode, fdom, id);
				if ( getDomain()->commandOutputs ) cout<<getDomain()->getDomainID()
						<<": INIT -> added "<<previousNode->toString()<<endl;
			}

		}

		lastReadLoc->setX(pos.getX());
		lastReadLoc->setY(pos.getY());

		return normal;
	} else {
		throw MissingOption(3-n);
	}
}

void Command::completeFront(FireFront* ff){
	if ( ff->getHead() == 0 ) return;
	FireNode* firstNode = ff->getHead();
	FireNode* lastNode = firstNode->getPrev();
	if ( firstNode->getState() == FireNode::init and lastNode->getState() == FireNode::init ){
		double perimRes = getDomain()->getPerimeterResolution();
		double distanceBetweenNodes = lastNode->distance2D(firstNode);
		double fdepth = 0.5*(firstNode->getFrontDepth()+lastNode->getFrontDepth());
		double kappa = 0.5*(firstNode->getCurvature()+lastNode->getCurvature());
		int numNewNodes = (int) floor(distanceBetweenNodes/(2.*perimRes));
		if ( numNewNodes > 0 ){
			FFPoint posinc = (1./(numNewNodes+1))
					*(firstNode->getLoc()-lastNode->getLoc());
			FFVector velinc = (1./(numNewNodes+1))
					*(firstNode->getVel()-lastNode->getVel());
			double timeinc = (1./(numNewNodes+1))
					*(firstNode->getTime()-lastNode->getTime());
			FireNode* prevNode = lastNode;
			FFPoint pos;
			FFVector vel;
			double t;
			for ( int k = 0; k < numNewNodes; k++ ){
				pos = lastNode->getLoc() + (k+1)*posinc;
				vel = lastNode->getVel() + (k+1)*velinc;
				t = lastNode->getTime() + (k+1)*timeinc;
				prevNode = getDomain()->addFireNode(pos, vel, t
						, fdepth, kappa, currentSession.ff, prevNode);
				prevNode->computeNormal();
			}
		}
	}
	if ( getDomain()->commandOutputs ){
		cout<<getDomain()->getDomainID()<<": "
			<< "****************************************" << endl;
		cout<<getDomain()->getDomainID()<<": "
			<< " BEFORE THE BEGINNING OF THE SIMULATION: " << endl
			<< ff->print(1);
		cout<<getDomain()->getDomainID()<<": "
			<< "****************************************" << endl;
	}

	if ( !currentSession.params->isValued("BMapFiles") ){
		/* testing the domain for burning matrix */
		double iniFrontDepth = currentSession.params->getDouble("initialFrontDepth");
		double iniBurningTime = currentSession.params->getDouble("initialBurningDuration");

		if ( !currentSession.params->isValued("noInitialScan"))
			{

				getDomain()->frontInitialBurningScan(ff->getTime(), ff, iniFrontDepth, iniBurningTime);
			}
	}

	currentFrontCompleted = true;
	delete lastReadLoc;
	lastReadLoc = 0;

}

int Command::stepSimulation(const string& arg, size_t& numTabs){
    if (getDomain() == 0) return normal;
    
	double dt = getFloat("dt",arg);
	endTime = startTime + dt;
	ostringstream etime;
	etime.precision(numeric_limits<double>::digits10);
	etime<<"t="<<endTime;
	goTo(etime.str().c_str(), numTabs);
    
    SimulationParameters *simParam = SimulationParameters::GetInstance();
    simParam->setParameter("ISOdate", simParam->FormatISODate(simParam->getInt("refTime") + getDomain()->getSimulationTime(), getDomain()->getReferenceYear(), getDomain()->getReferenceDay()));
    
	return normal;
}

int Command::goTo(const string& arg, size_t& numTabs){

	/* Advancing simulation to the prescribed time */
	endTime = getFloat("t",arg);

	if ( init ){
		/* finishing the initialization process */
		bmapOutputUpdate = currentSession.params->getDouble("bmapOutputUpdate");
		if ( !currentFrontCompleted ) completeFront(currentSession.ff);
		numAtmoIterations = currentSession.params->getInt("numAtmoIterations") - 1;
		init = false;
	}

	if ( endTime > startTime ){
		getDomain()->setTime(startTime);
		if ( getDomain()->commandOutputs ){
			cout.precision(numeric_limits<double>::digits10);
			cout<<getDomain()->getDomainID()<<": "
				<< "***************************************************" << endl;
			cout<<getDomain()->getDomainID()<<": "
				<< "   ADVANCING FOREFIRE SIMULATION FROM T="
				<< startTime << " to " << endTime << endl;
			cout<<getDomain()->getDomainID()<<": "
				<< "***************************************************" << endl;
		}


		try {

			if ( getDomain()->parallel ){
			//	cout<<getDomain()->getDomainID()<<" at time "<<startTime <<"to "<<endTime <<" writing in halo"<<endl;
				/* managing incoming information from other procs */
				/* ********************************************** */
				if ( getDomain()->commandOutputs ) cout<<getDomain()->getDomainID() << ": managing halo " << endl;
				//getDomain()->manageHaloFirenodes(startTime);
			}

			/* ******************************************** */
			/* Advancing the simulation to the desired time */
			/* ******************************************** */
			//cout<<getDomain()->getDomainID()<<"  iteration : "<<FireDomain::atmoIterNumber<<" and "<<getDomain()->getNumIterationAtmoModel()<<endl;
			FireDomain::atmoIterNumber = FireDomain::atmoIterNumber+1;
			currentSession.params->setInt("atmoIterNumber",FireDomain::atmoIterNumber);
			currentSession.fd->loadCellsInBinary();
			getDomain()->loadWindDataInBinary(endTime);
			currentSession.sim->goTo(endTime);
			getDomain()->dumpCellsInBinary();
			//currentSession.fd->dumpWindDataInBinary();
			startTime = endTime;
			getDomain()->increaseNumIterationAtmoModel();

			if ( getDomain()->parallel ){
				/* sending information to the other procs */
				/* ************************************** */
				//getDomain()->createFirenodesMatrices();
			}

			/* outputs */
			/* ******* */
			if ( getDomain()->commandOutputs ){
				cout<<getDomain()->getDomainID()<<": "
						<< "End of the step in domain "
						<< getDomain()->getDomainID() << endl
						<< currentSession.outStrRep->dumpStringRepresentation();
			}
			if ( currentSession.params->getInt("debugFronts") != 0 )
				printSimulation(currentSession.params->getParameter("ffOutputsPattern"), numTabs);
			/* burning map outputs */
			if ( bmapOutputUpdate > 0 ){
				if ( (int) ((endTime-refTime)/bmapOutputUpdate) > numBmapOutputs ){
					getDomain()->saveSimulation();
					numBmapOutputs++;
				}
			}

			/* backing up the state for future steps */
			/* ************************************* */
			getDomain()->validateTopology("advance");
			getDomain()->backupState();

			/* saving simulation if needed */
			/* *************************** */
			if ( getDomain()->getNumIterationAtmoModel() == numAtmoIterations )
				getDomain()->saveSimulation();

		} catch ( TopologicalException& e ) {
			cout<<getDomain()->getDomainID()<<": "<<e.what()<<endl;
			getDomain()->restoreValidState();
			getDomain()->setSafeTopologyMode(true);
			cout<<getDomain()->getDomainID()<<": "
					<< "**** MAKING THIS STEP IN SAFE TOPOLOGY MODE ****" << endl;
			try {

				/* ******************************************** */
				/* Advancing the simulation to the desired time */
				/* ******************************************** */
				currentSession.sim->goTo(endTime);
				// Getting out of safe topology mode
				getDomain()->setSafeTopologyMode(false);
				startTime = endTime;
				getDomain()->increaseNumIterationAtmoModel();

				if ( getDomain()->parallel ){
					/* sending information to the other procs */
					/* ************************************** */
					getDomain()->createFirenodesMatrices();
				}

				if ( getDomain()->commandOutputs ){
					cout<<getDomain()->getDomainID()<<": "
							<< "End of the step in domain "
							<< getDomain()->getDomainID() << endl
							<< currentSession.outStrRep->dumpStringRepresentation();
				}

				/* outputs */
				/* ******* */
				if ( getDomain()->commandOutputs ){
					cout<<getDomain()->getDomainID()<<": "
							<< "End of the step in domain "
							<< getDomain()->getDomainID() << endl
							<< currentSession.outStrRep->dumpStringRepresentation();
				}
				if ( currentSession.params->getInt("debugFronts") != 0 )
					printSimulation(currentSession.params->getParameter("ffOutputsPattern"), numTabs);
				/* burning map outputs */
				if ( bmapOutputUpdate > 0 ){
					if ( (int) ((endTime-refTime)/bmapOutputUpdate) > numBmapOutputs ){
						getDomain()->saveSimulation();
						numBmapOutputs++;
					}
				}

				/* backing up the state for future steps */
				/* ************************************* */
				getDomain()->validateTopology("advance");
				getDomain()->backupState();

				/* saving simulation if needed */
				/* *************************** */
				if ( getDomain()->getNumIterationAtmoModel() == numAtmoIterations )
					getDomain()->saveSimulation();
			} catch (...) {
				if ( getDomain()->commandOutputs ){
					cout<<getDomain()->getDomainID()<<": "
						<< "**** ERROR IN SAFE TOPOLOGY MODE, QUITING ****" << endl;
				}
				// TODO supersafe mode ?
				quit(arg, numTabs);
			}
		}
	}
	return normal;
}

int Command::printSimulation(const string& arg, size_t& numTabs){
	if ( getDomain() == 0 ) return normal;
    
    SimulationParameters *simParam = SimulationParameters::GetInstance();
    string finalStr = "";
    
    if (arg.size() > 0)
    {
        vector<string> parts;
        tokenize(arg, parts, "*");

        int partToEval = (arg.at(0) == '*') ? 0 : 1;

        for (int i = 0; i < parts.size(); i++) {
            finalStr += (i % 2 == partToEval) ? simParam->getParameter(parts[i]) : parts[i];
        }
    }
    
    replace(finalStr.begin(), finalStr.end(), ':', '-');
    
	ofstream outputfile(finalStr.c_str());
	if ( outputfile ) {
        simParam->setInt("count", simParam->getInt("count") + 1);
		outputfile<<currentSession.outStrRep->dumpStringRepresentation();
	} else {
		*currentSession.outStream<<currentSession.outStrRep->dumpStringRepresentation();
	}
	return normal;
}
int Command::systemExec(const string& arg, size_t& numTabs){


    SimulationParameters *simParam = SimulationParameters::GetInstance();
    string finalStr = "";

    if (arg.size() > 0)
    {
        vector<string> parts;
        tokenize(arg, parts, "*");

        int partToEval = (arg.at(0) == '*') ? 0 : 1;

        for (int i = 0; i < parts.size(); i++) {
            finalStr += (i % 2 == partToEval) ? simParam->getParameter(parts[i]) : parts[i];
        }
    }

    replace(finalStr.begin(), finalStr.end(), ':', '-');


	return std::system(finalStr.c_str());;
}

int Command::saveSimulation(const string& arg, size_t& numTabs){
	if ( getDomain() == 0 ) return normal;
	getDomain()->saveSimulation();
	return normal;
}

int Command::setParameters(const string& arg, size_t& numTabs){
	// Getting all the arguments, using the 'tokenize' function
	// Returns strings of the form 'arg=...'
	vector<string> tmpArgs;
	string delimiter = ";";
	tokenize(arg, tmpArgs, delimiter);
	for ( size_t i = 0, size = tmpArgs.size(); i < size; ++i ) {
		setParameter(tmpArgs[i], numTabs);
	}
	return normal;
}

int Command::setParameter(const string& arg, size_t& numTabs){
	// Getting all the arguments, using the 'tokenize' function
	// Returns strings of the form 'arg=...'
	vector<string> tmpArgs;
	string delimiter = "=";
	tokenize(arg, tmpArgs, delimiter);
	if ( tmpArgs.size() > 2 ){
		cout<<"problem in the number of arguments when setting "<<arg<<endl;
		return error;
	}
	currentSession.params->setParameter(tmpArgs[0], tmpArgs[1], true);
	return normal;
}

int Command::getParameter(const string& arg, size_t& numTabs){

    if (arg.size() == 0)
	{
		cout << "You have to specify one argument" << endl;
		return error;
	}
    
    string param = currentSession.params->getParameter(arg);
    
    if (param == "1234567890")
    {
		cout << "Parameter doesn't exist : " << arg << endl;
		return error;
	}
    
	*currentSession.outStream << param << endl;
	return normal;
}

int Command::triggerValue(const string& arg, size_t& numTabs){

	vector<string> tmpArgs;
	string delimiter = ";";
	tokenize(arg, tmpArgs, delimiter);
	if ( tmpArgs.size() < 2 ){
		return error;
	}


	if(tmpArgs[0]== "wind"){
				FFVector a = getVector("vel", arg);
				FFPoint b = getPoint("loc", arg);

				if(getDomain() != 0){
					if(getDomain()->getDataBroker() != 0){
						if(getDomain()->getDataBroker()->PwindULayer != 0){
							getDomain()->getDataBroker()->PwindULayer->setProjectionDirVector(a,b);
							}
						if(getDomain()->getDataBroker()->PwindVLayer != 0){
										getDomain()->getDataBroker()->PwindVLayer->setProjectionDirVector(a,b);
							}
						return normal;
						}

				}
		}

	if(tmpArgs[0]== "fuel"){

			if(getDomain()->getDataBroker() != 0){


	   			vector<string> tmpVal;
					string delimiter = "=";
					double val = 0;
					tokenize(tmpArgs[1], tmpVal, delimiter);
					istringstream iss(tmpVal[1]);
					if ( iss >> val ){
						getDomain()->updateFuelTable(tmpVal[0],val);

						return normal;
					}
					return error;
			}

	}
	if(tmpArgs[0]== "fuelIndice"){

			if(getDomain()->getDataBroker() != 0){

				FFPoint loc = getPoint("loc", arg);
				int fvalue = getInt("fuelType",arg);
			//	double oldval = currentSession.fd->getDataBroker()->getLayer("fuel")->getValueAt(loc,0);

				//cout <<"at "<< loc.x<< "we have fuel "<< oldval<<endl;
			//	cout <<"setting to "<< fvalue <<endl;

				getDomain()->getDataBroker()->getLayer("fuel")->setValueAt(loc,0.0,fvalue);
			//	double newval = currentSession.fd->getDataBroker()->getLayer("fuel")->getValueAt(loc,0);
			//	cout <<"at "<< loc.x<< " we have now fuel "<< fvalue <<endl;
				return normal;
			}

	}
    return error;

}
int Command::include(const string& arg, size_t& numTabs){
	/* the commands are read from a defined file
	 *  testing for the presence of the file */
	ifstream instream(arg.c_str());
	if ( instream ) {
		// reading all the commands (one command per line) of the input file
		string line;
		//size_t numLine = 0;
		while ( getline( instream, line ) ) {
			//    while ( getline( cin, line ) ) {
			//numLine++;
			// checking for comments or newline
			if((line[0] == '#')||(line[0] == '*')||(line[0] == '\n'))
				continue;
			// treating the command of the current line
			ExecuteCommand(line);
		}
	} else {
		cout << "wrong input file, check your settings..." << endl;
	}
	return normal;
}

int Command::help(const string& arg, size_t& numTabs){
	for(int i = 0; i < numberCommands ; i++) {
		cout << cmdDict[i].inlineCmd << " : " << cmdDict[i].usage << endl;
	}
	return normal;
}

int Command::man(const string& cmd, size_t& numTabs){
	commandMan::const_iterator cmdman = manpages.find( cmd );
	if ( cmdman == manpages.end() ) {
		if(cmd.at(0) != '!')
			cout << "unknown command, try 'help'." << endl;
	} else {
		// calling the right function
		cout << cmdman->second;
	}
	return normal;
}

int Command::loadData(const string& arg, size_t& numTabs){
    
    if (arg.size() == 0)
	{
		cout << "You have to specify the path to the data" << endl;
		return error;
	}
    
    vector<string> args;
	tokenize(arg, args, ";");
    
    if (args.size() > 2)
    {
        cout << "Error: Need 1 or 2 arguments, not more" << endl;
		return error;
    }
    
    SimulationParameters *simParam = SimulationParameters::GetInstance();
    string path = simParam->GetPath(args[0]);

    
	if (std::ifstream(path.c_str()).fail())
	{
		cout << "File "<< path<<" doesn't exist or no longer available" << endl;
		return error;
	}
#ifdef NETCDF_LEGACY 
	NcFile* ncFile = new NcFile(path.c_str(), NcFile::ReadOnly);

	if (!ncFile->is_valid())
	{
		cout << "File is not a valid NetCDF file" << endl;
		return error;
	}
    
    simParam->setParameter("NetCDFfile", args[0]);

	for (int layer = 0; layer < ncFile->num_vars(); layer++)
	{
		string varName(ncFile->get_var(layer)->name());

		if (varName == "domain")
		{

			NcVar* ncparams = ncFile->get_var(varName.c_str());
			int numParams = (size_t) ncparams->num_atts();
			NcAtt* ncparam;

			for (int i = 0; i < numParams; i++)
			{
				ncparam = ncparams->get_att(i);
				simParam->setParameter(ncparam->name(), ncparam->as_string(0));
			}

			if (ncparam != 0)
				delete ncparam;
		}
	}

#else
  simParam->setParameter("NetCDFfile", args[0]);
  try
   {
	NcFile dataFile(path.c_str(), NcFile::read);
	 	if (!dataFile.isNull()) {
			
			NcVar domVar = dataFile.getVar("domain");
			if (!domVar.isNull()) {
				  map<string,NcVarAtt> attributeList = domVar.getAtts();
				  map<string,NcVarAtt>::iterator myIter; 
					for(myIter=attributeList.begin();myIter !=attributeList.end();++myIter)
					{
					NcVarAtt att = myIter->second;
					NcType attValType = att.getType();
					string attsVal;
					int attiVal;
					float attfVal;
					switch ((int)attValType.getTypeClass()) {
							case NC_CHAR:
								
								att.getValues(attsVal); 
								simParam->setParameter(myIter->first, attsVal);
								break;
							case NC_INT:
								
								att.getValues(&attiVal); 
								simParam->setParameter(myIter->first, std::to_string(attiVal));
								break;
							case NC_FLOAT:
								 
								att.getValues(&attfVal); 
								simParam->setParameter(myIter->first, std::to_string(attfVal));
								break;
							default:
								std::cout << myIter->first << " attribute of unhandled type " <<attValType.getName() << endl;
								break;
						}
					}
			}
		}
	}catch(NcException& e)
	{
		e.what();
		cout<<"cannot read landscape file"<<endl;
	
	}
#endif
    if (args.size() > 1){
        double secs;
        int year, yday;
        if (simParam->ISODateDecomposition(args[1], secs, year, yday))
        {
            simParam->setInt("refYear", year);
            simParam->setInt("refDay", yday);
            simParam->setInt("refTime", secs);
            simParam->setParameter("ISOdate", args[1]);
        }
    }

    if (args.size() == 2){
		string com = "FireDomain[sw=("+simParam->getParameter("SWx")+".,"+simParam->getParameter("SWy")+".,"+simParam->getParameter("SWz")+".);ne=(";
		{
			std::ostringstream oss;
			oss << (simParam->getDouble("Lx") + simParam->getDouble("SWx"));
			oss << ".,";
			oss << (simParam->getDouble("Ly") + simParam->getDouble("SWy"));
			oss << ".,";
			oss << (simParam->getDouble("Lz") + simParam->getDouble("SWz"));
			oss << ".);t="+simParam->getParameter("t0")+".]";
			com += oss.str();
		}

		ExecuteCommand(com);

    }


    //getDomain()->addLayer("data", "windU", "windU");
    //getDomain()->addLayer("data", "windV", "windV");
	return normal;
}

int Command::clear(const string& arg, size_t& numTabs){
	delete currentSession.fd;
	delete currentSession.outStrRep;
	delete currentSession.sim;
	delete currentSession.params;
	return normal;
}

int Command::quit(const string& arg, size_t& numTabs){
	cout<<"*** QUITTING FOREFIRE, GOODBYE ***"<<endl;
	delete currentSession.fd;
	delete currentSession.outStrRep;
	delete currentSession.sim;
	delete currentSession.params;
	exit(0);
	return normal;
}


void Command::setOstringstream(ostringstream* oss){
	currentSession.outStream = oss;
}

void Command::tokenize(const string& str, vector<string>& tokens,
		const string& delimiter = " ") {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiter, lastPos);

	while (string::npos != pos || string::npos != lastPos){
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter, lastPos);
	}
}

string Command::getString(string opt, string arg){
	vector<string> tmpArgs;
	string tmparg;
	// Getting all the arguments, using the 'tokenize' function
	// Returns strings of the form 'arg=...'
	string delimiter1 = ";";
	tokenize(arg, tmpArgs, delimiter1);
	// Checking every 'tmpArgs' to see if it contains 'opt'
	// on the right-hand-side of 'arg=...' (i.e. arg=opt)
	vector<string> inArg;
	for ( size_t i = 0, size = tmpArgs.size(); i < size; ++i ) {
		string delimiter2 = "=";
		tokenize(tmpArgs[i], inArg, delimiter2);
		if ( opt.compare(inArg[0]) == 0 ) {
			return inArg[1];
		} else {
			inArg.clear();
		}
	}
	return stringError;
}

int Command::getInt(string opt, string arg){
	/* The 'getFloat()' method simply casts the results
	 * of the search for the option 'opt' given by 'getString()' */
	string tmpstr = getString(opt,arg);
	if ( tmpstr == stringError ) return INTERROR;
	string s;
	string::reverse_iterator it = tmpstr.rbegin();
	if ( *it == 's' ){
		// removing it
		size_t len = tmpstr.length();
		s = tmpstr.substr(0,len-1);
	} else {
		s = tmpstr;
	}
	istringstream iss(s);
	int d;
	if ( iss >> d ){
		return d;
	} else {
		return INTERROR;
	}
}

double Command::getFloat(string opt, string arg){
	/* The 'getFloat()' method simply casts the results
	 * of the search for the option 'opt' given by 'getString()' */
	string tmpstr = getString(opt,arg);
	if ( tmpstr == stringError ) return FLOATERROR;
	string s;
	string::reverse_iterator it = tmpstr.rbegin();
	if ( *it == 's' ){
		// removing it
		size_t len = tmpstr.length();
		s = tmpstr.substr(0,len-1);
	} else {
		s = tmpstr;
	}
	istringstream iss(s);
	double d;
	if ( iss >> d ){
		return d;
	} else {
		return FLOATERROR;
	}
}

FFPoint Command::getPoint(string opt, string arg){
	/* The 'getPoint()' method casts the results of the
	 * search for the option 'opt' given by 'getString()';
	 * as this result is in the form "(...,...,...)"
	 * this string has to be tokenized */

	string tmpstr = getString(opt, arg);
	if ( tmpstr == stringError ) {
		return pointError;
	}
	// droping the parenthesis
	string point = tmpstr.substr(1,getString(opt,arg).size()-2);

	vector<string> p;
	// Getting all the arguments, using the 'tokenize' function
	string delimiter = ",";
	tokenize(point, p, delimiter);
	// Passing the 'strings' into 'istringstream' for cast into 'double'
	double px, py, pz;
	istringstream psx(p[0]);
	istringstream psy(p[1]);
	istringstream psz(p[2]);
	if ( (psx>>px) && (psy>>py) && (psz>>pz) ) {
		// Possible shift at the initialization
		FFPoint shiftPos = FFPoint(currentSession.params->getDouble("SHIFT_ALL_POINT_ABSCISSA_BY")
				, currentSession.params->getDouble("SHIFT_ALL_POINT_ORDINATES_BY"));
		FFPoint relPos = FFPoint(px,py,pz);
		FFPoint pos = shiftPos + relPos;
		return pos;
	} else {
		return pointError;
	}
}

FFVector Command::getVector(string opt, string arg){
	/* The 'getPoint()' method casts the results of the
	 * search for the option 'opt' given by 'getString()';
	 * as this result is in the form "(...,...,...)"
	 * this string has to be tokenized */

	string tmpstr = getString(opt,arg);
	if ( tmpstr == stringError ) return vectorError;
	// droping the parenthesis
	string vec = tmpstr.substr(1,getString(opt,arg).size()-2);

	// Getting all the arguments, using the 'tokenize' function
	vector<string> v;
	string delimiter = ",";
	tokenize(vec, v, delimiter);
	// Passing the 'strings' into 'istringstream' for cast into 'double'
	double vx, vy, vz;
	istringstream vsx(v[0]);
	istringstream vsy(v[1]);
	istringstream vsz(v[2]);
	if ( (vsx>>vx) && (vsy>>vy) && (vsz>>vz) ) {
		return FFVector(vx,vy,vz);
	} else {
		return vectorError;
	}
}

size_t Command::argCount(string arg){
	vector<string> tmpArgs;
	string tmparg;
	// Getting all the arguments, using the 'tokenize' function
	// Returns strings of the form 'arg=...'
	string delimiter = ";";
	tokenize(arg, tmpArgs, delimiter);
	return tmpArgs.size();
}

void Command::setStartTime(const double& t){
	startTime = t;
}

void Command::setReferenceTime(const double& t){
	refTime = t;
}

double Command::getTime(){
	return startTime;
}

size_t Command::tabsCount(string cmdLine){
	size_t k = 0;
	while ( cmdLine[k] == '\t' ) {
		k++;
	}
	if ( firstCommand ){
		refTabs = k;
		firstCommand = false;
	}
	return k-refTabs;
}

string Command::removeTabs(string arg){
	size_t k = 0;
	while ( arg[k] == '\t' ){
		arg.erase(k,1);
	}
	return arg;
}

void Command::ExecuteCommand(string& line){
	/* The method 'ExecuteCommands' uses the 'tokenize' class
	 * to execute the command given by the user
	 * in the 'line'. If first looks at the possible options
	 * (listed in the tables 'short_options' and 'long_options')
	 * and then call the corresponding function */

	// Getting all the arguments, using the 'tokenize' function
	vector<string> command;
	size_t numTabs;
	string scmd;
	string postcmd;
	const string delimiter = "[";
	tokenize(line, command, delimiter);
	if ( command.size() > 1 ) {
		// counting and removing the tabs for the level of the command
		numTabs = tabsCount(command[0]);
		scmd = removeTabs(command[0]);
		if (command[1].size() > 1){
			// dropping the last parenthesis ']'
			postcmd= command[1].substr(command[1].rfind("]")+1,-1) ;
			command[1] = command[1].substr(0,command[1].rfind("]"));
		} else {
			command[1] = "";
		}
	} else {
		scmd = removeTabs(line);
	}

	// calling the right method using the 'translator'
	if( ((scmd)[0] == '#')||((scmd)[0] == '*')||((scmd)[0] == '\n')||((scmd)[0] == '\r') ){
		// this line is commented, nothing to do
		} else {
				if( postcmd.size() > 1){
					string atCommand=line.substr(0,line.rfind("]")+1);
					string whenCommand = line.substr(line.rfind("]")+1,-1);
					double whenDouble = getFloat("@t", whenCommand) ;
						if ((whenCommand)[0] == '@' && whenDouble != FLOATERROR){
						    currentSession.tt->insert(new FFEvent(new EventCommand(atCommand, whenDouble)));
						}else{
							cout <<whenDouble<< "unknown post operator  >"<< whenCommand << "< , try 'help[]'." << endl;
						}
				}
				else{
					commandMap::const_iterator curcmd = translator.find( scmd );
					if ( curcmd == translator.end() ) {
						if(scmd.at(0) != '!')
							cout << "unknown command  >"<< scmd << "< , try 'help[]'." << endl;
					} else {
						try {
							// calling the right function
							(curcmd->second)(command[1], numTabs);
						} catch( BadOption& ) {
							cout<<getDomain()->getDomainID()<<": "
									<< "argument(s) '" << command[1] << "' is (are) not fit for command '"
									<< command[0] << "'" << endl;
							cout << "type 'man[" << command[0] << "]' for more information" << endl;
						} catch( MissingOption& mo) {
							cout<<getDomain()->getDomainID()<<": "
									<< "command '" << command[0] << "' misses "
									<< mo.num << " options." << endl;
							cout << "type 'man[" << command[0] << "]' for more information" << endl;
						} catch( MissingTime& ) {
							cout<<getDomain()->getDomainID()<<": "
									<< "you have to specify a time for that command (as in 't=0.')" << endl;
						} catch (...) {
							cout<<getDomain()->getDomainID()<<": "
									<< " PROBLEM: Command associated to "<<line<<" ended in error"<< endl;
						}
					}
				}
		}
	command.clear();
}

FireDomain* Command::getDomain(){
	if (currentSession.fdp != 0) {
		return currentSession.fdp;
	}
	return currentSession.fd;
}

string Command::dumpString(){
	return currentSession.outStrRep->outputstr.str();
	return 0;
}

}
