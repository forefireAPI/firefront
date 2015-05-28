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

#include "StringRepresentation.h"

#define JSON_MODE  1
#define FF_MODE  0

using namespace std;

namespace libforefire {

// static variables initialization
FireDomain* StringRepresentation::domain = 0;
size_t StringRepresentation::currentLevel = 0;
ostringstream StringRepresentation::outputstr("");
bool writedOnce = false;

StringRepresentation::StringRepresentation(FireDomain* fdom) : Visitor() {
    lastLevel = -1;
	domain = fdom;
	setTime(domain->getTime());
	updateStep = SimulationParameters::GetInstance()->getDouble("outputsUpdate");
	setUpdateTime(domain->getTime());
	dumpMode = FF_MODE;
}

StringRepresentation::~StringRepresentation() {
}

void StringRepresentation::input(){
}

void StringRepresentation::update(){
	setTime(getUpdateTime());
}

void StringRepresentation::timeAdvance(){
	if ( updateStep < EPSILONT ){
		setUpdateTime(numeric_limits<double>::infinity());
	} else {
		setUpdateTime(getTime()+updateStep);
	}
}

void StringRepresentation::output(){

	if (writedOnce && (domain->getNumFF()==0) ) return;

	writedOnce= true;

	ostringstream oss;
	oss<<SimulationParameters::GetInstance()->getParameter("ffOutputsPattern")
			<<"."<<getTime();
	ofstream outputfile(oss.str().c_str());
	if ( outputfile ) {
		outputfile<<dumpStringRepresentation();
	} else {
		cout<<"could not open file "<<oss.str()<<" for writing domain file "<<endl;
	}
}

void StringRepresentation::visit(FireDomain* fd) {
    if (dumpMode == JSON_MODE)
    {
        outputstr << '{' << endl << '\t' << "\"fronts\": [";
        lastLevel = 0;
    }
    
    if (dumpMode == FF_MODE)
        outputstr << fd->toString() << endl;
}

void StringRepresentation::visit(FireFront* ff) {
    if (dumpMode == FF_MODE)
    {
        for ( size_t k=0; k < currentLevel; k++ ) outputstr << '\t';
        outputstr << ff->toString() << endl;
        return;
    }
        
    if (dumpMode == JSON_MODE)
    {
        SimulationParameters *simParam = SimulationParameters::GetInstance();
        
        if (ff->getDomain()->getSimulationTime() >= ff->getTime())
        {
            if (lastLevel >= 2)
                outputstr << '"';
            if (lastLevel >= 1)
                outputstr << endl << '\t' << "},";
            
            double t = simParam->getInt("refTime") + ff->getDomain()->getSimulationTime();
            int d = simParam->getInt("refDay");
            int y = simParam->getInt("refYear");

            outputstr.precision(3);
            outputstr << endl << '\t' << '{';
            outputstr << endl << "\t\t" << "\"area\":\"";
            outputstr << fixed << (ff->getArea() / 10000) << "ha\",";
            outputstr << endl << "\t\t" << "\"date\":\"";
            outputstr << SimulationParameters::FormatISODate(t, y, d) << "\",";
            outputstr << endl << "\t\t" << "\"projection\":\"";
            outputstr << SimulationParameters::GetInstance()->getParameter("projection") << "\",";
            outputstr << endl << "\t\t" << "\"coordinates\":\"";
            lastLevel = 1;
        }
    }
}

void StringRepresentation::visit(FireNode* fn) {
    if (dumpMode == JSON_MODE)
    {
        if (fn->getFront()->getDomain()->getSimulationTime() >= fn->getFront()->getTime())
        {
            if (lastLevel == 2)
                outputstr << ' ';

            outputstr.precision(3);
            outputstr << fixed << fn->getX() << ',' << fn->getY() << ',' << fn->getZ();
            lastLevel = 2;
        }
        return;
    }
    
    if (dumpMode == FF_MODE)
    {
        for ( size_t k=0; k < currentLevel; k++ ) outputstr << '\t';
        outputstr << fn->toString() << endl;
    }
}

string StringRepresentation::dumpStringRepresentation() {
    
    if (SimulationParameters::GetInstance()->getParameter("dumpMode") == "json")
        dumpMode = JSON_MODE;
    else if (SimulationParameters::GetInstance()->getParameter("dumpMode") == "ff")
        dumpMode = FF_MODE;
    
	currentLevel = 0;
    lastLevel = -1;
    
	outputstr.str("");
	if (domain != 0) domain->accept(this);
    
    if (dumpMode == JSON_MODE)
    {
        if (lastLevel >= 2)
            outputstr << '"';
        if (lastLevel >= 1)
            outputstr << endl << '\t' << '}';
        if (lastLevel >= 0)
            outputstr << ']' << endl << '}' << endl;
    }
    
	return outputstr.str();
}

void StringRepresentation::increaseLevel(){
	currentLevel++;
}

void StringRepresentation::decreaseLevel(){
	currentLevel--;
}

size_t StringRepresentation::getLevel(){
	return currentLevel;
}

string StringRepresentation::toString(){
	ostringstream oss;
	oss<<"string representation";
	return oss.str();
}

}
