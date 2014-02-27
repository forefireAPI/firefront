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

#ifndef PARAMETERDATALAYER_H_
#define PARAMETERDATALAYER_H_

#include "DataLayer.h"

using namespace std;

namespace libforefire {

// TODO X, XY and XYZ classes

/*! \class XYZTDataLayer
 * \brief Data Layer for variables stored in NetCDF format
 *
 *  XYZTDataLayer allows to create layers related to properties stored
 *  in a NetCDF format up to 4D. This class extracts one specific
 *  NetCDF variable of the file to store it in memory.
 */
template<typename T> class ParameterDataLayer : public DataLayer<T> {

	double SWCornerX; /*!< origin in the X direction */
	double SWCornerY; /*!< origin in the Y direction */
	double SWCornerZ; /*!< origin in the Z direction */
	double startTime; /*!< origin of time */

	double NECornerX; /*!< origin in the X direction */
	double NECornerY; /*!< origin in the Y direction */
	double NECornerZ; /*!< origin in the Z direction */
	double endTime; /*!< origin of time */

	string key;  /*!< key to read */

	size_t nx; /*!< size of the array in the X direction */
	size_t ny; /*!< size of the array in the X direction */
	size_t nz; /*!< size of the array in the Z direction */
	size_t nt; /*!< size of the array in the T direction */
	size_t size; /*!< size of the array */

	double dx; /*!< increment in the X direction */
	double dy; /*!< increment in the Y direction */
	double dz; /*!< increment in the Z direction */
	double dt; /*!< increment in the T direction */

	FFArray<T>* array; /*!< pointer to the FFArray containing the data */

	/*! \brief Interpolation method: lowest order */
	T getNearestData(FFPoint, const double&);

	/*! \brief Interpolation method: multilinear */
	FFPoint posToIndices(FFPoint&);
	T bilinearInterp(FFPoint, const double&);

public:
	/*! \brief Default constructor */
	ParameterDataLayer() : DataLayer<T>() {
		startTime = 0;
			nx = 1;
			ny = 1;
			nz = 1;
			nt = 1;
			array = NULL;
			size = 1;
			SWCornerX = 0.;
			SWCornerY = 0.;
			SWCornerZ = 0.;
			NECornerX = 0.;
			NECornerY = 0.;
			NECornerZ = 0.;
			endTime = 0.;
			dx = 0.;
			dy = 0.;
			dz = 0.;
			dt = 0.;
			interp = bilinear;

	};
	/*! \brief Constructor with a sole value */
	ParameterDataLayer(string name, string lookupkey) :
		DataLayer<T>(name) {
		startTime = 0;
		nx = 1;
		ny = 1;
		nz = 1;
		nt = 1;
		array = new FFArray<T>(name, val, nx, ny, nz, nt);
		size = 1;
		SWCornerX = 0.;
		SWCornerY = 0.;
		SWCornerZ = 0.;
		NECornerX = 0.;
		NECornerY = 0.;
		NECornerZ = 0.;
		endTime = 0.;
		dx = 0.;
		dy = 0.;
		dz = 0.;
		dt = 0.;
		interp = bilinear;
		key = lookupkey;
	}

	~ParameterDataLayer();

	/*! \brief interpolation method enum type */
	enum InterpolationMethod {
		nearestData = 0,
		bilinear = 1
	};

	/*! \brief interpolation method */
	InterpolationMethod interp;

	/*! \brief obtains the value at a given position in the array */
	T getVal(size_t = 0, size_t = 0, size_t = 0, size_t = 0);
	/*! \brief obtains the value at a given position in the array */
	T getVal(size_t = 0);
	/*! \brief obtains the position in the array for given location and time */
	size_t getPos(FFPoint&, const double&);

	/*! \brief computes the value at a given firenode */
	T getValueAt(FireNode*);
	/*! \brief computes the value at a given location and time */
	T getValueAt(FFPoint, const double&);
	/*! \brief sets the value at a given location and time */
	void setValueAt(FFPoint, const double&);
	/*! \brief directly stores the desired values in a given array */
	size_t getValuesAt(FireNode*, PropagationModel*, size_t);
	/*! \brief directly stores the desired values in a given array */
	size_t getValuesAt(FFPoint, const double&, FluxModel*, size_t);

	/*! \brief getter to the desired array at surface for a given time (should not be used) */
	void getMatrix(FFArray<T>**, const double&);
	/*! \brief stores array from a given array (should not be used) */
	void setMatrix(string&, double*, const size_t&, size_t&, const double&);

	/*! \brief print the related array (should not be used) */
	string print();
	string print2D(size_t, size_t);
	void dumpAsBinary(string, const double&
			, FFPoint&, FFPoint&, size_t&, size_t&);

};

template<typename T>
XYZTDataLayer<T>::~XYZTDataLayer() {
	delete array;
}

template<typename T>
T XYZTDataLayer<T>::getVal(size_t pos){
	return (*array)(pos);
}

template<typename T>
T XYZTDataLayer<T>::getVal(size_t i, size_t j, size_t k, size_t t){
	return (*array)(i, j, k, t);
}

template<typename T>
size_t XYZTDataLayer<T>::getPos(FFPoint& loc, const double& t){
	int i = nx>1 ? (int) (loc.getX()-SWCornerX)/dx : 0;
	int j = ny>1 ? (int) (loc.getY()-SWCornerY)/dy : 0;
	int k = nz>1 ? (int) (loc.getZ()-SWCornerZ)/dz : 0;
	int tt = nt>1 ? (int) (t-startTime)/dt : 0;

	if ( i < 0 or i > (int) nx-1 ){
		cout<<loc.getX()<<" is not within the 'x' domain"
				<<" of validity of array for "<<this->getKey()
				<<" ("<<SWCornerX<<"->"<<NECornerX<<")"<< endl;
		return 0;
	}
	if ( j < 0 or j > (int) ny-1 ){
		cout<<loc.getY()<<" is not within the 'y' domain"
				<<" of validity of array for "<<this->getKey()
				<<" ("<<SWCornerY<<"->"<<NECornerY<<")"<< endl;
		return 0;
	}
	if ( k < 0 or k > (int) nz-1 ){
		cout<<loc.getZ()<<" is not within the 'z' domain"
				<<" of validity of array for "<<this->getKey()
				<<" ("<<SWCornerZ<<"->"<<NECornerZ<<")"<< endl;
		return 0;
	}
	if ( tt < 0 or tt > (int) nt-1 ){
		cout<<t<<" is not within the 't' domain"
				<<" of validity of array for "<<this->getKey()
				<<" ("<<startTime<<"->"<<endTime<<")"<< endl;
		return 0;
	}
	return (size_t) i*ny*nz*nt + j*nz*nt + k*nt + tt;
}

template<typename T>
T XYZTDataLayer<T>::getValueAt(FireNode* fn){
	if ( size == 1 ) return (*array)(0, 0, 0, 0);
	FFPoint loc = fn->getLoc();
	double t = fn->getTime();
	if ( interp == nearestData ) return getNearestData(fn->getLoc(), t);
	if ( interp == bilinear ) return bilinearInterp(fn->getLoc(), t);
	cout<<"WARNING: unknown interpolation method in "
			<<"Array2DdataLayer<T>::getValueAt(FireNode*)"<<endl;
	return (T) 0;
}

template<typename T>
T XYZTDataLayer<T>::getValueAt(FFPoint loc, const double& t){
	if ( size == 1 ) return (*array)(0, 0, 0, 0);
	if ( interp == nearestData ) return getNearestData(loc, t);
	if ( interp == bilinear ) return bilinearInterp(loc, t);
	cout<<"WARNING: unknown interpolation method in "
			<<"Array2DdataLayer<T>::getValueAt(FireNode*)"<<endl;
	return (T) ;
}

template<typename T>
void XYZTDataLayer<T>::setValueAt(FFPoint loc, const double& t){
	if ( size == 1 ) return (*array)(0, 0, 0, 0);
	if ( interp == nearestData ) return getNearestData(loc, t);
	if ( interp == bilinear ) return bilinearInterp(loc, t);
	cout<<"WARNING: unknown interpolation method in "
			<<"Array2DdataLayer<T>::getValueAt(FireNode*)"<<endl;
	return (T) 0;
}


template<typename T>
size_t XYZTDataLayer<T>::getValuesAt(
		FireNode* fn, PropagationModel* model, size_t curItem){
	return 0;
}

template<typename T>
size_t XYZTDataLayer<T>::getValuesAt(
		FFPoint loc, const double& t, FluxModel* model, size_t curItem){
	return 0;
}

template<typename T>
T XYZTDataLayer<T>::getNearestData(FFPoint loc, const double& time){
	return (*array)(getPos(loc,time));
}

template<typename T>
FFPoint XYZTDataLayer<T>::posToIndices(FFPoint& loc){
	FFPoint ind;
	ind.setX((loc.getX()-SWCornerX)/dx-0.5);
	if ( ny > 1 ) ind.setY((loc.getY()-SWCornerY)/dy-0.5);
	if ( nz > 1 ) ind.setZ((loc.getZ()-SWCornerZ)/dz-0.5);
	return ind;
}

template<typename T>
T XYZTDataLayer<T>::bilinearInterp(FFPoint loc, const double& t){
	/* This method implements a multilinear interpolation */
	/* At this time no interpolation in 'Z' is done */

	/* searching the coordinates of the nodes around */
	FFPoint indices = posToIndices(loc);

	double ud = indices.getX() + FFConstants::epsilonx;
	double vd = indices.getY() + FFConstants::epsilonx;

	int uu = (int) ceil(ud-1);
	int vv = (int) ceil(vd-1);

	if ( uu < 0 ) uu = 0;
	if ( uu > ((int) nx) - 2 ) uu = ((int) nx) - 2;
	if ( vv < 0 ) vv = 0;
	if ( vv > ((int) ny) - 2 ) vv = ((int) ny) - 2;

	double udif = ud - ((double) uu);
	double vdif = vd - ((double) vv);

	double csw = (1.-udif) * (1.-vdif);
	double cse = udif * (1 - vdif);
	double cnw = (1 - udif) * vdif;
	double cne = udif * vdif;

	if ( nt>1 ) {
		int it = (int) (t-startTime)/dt;

		T tsw1 = getVal(uu,vv,0,it);
		T tnw1 = getVal(uu,vv+1,0,it);
		T tne1 = getVal(uu+1,vv+1,0,it);
		T tse1 = getVal(uu+1,vv,0,it);

		T tsw2 = getVal(uu,vv,0,it+1);
		T tnw2 = getVal(uu,vv+1,0,it+1);
		T tne2 = getVal(uu+1,vv+1,0,it+1);
		T tse2 = getVal(uu+1,vv,0,it+1);

		T val1 = csw*tsw1 + cse*tse1 + cnw*tnw1 + cne*tne1;
		T val2 = csw*tsw2 + cse*tse2 + cnw*tnw2 + cne*tne2;

		/* interpolation in time */
		double at = ( t - it*dt )/dt;
		return at*val1 + (1.-at)*val2;

	} else {
		T tsw = getVal(uu,vv);
		T tnw = getVal(uu,vv+1);
		T tne = getVal(uu+1,vv+1);
		T tse = getVal(uu+1,vv);
		T val = csw*tsw + cse*tse + cnw*tnw + cne*tne;
		return val;
	}
}

template<typename T>
void XYZTDataLayer<T>::getMatrix(
		FFArray<T>** matrix, const double& time){
	*matrix = array;
}

template<typename T>
void XYZTDataLayer<T>::setMatrix(string& mname, double* inMatrix
		, const size_t& sizein, size_t& sizeout, const double& time){
	// writing the incoming array in matrix
	// should not be done with this type of layer
}

template<typename T>
string XYZTDataLayer<T>::print(){
	return print2D(0,0);
}

template<typename T>
string XYZTDataLayer<T>::print2D(size_t k, size_t t){
	ostringstream oss;
	size_t jj;
	for ( int j = ny-1; j >= 0; j -= 10 ){
		jj = (size_t) j;
		for ( size_t i = 0; i < nx; i += 10 ){
			oss<<getVal(i, j, k, t)<<" ";
		}
		oss<<endl;
	}
	return oss.str();
}

template<typename T>
void XYZTDataLayer<T>::dumpAsBinary(string filename, const double& time
		, FFPoint& SWC, FFPoint& NEC, size_t& nnx, size_t& nny){

	T vals[nnx*nny];

	double ddx = (NEC.getX()-SWC.getX())/nnx;
	double ddy = (NEC.getY()-SWC.getY())/nny;

	FFPoint loc;
	loc.setX(SWC.getX());
	for ( size_t i = 0; i < nnx; i++ ) {
		loc.setY(SWC.getY());
		for ( size_t j = 0; j < nny; j++ ) {
			vals[i*nny+j] = getValueAt(loc, time);
			loc.setY(loc.getY()+ddy);
		}
		loc.setX(loc.getX()+ddx);
	}

	/* writing the interpolated matrix in a binary file */
	ostringstream outputfile;
	outputfile<<filename<<"."<<this->getKey();
	ofstream FileOut(outputfile.str().c_str(), ios_base::binary);
	FileOut.write(reinterpret_cast<const char*>(&nnx), sizeof(size_t));
	FileOut.write(reinterpret_cast<const char*>(&nny), sizeof(size_t));
	FileOut.write(reinterpret_cast<const char*>(&vals), sizeof(vals));
	FileOut.close();
}

}

#endif  /* NCXYZTDATALAYER_H_ */
