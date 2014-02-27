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

#ifndef FFCONSTANTS_H_
#define FFCONSTANTS_H_

#include<limits>

using namespace std;

namespace libforefire {

class FFConstants {
public:
	FFConstants(){}
	virtual ~FFConstants(){}

	static const double Pi = 3.141592653589793;

	static const double epsilont = 1.e-3;
	static const double epsilonx = 1.e-3;
	static const double epsilonv = 1.e-3;

	static const double infinity(){
		return numeric_limits<double>::infinity();
	}

	static const size_t infiniteLoop = 1000000;

	static const double no_future_event = -1.;

	static const unsigned long int fdid = 100000000;
	static const unsigned long int ffid = 1000000;

};

}

#endif /* FFCONSTANTS_H_ */
