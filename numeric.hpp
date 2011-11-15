/* 
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenVoronoi.
 *
 *  OpenVoronoi is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenVoronoi is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenVoronoi.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NUMERIC_HPP
#define NUMERIC_HPP

#include <qd/qd_real.h> 

namespace ovd {

namespace numeric {
    double chop8(double a);
    double chop(double val);
    qd_real chop(qd_real val);
    
    template<class Scalar>
    Scalar sq( Scalar x) {return x*x;}
} // numeric
} // ovd

#endif
