/*  
 *  Copyright 2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#pragma once

#include <boost/python.hpp>

#include "medial_axis_pocket.hpp"

namespace ovd {
namespace pyovd {
    
/// python wrapper for medial_axis_pocket
class medial_axis_pocket_py : public medial_axis_pocket {
public:
    /// create medial_axis_pocket object
    medial_axis_pocket_py(HEGraph& gi): medial_axis_pocket(gi) {}
    /// return list of MICs \todo replaced by get_mic_components()
    boost::python::list py_get_mic_list() {
        boost::python::list out;
        BOOST_FOREACH(MIC mic, mic_list) {
            boost::python::list m;
            m.append( mic.c2 );          //0
            m.append( mic.r2 );          //1
            m.append( mic.t1 );         //2
            m.append( mic.t2 );         //3
            m.append( mic.t3 );         //4
            m.append( mic.t4 );         //5
            m.append( mic.c1 );          //6
            m.append( mic.r1 );          //7
            m.append( mic.new_branch ); //8
            m.append( mic.c_prev );     //9
            m.append( mic.r_prev );     //10
            out.append(m);
        }
        return out;
    }
    
    /// return list of MICs. one list per connected component in the graph
    boost::python::list py_get_mic_components() {
        boost::python::list out;
        BOOST_FOREACH(MICList miclist, ma_components) {
            boost::python::list mic_out;
            BOOST_FOREACH(MIC mic, miclist) {
                boost::python::list m;
                m.append( mic.c2 );          //0
                m.append( mic.r2 );          //1
                m.append( mic.t1 );         //2
                m.append( mic.t2 );         //3
                m.append( mic.t3 );         //4
                m.append( mic.t4 );         //5
                m.append( mic.c1 );          //6
                m.append( mic.r1 );          //7
                m.append( mic.new_branch ); //8
                m.append( mic.c_prev );     //9
                m.append( mic.r_prev );     //10
                mic_out.append(m);
            }
            out.append(mic_out);
            //return out;
            
        }
        return out;
    }
};

} // pyovd
} // end ovd namespace
// end file 
