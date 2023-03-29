/*
Copyright 2018 Simon Weis and Philipp Schoenhoefer

This file is part of Pomelo.

Pomelo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pomelo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pomelo.  If not, see <http://www.gnu.org/licenses/>.

The development of Pomelo took place at the Friedrich-Alexander University of Erlangen and was funded by the German Research Foundation (DFG) Forschergruppe FOR1548 "Geometry and Physics of Spatial Random Systems" (GPSRS).
*/
#ifndef WRITERPOLY_H_GUARD_12345
#define WRITERPOLY_H_GUARD_12345

#include <iomanip>
#include <vector>
#include <map>
#include <cmath>

#include "IWriter.hpp"

class writerpoly : public IWriter
{
public:
    writerpoly()
    {
    };
    writerpoly(IWriter const& other)
    {
        faceCellMap = other.faceCellMap;
        p = other.p;
        faces = other.faces;
    }

    void print(std::ostream& f) const
    {
        f << "POINTS" << std::endl;
        f << std::fixed;
        for(auto it =  p.points.begin();
                it != p.points.end();
                ++it)
        {
            f << it->l << ":    " <<  std::setprecision(10) << it->x << " " << std::setprecision(10) << it-> y << " " << std::setprecision(10) << it->z<< "\n";
        }

        f << "POLYS" <<  std::endl;
        int removedFaces = 0;
        for (
            auto it = faces.begin();
            it != faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = faceCellMap.at(faceID);
            std::vector<unsigned int> testing;
            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
                bool doppelt = false;
                for(unsigned int kk = 0; kk < testing.size(); kk++ )
                {
                    if(testing[kk] == (*it2) )
                    {
                        doppelt = true;
                        break;
                    }
                }
                if(doppelt) continue;
                testing.push_back( (*it2) );
            }
            if(2 < testing.size())
            {
                f << it->first-removedFaces << ":    ";
                for(unsigned int kk = 0; kk < testing.size(); kk++ )
                {
                    f << testing[kk] << " ";
                }
                f << "< c(0, 0, 0, " << cellID << ")" << std::endl;
            }
            else removedFaces++;
  //          else std::cout << testing.size() << " " << cellID << std::endl;
        }

        f << "END";
    };


};

#endif
