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
#ifndef POINTPATTERN_GUARD
#define POINTPATTERN_GUARD

#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
struct point
{
    point(const double cx, const double cy, const double cz, const int cl, const long cf = -1, const long cC = -1): x(cx), y (cy), z(cz), l(cl), faceID(cf), cellID(cC) {};
    double x, y, z;
    int l;
    long faceID, cellID;
    point& operator= ( point const& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        l = rhs.l;
        cellID = rhs.cellID;
        faceID = rhs.faceID;
        return *this;
    }
    point operator+ ( point const& rhs) const
    {
        return point(point::x + rhs.x, point::y + rhs.y, point::z + rhs.z, point::l);
    }

    double length ( ) const
    {
        return std::sqrt(x*x + y*y + z*z);
    }
};

point operator/ (const point& p,  double const& f);
point operator* (double const& f, const point& p);
point operator* (const point& p,  double const& f);



inline bool checkdistancecloserthan (point const& a, point const& b, double e)
{
    const double dx =fabs(a.x - b.x);
    const double dy =fabs(a.y - b.y);
    const double dz =fabs(a.z - b.z);

    return (dx < e && dy < e && dz < e);
};

class pointpattern
{
public:
    void addpoint(const int l, const double x, const double y, const double z);
    void addpointForCell(const double x, const double y, const double z, const int l, const long cf, const long cC);
    void print() const;
    void removeduplicates ( const double epsilon);
    void removeduplicates ( const double epsilon, pointpattern& p);

    inline void clear()
    {
        points.clear();
        indexShift.clear();
    }

    std::vector<point> points;
    std::map<unsigned int, long> indexShift;    // first is particle ID, second is new Index or -1
    

    friend std::ostream& operator << (std::ostream &f, const pointpattern& p)
    {
        if(p.points.empty())
            return f;
        int oldl = p.points[0].l;
        f << std::fixed;
        f << std::setprecision(15);
        for (   auto it = p.points.begin();
                it != p.points.end();
                ++it)
        {
            if(oldl != (*it).l )
            {
                f << "\n\n";
                oldl = (*it).l;
            }
            f << (*it).l << " " <<  std::setw(5)<< (*it).x << " " << std::setw(5) << (*it).y << " " << std::setw(5) << (*it).z << "\n";
        }

        return f;
    };

    friend std::ostream& operator >> (std::ostream &f, const pointpattern& p)
    {
        if(p.points.empty())
            return f;
    double xx = p.points[0].x;
    double yy = p.points[0].y;
    double zz = p.points[0].z;
    int oldf = p.points[0].faceID;
    int oldc = p.points[0].cellID;
        for (   auto it = p.points.begin();
                it != p.points.end();
                ++it)
        {
            if(oldf != (*it).faceID)
            {
               	f << oldc << " " <<  std::setw(8)<< xx << " " << std::setw(8) << yy << " " << std::setw(8) << zz << "\n\n\n";
                oldf = (*it).faceID;
                oldc = (*it).cellID;
                xx = (*it).x;
                yy = (*it).y;
                zz = (*it).z;
            }
            f << (*it).cellID << " " <<  std::setw(8) << (*it).x << " " << std::setw(8) << (*it).y << " " << std::setw(8) << (*it).z << "\n";
        }
        f << oldc << " " <<  std::setw(8)<< xx << " " << std::setw(8) << yy << " " << std::setw(8) << zz;

        return f;
    };
};

#endif
