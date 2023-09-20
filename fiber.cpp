#include "fiber.h"
#include <math.h>
#include <iostream>
#include "gmsh.h"

Fiber::Fiber(double _diam, V3d _center, double _len, double _phi, double _theta)
{
    double dx = 0.5 * _len * sin(_theta) * cos(_phi);
    double dy = 0.5 * _len * sin(_theta) * sin(_phi);
    double dz = 0.5 * _len * cos(_theta);

    m_pstart(0) = _center(0) - dx;
    m_pstart(1) = _center(1) - dy;
    m_pstart(2) = _center(2) - dz;

    m_pend(0) = _center(0) + dx;
    m_pend(1) = _center(1) + dy;
    m_pend(2) = _center(2) + dz;

    m_diam = _diam;



    m_body_index = gmsh::model::occ::addCylinder(  m_pstart[0],
                                  m_pstart[1],
                                  m_pstart[2],
                                  m_pend[0] - m_pstart[0],
                                  m_pend[1] - m_pstart[1],
                                  m_pend[2] - m_pstart[2],
                                  m_diam);


}

V3d Fiber::getDirVec() const
{
    return (m_pend - m_pstart);
}

