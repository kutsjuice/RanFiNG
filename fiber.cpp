#include "fiber.h"
#include <math.h>
#include <iostream>
#include <vector>
#include "gmsh.h"

Fiber::Fiber(double _diam, V3d _center, double _len, double _phi, double _theta, long int _index)
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



    gmsh::model::occ::addCylinder(  m_pstart[0],
                                  m_pstart[1],
                                  m_pstart[2],
                                  m_pend[0] - m_pstart[0],
                                  m_pend[1] - m_pstart[1],
                                  m_pend[2] - m_pstart[2],
                                  m_diam/2,
                                  _index+1);
    gmsh::model::occ::addSphere(m_pstart[0],
                                m_pstart[1],
                                m_pstart[2],
                                m_diam/2,
                                _index+2);
    gmsh::model::occ::addSphere(m_pend[0],
                                m_pend[1],
                                m_pend[2],
                                m_diam/2,
                                _index+3);

    gmsh::vectorpair ov;
    std::vector<gmsh::vectorpair> ovv;
//    gmsh::model::occ::fuse({{3, fibA->getBodyInd()}}, {{3, fibB->getBodyInd()}}, ov, ovv);
    gmsh::model::occ::fuse({{3, _index+2}, {3, _index+3}}, {{3, _index+1}}, ov, ovv, _index);
    m_body_index = _index;

}

V3d Fiber::getDirVec() const
{
    return (m_pend - m_pstart);
}

