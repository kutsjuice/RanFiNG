#include "fiber.h"
#include <math.h>
#include <iostream>
#include <vector>
#include "gmsh.h"

Fiber::Fiber(double _diam, V3d _center, double _len, double _phi, double _theta, int _index)
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



    auto cyl_ind = gmsh::model::occ::addCylinder(  m_pstart[0],
                                  m_pstart[1],
                                  m_pstart[2],
                                  m_pend[0] - m_pstart[0],
                                  m_pend[1] - m_pstart[1],
                                  m_pend[2] - m_pstart[2],
                                  m_diam/2);
    auto sph1_ind = gmsh::model::occ::addSphere(m_pstart[0],
                                m_pstart[1],
                                m_pstart[2],
                                m_diam/2);
    auto sph2_ind = gmsh::model::occ::addSphere(m_pend[0],
                                m_pend[1],
                                m_pend[2],
                                m_diam/2);

    gmsh::vectorpair ov;
    std::vector<gmsh::vectorpair> ovv;
//    gmsh::model::occ::fuse({{3, fibA->getBodyInd()}}, {{3, fibB->getBodyInd()}}, ov, ovv);
    gmsh::model::occ::fuse({{3, sph1_ind}, {3, sph2_ind}}, {{3, cyl_ind}}, ov, ovv, _index);
//    gmsh::model::occ::remove({{3,_index+1},{3,_index+2}, {3,_index+3}});

    gmsh::vectorpair s_3d_entities;
    gmsh::model::occ::getEntities(s_3d_entities, 3);

    m_body_index = s_3d_entities.back().second;

}

V3d Fiber::getDirVec() const
{
    return (m_pend - m_pstart);
}


