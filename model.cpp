#include "model.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <random>

#include "gmsh.h"
#include <Eigen/Dense>
using namespace Eigen;

//! \brief distancePoint2Line calculates distance from point to line;
//! \param point
//! \param line_start
//! \param line_dir
//! \param _check_in_bounds if set, then function checks if the projection of point lays in bounds of line; if not - then returns Inf
//! \return distance from point to line; if _check_in_bounds set, then may return Inf, if point out of bounds;
//!
double distancePoint2Line(V3d point, V3d line_start, V3d line_dir, bool _check_in_bounds = false){

    auto SA = point - line_start;
    auto SD = line_dir;


    auto proj = SA.dot(SD) / SD.norm();

    if( _check_in_bounds && ((proj < 0) || (proj > line_dir.norm()))){
        return std::numeric_limits<double>::infinity();
    }
    auto orth = (SA - SA*proj);
    return orth.norm();

}

//!
//! \brief Model::Model
//! \param _size
//!
Model::Model(std::vector<double> _size): m_size{_size}
{

}

//!
//! \brief Model::addFiber
//! \param _params
//! \return
//!
int Model::addFiber(FiberParam params_)
{
    auto s_ind = getAvailableFiberIndex();
    m_fibers.emplace(s_ind, Fiber(params_));

    m_bodies.insert(m_fibers.at(s_ind).getBodyInd());
    m_body_fibers[s_ind] = {s_ind};
    return s_ind;
}
//!
//! \brief Model::isOneBody
//! \param _fib1
//! \param _fib2
//! \return
//!
bool Model::isOneBody(int fib1, int fib2) const
{
    return m_fibers.at(fib1).getBodyInd() == m_fibers.at(fib2).getBodyInd();
}
//!
//! \brief Model::intersect
//! \param _fib1
//! \param _fib2
//! \return
//!
bool Model::intersect(int fib1, int fib2) const
{
    auto distance = calcDistance(fib1, fib2, true);
    auto& fibA = m_fibers.at(fib1);
    auto& fibB = m_fibers.at(fib2);


    auto res = (distance <= (fibA.getDiam()/2 + fibB.getDiam()/2));
    return res;

}
//!
//! \brief Model::calcDistance
//! \param _fib1
//! \param _fib2
//! \param _inbounds
//! \return
//!
double Model::calcDistance(int fib1, int fib2, bool _inbounds) const
{
    auto& fibA = m_fibers.at(fib1);
    auto& fibB = m_fibers.at(fib2);
    // check parallel

    double eps = 0.001;
    double dist;    

    Vector3d AB = fibA.getDirVec();
    Vector3d CD = fibB.getDirVec();
    Vector3d AC = fibB.getSPoint() - fibA.getSPoint();
    Vector3d CB = fibA.getEPoint() - fibB.getEPoint();
    Vector3d AD = fibB.getEPoint() - fibA.getSPoint();
    Vector3d DB = fibA.getEPoint() - fibB.getEPoint();

    Vector3d n = AB.cross(CD);
    bool cover = false;
    if(n.norm() < eps){
        // case of parellel fibers
        //
        //                     * D
        //                    /
        //          B *      /
        //           /      /
        //          /      /
        //         /      * C
        //        /
        //     A *
        //
        // Algrithm:
        // 1. calculate projection of AC on AB Pr(AC)_AB
        // 2. calculae orthogonal vector Ort(AC)_|_AB = AC - Pr(AC)_AB
        // 3. Lenth of Ort(AC)_AB equal to distnace between AB and CD
        // 4. The shortest distance will eqal to Ort(AC)_AB these vectors cover each other at list partially. For that At least one projection of a first vector's point  should lay on second fibers


        Vector3d proj = AC.dot(AB)/AB.squaredNorm() * AB;
        Vector3d orth = (AC - proj);
        dist = orth.norm();
        //TODO check bounds
        double t = AC.dot(AB)/AB.squaredNorm();

        if ((t < 0) || (t > 1)){

            t = AD.dot(AB)/AB.squaredNorm();

            if ((t < 0) || (t > 1)){

                t = CB.dot(CD)/CD.squaredNorm();
                if ((t >= 0) && (t <= 1)){
                    cover = true;
                }
            }else{
                cover = true;
            }
        }else{
            cover = true;
        }
        // case, when fibers do not intersect each other


    } else{
        //not parallel
//        n = n / n.norm();
        dist = abs(n.dot(AC)) / n.norm();
        Vector3d d2 = CD.cross(n);
        Vector3d d1 = AB.cross(n);
        auto t1 = AC.dot(d2) / AB.dot(d2);
        auto t2 = - AC.dot(d1) / CD.dot(d1);

        //bounds check
        if((t1 >= 0.0) && (t1 <= 1.0) && (t2 >= 0.0) && (t2 <= 1.0)){
            cover = true;
        }
    }
    if (!cover && _inbounds){
        // calc 4 distances (AC, AD, CB, DB)
        // shortest distance will be the lowest value
        Vector<double, 4> distances {
            AC.norm(),
            CB.norm(),
            AD.norm(),
            DB.norm(),
        };
        dist = distances.minCoeff();
    }
    return dist;

}

//!
//! \brief Model::combine
//! \param _fib1
//! \param _fib2
//! \param _check_intersection
//! \return
//!
COMBINATION_RESULT Model::combine(int fib1, int fib2, bool _check_intersection){
    auto& fibA = m_fibers.at(fib1);
    auto& fibB = m_fibers.at(fib2);


    if(_check_intersection && !intersect(fib1, fib2)){
//        std::cout << "Fibers (" << fibA->getBodyInd() << ";" << fibB->getBodyInd() << ") don't intersect!"<< std::endl;

        return COMBINATION_RESULT::FAILED_NOT_INTESECTED;
    }

    if(fibA.getBodyInd() == fibB.getBodyInd()){
        return COMBINATION_RESULT::FAILED_ALREADY_COMBINED;
    }
    std::cout << "Fibers (" << fib1 << ";" << fib2 << ") Intersect!"<< std::endl;
    try
    {

        std::vector<std::pair<int, int> > s_new_bodies;
        std::vector<std::vector<std::pair<int, int> > > ovv;
        gmsh::vectorpair s_3d_entities_before;

        int fib_a_body_ind = fibA.getBodyInd();
        int fib_b_body_ind = fibB.getBodyInd();

        gmsh::model::occ::getEntities(s_3d_entities_before, 3);

        int s_index = 0;
        for (auto s_dimtag = s_3d_entities_before.begin(); s_dimtag != s_3d_entities_before.end(); s_dimtag++){
            if (s_dimtag->second > s_index){
                s_index = s_dimtag->second;
            }
        }
        s_index++;

        gmsh::model::occ::fuse({{3, fib_a_body_ind}}, {{3, fib_b_body_ind}}, s_new_bodies, ovv, s_index);


        m_bodies.erase(fib_a_body_ind);
        m_bodies.erase(fib_b_body_ind);
        m_bodies.insert(s_index);

        gmsh::vectorpair s_3d_entities_after;

        gmsh::model::occ::getEntities(s_3d_entities_after, 3);

        gmsh::vectorpair s_diff;

        std::set_difference(s_3d_entities_after.begin(),
                            s_3d_entities_after.end(),
                            s_3d_entities_before.begin(),
                            s_3d_entities_before.end(),
                            std::inserter(s_diff, s_diff.begin()));


        if (s_new_bodies.size() != 1 ){
            std::cerr << "diff doesn't correspond to s_index" << std::endl;
            return COMBINATION_RESULT::FAILED_UNKNOWN;
        }

        if (s_new_bodies.at(0).second != s_index){
            std::cerr << "more than a single diff" << std::endl;
            return COMBINATION_RESULT::FAILED_UNKNOWN;
        }


        m_body_fibers[s_index] = m_body_fibers.at(fib_a_body_ind);
        m_body_fibers[s_index].insert(m_body_fibers[s_index].end(), m_body_fibers.at(fib_b_body_ind).begin(), m_body_fibers.at(fib_b_body_ind).end());

        m_body_fibers.erase(fib_a_body_ind);
        m_body_fibers.erase(fib_b_body_ind);

        for (auto it = m_body_fibers[s_index].begin(); it != m_body_fibers[s_index].end(); it ++ ){
            m_fibers.at(*it).setBodyInd(s_index);
        }



        return COMBINATION_RESULT::SUCCESS;
    }

    catch(...){
        std::cerr << "Some error occured!" << std::endl;
        return COMBINATION_RESULT::FAILED_UNKNOWN;
    }
}
//!
//! \brief Model::createGeometry
//! \return
//!
GENERATION_RESULT Model::createGeometry(){
    double diam{1.5};
    double length{4.5};
    int fiber_number{30};


    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(100); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> X(-m_size[0]/2, m_size[0]/2);
    std::uniform_real_distribution<> Y(-m_size[1]/2, m_size[1]/2);
    std::uniform_real_distribution<> Z(-m_size[2]/2, m_size[2]/2);
    std::uniform_real_distribution<> THETA(-M_PI/3, M_PI/3);
    std::normal_distribution<> PHI(0.0, M_PI/3);


    FiberParam fib_param;
    fib_param.diameter = diam;
    fib_param.length = length;



    for(int i = 0; i < fiber_number; i++){

        fib_param.phi = PHI(gen);
        fib_param.theta = THETA(gen);
        fib_param.center = {X(gen), Y(gen), Z(gen)};

        addFiber(fib_param);
    }

    for(auto it1 = m_fibers.begin(); it1 != m_fibers.end(); it1++){
        for( auto it2 = std::next(it1); it2 != m_fibers.end(); it2++){
            combine(it1->first, it2->first, true);
        }
    }

    gmsh::vectorpair s_entities;
    gmsh::model::occ::getEntities(s_entities, 3);


    int s_box_ind{0};
    for (auto it = s_entities.begin(); it != s_entities.end(); it++){
        if (it->second > s_box_ind) s_box_ind = it->second;
    }
    s_box_ind++;
    gmsh::vectorpair ov;
    std::vector<gmsh::vectorpair> ovv;
    gmsh::model::occ::addBox(-m_size[0]/2, -m_size[1]/2, -m_size[2]/2, m_size[0], m_size[1], m_size[2], s_box_ind);

//    gmsh::model::occ::
    gmsh::model::occ::intersect(s_entities, {{3, s_box_ind}},ov,ovv);
    std::cout << "intersection performed" << std::endl;

    return GENERATION_RESULT::SUCCESS;



}

int Model::getAvailableFiberIndex() const
{
    while(m_fibers.find(m_fiber_index) != m_fibers.end()){
        m_fiber_index++;
    }
    return m_fiber_index;
}
