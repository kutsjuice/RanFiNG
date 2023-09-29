#include "model.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <random>

#include "gmsh.h"

//!
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


Model::Model(std::vector<double> _size): m_size{_size}
{

}


std::shared_ptr<Fiber> Model::addFiber(FiberParam _params)
{
    m_fibers.push_back(std::make_shared<Fiber>(_params, m_fiber_index));
    m_fiber_index += 5;
    return m_fibers.back();
}

bool Model::isOneBody(std::shared_ptr<Fiber> _fib1, std::shared_ptr<Fiber> _fib2) const
{
    return _fib1->getBodyInd() == _fib2->getBodyInd();
}

bool Model::intersect(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2) const
{
    auto distance = calcDistance(_fib1, _fib2, true);
    auto fibA = _fib1.lock();
    auto fibB = _fib2.lock();


    auto res = (distance <= (fibA->getDiam()/2 + fibB->getDiam()/2));
    return res;

}

double Model::calcDistance(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2, bool _inbounds) const
{
    auto fibA = _fib1.lock();
    auto fibB = _fib2.lock();
    auto n = fibA->getDirVec().cross(fibB->getDirVec());

    auto d = abs(n.dot(fibA->getSPoint() - fibB->getSPoint())) / n.norm();
    auto tA = fibA->getDirVec().cross(n).dot(fibB->getSPoint() - fibA->getSPoint()) / n.dot(n);
    auto tB = fibB->getDirVec().cross(n).dot(fibB->getSPoint() - fibA->getSPoint()) / n.dot(n);

    if(_inbounds && ( (tA < 0.0) || (tA > 1.0) || (tB < 0.0) || (tB > 1.0) ) ){
        d = std::numeric_limits<double>::infinity();

        auto d1 = distancePoint2Line(fibA->getSPoint(), fibB->getSPoint(), fibB->getDirVec(), true);
        if(d1 < d) d = d1;
        d1 = distancePoint2Line(fibA->getEPoint(), fibB->getSPoint(), fibB->getDirVec(), true);
        if(d1 < d) d = d1;
        d1 = distancePoint2Line(fibB->getSPoint(), fibA->getSPoint(), fibA->getDirVec(), true);
        if(d1 < d) d = d1;
        d1 = distancePoint2Line(fibB->getEPoint(), fibA->getSPoint(), fibA->getDirVec(), true);
        if(d1 < d) d = d1;
    }

    return d;

}


COMBINATION_RESULT Model::combine(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2, bool _check_intersection){
    auto fibA = _fib1.lock();
    auto fibB = _fib2.lock();


    if(_check_intersection && !intersect(fibA, fibB)){
//        std::cout << "Fibers (" << fibA->getBodyInd() << ";" << fibB->getBodyInd() << ") don't intersect!"<< std::endl;

        return COMBINATION_RESULT::FAILED_NOT_INTESECTED;
    }

    if(fibA->getBodyInd() == fibB->getBodyInd()){
        return COMBINATION_RESULT::FAILED_ALREADY_COMBINED;
    }
    std::cout << "Fibers (" << fibA->getBodyInd() << ";" << fibB->getBodyInd() << ") Intersect!"<< std::endl;
    try
    {

        std::vector<std::pair<int, int> > ov;
        std::vector<std::vector<std::pair<int, int> > > ovv;
        gmsh::model::occ::fuse({{3, fibA->getBodyInd()}}, {{3, fibB->getBodyInd()}}, ov, ovv);

        fibA->setBodyInd(ov[0].second);
        fibB->setBodyInd(ov[0].second);


        return COMBINATION_RESULT::SUCCESS;

    }

    catch(...){
        std::cerr << "Some error occured!" << std::endl;
        return COMBINATION_RESULT::FAILED_UNKNOWN;
    }
}

GENERATION_RESULT Model::createGeometry(){
    double diam{1.5};
    double length{4.5};
    int fiber_number{30};


    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> X(-m_size[0]/2, m_size[0]/2);
    std::uniform_real_distribution<> Y(-m_size[1]/2, m_size[1]/2);
    std::uniform_real_distribution<> Z(-m_size[2]/2, m_size[2]/2);
    std::uniform_real_distribution<> THETA(-M_PI/3, M_PI/3);


    FiberParam fib_param;
    fib_param.diameter = diam;
    fib_param.length = length;



    for(int i = 0; i < fiber_number; i++){

        fib_param.phi = 0.0;
        fib_param.theta = THETA(gen);
        fib_param.center = {X(gen), Y(gen), Z(gen)};

        addFiber(fib_param);
    }

    for(auto it1 = m_fibers.begin(); it1 != m_fibers.end(); it1++){
        for( auto it2 = std::next(it1); it2 != m_fibers.end(); it2++){
            combine(*it1, *it2, true);
        }
    }

    return GENERATION_RESULT::SUCCESS;



}
