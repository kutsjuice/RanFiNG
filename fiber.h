#ifndef FIBER_H
#define FIBER_H

#include <memory>

#include <Eigen/Dense>

typedef Eigen::Vector3d V3d;

struct FiberParam
{
    double diameter;
    V3d center;
    double length;
    double phi;
    double theta;
};

class Fiber
{
public:
    Fiber(double _diam, V3d _center, double _len, double _phi = 0.0, double _theta = 0.0, int index = -1);
    Fiber(FiberParam _p, int _index = -1)
        : Fiber(_p.diameter,
                _p.center,
                _p.length,
                _p.phi,
                _p.theta,
                _index){}

    int getBodyInd() const { return m_body_index; }
    void setBodyInd(int _index) {m_body_index = _index;}

    V3d getDirVec() const;
    V3d getSPoint() const {return m_pstart;}
    V3d getEPoint() const {return m_pend;}


    double getDiam() const {return m_diam;}

private:
    V3d m_pstart;
    V3d m_pend;
    double m_diam;

    int m_body_index;

};

#endif // FIBER_H
