#include <iostream>
#include <memory>
#include <cmath>


#include "fiber.h"
#include "model.h"
#include "body.h"

#include "gmsh.h"

using namespace std;

int testFibersCombination();
int testSimpleModel();
int testParallel();

int main(){
//    testFibersCombination();
    testSimpleModel();
//    testParallel();

}

int testFibersCombination()
{
    auto model = std::make_shared<Model>(std::vector<double>({10, 10, 10}));



    gmsh::initialize();

    gmsh::model::add("model");
    gmsh::logger::start();
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.035);
    // We first create two fiber:



    FiberParam fib_param_1;
    {
        fib_param_1.diameter = 0.1;
        fib_param_1.center = {0.0, 0.0, 0.0};
        fib_param_1.length = 0.5;
        fib_param_1.phi = 0.0;
        fib_param_1.theta = 0.0;
    }

    FiberParam fib_param_2;
    {
        fib_param_2.diameter = 0.1;
        fib_param_2.center = {0.0, 0.0, 0.0};
        fib_param_2.length = 0.5;
        fib_param_2.phi = 0.0;
        fib_param_2.theta = M_PI_2;
    }

    FiberParam fib_param_3;
    {
        fib_param_3.diameter = 0.1;
        fib_param_3.center = {0.0, 0.0, 0.0};
        fib_param_3.length = 0.5;
        fib_param_3.phi = M_PI_2;
        fib_param_3.theta = M_PI_2;
    }

    auto fib1 = model->addFiber(fib_param_1);
//    fib1->setBodyInd(1);
    auto fib2 = model->addFiber(fib_param_2);

    auto fib3 = model->addFiber(fib_param_3);

    auto d1 = model->calcDistance(fib1, fib2);
    auto d2 = model->calcDistance(fib1, fib3);
//    fib2->setBodyInd(1);

    auto res = model->isOneBody(fib1, fib2);

    auto res1 = model->intersect(fib1, fib2);
    auto res2 = model->intersect(fib1, fib3);

    gmsh::model::occ::synchronize();

    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;

    auto comb_res = model->combine(fib1, fib2);
    comb_res = model->combine(fib1, fib3);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::write("model.msh");
    gmsh::write("model.vtk");

    gmsh::logger::stop();

    gmsh::fltk::run();


    gmsh::finalize();
    return 0;
}


int testSimpleModel()
{
    auto model = std::make_shared<Model>(std::vector<double>({10, 10, 10}));



    gmsh::initialize();

    gmsh::model::add("model");
    gmsh::option::setNumber("General.Verbosity", 2);
//    gmsh::option::setNumber("")
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.5);
    gmsh::logger::start();

    model->createGeometry();


    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::write("model.msh");
    gmsh::write("model.vtk");

    gmsh::logger::stop();

    gmsh::fltk::run();


    gmsh::finalize();

    return 0;
}

int testParallel(){
    auto model = std::make_shared<Model>(std::vector<double>({10, 10, 10}));



    gmsh::initialize();

    gmsh::model::add("model");
    gmsh::logger::start();
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.035);
    // We first create two fiber:



    FiberParam fib_param_1;
    {
        fib_param_1.diameter = 0.1;
        fib_param_1.center = {0.0, 0.0, 0.0};
        fib_param_1.length = 0.5;
        fib_param_1.phi = 0.0;
        fib_param_1.theta = 0.0;
    }

    FiberParam fib_param_2;
    {
        fib_param_2.diameter = 0.1;
        fib_param_2.center = {0.05, 0.0, -0.51};
        fib_param_2.length = 0.5;
        fib_param_2.phi = 0.0;
        fib_param_2.theta = 0.0;
    }


    auto fib1 = model->addFiber(fib_param_1);
    //    fib1->setBodyInd(1);
    auto fib2 = model->addFiber(fib_param_2);

    auto d1 = model->calcDistance(fib1, fib2);
    //    fib2->setBodyInd(1);

    auto res = model->isOneBody(fib1, fib2);

    auto res1 = model->intersect(fib1, fib2);

    gmsh::model::occ::synchronize();

    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    if (res1){
        auto comb_res = model->combine(fib1, fib2);
    }

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::write("model.msh");
    gmsh::write("model.vtk");

    gmsh::logger::stop();

    gmsh::fltk::run();


    gmsh::finalize();
    return 0;

}
