#ifndef MODEL_H
#define MODEL_H

#include <memory>
#include <vector>
#include <unordered_set>
#include "fiber.h"
#include "body.h"
//#define uint

enum class COMBINATION_RESULT
{
    SUCCESS,
    FAILED_NOT_INTESECTED,
    FAILED_ALREADY_COMBINED,
    FAILED_UNKNOWN
};

enum class GENERATION_RESULT
{
    SUCCESS,
    FAILED
};

class Model
{
public:
    Model(std::vector<double> _size);

    int addFiber(FiberParam);

    bool isOneBody(int fib_1_ind, int fib_2_ind) const;

    bool intersect(int fib_1_ind, int fib_2_ind) const;

    double calcDistance(int fib_1_ind, int fib_2_ind, bool _inbounds = true) const;

    COMBINATION_RESULT combine(int fib_1_ind, int fib_2_ind, bool _check_intersection = false);

    GENERATION_RESULT createGeometry();

    int getAvailableFiberIndex() const;
private:
    std::unordered_map<int, Fiber> m_fibers;
    std::unordered_set<int> m_bodies{};
    std::vector<double> m_size;
    std::unordered_map<int, std::vector<int>> m_body_fibers;

    mutable int m_fiber_index{1};
};

#endif // MODEL_H
