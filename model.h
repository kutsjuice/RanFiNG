#ifndef MODEL_H
#define MODEL_H

#include <memory>
#include <vector>

#include "fiber.h"
#include "body.h"

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

    std::shared_ptr<Fiber> addFiber(FiberParam);

    bool isOneBody(std::shared_ptr<Fiber> _fib1, std::shared_ptr<Fiber> _fib2) const;

    bool intersect(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2) const;

    double calcDistance(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2, bool _inbounds = true) const;

    COMBINATION_RESULT combine(std::weak_ptr<Fiber> _fib1, std::weak_ptr<Fiber> _fib2, bool _check_intersection = false);

    GENERATION_RESULT createGeometry();
private:
    std::vector<std::shared_ptr<Fiber>> m_fibers;

    std::vector<double> m_size;

    unsigned int m_body_index{1};
};

#endif // MODEL_H
