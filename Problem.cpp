#include "Problem.h"

double Problem::summ(double sum) {
    for (int i = 1; i <= 100000; i++) {
        sum += sin(sin(sin(i))) * sin(sin(sin(i))) + cos(sin(sin(i))) * cos(sin(sin(i)));
    }
    sum -= 100000;
    return sum;
}

void Problem::getInfo()
{
    std::cout << getNameFunc() << std::endl;
    std::cout << "Optimum arg min f = " << xopt << std::endl;
    std::cout << "Optimum min f = " << zopt << std::endl;
}

SchwefelProblem::SchwefelProblem(const double r_)
{
    a = -500.0;
    b = 500.0;
    xopt = 420.97;
    zopt = -418.9829;
    r = r_;
}

double SchwefelProblem::func(const double x)
{
    return -1.0 * x * sin(sqrt(fabs(x))) + summ(0.0);
}

std::string SchwefelProblem::getNameFunc()
{
    return "Schwefel problem dimension one";
}

Task::Task()
{
    data = new SchwefelProblem;
}

Task::~Task()
{
    delete[] data;
}