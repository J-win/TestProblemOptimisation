#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <string>

struct Problem
{
    double a;
    double b;
    double xopt;
    double zopt;
    double r;

    const double pi = 3.14159265358979323846;
    const double e = 2.71828182845904523536;

    double summ(double sum);

    virtual double func(const double x) = 0;
    virtual std::string getNameFunc() = 0;
    void getInfo();
};

struct SchwefelProblem : public Problem
{
    SchwefelProblem(const double r_ = 2.0);
    double func(const double x) override;
    std::string getNameFunc() override;
};

struct Task
{
    Problem* data;
    Task();
    ~Task();
};