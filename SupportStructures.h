#pragma once

#include <list>
#include <queue>
#include <vector>
#include "Problem.h"

struct point
{
    double x;
    double z;
    point(const double x_, const double z_) : x(x_), z(z_) {}
    point(const point& p) : x(p.x), z(p.z) {}
};

struct interval
{
    double R;
    point* lp;
    point* rp;
    interval(const double R_ = 0.0, point* lp_ = nullptr, point* rp_ = nullptr) : R(R_), lp(lp_), rp(rp_) {}
    interval(const interval& i)
    {
        R = i.R;
        lp = i.lp;
        rp = i.rp;
    }
    interval& operator=(const interval& i)
    {
        R = i.R;
        lp = i.lp;
        rp = i.rp;
        return *this;
    }
};

bool operator<(const interval& i1, const interval& i2)
{
    return (i1.R < i2.R) ? true : false;
}

double Rfunc(const point& lp_, const point& rp_, const double m)
{
    double dx = rp_.x - lp_.x;
    double dz = rp_.z - lp_.z;
    return (m * dx + dz * dz / (m * dx) - 2.0 * (rp_.z + lp_.z));
}

point* insertUpList(std::list<point>* p, point* xk)
{
    std::list<point>::iterator itl, itr;
    itl = itr = (*p).begin();
    while ((itr != (*p).end()) && (itr->x < (*xk).x))
    {
        itl = itr;
        itr++;
    }
    (*p).insert(itr, (*xk));
    itl++;
    return &(*itl);
}

double funcFindM(std::list<point>& listPoints, const double r)
{
    std::list<point>::iterator itl, itr;
    double mm = 0.0;
    double m;
    itr = itl = listPoints.begin();
    itr++;

    while (itr != listPoints.end())
    {
        double max = fabs((itr->z - itl->z) / (itr->x - itl->x));
        if (mm < max)
        {
            mm = max;
        }
        itr++;
        itl++;
    }

    if (mm > 0.0)
    {
        m = r * mm;
    }
    else
    {
        m = 1.0;
    }

    return m;
}

std::pair<double, double> funcFindMinInList(std::list<point>& listPoints)
{
    std::list<point>::iterator itl;
    itl = listPoints.begin();
    double minf = itl->z;
    double minx = itl->x;
    itl++;

    while (itl != listPoints.end())
    {
        if (minf > itl->z)
        {
            minf = itl->z;
            minx = itl->x;
        }
        itl++;
    }

    return std::pair<double, double>(minx, minf);
}

void refillingQueue(std::list<point>& listPoints, std::priority_queue<interval>& queueIntervals, const double m)
{
    std::list<point>::iterator itl, itr;
    while (!queueIntervals.empty())
    {
        queueIntervals.pop();
    }
    itr = itl = listPoints.begin();
    itr++;
    while (itr != listPoints.end())
    {
        queueIntervals.push(interval(Rfunc(*itl, *itr, m), &(*itl), &(*itr)));
        itl++;
        itr++;
    }
}

void findNewPoint(const Task& testTask, const std::vector<double>& transferInterval, std::vector<double>& transferResults)
{
    double xl = transferInterval[0];
    double zl = transferInterval[1];
    double xr = transferInterval[2];
    double zr = transferInterval[3];
    double mm = transferInterval[4];

    double xk = 0.5 * (xr + xl) - ((zr - zl) / (2.0 * mm));
    transferResults[0] = xk;
    double zk = testTask.data->func(xk);
    transferResults[1] = zk;
}