#include <queue>
#include <cstdlib>
#include <chrono>
#include <omp.h>
#include "Problem.h"
#include "SupportStructures.h"

int AGSLinearVersion(int argc, char* argv[])
{
    Task testTask;

    int maxNumIter = 2000;
    double eps = 0.0000001;
    int end = 1;

    std::list<point> listPoints;
    std::priority_queue<interval> queueIntervals;
    double x;

    std::vector<double> transferInterval(5);
    std::vector<double> transferResult(2);

    double m = -1.0;
    int numIter = 0;

    auto st = std::chrono::steady_clock::now();

    x = testTask.data->a;
    listPoints.push_back(point(x, testTask.data->func(x)));
    x = testTask.data->b;
    listPoints.push_back(point(x, testTask.data->func(x)));

    do
    {
        double mold = m;
        m = funcFindM(listPoints, testTask.data->r);
        if (mold != m)
        {
            refillingQueue(listPoints, queueIntervals, m);
        }

        interval priorityInterval = queueIntervals.top();
        queueIntervals.pop();

        transferInterval[0] = priorityInterval.lp->x;
        transferInterval[1] = priorityInterval.lp->z;
        transferInterval[2] = priorityInterval.rp->x;
        transferInterval[3] = priorityInterval.rp->z;
        transferInterval[4] = m;

        findNewPoint(testTask, transferInterval, transferResult);

        double xk = transferResult[0];
        double zk = transferResult[1];

        point pointK(xk, zk);
        point* indPointK = insertUpList(&listPoints, &pointK);
        queueIntervals.push(interval(Rfunc(*priorityInterval.lp, *indPointK, m), priorityInterval.lp, indPointK));
        queueIntervals.push(interval(Rfunc(*indPointK, *priorityInterval.rp, m), indPointK, priorityInterval.rp));

        numIter++;

        end = 1;
        if (priorityInterval.rp->x - priorityInterval.lp->x <= eps)
        {
            end = 0;
        }

        if (numIter >= maxNumIter)
        {
            end = 0;
        }
    } while (end != 0);

    std::pair<double, double> optimum = funcFindMinInList(listPoints);

    auto fi = std::chrono::steady_clock::now();

    std::cout << "Arg min f = " << optimum.first << std::endl;
    std::cout << "Min f = " << optimum.second << std::endl;
    std::cout << "Number iterations = " << numIter << std::endl;
    std::cout << "Linear time work = " << std::chrono::duration_cast<std::chrono::milliseconds>(fi - st).count() / 1000.0 << std::endl;
    std::cout << std::endl;

    testTask.data->getInfo();
    return 0;
}

int AGSParallelVersion(int argc, char* argv[])
{
    Task testTask;

    int maxNumIter = 2000;
    double eps = 0.0000001;
    int numThreads = 4;
    int end = 1;

    std::list<point> listPoints;
    std::priority_queue<interval> queueIntervals;
    double x;

    std::vector<std::vector<double>> transferIntervals(numThreads, std::vector<double>(5));
    std::vector<std::vector<double>> transferResults(numThreads, std::vector<double>(2));
    std::vector<interval> priorityIntervals(numThreads);

    double m = -1.0;
    int numIter = 0;

    auto st = std::chrono::steady_clock::now();

    double pr = (testTask.data->b - testTask.data->a) / numThreads;
    for (int i = 0; i < numThreads; i++) {
        x = testTask.data->a + pr * i;
        listPoints.push_back(point(x, testTask.data->func(x)));
    }
    x = testTask.data->b;
    listPoints.push_back(point(x, testTask.data->func(x)));

    do
    {
        double mold = m;
        m = funcFindM(listPoints, testTask.data->r);
        if (mold != m)
        {
            refillingQueue(listPoints, queueIntervals, m);
        }

        for (int numThread = 0; numThread < numThreads; numThread++)
        {
            priorityIntervals[numThread] = queueIntervals.top();
            queueIntervals.pop();

            transferIntervals[numThread][0] = priorityIntervals[numThread].lp->x;
            transferIntervals[numThread][1] = priorityIntervals[numThread].lp->z;
            transferIntervals[numThread][2] = priorityIntervals[numThread].rp->x;
            transferIntervals[numThread][3] = priorityIntervals[numThread].rp->z;
            transferIntervals[numThread][4] = m;
        }

        #pragma omp parallel num_threads(numThreads)
        {
            int numThread = omp_get_thread_num();
            findNewPoint(testTask, transferIntervals[numThread], transferResults[numThread]);
        }

        for (int numThread = 0; numThread < numThreads; numThread++)
        {
            double xk = transferResults[numThread][0];
            double zk = transferResults[numThread][1];

            point pointK(xk, zk);
            point* indPointK = insertUpList(&listPoints, &pointK);
            queueIntervals.push(interval(Rfunc(*priorityIntervals[numThread].lp, *indPointK, m), priorityIntervals[numThread].lp, indPointK));
            queueIntervals.push(interval(Rfunc(*indPointK, *priorityIntervals[numThread].rp, m), indPointK, priorityIntervals[numThread].rp));

            numIter++;
        }

        end = 1;
        for (int numThread = 0; numThread < numThreads; numThread++)
        {
            if (priorityIntervals[numThread].rp->x - priorityIntervals[numThread].lp->x <= eps)
            {
                end = 0;
            }
        }

        if (numIter >= maxNumIter)
        {
            end = 0;
        }
    } while (end != 0);

    std::pair<double, double> optimum = funcFindMinInList(listPoints);

    auto fi = std::chrono::steady_clock::now();

    std::cout << "Arg min f = " << optimum.first << std::endl;
    std::cout << "Min f = " << optimum.second << std::endl;
    std::cout << "Number iterations = " << numIter << std::endl;
    std::cout << "Parallel time work = " << std::chrono::duration_cast<std::chrono::milliseconds>(fi - st).count() / 1000.0 << std::endl;
    std::cout << std::endl;

    testTask.data->getInfo();
    return 0;
}

int main(int argc, char* argv[])
{
    return AGSParallelVersion(argc, argv);
}