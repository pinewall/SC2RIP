#include "kinds.h"
#include "timers.h"
#include <iostream>
#include <string.h>
#include <sys/time.h>

using namespace std;

// decleare timers
Timers *timers = new Timers();

// default construction function
Timers::Timers()
{
    max_timers = 99;
    cycles1 = new double[max_timers];
    cycles2 = new double[max_timers];
    cputime = new double[max_timers];
    status = new char*[max_timers];
    for (int i = 0; i < max_timers; i++)
    {
        cycles1[i] = 0;
        cycles2[i] = 0;
        cputime[i] = 0.0;
        status[i] = new char[8];
        strcpy(status[i], "stopped");
    }
}

//self-define timer function, using gettimeofday
double Timers::system_time()
{
    struct timeval tv;
    double clockcycle = 0.0;
    gettimeofday(&tv, NULL);
    clockcycle = tv.tv_sec + tv.tv_usec * 1e-6;
    return clockcycle;
}

// member functions
void Timers::check(int timer)
{
    if (strcmp(status[timer], "running") == 0)
    {
        stop(timer);
        start(timer);
    }
}

void Timers::clear(int timer)
{
    cputime[timer] = 0.0;
}

double Timers::get(int timer)
{
    double timer_get = 0.0;
    if (strcmp(status[timer], "stopped") == 0)
    {
        timer_get = cputime[timer];
    }
    else
    {
        stop(timer);
        timer_get = cputime[timer];
        start(timer);
    }
    return timer_get;
}

void Timers::print(int timer)
{
    if (strcmp(status[timer], "stopped") == 0)
    {
        cout << "( CPU time for timer " << timer << "\t: " << cputime[timer] << " )" << endl;
    }
    else
    {
        stop(timer);
        cout << "( CPU time for timer " << timer << "\t: " << cputime[timer] << " )" << endl;
        start(timer);
    }
}

void Timers::start(int timer)
{
    if (strcmp(status[timer], "stopped") == 0)
    {
        cycles1[timer] = system_time();
        strcpy(status[timer], "running");
    }
}

void Timers::stop(int timer)
{
    if (strcmp(status[timer], "running") == 0)
    {
        cycles2[timer] = system_time();
        cputime[timer] += (cycles2[timer] - cycles1[timer]);
        strcpy(status[timer], "stopped");
    }
}

int Timers::get_max_timers()
{
    return max_timers;
}
