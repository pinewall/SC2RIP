#ifndef _TIMERS_H
#define _TIMERS_H 1
class Timers
{
    private:
    /* Variables */
    int max_timers;     // max number of timers allowed
    int cycles_max;     // max value of clock allowed by system
    double * cycles1;      // cycle number at start for each timer
    double * cycles2;      // cycle number at stop for each timer
    double clock_rate;  // clock_rate in seconds for each cycle
    double * cputime;   // accumulated cpu time in each timer
    char **status;      // timer status string
    /* char (*status)[8];*/

    public:
    /* functions */
    Timers();
    double system_time();
    void check(int timer);
    void clear(int timer);
    double get(int timer);
    void print(int timer);
    void start(int timer);
    void stop(int timer);

    int get_max_timers();
};

extern Timers *timers;          // tic-tac clock

#endif
