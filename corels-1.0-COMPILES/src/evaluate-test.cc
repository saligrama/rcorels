#include "evaluate.hh"

int main(int argc, char ** argv)
{
    logger = new NullLogger();

    int r = run_random_tests(1000, 30, 701, 0.15, 0, 2, curious_cmp, "curious", false,
                             100000, 0.0000001, time(NULL), 3);

    delete logger;

    return r;
}
