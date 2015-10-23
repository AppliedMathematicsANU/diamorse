#ifndef ANU_AM_PERFORMANCE_H
#define ANU_AM_PERFORMANCE_H

#include <fstream>
#include <sstream>
#include <string>
#include <sys/times.h>
#include <unistd.h>

class Stopwatch
{
private:
    bool const useCpuTime_;

    long accumulated_;
    long start_;
    bool isRunning_;

    tms mutable tmsCurrent_;

    long time() const
    {
        clock_t res = times(&tmsCurrent_);
        if (res < 0)
            throw "Cannot determine user time.";

    	if (useCpuTime_)
            return tmsCurrent_.tms_utime;
        else
            return res;
    }

public:
    Stopwatch(bool useCpuTime = true)
        : useCpuTime_(useCpuTime),
          accumulated_(0),
          start_(0),
          isRunning_(false)
    {
    }
    
    std::string mode() const
    {
    	if (useCpuTime_)
            return "CPU";
        else
            return "Real";
    }
    
    void resume()
    {
        if (!isRunning_)
        {
            isRunning_ = true;
            start_ = time();
        }
    }

    void start()
    {
        accumulated_ = 0;
        isRunning_ = true;
        start_ = time();
    }
    
    void stop()
    {
        if (isRunning_)
        {
            accumulated_ += time() - start_;
            isRunning_ = false;
        }
    }
    
    /**
     * Reports the elapsed time on this timer in milliseconds.
     */
    long elapsed() const
    {
        static long clktck = 0;

        if (clktck == 0 and (clktck = sysconf(_SC_CLK_TCK)) <= 0)
            throw "Cannot determine system clock rate";

    	return (accumulated_ + (isRunning_ ? time() - start_ : 0))
            * 1000 / clktck;
    }
    
    std::string format() const
    {
        return format(elapsed());
    }
    
    static std::string format(long const milliseconds)
    {
        std::stringstream ss;
    	ss << milliseconds / 10 / 100.0 << " seconds";
        return ss.str();
    }
};

struct memInfo
{
    size_t total;
    size_t resident;
    size_t heap;
    size_t stack;
};

#ifdef __linux__
#define ANU_AM_MEMORY_USAGE_DEFINED

memInfo memoryUsageInKB()
{
    std::string line;
    memInfo m;

    std::ifstream fp("/proc/self/status");

    while (not fp.eof())
    {
        std::getline(fp, line);
        if (line.find("VmSize:") == 0)
            std::stringstream(line.substr(7)) >> m.total;
        if (line.find("VmRSS:") == 0)
            std::stringstream(line.substr(6)) >> m.resident;
        if (line.find("VmData:") == 0)
            std::stringstream(line.substr(7)) >> m.heap;
        if (line.find("VmStk:") == 0)
            std::stringstream(line.substr(6)) >> m.stack;
    }

    fp.close();

    return m;
}

#endif //__linux__

//TODO implement for other operating systems

#endif // ANU_AM_PERFORMANCE_H
