#ifndef COMMON_H
#define COMMON_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
//#include "rmq.c"

// stream for normal messages
class StdOut {
    private: std::ofstream log_file;

    public: StdOut(const bool &logging_) : logging(logging_) {
        quiet = false;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (logging) log_file.open("data_log.txt");
    }

    public: inline StdOut& operator<<(std::ostream& (*r)(std::ostream&))
    {
        if (!quiet) std::cout << r;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (logging) log_file << r;
        return *this;
    }

    public: template<class T> inline StdOut& operator<<(const T &r) {
        if (!quiet) std::cout << r;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (logging) log_file << r;
        return *this;
    }
    
    public: const bool logging;
    public: bool quiet;
};

StdOut stdmsg(true);


// stream for error messages
class StdErr {
    private: std::ofstream log_file;
    private: bool log_open;

    private: bool open_log() {
        if (!logging) return false;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (!log_open) log_file.open("data_log_error.txt");
        log_open = true;
        return true;
    }

    public: StdErr(const bool &logging_) : logging(logging_) {
        log_open = false;
    }

    public: inline StdErr& operator<<(std::ostream& (*r)(std::ostream&))
{
        std::cerr << r;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (open_log()) log_file << r;
        return *this;
    }

    public: template<class T> inline StdErr& operator<<(const T &r)
{
        std::cerr << r;
        #ifdef WITH_MPI
        if (mpi_rank == 0)
        #endif
        if (open_log()) log_file << r;
        return *this;
    }

    public: const bool logging;

};

StdErr errmsg(true);

#ifdef WITH_MPI
#define ERROR_EXIT_COMMAND MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
#else
#define ERROR_EXIT_COMMAND exit(1);
#endif

#define ERROR(msg_arg) { \
    std::ostringstream os23; \
    os23 << "ERROR(" << __FILE__ << ':' << __LINE__ << "): " << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
}
#define ERROR_return(msg_arg) { \
    std::ostringstream os23; \
    os23 << "ERROR(" << __FILE__ << ':' << __LINE__ << "): " << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
    return false; \
}
#define ERROR_exit(msg_arg) { \
    std::ostringstream os23; \
    os23 << "ERROR(" << __FILE__ << ':' << __LINE__ << "): " << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
    ERROR_EXIT_COMMAND \
}

#define WARNING(msg_arg) { \
    std::ostringstream os23; \
    os23 << "WARNING(" << __FILE__ << ':' << __LINE__ << "): " << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
}
#define WARNING_return(msg_arg) { \
    std::ostringstream os23; \
    os23 << "WARNING(" << __FILE__ << ':' << __LINE__ << "): " << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
    return false; \
}

#define MSG_exit(msg_arg) { \
    std::ostringstream os23; \
    os23 << msg_arg << std::endl; \
    errmsg << os23.str() << std::flush; \
    ERROR_EXIT_COMMAND \
}
#define MSG(msg_arg) { \
    std::ostringstream os23; \
    os23 << msg_arg << std::endl; \
    stdmsg << os23.str() << std::flush; \
}
#define MSG_nonewline(msg_arg) { \
    std::ostringstream os23; \
    os23 << msg_arg; \
    stdmsg << os23.str() << std::flush; \
}
#define MSG_all(msg_arg) { \
    std::ostringstream os23; \
    os23 << msg_arg << std::endl; \
    stdmsg << os23.str() << std::flush; \
}
#define P(EX) { \
    std::ostringstream os23; \
    os23 << #EX << ": " << EX << std::endl; \
    stdmsg << os23.str() << std::flush; \
}

#endif // COMMON_H
