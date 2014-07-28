#ifndef GAUGE_H
#define GAUGE_H

#include "common.h"
#include <iostream>

namespace aw {

    struct gauge_exp {
        unsigned int counter, modifier;
    };

    inline void gauge_init(gauge_exp * d) {
        d->counter = 0;
        d->modifier = 1;
    }

    inline void gauge_inc(gauge_exp * d) {
        ++d->counter;
        if ((d->counter % d->modifier) == 0) {
            MSG_nonewline(d->modifier << "..");
            d->modifier *= 10;
        }
    }

    inline void gauge_end(gauge_exp * d) {
        MSG(d->counter);
        gauge_init(d);
    }

    // -----------------------------------------------------------------------------------

    struct gauge_percent {
        unsigned int counter, last, max, steps;
    };

    inline void gauge_init(gauge_percent * d, unsigned int _max, unsigned int _steps = 10) {
        d->counter = 0;
        d->last = 0;
        d->max = _max;
        d->steps = _steps;
    }

    inline void gauge_inc(gauge_percent * d) {
        if (d->counter == 0) MSG_nonewline("0%");
        ++d->counter;
        if (d->last < d->counter) {
            unsigned int v1 = (d->last * d->steps) / d->max;
            unsigned int v2 = (d->counter * d->steps) / d->max;
            if (v1 != v2) {
                MSG_nonewline(".." << (v2 * 100 / d->steps) << '%');
            }
            d->last = d->counter;
        }
    }

    inline void gauge_end(gauge_percent * d) {
        MSG("");
        gauge_init(d,d->max,d->steps);
    }
} // end of namespace

#endif // GAUGE_H
