// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cassert>
#include <cmath>
#include <numeric>

#include "thrust_common.h"

namespace edda {

/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation

                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
              this software for any application provided this
              copyright notice is preserved.

*/
template <class T>
T box_muller(T m, T s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
    T x1, x2, w, y1;
    static T y2;
    static char use_last = 0;

    if (use_last)		        /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * (T)rand()/RAND_MAX - 1.0;
            x2 = 2.0 * (T)rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );
}

}  // namespace edda

#endif  // STATISTICS_H_
