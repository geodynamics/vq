// Copyright (c) 2012-2014 Eric M. Heien, Michael K. Sachs, John B. Rundle
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "VCSimulation.h"

#ifndef _BG_POISSON_EVENTS_H_
#define _BG_POISSON_EVENTS_H_

/*!
 Creates background events with the given parameters at a Poisson rate during the simulation.
 */
class BGPoissonEvents : public SimPlugin {
    private:
        //! Background event parameters
        float _l;       //! Mean interevent time
        float _mag_min; //! Minimum event magnitude
        float _mag_max; //! Maximum event magnitude
        float _d;       //! Distance of events from main fault (meters)
        float _q;       //! Decay rate of events (distance)

        double last_year;   //! Last year of simulation when background events were generated

    public:
        virtual std::string name(void) const {
            return "Poisson background events";
        }
        virtual void initDesc(const SimFramework *_sim) const;
        virtual bool needsTimer(void) const {
            return true;
        };

        virtual void init(SimFramework *_sim);
        virtual SimRequest run(SimFramework *_sim);
};

#endif
