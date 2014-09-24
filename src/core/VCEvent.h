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

#include "VCBlock.h"
#include "QuakeLibIO.h"
#include <map>
#include <set>

#ifndef _VCEVENT_H_
#define _VCEVENT_H_

typedef std::map<BlockID, int> BlockIDProcMapping;
typedef std::set<FaultID> FaultIDSet;
typedef std::set<SectionID> SectionIDSet;
typedef std::map<FaultID, quakelib::ElementIDSet> FaultBlockMapping;
typedef std::map<SectionID, quakelib::ElementIDSet> SectionBlockMapping;
typedef std::vector<FaultID> FaultIDList;
typedef std::vector<SectionID> SectionIDList;

/*!
 Generalized events are similar to aftershocks but don't have a generation.
 These represent background events in the simulation off of main faults.
 We use floats to save space and improve speed since there can be on the
 order of 1e6 or more general events per VCEvent.
 */
class VCGeneralEvent {
    public:
        //! Magnitude of earthquake
        float mag;
        //! Time of earthquake
        float t;
        //! "x"-position (km) of earthquake
        float x;
        //! "y"-position (km) of earthquake
        float y;
        //! "z"-position (km) of earthquake
        float z;

        VCGeneralEvent(float _m=0.0, float _t=0.0, float _x=0.0, float _y=0.0, float _z=0.0) : mag(_m), t(_t), x(_x), y(_y), z(_z) {};
        void clear(void) {
            mag = t = x = y = z = 0;
        };

        quakelib::Vec<3> loc(void) const {
            return quakelib::Vec<3>(x,y,z);
        };

        // Comparison operator (for sorting)
        bool operator < (const VCGeneralEvent &as) const {
            return (this->t < as.t);
        };
        friend std::ostream &operator<<(std::ostream &os, const VCGeneralEvent &e);
};

typedef std::vector<VCGeneralEvent> VCGeneralEventSet;

/*!
 Aftershocks are events caused by a rupture. Aftershocks have a generation to indicate
 how far removed they are from the original event (aftershocks caused by aftershocks and so on).
 */
class VCEventAftershock : public VCGeneralEvent {
    public:
        unsigned int gen;  // Generation number of earthquake (0=mainshock, N=aftershock generation N)

        // Constructor
        VCEventAftershock(float _m=0.0, float _t=0.0, float _x=0.0, float _y=0.0, unsigned int _g=0)
            : VCGeneralEvent(_m, _t, _x, _y), gen(_g) { };

        void clear(void) {
            VCGeneralEvent::clear();
            gen = 0;
        };
        friend std::ostream &operator<<(std::ostream &os, const VCEventAftershock &e);
};

typedef std::vector<VCEventAftershock> AftershockVector;
typedef std::set<VCEventAftershock> AftershockSet;

#endif
