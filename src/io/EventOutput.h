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
#include "HDF5Data.h"

#ifdef VC_HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifndef _HDF5_DATA_SHARE_H_
#define _HDF5_DATA_SHARE_H_

#define PAUSE_FILE_NAME     "pause_vc"

char **event_fields = { "# Event number", "# Event year", "# Event trigger element ID", "# Event magnitude", "# Total initial shear stress of elements involved in event (Pascals)", "# Total final shear stress of elements involved in event (Pascals)", "# Total initial normal stress of elements involved in event (Pascals)", "# Total final normal stress of elements involved in event (Pascals)"};

char **sweep_fields = { "# Event number corresponding to this sweep", "# Sweep number", "# Element ID", "# Slip on element in this sweep (meters)", "# Area of element (square meters)", "# Element Lame mu parameter (Pascals)", "# Shear stress of element before sweep (Pascals)", "# Shear stress of element after sweep (Pascals)", "# Normal stress of element before sweep (Pascals)", "# Normal stress of element after sweep (Pascals)" };

/*!
 Manages the HDF5 format data output for VC program.
 */
class EventOutput : public SimPlugin {
    private:
#ifdef HDF5_FOUND
        HDF5DataWriter      *h5_data;
#else
        void                *h5_data;
#endif
        unsigned int        next_pause_check;
        std::ofstream       event_outfile, sweep_outfile;

        bool pauseFileExists(void);
        void writeEventFileHeader(void);
        void writeSweepFileHeader(void);

    public:
        virtual std::string name(void) const {
            return "Data Writer";
        }
        void initDesc(const SimFramework *_sim) const;

        virtual bool needsTimer(void) const {
            return true;
        };
        virtual void init(SimFramework *_sim);
        virtual SimRequest run(SimFramework *_sim);
        virtual void finish(SimFramework *_sim);
};

#endif
