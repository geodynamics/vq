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

#include "VCCommSpecExec.h"

/*!
 Notify the simulation that an event is starting (for use in speculative execution).
 */
void VCCommSpecExec::startEvent(void) {
	last_prediction_failed = false;
}

/*!
 Notify the simulation that the speculation failed and the event is being rewound.
 */
void VCCommSpecExec::speculationFailed(VCEventSweep &fail_sweep) {
	EventBlockMap::const_iterator		it;
	
	last_prediction_failed = true;
	num_predictions_failed++;
	
	if (spec_exec_method == SPEC_EXEC_ADAPTIVE) {
		for (it=fail_sweep.begin();it!=fail_sweep.end();++it) {
			boundary_dist_map[it->first] *= 0.5;
		}
	}
}

/*!
 Notify the simulation that the speculation failed and the event is being rewound.
 */
void VCCommSpecExec::speculationSuccess(VCEventSweep &success_sweep) {
	EventBlockMap::const_iterator		it;
	
	num_predictions_success++;
	
	if (spec_exec_method == SPEC_EXEC_ADAPTIVE) {
		for (it=success_sweep.begin();it!=success_sweep.end();++it) {
			boundary_dist_map[it->first] *= 1.05;
		}
	}
}
