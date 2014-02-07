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

#include "VCParams.h"
#include "VCBlock.h"
#include "VCEvent.h"

#ifndef _VC_COMM_SPEC_EXEC_H_
#define _VC_COMM_SPEC_EXEC_H_

class VCCommSpecExec {
protected:
	SpecExecMethod		spec_exec_method;
	
	//! Map of blocks to boundary of node (for speculative execution)
	std::map<BlockID, double>	boundary_dist_map;
	
	bool						last_prediction_failed;
	
	//! Statistics for speculative execution
	int							num_predictions;
	int							num_predicted_local;
	int							num_predictions_failed;
	int							num_predictions_success;
	
public:
	VCCommSpecExec(void) : spec_exec_method(SPEC_EXEC_UNDEFINED), num_predictions(0), num_predicted_local(0), num_predictions_failed(0), num_predictions_success(0) {};
	void setSpecExecMethod(const SpecExecMethod &new_method) { spec_exec_method = new_method; };
	
	void startEvent(void);
	void speculationFailed(VCEventSweep &fail_sweep);
	void speculationSuccess(VCEventSweep &success_sweep);
};

#endif
