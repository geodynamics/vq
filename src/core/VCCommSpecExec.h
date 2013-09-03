// Copyright (c) 2010, John B. Rundle <rundle@cse.ucdavis.edu>, 
// All rights reserved.
// 
// Redistribution and use of this code or any derivative works are
// permitted provided that the following conditions are met:
// 
// * Redistributions may not be sold, nor may they be used in a
// commercial product or activity.
// 
// * Redistributions that are modified from the original source must
// include the complete source code, including the source code for all
// components used by a binary built from the modified
// sources. However, as a special exception, the source code
// distributed need not include anything that is normally distributed
// (in either source or binary form) with the major components
// (compiler, kernel, and so on) of the operating system on which the
// executable runs, unless that component itself accompanies the
// executable.
// 
// * Redistributions must reproduce the above copyright notice, this list
// of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
