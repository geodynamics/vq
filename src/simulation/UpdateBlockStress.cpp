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

#include "UpdateBlockStress.h"
#include "SimError.h"

/*!
 Initialize the stress calculation by setting initial block slip and stresses
 then calculating the stress in the whole system.
 */
void UpdateBlockStress::init(SimFramework *_sim) {
    BlockList::iterator it;
    BlockID             bid;

    sim = static_cast<VCSimulation *>(_sim);
    tmpBuffer = new double[sim->numGlobalBlocks()];

    for (it=sim->begin(); it!=sim->end(); ++it) {
        bid = it->getBlockID();

        // Set stresses to their specified initial values
        sim->setShearStress(bid, it->getInitShearStress());
        sim->setNormalStress(bid, it->getInitNormalStress());

        // noise in the range [1-a, 1+a)
        it->state.slipDeficit = 0;

        if (sim->isLocalBlockID(bid)) {
            sim->decompressNormalRow(bid);
            //sim->setGreenNormal(bid, bid, 0.0);  // Erase diagonal for normal Greens matrix
            sim->compressNormalRow(bid, 0.7);
        }
    }

    // Compute initial stress on all blocks
    stressRecompute();
}

/*!
 Calculate the stress on all blocks and determine which block will be
 the next to fail and when it will fail. Use this to determine the
 slipDeficit on all blocks at the time of failure. Finally, use the
 new slipDeficit to recalculate the stress on all blocks at the time of failure.
 */
SimRequest UpdateBlockStress::run(SimFramework *_sim) {
    // Put a stress load on all blocks and determine which block will fail first
    int                     lid;
    BlockVal                fail_time;
    quakelib::Conversion    convert;

    // Calculate the current rates of stress change on all blocks
    stressRecompute();

    // Given the rates of change, determine which block will fail next
    fail_time.block_id = UNDEFINED_BLOCK_ID;
    nextTimeStep(fail_time);

    // Increment the simulation year to the next failure time and
    // update the slip on all other blocks
    sim->incrementYear(fail_time.val);

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        Block &local_block = sim->getBlock(sim->getGlobalBID(lid));
        local_block.state.slipDeficit -= local_block.slip_rate()*convert.year2sec(fail_time.val)*(1.0-local_block.aseismic());
    }

    // Recalculate the stress on all blocks with the new slip deficits
    stressRecompute();

    if (sim->getYear() > sim->getSimDuration()) return SIM_STOP_REQUIRED;
    else return SIM_CONTINUE;
}

/*!
 Determine the next time step in the simulation when a failure occurs.
 Return the block ID of the block responsible for the failure and the timestep until the failure.
 */
void UpdateBlockStress::nextTimeStep(BlockVal &fail_time) {
    BlockList::iterator     it;
    BlockVal                temp_block_fail;
    double                  ts;
    VCEvent                 new_event;
    BlockID                 gid;
    int                     lid;
    quakelib::Conversion    convert;

    // Set up the temporary buffer and update field
    for (it=sim->begin(); it!=sim->end(); ++it) {
        tmpBuffer[it->getBlockID()] = 0.0;

        // Set the update field to be the slip rate of each block
        sim->setUpdateField(it->getBlockID(), it->slip_rate());
    }

    // update the temporary buffer with the Greens function applied to the block slip rates
    sim->matrixVectorMultiplyAccum(tmpBuffer,
                                   sim->greenShear(),
                                   sim->getUpdateFieldPtr(),
                                   true);

    if (sim->doNormalStress()) {
        for (it=sim->begin(); it!=sim->end(); ++it) {
            sim->setUpdateField(it->getBlockID(), -it->friction()*it->slip_rate());
        }

        sim->matrixVectorMultiplyAccum(tmpBuffer,
                                       sim->greenNormal(),
                                       sim->getUpdateFieldPtr(),
                                       true);
    }

    // Go through the blocks and find which one will fail first
    temp_block_fail.val = DBL_MAX;
    temp_block_fail.block_id = UNDEFINED_BLOCK_ID;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);
        Block &block = sim->getBlock(gid);

        // Calculate the time until this block will fail
        // If the block has aseismic slip, the calculation is the exact solution of the
        // differential equation d(cff)/dt = rate_of_stress_change + cff*aseismic_frac*self_shear/recurrence
        if (block.aseismic() > 0) {
            double      A, B, K;
            A = -tmpBuffer[gid];
            B = -block.aseismic()*block.getSelfStresses()/block.getRecurrence();
            K = -log(A+B*block.getCFF())/B;
            ts = K + log(A)/B;
        } else {
            ts = convert.sec2year(block.getCFF()/tmpBuffer[gid]);
        }

        // Blocks with negative timesteps are skipped. These effectively mean the block
        // should have failed already but didn't. This can happen in the starting phase
        // of a simulation due to initial stresses. These blocks will fail anyway in the
        // event propagation section, so we can safely ignore them here. However, negative
        // timesteps should not happen after the simulation has progressed for a while.
        if (ts <= 0) continue;

        // If the time to slip is less than the current shortest time, record the block
        if (ts < temp_block_fail.val) {
            temp_block_fail.block_id = gid;
            temp_block_fail.val = ts;
        }
    }

    // Each node now has the time before the first failure among its blocks
    // Determine the time to first failure over all nodes
    sim->allReduceBlockVal(temp_block_fail, fail_time, BLOCK_VAL_MIN);

    // If we didn't find any blocks that slipped, abort the simulation
    assertThrow(fail_time.val < DBL_MAX, "System stuck, no blocks to move.");

    // Setup the current event to have the trigger block ID and time
    new_event.setEventTriggerOnThisNode(fail_time.block_id==temp_block_fail.block_id);
    new_event.setEventTrigger(fail_time.block_id);
    new_event.setEventYear(sim->getYear()+fail_time.val);
    new_event.setEventNumber(sim->getEventCount());
    sim->addEvent(new_event);
}

/*!
 Recompute the stress on each block based on the slip deficits of all the other blocks.
 */
void UpdateBlockStress::stressRecompute(void) {
    BlockList::iterator it;

    // Reset the shear and normal stress, and set the update field to be the slip deficit
    for (it=sim->begin(); it!=sim->end(); ++it) {
        sim->setShearStress(it->getBlockID(), 0.0);
        sim->setNormalStress(it->getBlockID(), it->getRhogd());
        sim->setUpdateField(it->getBlockID(), it->state.slipDeficit);
    }

    // Distribute the new update field over all nodes
    sim->distributeUpdateField();

    // Multiply the Greens shear function by the slipDeficit vector
    // to get shear on local blocks
    sim->matrixVectorMultiplyAccum(sim->getShearStressPtr(),
                                   sim->greenShear(),
                                   sim->getUpdateFieldPtr(),
                                   true);

    if (sim->doNormalStress()) {
        sim->matrixVectorMultiplyAccum(sim->getNormalStressPtr(),
                                       sim->greenNormal(),
                                       sim->getUpdateFieldPtr(),
                                       true);
    }

    // Recompute the CFF on blocks based on the new shear/normal stresses
    sim->computeCFFs();

}

void UpdateBlockStress::finish(SimFramework *_sim) {
    delete tmpBuffer;
}
