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

#include "RunEvent.h"


/*!
 At the end of each sweep after we have recalculated block CFF, we determine
 which blocks will have a failure due to dynamic or static stress changes.
 */
void RunEvent::markBlocks2Fail(Simulation *sim, const FaultID &trigger_fault) {
    int         lid;
    BlockID     gid;
    bool        add;

    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        gid = sim->getGlobalBID(lid);

        // Add this block if it has a static or dynamic CFF failure
        add = sim->cffFailure(gid) || sim->dynamicFailure(gid, trigger_fault);
        
        // Sachs allowed "additional failure"?
//        bool Block::additionalFailure(void) const {
//        return (state.stressS[0] > 0.0 && (state.cff-state.Fcff)/state.cff > dynamic_val);
//        }

		// Let each block fail at most 10 times (not per sweep but per event). This is somewhat arbitrary, but it prevents runaway ruptures.
        if(add && num_failures[gid] < 10) {
            num_failures[gid] += 1;
            sim->setFailed(gid, true);
            local_failed_elements.insert(gid);
        }
    }
}

/*!
 Process the list of blocks that failed on this node for this sweep.
 */
void RunEvent::processFailedBlocks(Simulation *sim, quakelib::ModelSweeps &sweeps) {
    quakelib::ElementIDSet::iterator    fit;
    double                              slip, stress_drop;

    // For each block that fails in this sweep, calculate how much it slips
    for (fit=local_failed_elements.begin(); fit!=local_failed_elements.end(); ++fit) {
        if (sim->isLocalBlockID(*fit)) {
            BlockID gid = *fit;
            Block &b = sim->getBlock(*fit);
            //
            // calculate the drop in stress from the failure

            // Heien method: 
            //stress_drop = sim->getCFF0(gid) - sim->getCFF(gid);
            // If this is the initial failure, use the stress drop
            //if (!stress_drop) stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);

            // Sachs method:
            stress_drop = sim->getStressDrop(gid) - sim->getCFF(gid);

            // Slip is in m
            slip = (stress_drop/sim->getSelfStresses(gid));

            if (slip < 0) slip = 0;

            // Record how much the block slipped in this sweep and initial stresses
            sweeps.setSlipAndArea(sweep_num,
                                  b.getBlockID(),
                                  slip,
                                  b.area(),
                                  b.lame_mu());
            sweeps.setInitStresses(sweep_num,
                                   b.getBlockID(),
                                   sim->getShearStress(gid),
                                   sim->getNormalStress(gid));

            sim->setSlipDeficit(gid, sim->getSlipDeficit(gid)+slip);
        }
    }
}


/*!
 Process the next aftershock. This involves determining a suitable rupture area from an empirical
 relationship, finding the nearest elements, choosing enough elements to match the empirical
 relationship, calculating the slip needed to generate the aftershock, and updating the stress
 field appropriately.
 */
void RunEvent::processAftershock(Simulation *sim) {
    std::map<double, BlockID>                   as_elem_dists;
    std::map<double, BlockID>::const_iterator   it;
    std::map<BlockID, double>                   elem_slips;
    EventAftershock                             as;
    BlockID                                     gid;
    quakelib::ElementIDSet                      id_set;
    quakelib::ElementIDSet::const_iterator      bit;
    quakelib::Conversion                        convert;
    quakelib::ModelSweeps                       event_sweeps;

    // Set the update field to the slip on each block
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        sim->setShearStress(gid, 0.0);
        sim->setNormalStress(gid, sim->getRhogd(gid));
        sim->setUpdateField(gid, sim->getSlipDeficit(gid));
    }

    // And distribute this around
    sim->distributeUpdateField();

    // Only process the aftershock stress effects on the root node
    if (sim->isRootNode()) {
        // Pop the next aftershock off the list
        as = sim->popAftershock();

        // Calculate the distance from the aftershock to all elements
        for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
            double as_to_elem_dist = sim->getBlock(gid).center().dist(as.loc());
            as_elem_dists.insert(std::make_pair(as_to_elem_dist, gid));
        }

        // Determine the target rupture area given the aftershock magnitude
        // TODO:
        // user_defined_constants (flag this for later revisions in which we move these contant definitions to a parameters file).
        double rupture_area = pow(10, -3.49+0.91*as.mag);
        double selected_rupture_area = 0;
        double selected_rupture_area_mu = 0;

        // Go through the elements, closest first, until we find enough to match the rupture area
        for (it=as_elem_dists.begin(); it!=as_elem_dists.end(); ++it) {
            Block &b=sim->getBlock(it->second);
            selected_rupture_area += convert.sqm2sqkm(b.area());
            selected_rupture_area_mu += b.area()*b.lame_mu();
            id_set.insert(it->second);

            if (selected_rupture_area > rupture_area) break;
        }

        // Determine the amount of slip needed to match the aftershock magnitude
        // The contribution of each block to moment is based on its fraction of total area*mu
        double total_moment = pow(10, (as.mag + 10.7)*(3.0/2.0))/1e7;

        for (bit=id_set.begin(); bit!=id_set.end(); ++bit) {
            Block &b=sim->getBlock(*bit);

            // Calculate the slip based on the earthquake moment
            double element_moment = total_moment*(b.area()*b.lame_mu()/selected_rupture_area_mu);
            double element_slip = element_moment/(b.lame_mu()*b.area());

            // Adjust the slip deficit appropriately
            sim->setUpdateField(*bit, sim->getUpdateField(*bit)+element_slip);

            // Create the sweep describing this aftershock
            // Since we don't distinguish sweeps, every slip occurs in sweep 0
            event_sweeps.setSlipAndArea(0, *bit, element_slip, b.area(), b.lame_mu());
            event_sweeps.setInitStresses(0, *bit, sim->getShearStress(*bit), sim->getNormalStress(*bit));
        }

    }

    // Broadcast the new slip deficit from the root node to all processes
    sim->broadcastUpdateField();

    // And update the slip deficit on each process to take this into account
    for (gid=0; gid<sim->numGlobalBlocks(); ++gid) {
        sim->setSlipDeficit(gid, sim->getUpdateField(gid));
    }

    // Calculate the new shear stresses and CFFs given the new update field values
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

    sim->computeCFFs();

    // Record final stresses on each block involved in the aftershock
    for (bit=id_set.begin(); bit!=id_set.end(); ++bit) {
        event_sweeps.setFinalStresses(0, *bit, sim->getShearStress(*bit), sim->getNormalStress(*bit));
    }

    // Add sweeps to list
    sim->getCurrentEvent().setSweeps(event_sweeps);
}


/*!
 Given an initial failed block, propagates the failure throughout the system
 by calculating changes in stress and using static and dynamic stress
 failure functions. A single step in the failure propagation is called a sweep
 and multiple sweeps comprise an entire event.
 */
void RunEvent::processStaticFailure(Simulation *sim) {
    BlockList::iterator     it;
    quakelib::ModelSweeps   event_sweeps;
    //BlockID                 triggerID, gid;
    BlockID                 triggerID;          // limit variable gid to local loop scopes.
    FaultID                 trigger_fault;
    int                     more_blocks_to_fail;
    bool                    final_sweep = false;
    std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> nbr_start_end;
    quakelib::ElementIDSet::const_iterator  nit;
    // Get the event trigger block and fault
    triggerID = sim->getCurrentEvent().getEventTrigger();
    trigger_fault = sim->getBlock(triggerID).getFaultID();
    sweep_num = 0;
    
    // schultz: Sachs' uses speculative execution (SpecExec.cpp) to guess whether the rupture will be local to one processor or distributed.
    // Here I am assuming that Heien's global/local scoping handles this.
    
    // Clear the list of failed blocks, and add the trigger block
    local_failed_elements.clear();
    num_failures.clear();
    
    // Clear the list of "loose" (can dynamically fail) blocks
    loose_elements.clear();

    if (sim->getCurrentEvent().getEventTriggerOnThisNode()) {
        local_failed_elements.insert(triggerID);
        sim->setFailed(triggerID, true);
        num_failures[triggerID] += 1;
    }
    
    // yoder (note): Comm::blocksToFail() executes a single MPI_Allreduce() (when MPI is present),
    // ... but since it's only the one "all-hands" call, it is not a likely candidate for heisen_hang.
    // also note that blocksToFail() involves some bool/int conversions that might warrant further inspection and clarification.
    more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());

    //
    // While there are still failed blocks to handle
    //sim->barrier();
    while (more_blocks_to_fail || final_sweep) {
        // write stress, slip, etc. to events and sweeps output (text or hdf5).
        sim->output_stress(sim->getCurrentEvent().getEventNumber(), sweep_num);

        // Share the failed blocks with other processors to correctly handle
        // faults that are split among different processors
        //sim->barrier();    // yoder: (debug)   (we're probably safe without this barrier() )... but at some point, i was able to generate a hang during distributeBlocks()
        // so let's try it with this in place...
        sim->distributeBlocks(local_failed_elements, global_failed_elements);
        //
        // Process the blocks that failed.
        // note: setInitStresses() called in processBlocksOrigFail().
        // note: processBlocksOrigFail() is entirely local (no MPI).
        processFailedBlocks(sim, event_sweeps);

        // Recalculate CFF for all blocks where slipped blocks don't contribute
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();
            sim->setShearStress(gid, 0.0);
            sim->setNormalStress(gid, sim->getRhogd(gid));
            sim->setUpdateField(gid, (sim->getFailed(gid) ? 0 : std::isnan(sim->getSlipDeficit(gid)) ? 0 :sim->getSlipDeficit(gid) )); // ... also check for nan values
            //sim->setUpdateField(gid, (sim->getFailed(gid) ? 0 : sim->getSlipDeficit(gid)));
        }

        // Distribute the update field values to other processors
        // (possible) MPI operations:
        //sim->barrier();    // yoder: (debug)
        sim->distributeUpdateField();
        //sim->barrier();    // yoder: (debug)


        // Set dynamic triggering on for any blocks neighboring blocks that slipped in the last sweep
        for (it=sim->begin(); it!=sim->end(); ++it) {
            BlockID gid = it->getBlockID();

            // Add block neighbors if the block has slipped
            if (sim->getUpdateField(gid) > 0) {
                nbr_start_end = sim->getNeighbors(gid);

                for (nit=nbr_start_end.first; nit!=nbr_start_end.second; ++nit) {
                    loose_elements.insert(*nit);
                }
            }
        }
    

        //
        // Calculate the CFFs based on the stuck blocks
        // multiply greenSchear() x getUpdateFieldPtr() --> getShearStressPtr() ... right?
        // assign stress values (shear stresses at this stage are all set to 0; normal stresses are set to (i think) sim->getRhogd(gid) -- see code a couple paragraphs above.
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

        sim->computeCFFs();

        // schultz: Skipping the secondary ruptures

        
        // use this iterator/counter to efficiently walk through the event_sweeps list when we update stresses:
        unsigned int event_sweeps_pos = 0;

        for (quakelib::ModelSweeps::iterator s_it=event_sweeps.begin(); s_it!=event_sweeps.end(); ++s_it, ++event_sweeps_pos) {
            //
            // yoder: as per request by KS, change isnan() --> std::isnan(); isnan() appears to throw an error on some platforms.
            if (std::isnan(s_it->_shear_final) and std::isnan(s_it->_normal_final)) {
                // note: the stress entries are initialized with nan values, but if there are cases where non nan values need to be updated,
                // this logic should be revisited.
                event_sweeps.setFinalStresses(sweep_num,
                                              s_it->_element_id,
                                              sim->getShearStress(s_it->_element_id),
                                              sim->getNormalStress(s_it->_element_id));

            }
        }

        //
        global_failed_elements.clear(); // we are done with these blocks
        local_failed_elements.clear();  // we are done with these blocks
        //
        // Find any blocks that fail because of the new stresses (all local; no MPI).
        markBlocks2Fail(sim, trigger_fault);

        if (final_sweep) {
            final_sweep = false;
        } else {
            more_blocks_to_fail = sim->blocksToFail(!local_failed_elements.empty());

            if (!more_blocks_to_fail) final_sweep = true;
        }

        sweep_num++;
    }

    //
    // output_stress() for final item in list.
    sim->output_stress(sim->getCurrentEvent().getEventNumber(), sweep_num);

    // Set the completed list as the sweep list for the entire event
    sim->collectEventSweep(event_sweeps);
    sim->getCurrentEvent().setSweeps(event_sweeps);

}


// TODO: Add rupture model version as a parameter



SimRequest RunEvent::run(SimFramework *_sim) {
    Simulation            *sim = static_cast<Simulation *>(_sim);
    int                     lid;
    BlockList::iterator     it;
    quakelib::ModelSweeps   event_sweeps;
    std::pair<quakelib::ElementIDSet::const_iterator, quakelib::ElementIDSet::const_iterator> nbr_start_end;
    quakelib::ElementIDSet::const_iterator  nit;

    // Save stress information at the beginning of the event
    // This is used to determine dynamic block failure
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) sim->saveStresses(sim->getGlobalBID(lid));

    // If there's a specific block that triggered the event, it's a static stress failure type event
    if (sim->getCurrentEvent().getEventTrigger() != UNDEFINED_ELEMENT_ID) {
        processStaticFailure(sim);
    } else {
        // Otherwise it's an aftershock
        processAftershock(sim);
    }

    // Record the stress in the system before and after the event.
    // yoder: note that recordEventStresses() emloyes a bit of MPI action, so it might be advisable to
    // add some barrier() blocking (which might have been added above to barrier()-wrap the process_{earthquake type}() call).
    recordEventStresses(sim);

    // Update the cumulative slip for this fault
    for (lid=0; lid<sim->numLocalBlocks(); ++lid) {
        BlockID gid = sim->getGlobalBID(lid);
        sim->setFailed(gid, false);
    }

    // TODO: reinstate this check
    // TODO: currently fails because processors may locally have no failures
    // TODO: notes to self(s) then: this single line works in SPP mode, or more specifically for a single node, so we can use an MPI_reduce() call to get the max or sum
    // ... and could this be causing heisen_hang? would this create a scenario where a child node would send/receive block data to root
    // but the root node would not send back anything (block not in list of failed blocks), so that child node would just hang there?
    // maybe, but i think that by this time, we'd be long since past that point.
    //       of all local sim->getCurrentEvent().size() values.
    //
    //assertThrow(sim->getCurrentEvent().size() > 0, "There was a trigger but no failed blocks.");
    //
#ifdef MPI_C_FOUND
    // so let's get the total event size by doing an MPI sum. note this can be an MPI_Reduce() to the root node only, or we can count on all processors.
    // this should be an equivalent operation (though slower in the latter case); i guess then, if we gather only to the root noode, we do the assert/throw
    // only on root node.:
    int local_event_size = sim->getCurrentEvent().size();
    int global_event_size = 0;
    //
    // aggregate and assert on root node:
    MPI_Reduce(&local_event_size, &global_event_size, 1, MPI_INT, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);

    if (sim->isRootNode()) {
        assertThrow(global_event_size > 0, "There was a trigger but no failed blocks.");
    };

    //
    //// aggregate and assert on all nodes:
    //MPI_Allreduce(&local_event_size, &global_event_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //assertThrow(sim->getCurrentEvent().size() > 0, "There was a trigger but no failed blocks.");
#else
    //global_event_size=local_event_size;
    assertThrow(global_event_size > 0, "There was a trigger but no failed blocks. (" << getpid() << "/" << sim->getNodeRank() << ")");

#endif
    return SIM_STOP_OK;
}

void RunEvent::recordEventStresses(Simulation *sim) {
    quakelib::ElementIDSet involved_blocks;
    double shear_init, shear_final, normal_init, normal_final;
    double total_shear_init, total_shear_final, total_normal_init, total_normal_final;

    involved_blocks = sim->getCurrentEvent().getInvolvedElements();

    sim->getInitialFinalStresses(involved_blocks, shear_init, shear_final, normal_init, normal_final);

    // If we have multiple processors, sum the values and store on the root node
#ifdef MPI_C_FOUND
    MPI_Reduce(&shear_init, &total_shear_init, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&shear_final, &total_shear_final, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&normal_init, &total_normal_init, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
    MPI_Reduce(&normal_final, &total_normal_final, 1, MPI_DOUBLE, MPI_SUM, ROOT_NODE_RANK, MPI_COMM_WORLD);
#else
    total_shear_init = shear_init;
    total_shear_final = shear_final;
    total_normal_init = normal_init;
    total_normal_final = normal_final;
#endif

    sim->getCurrentEvent().setEventStresses(total_shear_init, total_shear_final, total_normal_init, total_normal_final);
}
