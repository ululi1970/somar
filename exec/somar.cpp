/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/


// This is the driver of the Navier-Stokes solver. It contains the main() C++
// function. This file was not written to be user-friendly and should not need
// to be altered. -ES


// Standard headers
#include <iostream>
#include <sys/utsname.h>
using std::cout;
using std::endl;

#ifdef CH_MPI
#   include "mpi.h"
#endif

// Chombo headers
#include "SPMD.H"
#include "SPACE.H"
#include "parstream.H"
#include "ParmParse.H"
#include "memusage.H"

#ifndef NDEBUG
#   include "CH_Attach.H"
#endif

#ifdef CH_USE_MEMORY_TRACKING
#   include "Pool.H"
#endif

#ifndef CH_NTIMER
#   include "OldTimer.H"
#endif

// Project headers
#include "AMRLESMeta.H"
#include "Printing.H"
#include "LevelGeometry.H"
#include "AMRNavierStokesFactory.H"
#include "LepticAMR.H"
#include "FORT_PROTO.H"
#include "ProblemContext.H"


// Function prototypes
void setupMPIComms (const Real a_lesProcFrac);
void testMPIComms ();
void nsrun ();

// The LES code is not for public use.
// #define USE_LES
#ifdef USE_LES
extern "C" {
    // void FORTRAN_NAME(DIABLO,diablo) (
    //     int*, int*, int*, int*
    // );
    void FORTRAN_NAME(EDDY,eddy) (
        int*, int*, int*, int*
    );
}
#endif //USE_LES


// #define TRAP_FPE  //(should be off by default)
#ifdef TRAP_FPE
    // Previous versions of glibc require the following code:
#   include "parstream.H"
    extern "C"
    {
#       include <fpu_control.h>
    }

    /* IM: Invalid operation mask
     * DM: Denormalized operand mask
     * ZM: Zero-divide mask
     * OM: Overflow mask
     * UM: Underflow mask
     * PM: Precision (inexact result) mask */
    static void __attribute__ ((constructor)) trapfpe(void)
    {
        pout() << " Turning on floating-point traps! " << std::endl;
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_OM | _FPU_MASK_UM);
        fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_DM | _FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT;
        _FPU_SETCW(cw);
        /* On x86, this expands to: */
        /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
        /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
    }
#endif


// -----------------------------------------------------------------------------
// Setup and shutdown routines. This is the main()
// function, but the real work is done in mainLoop().
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
#   ifdef CH_MPI
        MPI_Init(&argc, &argv);
#   endif

    // Reset the stdout color
    cout << color::none << std::flush;

    // Open the input file
    char* in_file = argv[1];
    ParmParse pp(argc-2, argv+2, NULL, in_file);
    pout() << "\nReading from file: " << in_file << std::endl;

    ParmParse ppLES("les");

    Real lesProcFrac = 0.0;
#   ifdef USE_LES
        ppLES.query("proc_frac", lesProcFrac);
        pout() << "\tles.proc_frac = " << lesProcFrac << endl;
#   endif

#   ifndef CH_MPI
        if (lesProcFrac != 0.0) {
            MayDay::Warning("Cannot run LES without MPI. Continuing les.proc_frac = 0.0");
            lesProcFrac = 0.0;
        }
#   endif

    // This creates groups and comms for communication among amr and les.
    // Call this even if we are not using MPI.
    setupMPIComms(lesProcFrac);
    testMPIComms();

#   ifdef CH_MPI
#       ifdef CH_AIX
            H5dont_atexit();
#       endif
#   endif

#   ifdef TRAP_FPE
        trapfpe();
#   endif

#   ifndef CH_NTIMER
        // Start yer timers.
        OldTimer Everything;
        Everything.start();
#   endif

// #   ifndef NDEBUG
//     {
//         // Get the host name
//         utsname hostInfo;
//         int errCode = uname(&hostInfo);
//         if (errCode != 0) {
//             pout() << "\nerrCode = " << errCode << " recieved from uname function." << endl;
//         } else {
//             pout() << "\nHost info:"
//                    << "\n  sysname    = " << hostInfo.sysname
//                    << "\n  nodename   = " << hostInfo.nodename
//                    << "\n  release    = " << hostInfo.release
//                    << "\n  version    = " << hostInfo.version
//                    << "\n  machine    = " << hostInfo.machine
//                    << endl;
// #           ifdef _GNU_SOURCE
//                 pout() << "  domainname = " << hostInfo.domainname << endl;
// #           endif
//         }

//         // Only register the debugger on certain systems.
//         bool useDebugger = true;
//         useDebugger |= std::string(hostInfo.nodename) == std::string("iceman");
//         useDebugger |= (std::string(hostInfo.nodename) == std::string("scorpion") && numProc() == 1);
//         if (useDebugger) {
//             pout() << "Registering debugger..." << endl;
//             registerDebugger();
//         }
//     }
// #   endif
//     pout() << endl;


    // BEGIN: Chombo-only code.
    if (AMRLESMeta::amrSize > 0) {
        if (AMRLESMeta::isGroupMember(AMRLESMeta::amrGroup)) {
#           ifdef CH_MPI
                pout() << "Using MPI." << endl;

                int chomboRank = procID();
                int chomboNumRanks = numProc();
                pout() << "\nchomboRank = " << chomboRank
                       << "\nchomboNumRanks = " << chomboNumRanks
                       << endl;
#           else
                pout() << "Not using MPI." << endl;
#           endif

            // Make sure the user knows if we are in debug mode.
#           ifndef NDEBUG
                pout() << "Proc #" << AMRLESMeta::worldRank
                       << "/" << AMRLESMeta::worldSize-1
                       << ": *** DEBUG MODE ***"
                       << endl;
#           else
                pout() << "Proc #" << AMRLESMeta::worldRank
                       << "/" << AMRLESMeta::worldSize-1
                       << ": *** RELEASE MODE ***"
                       << endl;
#           endif

            pout() << "SpaceDim = " << SpaceDim << endl;

#           ifdef  CH_USE_MEMORY_TRACKING
                pout() << "Using memory tracking." << endl;
#           endif

            // Get the host name
            utsname hostInfo;
            int errCode = uname(&hostInfo);
            if (errCode != 0) {
                pout() << "\nerrCode = " << errCode << " recieved from uname function." << endl;
            } else {
                pout() << "\nHost info:"
                       << "\n  sysname    = " << hostInfo.sysname
                       << "\n  nodename   = " << hostInfo.nodename
                       << "\n  release    = " << hostInfo.release
                       << "\n  version    = " << hostInfo.version
                       << "\n  machine    = " << hostInfo.machine
                       << endl;
#               ifdef _GNU_SOURCE
                    pout() << "  domainname = " << hostInfo.domainname << endl;
#               endif
            }

#           ifndef NDEBUG
                // Only register the debugger on certain systems.
                bool useDebugger = true;
                useDebugger |= std::string(hostInfo.nodename) == std::string("iceman");
                useDebugger |= (std::string(hostInfo.nodename) == std::string("scorpion") && numProc() == 1);
                if (useDebugger) {
                    pout() << "Registering debugger..." << endl;
                    registerDebugger();
                }
#           endif
            pout() << endl;

            // Setup AMR and run the simulation
            nsrun();

            // MPI_Barrier(AMRLESMeta::amrComm);
            printf("worldRank %d: AMR code finished.\n", AMRLESMeta::worldRank);
        }
    }
    // END: Chombo-only code.

#ifdef USE_LES
    // BEGIN: LES-only code.
    if (AMRLESMeta::lesSize > 0) {
        if (AMRLESMeta::isGroupMember(AMRLESMeta::lesGroup)) {

            // Register the debugger?
#           ifndef NDEBUG
                if (false) {
                    registerDebugger();
                }
#           endif

             // FORTRAN_NAME(DIABLO,diablo)(
             //     &AMRLESMeta::lesComm,
             //     &AMRLESMeta::interComm,
             //     &AMRLESMeta::amr2lesLeader,
             //     &AMRLESMeta::les2amrLeader
             // );

             FORTRAN_NAME(EDDY,eddy)(
                 &AMRLESMeta::lesComm,
                 &AMRLESMeta::interComm,
                 &AMRLESMeta::amr2lesLeader,
                 &AMRLESMeta::les2amrLeader
             );

            // MPI_Barrier(AMRLESMeta::lesComm);
            printf("worldRank %d: LES code finished.\n", AMRLESMeta::worldRank);
        }
    }
    // End LES-only code
#endif //USE_LES

    cout << flush;
    MPI_Barrier(MPI_COMM_WORLD);

#   ifndef CH_NTIMER
        // Stop timing.
        Everything.stop();

        Real end_memory = get_memory_usage_from_OS();
        std::string end_time = formatTime(Everything.wc_time());

        pout()  << "\nEverything done.\n"
                << "mem usage: " << end_memory << "MB\n"
                << "elapsed time: " << end_time << endl;
        tout(0) << "Elapsed time = " << end_time << "." << endl;
#   endif

#   ifdef CH_USE_MEMORY_TRACKING
        dumpmemoryatexit();
#   endif

#   ifndef CH_NTIMER
        CH_TIMER_REPORT();
#   endif

#   ifdef CH_MPI
        MPI_Finalize();
#   endif

    return 0;
}


// -----------------------------------------------------------------------------
// This creates groups and comms for communication among amr and les.
// -----------------------------------------------------------------------------
void setupMPIComms (const Real a_lesProcFrac)
{
    using namespace AMRLESMeta;

    // Gather world info
    worldComm = MPI_COMM_WORLD;
    worldGroup = getCommGroup(MPI_COMM_WORLD);
    worldSize = getGroupSize(worldGroup);
    worldRank = getGroupRank(worldGroup);


    // Compute the sizes of the groups.
    amrSize = int(round( Real(worldSize)*(1.0-a_lesProcFrac) ));
    if (amrSize < 0 || worldSize < amrSize) {
        MayDay::Error("setupMPIComms produced an amrSize that is out of range");
    }

    lesSize = worldSize - amrSize;
    if (lesSize < 0 || worldSize < lesSize) {
        MayDay::Error("setupMPIComms produced an lesSize that is out of range");
    }

    // Now, we have to split worldComm.
    int membershipKey = ((worldRank < amrSize)? 0: 1);
    MPI_Comm splitComm;
    if (MPI_Comm_split(worldComm, membershipKey, worldRank, &splitComm) != MPI_SUCCESS) {
        MayDay::Error("setupMPIComms could not split MPI_COMM_WORLD");
    }
    if (worldRank < amrSize) {
        amrComm = splitComm;
        amrGroup = getCommGroup(amrComm);
        amrRank = getGroupRank(amrGroup);

        lesComm = MPI_COMM_NULL;
        lesGroup = MPI_GROUP_NULL;
        lesRank = MPI_UNDEFINED_RANK;

    } else {
        lesComm = splitComm;
        lesGroup = getCommGroup(lesComm);
        lesRank = getGroupRank(lesGroup);

        amrComm = MPI_COMM_NULL;
        amrGroup = MPI_GROUP_NULL;
        amrRank = MPI_UNDEFINED_RANK;
    }

    // Set Chombo's communicator
    if (isGroupMember(amrGroup)) {
        Chombo_MPI::comm = amrComm;
    }

    // Create the intercommunicator for the AMR and LES groups.
    if (amrSize > 0 && lesSize > 0) {
        amr2lesLeader = 0; // This is the rank in amrComm.
        les2amrLeader = 0; // This is the rank in lesComm.
        int tag = 435;

        // Use a dedicated peer communicator
        if (MPI_Comm_dup(MPI_COMM_WORLD, &amrlesPeerComm) != MPI_SUCCESS) {
            MayDay::Error("setupMPIComms had an error creating amrlesPeerComm");
        }

        if (isGroupMember(amrGroup)) {
            int remoteLeader = amrSize; // The local group specifies the worldRank.
            if (MPI_Intercomm_create(splitComm, amr2lesLeader, amrlesPeerComm, remoteLeader, tag, &interComm) != MPI_SUCCESS) {
                MayDay::Error("setupMPIComms had an error creating the AMR to LES intercommunicator");
            }

            // Collective communications never want a local rank, they want a special macro.
            amr2lesLeader = ((amrRank == amr2lesLeader)? MPI_ROOT: MPI_PROC_NULL);

        } else {
            int remoteLeader = 0; // The remote group specifies the remoteRank  (remote = les).
            if (MPI_Intercomm_create(splitComm, amr2lesLeader, amrlesPeerComm, remoteLeader, tag, &interComm) != MPI_SUCCESS) {
                MayDay::Error("setupMPIComms had an error creating the AMR to LES intercommunicator");
            }

            // Collective communications never want a local rank, they want a special macro.
            les2amrLeader = ((lesRank == les2amrLeader)? MPI_ROOT: MPI_PROC_NULL);
        }
    } else {
        amr2lesLeader = MPI_UNDEFINED_RANK;
        amrlesPeerComm = MPI_COMM_NULL;
        interComm = MPI_COMM_NULL;
    }

    // Name the communicators to help debugging along.
    if (isGroupMember(amrGroup)) {
        MPI_Comm_set_name(amrComm, "amrComm");
    }
    if (isGroupMember(lesGroup)) {
        MPI_Comm_set_name(lesComm, "lesComm");
    }
    if (amrSize > 0 && lesSize > 0) {
        MPI_Comm_set_name(interComm, "interComm");
    }
}


// -----------------------------------------------------------------------------
// Passes magic values among the AMR and LES groups to test the
// intercommunicator.
// -----------------------------------------------------------------------------
void testMPIComms () {
    using namespace AMRLESMeta;

    // Do we have an intercommunicator?
    if (amrSize <= 0 || lesSize <= 0) return;

    // Make sure Chombo is using amrComm.
    if (isGroupMember(amrGroup)) {
        if (Chombo_MPI::comm != amrComm) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm");
        }
        if (procID() != amrRank) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm to compute ranks");
        }
        if (numProc() != amrSize) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm to compute comm size");
        }
    }

    // AMR to LES intercommunication
    {
        const double targetMagic = 12345.09876;

        double localMagic = -1.0;
        if (amr2lesLeader == MPI_ROOT) {
            localMagic = targetMagic;
        }

        if (MPI_Bcast(&localMagic, 1, MPI_DOUBLE, amr2lesLeader, interComm) != MPI_SUCCESS) {
            MayDay::Error("testMPIComms could not broadcast a magic number from AMR to LES");
        }

        if (isGroupMember(lesGroup)) {
            if (localMagic != targetMagic) {
                MayDay::Error("testMPIComms received an invalid magic number while testing AMR to LES intercommunication");
            }
        } else {
            if (amr2lesLeader == MPI_ROOT) {
                if (localMagic != targetMagic) {
                    MayDay::Error("testMPIComms clobbered the magic number while testing AMR to LES intercommunication");
                }
            } else {
                if (localMagic != -1.0) {
                    MayDay::Error("testMPIComms broadcasted the magic number to too many ranks during AMR to LES intercommunication");
                }
            }
        }
    }

    // AMR to LES intercommunication
    {
        const double targetMagic = 92753.41608;

        double localMagic = -1.0;
        if (les2amrLeader == MPI_ROOT) {
            localMagic = targetMagic;
        }

        if (MPI_Bcast(&localMagic, 1, MPI_DOUBLE, les2amrLeader, interComm) != MPI_SUCCESS) {
            MayDay::Error("testMPIComms could not broadcast a magic number from LES to AMR");
        }

        if (isGroupMember(amrGroup)) {
            if (localMagic != targetMagic) {
                MayDay::Error("testMPIComms received an invalid magic number while testing LES to AMR intercommunication");
            }
        } else {
            if (les2amrLeader == MPI_ROOT) {
                if (localMagic != targetMagic) {
                    MayDay::Error("testMPIComms clobbered the magic number while testing LES to AMR intercommunication");
                }
            } else {
                if (localMagic != -1.0) {
                    MayDay::Error("testMPIComms broadcasted the magic number to too many ranks during LES to AMR intercommunication");
                }
            }
        }
    }
}


// // -----------------------------------------------------------------------------
// void nsrun ()
// {
// #ifndef CH_NTIMER
//     // Time the setup routines
//     OldTimer setupTmr;
//     setupTmr.start();
// #endif

//     // This is probably the first request for a ProblemContext object.
//     // After this, the input file will be read and finished with.
//     const ProblemContext* ctx = ProblemContext::getInstance();

//     // Create AMRFactory object
//     AMRNavierStokesFactory amrns_fact;

//     // Create AMR object
//     LepticAMR thisAMR;
//     thisAMR.maxGridSize(ctx->maxGridSize);
//     thisAMR.maxBaseGridSize(ctx->maxBaseGridSize);
//     thisAMR.define(ctx->max_level, ctx->refRatios, ctx->domain, &amrns_fact);
//     thisAMR.verbosity(ctx->verbosity);
//     MappedAMRLevel::verbosity(ctx->verbosity);

//     thisAMR.useSubcyclingInTime(ctx->useSubcycling);

//     thisAMR.plotInterval(ctx->plot_interval);
//     thisAMR.plotPeriod(ctx->plot_period);
//     thisAMR.plotPrefix(ctx->plot_prefix);
//     thisAMR.checkpointInterval(ctx->checkpoint_interval);
//     thisAMR.checkpointPrefix(ctx->check_prefix);
//     thisAMR.gridBufferSize(ctx->bufferSize);

//     thisAMR.maxGridSize(ctx->maxGridSize);
//     thisAMR.maxBaseGridSize(ctx->maxBaseGridSize);
//     thisAMR.splitDirs(ctx->splitDirs);
//     thisAMR.fillRatio(ctx->fill_ratio);
//     thisAMR.blockFactor(ctx->block_factor);
//     thisAMR.regridIntervals(ctx->regrid_intervals);

//     if (ctx->fixed_dt > 0) {
//         thisAMR.fixedDt(ctx->fixed_dt);
//     }
//     thisAMR.maxDtGrow(ctx->max_dt_grow);

//     if (ctx->isRestart) {
//         // Initialize from restart file
// #ifdef CH_USE_HDF5
//         HDF5Handle handle(ctx->restart_file, HDF5Handle::OPEN_RDONLY);
//         thisAMR.setupForRestart(handle);
//         handle.close();
// #else
//         MayDay::Error("AMRNavierStokes restart only defined with HDF5");
// #endif
//     } else {
//         // New run
//         if (ctx->hasPredefinedGrids) {
//             thisAMR.setupForFixedHierarchyRun(ctx->predefinedGrids, 1);
//         } else {
//             thisAMR.setupForNewAMRRun();
//         }
//     }

// #ifndef CH_NTIMER
//     setupTmr.stop();
//     pout() << "Total setup time = " << formatTime(setupTmr.wc_time()) << endl;
// #endif

//     // Run
//     thisAMR.run(ctx->stopTime, ctx->maxsteps);

//     // Conclude
//     thisAMR.conclude();

//     // Free statically allocated memory.
//     LevelGeometry::staticUndefine();
//     LepticMeshRefine::deleteBuffer();
//     ProblemContext::freeMemory();
// }




// #include "Constants.H"
// #include "LepticMeshRefine.H"
// #include "AnisotropicRefinementTools.H"
// #include "Printing.H"
// #include "MappedFineInterp.H"
// // -----------------------------------------------------------------------------
// void nsrun ()
// {
//     // This is probably the first request for a ProblemContext object.
//     // After this, the input file will be read and finished with.
//     const ProblemContext* ctx = ProblemContext::getInstance();

//     const char* infile  = "/home/eds/research/docs/BVSolver/img/3DDJL/plot_000000.3d.hdf5";
//     const char* outfile = "/home/eds/research/docs/BVSolver/img/3DDJL/flat_000000.hdf5";

//     pout() << "\nProcessing " << infile << ":\n";
//     if (procID() == 0) {
//         std::cout << infile << ": opening..." << std::flush;
//     }

//     // Open the current input HDF5 file.
//     pout() << "\tOpening file..." << std::flush;
//     HDF5Handle handle(infile, HDF5Handle::OPEN_RDONLY);
//     pout() << "done." << endl;

//     HDF5HeaderData header;
//     header.readFromFile(handle);

//     // maxLevel
//     if (header.m_int.find("max_level") == header.m_int.end()) {
//         MayDay::Error("File does not contain max_level");
//     }
//     const int maxLevel = header.m_int["max_level"];
//     pout() << "\tmax_level = " << maxLevel << "\n";

//     // numLevels
//     if (header.m_int.find("num_levels") == header.m_int.end()) {
//         MayDay::Error("File does not contain num_levels");
//     }
//     const int numLevels = header.m_int["num_levels"];
//     pout() << "\tnum_levels = " << numLevels << "\n";

//     // iter
//     if (header.m_int.find("iteration") == header.m_int.end()) {
//         MayDay::Error("File does not contain iteration");
//     }
//     const int iter = header.m_int["iteration"];
//     pout() << "\titeration = " << iter << "\n";

//     // curTime
//     if (header.m_real.find("time") == header.m_real.end()) {
//         MayDay::Error("File does not contain time");
//     }
//     const Real curTime = header.m_real["time"];
//     pout() << "\tcomposite time = " << curTime << "\n";

//     // numComps
//     if (header.m_int.find("num_components") == header.m_int.end()) {
//         MayDay::Error("File does not contain num_components");
//     }
//     const int numComps = header.m_int["num_components"];
//     pout() << "\tnum_components = " << numComps << "\n";

//     // Get the components we are after.
//     std::map<string,string> compNames = header.m_string;

//     // Read data at each level
//     Vector<LevelData<FArrayBox>*> vData(numLevels, NULL);
//     Vector<RealVect> vdx(numLevels, RealVect::Zero);
//     Vector<IntVect> vRefRatio(numLevels, IntVect::Zero);
//     Real dt = quietNAN;

//     for (int lev = 0; lev < numLevels; ++lev) {
//         char levelStr[20];
//         sprintf(levelStr, "%d", lev);
//         string label = string("level_") + levelStr;
//         handle.setGroup(label);
//         header.readFromFile(handle);
//         pout() << "\tLevel " << lev << endl;

//         // refRatio
//         if (header.m_intvect.find("ref_ratio") == header.m_intvect.end()) {
//             MayDay::Error("File does not contain ref_ratio");
//         }
//         vRefRatio[lev] = header.m_intvect["ref_ratio"];
//         pout() << "\t\tref_ratio = " << vRefRatio[lev] << "\n";

//         // dx
//         if (header.m_realvect.find("vec_dx") == header.m_realvect.end()) {
//             MayDay::Error("File does not contain vec_dx");
//         }
//         vdx[lev] = header.m_realvect["vec_dx"];
//         pout() << "\t\tvec_dx = " << vdx[lev] << "\n";

//         // dt
//         if (header.m_real.find("dt") == header.m_real.end()) {
//             MayDay::Error("File does not contain dt");
//         }
//         dt = header.m_real["dt"];
//         pout() << "\t\tdt = " << dt << "\n";

//         // levelTime
//         if (header.m_real.find("time") == header.m_real.end()) {
//             MayDay::Error("File does not contain time");
//         }
//         Real levelTime = header.m_real["time"];
//         pout() << "\t\tlevel time = " << levelTime << "\n";

//         // domBox
//         if (header.m_box.find("prob_domain") == header.m_box.end()) {
//             MayDay::Error("File does not contain prob_domain");
//         }
//         const Box domBox = header.m_box["prob_domain"];
//         pout() << "\t\tprob_domain = " << domBox << "\n";

//         // grids
//         Vector<Box> boxArray;
//         read(handle, boxArray, "boxes");

//         DisjointBoxLayout grids;
//         grids.defineAndLoadBalance(boxArray, NULL, ProblemDomain(domBox));

//         // data
//         vData[lev] = new LevelData<FArrayBox>;
//         int err = read(handle, *vData[lev], "data", grids);
//         if (err == 1) {
//             // Bad location
//             ostringstream msg;
//             msg << "Bad location error when reading data from file " << infile << " level " << lev;
//             MayDay::Warning(msg.str().c_str());

//         } else if (err < 0) {
//             // HDF5 error
//             ostringstream msg;
//             msg << "HDF5 error " << err << " when reading data from file " << infile << " level " << lev;
//             MayDay::Warning(msg.str().c_str());
//         }
//     } // end loop over levels (lev)

//     { // Resize vectors
//         Vector<LevelData<FArrayBox>*> vTmpLD(numLevels, NULL);
//         Vector<RealVect> vTmpRV(numLevels);
//         Vector<IntVect> vTmpIV(numLevels);
//         for (int l = 0; l < numLevels; ++l) {
//             vTmpLD[l] = vData[l];
//             vTmpRV[l] = vdx[l];
//             vTmpIV[l] = vRefRatio[l];
//         }
//         vData = vTmpLD;
//         vdx = vTmpRV;
//         vRefRatio = vTmpIV;
//     }

//     // Create geometry (needs nx and L from input file)
//     Vector<LevelGeometry*> vLevGeo(numLevels, NULL);
//     // for (int lev = 0; lev < numLevels; ++lev) {
//     //     vLevGeo[lev] = new LevelGeometry(vdx[lev]);
//     //     vLevGeo[lev]->regrid(vData[lev]->getBoxes());
//     //     if (lev > 0) {
//     //         vLevGeo[lev]->setCoarserPtr(vLevGeo[lev-1]);
//     //     }
//     // }

//     // We are done with this. Clean it up.
//     handle.close();


//     Box newDomBox(IntVect::Zero, ctx->nx - IntVect::Unit);
//     newDomBox.setBig(1,0);

//     Vector<LevelData<FArrayBox>*> vSliceData(numLevels, NULL);

//     for (int lev = 0; lev <= maxLevel; ++lev) {
//         // Put old data on 1 proc
//         const DisjointBoxLayout& oldGrids = vData[lev]->getBoxes();
//         const ProblemDomain& oldDomain = oldGrids.physDomain();

//         DisjointBoxLayout oldGrids0;
//         const Vector<Box>& oldBoxArray = oldGrids.boxArray();
//         const Vector<int> vProc0(oldBoxArray.size(), 0);
//         oldGrids0.define(oldBoxArray, vProc0, oldDomain);

//         LevelData<FArrayBox> oldData0(oldGrids0, numComps, IntVect::Zero);
//         vData[lev]->copyTo(oldData0);

//         // Set up this level's flattened grids
//         pout() << "Level " << lev << " newDomBox = " << newDomBox << endl;
//         ProblemDomain newDomain(newDomBox);
//         DisjointBoxLayout newGrids;
//         newGrids.define(Vector<Box>(1,newDomBox), Vector<int>(1,0), newDomain);

//         vdx[lev][1] = 1.0;
//         vLevGeo[lev] = new LevelGeometry(vdx[lev]);
//         vLevGeo[lev]->regrid(newGrids);
//         if (lev > 0) {
//             vLevGeo[lev]->setCoarserPtr(vLevGeo[lev-1]);
//         }

//         // Copy data to flat holder.
//         vSliceData[lev] = new LevelData<FArrayBox>(newGrids, numComps, IntVect::Zero);

//         // Interp from coarser
//         if (lev > 0) {
//             MappedFineInterp interpObj(newGrids, numComps, vRefRatio[lev-1], newDomain);
//             interpObj.pwcinterpToFine(*vSliceData[lev], *vSliceData[lev-1]);
//         }

//         DataIterator oldDit = oldGrids0.dataIterator();
//         DataIterator newDit = newGrids.dataIterator();
//         for (oldDit.reset(), newDit.reset(); oldDit.ok(); ++oldDit) {
//             FArrayBox& newFAB = (*vSliceData[lev])[newDit];
//             const Box& newValid = newGrids[newDit];

//             const FArrayBox& oldFAB = oldData0[oldDit];
//             const Box& oldValid = oldGrids0[oldDit];

//             for (int c = 0; c < numComps; ++c) {
//                 IntVect oldCC = oldValid.smallEnd();
//                 IntVect newCC = oldCC;
//                 newCC[1] = 0;
//                 for (int k = oldValid.smallEnd(2); k <= oldValid.bigEnd(2); ++k) {
//                     oldCC[2] = k;
//                     newCC[2] = k;
//                     for (int i = oldValid.smallEnd(0); i <= oldValid.bigEnd(0); ++i) {
//                         oldCC[1] = i;
//                         oldCC[0] = i;
//                         newCC[0] = i;

//                         if (oldFAB.box().contains(oldCC)) {
//                             newFAB(newCC,c) = oldFAB(oldCC,c);
//                         }
//                     }
//                 }
//             }
//         }

//         // Refine to next level
//         if (lev < maxLevel) {
//             pout() << "Level " << lev << " refRatio = " << vRefRatio[lev] << endl;
//             vRefRatio[lev][1] = 1;
//             newDomBox.refine(vRefRatio[lev]);
//             // newDomBox.setBig(1,0);
//         }
//     }

//     Vector<string> vNames(numComps);
//     vNames[0] = "x_Vel";
//     vNames[1] = "y_Vel";
//     vNames[2] = "z_Vel";
//     vNames[3] = "mag_vel";
//     vNames[4] = "divergence";
//     vNames[5] = "lambda-1";
//     vNames[6] = "pressure";
//     vNames[7] = "x_Vort";
//     vNames[8] = "y_Vort";
//     vNames[9] = "z_Vort";
//     vNames[10] = "mag_vort";
//     vNames[11] = "scalar_0";
//     vNames[12] = "scalar_0_pert";

//     _writeHDF5(vSliceData, *vLevGeo[0], 0, maxLevel, outfile, vNames);


//     for (int lev = 0; lev < vSliceData.size(); ++lev) {
//         delete vSliceData[lev];
//         vSliceData[lev] = NULL;
//     }
//     vSliceData.resize(0);
//     for (int lev = 0; lev < vData.size(); ++lev) {
//         delete vData[lev];
//         vData[lev] = NULL;
//     }
//     vData.resize(0);
//     for (int lev = 0; lev < vLevGeo.size(); ++lev) {
//         if (vLevGeo[lev] != NULL) {
//             vLevGeo[lev]->reset();
//             delete vLevGeo[lev];
//             vLevGeo[lev] = NULL;
//         }
//     }
//     vLevGeo.resize(0);


//     // Free statically allocated memory.
//     LevelGeometry::staticUndefine();
//     LepticMeshRefine::deleteBuffer();
//     ProblemContext::freeMemory();
// }


#include "Constants.H"
#include "LepticMeshRefine.H"
#include "AnisotropicRefinementTools.H"
#include "Printing.H"
#include "MappedFineInterp.H"
#include "MappedCoarseAverage.H"
#include "computeMappedSum.H"
#include "AMRNSF_F.H"
#include "BoxIterator.H"
#include <dirent.h>
#include <fstream>
#include <utility>
#include "SPMDI.H"
#include "PhysBCUtil.H"
// -----------------------------------------------------------------------------
void nsrun ()
{
    // This is probably the first request for a ProblemContext object.
    // After this, the input file will be read and finished with.
    const ProblemContext* ctx = ProblemContext::getInstance();

#define DO_AUTO_FILES
#ifdef DO_AUTO_FILES

    // Choose which files to process.
    Vector<string> fileNames;
    {
        // Input parameters...
        const string hdf5Path("/home/eds/research/docs/BVSolver/img/3DDJL");
        const string hdf5Prefix("plot_");
        const string hdf5Suffix(".3d.hdf5");

        // Create a sorted set of files to read.
        std::set<string> files;
        std::set<string>::const_iterator it;
        {
            const int prefixLength = hdf5Prefix.length();
            const int suffixLength = hdf5Suffix.length();
            const int minLength = prefixLength + suffixLength;

            DIR *dir;
            struct dirent *ent;
            if ((dir = opendir(hdf5Path.c_str())) != NULL) {
                // Found the directory. Now collect files.
                while ((ent = readdir(dir)) != NULL) {
                    string thisFile(ent->d_name);

                    // The file name must at least have enough characters to
                    // include the prefix and suffix.
                    if (thisFile.length() < minLength) continue;

                    // The filename must have the correct prefix.
                    string thisFilePrefix = thisFile.substr(0, prefixLength);
                    if (thisFilePrefix.compare(hdf5Prefix) != 0) continue;

                    // The filename must have the correct suffix.
                    string thisFileSuffix = thisFile.substr(thisFile.length()-suffixLength);
                    if (thisFileSuffix.compare(hdf5Suffix) != 0) continue;

                    // Looks like we have a match. Throw it on the set.
                    files.insert(thisFile);
                }
                closedir(dir);
            } else {
                // Did not find directory.
                MayDay::Error("Could not open HDF5 input directory");
            }

            // Print the results
            pout() << "HDF5 input...\n"
                   << "\tpath = " << hdf5Path << "\n"
                   << "\tfiles = \n";
            for (it = files.begin(); it != files.end(); ++it) {
                pout() << "\t\t" << *it << "\n";
            }
            pout() << "\n" << endl;
        }

        // Loop over the files
        int idx = 0;
        it = files.begin();
        for (; it != files.end(); ++it, ++idx) {
            // Compute input path+filename
            const string& curFile = *it;
            string curFileAndPath = hdf5Path;
            if (curFileAndPath.substr(curFileAndPath.length()-1,1).compare(string("/")) != 0) {
                curFileAndPath += "/";
            }
            curFileAndPath += curFile;
            fileNames.push_back(curFileAndPath);
        }
    }
    const int numFiles = fileNames.size();

#else

    // Choose which files to process.
    Vector<int> fileNums;
    fileNums.push_back(0);
    // fileNums.push_back(382);
    // fileNums.push_back(766);
    // fileNums.push_back(1150);
    const string filePath = "/home/eds/research/docs/BVSolver/img/3DDJL";

    // Construct full file path+names.
    const int numFiles = fileNums.size();
    Vector<string> fileNames(numFiles);
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx) {
        char thisFileName[128];
        sprintf(thisFileName, "%s/plot_%06d.3d.hdf5", filePath.c_str(), fileNums[fileIdx]);
        fileNames[fileIdx] = thisFileName;
        pout() << "fileNames[" << fileIdx << "] = " << fileNames[fileIdx] << endl;
    }
#endif


#define WRITE_CSV_DATA
#ifdef WRITE_CSV_DATA
    std::ofstream csvFile;
    csvFile.open("TE_slice.csv");
    csvFile << "time,x0,intKE,intPE,intAPE,intAPEt" << endl;
#endif

    // Loop over files and process
    for (int fileIdx = 0; fileIdx < numFiles; ++fileIdx) {
        const char* infile = fileNames[fileIdx].c_str();

        // Tell the user which file we are processing.
        pout() << "\nProcessing " << infile << ":" << endl;
        if (procID() == 0) {
            std::cout << infile << ": opening..." << std::flush;
        }


        // Open the current input HDF5 file.
        pout() << "\tOpening file..." << std::flush;
        HDF5Handle handle(infile, HDF5Handle::OPEN_RDONLY);
        pout() << "done." << endl;

        HDF5HeaderData header;
        header.readFromFile(handle);

        // maxLevel
        if (header.m_int.find("max_level") == header.m_int.end()) {
            MayDay::Error("File does not contain max_level");
        }
        const int maxLevel = header.m_int["max_level"];
        pout() << "\tmax_level = " << maxLevel << "\n";

        // numLevels
        if (header.m_int.find("num_levels") == header.m_int.end()) {
            MayDay::Error("File does not contain num_levels");
        }
        const int numLevels = header.m_int["num_levels"];
        pout() << "\tnum_levels = " << numLevels << "\n";

        // iter
        if (header.m_int.find("iteration") == header.m_int.end()) {
            MayDay::Error("File does not contain iteration");
        }
        const int iter = header.m_int["iteration"];
        pout() << "\titeration = " << iter << "\n";

        // curTime
        if (header.m_real.find("time") == header.m_real.end()) {
            MayDay::Error("File does not contain time");
        }
        const Real curTime = header.m_real["time"];
        pout() << "\tcomposite time = " << curTime << "\n";

        // numComps
        if (header.m_int.find("num_components") == header.m_int.end()) {
            MayDay::Error("File does not contain num_components");
        }
        const int numComps = header.m_int["num_components"];
        pout() << "\tnum_components = " << numComps << "\n";

        // Get the components we are after.
        std::map<string,string> compNames = header.m_string;
        {
            pout() << "\tComponents:\n";
            char comp_str[30];
            for (int comp = 0; comp < numComps; ++comp) {
                sprintf (comp_str, "component_%d", comp);
                pout() << "\t\t" << comp_str << " --> " << header.m_string[comp_str] << "\n";
            }
            pout() << endl;
        }

        // Read data at each level
        Vector<LevelData<FArrayBox>*> vData(numLevels, NULL);
        Vector<RealVect> vdx(numLevels, RealVect::Zero);
        Vector<IntVect> vRefRatio(numLevels, IntVect::Zero);
        Real dt = quietNAN;

        for (int lev = 0; lev < numLevels; ++lev) {
            char levelStr[20];
            sprintf(levelStr, "%d", lev);
            string label = string("level_") + levelStr;
            handle.setGroup(label);
            header.readFromFile(handle);
            pout() << "\tLevel " << lev << endl;

            // refRatio
            if (header.m_intvect.find("ref_ratio") == header.m_intvect.end()) {
                MayDay::Error("File does not contain ref_ratio");
            }
            vRefRatio[lev] = header.m_intvect["ref_ratio"];
            pout() << "\t\tref_ratio = " << vRefRatio[lev] << "\n";

            // dx
            if (header.m_realvect.find("vec_dx") == header.m_realvect.end()) {
                MayDay::Error("File does not contain vec_dx");
            }
            vdx[lev] = header.m_realvect["vec_dx"];
            pout() << "\t\tvec_dx = " << vdx[lev] << "\n";

            // dt
            if (header.m_real.find("dt") == header.m_real.end()) {
                MayDay::Error("File does not contain dt");
            }
            dt = header.m_real["dt"];
            pout() << "\t\tdt = " << dt << "\n";

            // levelTime
            if (header.m_real.find("time") == header.m_real.end()) {
                MayDay::Error("File does not contain time");
            }
            Real levelTime = header.m_real["time"];
            pout() << "\t\tlevel time = " << levelTime << "\n";

            // domBox
            if (header.m_box.find("prob_domain") == header.m_box.end()) {
                MayDay::Error("File does not contain prob_domain");
            }
            const Box domBox = header.m_box["prob_domain"];
            pout() << "\t\tprob_domain = " << domBox << "\n";

            // grids
            Vector<Box> boxArray;
            read(handle, boxArray, "boxes");

            DisjointBoxLayout grids;
            grids.defineAndLoadBalance(boxArray, NULL, ProblemDomain(domBox));

            // data
            vData[lev] = new LevelData<FArrayBox>;
            int err = read(handle, *vData[lev], "data", grids);
            if (err == 1) {
                // Bad location
                ostringstream msg;
                msg << "Bad location error when reading data from file " << infile << " level " << lev;
                MayDay::Warning(msg.str().c_str());

            } else if (err < 0) {
                // HDF5 error
                ostringstream msg;
                msg << "HDF5 error " << err << " when reading data from file " << infile << " level " << lev;
                MayDay::Warning(msg.str().c_str());
            }
        } // end loop over levels (lev)

        { // Resize vectors to only include extant levels
            Vector<LevelData<FArrayBox>*> vTmpLD(numLevels, NULL);
            Vector<RealVect> vTmpRV(numLevels);
            Vector<IntVect> vTmpIV(numLevels);
            for (int l = 0; l < numLevels; ++l) {
                vTmpLD[l] = vData[l];
                vTmpRV[l] = vdx[l];
                vTmpIV[l] = vRefRatio[l];
            }
            vData = vTmpLD;
            vdx = vTmpRV;
            vRefRatio = vTmpIV;
        }

        // We are done with this. Clean it up.
        handle.close();


        // Create geometry (needs nx and L from input file)
        Vector<LevelGeometry*> vLevGeo(numLevels, NULL);
        for (int lev = 0; lev < numLevels; ++lev) {
            vLevGeo[lev] = new LevelGeometry(vdx[lev]);
            vLevGeo[lev]->regrid(vData[lev]->getBoxes());
            if (lev > 0) {
                vLevGeo[lev]->setCoarserPtr(vLevGeo[lev-1]);
            }
        }

        // Coarsen the data down
        for (int lev = maxLevel; lev > 0; --lev) {
            MappedCoarseAverage avgObj(vLevGeo[lev]->getBoxes(),
                                       numComps,
                                       vLevGeo[lev]->getCrseRefRatio());
            avgObj.averageToCoarse(*vData[lev-1], *vData[lev]);
        }


        // 1. Count cells
        if (procID() == 0) std::cout << "counting cells..." << std::flush;

        pout() << endl;
        long long amrNumPts = 0;
        for (int lev = 0; lev < numLevels; ++lev) {
            long long levNumPts = 0;
            const Vector<Box>& levelBoxes = vData[lev]->getBoxes().boxArray();
            for (int ll = 0; ll < levelBoxes.size(); ll++) {
                levNumPts += levelBoxes[ll].numPts();
            }
            pout() << "\tNumber of cells on level " << lev << " = " << levNumPts << endl;

            amrNumPts += levNumPts;
        }
        pout() << "\tTotal number of cells in AMR grids = " << amrNumPts << endl;
        const long long nonAmrNumPts = vData[maxLevel]->getBoxes().physDomain().domainBox().numPts();
        pout() << "\tNumber of cells needed without AMR = " << nonAmrNumPts << endl;
        const double percentCells = 100.0 * double(amrNumPts) / double(nonAmrNumPts);
        pout() << "\tPercent AMR cells / non-AMR cells = " << percentCells << endl;


        // 2. Compute domain energy integrals
        if (procID() == 0) std::cout << "integrating..." << std::flush;

        // Construct kinetic energy
        Vector<LevelData<FArrayBox>*> vKE(numLevels, NULL);
        Vector<LevelData<FArrayBox>*> vPE(numLevels, NULL);
        Vector<LevelData<FArrayBox>*> vAPE(numLevels, NULL);
        Vector<LevelData<FArrayBox>*> vAPEt(numLevels, NULL);
        for (int lev = 0; lev < numLevels; ++lev) {
            const DisjointBoxLayout& grids = vData[lev]->getBoxes();
            DataIterator dit = grids.dataIterator();

            vKE[lev] = new LevelData<FArrayBox>(grids, 1);
            vPE[lev] = new LevelData<FArrayBox>(grids, 1);
            vAPE[lev] = new LevelData<FArrayBox>(grids, 1);
            vAPEt[lev] = new LevelData<FArrayBox>(grids, 1);

            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& KEFAB = (*vKE[lev])[dit];
                FArrayBox& PEFAB = (*vPE[lev])[dit];
                FArrayBox& APEFAB = (*vAPE[lev])[dit];
                FArrayBox& APEtFAB = (*vAPEt[lev])[dit];
                const FArrayBox velFAB(Interval(0,2), (*vData[lev])[dit]);
                const FArrayBox bFAB(Interval(11,11), (*vData[lev])[dit]);
                const FArrayBox bprimeFAB(Interval(12,12), (*vData[lev])[dit]);
                const FArrayBox& gdnFAB = vLevGeo[lev]->getCCgdn()[dit];
                const Box& valid = grids[dit];
                const RealVect dx = vdx[lev];

                // KE
                FORT_COMPUTEKINETICENERGY(
                    CHF_FRA1(KEFAB,0),
                    CHF_CONST_FRA(velFAB),
                    CHF_BOX(valid),
                    CHF_CONST_FRA(gdnFAB),
                    CHF_CONST_REALVECT(dx));

                // PE
                FArrayBox cartPosFAB(valid, SpaceDim);
                vLevGeo[lev]->fill_physCoor(cartPosFAB);

                FORT_COMPUTEPOTENTIALENERGY(
                    CHF_FRA1(PEFAB,0),
                    CHF_CONST_FRA1(bFAB,0),
                    CHF_CONST_FRA1(cartPosFAB,CH_SPACEDIM-1),
                    CHF_BOX(valid));

                // APE
                FArrayBox bbarFAB(grow(valid,1), 1);
                bbarFAB.copy(bFAB);
                bbarFAB.plus(bprimeFAB, -1.0);

                FORT_COMPUTEAPE(
                    CHF_FRA1(APEFAB,0),
                    CHF_CONST_FRA1(bFAB,0),
                    CHF_CONST_FRA1(bbarFAB,0),
                    CHF_BOX(valid),
                    CHF_CONST_REALVECT(dx));

                // APEt
                FORT_COMPUTEPOTENTIALENERGY(
                    CHF_FRA1(APEtFAB,0),
                    CHF_CONST_FRA1(bprimeFAB,0),
                    CHF_CONST_FRA1(cartPosFAB,CH_SPACEDIM-1),
                    CHF_BOX(valid));
            }
        }

        // Coarsen the KE down
        for (int lev = maxLevel; lev > 0; --lev) {
            MappedCoarseAverage avgObj(vLevGeo[lev]->getBoxes(),
                                       1,
                                       vLevGeo[lev]->getCrseRefRatio());
            avgObj.averageToCoarse(*vKE[lev-1], *vKE[lev]);
            avgObj.averageToCoarse(*vPE[lev-1], *vPE[lev]);
            avgObj.averageToCoarse(*vAPE[lev-1], *vAPE[lev]);
            avgObj.averageToCoarse(*vAPEt[lev-1], *vAPEt[lev]);
        }

        // // Compute integrals
        // Real vol = 0.0;
        // Real sumKE = computeMappedSum(vol, vKE, *vLevGeo[0]);
        // pout() << "\tAMR sum KE = " << sumKE << endl;
        // pout() << "\tAMR vol = " << vol << endl;

        // for (int lev = 0; lev < numLevels; ++lev) {
        //     vol = 0.0;
        //     sumKE = computeMappedSum(vol, *vKE[lev], NULL, *vLevGeo[lev]);
        //     pout() << "\tlevel " << lev << " sum KE = " << sumKE << endl;
        //     pout() << "\tlevel " << lev << " vol = " << vol << endl;
        // }



        // Find the center of the wave
        Real xp0 = 0.0;
        Real amrMaxKE = -123456.0;
        const Real rotAngle = 45.0 * Pi / 180.0;
        const Real sinA = sin(rotAngle);
        const Real cosA = cos(rotAngle);

        for (int lev = 0; lev < numLevels; ++lev) {
            const DisjointBoxLayout& grids = vLevGeo[lev]->getBoxes();
            DataIterator dit = grids.dataIterator();

            Real localXp0 = 0.0;
            Real localMaxKE = -123456.0;
            for (dit.reset(); dit.ok(); ++dit) {
                const FArrayBox velFAB(Interval(0,2), (*vData[lev])[dit]);
                const FArrayBox& KEFAB = (*vKE[lev])[dit];
                const Box& valid = grids[dit];
                const RealVect dx = vdx[lev];
                BoxIterator bit(valid);

                {
                    IntVect cc = valid.smallEnd();
                    for (cc[1] = valid.smallEnd(1); cc[1] <= valid.bigEnd(1); ++cc[1]) {
                        for (cc[0] = valid.smallEnd(0); cc[0] <= valid.bigEnd(0); ++cc[0]) {
                            // Real thisKE = 0.0;
                            // for (cc[2] = valid.smallEnd(2); cc[2] <= valid.bigEnd(2); ++cc[2]) {
                            //     thisKE += KEFAB(cc,0);
                            // }
                            // thisKE *= dx[2];

                            Box intBox(cc,cc);
                            intBox.setBig(2,valid.bigEnd(2));
                            Real thisKE = KEFAB.sum(intBox,0) * dx[2];

                            if (thisKE < localMaxKE) continue;

                            Real x = (Real(cc[0]) + 0.5) * dx[0];// - 128.0;// - 20.0;
                            Real y = (Real(cc[1]) + 0.5) * dx[1];// - 128.0;
                            localXp0 =  x*cosA + y*sinA;

                            localMaxKE = thisKE;
                        }
                    }
                }
            }

            Vector<Real> thisMaxKE(2);
            thisMaxKE[0] = localMaxKE;
            thisMaxKE[1] = localXp0;

            Vector<Vector<Real> > procMaxKE(numProc());
            gather(procMaxKE, thisMaxKE, 0);
            if (procID() == 0) {
                for (int idx = 0; idx < procMaxKE.size(); ++idx) {
                    if (procMaxKE[idx][0] < thisMaxKE[0]) continue;
                    thisMaxKE = procMaxKE[idx];
                }
            }
            broadcast(thisMaxKE, 0);
            if (thisMaxKE[0] >= amrMaxKE) {
                amrMaxKE = thisMaxKE[0];
                xp0 = thisMaxKE[1];
            }
        }
        pout() << "\n\txp0 = " << xp0 << endl;

        // // Now, apply the mask
        // for (int lev = 0; lev < numLevels; ++lev) {
        //     const DisjointBoxLayout& grids = vLevGeo[lev]->getBoxes();
        //     DataIterator dit = grids.dataIterator();

        //     for (dit.reset(); dit.ok(); ++dit) {
        //         FArrayBox& KEFAB = (*vKE[lev])[dit];
        //         FArrayBox& PEFAB = (*vPE[lev])[dit];
        //         FArrayBox& APEFAB = (*vAPE[lev])[dit];
        //         FArrayBox& APEtFAB = (*vAPEt[lev])[dit];
        //         const FArrayBox velFAB(Interval(0,2), (*vData[lev])[dit]);
        //         const Box& valid = grids[dit];
        //         const RealVect dx = vdx[lev];
        //         BoxIterator bit(valid);

        //         // Mask the KE
        //         for (bit.reset(); bit.ok(); ++bit) {
        //             const IntVect cc = bit();
        //             Real x = (Real(cc[0]) + 0.5) * dx[0];// - 128.0;// - 20.0;
        //             Real y = (Real(cc[1]) + 0.5) * dx[1];// - 128.0;
        //             Real xp =  x*cosA + y*sinA;
        //             Real yp = -x*sinA + y*cosA;

        //             // if (cc[0] == cc[1]) {
        //             if(-14.0 < yp && yp < 14.0) {
        //                 if ((xp0-14.0) < xp && xp < (xp0+14.0)) {
        //                     // KEFAB(cc) = 1.0;
        //                 } else {
        //                     KEFAB(cc) = 0.0;
        //                     PEFAB(cc) = 0.0;
        //                     APEFAB(cc) = 0.0;
        //                     APEtFAB(cc) = 0.0;
        //                 }
        //             } else {
        //                 KEFAB(cc) = 0.0;
        //                 PEFAB(cc) = 0.0;
        //                 APEFAB(cc) = 0.0;
        //                 APEtFAB(cc) = 0.0;
        //             }
        //         }
        //     }
        // }
        {
            Vector<string> vNames(1,"comp 0");
            ostringstream fname;
            fname << "KE_" << fileIdx << ".hdf5"; \
            _writeHDF5(vKE,
                       *vLevGeo[0],
                       0, maxLevel,
                       fname.str().c_str(),
                       vNames);

            // Compute integrals
            Real sumSliceKE = computeMappedSum(vKE, *vLevGeo[0]);
            pout() << "\tAMR sum slice KE = " << sumSliceKE << endl;

            Real sumSlicePE = computeMappedSum(vPE, *vLevGeo[0]);
            pout() << "\tAMR sum slice PE = " << sumSlicePE << endl;

            Real sumSliceAPE = computeMappedSum(vAPE, *vLevGeo[0]);
            pout() << "\tAMR sum slice APE = " << sumSliceAPE << endl;

            Real sumSliceAPEt = computeMappedSum(vAPEt, *vLevGeo[0]);
            pout() << "\tAMR sum slice APEt = " << sumSliceAPEt << endl;

#ifdef WRITE_CSV_DATA
            csvFile << curTime      << ","
                    << xp0          << ","
                    << sumSliceKE   << ","
                    << sumSlicePE   << ","
                    << sumSliceAPE  << ","
                    << sumSliceAPEt << endl;
#endif
        }


        // Free memory
        for (int lev = 0; lev < vKE.size(); ++lev) {
            delete vKE[lev];
            vKE[lev] = NULL;
        }
        vKE.resize(0);

        for (int lev = 0; lev < vPE.size(); ++lev) {
            delete vPE[lev];
            vPE[lev] = NULL;
        }
        vPE.resize(0);

        for (int lev = 0; lev < vAPE.size(); ++lev) {
            delete vAPE[lev];
            vAPE[lev] = NULL;
        }
        vAPE.resize(0);

        for (int lev = 0; lev < vAPEt.size(); ++lev) {
            delete vAPEt[lev];
            vAPEt[lev] = NULL;
        }
        vAPEt.resize(0);

        for (int lev = 0; lev < vData.size(); ++lev) {
            delete vData[lev];
            vData[lev] = NULL;
        }
        vData.resize(0);

        for (int lev = 0; lev < vLevGeo.size(); ++lev) {
            if (vLevGeo[lev] != NULL) {
                vLevGeo[lev]->reset();
                delete vLevGeo[lev];
                vLevGeo[lev] = NULL;
            }
        }
        vLevGeo.resize(0);

        // Free statically allocated memory.
        LevelGeometry::staticUndefine();
        LepticMeshRefine::deleteBuffer();


        // Tell the user that we are done with this file.
        if (procID() == 0) {
            std::cout << "done." << std::endl;
        }

    } // end loop over files (fileIdx)

#ifdef WRITE_CSV_DATA
    csvFile.close();
#endif

    // Free statically allocated memory.
    LevelGeometry::staticUndefine();
    LepticMeshRefine::deleteBuffer();
    ProblemContext::freeMemory();

    barrier();




    // Vector<string> vNames(numComps);
    // vNames[0] = "x_Vel";
    // vNames[1] = "y_Vel";
    // vNames[2] = "z_Vel";
    // vNames[3] = "mag_vel";
    // vNames[4] = "divergence";
    // vNames[5] = "lambda-1";
    // vNames[6] = "pressure";
    // vNames[7] = "x_Vort";
    // vNames[8] = "y_Vort";
    // vNames[9] = "z_Vort";
    // vNames[10] = "mag_vort";
    // vNames[11] = "scalar_0";
    // vNames[12] = "scalar_0_pert";
}
