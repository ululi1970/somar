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
#include "PhysBCUtil.H"
#include "PhysBCUtilF_F.H"
#include "ProblemContext.H"
#include "EllipticBCUtils.H"
#include "LevelGeometry.H"
#include "BoxIterator.H"


// Declare static members
bool     PhysBCUtil::s_staticMembersSet = false;
RealVect PhysBCUtil::s_domLength = RealVect::Zero;
bool     PhysBCUtil::s_useBackgroundScalar = false;

bool     PhysBCUtil::s_useSpongeLayer = false;
Real     PhysBCUtil::s_spongeWidth[CH_SPACEDIM][2];
Real     PhysBCUtil::s_spongeDtMult[CH_SPACEDIM][2];


// ************************ Constructors / destructors *************************

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PhysBCUtil::PhysBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
PhysBCUtil::~PhysBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Reads BC data from file
// -----------------------------------------------------------------------------
void PhysBCUtil::define ()
{
    if (s_staticMembersSet) return;

    const ProblemContext* ctx = ProblemContext::getInstance();

    s_domLength = ctx->domainLength;
    s_useBackgroundScalar = ctx->useBackgroundScalar;
    s_useSpongeLayer = ctx->useSpongeLayer;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        s_spongeWidth[dir][0] = ctx->spongeWidth[dir][0];
        s_spongeWidth[dir][1] = ctx->spongeWidth[dir][1];
        s_spongeDtMult[dir][0] = ctx->spongeDtMult[dir][0];
        s_spongeDtMult[dir][1] = ctx->spongeDtMult[dir][1];
    }

    s_staticMembersSet = true;
}


// ************************ ICs / background fields ****************************

// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void PhysBCUtil::setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    a_velFAB.setVal(0.0, a_velComp);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void PhysBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const
{
    a_scalarFAB.setVal(0.0);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void PhysBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar) {
        MayDay::Error("PhysBCUtil::setBackgroundScalar needs to be overridden if you plan to use a background scalar");
    } else {
        a_scalarFAB.setVal(0.0);
    }
}


// ******************** Miscellaneous utility functions ************************

// -----------------------------------------------------------------------------
// Sets the velocity ICs over an entire level.
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void PhysBCUtil::setVelIC (LevelData<FArrayBox>& a_vel,
                           const int             a_velComp,
                           const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_vel.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setVelIC(a_vel[dit], a_velComp, a_levGeo, dit());
    }
}


// -----------------------------------------------------------------------------
// Sets the scalar ICs over an entire level
// -----------------------------------------------------------------------------
void PhysBCUtil::setScalarIC (LevelData<FArrayBox>& a_scalar,
                              const int             a_scalarComp,
                              const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_scalar.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setScalarIC(a_scalar[dit], a_scalarComp, a_levGeo, dit());
    }
}


// -----------------------------------------------------------------------------
// Sets the background scalars over an entire level
// -----------------------------------------------------------------------------
void PhysBCUtil::setBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_scalar.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setBackgroundScalar(a_scalar[dit], a_scalarComp, a_levGeo, dit(), a_time);
    }
}


// -----------------------------------------------------------------------------
// Adds the background scalar to an existing LevelData
// -----------------------------------------------------------------------------
void PhysBCUtil::addBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const
{
    CH_assert(0 <= a_scalarComp);

    if (s_useBackgroundScalar) {
        const DisjointBoxLayout& grids = a_scalar.getBoxes();
        DataIterator dit = grids.dataIterator();

        const int numComps = a_scalar.nComp();
        LevelData<FArrayBox> backgroundScal(grids, numComps, a_scalar.ghostVect());
        setBackgroundScalar(backgroundScal, a_scalarComp, a_time, a_levGeo);

        for (dit.reset(); dit.ok(); ++dit) {
            a_scalar[dit].plus(backgroundScal[dit], 0, 0, numComps);
        }
    }
}


// -----------------------------------------------------------------------------
// Subtracts the background scalar from an existing LevelData
// -----------------------------------------------------------------------------
void PhysBCUtil::subtractBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                           const int             a_scalarComp,
                                           const Real            a_time,
                                           const LevelGeometry&  a_levGeo) const
{
    CH_assert(0 <= a_scalarComp);

    if (s_useBackgroundScalar) {
        const DisjointBoxLayout& grids = a_scalar.getBoxes();
        DataIterator dit = grids.dataIterator();

        const int numComps = a_scalar.nComp();
        LevelData<FArrayBox> backgroundScal(grids, numComps, a_scalar.ghostVect());
        setBackgroundScalar(backgroundScal, a_scalarComp, a_time, a_levGeo);

        for (dit.reset(); dit.ok(); ++dit) {
            a_scalar[dit].minus(backgroundScal[dit], 0, 0, numComps);
        }
    }
}


// -----------------------------------------------------------------------------
// Computes the Brunt–Väisälä frequency on a single grid.
// All ins and outs must be CC.
// a_bFAB is the background buoyancy field.
// a_dXidz needs SpaceDim comps...(dXi/dz, dNu/dz, dZeta/dz).
// This function is a bit raw, but dXidz is expensive to compute and should
// be cached by the user.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeNSq (FArrayBox&       a_NsqFAB,
                             const FArrayBox& a_bFAB,
                             const FArrayBox& a_dXidzFAB,
                             const Box&       a_destBox,
                             const RealVect&  a_dx) const
{
    // Sanity checks
    CH_assert(a_NsqFAB  .nComp() == 1);
    CH_assert(a_bFAB    .nComp() == 1);
    CH_assert(a_dXidzFAB.nComp() == SpaceDim);

    CH_assert(a_NsqFAB  .box().type() == a_destBox.type());
    CH_assert(a_bFAB    .box().type() == a_destBox.type());
    CH_assert(a_dXidzFAB.box().type() == a_destBox.type());
    CH_assert(a_destBox       .type() == IntVect::Zero);

    CH_assert(a_NsqFAB  .box().contains(a_destBox));
    CH_assert(a_bFAB    .box().contains(grow(a_destBox,1)));
    CH_assert(a_dXidzFAB.box().contains(a_destBox));

    // Calculate -dXi^i/dz * dB/dXi^i
    FORT_COMPUTE_CCNSQ(CHF_FRA1(a_NsqFAB,0),
                       CHF_CONST_FRA1(a_bFAB,0),
                       CHF_CONST_FRA(a_dXidzFAB),
                       CHF_CONST_REALVECT(a_dx),
                       CHF_BOX(a_destBox));
}


// -----------------------------------------------------------------------------
// Computes the Brunt–Väisälä frequency over an entire level.
// a_Nsq must be CC.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeNSq (LevelData<FArrayBox>& a_Nsq,
                             const LevelGeometry&  a_levGeo,
                             const Real            a_time) const
{
    // Sanity checks
    CH_assert(a_levGeo.getBoxes().physDomain() == a_Nsq.getBoxes().physDomain());
    CH_assert(a_levGeo.getBoxes().compatible(a_Nsq.getBoxes()));

    // Declare variables
    const RealVect& dx = a_levGeo.getDx();
    const GeoSourceInterface& geoSource = *(a_levGeo.getGeoSourcePtr());
    DataIterator dit = a_Nsq.dataIterator();

    // Loop over grids and compute Nsq.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& NsqFAB = a_Nsq[dit];
        const Box& destBox = NsqFAB.box();
        CH_assert(destBox.type() == IntVect::Zero);

        // Fill dXi^i/dz vector.
        FArrayBox dXidzFAB(destBox, SpaceDim);
        for (int dir = 0; dir < SpaceDim; ++dir) {
            geoSource.fill_dXidx(dXidzFAB,
                                 dir,           // fab comp
                                 dir,           // Xi index
                                 SpaceDim-1,    // z index
                                 dx);
        }

        // Fill background bouyancy field.
        FArrayBox BFAB(grow(destBox,1), 1);
        this->setBackgroundScalar(BFAB,
                                  0,        // scalar comp
                                  a_levGeo,
                                  dit(),
                                  a_time);

        // Calculate -dXi^i/dz * dB/dXi^i
        this->computeNSq(NsqFAB,
                         BFAB,
                         dXidzFAB,
                         destBox,
                         dx);
    }
}


// -----------------------------------------------------------------------------
// This is in case the BC's have an effect on the timestep.
// Pass in currently computed dt, along with the cfl and dx. If the effect
// of the BCs requires a decreased timestep, then the newly reduced timestep
// is returned. In the default case, this just returns a_dt back; however,
// derived classes may actually have an effect.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeBoundaryDt (Real&                a_dt,
                                    const Real           a_cfl,
                                    const LevelGeometry& a_levGeo) const
{
    // In the default case, nothing is done here.
}


// ************************* Sponge layer functions ****************************


// -----------------------------------------------------------------------------
// Splits the domain into its sponge layers and the interior. The locations of
// the splitting are given as face indices. If the domain is not split, these
// indices will lie outside of the domain.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeSpongeRegions (Tuple<Box, 2>&       a_spongeBox,
                                       Tuple<int, 2>&       a_splitFaceIndex,
                                       Box&                 a_interior,
                                       const LevelGeometry& a_levGeo,
                                       const int            a_dir)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // Gather some needed info
    const Box domBox = a_levGeo.getDomain().domainBox();
    const RealVect& dx = a_levGeo.getDx();
    Box interiorBox[2];

    // Lo side
    int s = 0;
    if (s_spongeWidth[a_dir][s] > 0.0) {
        const int numCells = ceil(s_spongeWidth[a_dir][s] / dx[a_dir]);
        a_splitFaceIndex[s] = domBox.smallEnd(a_dir) + numCells;
        a_spongeBox[s] = domBox;
        interiorBox[s] = a_spongeBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Hi side
    s = 1;
    if (s_spongeWidth[a_dir][s] > 0.0) {
        const int numCells = ceil(s_spongeWidth[a_dir][s] / dx[a_dir]);
        a_splitFaceIndex[s] = domBox.bigEnd(a_dir) - numCells + 1;
        interiorBox[s] = domBox;
        a_spongeBox[s] = interiorBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Compute the interior region
    a_interior = interiorBox[0] & interiorBox[1];
}


// -----------------------------------------------------------------------------
// Sets a_srcTerm = (target - state) / (layer time scale).
// To fill velocity sponge layer sources, just use a_comp's default value.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillSpongeLayerSrcTerm (LevelData<FArrayBox>&       a_srcTerm,
                                         const LevelData<FArrayBox>& a_state,
                                         const Real                  a_time,
                                         const Real                  a_dt,
                                         const LevelGeometry&        a_levGeo,
                                         const int                   a_comp)
{
    CH_TIME("PhysBCUtil::fillSpongeLayerSrcTerm");

    // Gather some needed information
    const int ncomp = a_srcTerm.nComp();
    const RealVect& dx = a_levGeo.getDx();
    const DisjointBoxLayout& grids = a_srcTerm.disjointBoxLayout();
    DataIterator dit = grids.dataIterator();

// LevelData<FArrayBox> spongeProfile(grids, 1);
// setValLevel(spongeProfile, 0.0);

    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(ncomp == a_state.nComp());
    CH_assert(a_comp >=  0 || ncomp == SpaceDim);
    CH_assert(a_comp == -1 || ncomp == 1);
    CH_assert(a_state.getBoxes() == a_srcTerm.getBoxes());

    // Start with an inactive sponge layer, then activate the sides we need.
    for (dit.reset(); dit.ok(); ++dit) {
        a_srcTerm[dit].setVal(0.0);
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Calculate the sponge regions and the interior
        Tuple<Box, 2> spongeBox;
        Tuple<int, 2> splitFaceIndex;
        Box interior;
        this->computeSpongeRegions(spongeBox, splitFaceIndex, interior, a_levGeo, dir);

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = ((iside == Side::Lo)? 0: 1);

            // Do we have a sponge on this side?
            if (spongeBox[s].isEmpty()) continue;

            // Calculate the sponge region and it's complement
            const Real splitFaceLoc = Real(splitFaceIndex[s]) * dx[dir];

            // Set sponge layer term coefficient
            const Real coeff = 1.0 / (s_spongeDtMult[dir][s] * a_dt);

            // Loop over the grids.
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& srcTermFAB = a_srcTerm[dit];
                const FArrayBox& stateFAB = a_state[dit];

                const Box valid = grids[dit] & srcTermFAB.box();
                const Box validSpongeBox = valid & spongeBox[s];

                // If the sponge region is empty, just move on.
                if (validSpongeBox.isEmpty()) continue;

                // Compute the target locations
                const Box targetBox = bdryBox(validSpongeBox, dir, iside, 1);
                const int targetFaceIndex = targetBox.smallEnd(dir);
                CH_assert(targetFaceIndex == targetBox.bigEnd(dir));

                FArrayBox targetFAB(targetBox, 1);
                Real ratio, ramp, pos;

                // This is a bit strange. If a_comps is -1, then we use that as a
                // flag to fill the velocity sponge. comp will loop over all
                // SpaceDim elements of the velocity. If a_comps is anything else,
                // then ncomp = 1, comp = 0, and a_comp signals which scalar sponge
                // we need to fill.
                for (int comp = 0; comp < ncomp; ++comp) {
                    // Compute the target field.
                    // TODO: Should we pass dir in?
                    if (a_comp == -1) {
                        this->fillVelSpongeLayerTarget(targetFAB, comp, dir, iside, a_levGeo, dit(), a_time);
                    } else {
                        this->fillScalarSpongeLayerTarget(targetFAB, a_comp, dir, iside, a_levGeo, dit(), a_time);
                    }

                    // Loop over the grid and fill the sponge source term.
                    BoxIterator bit(validSpongeBox);
                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& cc = bit();

                        IntVect targetfc = cc;
                        targetfc[dir] = targetFaceIndex;

                        pos = (Real(cc[dir]) + 0.5) * dx[dir];
                        ratio = Abs(pos - splitFaceLoc) / s_spongeWidth[dir][s];
                        ramp = this->spongeLayerRamp(ratio);

                        srcTermFAB(cc,comp) = coeff * ramp * (targetFAB(targetfc) - stateFAB(cc,comp));

                    } // end loop over sponge box (bit)
                } // end loop over state components (comp)
            } // end loop over grids (dit)
        } // end loop over sides (sit)
    } // end loop over directions (dir)

// if (ncomp == 1) writeLevelHDF5(spongeProfile, a_dt, false);
}


// -----------------------------------------------------------------------------
// Sets the Cartesian-based target velocity for the sponge layer.
// By default, this function throws an error.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                           const int            a_velComp,
                                           const int            a_spongeDir,
                                           const Side::LoHiSide a_spongeSide,
                                           const LevelGeometry& a_levGeo,
                                           const DataIndex&     a_di,
                                           const Real           a_time)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(a_target.nComp() == 1);
    CH_assert(0 <= a_spongeDir);
    CH_assert(a_spongeDir < SpaceDim);

    // I don't know what values the velocity should approach.
    const char* msg = "If you plan to use a sponge layer, "
                      "then you need to override "
                      "PhysBCUtil::fillVelSpongeLayerTarget";
    MayDay::Error(msg);
}


// -----------------------------------------------------------------------------
// Sets the target values for the scalar sponge layer. If we are using a
// background scalar, then this function set the perturbation to zero.
// Otherwise, an error is thrown and this function will need to be overridden.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                              const int            a_scalarComp,
                                              const int            a_spongeDir,
                                              const Side::LoHiSide a_spongeSide,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(a_target.nComp() == 1);
    CH_assert(0 <= a_spongeDir);
    CH_assert(a_spongeDir < SpaceDim);

    if (s_useBackgroundScalar) {
        // The outgoing scalar should approach its unperturbed value.
        a_target.setVal(0.0, a_scalarComp);
    } else {
        // I don't know what values the scalar should approach.
        const char* msg = "If you plan to use a sponge layer without a background scalar, "
                          "then you need to override PhysBCUtil::fillScalarSpongeLayerTarget";
        MayDay::Error(msg);
    }
}


// ****************************** Velocity BCs *********************************

// -----------------------------------------------------------------------------
// uStarFuncBC
// Pre-projection velocity BC.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::uStarFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// viscousSourceFuncBC
// Sets ghosts needed to calculate the viscous source term nu.L[vel]
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::viscousSourceFuncBC () const
{
    const bool isViscous = true;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// viscousSolveFuncBC (Used in single-component velocity TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::viscousSolveFuncBC (int a_dir) const
{
    const bool isViscous = true;
    return this->basicVelFuncBC(a_dir, isViscous);
}


// -----------------------------------------------------------------------------
// viscousRefluxBC (Used in sync step by the implicit refluxing op)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::viscousRefluxBC (int a_dir) const
{
    const bool isViscous = true;
    return this->basicVelFuncBC(a_dir, isViscous);
}


// -----------------------------------------------------------------------------
// viscousVelFuncBC
// Sets BCs on a generic viscous velocity field.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::viscousVelFuncBC () const
{
    const bool isViscous = true;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// tracingVelFuncBC
// Sets BCs on old-time velocity that will be used to trace characteristics.
// (Used by fill functions)
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::tracingVelFuncBC () const
{
    const bool isViscous = false;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// uDelUFuncBC
// Sets BCs on the advection term
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::uDelUFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// advectingVelFuncBC
// Sets BCs on the advecting velocity
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::advectingVelFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// vortFuncBC
// Sets BCs on velocity when computing vorticity. We do first-order extrap here.
// The effect is that derivatives at the boundaries become 1-sided differences.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::vortFuncBC (bool a_isViscous) const
{
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// velRiemannBC
// Sets BCs on FC velocity in the Riemann solver.
// NOTE: This needs to set BCs on the velocity in the Cartesian basis.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::velRiemannBC (int  a_velComp,
                                                         bool a_isViscous) const
{
    // Really, we don't need to set fluxes because we never request
    // velocity fluxes from the Riemann solver. Also, the BCs for each
    // direction (idir) are the same since predictScalar works on only
    // one velocity component at a time.
    BCMethodHolder protoBC = this->basicVelFuncBC(a_velComp, false);
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = protoBC;
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// velSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// NOTE: This needs to set BCs on the velocity in the Cartesian basis.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::velSlopeBC (int a_velComp, bool a_isViscous) const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ******************************* Scalar BCs **********************************

// -----------------------------------------------------------------------------
// diffusiveSourceFuncBC (Used to calculate the diffusive term nu.L[scalar])
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::diffusiveSourceFuncBC () const
{
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// diffusiveSolveFuncBC (used in scalar TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::diffusiveSolveFuncBC () const
{
    // This works better than 2nd order extrap for at least Diri BCs.
    // extrap exposes the flaws at box boundaries.

    // NOTE: basicScalarFuncBC is 1st order extrap in this function.
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// scalarRefluxSolveBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::scalarRefluxSolveBC (int a_scalarType) const
{
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// scalarTraceFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::scalarTraceFuncBC (int a_scalarType) const
{
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// scalarRiemannBC
// Sets BCs on FC scalars in the Riemann solver.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::scalarRiemannBC (int a_scalarType) const
{
    // Start with standard BC
    BCMethodHolder protoBC = this->basicScalarFuncBC();

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// scalarSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::scalarSlopeBC (int a_scalarType) const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ********************** Freestream preservation BCs **************************

// -----------------------------------------------------------------------------
// lambdaFuncBC
// Chombo uses Diri 1.0 BCs at inflow and 1st order extrap elsewhere.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::lambdaFuncBC () const
{
    BCMethodHolder holder;

    // Solid wall and outflow BCs
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> otherBCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder,
                                       IntVect::Unit,
                                       IntVect::Unit)
    );
    holder.addBCMethod(otherBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// FreestreamCorrFuncBC
// Chombo uses 0 Neum
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::FreestreamCorrFuncBC () const
{
    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstNeumBCGhostClass(IntVect::Zero, IntVect::Zero)
    );
    holder.addBCMethod(BCPtr);

    RefCountedPtr<BCFluxClass> BCFluxPtr(
        new EllipticConstNeumBCFluxClass(IntVect::Zero, IntVect::Zero)
    );
    holder.addBCMethod(BCFluxPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// gradELambdaFuncBC
// Chombo uses basicGradPressureFuncBC (ie, 2nd order extrap)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradELambdaFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// lambdaRiemannBC
// Sets BCs on FC lambda in the Riemann solver.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::lambdaRiemannBC () const
{
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Unit,
                                          RealVect::Unit)
    );

    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// scalarSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::lambdaSlopeBC () const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ****************************** Pressure BCs *********************************

// -----------------------------------------------------------------------------
// MacPressureFuncBC
// Used in levelMacProject solver
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::MacPressureFuncBC () const
{
    bool isHomogeneous = false;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradMacPressureFuncBC
// Used to calculate Grad[phi] in MAC projection
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradMacPressureFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// LevelPressureFuncBC
// Used in LevelProject solver
// Chombo uses homog basicPressureFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::LevelPressureFuncBC () const
{
    bool isHomogeneous = false;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradPiFuncBC
// Used to calculate CCGrad[Pi] in LevelProject
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradPiFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// SyncProjFuncBC
// Used in sync projection solver
// Chombo uses inhomog basicPressureFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::SyncProjFuncBC () const
{
    bool isHomogeneous = true;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradESyncFuncBC
// Used to calculate CCGrad[eSync] in sync projection
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradESyncFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// ************************ Miscellaneous BC functions *************************

// -----------------------------------------------------------------------------
// smoothingSolverBC (used for post-redgrid smoothing)
// Chombo uses 0 Diri at solid walls
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::smoothingSolverBC () const
{
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );

    BCMethodHolder holder;
    holder.addBCMethod(BCPtr);

    return holder;
}

// -----------------------------------------------------------------------------
// BCs for the streamfunction solver.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::streamSolverBC (int comp) const
{
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero)
    );

    BCMethodHolder holder;
    holder.addBCMethod(BCPtr);

    return holder;
}



// ************************* The basic BC functions ****************************
// ****** These create BCMethodHolders from RefCountedPtr<BC*Class>es. *********

// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> velBCPtr = RefCountedPtr<BCGhostClass>(
        new BasicVelocityBCGhostClass(0.0,      //s_inflowVel,
                                      -1,       //s_inflowDir,
                                      Side::Lo, //s_inflowSide,
                                      -1,       //s_outflowDir,
                                      Side::Lo, //s_outflowSide,
                                      a_veldir,
                                      a_isViscous)
    );
    holder.addBCMethod(velBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicScalarFuncBC   (Extrapolate BCs)
// Sets physical BCs on a generic passive scalar.
// Chombo uses 1st order extrap
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicScalarFuncBC () const
{
    BCMethodHolder holder;

    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicPressureFuncBC   (Zero Neum BCs)
// Sets physical BCs on pressures (used by the Poisson solvers).
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicPressureFuncBC (bool a_isHomogeneous) const
{
    BCMethodHolder holder;

    // This sets ghosts so that Grad[CCstate] = Grad[pressure] = 0 at bdry.
    RefCountedPtr<BCGhostClass> ghostBCPtr(
        new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                          RealVect::Zero)
    );
    holder.addBCMethod(ghostBCPtr);

    // This sets face values so that FCstate = Grad[pressure] = 0 at bdry.
    RefCountedPtr<BCFluxClass> fluxBCPtr(
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero)
    );
    holder.addBCMethod(fluxBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicGradPressureFuncBC   (Extrap BCs)
// Sets physical BCs on pressures before taking gradients.
// NOTE: We don't use Neum BCs because many papers say to use Extrap BCs.
// In fact, changing these to extrap cleared up an error in the velocity field!
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicGradPressureFuncBC () const
{
    BCMethodHolder holder;

    int extrapOrder = 2;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}
