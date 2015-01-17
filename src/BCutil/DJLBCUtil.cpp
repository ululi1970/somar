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
#include "DJLBCUtil.H"
#include "ProblemContext.H"
#include "BoxIterator.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ConvertFABF_F.H"
#include "EllipticBCUtils.H"

#include "AMRIO.H"

RealVect DJLBCUtil::s_L = RealVect::Zero;
Real DJLBCUtil::s_d = 0.1;
Real DJLBCUtil::s_z0 = 0.8;


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
DJLBCUtil::DJLBCUtil ()
{
    static bool paramsRead = false;
    if (!paramsRead) {
        const ProblemContext* ctx = ProblemContext::getInstance();

        s_L = ctx->domainLength;
        // s_d, s_z0

        paramsRead = true;
    }
}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
DJLBCUtil::~DJLBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* DJLBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new DJLBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void DJLBCUtil::setVelIC (FArrayBox&           a_velFAB,
                          const int            a_velComp,
                          const LevelGeometry& a_levGeo,
                          const DataIndex&     a_di) const
{
#if CH_SPACEDIM == 2
    // Sanity checks
    CH_assert(SpaceDim == 2); // Streamfunction method is different for other dims
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);
    CH_assert(a_velFAB.box().type() == IntVect::Zero);

    // Gather domain data
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box domBox = domain.domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();

    // Open input data file
    char infileName[100];
    sprintf(infileName, "DJLIC_%dx%d.bin", Nx[0], Nx[SpaceDim-1]);
    std::ifstream infile;
    infile.open(infileName, ios::in | ios::binary);

    // Gather data from file
    if (infile.is_open()) {
        double nmax = 0.0;
        double c = 0.0;
        Vector<double> x(Nx[0]+1, 0.0);
        Vector<double> z(Nx[SpaceDim-1]+1, 0.0);
        FArrayBox etaFAB(surroundingNodes(domBox), 1);

        // Move to beginning of file.
        infile.seekg(0, ios::beg);

        // Read N^2 scaling. (This is not used)
        infile.seekg(4, ios::cur);
        infile.read((char*)&nmax, sizeof(double));
        infile.seekg(4, ios::cur);

        // Read c (long-wave speed)
        infile.seekg(4, ios::cur);
        infile.read((char*)&c, sizeof(double));
        infile.seekg(4, ios::cur);

        // Read x coordinates
        infile.seekg(4, ios::cur);
        infile.read((char*)&x[0], sizeof(double)*x.size());
        infile.seekg(4, ios::cur);

        // Read z coordinates
        infile.seekg(4, ios::cur);
        infile.read((char*)&z[0], sizeof(double)*z.size());
        infile.seekg(4, ios::cur);

        // Read eta
        const IntVect etaShift = etaFAB.box().smallEnd();
        etaFAB.shift(-etaShift);
        CH_assert(etaFAB.box().smallEnd() == IntVect::Zero);
        for (int k = 0; k < Nx[SpaceDim-1]+1; ++k) {
            Vector<double> dataVec(Nx[0]+1, 0.0);

            infile.seekg(4, ios::cur);
            infile.read((char*)&dataVec[0], sizeof(double)*(Nx[0]+1));
            infile.seekg(4, ios::cur);

            for (int i = 0; i < Nx[0]+1; ++i) {
                IntVect nc(D_DECL(i,k,0));
                etaFAB(nc) = dataVec[i];
            }
        }
        etaFAB.shift(etaShift);

        // We are done reading data from file.
        infile.close();

        // Construct the velocity field.
        if (a_velComp == 0) {
            // u = c * eta_z
            BoxIterator bit(valid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                const IntVect nclb = cc;
                const IntVect ncrb = cc + BASISV(SpaceDim-1);

                const IntVect nclt = nclb + BASISV(0);
                const IntVect ncrt = ncrb + BASISV(0);

                Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[SpaceDim-1];
                Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[SpaceDim-1];
                a_velFAB(cc,a_velComp) = 0.5 * (detat + detab);
            }
        } else {
            // w = -c * eta_x
            BoxIterator bit(valid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                const IntVect nclb = cc;
                const IntVect ncrb = cc + BASISV(0);

                const IntVect nclt = nclb + BASISV(SpaceDim-1);
                const IntVect ncrt = ncrb + BASISV(SpaceDim-1);

                Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[0];
                Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[0];
                a_velFAB(cc,a_velComp) = -0.5 * (detat + detab);
            }
        }

    } else {
        std::ostringstream errmsg;
        errmsg << "Could not open " << infileName;
        MayDay::Error(errmsg.str().c_str());
    }

#else //CH_SPACEDIM == 3

    // Sanity checks
    CH_assert(SpaceDim == 3);
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);
    CH_assert(a_velFAB.box().type() == IntVect::Zero);

    // y-vel is just 0
    if (a_velComp == 1) {
        a_velFAB.setVal(0.0, a_velComp);
        return;
    }


    // Gather domain data
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box domBox = domain.domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();

    const Box flatDomBox = flattenBox(domBox, 1);
    const Box flatValid = flattenBox(valid, 1);

    // Open the source data file
    char infileName[100];
    sprintf(infileName, "DJLIC_%dx%d.bin", Nx[0], Nx[SpaceDim-1]);
    pout() << "infileName = " << infileName << endl;
    std::ifstream infile;
    infile.open(infileName, ios::in | ios::binary);

    if (infile.is_open()) {
        double nmax = 0.0;
        double c = 0.0;
        Vector<double> x(Nx[0]+1, 0.0);
        Vector<double> z(Nx[SpaceDim-1]+1, 0.0);
        FArrayBox etaFAB(surroundingNodes(flatDomBox), 1);

        infile.seekg(0, ios::beg);

        infile.seekg(4, ios::cur);
        infile.read((char*)&nmax, sizeof(double));
        pout() << "nmax = " << nmax << endl;
        infile.seekg(4, ios::cur);

        infile.seekg(4, ios::cur);
        infile.read((char*)&c, sizeof(double));
        pout() << "c = " << c << endl;
        infile.seekg(4, ios::cur);

        infile.seekg(4, ios::cur);
        infile.read((char*)&x[0], sizeof(double)*x.size());
        // pout() << "x = " << x << endl;
        infile.seekg(4, ios::cur);

        infile.seekg(4, ios::cur);
        infile.read((char*)&z[0], sizeof(double)*z.size());
        // pout() << "z = " << z << endl;
        infile.seekg(4, ios::cur);

        const IntVect etaShift = etaFAB.box().smallEnd();
        etaFAB.shift(-etaShift);
        CH_assert(etaFAB.box().smallEnd() == IntVect::Zero);
        for (int k = 0; k < Nx[SpaceDim-1]+1; ++k) {
            Vector<double> dataVec(Nx[0]+1, 0.0);

            infile.seekg(4, ios::cur);
            infile.read((char*)&dataVec[0], sizeof(double)*(Nx[0]+1));
            infile.seekg(4, ios::cur);

            for (int i = 0; i < Nx[0]+1; ++i) {
                IntVect nc(i,0,k);
                etaFAB(nc) = dataVec[i];
            }
        }
        etaFAB.shift(etaShift);

        infile.close();

        // Compute velocity in 2D slice.
        FArrayBox flatVelFAB(flatValid, 1);

        if (a_velComp == 0) {
            // u = c * eta_z
            BoxIterator bit(flatValid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                const IntVect nclb = cc;
                const IntVect ncrb = cc + BASISV(SpaceDim-1);

                const IntVect nclt = nclb + BASISV(0);
                const IntVect ncrt = ncrb + BASISV(0);

                Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[SpaceDim-1];
                Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[SpaceDim-1];
                flatVelFAB(cc,0) = 0.5 * (detat + detab);
                // a_velFAB(cc,a_velComp) = 0.5 * (detat + detab);
            }

        } else if (a_velComp == 1) {
            MayDay::Error("Handle this separately.");

        } else {
            // w = -c * eta_x
            BoxIterator bit(flatValid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                const IntVect nclb = cc;
                const IntVect ncrb = cc + BASISV(0);

                const IntVect nclt = nclb + BASISV(SpaceDim-1);
                const IntVect ncrt = ncrb + BASISV(SpaceDim-1);

                Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[0];
                Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[0];
                flatVelFAB(cc,0) = -0.5 * (detat + detab);
                // a_velFAB(cc,a_velComp) = -0.5 * (detat + detab);
            }
        }

        // Extrude velocity to spanwise dir.
        IntVect cc;
        for (cc[2] = valid.smallEnd(2); cc[2] <= valid.bigEnd(2); ++cc[2]) {
            for (cc[0] = valid.smallEnd(0); cc[0] <= valid.bigEnd(0); ++cc[0]) {
                cc[1] = 0;
                const Real val = flatVelFAB(cc,0);

                const int jcenter = (domBox.bigEnd(1) + domBox.smallEnd(0)) / 2;
                const Real jwidth = Real(domBox.size(1)) / 10.0;
                int jdist;
                Real amp, arg;

                for (cc[1] = valid.smallEnd(1); cc[1] <= valid.bigEnd(1); ++cc[1]) {
                    jdist = cc[1] - jcenter;
                    arg = Real(jdist) / jwidth;
                    amp = exp(-0.5*arg*arg);
                    a_velFAB(cc,a_velComp) = amp*val;
                }
            }
        }

    } else {
        std::ostringstream errmsg;
        errmsg << "Could not open " << infileName;
        MayDay::Error(errmsg.str().c_str());
    }
#endif
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void DJLBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                             const int            a_scalarComp,
                             const LevelGeometry& a_levGeo,
                             const DataIndex&     a_di) const
{
#if CH_SPACEDIM == 2
    CH_assert(SpaceDim == 2); // Streamfunction method is different for other dims
    CH_assert(a_scalarFAB.nComp() == 1);

    if (a_scalarComp == 0) {
        // Gather domain data
        const ProblemDomain& domain = a_levGeo.getDomain();
        const Box domBox = domain.domainBox();
        const Box valid = a_levGeo.getBoxes()[a_di] & a_scalarFAB.box();
        const Real dz = a_levGeo.getDx()[SpaceDim-1];
        const IntVect& Nx = domBox.size();

        char infileName[100];
        sprintf(infileName, "DJLIC_%dx%d.bin", Nx[0], Nx[SpaceDim-1]);
        pout() << "infileName = " << infileName << endl;
        std::ifstream infile;
        infile.open(infileName, ios::in | ios::binary);

        if (infile.is_open()) {
            double nmax = 0.0;
            double c = 0.0;
            Vector<double> x(Nx[0]+1, 0.0);
            Vector<double> z(Nx[SpaceDim-1]+1, 0.0);
            FArrayBox etaFAB(surroundingNodes(domBox), 1);

            // Move to beginning of file.
            infile.seekg(0, ios::beg);

            // Read N^2 scaling. (This is not used)
            infile.seekg(4, ios::cur);
            infile.read((char*)&nmax, sizeof(double));
            infile.seekg(4, ios::cur);

            // Read c (long-wave speed)
            infile.seekg(4, ios::cur);
            infile.read((char*)&c, sizeof(double));
            infile.seekg(4, ios::cur);

            // Read x coordinates
            infile.seekg(4, ios::cur);
            infile.read((char*)&x[0], sizeof(double)*x.size());
            infile.seekg(4, ios::cur);

            // Read z coordinates
            infile.seekg(4, ios::cur);
            infile.read((char*)&z[0], sizeof(double)*z.size());
            infile.seekg(4, ios::cur);

            // Read eta
            const IntVect etaShift = etaFAB.box().smallEnd();
            etaFAB.shift(-etaShift);
            CH_assert(etaFAB.box().smallEnd() == IntVect::Zero);
            for (int k = 0; k < Nx[SpaceDim-1]+1; ++k) {
                Vector<double> dataVec(Nx[0]+1, 0.0);

                infile.seekg(4, ios::cur);
                infile.read((char*)&dataVec[0], sizeof(double)*(Nx[0]+1));
                infile.seekg(4, ios::cur);

                for (int i = 0; i < Nx[0]+1; ++i) {
                    IntVect nc(D_DECL(i,k,0));
                    etaFAB(nc) = dataVec[i];
                }
            }
            etaFAB.shift(etaShift);

            // We are done reading data from file.
            infile.close();

            // Convert etaFAB centering to match a_scalarFAB.
            CH_assert(a_scalarFAB.box().type() == IntVect::Zero);
            FArrayBox ccEtaFAB(valid, 1);
            FORT_CONVERTFAB(
                CHF_FRA1(ccEtaFAB,0),
                CHF_BOX(valid),
                CHF_CONST_INTVECT(IntVect::Zero),
                CHF_CONST_FRA1(etaFAB,0),
                CHF_CONST_INTVECT(IntVect::Unit));

            // Construct total buoyancy.
            BoxIterator bit(valid);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                Real z = (Real(cc[SpaceDim-1]) + 0.5) * dz - ccEtaFAB(cc,0)/c;
                Real rho_bottom = 0.5 * (1.0 - tanh((0.0 - s_z0) / s_d));
                Real rho        = 0.5 * (1.0 - tanh((z   - s_z0) / s_d));
                Real rho_top    = 0.5 * (1.0 - tanh((1.0 - s_z0) / s_d));

                a_scalarFAB(cc,0) = (rho - rho_top) / (rho_bottom - rho_top);
            }
        } else {
            std::ostringstream errmsg;
            errmsg << "Could not open " << infileName;
            MayDay::Error(errmsg.str().c_str());
        }

    } else {
        MayDay::Error("scalar IC not defined for comp > 0");
    }

#else //CH_SPACEDIM == 3

    CH_assert(SpaceDim == 3);
    CH_assert(a_scalarFAB.nComp() == 1);

    if (a_scalarComp == 0) {
        // Gather domain data
        const ProblemDomain& domain = a_levGeo.getDomain();
        const Box domBox = domain.domainBox();
        const Box valid = a_levGeo.getBoxes()[a_di] & a_scalarFAB.box();
        const Real dz = a_levGeo.getDx()[SpaceDim-1];
        const IntVect& Nx = domBox.size();

        const Box flatDomBox = flattenBox(domBox, 1);
        const Box flatValid = flattenBox(valid, 1);

        Box flatDomNodeBox = flatDomBox;
        flatDomNodeBox.surroundingNodes(0);
        flatDomNodeBox.surroundingNodes(SpaceDim-1);

        char infileName[100];
        sprintf(infileName, "DJLIC_%dx%d.bin", Nx[0], Nx[SpaceDim-1]);
        pout() << "infileName = " << infileName << endl;
        std::ifstream infile;
        infile.open(infileName, ios::in | ios::binary);

        if (infile.is_open()) {
            double nmax = 0.0;
            double c = 0.0;
            Vector<double> x(Nx[0]+1, 0.0);
            Vector<double> z(Nx[SpaceDim-1]+1, 0.0);
            FArrayBox etaFAB(flatDomNodeBox, 1);
            // FArrayBox rhoFAB(domBox, 1);

            infile.seekg(0, ios::beg);

            infile.seekg(4, ios::cur);
            infile.read((char*)&nmax, sizeof(double));
            pout() << "nmax = " << nmax << endl;
            infile.seekg(4, ios::cur);

            infile.seekg(4, ios::cur);
            infile.read((char*)&c, sizeof(double));
            pout() << "c = " << c << endl;
            infile.seekg(4, ios::cur);

            infile.seekg(4, ios::cur);
            infile.read((char*)&x[0], sizeof(double)*x.size());
            // pout() << "x = " << x << endl;
            infile.seekg(4, ios::cur);

            infile.seekg(4, ios::cur);
            infile.read((char*)&z[0], sizeof(double)*z.size());
            // pout() << "z = " << z << endl;
            infile.seekg(4, ios::cur);

            const IntVect etaShift = etaFAB.box().smallEnd();
            etaFAB.shift(-etaShift);
            CH_assert(etaFAB.box().smallEnd() == IntVect::Zero);
            for (int k = 0; k < Nx[SpaceDim-1]+1; ++k) {
                Vector<double> dataVec(Nx[0]+1, 0.0);

                infile.seekg(4, ios::cur);
                infile.read((char*)&dataVec[0], sizeof(double)*(Nx[0]+1));
                infile.seekg(4, ios::cur);

                for (int i = 0; i < Nx[0]+1; ++i) {
                    IntVect nc(i,0,k);
                    etaFAB(nc) = dataVec[i];
                }
            }
            etaFAB.shift(etaShift);
            // writeFABname(&etaFAB, "etaFAB.hdf5");

            // const IntVect rhoShift = rhoFAB.box().smallEnd();
            // rhoFAB.shift(-rhoShift);
            // CH_assert(rhoFAB.box().smallEnd() == IntVect::Zero);
            // for (int k = 0; k < Nx[SpaceDim-1]; ++k) {
            //     Vector<double> dataVec(Nx[0], 0.0);

            //     infile.seekg(4, ios::cur);
            //     infile.read((char*)&dataVec[0], sizeof(double)*Nx[0]);
            //     infile.seekg(4, ios::cur);

            //     for (int i = 0; i < Nx[0]; ++i) {
            //         IntVect cc(D_DECL(i,k,0));
            //         rhoFAB(cc) = dataVec[i];
            //     }
            // }
            // rhoFAB.shift(rhoShift);
            // writeFABname(&rhoFAB, "rhoFAB.hdf5");

            infile.close();


            // Convert etaFAB centering.
            FArrayBox ccFlatEtaFAB(flatValid, 1);
            FORT_CONVERTFAB(
                CHF_FRA1(ccFlatEtaFAB,0),
                CHF_BOX(flatValid),
                CHF_CONST_INTVECT(IntVect::Zero),
                CHF_CONST_FRA1(etaFAB,0),
                CHF_CONST_INTVECT(IntVect(1,0,1)));

            // Set total b.
            FArrayBox flatScalarFAB(flatValid, 1);
            FArrayBox flatBackGroundFAB(flatValid, 1);
            BoxIterator bit(flatValid);

            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                Real z = (Real(cc[SpaceDim-1]) + 0.5) * dz;
                Real rho_bottom = 0.5 * (1.0 - tanh((0.0 - s_z0) / s_d));
                Real rho        = 0.5 * (1.0 - tanh((z   - s_z0) / s_d));
                Real rho_top    = 0.5 * (1.0 - tanh((1.0 - s_z0) / s_d));
                flatBackGroundFAB(cc,0) = (rho - rho_top) / (rho_bottom - rho_top);

                z -= ccFlatEtaFAB(cc,0)/c;
                rho = 0.5 * (1.0 - tanh((z   - s_z0) / s_d));
                flatScalarFAB(cc,0) = (rho - rho_top) / (rho_bottom - rho_top);
                // a_scalarFAB(cc,0) = (rho - rho_top) / (rho_bottom - rho_top);

            }

            // Extrude scalar in spanwise dir.
            IntVect cc;
            for (cc[2] = valid.smallEnd(2); cc[2] <= valid.bigEnd(2); ++cc[2]) {
                for (cc[0] = valid.smallEnd(0); cc[0] <= valid.bigEnd(0); ++cc[0]) {
                    cc[1] = 0;
                    const Real val = flatScalarFAB(cc,0);
                    const Real bgval = flatBackGroundFAB(cc,0);

                    const int jcenter = (domBox.bigEnd(1) + domBox.smallEnd(0)) / 2;
                    const Real jwidth = Real(domBox.size(1)) / 10.0;
                    int jdist;
                    Real amp, arg;

                    for (cc[1] = valid.smallEnd(1); cc[1] <= valid.bigEnd(1); ++cc[1]) {
                        jdist = cc[1] - jcenter;
                        arg = Real(jdist) / jwidth;
                        amp = exp(-0.5*arg*arg);
                        a_scalarFAB(cc,a_scalarComp) = amp*val + (1.0-amp)*bgval;
                    }
                }
            }

        } else {
            std::ostringstream errmsg;
            errmsg << "Could not open " << infileName;
            MayDay::Error(errmsg.str().c_str());
        }

    } else {
        MayDay::Error("scalar IC not defined for comp > 0");
    }

#endif
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void DJLBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                     const int            a_scalarComp,
                                     const LevelGeometry& a_levGeo,
                                     const DataIndex&     a_di,
                                     const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar && a_scalarComp == 0) {
        // Gather Cartesian coordinates.
        FArrayBox posFAB(a_scalarFAB.box(), 1);
        const RealVect& dx = a_levGeo.getDx();
        a_levGeo.getGeoSourcePtr()->fill_physCoor(posFAB, 0, SpaceDim-1, dx);

        // Loop over a_scalarFAB and set background scalar values.
        // We get the values from the bscalBCValues function.
        BoxIterator bit(a_scalarFAB.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();

            Real pos[CH_SPACEDIM];
            pos[CH_SPACEDIM-1] = posFAB(iv);

            int dirDummy;
            Side::LoHiSide sideDummy;
            Real derivScaleDummy = 1.0e300;
            Real value[1];
            this->bscalBCValues(pos, &dirDummy, &sideDummy, value, derivScaleDummy, a_time);
            a_scalarFAB(iv,a_scalarComp) = value[0];
        }

    } else {
        // This scalar does not have a background state.
        a_scalarFAB.setVal(0.0, a_scalarComp);
    }
}


// -----------------------------------------------------------------------------
// Simply sets a_value to the background density at any given point, which
// is all we need for the boundary conditions.
// This function conforms to the EllipticBCValueFunc typedef.
// -----------------------------------------------------------------------------
void DJLBCUtil::bscalBCValues (Real*           a_pos,
                               int*            a_dir,
                               Side::LoHiSide* a_side,
                               Real*           a_value,
                               Real            a_derivScale,
                               Real            a_time)
{
    Real z = a_pos[CH_SPACEDIM-1];
    Real rho_bottom = 0.5 * (1.0 - tanh((0.0 - s_z0) / s_d));
    Real rho        = 0.5 * (1.0 - tanh((z   - s_z0) / s_d));
    Real rho_top    = 0.5 * (1.0 - tanh((1.0 - s_z0) / s_d));

    a_value[0] = (rho - rho_top) / (rho_bottom - rho_top);
}


// -----------------------------------------------------------------------------
// Sets the (vector scaled) target velocity for the sponge layer. By default,
// this function persuades the velocity field to approach its inflow value.
// -----------------------------------------------------------------------------
void DJLBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
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

    // Assume the boundaries are in the far-field where not much is happening.
    if (a_spongeDir < SpaceDim-1) {
        a_target.setVal(0.0);
    } else {
        MayDay::Error("DJLBCUtil::fillVelSpongeLayerTarget "
                      "can only set a sponge target when a_spongeDir < SpaceDim-1");
    }
}


// // -----------------------------------------------------------------------------
// // basicVelFuncBC
// // Sets physical BCs on velocities.
// // -----------------------------------------------------------------------------
// BCMethodHolder DJLBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
// {
//     const IntVect hUnit = IntVect::Unit - BASISV(CH_SPACEDIM-1);
//     const IntVect vUnit = BASISV(CH_SPACEDIM-1);

//     BCMethodHolder holder;

//     //             Freeslip
//     // u: Neum 0 |==========| Neum 0
//     //             Freeslip

//     // Low order extrap in horizontal (sponged) directions
//     int extrapOrder = 0;
//     RefCountedPtr<BCGhostClass> horizBCPtr(
//         new EllipticExtrapBCGhostClass(extrapOrder,
//                                        hUnit,
//                                        hUnit)
//     );
//     holder.addBCMethod(horizBCPtr);

//     RefCountedPtr<BCFluxClass> fluxBCPtr(
//         new EllipticConstNeumBCFluxClass(RealVect::Zero,
//                                          RealVect::Zero,
//                                          BASISV(0),
//                                          BASISV(0))
//     );
//     holder.addBCMethod(fluxBCPtr);

//     // Free slip in vertical dir
//     RefCountedPtr<BCGhostClass> hiVertBCPtr = RefCountedPtr<BCGhostClass>(
//         new BasicVelocityBCGhostClass(0.0,             // inflowVel
//                                       -1,              // inflowDir
//                                       Side::Lo,        // inflowSide
//                                       -1,              // outflowDir
//                                       Side::Hi,        // outflowSide
//                                       a_veldir,
//                                       false,           // isViscous
//                                       vUnit,
//                                       vUnit)
//     );
//     holder.addBCMethod(hiVertBCPtr);

//     return holder;
// }


// // -----------------------------------------------------------------------------
// // basicScalarFuncBC   (Extrapolate BCs)
// // Sets physical BCs on a generic passive scalar.
// // Chombo uses 1st order extrap
// // -----------------------------------------------------------------------------
// BCMethodHolder DJLBCUtil::basicScalarFuncBC () const
// {
//     BCMethodHolder holder;

//     RefCountedPtr<BCGhostClass> BCPtr(
//         new EllipticConstDiriBCGhostClass(RealVect::Zero,
//                                           RealVect::Zero,
//                                           IntVect::Unit,
//                                           IntVect::Unit)
//     );
//     holder.addBCMethod(BCPtr);

//     return holder;
// }


// // -----------------------------------------------------------------------------
// // basicPressureFuncBC
// // Sets physical BCs on pressures (used by the Poisson solvers).
// // -----------------------------------------------------------------------------
// BCMethodHolder DJLBCUtil::basicPressureFuncBC (bool a_isHomogeneous) const
// {
//     BCMethodHolder holder;

//     const IntVect vmask = IntVect::Unit;
//     const IntVect hmask = IntVect::Zero;

//     RefCountedPtr<BCGhostClass> diriBCPtr(
//         new EllipticConstDiriBCGhostClass(RealVect::Zero,
//                                           RealVect::Zero,
//                                           hmask,
//                                           hmask)
//     );
//     holder.addBCMethod(diriBCPtr);

//     // This sets ghosts so that Grad[CCstate] = Grad[pressure] = 0 at bdry.
//     RefCountedPtr<BCGhostClass> neumBCPtr(
//         new EllipticConstNeumBCGhostClass(RealVect::Zero,
//                                           RealVect::Zero,
//                                           vmask,
//                                           vmask)
//     );
//     holder.addBCMethod(neumBCPtr);

//     // This sets face values so that FCstate = Grad[pressure] = 0 at bdry.
//     RefCountedPtr<BCFluxClass> neumBCFluxPtr(
//         new EllipticConstNeumBCFluxClass(RealVect::Zero,
//                                          RealVect::Zero,
//                                          vmask,
//                                          vmask)
//     );
//     holder.addBCMethod(neumBCFluxPtr);

//     return holder;
// }
