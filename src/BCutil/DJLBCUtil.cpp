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
#include "CubicSpline.H"
#include "Constants.H"
#include "Debug.H"


RealVect DJLBCUtil::s_L = RealVect::Zero;
Real DJLBCUtil::s_d = 0.1;
Real DJLBCUtil::s_z0 = 0.8;


// static const Real m = 1.0 / 6.4; // Envelope slope
// static const Real sigma = 192.0; // Envelope width
// static const Real offsetx = 0.0;
// static const Real offsety = 128.0;
// static const Real rotAngle = 0.0 * Pi/180.0;
// static const Real sinA = sin(rotAngle);
// static const Real cosA = cos(rotAngle);
// static const Real fileScale = 1;

// Oblique problem
static const Real m = 1.0 / 6.4; // Envelope slope
static const Real sigma = 192.0; // Envelope width
static const Real offsetx = 128.0;
static const Real offsety = 128.0;
static const Real rotAngle = 45.0 * Pi/180.0;
static const Real sinA = sin(rotAngle);
static const Real cosA = cos(rotAngle);
static const Real fileScale = 4;


// -----------------------------------------------------------------------------
// Reads a_eta from the DJLIC_[a_nx]x[a_nz].bin file.
// a_eta[i] is a Vector<Real> containing eta(x) at z[i].
// (a_nx, a_nz) are the number of cell centers in the domain, not nodes!
// Returns c.
// -----------------------------------------------------------------------------
Real readDJLICFile (Vector<Vector<Real> >& a_eta,
                    const int              a_nx,
                    const int              a_nz)
{
    char infileName[100];
    sprintf(infileName, "DJLIC_%dx%d.bin", a_nx, a_nz);
    pout() << "infileName = " << infileName << endl;
    std::ifstream infile;
    infile.open(infileName, ios::in | ios::binary);

    if (!infile.is_open()) {
        std::ostringstream errmsg;
        errmsg << "Could not open " << infileName;
        MayDay::Error(errmsg.str().c_str());
    }

    infile.seekg(0, ios::beg);

    // nmax
    double nmax = 0.0;
    infile.seekg(4, ios::cur);
    infile.read((char*)&nmax, sizeof(double));
    infile.seekg(4, ios::cur);

    // c
    double c = 0.0;
    infile.seekg(4, ios::cur);
    infile.read((char*)&c, sizeof(double));
    pout() << "c = " << c << endl;
    infile.seekg(4, ios::cur);

    // x
    Vector<double> x(a_nx+1, 0.0);
    infile.seekg(4, ios::cur);
    infile.read((char*)&x[0], sizeof(double)*x.size());
    infile.seekg(4, ios::cur);

    // z
    Vector<double> z(a_nz+1, 0.0);
    infile.seekg(4, ios::cur);
    infile.read((char*)&z[0], sizeof(double)*z.size());
    infile.seekg(4, ios::cur);

    // eta
    CH_assert(a_eta.size() >= a_nz+1);
    for (int k = 0; k <= a_nz; ++k) {
        Vector<double> dblVec(a_nx+1, 0.0);
        infile.seekg(4, ios::cur);
        infile.read((char*)&dblVec[0], sizeof(double)*(a_nx+1));
        infile.seekg(4, ios::cur);

        CH_assert(a_eta[k].size() >= a_nx+1);
        for (int i = 0; i <= a_nx; ++i) {
            a_eta[k][i] = dblVec[i];
        }
    }

    infile.close();

    return ((Real)c);
}


// -----------------------------------------------------------------------------
// Extrude state variable in spanwise dir.
// -----------------------------------------------------------------------------
void envelopeExtrusionVel (FArrayBox&       a_destFAB,
                           const int        a_destComp,
                           const Box&       a_valid,
                           const Box&       a_domBox,
                           const FArrayBox& a_flatSrcFAB) // comp assumed to be zero
{
    const Real m = 40.0;    // Envelope slope
    const Real sigma = 0.75; // Envelope width

    // Loop over the flat source region.
    IntVect cc;
    for (cc[2] = a_valid.smallEnd(2); cc[2] <= a_valid.bigEnd(2); ++cc[2]) {
        for (cc[0] = a_valid.smallEnd(0); cc[0] <= a_valid.bigEnd(0); ++cc[0]) {
            cc[1] = 0;
            const Real val = a_flatSrcFAB(cc,0);
            const Real jcenter = 0.5 * Real(a_domBox.bigEnd(1) + a_domBox.smallEnd(0));
            const Real jsize = Real(a_domBox.size(1));
            Real jfrac; // Range = -0.5 to 0.5
            Real amp;   // Range = 0.0 to 1.0

            for (cc[1] = a_valid.smallEnd(1); cc[1] <= a_valid.bigEnd(1); ++cc[1]) {
                jfrac = (Real(cc[1]) - jcenter) / jsize;
                amp = 0.5*(tanh(m*(jfrac+0.5*sigma))-tanh(m*(jfrac-0.5*sigma)));
                a_destFAB(cc,a_destComp) = amp*val;
            }
        }
    }
}
void envelopeExtrusionScal (FArrayBox&       a_destFAB,
                            const int        a_destComp,
                            const Box&       a_valid,
                            const Box&       a_domBox,
                            const FArrayBox& a_flatSrcFAB,        // comp assumed to be zero
                            const FArrayBox& a_flatBackGroundFAB) // comp assumed to be zero
{
    const Real m = 40.0;    // Envelope slope
    const Real sigma = 0.75; // Envelope width

    // Loop over the flat source region.
    IntVect cc;
    for (cc[2] = a_valid.smallEnd(2); cc[2] <= a_valid.bigEnd(2); ++cc[2]) {
        for (cc[0] = a_valid.smallEnd(0); cc[0] <= a_valid.bigEnd(0); ++cc[0]) {
            cc[1] = 0;
            const Real val = a_flatSrcFAB(cc,0);
            const Real bgval = a_flatBackGroundFAB(cc,0);

            const Real jcenter = 0.5 * Real(a_domBox.bigEnd(1) + a_domBox.smallEnd(0));
            const Real jsize = Real(a_domBox.size(1));
            Real jfrac; // Range = -0.5 to 0.5
            Real amp;   // Range = 0.0 to 1.0

            for (cc[1] = a_valid.smallEnd(1); cc[1] <= a_valid.bigEnd(1); ++cc[1]) {
                jfrac = (Real(cc[1]) - jcenter) / jsize;
                amp = 0.5*(tanh(m*(jfrac+0.5*sigma))-tanh(m*(jfrac-0.5*sigma)));
                a_destFAB(cc,a_destComp) = amp*val + (1.0-amp)*bgval;
            }
        }
    }
}




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

    // // Sanity checks
    // CH_assert(SpaceDim == 3);
    // CH_assert(a_velFAB.nComp() == SpaceDim);
    // CH_assert(0 <= a_velComp);
    // CH_assert(a_velComp < SpaceDim);
    // CH_assert(a_velFAB.box().type() == IntVect::Zero);

    // // y-vel is just 0
    // if (a_velComp == 1) {
    //     a_velFAB.setVal(0.0, a_velComp);
    //     return;
    // }


    // // Gather domain data
    // const ProblemDomain& domain = a_levGeo.getDomain();
    // const Box domBox = domain.domainBox();
    // const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    // const RealVect physDx = a_levGeo.getDx();
    // const IntVect& Nx = domBox.size();

    // const Box flatDomBox = flattenBox(domBox, 1);
    // const Box flatValid = flattenBox(valid, 1);

    // // Open the source data file
    // char infileName[100];
    // sprintf(infileName, "DJLIC_%dx%d.bin", Nx[0], Nx[SpaceDim-1]);
    // pout() << "infileName = " << infileName << endl;
    // std::ifstream infile;
    // infile.open(infileName, ios::in | ios::binary);

    // if (infile.is_open()) {
    //     double nmax = 0.0;
    //     double c = 0.0;
    //     Vector<double> x(Nx[0]+1, 0.0);
    //     Vector<double> z(Nx[SpaceDim-1]+1, 0.0);
    //     FArrayBox etaFAB(surroundingNodes(flatDomBox), 1);

    //     infile.seekg(0, ios::beg);

    //     infile.seekg(4, ios::cur);
    //     infile.read((char*)&nmax, sizeof(double));
    //     pout() << "nmax = " << nmax << endl;
    //     infile.seekg(4, ios::cur);

    //     infile.seekg(4, ios::cur);
    //     infile.read((char*)&c, sizeof(double));
    //     pout() << "c = " << c << endl;
    //     infile.seekg(4, ios::cur);

    //     infile.seekg(4, ios::cur);
    //     infile.read((char*)&x[0], sizeof(double)*x.size());
    //     infile.seekg(4, ios::cur);

    //     infile.seekg(4, ios::cur);
    //     infile.read((char*)&z[0], sizeof(double)*z.size());
    //     infile.seekg(4, ios::cur);

    //     const IntVect etaShift = etaFAB.box().smallEnd();
    //     etaFAB.shift(-etaShift);
    //     CH_assert(etaFAB.box().smallEnd() == IntVect::Zero);
    //     for (int k = 0; k < Nx[SpaceDim-1]+1; ++k) {
    //         Vector<double> dataVec(Nx[0]+1, 0.0);

    //         infile.seekg(4, ios::cur);
    //         infile.read((char*)&dataVec[0], sizeof(double)*(Nx[0]+1));
    //         infile.seekg(4, ios::cur);

    //         for (int i = 0; i < Nx[0]+1; ++i) {
    //             IntVect nc(i,0,k);
    //             etaFAB(nc) = dataVec[i];
    //         }
    //     }
    //     etaFAB.shift(etaShift);

    //     infile.close();

    //     // Compute velocity in 2D slice.
    //     FArrayBox flatVelFAB(flatValid, 1);

    //     if (a_velComp == 0) {
    //         // u = c * eta_z
    //         BoxIterator bit(flatValid);
    //         for (bit.reset(); bit.ok(); ++bit) {
    //             const IntVect& cc = bit();

    //             const IntVect nclb = cc;
    //             const IntVect ncrb = cc + BASISV(SpaceDim-1);

    //             const IntVect nclt = nclb + BASISV(0);
    //             const IntVect ncrt = ncrb + BASISV(0);

    //             Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[SpaceDim-1];
    //             Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[SpaceDim-1];
    //             flatVelFAB(cc,0) = 0.5 * (detat + detab);
    //         }

    //     } else if (a_velComp == 1) {
    //         MayDay::Error("Handle this separately.");

    //     } else {
    //         // w = -c * eta_x
    //         BoxIterator bit(flatValid);
    //         for (bit.reset(); bit.ok(); ++bit) {
    //             const IntVect& cc = bit();

    //             const IntVect nclb = cc;
    //             const IntVect ncrb = cc + BASISV(0);

    //             const IntVect nclt = nclb + BASISV(SpaceDim-1);
    //             const IntVect ncrt = ncrb + BASISV(SpaceDim-1);

    //             Real detat = (etaFAB(ncrt) - etaFAB(nclt)) / physDx[0];
    //             Real detab = (etaFAB(ncrb) - etaFAB(nclb)) / physDx[0];
    //             flatVelFAB(cc,0) = -0.5 * (detat + detab);
    //         }
    //     }

    //     // Extrude velocity to spanwise dir.
    //     envelopeExtrusionVel(a_velFAB, a_velComp, valid, domBox, flatVelFAB);

    // } else {
    //     std::ostringstream errmsg;
    //     errmsg << "Could not open " << infileName;
    //     MayDay::Error(errmsg.str().c_str());
    // }



    // Sanity checks
    CH_assert(SpaceDim == 3);
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);
    CH_assert(a_velFAB.box().type() == IntVect::Zero);

    // Gather domain data
    const Box domBox = a_levGeo.getDomain().domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();
    const RealVect L = a_levGeo.getDomainLength();

    // Read eta from file.
    Vector<Vector<Real> > eta(Nx[2]+1, Vector<Real>(Nx[0]+1, 0.0));
    readDJLICFile(eta, Nx[0]/fileScale, Nx[2]);

    // Compute locations of cell centers.
    Vector<Real> x(Nx[0]);
    for (int i = 0; i < Nx[0]; ++i) {
        x[i] = (Real(domBox.smallEnd(0) + i) + 0.5) * physDx[0];
    }
    const Real minX = x[0];
    const Real maxX = x[Nx[0]-1];


    // Loop over horizontal slices
    IntVect cc = valid.smallEnd();
    int k = cc[2] - domBox.smallEnd(2);
    for (; cc[2] <= valid.bigEnd(2); ++cc[2], ++k) {

        // Compute DJL velocity on this slice at cell centers...
        // u = c * eta_z
        Vector<Real> uDJL(Nx[0], 0.0);
        for (int i = 0; i < Nx[0]; ++i) {
            Real detar = eta[k+1][i+1] - eta[k][i+1];
            Real detal = eta[k+1][i  ] - eta[k][i  ];
            uDJL[i] = 0.5 * (detar + detal) / physDx[2];
        }
        // w = -c * eta_x
        Vector<Real> wDJL(Nx[0], 0.0);
        for (int i = 0; i < Nx[0]; ++i) {
            Real detat = eta[k+1][i+1] - eta[k+1][i];
            Real detab = eta[k  ][i+1] - eta[k  ][i];
            wDJL[i] = -0.5 * (detat + detab) / physDx[0];
        }

        // Construct splines of DJL velocity
        CubicSpline uDJLSpline, wDJLSpline;
        uDJLSpline.solve(uDJL, x);
        wDJLSpline.solve(wDJL, x);

        // Loop over this horizontal slice
        cc[0] = valid.smallEnd(0);
        int i = cc[0] - domBox.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0], ++i) {

            cc[1] = valid.smallEnd(1);
            int j = cc[1] - domBox.smallEnd(1);
            for (; cc[1] <= valid.bigEnd(1); ++cc[1], ++j) {

                // Compute this cell's location.
                Real thisX = (Real(cc[0]) + 0.5) * physDx[0] - offsetx;
                Real thisY = (Real(cc[1]) + 0.5) * physDx[1] - offsety;

                // Compute the location of the "source" DJL solution
                // and our distance away from the center of extrusion.
                Real srcX = thisY*sinA + thisX*cosA;
                Real dist = thisY*cosA - thisX*sinA;

                // Interpolate the DJL solution at the source location.
                Real srcU = ((minX <= srcX && srcX <= maxX)? uDJLSpline.interp(srcX): 0.0);
                Real srcW = ((minX <= srcX && srcX <= maxX)? wDJLSpline.interp(srcX): 0.0);

                // Now, rotate the source solution into position
                Real rotVel;
                if (a_velComp == 0) {
                    rotVel = srcU*cosA;
                } else if (a_velComp == 1) {
                    rotVel = -srcU*sinA;
                } else {
                    rotVel = srcW;
                }

                // Compute envelope at this distance from the center of extrusion.
                Real envelope = 0.5*(tanh(m*(dist+0.5*sigma))-tanh(m*(dist-0.5*sigma)));

                // Set 3D field.
                a_velFAB(cc, a_velComp) = envelope * rotVel;

            } // end loop over y (cc[1] and j)
        } // end loop over x (cc[0] and i)
    } // end loop over z (cc[2] and k)

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

#if 0
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
            infile.seekg(4, ios::cur);

            infile.seekg(4, ios::cur);
            infile.read((char*)&z[0], sizeof(double)*z.size());
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

            // We are done reading data from file.
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
            }

            // Extrude scalar in spanwise dir.
            envelopeExtrusionScal(a_scalarFAB, a_scalarComp, valid, domBox, flatScalarFAB, flatBackGroundFAB);

        } else {
            std::ostringstream errmsg;
            errmsg << "Could not open " << infileName;
            MayDay::Error(errmsg.str().c_str());
        }
    } else {
        MayDay::Error("scalar IC not defined for comp > 0");
    }

#else
    // Sanity checks
    CH_assert(SpaceDim == 3);
    CH_assert(a_scalarFAB.nComp() == 1);
    CH_assert(a_scalarComp == 0);
    CH_assert(a_scalarFAB.box().type() == IntVect::Zero);

    // Gather domain data
    const Box domBox = a_levGeo.getDomain().domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_scalarFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();
    const RealVect L = a_levGeo.getDomainLength();

    // Read eta from file.
    Vector<Vector<Real> > eta(Nx[2]+1, Vector<Real>(Nx[0]+1, 0.0));
    const Real c = readDJLICFile(eta, Nx[0]/fileScale, Nx[2]);

    // Compute locations of cell centers.
    Vector<Real> x(Nx[0]);
    for (int i = 0; i < Nx[0]; ++i) {
        x[i] = (Real(domBox.smallEnd(0) + i) + 0.5) * physDx[0];
    }
    const Real minX = x[0];
    const Real maxX = x[Nx[0]-1];


    // Loop over horizontal slices
    IntVect cc = valid.smallEnd();
    int k = cc[2] - domBox.smallEnd(2);
    for (; cc[2] <= valid.bigEnd(2); ++cc[2], ++k) {

        // Compute the background scalar for this slice.
        Real thisZ = (Real(cc[2]) + 0.5) * physDx[2];
        Real rho_bottom = 0.5 * (1.0 - tanh((0.0   - s_z0) / s_d));
        Real rho        = 0.5 * (1.0 - tanh((thisZ - s_z0) / s_d));
        Real rho_top    = 0.5 * (1.0 - tanh((1.0   - s_z0) / s_d));
        Real bgScalar   = (rho - rho_top) / (rho_bottom - rho_top);

        // Compute the DJL scalar for this slice.
        Vector<Real> bDJL(Nx[0], 0.0);
        for (int i = 0; i < Nx[0]; ++i) {
            Real ccEta = 0.25 * (eta[k][i] + eta[k+1][i] + eta[k][i+1] + eta[k+1][i+1]);
            Real z = thisZ - ccEta / c;
            rho = 0.5 * (1.0 - tanh((z   - s_z0) / s_d));
            bDJL[i] = (rho - rho_top) / (rho_bottom - rho_top);
        }

        // Construct splines of DJL buoyancy
        CubicSpline bDJLSpline;
        bDJLSpline.solve(bDJL, x);

        // Loop over this horizontal slice
        cc[0] = valid.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0]) {

            cc[1] = valid.smallEnd(1);
            for (; cc[1] <= valid.bigEnd(1); ++cc[1]) {

                // Compute this cell's location.
                Real thisX = (Real(cc[0]) + 0.5) * physDx[0] - offsetx;
                Real thisY = (Real(cc[1]) + 0.5) * physDx[1] - offsety;

                // Compute the location of the "source" DJL solution
                // and our distance away from the center of extrusion.
                Real srcX = thisY*sinA + thisX*cosA;
                Real dist = thisY*cosA - thisX*sinA;

                // Interpolate the DJL solution at the source location.
                Real srcB = ((minX <= srcX && srcX <= maxX)? bDJLSpline.interp(srcX): bgScalar);

                // Compute envelope at this distance from the center of extrusion.
                Real envelope = 0.5*(tanh(m*(dist+0.5*sigma))-tanh(m*(dist-0.5*sigma)));

                // Set 3D field.
                a_scalarFAB(cc,0) = envelope*srcB + (1.0-envelope)*bgScalar;

            } // end loop over y (cc[1])
        } // end loop over x (cc[0])
    } // end loop over z (cc[2] and k)

#endif // old vs new code

#endif // 2 or 3 dims
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
// // basicScalarFuncBC
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
