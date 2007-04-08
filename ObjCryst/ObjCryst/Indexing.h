/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2006- Vincent Favre-Nicolin vincefn@users.sourceforge.net

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*
*  source file for Indexing classes & functions
*
*/
#ifndef _OBJCRYST_INDEXING_H_
#define _OBJCRYST_INDEXING_H_

#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <time.h>
#include "RefinableObj/RefinableObj.h"
#include "RefinableObj/LSQNumObj.h"
namespace ObjCryst
{
/** Different lattice types.
*
*/
enum CrystalSystem
{ TRICLINIC, MONOCLINIC, ORTHOROMBIC, HEXAGONAL, RHOMBOEDRAL, TETRAGONAL, CUBIC};

enum CrystalCentering
{ LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F};

/** Estimate volume from number of peaks at a given dmin
* See J. Appl. Cryst. 20 (1987), 161
*
* \param dmin,dmax: 1/d limits between which the number of reflections has been observed
* \param nbrefl: number of observed reflections
* \param kappa: estimated percentage of reflections actually observed (default=1.0)
*/
float EstimateCellVolume(const float dmin, const float dmax, const float nbrefl, 
                         const CrystalSystem system,const CrystalCentering centering,const float kappa=1);

/** Lightweight class describing the reciprocal unit cell, for the fast computation of d*_hkl^2.
*
*
*/
class RecUnitCell
{
   public:
      RecUnitCell(const float zero=0,const float par0=0,const float par1=0,const float par2=0,
                  const float par3=0,const float par4=0,const float par5=0,CrystalSystem lattice=CUBIC);
      RecUnitCell(const RecUnitCell &old);
      void operator=(const RecUnitCell &rhs);
      // access to ith parameter
      //float operator[](const unsigned int i);
      /// Compute d*^2 for hkl reflection
      /// if deriv != -1, compute derivate versus the corresponding parameter
      ///
      /// If derivhkl=1,2,3, compute derivative versus h,k or l.
      float hkl2d(const float h,const float k,const float l,REAL *derivpar=NULL,const unsigned int derivhkl=0) const;
      /// Compute d*^2 for one hkl reflection: this functions computes a d*^2 range (min,max)
      /// for a given range of unit cell parameter (given in the delta parameter) around the
      /// current parameters.
      ///
      /// Used for DicVol algorithm
      void hkl2d_delta(const float h,const float k,const float l,const RecUnitCell &delta, float & dmin, float &dmax) const;
      std::vector<float> DirectUnitCell()const;
      /** The 6 parameters defining 1/d_hkl^2 = d*_hkl^2, for different crystal classes, from:
      * d*_hkl^2 = zero + a*^2 h^2 + b*^2 k^2 + c*^2 l^2 + 2 a*.b* hk + 2 b*.c* kl + 2 a*.c* hl
      *
      * for triclinic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2 + par[0]*par[1]*par[3] hk + par[1]*par[2]*par[4] kl + par[0]*par[2]*par[5] hl
      * for monoclinic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2 + par[0]*par[2]*par[3] hl
      * for orthorombic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2
      * for hexagonal:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[0]^2 k^2 + par[2]^2 l^2 + sqrt(3)/2*par[0]^2 hk
      * for rhomboedral:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[1]^2 k^2 + par[2]^2 l^2 + par[3] (hk + kl + hl)
      * for quadratic:
      *   d*_hkl^2 = zero + par[0]^2 h^2 + par[0]^2 k^2 + par[1]^2 l^2
      * for cubic
      *   d*_hkl^2 = zero + par[0]^2 (h^2 + k^2 + l^2)
      */
      float par[6];
      float zero;
      CrystalSystem mlattice;
};

/** Class to store positions of observed reflections.
*
*
*/
class PeakList
{
   public:
      PeakList();
      PeakList(const PeakList &old);
      void operator=(const PeakList &rhs);
      ~PeakList();
      void ImportDhklIntensity(std::istream &is);
      void ImportDhkl(std::istream &is);
      void Import2ThetaIntensity(std::istream &is, const float wavelength=1.5418);
      /// Add one peak
      ///\param d: 1/d for this peak (Angstroem)
      void AddPeak(const float d, const float iobs=1.0,const float dobssigma=0.0,const float iobssigma=0.0,
                   const int h=0,const int k=0, const int l=0,const float d2calc=0);
      void RemovePeak(unsigned int i);
      void Print(std::ostream &os) const;
      struct hkl
      {
         hkl(const float dobs=1.0,const float iobs=0.0,const float dobssigma=0.0,const float iobssigma=0.0,
             const int h=0,const int k=0, const int l=0,const float d2calc=0);
         /// Observed peak position 1/d
         float dobs;
         /// Error on peak position
         float dobssigma;
         /// Observed peak position 1/d^2
         float d2obs;
         /// Min value for observed peak position 1/(d+disgma/2)^2
         float d2obsmin;
         /// Min value for observed peak position 1/(d-disgma/2)^2
         float d2obsmax;
         /// Observed peak intensity
         float iobs;
         /// Error on observed peak intensity
         float iobssigma;
         /// Miller indices, after line is indexed
         mutable int h,k,l;
         /// Is this line indexed ?
         mutable bool isIndexed;
         /// Is this an impurity line ?
         mutable bool isSpurious;
         /// Indexing statistics
         mutable unsigned long stats;
         /// Calculated position, 1/d^2
         mutable float d2calc;
         /// 1/d^2 difference, obs-calc
         mutable float d2diff;
      };
      /// Get peak list
      vector<hkl> & GetPeakList();
      /// Get peak list
      const vector<hkl> & GetPeakList()const ;
      /// Predict peak positions
      /// Best h,k,l for each observed peak (for least-squares refinement)
      /// This is stored by the Score function, optionnally.
      mutable vector<hkl>mvHKL;
      /// Full list of calculated HKL positions for a given solution, up to a given resolution
      /// After finding a candidate solution, use score with pPredictedHKL=&mvPredictedHKL
      mutable list<hkl> mvPredictedHKL;
};

/// Compute score for a candidate RecUnitCell and a PeakList
float Score(const PeakList &dhkl, const RecUnitCell &ruc, const unsigned int nbSpurious=0,
            const bool verbose=false,const bool storehkl=false,
            const bool storePredictedHKL=false);

/** Algorithm class to find the correct indexing from observed peak positions.
*
*/
class CellExplorer:public RefinableObj
{
   public:
      CellExplorer(const PeakList &dhkl, const CrystalSystem lattice, const unsigned int nbSpurious);
      void Evolution(unsigned int ng,const bool randomize=true,const float f=0.7,const float cr=0.5,unsigned int np=100);
      void SetLengthMinMax(const float min,const float max);
      void SetAngleMinMax(const float min,const float max);
      void SetVolumeMinMax(const float min,const float max);
      void SetNbSpurious(const unsigned int nb);
      /// Allowed error on 1/d (squared!), used for dicvol
      void SetD2Error(const float err);
      void SetMinMaxZeroShift(const float min,const float max);
      virtual const string& GetClassName() const;
      virtual const string& GetName() const;
      virtual void Print() const;
      virtual unsigned int GetNbLSQFunction() const;
      virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQObs(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQWeight(const unsigned int) const;
      virtual const CrystVector_REAL & GetLSQDeriv(const unsigned int, RefinablePar &);
      virtual void BeginOptimization(const bool allowApproximations=false, const bool enableRestraints=false);
      void LSQRefine(int nbCycle=1, bool useLevenbergMarquardt=true, const bool silent=false);
      /// Run DicVOl algorithm, store only solutions with score >minScore or depth>=minDepth, 
      /// stop at the end of one volume interval (~400 A^3) if best score>stopOnScore, 
      /// or if one solution was found at depth>=stopOnDepth
      ///
      /// If stopOnDepth==0, do not stop for any depth
      void DicVol(const float minScore=10,const unsigned int minDepth=3,const float stopOnScore=50.0,const unsigned int stopOnDepth=6);
      /// Sort all solutions by score, remove duplicates
      void ReduceSolutions();
      float GetBestScore()const;
      const std::list<std::pair<RecUnitCell,float> >& GetSolutions()const;
      std::list<std::pair<RecUnitCell,float> >& GetSolutions();
   private:
      unsigned int RDicVol(RecUnitCell uc0, RecUnitCell uc1, unsigned int depth,unsigned long &nbCalc,const float minV,const float maxV);
      void Init();
      /// Max number of obs reflections to use
      std::list<std::pair<RecUnitCell,float> > mvSolution;
      unsigned int mnpar;
      const PeakList *mpPeakList;
      float mLengthMin,mLengthMax;
      float mAngleMin,mAngleMax;
      float mVolumeMin,mVolumeMax;
      float mZeroShiftMin,mZeroShiftMax;
      /// Min values for all parameters (7=unit cell +zero)
      float mMin[7];
      /// Max amplitude (max=min+amplitude) for all parameters
      /// All parameters are treated as periodic for DE (??)
      float mAmp[7];
      /// Lattice type for which we search
      CrystalSystem mlattice;
      unsigned int mNbSpurious;
      float mD2Error;
      LSQNumObj mLSQObj;
      mutable CrystVector_REAL mObs;
      mutable CrystVector_REAL mCalc;
      mutable CrystVector_REAL mWeight;
      mutable CrystVector_REAL mDeriv;
      /// Reciprocal unit cell used for least squares refinement
      RecUnitCell mRecUnitCell;
      /// Current best score
      float mBestScore;
      /// Number of solutions found during dicvol search, at each depth.
      std::vector<unsigned int> mvNbSolutionDepth;
      float mMinScoreReport;
      unsigned int mMaxDicVolDepth,mDicVolDepthReport;
};


}//namespace
#endif
