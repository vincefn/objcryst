#include <ctype.h>
#include <cmath>
#include <boost/format.hpp>

#include "cctbx/sgtbx/space_group.h"
#include "cctbx/sgtbx/space_group_type.h"
#include "cctbx/miller/sym_equiv.h"
#include "cctbx/sgtbx/brick.h"

#include "ObjCryst/ObjCryst/CIF.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/Quirks/Chronometer.h"

#define POSSIBLY_UNUSED(expr) (void)(expr)

using namespace std;

namespace ObjCryst
{
CIFData::CIFAtom::CIFAtom():
mLabel(""),mSymbol(""),mOccupancy(1.0),mBiso(0.0)
{}

CIFData::CIFData()
{}

void CIFData::ExtractAll(const bool verbose)
{
   (*fpObjCrystInformUser)("CIF: Extract Data...");
   // :TODO: convert cartesian to fractional coordinates and vice-versa, if unit cell is known
   // :TODO: Take care of values listed as "." and "?" instead of a real value.
   this->ExtractName(verbose);
   this->ExtractUnitCell(verbose);
   this->ExtractSpacegroup(verbose);
   this->ExtractAtomicPositions(verbose);
   this->ExtractAnisotropicADPs(verbose);
   this->ExtractPowderPattern(verbose);
   this->ExtractSingleCrystalData(verbose);
   (*fpObjCrystInformUser)("CIF: Finished Extracting Data...");
}

void CIFData::ExtractUnitCell(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_cell_length_a");
   if(positem!=mvItem.end())
   {
      (*fpObjCrystInformUser)("CIF: Extract Unit Cell...");
      mvLatticePar.resize(6);
      mvLatticePar[0]=CIFNumeric2REAL(positem->second);
      positem=mvItem.find("_cell_length_b");
      if(positem!=mvItem.end())
         mvLatticePar[1]=CIFNumeric2REAL(positem->second);
      positem=mvItem.find("_cell_length_c");
      if(positem!=mvItem.end())
         mvLatticePar[2]=CIFNumeric2REAL(positem->second);
      positem=mvItem.find("_cell_angle_alpha");
      if(positem!=mvItem.end())
         mvLatticePar[3]=CIFNumeric2REAL(positem->second);
      positem=mvItem.find("_cell_angle_beta");
      if(positem!=mvItem.end())
         mvLatticePar[4]=CIFNumeric2REAL(positem->second);
      positem=mvItem.find("_cell_angle_gamma");
      if(positem!=mvItem.end())
         mvLatticePar[5]=CIFNumeric2REAL(positem->second);
      if(verbose) cout<<"Found Lattice parameters:" <<mvLatticePar[0]<<" , "<<mvLatticePar[1]<<" , "<<mvLatticePar[2]
                      <<" , "<<mvLatticePar[3]<<" , "<<mvLatticePar[4]<<" , "<<mvLatticePar[5]<<endl;
      mvLatticePar[3]*=0.017453292519943295;// pi/180
      mvLatticePar[4]*=0.017453292519943295;
      mvLatticePar[5]*=0.017453292519943295;
      this->CalcMatrices();
   }
}

void CIFData::ExtractSpacegroup(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_space_group_IT_number");
   if(positem!=mvItem.end())
   {
      mSpacegroupNumberIT=positem->second;//CIFNumeric2Int()
      if(verbose) cout<<"Found spacegroup IT number:"<<mSpacegroupNumberIT<<endl;
   }
   else
   {
      positem=mvItem.find("_symmetry_Int_Tables_number");
      if(positem!=mvItem.end())
      {
         mSpacegroupNumberIT=positem->second;//CIFNumeric2Int()
         if(verbose) cout<<"Found spacegroup IT number (with OBSOLETE CIF #1.0 TAG):"<<mSpacegroupNumberIT<<endl;
      }
   }

   positem=mvItem.find("_space_group_name_Hall");
   if(positem!=mvItem.end())
   {
      mSpacegroupSymbolHall=positem->second;
      if(verbose) cout<<"Found spacegroup Hall symbol:"<<mSpacegroupSymbolHall<<endl;
   }
   else
   {
      positem=mvItem.find("_symmetry_space_group_name_Hall");
      if(positem!=mvItem.end())
      {
         mSpacegroupSymbolHall=positem->second;
         if(verbose) cout<<"Found spacegroup Hall symbol (with OBSOLETE CIF #1.0 TAG):"<<mSpacegroupSymbolHall<<endl;
      }
   }

   positem=mvItem.find("_space_group_name_H-M_alt");
   if(positem!=mvItem.end())
   {
      mSpacegroupHermannMauguin=positem->second;
      if(verbose) cout<<"Found spacegroup Hermann-Mauguin symbol:"<<mSpacegroupHermannMauguin<<endl;
   }
   else
   {
      positem=mvItem.find("_symmetry_space_group_name_H-M");
      if(positem!=mvItem.end())
      {
         mSpacegroupHermannMauguin=positem->second;
         if(verbose) cout<<"Found spacegroup Hall Hermann-Mauguin (with OBSOLETE CIF #1.0 TAG):"<<mSpacegroupHermannMauguin<<endl;
      }
   }
   // Try to extract symmetry_as_xyz
   for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
       loop!=mvLoop.end();++loop)
   {
      if(mvSymmetry_equiv_pos_as_xyz.size()>0) break;// only extract ONE list of symmetry strings
      map<ci_string,vector<string> >::const_iterator pos;
      pos=loop->second.find("_symmetry_equiv_pos_as_xyz");
      if(pos!=loop->second.end())
      {
         if(verbose) cout<<"Found list of _symmetry_equiv_pos_as_xyz:"<<endl;
         for(unsigned int i=0;i<pos->second.size();++i)
         {
            if(verbose) cout<<"   "<<pos->second[i]<<endl;
            mvSymmetry_equiv_pos_as_xyz.insert(pos->second[i]);
         }
      }
   }
}

void CIFData::ExtractName(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_chemical_name_systematic");
   if(positem!=mvItem.end())
   {
      mName=positem->second;
      if(verbose) cout<<"Found chemical name:"<<mName<<endl;
   }
   else
   {
      positem=mvItem.find("_chemical_name_mineral");
      if(positem!=mvItem.end())
      {
         mName=positem->second;
         if(verbose) cout<<"Found chemical name:"<<mName<<endl;
      }
      else
      {
         positem=mvItem.find("_chemical_name_structure_type");
         if(positem!=mvItem.end())
         {
            mName=positem->second;
            if(verbose) cout<<"Found chemical name:"<<mName<<endl;
         }
         else
         {
            positem=mvItem.find("_chemical_name_common");
            if(positem!=mvItem.end())
            {
               mName=positem->second;
               if(verbose) cout<<"Found chemical name:"<<mName<<endl;
            }
            else
            {
               positem=mvItem.find("_chemical_formula_moiety");
               if(positem!=mvItem.end())
               {
                  mName=positem->second;
                  if(verbose) cout<<"Found chemical name:"<<mName<<endl;
               }
               else
               {
                  positem=mvItem.find("_chemical_formula_sum");
                  if(positem!=mvItem.end())
                  {
                     mName=positem->second;
                     if(verbose) cout<<"Found chemical name:"<<mName<<endl;
                  }
               }
            }
         }
      }
   }
   /// Crystal formula
   positem=mvItem.find("_chemical_formula_analytical");
   if(positem!=mvItem.end())
   {
      mFormula=positem->second;
      if(verbose) cout<<"Found chemical formula:"<<mFormula<<endl;
   }
   else
   {
      positem=mvItem.find("_chemical_formula_structural");
      if(positem!=mvItem.end())
      {
         mFormula=positem->second;
         if(verbose) cout<<"Found chemical formula:"<<mFormula<<endl;
      }
      else
      {
         positem=mvItem.find("_chemical_formula_iupac");
         if(positem!=mvItem.end())
         {
            mFormula=positem->second;
            if(verbose) cout<<"Found chemical formula:"<<mFormula<<endl;
         }
         else
         {
            positem=mvItem.find("_chemical_formula_moiety");
            if(positem!=mvItem.end())
            {
               mFormula=positem->second;
               if(verbose) cout<<"Found chemical formula:"<<mFormula<<endl;
            }
         }
      }
   }
}

void CIFData::ExtractAtomicPositions(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
       loop!=mvLoop.end();++loop)
   {
      if(mvAtom.size()>0) break;// only extract ONE list of atoms, preferably fractional coordinates
      map<ci_string,vector<string> >::const_iterator posx,posy,posz,poslabel,possymbol,posoccup,posadp;
      posx=loop->second.find("_atom_site_fract_x");
      posy=loop->second.find("_atom_site_fract_y");
      posz=loop->second.find("_atom_site_fract_z");
      unsigned int nb = 0;
      if( (posx!=loop->second.end()) && (posy!=loop->second.end()) && (posz!=loop->second.end()))
      {
         nb=posx->second.size();
         mvAtom.resize(nb);
         for(unsigned int i=0;i<nb;++i)
         {
            mvAtom[i].mCoordFrac.resize(3);
            mvAtom[i].mCoordFrac[0]=CIFNumeric2REAL(posx->second[i]);
            mvAtom[i].mCoordFrac[1]=CIFNumeric2REAL(posy->second[i]);
            mvAtom[i].mCoordFrac[2]=CIFNumeric2REAL(posz->second[i]);
         }
         this->Fractional2CartesianCoord();
      }
      else
      {
         posx=loop->second.find("_atom_site_Cartn_x");
         posy=loop->second.find("_atom_site_Cartn_y");
         posz=loop->second.find("_atom_site_Cartn_z");
         if( (posx!=loop->second.end()) && (posy!=loop->second.end()) && (posz!=loop->second.end()))
         {
            nb=posx->second.size();
            mvAtom.resize(nb);
            for(unsigned int i=0;i<nb;++i)
            {
               mvAtom[i].mCoordCart.resize(3);
               mvAtom[i].mCoordCart[0]=CIFNumeric2REAL(posx->second[i]);
               mvAtom[i].mCoordCart[1]=CIFNumeric2REAL(posy->second[i]);
               mvAtom[i].mCoordCart[2]=CIFNumeric2REAL(posz->second[i]);
            }
            this->Cartesian2FractionalCoord();
         }
      }
      if(mvAtom.size()>0)
      {// Got the atoms, get names, symbols and adps
         (*fpObjCrystInformUser)("CIF: Extract Atoms...");
         possymbol=loop->second.find("_atom_site_type_symbol");
         if(possymbol!=loop->second.end())
            for(unsigned int i=0;i<nb;++i)
               mvAtom[i].mSymbol=possymbol->second[i];
         poslabel=loop->second.find("_atom_site_label");
         if(poslabel!=loop->second.end())
            for(unsigned int i=0;i<nb;++i)
            {
               mvAtom[i].mLabel=poslabel->second[i];
               if(possymbol==loop->second.end())
               {// There was no symbol, use the labels to guess it
                  int nbc=0;
                  if(mvAtom[i].mLabel.size()==1)
                     if(isalpha(mvAtom[i].mLabel[0])) nbc=1;
                  if(mvAtom[i].mLabel.size()>=2)
                  {
                     if(isalpha(mvAtom[i].mLabel[0]) && isalpha(mvAtom[i].mLabel[1])) nbc=2;
                     else if(isalpha(mvAtom[i].mLabel[0])) nbc=1;
                  }
                  if(nbc>0) mvAtom[i].mSymbol=mvAtom[i].mLabel.substr(0,nbc);
                  else mvAtom[i].mSymbol="H";//Something wen wrong, no symbol !
               }
            }
         // Occupancy ?
         posoccup=loop->second.find("_atom_site_occupancy");
         if(posoccup!=loop->second.end())
            for(unsigned int i=0;i<nb;++i)
               mvAtom[i].mOccupancy=CIFNumeric2REAL(posoccup->second[i]);
         // ADPs - Record ani, ovl or mpl as iso.
         REAL mult = 1.0;
         posadp=loop->second.find("_atom_site_B_iso_or_equiv");
         if(posadp==loop->second.end())
         {
            mult = 8 * M_PI * M_PI;
            posadp=loop->second.find("_atom_site_U_iso_or_equiv");
         }
         if(posadp!=loop->second.end())
            for(unsigned int i=0;i<nb;++i)
               mvAtom[i].mBiso = mult*CIFNumeric2REAL(posadp->second[i]);
         // Now be somewhat verbose
         if(verbose)
         {
            cout << "Found "<<nb<<" atoms. Waouh !"<<endl;
            for(unsigned int i=0;i<nb;++i)
            {
               cout<<mvAtom[i].mLabel<<" "<<mvAtom[i].mSymbol;
               if(mvAtom[i].mCoordFrac.size()>0)
               {
                  cout<<" , Fractional: ";
                  for(unsigned int j=0;j<mvAtom[i].mCoordFrac.size();++j)
                     cout<<mvAtom[i].mCoordFrac[j]<<" ";
               }
               if(mvAtom[i].mCoordCart.size()>0)
               {
                  cout<<" , Cartesian: ";
                  for(unsigned int j=0;j<mvAtom[i].mCoordCart.size();++j)
                     cout<<mvAtom[i].mCoordCart[j]<<" ";
               }
               cout<<" , Occupancy= "<<mvAtom[i].mOccupancy<<endl;
               cout<<" , Biso= "<<mvAtom[i].mBiso<<endl;
            }
         }
      }
   }
}

void CIFData::ExtractAnisotropicADPs(const bool verbose)
{

   typedef map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator LoopIter;
   typedef map<ci_string,vector<string> >::const_iterator EntryIter;

   const REAL utob = 8 * M_PI * M_PI;

   const char* uijlabels[] = {
       "_atom_site_aniso_U_11",
       "_atom_site_aniso_U_22",
       "_atom_site_aniso_U_33",
       "_atom_site_aniso_U_12",
       "_atom_site_aniso_U_13",
       "_atom_site_aniso_U_23",
   };

   const char* bijlabels[] = {
       "_atom_site_aniso_B_11",
       "_atom_site_aniso_B_22",
       "_atom_site_aniso_B_33",
       "_atom_site_aniso_B_12",
       "_atom_site_aniso_B_13",
       "_atom_site_aniso_B_23"
   };

   REAL mult[6];

   EntryIter anisolabels, beta11, beta22, beta33, beta12, beta13, beta23;

   EntryIter* betaiters[] = {&beta11, &beta22, &beta33, &beta12, &beta13, &beta23};

   for(LoopIter loop=mvLoop.begin(); loop!=mvLoop.end();++loop)
   {

      // Start with anisotropic factors. If we can find the the
      // "_atom_site_aniso_label" tag, we then want to look for the aniso
      // information.
      anisolabels = loop->second.find("_atom_site_aniso_label");

      // Move to the next loop if we can't find it here.
      if (anisolabels == loop->second.end()) continue;
      if(verbose) cout << "Found labels!" << endl;

      // We have a list of labels. Position the iterators for each of the
      // adps.
      for (int idx = 0; idx < 6; ++idx)
      {
         EntryIter& betaiter = *betaiters[idx];
         betaiter = loop->second.find(bijlabels[idx]);
         mult[idx] = 1.0;
         if(betaiter == loop->second.end())
         {
            betaiter = loop->second.find(uijlabels[idx]);
            mult[idx] = utob;
         }
         if(betaiter == loop->second.end()) mult[idx] = 0.0;
      }

      // Check that we have values. If not, then we can get out of here.
      bool havedata = false;
      for (int i = 0; i < 6; ++i)
      {
         if( mult[i] != 0 ) havedata = true;
      }
      if (!havedata) return;

      // Now loop over the labels, find the corresponding CIFAtom, and fill in
      // its information.
      size_t nb = anisolabels->second.size();
      if(verbose) cout << "Have " << nb << " labels." << endl;
      for (size_t i = 0; i < nb; ++i)
      {
         string label = anisolabels->second[i];
         if(verbose) cout << label << endl;

         // See if we have a CIFAtom with this label. If so, initialize the mBeta
         // vector.
         vector<CIFAtom>::iterator atom = mvAtom.begin();
         for (; atom != mvAtom.end(); ++atom)
         {
            if (atom->mLabel == label)
            {
               atom->mBeta.resize(6, 0.0);
               break;
            }
         }
         // If we didn't find the mvAtom, then we should move on to the next
         // label.
         if (atom == mvAtom.end()) continue;

         // Fill in what we've got, one entry at a time.
         for (int idx=0; idx<6; ++idx)
         {
            if (mult[idx] == 0)
            {
               if(verbose) cout << "skipping index " << idx << endl;
               continue;
            }

            EntryIter& betaiter = *betaiters[idx];

            if (betaiter->second.size() <= i) continue;

            double beta = CIFNumeric2REAL(betaiter->second[i]);
            atom->mBeta[idx] = mult[idx] * beta;

            if(verbose) cout << "mBeta " << idx << " " << atom->mBeta[idx] << endl;
         }
      }
   }
   return;
}


/// This is the default wavelength - whenever a "_diffrn_radiation_wavelength" or
/// "_pd_proc_wavelength"entry is found,
/// it is used as a new value for the default wavelength. Since the powder CIFs do not
/// include the wavelength, this could be useful if the crystal structure CIF (including
/// the wavelength) is parsed right before the powder pattern one.
static REAL defaultWavelength=1.0;

void CIFData::ExtractPowderPattern(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_diffrn_radiation_wavelength");
   if(positem==mvItem.end()) positem=mvItem.find("_pd_proc_wavelength");
   if(positem!=mvItem.end())
   {
      mWavelength=CIFNumeric2REAL(positem->second);
      defaultWavelength=mWavelength;
      cout<<"Found wavelength:"<<defaultWavelength<<endl;
   }
   else mWavelength=defaultWavelength;

   /// Now find the data
   for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
       loop!=mvLoop.end();++loop)
   {
      mDataType=WAVELENGTH_MONOCHROMATIC;
      map<ci_string,vector<string> >::const_iterator pos_x,pos_iobs,pos_weight,pos_mon,pos_wavelength;
      pos_wavelength=loop->second.find("_diffrn_radiation_wavelength");
      if(pos_wavelength!=loop->second.end())
      {
         cout<<"Found wavelength (in loop):"<<pos_wavelength->second[0];
         mWavelength=CIFNumeric2REAL(pos_wavelength->second[0]);
         defaultWavelength=mWavelength;
         cout<<" -> "<<defaultWavelength<<endl;
      }

      pos_iobs=loop->second.find("_pd_meas_counts_total");
      if(pos_iobs==loop->second.end()) pos_iobs=loop->second.find("_pd_meas_intensity_total");
      if(pos_iobs==loop->second.end()) pos_iobs=loop->second.find("_pd_proc_intensity_total");
      if(pos_iobs==loop->second.end()) pos_iobs=loop->second.find("_pd_proc_intensity_net");
      if(pos_iobs==loop->second.end()) continue;//no observed powder data found
      pos_weight=loop->second.find("_pd_proc_ls_weight");
      pos_x=loop->second.find("_pd_proc_2theta_corrected");
      if(pos_x==loop->second.end()) pos_x=loop->second.find("_pd_meas_angle_2theta");
      if(pos_x==loop->second.end()) pos_x=loop->second.find("_pd_meas_2theta_scan");
      if(pos_x==loop->second.end())
      {
         pos_x=loop->second.find("_pd_meas_time_of_flight");
         if(pos_x!=loop->second.end()) mDataType=WAVELENGTH_TOF;
      }

      bool x_fixed_step=false;
      REAL xmin = 0, xmax = 0, xinc = 0;
      POSSIBLY_UNUSED(xmax);
      if(pos_x==loop->second.end())
      {
         map<ci_string,string>::const_iterator pos_min,pos_max,pos_inc;
         pos_min=mvItem.find("_pd_proc_2theta_range_min");
         if(pos_min==mvItem.end()) pos_min=mvItem.find("_pd_meas_2theta_range_min");
         pos_max=mvItem.find("_pd_proc_2theta_range_max");
         if(pos_max==mvItem.end()) pos_max=mvItem.find("_pd_meas_2theta_range_max");
         pos_inc=mvItem.find("_pd_proc_2theta_range_inc");
         if(pos_inc==mvItem.end()) pos_inc=mvItem.find("_pd_meas_2theta_range_inc");
         if((pos_min!=mvItem.end()) && (pos_max!=mvItem.end()) && (pos_inc!=mvItem.end()) )
         {
            x_fixed_step=true;
            xmin=CIFNumeric2REAL(pos_min->second);
            xmax=CIFNumeric2REAL(pos_max->second);
            xinc=CIFNumeric2REAL(pos_inc->second);
         }
      }
      pos_mon=loop->second.find("_pd_meas_intensity_monitor");
      if(pos_mon==loop->second.end()) pos_mon=loop->second.find("_pd_meas_step_count_time");

      if( (pos_iobs!=loop->second.end()) && ( (pos_x!=loop->second.end()) || (x_fixed_step)) )
      {// Found powder data !
         const long nb=pos_iobs->second.size();
         if(verbose) cout<<"Found powder data, with "<<nb<<" data points"<<endl;
         mPowderPatternObs.resize(nb);
         mPowderPatternX.resize(nb);
         mPowderPatternSigma.resize(nb);
         REAL mult=1.0;
         if(mDataType!=WAVELENGTH_TOF) mult=0.017453292519943295;
         for(long i=0;i<nb;++i)
         {
            mPowderPatternObs[i]=CIFNumeric2REAL(pos_iobs->second[i]);
            if(x_fixed_step) mPowderPatternX[i]=(xmin+i*xinc)*mult;
            else mPowderPatternX[i]=CIFNumeric2REAL(pos_x->second[i])*mult;
            // :TODO: use esd on observed intensity, if available.
            if(pos_weight!=loop->second.end())
            {
               mPowderPatternSigma[i]=CIFNumeric2REAL(pos_weight->second[i]);
               if(mPowderPatternSigma[i]>0) mPowderPatternSigma[i]=1/sqrt(fabs(mPowderPatternSigma[i]));
               else mPowderPatternSigma[i]=sqrt(fabs(mPowderPatternObs[i])); // :KLUDGE: ?
            }
            else mPowderPatternSigma[i]=sqrt(fabs(mPowderPatternObs[i]));
            if(pos_mon!=loop->second.end())
            {//VCT or monitor
               const REAL mon=CIFNumeric2REAL(pos_mon->second[i]);
               if(mon>0)
               {
                  mPowderPatternObs[i]/=mon;
                  mPowderPatternSigma[i]/=sqrt(mon);
               }
            }
            //if((i<10) && verbose) cout<<i<<" "<<mPowderPatternX[i]/mult<<" "<<mPowderPatternObs[i]<<" "<<mPowderPatternSigma[i]<<endl;
         }
      }
   }
}

void CIFData::ExtractSingleCrystalData(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_diffrn_radiation_wavelength");
   if(positem==mvItem.end()) positem=mvItem.find("_pd_proc_wavelength");
   if(positem!=mvItem.end())
   {
      mWavelength=CIFNumeric2REAL(positem->second);
      defaultWavelength=mWavelength;
      cout<<"Found wavelength:"<<defaultWavelength<<endl;
   }
   else mWavelength=defaultWavelength;

   /// Now find the data
   for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
       loop!=mvLoop.end();++loop)
   {
      mDataType=WAVELENGTH_MONOCHROMATIC;
      map<ci_string,vector<string> >::const_iterator pos_h,pos_k,pos_l,pos_iobs,pos_sigma,pos_wavelength;
      pos_wavelength=loop->second.find("_diffrn_radiation_wavelength");
      if(pos_wavelength!=loop->second.end())
      {
         cout<<"Found wavelength (in loop):"<<pos_wavelength->second[0];
         mWavelength=CIFNumeric2REAL(pos_wavelength->second[0]);
         defaultWavelength=mWavelength;
         cout<<" -> "<<defaultWavelength<<endl;
      }

      pos_iobs=loop->second.find("_refln_F_squared_meas");
      if(pos_iobs==loop->second.end()) continue;//no observed powder data found
      pos_sigma=loop->second.find("_refln_F_squared_sigma");
      pos_h=loop->second.find("_refln_index_h");
      pos_k=loop->second.find("_refln_index_k");
      pos_l=loop->second.find("_refln_index_l");

      if( (pos_iobs!=loop->second.end()) && (pos_h!=loop->second.end()) && (pos_k!=loop->second.end()) && (pos_l!=loop->second.end()))
      {// Found single crystal data !
         const long nb=pos_iobs->second.size();
         if(verbose) cout<<"Found single crystal data, with "<<nb<<" data points"<<endl;
         mIobs.resize(nb);
         mH.resize(nb);
         mK.resize(nb);
         mL.resize(nb);
         mSigma.resize(nb);
         for(long i=0;i<nb;++i)
         {
            mIobs(i)=CIFNumeric2REAL(pos_iobs->second[i]);
            mH(i)=CIFNumeric2Int(pos_h->second[i]);
            mK(i)=CIFNumeric2Int(pos_k->second[i]);
            mL(i)=CIFNumeric2Int(pos_l->second[i]);
            if(pos_sigma!=loop->second.end()) mSigma(i)=CIFNumeric2REAL(pos_sigma->second[i]);
            else mSigma(i)=sqrt(fabs(abs(mIobs(i))));
         }
      }
   }
}

void CIFData::CalcMatrices(const bool verbose)
{
   if(mvLatticePar.size()==0) return;//:TODO: throw error
   REAL a,b,c,alpha,beta,gamma;//direct space parameters
   REAL aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
   POSSIBLY_UNUSED(aa);
   POSSIBLY_UNUSED(bb);
   POSSIBLY_UNUSED(betaa);
   POSSIBLY_UNUSED(gammaa);
   REAL v;//volume of the unit cell
   a=mvLatticePar[0];
   b=mvLatticePar[1];
   c=mvLatticePar[2];
   alpha=mvLatticePar[3];
   beta=mvLatticePar[4];
   gamma=mvLatticePar[5];

   v=sqrt(fabs(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
               +2*cos(alpha)*cos(beta)*cos(gamma)));

   aa=sin(alpha)/a/v;
   bb=sin(beta )/b/v;
   cc=sin(gamma)/c/v;

   alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
   betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
   gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );

   mOrthMatrix[0][0]=a;
   mOrthMatrix[0][1]=b*cos(gamma);
   mOrthMatrix[0][2]=c*cos(beta);

   mOrthMatrix[1][0]=0;
   mOrthMatrix[1][1]=b*sin(gamma);
   mOrthMatrix[1][2]=-c*sin(beta)*cos(alphaa);

   mOrthMatrix[2][0]=0;
   mOrthMatrix[2][1]=0;
   mOrthMatrix[2][2]=1/cc;

   // Invert upper triangular matrix
   REAL cm[3][3];
   cm[0][0]=mOrthMatrix[0][0];
   cm[0][1]=mOrthMatrix[0][1];
   cm[0][2]=mOrthMatrix[0][2];

   cm[1][0]=mOrthMatrix[1][0];
   cm[1][1]=mOrthMatrix[1][1];
   cm[1][2]=mOrthMatrix[1][2];

   cm[2][0]=mOrthMatrix[2][0];
   cm[2][1]=mOrthMatrix[2][1];
   cm[2][2]=mOrthMatrix[2][2];
   for(long i=0;i<3;i++)
      for(long j=0;j<3;j++)
         if(i==j) mOrthMatrixInvert[i][j]=1;
         else mOrthMatrixInvert[i][j]=0;
   for(long i=0;i<3;i++)
   {
      REAL a;
      for(long j=i-1;j>=0;j--)
      {
         a=cm[j][i]/cm[i][i];
         for(long k=0;k<3;k++) mOrthMatrixInvert[j][k] -= mOrthMatrixInvert[i][k]*a;
         for(long k=0;k<3;k++) cm[j][k] -= cm[i][k]*a;
      }
      a=cm[i][i];
      for(long k=0;k<3;k++) mOrthMatrixInvert[i][k] /= a;
      for(long k=0;k<3;k++) cm[i][k] /= a;
   }
   if(verbose)
   {
      cout <<"Fractional2Cartesian matrix:"<<endl
           <<mOrthMatrix[0][0]<<" "<<mOrthMatrix[0][1]<<" "<<mOrthMatrix[0][2]<<endl
           <<mOrthMatrix[1][0]<<" "<<mOrthMatrix[1][1]<<" "<<mOrthMatrix[1][2]<<endl
           <<mOrthMatrix[2][0]<<" "<<mOrthMatrix[2][1]<<" "<<mOrthMatrix[2][2]<<endl<<endl;
      /*
      cout <<cm[0][0]<<" "<<cm[0][1]<<" "<<cm[0][2]<<endl
           <<cm[1][0]<<" "<<cm[1][1]<<" "<<cm[1][2]<<endl
           <<cm[2][0]<<" "<<cm[2][1]<<" "<<cm[2][2]<<endl<<endl;
      */
      cout <<"Cartesian2Fractional matrix:"<<endl
           <<mOrthMatrixInvert[0][0]<<" "<<mOrthMatrixInvert[0][1]<<" "<<mOrthMatrixInvert[0][2]<<endl
           <<mOrthMatrixInvert[1][0]<<" "<<mOrthMatrixInvert[1][1]<<" "<<mOrthMatrixInvert[1][2]<<endl
           <<mOrthMatrixInvert[2][0]<<" "<<mOrthMatrixInvert[2][1]<<" "<<mOrthMatrixInvert[2][2]<<endl<<endl;
   }
}

void CIFData::f2c(REAL &x,REAL &y, REAL &z)
{
   const REAL x0=x,y0=y,z0=z;
   x=mOrthMatrix[0][0]*x0+mOrthMatrix[0][1]*y0+mOrthMatrix[0][2]*z0;
   y=mOrthMatrix[1][0]*x0+mOrthMatrix[1][1]*y0+mOrthMatrix[1][2]*z0;
   z=mOrthMatrix[2][0]*x0+mOrthMatrix[2][1]*y0+mOrthMatrix[2][2]*z0;
}

void CIFData::c2f(REAL &x,REAL &y, REAL &z)
{
   const REAL x0=x,y0=y,z0=z;
   x=mOrthMatrixInvert[0][0]*x0+mOrthMatrixInvert[0][1]*y0+mOrthMatrixInvert[0][2]*z0;
   y=mOrthMatrixInvert[1][0]*x0+mOrthMatrixInvert[1][1]*y0+mOrthMatrixInvert[1][2]*z0;
   z=mOrthMatrixInvert[2][0]*x0+mOrthMatrixInvert[2][1]*y0+mOrthMatrixInvert[2][2]*z0;
}

void CIFData::Cartesian2FractionalCoord()
{
   for(vector<CIFAtom>::iterator pos=mvAtom.begin();pos!=mvAtom.end();++pos)
   {
      pos->mCoordFrac.resize(3);
      pos->mCoordFrac[0]=pos->mCoordCart.at(0);
      pos->mCoordFrac[1]=pos->mCoordCart.at(1);
      pos->mCoordFrac[2]=pos->mCoordCart.at(2);
      c2f(pos->mCoordFrac[0],pos->mCoordFrac[1],pos->mCoordFrac[2]);
   }
}

void CIFData::Fractional2CartesianCoord()
{
   for(vector<CIFAtom>::iterator pos=mvAtom.begin();pos!=mvAtom.end();++pos)
   {
      pos->mCoordCart.resize(3);
      pos->mCoordCart[0]=pos->mCoordFrac.at(0);
      pos->mCoordCart[1]=pos->mCoordFrac.at(1);
      pos->mCoordCart[2]=pos->mCoordFrac.at(2);
      f2c(pos->mCoordCart[0],pos->mCoordCart[1],pos->mCoordCart[2]);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CIF::CIF(istream &is, const bool interpret,const bool verbose)
{
   string s;
   (*fpObjCrystInformUser)("CIF: Opening CIF");
   Chronometer chrono;
   chrono.start();
   //Copy to an iostream so that we can put back characters if necessary
   stringstream in;
   char c;
   while(is.get(c))in.put(c);
   const float t0read=chrono.seconds();
   s=(boost::format("CIF: Parsing CIF (reading dt=%5.3fs)")%t0read).str();
   (*fpObjCrystInformUser)(s);
   this->Parse(in);
   const float t1parse=chrono.seconds();
   s=(boost::format("CIF: Finished Parsing, Extracting...(parsing dt=%5.3fs)") % (t1parse-t0read)).str();
   (*fpObjCrystInformUser)(s);
   // Extract structure from blocks
   if(interpret)
      for(map<string,CIFData>::iterator posd=mvData.begin();posd!=mvData.end();++posd)
         posd->second.ExtractAll(verbose);
   const float t2interpret=chrono.seconds();
   s=(boost::format("CIF: Finished Import...(interpret dt=%5.3fs, total CIF import=%5.3fs)")%(t2interpret-t1parse)%t2interpret).str();
   (*fpObjCrystInformUser)(s);
}

bool iseol(const char c) { return ((c=='\n')||(c=='\r'));}

std::string trimString(const std::string &s)
{
   const size_t i0 = s.find_first_not_of(" \t\r\n");
   if (i0 == std::string::npos) return "";
   const size_t i1 = s.find_last_not_of(" \t\r\n");
   return s.substr(i0, i1-i0+1);
}

/// Read one value, whether it is numeric, string or text
string CIFReadValue(stringstream &in,char &lastc)
{
   bool vv=false;//very verbose ?
   string value;
   while(!isgraph(in.peek())) in.get(lastc);
   while(in.peek()=='#')
   {//discard these comments for now
      string tmp;
      getline(in,tmp);
      lastc='\r';
      while(!isgraph(in.peek())) in.get(lastc);
   }
   if(in.peek()==';')
   {//SemiColonTextField
      bool warning=!iseol(lastc);
      if(warning)
         cout<<"WARNING: Trying to read a SemiColonTextField but last char is not an end-of-line char !"<<endl;
      value="";
      in.get(lastc);
      while(in.peek()!=';')
      {
         string tmp;
         getline(in,tmp);
         value+=tmp+" ";
      }
      in.get(lastc);
      if(vv) cout<<"SemiColonTextField:"<<value<<endl;
      if(warning && !vv) cout<<"SemiColonTextField:"<<value<<endl;
      return trimString(value);
   }
   if((in.peek()=='\'') || (in.peek()=='\"'))
   {//QuotedString
      char delim;
      in.get(delim);
      value="";
      while(!((lastc==delim)&&(!isgraph(in.peek()))) )
      {
         in.get(lastc);
         value+=lastc;
      }
      if(vv) cout<<"QuotedString:"<<value<<endl;
      return trimString(value.substr(0,value.size()-1));
   }
   // If we got here, we have an ordinary value, numeric or unquoted string
   in>>value;
   if(vv) cout<<"NormalValue:"<<value<<endl;
   return value;
}

void CIF::Parse(stringstream &in)
{
   bool vv=false;//very verbose ?
   char lastc=' ';
   string block="";// Current block data
   while(!in.eof())
   {
      //stringstream mess;
      //mess<<"CIF: Parsing:"<<in.tellg();
      //(*fpObjCrystInformUser)(mess.str());
      while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
      if(in.eof()) break;
      if(vv) cout<<endl;
      if(in.peek()=='#')
      {//Comment
         string tmp;
         getline(in,tmp);
         if(block=="") mvComment.push_back(tmp);
         else mvData[block].mvComment.push_back(tmp);
         lastc='\r';
         if(vv)cout<<"Comment:"<<tmp<<endl;
         continue;
      }
      if(in.peek()=='_')
      {//Tag
         string tag,value;
         in>>tag;
         // Convert all dots to underscores to cover much of DDL2 with this DDL1 parser.
         for (string::size_type pos = tag.find('.'); pos != string::npos; pos = tag.find('.', ++ pos))
            tag.replace(pos, 1, 1, '_');
         value=CIFReadValue(in,lastc);
         if(value==string("?")) continue;//useless
         mvData[block].mvItem[ci_string(tag.c_str())]=value;
         if(vv)cout<<"New Tag:"<<tag<<" ("<<value.size()<<"):"<<value<<endl;
         continue;
      }
      if((in.peek()=='d') || (in.peek()=='D'))
      {// Data
         string tmp;
         in>>tmp;
         block=tmp.substr(5);
         if(vv) cout<<endl<<endl<<"NEW BLOCK DATA: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ->"<<block<<endl<<endl<<endl;
         mvData[block]=CIFData();
         continue;
      }
      if((in.peek()=='l') || (in.peek()=='L'))
      {// loop_
         vector<ci_string> tit;
         string tmp;
         in>>tmp; //should be loop_
         if(vv) cout<<"LOOP : "<<tmp;
         while(!in.eof())
         {//read titles
            while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
            if(in.peek()=='#')
            {
               getline(in,tmp);
               if(block=="") mvComment.push_back(tmp);
               else mvData[block].mvComment.push_back(tmp);
               continue;
            }
            if(in.peek()!='_')
            {
               if(vv) cout<<endl<<"End of loop titles:"<<(char)in.peek()<<endl;
               break;
            }
            in>>tmp;
            // Convert all dots to underscores to cover much of DDL2 with this DDL1 parser.
            for (string::size_type pos = tmp.find('.'); pos != string::npos; pos = tmp.find('.', ++ pos))
               tmp.replace(pos, 1, 1, '_');
            tit.push_back(ci_string(tmp.c_str()));
            if(vv) cout<<" , "<<tmp;
         }
         if(vv) cout<<endl;
         map<ci_string,vector<string> > lp;
         while(true)
         {
            while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
            if(in.eof()) break;
            if(vv) cout<<"LOOP VALUES...: "<<(char)in.peek()<<" "<<endl;
            if(in.peek()=='_') break;
            if(in.peek()=='#')
            {// Comment (in a loop ??)
               string tmp;
               getline(in,tmp);
               if(block=="") mvComment.push_back(tmp);
               else mvData[block].mvComment.push_back(tmp);
               lastc='\r';
               if(vv) cout<<"Comment in a loop (?):"<<tmp<<endl;
               continue;
            };
            const std::ios::pos_type pos=in.tellg();
            in>>tmp;
            if(vv) cout<<"WHATNEXT? "<<tmp;
            if(ci_string(tmp.c_str())=="loop_")
            {//go back and continue
               if(vv) cout<<endl<<"END OF LOOP :"<<tmp<<endl;
               in.seekg(pos);
               break;
            }
            if(tmp.size()>=5)
               if(ci_string(tmp.substr(0,5).c_str())=="data_")
               {//go back and continue
                  if(vv) cout<<endl<<"END OF LOOP :"<<tmp<<endl;
                  in.seekg(pos);
                  break;
               }
            // go back
            in.seekg(pos);
            for(unsigned int i=0;i<tit.size();++i)
            {//Read all values
               const string value=CIFReadValue(in,lastc);
               lp[tit[i]].push_back(value);
               if(vv) cout<<"     #"<<i<<" :  "<<value<<endl;
            }
         }
         // The key to the mvLoop map is the set of column titles
         set<ci_string> stit;
         for(unsigned int i=0;i<tit.size();++i) stit.insert(tit[i]);
         mvData[block].mvLoop[stit]=lp;
         continue;
      }
      // If we get here, something went wrong ! Discard till end of line...
      string junk;
      getline(in,junk);
      cout<<"WARNING: did not understand : "<<junk<<endl;
   }
}

REAL CIFNumeric2REAL(const string &s)
{
   if((s==".") || (s=="?")) return 0.0;
   // Use stream rather than sscanf to rely on C++ locale to read C-locale values
   REAL v=0;
   stringstream ss(s);
   ss.imbue(std::locale::classic());
   ss>>v;
   return v;
}

int CIFNumeric2Int(const string &s)
{
   if((s==".") || (s=="?")) return 0;
   // Use stream rather than sscanf to rely on C++ locale to read C-locale values
   int v=0;
   stringstream ss(s);
   ss.imbue(std::locale::classic());
   ss>>v;
   return v;
}

Crystal* CreateCrystalFromCIF(CIF &cif,bool verbose,bool checkSymAsXYZ)
{
   return CreateCrystalFromCIF(cif,verbose,checkSymAsXYZ,false,false);
}

Crystal* CreateCrystalFromCIF(CIF &cif,const bool verbose,const bool checkSymAsXYZ, const bool oneScatteringPowerPerElement, const bool connectAtoms)
{
   (*fpObjCrystInformUser)("CIF: Opening CIF");
   Chronometer chrono;
   chrono.start();

   // If oneScatteringPowerPerElement==true, we hold this to compute the average Biso per element
   std::map<ScatteringPower*,std::pair<REAL,unsigned int> > vElementBiso;

   Crystal *pCryst=NULL;
   for(map<string,CIFData>::iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
      if(pos->second.mvLatticePar.size()==6)
      {
         // If no atoms are listed and we already have a crystal structure defined,
         //asssume we don't want this one - e.g. like some IuCr journals single crystal
         //data cif files including cell parameters
         if((pos->second.mvAtom.size()==0) && (gCrystalRegistry.GetNb()>0)) continue;
         // Use unambigous Hall symbol if present, otherwise try HM symbol or spg number
         string spg;
         if(pos->second.mSpacegroupSymbolHall!="") try
         {
            tmp_C_Numeric_locale tmploc;
            cctbx::sgtbx::space_group cctbxspg(pos->second.mSpacegroupSymbolHall);
            cctbxspg.t_den();
            cctbxspg.n_smx();
            cctbxspg.n_ltr();
            cctbxspg.type();
            cctbxspg.type().number();
            cctbxspg.type().hall_symbol();
            cctbxspg.type().lookup_symbol();
            cctbxspg.match_tabulated_settings().extension();
            cctbxspg.match_tabulated_settings().hermann_mauguin();
            cctbxspg.type().universal_hermann_mauguin_symbol();
            cctbx::sgtbx::brick b(cctbxspg.type());
            spg=pos->second.mSpacegroupSymbolHall;
         }
         catch(exception)
         {
            VFN_DEBUG_MESSAGE("CreateCrystalFromCIF(): could not interpret Hall symbol:"<<pos->second.mSpacegroupSymbolHall, 10)
         }
         if((spg=="") && (pos->second.mSpacegroupHermannMauguin!="")) try
         {
            tmp_C_Numeric_locale tmploc;
            cctbx::sgtbx::space_group cctbxspg(cctbx::sgtbx::space_group_symbols(pos->second.mSpacegroupHermannMauguin));
            cctbxspg.t_den();
            cctbxspg.n_smx();
            cctbxspg.n_ltr();
            cctbxspg.type();
            cctbxspg.type().number();
            cctbxspg.type().hall_symbol();
            cctbxspg.type().lookup_symbol();
            cctbxspg.type().universal_hermann_mauguin_symbol();
            cctbxspg.match_tabulated_settings().extension();
            cctbxspg.match_tabulated_settings().hermann_mauguin();
            cctbx::sgtbx::brick b(cctbxspg.type());
            spg=pos->second.mSpacegroupHermannMauguin;
         }
         catch(exception)
         {
            VFN_DEBUG_MESSAGE("CreateCrystalFromCIF(): could not interpret Hermann-Mauguin symbol:"<<pos->second.mSpacegroupHermannMauguin, 10)
         }
         if((spg=="") && (pos->second.mSpacegroupNumberIT!=""))
         try
         {
            tmp_C_Numeric_locale tmploc;
            cctbx::sgtbx::space_group cctbxspg(cctbx::sgtbx::space_group_symbols(pos->second.mSpacegroupNumberIT));
            cctbxspg.t_den();
            cctbxspg.n_smx();
            cctbxspg.n_ltr();
            cctbxspg.type();
            cctbxspg.type().number();
            cctbxspg.type().hall_symbol();
            cctbxspg.type().lookup_symbol();
            cctbxspg.type().universal_hermann_mauguin_symbol();
            cctbxspg.match_tabulated_settings().extension();
            cctbxspg.match_tabulated_settings().hermann_mauguin();
            cctbx::sgtbx::brick b(cctbxspg.type());
            spg=pos->second.mSpacegroupNumberIT;
         }
         catch(exception)
         {
            VFN_DEBUG_MESSAGE("CreateCrystalFromCIF(): could not interpret spacegroup number (!) :"<<pos->second.mSpacegroupNumberIT, 10)
         }
         if(spg=="") spg="P1";
         if(verbose) cout<<"Create crystal with spacegroup: "<<spg
             <<" / "<<pos->second.mSpacegroupHermannMauguin
             <<" / "<<pos->second.mSpacegroupSymbolHall
             <<" / "<<pos->second.mSpacegroupNumberIT
             <<"-> "<<spg
             <<endl;
         (*fpObjCrystInformUser)("CIF: Create Crystal=");
         pCryst=new Crystal(pos->second.mvLatticePar[0],pos->second.mvLatticePar[1],pos->second.mvLatticePar[2],
                                     pos->second.mvLatticePar[3],pos->second.mvLatticePar[4],pos->second.mvLatticePar[5],spg);
         if(  (pos->second.mSpacegroupSymbolHall=="")
            &&(pos->second.mvSymmetry_equiv_pos_as_xyz.size()>0)
            &&(pos->second.mSpacegroupHermannMauguin!="")
            &&checkSymAsXYZ)
         {// Could not use a Hall symbol, but we have a list of symmetry_equiv_pos_as_xyz,
          // so check we have used the best possible origin
            tmp_C_Numeric_locale tmploc;
            static vector<string> origin_list;
            if(origin_list.size()==0)
            {//ugly ?
               origin_list.resize(5);
               origin_list[0]="";
               origin_list[1]=":1";
               origin_list[2]=":2";
               origin_list[3]=":R";
               origin_list[4]=":H";
            }
            // If we do not have an HM symbol, then use the one generated by cctbx (normally from spg number)
            string hmorig=pos->second.mSpacegroupHermannMauguin;
            if(hmorig=="") hmorig=pCryst->GetSpaceGroup().GetCCTbxSpg().match_tabulated_settings().hermann_mauguin();

            if(verbose) cout<<" Symmetry checking using symmetry_equiv_pos_as_xyz:"<<endl;
            string bestsymbol=hmorig;
            unsigned int bestscore=0;
            for(vector<string>::const_iterator posOrig=origin_list.begin();posOrig!=origin_list.end();++posOrig)
            {
               // The origin extension may not make sense, so we need to watch for exception
               try
               {
                  pCryst->GetSpaceGroup().ChangeSpaceGroup(hmorig+*posOrig);
               }
               catch(invalid_argument)
               {
                  continue;
               }

               // If the symbol is the same as before, the origin probably was not understood - no need to test
               if((posOrig!=origin_list.begin())&&(pCryst->GetSpaceGroup().GetName()==bestsymbol)) continue;

               unsigned int nbSymSpg=pCryst->GetSpaceGroup().GetCCTbxSpg().all_ops().size();
               unsigned int nbSymCommon=0;
               try
               {
                  for(unsigned int i=0;i<nbSymSpg;i++)
                  {
                     for(set<string>::const_iterator posSymCIF=pos->second.mvSymmetry_equiv_pos_as_xyz.begin();
                        posSymCIF!=pos->second.mvSymmetry_equiv_pos_as_xyz.end();++posSymCIF)
                     {
                        cctbx::sgtbx::rt_mx mx1(*posSymCIF);
                        cctbx::sgtbx::rt_mx mx2(pCryst->GetSpaceGroup().GetCCTbxSpg().all_ops()[i]);
                        mx1.mod_positive_in_place();
                        mx2.mod_positive_in_place();
                        if(mx1==mx2)
                        {
                           nbSymCommon++;
                           break;
                        }
                     }
                  }
                  if(verbose) cout<<"   Trying: "<<pCryst->GetSpaceGroup().GetName()
                      <<" nbsym:"<<nbSymSpg<<"(cctbx), "
                      <<pos->second.mvSymmetry_equiv_pos_as_xyz.size()<<"(CIF)"
                      <<",common:"<<nbSymCommon<<endl;
                  if(bestscore<((nbSymSpg==pos->second.mvSymmetry_equiv_pos_as_xyz.size())*nbSymCommon))
                  {
                     bestscore=(nbSymSpg==pos->second.mvSymmetry_equiv_pos_as_xyz.size())*nbSymCommon;
                     bestsymbol=pCryst->GetSpaceGroup().GetName();
                  }
               }
               catch(cctbx::error)
               {
                  cout<<"WOOPS: cctbx error ! Wrong symmetry_equiv_pos_as_xyz strings ?"<<endl;
               }
            }
            if(verbose) cout<<endl<<"Finally using spacegroup name:"<<bestsymbol<<endl;
            pCryst->GetSpaceGroup().ChangeSpaceGroup(bestsymbol);
         }
         if(pos->second.mName!="") pCryst->SetName(pos->second.mName);
         else if(pos->second.mFormula!="") pCryst->SetName(pos->second.mFormula);
         const float t1=chrono.seconds();
         (*fpObjCrystInformUser)((boost::format("CIF: Create Crystal:%s(%s)(dt=%6.3fs)")%pCryst->GetName() % pCryst->GetSpaceGroup().GetName() % t1).str());

         bool doInformUserAllAtoms=true;
         if(pos->second.mvAtom.size()>30) doInformUserAllAtoms = false;
         unsigned int ctatom=0;
         for(vector<CIFData::CIFAtom>::const_iterator posat=pos->second.mvAtom.begin();posat!=pos->second.mvAtom.end();++posat)
         {
            if( (posat->mLabel==".") || (posat->mSymbol==".") || (posat->mLabel.find("dummy")!=std::string::npos) || (posat->mSymbol.find("dummy")!=std::string::npos) )
            {
               if(doInformUserAllAtoms) (*fpObjCrystInformUser)("CIF: Ignoring DUMMY Atom:"+posat->mLabel+"(symbol="+posat->mSymbol+")");
               continue;
            }
            ctatom++;

            //const float t20=chrono.seconds();
            // Try to find an existing scattering power with the same properties, or create a new one
            ScatteringPower* sp=NULL;
            if(oneScatteringPowerPerElement)
            {
               for(unsigned int i=0;i<pCryst->GetScatteringPowerRegistry().GetNb();++i)
               {
                  if(pCryst->GetScatteringPowerRegistry().GetObj(i).GetSymbol()!=posat->mSymbol) continue;
                  vElementBiso[&(pCryst->GetScatteringPowerRegistry().GetObj(i))].first+=posat->mBiso;
                  vElementBiso[&(pCryst->GetScatteringPowerRegistry().GetObj(i))].second+=1;
                  sp=&(pCryst->GetScatteringPowerRegistry().GetObj(i));
                  break;
               }
               if(sp==NULL)
               {
                  if(verbose) cout<<"Scattering power "<<posat->mSymbol<<" not found, creating it..."<<endl;
                  sp = new ScatteringPowerAtom(posat->mSymbol,posat->mSymbol);
                  // Always extract isotropic DP, even with ADPs present
                  // :TODO: if only ADP are listed, calculate isotropic DP
                  vElementBiso[sp].first+=posat->mBiso;
                  vElementBiso[sp].second=1;
                  pCryst->AddScatteringPower(sp);
                  //const float t21=chrono.seconds();
                  //(*fpObjCrystInformUser)((boost::format("CIF: Add scattering power: %s (dt=%6.3fsCrystal creation=%6.3fs total)")% posat->mSymbol % (t21-t20) % t21).str());
               }
            }
            else
            {
               #if 0
               for(unsigned int i=0;i<pCryst->GetScatteringPowerRegistry().GetNb();++i)
               {
                  if(pCryst->GetScatteringPowerRegistry().GetObj(i).GetSymbol()!=posat->mSymbol) continue;
                  if(posat->mBeta.size() == 6)
                  {
                     if(  (pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(0)!=posat->mBeta[0])
                        ||(pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(1)!=posat->mBeta[1])
                        ||(pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(2)!=posat->mBeta[2])
                        ||(pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(3)!=posat->mBeta[3])
                        ||(pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(4)!=posat->mBeta[4])
                        ||(pCryst->GetScatteringPowerRegistry().GetObj(i).GetBij(5)!=posat->mBeta[5])) continue;
                  }
                  else if(posat->mBiso!=pCryst->GetScatteringPowerRegistry().GetObj(i).GetBiso()) continue;
                  sp=&(pCryst->GetScatteringPowerRegistry().GetObj(i));
                  break;
               }
               if(sp==NULL)
               #endif
               {
                  if(verbose) cout<<"Creating new scattering power for:"<<posat->mLabel<<endl;
                  sp = new ScatteringPowerAtom(posat->mLabel,posat->mSymbol);
                  // Always extract isotropic DP, even with ADPs present
                  // :TODO: if only ADP are listed, calculate isotropic DP
                  sp->SetBiso(posat->mBiso);
                  // ADPs ?
                  if(posat->mBeta.size() == 6)
                  {
                     for (int idx=0; idx<6; ++idx) sp->SetBij(idx, posat->mBeta[idx]);
                  }
                  pCryst->AddScatteringPower(sp);
                  //const float t21=chrono.seconds();
                  //(*fpObjCrystInformUser)((boost::format("CIF: Add scattering power: %s (dt=%6.3fsCrystal creation=%6.3fs total)") % posat->mLabel % (t21-t20) % t21).str());
               }
            }
            // (*fpObjCrystInformUser)("CIF: Add Atom:"+posat->mLabel+"("+sp->GetName()+")");
            pCryst->AddScatterer(new Atom(posat->mCoordFrac[0],posat->mCoordFrac[1],posat->mCoordFrac[2],
                                          posat->mLabel,sp,posat->mOccupancy));
            const float t22=chrono.seconds();
            if(doInformUserAllAtoms) (*fpObjCrystInformUser)((boost::format("CIF: new Atom: %s (%s) (Crystal creation=%6.3fs total)") % posat->mLabel % sp->GetName() % t22).str());
            else if (ctatom%20 == 0)(*fpObjCrystInformUser)((boost::format("CIF: imported %u atoms (Crystal creation=%6.3fs total)") % ctatom % t22).str());
         }
         if(oneScatteringPowerPerElement)
         {
            for(std::map<ScatteringPower*,std::pair<REAL,unsigned int> >::iterator pos=vElementBiso.begin();pos!=vElementBiso.end();++pos)
               pos->first->SetBiso(pos->second.first/pos->second.second);
         }
         (*fpObjCrystInformUser)((boost::format("CIF: Finished importing %u atoms (Crystal creation=%6.3fs total)") % ctatom % chrono.seconds()).str());
         if(connectAtoms)
         {
            (*fpObjCrystInformUser)("CIF: connecting atoms");
            pCryst->ConnectAtoms();
            unsigned int ctat=0;
            unsigned int ctmol=0;
            for(int i=0;i<pCryst->GetNbScatterer();i++)
            {
               if(pCryst->GetScatt(i).GetClassName()=="Atom") ctat +=1;
               else if(pCryst->GetScatt(i).GetClassName()=="Molecule") ctmol +=1;
            }
            (*fpObjCrystInformUser)((boost::format("CIF: finished connecting atoms (%u isolated atoms, %u molecules) (Crystal creation=%6.3fs total)") % ctat % ctmol % chrono.seconds()).str());
         }
      }
   return pCryst;
}

PowderPattern* CreatePowderPatternFromCIF(CIF &cif)
{
   PowderPattern* pPow=NULL;
   for(map<string,CIFData>::iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
   {
      if(pos->second.mPowderPatternObs.size()>10)
      {
         pPow=new PowderPattern();
         pPow->ImportPowderPatternCIF(cif);
         (*fpObjCrystInformUser)((boost::format("CIF: Imported POWDER PATTERN, with %d points") % pPow->GetNbPoint()).str());
      }
   }
   return pPow;
}

DiffractionDataSingleCrystal* CreateSingleCrystalDataFromCIF(CIF &cif, Crystal *pcryst)
{
   DiffractionDataSingleCrystal* pData=NULL;
   for(map<string,CIFData>::iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
   {
      if(pos->second.mH.numElements()>0)
      {
         (*fpObjCrystInformUser)((boost::format("CIF: Importing SINGLE CRYSTAL DIFFRACTION data")).str());
         if(pcryst==0)
         {
            if(gCrystalRegistry.GetNb()>0)
            {  // Use last Crystal created
               pcryst=&(gCrystalRegistry.GetObj(gCrystalRegistry.GetNb()-1));
               (*fpObjCrystInformUser)((boost::format("CIF: Importing SINGLE CRYSTAL DIFFRACTION data: using last Crystal structure as corresponding crystal [%s]") % pcryst->GetName().c_str()).str());
            }
            else
            {
               pcryst=new Crystal;
               pcryst->SetName("Crystal data for Single Crystal Diffraction data imported from CIF");
               (*fpObjCrystInformUser)((boost::format("CIF: Importing SINGLE CRYSTAL DIFFRACTION data: creating new empty Crystal structure")).str());
            }
         }
         pData=new DiffractionDataSingleCrystal(*pcryst);
         pData->SetHklIobs(pos->second.mH,pos->second.mK,pos->second.mL,pos->second.mIobs,pos->second.mSigma);
         (*fpObjCrystInformUser)((boost::format("CIF: Imported SINGLE CRYSTAL DIFFRACTION data, with %d reflections") % pData->GetNbRefl()).str());
      }
   }
   return pData;
}

}//namespace
