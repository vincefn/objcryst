#include <ctype.h>
#include <cmath>

#include "ObjCryst/CIF.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/Atom.h"
#include "ObjCryst/PowderPattern.h"

using namespace std;

namespace ObjCryst
{

CIFData::CIFAtom::CIFAtom():
mLabel(""),mSymbol(""),mOccupancy(1.0)
{}

CIFData::CIFData()
{}

void CIFData::ExtractAll(const bool verbose)
{
   // :TODO: convert cartesian to fractional coordinates and vice-versa, if unit cell is known
   // :TODO: Take care of values listed as "." and "?" instead of a real value.
   this->ExtractName(verbose);
   this->ExtractUnitCell(verbose);
   this->ExtractSpacegroup(verbose);
   this->ExtractAtomicPositions(verbose);
   this->ExtractPowderPattern(verbose);
}

void CIFData::ExtractUnitCell(const bool verbose)
{
   map<ci_string,string>::const_iterator positem;
   positem=mvItem.find("_cell_length_a");
   if(positem!=mvItem.end())
   {
      mvLatticePar.resize(6);
      mvLatticePar[0]=CIFNumeric2Float(positem->second);
      positem=mvItem.find("_cell_length_b");
      if(positem!=mvItem.end())
         mvLatticePar[1]=CIFNumeric2Float(positem->second);
      positem=mvItem.find("_cell_length_c");
      if(positem!=mvItem.end())
         mvLatticePar[2]=CIFNumeric2Float(positem->second);
      positem=mvItem.find("_cell_angle_alpha");
      if(positem!=mvItem.end())
         mvLatticePar[3]=CIFNumeric2Float(positem->second);
      positem=mvItem.find("_cell_angle_beta");
      if(positem!=mvItem.end())
         mvLatticePar[4]=CIFNumeric2Float(positem->second);
      positem=mvItem.find("_cell_angle_gamma");
      if(positem!=mvItem.end())
         mvLatticePar[5]=CIFNumeric2Float(positem->second);
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
      mSpacegroupNumberIT=CIFNumeric2Int(positem->second);
      if(verbose) cout<<"Found spacegroup IT number:"<<mSpacegroupNumberIT<<endl;
   }
   else
   {
      positem=mvItem.find("_symmetry_Int_Tables_number");
      if(positem!=mvItem.end())
      {
         mSpacegroupNumberIT=CIFNumeric2Int(positem->second);
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
      map<ci_string,vector<string> >::const_iterator posx,posy,posz,poslabel,possymbol,posoccup;
      posx=loop->second.find("_atom_site_fract_x");
      posy=loop->second.find("_atom_site_fract_y");
      posz=loop->second.find("_atom_site_fract_z");
      unsigned int nb;
      if( (posx!=loop->second.end()) && (posy!=loop->second.end()) && (posz!=loop->second.end()))
      {
         nb=posx->second.size();
         mvAtom.resize(nb);
         for(unsigned int i=0;i<nb;++i)
         {
            mvAtom[i].mCoordFrac.resize(3);
            mvAtom[i].mCoordFrac[0]=CIFNumeric2Float(posx->second[i]);
            mvAtom[i].mCoordFrac[1]=CIFNumeric2Float(posy->second[i]);
            mvAtom[i].mCoordFrac[2]=CIFNumeric2Float(posz->second[i]);
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
               mvAtom[i].mCoordCart[0]=CIFNumeric2Float(posx->second[i]);
               mvAtom[i].mCoordCart[1]=CIFNumeric2Float(posy->second[i]);
               mvAtom[i].mCoordCart[2]=CIFNumeric2Float(posz->second[i]);
            }
            this->Cartesian2FractionalCoord();
         }
      }
      if(mvAtom.size()>0)
      {// Got the atoms, get names and symbols
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
         posoccup=loop->second.find("atom_site_occupancy");
         if(posoccup!=loop->second.end())
            for(unsigned int i=0;i<nb;++i)
               mvAtom[i].mOccupancy=CIFNumeric2Float(posoccup->second[i]);
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
            }
         }
      }
   }
}

void CIFData::ExtractPowderPattern(const bool verbose)
{
   for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
       loop!=mvLoop.end();++loop)
   {
      mDataType=WAVELENGTH_MONOCHROMATIC;
      map<ci_string,vector<string> >::const_iterator pos_x,pos_iobs,pos_weight,pos_mon;
      pos_iobs=loop->second.find("_pd_meas_counts_total");
      if(pos_iobs==loop->second.end())
      {
         pos_iobs=loop->second.find("_pd_meas_intensity_total");
         if(pos_iobs==loop->second.end()) pos_iobs=loop->second.find("_pd_proc_intensity_total");
         if(pos_iobs==loop->second.end()) pos_iobs=loop->second.find("_pd_proc_intensity_net");
         if(pos_iobs!=loop->second.end()) pos_weight=loop->second.find("_pd_proc_ls_weight");
         else continue;//no powder data found
      }
      pos_x=loop->second.find("_pd_proc_2theta_corrected");
      if(pos_x==loop->second.end()) pos_x=loop->second.find("_pd_meas_angle_2theta");
      if(pos_x==loop->second.end())
      {
         pos_x=loop->second.find("_pd_meas_time_of_flight");
         if(pos_x!=loop->second.end()) mDataType=WAVELENGTH_TOF;
      }
      
      bool x_fixed_step=false;
      float xmin,xmax,xinc;
      if(pos_x==loop->second.end())
      {
         map<ci_string,string>::const_iterator pos_min,pos_max,pos_inc;
         pos_min=mvItem.find("_pd_proc_2theta_range_min");
         pos_max=mvItem.find("_pd_proc_2theta_range_max");
         pos_inc=mvItem.find("_pd_proc_2theta_range_inc");
         if((pos_min!=mvItem.end()) && (pos_max!=mvItem.end()) && (pos_inc!=mvItem.end()) )
         {
            x_fixed_step=true;
            xmin=CIFNumeric2Float(pos_min->second);
            xmax=CIFNumeric2Float(pos_max->second);
            xinc=CIFNumeric2Float(pos_inc->second);
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
         float mult=1.0;
         if(mDataType!=WAVELENGTH_TOF) mult=0.017453292519943295;
         for(long i=0;i<nb;++i)
         {
            mPowderPatternObs[i]=CIFNumeric2Float(pos_iobs->second[i]);
            if(x_fixed_step) mPowderPatternX[i]=(xmin+i*xinc)*mult;
            else mPowderPatternX[i]=CIFNumeric2Float(pos_x->second[i])*mult;
            // :TODO: use esd on observed intensity, if available.
            if(pos_weight!=loop->second.end())
            {
               mPowderPatternSigma[i]=CIFNumeric2Float(pos_weight->second[i]);
               if(mPowderPatternSigma[i]>0) mPowderPatternSigma[i]=1/sqrt(mPowderPatternSigma[i]);
               else mPowderPatternSigma[i]=0; // :KLUDGE: ?
            }
            else mPowderPatternSigma[i]=sqrt(mPowderPatternObs[i]);
            if(pos_mon!=loop->second.end())
            {//VCT or monitor
               mPowderPatternObs[i]/=CIFNumeric2Float(pos_mon->second[i]);
               mPowderPatternSigma[i]/=sqrt(CIFNumeric2Float(pos_mon->second[i]));
            }
            //if((i<10) && verbose) cout<<i<<" "<<mPowderPatternX[i]/mult<<" "<<mPowderPatternObs[i]<<" "<<mPowderPatternSigma[i]<<endl;
         }
      }
   }
}

void CIFData::CalcMatrices(const bool verbose)
{
   if(mvLatticePar.size()==0) return;//:TODO: throw error
   float a,b,c,alpha,beta,gamma;//direct space parameters
   float aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
   float v;//volume of the unit cell
   a=mvLatticePar[0];
   b=mvLatticePar[1];
   c=mvLatticePar[2];
   alpha=mvLatticePar[3];
   beta=mvLatticePar[4];
   gamma=mvLatticePar[5];
   
   v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
            +2*cos(alpha)*cos(beta)*cos(gamma));
   
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
   float cm[3][3];
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
      float a;
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

void CIFData::f2c(float &x,float &y, float &z)
{
   const float x0=x,y0=y,z0=z;
   x=mOrthMatrix[0][0]*x0+mOrthMatrix[0][1]*y0+mOrthMatrix[0][2]*z0;
   y=mOrthMatrix[1][0]*x0+mOrthMatrix[1][1]*y0+mOrthMatrix[1][2]*z0;
   z=mOrthMatrix[2][0]*x0+mOrthMatrix[2][1]*y0+mOrthMatrix[2][2]*z0;
}

void CIFData::c2f(float &x,float &y, float &z)
{
   const float x0=x,y0=y,z0=z;
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
   //Copy to an iostream so that we can put back characters if necessary
   stringstream in;
   char c;
   while(is.get(c))in.put(c);
   this->Parse(in);
   // Extract structure from blocks
   if(interpret)
      for(map<string,CIFData>::iterator posd=mvData.begin();posd!=mvData.end();++posd)
         posd->second.ExtractAll(verbose);
}

bool iseol(const char c) { return ((c=='\n')||(c=='\r'));}

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
      return value;
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
      return value.substr(0,value.size()-1);
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
      while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
      if(vv) cout<<endl;
      if(in.peek()=='#')
      {//Comment
         string tmp;
         getline(in,tmp);
         if(block=="") mvComment.push_back(tmp);
         else mvData[block].mvComment.push_back(tmp);
         lastc='\r';
         continue;
      }
      if(in.peek()=='_')
      {//Tag
         string tag,value;
         in>>tag;
         value=CIFReadValue(in,lastc);
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
         while(true)
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

float CIFNumeric2Float(const string &s)
{
   if((s==".") || (s=="?")) return 0.0;
   float v;
   const int n=sscanf(s.c_str(),"%f",&v);
   if(n!=1) return 0.0;
   return v;
}

int CIFNumeric2Int(const string &s)
{
   if((s==".") || (s=="?")) return 0;
   int v;
   const int n=sscanf(s.c_str(),"%d",&v);
   if(n!=1) return 0;
   return v;
}

Crystal* CreateCrystalFromCIF(CIF &cif)
{
   Crystal *pCryst=NULL;
   for(map<string,CIFData>::iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
      if(pos->second.mvLatticePar.size()==6)
      {
         string spg=pos->second.mSpacegroupSymbolHall;
         if(spg=="") spg=pos->second.mSpacegroupHermannMauguin;
         if(spg=="") spg=pos->second.mSpacegroupNumberIT;
         if(spg=="") spg="P1";
         Crystal *pCryst=new Crystal(pos->second.mvLatticePar[0],pos->second.mvLatticePar[1],pos->second.mvLatticePar[2],
                                     pos->second.mvLatticePar[3],pos->second.mvLatticePar[4],pos->second.mvLatticePar[5],spg);
         if(pos->second.mName!="") pCryst->SetName(pos->second.mName);
         else if(pos->second.mFormula!="") pCryst->SetName(pos->second.mFormula);
         
         for(vector<CIFData::CIFAtom>::const_iterator posat=pos->second.mvAtom.begin();posat!=pos->second.mvAtom.end();++posat)
         {
            if(pCryst->GetScatteringPowerRegistry().Find(posat->mSymbol,"ScatteringPowerAtom",true)<0)
            {
               cout<<"Scattering power "<<posat->mSymbol<<" not found, creating it..."<<endl;
               pCryst->AddScatteringPower(new ScatteringPowerAtom(posat->mSymbol,posat->mSymbol));
            }
            pCryst->AddScatterer(new Atom(posat->mCoordFrac[0],posat->mCoordFrac[1],posat->mCoordFrac[2],
                                          posat->mLabel,&(pCryst->GetScatteringPower(posat->mSymbol)),
                                          posat->mOccupancy));
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
         PowderPattern* pPow=new PowderPattern();
         pPow->ImportPowderPatternCIF(cif);
      }
   }
   return pPow;
}

}//namespace
