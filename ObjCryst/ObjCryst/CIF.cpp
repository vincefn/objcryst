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
   this->CalcMatrices(verbose);
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
               mvAtom[i].mLabel=poslabel->second[i];
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


CIF::CIF(istream &in, const bool interpret,const bool verbose)
{
   string tmp;
   do
   {
      getline(in,tmp);
      mvLine.push_back(tmp);
   }
   while(!in.eof());
   
   long i=0;
   //for(list<string>::const_iterator pos=mvLine.begin();pos!=mvLine.end();++pos)
   //   cout<< i++<<" : " <<*pos <<endl;
   this->Parse();
   
   // Extract structure from blocks
   if(interpret)
      for(map<string,CIFData>::iterator posd=mvData.begin();posd!=mvData.end();++posd)
         posd->second.ExtractAll(verbose);
}

/// Read one value, whether it is numeric, string or text
/// If last==true, this is the last value on the line,
/// so it may be a quoted string with the delimiter character allowe dinside the string... (what idiot supported _that_ in the specs ?)
string CIFReadValue(stringstream &sst,list<string>::const_iterator &pos, bool last=false)
{
   string value,tmp;
   if(sst.eof())
   {
      sst.clear();
      sst<<*pos++;//one row of data on multiple lines. STUPID FORMAT !
   }

   value="";
   while(value.size()==0)
   {
      if(last)
      {
         getline(sst,value);
         // remove leading spaces - getline reads everything
         while(true)
         {
            if(value.size()==0) break;
            if(value.at(0)==' ')
               value=value.substr(1);
            else break;
         }
      }
      else sst>>value;
      if(value.size()==0)
      {
         sst.clear();
         sst<<*pos++;
      }
   }
   //cout<<__FILE__<<":"<<__LINE__<<":"<<value<<" , "<<value.at(0)<<endl;
   if(value.at(0)==';')
   {//SemiColonTextField
      value="";
      for(;;)
      {
         tmp=*pos++;
         if(tmp.at(0)==';') break;
         value+=tmp+" ";
      }
      sst.clear();
      sst<<*pos++;
   }
   else
      if((value.at(0)=='\'') || (value.at(0)=='\"'))
      {//QuotedString => remove quotes and trailing spaces
         const char delim=value.at(0);
         value=value.substr(1);//remove leading ' or "
         if(!last)
         {// Make sure we read till the end of the quoted string
            const string::size_type loc = value.find(delim, 0 );
            if( loc == string::npos )
            {
               char c;
               while(sst.peek()!=delim){sst.get(c); value+=c;}
               sst.get(c);
            }
         }
         // Remove trailing spaces
         while(isgraph(value.at(value.size()-1))==0)
            value=value.substr(0,value.size()-1);
         value=value.substr(0,value.size()-1);
         //cout<<"QuotedString:"<<value<<endl;
      }
   return value;
}

void CIF::Parse()
{
   list<string>::const_iterator pos=mvLine.begin();
   string block="";
   while(pos!=mvLine.end())
   {
      //cout<<*pos<<endl;
      if(pos->size()==0)
      {
         ++pos;
         continue;
      }
      if(pos->at(0)=='#')
      {
         //cout<<"Comment :"<<*pos<<endl;
         if(block=="") mvComment.push_back(*pos++);
         else mvData[block].mvComment.push_back(*pos++);
         continue;
      }
      if(ci_string(pos->substr(0,5).c_str())==ci_string("data_"))
      {
         block=pos++->substr(5);
         //cout<<"NEW BLOCK DATA: "<<block<<endl;
         mvData[block]=CIFData();
         continue;
      }
      if(pos->substr(0,1)=="_")
      {
         stringstream sst;
         sst<<*pos++;
         string tag,value;
         sst>>tag;
         value=CIFReadValue(sst,pos,true);
         mvData[block].mvItem[ci_string(tag.c_str())]=value;
         //cout<<"New Tag:"<<tag<<" ("<<value.size()<<"):"<<value<<endl;
         continue;
      }
      if(ci_string(pos->substr(0,5).c_str())=="loop_")
      {
         vector<ci_string> tit;
         string tmp;
         //cout<<"LOOP: ";
         ++pos;
         {//extract string while ignoring leading and trailing spaces
            stringstream sst;
            sst<<*pos;
            sst>>tmp;
         }
         while(tmp.at(0)=='_')
         {
            tit.push_back(ci_string(tmp.c_str()));
            //cout<<" : "<<tmp;
            stringstream sst;
            sst<<*++pos;//increment before so that the first value after the titles is not missed
            sst>>tmp;
         }
         //cout<<endl;
         map<ci_string,vector<string> > lp;
         while(true)
         {
            stringstream sst;
            //cout<<"LOOP VALUES...: ";
            sst<<*pos++;
            for(unsigned int i=0;i<tit.size();++i)
            {
               const string value=CIFReadValue(sst,pos,i==(tit.size()-1));
               lp[tit[i]].push_back(value);
               //cout<<value<<"  ";
            }
            //cout<<endl<<*pos<<endl;
            if(pos==mvLine.end()) break;
            if(pos->size()==0) break;
            if(pos->at(0)!=' ') break;
         }
         // The key to the mvLoop map is the set of column titles
         set<ci_string> stit;
         for(unsigned int i=0;i<tit.size();++i) stit.insert(tit[i]);
         mvData[block].mvLoop[stit]=lp;
         continue;
      }
      // We should never get here !
      ++pos;
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

bool CreateCrystalFromCIF(std::istream &in)
{
   ObjCryst::CIF cif(in,true,true);
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
         return true;
      }
   return false;
}

}//namespace
