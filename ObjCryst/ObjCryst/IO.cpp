/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2006 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

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
*  source file for XMLInput/XMLOutput in ObjCryst++
*
*/

#include <stdio.h>//for sprintf

#include "ObjCryst/General.h"
#include "ObjCryst/IO.h"
#include "RefinableObj/IO.h"
#include "RefinableObj/GlobalOptimObj.h"
//#include "ObjCryst/SpaceGroup.h"
#include "ObjCryst/Scatterer.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/ZScatterer.h"
//#include "ObjCryst/ScatteringData.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/ScatteringPowerSphere.h"
#include "ObjCryst/Atom.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/PowderPattern.h"
#include "Quirks/VFNStreamFormat.h"
#include "ObjCryst/Molecule.h"

#include <iostream>
#include <fstream>
#include <sstream>

//#define USE_BACKGROUND_MAXLIKE_ERROR

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    Global functions
//
////////////////////////////////////////////////////////////////////////

void XMLCrystFileSaveGlobal(const string & filename)
{
   VFN_DEBUG_ENTRY("XMLCrystFileSaveGlobal(filename)",5)
   
   ofstream out(filename.c_str());
   if(!out){};//:TODO:
   XMLCrystFileSaveGlobal(out);
   out.close();
   VFN_DEBUG_EXIT("XMLCrystFileSaveGlobal(filename):End",5)
}

void XMLCrystFileSaveGlobal(ostream &out)
{
   VFN_DEBUG_ENTRY("XMLCrystFileSaveGlobal(ostream)",5)  
   XMLCrystTag tag("ObjCryst");
   time_t date=time(0);
   char strDate[40];
   strftime(strDate,sizeof(strDate),"%Y-%m-%dT%H:%M:%S%Z",gmtime(&date));//%Y-%m-%dT%H:%M:%S%Z
   tag.AddAttribute("Date",strDate);
   tag.AddAttribute("Revision","918");
   out<<tag<<endl;
   
   for(int i=0;i<gCrystalRegistry.GetNb();i++)
      gCrystalRegistry.GetObj(i).XMLOutput(out,1);
   
   for(int i=0;i<gDiffractionDataSingleCrystalRegistry.GetNb();i++)
      gDiffractionDataSingleCrystalRegistry.GetObj(i).XMLOutput(out,1);
   
   for(int i=0;i<gPowderPatternRegistry.GetNb();i++)
      gPowderPatternRegistry.GetObj(i).XMLOutput(out,1);
   
   for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
      gOptimizationObjRegistry.GetObj(i).XMLOutput(out,1);
   
   tag.SetIsEndTag(true);
   out<<tag;
   VFN_DEBUG_EXIT("XMLCrystFileSaveGlobal(ostream)",5)
}

ObjRegistry<XMLCrystTag> XMLCrystFileLoadObjectList(const string & filename)
{
   VFN_DEBUG_ENTRY("XMLCrystFileLoadObjectList(filename)",5)

   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   ObjRegistry<XMLCrystTag> reg;
   for(;;)
   {
      XMLCrystTag *pTag =new XMLCrystTag (is);
      if(true==is.eof())
      {
         VFN_DEBUG_EXIT("XMLCrystFileLoadObjectList(filename):End",5)
         for(int i=0;i<reg.GetNb();i++) reg.GetObj(i).Print();
         is.close();
         return reg;
      }
      //pTag->Print();
      if(("Crystal"==pTag->GetName()||
           "DiffractionDataSingleCrystal"==pTag->GetName()||
           "PowderPattern"==pTag->GetName()||
           "GlobalOptimObj"==pTag->GetName())
          && !(pTag->IsEndTag())) reg.Register(*pTag);
      else delete pTag;
   }
   return reg;
}

template<class T> void XMLCrystFileLoadObject(const string & filename,
                                              const string &tagName,
                                              const string &name, T*obj)
{
   VFN_DEBUG_ENTRY("XMLCrystFileLoadObject(filename,IOCrystTag,T&)",5)

   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   XMLCrystTag tag;
   while(true)
   {
      is>>tag;
      if(true==is.eof())
      {
         cout<<"XMLCrystFileLoadObject(filename,IOCrystTag,T&):Not Found !"<<endl;
         is.close();
         obj=0;
         VFN_DEBUG_EXIT("XMLCrystFileLoadObject(filename,IOCrystTag,T&)",5)
         return;
      }
      if(tagName!=tag.GetName())continue;
      for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         if("Name"==tag.GetAttributeName(i))
            if(name==tag.GetAttributeValue(i)) break;
   }
   VFN_DEBUG_MESSAGE("XMLCrystFileLoadObject(filename,IOCrystTag,T&):Found"<<tag,5)
   obj = new T;
   obj->XMLInput(is,tag);
   is.close();
   VFN_DEBUG_EXIT("XMLCrystFileLoadObject(filename,IOCrystTag,T&)",5)
}

template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     Crystal*);
template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     PowderPattern*);
template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     DiffractionDataSingleCrystal*);
//template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
//                                    ZScatterer*);
template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     PowderPatternBackground*);
template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     PowderPatternDiffraction*);
template void XMLCrystFileLoadObject(const string & ,const string &,const string &,
                                     MonteCarloObj*);

void XMLCrystFileLoadAllObject(const string & filename)
{
   VFN_DEBUG_ENTRY("XMLCrystFileLoadAllObject(filename,)",5)
   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   XMLCrystFileLoadAllObject(is);
   VFN_DEBUG_EXIT("XMLCrystFileLoadAllObject(filename,)",5)
}
void XMLCrystFileLoadAllObject(istream &is)
{
   VFN_DEBUG_ENTRY("XMLCrystFileLoadAllObject(istream)",5)
   XMLCrystTag tag;
   do {is>>tag;} while("ObjCryst"!=tag.GetName());
   
   while(true)
   {
      XMLCrystTag tag(is);
      if(true==is.eof()) break;
      if(tag.GetName()=="Crystal")
      {
         Crystal* obj = new Crystal;
         obj->XMLInput(is,tag);
      }
      if(tag.GetName()=="PowderPattern")
      {
         PowderPattern* obj = new PowderPattern;
         obj->XMLInput(is,tag);
      }
      if(tag.GetName()=="DiffractionDataSingleCrystal")
      {
         DiffractionDataSingleCrystal* obj = new DiffractionDataSingleCrystal;
         obj->XMLInput(is,tag);
      }
      if(tag.GetName()=="GlobalOptimObj")
      {
         MonteCarloObj* obj = new MonteCarloObj;
         obj->XMLInput(is,tag);
      }
   }
   VFN_DEBUG_EXIT("XMLCrystFileLoadAllObject(istream)",5)
}
////////////////////////////////////////////////////////////////////////
//
//    I/O ScatteringPowerAtom
//
////////////////////////////////////////////////////////////////////////
void ScatteringPowerAtom::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ScatteringPowerAtom::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ScatteringPowerAtom");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("Symbol",mSymbol);
   os <<tag<<endl;
   
   if(true==this->mIsIsotropic)
      this->GetPar(&mBiso).XMLOutput(os,"Biso",indent+1);
   os<<endl;

   this->GetPar("ML Error").XMLOutput(os,"ML Error",indent+1);
   os <<endl;

   this->GetPar("ML-Nb Ghost Atoms").XMLOutput(os,"ML-NbGhost",indent+1);
   os <<endl;
   
   this->GetPar("Formal Charge").XMLOutput(os,"Formal Charge",indent+1);
   os <<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
   XMLCrystTag tag2("RGBColour");
   os << tag2
      << mColourRGB[0]<<" "
      << mColourRGB[1]<<" "
      << mColourRGB[2];
   tag2.SetIsEndTag(true);
   os << tag2<<endl;
   
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ScatteringPowerAtom::XMLOutput():"<<this->GetName(),5)
}

void ScatteringPowerAtom::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ScatteringPowerAtom::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Symbol"==tagg.GetAttributeName(i)) mSymbol=tagg.GetAttributeValue(i);
   }
   this->Init(mName,mSymbol,mBiso);
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ScatteringPowerAtom"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("ScatteringPowerAtom::Exit():"<<this->GetName(),5)
         return;
      }
      if("RGBColour"==tag.GetName())
      {
         float r,g,b;
         is>>r>>g>>b;
         this->SetColour(r,g,b);
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("Biso"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mBiso).XMLInput(is,tag);
                  break;
               }
               if("ML Error"==tag.GetAttributeValue(i))
               {
                  this->GetPar("ML Error").XMLInput(is,tag);
                  break;
               }
               if("ML-NbGhost"==tag.GetAttributeValue(i))
               {
                  this->GetPar("ML-Nb Ghost Atoms").XMLInput(is,tag);
                  break;
               }
               if("Formal Charge"==tag.GetAttributeValue(i))
               {
                  this->GetPar("Formal Charge").XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O Atom
//
////////////////////////////////////////////////////////////////////////
void Atom::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Atom::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Atom");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("ScattPow",mpScattPowAtom->GetName());
   os <<tag;
   os <<endl;
   indent++;
   
   this->GetPar(mXYZ.data()+0).XMLOutput(os,"x",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+1).XMLOutput(os,"y",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+2).XMLOutput(os,"z",indent);
   os <<endl;
   
   this->GetPar(&mOccupancy).XMLOutput(os,"Occup",indent);
   os <<endl;
   
   tag.SetIsEndTag(true);
   indent--;
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   
   VFN_DEBUG_EXIT("Atom::XMLOutput():"<<this->GetName(),5)
}

void Atom::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("Atom::XMLInput():"<<this->GetName(),5)
   string scattPowName;
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("ScattPow"==tagg.GetAttributeName(i)) scattPowName=tagg.GetAttributeValue(i);
   }
   const ScatteringPower* scattPow=
         &(this->GetCrystal().GetScatteringPowerRegistry().GetObj(scattPowName));
   VFN_DEBUG_MESSAGE("Found Scattering Power:"<< scattPowName<<" at "<<scattPow,4);
   this->Init(0,0,0,mName,scattPow,1);
   while(true)
   {
      XMLCrystTag tag(is);
      if(("Atom"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("Atom::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("x"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+0).XMLInput(is,tag);
                  break;
               }
               if("y"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+1).XMLInput(is,tag);
                  break;
               }
               if("z"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+2).XMLInput(is,tag);
                  break;
               }
               if("Occup"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mOccupancy).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O ZAtom
//
////////////////////////////////////////////////////////////////////////
void ZAtom::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ZAtom::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ZAtom");
   tag.AddAttribute("Name",mName);
   if(0!=this->GetScatteringPower())//else it is a dummy atom
      tag.AddAttribute("ScattPow",this->GetScatteringPower()->GetName());
      
   tag.AddAttribute("BondAtom",this->GetZScatterer()
                                       .GetZAtomRegistry()
                                          .GetObj(this->GetZBondAtom())
                                             .GetName());
   tag.AddAttribute("AngleAtom",this->GetZScatterer()
                                       .GetZAtomRegistry()
                                          .GetObj(this->GetZAngleAtom())
                                             .GetName());
   tag.AddAttribute("DihedAtom",this->GetZScatterer()
                                       .GetZAtomRegistry()
                                          .GetObj(this->GetZDihedralAngleAtom())
                                             .GetName());
   os <<tag<<endl;
   indent++;
   
   
   this->GetZScatterer().GetPar(&mBondLength).XMLOutput(os,"BondLength",indent);
   os <<endl;
   
   this->GetZScatterer().GetPar(&mAngle).XMLOutput(os,"Angle",indent);
   os <<endl;
   
   this->GetZScatterer().GetPar(&mDihed).XMLOutput(os,"DihedAng",indent);
   os <<endl;
   
   this->GetZScatterer().GetPar(&mOccupancy).XMLOutput(os,"Occup",indent);
   os <<endl;
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ZAtom::XMLOutput():"<<this->GetName(),5)
}

void ZAtom::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ZAtom::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i))
      {
         this->SetName(tagg.GetAttributeValue(i));
         continue;
      }
      if("ScattPow"==tagg.GetAttributeName(i))
      {
         const ScatteringPower* scattPow=&(this->GetZScatterer()
                                                   .GetCrystal()
                                                      .GetScatteringPowerRegistry()
                                                         .GetObj(tagg.GetAttributeValue(i)));
         this->SetScatteringPower(scattPow);
         continue;
      }
      if("BondAtom"==tagg.GetAttributeName(i))
      {
         mAtomBond=this->GetZScatterer().GetZAtomRegistry().Find(tagg.GetAttributeValue(i));
         continue;
      }
      if("AngleAtom"==tagg.GetAttributeName(i))
      {
         mAtomAngle=this->GetZScatterer().GetZAtomRegistry().Find(tagg.GetAttributeValue(i));
         continue;
      }
      if("DihedAtom"==tagg.GetAttributeName(i))
      {
         mAtomDihed=this->GetZScatterer().GetZAtomRegistry().Find(tagg.GetAttributeValue(i));
         continue;
      }
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ZAtom"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("ZAtom::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("BondLength"==tag.GetAttributeValue(i))
               {
                  this->GetZScatterer().GetPar(&mBondLength).XMLInput(is,tag);
                  break;
               }
               if("Angle"==tag.GetAttributeValue(i))
               {
                  this->GetZScatterer().GetPar(&mAngle).XMLInput(is,tag);
                  break;
               }
               if("DihedAng"==tag.GetAttributeValue(i))
               {
                  this->GetZScatterer().GetPar(&mDihed).XMLInput(is,tag);
                  break;
               }
               if("Occup"==tag.GetAttributeValue(i))
               {
                  this->GetZScatterer().GetPar(&mOccupancy).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
   }
}

////////////////////////////////////////////////////////////////////////
//
//    I/O ZScatterer
//
////////////////////////////////////////////////////////////////////////
void ZScatterer::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ZScatterer::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ZScatterer");
   tag.AddAttribute("Name",mName);
   os <<tag<<endl;
   indent++;
   
   this->GetPar(mXYZ.data()+0).XMLOutput(os,"x",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+1).XMLOutput(os,"y",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+2).XMLOutput(os,"z",indent);
   os <<endl;
   
   this->GetPar(&mOccupancy).XMLOutput(os,"Occup",indent);
   os <<endl;
   
   this->GetPar(&mPhi).XMLOutput(os,"Phi",indent);
   os <<endl;
   
   this->GetPar(&mChi).XMLOutput(os,"Chi",indent);
   os <<endl;
   
   this->GetPar(&mPsi).XMLOutput(os,"Psi",indent);
   os <<endl;
   
   for(int i=0;i<mZAtomRegistry.GetNb();i++) mZAtomRegistry.GetObj(i).XMLOutput(os,indent);
   
   if(mZAtomRegistry.GetNb()>0)
   {
      for(int i=0;i<=indent;i++) os << "  " ;
      XMLCrystTag tag2("PivotAtom",false,true);
      tag2.AddAttribute("Name",this->GetZAtomRegistry().GetObj(mCenterAtomIndex).GetName());
      os <<tag2<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ZScatterer::XMLOutput():"<<this->GetName(),5)
}

void ZScatterer::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ZScatterer::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ZScatterer"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("ZScatterer::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("x"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+0).XMLInput(is,tag);
                  break;
               }
               if("y"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+1).XMLInput(is,tag);
                  break;
               }
               if("z"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+2).XMLInput(is,tag);
                  break;
               }
               if("Occup"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mOccupancy).XMLInput(is,tag);
                  break;
               }
               if("Phi"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPhi).XMLInput(is,tag);
                  break;
               }
               if("Chi"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mChi).XMLInput(is,tag);
                  break;
               }
               if("Psi"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPsi).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
      if("ZAtom"==tag.GetName())
      {
         //we must take care of possible dummy atoms
         const ScatteringPower* scattPow=0;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("ScattPow"==tag.GetAttributeName(i))
               scattPow=&(this->GetCrystal().GetScatteringPowerRegistry()
                                                .GetObj(tag.GetAttributeValue(i)));
         const long nb=mZAtomRegistry.GetNb();
         this->AddAtom("",scattPow,0,0,0,0,0,1);
         mZAtomRegistry.GetObj(nb).XMLInput(is,tag);
         // Update the name of refinable parameters
         {
            char buf [20];
            sprintf(buf,"%d-%d",(int)nb,(int)(mZAtomRegistry.GetObj(nb).GetZBondAtom()));
            this->GetPar(&(mZAtomRegistry.GetObj(nb).mBondLength))
               .SetName("Length"+(string)buf);
               
            sprintf(buf,"%d-%d-%d",(int)nb,(int)(mZAtomRegistry.GetObj(nb).GetZBondAtom()),
                                   (int)(mZAtomRegistry.GetObj(nb).GetZAngleAtom()));
            this->GetPar(&(mZAtomRegistry.GetObj(nb).mAngle))
               .SetName("Angle"+(string)buf);
               
            sprintf(buf,"%d-%d-%d-%d",(int)nb,(int)(mZAtomRegistry.GetObj(nb).GetZBondAtom()),
                                      (int)(mZAtomRegistry.GetObj(nb).GetZAngleAtom()),
                                      (int)(mZAtomRegistry.GetObj(nb).GetZDihedralAngleAtom()));
            this->GetPar(&(mZAtomRegistry.GetObj(nb).mDihed))
               .SetName("Dihed"+(string)buf);
         }
      }
      if("PivotAtom"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
            {
               mCenterAtomIndex=this->GetZAtomRegistry().Find(tag.GetAttributeValue(i));
            }
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O Crystal
//
////////////////////////////////////////////////////////////////////////
void Crystal::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Crystal::XMLOutput():"<<this->GetName(),5)
   
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Crystal");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("SpaceGroup",this->GetSpaceGroup().GetName());
   os <<tag<<endl;
   indent++;
   
   // :TODO: 
   this->GetPar("a").XMLOutput(os,"a",indent);
   os <<endl;
   
   this->GetPar("b").XMLOutput(os,"b",indent);
   os <<endl;
   
   this->GetPar("c").XMLOutput(os,"c",indent);
   os <<endl;
   
   this->GetPar("alpha").XMLOutput(os,"alpha",indent);
   os <<endl;
   
   this->GetPar("beta").XMLOutput(os,"beta",indent);
   os <<endl;
   
   this->GetPar("gamma").XMLOutput(os,"gamma",indent);
   os <<endl;
   
   for(unsigned int i=0;i<this->GetNbOption();i++)
   {
      this->GetOption(i).XMLOutput(os,indent);
      os <<endl<<endl;
   }
   
   for(int i=0;i<mScatteringPowerRegistry.GetNb();i++) 
      mScatteringPowerRegistry.GetObj(i).XMLOutput(os,indent);
   os <<endl;
   for(int i=0;i<mScattererRegistry.GetNb();i++) 
      mScattererRegistry.GetObj(i).XMLOutput(os,indent);
   os <<endl;

   if(mvBumpMergePar.size()>0)
   {
      VBumpMergePar::const_iterator pos;
      for(pos=mvBumpMergePar.begin();pos!=mvBumpMergePar.end();pos++)
      {
         for(int k=0;k<=indent;k++) os << "  " ;
         XMLCrystTag tagBump("AntiBumpDistance");
         tagBump.AddAttribute("ScattPow1",pos->first.first->GetName());
         tagBump.AddAttribute("ScattPow2",pos->first.second->GetName());
         {
            stringstream ss;
            ss << pos->second.mCanOverlap;
            tagBump.AddAttribute("AllowMerge",ss.str());
         }
         os<<tagBump;
         tagBump.SetIsEndTag(true);
         os<<sqrt(pos->second.mDist2)<<tagBump<<endl;
      }
      for(int k=0;k<=indent;k++) os << "  " ;
      XMLCrystTag tag2("AntiBumpScale");
      os << tag2<< mBumpMergeScale;
      tag2.SetIsEndTag(true);
      os << tag2<<endl;
   }
   if(mvBondValenceRo.size()>0)
   {
      map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>::const_iterator pos;
      for(pos=mvBondValenceRo.begin();pos!=mvBondValenceRo.end();pos++)
      {
         for(int k=0;k<=indent;k++) os << "  " ;
         XMLCrystTag tagBVRo("BondValenceRo");
         tagBVRo.AddAttribute("ScattPow1",pos->first.first->GetName());
         tagBVRo.AddAttribute("ScattPow2",pos->first.second->GetName());
         os<<tagBVRo;
         tagBVRo.SetIsEndTag(true);
         os<<pos->second<<tagBVRo<<endl;
      }
      for(int k=0;k<=indent;k++) os << "  " ;
      XMLCrystTag tag2("BondValenceCostScale");
      os << tag2<< mBondValenceCostScale;
      tag2.SetIsEndTag(true);
      os << tag2<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("Crystal::XMLOutput():"<<this->GetName(),5)
}

void Crystal::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("Crystal::XMLInput():"<<this->GetName(),5)
   //Remove Scatterers and Scattering Powers
      for(long i=0;i<mScatteringPowerRegistry.GetNb();i++)
      {
         this->RemoveSubRefObj(mScatteringPowerRegistry.GetObj(i));
         mScatteringPowerRegistry.GetObj(i).DeRegisterClient(*this);
      }
      mScatteringPowerRegistry.DeleteAll();
      for(long i=0;i<mScattererRegistry.GetNb();i++)
      {
         this->RemoveSubRefObj(mScattererRegistry.GetObj(i));
         mScattererRegistry.GetObj(i).DeRegisterClient(*this);
      }
      mScattererRegistry.DeleteAll();

   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("SpaceGroup"==tagg.GetAttributeName(i))
          this->Init(1,2,3,M_PI/2,M_PI/2,M_PI/2,tagg.GetAttributeValue(i),this->GetName());
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("Crystal"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("Crystal::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               this->GetPar(tag.GetAttributeValue(i)).XMLInput(is,tag);
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         this->InitRefParList();// Fix the "used" tag of refinable par after options
         continue;
      }
      if("AntiBumpDistance"==tag.GetName())
      {
         float dist;
         bool useMerge=false;
         bool allowMerge;
         string scattPow1;
         string scattPow2;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("AllowMerge"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss >>allowMerge;
               useMerge=true;
               continue;
            }
            if("ScattPow1"==tag.GetAttributeName(i)) scattPow1=tag.GetAttributeValue(i);
            if("ScattPow2"==tag.GetAttributeName(i)) scattPow2=tag.GetAttributeValue(i);
         }
         is>>dist;
         XMLCrystTag junk(is);//end tag
         if(useMerge)
            this->SetBumpMergeDistance(mScatteringPowerRegistry.GetObj(scattPow1),
                                    mScatteringPowerRegistry.GetObj(scattPow2),
                                    dist,allowMerge);
         else this->SetBumpMergeDistance(mScatteringPowerRegistry.GetObj(scattPow1),
                                    mScatteringPowerRegistry.GetObj(scattPow2),
                                    dist);
         continue;
      }
      if("BondValenceRo"==tag.GetName())
      {
         float ro;
         string scattPow1;
         string scattPow2;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("ScattPow1"==tag.GetAttributeName(i)) scattPow1=tag.GetAttributeValue(i);
            if("ScattPow2"==tag.GetAttributeName(i)) scattPow2=tag.GetAttributeValue(i);
         }
         is>>ro;
         XMLCrystTag junk(is);//end tag
         this->AddBondValenceRo(mScatteringPowerRegistry.GetObj(scattPow1),
                                mScatteringPowerRegistry.GetObj(scattPow2),ro);
         continue;
      }
      if("AntiBumpScale"==tag.GetName())
      {
         is>>mBumpMergeScale;
         XMLCrystTag junk(is);
      }
      if("BondValenceCostScale"==tag.GetName())
      {
         is>>mBondValenceCostScale;
         XMLCrystTag junk(is);
      }
      if("Atom"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Crystal::XMLInput():reading an Atom",5)
         Atom *at=new Atom;
         at->SetCrystal(*this);
         at->XMLInput(is,tag);
         this->AddScatterer(at);
         VFN_DEBUG_EXIT("Crystal::XMLInput():reading an Atom",5)
         continue;
      }
      if("ScatteringPowerAtom"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Crystal::XMLInput():reading a ScatteringPowerAtom",5)
         VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a ScatteringPowerAtom",5)
         ScatteringPowerAtom *sc=new ScatteringPowerAtom;
         sc->XMLInput(is,tag);
         this->AddScatteringPower(sc);
         VFN_DEBUG_EXIT("Crystal::XMLInput():reading a ScatteringPowerAtom",5)
         continue;
      }
      if("ScatteringPowerSphere"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Crystal::XMLInput():reading a ScatteringPowerSphere",5)
         VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a ScatteringPowerSphere",5)
         ScatteringPowerSphere *sc=new ScatteringPowerSphere;
         sc->XMLInput(is,tag);
         this->AddScatteringPower(sc);
         VFN_DEBUG_EXIT("Crystal::XMLInput():reading a ScatteringPowerSphere",5)
         continue;
      }
      if("ZScatterer"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Crystal::XMLInput():reading a ZScatterer",5)
         VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a ZScatterer",5)
         ZScatterer *z=new ZScatterer("",*this);
         z->XMLInput(is,tag);
         this->AddScatterer(z);
         VFN_DEBUG_EXIT("Crystal::XMLInput():reading a ZScatterer",5)
         continue;
      }
      if("Molecule"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Crystal::XMLInput():reading a Molecule",5)
         VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a Molecule",5)
         Molecule *z=new Molecule(*this,"");
         z->XMLInput(is,tag);
         this->AddScatterer(z);
         VFN_DEBUG_EXIT("Crystal::XMLInput():reading a Molecule",5)
         continue;
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O Radiation
//
////////////////////////////////////////////////////////////////////////
void Radiation::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Radiation::XMLOutput():"<<this->GetName(),5)
   XMLCrystTag tag("Radiation");
   if(WAVELENGTH_ALPHA12==this->GetWavelengthType())
      tag.AddAttribute("XRayTube",mXRayTubeName);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   indent++;

   mRadiationType.XMLOutput(os,indent);
   os<<endl;
   
   mWavelengthType.XMLOutput(os,indent);
   os<<endl;
   
   for(int i=0;i<indent;i++) os << "  " ;
   {
      XMLCrystTag tag2("LinearPolarRate");
      os << tag2<< mLinearPolarRate;
      tag2.SetIsEndTag(true);
      os << tag2<<endl;
   }
   
   if(WAVELENGTH_ALPHA12==this->GetWavelengthType())
   {
      for(int i=0;i<indent;i++) os << "  " ;
      {
         XMLCrystTag tag2("XRayTubeDeltaLambda");
         os << tag2<< mXRayTubeDeltaLambda;
         tag2.SetIsEndTag(true);
         os << tag2<<endl;
      }
      for(int i=0;i<indent;i++) os << "  " ;
      {
         XMLCrystTag tag2("XRayTubeAlpha2Alpha1Ratio");
         os << tag2<< mXRayTubeAlpha2Alpha1Ratio;
         tag2.SetIsEndTag(true);
         os << tag2<<endl;
      }
   }

   switch(this->GetWavelengthType())
   {
      case WAVELENGTH_MONOCHROMATIC: this->GetPar(mWavelength.data()).XMLOutput(os,indent);break;
      case WAVELENGTH_ALPHA12:
      {
         this->GetPar(mWavelength.data()).XMLOutput(os,"Wavelength",indent);
         break;
      }
      case WAVELENGTH_TOF:break;
      default: throw ObjCrystException("This radiation is not implemented !!");
   }
   os<<endl;
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag;
   
   VFN_DEBUG_EXIT("Radiation::XMLOutput():"<<this->GetName(),5)
}

void Radiation::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("Radiation::XMLInput():"<<this->GetName(),5)
   string scattPowName;
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("XRayTube"==tagg.GetAttributeName(i))
         if(tagg.GetAttributeValue(i)!="") // Something went wrong !
            this->SetWavelength(tagg.GetAttributeValue(i));
   }
   
   while(true)
   {
      XMLCrystTag tag(is);
      if(("Radiation"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("Radiation::Exit():"<<this->GetName(),5)
         return;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
            {
               if("Radiation"==tag.GetAttributeValue(i)) mRadiationType.XMLInput(is,tag);
               if("Spectrum"==tag.GetAttributeValue(i)) mWavelengthType.XMLInput(is,tag);
            }
      }
      if("LinearPolarRate"==tag.GetName())
      {
         is>>mLinearPolarRate;
         XMLCrystTag junk(is);
      }
      if("XRayTubeDeltaLambda"==tag.GetName())
      {
         is>>mXRayTubeDeltaLambda;
         XMLCrystTag junk(is);
      }
      if("XRayTubeAlpha2Alpha1Ratio"==tag.GetName())
      {
         is>>mXRayTubeAlpha2Alpha1Ratio;
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("Wavelength"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mWavelength.data()).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O DiffractionDataSingleCrystal
//
////////////////////////////////////////////////////////////////////////
void DiffractionDataSingleCrystal::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("DiffractionDataSingleCrystal");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("Crystal",this->GetCrystal().GetName());
   os <<tag<<endl;
   indent++;
   
   this->GetPar("Scale factor").XMLOutput(os,"Scale factor",indent);
   os <<endl;
   
   mRadiation.XMLOutput(os,indent);
   os <<endl;

   this->GetPar(&mGlobalBiso).XMLOutput(os,"globalBiso",indent);
   os <<endl;
   
   mGroupOption.XMLOutput(os,indent);
   os <<endl;
   
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag2("MaxSinThetaOvLambda");
   os << tag2<< mMaxSinThetaOvLambda;
   tag2.SetIsEndTag(true);
   os << tag2<<endl<<endl;
   
   if(mGroupOption.GetChoice()!=2)
   {
      XMLCrystTag tag3("HKLIobsSigmaWeightList");
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag3<<endl;
      
      for(long j=0;j<this->GetNbRefl();j++)
      {
         for(int i=0;i<=indent;i++) os << "  " ;
         os << mIntH(j) <<" "
            << mIntK(j) <<" "
            << mIntL(j) <<" "
            << mObsIntensity(j) <<" "
            << mObsSigma(j) <<" "
            << mWeight(j) <<" "
            <<endl;
      }
      
      tag3.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag3<<endl;
   }
   else
   {
      XMLCrystTag tag3("HKLIobsSigmaWeightGROUPList");
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag3<<endl;
      
      long first=0;
      for(long j=0;j<mNbGroup;j++)
      {
         XMLCrystTag tag4("HKLGroup");
         {
            stringstream s;
            s<<mGroupIobs(j);
            tag4.AddAttribute("Iobs",s.str());
         }
         {
            stringstream s;
            s<<mGroupSigma(j);
            tag4.AddAttribute("IobsSigma",s.str());
         }
         {
            stringstream s;
            s<<mGroupWeight(j);
            tag4.AddAttribute("Weight",s.str());
         }
         for(int i=0;i<=indent;i++) os << "  " ;
         os<<tag4<<endl;
         for(long k=first;k<mGroupIndex(j);k++)
         {
            for(int i=0;i<=indent;i++) os << "  " ;
            os << mIntH(k) <<" "<< mIntK(k) <<" "<< mIntL(k) <<" "<<endl;
         }
         for(int i=0;i<=indent;i++) os << "  " ;
         tag4.SetIsEndTag(true);
         os<<tag4<<endl;
         first=mGroupIndex(j);
      }
      
      tag3.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag3<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::XMLOutput():"<<this->GetName(),5)
}

void DiffractionDataSingleCrystal::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Crystal"==tagg.GetAttributeName(i)) 
         this->SetCrystal(gCrystalRegistry.GetObj(tagg.GetAttributeValue(i)));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("DiffractionDataSingleCrystal"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::XMLInput():"<<this->GetName(),5)
         return;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
            {
               string name=tag.GetAttributeValue(i);
               if(name=="Twinning correction") name="Group Reflections";
               mOptionRegistry.GetObj(name).XMLInput(is,tag);
            }
         continue;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("Scale factor"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mScaleFactor).XMLInput(is,tag);
                  break;
               }
            }
         }
      }
      if("Radiation"==tag.GetName()) mRadiation.XMLInput(is,tag);
      if("MaxSinThetaOvLambda"==tag.GetName())
      {
         is>>mMaxSinThetaOvLambda;
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("globalBiso"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mGlobalBiso).XMLInput(is,tag);
                  break;
               }
            }
         }
      }
      if("HKLIobsSigmaWeightList"==tag.GetName())
      {
         long nbrefl=0;
         CrystVector_long h(100),k(100),l(100);
         CrystVector_REAL iobs(100),sigma(100),weight(100);
         do
         {
            is >>h(nbrefl)>>k(nbrefl)>>l(nbrefl)
               >>iobs(nbrefl)>>sigma(nbrefl)>>weight(nbrefl);
            nbrefl++;
            if(nbrefl==iobs.numElements())
            {
               h.resizeAndPreserve(nbrefl+100);
               k.resizeAndPreserve(nbrefl+100);
               l.resizeAndPreserve(nbrefl+100);
               iobs.resizeAndPreserve(nbrefl+100);
               sigma.resizeAndPreserve(nbrefl+100);
               weight.resizeAndPreserve(nbrefl+100);
            }
            while(0==isgraph(is.peek())) is.get();
            //cout << is.peek()<<" "<<nbrefl<<endl;
         }
         while(is.peek()!='<');//until next tag
         XMLCrystTag junkEndTag(is);
         
         h.resizeAndPreserve(nbrefl);
         k.resizeAndPreserve(nbrefl);
         l.resizeAndPreserve(nbrefl);
         iobs.resizeAndPreserve(nbrefl);
         sigma.resizeAndPreserve(nbrefl);
         weight.resizeAndPreserve(nbrefl);
         this->SetHklIobs(h,k,l,iobs,sigma);
         this->SetWeight(weight);
         this->SortReflectionBySinThetaOverLambda();
         this->CalcIcalc();
         this->FitScaleFactorForRw();
      }
      if("HKLIobsSigmaWeightGROUPList"==tag.GetName())
      {
         mNbRefl=0;
         mNbGroup=0;
         // This must NOT be changed with this kind of data.
         mGroupOption.SetChoice(2);
         // So de-register the option so that it is hidden from the user's view
         mOptionRegistry.DeRegister(mGroupOption);
         mClockMaster.RemoveChild(mGroupOption.GetClock());
         mH.resize(500);
         mK.resize(500);
         mL.resize(500);
         mObsIntensity.resize(500);
         mObsSigma.resize(500);
         mGroupIndex.resize(500);
         mGroupIobs.resize(500);
         mGroupSigma.resize(500);
         mGroupWeight.resize(500);
         while(true)
         {
            XMLCrystTag grouptag(is);
            if(grouptag.GetName()=="HKLIobsSigmaWeightGROUPList") break;
            if(grouptag.GetName()=="HKLGroup")
            {
               for(unsigned int i=0;i<grouptag.GetNbAttribute();++i)
               {
                  if(grouptag.GetAttributeName(i)=="Iobs")
                  {
                     stringstream sst;
                     sst<<grouptag.GetAttributeValue(i);
                     sst>>mGroupIobs(mNbGroup);
                     continue;
                  }
                  if(grouptag.GetAttributeName(i)=="IobsSigma")
                  {
                     stringstream sst;
                     sst<<grouptag.GetAttributeValue(i);
                     sst>>mGroupSigma(mNbGroup);
                     continue;
                  }
                  if(grouptag.GetAttributeName(i)=="Weight")
                  {
                     stringstream sst;
                     sst<<grouptag.GetAttributeValue(i);
                     sst>>mGroupWeight(mNbGroup);
                     continue;
                  }
               }
               VFN_DEBUG_MESSAGE("Group #"<<mNbGroup<<" ,Iobs="<<mGroupIobs(mNbGroup)<<" ,Sigma="<<mGroupSigma(mNbGroup)<<" ,Weight="<<mGroupWeight(mNbGroup),2)
               do
               {
                  is >>mH(mNbRefl)>>mK(mNbRefl)>>mL(mNbRefl);
                  VFN_DEBUG_MESSAGE("         "<<mH(mNbRefl)<<" "<<mK(mNbRefl)<<" "<<mL(mNbRefl),2)
                  mGroupIndex(mNbRefl)=mNbGroup;
                  mNbRefl++;
                  if(mNbRefl==mH.numElements())
                  {
                     mH.resizeAndPreserve(mNbRefl+500);
                     mK.resizeAndPreserve(mNbRefl+500);
                     mL.resizeAndPreserve(mNbRefl+500);
                     mObsIntensity.resizeAndPreserve(mNbRefl+500);
                     mObsSigma.resizeAndPreserve(mNbRefl+500);
                     mGroupIndex.resizeAndPreserve(mNbRefl+500);
                  }
                  while(0==isgraph(is.peek())) is.get();
               }
               while(is.peek()!='<');//until end tag
               XMLCrystTag junkEndTag(is);
               if(++mNbGroup==mGroupIobs.numElements())
               {
                  mGroupIobs.resizeAndPreserve(mNbGroup+500);
                  mGroupSigma.resizeAndPreserve(mNbGroup+500);
                  mGroupWeight.resizeAndPreserve(mNbGroup+500);
               }
            }
         }
         mH.resizeAndPreserve(mNbRefl);
         mK.resizeAndPreserve(mNbRefl);
         mL.resizeAndPreserve(mNbRefl);
         mObsIntensity.resizeAndPreserve(mNbRefl);
         mObsSigma.resizeAndPreserve(mNbRefl);
         mWeight.resizeAndPreserve(mNbRefl);
         mGroupIndex.resizeAndPreserve(mNbRefl);
         
         mGroupIobs.resizeAndPreserve(mNbGroup);
         mGroupWeight.resizeAndPreserve(mNbGroup);
         mGroupSigma.resizeAndPreserve(mNbGroup);
      
         mHasObservedData=true;
         
         mMultiplicity.resize(mNbRefl);
         mMultiplicity=1;
         
         this->PrepareHKLarrays();
         this->SortReflectionBySinThetaOverLambda();
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPatternBackground
//
////////////////////////////////////////////////////////////////////////
void PowderPatternBackground::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("PowderPatternBackground::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("PowderPatternBackground");
   tag.AddAttribute("Name",this->GetName());
   os <<tag<<endl;
   indent++;
   
   mInterpolationModel.XMLOutput(os,indent);
   os<<endl;

   XMLCrystTag tag2("XIntensityList");
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag2<<endl;
   
   REAL scale=1.0;
   if(this->GetParentPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
      scale=RAD2DEG;

   for(long j=0;j<mBackgroundNbPoint;j++)
   {
   
      for(int i=0;i<=indent;i++) os << "  " ;
      os << mBackgroundInterpPointX(j)*scale <<" "
         << mBackgroundInterpPointIntensity(j) <<" "
         << !this->GetPar(mBackgroundInterpPointIntensity.data()+j).IsFixed()<<" "
         <<endl;
   }
   
   tag2.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag2<<endl;
   
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   this->GetPar("ML Model Error").XMLOutput(os,"ML Model Error",indent);
   os <<endl;
   #endif

   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("PowderPatternBackground::XMLOutput():"<<this->GetName(),5)
}

void PowderPatternBackground::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("PowderPatternBackground::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Interpolation"==tagg.GetAttributeName(i))
      {// Obsolete, but we must still read this
         if("Linear"==tagg.GetAttributeValue(i)) mInterpolationModel.SetChoice(0);
         if("Spline"==tagg.GetAttributeValue(i)) mInterpolationModel.SetChoice(1);
      }
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("PowderPatternBackground"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("PowderPatternBackground::Exit():"<<this->GetName(),5)
         return;
      }
      if(("TThetaIntensityList"==tag.GetName())||("XIntensityList"==tag.GetName()))
      {
         long nbPoint=0;
         CrystVector_REAL bckgd2Theta(100);
         CrystVector_REAL bckgd(100);
         CrystVector_bool fix(100);
         do
         {
            VFN_DEBUG_MESSAGE("PowderPatternBackground::XMLInput():"<<mBackgroundNbPoint,1)
            is >>bckgd2Theta(nbPoint)
               >>bckgd(nbPoint)
               >>fix(nbPoint);
            nbPoint++;
            if(nbPoint==bckgd2Theta.numElements())
            {
               bckgd2Theta.resizeAndPreserve(nbPoint+100);
               bckgd.resizeAndPreserve(nbPoint+100);
               fix.resizeAndPreserve(nbPoint+100);
            }
            while(0==isgraph(is.peek())) is.get();//Why do I need that ?
            //cout << is.peek()<<" "<<nbrefl<<endl;
         }
         while(is.peek()!='<');//until next tag
         bckgd2Theta.resizeAndPreserve(nbPoint);
         bckgd.resizeAndPreserve(nbPoint);
         if(this->GetParentPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) 
            bckgd2Theta*= DEG2RAD;
         this->SetInterpPoints(bckgd2Theta,bckgd);
         this->InitRefParList();
         //read closing tag
         XMLCrystTag junkEndTag(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("ML Model Error"==tag.GetAttributeValue(i))
               {
                  #ifdef USE_BACKGROUND_MAXLIKE_ERROR
                  this->GetPar("ML Model Error").XMLInput(is,tag);
                  break;
                  #endif
               }
            }
         }
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////
void PowderPatternDiffraction::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("PowderPatternCrystal");
   tag.AddAttribute("Name",this->GetName());
   tag.AddAttribute("Crystal",this->GetCrystal().GetName());
   {
      stringstream ss;
      ss<<this->IsIgnoringImagScattFact();
      tag.AddAttribute("IgnoreImagScattFact",ss.str());
   }
   os <<tag<<endl;
   indent++;

   if(mpReflectionProfile!=0) mpReflectionProfile->XMLOutput(os,indent);

   this->GetPar(&mGlobalBiso).XMLOutput(os,"globalBiso",indent);
   os <<endl;
   
   if(mCorrTextureMarchDollase.GetNbPhase()>0)
   {
      mCorrTextureMarchDollase.XMLOutput(os,indent);
   }
   
   if(mFhklObsSq.numElements()>0)
   {
      XMLCrystTag tag2("FhklObsSq");
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag2<<endl;
      
      for(long j=0;j<this->GetNbRefl();j++)
      {
         for(int i=0;i<=indent;i++) os << "  " ;
         os << mIntH(j) <<" "
            << mIntK(j) <<" "
            << mIntL(j) <<" "
            << mFhklObsSq(j) <<endl;
      }
      
      tag2.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag2<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("PowderPatternDiffraction::XMLOutput():"<<this->GetName(),5)
}

void PowderPatternDiffraction::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Crystal"==tagg.GetAttributeName(i)) 
         this->SetCrystal(gCrystalRegistry.GetObj(tagg.GetAttributeValue(i)));
      if("NeedLorentzCorr"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         bool b;
         ss>>b;
         //mNeedLorentzCorr=b;
         //mClockLorentzPolarSlitCorrPar.Reset();
      }
      if("NeedPolarCorr"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         bool b;
         ss>>b;
         //mNeedPolarCorr=b;
         //mClockLorentzPolarSlitCorrPar.Reset();
      }
      if("Polar_AFactor"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         float b;
         ss>>b;
         //mPolarAfactor=b;
         //mClockLorentzPolarSlitCorrPar.Reset();
      }
      if("NeedSlitApertureCorr"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         bool b;
         ss>>b;
         //mNeedSlitApertureCorr=b;
         //mClockLorentzPolarSlitCorrPar.Reset();
      }
      if("IgnoreImagScattFact"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         bool b;
         ss>>b;
         this->SetIsIgnoringImagScattFact(b);
         mClockLorentzPolarSlitCorrPar.Reset();
      }
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("PowderPatternCrystal"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("PowderPatternDiffraction::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("globalBiso"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mGlobalBiso).XMLInput(is,tag);
                  break;
               }
               if("U"==tag.GetAttributeValue(i))
               {
                  mpReflectionProfile->GetPar("U").XMLInput(is,tag);
                  break;
               }
               if("V"==tag.GetAttributeValue(i))
               {
                  mpReflectionProfile->GetPar("V").XMLInput(is,tag);
                  break;
               }
               if("W"==tag.GetAttributeValue(i))
               {
                  mpReflectionProfile->GetPar("W").XMLInput(is,tag);
                  break;
               }
               if("Eta0"==tag.GetAttributeValue(i))
               {
                  mpReflectionProfile->GetPar("Eta0").XMLInput(is,tag);
                  break;
               }
               if("Eta1"==tag.GetAttributeValue(i))
               {
                  mpReflectionProfile->GetPar("Eta1").XMLInput(is,tag);
                  break;
               }
               if("W0"==tag.GetAttributeValue(i))
               {
                  //:TODO: mpReflectionProfile->GetPar("Eta0").XMLInput(is,tag);
                  break;
               }
               if("W1"==tag.GetAttributeValue(i))
               {
                  //:TODO: mpReflectionProfile->GetPar("Eta1").XMLInput(is,tag);
                  break;
               }
               if("W2"==tag.GetAttributeValue(i))
               {
                  //:TODO: mpReflectionProfile->GetPar("Eta2").XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i)) 
            {
               if("Profile Type"!=tag.GetAttributeValue(i))
                  mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
            }
         }
         continue;
      }
      if("TextureMarchDollase"==tag.GetName())
      {
         mCorrTextureMarchDollase.XMLInput(is,tag);
         continue;
      }
      if("ReflectionProfilePseudoVoigt"==tag.GetName())
      {
         if(mpReflectionProfile==0)
         {
            mpReflectionProfile=new ReflectionProfilePseudoVoigt;
         }
         else
            if(mpReflectionProfile->GetClassName()!="ReflectionProfilePseudoVoigt")
            {
               delete mpReflectionProfile;
               mpReflectionProfile=new ReflectionProfilePseudoVoigt;
            }
         mpReflectionProfile->XMLInput(is,tag);
         continue;
      }
      if("ReflectionProfileDoubleExponentialPseudoVoigt"==tag.GetName())
      {
         if(mpReflectionProfile==0)
         {
            mpReflectionProfile
               =new ReflectionProfileDoubleExponentialPseudoVoigt(this->GetCrystal());
         }
         else
            if(mpReflectionProfile->GetClassName()!="ReflectionProfileDoubleExponentialPseudoVoigt")
            {
               delete mpReflectionProfile;
               mpReflectionProfile
                  =new ReflectionProfileDoubleExponentialPseudoVoigt(this->GetCrystal());
            }
         mpReflectionProfile->XMLInput(is,tag);
         continue;
      }
      if("FhklObsSq"==tag.GetName())
      {
         // We ignore the h,k,l arrays - just assume for now the order
         // is unchanged compared to when the extraction was made...
         long nbrefl=0;
         long junk;
         mFhklObsSq.resize(100);
         do
         {
            is >>junk>>junk>>junk
               >>mFhklObsSq(nbrefl);
            nbrefl++;
            if(nbrefl==mFhklObsSq.numElements()) mFhklObsSq.resizeAndPreserve(nbrefl+100);
            while(0==isgraph(is.peek())) is.get();
         }
         while(is.peek()!='<');//until next tag
         XMLCrystTag junkEndTag(is);
         
         mFhklObsSq.resizeAndPreserve(nbrefl);
         mClockGetFhklObsSq.Click();
      }
   }
}
////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPattern
//
////////////////////////////////////////////////////////////////////////
void PowderPattern::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("PowderPattern::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("PowderPattern");
   tag.AddAttribute("Name",mName);
   os <<tag<<endl;
   indent++;
   
   this->GetPar(&mXZero).XMLOutput(os,"Zero",indent);
   os <<endl;
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      this->GetPar(&mDIFC).XMLOutput(os,"TOF-DIFC",indent);
      os <<endl;

      this->GetPar(&mDIFA).XMLOutput(os,"TOF-DIFA",indent);
      os <<endl;
   }
   else
   {
      this->GetPar(&m2ThetaDisplacement).XMLOutput(os,"2ThetaDisplacement",indent);
      os <<endl;

      this->GetPar(&m2ThetaTransparency).XMLOutput(os,"2ThetaTransparency",indent);
      os <<endl;
   }
   
   for(unsigned int i=0;i<this->GetNbOption();i++)
   {
      this->GetOption(i).XMLOutput(os,indent);
      os <<endl<<endl;
   }
   
   mRadiation.XMLOutput(os,indent);
   os <<endl;
   {
      for(int i=0;i<indent;i++) os << "  " ;
      XMLCrystTag tag2("MaxSinThetaOvLambda");
      os << tag2<< mMaxSinThetaOvLambda;
      tag2.SetIsEndTag(true);
      os << tag2<<endl<<endl;
   }
   
   for(int j=0;j<mPowderPatternComponentRegistry.GetNb();j++)
   {
      mPowderPatternComponentRegistry.GetObj(j).XMLOutput(os,indent);
      XMLCrystTag tagg("PowderPatternComponent",false,true);
      {
         stringstream ss;
         ss<<mScaleFactor(j);
         tagg.AddAttribute("Scale",ss.str());
      }
      tagg.AddAttribute("Name",mPowderPatternComponentRegistry.GetObj(j).GetName());
      os<<endl;
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tagg<<endl<<endl;
   }
   XMLCrystTag tag2("XIobsSigmaWeightList");
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;

      REAL scale=1.0;
      if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) 
         scale=RAD2DEG;

      for(unsigned long j=0;j<this->GetNbPoint();j++)
      {
         for(int i=0;i<=indent;i++) os << "  " ;
         os << scale*mX(j) <<" "
            << mPowderPatternObs(j) <<" "
            << mPowderPatternObsSigma(j) <<" "
            << mPowderPatternWeight(j) <<" "
            <<endl;
      }
      tag2.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;
   
   for(int j=0;j<mExcludedRegionMinX.numElements();j++)
   {
      XMLCrystTag tag3("ExcludeX");
      for(int i=0;i<indent;i++) os << "  " ;
      if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
      {
         os << tag3 
            << mExcludedRegionMinX(j) <<" "
            << mExcludedRegionMaxX(j) ;
      }
      else
      {
         os << tag3 
            << mExcludedRegionMinX(j)*RAD2DEG <<" "
            << mExcludedRegionMaxX(j)*RAD2DEG ;
      }
      tag3.SetIsEndTag(true);
      os<<tag3<<endl;
   }
   
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("PowderPattern::XMLOutput():"<<this->GetName(),5)
}

void PowderPattern::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("PowderPattern::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("PowderPattern"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("PowderPattern::Exit():"<<this->GetName(),5)
         return;
      }
      if("Radiation"==tag.GetName()) mRadiation.XMLInput(is,tag);
      if("MaxSinThetaOvLambda"==tag.GetName())
      {
         is>>mMaxSinThetaOvLambda;
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if(("2ThetaZero"==tag.GetAttributeValue(i)) ||("Zero"==tag.GetAttributeValue(i)))
               {
                  this->GetPar(&mXZero).XMLInput(is,tag);
                  break;
               }
               if("2ThetaDisplacement"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&m2ThetaDisplacement).XMLInput(is,tag);
                  break;
               }
               if("2ThetaTransparency"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&m2ThetaTransparency).XMLInput(is,tag);
                  break;
               }
               if("TOF-DIFC"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mDIFC).XMLInput(is,tag);
                  break;
               }
               if("TOF-DIFA"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mDIFA).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
      if("PowderPatternBackground"==tag.GetName())
      {
         PowderPatternBackground *comp=new PowderPatternBackground;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternCrystal"==tag.GetName())
      {
         PowderPatternDiffraction *comp=new PowderPatternDiffraction;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternComponent"==tag.GetName())
      {
         REAL scale=1.0;
         string name;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Scale"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>scale;
               continue;
            }
            if("Name"==tag.GetAttributeName(i)) name=tag.GetAttributeValue(i);
         }
         this->AddPowderPatternComponent(gPowderPatternComponentRegistry.GetObj(name));
         mScaleFactor(mPowderPatternComponentRegistry.GetNb()-1)=scale;
         VFN_DEBUG_MESSAGE("->Adding Component :"<<name<<"with scale="<<scale,8);
         continue;
      }
      if("ExcludeX"==tag.GetName())
      {
         float min,max;
         is>>min>>max;
         if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
            this->AddExcludedRegion(min,max);
         else this->AddExcludedRegion(min*DEG2RAD,max*DEG2RAD);
         XMLCrystTag end(is);
         continue;
      }
      if("IobsSigmaWeightList"==tag.GetName())
      {
         // Old version, just for 2theta-intensity pattern (no TOF)
         VFN_DEBUG_ENTRY("Loading Iobs-Sigma-Weight List...",8);
         REAL min,step;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("TThetaMin"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>min;
               VFN_DEBUG_MESSAGE("2Theta min="<<min,8);
               min*=DEG2RAD;
            }
            if("TThetaStep"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>step;
               VFN_DEBUG_MESSAGE("2Theta step="<<step<<tag.GetAttributeValue(i),8);
               step*=DEG2RAD;
            }
         }
         while(0==isgraph(is.peek())) is.get();
         if(is.peek()=='<')
         {
            cout <<"PowderPattern::XMLInput(): no data point in the powder pattern !"<<endl;
            XMLCrystTag junk(is);
            VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
            continue;
         }
         mNbPoint=0;
         mPowderPatternObs.resize(500);
         mPowderPatternObsSigma.resize(500);
         mPowderPatternWeight.resize(500);
         do
         {
            is >>mPowderPatternObs(mNbPoint)
               >>mPowderPatternObsSigma(mNbPoint)
               >>mPowderPatternWeight(mNbPoint);
            mNbPoint++;
            VFN_DEBUG_MESSAGE("Point #"<<mNbPoint,5);
            if(mNbPoint==(unsigned long)mPowderPatternObs.numElements())
            {
               mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
               mPowderPatternWeight.resizeAndPreserve(mNbPoint+500);
            }
            while(0==isgraph(is.peek())) is.get();
         }
         while(is.peek()!='<');//until next tag
         this->SetPowderPatternPar(min,step,mNbPoint);
         mClockPowderPatternPar.Click();
         
         XMLCrystTag junk(is);
         VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
         continue;
      }
      if("XIobsSigmaWeightList"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Loading X-Iobs-Sigma-Weight List...",8);
         while(0==isgraph(is.peek())) is.get();
         if(is.peek()=='<')
         {
            cout <<"PowderPattern::XMLInput(): no data point in the powder pattern !"<<endl;
            XMLCrystTag junk(is);
            VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
            continue;
         }
         mNbPoint=0;
         mX.resize(500);
         mPowderPatternObs.resize(500);
         mPowderPatternObsSigma.resize(500);
         mPowderPatternWeight.resize(500);
         do
         {
            is >>mX(mNbPoint)
               >>mPowderPatternObs(mNbPoint)
               >>mPowderPatternObsSigma(mNbPoint)
               >>mPowderPatternWeight(mNbPoint);
            mNbPoint++;
            VFN_DEBUG_MESSAGE("Point #"<<mNbPoint,5);
            if(mNbPoint==(unsigned long)mPowderPatternObs.numElements())
            {
               mX.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
               mPowderPatternWeight.resizeAndPreserve(mNbPoint+500);
            }
            while(0==isgraph(is.peek())) is.get();
         }
         while(is.peek()!='<');//until next tag
         mX.resizeAndPreserve(mNbPoint);
         if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) 
            mX*=DEG2RAD;
         this->SetPowderPatternX(mX);
         
         XMLCrystTag junk(is);
         VFN_DEBUG_EXIT("Loading X-Iobs-Sigma-Weight List...",8);
         continue;
      }
   }
}
} //namespace
