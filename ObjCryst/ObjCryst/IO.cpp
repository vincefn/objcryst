/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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

#include <iostream>
#include <fstream>
#include <sstream>

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
   if(!out) {};//:TODO:
   
   XMLCrystTag tag("ObjCryst");
   time_t date=time(0);
   char strDate[40];
   strftime(strDate,sizeof(strDate),"%Y-%m-%dT%H:%M:%S%Z",gmtime(&date));//%Y-%m-%dT%H:%M:%S%Z
   tag.AddAttribute("Date",strDate);
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
   
   out.close();
   
   VFN_DEBUG_EXIT("XMLCrystFileSaveGlobal(filename):End",5)
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
   VFN_DEBUG_ENTRY("XMLCrystFileLoadAllObject(filename,IOCrystTag,T&)",5)
   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   
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
   VFN_DEBUG_EXIT("XMLCrystFileLoadAllObject(filename,IOCrystTag,T&)",5)
}
#if 0
////////////////////////////////////////////////////////////////////////
//
//    Global functions -OLD
//
////////////////////////////////////////////////////////////////////////

void IOCrystFileSaveGlobal(const string & filename)
{
   VFN_DEBUG_MESSAGE("IOCrystFileSaveGlobal(filename)",5)
   
   ofstream out(filename.c_str());
   if(!out) {};//:TODO:
   
   for(int i=0;i<gCrystalRegistry.GetNb();i++)
      gCrystalRegistry.GetObj(i).XMLOutput(out);
   
   for(int i=0;i<gDiffractionDataSingleCrystalRegistry.GetNb();i++)
      gDiffractionDataSingleCrystalRegistry.GetObj(i).XMLOutput(out);
   
   for(int i=0;i<gPowderPatternRegistry.GetNb();i++)
      gPowderPatternRegistry.GetObj(i).XMLOutput(out);
   
   for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
      gOptimizationObjRegistry.GetObj(i).XMLOutput(out);
   
   out.close();
   
   VFN_DEBUG_MESSAGE("IOCrystFileSaveGlobal(filename):End",5)
}

ObjRegistry<IOCrystTag> IOCrystFileLoadObjectList(const string & filename)
{
   VFN_DEBUG_MESSAGE("IOCrystFileLoadObjectList(filename)",5)

   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   ObjRegistry<IOCrystTag> reg;
   for(;;)
   {
      IOCrystTag *pTag =new IOCrystTag (is);
      if(true==is.eof())
      {
         VFN_DEBUG_MESSAGE("IOCrystFileLoadObjectList(filename):End",5)
         for(int i=0;i<reg.GetNb();i++) reg.GetObj(i).Print();
         is.close();
         return reg;
      }
      //pTag->Print();
      if(  (pTag->GetType()!="Param")
         &&(pTag->GetType()!="ParList")
         &&(pTag->GetType()!="Radiation")
         &&(pTag->GetType()!="Atom")
         &&(pTag->GetType()!="ScatteringPowerAtom")
         &&(pTag->GetType()!="ZScatterer")
         &&(pTag->IsClosingTag()==false)
         ) reg.Register(*pTag);
      else delete pTag;
   }
   
}

template<class T> void IOCrystFileLoadObject(const string & filename,
                                             const IOCrystTag &tag,
                                             T* obj)
{
   VFN_DEBUG_ENTRY("IOCrystFileLoadObject(filename,IOCrystTag,T&)",5)

   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   IOCrystTag ttag(is);;
   for(;;)
   {
      if(true==is.eof())
      {
         cout<<"IOCrystFileLoadObject(filename,IOCrystTag,T&):Not Found !"<<endl;
         is.close();
         return;
      }
      if(tag==ttag) break;
      ttag.XMLInput(is);
   }
   VFN_DEBUG_MESSAGE("IOCrystFileLoadObject(filename,IOCrystTag,T&):Found"<<ttag.GetName(),5)
   obj = new T;
   obj->XMLInputOld(is,ttag);
   is.close();
   VFN_DEBUG_EXIT("IOCrystFileLoadObject(filename,IOCrystTag,T&)",5)
}

template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    Crystal*);
template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    PowderPattern*);
template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    DiffractionDataSingleCrystal*);
//template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
//                                    ZScatterer*);
template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    PowderPatternBackground*);
template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    PowderPatternDiffraction*);
template void IOCrystFileLoadObject(const string &,const IOCrystTag &,
                                    OptimizationObj*);
                                    
void IOCrystFileLoadAllObject(const string & filename)
{
   ifstream is(filename.c_str());
   if(!is){};//:TODO:
   IOCrystTag tag(is);
   for(;;)
   {
      if(true==is.eof()) break;
      if(tag.GetType()=="Crystal")
      {
         Crystal* obj = new Crystal;
         obj->XMLInputOld(is,tag);
      }
      if(tag.GetType()=="PowderPattern")
      {
         PowderPattern* obj = new PowderPattern;
         obj->XMLInputOld(is,tag);
      }
      if(tag.GetType()=="DiffractionDataSingleCrystal")
      {
         DiffractionDataSingleCrystal* obj = new DiffractionDataSingleCrystal;
         obj->XMLInputOld(is,tag);
      }
      if(tag.GetType()=="GlobalOptimObj")
      {
         GlobalOptimObj* obj = new GlobalOptimObj;
         obj->XMLInputOld(is,tag);
      }
      tag.XMLInput(is);
   }
}
#endif
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
   
   for(int i=0;i<=indent;i++) os << "  " ;
   if(true==this->mIsIsotropic)
      this->GetPar(&mBiso).XMLOutput(os,"Biso",0);
   os<<endl;
   
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
               if("Biso"==tag.GetAttributeValue(i)) this->GetPar(&mBiso).XMLInput(is,tag);
               break;
            }
         }
         continue;
      }
   }
}
#if 0
void ScatteringPowerAtom::XMLInputOld(istream &is,const IOCrystTag &tag)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::XMLInput():"<<this->GetName(),5)
   switch(tag.GetVersion())
   {
      case 0:
      {
         string symbol;
         bool isIso;
         REAL bIso;
         //First read name of the object
         IOCrystExtractNameQuoted(is,symbol);
         is>> isIso;
         if(true==isIso)
         {
            is>>bIso;
            VFN_DEBUG_MESSAGE("->"<<tag.GetName()<<" "<<symbol<<" "<<isIso,5)
            this->Init(tag.GetName(),symbol,bIso);
            bool fix;
            is>>fix;
            this->GetPar(&mBiso).SetIsFixed(!fix);
         }
         this->UpdateDisplay();
         //else :TODO:
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::XMLInput():End",5)
}
#endif
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
#if 0
void Atom::XMLInputOld(istream &is,const IOCrystTag &tag)
{
   VFN_DEBUG_MESSAGE("Atom::XMLInput():"<<this->GetName(),5)
   switch(tag.GetVersion())
   {
      case 0:
      {
         REAL x,y,z,p;
         string scattPow;
         is >> x >> y>> z>> p;
         IOCrystExtractNameQuoted(is,scattPow);
         //cout <<name<<" "<< x<<" "<<y<<" "<<z<<" "<<p<<" "<<scattPow<<endl;
         //gScatteringPowerAtomRegistry.GetObj(scattPow).Print();
         this->Init(x,y,z,tag.GetName(),&gScatteringPowerAtomRegistry.GetObj(scattPow),p);
         bool fix;
         is>>fix;
         this->GetPar(mXYZ.data()+0).SetIsFixed(!fix);
         is>>fix;
         this->GetPar(mXYZ.data()+1).SetIsFixed(!fix);
         is>>fix;
         this->GetPar(mXYZ.data()+2).SetIsFixed(!fix);
         is>>fix;
         this->GetPar(&mOccupancy).SetIsFixed(!fix);
         this->UpdateDisplay();
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
   VFN_DEBUG_MESSAGE("Atom::XMLInput():End",5)
}
#endif
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
#if 0
void ZScatterer::XMLInputOld(istream &is,const IOCrystTag &tag)
{
   VFN_DEBUG_ENTRY("ZScatterer::XMLInput():"<<this->GetName(),5)
   switch(tag.GetVersion())
   {
      case 0:
      {
         this->SetName(tag.GetName());
         do
         {
            VFN_DEBUG_MESSAGE("ZScatterer::XMLInput():Reading new tag..",5)
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.IsClosingTag()==true)
            {
               VFN_DEBUG_EXIT("ZScatterer::XMLInput():"<<this->GetName()<<"->UpdateDisplay",5)
               this->UpdateDisplay();
               VFN_DEBUG_EXIT("ZScatterer::XMLInput():"<<this->GetName(),5)
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="XYZPhiChiPsiPopu")
               {
                  VFN_DEBUG_MESSAGE("ZScatterer::XMLInput():XMLInput position&popu",5)
                  REAL x,y,z,phi,chi,psi,p;
                  is >>x>>y>>z>>phi>>chi>>psi>>p;
                  this->SetX(x);
                  this->SetY(y);
                  this->SetZ(z);
                  this->SetPhi(phi*DEG2RAD);
                  this->SetChi(chi*DEG2RAD);
                  this->SetPsi(psi*DEG2RAD);
                  this->SetOccupancy(p);
                  bool fix;
                  is>>fix;
                  this->GetPar(mXYZ.data()+0).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(mXYZ.data()+1).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(mXYZ.data()+2).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mPhi).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mChi).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mPsi).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mOccupancy).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="Atom")
               {
                  VFN_DEBUG_MESSAGE("ZScatterer::XMLInput():XMLInput an atom",5)
                  string scattPow;
                  string name;
                  int abond,aangle,adihed;
                  REAL bond,angle,dihed,p;
                  IOCrystExtractNameQuoted(is,name);
                  is >>abond>>bond>>aangle>>angle>>adihed>>dihed>>p;
                  IOCrystExtractNameQuoted(is,scattPow);
                  //cout <<abond<<" "<<bond<<" "<<aangle<<" "<<angle
                  //     <<" "<<adihed<<" "<<dihed<<" "<<p<<" "<<scattPow<<endl;
                  this->AddAtom(name,&gScatteringPowerAtomRegistry.GetObj(scattPow),
                                abond,bond,aangle,angle*DEG2RAD,adihed,dihed*DEG2RAD,p);
                  const int j=mNbAtom-1;
                  bool fix;
                  is>>fix;
                  this->GetPar(&(mZAtomRegistry.GetObj(j).mBondLength)).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&(mZAtomRegistry.GetObj(j).mAngle)).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&(mZAtomRegistry.GetObj(j).mDihed)).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&(mZAtomRegistry.GetObj(j).mOccupancy)).SetIsFixed(!fix);
                  VFN_DEBUG_MESSAGE("ZScatterer::XMLInput():XMLInput an atom:Finished",5)
                  continue;
               }
            }
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
   VFN_DEBUG_EXIT("ZScatterer::XMLInput():"<<this->GetName(),5)
}
#endif
////////////////////////////////////////////////////////////////////////
//
//    I/O Crystal
//
////////////////////////////////////////////////////////////////////////
void Crystal::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Crystal::XMLOutput():"<<this->GetName(),5)
   using namespace std;
   // set locale settings: not available yet in GNU library ?
   //os.imbue(locale::classic());
   
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Crystal");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("SpaceGroup",mSpaceGroup.GetName());
   os <<tag<<endl;
   indent++;
   
   this->GetPar(mCellDim.data()+0).XMLOutput(os,"a",indent);
   os <<endl;
   
   this->GetPar(mCellDim.data()+1).XMLOutput(os,"b",indent);
   os <<endl;
   
   this->GetPar(mCellDim.data()+2).XMLOutput(os,"c",indent);
   os <<endl;
   
   this->GetPar(mCellDim.data()+3).XMLOutput(os,"alpha",indent);
   os <<endl;
   
   this->GetPar(mCellDim.data()+4).XMLOutput(os,"beta",indent);
   os <<endl;
   
   this->GetPar(mCellDim.data()+5).XMLOutput(os,"gamma",indent);
   os <<endl;
   
   mUseDynPopCorr.XMLOutput(os,indent);
   os <<endl<<endl;

   mConstrainLatticeToSpaceGroup.XMLOutput(os,indent);
   os <<endl<<endl;
   
   for(int i=0;i<mScatteringPowerRegistry.GetNb();i++) 
      mScatteringPowerRegistry.GetObj(i).XMLOutput(os,indent);
   os <<endl;
   for(int i=0;i<mScattererRegistry.GetNb();i++) 
      mScattererRegistry.GetObj(i).XMLOutput(os,indent);
   os <<endl;

   if(mBumpDistanceMatrix.numElements()>0)
   {
      const int num=mBumpDistanceMatrix.rows();
      for(int i=0;i<num;i++)
         for(int j=i;j<num;j++)
            if(.1<mBumpDistanceMatrix(i,j))
            {
               for(int k=0;k<=indent;k++) os << "  " ;
               XMLCrystTag tagBump("AntiBumpDistance");
               const ScatteringPower *scatt1=0;
               for(int k=0;k<mScatteringPowerRegistry.GetNb();k++)
                  if(i==mScatteringPowerRegistry.GetObj(k).GetDynPopCorrIndex())
                  {
                     scatt1=&(mScatteringPowerRegistry.GetObj(k));
                     break;
                  }
               const ScatteringPower *scatt2=0;
               for(int k=0;k<mScatteringPowerRegistry.GetNb();k++)
                  if(j==mScatteringPowerRegistry.GetObj(k).GetDynPopCorrIndex())
                  {
                     scatt2=&(mScatteringPowerRegistry.GetObj(k));
                     break;
                  }
               tagBump.AddAttribute("ScattPow1",scatt1->GetName());
               tagBump.AddAttribute("ScattPow2",scatt2->GetName());
               os<<tagBump;
               tagBump.SetIsEndTag(true);
               os<<sqrt(mBumpDistanceMatrix(i,j))<<tagBump<<endl;
            }
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
               if("a"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+0).XMLInput(is,tag);
                  break;
               }
               if("b"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+1).XMLInput(is,tag);
                  break;
               }
               if("c"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+2).XMLInput(is,tag);
                  break;
               }
               if("alpha"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+3).XMLInput(is,tag);
                  break;
               }
               if("beta"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+4).XMLInput(is,tag);
                  break;
               }
               if("gamma"==tag.GetAttributeValue(i))
               {
                  this->GetPar(mCellDim.data()+5).XMLInput(is,tag);
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
   }
}
#if 0
void Crystal::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("Crystal::XMLInput():"<<this->GetName(),5)
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
   
   
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.GetType()=="Param")
            {
               VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a Param",5)
               if(tag.GetName()=="Lattice")
               {
                  VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading Lattice",5)
                  REAL a,b,c,alpha,beta,gamma;
                  string spg;
                  is >>a>>b>>c>>alpha>>beta>>gamma;
                  IOCrystExtractNameQuoted(is,spg);
                  this->Init(a,b,c,alpha*DEG2RAD,beta*DEG2RAD,gamma*DEG2RAD,
                             spg,tagg.GetName());
                  bool fix;
                  is >>fix;
                  this->GetPar(mCellDim.data()+0).SetIsFixed(!fix);
                  is >>fix;
                  this->GetPar(mCellDim.data()+1).SetIsFixed(!fix);
                  is >>fix;
                  this->GetPar(mCellDim.data()+2).SetIsFixed(!fix);
                  is >>fix;
                  this->GetPar(mCellDim.data()+3).SetIsFixed(!fix);
                  is >>fix;
                  this->GetPar(mCellDim.data()+4).SetIsFixed(!fix);
                  is >>fix;
                  this->GetPar(mCellDim.data()+5).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="UseDynamicOccupancyCorr")
               {
                  VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading UseDynamicOccupancyCorr",5)
                  bool useDynPopCorr;
                  is >>useDynPopCorr;
                  this->SetUseDynPopCorr(useDynPopCorr);
                  continue;
               }
               if(tag.GetName()=="BumpMergeDistance")
               {
                  VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a Bump/Merge distance",5)
                  float dist;
                  bool allowMerge;
                  string scattPow1,scattPow2;
                  IOCrystExtractNameQuoted(is,scattPow1);
                  IOCrystExtractNameQuoted(is,scattPow2);
                  is >> dist >>allowMerge;
                  this->SetBumpMergeDistance(mScatteringPowerRegistry.GetObj(scattPow1),
                                             mScatteringPowerRegistry.GetObj(scattPow2),
                                             dist,allowMerge);
                  continue;
               }
               cout <<"Cannot recognize tag name:";
               IOCrystXMLOutputNameQuoted(cout,tag.GetName());
            }
            if(tag.IsClosingTag()==true)
            {
               VFN_DEBUG_MESSAGE("Crystal::XMLInput():End",5)
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Atom")
            {
               VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading an Atom",5)
               Atom *at=new Atom;
               at->XMLInputOld(is,tag);
               this->AddScatterer(at);
               continue;
            }
            if(tag.GetType()=="ScatteringPowerAtom")
            {
               VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a ScatteringPowerAtom",5)
               ScatteringPowerAtom *sc=new ScatteringPowerAtom;
               sc->XMLInputOld(is,tag);
               this->AddScatteringPower(sc);
               continue;
            }
            if(tag.GetType()=="ZScatterer")
            {
               VFN_DEBUG_MESSAGE("Crystal::XMLInput():reading a ZScatterer",5)
               ZScatterer *z=new ZScatterer("",*this);
               z->XMLInputOld(is,tag);
               this->AddScatterer(z);
               continue;
            }
            cout <<"Cannot recognize tag type:";
            IOCrystXMLOutputNameQuoted(cout,tag.GetType());
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
}
#endif
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
   XMLCrystTag tag2("LinearPolarRate");
   os << tag2<< mLinearPolarRate;
   tag2.SetIsEndTag(true);
   os << tag2<<endl;
   
   switch(this->GetWavelengthType())
   {
      case WAVELENGTH_MONOCHROMATIC: this->GetPar(mWavelength.data()).XMLOutput(os,indent);break;
      case WAVELENGTH_ALPHA12:
      {
         this->GetPar(mWavelength.data()).XMLOutput(os,"Wavelength",indent);
         break;
      }
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
      if("XRayTube"==tagg.GetAttributeName(i)) this->SetWavelength(tagg.GetAttributeValue(i));
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
#if 0
void Radiation::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("Radiation::XMLInput():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         do
         {
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.IsClosingTag()==true)
            {
               VFN_DEBUG_MESSAGE("Radiation::XMLInput():End",5)
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="RadiationType")
               {
                  string rad;
                  IOCrystExtractNameSpace(is,rad);
                  if(rad=="neutron")
                  {
                     this->SetRadiationType(RAD_NEUTRON);
                     continue;
                  }
                  if(rad=="xray")
                  {
                     this->SetRadiationType(RAD_XRAY);
                     continue;
                  }
                  if(rad=="RAD_ELECTRON")
                  {
                     this->SetRadiationType(RAD_NEUTRON);
                     continue;
                  }
               }
               if(tag.GetName()=="Wavelength")
               {
                  string rad;
                  IOCrystExtractNameSpace(is,rad);
                  if(rad=="monochromatic")
                  {
                     REAL l;
                     is>>l;
                     this->SetWavelength(l);
                     bool fix;
                     is>>fix;
                     this->GetPar(mWavelength.data()+0).SetIsFixed(!fix);
                     continue;
                  }
                  if(rad=="tube")
                  {
                     string tube;
                     IOCrystExtractNameSpace(is,tube);
                     REAL ratio;
                     is>>ratio;
                     this->SetWavelength(tube,ratio);
                     continue;
                  }
               }
            }
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
}
#endif
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
   
   mRadiation.XMLOutput(os,indent);
   os <<endl;

   this->GetPar(&mGlobalBiso).XMLOutput(os,"globalBiso",indent);
   os <<endl;
   
   mTwinningOption.XMLOutput(os,indent);
   os <<endl;
   
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag2("MaxSinThetaOvLambda");
   os << tag2<< mMaxSinThetaOvLambda;
   tag2.SetIsEndTag(true);
   os << tag2<<endl<<endl;
   
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
         VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::Exit():"<<this->GetName(),5)
         return;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
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
         this->SortReflectionByTheta();
         this->CalcIcalc();
         this->FitScaleFactorForRw();
      }
   }
}
#if 0
void DiffractionDataSingleCrystal::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::XMLInput():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.IsClosingTag()==true)
            {
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="Crystal")
               {
                  string cryst;
                  IOCrystExtractNameQuoted(is,cryst);
                  //:TODO: Check we found the Crystal !
                  this->SetCrystal(gCrystalRegistry.GetObj(cryst));
                  continue;
               }
               if(tag.GetName()=="IgnoreImagScattFact")
               {
                  bool b;
                  is>>b;
                  this->SetIsIgnoringImagScattFact(b);
                  continue;
               }
               if(tag.GetName()=="UseFastLessPreciseFunc")
               {
                  bool b;
                  is>>b;
                  this->SetUseFastLessPreciseFunc(b);
                  continue;
               }
            }
            if(tag.GetType()=="Radiation")
            {
               mRadiation.XMLInputOld(is,tag);
               continue;
            }
            if(tag.GetType()=="ParList")
            {
               if(tag.GetName()=="HKLIobsSigmaWeight")
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
                  
                  h.resizeAndPreserve(nbrefl);
                  k.resizeAndPreserve(nbrefl);
                  l.resizeAndPreserve(nbrefl);
                  iobs.resizeAndPreserve(nbrefl);
                  sigma.resizeAndPreserve(nbrefl);
                  weight.resizeAndPreserve(nbrefl);
                  this->SetHklIobs(h,k,l,iobs,sigma);
                  this->SetWeight(weight);
               }//HKLIobsSigmaWeight
            }//ParList
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
}
#endif
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
   switch(mBackgroundType)
   {
      case POWDER_BACKGROUND_CUBIC_SPLINE: tag.AddAttribute("Interpolation","Spline");break;
      case POWDER_BACKGROUND_LINEAR: tag.AddAttribute("Interpolation","Linear");break;
   }
   os <<tag<<endl;
   indent++;
   
   XMLCrystTag tag2("TThetaIntensityList");
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag2<<endl;
   
   for(long j=0;j<mBackgroundNbPoint;j++)
   {
      for(int i=0;i<=indent;i++) os << "  " ;
      os << mBackgroundInterpPoint2Theta(j)*RAD2DEG <<" "
         << mBackgroundInterpPointIntensity(j) <<" "
         << !this->GetPar(mBackgroundInterpPointIntensity.data()+j).IsFixed()<<" "
         <<endl;
   }
   
   tag2.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag2<<endl;
   
   
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
      {
         if("Spline"==tagg.GetAttributeValue(i)) mBackgroundType=POWDER_BACKGROUND_CUBIC_SPLINE;
         if("Linear"==tagg.GetAttributeValue(i)) mBackgroundType=POWDER_BACKGROUND_LINEAR;
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
      if("TThetaIntensityList"==tag.GetName())
      {
         mBackgroundNbPoint=0;
         mBackgroundInterpPoint2Theta.resize(100);
         mBackgroundInterpPointIntensity.resize(100);
         CrystVector_bool fix(100);
         do
         {
            is >>mBackgroundInterpPoint2Theta(mBackgroundNbPoint)
               >>mBackgroundInterpPointIntensity(mBackgroundNbPoint)
               >>fix(mBackgroundNbPoint);
            mBackgroundNbPoint++;
            if(mBackgroundNbPoint==mBackgroundInterpPoint2Theta.numElements())
            {
               mBackgroundInterpPoint2Theta.
                        resizeAndPreserve(mBackgroundNbPoint+100);
               mBackgroundInterpPointIntensity.
                        resizeAndPreserve(mBackgroundNbPoint+100);
               fix.resizeAndPreserve(mBackgroundNbPoint+100);
            }
            while(0==isgraph(is.peek())) is.get();//Why do I need that ?
            //cout << is.peek()<<" "<<nbrefl<<endl;
         }
         while(is.peek()!='<');//until next tag
         mBackgroundInterpPoint2Theta.resizeAndPreserve(mBackgroundNbPoint);
         mBackgroundInterpPointIntensity.resizeAndPreserve(mBackgroundNbPoint);
         mBackgroundInterpPoint2Theta *= DEG2RAD;
         mClockBackgroundPoint.Click();

         //Rebuild parameter list
         this->ResetParList();
         char buf[25];
         for(int j=0;j<mBackgroundNbPoint;j++)
         {
            sprintf(buf,"Background_Point_%d",j);
            RefinablePar tmp(buf,
                             mBackgroundInterpPointIntensity.data()+j,0.,1000.,
                             gpRefParTypeScattDataBackground,REFPAR_DERIV_STEP_RELATIVE,
                             false,true,true,false,1.);
            tmp.AssignClock(mClockBackgroundPoint);
            tmp.SetDerivStep(1e-3);
            tmp.SetIsFixed(!fix(j));
            this->AddPar(tmp);
         }
         //read closing tag
         XMLCrystTag junkEndTag(is);
      }
   }
}
#if 0
void PowderPatternBackground::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("PowderPatternBackground::XMLInput():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            tag.Print();
            if(tag.GetType()=="\\PowderPatternBackground")
            {
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="BackgroundType")
               {
                  string type;
                  IOCrystExtractNameSpace(is,type);
                  if(type=="Linear")
                  {
                     mBackgroundType=POWDER_BACKGROUND_LINEAR;
                     continue;
                  }
                  if(type=="Spline")
                  {
                     mBackgroundType=POWDER_BACKGROUND_CUBIC_SPLINE;
                     continue;
                  }
                  continue;
               }
            }
            if(tag.GetType()=="ParList")
            {
               if(tag.GetName()=="BackgroundPoints2ThetaIntensity")
               {
                  mBackgroundNbPoint=0;
                  mBackgroundInterpPoint2Theta.resize(100);
                  mBackgroundInterpPointIntensity.resize(100);
                  CrystVector_bool fix(100);
                  do
                  {
                     is >>mBackgroundInterpPoint2Theta(mBackgroundNbPoint)
                        >>mBackgroundInterpPointIntensity(mBackgroundNbPoint)
                        >>fix(mBackgroundNbPoint);
                     mBackgroundNbPoint++;
                     if(mBackgroundNbPoint==mBackgroundInterpPoint2Theta.numElements())
                     {
                        mBackgroundInterpPoint2Theta.
                                 resizeAndPreserve(mBackgroundNbPoint+100);
                        mBackgroundInterpPointIntensity.
                                 resizeAndPreserve(mBackgroundNbPoint+100);
                        fix.resizeAndPreserve(mBackgroundNbPoint+100);
                     }
                     while(0==isgraph(is.peek())) is.get();
                     //cout << is.peek()<<" "<<nbrefl<<endl;
                  }
                  while(is.peek()!='<');//until next tag
                  mBackgroundInterpPoint2Theta.resizeAndPreserve(mBackgroundNbPoint);
                  mBackgroundInterpPointIntensity.resizeAndPreserve(mBackgroundNbPoint);
                  mBackgroundInterpPoint2Theta *= DEG2RAD;
                  mClockBackgroundPoint.Click();
                  
                  //Rebuild parameter list
                  this->ResetParList();
                  char buf[25];
                  for(int j=0;j<mBackgroundNbPoint;j++)
                  {
                     sprintf(buf,"Background_Point_%d",j);
                     RefinablePar tmp(buf,
                                      mBackgroundInterpPointIntensity.data()+j,0.,1000.,
                                      gpRefParTypeScattDataBackground,REFPAR_DERIV_STEP_RELATIVE,
                                      false,true,true,false,1.);
                     tmp.AssignClock(mClockBackgroundPoint);
                     tmp.SetDerivStep(1e-3);
                     tmp.SetIsFixed(!fix(j));
                     this->AddPar(tmp);
                  }

                  //cout<<FormatVertVector<REAL>(mBackgroundInterpPoint2Theta,
                  //                              mBackgroundInterpPointIntensity);
               }//BackgroundPoints2ThetaIntensity
            }//ParList
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
   VFN_DEBUG_MESSAGE("PowderPatternBackground::XMLInput():End",5)
}
#endif
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
   //:TODO: Asymmetry parameter, and put the following as options...
   tag.AddAttribute("Name",this->GetName());
   tag.AddAttribute("Crystal",this->GetCrystal().GetName());
   {
      stringstream ss;
      ss<<this->IsIgnoringImagScattFact();
      tag.AddAttribute("IgnoreImagScattFact",ss.str());
   }
  os <<tag<<endl;
   indent++;
   
   mReflectionProfileType.XMLOutput(os,indent);
   os<<endl;
   
   this->GetPar(&mCagliotiU).XMLOutput(os,"U",indent);
   os<<endl;
   
   this->GetPar(&mCagliotiV).XMLOutput(os,"V",indent);
   os<<endl;
   
   this->GetPar(&mCagliotiW).XMLOutput(os,"W",indent);
   os<<endl;
   
   this->GetPar(&mPseudoVoigtEta0).XMLOutput(os,"Eta0",indent);
   os<<endl;
   
   this->GetPar(&mPseudoVoigtEta1).XMLOutput(os,"Eta1",indent);
   os<<endl;

   this->GetPar(&mGlobalBiso).XMLOutput(os,"globalBiso",indent);
   os <<endl;
   
   if(mCorrTextureMarchDollase.GetNbPhase()>0)
   {
      mCorrTextureMarchDollase.XMLOutput(os,indent);
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
               if("U"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiU).XMLInput(is,tag);
                  break;
               }
               if("V"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiV).XMLInput(is,tag);
                  break;
               }
               if("W"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiW).XMLInput(is,tag);
                  break;
               }
               if("Eta0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta0).XMLInput(is,tag);
                  break;
               }
               if("Eta1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta1).XMLInput(is,tag);
                  break;
               }
               if("globalBiso"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mGlobalBiso).XMLInput(is,tag);
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
      if("TextureMarchDollase"==tag.GetName())
      {
         mCorrTextureMarchDollase.XMLInput(is,tag);
         continue;
      }
   }
}
#if 0
void PowderPatternDiffraction::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::XMLInput():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.IsClosingTag()==true)
            {
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="Crystal")
               {
                  string cryst;
                  IOCrystExtractNameQuoted(is,cryst);
                  //:TODO: check crystal was found
                  this->SetCrystal(gCrystalRegistry.GetObj(cryst));
                  continue;
               }
               if(tag.GetName()=="ProfileType")
               {
                  string type;
                  IOCrystExtractNameSpace(is,type);
                  if(type=="pearson-VII")
                  {
                     mReflectionProfileType.SetChoice(PROFILE_PEARSON_VII);
                     mClockProfilePar.Click();
                     continue;
                  }
                  if(type=="gaussian")
                  {
                     mReflectionProfileType.SetChoice(PROFILE_GAUSSIAN);
                     mClockProfilePar.Click();
                     continue;
                  }
                  if(type=="lorentzian")
                  {
                     mReflectionProfileType.SetChoice(PROFILE_LORENTZIAN);
                     mClockProfilePar.Click();
                     continue;
                  }
                  if(type=="pseudo-voigt")
                  {
                     mReflectionProfileType.SetChoice(PROFILE_PSEUDO_VOIGT);
                     mClockProfilePar.Click();
                     continue;
                  }
                  if(type=="pseudo-voigt-fcj")
                  {
                     mReflectionProfileType.SetChoice(PROFILE_PSEUDO_VOIGT_FINGER_COX_JEPHCOAT);
                     mClockProfilePar.Click();
                     continue;
                  }
                  continue;
               }
               if(tag.GetName()=="CagliotiUVW")
               {
                  is>>mCagliotiU>>mCagliotiV>>mCagliotiW;
                  mCagliotiU*=DEG2RAD*DEG2RAD;
                  mCagliotiV*=DEG2RAD*DEG2RAD;
                  mCagliotiW*=DEG2RAD*DEG2RAD;
                  mClockProfilePar.Click();
                  bool fix;
                  is>>fix;
                  this->GetPar(&mCagliotiU).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mCagliotiV).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mCagliotiW).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="PseudoVoigtEta0Eta1")
               {
                  is>>mPseudoVoigtEta0>>mPseudoVoigtEta1;
                  mClockProfilePar.Click();
                  bool fix;
                  is>>fix;
                  this->GetPar(&mPseudoVoigtEta0).SetIsFixed(!fix);
                  is>>fix;
                  this->GetPar(&mPseudoVoigtEta1).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="NeedLorentzCorr")
               {
                  is>>mNeedLorentzCorr;
                  mClockLorentzPolarSlitCorrPar.Reset();
                  continue;
               }
               if(tag.GetName()=="NeedPolarCorr-AFactor")
               {
                  is>>mNeedPolarCorr>>mPolarAfactor;
                  mClockLorentzPolarSlitCorrPar.Reset();
                  continue;
               }
               if(tag.GetName()=="NeedSlitApertureCorr")
               {
                  is>>mNeedSlitApertureCorr;
                  mClockLorentzPolarSlitCorrPar.Reset();
                  continue;
               }
               if(tag.GetName()=="IgnoreImagScattFact")
               {
                  bool b;
                  is>>b;
                  this->SetIsIgnoringImagScattFact(b);
                  continue;
               }
            }//Param
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::XMLInput():End",5)
}
#endif
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
   
   this->GetPar(&m2ThetaZero).XMLOutput(os,"2ThetaZero",indent);
   os <<endl;
   
   this->GetPar(&m2ThetaDisplacement).XMLOutput(os,"2ThetaDisplacement",indent);
   os <<endl;
   
   this->GetPar(&m2ThetaTransparency).XMLOutput(os,"2ThetaTransparency",indent);
   os <<endl;
   
   
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
   
   XMLCrystTag tag2("IobsSigmaWeightList");
      {
         stringstream ss;
         ss<<m2ThetaMin*RAD2DEG;
         tag2.AddAttribute("TThetaMin",ss.str());
      }
      {
         stringstream ss;
         ss<<m2ThetaStep*RAD2DEG;
         tag2.AddAttribute("TThetaStep",ss.str());
      }
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;


      for(unsigned long j=0;j<this->GetNbPoint();j++)
      {
         for(int i=0;i<=indent;i++) os << "  " ;
         os << mPowderPatternObs(j) <<" "
            << mPowderPatternObsSigma(j) <<" "
            << mPowderPatternWeight(j) <<" "
            <<endl;
      }
      tag2.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;
   
   for(int j=0;j<mExcludedRegionMin2Theta.numElements();j++)
   {
      XMLCrystTag tag3("Exclude2Theta");
      for(int i=0;i<indent;i++) os << "  " ;
      os << tag3 
         << mExcludedRegionMin2Theta(j)*RAD2DEG <<" "
         << mExcludedRegionMax2Theta(j)*RAD2DEG ;
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
               if("2ThetaZero"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&m2ThetaZero).XMLInput(is,tag);
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
            }
         }
         continue;
      }
      if("PowderPatternBackground"==tag.GetName())
      {
         PowderPatternBackground *comp=new PowderPatternBackground;
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternCrystal"==tag.GetName())
      {
         PowderPatternDiffraction *comp=new PowderPatternDiffraction;
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
      if("Exclude2Theta"==tag.GetName())
      {
         float min,max;
         is>>min>>max;
         this->Add2ThetaExcludedRegion(min*DEG2RAD,max*DEG2RAD);
         XMLCrystTag end(is);
         continue;
      }
      if("IobsSigmaWeightList"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Loading Iobs-Sigma-Weight List...",8);
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("TThetaMin"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>m2ThetaMin;
               VFN_DEBUG_MESSAGE("2Theta min="<<m2ThetaMin,8);
               m2ThetaMin*=DEG2RAD;
            }
            if("TThetaStep"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>m2ThetaStep;
               VFN_DEBUG_MESSAGE("2Theta step="<<m2ThetaStep<<tag.GetAttributeValue(i),8);
               m2ThetaStep*=DEG2RAD;
            }
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

         mPowderPatternObs.resizeAndPreserve(mNbPoint);
         //cout << mPowderPatternObs.numElements()<<" "<<mNbPoint<<" "<<this<<endl;
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
         mPowderPatternWeight.resizeAndPreserve(mNbPoint);
         
         XMLCrystTag junk(is);
         VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
         continue;
      }
   }
}
#if 0
void PowderPattern::XMLInputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("PowderPattern::XMLInput():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            //tag.Print();
            if(tag.IsClosingTag()==true)
            {
               VFN_DEBUG_MESSAGE("PowderPattern::XMLInput():End",5)
               this->UpdateDisplay();
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="PowderPatternComponent")
               {
                  string comp;
                  IOCrystExtractNameQuoted(is,comp);
                  REAL scale;
                  is >> scale;
                  //:TODO: Check we found the component
                  this->AddPowderPatternComponent(
                           gPowderPatternComponentRegistry.GetObj(comp));
                  mScaleFactor(mPowderPatternComponentRegistry.GetNb()-1)=scale;
                  continue;
               }
               if(tag.GetName()=="2thetaMin-Step")
               {
                  is>>m2ThetaMin>>m2ThetaStep;
                  m2ThetaMin*=DEG2RAD;
                  m2ThetaStep*=DEG2RAD;
                  mClockPowderPatternPar.Click();
                  continue;
               }
               if(tag.GetName()=="UseFastLessPreciseFunc")
               {
                  bool b;
                  is>>b;
                  this->SetUseFastLessPreciseFunc(b);
                  continue;
               }
               if(tag.GetName()=="2ThetaZero")
               {
                  is>>m2ThetaZero;
                  m2ThetaZero*=DEG2RAD;
                  bool fix;
                  is>>fix;
                  this->GetPar(&m2ThetaZero).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="2ThetaDisplacement")
               {
                  is>>m2ThetaDisplacement;
                  m2ThetaDisplacement*=DEG2RAD;
                  bool fix;
                  is>>fix;
                  this->GetPar(&m2ThetaDisplacement).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="2ThetaTransparency")
               {
                  is>>m2ThetaTransparency;
                  m2ThetaTransparency*=DEG2RAD;
                  bool fix;
                  is>>fix;
                  this->GetPar(&m2ThetaTransparency).SetIsFixed(!fix);
                  continue;
               }
               if(tag.GetName()=="StatisticsExcludeBackground")
               {
                  is>>mStatisticsExcludeBackground;
                  continue;
               }
            }
            if(tag.GetType()=="PowderPatternBackground")
            {
               PowderPatternBackground* tmp=new PowderPatternBackground;
               tmp->XMLInputOld(is,tag);
            }
            if(tag.GetType()=="PowderPatternDiffraction")
            {
               PowderPatternDiffraction* tmp=new PowderPatternDiffraction;
               tmp->XMLInputOld(is,tag);
            }
            if(tag.GetType()=="Radiation")
            {
               Radiation rad;
               rad.XMLInputOld(is,tag);
               this->SetRadiation(rad);
               continue;
            }
            if(tag.GetType()=="ParList")
            {
               if(tag.GetName()=="IobsSigmaWeight")
               {
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
                     if(mNbPoint==(unsigned long)mPowderPatternObs.numElements())
                     {
                        mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
                        mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
                        mPowderPatternWeight.resizeAndPreserve(mNbPoint+500);
                     }
                     while(0==isgraph(is.peek())) is.get();
                     //cout << is.peek()<<" "<<nbrefl<<endl;
                  }
                  while(is.peek()!='<');//until next tag
                  
                  mPowderPatternObs.resizeAndPreserve(mNbPoint);
                  //cout << mPowderPatternObs.numElements()<<" "<<mNbPoint<<" "<<this<<endl;
                  mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
                  mPowderPatternWeight.resizeAndPreserve(mNbPoint);
               }//IobsSigmaWeight
            }//ParList
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
}
#endif
} //namespace
