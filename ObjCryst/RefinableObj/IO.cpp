/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: vincefn@users.sourceforge.net
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*
*  source file for generic Input/Output in 
*
*/


#include "RefinableObj/RefinableObj.h"
#include "RefinableObj/IO.h"
#include <sstream>

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    XMLCrystTag
//
////////////////////////////////////////////////////////////////////////
XMLCrystTag::XMLCrystTag():
mIsEndTag(false),mIsEmptyTag(false),mNbAttribute(0)
{}

XMLCrystTag::XMLCrystTag(istream &is):
mNbAttribute(0)
{
   is >> *this;
}

XMLCrystTag::XMLCrystTag(const string &tagName,
                         const bool isEndTag, 
                         const bool isEmptyTag):
mName(tagName),mIsEndTag(isEndTag),mIsEmptyTag(isEmptyTag),mNbAttribute(0)
{}

XMLCrystTag::~XMLCrystTag(){}

const string& XMLCrystTag::GetName()const{return mName;}
const string& XMLCrystTag::GetClassName()const{static string str="XMLCrystTag";return str;}

unsigned int XMLCrystTag::GetNbAttribute()const{return mNbAttribute;}

void XMLCrystTag::AddAttribute(const string &attName,const string &attValue)
{
   mAttributeName [mNbAttribute  ]=attName;
   mAttributeValue[mNbAttribute++]=attValue;
}
void XMLCrystTag::GetAttribute(const int attNum,string &attName,string &attValue)
{
   attName=mAttributeName  [attNum];
   attValue=mAttributeValue[attNum];
}

const string& XMLCrystTag::GetAttributeName(const int attNum)const
{return mAttributeName[attNum];}

const string& XMLCrystTag::GetAttributeValue(const int attNum)const
{return mAttributeValue[attNum];}

void XMLCrystTag::SetIsEndTag(const bool isEndTag){mIsEndTag=isEndTag;}
bool XMLCrystTag::IsEndTag()const{return mIsEndTag;}
void XMLCrystTag::SetIsEmptyTag(const bool isEmptyTag){mIsEmptyTag=isEmptyTag;}
bool XMLCrystTag::IsEmptyTag()const{return mIsEmptyTag;}
void XMLCrystTag::Print()const{cout<<*this;}
#ifdef __WX__CRYST__
WXCrystObj* XMLCrystTag::WXCreate(wxWindow *parent)
{
   //:TODO:
   //mpWXXMLCrystTag=new WXXMLCrystTag (parent,this);
   return mpWXXMLCrystTag;
}
WXCrystObj* XMLCrystTag::WXGet()
{
   return mpWXXMLCrystTag;
}
void XMLCrystTag::WXDelete()
{
   if(0!=mpWXXMLCrystTag) delete mpWXXMLCrystTag;
   mpWXXMLCrystTag=0;
}
void XMLCrystTag::WXNotifyDelete()
{
   mpWXXMLCrystTag=0;
}
#endif

ostream& operator<< (ostream& os, const XMLCrystTag&tag)
{
   if(true==tag.mIsEndTag)
   {
      os <<"</"<<tag.mName<<">";
      return os;
   }
   os <<"<"<<tag.mName;
   for(unsigned int i=0;i<tag.mNbAttribute;i++)
   {
      os<<" "<<tag.mAttributeName[i]<<"=\""<<tag.mAttributeValue[i]<<"\"";
   }
   if(true==tag.mIsEmptyTag) os <<"/>";
   else os <<">";
   return os;
}
istream& operator>> (istream& is, XMLCrystTag &tag)
{
   ios::fmtflags f=is.flags();
   is.unsetf(ios::skipws);//skip whitespaces
   tag.mIsEmptyTag=false;
   tag.mIsEndTag=false;
   char tmp;
   is>>tmp;
   while ((tmp!='<') && !(is.eof()) )is>>tmp;
   if(is.eof()) return is;//:TODO: throw exception ?
   while ((tmp==' ')||(tmp=='<'))is>>tmp;
   
   if('/'==tmp)
   {
      tag.mIsEndTag=true;
      while ((tmp==' ')||(tmp=='/'))is>>tmp;
   }
   
   string str="";
   do {str+=tmp;is>>tmp;VFN_DEBUG_MESSAGE(str,1)} while ((tmp!=' ')&&(tmp!='>')&&(tmp!='/'));
   tag.mName=str;
   
   string str2;
   while(true)
   {
      while(tmp==' ')is>>tmp;
      if(tmp=='>')
      {
         is.setf(f);
         return is;
      }
      if(tmp=='/')
      {
         is>>tmp;
         //if(tmp!='>') ; :TODO: 
         tag.mIsEmptyTag=true;
         is.setf(f);
         return is;
      }
      str="";
      do {str+=tmp;is>>tmp;VFN_DEBUG_MESSAGE(str,1)} while ((tmp!=' ')&&(tmp!='='));
      while(tmp!='"')is>>tmp;
      str2="";
      is>>tmp;
      if(tmp!='"') do {str2+=tmp;is>>tmp;VFN_DEBUG_MESSAGE(str2,1)} while (tmp!='"');
      is>>tmp;
      
      tag.AddAttribute(str,str2);
   }
   
   is.setf(f);
   return is;
}
////////////////////////////////////////////////////////////////////////
//
//    I/O RefinablePar
//
////////////////////////////////////////////////////////////////////////
void RefinablePar::Output(ostream &os,const string &name,int indent)const
{
   VFN_DEBUG_ENTRY("RefinablePar::Output():"<<this->GetName(),5)
   XMLCrystTag tag("Par");
   {
      stringstream ss;
      ss <<!(this->IsFixed());
      tag.AddAttribute("Refined",ss.str());
   }
   {
      stringstream ss;
      ss <<this->IsLimited();
      tag.AddAttribute("Limited",ss.str());
   }
   {
      stringstream ss;
      ss <<this->GetHumanMin();
      tag.AddAttribute("Min",ss.str());
   }
   {
      stringstream ss;
      ss <<this->GetHumanMax();
      tag.AddAttribute("Max",ss.str());
   }
   #if 0
   {//Useless. Periodicity cannot be changed for a given parameter
      stringstream ss;
      ss <<this->IsPeriodic();
      tag.AddAttribute("Periodic",ss.str());
   }
   #endif
	//the name of the parameter is saved last to enhance readability of saved files
   tag.AddAttribute("Name",name);
   
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag;
   tag.SetIsEndTag(true);
   os <<this->GetHumanValue()<<tag;
   VFN_DEBUG_EXIT("RefinablePar::Output():"<<this->GetName(),5)
}
void RefinablePar::Output(ostream &os,int indent)const
{
   this->Output(os,mName,indent);
}

void RefinablePar::Input(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("RefinablePar::Input():"<<this->GetName(),5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Name"==tag.GetAttributeName(i)) continue;//names must be set by the object
      if("Refined"==tag.GetAttributeName(i))
      {
         bool b;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>b;
         this->SetIsFixed(!b);
         continue;
      }
      if("Limited"==tag.GetAttributeName(i))
      {
         bool b;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>b;
         this->SetIsLimited(b);
         continue;
      }
      if("Min"==tag.GetAttributeName(i))
      {
         double f;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>f;
         this->SetHumanMin(f);
         continue;
      }
      if("Max"==tag.GetAttributeName(i))
      {
         double f;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>f;
         this->SetHumanMax(f);
         continue;
      }
      if("Periodic"==tag.GetAttributeName(i))
      {
         bool b;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>b;
         this->SetIsPeriodic(b);
         continue;
      }
   }
   double f;
   is>>f;
   this->SetHumanValue(f);
   XMLCrystTag junk(is);//read end tag
   VFN_DEBUG_EXIT("RefinablePar::Input():"<<this->GetName(),5)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O RefObjOpt
//
////////////////////////////////////////////////////////////////////////
void RefObjOpt::Output(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("RefObjOpt::Output():"<<this->GetName(),5)
   XMLCrystTag tag("Option",false,true);
   tag.AddAttribute("Name",this->GetName());
   {
      stringstream ss;
      ss <<this->GetChoice();
      tag.AddAttribute("Choice",ss.str());
   }
   tag.AddAttribute("ChoiceName",this->GetChoiceName(this->GetChoice()));

   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag;
   
   VFN_DEBUG_EXIT("RefObjOpt::Output():"<<this->GetName(),5)
}

void RefObjOpt::Input(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("RefObjOpt::Input():"<<this->GetName(),5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Name"==tag.GetAttributeName(i)) continue;//names must be set by the object
      if("ChoiceName"==tag.GetAttributeName(i)) continue;//just for info
      if("Choice"==tag.GetAttributeName(i))
      {
         int b;
         stringstream ss(tag.GetAttributeValue(i));
         ss >>b;
         this->SetChoice(b);
         continue;
      }
   }
   VFN_DEBUG_EXIT("RefObjOpt::Input():"<<this->GetName(),5)
}
////////////////////////////////////////////////////////////////////////
//
//    I/O RefinableObj // Does nothing ! Should be purely virtual...
//
////////////////////////////////////////////////////////////////////////
void RefinableObj::Output(ostream &os,int indent)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::Output():"<<this->GetName(),5)
}

void RefinableObj::Input(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_MESSAGE("RefinableObj::Input():"<<this->GetName(),5)
}
#if 0
void RefinableObj::InputOld(istream &is,const IOCrystTag &tag)
{
   VFN_DEBUG_MESSAGE("RefinableObj::Input():"<<this->GetName(),5)
}
#endif
////////////////////////////////////////////////////////////////////////
//
//    OLD Functions & objects
//
////////////////////////////////////////////////////////////////////////
#if 0
void IOCrystExtractNameSpace(istream &is,string &str)
{
   ios::fmtflags f=is.flags();
   is.unsetf(ios::skipws);//skip whitespaces
   char tmp;
   do {is>>tmp;} while (tmp==' ');
   str="";
   do {str+=tmp;is>>tmp;VFN_DEBUG_MESSAGE(str,1)} while ((tmp!=' ')&&(tmp!='>'));
   if(tmp=='>') is.putback(tmp);
   is.setf(f);
}

void IOCrystExtractNameLine(istream &is,string &str)
{
   ios::fmtflags f=is.flags();
   is.setf(ios::skipws);//skip leading whitespaces
   getline(is,str);
   is.setf(f);
}

void IOCrystExtractNameQuoted(istream &is,string &str)
{
   char tmp;
   do {is>>tmp;} while ((tmp!='\'')&&(tmp!='\"'));
   str="";
   is>>tmp;
   if((tmp=='\'')||(tmp=='\"')) return;
   ios::fmtflags f=is.flags();
   is.unsetf(ios::skipws);//do not skip whitespaces
   do {str+=tmp;is>>tmp;VFN_DEBUG_MESSAGE(str,1)} while ((tmp!='\'')&&(tmp!='\"'));
   is.setf(f);
}

void IOCrystOutputNameQuoted(ostream &os,const string &str)
{
   os << '\"' << str << '\"';
}

////////////////////////////////////////////////////////////////////////
//
//    IOCrystTag
//
////////////////////////////////////////////////////////////////////////

IOCrystTag::IOCrystTag(const string& type,const string& name, const unsigned long version):
mTagType(type),mTagName(name),mTagVersion(version),mIsClosingTag(false)
{}
IOCrystTag::IOCrystTag(istream &is)
{
   VFN_DEBUG_MESSAGE("IOCrystTag::IOCrystTag(istream &is)",2)
   char tmp;
   do
   {
      is>>tmp;
      if(true==is.eof())
      {
         mTagVersion=0;//closing tag
         mTagName="";
         mTagType="";
         mIsClosingTag=true;
         return;
      }
   } while (tmp!='<');
   IOCrystExtractNameSpace(is,mTagType);
   VFN_DEBUG_MESSAGE("IOCrystTag::IOCrystTag(istream &is):TagType:"<<mTagType,1)
   if(*(mTagType.c_str())=='\\')
   {
      mTagVersion=0;//closing tag
      mTagName="";
      mIsClosingTag=true;
   }
   else
   {  
      IOCrystExtractNameQuoted(is,mTagName);
      is >> mTagVersion;
      mIsClosingTag=false;
   }
   do {is>>tmp;} while (tmp!='>');
   VFN_DEBUG_MESSAGE("IOCrystTag::IOCrystTag(istream &is):End",2)
}
IOCrystTag::~IOCrystTag(){}

bool IOCrystTag::operator==(const IOCrystTag& rhs)const
{
   if( (rhs.GetType()==this->GetType()) && (rhs.GetName()==this->GetName())) return true;
   return false;
}
void IOCrystTag::Input(istream &is)
{
   VFN_DEBUG_MESSAGE("IOCrystTag::Input(istream &is)",2)
   char tmp;
   do {is>>tmp;} while ((tmp!='<') && (!is.eof()) );
   if(is.eof()) return;
   IOCrystExtractNameSpace(is,mTagType);
   VFN_DEBUG_MESSAGE("IOCrystTag::Input(istream &is):TagType:"<<mTagType,1)
   if(*(mTagType.c_str())=='\\')
   {
      mTagVersion=0;//closing tag
      mTagName="";
      mIsClosingTag=true;
   }
   else
   {  
      IOCrystExtractNameQuoted(is,mTagName);
      is >> mTagVersion;
      mIsClosingTag=false;
   }
   do {is>>tmp;} while (tmp!='>');
   VFN_DEBUG_MESSAGE("IOCrystTag::Input(istream &is):End",2)
}
const string &IOCrystTag::GetType() const{return mTagType;}
const string &IOCrystTag::GetName() const{return mTagName;}
unsigned long IOCrystTag::GetVersion()const{return mTagVersion;}
bool IOCrystTag::IsClosingTag()const{return mIsClosingTag;}
void IOCrystTag::Print() const
{
   if(mIsClosingTag==true) cout <<"<"<<mTagType<<">"<<endl;
   else
      cout <<"<"<<mTagType<<" \""<<mTagName<<"\" "<<mTagVersion<<">"<<endl;
}
const string &IOCrystTag::GetClassName() const{return mTagType;}
#ifdef __WX__CRYST__
WXCrystObj* IOCrystTag::WXCreate(wxWindow *parent)
{
   //mpWXIOCrystTag=new WXIOCrystTag (parent,this);
   return mpWXIOCrystTag;
}
WXCrystObj* IOCrystTag::WXGet()
{
   return mpWXIOCrystTag;
}
void IOCrystTag::WXDelete()
{
   if(0!=mpWXIOCrystTag) delete mpWXIOCrystTag;
   mpWXIOCrystTag=0;
}
void IOCrystTag::WXNotifyDelete()
{
   mpWXIOCrystTag=0;
}
#endif
#endif

} //namespace
