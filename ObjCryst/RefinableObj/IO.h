#ifndef _REFOBJ_IO_H_
#define _REFOBJ_IO_H_

#include <iostream>
#include <string>

namespace ObjCryst
{
class XMLCrystTag;
}
#ifdef __WX__CRYST__
   #include "wxCryst/wxCryst.h"
#endif

namespace ObjCryst
{
#ifdef __WX__CRYST__
/** \brief wxWindows representation of a XMLCrystTag (not implemented yet !)
*
* This will be used to choose objects to import from a save file.
*/
class WXXMLCrystTag: public WXCrystObj
{
   public:
      WXXMLCrystTag(wxWindow *parent, XMLCrystTag*);
      virtual void CrystUpdate();
      virtual void SetObjName(const string&);
      virtual string GetObjName()const;
      virtual bool Show(const bool);
   private:
      XMLCrystTag* mpTag;
};
#endif
/** \brief class to input or output a well-formatted xml beginning or ending tag.
* 
*/
class XMLCrystTag
{
   public:
      XMLCrystTag();
      XMLCrystTag(istream &is);
      XMLCrystTag(const string &tagName,const bool isEndTag=false, const bool isEmptyTag=false);
      ~XMLCrystTag();
      const string& GetName()const;
      const string& GetClassName()const;
      unsigned int GetNbAttribute()const;
      void AddAttribute(const string &attName,const string &attValue);
      void GetAttribute(const int attNum,string &attName,string &attValue);
      const string& GetAttributeName(const int attNum)const;
      const string& GetAttributeValue(const int attNum)const;
      void SetIsEndTag(const bool isEndTag);
      bool IsEndTag()const;
      void SetIsEmptyTag(const bool isEmptyTag);
      bool IsEmptyTag()const;
      void Print()const;
   private:
      string mName;
      bool mIsEndTag;
      bool mIsEmptyTag;
      unsigned int mNbAttribute;
      string mAttributeName [20];
      string mAttributeValue[20];
      friend ostream& operator<< (ostream&, const XMLCrystTag&);
      friend istream& operator>> (istream&, XMLCrystTag&);
   #ifdef __WX__CRYST__
   public:
      /// Create a WXCrystObj for this object. (not implemented yet)
      WXCrystObj* WXCreate(wxWindow*);
      WXCrystObj* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXXMLCrystTag *mpWXXMLCrystTag;
   #endif
};

/// Output an XMLCrystTag to a stream
ostream& operator<< (ostream&, const XMLCrystTag&);
/// Input an XMLCrystTag from a stream
istream& operator>> (istream&, XMLCrystTag&);

#if 0
//OLD

void IOCrystExtractNameSpace(istream &is,string &str);
void IOCrystExtractNameLine(istream &is,string &str);
void IOCrystExtractNameQuoted(istream &is,string &str);
void IOCrystXMLOutputNameQuoted(ostream &os,const string &str);

#ifdef __WX__CRYST__
class IOCrystTag;

class WXIOCrystTag: public WXCrystObj
{
   public:
      WXIOCrystTag(wxWindow *parent, IOCrystTag*);
      virtual void CrystUpdate();
      virtual void SetObjName(const string&);
      virtual string GetObjName()const;
      virtual bool Show(const bool);
   private:
      IOCrystTag* mpTag;
};
#endif
/// OLD
/// \internal Tag used to delimitate objects (example: "<Crystal 0>")
/// This includes the name of the corresponding ObjCryst class, as well
/// as a version number for the description.
class IOCrystTag
{
   public:
      IOCrystTag(const string& type,const string& name, const unsigned long version=0);
      IOCrystTag(istream &is);
      virtual ~IOCrystTag();
      void XMLInput(istream &is);
      bool operator==(const IOCrystTag&)const;
      const string &GetType()const;
      const string &GetName()const;
      unsigned long GetVersion()const;
      bool IsClosingTag()const;
      void Print()const;
      /// This is the same as GetType(), but allows using a registry
      const string &GetClassName()const;
   private:
      string mTagType;
      string mTagName;
      unsigned long mTagVersion;
      bool mIsClosingTag;
   #ifdef __WX__CRYST__
   public:
      /// Create a WXCrystObj for this object.
      virtual WXCrystObj* WXCreate(wxWindow*);
      WXCrystObj* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXIOCrystTag *mpWXIOCrystTag;
   #endif
};
#endif 

}//namespace ObjCryst

#endif //_REFOBJ_IO_H_
