/*  ObjCryst++ Object-Oriented Crystallographic Library
 (c) 2000- Vincent Favre-Nicolin vincefn@gmail.com
 
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

#include <cstdlib>
#include <boost/date_time.hpp>
#include <sstream>
#include "ObjCryst/ObjCryst/Undo.h"

namespace ObjCryst
{
XMLConfigHistory gConfigHistory;

////////////////////////////////////////////////////////////////////////
//
//    XMLConfig
//
////////////////////////////////////////////////////////////////////////
XMLConfig::XMLConfig()
{
   mClock.Click();
   mTime = boost::posix_time::second_clock::local_time();
   
   for(std::vector<Crystal*>::const_iterator pos=gCrystalRegistry.begin();pos!=gCrystalRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      stringstream ss;
      (*pos)->XMLOutput(ss,0);
      mvpCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
   }
   
   for(std::vector<DiffractionDataSingleCrystal*>::const_iterator pos=gDiffractionDataSingleCrystalRegistry.begin();pos!=gDiffractionDataSingleCrystalRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      stringstream ss;
      (*pos)->XMLOutput(ss,0);
      mvpDiffractionDataSingleCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
   }
   
   for(std::vector<PowderPattern*>::const_iterator pos=gPowderPatternRegistry.begin();pos!=gPowderPatternRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      stringstream ss;
      (*pos)->XMLOutput(ss,0);
      mvpPowderPatternXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
   }
   
   for(std::vector<OptimizationObj*>::const_iterator pos=gOptimizationObjRegistry.begin();pos!=gOptimizationObjRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      stringstream ss;
      (*pos)->XMLOutput(ss,0);
      mvpOptimizationObjXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
   }
   mClock.Click();
}

XMLConfig::XMLConfig(const XMLConfig &old)
{
   mTime = boost::posix_time::second_clock::local_time();
   
   for(std::vector<Crystal*>::const_iterator pos=gCrystalRegistry.begin();pos!=gCrystalRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(&old): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      std::map<Crystal*, boost::shared_ptr<std::string> >::const_iterator pos2 = old.mvpCrystalXML.find(*pos);
      if(pos2 != old.mvpCrystalXML.end())
      {
         if((*pos)->GetClockMaster() <= old.GetClock())
            mvpCrystalXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
         else
         {
            stringstream ss;
            (*pos)->XMLOutput(ss,0);
            if(ss.str() == *(pos2->second))
               mvpCrystalXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
            else
               mvpCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
         }
      }
      else
      {
         stringstream ss;
         (*pos)->XMLOutput(ss,0);
         mvpCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
      }
   }
   
   for(std::vector<DiffractionDataSingleCrystal*>::const_iterator pos=gDiffractionDataSingleCrystalRegistry.begin();pos!=gDiffractionDataSingleCrystalRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(&old): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      std::map<DiffractionDataSingleCrystal*, boost::shared_ptr<std::string> >::const_iterator pos2 = old.mvpDiffractionDataSingleCrystalXML.find(*pos);
      if(pos2 != old.mvpDiffractionDataSingleCrystalXML.end())
      {
         if((*pos)->GetClockMaster() <= old.GetClock())
            mvpDiffractionDataSingleCrystalXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
         else
         {
            stringstream ss;
            (*pos)->XMLOutput(ss,0);
            if(ss.str() == *(pos2->second))
               mvpDiffractionDataSingleCrystalXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
            else
               mvpDiffractionDataSingleCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
         }
      }
      else
      {
         stringstream ss;
         (*pos)->XMLOutput(ss,0);
         mvpDiffractionDataSingleCrystalXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
      }
   }

   for(std::vector<PowderPattern*>::const_iterator pos=gPowderPatternRegistry.begin();pos!=gPowderPatternRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(&old): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      std::map<PowderPattern*, boost::shared_ptr<std::string> >::const_iterator pos2 = old.mvpPowderPatternXML.find(*pos);
      if(pos2 != old.mvpPowderPatternXML.end())
      {
         if((*pos)->GetClockMaster() <= old.GetClock())
            mvpPowderPatternXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
         else
         {
            stringstream ss;
            (*pos)->XMLOutput(ss,0);
            if(ss.str() == *(pos2->second))
               mvpPowderPatternXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
            else
               mvpPowderPatternXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
         }
      }
      else
      {
         stringstream ss;
         (*pos)->XMLOutput(ss,0);
         mvpPowderPatternXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
      }
   }

   for(std::vector<OptimizationObj*>::const_iterator pos=gOptimizationObjRegistry.begin();pos!=gOptimizationObjRegistry.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::XMLConfig(&old): storing "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      std::map<OptimizationObj*, boost::shared_ptr<std::string> >::const_iterator pos2 = old.mvpOptimizationObjXML.find(*pos);
      if(pos2 != old.mvpOptimizationObjXML.end())
      {
         stringstream ss;
         (*pos)->XMLOutput(ss,0);
         if(ss.str() == *(pos2->second))
            mvpOptimizationObjXML[*pos] = boost::shared_ptr<std::string>(pos2->second);
         else
            mvpOptimizationObjXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
      }
      else
      {
         stringstream ss;
         (*pos)->XMLOutput(ss,0);
         mvpOptimizationObjXML[*pos] = boost::shared_ptr<std::string>(new std::string(ss.str()));
      }
   }
   mClock.Click();
}
   
void XMLConfig::Restore()
{
   VFN_DEBUG_ENTRY("XMLConfig::Restore()", 10)
   std::list<RefinableObj*> vpRefObjUpdated;
   std::list<OptimizationObj*> vpOptimObjUpdated;
   std::list<RefinableObj*> vpRefObjDeleted;
   std::list<OptimizationObj*> vpOptimObjDeleted;
   // Restore objects configuration and load new ones
   for(std::map<Crystal*, boost::shared_ptr<std::string> >::const_iterator pos= mvpCrystalXML.begin(); pos!=mvpCrystalXML.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::Restore(): found "<< pos->first->GetClassName() << ":" << pos->first->GetName(), 10)
      Crystal *p;
      if(gCrystalRegistry.Find(pos->first) >= 0)
      {
         p = pos->first;
         stringstream ss;
         p->XMLOutput(ss,0);
         if(ss.str()== *(pos->second)) continue;
      }
      else
         p = new Crystal;
      stringstream ss;
      ss << *(pos->second);
      VFN_DEBUG_MESSAGE("XMLConfig::Restore():Crystal:xml="<<*(pos->second), 5)
      VFN_DEBUG_MESSAGE("XMLConfig::Restore():Crystal:xml="<<ss.str(), 5)
      XMLCrystTag tag(ss);
      VFN_DEBUG_MESSAGE("XMLConfig::Restore():Crystal:tag="<<tag, 5)
      VFN_DEBUG_ENTRY("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      p->XMLInput(ss,tag);
      VFN_DEBUG_EXIT("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      vpRefObjUpdated.push_back(p);
   }

   for(std::map<DiffractionDataSingleCrystal*, boost::shared_ptr<std::string> >::const_iterator pos= mvpDiffractionDataSingleCrystalXML.begin(); pos!=mvpDiffractionDataSingleCrystalXML.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::Restore(): found "<< pos->first->GetClassName() << ":" << pos->first->GetName(), 10)
      DiffractionDataSingleCrystal *p;
      if(gDiffractionDataSingleCrystalRegistry.Find(pos->first) >= 0)
      {
         p = pos->first;
         stringstream ss;
         p->XMLOutput(ss,0);
         if(ss.str()== *(pos->second)) continue;
      }
      else
         p = new DiffractionDataSingleCrystal;
      stringstream ss;
      ss << pos->second;
      XMLCrystTag tag(ss);
      VFN_DEBUG_ENTRY("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      p->XMLInput(ss,tag);
      VFN_DEBUG_EXIT("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      vpRefObjUpdated.push_back(p);
   }
   
   for(std::map<PowderPattern*, boost::shared_ptr<std::string> >::const_iterator pos= mvpPowderPatternXML.begin(); pos!=mvpPowderPatternXML.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::Restore(): found "<< pos->first->GetClassName() << ":" << pos->first->GetName(), 10)
      PowderPattern *p;
      if(gPowderPatternRegistry.Find(pos->first) >= 0)
      {
         p = pos->first;
         stringstream ss;
         p->XMLOutput(ss,0);
         if(ss.str()== *(pos->second)) continue;
      }
      else
         p = new PowderPattern;
      stringstream ss;
      ss << pos->second;
      XMLCrystTag tag(ss);
      VFN_DEBUG_ENTRY("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      p->XMLInput(ss,tag);
      VFN_DEBUG_EXIT("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      vpRefObjUpdated.push_back(p);
   }
   
   for(std::map<OptimizationObj*, boost::shared_ptr<std::string> >::const_iterator pos= mvpOptimizationObjXML.begin(); pos!=mvpOptimizationObjXML.end();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::Restore(): found "<< pos->first->GetClassName() << ":" << pos->first->GetName(), 10)
      OptimizationObj *p;
      if(gOptimizationObjRegistry.Find(pos->first) >= 0)
      {
         p = pos->first;
         stringstream ss;
         p->XMLOutput(ss,0);
         if(ss.str()== *(pos->second)) continue;
      }
      else
         p = new MonteCarloObj;
      stringstream ss;
      ss << pos->second;
      XMLCrystTag tag(ss);
      VFN_DEBUG_ENTRY("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      p->XMLInput(ss,tag);
      VFN_DEBUG_EXIT("XMLConfig::Restore()"<<p->GetClassName()<<":"<<p->GetName(), 10)
      vpOptimObjUpdated.push_back(p);
   }

   // List objects which don't exist in the loaded configuration and which should be deleted
   for(std::vector<Crystal*>::const_iterator pos=gCrystalRegistry.begin();pos!=gCrystalRegistry.end();++pos)
      if(mvpCrystalXML.find(*pos) == mvpCrystalXML.end()) vpRefObjDeleted.push_back(*pos);

   for(std::vector<DiffractionDataSingleCrystal*>::const_iterator pos=gDiffractionDataSingleCrystalRegistry.begin();pos!=gDiffractionDataSingleCrystalRegistry.end();++pos)
      if(mvpDiffractionDataSingleCrystalXML.find(*pos) == mvpDiffractionDataSingleCrystalXML.end()) vpRefObjDeleted.push_back(*pos);

   for(std::vector<PowderPattern*>::const_iterator pos=gPowderPatternRegistry.begin();pos!=gPowderPatternRegistry.end();++pos)
      if(mvpPowderPatternXML.find(*pos) == mvpPowderPatternXML.end()) vpRefObjDeleted.push_back(*pos);

   for(std::vector<OptimizationObj*>::const_iterator pos=gOptimizationObjRegistry.begin();pos!=gOptimizationObjRegistry.end();++pos)
      if(mvpOptimizationObjXML.find(*pos) == mvpOptimizationObjXML.end()) vpOptimObjDeleted.push_back(*pos);

   // Objects are destroyed after listing all, to avoid unvalidating the ObjRegistry containers (std::vector)
   // First optimisation obj, then diffraction data, then crystals
   for(std::list<OptimizationObj*>::const_iterator pos = vpOptimObjDeleted.begin();pos!=vpOptimObjDeleted.end();++pos)
      {
         VFN_DEBUG_MESSAGE("XMLConfig::Restore(): deleting "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
         delete *pos;
      }
   for(std::list<RefinableObj*>::const_reverse_iterator pos = vpRefObjDeleted.rbegin();pos!=vpRefObjDeleted.rend();++pos)
   {
      VFN_DEBUG_MESSAGE("XMLConfig::Restore(): deleting "<< (*pos)->GetClassName() << ":" << (*pos)->GetName(), 10)
      delete *pos;
   }

   // Update display. Order should not matter
   for(std::list<OptimizationObj*>::const_iterator pos = vpOptimObjUpdated.begin() ; pos!=vpOptimObjUpdated.end() ; ++pos) (*pos)->UpdateDisplay();
   for(std::list<RefinableObj*>::const_iterator pos = vpRefObjUpdated.begin() ; pos!=vpRefObjUpdated.end() ; ++pos) (*pos)->UpdateDisplay();
   VFN_DEBUG_EXIT("XMLConfig::Restore()", 10)
}

bool XMLConfig::operator==(const ObjCryst::XMLConfig &rhs) const
{
   if(mvpCrystalXML.size() != rhs.mvpCrystalXML.size()) return false;
   if(mvpDiffractionDataSingleCrystalXML.size() != rhs.mvpDiffractionDataSingleCrystalXML.size()) return false;
   if(mvpPowderPatternXML.size() != rhs.mvpPowderPatternXML.size()) return false;
   if(mvpOptimizationObjXML.size() != rhs.mvpOptimizationObjXML.size()) return false;
   
   for(std::map<Crystal*, boost::shared_ptr<std::string> >::const_iterator pos = mvpCrystalXML.begin(); pos!=mvpCrystalXML.end(); pos++)
   {
      std::map<Crystal*, boost::shared_ptr<std::string> >::const_iterator pos2 = rhs.mvpCrystalXML.find(pos->first);
      if(pos2==rhs.mvpCrystalXML.end()) return false;
      if(*(pos->second) != *(pos2->second)) return false;
   }

   for(std::map<DiffractionDataSingleCrystal*, boost::shared_ptr<std::string> >::const_iterator pos = mvpDiffractionDataSingleCrystalXML.begin(); pos!=mvpDiffractionDataSingleCrystalXML.end(); pos++)
   {
      std::map<DiffractionDataSingleCrystal*, boost::shared_ptr<std::string> >::const_iterator pos2 = rhs.mvpDiffractionDataSingleCrystalXML.find(pos->first);
      if(pos2==rhs.mvpDiffractionDataSingleCrystalXML.end()) return false;
      if(*(pos->second) != *(pos2->second)) return false;
   }

   for(std::map<PowderPattern*, boost::shared_ptr<std::string> >::const_iterator pos = mvpPowderPatternXML.begin(); pos!=mvpPowderPatternXML.end(); pos++)
   {
      std::map<PowderPattern*, boost::shared_ptr<std::string> >::const_iterator pos2 = rhs.mvpPowderPatternXML.find(pos->first);
      if(pos2==rhs.mvpPowderPatternXML.end()) return false;
      if(*(pos->second) != *(pos2->second)) return false;
   }

   for(std::map<OptimizationObj*, boost::shared_ptr<std::string> >::const_iterator pos = mvpOptimizationObjXML.begin(); pos!=mvpOptimizationObjXML.end(); pos++)
   {
      std::map<OptimizationObj*, boost::shared_ptr<std::string> >::const_iterator pos2 = rhs.mvpOptimizationObjXML.find(pos->first);
      if(pos2==rhs.mvpOptimizationObjXML.end()) return false;
      if(*(pos->second) != *(pos2->second)) return false;
   }
   return true;
}

bool XMLConfig::operator!=(const ObjCryst::XMLConfig &rhs) const
{
   return !(this->operator==(rhs));
}

const RefinableObjClock& XMLConfig::GetClock() const { return mClock;}


////////////////////////////////////////////////////////////////////////
//
//    XMLConfigHistory
//
////////////////////////////////////////////////////////////////////////
XMLConfigHistory::XMLConfigHistory(const int maxnb):
mMaxNb(maxnb), mpCurrentConfig(), mLock(false)
{}

bool XMLConfigHistory::Store()
{
   if(mLock) return false;
   mLock=true;
   VFN_DEBUG_ENTRY("XMLConfigHistory::Store()", 10)
   if(mvpConfig.size()==0) mvpConfig.push_front(boost::shared_ptr<XMLConfig>(new XMLConfig));
   else
   {
      boost::shared_ptr<XMLConfig> previous = mpCurrentConfig;
      if(mpCurrentConfig==NULL) // C++11 could use nullptr
         previous = mvpConfig.front();
      boost::shared_ptr<XMLConfig> newconf(new XMLConfig(*previous));
      if(*newconf == *previous)
      {
         VFN_DEBUG_EXIT("XMLConfigHistory::Store(): nothing new", 10)
         mLock=false;
         return false;
      }
      mvpConfig.push_front(newconf);
   }
   mpCurrentConfig = boost::shared_ptr<XMLConfig>();
   if(mvpConfig.size()>mMaxNb) mvpConfig.pop_back();
   (*fpObjCrystInformUser)("XMLConfigHistory::Store()");
   
   VFN_DEBUG_EXIT("XMLConfigHistory::Store(): new config", 10)
   mLock=false;
   return true;
}
   
void XMLConfigHistory::Restore(boost::shared_ptr<XMLConfig> &p)
{
   mLock=true;
   p->Restore();
   mpCurrentConfig = p;
   mLock=false;
}

std::list<boost::shared_ptr<XMLConfig> > XMLConfigHistory::GetList() { return mvpConfig;}

const std::list<boost::shared_ptr<XMLConfig> > XMLConfigHistory::GetList() const { return mvpConfig;}

bool XMLConfigHistory::Previous()
{
   if(mvpConfig.size()==0) return false;
   if(mpCurrentConfig == mvpConfig.back())
   {
      fpObjCrystInformUser("XMLConfigHistory::Previous(): already at last stored configuration");
      return false;
   }
   mLock=true;
   boost::shared_ptr<XMLConfig> p;
   if(mpCurrentConfig==NULL) p = mvpConfig.front();
   else
   {
      std::list<boost::shared_ptr<XMLConfig> >::const_iterator pos;
      for(pos=mvpConfig.begin();pos!=mvpConfig.end();++pos)
         if(*pos == mpCurrentConfig) {++pos ; break;}
      p = *pos;
   }
   (*fpObjCrystInformUser)("XMLConfigHistory::Previous()");
   mpCurrentConfig = p;
   p->Restore();
   mLock=false;
   return true;
}

bool XMLConfigHistory::Next()
{
   if((mvpConfig.size()==0) || (mpCurrentConfig==NULL)) return false;
   if(mpCurrentConfig == mvpConfig.front())
   {
      fpObjCrystInformUser("XMLConfigHistory::Next(): already at first stored configuration");
      return false;
   }
   mLock=true;
   (*fpObjCrystInformUser)("XMLConfigHistory::Next()");

   std::list<boost::shared_ptr<XMLConfig> >::const_reverse_iterator pos;
   for(pos=mvpConfig.rbegin();pos!=mvpConfig.rend();++pos)
      if(*pos == mpCurrentConfig) {++pos ; break;}
   (*pos)->Restore();
   mpCurrentConfig = *pos;
   mLock=false;
   return true;
}


} // namespace ObjCryst
