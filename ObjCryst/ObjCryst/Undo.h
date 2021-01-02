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

#ifndef _OBJCRYST_UNDO_H_
#define _OBJCRYST_UNDO_H_

#include <string>
#include <map>
// With C++11, could used std shared_ptr instead of boost
#include <boost/shared_ptr.hpp>
#include <boost/date_time.hpp>

#include "ObjCryst/RefinableObj/IO.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"

namespace ObjCryst
{
/** Class to store & restore all top-level ObjCryst++ objects as their XML description.
 * Storage uses a shared_ptr to avoid useless copies.
 *
 */
class XMLConfig
{
   public:
      /// Default cronstructor
      XMLConfig();
      /// Constructor based on a previous instance
      XMLConfig(const XMLConfig &old);
      /// Restore the saved XML configuration
      void Restore();
      /// Equality operator
      bool operator==(const XMLConfig &rhs) const;
      /// Non-equality operator
      bool operator!=(const XMLConfig &rhs) const;
      /// Access the clock corresponding to the creation of this config
      const RefinableObjClock& GetClock() const;
   protected:
      std::map<Crystal*, boost::shared_ptr<std::string> > mvpCrystalXML;
      std::map<PowderPattern*, boost::shared_ptr<std::string> > mvpPowderPatternXML;
      std::map<DiffractionDataSingleCrystal*, boost::shared_ptr<std::string> > mvpDiffractionDataSingleCrystalXML;
      std::map<OptimizationObj*, boost::shared_ptr<std::string> > mvpOptimizationObjXML;
      boost::posix_time::ptime mTime;
      /// Clock recording the creation of this XMLConfig
      RefinableObjClock mClock;
};
   
/** Class to store XMLConfig configurations
 *
 */
class XMLConfigHistory
{
   public:
      /** Constructor
       * \param maxnb: maximum size of history
       */
      XMLConfigHistory(const int maxnb=20);
      /// Store a new configuration. If there is no difference with the previous one,
      /// nothing is done.
      /// \return: true if a new configuration was saved, false otherwise
      bool Store();
      /// Restore a configuration
      void Restore(boost::shared_ptr<XMLConfig> &p);
      /// Access the list of configurations
      std::list<boost::shared_ptr<XMLConfig> > GetList();
      /// Access the list of configurations
      const std::list<boost::shared_ptr<XMLConfig> > GetList() const;
      /** Reload the previous (older) configuration in the list
       * \return: true if the loading was done, false otherwise (no available configuration or last)
       */
      bool Previous();
      /** Reload the next (more recent) configuration in the list
       * \return: true if the loading was done, false otherwise (no available configuration or last)
       */
      bool Next();
   private:
      /// Maximum number of stored configurations
      int mMaxNb;
      /// List of configurations
      std::list<boost::shared_ptr<XMLConfig> > mvpConfig;
      /// The current configuration loaded (can be null if none is selected)
      /// This is reset every time a
      boost::shared_ptr<XMLConfig> mpCurrentConfig;
};

/// Global object to hold configurations history
extern XMLConfigHistory gConfigHistory;

}
#endif //_OBJCRYST_UNDO_H_
