#ifndef _OBJCRYST_CIF_H
#define _OBJCRYST_CIF_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>

#include "ObjCryst/Quirks/ci_string.h"
namespace ObjCryst
{// Forward declaration
   class CIF;
}
#include "ObjCryst/ObjCryst/PowderPattern.h" // For CreatePowderPatternFromCIF only.
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h" // For CreateSingleCrystalDataFromCIF only.
#include "ObjCryst/ObjCryst/Crystal.h" // For CreateCrystalFromCIF only.
#include "ObjCryst/ObjCryst/General.h" // TO identify wavelength type in CIFData::ExtractPowderPattern.

namespace ObjCryst
{
/// Convert one CIF value to a floating-point value
/// Return 0 if no value can be converted (e.g. if '.' or '?' is encountered)
REAL CIFNumeric2REAL(const std::string &s);
/// Convert one CIF value to a floating-point value
/// Return 0 if no value can be converted (e.g. if '.' or '?' is encountered)
int CIFNumeric2Int(const std::string &s);

/** The CIFData class holds all the information from a \e single data_ block from a cif file.
* 
* It is a placeholder for all comments, item and loop data, as raw strings copied from
* a cif file.
*
* It is also used to interpret this data to extract parts of the cif data, i.e.
* only part of the core cif dictionnary are recognized. CIF tags currently recognized
* include ("tag1 > tag2" means tag1 is preferred to tag2 when extracting the info, only one is reported):
*  - crystal name: _chemical_name_systematic > _chemical_name_mineral > _chemical_name_structure_type > _chemical_name_common
*  - crystal formula: _chemical_formula_analytical > _chemical_formula_structural > _chemical_formula_iupac > _chemical_formula_moiety
*  - unit cell:  _cell_length_{a,b,c} ; _cell_angle_{alpha,beta,gamma}
*  - spacegroup number: _space_group_IT_number > _symmetry_Int_Tables_number
*  - spacegroup Hall symbol: _space_group_name_Hall > _symmetry_space_group_name_Hall
*  - spacegroup Hermann-Mauguin symbol:_space_group_name_H-M_alt > _symmetry_space_group_name_H-M
*  - atom coordinates: _atom_site_fract_{x} ; _atom_site_Cartn_{x,y,z}
*  - atom occupancy: _atom_site_occupancy
*  - atom label & symbol: _atom_site_type_symbol ; _atom_site_label
*  - atom adps: _atom_site_aniso_{U,B}_{11,22,33,12,13,23} > _atom_site_{U,B}_iso_or_equiv
*
*
* Cartesian coordinates are stored in Angstroems, angles in radians.
*
* To import PowderPattern data, the following tags are used:
* - observed intensity: _pd_meas_counts_total > _pd_meas_intensity_total > _pd_proc_intensity_total > _pd_proc_intensity_net
* - uncertainty on intensity: deducted from _pd_proc_ls_weight or square root of intensity
* - coordinates: _pd_proc_2theta_corrected > _pd_meas_angle_2theta > _pd_meas_time_of_flight > _pd_proc_2theta_range_{min,max,inc}
* - intensity normalizer (optional): _pd_meas_intensity_monitor > _pd_meas_step_count_time
*
* If another data field is needed, it is possible to directly access the string data 
* (CIFData::mvComment , CIFData::mvItem and CIFData::mvLoop) to search for the correct tags.
*/
class CIFData
{
   public:
      CIFData();
      
      /// Extract lattice parameters, spacegroup (symbol or number), atomic positions,
      /// chemical name and formula if available.
      /// All other data is ignored
      void ExtractAll(const bool verbose=false);
      /// Extract name & formula for the crystal
      void ExtractName(const bool verbose=false);
      /// Extract unit cell
      void ExtractUnitCell(const bool verbose=false);
      /// Extract spacegroup number or symbol
      void ExtractSpacegroup(const bool verbose=false);
      /// Extract all atomic positions. Will generate cartesian from fractional
      /// coordinates or vice-versa if only cartesian coordinates are available.
      void ExtractAtomicPositions(const bool verbose=false);
      /// Extract anisotropic atomic displacement parameters. Isotropic ADPs
      //  may be extracted in ExtractAtomicPositions. This takes prescedent.
      void ExtractAnisotropicADPs(const bool verbose=false);
      /// Extract Powder Diffraction data, with Iobs, sigma(Iobs) and either 2theta
      /// or time-of-flight position.
      void ExtractPowderPattern(const bool verbose=false);
      /// Extract single crystal data, with Iobs, sigma(Iobs) and h,k,l
      void ExtractSingleCrystalData(const bool verbose=false);
      /// Generate fractional coordinates from cartesian ones for all atoms
      /// CIFData::CalcMatrices() must be called first
      void Cartesian2FractionalCoord();
      /// Generate cartesian coordinates from fractional ones for all atoms
      /// CIFData::CalcMatrices() must be called first
      void Fractional2CartesianCoord();
      /// Convert from fractional to cartesian coordinates
      /// CIFData::CalcMatrices() must be called first
      void f2c(REAL &x,REAL &y, REAL &z);
      /// Convert from cartesia to fractional coordinates
      /// CIFData::CalcMatrices() must be called first
      void c2f(REAL &x,REAL &y, REAL &z);
      /// Calculate real space transformation matrices
      /// requires unit cell parameters
      void CalcMatrices(const bool verbose=false);
      /// Comments from CIF file, in the order they were read
      std::list<std::string> mvComment;
      /// Individual CIF items
      std::map<ci_string,std::string> mvItem;
      /// CIF Loop data
      std::map<std::set<ci_string>,std::map<ci_string,std::vector<std::string> > > mvLoop;
      /// Lattice parameters, in ansgtroem and degrees - vector size is 0 if no
      /// parameters have been obtained yet.
      std::vector<REAL> mvLatticePar;
      /// Spacegroup number from International Tables (_space_group_IT_number), or -1.
      std::string mSpacegroupNumberIT;
      /// Spacegroup Hall symbol (or empty string) (_space_group_name_Hall)
      std::string mSpacegroupSymbolHall;
      /// Spacegroup Hermann-Mauguin symbol (or empty string) (_space_group_name_H-M_alt)
      std::string mSpacegroupHermannMauguin;
      /// Map of _symmetry_equiv_pos_as_xyz strings
      std::set<string> mvSymmetry_equiv_pos_as_xyz;
      /// Crystal name. Or empty string if none is available.
      std::string mName;
      /// Formula. Or empty string if none is available.
      std::string mFormula;
      /// Atom record 
      struct CIFAtom
      {
         CIFAtom();
         /// Label of the atom, or empty string (_atom_site_label).
         std::string mLabel;
         /// Symbol of the atom, or empty string (_atom_type_symbol or _atom_site_type_symbol).
         std::string mSymbol;
         /// Fractionnal coordinates (_atom_site_fract_{x,y,z}) or empty vector.
         std::vector<REAL> mCoordFrac;
         /// Cartesian coordinates in Angstroem (_atom_site_Cartn_{x,y,z}) or empty vector.
         /// Transformation to fractionnal coordinates currently assumes 
         /// "a parallel to x; b in the plane of y and z" (see _atom_sites_Cartn_transform_axes)
         std::vector<REAL> mCoordCart;
         /// Site occupancy, or -1
         REAL mOccupancy;

         /// ADP tensor
         std::vector<REAL> mBeta;

         /// Biso
         REAL mBiso;
      };
      /// Atoms, if any are found
      std::vector<CIFAtom> mvAtom;
      /// Fractionnal2Cartesian matrix
      REAL mOrthMatrix[3][3];
      /// Cartesian2Fractionnal matrix
      REAL mOrthMatrixInvert[3][3];
      /// Powder pattern data
      std::vector<REAL> mPowderPatternObs,mPowderPatternX,mPowderPatternSigma;
      /// Single crystal data
      CrystVector_long mH,mK,mL;
      /// Single crystal data
      CrystVector_REAL mIobs,mSigma;
      /// Is this X-Ray 2theta, time-of-flight ?
      WavelengthType mDataType;
      /// Wavelength
      REAL mWavelength;
};

/** Main CIF class - parses the stream and separates data blocks, comments, items, loops.
* All values are stored as string, and Each CIF block is stored in a separate CIFData object.
* No interpretaion is made here - this must be done from all CIFData objects.
*/
class CIF
{
   public:
      /// Creates the CIF object from a stream
      ///
      /// \param interpret: if true, interpret all data blocks. See CIFData::ExtractAll()
      CIF(std::istream &in, const bool interpret=true,const bool verbose=false);
   //private:
      /// Separate the file in data blocks and parse them to sort tags, loops and comments.
      /// All is stored in the original strings.
      void Parse(std::stringstream &in);
      /// The data blocks, after parsing. The key is the name of the data block
      std::map<std::string,CIFData> mvData;
      /// Global comments, outside and data block
      std::list<std::string> mvComment;
};

// Forward declarations
class Crystal;
class PowderPattern;

/** Extract Crystal object(s) from a CIF, if possible.
* Returns a null pointer if no crystal structure could be extracted
* (the minimum data is the unit cell parameters).
*
* \param checkSymAsXYZ: if true, and the CIF file does not have a Hall symbol
* but has a list of symmetry_equiv_pos_as_xyz, check we have the correct
* setting by trying different ones using cctbx
*/
Crystal* CreateCrystalFromCIF(CIF &cif,const bool verbose=true,const bool checkSymAsXYZ=true);

/// Create PowderPattern object(s) from a CIF, if possible.
/// Returns a null pointer if no pattern could be extracted.
/// No components (background, crystal data) are created.
PowderPattern* CreatePowderPatternFromCIF(CIF &cif);

/// Create DiffractionDataSingleCrystal object(s) from a CIF, if possible.
/// Returns a null pointer if no data could be extracted.
/// A Crystal object must be supplied - if none is given, the last Crystal
/// object will be used. If no Cyrstal data exists, a new one will be created.
DiffractionDataSingleCrystal* CreateSingleCrystalDataFromCIF(CIF &cif, Crystal *pryst=0);

}

#endif
