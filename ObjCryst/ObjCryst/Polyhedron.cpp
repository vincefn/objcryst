#include "ObjCryst/Polyhedron.h"

namespace ObjCryst
{

Molecule* MakeTetrahedron(Crystal &cryst,const string &name,
                          const ScatteringPower *centralAtom,
                          const ScatteringPower *peripheralAtom,
                          const REAL dist)
{
   Molecule *mol=new Molecule(cryst,name);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom(0.,0.,dist,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(0.,dist*0.943,-dist*0.333,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(dist*0.817,-dist*0.472,-dist*0.333,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(-dist*0.817,-dist*0.472,-dist*0.333,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),dist,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),109.5*DEG2RAD,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(3),109.5*DEG2RAD,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(4),109.5*DEG2RAD,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),109.5*DEG2RAD,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(4),109.5*DEG2RAD,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(4),109.5*DEG2RAD,0.01,0.05);
   return mol;
}

Molecule* MakeOctahedron(Crystal &cryst,const string &name,
                         const ScatteringPower *centralAtom,
                         const ScatteringPower *peripheralAtom,
                         const REAL dist)
{
   Molecule *mol=new Molecule(cryst,name);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom(0.,0.,dist,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(dist,0.,0.,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(0.,dist,0.,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(-dist,0.,0.,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddAtom(0.,-dist,0.,peripheralAtom,peripheralAtom->GetName()+"5");
   mol->AddAtom(0.,0.,-dist,peripheralAtom,peripheralAtom->GetName()+"6");
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(5),dist,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(6),dist,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(3),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(4),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(5),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(2),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(3),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(4),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(5),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(4),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(5),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(5),mol->GetAtom(0),mol->GetAtom(2),M_PI/2.,0.01,0.05);
   return mol;
}

Molecule* MakeSquarePlane(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL d)
{
   Molecule *mol=new Molecule(cryst,name);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom( d,0.,0.,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(0., d,0.,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(-d,0.,0.,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(0.,-d,0.,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),d,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(4),M_PI/2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(1),M_PI/2.,0.01,0.05);
   return mol;
}

Molecule* MakeCube(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL d)
{
   Molecule *mol=new Molecule(cryst,name);
   const REAL d0=d/sqrt(3.);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom( d0,-d0, d0,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom( d0, d0, d0,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(-d0, d0, d0,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(-d0,-d0, d0,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddAtom( d0,-d0,-d0,peripheralAtom,peripheralAtom->GetName()+"5");
   mol->AddAtom( d0, d0,-d0,peripheralAtom,peripheralAtom->GetName()+"6");
   mol->AddAtom(-d0, d0,-d0,peripheralAtom,peripheralAtom->GetName()+"7");
   mol->AddAtom(-d0,-d0,-d0,peripheralAtom,peripheralAtom->GetName()+"8");
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(5),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(6),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(7),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(8),d,0.01,0.05);
   
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(2),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(3),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(4),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(1),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(6),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(6),mol->GetAtom(7),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(8),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(5),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(5),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(6),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(7),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(8),2*d0,0.01,0.05);
   #if 0
   const REAL a=2*atan(1/sqrt(2.));
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(4),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(1),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(5),mol->GetAtom(0),mol->GetAtom(6),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(7),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(7),mol->GetAtom(0),mol->GetAtom(8),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(8),mol->GetAtom(0),mol->GetAtom(1),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(5),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(6),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(7),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(8),a,0.01,0.05);
   #endif
   return mol;
}

Molecule* MakeAntiPrismTetragonal(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL d)
{
   Molecule *mol=new Molecule(cryst,name);
   const REAL d0=d/sqrt(3.);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom( d0,-d0, d0,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom( d0, d0, d0,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(-d0, d0, d0,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(-d0,-d0, d0,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddAtom( d0*sqrt(2.),0.,-d0,peripheralAtom,peripheralAtom->GetName()+"5");
   mol->AddAtom(0., d0*sqrt(2.),-d0,peripheralAtom,peripheralAtom->GetName()+"6");
   mol->AddAtom(-d0*sqrt(2.),0.,-d0,peripheralAtom,peripheralAtom->GetName()+"7");
   mol->AddAtom(0.,-d0*sqrt(2.),-d0,peripheralAtom,peripheralAtom->GetName()+"8");
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(5),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(6),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(7),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(8),d,0.01,0.05);
   
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(2),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(3),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(4),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(1),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(6),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(6),mol->GetAtom(7),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(8),2*d0,0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(5),2*d0,0.01,0.05);
   
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(1),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(2),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(6),mol->GetAtom(2),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(6),mol->GetAtom(3),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(3),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(4),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(4),d0*sqrt(5.),0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(1),d0*sqrt(5.),0.01,0.05);
   #if 0
   const REAL a=2*atan(1/sqrt(2.));
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(4),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(1),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(5),mol->GetAtom(0),mol->GetAtom(6),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(6),mol->GetAtom(0),mol->GetAtom(7),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(7),mol->GetAtom(0),mol->GetAtom(8),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(8),mol->GetAtom(0),mol->GetAtom(1),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(5),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(6),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(3),mol->GetAtom(0),mol->GetAtom(7),a,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(4),mol->GetAtom(0),mol->GetAtom(8),a,0.01,0.05);
   #endif
   return mol;
}

Molecule* MakePrismTrigonal(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL d)
{
   const REAL a=sqrt(3./7.)*d;
   const REAL a1=a*2./sqrt(3.);
   Molecule *mol=new Molecule(cryst,name);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom( a1   ,0.             , a,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(-a1/2., a1*sqrt(3.)/2., a,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(-a1/2.,-a1*sqrt(3.)/2., a,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom( a1   ,0.             ,-a,peripheralAtom,peripheralAtom->GetName()+"4");
   mol->AddAtom(-a1/2., a1*sqrt(3.)/2.,-a,peripheralAtom,peripheralAtom->GetName()+"5");
   mol->AddAtom(-a1/2.,-a1*sqrt(3.)/2.,-a,peripheralAtom,peripheralAtom->GetName()+"6");

   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(5),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(6),d,0.01,0.05);

   mol->AddBond(mol->GetAtom(1),mol->GetAtom(2),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(3),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(3),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(6),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(6),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(4),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(6),2*a,0.01,0.05);
   return mol;
}

Molecule* MakeIcosahedron(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL d)
{
   Molecule *mol=new Molecule(cryst,name);
   
   const REAL g0=(1.+sqrt(5.))/2.;
   const REAL a=d/sqrt(1.+g0*g0);
   const REAL g=g0*a;
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());

   mol->AddAtom(0., g, a,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(0., g,-a,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(0.,-g, a,peripheralAtom,peripheralAtom->GetName()+"3");
   mol->AddAtom(0.,-g,-a,peripheralAtom,peripheralAtom->GetName()+"4");

   mol->AddAtom( a,0., g,peripheralAtom,peripheralAtom->GetName()+"5");
   mol->AddAtom(-a,0., g,peripheralAtom,peripheralAtom->GetName()+"6");
   mol->AddAtom( a,0.,-g,peripheralAtom,peripheralAtom->GetName()+"7");
   mol->AddAtom(-a,0.,-g,peripheralAtom,peripheralAtom->GetName()+"8");
   
   mol->AddAtom( g, a,0.,peripheralAtom,peripheralAtom->GetName()+"9");
   mol->AddAtom( g,-a,0.,peripheralAtom,peripheralAtom->GetName()+"10");
   mol->AddAtom(-g, a,0.,peripheralAtom,peripheralAtom->GetName()+"11");
   mol->AddAtom(-g,-a,0.,peripheralAtom,peripheralAtom->GetName()+"12");

   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(4),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(5),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(6),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(7),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(8),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(9),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(10),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(11),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(12),d,0.01,0.05);

   mol->AddBond(mol->GetAtom(1),mol->GetAtom(2),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(11),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(6),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(1),mol->GetAtom(9),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(2),mol->GetAtom(11),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(11),mol->GetAtom(6),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(6),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(5),mol->GetAtom(9),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(9),mol->GetAtom(2),2*a,0.01,0.05);

   mol->AddBond(mol->GetAtom(4),mol->GetAtom(3),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(12),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(8),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(7),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(4),mol->GetAtom(10),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(12),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(12),mol->GetAtom(8),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(7),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(10),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(10),mol->GetAtom(3),2*a,0.01,0.05);

   mol->AddBond(mol->GetAtom(3),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(6),2*a,0.01,0.05);
   
   mol->AddBond(mol->GetAtom(12),mol->GetAtom(6),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(12),mol->GetAtom(11),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(11),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(8),mol->GetAtom(2),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(2),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(7),mol->GetAtom(9),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(10),mol->GetAtom(9),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(10),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(5),2*a,0.01,0.05);
   mol->AddBond(mol->GetAtom(3),mol->GetAtom(6),2*a,0.01,0.05);
   return mol;
}

Molecule* MakeTriangle(Crystal &cryst,const string &name,
                       const ScatteringPower *centralAtom,
                       const ScatteringPower *peripheralAtom,
                       const REAL d)
{
   Molecule *mol=new Molecule(cryst,name);
   mol->AddAtom(0.,0.,0.,centralAtom,centralAtom->GetName());
   mol->AddAtom(d,0.,0.,peripheralAtom,peripheralAtom->GetName()+"1");
   mol->AddAtom(-d/2., d*sqrt(3.)/2,0.,peripheralAtom,peripheralAtom->GetName()+"2");
   mol->AddAtom(-d/2.,-d*sqrt(3.)/2,0.,peripheralAtom,peripheralAtom->GetName()+"3");

   mol->AddBond(mol->GetAtom(0),mol->GetAtom(1),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(2),d,0.01,0.05);
   mol->AddBond(mol->GetAtom(0),mol->GetAtom(3),d,0.01,0.05);

   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(2),M_PI/3.*2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(1),mol->GetAtom(0),mol->GetAtom(3),M_PI/3.*2.,0.01,0.05);
   mol->AddBondAngle(mol->GetAtom(2),mol->GetAtom(0),mol->GetAtom(3),M_PI/3.*2.,0.01,0.05);
   return mol;
}

}//namespace
