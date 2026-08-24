// Unit tests for RefinablePar/RefinableObj proxy behavior.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "ObjCryst/RefinableObj/RefinableObj.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

class TestRefObj: public ObjCryst::RefinableObj
{
   public:
      TestRefObj(const std::string &className, const std::string &name)
         : mClassName(className)
      {
         this->SetName(name);
      }
      const std::string& GetClassName() const {return mClassName;}
   private:
      std::string mClassName;
};

void TestProxyWithOneParentPrefix()
{
   using namespace ObjCryst;
   REAL x=0.1;
   TestRefObj source("ChildClass","child");
   source.AddPar(new RefinablePar("x",&x,-1.,1.,gpRefParTypeObjCryst));

   TestRefObj parent("ParentClass","parent");
   source.RegisterClient(parent);

   RefinableObj holder;
   holder.AddPar(source,false);

   Check(holder.GetNbPar()==1,"Proxy parameter should be added");
   Check(source.GetPar(0L).GetName()=="x","Original parameter name should stay unchanged");
   Check(source.GetParNameHierarchy(source.GetPar(0L),0)=="x",
         "Hierarchy helper with zero parent should return plain parameter name");
   Check(source.GetParNameHierarchy(source.GetPar(0L),1)=="ParentClass:parent:x",
         "Hierarchy helper should build one-parent prefix");
   Check(holder.GetPar(0L).GetName()=="ParentClass:parent:x",
         "Proxy parameter should use one-parent prefix");
}

void TestProxyWithTwoParentPrefix()
{
   using namespace ObjCryst;
   REAL y=0.2;
   TestRefObj source("ChildClass","child");
   source.AddPar(new RefinablePar("y",&y,-1.,1.,gpRefParTypeObjCryst));

   TestRefObj parent("ParentClass","parent");
   TestRefObj grandParent("GrandParentClass","grand");
   source.RegisterClient(parent);
   parent.RegisterClient(grandParent);

   Check(source.GetParNameHierarchy(source.GetPar(0L),2)=="GrandParentClass:grand:ParentClass:parent:y",
         "Hierarchy helper should build two-parent prefix");

   RefinableObj holder;
   holder.AddPar(source,false);

   Check(holder.GetNbPar()==1,"Proxy parameter should be added");
   Check(holder.GetPar(0L).GetName()=="GrandParentClass:grand:ParentClass:parent:y",
         "Proxy parameter should use two-parent prefix");
}

void TestProxyDelegatesToOriginal()
{
   using namespace ObjCryst;
   REAL v=0.25;
   TestRefObj source("SourceClass","source");
   source.AddPar(new RefinablePar("v",&v,-1.,1.,gpRefParTypeObjCryst));

   RefinableObj holder;
   holder.AddPar(source,false);
   RefinablePar &proxy=holder.GetPar(0L);
   RefinablePar &original=source.GetPar(0L);

   proxy.SetMin(-0.3);
   proxy.SetMax(0.4);
   CheckNearAbs(original.GetMin(), -0.3, 1e-12, "Proxy minimum should update original");
   CheckNearAbs(original.GetMax(), 0.4, 1e-12, "Proxy maximum should update original");

   original.SetIsUsed(false);
   Check(!proxy.IsUsed(), "Original used flag should be visible through proxy");

   proxy.SetValue(0.35);
   CheckNearAbs(original.GetValue(), 0.35, 1e-12, "Proxy value update should mutate original");
   CheckNearAbs(v, 0.35, 1e-12, "Proxy value update should mutate underlying scalar");

   const RefParType delegatedType(gpRefParTypeObjCryst,"DelegatedType");
   proxy.SetType(&delegatedType);
   Check(original.GetType()==&delegatedType, "Proxy type update should propagate to original");
}

void TestProxyDeepHierarchy()
{
   // Verifies that GetParNameHierarchy follows the ancestry chain recursively
   // through three levels of the client registry.
   using namespace ObjCryst;
   REAL z=0.3;
   TestRefObj source("ChildClass","child");
   source.AddPar(new RefinablePar("z",&z,-1.,1.,gpRefParTypeObjCryst));

   TestRefObj parent("ParentClass","parent");
   TestRefObj grandParent("GrandParentClass","grand");
   TestRefObj greatGrandParent("GreatGrandParentClass","great");
   source.RegisterClient(parent);
   parent.RegisterClient(grandParent);
   grandParent.RegisterClient(greatGrandParent);

   Check(source.GetParNameHierarchy(source.GetPar(0L),0)=="z",
         "Depth 0: plain name");
   Check(source.GetParNameHierarchy(source.GetPar(0L),1)=="ParentClass:parent:z",
         "Depth 1: one prefix");
   Check(source.GetParNameHierarchy(source.GetPar(0L),2)=="GrandParentClass:grand:ParentClass:parent:z",
         "Depth 2: two prefixes");
   Check(source.GetParNameHierarchy(source.GetPar(0L),3)=="GreatGrandParentClass:great:GrandParentClass:grand:ParentClass:parent:z",
         "Depth 3: three prefixes (recursive traversal)");
}

void TestProxyCollisionFallback()
{
   // When two source objects produce the same hierarchy name for a parameter,
   // the second proxy should get an address suffix and a warning is printed.
   using namespace ObjCryst;
   REAL x1=0.1, x2=0.2;
   TestRefObj source1("ChildClass","child");
   source1.AddPar(new RefinablePar("x",&x1,-1.,1.,gpRefParTypeObjCryst));

   TestRefObj source2("ChildClass","child");
   source2.AddPar(new RefinablePar("x",&x2,-1.,1.,gpRefParTypeObjCryst));

   TestRefObj parent("ParentClass","parent");
   source1.RegisterClient(parent);
   source2.RegisterClient(parent);

   RefinableObj holder;
   holder.AddPar(source1, false);
   // source2 generates the same hierarchy name "ParentClass:parent:x" -> collision
   holder.AddPar(source2, false);

   Check(holder.GetNbPar()==2, "Both proxy parameters should be added");
   Check(holder.GetPar(0L).GetName()=="ParentClass:parent:x",
         "First proxy should have plain hierarchy name");
   const std::string name2=holder.GetPar(1L).GetName();
   const std::string prefix="ParentClass:parent:x:";
   Check(name2.size()>prefix.size() && name2.substr(0,prefix.size())==prefix,
         "Second proxy should have address suffix appended after hierarchy name");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_refinableobj <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "proxy-one-parent-prefix") TestProxyWithOneParentPrefix();
   else if(testName == "proxy-two-parent-prefix") TestProxyWithTwoParentPrefix();
   else if(testName == "proxy-delegation") TestProxyDelegatesToOriginal();
   else if(testName == "proxy-deep-hierarchy") TestProxyDeepHierarchy();
   else if(testName == "proxy-collision-fallback") TestProxyCollisionFallback();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
