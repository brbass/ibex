#include "Manufactured_Cross_Sections.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Cross_Sections::
Manufactured_Cross_Sections(shared_ptr<Angular_Discretization> angular,
                            shared_ptr<Energy_Discretization> energy):
    
    angular_(angular),
    energy_(energy)
{
    Assert(angular_);
    Assert(energy_);
}
