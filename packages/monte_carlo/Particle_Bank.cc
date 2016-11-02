#include "Particle_Bank.hh"

#include "Check.hh"

Particle_Bank::
Particle_Bank():
    number_of_particles_(0)
{
}

void Particle_Bank::
check_class_invariants()
{
    Assert(particle_bank_.size() == number_of_particles_);
}
