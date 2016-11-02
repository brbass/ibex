#include "Particle.hh"

#include "Check.hh"

Particle::
Particle(int dimension,
         int group,
         double weight,
         vector<double> position,
         vector<double> direction):
    dimension_(dimension),
    group_(group),
    weight_(weight),
    position_(position),
    direction_(direction)
{
}

void Particle::
check_class_invariants() const
{
    Assert(dimension >= 1);
    Assert(weight > 0);
    Assert(position.size() == dimension);
    Assert(direction.size() == dimension);
}
