#ifndef Particle_Bank_hh
#define Particle_Bank_hh

#include <stack>

#include "Particle.hh"

class Particle_Bank
{
public:

    Particle_Bank();
    
    int number_of_particles() const
    {
        return number_of_particles_;
    }
    
    void add_particle(Particle &particle)
    {
        particle_bank_.push(particle);
        number_of_particles_ += 1;
    }
    
    Particle get_particle()
    {
        Particle particle(particle_bank_.top());
        particle_bank_.pop();
        number_of_particles_ -= 1;
        
        return particle;
    }
    
    void check_class_invariants() const;
    
private:
    
    int number_of_particles_;
    std::stack<Particle> particle_bank_;
};

#endif
