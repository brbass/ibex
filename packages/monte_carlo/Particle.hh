#ifndef Particle_hh
#define Particle_hh

#include <vector>

class Particle
{
public:

    Particle(int dimension,
             int group,
             double weight,
             std::vector<double> position,
             std::vector<double> direction);
    
    // View data
    
    int dimension() const
    {
        return dimension_;
    }
    
    int group() const
    {
        return group_;
    }
    
    double weight() const
    {
        return weight_;
    }
    
    std::vector<double> const &position() const
    {
        return position_;
    }
    
    std::vector<double> const &angle() const
    {
        return angle_;
    }
    
    // Change data
    
    void set_group(int group)
    {
        group_ = group;
    }
    
    void set_weight(double weight)
    {
        weight_ = weight;
    }
    
    void set_position(std::vector<double> &position)
    {
        position_ = position;
    }
    
    void set_angle(std::vector<double> &angle)
    {
        angle_ = angle;
    }

    void check_class_invariants() const;
    
private:
    
    int dimension_;
    int group_;
    double weight_;
    std::vector<double> position_;
    std::vector<double> angle_;
};

#endif
