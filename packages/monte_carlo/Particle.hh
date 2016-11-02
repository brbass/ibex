#ifndef Particle_hh
#define Particle_hh

#include <vector>

using std::vector;

class Particle
{
public:

    Particle(int dimension,
             int group,
             double weight,
             vector<double> position,
             vector<double> direction);
    
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
    
    vector<double> const &position() const
    {
        return position_;
    }
    
    vector<double> const &angle() const
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
    
    void set_position(vector<double> &position)
    {
        position_ = position;
    }
    
    void set_angle(vector<double> &angle)
    {
        angle_ = angle;
    }

    void check_class_invariants() const;
    
private:
    
    int dimension_;
    int group_;
    double weight_;
    vector<double> position_;
    vector<double> angle_;
};

#endif
