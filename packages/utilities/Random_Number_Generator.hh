#ifndef Random_Number_Generator_hh
#define Random_Number_Generator_hh

#include <chrono>
#include <random>
#include <vector>

/*
  Simple class to generate several types or sequences of random numbers
*/
template<class T>
class Random_Number_Generator
{
public:

    // Constructor
    Random_Number_Generator(T lower_bound,
                            T upper_bound);

    // Constructor
    Random_Number_Generator(T lower_bound,
                            T upper_bound,
                            int seed);
    
    
    // Returns a random scalar
    T scalar();
    
    // Returns a sequence of random scalars
    std::vector<T> vector(int number_of_elements);

    void seed(T seed);
    
private:
    
    // Get seed for random number generator
    // Uses time since epoch
    T get_seed();
};

template<>
class Random_Number_Generator<int>
{
public:
    
    Random_Number_Generator(int lower_bound,
                            int upper_bound):
        generator(get_seed()),
        distribution(lower_bound, upper_bound)
    {
    }

    Random_Number_Generator(int lower_bound,
                            int upper_bound,
                            int seed):
        generator(seed),
        distribution(lower_bound, upper_bound)
    {
    }
    
    int scalar()
    {
        return distribution(generator);
    }

    std::vector<int> vector(int number_of_elements)
    {
        std::vector<int> vec(number_of_elements);

        for (unsigned i = 0; i < number_of_elements; ++i)
        {
            vec[i] = distribution(generator);
        }

        return vec;
    }
    
    void seed(int seed)
    {
        generator.seed(seed);
    }
    
private:

    int get_seed()
    {
        using namespace std::chrono;
        
        duration<int> time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch());
        
        return time.count();
    }
    
    std::default_random_engine generator;
    std:: uniform_int_distribution<int> distribution;
};

template<>
class Random_Number_Generator<double>
{
public:
    
    Random_Number_Generator(double lower_bound,
                            double upper_bound):
        generator(get_seed()),
        distribution(lower_bound, upper_bound)
    {
    }

    Random_Number_Generator(double lower_bound,
                            double upper_bound,
                            int seed):
        generator(seed),
        distribution(lower_bound, upper_bound)
    {
    }
    double scalar()
    {
        return distribution(generator);
    }

    std::vector<double> vector(int number_of_elements)
    {
        std::vector<double> vec(number_of_elements);

        for (unsigned i=0; i<number_of_elements; ++i)
        {
            vec[i] = distribution(generator);
        }

        return vec;
    }
    
    void seed(double seed)
    {
        generator.seed(seed);
    }
    
private:

    double get_seed()
    {
        using namespace std::chrono;
        
        duration<double> time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch());
        
        return time.count();
    }
    
    std::default_random_engine generator;
    std:: uniform_real_distribution<double> distribution;
};


#endif
