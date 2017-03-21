#ifndef RBF_Factory_hh
#define RBF_Factory_hh

#include <memory>
#include <string>

class RBF;

class RBF_Factory
{
public:

    RBF_Factory();
    
    std::shared_ptr<RBF> get_rbf(std::string rbf_type);
};

#endif
