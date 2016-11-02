#ifndef Vector_Operator_Functions_hh
#define Vector_Operator_Functions_hh

#include <memory>

#include "Vector_Operator.hh"

using std::shared_ptr;

// Add two operators
shared_ptr<Vector_Operator> operator+(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2);

// Subtract two operators
shared_ptr<Vector_Operator> operator-(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2);

// Product of two operators
shared_ptr<Vector_Operator> operator*(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2);

#endif
