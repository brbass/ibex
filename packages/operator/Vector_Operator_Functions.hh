#ifndef Vector_Operator_Functions_hh
#define Vector_Operator_Functions_hh

#include <memory>

#include "Vector_Operator.hh"

// Add two operators
std::shared_ptr<Vector_Operator> operator+(std::shared_ptr<Vector_Operator> const op1,
                                           std::shared_ptr<Vector_Operator> const op2);

// Subtract two operators
std::shared_ptr<Vector_Operator> operator-(std::shared_ptr<Vector_Operator> const op1,
                                           std::shared_ptr<Vector_Operator> const op2);

// Product of two operators
std::shared_ptr<Vector_Operator> operator*(std::shared_ptr<Vector_Operator> const op1,
                                           std::shared_ptr<Vector_Operator> const op2);

#endif
