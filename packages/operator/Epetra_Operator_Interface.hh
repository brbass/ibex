#ifndef Epetra_Operator_Interface_hh
#define Epetra_Operator_Interface_hh

#include <memory>

#ifdef EPETRA_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include "Vector_Operator.hh"

using std::shared_ptr;

/*
  Wraps a Vector_Operator to create an Epetra_Operator object
  
  Useful for Krylov solves, in which the matrix does not need to be explicit
*/
class Epetra_Operator_Interface: public Epetra_Operator
{
public:

    // Creator
    Epetra_Operator_Interface(shared_ptr<Epetra_Comm> const &comm,
                              shared_ptr<Epetra_Map> const &map,
                              shared_ptr<Vector_Operator> const &oper);

    // Destructor
    ~Epetra_Operator_Interface();

    // Cannot use transpose
    virtual int SetUseTranspose(bool UseTranspose) override
    {
        return -1;
    }
    
    // Apply the Vector_Operator
    virtual int Apply(Epetra_MultiVector const &X,
                      Epetra_MultiVector &Y) const override;
    
    // Do nothing, as explicit inverse is not available
    virtual int ApplyInverse(Epetra_MultiVector const &X,
                             Epetra_MultiVector &Y) const override
    {
        return 1;
    }
    
    // Cannot provide inf norm
    virtual double NormInf() const override
    {
        return 0.;
    }
    
    // Label for object
    virtual const char *Label() const override
    {
        return "Epetra_Interface";
    }

    // Cannot use transpose
    virtual bool UseTranspose() const override
    {
        return false;
    }

    // Cannot provide inf norm
    virtual bool HasNormInf() const override
    {
        return false;
    }

    // Return associated Epetra_Comm
    virtual const Epetra_Comm &Comm() const override
    {
        return *comm_;
    }

    // Return associated Epetra_Map
    virtual const Epetra_Map &OperatorDomainMap() const override
    {
        return *map_;
    }

    // Return associated Epetra_Map
    virtual const Epetra_Map &OperatorRangeMap() const override
    {
        return *map_;
    }
    

private:
    
    shared_ptr<Epetra_Comm> comm_;
    shared_ptr<Epetra_Map> map_;
    shared_ptr<Vector_Operator> oper_;
};

#endif
