#include "Check.hh"
#include "Epetra_Operator_Interface.hh"

using std::shared_ptr;
using std::vector;

Epetra_Operator_Interface::
Epetra_Operator_Interface(shared_ptr<Epetra_Comm> const &comm,
                          shared_ptr<Epetra_Map> const &map,
                          shared_ptr<Vector_Operator> const &oper):
    comm_(comm),
    map_(map),
    oper_(oper)
{
    Check(oper_->square());
}

Epetra_Operator_Interface::
~Epetra_Operator_Interface()
{
}

int Epetra_Operator_Interface::
Apply(Epetra_MultiVector const &X,
      Epetra_MultiVector &Y) const
{
    Assert(X.NumVectors() == Y.NumVectors());

    int const size = oper_->row_size();
    int const number_of_vectors = X.NumVectors();
    
    vector<double> x(size * number_of_vectors);
    X.ExtractCopy(&x[0], X.Stride());
    
    if (number_of_vectors > 1)
    {
        for (int i = 0; i < number_of_vectors; ++i)
        {
            vector<double> x1(x.begin() + i * size, x.begin() + (i + 1) * size);
            
            (*oper_)(x1);
            
            for (int j = 0; j < size; ++j)
            {
                Y.ReplaceGlobalValue(j, i, x1[j]);
            }
        }
    }
    else
    {
        (*oper_)(x);
        
        for (int i = 0; i < size; ++i)
        {
            Y.ReplaceGlobalValue(i, 0, x[i]);
        }
    }
    
    return 0;
}
