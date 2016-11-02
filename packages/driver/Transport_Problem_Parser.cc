#include "Transport_Problem_Parser.hh"

using namespace std;

Transport_Problem_Parser::
Transport_Problem_Parser(pugi::xml_node &input_file,
                         shared_ptr<Solver> solver):
    Parser(input_file)
{
    string problem_type_str = XML_Functions::child_value<string>(input_file_, "problem_type");
    
    Transport_Problem::Problem_Type problem_type;
    
    if (problem_type_str == "steady_state")
    {
        problem_type = Transport_Problem::Problem_Type::STEADY_STATE;
    }
    else if (problem_type_str == "k_eigenvalue")
    {
        problem_type = Transport_Problem::Problem_Type::K_EIGENVALUE;
    }
    else if (problem_type_str == "time_dependent")
    {
        problem_type = Transport_Problem::Problem_Type::TIME_DEPENDENT;
    }
    else
    {
        AssertMsg(false, "transport type \"" + problem_type_str +  "\" not found");
    }

    transport_ = make_shared<Transport_Problem>(problem_type,
                                                solver);
}

