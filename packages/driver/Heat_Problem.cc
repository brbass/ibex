#include "Heat_Problem.hh"

#include <iostream>

#include "Timer.hh"

using namespace std;

Heat_Problem::
Heat_Problem(XML_Node input_node,
                  XML_Node output_node,
                  bool print):
    input_node_(input_node),
    output_node_(output_node),
    print_(print)
{
}

void Heat_Problem::
solve()
{
    Timer timer;
    XML_Node problem_node = input_node_.get_child("problem");
    string type = problem_node.get_attribute<string>("type");
    
    timer.stop();
    times_.emplace_back(timer.time(), "total");
    output_timing();
}

void Heat_Problem::
output_timing()
{
    XML_Node timing_node = output_node_.append_child("timing");
    for (pair<double, string> time : times_)
    {
        timing_node.set_child_value(time.first, time.second);
    }
}

void Heat_Problem::
print_message(string message) const
{
    if (print_)
    {
        cout << endl;
        cout << "Heat Problem:  ";
        cout << message << endl;
    }
}
