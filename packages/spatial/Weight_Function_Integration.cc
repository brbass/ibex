#include "Weight_Function_Integration.hh"

using std::shared_ptr;
using std::vector;

Weight_Function_Integration::
Weight_Function_Integration(Weight_Function::Options options,
                            vector<shared_ptr<Basis_Function> > const &bases,
                            vector<shared_ptr<Weight_Function> > const &weights,
                            shared_ptr<Solid_Geometry> solid,
                            vector<vector<double> > limits,
                            vector<int> dimensional_cells):
    options_(options),
    bases_(bases),
    weights_(weights),
    solid_(solid)
{
    mesh_ = make_shared<Mesh>(*this,
                              solid->dimension(),
                              limits,
                              dimensional_cells);
}

Weight_Function_Integration::Mesh
Mesh(Weight_Function_Integration const &wfi,
     int dimension,
     vector<vector<double> > limits,
     vector<int> dimensional_cells):
    wfi_(wfi),
    dimension_(dimension),
    limits_(limits),
    dimensional_cells_(dimensional_cells)
{
    initialize();
    kd_tree_ = get_kd_tree();
}

void Weight_Function_Integration::Mesh
initialize()
{
    vector<vector<double> > positions;
    
    // Check sizes
    Assert(dimensional_cells_.size() == dimension);
    Assert(limits_.size() == dimension);
    
    // Get total number of nodes
    dimensional_nodes_.resize(dimension_);
    number_of_global_nodes_ = 1;
    number_of_global_cells_ = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        Assert(dimensional_cells_[d] >= 1);
        dimensional_nodes_[d] = dimensional_cells_[d] + 1;
        number_of_global_nodes_ *= dimensional_nodes_[d];
        number_of_global_cells_ *= dimensional_cells_[d];
    }
    
    // Get intervals between cells
    intervals_.resize(dimension);
    for (int d = 0; d < dimension; ++d)
    {
        intervals_[d] = (limits_[d][1] - limits_[d][0]) / static_cast<double>(dimensional_cells_[d]);
    }
    

    // Initialize nodes
    nodes_.resize(number_of_global_nodes_);
    switch (dimension)
    {
    case 1:
    {
        int d_i = 0;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            double index = i;
            Node &node = nodes_[i];
            node.position = {limits_[d_i][0] + intervals_[d_i] * i};
            
            if (i != 0)
            {
                node.neighboring_cells.push_back(i - 1);
            }
            if (i != dimensional_nodes_[d] - 1)
            {
                node.neighboring_cells.push_back(i);
            }
        }
        break;
    }
    case 2:
    {
        int d_i = 0;
        int d_j = 0;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++i)
            {
                int index = j + dimensional_nodes_[1] * i;
                Node &node = nodes_[index];

                // Get node position
                node.position = {limits_[d_i][0] + intervals_[d_i] * i,
                                 limits_[d_j][0] + intervals_[d_j] * j};

                // Get neighboring cells
                for (int ci = i - 1; ci <= i; ++ci)
                {
                    if (ci < 0 || ci >= dimensional_cells_[0])
                    {
                        continue;
                    }
                    for (int cj = i - 1; cj <= i; ++cj)
                    {
                        if (cj < 0 || cj >= dimensional_cells_[1])
                        {
                            continue;
                        }
                        int c_index = cj + dimensional_cells_[1] * ci;
                        node.neighboring_cells.push_back(c_index);
                    }
                }
            }
        }
        break;
    }
    case 3:
    {
        int d_i = 0;
        int d_j = 0;
        int d_k = 0;
        for (int i = 0; i < dimensional_nodes_[0]; ++i)
        {
            for (int j = 0; j < dimensional_nodes_[1]; ++i)
            {
                for (int k = 0; k < dimensional_nodes_[2]; ++k)
                {
                    int index = k + dimensional_nodes_[2] * (j + dimensional_nodes_[1] * i);
                    Node &node = nodes_[index];

                    // Get node position
                    node.position = {limits_[d_i][0] + intervals_[d_i] * i,
                                     limits_[d_j][0] + intervals_[d_j] * j,
                                     limits_[d_k][0] + intervals_[d_k] * k};
                    
                    // Get neighboring cells
                    for (int ci = i - 1; ci <= i; ++ci)
                    {
                        if (ci < 0 || ci >= dimensional_cells_[0])
                        {
                            continue;
                        }
                        for (int cj = i - 1; cj <= i; ++cj)
                        {
                            if (cj < 0 || cj >= dimensional_cells_[1])
                            {
                                continue;
                            }
                            for (int ck = i - 1; ck <= i; ++ck)
                            {
                                if (ck < 0 || ck >= dimensional_cells_[2])
                                {
                                    continue;
                                }
                                
                                int c_index = ck + dimensional_cells_[2] * (cj + dimensional_cells_[1] * ci);
                                node.neighboring_cells.push_back(c_index);
                            }
                        }
                    }
                }
            }
        }
        break;
    }
    default:
        AssertMsg(false, "dimension (" << dimension_ << ") not found");
    }

    // Get cells
    switch (dimension_)
    {
    case 1:
    {
        break;
    }
    case 2:
    {
        for (int i = 0; i < dimensional_cells_[0]; ++i)
        {
            for (int j = 0; j < dimensional_cells_[1]; ++j)
            {
                int index = j + dimensional_cells_[1] * i;
                vector<int> indices = {i, j};
                Cell &cell = cells[index];
                
                // Set neighboring nodes and cells
                for (int ni = i; ni <= i + 1; ++ni)
                {
                    for (int nj = j; nj <= j + 1; ++nj)
                    {
                        int n_index = nj + dimensional_nodes_[1] * ni;
                        Node &node = nodes_[n_index];
                        node.neighboring_cells.push_back(index);
                        cell.neighboring_nodes.push_back(n_index);
                    }
                }
                
                // Set upper and lower limits
                cell.limits.assign(dimension_, vector<double>(2));
                for (int d = 0; d < dimension_; ++d)
                {
                    cell.limits[d][0] = indices[d]
                }
            }
        }
        break;
    }
    case 3:
    {
        break;
    }
    default:
        AssertMsg(false, "dimension (" << dimension_ << ") not found");
    }
    
    // Get KD tree
    vector<vector<double> > kd_positions(number_of_global_nodes_, vector<double>(dimension_));
    for (int i = 0; i < number_of_global_nodes_; ++i)
    {
        kd_positions[i] = nodes_[i].position;
    }
    kd_tree_ = make_shared<KD_Tree>(dimension_,
                                    number_of_global_nodes_,
                                    kd_positions);
    
}

