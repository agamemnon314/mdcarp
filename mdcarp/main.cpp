#include <filesystem>
#include <fstream>
#include <iostream>
#include "mdcarp.h"

int main()
{
    ListGraph g;
    auto u = g.addNode();
    auto v = g.addNode();
    g.addEdge(u, v);
    ListGraph::EdgeMap<int> edge_cost(g,0);
    ListGraph::EdgeMap<bool> is_required(g, false);
    ListGraph::NodeMap<bool> is_depot(g);

    int n_repeat_times = 100;
    std::vector<int> n_depots = { 2,5,10,20 };
    std::vector<int> n_required_edges_multiplier = { 5,10,20 };


    std::filesystem::create_directory("./output");
    std::ofstream f_out;
    f_out.open("./output/solution_kpi.csv");
    if (!f_out.is_open())
    {
        std::cerr << "cannot open the file";
        return -1;
    }
    f_out << "name,n_depots,n_required_edges,n_nodes,";
    f_out << "lb_opt_1,lb_opt_2,lb_opt,sol_1,sol2,sol,approx,ub_approx" << std::endl;

    for (int nd : n_depots)
    {
        for (int multiplier : n_required_edges_multiplier)
        {
            for (int iter = 0; iter < n_repeat_times; ++iter)
            {
                int n_r = multiplier * nd;
                MDRPP mdrpp(nd, n_r, g, edge_cost, is_required, is_depot);

                std::cout << "Iteration " << iter << " ==============" << std::endl;
                std::cout << "n depots " << nd << ", " << "n required edges " << n_r << std::endl;

                mdrpp.solve();
                std::string instance_name = std::to_string(nd) + "_" + std::to_string(n_r) + "_" + std::to_string(iter);
                mdrpp.write(instance_name);

                f_out << instance_name << ",";
                f_out << mdrpp.n_depots << ",";
                f_out << mdrpp.n_required_edges << ",";
                f_out << countNodes(mdrpp.g) << ",";
                f_out << mdrpp.lb_opt_1 << ",";
                f_out << mdrpp.lb_opt_2 << ",";
                f_out << mdrpp.lb_opt << ",";
                f_out << mdrpp.sol1.first << ",";
                f_out << mdrpp.sol2.first << ",";
                f_out << mdrpp.sol.first << ",";
                f_out << mdrpp.approx_rate << ",";
                f_out << 2.0 - 1.0 / static_cast<double>(2 * mdrpp.n_depots + 1) << "," << std::endl;
                
            }
        }
    }

    f_out.close();
    
    
    return 0;
}