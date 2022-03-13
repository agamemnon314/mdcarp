#include "mdcarp.h"
#include <random>
#include <vector>
#include <filesystem>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/lgf_writer.h>

void MDRPP::generate_complete_graph(int n_nodes) const
{
	std::random_device rd;
	std::mt19937 rng(rd()); 
	std::uniform_int_distribution<int> uid(1, 100);

	//add nodes
	std::vector<ListGraph::Node> nodes;
	for (int i = 0; i < n_nodes; ++i)
	{
		auto u = g.addNode();
		is_depot[u] = false;
		nodes.push_back(u);
	}
	std::shuffle(nodes.begin(), nodes.end(), rng);
	for (int i = 0; i < n_depots; ++i)
	{
		is_depot[nodes[i]] = true;
	}

	//add edges
	std::vector<ListGraph::Edge> edges;
	for (int i = 0; i < n_nodes -1; ++i)
	{
		auto u = g.nodeFromId(i);
		for (int j = i+1; j < n_nodes; ++j)
		{
			auto v = g.nodeFromId(j);
			auto e = g.addEdge(u, v);
			edge_cost[e] = uid(rng);
			is_required[e] = false;
			edges.push_back(e);
		}
	}

	std::shuffle(edges.begin(), edges.end(), rng);
	for (int i = 0; i < n_required_edges; ++i)
	{
		is_required[edges[i]] = true;
	}
}

void MDRPP::get_reduced_graph() const
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uid(1, 100);
	int n_nodes = countNodes(g);
	std::vector<ListGraph::Edge> edges_to_remove;
	ListGraph::NodeMap<bool> is_node_required(g, false);

	for (ListGraph::EdgeIt e_it(g); e_it != INVALID; ++e_it)
	{
		if(is_required[e_it])
		{
			auto u = g.u(e_it);
			auto v = g.v(e_it);
			auto u1 = g.addNode();
			is_depot[u1] = false;
			is_node_required[u1] = true;
			auto v1 = g.addNode();
			is_depot[v1] = false;
			is_node_required[v1] = true;
			auto uu1 = g.addEdge(u, u1);
			edge_cost[uu1] = uid(rng);
			auto vv1 = g.addEdge(v, v1);
			edge_cost[vv1] = uid(rng);
			auto u1v1 = g.addEdge(u1, v1);
			edge_cost[u1v1] = edge_cost[e_it];
			is_required[e_it] = false;
			is_required[u1v1] = true;
		}
		else
		{
			edges_to_remove.push_back(e_it);
		}
		
	}

	std::vector<ListGraph::Node> reduced_nodes;
	std::vector<ListGraph::Node> nodes_to_remove;
	for (ListGraph::NodeIt u_it(g); u_it != INVALID; ++u_it)
	{
		if (is_depot[u_it] || is_node_required[u_it])
		{
			reduced_nodes.push_back(u_it);
		}
		else
		{
			nodes_to_remove.push_back(u_it);
		}
	}


	std::map<int, std::map<int, int>> spl_map;
	Dijkstra<ListGraph, ListGraph::EdgeMap<int>> dijkstra(g, edge_cost);
	Dijkstra<ListGraph, ListGraph::EdgeMap<int>>::DistMap dist_map(g);
	dijkstra.distMap(dist_map);
	dijkstra.init();

	for(auto u: reduced_nodes)
	{
		dijkstra.run(u);
		for (auto v : reduced_nodes)
		{
			if (g.id(u) < g.id(v))
			{
				spl_map[g.id(u)][g.id(v)] = dist_map[v];
			}
		}
	}

	for(auto e:edges_to_remove)
	{
		g.erase(e);
	}
	for (auto u : nodes_to_remove)
	{
		g.erase(u);
	}
	for (auto u : reduced_nodes)
	{
		for (auto v : reduced_nodes)
		{
			if (g.id(u) < g.id(v))
			{
				auto e = g.addEdge(u, v);
				edge_cost[e] = spl_map[g.id(u)][g.id(v)];
			}
		}
	}
}

bool MDRPP::check_graph_property() const
{
	for (ListGraph::EdgeIt e_it(g); e_it != INVALID; ++e_it)
	{
		if (is_required[e_it]) 
		{
			if (is_depot[g.u(e_it)] || is_depot[g.v(e_it)])
			{
				std::cout << "required edge adjacent to depot!" << std::endl;
				return false;
			}
		}
	}
	for (ListGraph::NodeIt u_it(g); u_it != INVALID; ++u_it)
	{
		int required_count = 0;
		for (ListGraph::IncEdgeIt e_it(g,u_it); e_it != INVALID; ++e_it)
		{
			if(is_required[e_it])
			{
				required_count++;
				if (required_count>1)
				{
					std::cout << "adjacent required edges!" << std::endl;
					return false;
				}
			}
		}
	}
	return true;
}

bool MDRPP::check_solution_legality(std::vector<std::vector<ListGraph::Edge>>& walks) const
{
	ListGraph::EdgeMap<bool> edge_required_not_in_walks(g, false);
	for (ListGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		if(is_required[e])
		{
			edge_required_not_in_walks[e] = true;
		}
	}

	for (auto& walk : walks)
	{
		if (walk.empty())
		{
			continue;
		}
		bool has_depot = false;
		for (auto e : walk)
		{
			edge_required_not_in_walks[e] = false;
			if (is_depot[g.u(e)] || is_depot[g.v(e)])
			{
				has_depot = true;
			}
		}
		if (!has_depot)
		{
			return false;
		}
	}
	if (countEdges(filterEdges(g, edge_required_not_in_walks)))
	{
		return false;
	}
	return true;
}

ListGraph::Edge MDRPP::find_shortest_path(ListGraph::Edge e1, ListGraph::Edge e2) const
{
	auto u1 = g.u(e1);
	auto v1 = g.v(e1);
	auto u2 = g.u(e2);
	auto v2 = g.v(e2);

	ListGraph::Edge shortest_edge = INVALID;

	for (auto e = findEdge(g, u1, u2); e != INVALID; e = findEdge(g, u1, u2, e))
	{
		if (shortest_edge==INVALID)
		{
			shortest_edge = e;
		}
		else if(e != INVALID && edge_cost[e] < edge_cost[shortest_edge])
		{
			shortest_edge = e;
		}
	}
	

	for (auto e = findEdge(g, v1, v2); e != INVALID; e = findEdge(g, v1, v2, e))
	{
		if (shortest_edge == INVALID)
		{
			shortest_edge = e;
		}
		else if (e != INVALID && edge_cost[e] < edge_cost[shortest_edge])
		{
			shortest_edge = e;
		}
	}

	for (auto e = findEdge(g, u1, v2); e != INVALID; e = findEdge(g, u1, v2, e))
	{
		if (shortest_edge == INVALID)
		{
			shortest_edge = e;
		}
		else if (e != INVALID && edge_cost[e] < edge_cost[shortest_edge])
		{
			shortest_edge = e;
		}
	}

	for (auto e = findEdge(g, v1, u2); e != INVALID; e = findEdge(g, v1, u2, e))
	{
		if (shortest_edge == INVALID)
		{
			shortest_edge = e;
		}
		else if (e != INVALID && edge_cost[e] < edge_cost[shortest_edge])
		{
			shortest_edge = e;
		}
	}
	return shortest_edge;
}

ListGraph::Edge MDRPP::find_shortest_path(ListGraph::Node u, ListGraph::Edge e) const
{
	auto v = g.u(e);
	auto w = g.v(e);

	ListGraph::Edge shortest_edge = INVALID;

	for (auto e1 = findEdge(g, u, v); e1 != INVALID; e1 = findEdge(g, u, v, e1))
	{
		if (shortest_edge == INVALID)
		{
			shortest_edge = e1;
		}
		else if (e1 != INVALID && edge_cost[e1] < edge_cost[shortest_edge])
		{
			shortest_edge = e1;
		}
	}

	for (auto e1 = findEdge(g, u, w); e1 != INVALID; e1 = findEdge(g, u, w, e1))
	{
		if (shortest_edge == INVALID)
		{
			shortest_edge = e1;
		}
		else if (e1 != INVALID && edge_cost[e1] < edge_cost[shortest_edge])
		{
			shortest_edge = e1;
		}
	}
	return shortest_edge;
}

ListGraph::Edge MDRPP::find_orig_edge(ListGraph::Node u, ListGraph::Node v,
	ListGraph::NodeMap<bool>& is_depot_g1,
	ListGraph::NodeMap<ListGraph::Node>& orig_node,
	ListGraph::NodeMap<ListGraph::Edge>& orig_edge) const
{
	ListGraph::Edge e_g = INVALID;
	if (is_depot_g1[u] && is_depot_g1[v])
	{
		for (auto e1 = findEdge(g, orig_node[u], orig_node[v]); e1 != INVALID; e1 = findEdge(g, orig_node[u], orig_node[v], e1))
		{
			if (e_g == INVALID)
			{
				e_g = e1;
			}
			else if (e1 != INVALID && edge_cost[e1] < edge_cost[e_g])
			{
				e_g = e1;
			}
		}
	}
	else if (is_depot_g1[u] && !is_depot_g1[v])
	{
		e_g = find_shortest_path(orig_node[u], orig_edge[v]);
	}
	else if (!is_depot_g1[u] && is_depot_g1[v])
	{
		e_g = find_shortest_path(orig_node[v], orig_edge[u]);
	}
	else
	{
		e_g = find_shortest_path(orig_edge[u], orig_edge[v]);
	}
	return e_g;
}

void MDRPP::export_graph(const std::string& file_name) const
{
	std::ofstream f_out;
	f_out.open(file_name);
	if (f_out.is_open())
	{

		graphWriter(g, f_out).nodeMap("is_depot", is_depot)
			.edgeMap("is_required", is_required)
			.edgeMap("edge_cost", edge_cost)
			.run();
		f_out.close();
	}
	else
	{
		std::cerr << "cannot open the file";
	}
}

void MDRPP::export_sol(const std::string& file_name, std::vector<std::vector<ListGraph::Edge>>& walks) const
{
	std::ofstream f_out;
	f_out.open(file_name);
	if (f_out.is_open())
	{
		for (auto& walk : walks)
		{
			f_out << "walk ===========" << std::endl;
			for (auto e : walk)
			{
				f_out << "(" << g.id(g.u(e)) << ", " << g.id(g.v(e)) << ")"<<std::endl;
			}
			f_out << std::endl;
		}

	}
	else
	{
		std::cerr << "cannot open the file";
	}
}


void MDRPP::solve()
{
	solve_by_mrpp_1();
	if (!check_solution_legality(sol1.second))
	{
		std::cout << "Illegal Solution found !!!!!!!!!!!!!!!!!!!!!==============" << std::endl;
	}
	else
	{
		std::cout << "Solution cost of MRPP1 : " << sol1.first << std::endl;
	}

	solve_by_mrpp_2();
	if (!check_solution_legality(sol2.second))
	{
		std::cout << "Illegal Solution found !!!!!!!!!!!!!!!!!!!!!==============" << std::endl;
	}
	else
	{
		std::cout << "Solution cost of MRPP2 : " << sol2.first << std::endl;
	}

	double alpha_1 = 2 - 1.0 / static_cast<double>(n_depots);
	double beta_1 = 1.0 / static_cast<double>(n_depots);
	double alpha_2 = 3.0;
	double beta_2 = 2.0;

	double x = (alpha_2 * static_cast<double>(sol1.first) - alpha_1 * static_cast<double>(sol2.first)) / (alpha_1 * beta_2 + alpha_2 * beta_1);
	lb_opt_1 = (static_cast<double>(sol1.first) - beta_1*x) / alpha_1;
	lb_opt_2 = (static_cast<double>(sol2.first) + beta_2*x) / alpha_2;
	lb_opt = std::max(lb_opt_1, lb_opt_2);

	sol.first = std::min(sol1.first, sol2.first);
	approx_rate = sol.first / lb_opt;

	std::cout << "Solution cost of MRPP : " << sol.first << std::endl;
	std::cout << "OPT lower bound: " << lb_opt << std::endl;
	std::cout << "Approximation ratio of MRPP : " << approx_rate << std::endl;
	std::cout << "Theoretical Approximation ratio upper bound : " << 2.0 - 1.0 / static_cast<double>(2 * n_depots + 1) << std::endl;
}

void MDRPP::write(const std::string& instance_name)
{
	std::filesystem::create_directory("./output/instance/");
	std::filesystem::create_directory("./output/solution/");

	const std::string graph_file_name = "./output/instance/mdrpp_" + instance_name + ".lgf";
	export_graph(graph_file_name);
	const std::string sol_file_name_1 = "./output/solution/mdrpp_" + instance_name + "_1.sol";
	const std::string sol_file_name_2 = "./output/solution/mdrpp_" + instance_name + "_2.sol";
	export_sol(sol_file_name_1, sol1.second);
	export_sol(sol_file_name_2, sol2.second);
}
