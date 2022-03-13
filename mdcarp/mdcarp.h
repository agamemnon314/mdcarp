#pragma once
#include <lemon/list_graph.h>

using namespace lemon;

class MDRPP
{
public:
	int n_depots;
	int n_required_edges;
	ListGraph& g;
	ListGraph::EdgeMap<int>& edge_cost;
	ListGraph::EdgeMap<bool>& is_required;
	ListGraph::NodeMap<bool>& is_depot;
	std::pair<int, std::vector<std::vector<ListGraph::Edge>>> sol;
	std::pair<int, std::vector<std::vector<ListGraph::Edge>>> sol1;
	std::pair<int, std::vector<std::vector<ListGraph::Edge>>> sol2;
	double approx_rate;
	double approx_rate_improved_1;
	double approx_rate_improved_2;
	double lb_opt;
	double lb_opt_1;
	double lb_opt_1_1;
	double lb_opt_2;

	bool _debug = true;

	MDRPP(int n_depots, int n_required_edges,
		ListGraph& g,
		ListGraph::EdgeMap<int>& edge_cost,
		ListGraph::EdgeMap<bool>& is_required,
		ListGraph::NodeMap<bool>& is_depot) :
		n_depots(n_depots), n_required_edges(n_required_edges),
		g(g),
		edge_cost(edge_cost),
		is_required(is_required),
		is_depot(is_depot)
	{
		g.clear();
		const int graph_size = 2 * (n_depots + n_required_edges);

		generate_complete_graph(graph_size);
		get_reduced_graph();
		if (!check_graph_property())
		{
			std::cout << "Invalid instance!!!!!" << std::endl;
		}
	}

	void generate_complete_graph(int n_nodes) const;
	void get_reduced_graph() const;
	bool check_graph_property() const;
	bool check_solution_legality(std::vector<std::vector<ListGraph::Edge>>& walks) const;


	ListGraph::Edge find_shortest_path(ListGraph::Edge e1, ListGraph::Edge e2) const;
	ListGraph::Edge find_shortest_path(ListGraph::Node u, ListGraph::Edge e) const;
	ListGraph::Edge find_orig_edge(ListGraph::Node u, ListGraph::Node v,
		ListGraph::NodeMap<bool>& is_depot_g1,
		ListGraph::NodeMap<ListGraph::Node>& orig_node,
		ListGraph::NodeMap<ListGraph::Edge>& orig_edge) const;

	void solve_by_mrpp_1();
	void solve_by_mrpp_2();
	void short_cut_walks(std::vector<std::vector<ListGraph::Edge>>& walks);
	void solve();
	void write(const std::string& instance_name);

	void export_graph(const std::string& file_name) const;
	void export_sol(const std::string& file_name, std::vector<std::vector<ListGraph::Edge>>& walks) const;
};
