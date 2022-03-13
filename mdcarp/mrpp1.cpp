#include "mdcarp.h"
#include <lemon/kruskal.h>
#include <lemon/matching.h>
#include <lemon/euler.h>
#include <lemon/adaptors.h>

void MDRPP::solve_by_mrpp_1()
{
	ListGraph g1;
	auto u = g1.addNode();
	auto v = g1.addNode();
	g1.addEdge(u, v);
	ListGraph::EdgeMap<int> edge_cost_for_matching(g1);
	ListGraph::EdgeMap<int> edge_cost_for_forest(g1);
	ListGraph::NodeMap<bool> is_depot_g1(g1, false);
	ListGraph::NodeMap<ListGraph::Node> orig_node(g1, INVALID);
	ListGraph::NodeMap<ListGraph::Edge> orig_edge(g1, INVALID);
	g1.clear();

	// add node to g1 from depots of g
	for (ListGraph::NodeIt u(g); u != INVALID; ++u)
	{
		if (is_depot[u])
		{
			auto u1 = g1.addNode();
			is_depot_g1[u1] = true;
			orig_node[u1] = u;
		}
	}

	// add node to g1 from required edges of g
	for (ListGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		if (is_required[e])
		{
			auto u = g1.addNode();
			is_depot_g1[u] = false;
			orig_edge[u] = e;
		}
	}

	// add edges to g1
	for (ListGraph::NodeIt u(g1); u != INVALID; ++u)
	{
		for (ListGraph::NodeIt v(g1); v != INVALID; ++v)
		{
			if (g1.id(u) < g1.id(v))
			{
				auto e = g1.addEdge(u, v);
				ListGraph::Edge orig_e = find_orig_edge(u, v, is_depot_g1, orig_node, orig_edge);
				edge_cost_for_matching[e] = -edge_cost[orig_e];
				edge_cost_for_forest[e] = -edge_cost_for_matching[e];

				if (is_depot_g1[u] && is_depot_g1[v])
				{
					edge_cost_for_forest[e] = 0;
				}
			}
		}
	}

	// find constrained forest F
	std::vector<ListGraph::Edge> tree_edges(countNodes(g1) - 1);
	kruskal(g1, edge_cost_for_forest, tree_edges.begin());

	ListGraph::EdgeMap<bool> forest_edge_map(g1, false);
	for (auto e : tree_edges)
	{
		if (is_depot_g1[g1.u(e)] && is_depot_g1[g1.v(e)])
		{
			continue;
		}
		forest_edge_map[e] = true;
	}

	FilterEdges<ListGraph> sub_graph_f(g1, forest_edge_map);

	ListGraph::NodeMap<bool> odd_degree_node_filter(g1, false);
	for (FilterEdges<ListGraph>::NodeIt u(sub_graph_f); u != INVALID; ++u)
	{
		int degree = 0;
		for (FilterEdges<ListGraph>::IncEdgeIt e(sub_graph_f, u); e != INVALID; ++e)
		{
			degree++;
		}
		if (degree % 2 == 1)
		{
			odd_degree_node_filter[u] = true;
		}
	}

	FilterNodes<ListGraph> sub_graph_odd_degree(g1, odd_degree_node_filter);

	MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<int>> max_weighted_perfect_matching(sub_graph_odd_degree, edge_cost_for_matching);
	bool result = max_weighted_perfect_matching.run();
	if (!result)
	{
		std::cout << "fatal error!!!!!!!!!!!!!" << std::endl;
	}

	ListGraph solution_graph;
	ListGraph::EdgeMap<bool> disable_all_edges(g, false);
	ListGraph::NodeMap<ListGraph::Node> ncr(solution_graph);
	ListGraph::EdgeMap<ListGraph::Edge> ecr(solution_graph);
	ListGraph::NodeMap<ListGraph::Node> nr(g);
	graphCopy(filterEdges(g, disable_all_edges), solution_graph).nodeRef(nr).nodeCrossRef(ncr).run();


	double cost_F_M = 0;
	double cost_F = 0;
	double cost_M = 0;

	// add orig edge from F+M
	for (ListGraph::EdgeIt e(g1); e != INVALID; ++e)
	{
		auto u1 = g1.u(e);
		auto v1 = g1.v(e);

		auto orig_e = find_orig_edge(u1, v1, is_depot_g1, orig_node, orig_edge);
		auto orig_u = g.u(orig_e);
		auto orig_v = g.v(orig_e);
		if (forest_edge_map[e])
		{
			auto sol_e = solution_graph.addEdge(nr[orig_u], nr[orig_v]);
			ecr[sol_e] = orig_e;
			cost_F = cost_F + edge_cost[orig_e];
		}
		if (max_weighted_perfect_matching.matching(e) && odd_degree_node_filter[u1] && odd_degree_node_filter[v1])
		{
			auto sol_e = solution_graph.addEdge(nr[orig_u], nr[orig_v]);
			ecr[sol_e] = orig_e;
			cost_M = cost_M + edge_cost[orig_e];
		}
	}

	cost_F_M = cost_F + cost_M;


	int cost_R = 0;
	// add orig edge from required edges
	for (ListGraph::EdgeIt orig_e(g); orig_e != INVALID; ++orig_e)
	{
		if (is_required[orig_e])
		{
			auto orig_u = g.u(orig_e);
			auto orig_v = g.v(orig_e);

			auto sol_u = nr[orig_u];
			auto sol_v = nr[orig_v];

			int degree_sol_u = 0;
			int degree_sol_v = 0;
			for (ListGraph::IncEdgeIt sol_e(solution_graph, sol_u); sol_e != INVALID; ++sol_e)
			{
				degree_sol_u++;
			}
			for (ListGraph::IncEdgeIt sol_e(solution_graph, sol_v); sol_e != INVALID; ++sol_e)
			{
				degree_sol_v++;
			}

			if (degree_sol_u % 2 == 1 && degree_sol_v % 2 == 1)
			{
				auto sol_e = solution_graph.addEdge(sol_u, sol_v);
				ecr[sol_e] = orig_e;
				cost_R = cost_R + edge_cost[orig_e];
			}
			else if (degree_sol_u % 2 == 0 && degree_sol_v % 2 == 0)
			{
				auto sol_e = solution_graph.addEdge(sol_u, sol_v);
				ecr[sol_e] = orig_e;
				sol_e = solution_graph.addEdge(sol_u, sol_v);
				ecr[sol_e] = orig_e;
				cost_R = cost_R + 2*edge_cost[orig_e];
			}
			else
			{
				std::cout << "fatal error! F+M is not not Eulerian graph!" << std::endl;
			}
		}
	}

	int cost = 0;
	sol1.second.clear();
	ListGraph::NodeMap<int>component_label(solution_graph);
	int n_component = connectedComponents(solution_graph, component_label);
	for (int i = 0; i < n_component; ++i)
	{
		bool has_depot = false;
		std::vector<ListGraph::Edge> walk;
		ListGraph::NodeMap<bool> component_filter(solution_graph, false);
		for (ListGraph::NodeIt u(solution_graph); u != INVALID; ++u)
		{
			if (component_label[u] == i)
			{
				component_filter[u] = true;
				if (is_depot[ncr[u]])
				{
					has_depot = true;
				}
			}
		}

		if (!has_depot)
		{
			std::cout << "fatal error! no depot!" << std::endl;
		}

		FilterNodes<ListGraph> component(solution_graph, component_filter);



		if (!eulerian(component))
		{
			std::cout << "fatal error! not Eulerian graph!" << std::endl;
		}
		for (EulerIt<FilterNodes<ListGraph>> e(component); e != INVALID; ++e)
		{
			auto orig_e = ecr[e];
			walk.push_back(orig_e);
			cost = cost + edge_cost[orig_e];
		}
		sol1.second.push_back(walk);
	}
	sol1.first = cost;


	double lb_C_R_star = 0;
	for (ListGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		if (is_required[e])
		{
			lb_C_R_star = lb_C_R_star + edge_cost[e];
		}
	}

	if(_debug)
	{
		std::cout << "F: " << cost_F << std::endl;
		std::cout << "M: " << cost_M << std::endl;
		std::cout << "c_added_R: " << cost_R << std::endl;
	}

	double alpha = 2.0 - 1.0 / static_cast<double>(n_depots);
	lb_opt_1 = (cost_F + cost_M) / alpha + lb_C_R_star;

	double beta = 2.0 * static_cast<double>(n_depots) - 2.0;

	double l_max = (2.0 * cost_M - cost_F) / (2.0 * static_cast<double>(n_depots) - 1.0);

	alpha = (cost_F + cost_M) / std::max(cost_F + l_max, 2.0 * cost_M - beta * l_max);
	lb_opt_1_1 = (cost_F + cost_M) / alpha + lb_C_R_star;


}
