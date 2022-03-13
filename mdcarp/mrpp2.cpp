#include "mdcarp.h"
#include <lemon/matching.h>
#include <lemon/euler.h>
#include <lemon/adaptors.h>
#include <lemon/full_graph.h>
#include <lemon/kruskal.h>


void MDRPP::solve_by_mrpp_2()
{
	ListGraph ddg;
	auto u = ddg.addNode();
	auto v = ddg.addNode();
	ddg.addEdge(u, v);
	ListGraph::NodeMap<ListGraph::Node> ncr_ddg(ddg);
	ListGraph::NodeMap<ListGraph::Node> nr_dgg(g);
	ListGraph::EdgeMap<ListGraph::Edge> ecr_dgg(ddg);
	ListGraph::EdgeMap<ListGraph::Edge> er_dgg(g);
	ListGraph::EdgeMap<int> edge_cost_ddg(ddg, 0);
	ListGraph::EdgeMap<bool> is_in_E1(ddg, false);
	ListGraph::EdgeMap<bool> is_in_E2(ddg, false);
	ddg.clear();

	ListGraph::EdgeMap<bool> disable_all_edges(g);
	graphCopy(filterEdges(g, disable_all_edges), ddg).nodeCrossRef(ncr_ddg).edgeCrossRef(ecr_dgg).nodeRef(nr_dgg).edgeRef(er_dgg).run();

	for (ListGraph::NodeIt u(ddg); u != INVALID; ++u)
	{
		if (is_depot[ncr_ddg[u]])
		{
			auto d = ddg.addNode();
			ncr_ddg[d] = ncr_ddg[u];
			auto e_dd = ddg.addEdge(d, u);
			is_in_E2[e_dd] = true;
			ecr_dgg[e_dd] = INVALID;
			for (ListGraph::IncEdgeIt e(g, ncr_ddg[u]); e != INVALID; ++e)
			{
				if (is_depot[g.oppositeNode(ncr_ddg[u], e)])
				{
					continue;
				}
				auto v = nr_dgg[g.oppositeNode(ncr_ddg[u], e)];
				auto e1 = ddg.addEdge(u, v);
				is_in_E1[e1] = true;
				edge_cost_ddg[e1] = -edge_cost[e];
				ecr_dgg[e1] = e;
				auto e2 = ddg.addEdge(d, v);
				is_in_E1[e2] = true;
				edge_cost_ddg[e2] = -edge_cost[e];
				ecr_dgg[e2] = e;
			}
		}
		else
		{
			for (ListGraph::IncEdgeIt e(g, ncr_ddg[u]); e != INVALID; ++e)
			{
				if (is_depot[g.oppositeNode(ncr_ddg[u], e)])
				{
					continue;
				}
				auto v = nr_dgg[g.oppositeNode(ncr_ddg[u], e)];
				auto e1 = ddg.addEdge(u, v);
				is_in_E1[e1] = true;
				edge_cost_ddg[e1] = -edge_cost[e];
				ecr_dgg[e1] = e;
			}
		}
	}

	MaxWeightedPerfectMatching<ListGraph> mwpm(ddg, edge_cost_ddg);
	bool result = mwpm.run();
	if (!result)
	{
		std::cout << "fatal error!!!!!!!!!!!!!" << std::endl;
	}


	ListGraph h;
	ListGraph::NodeMap<ListGraph::Node> ncr_h(h);
	ListGraph::NodeMap<ListGraph::Node> nr_h(g);
	ListGraph::EdgeMap<ListGraph::Edge> ecr_h(h);
	ListGraph::EdgeMap<ListGraph::Edge> er_h(g);

	graphCopy(filterEdges(g, is_required), h).nodeCrossRef(ncr_h).edgeCrossRef(ecr_h).nodeRef(nr_h).edgeRef(er_h).run();

	double cost_R = 0;

	for (ListGraph::EdgeIt e(h); e != INVALID; ++e)
	{
		cost_R = cost_R + edge_cost[ecr_h[e]];
	}

	double cost_M = 0;
	for (ListGraph::EdgeIt e(ddg); e != INVALID; ++e)
	{
		if (mwpm.matching(e) && !is_in_E2[e])
		{
			auto u_h = nr_h[g.u(ecr_dgg[e])];
			auto v_h = nr_h[g.v(ecr_dgg[e])];
			auto new_edge = h.addEdge(u_h, v_h);
			ecr_h[new_edge] = ecr_dgg[e];
			er_h[ecr_dgg[e]] = new_edge;
			cost_M = cost_M + edge_cost[ecr_h[new_edge]];
		}
	}

	ListGraph::NodeMap<int>component_label(h);
	int n_component = connectedComponents(h, component_label);
	int n_component_with_depot = 0;
	std::vector<bool> has_depot(n_component, false);

	for (ListGraph::NodeIt u(h); u != INVALID; ++u)
	{
		if (is_depot[ncr_h[u]])
		{
			has_depot[component_label[u]] = true;
		}
	}
	for (int i = 0; i < n_component; ++i)
	{
		if(has_depot[i])
		{
			n_component_with_depot++;
		}
	}

	FullGraph contract_graph(n_component - n_component_with_depot + 1);
	FullGraph::EdgeMap<ListGraph::Edge>shortest_cut_edge(contract_graph, INVALID);
	FullGraph::EdgeMap<int>shortest_cut_edge_cost(contract_graph, 0);
	std::vector<FullGraph::Node> contract_nodes;
	for (FullGraph::NodeIt u(contract_graph); u != INVALID; ++u)
	{
		contract_nodes.push_back(u);
	}

	std::map<int, FullGraph::Node> component_label_to_node;
	int idx_u = 1;
	for (int i = 0; i < n_component; ++i)
	{
		if (has_depot[i])
		{
			component_label_to_node[i] = contract_nodes[0];
		}
		else
		{
			component_label_to_node[i] = contract_nodes[idx_u];
			idx_u++;
		}
	}

	for (ListGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		auto u = g.u(e);
		auto v = g.v(e);
		auto u_h = nr_h[u];
		auto v_h = nr_h[v];
		auto l_u = component_label[u_h];
		auto l_v = component_label[v_h];
		auto u_c = component_label_to_node[l_u];
		auto v_c = component_label_to_node[l_v];
		if ( u_c!= v_c)
		{
			int cost = edge_cost[e];
			auto contract_edge = findEdge(contract_graph, u_c, v_c);
			if (shortest_cut_edge[contract_edge] == INVALID || shortest_cut_edge_cost[contract_edge] > cost)
			{
				shortest_cut_edge[contract_edge] = e;
				shortest_cut_edge_cost[contract_edge] = cost;
			}
		}
	}

	std::vector<FullGraph::Edge> tree_edges(countNodes(contract_graph) - 1);
	kruskal(contract_graph, shortest_cut_edge_cost, tree_edges.begin());

	double cost_E = 0;
	for (auto e : tree_edges)
	{
		auto e_g = shortest_cut_edge[e];
		auto u_g = g.u(e_g);
		auto v_g = g.v(e_g);
		auto u_h = nr_h[u_g];
		auto v_h = nr_h[v_g];
		auto e_h = h.addEdge(v_h, u_h);
		ecr_h[e_h] = e_g;
		er_h[e_g] = e_h;
		e_h = h.addEdge(v_h, u_h);
		ecr_h[e_h] = e_g;
		er_h[e_g] = e_h;
		cost_E = cost_E + edge_cost[e_g];
	}

	int cost = 0;
	sol2.second.clear();
	n_component = connectedComponents(h, component_label);
	for (int i = 0; i < n_component; ++i)
	{
		std::vector<ListGraph::Edge> walk;
		ListGraph::NodeMap<bool> component_filter(h, false);
		for (ListGraph::NodeIt u(h); u != INVALID; ++u)
		{
			if (component_label[u] == i)
			{
				component_filter[u] = true;
			}
		}
		FilterNodes<ListGraph> component(h, component_filter);

		if (!eulerian(component))
		{
			for (FilterNodes<ListGraph>::EdgeIt e(component); e != INVALID; ++e)
			{
				std::cout << component.id(component.u(e)) << "--" << component.id(component.v(e)) << std::endl;

			}
			std::cout << "fatal error! not Eulerian graph!" << std::endl;
		}
		for (EulerIt<FilterNodes<ListGraph>> e(component); e != INVALID; ++e)
		{
			auto orig_e = ecr_h[e];
			walk.push_back(orig_e);
			cost = cost + edge_cost[orig_e];
		}
		sol2.second.push_back(walk);
	}
	sol2.first = cost;

	lb_opt_2 = std::max(cost_M, cost_E) + cost_R;
}

