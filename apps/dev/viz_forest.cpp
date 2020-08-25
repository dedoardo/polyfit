// Polyvec
#include <polyvec/curve-tracer/fit_classifier.hpp>

// libc++
#include <fstream>
#include <iostream>
#include <map>

using namespace polyfit;
using namespace polyvec;

struct TreeIterationInfo
{
	const andres::ml::DecisionTree<>* tree;
	std::vector<const andres::ml::DecisionNode<double, unsigned char>*> current_level_nodes;
};

void group_trees_by_structure(const std::vector<TreeIterationInfo>& trees, std::vector<std::vector<const andres::ml::DecisionTree<>*>>& out_groups)
{
	std::map<std::vector<unsigned char>, std::vector<TreeIterationInfo>> grouping; //group trees by feature indices of the current level, leaves will be -1
	//group all trees by the structure of the current level
	for (auto& tree : trees)
	{
		std::vector<uint8_t> characterization;
		for (auto node : tree.current_level_nodes)
		{
			if (node->isLeaf())
				characterization.push_back(-1);
			else
				characterization.push_back(node->featureIndex());
		}
		grouping[characterization].push_back(tree);
	}

	//now generate for every group the next levels
	for (auto& group : grouping)
	{
		std::vector<TreeIterationInfo> next_levels(group.second.size());
		bool has_next_level = false;
		//for every tree in the current group
		for (int i = 0; i < group.second.size(); ++i)
		{
			auto tree = group.second[i].tree;
			next_levels[i].tree = tree;
			for (auto& node : group.second[i].current_level_nodes)
			{
				if (node->isLeaf())
					continue;
				//add the child nodes to the next level
				next_levels[i].current_level_nodes.push_back(&tree->decisionNode(node->childNodeIndex(0)));
				next_levels[i].current_level_nodes.push_back(&tree->decisionNode(node->childNodeIndex(1)));
				has_next_level = true;
			}
		}

		if (has_next_level)
			group_trees_by_structure(next_levels, out_groups);
		else
		{
			//we finished analyzing this group, report the result
			out_groups.emplace_back();
			for (auto& tree : group.second)
			{
				out_groups.back().push_back(tree.tree);
			}
		}

	}
}

std::vector<std::vector<const andres::ml::DecisionTree<>*>> group_trees_by_structure(const andres::ml::DecisionForest<>& forest)
{	
	std::vector<TreeIterationInfo> trees(forest.size());
	for (int i = 0; i < forest.size(); ++i)
	{
		auto& tree = forest.decisionTree(i);
		trees[i].tree = &tree;
		trees[i].current_level_nodes.push_back(&tree.decisionNode(0));
	}
	std::vector<std::vector<const andres::ml::DecisionTree<>*>> groups;
	group_trees_by_structure(trees, groups);
	return groups;
}

void print_subtree_group(const std::vector<const andres::ml::DecisionNode<double, unsigned char>*>& roots, const std::vector<const andres::ml::DecisionTree<>*>& trees, int indent = 0)
{
	double feature_sum = 0;
	double square_sum = 0;
	size_t label_0 = 0;
	for (auto node : roots)
	{
		if (node->isLeaf())
		{
			if (node->label() == 0)
				++label_0;
		}
		else
		{
			feature_sum += node->threshold();
			square_sum += node->threshold() * node->threshold();
		}
	}
	for (int i = 0; i < indent; ++i)
		if (i == indent - 1)
			std::cout << "+-";
		else
			std::cout << "| ";
	if (roots.front()->isLeaf())
		std::cout << "Leaf (" << label_0 << " out of " << roots.size() << " with label 0)" << std::endl;
	else
	{
		double average = feature_sum / roots.size();
		double variance = square_sum / roots.size() - average * average;
		std::cout << "Feature " << roots.front()->featureIndex() << " (avg threshold: " << average << ", std dev: " << std::sqrt(variance) << ")" << std::endl;
		std::vector < const andres::ml::DecisionNode<double, unsigned char>*> left_roots(roots.size());
		std::vector < const andres::ml::DecisionNode<double, unsigned char>*> right_roots(roots.size());
		for(int i = 0; i < roots.size(); ++i)
		{
			auto r = roots[i];
			left_roots[i] = &trees[i]->decisionNode(r->childNodeIndex(0));
			right_roots[i] = &trees[i]->decisionNode(r->childNodeIndex(1));
		}
		print_subtree_group(left_roots, trees, indent + 1);
		print_subtree_group(right_roots, trees, indent + 1);
	}
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Usage: <forest:file>" << std::endl;
		return EXIT_FAILURE;
	}

	const std::string model_uri = argv[1];
	FitClassifierRandomForest classifier;
	classifier.load_from_file(model_uri);	

	auto groups = group_trees_by_structure(classifier.get_forest());
	std::cout << "Forest has trees with " << groups.size() << " distinct structures." << std::endl;

	std::sort(groups.begin(), groups.end(), 
		[](const std::vector<const andres::ml::DecisionTree<>*>& first, const std::vector<const andres::ml::DecisionTree<>*>& second) { 
			return first.size() > second.size(); 
	});

	for (int i = 0; i < std::min<size_t>(10, groups.size()); ++i)
	{
		std::cout << "Structure " << i << " (" << groups[i].size() << " trees):" << std::endl;
		std::vector<const andres::ml::DecisionNode<double, unsigned char>*> roots;
		for (auto& t : groups[i])
			roots.push_back(&t->decisionNode(0));
		print_subtree_group(roots, groups[i]);
		std::cout << std::endl;
	}	

	int leaves = 0;
	for (int i_tree = 0; i_tree < classifier.get_forest().size(); ++i_tree)
	{
		auto& tree = classifier.get_forest().decisionTree(i_tree);
		for (int i_node = 0; i_node < tree.size(); ++i_node)
		{
			auto& node = tree.decisionNode(i_node);
			if (node.isLeaf())
				++leaves;
		}
	}
	std::cout << "Forest has " << ((double)leaves / classifier.get_forest().size()) << " leaves per tree in average." << std::endl;
}