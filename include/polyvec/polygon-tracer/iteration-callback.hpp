#pragma once

#include <functional>

using IterationCallback = std::function<void(
	const polyfit::mat2x& B,
	const polyfit::vecXi& P,
	const std::vector<polyfit::BoundaryGraph::Edge>& E
	)>;