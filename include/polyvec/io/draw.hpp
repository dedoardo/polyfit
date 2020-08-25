#pragma once

// polyvec
#include <polyvec/api.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/curve-tracer/curve.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Draw)

void indices(const mat2x& B);
void indices(const mat2x& B, const vecXi& P);
void boundary ( const mat2x& B, bool circular );
void graph(const mat2x& B, const std::vector<BoundaryGraph::Edge>& E);
void polygon(const mat2x& B, const vecXi& P, bool circular);
void polygon ( const mat2x& P, bool circular );
void polygon_fill ( const mat2x& P, bool circular );
void curve(polyvec::GlobFitCurve* curve);

void closure_same_line(const mat2x& B, const vecXi& P, Vertex i, Vertex j, bool circular);
void closure_same_circle(const mat2x& B, const vecXi& P, Vertex i, Vertex j, bool circular);

NAMESPACE_END ( Draw )
NAMESPACE_END ( polyfit )