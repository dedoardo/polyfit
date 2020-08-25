#include <polyvec/tracer-public.hpp>
#include <polyvec/curve-tracer/find-fitting-points.hpp>
#include <dev/pdf.hpp>
#include "../../apps/apps_polyvec/drawing.hpp"
#include <polyvec/polygon-tracer/minimum.hpp>

#include <stdlib.h>
#include <time.h>

using namespace polyfit;
using namespace polyvec;

NAMESPACE_BEGIN(polyvectest)
NAMESPACE_BEGIN(CurveTracer)

int fit_points(int argc, char ** argv)
{
    const char* image_uri = "/home/sparkon/code/polyvec/baseline-images/binary-color/axe-32.png";

    Log::open((char*)stdout);

    IO::Image I;

    if (!IO::read_image(image_uri, I))
    {
        PF_LOGF("image not found %s", image_uri);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "image %s", image_uri);
    PF_LOGF("read image %d %d from %s", I.width, I.height, image_uri);

    // extracts overlapping closed boundaries
    std::vector<mat2x> boundaries;
    std::vector<vec4> colors;
    ImageSegment::expand_and_cleanup(I);
    ImageSegment::extract_closed_regions(I, boundaries, colors);

    if (!boundaries.size())
    {
        PF_LOGS("boundary extraction failed");
        return EXIT_FAILURE;
    }

    // setting default trace options
    VectorOptions opt;
    opt.error_continuity_limit = PF_RAD(60);
    opt.error_continuity_weight = .25;
    opt.error_inflection_limit = PF_RAD(90);
    opt.error_inflection_weight = .25;
    opt.error_smoothness_limit = PF_RAD(135);
    opt.error_smoothness_weight = 1.;
    opt.error_accuracy_weight = 1.;
    opt.accuracy_thr = 1.0;
    opt.accuracy_metric = AccuracyMetric::Implicit;
    opt.make_global();

    // finds the closed boundary most likely to be a binary image
    mat2x &R = boundaries[ImageSegment::find_binary_color_region(boundaries)];
    PF_LOGF("boundary has %lld points", R.cols());

    // polygon vertices
    vecXi PV;

    // fitting polyline
    mat2x B;
    std::vector<BoundaryGraph::Edge> E;

    if (!PolygonTracer::minimum(R, PV, B, E, true))
    {
        PF_LOGS("failed to trace polygon");
        return EXIT_FAILURE;
    }

    PF_LOGF("traced %d vertices", PV.size());

    // polygon points
    mat2x P;

    // converts the shortest cycle indices into primary graph components
    BoundaryGraph::trace_to_points(B, PV, P);

    srand(1519251251);

    int e_start = rand() % PV.size();
    int e_end = rand() % PV.size();
    double t_start = (double)rand() / RAND_MAX;
    double t_end = (double)rand() / RAND_MAX;

    PF_LOGF("e_start %d t_start %f e_end %d t_end %f", e_start, t_start, e_end, t_end);

    mat2xi vertices;
    mat2x midpoints;
    PF_ASSERT(::polyfit::CurveTracer::find_pixel_centers_for_subpath(B, PV, e_start, t_start, e_end, t_end, vertices, &midpoints, true));
    PF_LOGF("found %d vertices", vertices.cols());

#if 1
    DevicePDF* pdf = new DevicePDF(misc::sfmt("/home/sparkon/data/polyvec/out/fit-points.pdf").c_str(), 1, 1);
    draw_raster(B);
    draw_raster_indices(B);

    vec2 src0 = B.col(PV(e_start));
    vec2 src1 = B.col(PV((e_start + 1) %PV.size()));

    vec2 dst0 = B.col(PV(e_end));
    vec2 dst1 = B.col(PV((e_end + 1) % PV.size()));

    int e = e_start;
    while (e != e_end) {
        int enext = (e + 1) % PV.size();
        vec2 p0 = B.col(PV(e));
        vec2 p1 = B.col(PV(enext));

        draw::line(p0, p1, Style::outline(colors::calm_blue, 2.));
        e = enext;
    }

    draw::point(misc::lerp(src0, src1, t_start), 10., Style::fill(colors::forest_green));
    draw::point(misc::lerp(dst0, dst1, t_end), 10., Style::fill(colors::green));

    for (int j = 0; j < vertices.cols(); ++j) {
        draw::point(midpoints.col(j), 5., Style::fill(colors::red));
    }
    pdf->draw(0, 0);
    delete pdf;

#endif

//    std::vector<int> expected;
//    for (int i = 0; i < vertices.size(); ++i) {
//        PF_ASSERT(vertices(i) == expected[i]);
//    }

    return EXIT_SUCCESS;
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvectest)