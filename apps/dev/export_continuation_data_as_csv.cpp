/*
    Exports all the error data for the continuation in a column. If the file already exists,
    it won't be overwritten.   
*/
const char* usage =
"<image_uri>      Path to the image file to be vectorized\n"
"<classifier_uri> Trained corner classification model\n"
"<csv_file>       Path to the file where the CSV data will be written.";

// Polyvec
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/system.hpp>

// libc++
#include <cstdlib>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace polyfit;
using namespace polyvec;

struct CSVColumn {
    using EvalFn = std::function<string(const Regularity::EdgeData::Continuation& r)>;
    CSVColumn(const string& name, const EvalFn& eval) :
        name(name), eval(eval) { }

    string name;
    EvalFn eval;
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        return EXIT_FAILURE;
    }

    const char* image_uri = argv[1];
    const char* classifier_uri = argv[2];
    const char* csv_uri = argv[3];

    init_pipeline();

    InputData input_data = load_input(image_uri);
    PolygonData polygon_data = extract_polygon(input_data);
    std::vector<polyvec::AttemptInfo> attempts;
    SplineData spline_data = extract_curves_with_classifier(polygon_data, classifier_uri, &attempts);

    std::vector<CSVColumn> csv_columns;
    csv_columns.emplace_back("id", [image_uri](const Regularity::EdgeData::Continuation& r) -> string {
        return StringUtils::fmt("%s-%d-%d-%d-%d", image_uri, r.v0_prev, r.v0, r.v1, r.v1_next).c_str();
    });

    csv_columns.emplace_back("angle_polygon_0", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(r.angle_polygon_0);
    });

    csv_columns.emplace_back("angle_polygon_1", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(r.angle_polygon_1);
    });

    csv_columns.emplace_back("angle_polygon_max", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(max(r.angle_polygon_0, r.angle_polygon_1));
    });

    csv_columns.emplace_back("angle_polygon_min", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(min(r.angle_polygon_0, r.angle_polygon_1));
    });

    csv_columns.emplace_back("angle_polygon_avg", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(.5 * (r.angle_polygon_0 + r.angle_polygon_1));
    });

    csv_columns.emplace_back("angle_continuation_0", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(r.angle_continuation_0);
    });

    csv_columns.emplace_back("angle_continuation_1", [](const Regularity::EdgeData::Continuation& r) -> string {
        return to_string(r.angle_continuation_1);
    });

    csv_columns.emplace_back("continuation_angle_0", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(r.angle_continuation_0);
    });

    csv_columns.emplace_back("continuation_angle_1", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(r.angle_continuation_1);
    });

    csv_columns.emplace_back("continuation_angle_max", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(max(r.angle_continuation_0, r.angle_continuation_1));
    });

    csv_columns.emplace_back("continuation_angle_min", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(min(r.angle_continuation_0, r.angle_continuation_1));
    });

    csv_columns.emplace_back("continuation_angle_avg", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(.5 * (r.angle_continuation_0 + r.angle_continuation_1));
    });

    csv_columns.emplace_back("separation_angle_0", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(r.angle_separation_0);
    });

    csv_columns.emplace_back("separation_angle_1", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(r.angle_separation_0);
    });

    csv_columns.emplace_back("separation_angle_min", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(min(r.angle_separation_0, r.angle_separation_1));
    });

    csv_columns.emplace_back("midpoint_distance", [](const Regularity::EdgeData::Continuation& r)->string {
        return to_string(sqrt(r.distance_midpoints_sq));
    });

    bool write_column_names = !os::file_exists(csv_uri);

    std::fstream csv(csv_uri, std::fstream::out | std::fstream::app);
    if (!csv.good()) {
        fprintf(stderr, "Failed to open csv file %s\n", csv_uri);
        return EXIT_FAILURE;
    }

    // Writing column names
    if (write_column_names) {
        for (size_t i = 0; i < csv_columns.size(); ++i) {
            csv << csv_columns[i].name;
            if (i != csv_columns.size() - 1) {
                csv << ", ";
            }
        }
        csv << std::endl;
    }

    // Writing rows
    for (const auto& e : polygon_data.RE) {
        if (e.type != Regularity::Type::CONTINUATION) {
            continue;
        }

        for (size_t i = 0; i < csv_columns.size(); ++i) {
            csv << csv_columns[i].eval(e.data.continuation);
            if (i != csv_columns.size() - 1) {
                csv << ", ";
            }
        }
        csv << std::endl;
    }

    csv.close();
    return EXIT_SUCCESS;
}
