// Header
#include <polyvec/io/ascii.hpp>

#include <polyvec/api.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/misc.hpp>

#include <fstream>

using namespace polyvec;
using namespace std;

namespace polyfit {
	namespace IO {
		void read_raster(const std::string& addr, Eigen::Matrix2Xd& points) {
			std::ifstream ifs(addr);
			assert_break(ifs.good());

			std::string line;

			while (std::getline(ifs, line)) {
				float x, y;
				sscanf(line.c_str(), "%f %f", &x, &y);

				points.conservativeResize(2, points.cols() + 1);
				points(0, points.cols() - 1) = (double)x;
				points(1, points.cols() - 1) = (double)y;
			}
		}

		bool read_polygon(const std::string & addr, std::vector<Eigen::Index> & indices) {
			std::ifstream ifs(addr);
			if (!ifs.good()) {
				return false;
			}

			indices.clear();
			std::string line;

			while (std::getline(ifs, line)) {
				int idx;
				sscanf(line.c_str(), "%d", &idx);
				indices.emplace_back(idx);
			}

			return true;
		}

		void read_polygon_g0(const std::string& addr, std::vector<Eigen::Index>& indices) {
			std::ifstream ifs(addr);
			assert_break(ifs.good());

			indices.clear();
			std::string line;

			while (std::getline(ifs, line)) {
				int idx;
				sscanf(line.c_str(), "%d", &idx);
				indices.emplace_back(idx);
			}
		}

		void read_manual_polygon(const std::string& raster_addr, const std::string& polygon_addr, Eigen::Matrix2Xd& raster, Eigen::Matrix2Xd& polygon, std::vector<Eigen::Index>& polygon_vertices) {
			read_raster(raster_addr, raster);

			std::vector<Eigen::Index> indices;
			read_polygon(polygon_addr, indices);

			polygon.resize(2, indices.size());
			polygon_vertices.resize(indices.size());

			int i = 0;
			for (auto idx : indices) {
				vec2 point;
				if (idx >= raster.cols()) {
					idx -= raster.cols();
					const vec2 pt_prev = raster.col(idx);
					const vec2 pt_next = CircularAt(raster, idx + 1);

					point = misc::lerp(pt_prev, pt_next, .5);
				}
				else {
					point = raster.col(idx);
				}

				polygon.col(i) = point;
				polygon_vertices[i] = idx;
				++i;
			}
		}

		void read_baseline_boundary(const std::string& uri, mat2x& B) {
			read_raster(uri, B);
		}

		void read_baseline_polygon(const std::string& uri, vecXi& P) {
			std::ifstream ifs(uri);
			assert_break(ifs.good());

			P.resize(0);

			std::string line;

			while (std::getline(ifs, line)) {
				int idx;
				sscanf(line.c_str(), "%d", &idx);
				MatrixUtils::append(P, idx);
			}
		}

		bool write_baseline_boundary(const std::string& uri, const mat2x& B) {
			FILE* fp = fopen(uri.c_str(), "w");
			if (!fp) {
				return false;
			}

			for (Vertex i = 0; i < B.cols(); ++i) {
				fprintf(fp, "%.3f %.3f\n", B(0, i), B(1, i));
			}

			return true;
		}

		bool write_baseline_polygon(const std::string& uri, const vecXi& P) {
			FILE* fp = fopen(uri.c_str(), "w");
			if (!fp) {
				return false;
			}

			for (Vertex i = 0; i < P.size(); ++i) {
				fprintf(fp, "%d\n", P(i));
			}

			return true;
		}

		std::string model_id(const std::string& image_uri) {
			const char* b = &image_uri[0];
			const char* r = &image_uri[image_uri.size() - 1];
			string id;

			// skip extension
			while (r != b && *r != '.') {
				--r;
			}

			--r;

			// filename to directory
			while (r != b && *r != '\\'  && *r != '/') {
				id += *r--;
			}

			// replace directory separator
			id += '-';
			
			// skipping separator
			while (r != b && (*r == '\\' || *r == '/')) {
				--r;
			}

			// directory name
			while (r != b && *r != '\\' && *r != '/') {
				id += *r--;
			}

			reverse(id.begin(), id.end());
			return id;
		}

		bool write_points(const std::string& uri, const Eigen::MatrixXd& P) {
			FILE* fp = fopen(uri.c_str(), "w");
			if (!fp) {
				return false;
			}

			for (Vertex i = 0; i < P.cols(); ++i) {
				for (Vertex j = 0; j < P.rows(); ++j) {
					fprintf(fp, "%f", P(j, i));

					if (j < P.rows() - 1) {
						fprintf(fp, " ");
					}
				}

				fprintf(fp, "\n");
			}

			return true;
		}


		bool read_parameter_limits(const std::string& uri, double& lim_accuracy, double& lim_smoothness, double& lim_continuity, double& lim_inflection) {
			std::ifstream ifs(uri);
			if (!ifs.good()) {
				return false;
			}

			string line;
			std::getline(ifs, line);


			return true;
		}
	}
}