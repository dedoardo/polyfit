#include <polyvec/image-segment/image_segment.hpp>

#include <polyvec/core/log.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/geometry/path.hpp>

// libc++
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <vector>
#include <memory>

#define array_len(x) (sizeof(x)/sizeof((x)[0]))
#define index_null ((Vertex)~0)
#define int64_max std::numeric_limits<int64_t>::max()

using namespace std;

namespace polyfit {
    uint32_t pack ( const uint8_t* xyzw ) {
        uint32_t ret = 0x0;
        ret |= ( uint32_t ) xyzw[0] << 0;
        ret |= ( uint32_t ) xyzw[1] << 8;
        ret |= ( uint32_t ) xyzw[2] << 16;
        ret |= ( uint32_t ) xyzw[3] << 24;
        return ret;
    }

    void unpack ( uint32_t xyzw, uint8_t& x, uint8_t& y, uint8_t& z, uint8_t& w ) {
        x = ( uint8_t ) ( ( xyzw >> 0 ) & 0xff );
        y = ( uint8_t ) ( ( xyzw >> 8 ) & 0xff );
        z = ( uint8_t ) ( ( xyzw >> 16 ) & 0xff );
        w = ( uint8_t ) ( ( xyzw >> 24 ) & 0xff );
    }

    vec4i unpack ( uint32_t xyzw ) {
        return {
            ( int ) ( ( xyzw >> 0 ) & 0xff ),
            ( int ) ( ( xyzw >> 8 ) & 0xff ),
            ( int ) ( ( xyzw >> 16 ) & 0xff ),
            ( int ) ( ( xyzw >> 24 ) & 0xff )
        };
    }

    enum PixelAdjBit { // There is a little redundancy in `Nhb` and `Adj`, just stick to one (todo)
        None = 0,
        Top = 1 << 0,
        TopRight = 1 << 1,
        Right = 1 << 2,
        BottomRight = 1 << 3,
        Bottom = 1 << 4,
        BottomLeft = 1 << 5,
        Left = 1 << 6,
        TopLeft = 1 << 7
    };

    struct PixelGraphNode {
        static const PixelAdjBit neighbors[8];
        static const PixelAdjBit aa_neighbors[4];

        void        init_from_rgba ( int64_t xy, uint8_t r, uint8_t g, uint8_t b, uint8_t a);
        bool        connected ( PixelAdjBit adj ) const;
        void        disconnect ( int bits );
        int32_t     valence() const;
        PixelAdjBit next ( PixelAdjBit adj_ignore ) const;

        // Color information
        uint8_t& y();
        uint8_t& u();
        uint8_t& v();
        uint8_t y() const;
        uint8_t u() const;
        uint8_t v() const;

        // todo: convert to single neighbor(PixelAdjBit)...
        // careful! diagonal accessor don't do what you think..
        int64_t tl() const;
        int64_t tr() const;
        int64_t br() const;
        int64_t bl() const;

        vec2i point() const;

        // It could be a little more compact.. (xy, yuv, alpha) should be known by the processing class..
        int64_t xy; // packed xy
        uint8_t yuv[3];
        uint8_t alpha;
        uint8_t adj;
    };

    struct PixelRegion;
    struct PixelBoundaryNode;

    // A pixel boundary is closed and every point is overlapping with other boundaries,
    // it can be shared between multiple pixel regions
    struct PixelRegionBoundary {
        static int compare ( const PixelRegionBoundary& lhs, const PixelRegionBoundary& rhs );
        std::vector<int64_t> points;
    };

    // Unique boundary segments, shared by only two pixel segments
    struct PixelBoundary {
        PixelRegion* region_0;
        PixelRegion* region_1;
        PixelBoundaryNode* src;
        PixelBoundaryNode* dst;
        std::vector<int64_t>       points; // includes `src` and `dst`
    };

    enum BoundaryContinuity {
        G0,
        G1
    };

    // Valence 3/4 pixels between multiple colored regions
    struct PixelBoundaryNode {
        int64_t             xy;         // Location
        std::vector<PixelBoundary*> boundaries; // *Unique* outgoing boundaries from the node
        std::vector<PixelRegion*>   regions;    // Regions overlapping the node

        Vertex count_unique_regions() const;
        vec2i  point() const;
    };

    // Hierarchical pixel region which is delimited by a closed boundary
    // and optionally encompasses. Hierarchy is necessary when exporting the
    // vectorized result to make sure that regions are drawn in the correct order
    // The boundaries stored in the pixel region are the duplicate originals
    struct PixelRegion {
        // Packs region Vertex in 64 bit integer
        static void init ( uint64_t* bits );
        static bool add ( uint64_t* bits, uint8_t region );
        static int  count ( uint64_t bits );
        static int  count_unique ( uint64_t bits );
        static int  at ( uint64_t bits, int idx );

        uint8_t index;
        uint8_t yuv[3];
        uint8_t alpha;
        std::vector<PixelRegionBoundary*> boundaries;
    };

    // Origin at the topleft of the image
    struct PixelGraph {
        static int64_t pack ( vec2i xy );
        static vec2i unpack ( int64_t xy );
        static vec2 unpack_edge ( int64_t xy );

        // Opposite adjacent bit
        static PixelAdjBit opposite ( PixelAdjBit adj );

        // Access
        const PixelGraphNode* at ( int32_t x, int32_t y ) const;
        PixelGraphNode* at ( int32_t x, int32_t y );
        const PixelGraphNode* at ( vec2i xy ) const;
        PixelGraphNode* at ( vec2i xy );
        const PixelGraphNode* at ( int64_t xy ) const;
        PixelGraphNode* at ( int64_t xy );

        // Query neighbor
        PixelGraphNode* nhb ( PixelGraphNode* node, PixelAdjBit adj );
        const PixelGraphNode* nhb ( const PixelGraphNode* node, PixelAdjBit adj ) const;

        int32_t point2idx ( int32_t x, int32_t y ) const;
        vec2i    idx2point ( int32_t idx ) const;

        vec2i    nhb_point ( int32_t x, int32_t y, PixelAdjBit adj ) const;

        Vertex find_outgoing_boundaries ( PixelBoundaryNode* node, std::vector<std::unique_ptr<PixelBoundary>>& boundaries_out );

        int32_t heuristic_sparse ( int32_t x, int32_t y );
        int32_t heuristic_curves ( PixelGraphNode* edge_top, PixelGraphNode* edge_bottom );

        void save ( const std::string& addr ) const;
        void draw() const;
        void draw_boundaries() const;
        void draw_regions ( Vertex idx = -1 ) const;

        bool validate() const;

        // Region hierarchy
        vector<unique_ptr<PixelRegion>>             regions;
        vector<unique_ptr<PixelRegionBoundary>>     region_boundaries;
        vector<unique_ptr<PixelBoundary>>           boundaries;
        vector<unique_ptr<PixelBoundaryNode>>       boundary_nodes;
        unordered_map<int64_t, PixelBoundaryNode*> boundary_nodes_map;

        // Raster data
        unique_ptr<PixelGraphNode[]>    nodes = nullptr;
        int32_t                  width = 0;
        int32_t                  height = 0;
    };


    namespace {
        using namespace std;
        bool read_next_line ( ifstream& filestream, stringstream& ss ) {
            PF_ASSERT ( filestream.is_open() );
            string segment;
            bool flag;

            // Read an actual segment
            do {
                if ( filestream.eof() ) {
                    return false;
                }

                std::getline ( filestream, segment );
                ss.clear();
                ss.str ( segment );

                if ( segment.empty() || segment[0] == '#' ) {
                    flag = true;
                } else {
                    flag = false;
                }
            } while ( flag );

            return true;
        }

        template <typename T>
        T read_next_token ( stringstream& ss ) {
            T token;
            ss >> token;
            return token;
        }

        bool extract_name_value ( const string& segment, string& name, string& val ) {
            name.clear();
            val.clear();

            const char* r = segment.data();

            while ( *r != '\0' && *r == ' ' ) {
                ++r;
            }

            const char* name_beg = r;

            while ( *r != '\0' && *r != ' ' ) {
                ++r;
            }

            name = string ( name_beg, r - name_beg );
            val = r + 1;

            return !name.empty() && !val.empty();
        }
    }

    const PixelAdjBit PixelGraphNode::neighbors[8] = {
        PixelAdjBit::Top,
        PixelAdjBit::TopRight,
        PixelAdjBit::Right,
        PixelAdjBit::BottomRight,
        PixelAdjBit::Bottom,
        PixelAdjBit::BottomLeft,
        PixelAdjBit::Left,
        PixelAdjBit::TopLeft,
    };

    const PixelAdjBit PixelGraphNode::aa_neighbors[4] = {
        PixelAdjBit::Top,
        PixelAdjBit::Right,
        PixelAdjBit::Bottom,
        PixelAdjBit::Left
    };

    uint8_t& PixelGraphNode::y() {
        return yuv[0];
    }

    uint8_t& PixelGraphNode::u() {
        return yuv[1];
    }

    uint8_t& PixelGraphNode::v() {
        return yuv[2];
    }

    uint8_t PixelGraphNode::y() const {
        return yuv[0];
    }

    uint8_t PixelGraphNode::u() const {
        return yuv[1];
    }

    uint8_t PixelGraphNode::v() const {
        return yuv[2];
    }

    void PixelGraphNode::init_from_rgba ( int64_t _xy, uint8_t r, uint8_t g, uint8_t b, uint8_t _alpha ) {
        memset ( this, 0xffffffff, sizeof ( PixelGraphNode ) );
        y() = ( uint8_t ) ( 0.299 * r + 0.587 * g + 0.114 * b );
        //u() = ( uint8_t ) ( ( uint16_t ) ( -0.169 * r - 0.331 * g + 0.5 * b ) + 128 );
        //v() = ( uint8_t ) ( ( uint16_t ) ( 0.5 * r - 0.419 * g - 0.081 * b ) + 128 );
        u() = ( uint8_t ) ( ( uint16_t ) ( -0.14713 * r - 0.28886 * g + 0.436 * b ) + 128 );
        v() = ( uint8_t ) ( ( uint16_t ) ( 0.615 * r - 0.51499 * g - 0.10001 * b ) + 128 );
        alpha = _alpha;
        xy = _xy;
    }

    bool PixelGraphNode::connected ( PixelAdjBit adj_bit ) const {
        return adj & adj_bit;
    }

    void PixelGraphNode::disconnect ( int remove_adj ) {
        for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
            if ( remove_adj & PixelGraphNode::neighbors[inhb] ) {
                adj = adj & ~ ( PixelGraphNode::neighbors[inhb] );
            }
        }
    }

    int32_t PixelGraphNode::valence() const {
        int32_t valence = 0;

        for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
            valence += ( adj & PixelGraphNode::neighbors[inhb] ) ? 1 : 0;
        }

        return valence;
    }

    PixelAdjBit PixelGraphNode::next ( PixelAdjBit adj_ignore ) const {
        PF_ASSERT ( valence() <= 2 );

        for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
            PixelAdjBit adj_next = PixelGraphNode::neighbors[inhb];

            if ( adj & adj_next && adj_next != adj_ignore ) {
                return adj_next;
            }
        }

        return PixelAdjBit::None;
    }

    int64_t PixelGraphNode::tl() const {
        return xy;
    }

    int64_t PixelGraphNode::tr() const {
        return PixelGraph::pack ( PixelGraph::unpack ( xy ) + vec2i ( 1, 0 ) );
    }

    int64_t PixelGraphNode::br() const {
        return PixelGraph::pack ( PixelGraph::unpack ( xy ) + vec2i ( 1, 1 ) );
    }

    int64_t PixelGraphNode::bl() const {
        return PixelGraph::pack ( PixelGraph::unpack ( xy ) + vec2i ( 0, 1 ) );
    }

    vec2i PixelGraphNode::point() const {
        return PixelGraph::unpack ( xy );
    }

    // Compares the two sequences sequentially (unordered start) regardless of orientatation
    // -1 if the two paths have different lengths
    int PixelRegionBoundary::compare ( const PixelRegionBoundary& lhs, const PixelRegionBoundary& rhs ) {
        if ( lhs.points.size() != rhs.points.size() ) {
            return -1;
        }

        // Finding an identical point to start matching from
        Vertex n_points = lhs.points.size();

        for ( Vertex ipt_off = 0; ipt_off < n_points; ++ipt_off ) {
            if ( lhs.points[ipt_off] == rhs.points[0] ) {

                // Left to right or right to left ?
                if ( lhs.points[ipt_off + 1] == rhs.points[1] ) {
                    for ( Vertex ipt = 0; ipt < n_points; ++ipt ) {
                        if ( lhs.points[ipt_off + ipt] != rhs.points[ipt] ) {
                            return 1;
                        }
                    }
                } else if ( lhs.points[ipt_off - 1] == rhs.points[1] ) {
                    for ( Vertex ipt = 0; ipt < n_points; ++ipt ) {
                        if ( lhs.points[ipt_off - ipt] != rhs.points[ipt] ) {
                            return 1;
                        }
                    }
                } else {
                    return ( int ) n_points;
                }

                return 0;
            }
        }

        return ( int ) n_points;
    }

#define no_region 0xff

    Vertex PixelBoundaryNode::count_unique_regions() const {
        Vertex count = regions.size();

        for ( Vertex iregion = 0; iregion < regions.size(); ++iregion ) {
            for ( Vertex iregion_alt = iregion + 1; iregion_alt < regions.size(); ++iregion_alt ) {
                if ( regions[iregion] == regions[iregion_alt] ) {
                    PF_ASSERT ( regions[iregion]->index == regions[iregion_alt]->index );
                    --count;
                }
            }
        }

        return count;
    }

    vec2i PixelBoundaryNode::point() const {
        return PixelGraph::unpack ( xy );
    }

    void PixelRegion::init ( uint64_t* bits ) {
        uint8_t* b = ( uint8_t* ) bits;

        for ( int iregion = 0; iregion < sizeof ( uint64_t ) / sizeof ( uint8_t ); ++iregion ) {
            b[iregion] = no_region;
        }
    }

    bool PixelRegion::add ( uint64_t* bits, uint8_t region ) {
        uint8_t* b = ( uint8_t* ) bits;

        for ( int iregion = 0; iregion < sizeof ( uint64_t ) / sizeof ( uint8_t ); ++iregion ) {
            if ( b[iregion] == no_region ) {
                b[iregion] = region;
                return true;
            }
        }

        return false;
    }

    int PixelRegion::count ( uint64_t bits ) {
        int ret = 0;
        uint8_t* b = ( uint8_t* ) & bits;

        for ( int iregion = 0; iregion < sizeof ( uint64_t ) / sizeof ( uint8_t ); ++iregion ) {
            if ( b[iregion] == no_region ) {
                break;
            }

            ++ret;
        }

        return ret;
    }

    int PixelRegion::count_unique ( uint64_t bits ) {
        int ret = 0;
        uint8_t* b = ( uint8_t* ) & bits;

        for ( int iregion = 0; iregion < sizeof ( uint64_t ) / sizeof ( uint8_t ); ++iregion ) {
            if ( b[iregion] == no_region ) {
                break;
            }

            bool skip_count = false;

            for ( int iregion_alt = iregion - 1; iregion_alt >= 0; --iregion_alt ) {
                if ( b[iregion] == b[iregion_alt] ) {
                    skip_count = true;
                    break;
                }
            }

            ret += skip_count ? 0 : 1;
        }

        return ret;
    }

    int PixelRegion::at ( uint64_t bits, int idx ) {
        uint8_t* b = ( uint8_t* ) & bits;
        return ( int ) b[idx];
    }

    int64_t PixelGraph::pack ( vec2i xy ) {
        int64_t ret;
        ret = ( int64_t ) xy.x();
        ret |= ( ( int64_t ) xy.y() ) << 32;
        return ret;
    }

    vec2i PixelGraph::unpack ( int64_t xy ) {
        vec2i ret;
        ret.x() = ( int32_t ) ( xy & 0xffffffff );
        ret.y() = ( int32_t ) ( xy >> 32 );
        return ret;
    }

    vec2 PixelGraph::unpack_edge ( int64_t xy ) {
        return unpack ( xy ).cast<double>() - vec2 ( .5, .5 );
    }

    PixelAdjBit PixelGraph::opposite ( PixelAdjBit adj ) {
        switch ( adj ) {
        case PixelAdjBit::Top:
            return PixelAdjBit::Bottom;

        case PixelAdjBit::TopRight:
            return PixelAdjBit::BottomLeft;

        case PixelAdjBit::Right:
            return PixelAdjBit::Left;

        case PixelAdjBit::BottomRight:
            return PixelAdjBit::TopLeft;

        case PixelAdjBit::Bottom:
            return PixelAdjBit::Top;

        case PixelAdjBit::BottomLeft:
            return PixelAdjBit::TopRight;

        case PixelAdjBit::Left:
            return PixelAdjBit::Right;

        case PixelAdjBit::TopLeft:
            return PixelAdjBit::BottomRight;

        default:
            return PixelAdjBit::None;
        }
    }

    const PixelGraphNode* PixelGraph::at ( int32_t x, int32_t y ) const {
        return at ( vec2i ( x, y ) );
    }

    PixelGraphNode* PixelGraph::at ( int32_t x, int32_t y ) {
        return at ( vec2i ( x, y ) );
    }

    const PixelGraphNode* PixelGraph::at ( vec2i xy ) const {
        return const_cast<PixelGraph&> ( *this ).at ( xy );
    }

    PixelGraphNode* PixelGraph::at ( vec2i p ) {
        PF_ASSERT ( nodes );

        if ( p.x() < 0 || p.x() >= width ||
                p.y() < 0 || p.y() >= height ) {
            return nullptr;
        }

        return nodes.get() + ( p.y() * width + p.x() );
    }

    const PixelGraphNode* PixelGraph::at ( int64_t xy ) const {
        return at ( PixelGraph::unpack ( xy ) );
    }

    PixelGraphNode* PixelGraph::at ( int64_t xy ) {
        return at ( PixelGraph::unpack ( xy ) );
    }

    PixelGraphNode* PixelGraph::nhb ( PixelGraphNode* node, PixelAdjBit adj ) {
        PF_ASSERT ( node );

        int32_t x = PixelGraph::unpack ( node->xy ).x();
        int32_t y = PixelGraph::unpack ( node->xy ).y();

        if ( node->adj & adj ) {
            switch ( adj ) {
            case PixelAdjBit::Top:
                return at ( x, y - 1 );

            case PixelAdjBit::TopRight:
                return at ( x + 1, y - 1 );

            case PixelAdjBit::Right:
                return at ( x + 1, y );

            case PixelAdjBit::BottomRight:
                return at ( x + 1, y + 1 );

            case PixelAdjBit::Bottom:
                return at ( x, y + 1 );

            case PixelAdjBit::BottomLeft:
                return at ( x - 1, y + 1 );

            case PixelAdjBit::Left:
                return at ( x - 1, y );

            case PixelAdjBit::TopLeft:
                return at ( x - 1, y - 1 );

            default:
                PF_ABORT;
            }
        }

        return nullptr;
    }

    const PixelGraphNode* PixelGraph::nhb ( const PixelGraphNode* node, PixelAdjBit adj ) const {
        return const_cast<PixelGraph&> ( *this ).nhb ( const_cast<PixelGraphNode*> ( node ), adj );
    }

    int32_t PixelGraph::point2idx ( int32_t x, int32_t y ) const {
        return y * width + x;
    }

    vec2i PixelGraph::idx2point ( int32_t idx ) const {
        return vec2i ( idx % width, idx / width );
    }

    vec2i PixelGraph::nhb_point ( int32_t x, int32_t y, PixelAdjBit adj ) const {
        switch ( adj ) {
        case PixelAdjBit::Top:
            return vec2i ( x, y - 1 );

        case PixelAdjBit::TopRight:
            return vec2i ( x + 1, y - 1 );

        case PixelAdjBit::Right:
            return vec2i ( x + 1, y );

        case PixelAdjBit::BottomRight:
            return vec2i ( x + 1, y + 1 );

        case PixelAdjBit::Bottom:
            return vec2i ( x, y + 1 );

        case PixelAdjBit::BottomLeft:
            return vec2i ( x - 1, y + 1 );

        case PixelAdjBit::Left:
            return vec2i ( x - 1, y );

        case PixelAdjBit::TopLeft:
            return vec2i ( x - 1, y - 1 );

        default:
            PF_ABORT;
            return vec2i ( -1, -1 );
        }
    }

    // Auxiliary data structure for the unique boundary extraction given regions with multiple boundaries and
    // different ordering.
    struct PixelBoundaryRef {
        void add ( PixelRegion* region, PixelRegionBoundary* boundary, Vertex ipoint, Vertex direction ) {
            if ( point == index_null ) {
                point = boundary->points[ipoint];
            }

            if ( region_0 == nullptr && boundary_0 == nullptr && ipoint_0 == index_null && direction_0 == index_null ) {
                region_0 = region;
                boundary_0 = boundary;
                ipoint_0 = ipoint;
                direction_0 = direction;
            } else if ( region_1 == nullptr && boundary_1 == nullptr && ipoint_1 == index_null && direction_1 == index_null ) {
                region_1 = region;
                boundary_1 = boundary;
                ipoint_1 = ipoint;
                direction_1 = direction;
            } else {
                PF_ASSERT ( false );
            }
        }

        int64_t point = ( int64_t ) ~0;

        PixelRegion* region_0 = nullptr;
        PixelRegionBoundary* boundary_0 = nullptr;
        Vertex          ipoint_0 = ( Vertex ) ~0;
        Vertex          direction_0 = index_null;

        PixelRegion* region_1 = nullptr;
        PixelRegionBoundary* boundary_1 = nullptr;
        Vertex          ipoint_1 = index_null;
        Vertex          direction_1 = index_null;
    };

    // Given a valence 3/4 corner, it extracts a list of boundary pieces up to the next shared node
    Vertex PixelGraph::find_outgoing_boundaries ( PixelBoundaryNode* node, vector<unique_ptr<PixelBoundary>>& boundaries_out ) {
        PF_ASSERT ( node );

        Vertex n_pieces = node->regions.size();
        PF_ASSERT ( n_pieces > 2 );

        // Initializing references to outgoing boundaries
        PixelRegionBoundary* out_boundaries[4] = { nullptr, nullptr, nullptr, nullptr };

        // Since regions can have multiple boundaries, finding the specific boundary which overlaps the shared node
        for ( Vertex iregion = 0; iregion < n_pieces; ++iregion ) {
            PixelRegion* region = node->regions[iregion];

            for ( Vertex iboundary = 0; iboundary < region->boundaries.size(); ++iboundary ) {
                PixelRegionBoundary* boundary = region->boundaries[iboundary];
                Vertex ipoint = distance ( boundary->points.begin(), std::find ( boundary->points.begin(), boundary->points.begin(), node->xy ) );

                if ( ipoint != index_null ) {
                    out_boundaries[iregion] = boundary;
                    break;
                }
            }

            // If the shared node has those regions registered, there should be a point along any of the boundaries
            // of the region matching the node.
            PF_ASSERT ( out_boundaries[iregion] != nullptr );
        }

        int64_t          out_points[4] = { index_null };
        PixelBoundaryRef out_refs[4] = { PixelBoundaryRef() };
        Vertex            n_out_refs = 0;

        PF_VERBOSE_F( "Extracting boundaries from node %d %d", node->point().x(), node->point().y() );

        // Finding the pair of regions (boundaries) which specifically match the boundary
        for ( Vertex ipiece = 0; ipiece < n_pieces; ++ipiece ) {
            PixelRegion* region = node->regions[ipiece];
            PixelRegionBoundary* boundary = out_boundaries[ipiece];

            // Generally, we expect that a boundary has a single point touching the shared node, but it can happen
            // that two of the regions can be the same, that is why we need to try both
            Vertex ipoints[2] = { index_null, index_null };
            ipoints[0] = distance ( boundary->points.begin(), find ( boundary->points.begin(), boundary->points.end(), node->xy ) );
            ipoints[1] = distance ( boundary->points.begin(), find ( boundary->points.begin(), boundary->points.end(), node->xy ) );
            PF_ASSERT ( ipoints[0] != index_null );

            if ( ipoints[1] != index_null ) {
                // And let's avoid reprocessing the same boundary
                bool skip_piece = false;

                for ( Vertex ipiece_prev = ipiece - 1; ipiece_prev >= 0; --ipiece_prev ) {
                    if ( node->regions[ipiece_prev] == node->regions[ipiece] ) {
                        skip_piece = true;
                        break;
                    }
                }

                if ( skip_piece ) {
                    continue;
                }
            }

            Vertex n_ipoints = 1 + ( ( ipoints[1] != index_null ) ? 1 : 0 );

            for ( Vertex i = 0; i < n_ipoints; ++i ) {
                Vertex ipoint = ipoints[i];
                Vertex ipoint_prev = ipoint - 1;
                Vertex ipoint_next = ipoint + 1;
                int64_t point_prev = boundary->points[ipoint_prev];
                int64_t point_next = boundary->points[ipoint_next];

                PF_VERBOSE_F( "Region %d: %lld - %lld - %lld", region->index, ipoint_prev, ipoint, ipoint_next );

                // Previous point
                {
                    Vertex idx = -1;

                    for ( Vertex iout = 0; iout < n_out_refs; ++iout ) {
                        if ( out_points[iout] == point_prev ) {
                            idx = iout;
                            break;
                        }
                    }

                    if ( idx == -1 ) {
                        idx = n_out_refs++;
                        out_points[idx] = point_prev;
                    }

                    out_refs[idx].add ( region, boundary, ipoint_prev, -1 );
                }

                // Next point
                {
                    Vertex idx = -1;

                    for ( Vertex iout = 0; iout < n_out_refs; ++iout ) {
                        if ( out_points[iout] == point_next ) {
                            idx = iout;
                            break;
                        }
                    }

                    if ( idx == -1 ) {
                        idx = n_out_refs++;
                        out_points[idx] = point_next;
                    }

                    out_refs[idx].add ( region, boundary, ipoint_next, +1 );
                }
            }
        }

        // If the above loop worked correctly, all outgoing boundaries should be shared by two pixel regions
        for ( Vertex iout = 0; iout < n_out_refs; ++iout ) {
            PF_ASSERT ( out_refs[iout].boundary_0 != nullptr && out_refs[iout].ipoint_0 != index_null
                        && out_refs[iout].direction_0 != index_null );
            PF_ASSERT ( out_refs[iout].boundary_1 != nullptr && out_refs[iout].ipoint_1 != index_null
                        && out_refs[iout].direction_1 != index_null );
        }

        // Finally walking along each of the boundaries/directions obtained in the previous steps
        // up until the next shared node, which determines the end of the boundary piece
        for ( Vertex iout = 0; iout < n_out_refs; ++iout ) {
            const PixelBoundaryRef& ref = out_refs[iout];
            Vertex   inext_0 = ref.ipoint_0;
            Vertex   inext_1 = ref.ipoint_1;
            int64_t next = out_points[iout];

            PixelBoundary* boundary = new PixelBoundary;
            boundary->points.emplace_back ( node->xy );
            boundary->points.emplace_back ( next );

            while ( boundary_nodes_map.find ( next ) == boundary_nodes_map.end() ) {
                inext_0 += ref.direction_0;
                inext_1 += ref.direction_1;
                int64_t next_0 = ref.boundary_0->points[inext_0];
                int64_t next_1 = ref.boundary_1->points[inext_1];
                PF_ASSERT ( next_0 == next_1 ); // We are walking along the same unique boundary
                next = next_0;
                boundary->points.emplace_back ( next );
            }

            boundary->region_0 = ref.region_0;
            boundary->region_1 = ref.region_1;
            boundary->src = node;
            boundary->dst = boundary_nodes_map[next];
            boundaries_out.emplace_back ( unique_ptr<PixelBoundary> ( boundary ) );
        }

        return n_pieces;
    }

    // Counts the connected component in a 8x8 window around the edge
    int32_t PixelGraph::heuristic_sparse ( int32_t src_x, int32_t src_y ) {
        int32_t x_min = max ( src_x - 3, 0 );
        int32_t y_min = max ( src_y - 3, 0 );
        int32_t x_max = min ( src_x + 4, width - 1 );
        int32_t y_max = min ( src_y + 4, height - 1 );

        // todo: this can be done without any extra memory...
        static thread_local vector<int64_t>            traverse_stack;
        static thread_local unordered_set<int64_t>  traverse_visited;
        traverse_stack.clear();
        traverse_visited.clear();

        int32_t count = 0;
        int64_t src = pack ( { src_x, src_y } );
        traverse_stack.push_back ( src );
        traverse_visited.insert ( src );

        while ( !traverse_stack.empty() ) {
            vec2i xy = unpack ( traverse_stack.back() );
            PixelGraphNode* node = at ( xy );
            ++count;
            traverse_stack.pop_back();

            for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
                vec2i nhb_xy = nhb_point ( xy.x(), xy.y(), PixelGraphNode::neighbors[inhb] );
                int64_t nhb = pack ( nhb_xy );

                const bool connected = node->adj & PixelGraphNode::neighbors[inhb];
                const bool unvisited = traverse_visited.find ( nhb ) == traverse_visited.end();
                const bool in_window = nhb_xy.x() >= x_min && nhb_xy.x() <= x_max &&
                                       nhb_xy.y() >= y_min && nhb_xy.y() < y_max;

                if ( connected && unvisited && in_window ) {
                    traverse_stack.push_back ( nhb );
                    traverse_visited.insert ( nhb );
                }
            }
        }

        return count;
    }

    // Counts the length of the curve (valence 2) the two pixels are part of
    int32_t PixelGraph::heuristic_curves ( PixelGraphNode* edge_top, PixelGraphNode* edge_bottom ) {
        PixelGraphNode* nhb_bl = nhb ( edge_top, PixelAdjBit::BottomLeft );
        PixelGraphNode* nhb_br = nhb ( edge_top, PixelAdjBit::BottomRight );
        PF_ASSERT ( logic_xor ( nhb_bl == edge_bottom, nhb_br == edge_bottom ) );

        PixelAdjBit top_nhb_bit = nhb_bl == edge_bottom ? PixelAdjBit::BottomLeft : PixelAdjBit::BottomRight;
        PixelAdjBit bottom_nhb_bit = PixelGraph::opposite ( top_nhb_bit );

        int32_t length = 1; // The edge from src to bottom is not counted in the loops below
        // From top node
        {
            PixelAdjBit ignore_adj = top_nhb_bit;
            PixelGraphNode* next = edge_top;

            while ( next->valence() == 2 ) {
                ++length;
                PixelAdjBit adj_next = next->next ( ignore_adj );
                next = nhb ( next, adj_next );
                ignore_adj = PixelGraph::opposite ( adj_next );

                if ( next == edge_top ) {
                    break;
                }
            }
        }

        // From bottom node
        {
            PixelAdjBit ignore_adj = bottom_nhb_bit;
            PixelGraphNode* next = edge_bottom;

            while ( next->valence() == 2 ) {
                ++length;
                PixelAdjBit adj_next = next->next ( ignore_adj );
                next = nhb ( next, adj_next );
                ignore_adj = PixelGraph::opposite ( adj_next );

                if ( next == edge_bottom ) {
                    break;
                }
            }
        }

        PF_ASSERT ( length > 0 ); // min should be one if the above code worked correctly
        return length;
    }

    // Per-pixel storage for boundary extraction
    struct PointNhb {
        int64_t  next_1 = int64_max;
        int64_t  prev_1 = int64_max;
        int64_t  next_2 = int64_max;
        int64_t  prev_2 = int64_max;

        // Handling overlapping/partial insertions caused by unordered
        // processing of pixels
        void add ( int64_t prev, int64_t next ) {
            PF_ASSERT ( next != int64_max );

            if ( next == next_1 ) {
                prev_1 = prev_1 == int64_max ? prev : prev_1;
            } else if ( next == next_2 ) {
                prev_2 = prev_2 == int64_max ? prev : prev_2;
            } else if ( next_1 == int64_max ) {
                next_1 = next;
                prev_1 = prev;
            } else if ( next_2 == int64_max ) {
                next_2 = next;
                prev_2 = prev;
            } else {
                PF_ASSERT ( false );
            }
        }

        int64_t next ( int64_t current, int64_t prev ) const {
            PF_ASSERT ( prev != int64_max );

            if ( prev == prev_1 ) {
                return next_1;
            }

            if ( prev == prev_2 ) {
                return next_2;
            }

            PF_ASSERT ( next_1 != int64_max );

            if ( next_2 == int64_max ) {
                return next_1;
            }

            if ( next_1 == prev ) {
                return next_2;
            }

            if ( next_2 == prev ) {
                return next_1;
            }

            // Consistently ""guessing"" the next node using the same clockwise convention
            // todo: maybe handle this in the main adjacency code?
            int in_dx = PixelGraph::unpack ( current ).x() - PixelGraph::unpack ( prev ).x();
            int in_dy = PixelGraph::unpack ( current ).y() - PixelGraph::unpack ( prev ).y();

            // incoming direction
            int from = 0;

            if ( in_dx == 1 && in_dy == 0 ) { // Top
                from = 0;
            } else if ( in_dx == 0 && in_dy == 1 ) { // Right
                from = 1;
            } else if ( in_dx == -1 && in_dy == 0 ) { // Bottom
                from = 2;
            } else if ( in_dx == 0 && in_dy == -1 ) { // Left
                from = 3;
            } else {
                PF_ASSERT (
                    false ); // Is this happening because `vertex_first` calls next()? It can be fixed with a little thinking..
            }

            int out_dx_1, out_dy_1, out_dx_2, out_dy_2;
            out_dx_1 = PixelGraph::unpack ( next_1 ).x() - PixelGraph::unpack ( current ).x();
            out_dy_1 = PixelGraph::unpack ( next_1 ).y() - PixelGraph::unpack ( current ).y();
            out_dx_2 = PixelGraph::unpack ( next_2 ).x() - PixelGraph::unpack ( current ).x();
            out_dy_2 = PixelGraph::unpack ( next_2 ).y() - PixelGraph::unpack ( current ).y();

            // outgoing direction
            for ( int dir_off = 0; dir_off < 4; ++dir_off ) {
                int dir = ( from + dir_off ) % 4;

                if ( dir == 0 ) { // Top -> Right
                    if ( out_dx_1 == 0 && out_dy_1 == 1 ) {
                        return next_1;
                    }

                    if ( out_dx_2 == 0 && out_dy_2 == 1 ) {
                        return next_2;
                    }
                } else if ( dir == 1 ) { // Right -> Bottom
                    if ( out_dx_1 == -1 && out_dy_1 == 0 ) {
                        return next_1;
                    }

                    if ( out_dx_2 == -1 && out_dy_2 == 0 ) {
                        return next_2;
                    }
                } else if ( dir == 2 ) { // Bottom -> Left
                    if ( out_dx_1 == 0 && out_dy_1 == -1 ) {
                        return next_1;
                    }

                    if ( out_dx_2 == 0 && out_dy_2 == -1 ) {
                        return next_2;
                    }
                } else if ( dir == 3 ) { // Left -> Top
                    if ( out_dx_1 == 1 && out_dy_1 == 0 ) {
                        return next_1;
                    }

                    if ( out_dx_2 == 1 && out_dy_2 == 0 ) {
                        return next_2;
                    }
                }
            }

            PF_ASSERT ( false );
            return index_null;
        }

        bool is_shared() const {
            return next_1 != int64_max && next_2 != int64_max;
        }
    };

    void pixel_graph_init ( PixelGraph& G, const uint8_t* rgba_pixels, int32_t width, int32_t height, int channels ) {
        // Initializing pixel graph
        G.width = ( uint32_t ) width;
        G.height = ( uint32_t ) height;

        uint32_t n_pixels = ( uint32_t ) ( G.width * G.height );
        G.nodes = unique_ptr<PixelGraphNode[]> ( new PixelGraphNode[n_pixels] );

        // Creating graph nodes and converting to YUV
        for ( int32_t y = 0; y < G.height; ++y ) {
            for ( int32_t x = 0; x < G.width; ++x ) {
                const uint8_t* rgba = rgba_pixels + ( y * width + x ) * channels;
                G.at ( x, y )->init_from_rgba ( PixelGraph::pack ( vec2i ( x, y ) ), rgba[0], rgba[1], rgba[2], rgba[3] );

                int remove = 0x0;

                if ( x == 0 ) {
                    remove |= PixelAdjBit::TopLeft | PixelAdjBit::Left | PixelAdjBit::BottomLeft;
                }

                if ( y == 0 ) {
                    remove |= PixelAdjBit::TopLeft | PixelAdjBit::Top | PixelAdjBit::TopRight;
                }

                if ( x == G.width - 1 ) {
                    remove |= PixelAdjBit::TopRight | PixelAdjBit::Right | PixelAdjBit::BottomRight;
                }

                if ( y == G.height - 1 ) {
                    remove |= PixelAdjBit::BottomRight | PixelAdjBit::Bottom | PixelAdjBit::BottomLeft;
                }

                G.at ( x, y )->disconnect ( remove );
            }
        }

        // Disconnecting different pixels
        for ( int32_t y = 0; y < ( int32_t ) G.height; ++y ) {
            for ( int32_t x = 0; x < G.width; ++x ) {
                PixelGraphNode* node = G.at ( x, y );

                for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
                    PixelGraphNode* adj = G.nhb ( node, PixelGraphNode::neighbors[inhb] );

                    if ( !adj ) {
                        continue;
                    }

                    // Hqx thresholds https://web.archive.org/web/20070703061942/http://www.hiend3d.com/hq3x.html
                    // Inkscape's implementation uses different ones for uv (10) (9) ...
                    bool different_colors =
                        std::abs ( node->y() - adj->y() ) > 0x30 ||
                        std::abs ( node->u() - adj->u() ) > 7 ||
                        std::abs ( node->v() - adj->v() ) > 6 ||
                        node->alpha != adj->alpha;

                    if ( different_colors ) {
                        node->disconnect ( PixelGraphNode::neighbors[inhb] );
                    }
                }
            }
        }
    }

    void pixel_graph_remove_crossing_edges ( PixelGraph& G ) {
        // Pruning crossing edges using depixelize's pixel heuristics
        // Keeping the diagonal with the bigger number of connected components
        for ( int32_t y = 0; y < G.height; ++y ) {
            for ( int32_t x = 0; x < G.width; ++x ) {
                PixelGraphNode* tl = G.at ( x, y );

                if ( ( tl->adj & PixelAdjBit::BottomRight ) &&
                        ! ( tl->adj & PixelAdjBit::Right ) ) {
                    PixelGraphNode* tr = G.at ( x + 1, y );
                    PixelGraphNode* bl = G.at ( x, y + 1 );
                    PixelGraphNode* br = G.at ( x + 1, y + 1 );

                    if ( tr->adj & PixelAdjBit::BottomLeft ) {
                        int32_t weight_left2right = 0;
                        int32_t weight_right2left = 0;

                        // Sparse pixels
                        int32_t weight_left2right_sparse = G.heuristic_sparse ( x, y );
                        int32_t weight_right2left_sparse = G.heuristic_sparse ( x + 1, y );

                        if ( weight_left2right_sparse >= weight_right2left_sparse ) {
                            weight_right2left += weight_left2right_sparse - weight_right2left_sparse;
                        } else {
                            weight_left2right += weight_right2left_sparse - weight_left2right_sparse;
                        }

                        // Curves
                        int32_t weight_left2right_curves = 4 * G.heuristic_curves ( tl, br );
                        int32_t weight_right2left_curves = 4 * G.heuristic_curves ( tr, bl );

                        weight_left2right += weight_left2right_curves;
                        weight_right2left += weight_right2left_curves;

                        if ( weight_left2right <= weight_right2left ) {
                            tl->disconnect ( PixelAdjBit::BottomRight );
                            br->disconnect ( PixelAdjBit::TopLeft );
                            //draw::line(tl->point().cast<double>(), br->point().cast<double>(), Style::outline(colors::red, 1.));
                        }

                        if ( weight_left2right >= weight_right2left ) {
                            tr->disconnect ( PixelAdjBit::BottomLeft );
                            bl->disconnect ( PixelAdjBit::TopRight );
                            //draw::line(tr->point().cast<double>(), bl->point().cast<double>(), Style::outline(colors::red, 1.));
                        }
                    }
                }
            }
        }
    }

    // Finding boundary edges between the regions.For each new pixel encountered its region is identified by flooding.
    // The outgoing edges which have been removed are then put into a neighbor map storing the next (ccw) point for each edge.
    // The map is traversed starting from any point to obtain the boundary polyline
    // I can't use the int32_t as they are for the edge_nhb_map because edges in width + 1, height + 1
    // and the index wraps around.
    // The PointNhb map is technically not complete (`prev_` can miss) as we cannot guarantee the order
    // of insertion of the points. Nonetheless this is enough information to resolve the boundaries as
    // the `prev` nodes which are missing are for non-shared nodes and not required to extract the boundary
    void pixel_graph_flood_regions ( PixelGraph& G ) {
        std::unordered_map<int64_t, PointNhb> edge_nhb_map;
        std::unordered_set<int64_t>          edge_nhb_visit_map;
        std::unordered_set<int32_t>          flood_visit_map;
        std::vector<int32_t>                         flood_visit_stack;
        std::vector<std::vector<vec2>>                      region_boundaries;
        uint8_t                              region_index = 0; // Region index counter

        for ( int32_t floody = 0; floody < G.height; ++floody ) {
            for ( int32_t floodx = 0; floodx < G.width; ++floodx ) {
                int32_t floodidx = G.point2idx ( floodx, floody );

                // Pixel already processed
                if ( flood_visit_map.find ( floodidx ) != flood_visit_map.end() ) {
                    continue;
                }

                PF_ASSERT ( flood_visit_stack.empty() );
                flood_visit_stack.push_back ( floodidx );
                flood_visit_map.insert ( floodidx );
                int32_t pixels_in_region = 0;

                while ( !flood_visit_stack.empty() ) {
                    ++pixels_in_region;
                    int32_t idx = flood_visit_stack.back();
                    int32_t x = G.idx2point ( idx ).x();
                    int32_t y = G.idx2point ( idx ).y();
                    flood_visit_stack.pop_back();

                    PixelGraphNode* node = G.nodes.get() + idx;

                    // todo: the two cases (diagonal|aa) below can be put in their own loops, it's not as immediate to debug though
                    //
                    // First resolving diagonal neighbors for sparse pixel regions. If the diagonal
                    // dependency has been resolved we need to mask the bits to avoid processing the
                    // respective axis aligned bits. This increases the complexity of the neighbor map
                    // a point can have multiple neighbors depending on the incoming direction
                    int skip_mask = 0x0;

                    // Left to Right diagonal
                    if ( node->connected ( PixelAdjBit::BottomRight ) && !node->connected ( PixelAdjBit::Right )
                            && !node->connected ( PixelAdjBit::Bottom ) ) {
                        int64_t shared = node->br();
                        int64_t src_to_shared = node->tr();
                        int64_t shared_to_dst = PixelGraph::pack ( vec2i ( x + 2, y + 1 ) );

                        // Right edge of current pixel
                        edge_nhb_map[src_to_shared].add ( int64_max, shared );

                        // Top edge of bottomleft pixel
                        edge_nhb_map[shared].add ( src_to_shared, shared_to_dst );
                        skip_mask |= ( PixelAdjBit::Right | PixelAdjBit::Bottom );
                    }

                    if ( node->connected ( PixelAdjBit::TopLeft ) && !node->connected ( PixelAdjBit::Top )
                            && !node->connected ( PixelAdjBit::Left ) ) {
                        int64_t shared = node->tl();
                        int64_t src_to_shared = node->bl();
                        int64_t shared_to_dst = PixelGraph::pack ( vec2i ( x - 1, y ) );

                        // Left edge of current pixel
                        edge_nhb_map[src_to_shared].add ( int64_max, shared );

                        // Bottom edge of topleft pixel
                        edge_nhb_map[shared].add ( src_to_shared, shared_to_dst );
                        skip_mask |= ( PixelAdjBit::Top | PixelAdjBit::Left );
                    }

                    // Right to left diagonal
                    if ( node->connected ( PixelAdjBit::TopRight ) && !node->connected ( PixelAdjBit::Top )
                            && !node->connected ( PixelAdjBit::Right ) ) {
                        int64_t shared = node->tr();
                        int64_t src_to_shared = node->tl();
                        int64_t shared_to_dst = PixelGraph::pack ( vec2i ( x + 1, y - 1 ) );

                        // Top edge of current pixel
                        edge_nhb_map[src_to_shared].add ( int64_max, shared );

                        // Left edge of topright pixel
                        edge_nhb_map[shared].add ( src_to_shared, shared_to_dst );
                        skip_mask |= ( PixelAdjBit::Top | PixelAdjBit::Right );
                    }

                    if ( node->connected ( PixelAdjBit::BottomLeft ) && !node->connected ( PixelAdjBit::Bottom )
                            && !node->connected ( PixelAdjBit::Left ) ) {
                        int64_t shared = node->bl();
                        int64_t src_to_shared = node->br();
                        int64_t shared_to_dst = PixelGraph::pack ( vec2i ( x, y + 2 ) );

                        // Bottom edge of current pixel
                        edge_nhb_map[src_to_shared].add ( int64_max, shared );

                        // Right edge of bottomleft pixel
                        edge_nhb_map[shared].add ( src_to_shared, shared_to_dst );
                        skip_mask |= ( PixelAdjBit::Bottom | PixelAdjBit::Left );
                    }

                    // Axis-aligned adjacencies, the previous point is unknown
                    if ( !node->connected ( PixelAdjBit::Left ) && ! ( skip_mask & PixelAdjBit::Left ) ) {
                        int64_t src = PixelGraph::pack ( vec2i ( x, y + 1 ) );
                        int64_t dst = PixelGraph::pack ( vec2i ( x, y ) );
                        edge_nhb_map[src].add ( int64_max, dst );
                    }

                    if ( !node->connected ( PixelAdjBit::Bottom ) && ! ( skip_mask & PixelAdjBit::Bottom ) ) {
                        int64_t src = PixelGraph::pack ( vec2i ( x + 1, y + 1 ) );
                        int64_t dst = PixelGraph::pack ( vec2i ( x, y + 1 ) );
                        edge_nhb_map[src].add ( int64_max, dst );
                    }

                    if ( !node->connected ( PixelAdjBit::Right ) && ! ( skip_mask & PixelAdjBit::Right ) ) {
                        int64_t src = PixelGraph::pack ( vec2i ( x + 1, y ) );
                        int64_t dst = PixelGraph::pack ( vec2i ( x + 1, y + 1 ) );
                        edge_nhb_map[src].add ( int64_max, dst );
                    }

                    if ( !node->connected ( PixelAdjBit::Top ) && ! ( skip_mask & PixelAdjBit::Top ) ) {
                        int64_t src = PixelGraph::pack ( vec2i ( x, y ) );
                        int64_t dst = PixelGraph::pack ( vec2i ( x + 1, y ) );
                        edge_nhb_map[src].add ( int64_max, dst );
                    }

                    // Adding neighboring nodes
                    for ( int inhb = 0; inhb < array_len ( PixelGraphNode::neighbors ); ++inhb ) {
                        vec2i node_nhb_point = G.nhb_point ( x, y, PixelGraphNode::neighbors[inhb] );
                        int32_t node_nhb_idx = G.point2idx ( node_nhb_point.x(), node_nhb_point.y() );
                        PixelGraphNode* node_nhb = G.nhb ( node, PixelGraphNode::neighbors[inhb] );

                        if ( !node_nhb || flood_visit_map.find ( node_nhb_idx ) != flood_visit_map.end() ) {
                            continue;
                        }

                        flood_visit_map.insert ( node_nhb_idx );
                        flood_visit_stack.push_back ( node_nhb_idx );
                    }
                }

                PF_ASSERT ( !edge_nhb_map.empty() );
                edge_nhb_visit_map.clear();

                PixelRegion* region = new PixelRegion;
                region->index = region_index;
                std::memcpy ( &region->yuv, G.at ( floodx, floody )->yuv, 3 );
                region->alpha = G.at(floodx, floody)->alpha;

                // Getting the boundary points for this region from the edge neighbor map
                // any starting edge (except shared nodes) should work and bring us to the origin model.
                // This has to be done once for each region
                auto edge_it = edge_nhb_map.begin();

                for ( ; edge_it != edge_nhb_map.end(); ++edge_it ) {
                    if ( edge_nhb_visit_map.find ( edge_it->first ) != edge_nhb_visit_map.end() ) {
                        continue;
                    }

                    int64_t vertex_first = edge_it->first;
                    const PointNhb* nhb = &edge_it->second;

                    if ( nhb->is_shared() ) { // Starting from valence 2 nodes
                        continue;
                    }

                    int64_t vertex_prev = vertex_first;
                    int64_t vertex_next = nhb->next ( vertex_prev, vertex_first );

                    PixelRegionBoundary* boundary = new PixelRegionBoundary;

                    // Adding first point
                    boundary->points.push_back ( vertex_first );
                    edge_nhb_visit_map.insert ( vertex_first );

                    //dbg::info(FMT("Vertex first %d %d",
                    //PixelGraph::unpack(vertex_first).x(), PixelGraph::unpack(vertex_first).y()));

                    // Adding remaining points up to last to end
                    while ( vertex_next != vertex_first ) {
                        boundary->points.push_back ( vertex_next );
                        edge_nhb_visit_map.insert ( vertex_next );

                        vec2i xy = PixelGraph::unpack ( vertex_next );
                        //dbg::info(FMT("Vertex next %d %d",xy.x(), xy.y()));

                        nhb = &edge_nhb_map[vertex_next];
                        int64_t vertex_next_next = edge_nhb_map[vertex_next].next ( vertex_next, vertex_prev );
                        PF_ASSERT ( vertex_next_next != int64_max );

                        vertex_prev = vertex_next;
                        vertex_next = vertex_next_next;
                    }

                    // This is a leak -- fix it.
                    region->boundaries.push_back ( boundary );
                }

                G.regions.push_back ( unique_ptr<PixelRegion> ( region ) );
                edge_nhb_map.clear();

                PF_ASSERT ( region_index <= 254 ); // More than 255 regions? just need a small modification
                ++region_index;
            }
        }

        PF_VERBOSE_F( "Extracted %d regions", region_index );
    }

    namespace ImageSegment {
        void extract_closed_regions ( const IO::Image& I, std::vector<mat2x>& boundaries, std::vector<vec4>& colors ) {
            PF_ASSERT ( I.pixels && I.channels );
            boundaries.clear();

            PixelGraph G;
            pixel_graph_init ( G, I.pixels, I.width, I.height, I.channels );
            pixel_graph_remove_crossing_edges ( G );
            pixel_graph_flood_regions ( G );

            for ( int iregion = 0; iregion < G.regions.size(); ++iregion ) {
                PixelRegion* region = G.regions[iregion].get();
                PF_VERBOSE_F( "Region %d YUV %d %d %d", iregion, region->yuv[0], region->yuv[1], region->yuv[2] );

                for ( int iboundary = 0; iboundary < region->boundaries.size(); ++iboundary ) {
                    const int n_points = region->boundaries[iboundary]->points.size();
                    PF_VERBOSE_F ( "Boundary %d points %d", iboundary, n_points );

                    mat2x points ( 2, n_points );

                    for ( int p = 0; p < points.cols(); ++p ) {
                        points.col ( p ) = PixelGraph::unpack_edge ( region->boundaries[iboundary]->points[p] );
                    }

                    vec2i coord = PixelGraph::unpack ( region->boundaries[iboundary]->points[0] );
                    const double u = (double)(region->yuv[1] - 128);
                    const double v = (double)(region->yuv[2] - 128);

                    vec4 color;
                    color.x() = (1. * region->yuv[0] + 0. * u + 1.13983 * v) / 255;
                    color.y() = (1. * region->yuv[0] + -0.39465 * u + -0.58060 * v) / 255;
                    color.z() = (1. * region->yuv[0] + 2.03211 * u + 0. * v) / 255;
                    color.w() = (double)(region->alpha) / 255;

                    boundaries.emplace_back ( std::move ( points ) );
                    colors.emplace_back ( color );
                }
            }
        }

        void expand_and_cleanup ( IO::Image& I ) {
            PF_ASSERT ( I.channels == 4 && I.pixels );

            // expand by 1 pixel
            // ---
            int width_exp = I.width + 2;
            int height_exp = I.height + 2;
            uint8_t* pixels_exp = ( uint8_t* ) malloc ( I.channels * width_exp * height_exp );

            uint32_t bg = find_background_color ( I );

            // Initializing new image w/ top-left corner of old
            for ( int i = 0; i < width_exp * height_exp * I.channels; i += I.channels ) {
                std::memcpy ( pixels_exp + i, &bg, I.channels );
            }

            for ( int y = 0; y < I.height; ++y ) {
                uint8_t* dst = pixels_exp + I.channels * ( ( y + 1 ) * width_exp + 1 );
                uint8_t* src = I.pixels + I.channels * ( y * I.width );
                std::memcpy ( dst, src, I.channels * I.width );
            }

            // Setting all pixels with alpha 0 to the same color to avoid island regions
            for ( int ipixel = 0; ipixel < width_exp * height_exp * I.channels; ipixel += I.channels ) {
                uint8_t* rgba = pixels_exp + ipixel;

                if ( rgba[3] == 0 ) {
                    rgba[0] = rgba[1] = rgba[2] = 0;
                }
            }

            I.reset ( pixels_exp, width_exp, height_exp );

            // recolor tiny regions
            // ---
            std::unordered_map<uint32_t, Vertex> histogram;

            for ( int ipixel = 0; ipixel < I.width * I.height * I.channels; ipixel += I.channels ) {
                const uint8_t* rgba = I.pixels + ipixel;
                uint32_t key = pack ( rgba );

                if ( histogram.find ( key ) == histogram.end() ) {
                    histogram[key] = 0;
                }

                ++histogram[key];
            }

            PF_VERBOSE_F ( "Color histogram entries %zu", histogram.size() );

            // Finding new color mappings for small regions
            const Vertex  small_component_thr = ( Vertex ) round ( ( double ) max ( I.width, I.height ) / 10 );
            const uint8_t channel_similarity_thr = 64;
            int n_quantized_colors = 0;

            std::unordered_map<uint32_t, uint32_t> quantized_colors_map;

            for ( auto histogram_it : histogram ) {
                uint32_t rgba_packed = histogram_it.first;
                int count = histogram_it.second;

                if ( count < small_component_thr ) {
                    vec4i rgba = unpack ( rgba_packed );
                    uint32_t rgba_alt_best = 0x0;
                    int alt_count = 0;

                    // Empirically determining whether it's an image feature we should keep by checking
                    // how close it is to other color
                    for ( auto histogram_alt_it = histogram.begin(); histogram_alt_it != histogram.end(); ++histogram_alt_it ) {
                        uint32_t rgba_packed_alt = histogram_alt_it->first;
                        int count_alt = histogram_alt_it->second;

                        if ( count_alt < small_component_thr ) {
                            continue;
                        }

                        const vec4i rgba_alt = unpack ( rgba_packed_alt );
                        vec4i rgba_delta = ( rgba - rgba_alt );

                        if ( rgba_delta.cwiseAbs().maxCoeff() <= channel_similarity_thr &&
                                histogram_alt_it->second > alt_count ) {
                            alt_count = histogram_alt_it->second;
                            rgba_alt_best = rgba_packed_alt;
                        }
                    }

                    if ( alt_count == 0 ) { // No replacement
                        quantized_colors_map[rgba_packed] = 0x0;
                        continue;
                    }

                    quantized_colors_map[rgba_packed] = rgba_alt_best;
                    ++n_quantized_colors;
                    continue;
                }

                quantized_colors_map[rgba_packed] = rgba_packed;
            }

            // Recoloring image
            for ( int ipixel = 0; ipixel < I.width * I.height * I.channels; ipixel += I.channels ) {
                uint8_t* rgba = I.pixels + ipixel;
                uint32_t rgba_packed = pack ( rgba );
                uint32_t new_rgba_packed = quantized_colors_map[rgba_packed];

                unpack ( new_rgba_packed, rgba[0], rgba[1], rgba[2], rgba[3] );
            }

            PF_VERBOSE_F ( "Color reduction %zu -> %lld", histogram.size(), histogram.size() - n_quantized_colors );
        }

        uint32_t find_background_color ( IO::Image& I ) {
            unordered_map<uint32_t, Vertex> histogram;
            histogram.reserve ( I.width * 2 + I.height * 2 );

            // Top/Bottom sides
            for ( int x = 0; x < I.width; ++x ) {
                const uint8_t* rgba_top = I.pixels + I.channels * x;
                const uint8_t* rgba_bottom = I.pixels + I.channels * ( ( I.height - 1 ) * I.width + x );
                const uint32_t key_top = pack ( rgba_top );
                const uint32_t key_bottom = pack ( rgba_bottom );

                if ( histogram.find ( key_top ) == histogram.end() ) {
                    histogram[key_top] = 0;
                }

                if ( histogram.find ( key_bottom ) == histogram.end() ) {
                    histogram[key_bottom] = 0;
                }

                ++histogram[key_top];
                ++histogram[key_bottom];
            }

            // Left/Right sides
            for ( int y = 0; y < I.height; ++y ) {
                const uint8_t* rgba_left = I.pixels + I.channels * ( I.width * y );
                const uint8_t* rgba_right = I.pixels + I.channels * ( I.width * y + ( I.width - 1 ) );
                const uint32_t key_left = pack ( rgba_left );
                const uint32_t key_right = pack ( rgba_right );

                if ( histogram.find ( key_left ) == histogram.end() ) {
                    histogram[key_left] = 0;
                }

                if ( histogram.find ( key_right ) == histogram.end() ) {
                    histogram[key_right] = 0;
                }

                ++histogram[key_left];
                ++histogram[key_right];
            }

            int    max_count = 0;
            uint32_t bg_rgba = 0;

            for ( auto histogram_it : histogram ) {
                if ( histogram_it.second > max_count ) {
                    max_count = histogram_it.second;
                    bg_rgba = histogram_it.first;
                }
            }

            PF_ASSERT ( max_count > 0 );
            return bg_rgba;
        }

        int find_binary_color_region ( const std::vector<mat2x>& boundaries ) {
            int iregion = -1;
            int iregion_corners = 0;

            for ( int i = 0; i < boundaries.size(); ++i ) {
                std::vector<int> convexity;
                PathUtils::compute_convexities ( boundaries[i], convexity );

                // skipping squared contour
                int corners = 0;;

                for ( int j = 0; j < boundaries[i].cols(); ++j ) {
                    if ( convexity[j] != 0 ) {
                        ++corners;
                    }
                }

                PF_VERBOSE_F ( "boundary %d corners %d", i, corners );

                if ( corners != 4 && corners > iregion_corners ) {
                    iregion = i;
                    iregion_corners = corners;
                }
            }

            // are we actually rasterizing a square?
            if (iregion == -1) {
                return boundaries.size() - 1;
            }

            PF_ASSERT ( iregion >= 0 );
            return iregion;
        }
    }
}