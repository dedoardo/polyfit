// polyvec
#include <polyvec/mc/mesh_connectivity.hpp>
#include <polyvec/mc/bit_vector.hpp>
#include <polyvec/utils/matrix.hpp>

// libc++
#include <iostream>
#include <map>
#include <utility>
#include <sstream>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)


// ========================================================
//     STUFF NEEDED TO MAKE THE CODE BUILD
// ========================================================

NAMESPACE_BEGIN(consolecolors)
static inline const char * const black() { return  "\033[0;30m";}
static inline const char * const darkblue () { return "\033[0;34m";}
static inline const char * const darkgreen() { return "\033[0;32m";}
static inline const char * const darkteal() { return "\033[0;36m";}
static inline const char * const darkred() { return "\033[0;31m";}
static inline const char * const darkpink() { return "\033[0;35m";}
static inline const char * const darkyellow() { return "\033[0;33m";}
static inline const char * const gray() { return "\033[0;37m";}
static inline const char * const darkgray() { return "\033[1;30m";}
static inline const char * const blue() { return "\033[1;34m";}
static inline const char * const green() { return "\033[1;32m";}
static inline const char * const teal() { return "\033[1;36m";}
static inline const char * const red() { return "\033[1;31m";}
static inline const char * const pink() { return "\033[1;35m";}
static inline const char * const yellow() { return "\033[1;33m";}
static inline const char * const white() { return "\033[1;37m";}
static inline const char * const reset() { return "\033[0m";}
NAMESPACE_END(consolecolors)

#define MINIMESH_UNUSED(x) (void)(x)

#define MINIMESH_STR(X) static_cast<std::ostringstream&>(std::ostringstream().flush() << X).str()

inline void force_assert_checkpoint()
{
  const double x = 101010101;
  MINIMESH_UNUSED(x);
}

#define force_assert_msg(EXPR, MSG)                                                                               \
  if(!(EXPR))                                                                                                     \
  {                                                                                                               \
    std::ostringstream ss;                                                                                        \
    ss << std::scientific;                                                                                        \
    printf("Assertion in %s, %d \n", __func__, __LINE__);                                                         \
    ss.str("");                                                                                                   \
    ss << #EXPR << "\n";                                                                                          \
    ss << MSG;                                                                                                    \
    printf("%s %s %s\n", consolecolors::red(), ss.str().c_str(), consolecolors::reset()); \
    force_assert_checkpoint();                                                                                    \
    throw int(0);                                                                                                 \
  }

#define force_assert(EXPR) force_assert_msg(EXPR, "")

// ========================================================
//                      ITERATORS
// ========================================================

// ===================== VERTEX ===========================

Mesh_connectivity::Vertex_data &
Mesh_connectivity::Vertex_iterator::data()
{
  force_assert(_index >= 0);
  force_assert(_index < _parent->n_total_vertices());

  return _parent->_vertices[index()];
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Vertex_iterator::half_edge()
{
  return _parent->half_edge_at(data().half_edge);
}

int
Mesh_connectivity::Vertex_iterator::index()
{
  return _index;
}

bool
Mesh_connectivity::Vertex_iterator::is_active()
{
  return data().is_active;
}

void
Mesh_connectivity::Vertex_iterator::deactivate()
{
  force_assert(is_active() && "Attempting to delete already deleted vertex");

  data().is_active = false;
  _parent->_inactive_vertices.push(index());
}

bool
Mesh_connectivity::Vertex_iterator::is_equal(Vertex_iterator o)
{
  return (_parent == o._parent) && (_index == o._index);
}

// ========================= HALF EDGE ===========================
Mesh_connectivity::Half_edge_data &
Mesh_connectivity::Half_edge_iterator::data()
{
  force_assert(_index >= 0);
  force_assert(_index < _parent->n_total_half_edges());

  return _parent->_half_edges[index()];
}

int
Mesh_connectivity::Half_edge_iterator::index()
{
  return _index;
}

bool
Mesh_connectivity::Half_edge_iterator::is_active()
{
  return data().is_active;
}

void
Mesh_connectivity::Half_edge_iterator::deactivate()
{
  force_assert(is_active() && "Attempting to delete already deleted half edge");

  data().is_active = false;
  _parent->_inactive_half_edges.push(index());
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Half_edge_iterator::next()
{
  return _parent->half_edge_at(data().next);
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Half_edge_iterator::prev()
{
  return _parent->half_edge_at(data().prev);
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Half_edge_iterator::twin()
{
  return _parent->half_edge_at(data().twin);
}

Mesh_connectivity::Vertex_iterator
Mesh_connectivity::Half_edge_iterator::dest()
{
  return twin().origin();
}

Mesh_connectivity::Vertex_iterator
Mesh_connectivity::Half_edge_iterator::origin()
{
  return _parent->vertex_at(data().origin);
}

Mesh_connectivity::Face_iterator
Mesh_connectivity::Half_edge_iterator::face()
{
  return _parent->face_at(data().face);
}

bool
Mesh_connectivity::Half_edge_iterator::is_equal(Half_edge_iterator o)
{
  return (_parent == o._parent) && (_index == o._index);
}

// =========================== FACE ==============================

Mesh_connectivity::Face_data &
Mesh_connectivity::Face_iterator::data()
{
  force_assert(_index >= 0);
  force_assert(_index < _parent->n_total_faces());

  return _parent->_faces[_index];
}

int
Mesh_connectivity::Face_iterator::index()
{
  return _index;
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Face_iterator::half_edge()
{
  return _parent->half_edge_at(data().half_edge);
}

int
Mesh_connectivity::Face_iterator::n_vertices()
{
  force_assert(this->is_active());

  if(this->is_equal(_parent->hole()))
  {
    return 0;
  }
  else
  {
    Mesh_connectivity::Half_edge_iterator he = this->half_edge();
    Mesh_connectivity::Half_edge_iterator he_end = this->half_edge();
    int answer = 0;
    do
    {
      ++answer;
      he = he.next();
    } while(!he.is_equal(he_end));
    return answer;
  }
}

bool
Mesh_connectivity::Face_iterator::is_active()
{
  return data().is_active;
}

void
Mesh_connectivity::Face_iterator::deactivate()
{
  force_assert(is_active() && "Attempting to delete already deleted face");

  data().is_active = false;
  _parent->_inactive_faces.push(index());
}

bool
Mesh_connectivity::Face_iterator::is_equal(Face_iterator other)
{
  return (other._parent == this->_parent) && (other._index == this->_index);
}

// ===================== RING ITERATOR ===========================

Mesh_connectivity::Vertex_ring_iterator::Vertex_ring_iterator(Half_edge_iterator half_edge)
: _cur(half_edge)
, _end(half_edge)
{
  force_assert(_cur.is_active());
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::Vertex_ring_iterator::half_edge()
{
  return _cur;
}


bool
Mesh_connectivity::Vertex_ring_iterator::reset_boundary()
{
  Half_edge_iterator snapshot = _cur;
  do
  {
    // start from the boundary edge located to the right of the vertex.
    if(_cur.face().is_equal(_cur._parent->hole()) && _cur.next().face().is_equal(_cur._parent->hole()))
    {
      _end = _cur;
      return true;
    }
    _cur = _cur.next().twin();
  } while(!_cur.is_equal(snapshot));
  return false;
}

bool
Mesh_connectivity::Vertex_ring_iterator::advance()
{
  _cur = _cur.next().twin();
  return (!_cur.is_equal(_end));
}


// ===================== CREATING ITERATORS ===========================

Mesh_connectivity::Vertex_iterator
Mesh_connectivity::vertex_at(const int vertex_id)
{
  force_assert(vertex_id < n_total_vertices());
  return Vertex_iterator(vertex_id, this);
}

Mesh_connectivity::Vertex_ring_iterator
Mesh_connectivity::vertex_ring_at(const int vertex_id)
{
  force_assert(vertex_id < n_total_vertices());
  // cannot traverse around a dead vertex
  force_assert(vertex_at(vertex_id).is_active());
  // we must initialize the iterator with a half edge ending
  // at this vertex.
  return Vertex_ring_iterator(vertex_at(vertex_id).half_edge().twin());
}

Mesh_connectivity::Face_iterator
Mesh_connectivity::face_at(const int face_id)
{
  if(face_id == hole_index)
  {
    return hole();
  }
  else
  {
    force_assert(face_id < n_total_faces());
    return Face_iterator(face_id, this);
  }
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::half_edge_at(const int he_id)
{
  force_assert(he_id < n_total_half_edges());
  if( he_id != invalid_index )
    return Half_edge_iterator(he_id, this);
  else 
    return Half_edge_iterator(invalid_index, this);
}


Mesh_connectivity::Face_iterator
Mesh_connectivity::hole()
{
  return Face_iterator(hole_index, this);
}

// ========================================================
// Number of current entities
// ========================================================

int
Mesh_connectivity::n_active_vertices()
{
  int ans = (int)_vertices.size() - (int)_inactive_vertices.size();
  force_assert(ans >= 0);
  return ans;
}

int
Mesh_connectivity::n_active_half_edges()
{
  int ans = (int)_half_edges.size() - (int)_inactive_half_edges.size();
  force_assert(ans >= 0);
  return ans;
}

int
Mesh_connectivity::n_active_faces()
{
  int ans = (int)_faces.size() - (int)_inactive_faces.size();
  force_assert(ans >= 0);
  return ans;
}

int
Mesh_connectivity::n_total_vertices()
{
  return (int)_vertices.size();
}

int
Mesh_connectivity::n_total_half_edges()
{
  return (int)_half_edges.size();
}

int
Mesh_connectivity::n_total_faces()
{
  return (int)_faces.size();
}

// ========================================================
//              Adding entities to the mesh
// ========================================================

Mesh_connectivity::Vertex_iterator
Mesh_connectivity::add_vertex(const bool allow_recycling)
{
  int answer_id;
  if(allow_recycling && _inactive_vertices.size())
  {
    answer_id = _inactive_vertices.top();
    _inactive_vertices.pop();
  }
  else
  {
    answer_id = (int)_vertices.size();
    _vertices.push_back(Vertex_data());
  }
  Vertex_iterator answer = vertex_at(answer_id);
  answer.data().is_active = true;
  return answer;
}

Mesh_connectivity::Half_edge_iterator
Mesh_connectivity::add_half_edge(const bool allow_recycling)
{
  int answer_id;
  if(allow_recycling && _inactive_half_edges.size())
  {
    answer_id = _inactive_half_edges.top();
    _inactive_half_edges.pop();
  }
  else
  {
    answer_id = (int)_half_edges.size();
    _half_edges.push_back(Half_edge_data());
  }
  Half_edge_iterator answer = half_edge_at(answer_id);
  answer.data().is_active = true;
  return answer;
}

Mesh_connectivity::Face_iterator
Mesh_connectivity::add_face(const bool allow_recycling)
{
  int answer_id;
  if(allow_recycling && _inactive_faces.size())
  {
    answer_id = _inactive_faces.top();
    _inactive_faces.pop();
  }
  else
  {
    answer_id = (int)_faces.size();
    _faces.push_back(Face_data());
  }
  Face_iterator answer = face_at(answer_id);
  answer.data().is_active = true;
  return answer;
}

// ========================================================
//         Continuous indices map and defragmenting
// ========================================================

void
Mesh_connectivity::compute_defragmention_maps(Defragmentation_maps & m)
{
  m.old2new_faces.resize(n_total_faces());
  m.old2new_half_edges.resize(n_total_half_edges());
  m.old2new_vertices.resize(n_total_vertices());
  //
  m.new2old_faces.resize(n_active_faces());
  m.new2old_half_edges.resize(n_active_half_edges());
  m.new2old_vertices.resize(n_active_vertices());

  int newid;

  //
  // Continuous indices for vertices
  //
  newid = 0;
  for(int voldid = 0; voldid < n_total_vertices(); ++voldid)
  {
    Vertex_iterator vertex = vertex_at(voldid);
    if(vertex.is_active())
    {
      m.old2new_vertices[voldid] = newid;
      m.new2old_vertices[newid] = voldid;
      ++newid;
    }
    else
    {
      m.old2new_vertices[voldid] = invalid_index;
    }
  }

  //
  // Continuous indices for half edges
  //
  newid = 0;
  for(int heoldid = 0; heoldid < n_total_half_edges(); ++heoldid)
  {
    Half_edge_iterator he = half_edge_at(heoldid);
    if(he.is_active())
    {
      m.old2new_half_edges[heoldid] = newid;
      m.new2old_half_edges[newid] = heoldid;
      ++newid;
    }
    else
    {
      m.old2new_half_edges[heoldid] = invalid_index;
    }
  }

  //
  // Continuous indices for faces
  //
  newid = 0;
  for(int foldid = 0; foldid < n_total_faces(); ++foldid)
  {
    Face_iterator face = face_at(foldid);
    if(face.is_active())
    {
      m.old2new_faces[foldid] = newid;
      m.new2old_faces[newid] = foldid;
      ++newid;
    }
    else
    {
      m.old2new_faces[foldid] = invalid_index;
    }
  }
}

void
Mesh_connectivity::defragment(Defragmentation_maps & m, Mesh_connectivity & defragmented_mesh)
{
  compute_defragmention_maps(m);

  defragmented_mesh.clear();
  defragmented_mesh._faces.reserve(n_active_faces());
  defragmented_mesh._half_edges.reserve(n_active_half_edges());
  defragmented_mesh._vertices.reserve(n_active_vertices());

  // Copy vertices
  for(int vnewid = 0; vnewid < n_active_vertices(); ++vnewid)
  {
    const int voldid = m.new2old_vertices[vnewid];

    Vertex_iterator vertex_new = defragmented_mesh.add_vertex();
    Vertex_iterator vertex_old = this->vertex_at(voldid);
    vertex_new.data().half_edge = m.old2new_half_edges[vertex_old.data().half_edge];
  }

  // Copy faces
  for(int fnewid = 0; fnewid < n_active_faces(); ++fnewid)
  {
    const int foldid = m.new2old_faces[fnewid];

    Face_iterator face_new = defragmented_mesh.add_face();
    Face_iterator face_old = this->face_at(foldid);
    face_new.data().half_edge = m.old2new_half_edges[face_old.data().half_edge];
  }

  // Copy half edges
  for(int henewid = 0; henewid < n_active_half_edges(); ++henewid)
  {
    const int heoldid = m.new2old_half_edges[henewid];
    Half_edge_iterator he_new = defragmented_mesh.add_half_edge();
    Half_edge_iterator he_old = this->half_edge_at(heoldid);

    he_new.data().face = (he_old.data().face == hole_index ? hole_index : m.old2new_faces[he_old.data().face]);
    he_new.data().next = m.old2new_half_edges[he_old.data().next];
    he_new.data().prev = m.old2new_half_edges[he_old.data().prev];
    he_new.data().twin = m.old2new_half_edges[he_old.data().twin];
    he_new.data().origin = m.old2new_vertices[he_old.data().origin];
  }
}

void
Mesh_connectivity::defragment_in_place(Defragmentation_maps & m)
{
  Mesh_connectivity defragmented_mesh;
  defragment(m, defragmented_mesh);
  this->swap(defragmented_mesh);
}

// ========================================================
//              Building mesh connectivity
// ========================================================

void
Mesh_connectivity::build_from_triangles(const int n_vertices, const std::vector<int> & triangle_verts)
{
  force_assert_msg((int)triangle_verts.size() % 3 == 0, "Triangle connectivity list is corrupted");
  const int n_triangles = (int)triangle_verts.size() / 3;

  std::vector<int> polygon_adj(n_triangles + 1);
  polygon_adj[0] = 0;
  for(int i = 0; i < n_triangles; ++i)
  {
    polygon_adj[i + 1] = polygon_adj[i] + 3;
  }

  this->build_from_polygons(n_vertices, triangle_verts, polygon_adj);
}

void
Mesh_connectivity::build_from_polygons(const int n_vertices,
    const std::vector<int> & polygon_verts,
    const std::vector<int> & polygon_adj)
{

  // Clear current data
  this->clear();

  force_assert_msg(polygon_adj.back() == (int)polygon_verts.size(), "polygon vertices corrupted");
  const int n_polygons = (int)polygon_adj.size() - 1;

  //
  // Add the vertices
  //
  _vertices.reserve(n_vertices);
  for(int vid = 0; vid < n_vertices; vid++)
  {
    Vertex_iterator vertex = add_vertex();
    force_assert(vertex.index() == vid);
  }

  //
  // Add the half edges and build a map to keep track of the
  //  connectivity
  //
  typedef std::map<std::pair<int, int>, int> Edge_map;
  auto get_origin = [](Edge_map::const_iterator in) -> int { return in->first.first; };
  auto get_dest = [](Edge_map::const_iterator in) -> int { return in->first.second; };
  auto get_half_edge = [](Edge_map::const_iterator in) -> int { return in->second; };
  auto create_edge_map_member = [](int beg, int end, int he) -> Edge_map::value_type {
    return std::make_pair(std::make_pair(beg, end), he);
  };
  MINIMESH_UNUSED(get_origin); 
  Edge_map edge_map;

  _faces.reserve(n_polygons);
  _half_edges.reserve(polygon_adj.back()); // this won't be enough if the mesh is open
  for(int polyid = 0; polyid < n_polygons; ++polyid)
  {
    int n_polygon_vertices = polygon_adj[polyid + 1] - polygon_adj[polyid];
    std::vector<Half_edge_iterator> half_edges(n_polygon_vertices);
    std::vector<Vertex_iterator> vertices(n_polygon_vertices);

    // create half edges and the faces. Also get handles to vertices
    Face_iterator face = add_face();

    for(int voffset = 0; voffset < n_polygon_vertices; ++voffset)
    {
      half_edges[voffset] = add_half_edge();
      vertices[voffset] = vertex_at(polygon_verts[polygon_adj[polyid] + voffset]);
    }

    //  set the data as much as you can
    for(int voffset = 0; voffset < n_polygon_vertices; ++voffset)
    {
      int voffsetp1 = (voffset + 1) % n_polygon_vertices;
      int voffsetm1 = (voffset - 1 + n_polygon_vertices) % n_polygon_vertices;

      half_edges[voffset].data().next = half_edges[voffsetp1].index();
      half_edges[voffset].data().prev = half_edges[voffsetm1].index();
      half_edges[voffset].data().face = face.index();
      half_edges[voffset].data().origin = vertices[voffset].index();

      vertices[voffset].data().half_edge = half_edges[voffset].index();
    }
    face.data().half_edge = half_edges[0].index();

    //
    // Now try to set the twin for the half edge
    //
    for(int voffset = 0; voffset < n_polygon_vertices; ++voffset)
    {
      int voffsetp1 = (voffset + 1) % n_polygon_vertices;

      // try to find your twin
      Edge_map::iterator it = edge_map.find(std::make_pair(vertices[voffsetp1].index(), vertices[voffset].index()));


      if((it != edge_map.end()))
      {
        force_assert(get_origin(it) == vertices[voffsetp1].index());
        force_assert(get_dest(it) == vertices[voffset].index());
        half_edges[voffset].data().twin = get_half_edge(it);
        half_edge_at(get_half_edge(it)).data().twin = half_edges[voffset].index();
        edge_map.erase(it);
      }
      else
      {
        std::pair<Edge_map::iterator, bool> scss = edge_map.insert(create_edge_map_member(
            vertices[voffset].index(), vertices[voffsetp1].index(), half_edges[voffset].index()));
        force_assert_msg(scss.second == true,
            " Mesh non manifold, half edge duplicated, verts:  " << vertices[voffset].index() << " "
                                                                 << vertices[voffsetp1].index());
      }
    }

  } // End of polygons

  force_assert(n_polygons == n_total_faces());

  //
  // Deativate unreferenced vertices
  //
  for(int vid = 0; vid < n_total_vertices(); ++vid)
  {
    Mesh_connectivity::Vertex_iterator viter = vertex_at(vid);
    if(viter.data().half_edge == invalid_index)
      viter.deactivate();
  }

  //
  // Add half edges for the holes
  //
  for(Edge_map::iterator it = edge_map.begin(); it != edge_map.end();)
  {
    Half_edge_iterator he_interior = half_edge_at(get_half_edge(it));
    Half_edge_iterator he_bdry = add_half_edge();
    he_bdry.data().face = hole_index;
    he_bdry.data().twin = he_interior.index();
    he_bdry.data().origin = get_dest(it);
    he_interior.data().twin = he_bdry.index();

    it = edge_map.erase(it);
  } // End of residue iterators
  force_assert_msg(edge_map.size() == 0, "[Error] Corrupted mesh, n_mems left: " << edge_map.size());

  //
  // Set he->next for hole half edges
  //
  for(int heid = 0; heid < n_total_half_edges(); heid++)
  {
    Half_edge_iterator he = half_edge_at(heid);

    if(!he.face().is_equal(hole()))
    {
      force_assert(he.data().next != invalid_index);
      force_assert(he.data().prev != invalid_index);
    }
    else
    {
      force_assert(he.data().next == invalid_index);
      Half_edge_iterator he_next = he.twin();
      while(!he_next.face().is_equal(hole()))
      {
        he_next = he_next.prev().twin();
      }
      he.data().next = he_next.index();
      he_next.data().prev = he.index();
    }

  } // End of hald edges
} // All done with connectivity

// ========================================================
//              Check sanity slowly
// ========================================================

bool
Mesh_connectivity::check_sanity_slowly(bool verbose)
{
  //
  // Let's not use force_assert
  //
  auto checkpoint = []() {
    const int x = 12302;
    MINIMESH_UNUSED(x);
  };

#define soft_assert_msg(CND, MSG)                      \
  if(!(CND))                                           \
  {                                                    \
    if(verbose){                                       \
    std::cout << consolecolors::red();                 \
    std::cout << "Mesh insane: " << #CND << std::endl; \
    std::cout << MSG << std::endl;                     \
    std::cout << consolecolors::reset();               \
    }                                                  \
    checkpoint();                                      \
    return false;                                      \
  }

  //
  // First check that each face.he.face == face
  //
  for(int fid = 0; fid < n_total_faces(); ++fid)
  {
    Face_iterator face = face_at(fid);
    if(face.is_active())
    {
      Half_edge_iterator he = face.half_edge();
      soft_assert_msg(he.is_active(), "");
      soft_assert_msg(face.is_equal(he.face()), "face.he.face == face failed");
    }
  }

  //
  // Now check that the half edges of each face construct a loop
  // Also the half-edges should be active
  //
  const int max_face_half_edges = 100;
  for(int fid = 0; fid < n_total_faces(); ++fid)
  {
    int num_face_half_edges = 0;
    Face_iterator face = face_at(fid);
    if(face.is_active())
    {
      Half_edge_iterator he = face.half_edge();
      Half_edge_iterator he_end = face.half_edge();
      do
      {
        ++num_face_half_edges;
        soft_assert_msg(he.is_active(), "active face has inactive half edge");
        soft_assert_msg(face.is_equal(he.face()), "active face has inactive half edge");
        soft_assert_msg(num_face_half_edges < max_face_half_edges, "face has too many half edges");
      } while(!he.is_equal(he_end));
    }
  }


  //
  // Check mutual relationship between half edges
  //
  for(int heid = 0; heid < n_total_half_edges(); ++heid)
  {
    Half_edge_iterator he = half_edge_at(heid);
    if(he.is_active())
    {
      soft_assert_msg(he.twin().is_active(), "");
      soft_assert_msg(he.twin().twin().is_equal(he), "");
      //
      soft_assert_msg(he.next().is_active(), "");
      soft_assert_msg(he.next().prev().is_equal(he), "");
      //
      soft_assert_msg(he.prev().is_active(), "");
      soft_assert_msg(he.prev().next().is_equal(he), "");
      //
      soft_assert_msg(!he.prev().is_equal(he), "");
      soft_assert_msg(!he.next().is_equal(he), "");
      soft_assert_msg(!he.prev().is_equal(he.next()), "");
      //
      soft_assert_msg(he.origin().is_active(), "");
      soft_assert_msg(he.face().is_equal(hole()) ||  he.face().is_active(), "");
      // Both the he and twin should not be on the boundary (a tiny line!)
      soft_assert_msg( (!he.face().is_equal(hole())) || (!he.twin().face().is_equal(hole())), "");

    }
  }

  //
  // Vertices
  //
  const int max_vertex_outgoing_half_edges = 100;
  Bit_vector is_vertex_visited(n_total_vertices());

  for(int vid = 0; vid < n_total_vertices(); ++vid)
  {
    Vertex_iterator viter = vertex_at(vid);

    if(viter.is_active())
    {
      is_vertex_visited.reset_all();
      int n_outgoing_boundary_half_edges = 0;
      int n_incoming_boundary_half_edges = 0;
      int n_traversed_edges = 0;

      // The first one is force because it means will crash the code
      force_assert(viter.data().half_edge != invalid_index);
      //
      soft_assert_msg(viter.half_edge().is_active(), "");

      //
      Vertex_ring_iterator ring_iter = vertex_ring_at(vid);
      do
      {      
        const int other_vertex = ring_iter.half_edge().origin().index();
        soft_assert_msg(ring_iter.half_edge().is_active(), "");
        soft_assert_msg(ring_iter.half_edge().twin().is_active(), "");
        soft_assert_msg(ring_iter.half_edge().dest().index() == vid, "");

        // The number of outgoing half-edges should not be unbounded
        ++n_traversed_edges;
        soft_assert_msg(n_traversed_edges < max_vertex_outgoing_half_edges, "");

        // There should not be more than one edge between two vertices
        soft_assert_msg(!is_vertex_visited[other_vertex], "More than one edge between two");
        is_vertex_visited.set( other_vertex );

        // Not more than one boundary half-edge can go out of the vertex
        if(ring_iter.half_edge().twin().face().is_equal( hole() ) ) ++n_outgoing_boundary_half_edges;
        if(ring_iter.half_edge().face().is_equal( hole() ) ) ++n_incoming_boundary_half_edges;
        soft_assert_msg(n_outgoing_boundary_half_edges <= 1, "");
        soft_assert_msg(n_outgoing_boundary_half_edges <= 1, "");

      } while(ring_iter.advance());
    } // end of active vertex
  } // End of vertices

  return true;
#undef soft_assert_msg
}

// ========================================================
//              COPYING THE MESH
// ========================================================

void
Mesh_connectivity::copy(const Mesh_connectivity & other)
{
  this->_vertices = other._vertices;
  this->_faces = other._faces;
  this->_half_edges = other._half_edges;

  this->_inactive_vertices = other._inactive_vertices;
  this->_inactive_faces = other._inactive_faces;
  this->_inactive_half_edges = other._inactive_half_edges;
}

void
Mesh_connectivity::swap(Mesh_connectivity & other)
{
  this->_vertices.swap(other._vertices);
  this->_faces.swap(other._faces);
  this->_half_edges.swap(other._half_edges);

  this->_inactive_vertices.swap(other._inactive_vertices);
  this->_inactive_faces.swap(other._inactive_faces);
  this->_inactive_half_edges.swap(other._inactive_half_edges);
}

void
Mesh_connectivity::clear()
{
  this->_vertices.clear();
  this->_faces.clear();
  this->_half_edges.clear();

  this->_inactive_vertices = std::stack<int>();
  this->_inactive_faces = std::stack<int>();
  this->_inactive_half_edges = std::stack<int>();
}

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

// *INDENT-ON*