// Copied from minimesh.
// https://gitlab.com/hooshi/minimesh-public
//
// mesh_connectivity.hpp
//
// Author: Shayan Hoshyari
//

#ifndef MINIMESH_MESH_CONNECTIVITY_IS_INCLUDED
#define MINIMESH_MESH_CONNECTIVITY_IS_INCLUDED

//
#include <stack>
#include <vector>

//
#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)


//
// Surface_mesh_connectivity
//
// A class for storing and traversing a mesh using the half-edge
// data structure.
//

class Mesh_connectivity
{

public:
  // ========================================================
  // The connectivity data
  // ========================================================

  constexpr static int hole_index = -2;
  constexpr static int invalid_index = -3;

  //
  // Vertex_data: The connectivity data for a vertex.
  //
  class Vertex_data
  {
  public:
    // Index of a half-edge going out of this vertex.
    int half_edge = invalid_index;

  private:
    // Is this vertex present in the mesh, or is it deleted, i.e.,
    // it is only an empty space in the vector of vertices data.
    //
    // You should not access this member directly. Use
    // Vertex_iterator:: is_active() and deactivate()
    // Note that once you deactivate a vertex, you cannot activate
    // it again. You can just add new vertices.
    bool is_active = false;

    // Please do not create this class yourself.
    // Use Surface_mesh_connectivity::add_vertex()
    // And then access data via Vertex_iterator::data()
    Vertex_data() = default;
    friend class Mesh_connectivity;
  };

  //
  // Half_edge_data: Connectivity data for a half-edge
  //
  struct Half_edge_data
  {
    // The index of the half-edge after
    int next = invalid_index;
    
    // The index of the half-edge before
    int prev = invalid_index;

    // The index of the twin half edge
    int twin = invalid_index;

    // The index of the face, that the half edge belongs to
    int face = invalid_index;

    // The index of the vertex that this half edge goes out of.
    int origin = invalid_index;


  private:
    // Is this half-edge present in the mesh, or is it deleted, i.e.,
    // it is only an empty space in the vector of half-edges data.
    //
    // You should not access this member directly. Use
    // Half_edge_iterator:: is_active() and deactivate()
    // Note that once you deactivate a Half_edge, you cannot activate
    // it again. You can just add new half edges.
    bool is_active = false;

    // Please do not create this class yourself.
    // Use Surface_mesh_connectivity::add_half_edge()
    // And then access data via Face_iterator::data()
    Half_edge_data() = default;
    friend class Mesh_connectivity;
  };

  //
  // Face_data: Connectivity data for a face
  //
  struct Face_data
  {

    // The index of a half-edge belonging to this face.
    int half_edge = invalid_index;

  private:
    // Is this face present in the mesh, or is it deleted, i.e.,
    // it is only an empty space in the vector of faces data.
    //
    // You should not access this member directly. Use
    // Face_iterator:: is_active() and deactivate()
    // Note that once you deactivate a face, you cannot activate
    // it again. You can just add new faces.
    bool is_active = false;

    // Please do not create this class yourself.
    // Use Surface_mesh_connectivity::add_face()
    // And then access data via Face_iterator::data()
    Face_data() = default;
    friend class Mesh_connectivity;
  };


  // ========================================================
  // Half edge iterators for traversing the mesh
  //  and removing entities from the mesh
  // ========================================================

  // Forward declarations are necessary here.
  class Vertex_iterator;
  class Half_edge_iterator;
  class Face_iterator;

  //
  // Vertex_iterator: An iterator to traverse mesh vertices
  //
  class Vertex_iterator
  {
  public:
    // Returns a reference to the data of this vertex.
    //  The data object can be used to modify connectivity.
    Vertex_data & data();

    // Get the index of the vertex.
    int index();

    // Is this vertex active or deleted
    bool is_active();

    // Deletes the vertex
    void deactivate();

    // Returns an iterator to a half edge going *out* of this vertex.
    Half_edge_iterator half_edge();

    bool is_equal(Vertex_iterator);

    Vertex_iterator() = default;

  private:
    // The mesh that the iterator belongs to
    Mesh_connectivity * _parent;

    // Index of the vertex
    int _index;

    // clang-format off
    Vertex_iterator(const int index, Mesh_connectivity *parent): _parent(parent), _index(index){}
    // clang-format on

    // Don't create this class yourself. Use
    // Surface_mesh_connectivity::vertex_at()
    friend class Mesh_connectivity;
  };

  //
  // Vertex_iterator: An iterator to traverse mesh half edges
  //
  class Half_edge_iterator
  {
  public:
    // Returns a reference to the data of this half edge.
    //  The data object can be used to modify connectivity.
    Half_edge_data & data();

    // Index of half edge
    int index();

    // Is half edge active or deleted
    bool is_active();

    // Delete the half edge
    void deactivate();

    // Iterator to the next half edge
    Half_edge_iterator next();

    // Iterator to the previous half edge
    Half_edge_iterator prev();

    // Iterator to the twin half edge
    Half_edge_iterator twin();

    // Iterator to the vertex at the end of half edge
    Vertex_iterator dest();

    // Iterator to the vertex at the beginning of half edge
    Vertex_iterator origin();

    // The face that the half edge belongs to
    Face_iterator face();


    bool is_equal(Half_edge_iterator);

    Half_edge_iterator() = default;

  private:
    // The mesh that the iterator belongs to
    Mesh_connectivity * _parent;

    // Index of the half edge
    int _index;

    // clang-format off
    Half_edge_iterator(const int index, Mesh_connectivity *parent): _parent(parent), _index(index){}
    // clang-format on

    // Don't create this class yourself. Use
    // Surface_mesh_connectivity::half_edge_at()
    friend class Mesh_connectivity;
  };

  //
  // Face_iterator: An iterator to traverse mesh faces
  //
  class Face_iterator
  {
  public:
    // Returns a reference to the data of this half edge.
    //  The data object can be used to modify connectivity.
    Face_data & data();

    // Index of the face
    int index();

    // Iterator to a half edge on the face
    Half_edge_iterator half_edge();

    // get number of vertices (or edges)
    int n_vertices();

    // Is the face that this iterator referes to active or deleted
    bool is_active();

    // Delete the face that this iterator represents
    void deactivate();

    // Are two iterators equal
    bool is_equal(Face_iterator);

    Face_iterator() = default;

  private:
    // The mesh that the iterator belongs to
    Mesh_connectivity * _parent;

    // Index of the face
    int _index;

    // clang-format off
    Face_iterator(const int index, Mesh_connectivity *parent): _parent(parent), _index(index){}
    // clang-format on

    // Don't create this class yourself. Use
    // Surface_mesh_connectivity::half_edge_at()
    friend class Mesh_connectivity;
  };

  //
  // Vertex_ring_iterator: An iterator to traverse the half edges
  //  and vertices around a vertex.
  //

  // Broken for Kopf Junctions
  
  class Vertex_ring_iterator
  {
  public:
    // An iterator to one of the half edges surrounding the vertex
    // that *points to* the vertex of interest
    Half_edge_iterator half_edge();

    // Attempts to see if the vertex is on the boundary. If so returns
    // true, and false otherwise.
    // In the event of success puts the half_edge to one on the boundary
    // that points to the vertex of interest.
    bool reset_boundary();

    // Updates half_edge to present the next half-edge surrounding
    // the vertex
    // Returns true if the half edge is not seen before
    // when all the half edges are traversed, returns false.
    bool advance();

    Vertex_ring_iterator() = default;

  private:
    // Current half edge being traversed
    Half_edge_iterator _cur;

    // The first half edge that was seen
    Half_edge_iterator _end;

    // Half_edge is a half edge that points to the vertex of interets
    Vertex_ring_iterator(Half_edge_iterator half_edge);

    // Don't create this class yourself. Use
    // Surface_mesh_connectivity::half_edge_at()
    friend class Mesh_connectivity;
  };

  // Returns an iterator for the vertex with id of vertex_id
  Vertex_iterator vertex_at(const int vertex_id);

  // Returns a ring iterator for vertex with vertex_id
  // The ring iterator can be used to traverse all half edges
  // surrounding vertex_id
  Vertex_ring_iterator vertex_ring_at(const int vertex_id);

  // Returns an iterator for the face with id of face_id
  Face_iterator face_at(const int face_id);

  // Returns an iterator for the half_edge with id of he_id
  Half_edge_iterator half_edge_at(const int he_id);

  // A dummy face to represent a hole
  Face_iterator hole();


  // ========================================================
  // Number of current entities
  // ========================================================

  // Returns
  // # of active vertices = # of total vertices - # of deleted vertices
  int n_active_vertices();

  // Returns
  // # of active half edges = # of total vertices - # of deleted half edges
  int n_active_half_edges();

  // Returns
  // # of active faces = # of total vertices - # of deleted faces
  int n_active_faces();

  // Returns number of total vertices (including the deleted ones)
  int n_total_vertices();

  // Returns number of total half edges (including the deleted ones)
  int n_total_half_edges();

  // Returns number of total faces (including the deleted ones)
  int n_total_faces();

  // ========================================================
  // Adding entities to the mesh
  // ========================================================

  // Adds a vertex to the mesh, and returns an iterator to it.
  // if allow_recycling is true, then it tries to revive a deleted vertex
  // otherwise, it adds a totally new vertex and does not touch the dead ones.
  Vertex_iterator add_vertex(const bool allow_recycling = true);

  // Adds a face to the mesh, and returns an iterator to it.
  // if allow_recycling is true, then it tries to revive a deleted face,
  // otherwise, it adds a totally new face and does not touch the dead ones.
  Face_iterator add_face(const bool allow_recycling = true);

  // Adds a half edge to the mesh, and returns an iterator to it.
  // if allow_recycling is true, then it tries to revive a deleted half edge,
  // otherwise, it adds a totally new half edge and does not touch the dead ones.
  Half_edge_iterator add_half_edge(const bool allow_recycling = true);

  // ========================================================
  // Defragmenting the mesh
  // ========================================================

  struct Defragmentation_maps
  {
    std::vector<int> old2new_faces;
    std::vector<int> new2old_faces;
    //
    std::vector<int> old2new_half_edges;
    std::vector<int> new2old_half_edges;
    //
    std::vector<int> old2new_vertices;
    std::vector<int> new2old_vertices;
  };
  void compute_defragmention_maps(Defragmentation_maps &);
  void defragment(Defragmentation_maps & m, Mesh_connectivity & out);
  void defragment_in_place(Defragmentation_maps & m);

  // ========================================================
  // Build the mesh from a given connectivity
  // ========================================================

  // Default constructor
  // does not do anything
  Mesh_connectivity() = default;

  // Build the connectivity from a list of vertices,
  // and vertices of triangles
  //
  // Input:
  //
  //   const int n_vertices -> number of points (can include inactive ones)
  //
  //   triangles_vertices, list of vertex number belonging to triangles in the
  //   order triangle0_v0, triangle0_v1, triangle0_v2, ..., trianglen_v0, trianglen_v1, trianglen_v2
  //   such that the number of triangles if triangle_verts.size()/3.
  void build_from_triangles(const int n_vertices, const std::vector<int> & triangle_verts);

  // Build the connectivity from a list of vertices,
  // and vertices of polygons with arbitrary number of vertices
  //
  //   const int n_vertices -> number of points (can include inactive ones)
  //
  //   polygon_verts, list of vertex number belonging to polygons in the
  //   order polygon0_v0, polygon0_v1, ..., polygon0_vend, ..., polygonn_v0, polygonn_v1, ..., polygonn_vend
  //
  //   polygon_adj, a list polygon_adj.size() == number of polygons +1
  //   It should be organized such that polygon_verts[ polygon_adj[i] ], ... ,polygon_verts[ polygon_adj[i+1]-1 ] are
  //   the vertices of polygon i.
  void build_from_polygons(const int n_vertices,
      const std::vector<int> & polygon_verts,
      const std::vector<int> & polygon_adj);



  // ========================================================
  // Check if the mesh is not corrupted
  // ========================================================
  bool check_sanity_slowly(bool verbose=false);

  // ========================================================
  // Copying the mesh
  // ========================================================

  // copy another mesh
  void copy(const Mesh_connectivity &);

  // swap the data with another mesh
  void swap(Mesh_connectivity &);

  // clear the mesh data
  void clear();

private:
  // List of all vertices (deleted or alive)
  std::vector<Vertex_data> _vertices;
  // A stack to deleted vertices, so that we can revive
  // one in O(1).
  std::stack<int> _inactive_vertices;

  // List of all half edges (deleted or alive)
  std::vector<Half_edge_data> _half_edges;
  // A stack to deleted half edges, so that we can revive
  // one in O(1).
  std::stack<int> _inactive_half_edges;

  // List of all face (deleted or alive)
  std::vector<Face_data> _faces;
  // A stack to deleted faces, so that we can revive
  // one in O(1).
  std::stack<int> _inactive_faces;

  // Shayan: Let them be for now.
  //
  // To prevent accidental copying of the mesh (which is a considerable
  // amount of data, the = operator is disables
  // Mesh_connectivity & operator=(const Mesh_connectivity &) = delete;
  // Mesh_connectivity(const Mesh_connectivity &) = delete;
};

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)



#endif /* MINIMESH_MESH_CONNECTIVITY_IS_INCLUDED */
