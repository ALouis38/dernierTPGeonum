#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <algorithm>

#include <glm/glm.hpp>

using namespace glm;
using namespace std;


class Edge
{

public:

    Edge(unsigned int i0, unsigned int i1): m_i0(std::min(i0, i1)), m_i1(std::max(i0, i1)){}    // Direct ctor, stores sorted attributes
    
    unsigned int m_i0, m_i1;    // Indices of the extremities of the edge
    
    bool operator ==(const Edge& e) const { return m_i0 == e.m_i0 && m_i1 == e.m_i1; }  // equality operator
};


class Mesh
{

public:

  // Ctors
  Mesh(){}                      // Empty constructor
  Mesh(const char* filename);   // Constructs a mesh from an OFF file

  // Data
  vector< vec3 > vertices;      // Vertices positions
  vector< vector< unsigned int > > faces;   // Polygon faces indice list

  // Attribute acces
  vector< unsigned int > get_face(const unsigned int i) const;  // Returns a position
  vec3 get_vertex(const unsigned int i) const;  // Return an indice list defining a face

  // Output
  void write(const char* filename) const;       // Outputs an OFF file
  void write_obj(const char* filename) const;   // Outputs an OBJ file
  
  
  // Neighborhood computation
  vector< vector< unsigned int > > get_neighborhoods() const;   // Computes adjacent vertices of each vertex of the mesh (ordered)
  vector< Edge > get_edges() const;     // Computes edges of the mesh
  
  vector< vector< unsigned int > > get_vertex_faces() const;    // Computes adjacent faces of each vertex of the mesh
  vector< vector< unsigned int > > get_edge_faces(vector< Edge > edges) const;      // Computes adjacent faces of each edge
  
  vector< vector< unsigned int > > get_vertex_edges(vector< Edge > edges) const;    // Computes adjacent edges of each vertex
  vector< vector< unsigned int > > get_vertex_edges(vector< vector< unsigned int > > vert_neib, vector< Edge > edges) const;    // Computes adjacent edges of each vertex (ordered)
  
  unsigned int get_oriented_face(unsigned int i_e0, unsigned int i_e1) const;   // Computes the face containing the oriented edge i_e0 -> i_e1
  unsigned int get_next_vert(unsigned int i_e0, unsigned int i_e1) const;       // Computes the vertex following the sequence i_e0, i_e1 in the corresponding face

  
  // Subdivision
  Mesh subdivide() const;   // Computes the Catmull-Clark subdivision of the mesh

};




#endif // MESH_H
