#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <glm/glm.hpp>  // OpenGL Mathematics : http://glm.g-truc.net/0.9.5/index.html
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>

#include <Mesh.h>

using namespace glm;
using namespace std;


int main()
{
    string input_file_name = "/home/henry/Documents/Qt_Project/dernierTPGeonum/data/cube_quad.off";
    
    cout << "Reading file '" << input_file_name << "'..." << endl;
    Mesh m(input_file_name.c_str());
    cout << "... done" << endl;
    
    cout << "Mesh extracted." << endl;
    cout << "Number of vertices : " << m.vertices.size() << endl;
    for(unsigned int i = 0; i < m.vertices.size(); i++)
    {
        vec3 p = m.get_vertex(i);
        cout << i << " : " << p[0] << " " << p[1] << " " << p[2] << endl;
    }
    
    cout << "Number of faces : " << m.faces.size() << endl;
    for(unsigned int i = 0; i < m.faces.size(); i++)
    {
        vector<unsigned int> f = m.get_face(i);
        cout << i << " :";

        for(unsigned int j = 0; j < f.size(); j++)
        {
            cout << " " << f[j];
        }
        cout << endl;
    }
    
    
    cout << endl << endl;
    
    
    vector< vector< unsigned int > > vert_neib = m.get_neighborhoods();
    cout << "Vertices->vertices adjacency: " << endl;
    for(unsigned int i = 0; i < vert_neib.size(); i++)
    {
        cout << i << " : ";
        for(unsigned int j = 0; j < vert_neib[i].size(); j++)
        {
            cout << vert_neib[i][j] << ", ";
        }
        cout << endl;
    }


    cout << endl << endl;


    vector< vector< unsigned int > > vert_faces =  m.get_vertex_faces();
    cout << "Vertices->faces adjacency: " << endl;
    for(unsigned int i = 0; i < vert_faces.size(); i++)
    {
        cout << i << " : ";
        for(unsigned int j = 0; j < vert_faces[i].size(); j++)
        {
            cout << vert_faces[i][j] << ", ";
        }
        cout << endl;
    }
    
    
    vector< Edge > edges = m.get_edges();
    cout << "Edges: " << endl;
    for(unsigned int i = 0; i < edges.size(); i++)
    {
        cout << i << " : " << edges[i].m_i0 << " - " << edges[i].m_i1 << endl;
    }
    
    
    cout << endl << endl;
    
    
    vector< vector< unsigned int > > edge_faces = m.get_edge_faces(edges);
    cout << "Edges->faces adjacency: " << endl;
    for(unsigned int i = 0; i < edge_faces.size(); i++)
    {
        cout << i << " : ";
        for(unsigned int j = 0; j < edge_faces[i].size(); j++)
        {
            cout << edge_faces[i][j] << ", ";
        }
        cout << endl;
    }
    
    
    
    cout << endl << endl;
    

    vector< vector< unsigned int > > vert_edges = m.get_vertex_edges(vert_neib, edges);
    cout << "Vertices->edges adjacency: " << endl;
    for(unsigned int i = 0; i < vert_edges.size(); i++)
    {
        cout << i << " : ";
        for(unsigned int j = 0; j < vert_edges[i].size(); j++)
        {
            cout << vert_edges[i][j] << ", ";
        }
        cout << endl;
    }
    
    
    cout << endl << endl;


    cout << endl << endl;

    Mesh m_subdiv = m;
  
    string output_file_name = "/home/henry/Documents/Qt_Project/dernierTPGeonum/output/cube_subdiv0.obj";
    
    cout << "Writing file '" << output_file_name << "'..." << endl;
    m_subdiv.write_obj(output_file_name.c_str());
    cout << "... done" << endl;


    for(int i = 1; i <= 3; i++)
    {
    
        cout << endl << endl;

        cout << "Computing Subdivision #" << i << "..." << endl;
        
        m_subdiv = m_subdiv.subdivide();
        
        cout << "... done" << endl;
        
        
        cout << endl << endl;
        
        
        stringstream ss;
        ss << i;
        output_file_name = "/home/henry/Documents/Qt_Project/dernierTPGeonum/output/cube_subdiv" + ss.str() + ".obj";
        
        cout << "Writing file '" << output_file_name << "'..." << endl;
        m_subdiv.write_obj(output_file_name.c_str());
        cout << "... done" << endl;
    }
    

	return EXIT_SUCCESS;
}

