#include <cstdlib>
#include <iostream>
#include <math.h>

#include <Mesh.h>

using namespace glm;
using namespace std;


Mesh::Mesh(const char* filename) 
{
	int j = 0;
    unsigned int tmp;
	FILE *file;
    int   error;

    if((file=fopen(filename,"r"))==NULL)
	{
		std::cerr << "Unable to read : " << filename << std::endl;
	}

	// create mesh
    vertices = vector<vec3>();
    faces    = vector<vector< unsigned int > >();

    unsigned int nb_vertices, nb_faces;

	error = fscanf(file,"OFF\n%d %d %d\n",&(nb_vertices),&(nb_faces),&tmp);
	if(error==EOF) 
	{
		std::cerr << "Unable to read : " << filename << std::endl;
	}

    vertices.resize(nb_vertices);
    faces.resize(nb_faces);

	// reading vertices
	for(int i=0;i<nb_vertices;++i) 
	{
		error = fscanf(file,"%f %f %f\n",&(vertices[i][0]),&(vertices[i][1]),&(vertices[i][2]));
		if(error==EOF) 
		{
			std::cerr << "Unable to read vertices of : " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	// reading faces
	j = 0;
    for(int i = 0; i < nb_faces; ++i)
	{
        error = fscanf(file,"%d ",&tmp);
		
		if(error==EOF) 
		{
			std::cerr << "Unable to read faces of : " << filename << std::endl;
			exit(EXIT_FAILURE);
        }

        unsigned int N = tmp;

        for(unsigned int j = 0; j < N; j++)
        {
            error = fscanf(file,"%d ",&tmp);

            if(error==EOF)
            {
                std::cerr << "Unable to read faces of : " << filename << std::endl;
                exit(EXIT_FAILURE);
            }

            faces[i].push_back(tmp);
        }
	}

	fclose(file);
}


vec3 Mesh::get_vertex(const unsigned int i) const
{
    return vertices[i];
}


vector< unsigned int > Mesh::get_face(const unsigned int i) const
{
    return faces[i];
}


void Mesh::write(const char* filename) const 
{
	FILE *file;

    // file opening
	if((file=fopen(filename,"w"))==NULL) 
	{
		std::cerr << "Unable to open : " << filename << std::endl;
	}
	
	
	// header
	fprintf(file, "OFF\n");
    fprintf(file,"%lu %lu\n", vertices.size(), faces.size());
	
	// vertices positions
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        vec3 p = get_vertex(i);
        fprintf(file,"%f %f %f\n", p[0], p[1], p[2]);
    }
    
    // faces indices
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);
        fprintf(file,"3 %i %i %i\n", f[0], f[1], f[2]);
    }
}


void Mesh::write_obj(const char* filename) const 
{
	FILE *file;

    // file opening
	if((file=fopen(filename,"w"))==NULL) 
	{
		std::cerr << "Unable to open : " << filename << std::endl;
	}
	
	
	// vertices positions
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        vec3 p = get_vertex(i);
        fprintf(file,"v %f %f %f\n", p[0], p[1], p[2]);
    }
    
    // faces indices
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);
        fprintf(file,"f");

        for(unsigned int j = 0; j < f.size(); j++)
        {
            fprintf(file," %i", f[j]+1);
        }

        fprintf(file,"\n");
    }
}



unsigned int edges_common_face(unsigned int i_e0, unsigned int i_e1, vector< vector<unsigned int> > edge_faces)
{
    // for all faces adjacent to edge i_e0
    for(unsigned int i = 0; i < edge_faces[i_e0].size(); i++)
    {
        // and for all faces adjacent to edge i_e1
        for(unsigned int j = 0; j < edge_faces[i_e1].size(); j++)
        {
            // return identical face index
            if(edge_faces[i_e0][i] == edge_faces[i_e1][j])
                return edge_faces[i_e0][i];
        }
    }
    
    // no common face has been found -> error display
    cerr << "common face not found between edges " << i_e0 << " and " << i_e1 << "." << endl;
    
    return 0;
}

unsigned int vertices_common_face(unsigned int i0, unsigned int i1, unsigned int i2, vector< vector<unsigned int> > vert_faces)
{
    // for all faces adjacent to vertex i0
    for(unsigned int i = 0; i < vert_faces[i0].size(); i++)
    {
        // and for all faces adjacent to vertex i1
        for(unsigned int j = 0; j < vert_faces[i1].size(); j++)
        {
            // if thay have a common face
            if(vert_faces[i0][i] == vert_faces[i1][j])
            {
                // for all faces adjacent to vertex i2
                for(unsigned int k = 0; k < vert_faces[i2].size(); k++)
                {
                    // returns the common face
                    if(vert_faces[i0][i] == vert_faces[i2][k])
                        return vert_faces[i0][i];
                }
            }
        }
    }
    
    // no common face has been found -> error display
    cerr << "common face not found between vertices " << i0 << ", " << i1 << ", and " << i2 << "." << endl;
    
    return 0;
}

unsigned int trouverFaceCommune(vector< vector< unsigned int > > &facesOfEdge, unsigned int e1, unsigned int e2) {

    if (facesOfEdge[e1][0] == facesOfEdge[e2][0]) {
        return facesOfEdge[e1][0];

    } else if (facesOfEdge[e1][0] == facesOfEdge[e2][1]) {
        return facesOfEdge[e1][0];

    } else if (facesOfEdge[e1][1] == facesOfEdge[e2][0]) {
        return facesOfEdge[e1][1];

    } else if (facesOfEdge[e1][1] == facesOfEdge[e2][1]) {
        return facesOfEdge[e1][1];
    }

    return -1;
}

int trouverDoublon(vector< vec3 > &vertices, vec3 vertex) {
    for (int i = 0; i < vertices.size(); ++i) {
        vec3 p = vertices[i];
        if (p.x == vertex.x && p.y == vertex.y && p.z == vertex.z) {
            return i;
        }
    }
    return -1;
}

Mesh Mesh::subdivide() const
{
    Mesh output;
    

    
    //=======================================================
    //
    // TODO : implémenter le schema de subdivision de Catmull-Clark
    //
    //=======================================================


    //1. Ajouter le barycentre de chaque face du maillage

    vector<vec3> listBFaces; // liste des barycentre de chaque face. listBFaces[i] = indince barycentre de faces[i] dans la liste des vertices
    vector<vec3> listBEdges; // liste des barycentre de chaque face. listBEdges[i] = indince barycentre de get_edges()[i] dans la liste des vertices
    for (int i = 0; i < faces.size(); ++i) {

        vec3 somme(0, 0, 0);

        //Recherche du barycentre b de la face f
        for (int j = 0; j < get_face(i).size(); ++j) {
            somme += get_vertex(get_face(i)[j]);
        }

        somme /= (float) get_face(i).size();
        listBFaces.push_back(somme);
    }

    //1-2. Ajouter le barycentre de chaque arête

//    for (int i = 0; i < get_edges().size(); ++i) {
//        Edge e = get_edges()[i];

//        vec3 somme;
//        somme += get_vertex(e.m_i0);
//        somme += get_vertex(e.m_i1);
//        somme /= 2.0;
//        listBEdges.push_back(somme);
//    }


    //2. Ajouter le barycentre du tetraedre

    cout <<"1" << endl;
    vector< Edge > edges = get_edges();
    cout <<"2" << endl;
    vector< vector< unsigned int > > facesOfEdge = get_edge_faces(edges);

    vector<vec3> listBTetraedre;

    for (int i = 0; i < edges.size(); ++i) {
        Edge e = edges[i];
        unsigned int s0 = e.m_i0;
        unsigned int s1 = e.m_i1;
        unsigned int f0 = facesOfEdge[i][0];
        unsigned int f1 = facesOfEdge[i][1];

        // Barycentre situe au meme indice que l'arete correspondante
        vec3 barycentre;
        barycentre += get_vertex(s0);
        barycentre += get_vertex(s1);
        barycentre += listBFaces[f0];
        barycentre += listBFaces[f1];
        barycentre /= 4.0;

        listBTetraedre.push_back(barycentre);
    }


    //3. Déplacer les points du maillage

    vector< vector< unsigned int > > listVertexFaces = get_vertex_faces();
    vector< vector< unsigned int > > listVertexEdges = get_vertex_edges(get_edges());

    vector< vec3 > listVertexTmp;

    for (int i = 0; i < vertices.size(); ++i) {
        vector< unsigned int > listFaces = listVertexFaces[i];
        vector< unsigned int > listEdges = listVertexEdges[i];

        // Barycentre des sfi
        vec3 F;
        for (int k = 0; k < listFaces.size(); ++k) {
            F += listBFaces[listFaces[k]];
        }
        F /= (float) listFaces.size();

        // Barycentre des sai
        vec3 A;
        for (int k = 0; k < listEdges.size(); ++k) {
            //A += listBEdges[k];
            A += listBTetraedre[listEdges[k]];
        }
        A /= (float) listEdges.size();

        vec3 s = get_vertex(i);
        int n = vertices.size();
        float x = (F[0] + 2*A[0] + (n - 3) * s[0]) / n;
        float y = (F[1] + 2*A[1] + (n - 3) * s[1]) / n;
        float z = (F[2] + 2*A[2] + (n - 3) * s[2]) / n;

        listVertexTmp.push_back(vec3(x, y, z));
    }


    //4. Former les faces avec les quadruplets

    vector< unsigned int > faces;
    vector< vec3 > vertices_mesh;
    vector< vector< unsigned int > > faces_mesh;

    for (int i = 0; i < listVertexTmp.size(); ++i) {
        vec3 s = listVertexTmp[i];
        //cout << "x " << s.x << ", y " << s.y << ", z " << s.z << endl;

        int nbArete = listVertexEdges[i].size();
//        int w = 0;
//        int e1;
//        int e2 = listVertexEdges[i][0];
//        int k = 0;
//        do {
//            e1 = e2;
//            int f = facesOfEdge[e1][k];
//            e2 = -1;
//            for (int j = 0; j < listVertexEdges[i].size() && e2 == -1; ++j) {
//                if (j != e1) {
//                    if (facesOfEdge[j][0] == f) {
//                        e2 = j;
//                        k = 1;
//                    } else if (facesOfEdge[j][1] == f) {
//                        e2 = j;
//                        k = 0;
//                    }

//                }
//            }

//            vec3 s1 = listBTetraedre[e1];
//            vec3 s2 = listBFaces[f];
//            vec3 s3 = listBTetraedre[e2];

//            int ts = trouverDoublon(vertices_mesh, s);
//            if (ts == -1) {
//                vertices_mesh.push_back(s);
//                ts = vertices_mesh.size()-1;
//            }
//            int ts1 = trouverDoublon(vertices_mesh, s1);
//            if (ts1 == -1) {
//                vertices_mesh.push_back(s1);
//                ts1 = vertices_mesh.size()-1;
//            }
//            int ts2 = trouverDoublon(vertices_mesh, s2);
//            if (ts2 == -1) {
//                vertices_mesh.push_back(s2);
//                ts2 = vertices_mesh.size()-1;
//            }
//            int ts3 = trouverDoublon(vertices_mesh, s3);
//            if (ts3 == -1) {
//                vertices_mesh.push_back(s3);
//                ts3 = vertices_mesh.size()-1;
//            }

//            faces.clear();
//            faces.push_back(ts);
//            faces.push_back(ts1);
//            faces.push_back(ts2);
//            faces.push_back(ts3);

//            faces_mesh.push_back(faces);

//            w++;
//        } while (w < nbArete);


        for (int j = 0; j < nbArete; ++j) {
            int e1 = listVertexEdges[i][j];
            int e2 = listVertexEdges[i][(j+1)%nbArete];
            int f = trouverFaceCommune(facesOfEdge, e1, e2);

            vec3 s1 = listBTetraedre[e1];
            vec3 s2 = listBFaces[f];
            vec3 s3 = listBTetraedre[e2];

            int ts = trouverDoublon(vertices_mesh, s);
            if (ts == -1) {
                vertices_mesh.push_back(s);
                ts = vertices_mesh.size()-1;
            }
            int ts1 = trouverDoublon(vertices_mesh, s1);
            if (ts1 == -1) {
                vertices_mesh.push_back(s1);
                ts1 = vertices_mesh.size()-1;
            }
            int ts2 = trouverDoublon(vertices_mesh, s2);
            if (ts2 == -1) {
                vertices_mesh.push_back(s2);
                ts2 = vertices_mesh.size()-1;
            }
            int ts3 = trouverDoublon(vertices_mesh, s3);
            if (ts3 == -1) {
                vertices_mesh.push_back(s3);
                ts3 = vertices_mesh.size()-1;
            }

            faces.clear();
            faces.push_back(ts);
            faces.push_back(ts3);
            faces.push_back(ts2);
            faces.push_back(ts1);

            faces_mesh.push_back(faces);
        }
    }


    output.vertices = vertices_mesh;
    output.faces = faces_mesh;


    //output = *this;     // place holder : current mesh copy


    return output;
}


vector< vector< unsigned int > > Mesh::get_neighborhoods() const
{
    vector< vector< unsigned int > > output;
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        output.push_back(vector< unsigned int >());
    }
    
    // unordered neighborhood computation
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);
        for(unsigned int j=0; j < f.size(); j++)
        {
            unsigned int next = f[(j+1)%f.size()];
            if(std::find(output[f[j]].begin(), output[f[j]].end(), next) == output[f[j]].end())
                output[f[j]].push_back(next);

            next = f[(j+f.size()-1)%f.size()];
            if(std::find(output[f[j]].begin(), output[f[j]].end(), next) == output[f[j]].end())
                output[f[j]].push_back(next);
        }
    }

    // neighborhood ordering
    vector< vector< unsigned int > > output2;
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        output2.push_back(vector< unsigned int >());
        
        output2[i].push_back(output[i][0]);
        for(unsigned int j = 1; j < output[i].size(); j++)
        {
            unsigned int i_prev = i;
            unsigned int i_next = output2[i][j-1];

            while(i_next != i)
            {
                unsigned int temp = i_next;
                i_next = get_next_vert(i_prev, i_next);
                i_prev = temp;
            }


            output2[i].push_back(i_prev);
        }
    }

    return output2;
}


std::vector< Edge > Mesh::get_edges() const
{
    vector< vector< unsigned int > > neib = get_neighborhoods();
    vector< Edge> output;
    
    // for all vertices
    for(unsigned int i = 0; i < neib.size(); i++)
    {
        // for all neighbors
        for(unsigned int j=0; j < neib[i].size(); j++)
        {
            Edge e(i, neib[i][j]);
            
            // add the edge if it is not already created
            if(std::find(output.begin(), output.end(), e) == output.end())
                output.push_back(e);
        }
    }
    
    return output;
}


vector< vector< unsigned int > > Mesh::get_vertex_faces() const
{
    vector< vector< unsigned int > > output;
    
    // for all vertices
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        // create an empty incident face set
        output.push_back(vector< unsigned int >());
    }
    
    // for all faces
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);
        
        // for all composing vertices
        for(unsigned int j=0; j < f.size(); j++)
        {
            // add self as incident
            output[f[j]].push_back(i);
        }
    }
    
    return output;
}


vector< vector< unsigned int > > Mesh::get_vertex_edges(vector< Edge > edges) const
{
    vector< vector< unsigned int > > output;
    
    // for all vertices
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        // create an empty incident edge set
        output.push_back(vector< unsigned int >());
    }
    
    // for all edges
    for(unsigned int i = 0; i < edges.size(); i++)
    {
        // add self as incident to both extremities
        output[edges[i].m_i0].push_back(i);
        output[edges[i].m_i1].push_back(i);
    }
    
    return output;
}


vector< vector< unsigned int > > Mesh::get_vertex_edges(vector< vector< unsigned int > > vert_neib, vector< Edge > edges) const
{
    vector< vector< unsigned int > > output;
    
    // for all vertices
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        // add an empty incident edges set
        output.push_back(vector< unsigned int >());
        
        // for all incident edges
        for(unsigned int j = 0; j < vert_neib[i].size(); j++)
        {
            // add corresponding edge index
            unsigned int k = std::find(edges.begin(), edges.end(), Edge(i, vert_neib[i][j])) - edges.begin();
            output[i].push_back(k);
        }
    }
    
    return output;
}


vector< vector< unsigned int > > Mesh::get_edge_faces(vector< Edge > edges) const
{
    vector< vector< unsigned int > > vf = get_vertex_faces();
    
    vector< vector< unsigned int > > output;
    
    // for all edges
    for(unsigned int i = 0; i < edges.size(); i++)
    {
        // add an empty incident face set
        output.push_back(vector< unsigned int >());
        
        // retrieve extremities indices
        unsigned int p0 = edges[i].m_i0;
        unsigned int p1 = edges[i].m_i1;
        
        // for all faces adjacent to one extremity
        for(unsigned int j0 = 0; j0 < vf[p0].size(); j0++)
        {
            // and for all faces adjacent to the other extremity
            for(unsigned int j1 = 0; j1 < vf[p1].size(); j1++)
            {
                // add identical faces
                if(vf[p0][j0] == vf[p1][j1])
                    output[i].push_back(vf[p0][j0]);
            }
        }
    }
    
    return output;
}


unsigned int Mesh::get_oriented_face(unsigned int i_e0, unsigned int i_e1) const
{
    // for all faces
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);
        
        // for all composing vertices
        for(unsigned int j=0; j < f.size(); j++)
        {
            // return the right face
            if(f[j] == i_e0 && f[(j+1)%f.size()] == i_e1)
                return i;
        }
    }
    
    // no common face has been found -> error display
    cerr << "oriented face not found for edge " << i_e0 << " -> " << i_e1 << "." << endl;
    
    return 0;
}

unsigned int Mesh::get_next_vert(unsigned int i_e0, unsigned int i_e1) const
{
    // for all faces
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        vector< unsigned int > f = get_face(i);

        // for all composing vertices
        for(unsigned int j = 0; j < f.size(); j++)
        {
            // return the following vertex index in the right face
            if(f[j] == i_e0 && f[(j+1)%f.size()] == i_e1)
                return f[(j+2)%f.size()];
        }
    }
    
    // no common face has been found -> error display
    cerr << "Next vertex not found for edge " << i_e0 << " -> " << i_e1 << "." << endl;
    
    return 0;
}


