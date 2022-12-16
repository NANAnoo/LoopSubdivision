///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity({0.0,0.0,0.0})
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	first_directed_edge_of_vertex.resize(0);
	faces.resize(0);
	other_half_of_edge.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];

    // reading lines & adding them in appropriate places
    while (true) { // not eof
        // token for identifying meaning of line
        std::string token;

        // character to read
        geometryStream >> token;

        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
        if (token == "#") { // comment
            // read and discard the line
            geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
        } // comment
        else if (token == "Vertex") { // vertex
            // variables for the read
            unsigned int vertexID;
            geometryStream >> vertexID;
            // it has to be next valid 0-based ID, so
            // reject line if it isn't
            if (vertexID != vertices.size()) { // bad vertex ID
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad vertex ID

            // read in the new vertex position
            float x, y, z;
            geometryStream >> x >> y >> z;
            vertices.push_back({x, y, z});
        } // vertex
        else if (token == "Normal") { // normal
            // variables for the read
            unsigned int normalID;
            geometryStream >> normalID;
            // it has to be next valid 0-based ID, so
            // reject line if it isn't
            if (normalID != normals.size()) { // bad ID
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad ID

            // read in the new normal
            float x, y, z;
            geometryStream >> x >> y >> z;
            normals.push_back({x, y, z});
        } // normal
        else if (token == "FirstDirectedEdge") { // first directed edge
            // variables for the read
            unsigned int FDE_ID;
            geometryStream >> FDE_ID;
            // it has to be next valid 0-based ID, so
            // reject line if it isn't
            if (FDE_ID != first_directed_edge_of_vertex.size()) { // bad ID
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad ID

            // read in the new FDE
            HalfEdgeRef value;
            geometryStream >> value;
            first_directed_edge_of_vertex.push_back(value);
        } // first directed edge
        else if (token == "Face") { // face
            // variables for the read
            unsigned int faceID;
            geometryStream >> faceID;
            // it has to be next valid 0-based ID, so
            // reject line if it isn't
            if (faceID != faces.size()) { // bad face ID
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad face ID

            // read in the new face vertex
            // read in the new normal
            unsigned int v0, v1, v2;
            geometryStream >> v0 >> v1 >> v2;
            faces.push_back({v0, v1, v2});
        } // face
        else if (token == "OtherHalf") { // other half
            // variables for the read
            unsigned int otherHalfID;
            geometryStream >> otherHalfID;
            // it has to be next valid 0-based ID, so
            // reject line if it isn't
            if (otherHalfID != other_half_of_edge.size()) { // bad ID
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad ID

            // read in the new other half edge index
            HalfEdgeRef value;
            geometryStream >> value;
            other_half_of_edge.push_back(value);
        } // other half
    } // not eof
	// compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = {0.0, 0.0, 0.0};

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set
	return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    // First Block, information
    geometryStream << "# University of Leeds 2022-2023\n";
    geometryStream << "# Vertices=" << getVertexCount() << " Faces=" << getFaceCount() << std::endl;
    geometryStream << "# \n";

    // Second Block, Vertices
    geometryStream.precision(6);
    for (int i = 0; i < getVertexCount(); i++) {
        // Vertex $(index) $(x) $(y) $(z)
        geometryStream.setf(std::ios::showbase);
        geometryStream << "Vertex " << i;
        geometryStream.setf(std::ios::showpoint);
        geometryStream << " " << vertices[i].x
               << " " << vertices[i].y
               << " " << vertices[i].z
               << std::endl;
    }

    // Third Block, normals
    for (int i = 0; i < normals.size(); i++) {
        // Normal $(index) $(x) $(y) $(z)
        geometryStream.setf(std::ios::showbase);
        geometryStream << "Normal " << i;
        geometryStream.setf(std::ios::showpoint);
        geometryStream << " " << normals[i].x
               << " " << normals[i].y
               << " " << normals[i].z
               << std::endl;
    }

    // Fourth Block, first directed edge
    for (int i = 0; i < first_directed_edge_of_vertex.size(); i++) {
        // FirstDirectedEdge $(index) $(FirstDirectedEdge)
        geometryStream.setf(std::ios::showbase);
        geometryStream << "FirstDirectedEdge " << i;
        geometryStream << " " << first_directed_edge_of_vertex[i]
               << std::endl;
    }

    // Fifth Block, Faces
    geometryStream.setf(std::ios::showbase);
    for (int i = 0; i < getFaceCount(); i++) {
        // Face $(index) $(v0_index) $(v1_index) $(v2_index)
        geometryStream << "Face " << i
               << " " << faces[i].vertex_index[0]
               << " " << faces[i].vertex_index[1]
               << " " << faces[i].vertex_index[2]
               << std::endl;
    }

    // Last block, other half edge
    for (int i = 0; i < other_half_of_edge.size(); i++) {
        // OtherHalf $(index) $(OtherHalf)
        geometryStream.setf(std::ios::showbase);
        geometryStream << "OtherHalf " << i;
        geometryStream << " " << (int) ((other_half_of_edge[i] == UNKNOWN_HALF_EDGE) ? -1 : other_half_of_edge[i])
               << std::endl;
    }
    } // WriteObjectStream()

void DirectedEdgeSurface::findNeighbours(unsigned int vertex, std::vector<unsigned int> &ring) {
    HalfEdgeRef first_edge = firstDirectedHalfEdgeOnVertex(vertex);
    HalfEdgeRef next_edge = first_edge;
    do {
        next_edge = nextHalfEdge(otherHalfEdge(next_edge));
        ring.push_back(toVertexIndexOfHalfEdge(next_edge));
    } while (next_edge != first_edge);
}

void DirectedEdgeSurface::loopSubdivisionLocally() {
    // first_half_directed_edge and second_half_directed_edge of each old half edge
    // init as UNKNOWN_HALF_EDGE to indicate "a vertex that hasn't been accessed"
    std::vector<HalfEdgeRef> first_half_directed_edge(other_half_of_edge.size(), UNKNOWN_HALF_EDGE);
    std::vector<HalfEdgeRef> second_half_directed_edge(other_half_of_edge.size(), UNKNOWN_HALF_EDGE);

    std::vector<Vec3> new_pos_of_old_vertices(vertices.size());
    // update position of old vertices in place
    std::vector<unsigned int> one_ring;
    for (unsigned int v = 0; v < vertices.size(); v++) {
        one_ring.clear();
        findNeighbours(v, one_ring);
        unsigned int ring_size = one_ring.size();
        float miu = ring_size == 3 ? 3.f / 16 : 3.f / (8 * float(ring_size));
        float center_w = 1.f - float(ring_size) * miu;
        new_pos_of_old_vertices[v].x = vertices[v].x * center_w;
        new_pos_of_old_vertices[v].y = vertices[v].y * center_w;
        new_pos_of_old_vertices[v].z = vertices[v].z * center_w;
        for (auto n: one_ring) {
            new_pos_of_old_vertices[v].x += vertices[n].x * miu;
            new_pos_of_old_vertices[v].y += vertices[n].y * miu;
            new_pos_of_old_vertices[v].z += vertices[n].z * miu;
        }
    }

    // a map to store where there is a middle point on an edge
    struct edge {
        unsigned int v1;
        unsigned int v2;
    };
    struct edge_hash {
        unsigned int operator()(const edge &e) const {
            return (e.v1 + 1) * (e.v2 + 1) * 157;
        }
    };
    struct edge_eq {
        bool operator()(const edge &e1, const edge &e2) const {
            return (e1.v1 == e2.v1 && e1.v2 == e2.v2)
                   || (e1.v1 == e2.v2 && e1.v2 == e2.v1);
        }
    };
    std::unordered_map<edge, unsigned int, edge_hash, edge_eq> mid_vid_on_edge;

    // resize first directed edge and set all value as UNKNOWN_HALF_EDGE
    unsigned int result_vertices_count = getVertexCount() + other_half_of_edge.size() / 2;
    unsigned int result_edge_count = 4 * other_half_of_edge.size();
    first_directed_edge_of_vertex = std::vector<HalfEdgeRef>(result_vertices_count, UNKNOWN_HALF_EDGE);
    other_half_of_edge.resize(result_edge_count, UNKNOWN_HALF_EDGE);

    // main loop is based the count of old face number
    unsigned int old_face_cnt = getFaceCount();
    unsigned int old_vertex_cnt = getVertexCount();
    // Start loop subdivision, iterate all faces
    // insert vertices and faces
    for (unsigned int f_index = 0; f_index < old_face_cnt; f_index++) {
        // every old face are divided into four
        // keep the old face index and set it as the center triangle
        Face old_f = faces[f_index];
        Face newFaces[3];
        unsigned int next_face_id = faces.size();
        unsigned int new_face_ids[3] = {next_face_id, next_face_id + 1, next_face_id + 2};
        // First, get midpoint ids on three edges
        unsigned int mid_vertex_ids[3];
        // iterate three old edges in CCW
        for (unsigned int i = 0; i < 3; i++) {
            unsigned int old_vertex_from = old_f.vertex_index[i];
            unsigned int old_vertex_to = old_f.vertex_index[(i + 1) % 3];
            // check if the midpoint of current edge has an id
            edge e = {old_vertex_from, old_vertex_to};
            auto ptr = mid_vid_on_edge.find(e);
            if (ptr != mid_vid_on_edge.end()) {
                // use the existed id
                mid_vertex_ids[i] = ptr->second;
            } else {
                // alloc a new id
                mid_vertex_ids[i] = vertices.size();
                mid_vid_on_edge[e] = mid_vertex_ids[i];
                // generate a new vertex
                HalfEdgeRef from_to = nextHalfEdge(3 * f_index + i);
                unsigned int vertex_right = old_f.vertex_index[(i + 2) % 3];
                unsigned int vertex_left =
                        toVertexIndexOfHalfEdge(nextHalfEdge(otherHalfEdge(from_to)));
                // add a new vertex
                vertices.push_back({3.f / 8 * (vertices[old_vertex_from].x + vertices[old_vertex_to].x) +
                                 1.f / 8 * (vertices[vertex_left].x + vertices[vertex_right].x),
                                 3.f / 8 * (vertices[old_vertex_from].y + vertices[old_vertex_to].y) +
                                 1.f / 8 * (vertices[vertex_left].y + vertices[vertex_right].y),
                                 3.f / 8 * (vertices[old_vertex_from].z + vertices[old_vertex_to].z) +
                                 1.f / 8 * (vertices[vertex_left].z + vertices[vertex_right].z)
                                });
                // add a new normal
                normals.push_back({
                    0.5f * (normals[old_vertex_from].x + normals[old_vertex_to].x),
                    0.5f * (normals[old_vertex_from].y + normals[old_vertex_to].y),
                    0.5f * (normals[old_vertex_from].z + normals[old_vertex_to].z),
                });
            }
        }
        // update first directed edge of the old and new vertex
        for (int i = 0; i < 3; i++) {
            // old vertex
            unsigned int v = old_f.vertex_index[i];
            if (first_directed_edge_of_vertex[v] == UNKNOWN_HALF_EDGE) {
                first_directed_edge_of_vertex[v] = 3 * new_face_ids[i] + 1;
            }
            // new vertex
            v = mid_vertex_ids[i];
            HalfEdgeRef e = 3 * f_index + (i + 1) % 3;
            if (first_directed_edge_of_vertex[v] == UNKNOWN_HALF_EDGE) {
                first_directed_edge_of_vertex[v] = e;
            }
        }
        // update first & second directed edges
        for (int i = 0; i < 3; i++) {
            // old edge from the vertex
            HalfEdgeRef old_edge = nextHalfEdge(3 * f_index + i);
            // update first & second half edge
            first_half_directed_edge[old_edge] = 3 * new_face_ids[i] + 1;
            second_half_directed_edge[old_edge] = 3 * new_face_ids[(i + 1) % 3];
        }
        // update vertices id of new faces
        // then push new faces into mesh
        for (unsigned int i = 0; i < 3; i++) {
            // use the old vertex as the first vertex in the new face
            newFaces[i].vertex_index[0] = old_f.vertex_index[i];
            // use the middle vertex as the new id in CCW
            newFaces[i].vertex_index[1] = mid_vertex_ids[i];
            newFaces[i].vertex_index[2] = mid_vertex_ids[(i + 2) % 3];
            // push it into mesh
            faces.push_back(newFaces[i]);
        }
        // update the vertex id of old face
        for (unsigned int i = 0; i < 3; i++) {
            faces[f_index].vertex_index[i] = mid_vertex_ids[i];
        }
    }
    // update other half
    for (unsigned int old_edge = 0; old_edge < first_half_directed_edge.size(); old_edge++) {
        // update half edge on new faces
        HalfEdgeRef first_half = first_half_directed_edge[old_edge];
        HalfEdgeRef other_of_first_half = second_half_directed_edge[other_half_of_edge[old_edge]];
        other_half_of_edge[first_half] = other_of_first_half;
        other_half_of_edge[other_of_first_half] = first_half;
        // update half edge on old faces
        HalfEdgeRef second_half = second_half_directed_edge[old_edge];
        other_half_of_edge[old_edge] = second_half + 2;
        other_half_of_edge[second_half + 2] = old_edge;
    }

    // update old vertex positions
    for (unsigned int v = 0; v < old_vertex_cnt; v++) {
        vertices[v].x = new_pos_of_old_vertices[v].x;
        vertices[v].y = new_pos_of_old_vertices[v].y;
        vertices[v].z = new_pos_of_old_vertices[v].z;
    }
}

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faces.size(); face ++)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Vec3 pq = vertices[faces[face].vertex_index[1]] - vertices[faces[face].vertex_index[0]];
			Vec3 pr = vertices[faces[face].vertex_index[2]] - vertices[faces[face].vertex_index[0]];

			// take their cross product and normalise
			Vec3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex : faces[face].vertex_index)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[vertex].x * scale,
					normals[vertex].y * scale,
					normals[vertex].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[vertex].x,
				vertices[vertex].y,
				vertices[vertex].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

