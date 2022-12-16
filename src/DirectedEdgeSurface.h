///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.h
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

// include guard for DirectedEdgeSurface
#ifndef _DIRECTED_EDGE_SURFACE_H
#define _DIRECTED_EDGE_SURFACE_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#include <limits>
#include "math.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// the render parameters
#include "RenderParameters.h"
// the image class for a texture
#include "RGBAImage.h" 

#define UNKNOWN_HALF_EDGE std::numeric_limits<unsigned int>::max()

class DirectedEdgeSurface
    { // class DirectedEdgeSurface
    public:
    // structures in directed edge
    // record xyz value

    // I want to use the implementation from Assignment 1 
    // so I move some codes from Assignment 1, to make things easier
    typedef unsigned int HalfEdgeRef;
    struct Vec3
    {
        float x;
        float y;
        float z;

        // only impl operators that used in render,
        Vec3 operator+(const Vec3 &other) {
            return {x + other.x, y + other.y, z + other.z};
        }
        Vec3 operator-(const Vec3 &other) {
            return {x - other.x, y - other.y, z - other.z};
        }
        Vec3 operator/(float value) {
            return {x / value, y / value, z / value};
        }

        Vec3 cross(const Vec3 &other) {
            return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
        }

        Vec3 unit() {
            float l = length();
            return {x / l, y / l, z / l };
        }

        float length() {
            return sqrt(x * x + y * y + z * z);
        }
    };
    // record vertices index in CCW
    struct Face
    {
        unsigned int vertex_index[3];
    };
    
    // getters
    unsigned int getFaceCount() const { return faces.size(); }
    unsigned int getVertexCount() const { return vertices.size(); }
    const DirectedEdgeSurface::Face& getFaceAt(unsigned int index) { return faces[index]; }
    const DirectedEdgeSurface::Vec3& getVertexAt(unsigned int index) { return vertices[index]; }

    // ------------------------------all methods below this line are O(1)------------------------------------------------------
    // get index of the face that the edge is on
    static inline unsigned int faceIndexOfHalfEdge(HalfEdgeRef edge) { return edge / 3;}

    // get the from_vertex index of the edge
    // edge : [from vertex] -> to vertex
    inline unsigned int fromVertexIndexOfHalfEdge(HalfEdgeRef edge) {
        return getFaceAt(faceIndexOfHalfEdge(edge)).vertex_index[(edge + 2) % 3];
    }

    // get the to_vertex index of an edge
    // edge : from vertex -> [to vertex]
    inline unsigned int toVertexIndexOfHalfEdge(HalfEdgeRef edge) {
        return this->faces[faceIndexOfHalfEdge(edge)].vertex_index[edge % 3];
    }

    // first half edge of the vertex
    inline HalfEdgeRef firstDirectedHalfEdgeOnVertex(unsigned int vertex_index) {
        return first_directed_edge_of_vertex[vertex_index];
    }
    // next half edge of the edge on the same face
    static inline HalfEdgeRef nextHalfEdge(HalfEdgeRef edge) {
        return (edge + 1) % 3 + 3 * faceIndexOfHalfEdge(edge);
    }
    // previous half edge of the edge on the same face
    static inline HalfEdgeRef prevHalfEdge(HalfEdgeRef edge) {
        return (edge + 2) % 3 + 3 * faceIndexOfHalfEdge(edge);
    }
    // the opposite edge of the edge
    inline HalfEdgeRef otherHalfEdge(HalfEdgeRef edge) {
        return other_half_of_edge[edge];
    }
    // ------------------------------all methods until this line are O(1)------------------------------------------------------
    
    /* store the vertex coords */
    std::vector<Vec3> vertices;

    /* store the vertex normals */
    std::vector<Vec3> normals;

    /* store the face vertices index in each level*/
    std::vector<Face> faces;

    /* store the first directed edge*/
    std::vector<HalfEdgeRef> first_directed_edge_of_vertex;

    /* store the opposite edge of each edge*/
    std::vector<HalfEdgeRef> other_half_of_edge;

    /* get one-ring of a vertex */
    void findNeighbours(unsigned int vertex, std::vector<unsigned int> &ring);

    // centre of gravity - computed after reading
    Vec3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    // constructor will initialise to safe values
    DirectedEdgeSurface();

    // Loop Subdivision, core part for assignment 2
    void loopSubdivisionLocally();
    
    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    // routine to render
    void Render(RenderParameters *renderParameters);

    }; // class DirectedEdgeSurface

// end of include guard for DirectedEdgeSurface
#endif
