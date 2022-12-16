Assignment 2
    - Loop Subdivision

Name : Hao Zhang
ID   : 201670059

Description:

    - [done] Implemented with the standard 3
        - Edit the triangles in place, there is only one DirectedEdgeSurface object during the whole time.
          ~ NO 2nd or more objects
          ~ NOT using the simple rebuild half-edge data
          ~ Generate new vertex, new face one by one with repairing the half-edge

        - Add a button in GUI to subdivide surface into a bigger level (0->1->2->...), 
          and two labels that show current level and vertices count in current surface

    - [done] Each vertex remains vertex index after subdivision 

    - [done] Each edge have a new vertex, and the new vertices are added at the end of the array sequentially.

    - [done] Each triangle is devided into four, and the centre one remains the same index, 
             while the other three are added at the end of the array sequentially.

Build:

- feng-linux

    module load legacy-eng
    module add gcc
    module add qt/5.13.0

    mkdir build && cd build
    qmake ../
    make -j20

Run:

./LoopSubdivision ${path_to_models}/model_name.diredgenormal
