Some General conventions used throughout the code:
1. point - is a tuple of (x,y,z)
2. Points - is a list of points in the form [(x1,y1,z1), (x2,y2,z2)]
3. triangle - a tuple of (p1,p2,p3) which defines the three points that 
   form the triangle through their index in the Points list
4. triangles - a list of triangles in the form [(p1,p2,p3), (p4,p5,p6)]

There may be confusion as point may refer either to (x,y,z) or its index in Points
Same thing for triangles: triangle may refer either to (p1,p2,p3) or its index in Triangles
Sorry for any confusion

There are two important files: findtri.py and fixtri.py

findtri.py
    To use: run init(fileName)
        - fileName must have either .stl or .vtk extension
        - fileName should refer to a file in the same level as the program
        - IMPORTANT: first set runlength to a reasonable amount as to not overload your computer;
          to find all intersecting triangles, reset runlength to len(tris)
        - will return a vtk of intersecting triangles
        - Saves three lists:
            - intersectingTriangles.npy
            - triangles.npy
            - points.npy
            - these are important for fixtri.py
    You can make changes in init(fileName) depending on what you want
        - intersectingTriangles is the array containing the intersecting triangles in the form
          [(t1,t2), (t3,t4)] where t1 and t2 intersect and t3 and t4 intersect
          t1 is the index of a triangle in the list "triangles"
        - if you don't want to search over the entire set of trianlges, set tris equal to the list
          of triangles you want to search over. Set tris = triangles to search over entire file
          - ex: find intersections in a certain bounding box using the list returned from findTrisInBox()
        - writeFileWithTriangles will write a vtk file from a list of triangles
        - writeFileWithTwoTris will write a vtk file for 2 triangles (meant to be called iteratively)
            ex:
            for pair in intersectingTriangles:
                writeFileWithTwoTris(pair[0], pair[1], fileName)

fixtri.py
This is a multistep process
1. Make sure intersectingTriangles.npy, triangles.npy, and points.npy have been created by findtri.py.
   Load them in.
2. Run:
     tri_groups = getTriangleGroups(triint, intersecting)
     np.save("triangle_groups", tri_groups)
   This will group intersecting triangles into "chains" or "groups" that are connected by intersections. 
   This is necessary because different size chains may need to be dealt with differently.
3. Load in tri_groups if necesasry and run sortGroups(tri_groups)
   Sort groups into different sizes, .npy files will be saved into new directory fixtriNumpyArrays
   - error_groups - groups with 2 triangles; called error groups because this should very likely be empty.
     However, this may differ with your project
   - four_groups - groups with 4 triangles (special case)
   - multi_groups - groups with 3 or 5 or more triangles (everything else)
   - Theoretically there should ONLY be intersections in four_groups. If there are any in error_groups,
     your triangles in your mesh likely aren't connected with each other. If there are any in multi_groups,
     two entire sheets of triangles are likely overlapping.
4. Solve four_groups
    4a) Load in four_groups if necessary and run
          getFourList(four_groups)
          fourList = np.load("4list.npy")
        fourList is a list of all triangles in four_groups with no duplicates
    4b) [OPTIONAL] Run writeFileWithTriangles(fourList, "intersectingQuads.vtk")
        Writes a vtk file with the intersecting triangles to visually see where they are.
    4c) Run fixFourGroups(four_groups)
        Name self-explanatory; fixes fourGroups and returns a list of new triangles.
        Note: fixFourGroups only fixes quads; four groups that don't intersect like quads will be ignored
        and not added to list of new triangles
    4d) Run finalWriter() to create vtk with the new fixed four_groups



   
