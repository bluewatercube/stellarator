import numpy as np
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter

'''
***********************************************************************
                                FUNCTIONS
***********************************************************************
'''

def getPoints(data):
    """
    Parses STL or VTK data for points
    Parameters:
    data: STL or VTK data from getOutput()
    Returns:
    ndarr: Points as a numpy array of tuples
    """
    arr = []
    #print(data.GetNumberOfPoints())
    for i in range(data.GetNumberOfPoints()):
        #print(i)
        tupledata = data.GetPoint(i)
        #print(tupledata)
        arr.append(tupledata)

    np_arr = np.array(arr)
    return np_arr

def getTriangles(data):
    """
    Parses STL or VTK data for triangles
    Parameters:
    data: STL or VTK data from getOutput()
    Returns:
    list: Triangles as a list of point indeces
    """
    arr = []
    cellArray = data.GetPolys() #this is vtkCellArray
    cellArray.InitTraversal()
    idList = vtk.vtkIdList()
    
    #print(cellArray.GetNumberOfCells())

    for i in range(0,cellArray.GetNumberOfCells()):
        cellArray.GetNextCell(idList)
        #print(type(idList))
        if idList.GetNumberOfIds()==3:
            point1 = idList.GetId(0)
            point2 = idList.GetId(1)
            point3 = idList.GetId(2)
            # print("point 1,2,3")
            # print(point1, point2, point3)
            arr.append((int(point1), int(point2), int(point3)))
    
    return arr

def pointSameCheck(triangle1, triangle2):
    """
    Checks if two triangles share the same point
    Returns false if two triangles have 2 same points
    """
    if len(list(set(triangle1).intersection(set(triangle2)))) >= 1:
        return False
    return True

def pointNearCheck(triangle1, triangle2, rf):
    """
    Checks if points of 2 triangles are within 0.000001 of each other
    Parameters:
    triangle1: tuple of triangle 1's point indeces
    triangle2: tuple of triangle 2's point indeces
    rf: round factor
    """
    #append 6 3-number tuples representing triangle points
    trianglepoints = []
    for ptIndex in triangle1:
        #trianglepoints.append(points[ptIndex])
        trianglepoints.append((
            round(points[ptIndex][0], rf),
            round(points[ptIndex][1], rf),
            round(points[ptIndex][2], rf),
        ))
        
    for ptIndex in triangle2:
        trianglepoints.append((
            round(points[ptIndex][0], rf),
            round(points[ptIndex][1], rf),
            round(points[ptIndex][2], rf),
        ))
    
    #print(set(trianglepoints))

    if len(set(trianglepoints)) != 6:
        return False
    return True

def boundingBox (triangle):
    """
    Returns bounding box of triangle in [(x1, x2), (y1, y2), (z1, z2)]
    Triangle should be tuple of point values eg (1,3,4)
    Parameters:
    triangle: tuple of triangle's point indeces
    """
    # triangle = np.array([0,2,3])
    # points = np.array([(1.1,2.3,3.2), (2.5,3.1,4.5), (5.6, 4.6, 5.8), (7.8,6.8,9.7)])
    box=[]
    box.append(getHighestandLowest(points[triangle[0]][0], points[triangle[1]][0], points[triangle[2]][0]))
    box.append(getHighestandLowest(points[triangle[0]][1], points[triangle[1]][1], points[triangle[2]][1]))
    box.append(getHighestandLowest(points[triangle[0]][2], points[triangle[1]][2], points[triangle[2]][2]))
    return box

def boundingBoxCheck(box1, box2):
    """
    Determines if two triangle's bounding boxes overlap
    Boxes should be in format [(x1, x2), (y1, y2), (z1, z2)]
    """
    if box1[0][0] > box2[0][1] or box1[0][1] < box2[0][0] or \
        box1[1][0] > box2[1][1] or box1[1][1] < box2[1][0] or \
        box1[2][0] > box2[2][1] or box1[2][1] < box2[2][0]:
        
        return False
    else:
        return True

def planeCheck(triangle1, triangle2):
    """
    Returns True if 2 points are on one side of the plane and 1 is on the other side
    Returns False if all 3 points are on one side of the plane or all three points are on the plane
    """
    planeForm1 = findPlane(triangle1)
    planeForm2 = findPlane(triangle2)
    #check triangle1 against planeForm2
    discriminant1 = 0
    for point in triangle1:
        discriminant1 += pointPlane(planeForm2, points[point])
    #cases: all 3 on one side: 15, -3
    #       2 on one side, 1 on plane: -2, 10
    #       1 on one side, 2 on plane: -1, 5
    #       all on plane: 0
    #       2 on one side one on other: 9, 3
    if discriminant1 in [15, 0, -3]:
        return False
   
    #check triangle2 against planeForm1
    discriminant2 = 0
    for point in triangle2:
        discriminant2 += pointPlane(planeForm1, points[point])
    if discriminant2 in [15, 0, -3]:
        return False

    return True    
    
def triangleIntersect(triangle1, triangle2):
    """
    Applies triangle intersection check algorithm to determine if intersecting
    Parameters:
    triangle1: tuple of triangle 1's point indeces
    triangle2: tuple of triangle 2's point indeces
    Returns:
    bool: if triangles intersect
    """

    a1,b1,c1,d1 = findPlane(triangle1)
    a2,b2,c2,d2 = findPlane(triangle2)
    
    vector1 = (a1,b1,c1)
    vector2 = (a2,b2,c2)

    pointPairs1 = [(points[triangle1[0]], points[triangle1[1]]), 
                    (points[triangle1[1]], points[triangle1[2]]),
                    (points[triangle1[2]], points[triangle1[0]])]
    
    pointPairs2 = [(points[triangle2[0]], points[triangle2[1]]), 
                    (points[triangle2[1]], points[triangle2[2]]),
                    (points[triangle2[2]], points[triangle2[0]])]
    
    intersectionPoints1 = []
    intersectionPoints2 = []

    fac_list = []
    
    for pair in pointPairs1:
        intersect, fac = isect_line_plane_v3(pair[0], pair[1], points[triangle2[0]],vector2)
        if fac >= 0 and fac <= 1:
            fac_list.append(abs(0.5 - fac))
            intersectionPoints1.append(intersect)

    for pair in pointPairs2:
        intersect, fac = isect_line_plane_v3(pair[0], pair[1], points[triangle1[0]],vector1)
        if fac >= 0 and fac <= 1:
            fac_list.append(abs(0.5 - fac))
            intersectionPoints2.append(intersect)

    fac = 0

    if len(fac_list) > 0:
        fac = 1.0 * sum(fac_list)/len(fac_list)
    

    intersectionPoints1 = list(set(intersectionPoints1))
    intersectionPoints2 = list(set(intersectionPoints2))

    if len(intersectionPoints1) == 1:
        intersectionPoints1.append(intersectionPoints1[0])
    elif len(intersectionPoints2) == 1:
        intersectionPoints2.append(intersectionPoints2[0])

    if len(intersectionPoints1) != 2 or len(intersectionPoints2) != 2:

        planeForm1 = findPlane(triangle1)
        planeForm2 = findPlane(triangle2)
        
        totDis1,totDis2 = 0,0
        
        #check triangle1 against planeForm2
        for point in triangle1:
            # print(pointPlane(planeForm1, points[point]))
            discriminant1 = planeForm2[0]*points[point][0] + planeForm2[1]*points[point][1] + planeForm2[2]*points[point][2] + planeForm2[3]
            if round(discriminant1, 6) == 0:
                pass
            elif round(discriminant1, 6) < 0:
                totDis1 += -1
            else:
                totDis1 += 5
        
        #check triangle2 against planeForm1
        for point in triangle2:
            # print(pointPlane(planeForm2, points[point]))
            discriminant2 = planeForm1[0]*points[point][0] + planeForm1[1]*points[point][1] + planeForm1[2]*points[point][2] + planeForm1[3]
            if round(discriminant2, 6) == 0:
                pass
            elif round(discriminant2, 6) < 0:
                totDis2 += -1
            else:
                totDis2 += 5


        return (0, 0)

        
    return ( pointComparer(intersectionPoints1, intersectionPoints2), fac )
    

'''
***********************************************************************
                            HELPER FUNCTIONS
***********************************************************************
'''

def pointComparer(points1, points2):
    """
    Based upon ordering of the points, returns whether the triangles intersect or not
    points is of form [(0.0, 0.0, 1.0), (3, 0.0, 1.0), (4, 0.0, 1.0), (1.0, 0.0, 1.0)]
    points1 and points2 are int points
    """
    index = -1
    if points1[0][0] != points2[0][0] and \
        points1[1][0] != points2[0][0] and \
        points1[0][0] != points2[1][0] and \
        points2[1][0] != points2[1][0]:
        index = 0
    elif points1[0][1] != points2[0][1] and \
        points1[1][1] != points2[0][1] and \
        points1[0][1] != points2[1][1] and \
        points2[1][1] != points2[1][1]:
        index = 1
    else:
        index = 2
    
    comparater1 = [points1[0][index], points1[1][index]]
    comparater2 = [points2[0][index], points2[1][index]]

    low1, high1, low2, high2 = 0,0,0,0

    if comparater1[0] < comparater1[1]:
        low1 = comparater1[0]
        high1 = comparater1[1]
    else:
        low1 = comparater1[1]
        high1 = comparater1[0]

    if comparater2[0] < comparater2[1]:
        low2 = comparater2[0]
        high2 = comparater2[1]
    else:
        low2 = comparater2[1]
        high2 = comparater2[0]


    if low1-high2 >= 0 or low2-high1 >= 0:
        return 0
    
    #case 1: tri1 is in tri2
    if low1 >= low2 and high1 <= high2:
        return high1 - low1
    #case 2: tri2 is in tri1
    if low1 <= low2 and high1 >= high2:
        return high2 - low2
    #case 3: tri1 is linked with tri2 and tri 1 is higher
    if high1 > high2 and low1 > low2:
        return high2 - low1
    # case 4: tri2 is linked with tri 1 and tri2 is higher
    if high2 > high1 and low2 > low1:
        return high1 - low2 

def isect_line_plane_v3(p0, p1, p_co, p_no, epsilon=1e-6):
    """ 
    Helps determine if triangles intersect
    p0, p1: define the line
    p_co, p_no: define the plane:
        p_co is a point on the plane (plane coordinate).
        p_no is a normal vector defining the plane direction;
             (does not need to be normalized).

    return a Vector or None (when the intersection can't be found).
    """

    u = sub_v3v3(p1, p0)
    dot = dot_v3v3(p_no, u)

    if abs(dot) > epsilon:
        # the factor of the point between p0 -> p1 (0 - 1)
        # if 'fac' is between (0 - 1) the point intersects with the segment.
        # otherwise:
        #  < 0.0: behind p0.
        #  > 1.0: infront of p1.
        w = sub_v3v3(p0, p_co)
        fac = 1.0* -dot_v3v3(p_no, w) / dot
        u = mul_v3_fl(u, fac)
        return (add_v3v3(p0, u), round(fac,6))
    else:
        # The segment is parallel to plane
        return (None, None)

def add_v3v3(v0, v1):
    """ Performs vector addition """
    return ( v0[0] + v1[0],v0[1] + v1[1], v0[2] + v1[2])

def sub_v3v3(v0, v1):
    """ Performs vector subtraction """
    return (v0[0] - v1[0],v0[1] - v1[1],v0[2] - v1[2])

def dot_v3v3(v0, v1):
    """ Calculates dot product of vectors """
    return ( (v0[0] * v1[0]) + (v0[1] * v1[1]) + (v0[2] * v1[2]) )

def len_squared_v3(v0):
    """ Returns magnitude of vector """
    return dot_v3v3(v0, v0)

def mul_v3_fl(v0, f):
    """ Performs scalar multiplication """
    return ( v0[0] * f, v0[1] * f, v0[2] * f )

def getHighestandLowest(x, y, z):
    """ Helper function for boundingBox """
    highest = x
    lowest = x
    if x < y:
        highest = y
    else: 
        lowest = y
    
    if z > highest:
        highest = z

    if z < lowest:
        lowest = z

    return (lowest, highest)

def findPlane(triangle):
    """ Returns tuple corresponding to normal vector of plane """
    point1 = points[triangle[0]]
    point2 = points[triangle[1]]
    point3 = points[triangle[2]]
    v1 = np.subtract(point3,point2)
    v2 = np.subtract(point2,point1)
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = -1*(a*point1[0] + b*point1[1] + c*point1[2])
    return (a,b,c,d)

def pointPlane(planeFormula, point):
    """
    Tells which side of the plane the point is on
    1 for >, 0 for =, -1 for <
    """
    discriminant = planeFormula[0]*point[0] + planeFormula[1]*point[1] + planeFormula[2]*point[2] + planeFormula[3]
    if round(discriminant, 6) == 0:
        return 0
    elif round(discriminant, 6) < 0:
        return -1
    else:
        return 5

def getScaleFactor(length, fac):
    """ Returns the length and fac multiplied by weighted factors summed """
    return 1.0 * length / 2.75 / 2.843 + 0.5 - fac
    

'''
***********************************************************************
                                RUNTIME
***********************************************************************
'''

def getIntersectingTris(triangles,runSize):
    """
    Runner function that gets list of intersecting triangles
    Parameters:
    tris: list of triangles
    runSize: number of triangles to check against each other
    Returns:
    list of tuples of triangle pairs
    """
    trueCount = 0
    intersectingTriangles = []

    intersectingWithFactor = []
    boxes = []  #array of triangle bounding boxes
    for i in range(0, 100):
        print(triangles[i])
    for y in range(runSize):
        boxes.append(boundingBox(triangles[y]))
    for y in range(runSize):
        if y % 100 == 0:
            print str(y) + " / " + str(runSize)
        for z in range(y+1,runSize):
            tri1 = triangles[y]
            tri2 = triangles[z]
            if pointSameCheck(tri1, tri2) and \
            boundingBoxCheck(boxes[y], boxes[z]) and \
            pointNearCheck(tri1, tri2, 2) and \
            planeCheck(tri1, tri2):

                intersectlen, fac = triangleIntersect(tri1, tri2)

                if intersectlen > 0 and fac>0:
                    
                    intersectingWithFactor.append( [(tri1, tri2) , getScaleFactor(intersectlen, fac)] )
                    trueCount += 1
                
    print str(runSize) + " / " + str(runSize)
    print "Finished detecting intersections."
    #print ("trueCount",trueCount)

    #sort intersectingWithFactor by scale factor
    intersectingWithFactor = sorted(intersectingWithFactor, key=itemgetter(1), reverse = True)
    

    for thing in intersectingWithFactor:
        intersectingTriangles.append(thing[0])

    return (intersectingTriangles, trueCount)

def writeFileWithTwoTris(triangle1, triangle2, fileName):
    """
    Writes file containing two traingles
    filename of the form "file.vtk"
    """
    fileOut = open(fileName, "w+")
    fileOut.write("# vtk DataFile Version 3.0\n")
    fileOut.write("Intersecting Triangles\n")
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")
    fileOut.write("POINTS 6 float\n")

    #pt is int representing index in points array
    #point is a tuple with XYZ coord
    for pt in triangles[triangle1]:
        point = points[pt]
        fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")

    for pt in triangles[triangle2]:
        point = points[pt]
        fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")

    fileOut.write("POLYGONS 2 8\n")
    fileOut.write("3 0 1 2\n")
    fileOut.write("3 3 4 5\n")

def writeFileWithTriangleTuples(intersecting, fileName):
    """
    same as writeFileWithTriangle but with the triangles being tuples instead of indices
    """
    fileOut = open(fileName, "w")
    fileOut.write("# vtk DataFile Version 3.0\n")
    fileOut.write("Intersecting Triangles\n")
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")

    numPoints = len(points)
    numCells = 2*len(intersecting) #num triangles
    cellListSize = numCells * 4

    fileOut.write("POINTS " + str(numPoints) + " float\n")
    #print points
    for point in points:
        fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")

    fileOut.write("POLYGONS " + str(numCells) + " " + str(cellListSize)+ "\n")
    for trianglePair in intersecting:
        for triangle in trianglePair:
            fileOut.write("3 " + " ".join(map(str,triangle))+"\n")

def writeFileWithTriangles(triangleList, fileName):
    """
    fileName should be in the form "file.vtk"
    triangleList is the list of intersecting triangles
    writes out to intersect.vtk - puts all pairs in a single file
    """
    fileOut = open(fileName, "w+")
    fileOut.write("# vtk DataFile Version 3.0\n")
    fileOut.write("Intersecting Triangles\n")
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")
    fileOut.write("POINTS %d float\n" % (3*len(triangleList)))

    #pt is int representing index in points array
    #point is a tuple with XYZ coord
    for triangle in triangleList:
        for pt in triangles[triangle]:
            point = points[pt]
            fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")

    fileOut.write("POLYGONS %d %d\n" % (len(triangleList), 4*len(triangleList)) )
    for x in range(len(triangleList)):
        fileOut.write("3 %d %d %d\n" % (3*x, 3*x+1, 3*x+2))
    fileOut.close()

def findTrisInBox(x1,x2,y1,y2,z1,z2):
    """ Returns indexes of triangles in bounding box """
    tris_in_box = []
    for i in range(len(triangles)):
        tri = triangles[i]
        in_box = True
        for point in tri:
            pt = points[point] #access tuple from points arr
            x = pt[0]
            y = pt[1]
            z = pt[2]
            if not( x > x1 and x < x2 and \
                y > y1 and y < y2 and \
                z > z1 and z < z2 ):
                in_box = False
                break
        if in_box == True:
            tris_in_box.append(i)
    return tris_in_box

def createFixedVTK(new_triangles, filename):
    """ new_triangles is a list of tuples of point indices """
    fileOut = open(filename, "w+")
    fileOut.write("# vtk DataFile Version 3.0\n")
    group_string = ""
    fileOut.write("Fixed Triangles\n")
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")
    
    num_points = len(new_triangles) * 3

    fileOut.write("POINTS " + str(num_points) + " float\n")
    
    for tri in new_triangles:
        #pt is index of pt tuple in larger points list
        for pt in tri:
            point = points[pt]
            fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")
    
    #get num triangles
    num_tri = len(new_triangles)
    num_to_rep = num_tri * 4

    fileOut.write("POLYGONS %s %s\n" % (str(num_tri), str(num_to_rep)))
    
    count = 0
    for i in range(0,num_tri):
        fileOut.write("3 %s %s %s\n" % ( str(count), str(count+1), str(count+2) ))
        count += 3

def getVtkData(filename):
    """ Gets data from VTK """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    return data

def getStlData(filename):
    """ Gets data from STL """
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    return data

def init(filename):
    """
    Extracts points and triangles from file and assigns to global variables
    filename example  = "file.stl"
    file extensions supported: stl and vtk
    """
    global points
    global triangles

    #get file type
    fileType = filename[-3:]
    if fileType == "stl":
        data = getStlData(filename)
    elif fileType == "vtk":
        data = getVtkData(filename)
    else:
        raise Exception('File type not supported')
    
    #parse data for points and triangles
    points = getPoints(data)
    triangles = getTriangles(data)
    
    # #save to numpy arrays for use in fixtri.py
    # np.save("points", points)
    # np.save("triangles", triangles)
    
    # runLength = len(triangles)

    # #trueCount is number of intersections
    # #intersecting is list of intersecting triangle pairs
    # print("runlength: " + str(runLength))
    # intersecting, trueCount = getIntersectingTris(triangles, runLength)

    # print("Number of intersections: " + (str(trueCount) ))

    # # creates intersectingTriangles.npy; needed for fixtri.py
    # np.save("intersectingTriangles", intersecting)
    intersecting= np.load("intersectingTriangles.npy")
    # write file to intersecting.vtk
    writeFileWithTriangleTuples(intersecting, "intersecting.vtk")
    


''' *********** GLOBAL VARIABLES *********** '''
points = []
triangles = []

createFixedVTK(triangles, "fixedStellarator.vtk")

# mins = [-25.9, -19.7, 19.8]
# maxes = [24.2, 17.8, 56.5]
# box_tris = findTrisInBox(mins[0],maxes[0],mins[1],maxes[1],mins[2],maxes[2])
# print(box_tris)
# print(len(box_tris))
# print(len(triangles))

#total triangles: 15576 (for vessel_triangles.vtk)
#total for inout.vtk: 8947
# precision: 1390
#ves_lr1: 234160
# cut down: 865
# repaired = 1924
# runLength = 865
# print("runlength: " + str(runLength))
# 

# fileDirectory = "./resolved"
# if not os.path.exists(fileDirectory):
#     os.makedirs(fileDirectory)

# counter = 1
# for pair in intersecting:
    
#     writeMultipleFiles(pair[0], pair[1], fileDirectory + "/int{:04}".format(counter) + ".vtk")
#     counter += 1
# print(intersecting)
# print ("trueCount",trueCount)



    



