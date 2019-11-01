import numpy as np
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter


# ******************** FUNCTIONS *********************

#parses file for points
def getPoints(data):
    arr = []
    #print(data.GetNumberOfPoints())
    for i in range(data.GetNumberOfPoints()):
        #print(i)
        tupledata = data.GetPoint(i)
        #print(tupledata)
        arr.append(tupledata)

    np_arr = np.array(arr)
    return np_arr

#parses file for triangles
def getTriangles(data):
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

    #for i in range(0, numCells):
    #    numIds = 3
    #    print(tupledata.GetCell(cellLocation,numIds))
    #    cellLocation += 1 + numIds
    #tupledata.InitTraversal()

# checks if triangles has the same point
# return false if two triangles have 2 same points
def pointSameCheck(triangle1, triangle2):
    
    if len(list(set(triangle1).intersection(set(triangle2)))) >= 1:
        return False
    return True
    
def pointNearCheck(triangle1, triangle2, rf):

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

# returns bounding box of triangle in [(x1, x2), (y1, y2), (z1, z2)]
# triangle should be tuple of point values eg (1,3,4)
def boundingBox (triangle):
    
    # triangle = np.array([0,2,3])
    # points = np.array([(1.1,2.3,3.2), (2.5,3.1,4.5), (5.6, 4.6, 5.8), (7.8,6.8,9.7)])
    box=[]
    box.append(getHighestandLowest(points[triangle[0]][0], points[triangle[1]][0], points[triangle[2]][0]))
    box.append(getHighestandLowest(points[triangle[0]][1], points[triangle[1]][1], points[triangle[2]][1]))
    box.append(getHighestandLowest(points[triangle[0]][2], points[triangle[1]][2], points[triangle[2]][2]))
    return box

# boxes should be in format [(x1, x2), (y1, y2), (z1, z2)]
# returns true or false
def boundingBoxCheck(box1, box2):
    if box1[0][0] > box2[0][1] or box1[0][1] < box2[0][0] or \
        box1[1][0] > box2[1][1] or box1[1][1] < box2[1][0] or \
        box1[2][0] > box2[2][1] or box1[2][1] < box2[2][0]:
        
        return False
    else:
        return True

# planeFormula of form (a,b,c,d)
# return type of discriminant
def planeCheck(triangle1, triangle2):
    
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

    
    # return 0 for false
    # return 1 for 2 on one side one on other
    # return 2 for 1 point on plane with two points on either side
    
    
# will return true or false depending on whether the triangles intersect
# triangle should be tuple of point indeces referring to points array eg (1,3,4)
def triangleIntersect(triangle1, triangle2):

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
    

    
# ************************ HELPER FUNCTIONS ***********************

# points is of form [(0.0, 0.0, 1.0), (3, 0.0, 1.0), (4, 0.0, 1.0), (1.0, 0.0, 1.0)]
#points1 and points2 are int points
def pointComparer(points1, points2):
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


# intersection function
def isect_line_plane_v3(p0, p1, p_co, p_no, epsilon=1e-6):
    """ 
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
    return (
        v0[0] + v1[0],
        v0[1] + v1[1],
        v0[2] + v1[2],
        )


def sub_v3v3(v0, v1):
    return (
        v0[0] - v1[0],
        v0[1] - v1[1],
        v0[2] - v1[2],
        )


def dot_v3v3(v0, v1):
    return (
        (v0[0] * v1[0]) +
        (v0[1] * v1[1]) +
        (v0[2] * v1[2])
        )


def len_squared_v3(v0):
    return dot_v3v3(v0, v0)


def mul_v3_fl(v0, f):
    return (
        v0[0] * f,
        v0[1] * f,
        v0[2] * f,
        )

def getHighestandLowest(x, y, z):
    # helper function for boundingBox
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

#returns tuple corresponding to coeff of plane formula
def findPlane(triangle):
    point1 = points[triangle[0]]
    point2 = points[triangle[1]]
    point3 = points[triangle[2]]
    v1 = np.subtract(point3,point2)
    v2 = np.subtract(point2,point1)
    cp = np.cross(v1, v2)
    a, b, c = cp
    # d = -1*(np.dot(cp, point3))
    d = -1*(a*point1[0] + b*point1[1] + c*point1[2])
    return (a,b,c,d)

# tells which side of the plane point is on
# 1 for >, 0 for =, -1 for <
def pointPlane(planeFormula, point):
    discriminant = planeFormula[0]*point[0] + planeFormula[1]*point[1] + planeFormula[2]*point[2] + planeFormula[3]
    if round(discriminant, 6) == 0:
        return 0
    elif round(discriminant, 6) < 0:
        return -1
    else:
        return 5

#returns the length and fac multiplied by weighted factors summed
def getScaleFactor(length, fac):
    return 1.0 * length / 2.75 / 2.843 + 0.5 - fac
    
#returns intersectingTriangles (list of tuples of triangle pairs)
def runner(runSize):
    trueCount = 0
    intersectingTriangles = []

    intersectingWithFactor = []
    boxes = []
    for y in range(runSize):
        boxes.append(boundingBox(triangles[y]))
    for y in range(runSize):
        for z in range(y+1,runSize):
            if pointSameCheck(triangles[y], triangles[z]) and \
            boundingBoxCheck(boxes[y], boxes[z]) and \
            pointNearCheck(triangles[y], triangles[z], 2) and \
            planeCheck(triangles[y], triangles[z]):

                intersectlen, fac = triangleIntersect(triangles[y], triangles[z])

                if intersectlen > 0 and fac>0:
                    
                    intersectingWithFactor.append( [(y,z) , getScaleFactor(intersectlen, fac)] )
                    trueCount += 1
                
    print ("trueCount",trueCount)

    #sort intersectingWithFactor by scale factor
    intersectingWithFactor = sorted(intersectingWithFactor, key=itemgetter(1), reverse = True)
    

    for thing in intersectingWithFactor:
        intersectingTriangles.append(thing[0])

    return (intersectingTriangles, trueCount)

# writes out to intersect.vtk - puts all pairs in a single file
def writeOneFile(intersecting, fileName):
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

# list of intersecting triangles with an index 
def writeMultipleFiles(triangle1, triangle2, fileName):
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




# ************************ RUNNER STUFF ***********************
reader = vtk.vtkPolyDataReader()
reader.SetFileName('vessel_triangles.vtk')
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

data = reader.GetOutput()

points = getPoints(data)

triangles = getTriangles(data)

#total triangles: 15576
runLength = 15576
print("runlength: " + str(runLength))
intersecting, trueCount = runner(runLength)


fileDirectory = "./garbage"
if not os.path.exists(fileDirectory):
    os.makedirs(fileDirectory)

# counter = 1
# for pair in intersecting:
    
#     writeMultipleFiles(pair[0], pair[1], fileDirectory + "/int{:04}".format(counter) + ".vtk")
#     counter += 1

print ("trueCount",trueCount)


# np.save("intersectingTriangles", intersecting)

# writeOneFile(intersecting)
    



