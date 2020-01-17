import numpy as np
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter
import scipy
from scipy.spatial import Delaunay

"""
*********************************************************************
                            PARSE DATA FUNCTIONS
*********************************************************************
"""

#don't really need these anymore?
# def getPoints(data):
#     """
#     Parses STL or VTK data for points
#     Parameters:
#     data: STL or VTK data from getOutput()
#     Returns:
#     ndarr: Points as a numpy array of tuples
#     """
#     arr = []
#     #print(data.GetNumberOfPoints())
#     for i in range(data.GetNumberOfPoints()):
#         tupledata = data.GetPoint(i)
#         arr.append(tupledata)

#     np_arr = np.array(arr)
#     return np_arr

# def getTriangles(data):
#     """
#     Parses STL or VTK data for triangles
#     Parameters:
#     data: STL or VTK data from getOutput()
#     Returns:
#     list: Triangles as a list of point indeces
#     """
#     arr = []
    
#     cellArray = data.GetPolys() #this is vtkCellArray
#     cellArray.InitTraversal()
#     idList = vtk.vtkIdList()

#     for i in range(0,cellArray.GetNumberOfCells()):
#         cellArray.GetNextCell(idList)
#         if idList.GetNumberOfIds()==3:
#             point1 = idList.GetId(0)
#             point2 = idList.GetId(1)
#             point3 = idList.GetId(2)
#             arr.append((int(point1), int(point2), int(point3)))
    
#     return arr

def getTriangleGroups(data, intersecting):
    """
    Returns a list of groups of triangles that are intersecting together
    Format: [(group), (group), (group)] ex: [(1, 4, 5), (2, 3), (6, 7, 8, 9)]
    Each group is a tuple of indeces
    """
    triangles = list(data[:])
    groups = []
    visited = []
    while len(triangles):
        groups.append([triangles[0]])
        searchOrder = [triangles[0]]
        visited.append(triangles[0])
        while len(searchOrder) > 0:
            #print("sersearchOrder", searchOrder)
            searchObject = searchOrder[0]
            for pair in intersecting:
                if pair[0] == searchObject and pair[1] not in visited:
                    groups[len(groups)-1].append(pair[1]) #add tri2 to SO
                    searchOrder.insert(0, pair[1])
                    visited.append(pair[1])
                elif pair[1] == searchObject and pair[0] not in visited:
                    groups[len(groups)-1].append(pair[0])
                    searchOrder.insert(0, pair[0])
                    visited.append(pair[0])
                    
            searchOrder.remove(searchObject)

        # print("triangles", triangles) 
        # print("groups", groups)
        # print("triangles", triangles)
        for item in groups[len(groups)-1]:
            triangles.remove(item)
    
    return groups
    
def dfs(visited,triangles,node,search):
    if node not in visited:
        visited.append(node)
        for neighbor in search[node]:
            dfs(visited,triangles,neighbor,search)

#gets list of inter triangles
def getTriList(triangles):
    triabc = set()
    for pair in triangles:
        triabc.add(pair[0])
        triabc.add(pair[1])
    
    np.save("triint", list(triabc))
    
def getFourList(triangles):
    triabc = set()
    for pair in triangles:
        triabc.add(pair[0])
        triabc.add(pair[1])
        triabc.add(pair[2])
        triabc.add(pair[3])
    
    np.save("4list", list(triabc))

""" *********************************************************************
                            SEARCHING FUNCTIONS
    ********************************************************************* """

def identifyIntersectType(triarr):
    """
    Returns two lists as a tuple: (line over triangles, node over triangles)
    it's type "line over" if a triangle intersects with 2 triangles
    and an adjacent triangle intersects with the same 2 triangles
    """
          
    intersectionsList = triarr[:]
    
    #check line over
    for pair in triarr:
        inter1 = getInter(pair[0], triarr) #--> tri2 (other tri in pair) and tri3
        inter2 = getInter(pair[1], triarr) #--> tri1 (other tri in pair) and tri4
        
        tripair = TriPair(pair[0], pair[1], inter1, inter2)
        intersectionsList.append(tripair)
        
def vtk_tri_groups(tri_groups):  
    """
    Parameter: groups from getTriangleGroups (list of lists)
    sorts into different types of triangle groups:
        error groups (1 or 2 triangles)
        quads (4 triangles)
        multigroups (>4 triangles)
    """ 
    fileDirectory = "./vtktrigroups"
    direrror = "./vtktrigroups/groupserror"
    dir4 = "./vtktrigroups/groups4"
    dirmulti = "./vtktrigroups/groupsmulti"
    if not os.path.exists(fileDirectory):
        os.makedirs(fileDirectory)   
        os.makedirs(direrror)
        os.makedirs(dir4)
        os.makedirs(dirmulti)  
    
    counterror, count4, countmulti = 0, 0, 0
    
    for group in tri_groups:
        if len(group) in [1,2]:
            filename = direrror + "/int{:04}".format(counterror) + ".vtk"
            counterror += 1
        elif len(group) == 4:
            filename = dir4 + "/int{:04}".format(count4) + ".vtk"
            count4 += 1
        else:
            filename = dirmulti + "/int{:04}".format(countmulti) + ".vtk"
            countmulti += 1
        writeGroupFile(group, filename)
                
def writeGroupFile(group, filename):
    """
    Parameters: group is list of triangle indeces
    writes the file for a single group of triangles
    """
    fileOut = open(filename, "w+")
    fileOut.write("# vtk DataFile Version 3.0\n")
    group_string = ""
    for tri in group:
        group_string += str(tri) + " "
    fileOut.write("Intersecting Triangles %s\n" % (group_string))
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")

    #get number of points
    num_points = len(group) * 3
    fileOut.write("POINTS " + str(num_points) + " float\n")
    
    #add point data to file
    for tri in group:
        #pt is index of pt tuple in larger points list
        for pt in triangles[tri]:
            point = points[pt]
            fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")
    
    #get num triangles
    num_tri = len(group)
    num_to_rep = num_tri * 4

    fileOut.write("POLYGONS %s %s\n" % (str(num_tri), str(num_to_rep)))
    
    count = 0
    for i in range(0,num_tri):
        fileOut.write("3 %s %s %s\n" % ( str(count), str(count+1), str(count+2) ))
        count += 3

def findDoubles(group):
    """
    returns a list of triangle intersecting pairs
    """
    global triangles
    # print("points")
    # for tri in group:
    #     print(triangles[tri])
        
    doubles = []
    for x in range(len(group)):
        for y in range(x+1, len(group)):
            if len(list(set(triangles[group[x]]).intersection(set(triangles[group[y]])))) == 2:
                doubles.append([group[x], group[y]])
    return doubles

def lineErrorSearch(group):
    """
    returns groups of quads 
    """
    doubles = findDoubles(group)
    # print("doubles")
    # print(doubles)
    quads = []
    for x in range(len(doubles)):
        for y in range(x+1, len(doubles)):
            if doubles[x][0] != doubles[y][0] and \
                doubles[x][1] != doubles[y][0] and \
                doubles[x][0] != doubles[y][1] and \
                doubles[x][1] != doubles[y][1]:
                # [doubles[x][0], doubles[y][0]] in intersectinglist and \
                # [doubles[x][1], doubles[y][0]] in intersectinglist and \
                # [doubles[x][0], doubles[y][1]] in intersectinglist and \
                # [doubles[x][1], doubles[y][1]] in intersectinglist:
                quads.append([(doubles[x][0], doubles[x][1]), (doubles[y][0], doubles[y][1])])
    return quads
    
def getInter(tri, triarr):
    """
    returns list of triangle tuples that given tri intersects with
    """
    interarr = []
    for pair in triarr:
        if tri == pair[0]:
            interarr.append(pair[1])
        elif tri == pair[1]:
            interarr.append(pair[0])
    
    return interarr

def sortGroups(tri_groups):
    """
    saves different group types into np file
    used mainly for identifying four triangle groups 
    """
    error_groups = []
    four_groups = []
    multi_groups = []

    fileDirectory = "./fixtriNumpyArrays"
    if not os.path.exists(fileDirectory):
        os.makedirs(fileDirectory)

    for group in tri_groups:
        if len(group) == 1 or len(group) == 2:
            error_groups.append(group)
        elif len(group) == 4:
            four_groups.append(group)
        else:
            multi_groups.append(group)
    np.save("./fixtriNumpyArrays/error_groups", error_groups)
    np.save("./fixtriNumpyArrays/four_groups", four_groups)
    np.save("./fixtriNumpyArrays/multi_groups", multi_groups)

""" *********************************************************************
                            FIXING FUNCTIONS
    ********************************************************************* """
#delete
# def error_triangle_groups(tri_groups):
#     # trigroup2 = []
#     # for group in tri_groups:
#     #     if len(group)==2 or len(group)==1:
#     #         trigroup2.append(group)
#     # vtk_tri_groups(trigroup2)
#     pass

def four_triangle_groups(four_groups):
    """
    print which four-groups are quads
    """
    for group in four_groups:
        print(lineErrorSearch(group))

def multi_triangle_groups(multi_groups):    
    for group in multi_groups:
        print(lineErrorSearch(group))

def fixQuad(quad):
    """
    Parameter: quad is list of 2 tuples
    fixes quads by redrawing the triangles from the same set of 4 points
    ex: A[BC] and [BC]D become [AD]B and [AD]C
    determine which double is on inside surface --> indoub
    """
    index_avg_1 = sum(quad[0])/2.0
    index_avg_2 = sum(quad[1])/2.0

    if index_avg_1 > index_avg_2:
        indoub = quad[0]
        outdoub = quad[1]
    else:
        indoub = quad[1]
        outdoub = quad[0]
    
    #determine shared pts and indiv pts
    tri1 = triangles[indoub[0]]
    tri2 = triangles[indoub[1]]

    shared_pts = list(set(tri1).intersection(tri2))
    indiv_pts = list((set(tri1) - set(tri2) | set(tri2) - set(tri1)))
    #shared becomes indiv, indiv becomes shared
    new_tri1 = indiv_pts + [shared_pts[0]]
    new_tri2 = indiv_pts + [shared_pts[1]]
    
    return [tuple(new_tri1), tuple(new_tri2), triangles[outdoub[0]], triangles[outdoub[1]]]

def fixFourGroups(four_groups):
    """
    Goes through all four groups and fixes them with fixQuad()
    Returns all_new_tris: list of tuples of new triangle indeces
    indeces refer to points array
    [(p1, p2, p3),(p4, p5, p6), (p7, p8, p9)]
    """
    #print four_groups
    all_new_tris = []
    quads = []

    # get all quads from the four groups (because some fourgroups aren't quads)
    for group in four_groups:
        quad = lineErrorSearch(group)
        if quad != []:
            quads += quad
    
    # iterate over quads, combine each list of fixed triangles
    # with list of all fixed triangles
    for quad in quads:
        new_tris = fixQuad(quad)
        all_new_tris += new_tris
    
    return all_new_tris

def createFixedVTK(new_triangles, filename):
    """
    Creates a VTK with new fixed triangles
    new_triangles is a list of tuples of point indeces
    """
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

def writeFileWithTriangles(triangleList, fileName):
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
    

def finalwriter(fileName):
    global triangles
    global new_tris
    #edit triangles to remove those in four groups
    #fourList is list of indexes
    # for i in range(len(triangles)):
    #     if i in fourList:
    #         print(triangles[i])
    #         del triangles[i]
    
    value_list = []

    for index in fourList:
        value_list.append(triangles[index])
    for val in value_list:
        triangles.remove(val)
    
    triangles += new_tris
    print("length")
    print(len(triangles))

    createFixedVTK(triangles,fileName)
    
""" ************************ GET DATA *********************** """

# reader = vtk.vtkPolyDataReader()
# reader.SetFileName('ves_lr1.vtk')
# reader.ReadAllScalarsOn()
# reader.ReadAllVectorsOn()
# reader.Update()

# data = reader.GetOutput()


# points = getPoints(data)
# triangles = getTriangles(data) #list of tuples of point indexes
points = np.load("points.npy")
triangles = np.load("triangles.npy")
    print(len(triangles))

    createFixedVTK(triangles,fileName)
    
""" ************************ GET DATA *********************** """

# reader = vtk.vtkPolyDataReader()
# reader.SetFileName('ves_lr1.vtk')
# reader.ReadAllScalarsOn()
# reader.ReadAllVectorsOn()
# reader.Update()

# data = reader.GetOutput()


# points = getPoints(data)
# triangles = getTriangles(data) #list of tuples of point indexes
points = np.load("points.npy")
triangles = np.load("triangles.npy")
    print(len(triangles))

    createFixedVTK(triangles,fileName)
    
""" ************************ GET DATA *********************** """

# reader = vtk.vtkPolyDataReader()
# reader.SetFileName('ves_lr1.vtk')
# reader.ReadAllScalarsOn()
# reader.ReadAllVectorsOn()
# reader.Update()

# data = reader.GetOutput()


# points = getPoints(data)
# triangles = getTriangles(data) #list of tuples of point indexes
points = np.load("points.npy")
triangles = np.load("triangles.npy")
    print(len(triangles))

    createFixedVTK(triangles,fileName)
    
""" ************************ GET DATA *********************** """

# reader = vtk.vtkPolyDataReader()
# reader.SetFileName('ves_lr1.vtk')
# reader.ReadAllScalarsOn()
# reader.ReadAllVectorsOn()
# reader.Update()

# data = reader.GetOutput()


# points = getPoints(data)
# triangles = getTriangles(data) #list of tuples of point indexes
points = np.load("points.npy")
triangles = np.load("triangles.npy")
    print(len(triangles))

    createFixedVTK(triangles,fileName)
    
""" ************************ GET DATA *********************** """

# reader = vtk.vtkPolyDataReader()
# reader.SetFileName('ves_lr1.vtk')
# reader.ReadAllScalarsOn()
# reader.ReadAllVectorsOn()
# reader.Update()

# data = reader.GetOutput()


# points = getPoints(data)
# triangles = getTriangles(data) #list of tuples of point indexes
points = np.load("points.npy")
triangles = np.load("triangles.npy")

# get intersecting triangles
intersecting = np.load("intersectingTriangles.npy")
intersectinglist = intersecting.tolist() #list of 2-int lists of tri indexes
# print(intersectinglist)
triint = np.load("triint.npy")
# writeFileWithTriangles(triint, "trinint.vtk")
# print(intersecting)
""" ************************ RUNNER ************************** """
# print(triint)
# getTriList(intersectinglist)
# print(identifyIntersectType(triangles[:runLength])[:5])
print("\n\n\n********************************************")
# tri_groups = getTriangleGroups(triint, intersecting)
# np.save("triangle_groups", tri_groups)

tri_groups = np.load("triangle_groups.npy")
# print(tri_groups)
# vtk_tri_groups(tri_groups)
sortGroups(tri_groups)
error_groups = np.load("./fixtriNumpyArrays/error_groups.npy")
four_groups = np.load("./fixtriNumpyArrays/four_groups.npy")
# print(four_groups)
multi_groups = np.load("./fixtriNumpyArrays/multi_groups.npy")
# getFourList(four_groups)
fourList = np.load("4list.npy")
# print(fourList)
writeFileWithTriangles(fourList, "intersectingQuads.vtk")
# print(np.array([8185, 52344]) in intersecting[:, len(intersecting)])


#print(lineErrorSearch(multi_groups[1]))
#print(lineErrorSearch(multi_groups[2]))

#multi_triangle_groups(multi_groups)
# print("FOUR")
# print(four_groups)


#print(triint)
# print(intersecting)
# print("randomrandomrandom")
# print(triangles[406])
# print(triangles[407])


""" testing for fixing four groups """
# new_tris = fixFourGroups()

# new_tris = fixQuad([[406,407], [2539,2538]])
# print(triangles[406])
# print(triangles[407])
# print(triangles[2539])
# print(triangles[2538])
# print(new_tris)
# for group in four_groups:
#     stringer=""
#     for triangle in group:
#         for point in triangles[triangle]:
#             stringer += str(point) + " "
#         stringer += ", "
#     print(stringer)
# print(fourList)
# createFixedVTK(new_tris, "fixed.vtk")
#print(intersectinglist[0:100])
# finalwriter("final.vtk")