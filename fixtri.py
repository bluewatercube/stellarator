import numpy as np
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter

class TriPair:

    def __init__(self, tri1, tri2, inter1, inter2):
        self.tri1 = tri1 #tuple of point indexes(int)
        self.tri2 = tri2
        self.inter1 = inter1 #tuple of triangle indexes
        self.inter2 = inter2

# *************************** PARSE DATA FUNCTIONS ********************************

#parses file for points
def getPoints(data):
    arr = []
    #print(data.GetNumberOfPoints())
    for i in range(data.GetNumberOfPoints()):
        tupledata = data.GetPoint(i)
        arr.append(tupledata)

    np_arr = np.array(arr)
    return np_arr

#parses file for triangles
def getTriangles(data):
    arr = []
    
    cellArray = data.GetPolys() #this is vtkCellArray
    cellArray.InitTraversal()
    idList = vtk.vtkIdList()

    for i in range(0,cellArray.GetNumberOfCells()):
        cellArray.GetNextCell(idList)
        if idList.GetNumberOfIds()==3:
            point1 = idList.GetId(0)
            point2 = idList.GetId(1)
            point3 = idList.GetId(2)
            arr.append((int(point1), int(point2), int(point3)))
    
    return arr

# returns a list of groups of triangles that are intersecting together
# format: [(group), (group), (group)]
def getTriangleGroups(data, intersecting):
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

# ************************ TRIANGLE FIXING FUNCTIONS ***********************
#returns two lists as a tuple: (line over triangles, node over triangles)
def identifyIntersectType(triarr):
    """
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
        
#takes in groups from getTriangleGroups (list of lists)
def vtk_tri_groups(tri_groups):   
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
                

#group is list of triangle indeces
def writeGroupFile(group, filename):
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

def error_triangle_groups(tri_groups):
    # trigroup2 = []
    # for group in tri_groups:
    #     if len(group)==2 or len(group)==1:
    #         trigroup2.append(group)
    # vtk_tri_groups(trigroup2)
    pass

def four_triangle_groups(four_groups):

    for group in four_groups:
        print(lineErrorSearch(group))

def multi_triangle_groups(multi_groups):    
    for group in multi_groups:
        print(lineErrorSearch(group))
    
    #check
    # [(set),(set),(set)]

    # for x in range(len(triarr)):
    #     for y in range(x+1, len(triarr)):
    #         pair1 = triarr[x]
    #         pair2 = triarr[y]
    #         for intersection1 in intersections1:
    #             for intersection2 in intersections2:
    #                 if len(list(set(intersection1).intersection(set(intersection2)))) >= 2:

#quad is list of 2 tuples
def fixQuad(quad):
    #ex: A[BC] and [BC]D become [AD]B and [AD]C
    #determine which double is on inside surface --> indoub
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

def findDoubles(group):
    doubles = []
    for x in range(len(group)):
        for y in range(x+1, len(group)):
            if len(list(set(triangles[group[x]]).intersection(set(triangles[group[y]])))) == 2:
                doubles.append([group[x], group[y]])
    return doubles

def lineErrorSearch(group):
    doubles = findDoubles(group)
    quads = []
    counter = 0
    for x in range(len(doubles)):
        for y in range(x+1, len(doubles)):
            if doubles[x][0] != doubles[y][0] and \
                doubles[x][1] != doubles[y][0] and \
                doubles[x][0] != doubles[y][1] and \
                doubles[x][1] != doubles[y][1] and \
                [doubles[x][0], doubles[y][0]] in intersectinglist and \
                [doubles[x][1], doubles[y][0]] in intersectinglist and \
                [doubles[x][0], doubles[y][1]] in intersectinglist and \
                [doubles[x][1], doubles[y][1]] in intersectinglist:
                if counter < 100:
                    counter+=1
                quads.append([(doubles[x][0], doubles[x][1]), (doubles[y][0], doubles[y][1])])
    return quads
    
        
#returns list of triangle tuples that given tri intersects with
def getInter(tri, triarr):
    interarr = []
    for pair in triarr:
        if tri == pair[0]:
            interarr.append(pair[1])
        elif tri == pair[1]:
            interarr.append(pair[0])
    
    return interarr

def sortGroups(tri_groups):
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

def fixFourGroups():
    #print four_groups
    all_new_tris = []
    quads = []
    for group in four_groups:
        quad = lineErrorSearch(group)
        if quad != []:
            quads += quad
    print(quads)
    
    for quad in quads:
        new_tris = fixQuad(quad)
        all_new_tris += new_tris
    
    return all_new_tris


#new_triangles is a list of tuples of point indeces
def createFixedVTK(new_triangles, filename):
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

# ************************ GET DATA ***********************
reader = vtk.vtkPolyDataReader()
reader.SetFileName('vessel_triangles.vtk')
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

data = reader.GetOutput()

points = getPoints(data)
triangles = getTriangles(data)

# get intersecting triangles
intersecting = np.load("intersectingTriangles.npy")
intersectinglist = intersecting.tolist()
triint = np.load("triint.npy")
# ************************ RUNNER **************************

# getTriList(intersecting)
# print(identifyIntersectType(triangles[:runLength])[:5])
print("\n\n\n********************************************")
# tri_groups = getTriangleGroups(triint, intersecting)
# np.save("triangle_groups", tri_groups)

tri_groups = np.load("triangle_groups.npy")
#vtk_tri_groups(tri_groups)
# sortGroups(tri_groups)

error_groups = np.load("./fixtriNumpyArrays/error_groups.npy")
four_groups = np.load("./fixtriNumpyArrays/four_groups.npy")
multi_groups = np.load("./fixtriNumpyArrays/multi_groups.npy")



# print(np.array([8185, 52344]) in intersecting[:, len(intersecting)])


print(lineErrorSearch(multi_groups[1]))
print(lineErrorSearch(multi_groups[2]))

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
# print(new_tris)
# # new_tris = fixQuad([[406,407], [2539,2538]])
# createFixedVTK(new_tris, "fixed.vtk")