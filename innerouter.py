import numpy as np
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter

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


# searchList is a list initially containing original point
# and group is a list initially empty
def groupFinder(start):
    global group
    global triangleGroup

    explored = []
    queue = [start]
    while queue:
        node = queue.pop(0)
        explored.append(node)
        for index in range(len(triangles)):
            triangle = triangles[index]
            if node in triangle:
                if index not in triangleGroup:
                    triangleGroup.append(index)
                for point in triangle:
                    if point != node and (point not in explored):
                        explored.append(point)
                        queue.append(point)
    return explored
        

    # searchitem = searchList[0]
    # # print("group")
    # print(len(group))
    # # print("search")
    # # print(searchList)
    # # print()
    # trip = False
    # for index in range(len(triangles)):
    #     triangle = triangles[index]
    #     if searchitem in triangle:
    #         if index not in triangleGroup:
    #             triangleGroup.append(index)
    #         for point in triangle:
    #             if point != searchitem and (point not in group):
                    # trip = True
    #                 group.append(point)
    #                 searchList.append(point)
    # del searchList[0]
    
    # if trip:
    #     return
    # else: 
    #     return groupFinder(searchList)
        

def writeOneFile(group, fileName):
    fileOut = open(fileName, "w")
    fileOut.write("# vtk DataFile Version 3.0\n")
    fileOut.write("group Triangles\n")
    fileOut.write("ASCII\n")
    fileOut.write("DATASET POLYDATA\n")

    numPoints = len(points)
    numCells = len(group) #num triangles
    cellListSize = numCells * 4

    fileOut.write("POINTS " + str(numPoints) + " float\n")
    #print points
    for point in points:
        fileOut.write(str(point[0])+" "+str(point[1])+" "+str(point[2])+"\n")

    fileOut.write("POLYGONS " + str(numCells) + " " + str(cellListSize)+ "\n")
    for triangle in group:
        fileOut.write("3 " + " ".join(map(str,triangles[triangle]))+"\n")

                 
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
    



reader = vtk.vtkPolyDataReader()
reader.SetFileName('precisionCutTri.vtk')
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

data = reader.GetOutput()

points = getPoints(data)
triangles = getTriangles(data)



group = []
triangleGroup = []
group = groupFinder(493)
# print(triangles)
# print("group")
# print(group)
print("triangleGroup")
print(triangleGroup)
print(len(triangleGroup))

fileDirectory = "./triangleSearch"
if not os.path.exists(fileDirectory):
    os.makedirs(fileDirectory)
for x in range(int(len(triangleGroup)/75)):
    writeFileWithTriangles(triangleGroup[0 : 50 + 75*x], "./triangleSearch/spread_"+str(x)+".vtk")
