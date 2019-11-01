import pymesh2 as pymesh
import numpy as np
# from .meshio import form_mesh
import vtk
from vtk import vtkPolyDataReader
from vtk.util import numpy_support as VN    
import os
from operator import itemgetter

# ****************************** FUNCTIONS *****************************
def getPoints(data):
    arr = []
    #print(data.GetNumberOfPoints())
    for i in range(data.GetNumberOfPoints()):
        #print(i)
        tupledata = data.GetPoint(i)
        #print(tupledata)
        arr.append(list(tupledata))

    np_arr = np.array(arr)
    return np_arr

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
            arr.append([int(point1), int(point2), int(point3)])
    
    return np.array(arr)

# ************************** RUNTIME ************************************

reader = vtk.vtkPolyDataReader()
reader.SetFileName('vessel_triangles.vtk')
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

data = reader.GetOutput()


points = getPoints(data)

triangles = getTriangles(data)

mesh = pymesh.form_mesh(points, triangles)
cleanmesh = pymesh.detect_self_intersection(mesh)
pymesh.save_mesh(filename, mesh)

print(cleanmesh.vertices)
print(cleanmesh.faces)
