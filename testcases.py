# ************************ TEST CASES ***********************

# intersecting with piercing triangle
# points = np.array([(0,0,0), (0,0,5), (5,0,-1), (0,1,1), (5,1,1), (3,-2,1)])
# 
# intersecting with linking triangle
# points = np.array([(0,0,0), (5,-3,5), (5,3,5), (5,0,0), (0,0,5), (3,3,50)])
# 
# intersecting with piercing triangle
# points = np.array([(0,0,0), (0,5,0), (8,0,0), (6,8,3), (6,8,-2),(6,-4,-2)])
# 
# case 3:  intersecting
# points = np.array([(0,0,0), (2.1,0,0), (1.5, 48, 0), (1,1,1), (1.2,-1,3), (2,-1,-50)])
# 
# another non intersecting triangle
# points = np.array([(0,0,0), (10,0,0), (0,0,10), (9, -1, 9), (7,1,7), (9,1,9)])
#   
# another non intersecting triangle
# points = np.array([(1,1,1), (1,3,4), (1,5,2), (0,2,5)), (2,3,5), (3,5,5)])
#
# intersecting by 0.2
# points = np.array([(0,0,0), (2.1,0,1), (1.5, 0,4), (1,1,1), (1.2,3,-1), (1,-1,0)])

# triangle1 = np.array([0,1,2])
# triangle2 = np.array([3,4,5])
# triangles = [triangle1, triangle2]
# print(triangleIntersect(triangle1, triangle2))

# boxes = []
# for y in range(2):
#     boxes.append(boundingBox(triangles[y]))

# if twoPointcheck(triangles[0], triangles[1]) and \
#             boundingBoxCheck(boxes[0], boxes[1]) and \
#             planeCheck(triangles[0], triangles[1]) and \
#             triangleIntersect(triangles[0], triangles[1]):
#             print("yay")
# else:
#     print("sad")