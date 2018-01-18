import math

# input points/utools here
worldPnt = {'x': 535.533, 'y': 659.809, 'z': 800.271, 'w': 179.232,
            'p': 1.968, 'r': 136.681}
utoolPnt = {'x': 884.207, 'y': 610.492, 'z': 699.190, 'w': 179.084,
            'p': 0.077, 'r': 49.245}
uframePnt = {'x': -323.448, 'y': 492.870, 'z': 1252.489, 'w': -159.618,
            'p': -19.331, 'r': 133.385}
bothPnt = {'x': -70.911, 'y': 443.500, 'z': 991.701, 'w': -157.493,
            'p': 18.975, 'r': 53.217}

uframe = {'x': 1437.991, 'y': 167.013, 'z': -126.546, 'w': -0.006,
          'p': -29.753, 'r': 0.013}
utool = {'x': -283.869, 'y': 201.826, 'z': 113.61, 'w': 1.085,
         'p': -0.757, 'r': 87.435}

def eul2Mat(w, p, r):
        # works
        # converts euler angles to 3x3 rotation matrix
        
        # initialize matrix
        M = [[0, 0, 0], 
             [0, 0, 0], 
             [0, 0, 0]]
        
        # convert euler angles from degrees to radians
        w = math.radians(w)
        p = math.radians(p)
        r = math.radians(r)
        
        # calculate matrix elements
        M[0][0] = math.cos(p) * math.cos(r)
        M[0][1] = (math.sin(w) * math.sin(p) * math.cos(r)) - (math.cos(w) * math.sin(r))
        M[0][2] = (math.sin(w) * math.sin(r)) + (math.cos(w) * math.sin(p) * math.cos(r))
        M[1][0] = math.cos(p) * math.sin(r)
        M[1][1] = (math.cos(w) * math.cos(r)) + (math.sin(w) * math.sin(p) * math.sin(r))
        M[1][2] = (math.cos(w) * math.sin(p) * math.sin(r)) - (math.sin(w) * math.cos(r))
        M[2][0] = -math.sin(p)
        M[2][1] = math.sin(w) * math.cos(p)
        M[2][2] = math.cos(w) * math.cos(p)

        return M

def mat2Eul(M):
        """Converts 3x3 rotation matrix to euler angles.
        """
        r = math.degrees(math.atan2(M[1][0], M[0][0]))
        p = math.degrees(-math.asin(M[2][0]))
        w = math.degrees(math.atan2(M[2][1], M[2][2]))

        return w, p, r

def matMult(A, B):
        """Multiplies matricies A and B (A*B = C).
        Full disclosure, I did not write this function.
        Stolen from: https://stackoverflow.com/questions/10508021/matrix-multiplication-in-python
        """

        zip_b = zip(*B)
        zip_b = list(zip_b)
        return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) for col_b in zip_b] for row_a in A]

def rotMatInv(M):
        # inverts the 3x3 rotation matrix 'M'
        # the inverse of a rotation matrix is its transpose

        # initialize output matrix
        MI = [[0, 0, 0],
              [0, 0, 0],
              [0, 0, 0]]

        # calculate matrix elements
        MI[0][0] = M[0][0]
        MI[0][1] = M[1][0]
        MI[0][2] = M[2][0]
        MI[1][0] = M[0][1]
        MI[1][1] = M[1][1]
        MI[1][2] = M[2][1]
        MI[2][0] = M[0][2]
        MI[2][1] = M[1][2]
        MI[2][2] = M[2][2]

        return MI

def rotTransMat(rotM, translation):
        # turns a rotation matrix and a translation vector into
        # a 4x4 rotation and translation matrix
        
        # create output matrix
        M = [[rotM[0][0], rotM[0][1], rotM[0][2], translation['x']],
             [rotM[1][0], rotM[1][1], rotM[1][2], translation['y']],
             [rotM[2][0], rotM[2][1], rotM[2][2], translation['z']],
             [0, 0, 0, 1]]
        return M

def rot_matrix_4x4(rotM):
        # changes 3x3 rotation matrix to a 4x4 matrix
        
        # create output matrix
        M = [[rotM[0][0], rotM[0][1], rotM[0][2], 0],
             [rotM[1][0], rotM[1][1], rotM[1][2], 0],
             [rotM[2][0], rotM[2][1], rotM[2][2], 0],
             [0, 0, 0, 1]]
        return M

def multMatVec(M, vector):
        # multiplies a 4x4 matrix and a 4x1 vector

        # initialize output vector
        vector1 = {'x': 0, 'y': 0, 'z': 0}
        
        # calculate vector elements
        vector1['x'] = (M[0][0]*vector['x']) + (M[0][1]*vector['y']) + (M[0][2]*vector['z']) + M[0][3]
        vector1['y'] = (M[1][0]*vector['x']) + (M[1][1]*vector['y']) + (M[1][2]*vector['z']) + M[1][3]
        vector1['z'] = (M[2][0]*vector['x']) + (M[2][1]*vector['y']) + (M[2][2]*vector['z']) + M[2][3]
        return vector1

def printMat(M):
        # prints matrix in readable format
        
        for row in M:   
                print(row)


'''
**********

End of helper functions.
Below are point conversion functions that utilize the helper functions.

**********
'''

def world_to_utool(worldPnt, utool):
        # converts a world frame point to a utool point
        
        rotM = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        M4x4 = rotTransMat(rotM, worldPnt)
        gunPnt = multMatVec(M4x4, utool)

        # convert world-coord angles to utool angles
        rotA = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        rotB = eul2Mat(utool['w'], utool['p'], utool['r'])
        gunPnt['w'], gunPnt['p'], gunPnt['r'] = mat2Eul(matMult(rotA, rotB))
        return gunPnt

def utool_set_points(worldPnt, utoolPnt):
        # finds accurate utool set points given gun-tip point and faceplate point
        
        rotM = rotMatInv(eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r']))
        rotM2 = rot_matrix_4x4(rotM)
        transM = [[1, 0, 0, -worldPnt['x']], 
                  [0, 1, 0, -worldPnt['y']], 
                  [0, 0, 1, -worldPnt['z']],
                  [0, 0, 0, 1]]
        M = matMult(rotM2, transM)
        utool = multMatVec(M, utoolPnt)
        return utool

def utool_to_world(utoolPnt, utool):
        # converts a utool point to a world coordinate point
        
        utoolCopy = dict(utool)

        # change utool angles to world angles
        rotA = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        rotB = eul2Mat(utoolCopy['w'], utoolCopy['p'], utoolCopy['r'])
        utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(
                                                        matMult(
                                                        rotA, rotMatInv(rotB)))

        # set up point translation (utool setpoints)
        for key in utoolCopy:
                utoolCopy[key] *= -1
                
        # set up utool matrix
        rotM = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        M4x4 = rotTransMat(rotM, utoolPnt)

        worldPnt = multMatVec(M4x4, utoolCopy)
        worldPnt['w'] = utoolPnt['w']
        worldPnt['p'] = utoolPnt['p']
        worldPnt['r'] = utoolPnt['r']
        return worldPnt

def frame_to_world(uframePnt, uframe):
        # converts uframe point to a world frame point

        # convert xyz values to world frame
        rotA = eul2Mat(uframe['w'], uframe['p'], uframe['r'])
        four = rotTransMat(rotA, uframe)
        worldPnt = multMatVec(four, uframePnt)

        # convert wpr angles to world frame angles
        rotB = eul2Mat(uframePnt['w'], uframePnt['p'], uframePnt['r'])
        worldPnt['w'], worldPnt['p'], worldPnt['r'] = mat2Eul(matMult(rotA,
                                                                      rotB))
        return worldPnt

def world_to_frame(worldPnt, uframe):
        """Converts a world point to a user frame point.
        """
        # set up translation matrix (uframe setpoints)
        transM = [[1, 0, 0, -uframe['x']], 
                  [0, 1, 0, -uframe['y']], 
                  [0, 0, 1, -uframe['z']],
                  [0, 0, 0, 1]]
                
        # convert xyz values to uframe
        rotA3 = rotMatInv(eul2Mat(uframe['w'], uframe['p'], uframe['r']))
        rotA4 = rot_matrix_4x4(rotA3)
        M = matMult(rotA4, transM)
        ufPnt = multMatVec(M, worldPnt)

        # convert wpr angles to uframe
        rotB = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        ufPnt['w'], ufPnt['p'], ufPnt['r'] = mat2Eul(matMult(rotA3, rotB))
        return ufPnt

def convertAngles(oldPnt, newEuler, opposite=False):
        """Converts Euler angles from a point in one booth to another booth.
           Old point Euler angles have already been changed to UGO angles.
        """
        # convert old UGO angles to a rotation matrix
        C = eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r'])

        # for opposite booths, find the z-angle (pointing straight down x-axis
        # is 0 degrees) for the old booth, flip it to the other side (by
        # reversing its sign), and recalculate UGO angle rotation matrix
        if opposite:
                oldZAngle = getZAngle(C)
                zRotM = eul2Mat(0, 0, -2 * oldZAngle)
                C = matMult(zRotM, C)        
        
        # find new booth's faceplate-to-UGO rotation matrix
        newB = eul2Mat(newEuler['w'], newEuler['p'], newEuler['r'])

        # find new booth's origin-to-faceplate euler angles
        newA = matMult(C, rotMatInv(newB))
        w, p, r = mat2Eul(newA)
        return w, p, r

def normalize_Z(eulerMat):
        # finds the z-angle when the gun is pointing straight down the x-axis
        xMat = eul2Mat(0, 90, 0)

        rotM = matMult(xMat, rotMatInv(eulerMat))

        w, p, r = mat2Eul(rotM)
        print(mat2Eul(rotM))

        return r

def getZAngle(M):
        """Returns the z-angle of the projection of a 3D vector on the 2D X-Y
           plane. Rotation matrix "M" must be in UGO angles. Returns angle in
           radians.
        """
        # multiply unit vector [0, 0, 1] and rotation matrix to get X-Y
        # projection points
        x = M[0][2]
        y = M[1][2]

        # find tangent of X-Y point
        return math.atan2(y, x)

def UGO(pnt, Euler):
    """Convert a point's angles to Universal Gun Orientation (UGO) angles.
    """
    pnt1 = dict(pnt)
    rotM = matMult(eul2Mat(pnt1['w'], pnt1['p'], pnt1['r']),
            eul2Mat(Euler['w'], Euler['p'], Euler['r']))
    pnt1['w'], pnt1['p'], pnt1['r'] = mat2Eul(rotM)
    return pnt1

def convert25ttcp(pnt, oldTbl):
        """Converts the turntable centerpoint (TTCP) x-value, when the old
           booth is 25. The new TTCP x-value is based upon the E1-value of the
           current point.
        """
        oldTbl1 = dict(oldTbl)
        oldTbl1['x'] = 1012.072 - (2*(-582.815 - pnt['e1']))
        return oldTbl1

utPnt = world_to_utool(worldPnt, utool)
print(world_to_frame(utPnt, uframe))

#print(utool_set_points(worldPnt, utoolPnt))
