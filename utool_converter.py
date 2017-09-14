import math

# input points/utools here
worldPnt = {'x': 562.729, 'y': -909.742, 'z': -237.066, 'w': -88.849,
            'p': -5.162, 'r': 4.281}
uframePnt = {'x': 541.076, 'y': -926.536, 'z': 211.193, 'w': -110.060,
            'p': -17.876, 'r': 7.497}
utoolPnt = {'x': 562.727, 'y': -909.742, 'z': -237.069, 'w': -88.846,
            'p': -6.956, 'r': 4.245}

uframe = {'x': 1.861, 'y': 0, 'z': -5.478, 'w': 20, 'p': 12.369, 'r': 3.569}
utool = {'x': 0, 'y': 0, 'z': 0, 'w': 0, 'p': 0, 'r': -1.8}

def eul2Mat(w, p, r):
        """Converts euler angles to 3x3 rotation matrix.
        """
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
        M[0][1] = (math.sin(w) * math.sin(p) * math.cos(r)) - \
                  (math.cos(w) * math.sin(r))
        M[0][2] = (math.sin(w) * math.sin(r)) + \
                  (math.cos(w) * math.sin(p) * math.cos(r))
        M[1][0] = math.cos(p) * math.sin(r)
        M[1][1] = (math.cos(w) * math.cos(r)) + \
                  (math.sin(w) * math.sin(p) * math.sin(r))
        M[1][2] = (math.cos(w) * math.sin(p) * math.sin(r)) - \
                  (math.sin(w) * math.cos(r))
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
        Taken from: https://stackoverflow.com/questions/10508021/
        matrix-multiplication-in-python
        """
        zip_b = zip(*B)
        zip_b = list(zip_b)
        return [[sum(ele_a*ele_b for ele_a, ele_b in zip(
                row_a, col_b)) for col_b in zip_b] for row_a in A]

def rotMatInv(M):
        """Inverts the 3x3 rotation matrix 'M'.
           The inverse of a rotation matrix is its transpose.
        """
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
        """Turns a rotation matrix and a translation vector into
           a 4x4 rotation and translation matrix.
        """
        # create output matrix
        M = [[rotM[0][0], rotM[0][1], rotM[0][2], translation['x']],
             [rotM[1][0], rotM[1][1], rotM[1][2], translation['y']],
             [rotM[2][0], rotM[2][1], rotM[2][2], translation['z']],
             [0, 0, 0, 1]]
        return M

def rot_matrix_4x4(rotM):
        """Changes a 3x3 rotation matrix to a 4x4 matrix.
        """        
        # create output matrix
        M = [[rotM[0][0], rotM[0][1], rotM[0][2], 0],
             [rotM[1][0], rotM[1][1], rotM[1][2], 0],
             [rotM[2][0], rotM[2][1], rotM[2][2], 0],
             [0, 0, 0, 1]]
        return M

def multMatVec(M, vector):
        """Multiplies a 4x4 matrix and a 4x1 vector.
        """
        # initialize output vector
        vector1 = {'x': 0, 'y': 0, 'z': 0}
        
        # calculate vector elements
        vector1['x'] = (M[0][0]*vector['x']) + (M[0][1]*vector['y']) + \
                       (M[0][2]*vector['z']) + M[0][3]
        vector1['y'] = (M[1][0]*vector['x']) + (M[1][1]*vector['y']) + \
                       (M[1][2]*vector['z']) + M[1][3]
        vector1['z'] = (M[2][0]*vector['x']) + (M[2][1]*vector['y']) + \
                       (M[2][2]*vector['z']) + M[2][3]
        return vector1

def printMat(M):
        """Prints a matrix in readable format.
        """
        for row in M:   
                print(row)


'''
**********

End of helper functions.
Below are point conversion functions that utilize the helper functions.

**********
'''

def world_to_utool(worldPnt, utool):
        """Converts a world frame point to a utool point.
        """
        # set up world-coord matrix and translate by utool xyz setpoints
        rotM = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        M4x4 = rotTransMat(rotM, worldPnt)
        gunPnt = multMatVec(M4x4, utool)

        # convert world-coord angles to utool angles
        rotA = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        rotB = eul2Mat(utool['w'], utool['p'], utool['r'])
        gunPnt['w'], gunPnt['p'], gunPnt['r'] = mat2Eul(matMult(rotA, rotB))
        return gunPnt

def utool_set_points(worldPnt, utoolPnt):
        """Finds accurate utool set points given a gun-tip point and a
           faceplate point.
        """
        # set up rotation matrix and translation matrix
        rotM = rotMatInv(eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r']))
        rotM2 = rot_matrix_4x4(rotM)
        transM = [[1, 0, 0, -worldPnt['x']], 
                  [0, 1, 0, -worldPnt['y']], 
                  [0, 0, 1, -worldPnt['z']],
                  [0, 0, 0, 1]]

        # multiply matricies and translate utool point
        M = matMult(rotM2, transM)
        utool = multMatVec(M, utoolPnt)
        return utool

def utool_to_world(utoolPnt, utool):
        """Converts a utool point to a world coordinate point.
        """
        # reverse utool setpoint value signs
        utoolCopy = dict(utool)
        for key in utoolCopy:
                utoolCopy[key] *= -1

        # change utool angles to world angles
        rotA = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        rotB = eul2Mat(utoolCopy['w'], utoolCopy['p'], utoolCopy['r'])
        utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(matMult(rotA,
                                                                      rotB))

        # set up utool matrix
        rotM = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        M4x4 = rotTransMat(rotM, utoolPnt)

        worldPnt = multMatVec(M4x4, utoolCopy)
        worldPnt['w'] = utoolPnt['w']
        worldPnt['p'] = utoolPnt['p']
        worldPnt['r'] = utoolPnt['r']
        return worldPnt

def frame_to_world(uframePnt, uframe):
        """Converts uframe point to a world frame point.
        """
        # convert xyz values to world frame
        rotA = eul2Mat(uframe['w'], uframe['p'], uframe['r'])
        four = rotTransMat(rotA, uframe)
        worldPnt = multMatVec(four, uframePnt)

        # convert wpr angles to world frame angles
        rotB = eul2Mat(uframePnt['w'], uframePnt['p'], uframePnt['r'])
        worldPnt['w'], worldPnt['p'], worldPnt['r'] = mat2Eul(matMult(rotA,
                                                                      rotB))
        return worldPnt

print(frame_to_world(uframePnt, uframe))
