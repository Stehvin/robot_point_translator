import math
import boothSetupInfo

def execute(oldBooth, newBooth, newCurUtoolNum, oldPnts):
        """Execute the point conversion program.
        """
        # get parameters from excel document
        oldUtool, oldTbl, oldEuler, oldDustColLoc, oldDustColAng = \
                  boothSetupInfo.boothInfo(oldBooth)
        newUtool, newTbl, newEuler, newDustColLoc, newDustColAng = \
                  boothSetupInfo.boothInfo(newBooth)
        newCurUtool = boothSetupInfo.curUtool(newBooth, newCurUtoolNum)

        # iterate through points, translating each one
        print()
        newPnts = {}
        newOppPnts = {}
        for pntNum, pnt in oldPnts.items():

                # convert original point to gun-tip point
                oldCurUtool = boothSetupInfo.curUtool(oldBooth, pnt['ut'])
                oldCurUframe = boothSetupInfo.curUframe(oldBooth, pnt['uf'])
                if oldCurUframe != {'x':0, 'y':0, 'z':0, 'w':0, 'p':0, 'r':0}:

                        # convert point from uframe to world frame
                        pnt = frame_to_world(pnt, oldCurUframe)

                if oldCurUtool != oldUtool:

                        # convert point from original utool to nozzle utool
                        pnt = utool_to_world(pnt, oldCurUtool)
                        pnt = world_to_utool(pnt, oldUtool)
                
                # convert z-angle ("roll" angle)
                if oldDustColLoc == newDustColLoc:
                        opposite = False
                        zAngleDeg = pnt['r']
                elif (oldDustColLoc == 'R' and newDustColLoc == 'L') or \
                     (oldDustColLoc == 'L' and newDustColLoc == 'R'):
                        opposite = True
                        zAngleDeg = -pnt['r']
                elif newDustColLoc == 'C' or oldDustColLoc == 'C':
                        opposite = False
                        zAngleDeg = newDustColAng
                
                # find new point
                newPnts[pntNum] = spfTranslate(pnt, oldTbl, newTbl, oldEuler,
                                               newEuler, zAngleDeg, opposite)
                
                # convert new gun-tip point to faceplate point
                # (or different utool point)
                if newCurUtool != newUtool:
                        newPnts[pntNum] = utool_to_world(newPnts[pntNum],
                                                         newUtool)
                        newPnts[pntNum] = world_to_utool(newPnts[pntNum],
                                                         newCurUtool)
        return dict(newPnts)

def eul2Mat(w, p, r):
        """Converts Euler angles to 3x3 rotation matrix.
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
        """Converts 3x3 rotation matrix to Euler angles.
        """
        r = math.degrees(math.atan2(M[1][0], M[0][0]))
        p = math.degrees(-math.asin(M[2][0]))
        w = math.degrees(math.atan2(M[2][1], M[2][2]))

        return w, p, r

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

def matMult(A, B):
        """Multiplies matricies A and B (A*B = C).
        Full disclosure, I did not write this function.
        Taken from:
        https://stackoverflow.com/questions/10508021/
        matrix-multiplication-in-python
        """
        zip_b = zip(*B)
        zip_b = list(zip_b)
        return [[sum(ele_a*ele_b for ele_a, ele_b in zip(
                row_a, col_b)) for col_b in zip_b] for row_a in A]

def convertAngles(oldPnt, oldEuler, newEuler):
        """Converts Euler angles from a point in one booth to another booth.
        """
        # find old booth's origin-to-gun rotation matrix
        A = eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r'])
        B = eul2Mat(oldEuler['w'], oldEuler['p'], oldEuler['r'])
        C = matMult(A, B)

        # find new booth's faceplate-to-gun rotation matrix
        newB = eul2Mat(newEuler['w'], newEuler['p'], newEuler['r'])

        # find new booth's origin-to-faceplate euler angles     
        newA = matMult(C, rotMatInv(newB))
        w, p, r = mat2Eul(newA)
        return w, p, r

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

def world_to_utool(worldPnt, utool):
        """Converts a world frame point to a utool point.
        """
        # set up world-coord matrix and translate by utool xyz setpoints
        rotM = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        M4x4 = rotTransMat(rotM, worldPnt)
        utoolPnt = multMatVec(M4x4, utool)

        # change world-coord angles to utool angles
        rotA = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
        rotB = eul2Mat(utool['w'], utool['p'], utool['r'])
        utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(matMult(rotA,
                                                                      rotB))
        return utoolPnt

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

def quadratic(a, b, c):
        """Find the solutions to a quadratic formula,
           where 0 = (a)x^2 + (b)x + (c)
        """
        term1 = (-1) * b
        term2 = math.sqrt((b**2) - (4*a*c))
        term3 = 2 * a

        solutions = []  
        signs = [1, -1]
        for sign in signs:
                try:
                        x = (term1 + (sign * term2)) / term3
                except:
                        x = None
                solutions.append(x)
        return solutions

def spfTranslate(oldPnt, oldTbl, newTbl, oldEuler, newEuler, zAngleDeg,
                 opposite=False):
        """Translate point from old booth to new booth.
        """
        # get original booth center offset angle (alpha)
        oldFi = math.atan2(oldTbl['y'] - oldPnt['y'], oldTbl['x'] - oldPnt['x'])
        oldAlpha = math.radians(oldPnt['r']) - oldFi
        
        # change old z-Angle if opposite points are necessary
        if opposite:
                oldPnt['r'] = math.degrees(oldFi - oldAlpha)
                oldAlpha = -oldAlpha
        
        # initialize newPnt and set angles
        newPnt = {'x': 0, 'y': 0, 'z': 0, 'w': 0, 'p': 0, 'r': 0}
        newPnt['w'], newPnt['p'], newPnt['r'] = convertAngles(oldPnt, oldEuler,
                                                              newEuler)
        newPnt['r'] = zAngleDeg

        # get z-height and set E1 value
        newPnt['z'] = oldPnt['z'] - oldTbl['z'] + newTbl['z']
        try:
                newPnt['e1'] = oldPnt['e1']
        except:
                newPnt['e1'] = 91.588
        
        # x and y system of equations
        angVar = math.tan(math.atan2(oldTbl['y'] - oldPnt['y'], oldTbl['x'] - \
                                     oldPnt['x']) - \
                          math.radians(oldPnt['r']) + math.radians(zAngleDeg))
        distSq = ((oldPnt['x'] - oldTbl['x'])**2) + \
                 ((oldPnt['y'] - oldTbl['y'])**2)
        
        # execute quadratic equation solver
        a = 1 + (angVar**2)
        b = ((-2)*(newTbl['x'])) - (2*(angVar**2)*newTbl['x'])
        c = (newTbl['x']**2) + ((angVar**2)*(newTbl['x']**2)) - distSq
        xPnts = quadratic(a, b, c)
        
        for x in xPnts:
                if x:
                        # solve for y
                        y = newTbl['y'] - (angVar * newTbl['x']) + (angVar * x)

                        # compare new and old booth center offset angles
                        newAlpha = math.radians(newPnt['r']) - \
                                   math.atan2(newTbl['y'] - y, newTbl['x'] - x)
                        if math.isclose(newAlpha, oldAlpha, rel_tol=1e-06) or \
                           math.isclose(newAlpha - (2*math.pi), oldAlpha,
                                        rel_tol=1e-06):
                                newPnt['x'] = x
                                newPnt['y'] = y
                                break
        return newPnt
