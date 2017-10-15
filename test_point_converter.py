import math
import boothSetupInfo

def main():
        oldPnts = {2: {'uf': 1, 'ut': 1, 'x': -194.62, 'y': -194.41,
                       'z': 130.67, 'w': 179.08, 'p': 0.08, 'r': 54.64}}
        oldBooth = 6
        newBooth = 13
        newCurUtoolNum = 1
    
        newPnts = execute(oldBooth, newBooth, newCurUtoolNum, oldPnts)
        print(newPnts)

def execute(oldBooth, newBooth, newCurUtoolNum, oldPnts):
        """Execute the point conversion program.
        """
        # get parameters from excel document
        oldUtool, oldTbl, oldEuler, oldDustColLoc, oldDustColAng = \
                  boothSetupInfo.boothInfo(oldBooth)[0:5]
        newUtool, newTbl, newEuler, newDustColLoc, newDustColAng = \
                  boothSetupInfo.boothInfo(newBooth)[0:5]
        newCurUtool = boothSetupInfo.curUtool(newBooth, newCurUtoolNum)

        # iterate through points, translating each one
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
                        # and convert angles to UGO (Universal Gun Orientation)
                        # angles
                        pnt = utool_to_world(pnt, oldCurUtool)
                        pnt = world_to_utool(pnt, oldUtool)
                        print("Old UF0, Gun-tip UTool Coordinates:")
                        print((pnt['x'], pnt['y'], pnt['z'], pnt['w'],
                               pnt['p'], pnt['r']))
                        print()
                        pnt = UGO(pnt, oldEuler)

                # convert z-angle ("roll" angle)
                if oldDustColLoc == newDustColLoc:
                        center = False
                        opposite = False
                        zAngleDeg = getZAngle(eul2Mat(pnt['w'], pnt['p'],
                                                      pnt['r']))
                elif (oldDustColLoc == 'R' and newDustColLoc == 'L') or \
                     (oldDustColLoc == 'L' and newDustColLoc == 'R'):
                        center = False
                        opposite = True
                        zAngleDeg = -getZAngle(eul2Mat(pnt['w'], pnt['p'],
                                                       pnt['r']))
                elif newDustColLoc == 'C' or oldDustColLoc == 'C':
                        center = True
                        opposite = False
                        zAngleDeg = newDustColAng
                
                # find new point
                newPnts[pntNum] = spfTranslate(pnt, oldTbl, newTbl, oldEuler,
                                               newEuler, zAngleDeg, pntNum,
                                               center, opposite)
                
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
        https://stackoverflow.com/questions/10508021/matrix-
        multiplication-in-python
        """
        zip_b = zip(*B)
        zip_b = list(zip_b)
        return [[sum(ele_a*ele_b for ele_a, ele_b in zip(
                row_a, col_b)) for col_b in zip_b] for row_a in A]

def convertAngles(oldPnt, newEuler, zAngleDeg, center=False, opposite=False):
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
                
        # for center booths, take the same steps described above for opposite
        # booths, except the new z-angle is based on the new dust collector
        # location (zAngleDeg)
        if center:
                oldZAngle = getZAngle(C)
                zRotM = eul2Mat(0, 0, zAngleDeg - oldZAngle)
                C = matMult(zRotM, C)
        
        # find new booth's faceplate-to-UGO rotation matrix
        newB = eul2Mat(newEuler['w'], newEuler['p'], newEuler['r'])

        # find new booth's origin-to-faceplate euler angles
        newA = matMult(C, rotMatInv(newB))
        w, p, r = mat2Eul(newA)
        return w, p, r

def printMat(M):
        # prints matrix in readable format
        for row in M:   
                print(row)

def getZAngle(M):
        """Returns the z-angle of the projection of a 3D vector on the 2D X-Y
           plane. Rotation matrix "M" must be in UGO angles. Returns angle in
           degrees.
        """
        # multiply unit vector [0, 0, 1] and rotation matrix to get X-Y
        # projection points
        x = M[0][2]
        y = M[1][2]

        # find tangent of X-Y point
        return math.degrees(math.atan2(y, x))

def getZAngleWorld(pnt, Euler):
        """Performs same function as the getZAngle function, except the input
           point is a world-frame point.
        """
        # convert world-frame Euler angles to UGO rotation matrix
        pnt1 = UGO(pnt, Euler)
        M = eul2Mat(pnt1['w'], pnt1['p'], pnt1['r'])
        
        # multiply unit vector [0, 0, 1] and rotation matrix to get X-Y
        # projection points
        x = M[0][2]
        y = M[1][2]

        # find tangent of X-Y point
        return math.degrees(math.atan2(y, x))

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

        # assign E1 coordinate, if it exists
        try:
                utoolPnt['e1'] = worldPnt['e1']
        except:
                pass
        return utoolPnt

def utool_to_world(utoolPnt, utool):
        """Converts a utool point to a world coordinate point.
        """
        utoolCopy = dict(utool)

        # change utool angles to world angles
        rotA = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        rotB = rotMatInv(eul2Mat(utoolCopy['w'], utoolCopy['p'],
                                 utoolCopy['r']))
        utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(matMult(rotA,
                                                                      rotB))

        # reverse utool setpoint value signs
        for key in utoolCopy:
                utoolCopy[key] *= -1
        
        # set up utool matrix
        rotM = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
        M4x4 = rotTransMat(rotM, utoolPnt)

        worldPnt = multMatVec(M4x4, utoolCopy)
        worldPnt['w'] = utoolPnt['w']
        worldPnt['p'] = utoolPnt['p']
        worldPnt['r'] = utoolPnt['r']

        # assign E1 coordinate, if it exists
        try:
                worldPnt['e1'] = utoolPnt['e1']
        except:
                pass
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

def UGO(pnt, Euler):
        """Convert a point's angles to Universal Gun Orientation (UGO) angles.
        """
        pnt1 = dict(pnt)
        rotM = matMult(eul2Mat(pnt1['w'], pnt1['p'], pnt1['r']),
                eul2Mat(Euler['w'], Euler['p'], Euler['r']))
        pnt1['w'], pnt1['p'], pnt1['r'] = mat2Eul(rotM)
        return pnt1

def spfTranslate(oldPnt, oldTbl, newTbl, oldEuler, newEuler, zAngleDeg,
                 pntNum, center=False, opposite=False):
        """Translate point from old booth to new booth.
        """
        # get original booth center offset angle (alpha)
        oldFi = math.atan2(oldTbl['y'] - oldPnt['y'], oldTbl['x'] - oldPnt['x'])
        oldAlpha = math.radians(getZAngle(eul2Mat(oldPnt['w'], oldPnt['p'],
                                                  oldPnt['r']))) - oldFi

        # print old point values
        print("Old Point Values:")
        print("UGO Angles:")
        print(mat2Eul(eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r'])))
        print("Z-Height:")
        print(oldPnt['z'] - oldTbl['z'])
        print("X-Y Distance to Center of TT:")
        print(math.sqrt(((oldPnt['x'] - oldTbl['x'])**2) + \
                 ((oldPnt['y'] - oldTbl['y'])**2)))
        print("Offset Angle (Alpha):")
        print(oldAlpha)
        print("Coordinates:")
        print((oldPnt['x'], oldPnt['y'], oldPnt['z'], oldPnt['w'], oldPnt['p'],
               oldPnt['r']))
        print()
        
        # change old z-angle to point on other side of the table
        # if opposite points are necessary (reverse the sign of alpha)
        if opposite:
                oldAlpha = -oldAlpha
                print("Old Point Values (Opposite Change):")
                print("Offset Angle (Alpha):")
                print(oldAlpha)
                print()

        # initialize new point
        newPnt = {'x': 0, 'y': 0, 'z': 0, 'w': 0, 'p': 0, 'r': 0}

        # convert old UGO (Universal
        # Gun Orientation) angles to new UGO angles, then to new world angles
        print("Old Z-Angle:")
        print(getZAngle(eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r'])))
        newPnt['w'], newPnt['p'], newPnt['r'] = convertAngles(
            oldPnt, newEuler, zAngleDeg, center, opposite)
        print("New Z-Angle:")
        print(getZAngleWorld(newPnt, newEuler))
        print()

        # get z-height and set E1 value
        newPnt['z'] = oldPnt['z'] - oldTbl['z'] + newTbl['z']
        try:
                newPnt['e1'] = oldPnt['e1']
        except:
                newPnt['e1'] = 91.588

        # x and y system of equations
        angVar = math.tan(math.radians(zAngleDeg) - oldAlpha)
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
                        newAlpha = math.radians(zAngleDeg) - \
                                   math.atan2(newTbl['y'] - y, newTbl['x'] - x)
                        if math.isclose(newAlpha, oldAlpha, rel_tol=1e-06) or \
                           math.isclose(newAlpha - (2*math.pi), oldAlpha,
                                        rel_tol=1e-06) or \
                           math.isclose(newAlpha + (2*math.pi), oldAlpha,
                                        rel_tol=1e-06):
                                newPnt['x'] = x
                                newPnt['y'] = y
                                break

        # print new point values
        print("New Point Values:")
        print("UGO Angles:")
        print(mat2Eul(matMult(
            eul2Mat(newPnt['w'], newPnt['p'], newPnt['r']),
            eul2Mat(newEuler['w'], newEuler['p'], newEuler['r']))))
        print("Z-Height:")
        print(newPnt['z'] - newTbl['z'])
        print("X-Y Distance to Center of TT:")
        print(math.sqrt(((newPnt['x'] - newTbl['x'])**2) + \
                 ((newPnt['y'] - newTbl['y'])**2)))
        print("Offset Angle (Alpha):")
        print(newAlpha)
        print("New Coordinates:")
        print((newPnt['x'], newPnt['y'], newPnt['z']))
        print()
        return newPnt

if __name__ == "__main__":
        main()
