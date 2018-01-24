import math
import re

import boothSetupInfo


def execute(oldBooth, newBooth, newProgUtoolNum, newProgUframeNum, oldPnts):
    """Execute the point conversion program."""
    # Get parameters from excel document.
    oldNozUtool, oldTbl, oldEuler, oldDustColLoc, oldDustColAng = \
        boothSetupInfo.boothInfo(oldBooth)[0:5]
    newNozUtool, newTbl, newEuler, newDustColLoc, newDustColAng = \
        boothSetupInfo.boothInfo(newBooth)[0:5]

    # Convert booth strings into actual booth numbers (e.g. '2a' -> 2).
    regexNum = re.compile('[0-9]+')
    try:
        oldBooth = int(regexNum.search(oldBooth).group(0))
    except:
        pass
    try:
        newBooth = int(regexNum.search(newBooth).group(0))
    except:
        pass

    # Get the user-specified Utool and Uframe setpoints for the new
    # program.
    newProgUtool = boothSetupInfo.curUtool(newBooth, newProgUtoolNum)
    newProgUframe = boothSetupInfo.curUframe(newBooth, newProgUframeNum)

    # Initialize new points dictionary and incremental motion
    # Utool/Uframe "zero-point queue" set.
    newPnts = {}
    incUfUt = set()

    # Iterate through points, translating each one.
    for pntNum, pnt in oldPnts.items():

        # If point is a raster (incremental motion) point, add its
        # Uframe and Utool to the zero-point queue.
        if pnt['inc']:
            incUfUt.add((pnt['uf'], pnt['ut']))

        # Translate point from old booth to new booth.
        newPnts[pntNum] = convert(pnt, oldNozUtool, oldTbl, oldEuler,
                                  oldDustColLoc, oldBooth, newProgUtool,
                                  newProgUframe, newNozUtool, newTbl, newEuler,
                                  newDustColLoc, newDustColAng)

    # [raster points] Translate zero-points in the zero-point queue.
    newZeroPnts = {}
    for ufut in incUfUt:
        newZeroPnts[ufut[0]] = {}
        newZeroPnts[ufut[0]][ufut[1]] = {
            'uf': ufut[0], 'ut': ufut[1], 'x': 0.0, 'y': 0.0, 'z': 0.0,
            'w': 0.0, 'p': 0.0, 'r': 0.0, 'e1': 0.0}
        newZeroPnts[ufut[0]][ufut[1]] = convert(
            newZeroPnts[ufut[0]][ufut[1]], oldNozUtool, oldTbl,
            oldEuler, oldDustColLoc, oldBooth, newProgUtool,
            newProgUframe, newNozUtool, newTbl, newEuler, newDustColLoc,
            newDustColAng)

    # [raster points] Convert raster (incremental motion) points.
    for pntNum, pnt in newPnts.items():
        if pnt['inc']:
            incTranslate(pnt, newZeroPnts)
    return dict(newPnts)


def eul2Mat(w, p, r):
    """Converts Euler angles to 3x3 rotation matrix."""
    # Initialize matrix.
    M = [[0, 0, 0], 
         [0, 0, 0], 
         [0, 0, 0]]
    
    # Convert Euler angles from degrees to radians.
    w = math.radians(w)
    p = math.radians(p)
    r = math.radians(r)

    # Calculate matrix elements.
    M[0][0] = math.cos(p) * math.cos(r)
    M[0][1] = ((math.sin(w)*math.sin(p)*math.cos(r))
               - (math.cos(w)*math.sin(r)))
    M[0][2] = ((math.sin(w)*math.sin(r))
               + (math.cos(w)*math.sin(p)*math.cos(r)))
    M[1][0] = math.cos(p) * math.sin(r)
    M[1][1] = ((math.cos(w)*math.cos(r))
               + (math.sin(w)*math.sin(p)*math.sin(r)))
    M[1][2] = ((math.cos(w)*math.sin(p)*math.sin(r))
               - (math.sin(w)*math.cos(r)))
    M[2][0] = -math.sin(p)
    M[2][1] = math.sin(w) * math.cos(p)
    M[2][2] = math.cos(w) * math.cos(p)
    return M


def mat2Eul(M):
    """Converts 3x3 rotation matrix to Euler angles."""
    r = math.degrees(math.atan2(M[1][0], M[0][0]))
    p = math.degrees(-math.asin(M[2][0]))
    w = math.degrees(math.atan2(M[2][1], M[2][2]))
    return w, p, r


def rotMatInv(M):
    """Inverts the 3x3 rotation matrix 'M'.
    The inverse of a rotation matrix is its transpose.
    """
    # Initialize output matrix.
    MI = [[0, 0, 0],
          [0, 0, 0],
          [0, 0, 0]]

    # Calculate matrix elements.
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
    """Turns a rotation matrix and a translation vector into a 4x4
    rotation and translation matrix.
    """
    # Create output matrix.
    M = [[rotM[0][0], rotM[0][1], rotM[0][2], translation['x']],
         [rotM[1][0], rotM[1][1], rotM[1][2], translation['y']],
         [rotM[2][0], rotM[2][1], rotM[2][2], translation['z']],
         [0, 0, 0, 1]]
    return M


def rot_matrix_4x4(rotM):
    """Changes a 3x3 rotation matrix to a 4x4 matrix."""
    # Create output matrix
    M = [[rotM[0][0], rotM[0][1], rotM[0][2], 0],
         [rotM[1][0], rotM[1][1], rotM[1][2], 0],
         [rotM[2][0], rotM[2][1], rotM[2][2], 0],
         [0, 0, 0, 1]]
    return M


def multMatVec(M, vector):
    """Multiplies a 4x4 matrix and a 4x1 vector."""
    # Initialize output vector.
    vector1 = {'x': 0, 'y': 0, 'z': 0}
    
    # Calculate vector elements.
    vector1['x'] = ((M[0][0]*vector['x']) + (M[0][1]*vector['y'])
                    + (M[0][2]*vector['z']) + M[0][3])
    vector1['y'] = ((M[1][0]*vector['x']) + (M[1][1]*vector['y'])
                    + (M[1][2]*vector['z']) + M[1][3])
    vector1['z'] = ((M[2][0]*vector['x']) + (M[2][1]*vector['y'])
                    + (M[2][2]*vector['z']) + M[2][3])
    return vector1


def matMult(A, B):
    """Multiplies matricies A and B (A*B = C).
    Full disclosure, I did not write this function. Taken from:
    https://stackoverflow.com/questions/10508021/matrix-multiplication-
    in-python
    """
    zip_b = zip(*B)
    zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(
        row_a, col_b)) for col_b in zip_b] for row_a in A]


def convertAngles(oldPnt, newEuler, zAngleDeg, center=False, opposite=False):
    """Converts Euler angles from a point in one booth to another booth.
    Old point Euler angles have already been changed to UGO angles.
    """
    # Convert old UGO angles to a rotation matrix.
    C = eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r'])

    # For opposite booths, find the z-angle (pointing straight down
    # x-axis is 0 degrees) for the old booth, flip it to the other side
    # (by reversing its sign), and recalculate the UGO angle rotation
    # matrix.
    if opposite:
        oldZAngle = getZAngle(C)
        zRotM = eul2Mat(0, 0, -2 * oldZAngle)
        C = matMult(zRotM, C)

    # For center booths, take the same steps described above for
    # opposite booths, except the new z-angle is based on the new dust
    # collector location (zAngleDeg).
    if center:
        oldZAngle = getZAngle(C)
        zRotM = eul2Mat(0, 0, zAngleDeg - oldZAngle)
        C = matMult(zRotM, C)
    
    # Find the new booth's faceplate-to-UGO rotation matrix.
    newB = eul2Mat(newEuler['w'], newEuler['p'], newEuler['r'])

    # Find the new booth's origin-to-faceplate Euler angles.
    newA = matMult(C, rotMatInv(newB))
    w, p, r = mat2Eul(newA)
    return w, p, r


def getZAngle(M):
    """Returns the z-angle of the projection of a 3D vector on the 2D
    X-Y plane. Rotation matrix "M" must be in UGO angles. Returns angle
    in degrees.
    """
    # Multiply unit vector [0, 0, 1] and rotation matrix to get X-Y
    # projection points.
    x = M[0][2]
    y = M[1][2]

    # Find tangent of X-Y point.
    return math.degrees(math.atan2(y, x))


def frame_to_world(uframePnt, uframe):
    """Converts Uframe point to a world frame point."""
    worldPnt = dict(uframePnt)
    
    # Convert xyz values to world frame.
    rotA = eul2Mat(uframe['w'], uframe['p'], uframe['r'])
    four = rotTransMat(rotA, uframe)
    vec = multMatVec(four, uframePnt)
    worldPnt['x'], worldPnt['y'], worldPnt['z'] = vec['x'], vec['y'], vec['z']

    # Convert wpr angles to world frame angles.
    rotB = eul2Mat(uframePnt['w'], uframePnt['p'], uframePnt['r'])
    worldPnt['w'], worldPnt['p'], worldPnt['r'] = mat2Eul(matMult(rotA, rotB))
    return worldPnt


def world_to_frame(worldPnt, uframe):
    """Converts a world point to a user frame point."""
    ufPnt = dict(worldPnt)
    
    # Set up translation matrix (Uframe setpoints).
    transM = [[1, 0, 0, -uframe['x']],
              [0, 1, 0, -uframe['y']],
              [0, 0, 1, -uframe['z']],
              [0, 0, 0, 1]]
        
    # Convert xyz values to Uframe.
    rotA3 = rotMatInv(eul2Mat(uframe['w'], uframe['p'], uframe['r']))
    rotA4 = rot_matrix_4x4(rotA3)
    M = matMult(rotA4, transM)
    vec = multMatVec(M, worldPnt)
    ufPnt['x'], ufPnt['y'], ufPnt['z'] = vec['x'], vec['y'], vec['z']

    # Convert wpr angles to Uframe.
    rotB = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
    ufPnt['w'], ufPnt['p'], ufPnt['r'] = mat2Eul(matMult(rotA3, rotB))
    return ufPnt


def world_to_utool(worldPnt, utool):
    """Converts a world frame point to a Utool point."""
    utoolPnt = dict(worldPnt)
    
    # Set up world-coord matrix and translate by Utool xyz setpoints.
    rotM = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
    M4x4 = rotTransMat(rotM, worldPnt)
    vec = multMatVec(M4x4, utool)
    utoolPnt['x'], utoolPnt['y'], utoolPnt['z'] = vec['x'], vec['y'], vec['z']

    # Change world-coord angles to Utool angles.
    rotA = eul2Mat(worldPnt['w'], worldPnt['p'], worldPnt['r'])
    rotB = eul2Mat(utool['w'], utool['p'], utool['r'])
    utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(matMult(rotA, rotB))
    return utoolPnt


def utool_to_world(utoolPnt, utool):
    """Converts a utool point to a world coordinate point."""
    utoolCopy = dict(utool)
    worldPnt = dict(utoolPnt)

    # Change Utool angles to world angles.
    rotA = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
    rotB = rotMatInv(eul2Mat(utoolCopy['w'], utoolCopy['p'], utoolCopy['r']))
    utoolPnt['w'], utoolPnt['p'], utoolPnt['r'] = mat2Eul(matMult(rotA, rotB))

    # Reverse Utool setpoint value signs.
    for key in utoolCopy:
        utoolCopy[key] *= -1
    
    # Set up Utool matrix.
    rotM = eul2Mat(utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
    M4x4 = rotTransMat(rotM, utoolPnt)
    vec = multMatVec(M4x4, utoolCopy)
    worldPnt['x'], worldPnt['y'], worldPnt['z'] = vec['x'], vec['y'], vec['z']
    worldPnt['w'], worldPnt['p'], worldPnt['r'] = (
        utoolPnt['w'], utoolPnt['p'], utoolPnt['r'])
    return worldPnt


def quadratic(a, b, c):
    """Find the solutions to a quadratic formula, where:
    0 = (a)x^2 + (b)x + (c)
    """
    term1 = (-1) * b
    term2 = math.sqrt((b**2) - (4*a*c))
    term3 = 2 * a

    solutions = []  
    signs = [1, -1]
    for sign in signs:
        try:
            x = (term1+(sign*term2)) / term3
        except:
            x = None
        solutions.append(x)
    return solutions


def UGO(pnt, Euler):
    """Convert a point's angles to Universal Gun Orientation (UGO)
    angles.
    """
    pnt1 = dict(pnt)
    rotM = matMult(eul2Mat(pnt1['w'], pnt1['p'], pnt1['r']),
                   eul2Mat(Euler['w'], Euler['p'], Euler['r']))
    pnt1['w'], pnt1['p'], pnt1['r'] = mat2Eul(rotM)
    return pnt1


def convert25ttcp(pnt, oldTbl):
    """Converts the turntable centerpoint (TTCP) x-value, when the old
    booth is 25. The new TTCP x-value is based upon the E1-value of
    the current point.
    """
    oldTbl1 = dict(oldTbl)
    oldTbl1['x'] = 1012.072 - (2*(-582.815 - pnt['e1']))
    return oldTbl1


def incTranslate(pnt, newZeroPnts):
    """Translates raster (incremental motion) points using the old
    booth's zero-point, which was previously translated to the new
    booth.
    """
    pnt['x'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['x']
    pnt['y'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['y']
    pnt['z'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['z']
    pnt['w'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['w']
    pnt['p'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['p']
    pnt['r'] -= newZeroPnts[pnt['olduf']][pnt['oldut']]['r']
    pnt['e1'] = 0.00
    return None


def convert(pnt, oldNozUtool, oldTbl, oldEuler, oldDustColLoc, oldBooth,
            newProgUtool, newProgUframe, newNozUtool, newTbl, newEuler,
            newDustColLoc, newDustColAng):
    """Converts old booth point to UGO, nozzle utool, world frame point.
    Then, runs the spfTranslate function which translates the point to a
    new booth. After that, this function changes the point to a Utool,
    Uframe point.
    """
    # Get the original program's Utool and Uframe info.
    oldProgUtool = boothSetupInfo.curUtool(oldBooth, pnt['ut'])
    oldProgUframe = boothSetupInfo.curUframe(oldBooth, pnt['uf'])
    
    # Convert point from Uframe to world frame.
    pnt = frame_to_world(pnt, oldProgUframe)

    # Convert point from original Utool to nozzle Utool and convert
    # angles to UGO (Universal Gun Orientation) angles.
    pnt = utool_to_world(pnt, oldProgUtool)
    pnt = world_to_utool(pnt, oldNozUtool)
    pnt = UGO(pnt, oldEuler)

    # Convert z-angle ("roll" angle).
    if oldDustColLoc == newDustColLoc:
        center = False
        opposite = False
        zAngleDeg = getZAngle(eul2Mat(pnt['w'], pnt['p'], pnt['r']))
    elif ((oldDustColLoc == 'R' and newDustColLoc == 'L') or
          (oldDustColLoc == 'L' and newDustColLoc == 'R')):
        center = False
        opposite = True
        zAngleDeg = -getZAngle(eul2Mat(pnt['w'], pnt['p'], pnt['r']))
    elif newDustColLoc == 'C' or oldDustColLoc == 'C':
        center = True
        opposite = False
        zAngleDeg = newDustColAng

    # Convert turntable centerpoint x-value for booth 25.
    if oldBooth == 25:
        oldTbl = convert25ttcp(pnt, oldTbl)

    # Find new point.
    newPnt = spfTranslate(pnt, oldTbl, newTbl, oldEuler, newEuler, zAngleDeg,
                          center, opposite)
        
    # Convert new point from nozzle Utool to user-specified Utool.
    newPnt = utool_to_world(newPnt, newNozUtool)
    newPnt = world_to_utool(newPnt, newProgUtool)

    # Convert new point from world frame to user-specified frame.
    newPnt = world_to_frame(newPnt, newProgUframe)
    return dict(newPnt)


def spfTranslate(oldPnt, oldTbl, newTbl, oldEuler, newEuler, zAngleDeg,
                 center=False, opposite=False):
    """Translate UGO, nozzle utool, world frame point from old booth to
    new booth.
    """
    # Get original booth center offset angle (alpha).
    oldFi = math.atan2(oldTbl['y'] - oldPnt['y'], oldTbl['x'] - oldPnt['x'])
    oldAlpha = math.radians(
        getZAngle(eul2Mat(oldPnt['w'], oldPnt['p'], oldPnt['r']))) - oldFi
    
    # Change old z-angle to point on other side of the table if opposite
    # points are necessary (reverse the sign of alpha).
    if opposite:
        oldAlpha = -oldAlpha

    # Initialize new point.
    newPnt = dict(oldPnt)
    newPnt['olduf'], newPnt['oldut'] = newPnt['uf'], newPnt['ut']
    del newPnt['uf'], newPnt['ut']

    # Convert old UGO (Universal Gun Orientation) angles to new UGO
    # angles, then to new world angles.
    newPnt['w'], newPnt['p'], newPnt['r'] = convertAngles(
        oldPnt, newEuler, zAngleDeg, center, opposite)

    # Get z-height and set E1 value.
    newPnt['z'] = oldPnt['z'] - oldTbl['z'] + newTbl['z']
    try:
        newPnt['e1'] = oldPnt['e1']
    except:
        newPnt['e1'] = 91.588

    # Set up x and y system of equations.
    angVar = math.tan(math.radians(zAngleDeg) - oldAlpha)
    distSq = ((oldPnt['x']-oldTbl['x'])**2) + ((oldPnt['y']-oldTbl['y'])**2)
    
    # Execute quadratic equation solver.
    a = 1 + (angVar**2)
    b = ((-2)*(newTbl['x'])) - (2*(angVar**2)*newTbl['x'])
    c = (newTbl['x']**2) + ((angVar**2)*(newTbl['x']**2)) - distSq
    xPnts = quadratic(a, b, c)
    
    for x in xPnts:
        if x:
            # Solve for y.
            y = newTbl['y'] - (angVar*newTbl['x']) + (angVar*x)
            
            # Compare new and old booth center offset angles.
            newAlpha = (math.radians(zAngleDeg)
                        - math.atan2(newTbl['y'] - y, newTbl['x'] - x))
            if (math.isclose(newAlpha, oldAlpha, rel_tol=1e-06) or
                math.isclose(newAlpha - (2*math.pi), oldAlpha,
                             rel_tol=1e-06) or
                math.isclose(newAlpha + (2*math.pi), oldAlpha,
                             rel_tol=1e-06)):
                newPnt['x'] = x
                newPnt['y'] = y
                break
    return newPnt
