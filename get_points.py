import re

def makePntsDict(path):
        # set up regular expressions
        regex_pnt = re.compile('P\[(\d+)(\s?:\s?".*")?\]\s?\{\n')
        regex_UT = re.compile('UT\s:\s(\d{0,2}),')
        regex_UF = re.compile('UF\s:\s(\d{0,2}),')
        regex_x = re.compile('X\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_y = re.compile('Y\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_z = re.compile('Z\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_w = re.compile('W\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_p = re.compile('P\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_r = re.compile('R\s=\s+(-?\d{0,4}\.\d{2,3})')
        regex_e1 = re.compile('E1=\s+(-?\d{0,4}\.\d{2,3})')

        # initialize points dictionary
        pntsDict = {}

        # open input file (.ls format)
        with open(path, 'r') as f:
                lines = f.readlines()
                for lineIndex, line in enumerate(lines, 1):
                        
                        # if line is the beginning of a point, get the point
                        # number and put point values into a dictionary
                        pntMatch = regex_pnt.search(line)
                        if pntMatch:
                                
                                pntDict = {}

                                # get the point number
                                pntNum = pntMatch.group(1)

                                # get point values in dictionary
                                pntDict['x'] = float(regex_x.search(
                                        lines[lineIndex + 2]).group(1))
                                pntDict['y'] = float(regex_y.search(
                                        lines[lineIndex + 2]).group(1))		
                                pntDict['z'] = float(regex_z.search(
                                        lines[lineIndex + 2]).group(1))
                                pntDict['w'] = float(regex_w.search(
                                        lines[lineIndex + 3]).group(1))
                                pntDict['p'] = float(regex_p.search(
                                        lines[lineIndex + 3]).group(1))
                                pntDict['r'] = float(regex_r.search(
                                        lines[lineIndex + 3]).group(1))

                                if lines[lineIndex + 4] != "};\n":
                                        pntDict['e1'] = float(regex_e1.search(
                                                lines[lineIndex + 4]).group(1))

                                # get the uframe and utool numbers
                                pntDict['uf'] = int(regex_UF.search(
                                        lines[lineIndex + 1]).group(1))
                                pntDict['ut'] = int(regex_UT.search(
                                        lines[lineIndex + 1]).group(1))

                                # each point dictionary is an element of the
                                # points dictionary
                                pntsDict[int(pntNum)] = pntDict

        return pntsDict
