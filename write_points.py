import os
import re

def formatDict(pntsDict):
        # format points dictionary for proper output to .ls files
        for pntNum, pnt in pntsDict.items():
                for coord in pnt:
                        if pnt[coord] < 1 and pnt[coord] > 0:
                                pnt[coord] = format(pnt[coord], '.3f')[1:]
                        elif pnt[coord] < 0 and pnt[coord] > -1:
                                pnt[coord] = "-" + format(pnt[coord], '.3f')[2:]
                        else:
                                pnt[coord] = format(pnt[coord], '.3f')

def pointWrite(outfile, line, pntMatch, pntStatusCount, pntsDict,
               newCurUtoolNum, newBooth):
        """Properly format and write each line after a point has been
            discovered.
        """
        # set up uframe and utool regular expressions
        regex_UT = re.compile('UT\s:\s(\d{0,2}),')
        regex_UF = re.compile('UF\s:\s(\d{0,2}),')

        # first two lines can be printed verbatim
        if pntStatusCount == 0:
                outfile.write(line)
        elif pntStatusCount == 1:
                outfile.write(line)

        # third line - UT must change for different utools, UF = 0
        # e.g. utool number = 1 -> UT : 1
        elif pntStatusCount == 2:
                newLine = regex_UF.sub("UF : 0,", line)
                newLine = regex_UT.sub("UT : {},".format(newCurUtoolNum),
                                       newLine)
                outfile.write(newLine)
        
        # fourth and fifth lines must be reformated to paste new point
        # coordinates
        elif pntStatusCount == 3:
                pntNum = int(pntMatch.group(1))
                outfile.write("\t" + \
                              "X = {:>9}".format(pntsDict[pntNum]['x']) + \
                              "  mm,\t" + \
                              "Y = {:>9}".format(pntsDict[pntNum]['y']) + \
                              "  mm,\t" + \
                              "Z = {:>9}".format(pntsDict[pntNum]['z']) + \
                              "  mm,\n")
        elif pntStatusCount == 4:
                pntNum = int(pntMatch.group(1))
                outfile.write("\t" + \
                              "W = {:>9}".format(pntsDict[pntNum]['w']) + \
                              " deg,\t" + \
                              "P = {:>9}".format(pntsDict[pntNum]['p']) + \
                              " deg,\t" + \
                              "R = {:>9}".format(pntsDict[pntNum]['r']) + \
                              " deg")
                if newBooth in (1, 2):
                        outfile.write(",\n\t" + \
                                "E1 = {:>9}".format(pntsDict[pntNum]['e1']) + \
                                "deg\n")
                else:
                        outfile.write("\n")

def writeFile(inputPath, outputPath, pntsDict, newCurUtoolNum, newBooth):
        """Write old program to new program file path, using new points
           dictionary.
        """
        # format point dictionary decimals
        formatDict(pntsDict)
        
        # set up regular expression to match point beginning
        regex_pnt = re.compile('P\[(\d+)(:".*")?\]\{\n')

        # initialize beginning of point checker to false
        pntStatus = False

        # get output file name (without .ls)
        newProgName = os.path.basename(outputPath)[:-3]
        
        # copy input file to output file
        with open(inputPath, 'r') as infile:
                with open(outputPath, 'w') as outfile:
                        for line in infile:
                                
                                # search line for beginning of point
                                pntMatch = regex_pnt.search(line)
                                
                                # beginning of point found
                                if pntMatch:
                                        pntStatus = True
                                        pntStatusCount = 0
                                        pntMatchCopy = pntMatch
                                
                                # end of point found
                                elif line == "};\n":
                                        pntStatus = False
                                
                                # execute pointWrite function if line is
                                # within a point
                                if pntStatus:
                                        pointWrite(outfile, line, pntMatchCopy,
                                                   pntStatusCount, pntsDict,
                                                   newCurUtoolNum, newBooth)
                                        pntStatusCount += 1
                                elif line[0:5] == r"/PROG":
                                        
                                        # change program name
                                        outfile.write("/PROG  " + \
                                                      newProgName + "\n")
                                else:
                                        outfile.write(line)

