import os
import re
import boothSetupInfo

def formatDict(pntsDict):
        # format points dictionary for proper output to .ls files
        for pntNum, pnt in pntsDict.items():
                for coord in pnt:
                        pnt[coord] = round(pnt[coord], 3)
                        if pnt[coord] < 1 and pnt[coord] > 0:
                                pnt[coord] = format(pnt[coord], '.3f')[1:]
                        elif pnt[coord] < 0 and pnt[coord] > -1:
                                pnt[coord] = "-" + format(pnt[coord], '.3f')[2:]
                        else:
                                pnt[coord] = format(pnt[coord], '.3f')

def pointWrite(outfile, line, pntMatch, pntStatusCount, pntsDict,
               newProgUtoolNum, newProgUframeNum, newBooth, configStr):
        """Properly format and write each line after a point has been
           discovered.
        """     
        # set up uframe, utool, and configuration string regular expressions
        regex_UT = re.compile('UT\s:\s(\d{0,2}),')
        regex_UF = re.compile('UF\s:\s(\d{0,2}),')
        regex_ConStr = re.compile("CONFIG\s:\s'\w\s\w\s\w,")

        # first two lines can be printed verbatim
        if pntStatusCount == 0:
                outfile.write(line)
        elif pntStatusCount == 1:
                outfile.write(line)

        # third line - UT and UF must change for different utools/uframes,
        # configuration string must change to new booth's configuration string
        elif pntStatusCount == 2:
                newLine = regex_UF.sub("UF : {},".format(newProgUframeNum),
                                       line)
                newLine = regex_UT.sub("UT : {},".format(newProgUtoolNum),
                                       newLine)
                newLine = regex_ConStr.sub("CONFIG : '" + configStr + ",",
                                           newLine)
                outfile.write(newLine)
        
        # reformat fourth and fifth lines to paste new point coordinates
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
                if newBooth in (1, 2, 13, 16, 19):
                        outfile.write(",\n\t" + \
                                      "E1= {:>9}".format(
                                              pntsDict[pntNum]['e1']) + \
                                      " deg\n")
                else:
                        outfile.write("\n")

def writeFile(inputPath, outputPath, pntsDict, newProgUtoolNum,
              newProgUframeNum, newBooth):
        """Write old program to new program file path, using new points
           dictionary.
        """
        # get new booth's configuration string
        configStr = boothSetupInfo.boothInfo(newBooth)[5]
        
        # convert booth string into actual booth number (e.g. '2a' -> '2')
        regexNum = re.compile('[0-9]+')
        try:
                newBooth = int(regexNum.search(newBooth).group(0))
        except:
                pass
        
        # format point dictionary decimals
        formatDict(pntsDict)
        
        # set up regular expression to match point beginning
        regex_pnt = re.compile('P\[(\d+)(\s?:\s?".*")?\]\s?\{\n')

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
                                
                                # execute pointWrite function if line is within
                                # a point
                                if pntStatus:
                                        pointWrite(outfile, line, pntMatchCopy,
                                                   pntStatusCount, pntsDict,
                                                   newProgUtoolNum,
                                                   newProgUframeNum, newBooth,
                                                   configStr)
                                        pntStatusCount += 1
                                elif line[0:5] == r"/PROG":
                                        
                                        # change program name
                                        outfile.write("/PROG  " + \
                                                      newProgName + "\n")
                                else:
                                        outfile.write(line)

