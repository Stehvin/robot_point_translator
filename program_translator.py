import get_points
import point_converter
import write_points

'''
INPUT INFORMATION
'''

inputPath = r"C:\Users\Example\test.ls"
outputPath = r"C:\Users\Example\OutputTest.ls"
oldBooth = 6
newBooth = 5
newProgUtoolNum = 6
newProgUframeNum = 0

'''
END OF INPUT INFORMATION
'''

# get points from input file
pntsDict = get_points.makePntsDict(inputPath)

# modify points from old booth to new booth
pntsDict = point_converter.execute(oldBooth, newBooth, newProgUtoolNum,
                                   newProgUframeNum, pntsDict)

# output file with new points
write_points.writeFile(inputPath, outputPath, pntsDict, newProgUtoolNum,
                       newProgUframeNum, newBooth)
