import get_points
import point_converter
import write_points


'''
INPUT INFORMATION
'''

inputPath = r"C:\Users\Example\inTestProg.ls"
outputPath = r"C:\Users\Example\outTestProg.ls"
oldBooth = 6
newBooth = 5
newProgUtoolNum = 6
newProgUframeNum = 0

'''
END OF INPUT INFORMATION
'''

# Get points from input file.
pntsDict = get_points.makePntsDict(inputPath)

# Convert points from old booth to new booth.
pntsDict = point_converter.execute(oldBooth, newBooth, newProgUtoolNum,
                                   newProgUframeNum, pntsDict)

# Write output file with new points.
write_points.writeFile(inputPath, outputPath, pntsDict, newProgUtoolNum,
                       newProgUframeNum, newBooth)
