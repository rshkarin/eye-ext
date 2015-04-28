from ij import IJ, ImagePlus, ImageStack, Prefs, WindowManager
from ij.io import TiffDecoder, FileSaver
from ij.plugin.frame import RoiManager
from ij.process import ImageStatistics, ByteProcessor, FloatProcessor, ShortProcessor, ImageConverter, AutoThresholder
from ij.plugin.filter import ParticleAnalyzer, Binary
import ij.plugin.filter.PlugInFilter
from ij.plugin import Slicer
from ij.measure import ResultsTable, Measurements
from ij.gui import ImageWindow
from java.lang import Double
from java.lang import Character
from fiji.threshold import Auto_Threshold
import os, sys, datetime, math, csv, codecs, cStringIO

inputPath = "/Users/rshkarin/Documents/Eye centroids detection/Data"
outputPath = "/Users/rshkarin/Documents/Eye centroids detection/Eyes"

fishPrefix = "fish"
fileExt = ".tif"
fishNumbers = ["202"]
cacheFolderName = "Temp"

Prefs.blackBackground = True

def printLog(title, message):
	print "%s - %s (%s)" % (datetime.datetime.now(), title, message)

'''
Reslice from specific direction
'''
def performReslice(imp, direction):
	newStack = imp.duplicate()

	IJ.run(newStack, "Reslice [/]...", "output=1.000 start=" + direction + " avoid")

	return IJ.getImage()

'''
Filter random noise in the beginning and in the end

Input: ImagePlus, num of slices from the beginning and the end in percentage
Output: Filtered 8-bit ImagePlus
'''
def filterNoisedSices(imp, percentOfBoundSlices):
	stackSize = imp.getStackSize()
	numSlicesCheck = stackSize / 100. * percentOfBoundSlices

	firstBound, lastBound = range(numSlicesCheck), range(stackSize - numSlicesCheck, stackSize)

	filteredStak = ImageStack(imp.getWidth(), imp.getHeight())
	testingImage = imp.duplicate()

	for i in range(testingImage.getStackSize()):
		if i in firstBound or lastBound:
			ipt = testingImage.getImageStack().getProcessor(i + 1)

			hist = ipt.getHistogram()
			lowTH = Auto_Threshold.Otsu(hist)
			ipt.threshold(lowTH)

			binner = Binary()
			binner.setup("dilate", None)

			for _ in range(7):
				binner.run(ipt)

			pixels = ipt.convertToFloat().getPixels()

			totalCount = len(pixels)
			#print totalCount
			nonNull = 0
			for p in pixels:
				if p > 0:
					nonNull = nonNull + 1
			#print nonNull
			#arr = filter(lambda x: x > 0, list(pixels))
			#print arr
			#nonNull = len(arr)

			if float(nonNull)/totalCount < 0.5:
				filteredStak.addSlice(imp.getImageStack().getProcessor(i + 1))
		else:
			filteredStak.addSlice(imp.getImageStack().getProcessor(i + 1))

	outImage = ImagePlus("filtered_stack_" + imp.getTitle(), filteredStak)
	ImageConverter(outImage).convertToGray8()

	return outImage

'''
Determine stack Z boundaries
'''
def getStackBoundaries(imp):
	fisrtSlice, lastSlice = 0, 0

	for i in range(imp.getStackSize()):
		if sum(imp.getImageStack().getProcessor(i + 1).convertToShort(False).getPixels()) > 0:
			fisrtSlice = i + 1
			break

	for i in reversed(range(imp.getStackSize())):
		if sum(imp.getImageStack().getProcessor(i + 1).convertToShort(False).getPixels()) > 0:
			lastSlice = i + 1
			break
			
	return fisrtSlice, lastSlice

'''
Inverting binary stack
'''
def prepareInvertedBinaryStack(imp, lZBound, hZBound):
	invertedStak = ImageStack(imp.getWidth(), imp.getHeight())
	bounds = range(lZBound, hZBound + 1)

	for i in range(imp.getStackSize()):
		if (i + 1) in bounds:
			hist = imp.getImageStack().getProcessor(i + 1).getHistogram()
			lowTH = Auto_Threshold.Otsu(hist)

			imp.getImageStack().getProcessor(i + 1).threshold(lowTH) 
			imp.getImageStack().getProcessor(i + 1).invert()

			invertedStak.addSlice(imp.getImageStack().getProcessor(i + 1))

	outImage = ImagePlus("inverted_binary_" + imp.getTitle(), invertedStak)

	return outImage

def test():
	imp = IJ.openImage("/Users/rshkarin/Documents/Fishes/fish202_filtered/filtered_fish202_8bit_416x396x1857.tif")
	#impFront = IJ.openImage("/Users/rshkarin/Documents/Eyes/inverted_binary_fish202_8bit_416x396x1857.tif")
	#impRight = IJ.openImage("/Users/rshkarin/Documents/Eyes/inverted_binary_right_fish202_8bit_416x396x1857.tif")
	#impTop = IJ.openImage("/Users/rshkarin/Documents/Eyes/inverted_binary_top_fish202_8bit_416x396x1857.tif")

	#statisticsDict, headings = getCirlcesStatistics(impFront, "/Users/rshkarin/Documents/Eyes/fish202", "202", "front", 50, 0.85)
	#cirRight = getCirlcesStatistics(impRight, 50, 0.85)
	#cirTop = getCirlcesStatistics(impTop, 50, 0.85)
	#print headings
	#x,y = getStackCirlceCurve(impFront, statisticsDict, headings, 0)
	#filteredStack = filterNoisedSices(imp, 5)
	x, y, width, height = getEyesBoundingArea(imp, 1681)
	LEyeXY, REyeXY = getLeftRightEyesCentroids(x, y, width, height)
	print LEyeXY, REyeXY

'''
Write statistic
'''
def writeStatistics(outputPath, fileName, fileExt, statDict, headers):
	if not os.path.exists(outputPath):
		os.makedirs(outputPath)

	outFile = open(os.path.join(outputPath, fileName + fileExt), 'w')
	csvWriter = csv.writer(outFile, delimiter=';')
	csvWriter.writerow(headers)

	for row in range(len(statDict[headers[0]])):
		outRow = [statDict[header][row] for header in headers]
		csvWriter.writerow(outRow)

	outFile.close()

'''
Return specific statistic of the stack

Kind of statistic:
1 - Area
2 - The product of centroids (X*Y)
'''
def getStackCirlceCurve(imp, statistics, headings, kind=0):
	y = []
	
	for i in x:
		if i in statistics['Slice']:
			if kind == 0:
				y.append(statistics['Area'])
			elif kind == 1:
				y.append(statistics['X'] * statistics['Y'])
			else:
				y.append(statistics['Area'])

	return y

'''
Get the circle statistic of the stack along the specific direction

Returns: statistic, headings
'''
def getCirlcesStatistics(imp, outputPath, fishNumber, direction, min_area=0, min_circularity=0.0):
	cacheFilename = '_'.join(["fish" + str(fishNumber), str(direction), str(min_circularity).replace('.','p'), str(min_area)])

	outputPath = os.path.join(outputPath, cacheFolderName)

	statFileExt = '.csv'
	headings = []
	statisticsDict = {}

	if not os.path.isfile(os.path.join(outputPath, cacheFilename) + statFileExt):
		paOptions = ParticleAnalyzer.SHOW_PROGRESS + ParticleAnalyzer.CLEAR_WORKSHEET
		paMeasurements = Measurements.AREA + Measurements.CENTROID + Measurements.CIRCULARITY + Measurements.STACK_POSITION + Measurements.SLICE + Measurements.RECT + Measurements.SHAPE_DESCRIPTORS + Measurements.STACK_POSITION
		rt = ResultsTable()

		pa = ParticleAnalyzer(paOptions, paMeasurements, rt, min_area, Double.POSITIVE_INFINITY, min_circularity, 1.0)
		pa.setHideOutputImage(True)

		statisticsDict = {}

		for i in range(imp.getStackSize()):
			imp.setSliceWithoutUpdate(i + 1)
			pa.analyze(imp)

		#headings = [u'Index']
		#headings.extend(list(rt.getHeadings()))
		headings = list(rt.getHeadings())
		for header in headings:
			statisticsDict[header] = list(rt.getColumn(rt.getColumnIndex(header)))

		rt.reset()

		writeStatistics(outputPath, cacheFilename, statFileExt, statisticsDict, headings)
	else:
		csvfile = open(os.path.join(outputPath, cacheFilename + statFileExt), 'rb')
		csvReader = csv.reader(csvfile, delimiter=';')
		headings = csvReader.next()
		for header in headings:
			statisticsDict[header] = []

		for row in csvReader:
			for i in range(len(headings)):
				if row[i].isdigit():
					statisticsDict[headings[i]].append(int(row[i]))
				else:
					statisticsDict[headings[i]].append(float(row[i]))

	return statisticsDict, headings

'''
Get bounding box of eyes region at specific slice
'''
def getEyesBoundingArea(imp, sliceIdx, min_area=600, depth=60):
	paOptions = ParticleAnalyzer.SHOW_PROGRESS + ParticleAnalyzer.CLEAR_WORKSHEET
	paMeasurements = Measurements.AREA + Measurements.RECT
	rt = ResultsTable()

	sliceIdx = int(sliceIdx)

	pa = ParticleAnalyzer(paOptions, paMeasurements, rt, min_area, Double.POSITIVE_INFINITY, 0.0, 1.0)
	pa.setHideOutputImage(True)

	hist = imp.getImageStack().getProcessor(sliceIdx).getHistogram()
	lowTH = Auto_Threshold.Otsu(hist)
	imp.getImageStack().getProcessor(sliceIdx).threshold(lowTH)

	imp.setSliceWithoutUpdate(sliceIdx)
	binner = Binary()
	binner.setup("fill", None)
	binner.run(imp.getImageStack().getProcessor(sliceIdx))

	pa.analyze(imp)

	statisticsDict = {}
	headings = list(rt.getHeadings())
	for header in headings:
		statisticsDict[header] = list(rt.getColumn(rt.getColumnIndex(header)))

	rt.reset()

	maxAreaIdx = statisticsDict['Area'].index(max(statisticsDict['Area']))
	x, y = statisticsDict['BX'][maxAreaIdx], statisticsDict['BY'][maxAreaIdx]
	width, height = statisticsDict['Width'][maxAreaIdx], statisticsDict['Height'][maxAreaIdx]

	z = sliceIdx - int(depth/2)
	z = 0 if z < 0 else z

	return (x, y, z, width, height, depth)

'''
Get left and right eyes' bounding box centroids
'''
def getLeftRightEyesCentroids(x, y, width, height):
	width2 = width / 2
	LEyeX = width2 / 2
	REyeX = width2 / 2

	height2 = height /2
	LEyeY = height2
	REyeY = height2

	return (x + LEyeX, y + LEyeY), (x + width2 + REyeX, y + REyeY)

'''
Save as tiff stack
'''
def saveTiff(imp, outputPath, currentFileName):
	fileSaver = FileSaver(imp)
	currentFileName = 'inverted_binary_' + currentFileName
	fileSaver.saveAsTiffStack(os.path.join(outputPath, currentFileName))

'''
0 - Frontal (x, y, z)
1 - Top (x, z, y) with respect to frontal
2 - Right (y, z, x) with respect to frontal
'''
def getEyeCoordinate(orientFrontCoord=None, orientTopCoord=None, orientRightCoord=None):
	nCoords = 0
	nCoords = nCoords + (1 if orientFrontCoord else 0)
	nCoords = nCoords + (1 if orientTopCoord else 0)
	nCoords = nCoords + (1 if orientRightCoord else 0)

	xVal, yVal, zVal = 0., 0., 0.

	xVal = xVal + (orientFrontCoord[0] if orientFrontCoord else 0)
	xVal = xVal + (orientTopCoord[0] if orientTopCoord else 0)
	xVal = xVal + (orientRightCoord[1] if orientRightCoord else 0)

	yVal = yVal + (orientFrontCoord[1] if orientFrontCoord else 0)
	yVal = yVal + (orientTopCoord[2] if orientTopCoord else 0)
	yVal = yVal + (orientRightCoord[0] if orientRightCoord else 0)

	zVal = zVal + (orientFrontCoord[2] if orientFrontCoord else 0)
	zVal = zVal + (orientTopCoord[1] if orientTopCoord else 0)
	zVal = zVal + (orientRightCoord[1] if orientRightCoord else 0)

	#print xVal, yVal, zVal

	return (xVal / float(nCoords), yVal / float(nCoords), zVal / float(nCoords))

'''
0 - Frontal (x, y, z)
1 - Top (x, z, y) with respect to frontal
2 - Right (y, z, x) with respect to frontal
3 - Skip slice
CC - centroid coordinate

Return values
 1 - Left eye
 2 - Right eye
-1 - Unknown
'''
def eyesSide(coord, BBox=None, sliceIdx=0, orientation=0, restrictedByDepth=False):
	LEyeCC = ()
	REyeCC = ()

	if not BBox:
		return -1

	if restrictedByDepth:
		if orientation == 0:
			if sliceIdx not in range(int(BBox[2]), int(BBox[2] + BBox[5] + 1)):
				print "Frontal Value %f not in range from %d to %d" % (sliceIdx, int(BBox[2]), int(BBox[2] + BBox[5] + 1))
				return 3
		elif orientation == 1:
			if sliceIdx not in range(int(BBox[1]), int(BBox[1] + BBox[4] + 1)):
				print "Top Value %f not in range from %d to %d" % (sliceIdx, int(BBox[1]), int(BBox[1] + BBox[4] + 1))
				return 3
		elif orientation == 2:
			if sliceIdx not in range(int(BBox[0]), int(BBox[0] + BBox[3] + 1)):
				print "Right Value %f not in range from %d to %d" % (sliceIdx, int(BBox[0]), int(BBox[0] + BBox[3] + 1))
				return 3
		else:
			print "Unknown error of input coordinates and bounding box."

	if orientation == 0:
		LEyeCC, REyeCC = getLeftRightEyesCentroids(BBox[0], BBox[1], BBox[3], BBox[4])
	elif orientation == 1:
		LEyeCC, REyeCC = getLeftRightEyesCentroids(BBox[0], BBox[2], BBox[3], BBox[5])
	elif orientation == 2:
		LEyeCC, REyeCC = getLeftRightEyesCentroids(BBox[1], BBox[2], BBox[4], BBox[5])

	#print LEyeCC, REyeCC

	LDist = math.sqrt((coord[0] - LEyeCC[0]) ** 2 + (coord[1] - LEyeCC[1]) ** 2)
	RDist = math.sqrt((coord[0] - REyeCC[0]) ** 2 + (coord[1] - REyeCC[1]) ** 2)

	if LDist < RDist:
		return 1
	elif RDist < LDist:
		return 2
	else:
		return -1

def average(s):
	return sum(s) * 1.0 / len(s)

def variance(x):
	return map(lambda y: (y - average(x)) ** 2, x)

def stdev(x):
	return math.sqrt(average(variance(x)))

def calcStd(i, array, num_neighbours):
	L2 = (num_neighbours - 1) / 2

	values = array[i - L2:i + L2 + 1]

	if not sum(values):
		return 0

	num_neighbours_recalc = num_neighbours if len(values) == num_neighbours else len(values)

	return 1./stdev(values)

def gauss(n=11,sigma=1):
	r = range(-int(n/2),int(n/2)+1)
	return [1 / (sigma * math.sqrt(2*math.pi)) * math.exp(-float(x)**2/(2*sigma**2)) for x in r]

def calcGauss(i, array, num_neighbours, gaussVals):
	L2 = (num_neighbours - 1) / 2

	values = array[i - L2:i + L2 + 1]

	if not sum(values):
		return 0

	new_vals = values

	for i in range(len(values)):
		new_vals[i] = values[i] * gaussVals[i]

	return sum(new_vals)

'''
Obtain preliminary BBox position on frontal Z direction
'''
def getPreliminaryBBoxZPosition(dotXY, areas, sliceIdxs, maxPrelimBBoxZPosResidual=100):
	gaussVals = gauss(9, 9)

	smoothAreaValues = [calcGauss(i, areas, 9, gaussVals) for i in range(len(areas))]
	smoothXYValues = [calcStd(i, dotXY, 5) for i in range(len(dotXY))]

	maxSmoothXY = max(smoothXYValues)
	maxSmoothArea = max(smoothAreaValues)

	residualsXY = []
	residualsAreas = []

	for xyVal in dotXY:
		residualsXY.append(abs(xyVal - maxSmoothXY))

	for areaVal in areas:
		residualsAreas.append(abs(areaVal - maxSmoothArea))

	idxMinResidualSmoothXY = residualsXY.index(min(residualsXY))
	idxMinResidualSmoothArea = residualsAreas.index(min(residualsAreas))

	prelimBBoxZPosXY = sliceIdxs[idxMinResidualSmoothXY]
	prelimBBoxZPosArea = sliceIdxs[idxMinResidualSmoothArea]

	if abs(prelimBBoxZPosXY - prelimBBoxZPosArea) > maxPrelimBBoxZPosResidual:
		print "Warning: very big difference between XY and area Z positions, the biggest will selected."

	return prelimBBoxZPosXY

'''
Get eye coordinates by direction
'''
def getEyesCoordinatesByOrientation(centerX, centerY, indices, sliceIdxs, orientation, BBox, restrictedByDepth, neighbourSlices):
	leftEye = []
	rightEye = []

	print centerX
	
	for x, y, idx, sliceIdx in zip(centerX, centerY, indices, sliceIdxs):
		if eyesSide((x, y), BBox, sliceIdx, orientation, restrictedByDepth) == 1:
			leftEye.append(idx)
		elif eyesSide((x, y), BBox, sliceIdx, orientation, restrictedByDepth) == 2:
			rightEye.append(idx)
		elif eyesSide((x, y), BBox, sliceIdx, orientation, restrictedByDepth) == 3:
			print "Slice %d: skipped." % sliceIdx
		else:
			print "Slice %d: something is wrong." % sliceIdx

	lEyeValuesArea = [areas[lIdx] for lIdx in leftEye]
	rEyeValuesArea = [areas[rIdx] for rIdx in rightEye]

	#print "LEyeMaxAreaSlice = %d, REyeMaxAreaSlice = %d" % (sliceIdxs[leftEye[lEyeValuesArea.index(max(lEyeValuesArea))]], sliceIdxs[rightEye[rEyeValuesArea.index(max(rEyeValuesArea))]])
	
	lEyeValuesXY = [dotXY[lIdx] for lIdx in leftEye]
	rEyeValuesXY = [dotXY[rIdx] for rIdx in rightEye]

	lEyeValuesXYFilt = [calcStd(i, lEyeValuesXY, 5) for i in range(len(lEyeValuesXY))]
	rEyeValuesXYFilt = [calcStd(i, rEyeValuesXY, 5) for i in range(len(rEyeValuesXY))]

	#First pass - detect preliminary eyes position
	lEyeXYMaxVal = max(lEyeValuesXYFilt)
	rEyeXYMaxVal = max(rEyeValuesXYFilt)

	lEyeValuesXYFiltMaxIndex = lEyeValuesXYFilt.index(lEyeXYMaxVal)
	rEyeValuesXYFiltMaxIndex = rEyeValuesXYFilt.index(rEyeXYMaxVal)

	lEyePrelimZPos = sliceIdxs[leftEye[lEyeValuesXYFilt.index(lEyeXYMaxVal)]]
	rEyePrelimZPos = sliceIdxs[rightEye[rEyeValuesXYFilt.index(rEyeXYMaxVal)]]

	lEyeAllowedRange = leftEye[((lEyeValuesXYFiltMaxIndex - neighbourSlices) if (lEyeValuesXYFiltMaxIndex - neighbourSlices) >= 0 else 0) : 
								(lEyeValuesXYFiltMaxIndex + neighbourSlices) if (lEyeValuesXYFiltMaxIndex + neighbourSlices) < len(leftEye) else len(leftEye) - 1]

	rEyeAllowedRange = rightEye[((rEyeValuesXYFiltMaxIndex - neighbourSlices) if (rEyeValuesXYFiltMaxIndex - neighbourSlices) >= 0 else 0) : 
								(rEyeValuesXYFiltMaxIndex + neighbourSlices) if (rEyeValuesXYFiltMaxIndex + neighbourSlices) < len(rightEye) else len(rightEye) - 1]

	lEyeAreaAllowed = [areas[i] for i in lEyeAllowedRange]
	rEyeAreaAllowed = [areas[i] for i in rEyeAllowedRange]

	#Seond pass - detect end eyes position
	lEyeEndXPos = centerX[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndXPos = centerX[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	lEyeEndYPos = centerY[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndYPos = centerY[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	lEyeEndZPos = sliceIdxs[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndZPos = sliceIdxs[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	#print lEyePrelimZPos, rEyePrelimZPos
	#print lEyeEndZPos, rEyeEndZPos

	print "LEye coord " + str((lEyeEndXPos, lEyeEndYPos, lEyeEndZPos))
	print "REye coord " + str((rEyeEndXPos, rEyeEndYPos, rEyeEndZPos))

	return (lEyeEndXPos, lEyeEndYPos, lEyeEndZPos), (rEyeEndXPos, rEyeEndYPos, rEyeEndZPos)


'''
Obtain eyes' coordinates from specific direction
'''
def getEyesCoordinates(imp, statisticsDict, headings, orientation=0, BBox=None, restrictedByDepth=False, neighbourSlices=60, maxPrelimBBoxZPosResidual=100):
	centerX = statisticsDict['X']
	centerY = statisticsDict['Y']
	sliceIdxs = statisticsDict['Slice']
	areas = statisticsDict['Area']
	dotXY = map(lambda x, y: x * y, statisticsDict['X'], statisticsDict['Y'])
	indices = range(len(sliceIdxs))

	if not BBox:
		prelimBBoxSliceIdx = getPreliminaryBBoxZPosition(dotXY, areas, sliceIdxs, maxPrelimBBoxZPosResidual)
		BBox = getEyesBoundingArea(imp, prelimBBoxSliceIdx)

	LEyeCoord, REyeCoord = getEyesCoordinatesByOrientation(centerX, centerY, indices, sliceIdxs, orientation, BBox, restrictedByDepth, neighbourSlices)

	return LEyeCoord, REyeCoord, BBox

'''
Perform filtering, binarizing, inverting and reslicing of the stack from Front, Top and Right 
'''
def extractEyes(inputPath, outputPath, fishNumbers, fishPrefix, fileExt):
	for fishNumber in fishNumbers:
		currentPath = os.path.join(inputPath, fishPrefix + fishNumber)
		files = [f for f in os.listdir(currentPath) if os.path.isfile(os.path.join(currentPath,f)) and f.startswith(fishPrefix + fishNumber) and f.endswith(fileExt)]

		currentFileName = ''
		if files:
			currentFileName = files[0]

		pathToVolume = os.path.join(currentPath, currentFileName)

		printLog("Opening", pathToVolume)
		imp = IJ.openImage(pathToVolume)

		if imp is None:  
			print "Could not open image from file:", currentFileName  
			continue

		td = TiffDecoder(currentPath, currentFileName)

		stackSize = imp.getStackSize()
		stackDims = imp.getDimensions()

		stackInfo = td.getTiffInfo()

		printLog("Obtaining Z bounds", pathToVolume)
		lowZbound, highZbound = getStackBoundaries(imp)

		printLog("Filtring noisy slices", pathToVolume)
		filteredStack = filterNoisedSices(imp, 5)
		filteredTopReslice = performReslice(filteredStack, "Top")

		printLog("Binarizing and inverting", pathToVolume)
		invertedStack = prepareInvertedBinaryStack(filteredStack, lowZbound, highZbound)
		filteredStack.close()

		invertedTopReslice = prepareInvertedBinaryStack(filteredTopReslice, 0, filteredTopReslice.getStackSize())
		filteredTopReslice.close()

		printLog("Saving as tiff stack", pathToVolume)
		saveTiff(invertedStack, outputPath, currentFileName)
		saveTiff(invertedTopReslice, outputPath, 'top_' + currentFileName)

		printLog("Gathering the area and circularity statistic", pathToVolume)
		statisticsDictFrontal, headingsFrontal = getCirlcesStatistics(invertedStack, currentPath, fishNumber, "front", 50, 0.85)
		statisticsDictTop, headingsTop = getCirlcesStatistics(invertedTopReslice, currentPath, fishNumber, "top", 50, 0.85)

		printLog("Get eyes coordinates", pathToVolume)
		LEyeFrontalCoord, REyeFrontalCoord, BBox = getEyesCoordinates(invertedStack, statisticsDictFrontal, headingsFrontal)
		LEyeTopCoord, REyeTopCoord, = getEyesCoordinates(invertedTopReslice, statisticsDictTop, headingsTop, 1, BBox, True)

		printLog("Interpolate coordinates", pathToVolume)
		LEyeCoord = getEyeCoordinate(LEyeFrontalCoord, LEyeTopCoord)
		REyeCoord = getEyeCoordinate(REyeFrontalCoord, REyeTopCoord)

		print str(LEyeCoord), str(REyeCoord)

		imp.close()
		invertedStack.close()
		invertedTopReslice.close()

'''
Obtains the area statistic of the specified stack and saves into temp folder
'''
def getAreaStatistic(imp, outputPath, fishNumber, direction, min_area=0, min_circularity=0.0):
	cacheFilename = '_'.join(["area_fish" + str(fishNumber), str(direction), str(min_circularity).replace('.','p'), str(min_area)])

	outputPath = os.path.join(outputPath, cacheFolderName)

	statFileExt = '.csv'
	headings = []
	statisticsDict = {}

	if not os.path.isfile(os.path.join(outputPath, cacheFilename) + statFileExt):
		paOptions = ParticleAnalyzer.SHOW_PROGRESS + ParticleAnalyzer.CLEAR_WORKSHEET
		paMeasurements = Measurements.AREA + Measurements.STACK_POSITION + Measurements.SLICE + Measurements.CIRCULARITY
		rt = ResultsTable()

		pa = ParticleAnalyzer(paOptions, paMeasurements, rt, min_area, Double.POSITIVE_INFINITY, min_circularity, 1.0)
		pa.setHideOutputImage(True)

		statisticsDict = {}

		for i in range(imp.getStackSize()):
			imp.setSliceWithoutUpdate(i + 1)
			pa.analyze(imp)

		headings = list(rt.getHeadings())
		for header in headings:
			statisticsDict[header] = list(rt.getColumn(rt.getColumnIndex(header)))

		print headings

		rt.reset()

		writeStatistics(outputPath, cacheFilename, statFileExt, statisticsDict, headings)

def test2():
	impBin = IJ.openImage("/Users/rshkarin/Documents/Fishes/fish202_binary/binary_median5_fish202_8bit_416x396x1857.tif")
	getAreaStatistic(impBin, "/Users/rshkarin/Documents/Fishes/fish202_binary", 202, "frontal_meadian5", 1200)

extractEyes(inputPath, outputPath, fishNumbers, fishPrefix, fileExt)
#test()
#test2()





