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

inputPath = "/Users/rshkarin/Documents/Fishes"
outputPath = "/Users/rshkarin/Documents/Eyes"

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
			nonNull = len(filter(lambda x: x > 0, pixels))

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
	x = range(imp.getStackSize())
	y = []
	
	for i in x:
		if i in statistics['Slice']:
			if kind == 0:
				y.append(statistics['Area'])
			elif kind == 1:
				y.append(statistics['X'] * statistics['Y'])
			else:
				y.append(statistics['Area'])

	return x, y

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
def getEyesBoundingArea(imp, sliceIdx, min_area=600):
	paOptions = ParticleAnalyzer.SHOW_PROGRESS + ParticleAnalyzer.CLEAR_WORKSHEET
	paMeasurements = Measurements.AREA + Measurements.RECT
	rt = ResultsTable()

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

	return x, y, width, height

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
		filteredRightReslice = performReslice(filteredStack, "Right")

		printLog("Binarizing and inverting", pathToVolume)
		invertedStack = prepareInvertedBinaryStack(filteredStack, lowZbound, highZbound)
		filteredStack.close()

		invertedTopReslice = prepareInvertedBinaryStack(filteredTopReslice, 0, filteredTopReslice.getStackSize())
		filteredTopReslice.close()

		invertedRightReslice = prepareInvertedBinaryStack(filteredRightReslice, 0, filteredRightReslice.getStackSize())
		filteredRightReslice.close()

		printLog("Saving as tiff stack", pathToVolume)
		saveTiff(invertedStack, outputPath, currentFileName)
		saveTiff(invertedTopReslice, outputPath, 'top_' + currentFileName)
		saveTiff(invertedRightReslice, outputPath, 'right_' + currentFileName)

		imp.close()
		invertedStack.close()
		invertedTopReslice.close()
		invertedRightReslice.close()

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

#extractEyes(inputPath, outputPath, fishNumbers, fishPrefix, fileExt)
#test()
test2()





