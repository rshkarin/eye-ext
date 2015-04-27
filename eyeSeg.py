from ij import IJ, ImagePlus, ImageStack
from ij.io import TiffDecoder, FileSaver
from ij.plugin.frame import RoiManager
from ij.process import ImageStatistics, ByteProcessor, FloatProcessor, ShortProcessor, ImageConverter
from ij.plugin.filter import ParticleAnalyzer, Binary
import ij.plugin.filter.PlugInFilter;
from ij.measure import ResultsTable, Measurements
from java.lang import Double
from java.lang import Character
import os, sys, datetime, math, csv, codecs, cStringIO
from fiji.threshold import Auto_Local_Threshold as ALT
from ij import Prefs

Prefs.blackBackground = True

inputPath = "/Users/rshkarin/Documents/Fishes"
outputPath = "/Users/rshkarin/Documents/Eyes"

fishPrefix = "fish"
fileExt = ".tif"
fishNumbers = ["202"]

def performReslice(imp, direction):
	newStack = imp.duplicate()

	IJ.run(newStack, "Reslice [/]...", "output=1.000 start=" + direction + " avoid")

	return IJ.getImage()

def test(inputPath, outputPath, fishPrefix, fileExt, fishNumbers):
	for fishNumber in fishNumbers:
		currentPath = os.path.join(inputPath, fishPrefix + fishNumber)
		
		files = [f for f in os.listdir(currentPath) if os.path.isfile(os.path.join(currentPath,f)) and f.startswith(fishPrefix + fishNumber) and f.endswith(fileExt)]
		print files
		currentFileName = ''
		
		if files:
			currentFileName = files[0]
			pathToVolume = os.path.join(currentPath, currentFileName)

			print "Opening"

			imp = IJ.openImage(pathToVolume)

			if imp is None:
				print "Could not open image from file:", currentFileName
				continue

			td = TiffDecoder(currentPath, currentFileName)
			stackSize = imp.getStackSize()
			stackDims = imp.getDimensions()

			print "Eyes extracting"

			LEyeCentroidCoord = (55.,256.,1698.)
			#REyeCentroidCoord = (0.,0.,0.)

			LEyeSizes = calculateEyeSizeByPupilRadius(23.0)
			#REyeSizes = calculateEyeSizeByPupilRadius(23.0)

			extractEyeByCentroid(imp, LEyeCentroidCoord, LEyeSizes)
			#REye, REyeMask = extractEyeByCentroid(imp, LEyeCentroidCoord)

def calculateEyeSizeByPupilRadius(pupil_radius, borderEyeStd=(0.25, 0.25, 0.25)):
	width = pupil_radius * 4.5
	height = pupil_radius * 5.9
	depth = pupil_radius * 5.9

	width_std = pupil_radius * borderEyeStd[0]
	height_std = pupil_radius * borderEyeStd[1]
	depth_std = pupil_radius * borderEyeStd[2]

	return (int(width + 2 * width_std), int(height + 2 * height_std), int(depth + 2 * depth_std))

def extractEyeByCentroid(imp, coord, sizes):
	#dims
	stackDims = (imp.getImageStack().getWidth(), imp.getImageStack().getHeight(), imp.getImageStack().getSize())

	#get coords of roi of crop
	x0, y0, z0 = int(coord[0] - sizes[0]/2), int(coord[1] - sizes[1]/2), int(coord[2] - sizes[2]/2)

	print x0, y0, z0

	x0 = 0 if x0 < 0 else x0
	x0 = (stackDims[0] - 1) if x0 >= imp.getImageStack().getWidth() else x0

	y0 = 0 if y0 < 0 else y0
	y0 = (stackDims[1] - 1) if y0 >= stackDims[1] else y0

	z0 = 0 if z0 < 0 else z0
	z0 = (stackDims[2] - 1) if z0 >= stackDims[2] else z0

	#get cropped volume
	cropped_eye = ImagePlus("cropped_eye",imp.getImageStack().crop(x0, y0, z0, sizes[0], sizes[1], sizes[2]))

	resliced_cropped_eye = performReslice(cropped_eye, "Left")

	#eye mask
	mask_eye = resliced_cropped_eye.duplicate()
	resliced_cropped_eye.close()

	binnerFill = Binary()
	binnerClose = Binary()

	binnerFill.setup("fill", None)
	binnerClose.setup("close", None)


	for i in range(mask_eye.getStackSize()):
		mask_eye.setSliceWithoutUpdate(i + 1)
		ALT().exec(mask_eye, "Bernsen", 15, 0, 0, True)

		mask_eye.getImageStack().getProcessor(i + 1).invert()

		binnerFill.run(mask_eye.getImageStack().getProcessor(i + 1))

		for _ in range(5):	
			binnerClose.run(mask_eye.getImageStack().getProcessor(i + 1))

		mask_eye.getImageStack().getProcessor(i + 1).invert()

	saveTiff(mask_eye, outputPath, 'left_fish202.tif')

	#filteredNuc = ImageCalculator().run("AND create", impDilatedNuc, impNuc)

def saveTiff(imp, outputPath, currentFileName):
	fileSaver = FileSaver(imp)
	currentFileName = 'eye_' + currentFileName
	fileSaver.saveAsTiffStack(os.path.join(outputPath, currentFileName))

test(inputPath, outputPath, fishPrefix, fileExt, fishNumbers)