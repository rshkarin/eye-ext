import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import re
import math

areas = {}
axisDot = {}

sliceIdxs = []
dotXY = []
areas = []
circs = []
centerX = []
centerY = []
indices = []

'''
Estimate a fish orientation - plot slice of total area
Fit curve, the parameters of curve show how the fish oriented
'''

def average(s):
	return sum(s) * 1.0 / len(s)

def variance(x):
	return map(lambda y: (y - average(x)) ** 2, x)

def stdev(x):
	return math.sqrt(average(variance(x)))

def normpdf(x, mu, sigma):
	u = (x-mu)/abs(sigma)
	y = math.exp(-u*u/2.0)
	return y

def calcReciepError(i, array, num_neighbours):
	L2 = (num_neighbours - 1) / 2

	values = array[i - L2:i + L2 + 1]

	if not np.sum(values):
		return 0

	num_neighbours_recalc = num_neighbours if len(values) == num_neighbours else len(values)
	normal_dist = [v * normpdf(v, average(values), stdev(values)) + stdev(values) for v in values]

	err_sum = 0
	for (v,n) in zip(values, normal_dist):
		print "%d: %f - %f, mean: %f, std: %f" % (i, v, n, average(values), stdev(values))
		err_sum += abs(v - n) ** 2

	return 1./err_sum
	#return 1.0/np.sum(np.abs(np.array(values, dtype=np.float64) - np.random.normal(np.mean(values), np.std(values), num_neighbours_recalc)) ** 2)

def calcStd(i, array, num_neighbours):
	L2 = (num_neighbours - 1) / 2

	values = array[i - L2:i + L2 + 1]

	if not np.sum(values):
		return 0

	num_neighbours_recalc = num_neighbours if len(values) == num_neighbours else len(values)

	return 1./np.std(values)

def getLeftRightEyesCentroids(x, y, width, height):
	width2 = width / 2
	LEyeX = width2 / 2
	REyeX = width2 / 2

	height2 = height /2
	LEyeY = height2
	REyeY = height2

	return (x + LEyeX, y + LEyeY), (x + width2 + REyeX, y + REyeY)

'''
0 - Frontal (x, y, z)
1 - Top (x, z, y) with respect to frontal
2 - Right (y, z, x) with respect to frontal
3 - Skip slice
CC - centroid coordinate
'''

'''
Return values
 1 - Left eye
 2 - Right eye
-1 - Unknown
'''
def eyesSide(coord, BBox=(5.0, 158.0, 1601.0, 302.0, 186.0, 160.0), sliceIdx=0, orientation=0, restrictedByDepth=False):
	LEyeCC = ()
	REyeCC = ()

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

def smoothStep(x, stopX, k=20.0):
	return 0.5 + 0.5 * math.tanh((x - stopX)/k)

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

def test2():
	sliceIdxs = []
	dotXY = []
	areas = []
	circs = []
	centerX = []
	centerY = []
	indices = []

	i = -1
	with open('/Users/rshkarin/Documents/eyes.csv', 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter=',')
		for r in rd:
			if re.match("^\d+?\.\d+?$", r[2]) and re.match("^\d+?\.\d+?$", r[3]):
				centerX.append(float(r[2]))
				centerY.append(float(r[3]))
				indices.append(i)

			if r[1].isdigit() and r[8].isdigit():
				areas.append(int(r[1]))
				sliceIdxs.append(int(r[8]))
				dotXY.append(float(r[2]) * float(r[3]))

			i = i + 1

	leftEye = []
	rightEye = []

	sliceNum = 1850
	neighbourSlices = 60
	maxPrelimBBoxZPosResidual = 100

	#smooth all circle data
	kernel = signal.get_window(('gaussian', 9), 9)
	kernel = kernel / sum(kernel)

	smoothAreaValues = np.convolve(areas, kernel)
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

	#print "BBox ares slice(XY): %d" % prelimBBoxZPosXY
	#print "BBox ares slice(Area): %d" % prelimBBoxZPosArea
	#print "BBox ares slice average: %d" % int((prelimBBoxZPosXY + prelimBBoxZPosArea)/2.0)

	validSliceIdxs = range(sliceNum)

	for x, y, idx, sliceIdx in zip(centerX, centerY, indices, sliceIdxs):
		if sliceIdx in validSliceIdxs:
			if eyesSide((x, y)) == 1:
				leftEye.append(idx)
			elif eyesSide((x, y)) == 2:
				rightEye.append(idx)
			else:
				print "Lol wut?"

	lEyeValuesArea = [areas[lIdx] for lIdx in leftEye]
	rEyeValuesArea = [areas[rIdx] for rIdx in rightEye]

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

	#lEyeAllowedRange = range(lEyePrelimZPos - neighbourSlices : lEyePrelimZPos + neighbourSlices)
	#rEyeAllowedRange = range(rEyePrelimZPos - neighbourSlices : rEyePrelimZPos + neighbourSlices)

	lEyeAllowedRange = leftEye[((lEyeValuesXYFiltMaxIndex - neighbourSlices) if (lEyeValuesXYFiltMaxIndex - neighbourSlices) >= 0 else 0) : 
								(lEyeValuesXYFiltMaxIndex + neighbourSlices) if (lEyeValuesXYFiltMaxIndex + neighbourSlices) < len(leftEye) else len(leftEye) - 1]

	rEyeAllowedRange = rightEye[((rEyeValuesXYFiltMaxIndex - neighbourSlices) if (rEyeValuesXYFiltMaxIndex - neighbourSlices) >= 0 else 0) : 
								(rEyeValuesXYFiltMaxIndex + neighbourSlices) if (rEyeValuesXYFiltMaxIndex + neighbourSlices) < len(rightEye) else len(rightEye) - 1]

	lEyeAreaAllowed = [areas[i] for i in lEyeAllowedRange]
	rEyeAreaAllowed = [areas[i] for i in rEyeAllowedRange]

	#Second pass - detect end eyes position
	lEyeEndXPos = centerX[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndXPos = centerX[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	lEyeEndYPos = centerY[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndYPos = centerY[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	lEyeEndZPos = sliceIdxs[leftEye[lEyeValuesArea.index(max(lEyeAreaAllowed))]]
	rEyeEndZPos = sliceIdxs[rightEye[rEyeValuesArea.index(max(rEyeAreaAllowed))]]

	#print lEyePrelimZPos, rEyePrelimZPos
	#print lEyeEndZPos, rEyeEndZPos

	print "LEye frontal coord " + str((lEyeEndXPos, lEyeEndYPos, lEyeEndZPos))
	print "REye frontal coord " + str((rEyeEndXPos, rEyeEndYPos, rEyeEndZPos))

	plt.figure()
	plt.plot(lEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(lEyeValuesArea)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesArea)
	plt.show()
	
	'''
	plt.figure()
	plt.plot(lEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(lEyeValuesXY)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesXY)
	plt.show()
	'''

def test3():
	sliceIdxs = []
	dotXY = []
	areas = []
	circs = []
	centerX = []
	centerY = []
	indices = []

	i = -1
	with open('/Users/rshkarin/Documents/eyesTop.csv', 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter=',')
		for r in rd:
			if re.match("^\d+?\.\d+?$", r[2]) and re.match("^\d+?\.\d+?$", r[3]):
				centerX.append(float(r[2]))
				centerY.append(float(r[3]))
				indices.append(i)

			if r[1].isdigit() and r[8].isdigit():
				areas.append(int(r[1]))
				sliceIdxs.append(int(r[8]))
				dotXY.append(float(r[2]) * float(r[3]))

			i = i + 1

	leftEye = []
	rightEye = []

	BBox = (5.0, 158.0, 1631.0, 302.0, 186.0, 100.0)
	orientation = 1
	restrictedByDepth = True
	neighbourSlices = 60
	maxPrelimBBoxZPosResidual = 100

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

	print "LEye top coord " + str((lEyeEndXPos, lEyeEndYPos, lEyeEndZPos))
	print "REye top coord " + str((rEyeEndXPos, rEyeEndYPos, rEyeEndZPos))
	
	plt.figure()
	plt.plot(lEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesXYFilt)
	plt.show()

	plt.figure()
	plt.plot(lEyeValuesArea)
	plt.show()

	plt.figure()
	plt.plot(rEyeValuesArea)
	plt.show()

'''
Test of determining fish orientation by fitting polynom
'''
def test5():
	i = -1
	#with open('/Users/rshkarin/Documents/Fishes/fish202_binary/Temp/area_fish202_frontal_0p0_1500.csv', 'rb') as csvfile:
	with open('/Users/rshkarin/Documents/Fishes/fish202_binary/Temp/area_fish202_frontal_meadian5_0p0_1200.csv', 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter=';')
		for r in rd:
			if re.match("^\d+?\.\d+?$", r[0]) or re.match("^\d+?\.\d+?$", r[1]):
				areas.append(float(r[0]))
				circs.append(float(r[1]))

			i = i + 1

	kernel = signal.get_window(('gaussian', 20), 250)
	kernel = kernel / sum(kernel)
	smoothAreaValues = np.convolve(areas, kernel)



	x = range(len(smoothAreaValues))
	xp = np.linspace(0, len(smoothAreaValues), len(smoothAreaValues))
	z = np.polyfit(x, smoothAreaValues, 6)
	z2 = np.polyfit(x, smoothAreaValues[::-1], 6)
	print z
	print z2
	p = np.poly1d(z)
	rp = np.poly1d(z2)

	plt.figure()
	plt.plot(circs)
	plt.show()

	plt.figure()
	plt.plot(areas)
	plt.show()

	plt.figure()
	plt.plot(x, smoothAreaValues, xp, p(xp), '--')
	plt.show()

	plt.figure()
	plt.show()

'''
Test eyes' coordinates interpoaltion
'''
def test6():
	LEyeFrontalCoord = (56.004, 256.842, 1694)
	REyeFrontalCoord = (256.13, 256.433, 1688)

	LEyeTopCoord = (55.903, 1693.078, 257)
	REyeTopCoord = (256.173, 1687.101, 257)

	LEyeFinalCoord = getEyeCoordinate(LEyeFrontalCoord, LEyeTopCoord)
	REyeFinalCoord = getEyeCoordinate(REyeFrontalCoord, REyeTopCoord)

	print "LEyeFinalCoord = " + str(LEyeFinalCoord)
	print "REyeFinalCoord = " + str(REyeFinalCoord)

if __name__ == "__main__":
	#test2()
	#test3()
	test6()
	#test5()

