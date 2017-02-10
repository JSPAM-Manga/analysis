# Graham West

import numpy  as np
import math   as math
import sys

####### MAIN #######
####### MAIN #######

def main():
	
#	fileBase = sys.argv[1]
#	nSim = int(sys.argv[2])
#	print(fileBase)
	
	fileBase = '1237678620102623480'
        rundir = "/nfshome/gtw2i/Research/JSPAM-master_ForEditing/fortran/"

	fileWrite =  rundir + 'Runs_' + fileBase + '/' + fileBase
	nSim = 30
	nBin = 75
	distPerc = 1.1

        # bin all the files in the directory
        binC, binV, galaxyCenters = ReadAndBin(nSim, fileWrite, nBin, distPerc)

        # open files for the results
	fileSave = 'ErrorData_' + fileBase + '_75x75_Average_05.txt'	
	writeID = open(fileSave, 'w')

        # this opens up the run parameter file for access for future comparisions between the run data
        # and the fits
	paramFile = rundir + fileBase + "_combined.txt"
	paramID = open(paramFile,'r')

	
	for i in range(nSim):
		j = i + 1

                # parse the parameter file information to get the citizen science fitness
		line = paramID.readline()
		lineArr = line.split( ' ' )
		lineArr = lineArr[0].split( ',' )
		CitFit = lineArr[1]
                
                # zero the parameters before calling the subroutine - probably not necessary 
                MSE = 0
                AS = 0
                At = 0
                Ovr = 0
		
                MSE, As, At, Ovr = ErrorFunction( nBin, binC[i,:,:], binV[i,:,:], binV[0,:,:] )		
		writeID.write(str(j) +' '+ str(CitFit) +' '+ str(MSE) +' '+ str(As) +' '+ str(At) +' '+ str(Ovr) +'\n')
	
	writeID.close
        paramID.close



####### MAIN #######
####### MAIN #######

def ErrorFunction( nBin, Ci, Vi, V1 ):
	
	As = 0   # area simulation
	At = 0   # area target
	Ovr = 0  # area overlap
	MSE = 0
	morph = 0
	error = 0		
	W = np.zeros((nBin,nBin))
	
	for i in range(nBin):
		for j in range(nBin):
			
			isOvr = 0
                        targetPixel = V1[j,i]
                        simulationPixel = Vi[j,i]
                        tArea = 0
                        sArea = 0
                        oArea = 0

                        # determine if the pixels are nonzero, and set the increment the areas 
                        if targetPixel <> 0.0:
				tArea = 1

                        if simulationPixel <> 0.0:
                                sArea = 1

                        # if both pixels are non-zero, update the overlap area and the MSE
                        if sArea + tArea == 2:
                                oArea = 1
                                MSE = MSE + (simulationPixel - targetPixel)**2


                        As = As + sArea
                        At = At + tArea
                        Ovr = Ovr + oArea
	
	MSE = MSE/Ovr
	
	# stuff with overlap
	#morph = 1 - (Ovr/(As+At-Ovr))
	#error = morph*MSE
	
	return MSE, As, At, Ovr, 



def ReadAndBin(nFile, fileBase, nBin, distPerc):
	
	binC = np.zeros((nFile,nBin,nBin))
	binV = np.zeros((nFile,nBin,nBin))
	
	filename = fileBase + '_00' + str(1) + '.txt'
	RV = np.loadtxt(filename)
	RV_shape = RV.shape
	RV_All = np.zeros((nFile,RV_shape[0],RV_shape[1]))
	maxAll  = 0
	maxLine = 0
	galaxyCenters = []

	for j in range(nFile):
                i = j + 1
		if( i < 10 ):
			filename = fileBase + '_00' + str(i) + '.txt'
		elif( i < 100 ):
			filename = fileBase + '_0'  + str(i) + '.txt'
		else:
			filename = fileBase + '_'   + str(i) + '.txt'
	
                print filename
                #print(j, end=" ")
                particleData = np.loadtxt(filename)
                nparticles = len(particleData)
		RV_All[j,:,:] = particleData
		maxLine = np.amax((np.absolute(RV_All[j,:,0:2])))
		galaxyCenters.append(particleData[nparticles-1,:])

		if( maxLine > maxAll ):
			maxAll = maxLine

        print()
	
	dist = distPerc*maxAll	
	for j in range(nFile):
                print j
		#print(j, end=" ")
		binC[j,:,:], binV[j,:,:] = BinField(nBin, dist, RV_All[j,:,:])
	print()
	return binC, binV, galaxyCenters


def SaveBins(nFile,fileBase,BIN):
	
	for j in range(nFile):		
		i = j + 1
		if( i < 10 ):
			filename = fileBase + '_00' + str(i) + '.txt'
		elif( i < 100 ):
			filename = fileBase + '_0'  + str(i) + '.txt'
		else:
			filename = fileBase + '_'  + str(i) + '.txt'		
		np.savetxt(filename,BIN[j,:,:])




def BinField(nBin, dist, RV):

	nPts = np.size(RV,0)
	step = 2*dist/nBin
	
	binCnt = np.zeros((nBin,nBin))
	binVel = np.zeros((nBin,nBin))
	
	for i in range(nPts):		
                x = float(RV[i,0])
                y = float(RV[i,1])
                vz = float(RV[i,5])
                
                xmin = -dist
                xmax =  dist
                ymin = -dist
                ymax =  dist

                ii = (x - xmin) / (xmax - xmin) * nBin
                jj = (y - ymin) / (ymax - ymin) * nBin
                
                if ii > 0 and ii < nBin and jj > 0 and jj < nBin:
                        binCnt[jj,ii] = binCnt[jj,ii] + 1
                        binVel[jj,ii] = binVel[jj,ii] + vz
                        
	for i in range(nBin):
		for j in range(nBin):
			if( binCnt[i,j] > 1 ):
				binVel[i,j] = binVel[i,j]/binCnt[i,j]	
	return binCnt, binVel



def print2D(A):
	
	print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in A]))





main()
