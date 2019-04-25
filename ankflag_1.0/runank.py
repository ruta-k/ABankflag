import numpy as np
import os
import time as tm
from convertfits import *
from inputs import *

#	--------------		Convert inputs to numbers

if (CLEARSCRATCH):
	print('\nClearing scratch directory....\n')
	os.system('rm '+scratchdir+'*')


exmode		=	['baseline', 'uvbin']
flagwhat	=	['vis_ind', 'chan_ind', 'rec_ind', 'vis_block', 'chan_block', 'rec_block']
flagon		=	['mean', 'rms', 'mean_rms']
statused	=	['median', 'mean']
datype		=	['re', 'im', 'am', 'ph']
blkorder	=	['ascending', '', 'descending']

print('Total flagging rounds	=	%d'%len(FLAGPARS))

flagparfile	=	open(scratchdir+'flagpars.pars','w')
flagparfile.write('%d	%d	%d	%d	%d	'%(ANTS, exmode.index(FLAGMODE), ugrids, flagon.index(BASEFLAGMEAN[0])+1, vgrids))
flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(BASEFLAGMEAN[1], BASEFLAGMEAN[2], len(FLAGPARS), THREADS, WRITEOUT, BLOCKPOW, BASEFLAGMEAN[3]))

for flpar in FLAGPARS:
	flagparfile.write('%d	%d	%d	%d	%d	'%(flagwhat.index(flpar[0]), flagon.index(flpar[1])+1, statused.index(flpar[2]), datype.index(flpar[3]), -(blkorder.index(flpar[7])-1)))
	flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(flpar[4], flpar[5], flpar[6]+1, flpar[8], flpar[9], flpar[10], flpar[11]))

flagparfile.close()

#	----------------------------	Convert FITS to binary files
start0	=	tm.time()	

if (CONVERTFITS):

	infile		=	fits.open(infilename)
	data		=	infile[0].data
	
	if (FLAGMODE==exmode[1]):		
		uvfitstobinary(data,scratchdir,ugrids,vgrids,plotuv)

	elif (FLAGMODE==exmode[0]):		
		baselinestobinary(ANTS,data,scratchdir)
			
	else:
		print("Unknown flagging mode !!!!			Please tell me how to execute it ........")
	
	infile.close()	
	print("\nConvertion done in 		%d seconds\n"%(tm.time()-start0))
	
#	------------------------------		Flag data
start1	=	tm.time()

if (DOFLAG):
	status	=	os.system('./ankflag')	
	print("\nFlagging done in 		%d seconds\n"%(tm.time()-start1))
	
#	------------------------------		Convert back binary files to FITS	

if (READBACK):

	infile2		=	fits.open(infilename)		
	data2		=	infile2[0].data	

	if (FLAGMODE==exmode[1]):
		bintofits	=	uvfitsfrombinary(data2,scratchdir,ugrids,vgrids)

	elif (FLAGMODE==exmode[0]):
		baselinesfrombinary(ANTS,data2,scratchdir)

	#	---------------------------------------------------------------------
	#					Plot Baselines 
	#	---------------------------------------------------------------------

	if (SHOWBASE):

		infile	=	fits.open(infilename)
		data	=	infile[0].data

		blid	=	[]
		for a in range (1,ANTS):
			for b in range (a+1,ANTS+1):
				blid.append([a,b,256*a+b])
		blid	=	np.array(blid)
		nbase	=	len(blid)

		print('\nIdeally total baselines =	%d'%nbase)

		flaggingstatus	=	[]
		for i in range (0,nbase):	
			flaggingstatus.append(showbasecomparison(blid[i],data,data2,SHOWTF))

		flaggingstatus	=	np.array(flaggingstatus)
		avgflag			=	np.mean(flaggingstatus,axis=0)
		
		print('\n\nAverage flagging fraction	'),
		for p in range (0,npols):
			print('%.3f %.3f	'%(avgflag[p],avgflag[p+npols])),
		print('\n')
		
		print('\nFlagged data			'),
		for p in range (0,npols):
			print('%.3f	'%((avgflag[p+npols]-avgflag[p])/(1.0-avgflag[p]))),
		print('\n')

		infile.close()

	if (WRITEOUT):
		infile2.writeto(outfilename,output_verify='warn',overwrite=True)

	infile2.close()
	print("\nEverything done in 		%d seconds\n"%(tm.time()-start0))
	#	----------------------------------------------------------------------------

if (CLEARSCRATCH):
	print('\nClearing scratch directory....\n')
	os.system('rm '+scratchdir+'*')






























