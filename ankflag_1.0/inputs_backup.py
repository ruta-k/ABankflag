infilename		=	'/data/apurba/13064273_2nd/uvsub_avg.fits'		#	Name of the input FITS file
outfilename		=	'/data/apurba/13064273_2nd/uvsub_avg_out.fits'	#	Name of the output file
scratchdir		=	'scratch/'			#	Scratch directory with at least twice the size of the inputs FITS file of space

CLEARSCRATCH	=	True											#	Clear the scratch directory
ukey			=	'UU---SIN'
vkey			=	'VV---SIN'
wkey			=	'WW---SIN'
npols			=	2												#	Number of polarizations

ugrids			=	30
vgrids			=	30
plotuv			=	False

CONVERTFITS		=	1											#	Convert FITS to binary ?	
DOFLAG			=	1											#	Do flagging ?
ANTS			=	30											#	Total number of GMRT antennas
THREADS			=	60											#	Number of parallel threads
FLAGMODE		=	'baseline' 									#	'baseline' / 'uvbin' 
BASEFLAGMEAN	=	['mean_rms',	1.3,	1.2,	0.01]		#	[FLAGON,	tolerance_mean,	tolearnce_rms,	min fraction]]			ONLY for 'baseline'
READBACK		=	1											#	Read back baselines ?
SHOWBASE		=	1											#	SHOW baseline stats ?
SHOWTF			=	0											#	Show time-frequency plots ?
WRITEOUT		=	1											#	Write output ?
BLOCKPOW		=	0.8											#	Power low for Block non-Gaussianity (DON'T CHANGE UNLESS YOU KNOW WHAT IT IS !)

#	For bandpass calibrated data
'''
FLAGPARS		=	[	[ 'chan_ind',	'mean',		'median',	'am',	1.4,	0.2,	2,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'median',	'am',	1.3,	0.2,	2,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'median',	'am',	1.3,	0.1,	2,	'',				0,	0,	0.0,	0.0],
						[ 'vis_ind',	'mean',		'mean',		'am',	1.3,	0.0,	2,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'mean',		'am',	1.2,	0.2,	2,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'mean',		'am',	1.2,	0.2,	2,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'mean',		'am',	1.2,	0.1,	2,	'',				0,	0,	0.0,	0.0],
						[ 'chan_block',	'mean',		'mean',		'am',	1.2,	0.5,	2,	'ascending',	2,	0,	0.1,	0.0],
						[ 'rec_block',	'mean',		'mean',		'am',	1.3,	0.5,	2,	'ascending',	0,	2,	0.0,	0.05]	]

'''
#	For continuum subtracted data

FLAGPARS		=	[	[ 'chan_ind',	'mean',		'median',	'am',	1.4,	0.01,	1,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'median',	'am',	1.3,	0.01,	1,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'median',	'am',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'vis_ind',	'mean',		'median',	'am',	1.3,	0.0,	1,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'rms',		'median',	're',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'rms',		'median',	'im',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'vis_ind',	'mean',		'mean',		'am',	1.2,	0.0,	1,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'mean',		'am',	1.2,	0.02,	1,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'mean',		'am',	1.2,	0.02,	1,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'mean',		'am',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'rms',		'mean',		're',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'rms',		'mean',		'im',	1.2,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_block',	'mean',		'median',	'am',	1.2,	0.5,	1,	'ascending',	2,	0,	0.1,	0.0],
						[ 'rec_block',	'mean',		'median',	'am',	1.2,	0.5,	1,	'ascending',	0,	2,	0.0,	0.05],
						[ 'vis_block',	'mean',		'mean',		'am',	1.2,	0.5,	1,	'ascending',	3,	3,	0.05,	0.05],
						[ 'rec_block',	'mean',		'mean',		'am',	1.2,	0.5,	1,	'ascending',	0,	2,	0.0,	0.05],
						[ 'rec_block',	'rms',		'mean',		're',	1.2,	0.5,	1,	'ascending',	0,	2,	0.0,	0.05],
						[ 'rec_block',	'rms',		'mean',		'im',	1.2,	0.5,	1,	'ascending',	0,	2,	0.0,	0.05],
						[ 'vis_ind',	'mean',		'mean',		'am',	1.2,	0.0,	1,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'mean',		'am',	1.1,	0.02,	1,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'mean',		'am',	1.1,	0.01,	0,	'',				0,	0,	0.0,	0.0]	]


#
#	FLAGON			=	'mean'		/	'rms'		/	'mean_rms'
#	STATYPE			=	'mean'		/	'median'
#	DATATYPE		=	're'		/	'im'		/	'am'		(/	'ph'	---- NOT YET SUPPORTED)
#	ORDER			=	'ascending'			(/	'descending' --- NOT YET SUPPORTED)

#	Vis individual	=	[ 'vis_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	0, 				fit order,  '', 	0, 			0, 			0, 				0			 	]

#	Chan individual	=	[ 'chan_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	'',		0,			0,			0,				0			 	]

#	Rec individual	=	[ 'rec_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,			0,	'',		0,			0,			0,				0			 	]

#	Vis block		=	[ 'vis_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	chan_block,	rec_block,	chan_max_frac,	rec_max_frac 	]

#	Chan block		=	[ 'chan_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	chan_block,	0,			chan_max_frac,	0 				]

#	Rec block		=	[ 'rec_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	0,			rec_block,	0,				rec_max_frac 	]







































































