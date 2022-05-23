#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 10:00:35 2021

@author: katherinedelgado and Erin Barnhart
"""
import glob
import csv
import numpy
from pylab import *
import scipy.ndimage as ndimage
from PIL import Image, ImageSequence
import ResponseClassSimple_EB
import skimage.util as sku
import skimage.transform as skt
from pystackreg import StackReg
import alignment
from readlif.reader import LifFile

def get_file_names(parent_directory,file_type = 'all',label = ''):
	if file_type == 'csv':
		file_names = glob.glob(parent_directory+'/*'+label+'*.csv')
	elif file_type == 'all':
		file_names = glob.glob(parent_directory+'/*'+label+'*')
	return file_names


def read_csv_file(filename, header=True):
	data = []
	with open(filename, newline='') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			data.append(row)
	if header==True:
		out_header = data[0]
		out = data[1:]
		return out, out_header
	else:
		return out

def write_csv(data,header,filename):
	with open(filename, "w") as f:
		writer= csv.writer(f)
		writer.writerow(header)
		for row in data:
			writer.writerow(row)

def read_tif(filename):
    tiff = Image.open(filename)
    return numpy.asarray(tiff)

def read_tifs(filename):
	# Open image as PIL ImageSequence
	tiffs = Image.open(filename)
	# Convert each image page to a numpy array, convert array stack into 3D array
	return numpy.array([numpy.array(page) for page in ImageSequence.Iterator(tiffs)], dtype=np.uint8)

def save_tif(image_array,filename):
	out = Image.fromarray(image_array)
	out.save(filename)


def count_frames(filename,gt_index = 0,threshold=1):
	"""Reads in a stimulus output file and assigns an image frame number to each stimulus frame."""
	rows,header = read_csv_file(filename) #import stim file
	R = numpy.asarray(rows,dtype='float') # convert stim file list to an array
	#set up the output array
	output_array=numpy.zeros((R.shape[0],R.shape[1]+2))
	header.extend(['dt','frames']) 
	#calculate change in voltage signal
	vs = [0]
	vs.extend(R[1:,-1]-R[:-1,-1])
	#count image frames based on the change in voltage signal
	count_on = 0
	F_on = [0]
	count_off = 0
	F_off = [0]
	frame_labels = [0]
	n = 1
	while n<len(vs)-1:
		if vs[n]>vs[n-1] and vs[n]>vs[n+1] and vs[n] > threshold:
			count_on = count_on+1
			F_on.extend([count_on])
			F_off.extend([count_off])
		elif vs[n]<vs[n-1] and vs[n]<vs[n+1] and vs[n] < threshold*-1:
			count_off = count_off - 1
			F_off.extend([count_off])
			F_on.extend([F_on])
		else:
			F_on.extend([count_on])
			F_off.extend([count_off])
		frame_labels.extend([count_on*(count_on+count_off)])
		n=n+1
	frame_labels.extend([0])
	output_array[:,:R.shape[1]] = R
	output_array[:,-1] = frame_labels
	OAS = output_array[output_array[:,-1].argsort()]
	i1 = numpy.searchsorted(OAS[:,-1],1)
	OASc = OAS[i1:,:] #removes rows before image series start
	output_list = []
	n = 0
	frame=1
	same_frame = []
	while n<len(OASc):
		if int(OASc[n,-1])==frame:
			same_frame.append(list(OASc[n,:]))
		else:
			output_list.append(same_frame[0])
			same_frame = []
			same_frame.append(list(OASc[n,:]))
			frame=frame+1
		n=n+1
	output_list.append(same_frame[0])
	OL = numpy.asarray(output_list)
	gt = OL[:,gt_index]
	dt = gt[1:]-gt[:-1]
	OL[1:,-2]=dt
	return OL,R,header

def find_dropped_frames(frames,time_interval,stim_data,stim_data_OG,gt_index):
	stim_frames = stim_data[-1,-1]
	print('number of image frames = ' +str(frames))
	print('number of stimulus frames = '+str(int(stim_frames)))
	if stim_frames != frames:
		print('uh oh!')
		target_T = frames * time_interval
		stim_T = numpy.sum(stim_data[:,-2])
		print('total time should be '+str(target_T)+' seconds')
		print('total time from stim file is '+str(stim_T)+' seconds')
		max_t_step = numpy.max(stim_data[:,-2])
		if numpy.round(max_t_step/time_interval)>=2:
			print('stimulus dropped at least one frame!')
			OUT = []
			num_df = 0
			for row in stim_data:
				if numpy.round(row[-2]/time_interval)>=2:
					num_df=num_df+1
					gt_dropped = row[gt_index]-time_interval
					stim_frame = numpy.searchsorted(stim_data_OG[:,gt_index],gt_dropped)
					print('frame dropped in row '+str(stim_frame)+' of original stim file (maybe)')
					OUT.append(stim_data_OG[stim_frame])
					OUT.append(row[:-2])
				else:
					OUT.append(row[:-2])
			print('found '+str(num_df)+' potential dropped frames')
		else:
			print('stim frames and image frames do not match, but no dropped frames found...double check the stim file')
	else:
		print('looks fine!')


def parse_stim_file(stim_info_array,frame_index = -1,rt_index = 1,st_index = None):
	"""Get frame numbers, global time, relative time per epoch, and stim_state (if it's in the stim_file)"""
	frames = stim_info_array[:,frame_index]
	rel_time = stim_info_array[:,rt_index]
	if st_index == 'None':
		stim_type = list(numpy.ones(len(frames),dtype = int))
	else:
		stim_type = stim_info_array[:,int(st_index)]
	return frames, rel_time, stim_type

def define_stim_state(rel_time,on_time,off_time):
	"""Define stimulus state (1 = ON; 0 = OFF) based on relative stimulus time."""
	stim_state = []
	for t in rel_time:
		if t>on_time and t<off_time:
			stim_state.extend([1])
		else:
			stim_state.extend([0])
	return stim_state

def segment_ROIs(mask_image):
	"""convert binary mask to labeled image"""
	labels = ndimage.measurements.label(mask_image)
	return labels[0]

def generate_ROI_mask(labels_image, ROI_int):
	return labels_image == ROI_int

def measure_ROI_fluorescence(image,mask):
	"""measure average fluorescence in an ROI"""
	masked_ROI = image * mask
	return numpy.sum(masked_ROI) / numpy.sum(mask)

def measure_ROI_ts(images,mask):
	out = []
	for image in images:
		out.append(measure_ROI_fluorescence(image,mask))
	return out

def measure_multiple_ROIs(images,mask_image):
	labels = segment_ROIs(mask_image)
	out = []
	num = []
	n = 1
	while n<=numpy.max(labels):
		mask = generate_ROI_mask(labels,n)
		out.append(measure_ROI_ts(images,mask))
		num.append(n)
		n=n+1
	return out,num,labels

def get_input_dict(row,header):
	input_dict = {}
	for r,h in zip(row,header):
		input_dict[h]=r
	return input_dict

def extract_response_objects(image_file,mask_file,stim_file,input_dict):
	"""inputs are file names for aligned images, binary mask, and unprocessed stimulus file
	outputs a list of response objects"""
	#read files
	I = read_tifs(image_file)
    #add line to smooth images here?
	mask = read_tifs(mask_file)
	labels = segment_ROIs(mask)
	print('number of ROIs = '+ str(numpy.max(labels)))
	#process stimulus file
	stim_data,stim_data_OG,header = count_frames(stim_file)
	if (len(I))!=int(stim_data[-1][-1]):
		print("number of images does not match stimulus file")
		print('stimulus frames = ' + str(int(stim_data[-1][-1])))
		print('image frames = ' + str(len(I)))
		#stim_data = fix_dropped_frames(len(I),float(input_dict['time_interval']),stim_data,stim_data_OG,int(input_dict['gt_index']))
	#get frames, relative time, stimuulus type, and stimulus state from stim data
	fr,rt,st = parse_stim_file(stim_data,rt_index = int(input_dict['rt_index']),st_index = input_dict['st_index'])
	ss = define_stim_state(rt,float(input_dict['on_time']),float(input_dict['off_time']))
	#measure fluorscence intensities in each ROI
	responses,num,labels = measure_multiple_ROIs(I,mask)
	#load response objects
	response_objects = []
	for r,n in zip(responses,num):
		ro = ResponseClassSimple_EB.Response(F=r,stim_time = rt,stim_state = ss,ROI_num = n,stim_type = st)
		ro.sample_name = input_dict['sample_name']
		ro.reporter_name = input_dict['reporter_name']
		ro.driver_name = input_dict['driver_name']
		ro.stimulus_name = input_dict['stimulus_name']
		ro.time_interval = float(input_dict['time_interval'])
		response_objects.append(ro)
	return response_objects, stim_data, header, labels

def segment_individual_responses(response_objects,input_dict):
	for ro in response_objects:
		ro.segment_responses(int(input_dict['frames_before']),int(input_dict['frames_after']))
		ro.measure_dff(int(input_dict['baseline_start']),int(input_dict['baseline_stop']))

def measure_average_dff(response_objects,input_dict):
	for ro in response_objects:
		ro.measure_average_dff(int(input_dict['epoch_length']))

def save_raw_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num]
		out.extend(ro.F)
		OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(ro.F),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def save_individual_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		for st,dff in zip(ro.stim_type_ind,ro.dff):
			out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num,st]
			out.extend(dff)
			OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num','stim_type']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(dff),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def save_average_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		for st,a in zip(ro.stim_type_ind,ro.average_dff):
			out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num,st]
			out.extend(a)
			OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num','stim_type']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(a),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def plot_raw_responses(response_objects,filename):
	for ro in response_objects:
		plot(ro.F)
		savefig(filename+'-'+str(ro.ROI_num)+'-raw.png',dpi=300,bbox_inches = 'tight')
		clf()

def plot_average_responses(response_objects,filename):
	for ro in response_objects:
		for a in ro.average_dff:
			plot(a)
		savefig(filename+'-'+str(ro.ROI_num)+'-average.png',dpi=300,bbox_inches = 'tight')
		clf()

def alignMultiPageTiff(ref, img):
    tmat = []
    aligned_images = []
    for t in img:
        #print(t.shape)
        mat = alignment.registerImage(ref, t, mode="rigid") #transformation matrix
        a = alignment.transformImage(t,mat) #aligned image
        aligned_images.append(a)
        tmat.append(mat)
    A = numpy.asarray(aligned_images)
    return A, tmat

def alignFromMatrix(img, tmat):
    aligned_images = []
    for t,mat in zip(img, tmat):
        a = alignment.transformImage(t, mat)
        aligned_images.append(a)
    A = numpy.asarray(aligned_images)
    return A

def saveMultipageTif(numpyArray, saveFile):
   # use list comprehension to convert 3D numpyArray into 1D pilArray, which is a list of 2D PIL images (8-bit grayscale via mode="L")
   pilArray = [
       Image.fromarray(numpyArray[x], mode="L")
       for x in range(numpyArray.shape[0])]
   # saveFile is a complete string path in which to save your multipage image. Note, saveFile should end with ".tif"
   pilArray[0].save(
       saveFile, compression="tiff_deflate",
       save_all=True, append_images=pilArray[1:])

def loadLifFile(file):
    """
    Load entire lif file as an object from which to extract imaging samples
    @param file: path to .lif file to be loaded
    @type file: string
    @return: LifFile iterable that contains all imaged samples in lif file
    @rtype: readlif.reader.LifFile
    """
    lif = LifFile(file)
    return lif

def getLifImage(lif, idx, dtype=numpy.uint8):
    """
    Extract an imaged sample as a hyperstack from a pre-loaded .lif file
    @param lif: LifFile iterable that contains all imaged samples in lif file
    @type lif: readlif.reader.LifFile
    @param idx: index of desired image within lif file (0-indexed)
    @type idx: int
    @param dtype: data type of image array to be saved
    @type dtype: np.dtype, passed to np.ndarray.astype constructor
    @return: lif imaging sample converted to 5D hyperstack array
    @rtype: numpy.ndarray
    """
    image = lif.get_image(img_n=idx)
    stack = [[[numpy.array(image.get_frame(t=t, c=c))
               for c in range(image.channels)]
             for t in range(image.dims.t)]]
    stack = numpy.array(stack, dtype=dtype)
    return stack










