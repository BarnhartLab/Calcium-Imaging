import numpy
import ResponseClassSimple_EB
import ResponseTools_V4
import os
from pylab import *

parent_dir = '/Users/katherinedelgado/Desktop/Lab_2022/Mi1_Activity/'

input_csv = parent_dir+'inputs.csv'
rows,header = ResponseTools_V4.read_csv_file(input_csv)

for row in rows:
	input_dict = ResponseTools_V4.get_input_dict(row,header)
	#print(input_dict)
	print('checking stim file for sample ' + input_dict['sample_name'] + ' ' + input_dict['reporter_name'] + ' ' + input_dict['stimulus_name'])
	sample_dir = parent_dir+input_dict['sample_name']
	image_dir = sample_dir+'/aligned_images/'
	stim_dir = sample_dir+'/stim_files/'
	output_dir = stim_dir+'parsed/'
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	stim_file = ResponseTools_V4.get_file_names(stim_dir,file_type = 'csv',label = input_dict['stimulus_name'])[0]
	if input_dict['aligned']=='TRUE':
		image_file = image_dir + input_dict['ch1_name']+'-aligned.tif'
	else:
		image_file = image_dir + input_dict['ch1_name']+'.tif'

	I = ResponseTools_V4.read_tifs(image_file)

	frames = len(I)
	time_interval = float(input_dict['time_interval'])
	gt_index = int(input_dict['gt_index'])

	stim_data,stim_data_OG,dataheader = ResponseTools_V4.count_frames(stim_file,gt_index = gt_index)
	ResponseTools_V4.find_dropped_frames(frames,time_interval,stim_data,stim_data_OG,gt_index)

	if input_dict['verbose']=='TRUE':
		ResponseTools_V4.write_csv(stim_data,dataheader,stim_dir + 'parsed/parsed-'+input_dict['stimulus_name'] +'.csv')