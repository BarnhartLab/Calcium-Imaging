import numpy
import ResponseClassSimple_EB
import ResponseTools_V4
import os
from pylab import *

parent_dir = '/Users/katherinedelgado/Desktop/Baseline_Testing_2/'

input_csv = parent_dir+'inputs.csv'
rows,header = ResponseTools_V4.read_csv_file(input_csv)

for row in rows:
	input_dict = ResponseTools_V4.get_input_dict(row,header)
	print('analyzing sample ' + input_dict['sample_name'] + ' ' + input_dict['reporter_name'] + ' ' + input_dict['stimulus_name'])

	sample_dir = parent_dir+input_dict['sample_name']
	image_dir = sample_dir+'/aligned_images/'
	mask_dir = sample_dir+'/masks/'
	stim_dir = sample_dir+'/stim_files/'
	plot_dir = sample_dir+'/plots/'
	output_dir = sample_dir+'/measurements/'

	if input_dict['aligned']=='TRUE':
		image_file = image_dir + input_dict['ch1_name']+'-aligned.tif'
	else:
		image_file = image_dir + input_dict['ch1_name']+'.tif'
	mask_file = mask_dir + input_dict['mask_name']+'.tif'
	stim_file = ResponseTools_V4.get_file_names(stim_dir,file_type = 'csv',label = input_dict['stimulus_name'])[0]
	
	response_objects, stim_data, dataheader, labels = ResponseTools_V4.extract_response_objects(image_file,mask_file,stim_file,input_dict)
	ResponseTools_V4.save_raw_responses_csv(response_objects,output_dir+input_dict['ch1_name']+'-raw.csv')
	ResponseTools_V4.plot_raw_responses(response_objects,plot_dir+input_dict['ch1_name'])
	print('extracted '+str(len(response_objects))+' response objects')
	
	ResponseTools_V4.segment_individual_responses(response_objects,input_dict)
	ResponseTools_V4.measure_average_dff(response_objects,input_dict)
	ResponseTools_V4.save_individual_responses_csv(response_objects,output_dir+input_dict['ch1_name']+'-individual.csv')
	ResponseTools_V4.save_average_responses_csv(response_objects,output_dir+input_dict['ch1_name']+'-average.csv')
	ResponseTools_V4.plot_average_responses(response_objects,plot_dir+input_dict['ch1_name'])

	if input_dict['verbose']=='TRUE':
		ResponseTools_V4.write_csv(stim_data,dataheader,stim_dir + 'parsed/parsed-'+input_dict['stimulus_name'] +'.csv')