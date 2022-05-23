
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:39:08 2022

@author: katherinedelgado and Erin Barnhart
"""
import csv
import alignment 
import numpy as np
import ResponseTools_V4
import skimage.util as sku
import skimage.transform as skt
from pystackreg import StackReg
from PIL import Image, ImageSequence

'''Script to parse images from lif files, align them, and output aligned images and an average projection for masking'''

'''Specify parent directory. This should be a folder that will house all z-plane directories'''

parent_dir = '/Users/katherinedelgado/Desktop/Alignment/' #Change to your parent directory

'''Specify CSV file name. This CSV should be in parent directory'''

csv_filename = parent_dir+'inputs.csv'
data,header = ResponseTools_V4.read_csv_file(csv_filename)

'''Iterate over each line of the CSV where one line contains information from one image or job in the lif file'''
for d in data:
    '''Each bracketed number corresponds to an index (column) of the CSV. 
        If you make changes to the CSV structure change accordingly.''' #Check to make sure your csv columns match the indicies below
    sample = d[0]
    lif_name = d[1]
    job_index = d[2]
    ch1_name = d[3]
    ch1_index = int(d[4])
    use_ch2 = d[5]
    ch2_name = d[6]
    ch2_index = int(d[7])
    use_target = d[8]
    target_name = d[9]
    target_start = int(d[10]) #in frames
    target_stop = int(d[11])
    save_avg = d[12]
    
    print(sample)
    
    #Check all filepaths below to make sure they match yours
    
    '''Parse lif file'''
    lifFileName = parent_dir+'/lif_files/'+lif_name+'.lif'
    image = ResponseTools_V4.loadLifFile(lifFileName)
    job = ResponseTools_V4.getLifImage(image, job_index)
    
    '''Extract ch1 from hyperstack'''
    ch1 = job[0,:,ch1_index,:,:]
    '''Use target for alignment. Either use specified target or make target from ch1'''
    if use_target == 'TRUE':
        target_filename = parent_dir+sample+'/images/'+target_name+'.tif'
        target = ResponseTools_V4.read_tif(target_filename)
        print('using designated target to align')
    else:
        subset_target = ch1[target_start:target_stop]
        print(subset_target.shape)
        target = np.average(subset_target,axis=0)
        print(target.shape)
        ResponseTools_V4.save_tif(target, parent_dir+sample+'/images/'+target_name+'.tif')
        print('using automatically generated average projection to align')
    A1, tmat = ResponseTools_V4.alignMultiPageTiff(target, ch1)
    ResponseTools_V4.saveMultipageTif(ch1, parent_dir +sample+'/images/'+ch1_name+'.tif')
    ResponseTools_V4.saveMultipageTif(A1, parent_dir +sample+'/aligned_images/'+ch1_name+'-aligned.tif')
    
    '''If ch2 exists, align to specified target'''
    if use_ch2 == 'TRUE':
        print ('will do ch2')
        ch2 = job[0,:,ch2_index,:,:]
        A2 = ResponseTools_V4.alignFromMatrix(ch2, tmat)
        ResponseTools_V4.saveMultipageTif(ch2, parent_dir +sample+'/images/'+ch2_name+'.tif')
        ResponseTools_V4.saveMultipageTif(A2, parent_dir +sample+'/aligned_images/'+ch2_name+'-aligned.tif')
    
    '''Make and save aligned average projections'''
    if save_avg == 'TRUE':
        print ('saving average projection for masking')
        AVG1 = np.average(A1, axis=0)
        AVG2 = np.average(A2, axis=0)
        ResponseTools_V4.save_tif(AVG1, parent_dir+ sample+'/aligned_images/'+ch1_name+'-avg.tif')
        ResponseTools_V4.save_tif(AVG2, parent_dir+ sample+'/aligned_images/'+ch2_name+'-avg.tif')
        
    
    
    
    
    
        

    