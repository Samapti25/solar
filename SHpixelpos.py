import os
from astropy.io import fits
import sunpy.map
from sunpy.coordinates import NorthOffsetFrame
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import matplotlib.pyplot as plt
import matplotlib.colors
import glob
from scipy import ndimage
import numpy as np
import sunpy.timeseries
import scipy.ndimage as ndimage
from scipy.ndimage import label

folder_path = "/home/samapti-lakshan/carrington"
save_dir= "/home/samapti-lakshan/carrington/pixelpos"
files = os.listdir(folder_path)
fits_files = glob.glob(os.path.join(folder_path, '*.fits')) 
for file in fits_files:
	with fits.open(file)as f:
		header = f[0].header
		print (header)
		data=f[0].data
		nx=data.shape[1]
		ny=data.shape[0]
		print(nx)
		print(ny)
		longitude=np.linspace(0,360,nx)
		sin_lat=np.linspace(-1,1,ny)
		latitude=np.arcsin(sin_lat)*(180/np.pi)
		
		mask= (latitude<=0)  & (latitude>=-90)
		crop_data= data[mask, :]
		crop_latitude=latitude[mask]
		
		threshold=200
		binary_mask=np.abs(crop_data)>threshold
		labeled_array, num_features = label(binary_mask)
		
		
		sizes = ndimage.sum(binary_mask, labeled_array, range(1, num_features + 1))
		
		areathres=700
		
		mask_filtered = np.zeros_like(binary_mask)
		
		for i, size in enumerate(sizes):
			if size >=areathres:
				mask_filtered[labeled_array == (i + 1)] = 1
		
		
		
		
		y_indices, x_indices = np.where(mask_filtered==1)
		print (y_indices)
		lat=crop_latitude[y_indices]
		lon=longitude[x_indices]
		sorted_indices= np.argsort(lon)
		lon_sorted=lon[sorted_indices]
		lat_sorted= lat[sorted_indices]
		output_lat=np.column_stack((lat_sorted))
		output= np.column_stack((lat_sorted, lon_sorted))
		
		
		base_name = os.path.basename(file).replace('.fits', '')
		txt_name1= f" {base_name}_SHpixelpos_lat.txt"
		txt1= os.path.join(save_dir, txt_name1)
		np.savetxt(txt1,output_lat)
		txt_name2= f" {base_name}_SHpixelpos.txt"
		txt2= os.path.join(save_dir, txt_name2)
		np.savetxt(txt2,output)
		
		
		
		
		
