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
save_dir= "/home/samapti-lakshan/carrington/areathresimage"
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
		
		threshold=200
		binary_mask=np.abs(data)>threshold
		labeled_array, num_features = label(binary_mask)
		
		
		sizes = ndimage.sum(binary_mask, labeled_array, range(1, num_features + 1))
		
		areathres=700
		
		mask_filtered = np.zeros_like(binary_mask)
		
		for i, size in enumerate(sizes):
			if size >=areathres:
				mask_filtered[labeled_array == (i + 1)] = 1
		
		
		
		fig = plt.figure()
		plt.imshow(mask_filtered, origin="lower", cmap='gray',extent=[longitude.min(),longitude.max(),latitude.min(),latitude.max()] )
		plt.colorbar()
		plt.xticks(np.arange(longitude.min(),longitude.max()+1,60))
		plt.yticks(np.arange(latitude.min(),latitude.max()+1,30))
		base_name = os.path.basename(file).replace('.fits', '')
		plt.title(f"{base_name}")
		image_name= f" {base_name}_areathres.jpeg"
		mag_image= os.path.join(save_dir, image_name)
		plt.savefig(mag_image, dpi=300)
		plt.show()
