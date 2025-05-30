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

folder_path = "/home/anu/Desktop/samapti/carrington"
save_dir= "/home/anu/Desktop/samapti/carrington/pixelthresimage"
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
		
		fig = plt.figure()
		plt.imshow(binary_mask, origin="lower", cmap='gray',extent=[longitude.min(),longitude.max(),latitude.min(),latitude.max()] )
		plt.colorbar()
		plt.xticks(np.arange(longitude.min(),longitude.max()+1,60))
		plt.yticks(np.arange(latitude.min(),latitude.max()+1,30))
		base_name = os.path.basename(file).replace('.fits', '')
		plt.title(f"{base_name}")
		image_name= f" {base_name}_thres.jpeg"
		mag_image= os.path.join(save_dir, image_name)
		plt.savefig(mag_image, dpi=300)
		plt.show()
