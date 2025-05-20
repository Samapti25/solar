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
from sunpy.map.header_helper import make_heliographic_header   #  Creating Carrington Maps
import scipy.ndimage as ndimage
from scipy.ndimage import label

folder_path = "/home/samapti-lakshan/carringtonmap"
save_dir1= "/home/samapti-lakshan/carringtonmap/image_intensity"
save_dir2= "/home/samapti-lakshan/carringtonmap/binary_mask"
save_dir4= "/home/samapti-lakshan/carringtonmap/label_info"
save_dir5= "/home/samapti-lakshan/carringtonmap/cluster"
save_dir= "/home/samapti-lakshan/carringtonmap/threshold value(500)_SH"

#save_dir= "/home/samapti-lakshan/imshowimage"
files = os.listdir(folder_path)
fits_files = glob.glob(os.path.join(folder_path, '*.fits')) 
#my_timeseries= sunpy.timeseries.TimeSeries(fits_files, allow_errors=True)
#print (my_timeseries)
#print(plt.colormaps())

for file in fits_files:
	#print(file)
	#file_path= os.path.join(folder_path, file)
	with fits.open(file)as f:
		
		#print(data)
		header = f[0].header
		#header["CROTA2"]=0
		#file.flush()
		print (header)
		data=f[0].data
		#print(data)
		#print(data)
		#mask= np.ma.masked_invalid(data)
		#print(mask)
		#print(mask.max())
		#print(np.nanmax(data))
		#print(np.unravel_index(np.argmax(mask), mask.shape))
		
		# For Rotation🗨️
		#rotated_data=np.rot90(data, k=2)
		#print(f,{file})
		#print(f.read())
		# For plotting📈️
		
		#plots.append(data)
		#print(f[0].CROTA2)
		my_map= sunpy.map.Map(file)              #allow_errors=True
		#rotated_map= my_map.rotate()
		#shape = (720, 1440)         # for Carrington 
		
		#carr_header = make_heliographic_header(rotated_map.date, rotated_map.observer_coordinate, shape, frame='carrington')
		#outmap = rotated_map.reproject_to(carr_header)
		# Try for heliographic_stonyhurst😐️
		#all_hpc = sunpy.map.all_coordinates_from_map(rotated_map)
		#all_hgs = all_hpc.transform_to("heliographic_stonyhurst")
		# giving new dimensions
		#new_dimensions = [4000, 4000] * u.pixel
		#my_resampled_map = rotated_map.resample(new_dimensions)
		
		#print (my_map.coordinate_frame)
		#print(my_map.unit)
		#print(my_map.data.min())
		#contours = my_map.find_contours(4000 * u.G)
		#print(data.max())
		#pixel_pos = np.argwhere(my_map.data == my_map.data.min()) * u.pixel
		#hpc_max = my_map.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
		#world=my_map.wcs.pixel_to_world(0,0)
		#print(world)
		#top_right = SkyCoord(360 * u.deg, 60 * u.deg, frame=my_map.coordinate_frame)
		#bottom_left = SkyCoord(0 * u.deg, -60 * u.deg, frame=my_map.coordinate_frame)
		#crop_map = my_map.submap(bottom_left, top_right=top_right)
		
		# mask that is cut out the point
		#mask = my_map.data < my_map.max() * 0.1
		#print(mask)
		#mymap.mask=mask
		#data2 = ndimage.gaussian_filter(my_map.data * ~mask, 14)
		#print(data2)
		#data2[data2 < 100] = 0
		#my_map2 = sunpy.map.Map(data2, my_map.meta)
		#labels, n = ndimage.label(my_map2.data)
		wcs= my_map.wcs
		nx=my_map.data.shape[1]
		ny=my_map.data.shape[0]
		#print(nx)
		#print(ny)
		longitude=np.linspace(0,360,nx)
		#print(longitude)
		sin_lat=np.linspace(-1,1,ny)
		#print(sin_lat)
		latitude=np.arcsin(sin_lat)*(180/np.pi)
		#print(latitude)
		mask=(latitude >=-65) & (latitude<=65)
		mask2=(latitude >=0) & (latitude<=65)
		mask3=(latitude >=-65) & (latitude<=0)
		#print (mask)
		crop_data= data[mask, :]
		#print(crop_data)
		crop_latitude=latitude[mask]
		#print(crop_latitude)
		#print(crop_latitude.min())
		
		
		#print(data.max())
		
		#fig = plt.figure(figsize=(12,4))		
		#ax = fig.add_subplot(projection=my_map)
		
		#ax.set_axisbelow(True)
		#my_map.plot(axes=ax, origin="lower")
		#ax.imshow(data, origin="lower", cmap='Greys',vmax=100, vmin=-100) 
		#ax.plot_coord(world, 'wx', fillstyle='none', markersize=10)
		#ax.coords.grid(False)
		#ax.set_aspect("equal")
		
		
		#rotated_map.draw_contours([ 70] * u.percent, axes=ax)
		#plt.colorbar()
		
		#rotated_map.draw_limb(axes=ax, color="white")
		#rotated_map.draw_grid(axes=ax)
		#ax.set_aspect("equal")
		#ax.set_xlim((-1000*u.arcsec).value,(1000*u.arcsec).value)
		#ax.set_ylim((-1000*u.arcsec).value,(1000*u.arcsec).value)
		
		#ax.contour(labels)
		#ax.invert_xaxis()
		#ax.invert_yaxis()
		#for contour in contours:
			#ax.plot_coord(contour)
		
		
		#plt.imshow(data, origin="lower", cmap='Greys') #, vmax=100, vmin=-100)
		#nx=data.shape[1]
		#ny=data.shape[0]
		#print(ny)
		#longitude=np.linspace(0,360,nx)
		#latitude=np.linspace(-1,1,ny)
		#print(latitude)
		#print(longitude)
		threshold=200
		binary_mask=np.abs(data)>threshold
		#print(binary_mask)
		labeled_array, num_features = label(binary_mask)
		#print("Labeled array (region IDs):")
		#print(labeled_array)
		

		#print(f"Total number of regions: {num_features}")
		

		centers = ndimage.center_of_mass(binary_mask, labeled_array, range(1, num_features + 1))
		#print(centers)
		
		
		value=np.where((np.abs(data)>threshold),data,0)                        #|(data==0)
		
		#print(value)
		y, x =np.where((np.abs(data)>threshold))
		lon_val=longitude[x]
		lat_val=latitude[y]
		#print(lon_val)
		val_vals=data[y,x]
		output=np.column_stack((lon_val,lat_val,val_vals))
		#print(output)

		fig = plt.figure()
		plt.imshow(value, origin="lower",cmap='grey', aspect="auto", extent=[longitude.min(),longitude.max(),latitude.min(),latitude.max()] )      #extent=[longitude.min(),longitude.max(),crop_latitude.min(),crop_latitude.max()]
		plt.colorbar()
		label_info=[]
		cluster=[]
		for j in range(1, num_features + 1):
    			coords = np.column_stack(np.where(labeled_array == j))  # [ [y1, x1], [y2, x2], ... ]
    			for y, x in coords:
        			lon = longitude[x]
        			lat = latitude[y]
        			cluster.append([j, lon, lat])
		cluster = np.array(cluster)
        	
		for i, (y_idx, x_idx) in enumerate(centers):
    			x_idx = int(round(x_idx))
    			y_idx = int(round(y_idx))

    			if 0 <= x_idx < len(longitude) and 0 <= y_idx < len(crop_latitude):
        			x_val = longitude[x_idx]
        			y_val = crop_latitude[y_idx]
        			label_info.append([i + 1, x_val, y_val])


        			#plt.text(x_val, y_val, str(i + 1), color='black', fontsize=6, ha='center', va='center')
		label_info = np.array(label_info)
#extent=[longitude.min(),longitude.max(),crop_latitude.min(),crop_latitude.max()] ,aspect="auto" 
		
		#plt.xticks(np.arange(longitude.min(),longitude.max()+1,60))
		#plt.yticks(np.arange(-60,61,20))
		#plt.xlabel("carrington longitude ", color="green")
		#plt.ylabel("latitude", color="green")

		base_name = os.path.basename(file).replace('.fits', '')
		plt.title(f"{base_name}")
		txt_name= f" {base_name}_label(2000).txt"
		txt_name2= os.path.join(save_dir, txt_name)
		txt_name4= os.path.join(save_dir2, txt_name)
		txt_name5= os.path.join(save_dir4, txt_name)
		txt_name6= os.path.join(save_dir5, txt_name)
		
		#np.savetxt(txt_name2,output)
		#np.savetxt(txt_name4,binary_mask)
		#np.savetxt(txt_name5,label_info)
		#np.savetxt(txt_name6,cluster)
		
		image_name= f" {base_name}_intensity.jpeg"
		mag_image= os.path.join(save_dir1, image_name)
		plt.savefig(mag_image, dpi=300)
		
		
		#plt.figtext(0.6, 0.8, f'Number of regions = {n}', color='red')
		
		#plt.show()
# if used "Greys" in color then it rotated the colorbar i.e. +ve is black and -ve is white
