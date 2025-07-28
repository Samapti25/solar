#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from astropy.io import fits
import sunpy.map
#from sunpy.coordinates import NorthOffsetFrame
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
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator


# In[13]:



R_sun = 6.96e10

dtr=np.pi/180
rtd=180/np.pi
folder_path = "/home/anu/Desktop/samapti/cycle24(feb-july)/march"
save_dir1="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
save_dir2= "/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
save_dir3= "/home/anu/Desktop/samapti/cycle24(feb-july)/march/patch"
save_dir4="/home/anu/Desktop/samapti/2025/before shift"

files = os.listdir(folder_path)
fits_files = glob.glob(os.path.join(folder_path, '*.fits'))
for file in fits_files:
    longitude=[]
    longitude_rad=[]
    sinlat=[]
    latitude=[]
    print(fits_files)
    with fits.open(file)as f:
        header = f[0].header
        #print(header)
        data=f[0].data
        newdata=np.roll(data,-600)
        nx=data.shape[1]
        ny=data.shape[0]
        print(f"{nx}  {ny}")
        d_sinlat=2/ny
        dlon_rad=(360/nx)*dtr
        for i in range (ny):
            sin_lat= -1+i*(2/ny)
            Latitude =np.arcsin(sin_lat)*rtd
            Latitude_rad=  (Latitude*dtr)
            latitude.append(Latitude)
            sinlat.append(sin_lat)
        for j in range(nx):
            Longitude = 0+j*(360/nx)
            Longitude_rad=Longitude*dtr 
            longitude.append(Longitude)
            longitude_rad.append(Longitude_rad)
        Flux=(data)*(R_sun**2)*(d_sinlat)*(dlon_rad)
        S= (R_sun**2)*(d_sinlat)*(dlon_rad)
        print(np.nansum(Flux>0))
        print(np.nansum(Flux<0))
        print(S)

        
        
        lon=np.array(longitude)
        lat= np.array(latitude)
        Sinlat=np.array( sinlat)
        x,y=np.meshgrid(lon, lat)
        X,Y=np.meshgrid(lon, Sinlat )
     
        
        fluxthres1=2.3e18
        fluxthres2= -2.3e18
        
        binary_mask=(Flux>fluxthres1)|(Flux<fluxthres2)
        labeled_array, num_features = label(binary_mask)
        totflux=ndimage.sum(Flux,labeled_array, range(1, num_features + 1))
        Totflux= np.zeros_like(Flux)
        for i, val in enumerate(totflux):
            Totflux[labeled_array == (i + 1)] = val
        regionthres1 = 2e20
        regionthres2 = -2e20
        binary_mask2=(Totflux>regionthres1)| (Totflux<regionthres2)
        y_indices, x_indices = np.where(binary_mask2==1)

        Lat=lat[y_indices]
        Lon=lon[x_indices]
        mask1=(Lat<=0)& (Lat>=-90)
        SH_lat=Lat[mask1]
        SH_lon=Lon[mask1]
        mask2= (Lat<=90) & (Lat>=0)
        NH_lat = Lat[mask2]
        NH_lon = Lon[mask2]

        sorted_indices1= np.argsort(SH_lon)

        SHlon_sorted=SH_lon[sorted_indices1]
        SHlat_sorted=SH_lat[sorted_indices1]
        output1=np.column_stack((SHlon_sorted,SHlat_sorted))
        sorted_indices2 = np.argsort(NH_lon)
        NHlon_sorted= NH_lon[sorted_indices2]
        NHlat_sorted = NH_lat[sorted_indices2]
        output2 = np.column_stack((NHlon_sorted, NHlat_sorted))
        base_name = os.path.basename(file).replace('.fits', '')
        txt_name1= f" {base_name}_SH.txt"
        txt1= os.path.join(save_dir1, txt_name1)
        txt_name2= f" {base_name}_NH.txt"
        txt2= os.path.join(save_dir2, txt_name2)
        #np.savetxt(txt1,output1)
        #np.savetxt(txt2, output2)
        
        
        fig = plt.figure(figsize=(12,6))
        plt.imshow(data, origin="lower", cmap='gray', aspect=1.5)              #vmax=100, vmin=-100
        image_name= f" {base_name}.jpeg"
        mag_image= os.path.join(save_dir4, image_name)
        plt.savefig(mag_image, dpi=300)
        plt.close(fig)
        fig = plt.figure(figsize=(12,6))
        plt.imshow(data, origin="lower", cmap='gray', aspect=1.5)              #vmax=100, vmin=-100
        image_name= f" {base_name}.jpeg"
        mag_image= os.path.join(save_dir3, image_name)
        plt.savefig(mag_image, dpi=300)
        plt.close(fig)
        
        fig=plt.figure(figsize=(12,6))
        plt.pcolormesh(x,y,data, cmap='gray')
        image_name= f" {base_name}_flux.jpeg"
        mag_image= os.path.join(save_dir3, image_name)
        #plt.savefig(mag_image, dpi=300)
        plt.close(fig)
        
        
        fig=plt.figure(figsize=(12,6))
        ax=fig.add_subplot(111)                          #.add_subplot(nrows, ncolumns, index or position)
        #ax=fig.add_subplot(122)
        #plt.pcolormesh(X,Y,Flux, cmap='gray')
        pc=ax.pcolormesh(X,Y,Flux, cmap='gray')
        image_name= f" {base_name}_sinflux.jpeg"
        mag_image= os.path.join(save_dir3, image_name)
        #plt.savefig(mag_image, dpi=300)
        ax.grid(True)
        #pc=ax.pcolormesh(X,Y,Flux, cmap='gray')
        #ax.set_xticks(np.arange(0,361,30))
        #ax.set_yticks(np.arange(-1,1,0.5))
        #ax.grid(True)
        plt.close(fig)
        
        
        fig=plt.figure(figsize=(12,6))
        plt.pcolormesh(x,y, binary_mask2, cmap='gray')
        image_name= f" {base_name}_regionmask.jpeg"
        mag_image= os.path.join(save_dir3, image_name)
        #plt.savefig(mag_image, dpi=300)
        plt.close(fig)



# In[ ]:


LON_FRST=   826966.79999980924 / [degree] 
LON_LAST=   827326.69999980927 / [degree] Last Carrington Time
CRLN_OBS=           13.3944387 / [deg] #06.18
CRPIX1  =                1800.
   CRVAL1  =   827146.79999980924
LON_FRST=   826953.49999961851 / [degree]
LON_LAST=   827313.39999961853 / [degree] Last Carrington Time #06.17
CRLN_OBS=   26.631261800000001 / [deg]     


# In[25]:


print(826953.49999961851 -827313.39999961853)
print(-80+120)
print(2250.25%1)
print(  13.3944387 -60)
print(13.3944387-26.631261800000001)
print(26.631261800000001-39.868099200000003 )


# In[19]:


print((46.6055613+13.3944387)/0.10000000000000001)


# In[20]:


newdata=np.roll(data, -600)


# In[23]:


#2025.06.18


fig = plt.figure(figsize=(12,6))
plt.imshow(data, origin="lower", cmap='gray', aspect=1)              #vmax=100, vmin=-100
plt.colorbar()
fig = plt.figure(figsize=(12,6))
plt.imshow(newdata, origin="lower", cmap='gray', aspect=1)              #vmax=100, vmin=-100
plt.colorbar()


# In[22]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))
RMS=[]

for Sfile in txt_files:
    P=[]
    Q=[]
    print(Sfile)
    with open (Sfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                P.append(fhi)
                Q.append(lamda)
    M=np.array(P)
    def model(lon,p,q1,r1):
        S= (lon*np.pi)/180
        return p+q1*np.sin(S+r1)
    Params, Cov = curve_fit(model,M,Q,method = 'trf')
    print(f"p = {Params[0]}, q1 = {Params[1]}, r1={Params[2]}")
    Y=model(M, *Params)
#print(SY)
    SQ=np.subtract(Q,Y)
#print(SQ)
    mean=np.mean(SQ**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    #print(RMS)
#print(Srms)
    Error1=model(lon, *Params)+rms
    Error2=model(lon, *Params)-rms
    
    base_name = os.path.basename(Sfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    ax=fig.add_subplot(111) 
    plt.plot (lon,model(lon, *Params), color="Red" ,linewidth=1)
    plt.plot(lon,Error1, color="green",linewidth=1)
    plt.plot(lon,Error2, color="blue",linewidth=1)
    
    pc=ax.pcolormesh(x,y,Flux,cmap='gray')
    ax.grid(True)

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_STF.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    #plt.savefig(mag_image, dpi=300)
    #plt.close(fig)


# In[119]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Sfile in txt_files:
    P=[]
    Q=[]
    print(Sfile)
    with open (Sfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                P.append(fhi)
                Q.append(lamda)
    M=np.array(P)
    def model(lon,p,q1,r1,q2,r2):
        S= (lon*np.pi)/180
        return p+q1*np.sin(S+r1)+q2*np.sin(2*S+r2)
    Params2, Cov2 = curve_fit(model,M,Q,method = 'trf')
    print(f"p = {Params2[0]}, q1 = {Params2[1]}, r1={Params2[2]},q2={Params2[3]},r2={Params2[4]}")
    Y2=model(M, *Params2)
#print(SY)
    SQ=np.subtract(Q,Y2)
#print(SQ)
    mean=np.mean(SQ**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    #print(RMS)
#print(Srms)
    Error3=model(lon, *Params2)+rms
    Error4=model(lon, *Params2)-rms
    
    base_name = os.path.basename(Sfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *Params2), color="Red" ,linewidth=1)
    plt.plot(lon,Error3, color="green",linewidth=1)
    plt.plot(lon,Error4, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_STF2.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[120]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Sfile in txt_files:
    P=[]
    Q=[]
    print(Sfile)
    with open (Sfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                P.append(fhi)
                Q.append(lamda)
    M=np.array(P)
    def model(lon,p,q1,r1,q2,r2,q3,r3):
        S= (lon*np.pi)/180
        return p+q1*np.sin(S+r1)+q2*np.sin(2*S+r2)+q3*np.sin(3*S+r3)
    Params3, Cov3 = curve_fit(model,M,Q,method = 'trf')
    print(f"p = {Params3[0]}, q1 = {Params3[1]}, r1={Params3[2]},q2={Params3[3]},r2={Params3[4]},q3={Params3[5]},r3={Params3[6]}")
    Y3=model(M, *Params3)
#print(SY)
    SQ=np.subtract(Q,Y3)
#print(SQ)
    mean=np.mean(SQ**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    
#print(Srms)
    Error5=model(lon, *Params3)+rms
    Error6=model(lon, *Params3)-rms
    
    base_name = os.path.basename(Sfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *Params3), color="Red" ,linewidth=1)
    plt.plot(lon,Error5, color="green",linewidth=1)
    plt.plot(lon,Error6, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_STF3.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[121]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Sfile in txt_files:
    P=[]
    Q=[]
    print(Sfile)
    with open (Sfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                P.append(fhi)
                Q.append(lamda)
    M=np.array(P)
    def model(lon,p,q1,r1,q2,r2,q3,r3,q4,r4):
        S= (lon*np.pi)/180
        return p+q1*np.sin(S+r1)+q2*np.sin(2*S+r2)+q3*np.sin(3*S+r3)+q4*np.sin(4*S+r4)
    Params4, Cov4 = curve_fit(model,M,Q,method = 'trf')
    print(f"p = {Params4[0]}, q1 = {Params4[1]}, r1={Params4[2]},q2={Params4[3]},r2={Params4[4]},q3={Params4[5]},r3={Params4[6]},q4={Params4[7]},r4={Params4[8]}")
    Y4=model(M, *Params4)
#print(SY)
    SQ=np.subtract(Q,Y4)
#print(SQ)
    mean=np.mean(SQ**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
#print(Srms)
    Error7=model(lon, *Params4)+rms
    Error8=model(lon, *Params4)-rms
    
    base_name = os.path.basename(Sfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *Params4), color="Red" ,linewidth=1)
    plt.plot(lon,Error7, color="green",linewidth=1)
    plt.plot(lon,Error8, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_STF4.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[122]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/South hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Sfile in txt_files:
    P=[]
    Q=[]
    print(Sfile)
    with open (Sfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                P.append(fhi)
                Q.append(lamda)
    M=np.array(P)
    def model(lon,p,q1,r1,q2,r2,q3,r3,q4,r4,q5,r5):
        S= (lon*np.pi)/180
        return p+q1*np.sin(S+r1)+q2*np.sin(2*S+r2)+q3*np.sin(3*S+r3)+q4*np.sin(4*S+r4)+q5*np.sin(5*S+r5)
    Params5, Cov5 = curve_fit(model,M,Q,method = 'trf')
    print(f"p = {Params5[0]}, q1 = {Params5[1]}, r1={Params5[2]},q2={Params5[3]},r2={Params5[4]},q3={Params5[5]},r3={Params5[6]},q4={Params5[7]},r4={Params5[8]},q5={Params5[9]},r5={Params5[10]}")
    Y5=model(M, *Params5)
#print(SY)
    SQ=np.subtract(Q,Y5)
#print(SQ)
    mean=np.mean(SQ**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    #print(RMS)
#print(Srms)
    Error9=model(lon, *Params5)+rms
    Error10=model(lon, *Params5)-rms
    
    base_name = os.path.basename(Sfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *Params5), color="Red" ,linewidth=1)
    plt.plot(lon,Error9, color="green",linewidth=1)
    plt.plot(lon,Error10, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_STF5.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[123]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))
RMS=[]

for Nfile in txt_files:
    A=[]
    B=[]
    print(Nfile)
    with open (Nfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                A.append(fhi)
                B.append(lamda)
    N=np.array(A)
    def model(lon,a,b1,c1):
        S= (lon*np.pi)/180
        return a+b1*np.sin(S+c1)
    params1, Cov1 = curve_fit(model,N,B,method = 'trf')
    print(f"a= {params1[0]}, b= {params1[1]}, c={params1[2]}")
    Y1=model(N, *params1)
#print(SY)
    NB=np.subtract(B,Y1)
#print(SQ)
    mean=np.mean(NB**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    
    
#print(Srms)
    error1=model(lon, *params1)+rms
    error2=model(lon, *params1)-rms
    
    base_name = os.path.basename(Nfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *params1), color="Red" ,linewidth=1)
    plt.plot(lon,error1, color="green",linewidth=1)
    plt.plot(lon,error2, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_NTF.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[124]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Nfile in txt_files:
    A=[]
    B=[]
    print(Nfile)
    with open (Nfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                A.append(fhi)
                B.append(lamda)
    N=np.array(A)
    def model(lon,a,b1,c1,b2,c2):
        S= (lon*np.pi)/180
        return a+b1*np.sin(S+c1)+b2*np.sin(2*S+c2)
    params2, Cov2 = curve_fit(model,N,B,method = 'trf')
    print(f"a= {params2[0]}, b1= {params2[1]}, c1={params2[2]},b2={params2[3]},c2={params2[4]}")
    Y2=model(N, *params2)
#print(SY)
    NB=np.subtract(B,Y2)
#print(SQ)
    mean=np.mean(NB**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    error3=model(lon, *params2)+rms
    error4=model(lon, *params2)-rms
    
    base_name = os.path.basename(Nfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *params2), color="Red" ,linewidth=1)
    plt.plot(lon,error3, color="green",linewidth=1)
    plt.plot(lon,error4, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_NTF2.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[125]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Nfile in txt_files:
    A=[]
    B=[]
    print(Nfile)
    with open (Nfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                A.append(fhi)
                B.append(lamda)
    N=np.array(A)
    def model(lon,a,b1,c1,b2,c2,b3,c3):
        S= (lon*np.pi)/180
        return a+b1*np.sin(S+c1)+b2*np.sin(2*S+c2)+b3*np.sin(3*S+c3)
    params3, Cov3 = curve_fit(model,N,B,method = 'trf')
    print(f"a= {params3[0]}, b1= {params3[1]}, c1={params3[2]},b2={params3[3]},c2={params3[4]},b3={params3[5]},c3={params3[6]}")
    Y3=model(N, *params3)
#print(SY)
    NB=np.subtract(B,Y3)
#print(SQ)
    mean=np.mean(NB**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    error5=model(lon, *params3)+rms
    error6=model(lon, *params3)-rms
    
    base_name = os.path.basename(Nfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *params3), color="Red" ,linewidth=1)
    plt.plot(lon,error5, color="green",linewidth=1)
    plt.plot(lon,error6, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_NTF3.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[126]:


folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Nfile in txt_files:
    A=[]
    B=[]
    print(Nfile)
    with open (Nfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                A.append(fhi)
                B.append(lamda)
    N=np.array(A)
    def model(lon,a,b1,c1,b2,c2,b3,c3,b4,c4):
        S= (lon*np.pi)/180
        return a+b1*np.sin(S+c1)+b2*np.sin(2*S+c2)+b3*np.sin(3*S+c3)+b4*np.sin(4*S+c4)
    params4, Cov4 = curve_fit(model,N,B,method = 'trf')
    print(f"a= {params4[0]}, b1= {params4[1]}, c1={params4[2]},b2={params4[3]},c2={params4[4]},b3={params4[5]},c3={params4[6]},b4={params4[7]},c4={params4[8]}")
    Y4=model(N, *params4)
#print(SY)
    NB=np.subtract(B,Y4)
#print(SQ)
    mean=np.mean(NB**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    error7=model(lon, *params4)+rms
    error8=model(lon, *params4)-rms
    
    base_name = os.path.basename(Nfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *params4), color="Red" ,linewidth=1)
    plt.plot(lon,error7, color="green",linewidth=1)
    plt.plot(lon,error8, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_NTF4.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[127]:



folder_path3="/home/anu/Desktop/samapti/carrington/synoptic fits files/North hemisphere"
Files = os.listdir(folder_path3)
txt_files = glob.glob(os.path.join(folder_path3, '*.txt'))

for Nfile in txt_files:
    A=[]
    B=[]
 
    with open (Nfile, "r") as infile:
          for line in infile:
            parts=line.split()
            if len(parts)==2:
                fhi=float(parts[0])
                lamda=float(parts[1])
                A.append(fhi)
                B.append(lamda)
    N=np.array(A)
    def model(lon,a,b1,c1,b2,c2,b3,c3,b4,c4,b5,c5):
        S= (lon*np.pi)/180
        return a+b1*np.sin(S+c1)+b2*np.sin(2*S+c2)+b3*np.sin(3*S+c3)+b4*np.sin(4*S+c4)+b5*np.sin(5*S+c5)
    params5, Cov5 = curve_fit(model,N,B,method = 'trf')
    print(f"a= {params5[0]}, b1= {params5[1]}, c1={params5[2]},b2={params5[3]},c2={params5[4]},b3={params5[5]},c3={params5[6]},b4={params5[7]},c4={params5[8]},b5={params5[9]},c5={params5[10]}")
    Y5=model(N, *params5)
#print(SY)
    NB=np.subtract(B,Y5)
#print(SQ)
    mean=np.mean(NB**2)
#print(Smean)
    rms=np.sqrt(mean)
    print(rms)
    RMS.append(rms)
    error9=model(lon, *params5)+rms
    error10=model(lon, *params5)-rms
    
    base_name = os.path.basename(Nfile).replace('.txt', '')
    fig = plt.figure(figsize=(12,6))
    plt.plot (lon,model(lon, *params5), color="Red" ,linewidth=1)
    plt.plot(lon,error9, color="green",linewidth=1)
    plt.plot(lon,error10, color="blue",linewidth=1)
    plt.pcolormesh(x,y,Flux,cmap='gray')

    plt.xlabel("carington longitutde", color="black", fontsize=14, fontstyle="oblique")
    plt.ylabel(" latitude", color="black", fontsize=14, fontstyle="oblique")
    image_name= f" {base_name}_NTF5.jpeg"
    mag_image= os.path.join(save_dir4, image_name)
    plt.savefig(mag_image, dpi=300)
    plt.close(fig)


# In[ ]:




