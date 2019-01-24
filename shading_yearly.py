# -*- coding: utf-8 -*-
"""
Created on Wed May 23 09:21:41 2018

@author: gkhawam
See shading_generator for remarks about the solar model
"""
import numpy as np
import math as m
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import sys
import scipy.linalg

def two_intersection(x1serie,y1serie,x2serie,y2serie,num_points):
# This function find the two smallest distances between serie1 and serie2 (the intersection points hopefully)
# possible numerical problems especially if we have a coarse mesh
# sometimes it find 2 consective points to have the 2 lowest disatnces 
    
    m1, m2 = 100000, 100000
    for i in range(num_points):
        for j in range(num_points):
            Distance = ((((x2serie[j] - x1serie[i]) ** 2) + (y2serie[j] - y1serie[i]) ** 2)) ** 0.5
            if Distance <= m1:
                m1,m2 = Distance,m1
                idx_geo1,idx_shade1 = i,j
            elif Distance < m2:
                m2 = Distance
                idx_geo2,idx_shade2 = i,j
            
    m = m1
    if m2 > m:
        m = m2
    # m is the distance between the intersection point (the larger value)
    # it will be used as a numerical threshold 
    return idx_geo1,idx_shade1,idx_geo2,idx_shade2,m

Pi = np.pi
# if you want to change location you need to fill geometrical info in line 60
location = 'Tucson'
year = 2016
#the orientation of the raceway straight part
orientation_angle = 45   # in degree where 0 is east-west and 90 is North-south
# to control the mesh size
width_inc = 100
ang_inc = 100

# simulation period and calculating the number of time steps
start_doy = 1
end_doy = 365
Hour_inc = 1    # 1 inc per hour
t_inc = Hour_inc/24
time_steps = int((end_doy - start_doy)/t_inc)

# Geometry input parameters, you can use any unit you want, the units here are in inch
# the distance between the water surface and the top of the paddlewheel     
roughness = 7.5 
# the thickness of the middle bar
medium_thickness = 1.5  
# the length of the straight part of the paddlewheel racew
width = 77.75
# the inside radius of the semi circular part of the raceway
radius = 22.125 
# water surface area, square inch     
Area = (2*radius-medium_thickness)*width + Pi*radius**2    
# the length of one blade of the paddlewheel parallel to the straight part of the raceway
paddle_width = 13.875
# the depth of the the empty raceay, not used in this simulation
tub_depth = 13.5
max_mid_shade = width*radius

Num_points = 2 * ang_inc + 2 * width_inc
orientation_angle = m.radians(orientation_angle)

# numerical variable used to 2-D draw the geometry mesh
finite_angle = Pi / 2 / (ang_inc)
finite_base = 2 * radius * np.sin(finite_angle)

DOY_rec = np.array([])
delta_y_rec = np.array([])
delta_x_rec = np.array([])
wall_shaded_area = np.array([])


for i in range(time_steps):
    
    DOY = start_doy + i*t_inc
    
    # first value of solar_angles is zenith, 2nd is azimuth
    solar_angles = solar_model(DOY,year,location)
    Zenith = solar_angles[0]
    Azimuth_N = solar_angles[1]     # Azimuth with north is the reference going east
    
    Elev_angle = Pi/2 - Zenith
    Azimuth_S = Pi + Azimuth_N      # Azimuth with the south is the reference and going west
    
    
    delta_y = (roughness / np.tan(Elev_angle)) * np.cos(Azimuth_S)      #in inch
    delta_x = (roughness / np.tan(Elev_angle)) * np.sin(Azimuth_S)      #in inch
    
    # the absolute value for delta_x
    abs_delta_x = np.abs(delta_x)
    
    # this if statment restricts when the  result of the algorithm is recorded
    # 16 inch was an experimental number that the algorithm gives resonable results below it
    
    if Elev_angle > 0 and Elev_angle < Pi/2 and  abs_delta_x < 16:  
        # initializing variables
        x_geo = np.ones(Num_points)
        y_geo = np.ones(Num_points)
        x_shade = np.ones(Num_points)
        y_shade = np.ones(Num_points)
        x_medium = np.ones(width_inc)
        y_medium = np.ones(width_inc)
        
        # initial values
        x_geo[0] = 0
        y_geo[0] = 0
        x_shade[0] = delta_x
        y_shade[0] = delta_y
        x_medium[0] = x_geo[0] - (radius + medium_thickness/2)*np.cos((Pi/2)-orientation_angle)
        y_medium[0] = y_geo[0] + (radius + medium_thickness/2)*np.sin((Pi/2)-orientation_angle)
        
        j = 2 # this a counter used in the circular angles calculation
        
        for i in range(width_inc - 1):
            x_medium[i+1] = x_medium[i] + (width/(width_inc+1))*np.cos(orientation_angle)
            y_medium[i+1] = y_medium[i] + (width/(width_inc+1))*np.sin(orientation_angle)
            
        for ii in range(Num_points-1):
            if ii <= width_inc:
                y_geo[ii+1] = y_geo[ii] + (width / (width_inc - 1)) * np.sin(orientation_angle)
                x_geo[ii+1] = x_geo[ii] + (width / (width_inc - 1)) * np.cos(orientation_angle)
            elif ii > width_inc and ii <= (width_inc + ang_inc):
                y_geo[ii+1] = y_geo[ii] + finite_base * np.sin(j * finite_angle + orientation_angle)
                x_geo[ii+1] = x_geo[ii] + finite_base * np.cos(j * finite_angle + orientation_angle)
                j = j+2
            elif ii > (width_inc + ang_inc) and ii < (2*width_inc + ang_inc +1):
                y_geo[ii+1] = y_geo[ii] - (width / (width_inc - 1)) * np.sin(orientation_angle)
                x_geo[ii+1] = x_geo[ii] - (width / (width_inc - 1)) * np.cos(orientation_angle)
            else:
                y_geo[ii+1] = y_geo[ii] + finite_base * np.sin(j * finite_angle + orientation_angle)
                x_geo[ii+1] = x_geo[ii] + finite_base * np.cos(j * finite_angle + orientation_angle)
                j = j+2
           
            #x_geo[-1] = x_geo[0]
            #y_geo[-1] = y_geo[0]
            
            x_shade = x_geo + delta_x
            y_shade = y_geo + delta_y
            x_medium_shade = x_medium + delta_x
            y_medium_shade = y_medium + delta_y
        
        # finding the two intersection points using the function defined above   
        [idx_geo_int1,idx_shade_int1,idx_geo_int2,idx_shade_int2,distance] = two_intersection(x_geo,y_geo,x_shade,y_shade,Num_points)
        
        # rearranging so the first intersection point is always the first one in terms of geometry index order
        if idx_geo_int1 > idx_geo_int2:
            idx_geo_int1,idx_geo_int2 = idx_geo_int2, idx_geo_int1
            idx_shade_int1,idx_shade_int2 = idx_shade_int2,idx_shade_int1
            
        x_polygon = np.array([])
        y_polygon = np.array([])
        
        # this is a counter for the polygon index
        #ii = 0
        # this is a counter for the indices
        jj = idx_geo_int1
        while (jj != idx_geo_int2):
            x_polygon =np.append(x_polygon, x_geo[jj])
            y_polygon =np.append(y_polygon, y_geo[jj])
            jj = jj - 1
            
            if jj == -1:
                jj = Num_points -1
                
        jj = idx_shade_int2
        while (jj != idx_shade_int1):
        
            x_polygon =np.append(x_polygon, x_shade[jj])
            y_polygon =np.append(y_polygon, y_shade[jj])
            jj = jj + 1 
            
            if jj == Num_points:
                jj= 0
        
        # this part apply the sholace formula (explained in wikipedia)
        polygon_shade = 0
        for i in range(len(x_polygon)):
            #if i == len(x_polygon-1):
               # polygon_shade = polygon_shade + (x_polygon[i] * y_polygon[0]) - (y_polygon[i] * x_polygon[0])
           # else:
           polygon_shade = polygon_shade + (x_polygon[i-1] * y_polygon[i]) - (y_polygon[i-1] * x_polygon[i])
        
        # polygon shade is basically the shade resulting from the outer walls of the tub        
        polygon_shade = np.abs(polygon_shade)/2
        #middle shade is the shade resulting from the middle barrier in the raceway
        middle_shade = width *((delta_y*np.cos(orientation_angle))**2 + (delta_x*np.sin(orientation_angle))**2)**0.5
        #paddle shade is the shade resulting from the paddle wheel
        # this is a simplified calculoation assuming that the paddlewheel is not moving and overlay horizontally to the water surface (over prediction)
        # delta_y is always positvie during the day that's why we didn't use abs
        # we multiplied the paddlewidth by 2 because I measured the width of one blade of the paddlewheel and based on our flat assumption..
        # we will have 2 blades parallel to the surface
        # the name width could be misleading but since I named straight part of the raceway (width) that blade measurment was parallel to the straight part of the raceway
        
        paddle_shade =  2 * paddle_width * (radius - (delta_y*np.cos(orientation_angle) + np.abs(delta_x)*np.sin(orientation_angle) - medium_thickness / 2))
        
        total_shade = polygon_shade + middle_shade + paddle_shade
        
        # it was giving me wall shaded area less than 1 in this script
        
        if polygon_shade >= Area or polygon_shade < 100 or middle_shade < 1 or middle_shade > max_mid_shade:
            polygon_shade = 0
            total_shade = Area
        else:    
            DOY_rec = np.append(DOY_rec,DOY)
            delta_y_rec = np.append(delta_y_rec,delta_y)
            delta_x_rec = np.append(delta_x_rec,delta_x)
            wall_shaded_area = np.append(wall_shaded_area,polygon_shade)
    else:
        polygon_shade = 0
        total_shade = Area
        
    shading_percent = (total_shade/Area)*100
    shading_percent = round(shading_percent,2)

shading_script = pd.DataFrame(columns = ['DOY', 'delta_x', 'delta_y', 'wall shaded area'])
shading_script['DOY'] = DOY_rec
shading_script['delta_x'] = delta_x_rec
shading_script['delta_y'] = delta_y_rec
shading_script['wall shaded area'] = wall_shaded_area

x = shading_script['delta_x']
y = shading_script['delta_y']
z = shading_script['wall shaded area']

data=np.c_[x,y,z]

mn = np.min(data,axis=0)
mx = np.max(data,axis=0)
X,Y = np.meshgrid(np.linspace(mn[0],mx[0],100),np.linspace(mn[1],mx[1],100))
XX = X.flatten()
YY = Y.flatten()

# best-fit quadratic curve
A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])

# evaluate it on a grid
Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlabel('Z')
ax.axis('equal')
ax.axis('tight')
plt.show()

print(C)

"""
shading model for the wall
wall shaded area = C[0] + C[1]*delta_x  + C[2]*delta_y + C[3]*delta_x*delta_y + C[4]*(delta_x)^2 + C[5]*(delta_y)^2

middle_shade = width *((delta_y*np.cos(orientation_angle))**2 + (delta_x*np.sin(orientation_angle))**2)**0.5
paddle_shade =  2 * paddle_width * (radius - (delta_y*np.cos(orientation_angle) + np.abs(delta_x)*np.sin(orientation_angle) - medium_thickness / 2))

tot_shade(delta_x,delta_y) = wall + middle + paddle
and delta_x,delta_y are calculated and they are functions of solar angles (time) and water depth 

"""
