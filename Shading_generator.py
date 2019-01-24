# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 12:49:00 2018

@author: gkhawam
        gao, song
"""
import numpy as np
import math as m
import datetime as dt
import matplotlib.pyplot as plt
import sys


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


# the solar model is developed by Gao, song                   
def solar_model(time, year, location = None, longitude = None, latitude = None, elevation = None, timezone = None):
    '''
    The function calculates for zenith angle and azimuth angle in radians.
    time: (DOY)
    starting_DOY:starting day for the simulation
    year: the year of the simulation
    location: city, not required, but I set geographic info. for some cities we simulate regularly
    elevation (m):
    longitude (deg):
    latitude (deg):
    timezone: UTC time zone
    '''
    check_list = [location, longitude, latitude, elevation, timezone]
    if check_list.count(None) == len(check_list) :
        print('The function requires longitude, latitude, elevation, timezone \n')
        print('Or input location = Tucson, Mesa')
        sys.exit()
    
    geo_info = {#'city':[latitude, longitude, UCT-tz, elevation]
                'Tucson':[32.28, -110.95, -7, 800],
                'Mesa':[],
                'Pecos':[]}
    
    if location in geo_info:
        latitude = geo_info[location][0]
        longitude = geo_info[location][1]
        timezone = geo_info[location][2]
        elevation = geo_info[location][3]
    else:
        print('Please input latitude, longitude, timezone and elevation')
    lat_r = latitude/180*np.pi
    lon_r = longitude/180*np.pi
    
    n = 0
    for i in range(1900,year):
        if i%4 == 0:
            n += 366
        else:
            n+=365
    JulD = n + time + 2415018.5 - (timezone)/24
    LT = time - int(time)
    JC = (JulD - 2451545) / 36525
    x = 46.815 + JC * (0.00059 - JC * 0.001813)
    M_OE = 23 + (26 + (21.448 - JC * x) / 60) / 60
    EEO = 0.016708634 - JC * (0.000042037 + 0.0000001267 * JC)
    GMAS = 357.52911 + JC * (35999.05029 - 0.0001537 * JC)
    GMAS_r = m.radians(GMAS)
    GMLS = (280.46646 + JC * (36000.76983 + JC * 0.0003032))%360
    GMLS_r = m.radians(GMLS)
    Obliq_C = M_OE + 0.00256 * np.cos((125.04 - 1934.136 * JC) / 180 * np.pi)
    Obliq_C_r = m.radians(Obliq_C)
    
    SEC = np.sin(GMAS_r) * (1.914602 - JC * (0.004817 + 0.000014 * JC)) + np.sin(2 * GMAS_r) * (0.019993 - 0.000101 * JC) + np.sin(3 * GMAS_r) * 0.000289
    STL = GMLS + SEC
    SAL = STL - 0.00569 - 0.00478 * np.sin((125.04 - 1934.136 * JC) / 180 * np.pi)
    SAL_r = m.radians(SAL)
    sin_Delta = np.sin(Obliq_C_r) * np.sin(SAL_r)
    Delta_r = np.arcsin(sin_Delta)     #in radians   
    Var_y = np.tan((Obliq_C / 2) / 180 * np.pi) * np.tan((Obliq_C / 2) / 180 * np.pi)
    EOT_prime = Var_y * np.sin(2 * GMLS_r) - 2 * EEO * np.sin(GMAS_r) + 4 * EEO * Var_y * np.sin(GMAS_r) * np.cos(2 * GMLS_r) - 0.5 * Var_y * Var_y * np.sin(4 * GMLS_r) - 1.25 * EEO * EEO * np.sin(2 * GMAS_r)
    EOT = 4 * EOT_prime / np.pi * 180
            
    TST = (LT * 1440 + EOT + 4 * longitude - 60 * timezone)%1440
    if TST / 4 < 0:
        Omega = TST/4+180
    else:
        Omega = TST/4 - 180
        
    Omega_r = m.radians(Omega)
    cos_Zenith = np.sin(lat_r) * np.sin(Delta_r) + np.cos(lat_r) * np.cos(Delta_r) * np.cos(Omega_r)
    Zenith_r = np.arccos(cos_Zenith)             #in radians
    
    Aprime_r = np.arccos((np.sin(lat_r) * np.cos(Zenith_r) - np.sin(Delta_r)) / (np.cos(lat_r) * np.sin(Zenith_r)))
    Aprime = Aprime_r / np.pi * 180
    
    if Omega > 0:
        Azimuth = (Aprime + 180) % 360   #in degrees
    else:
        Azimuth = (540 - Aprime) % 360   #in degrees                
    Azimuth_r = Azimuth / 180 * np.pi
    
    
    return Zenith_r, Azimuth_r

# if you want to change location you need to fill geometrical info in line 60
location = 'Tucson'
year = 2016
#the orientation of the raceway straight part
orientation_angle = 0   # in degree where 0 is east-west and 90 is North-south
# to control the mesh size
width_inc = 100
ang_inc = 100

# time input parameters
snap_day = 1
snap_month = 1
snap_hour = 12
snap_minute = 0
snap_second = 0

Pi = np.pi
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

Num_points = 2 * ang_inc + 2 * width_inc
orientation_angle = m.radians(orientation_angle)

date=dt.datetime(year,snap_month,snap_day)
DOY= date.timetuple().tm_yday
DOY = DOY + snap_hour/24 + snap_minute/24/60 + snap_second/24/3600

# numerical variable used to 2-D draw the geometry mesh
finite_angle = Pi / 2 / (ang_inc)
finite_base = 2 * radius * np.sin(finite_angle)

# first value of solar_angles is zenith, 2nd is azimuth
solar_angles = solar_model(DOY,year,location)
Zenith = solar_angles[0]
Azimuth_N = solar_angles[1]     # Azimuth with north is the reference going east

Elev_angle = Pi/2 - Zenith
Azimuth_S = Pi + Azimuth_N      # Azimuth with the south is the reference and going west


delta_y = (roughness / np.tan(Elev_angle)) * np.cos(Azimuth_S)      #in inch
delta_x = (roughness / np.tan(Elev_angle)) * np.sin(Azimuth_S)      #in inch

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
shading_percent = (total_shade/Area)*100
shading_percent = round(shading_percent,2)

print ('shading percentage = ' +str(shading_percent) +' %')   
    
plt.figure(1)
plt.scatter(x_geo,y_geo, label = 'tub geometry') 
plt.scatter(x_shade,y_shade, label = 'shade')
plt.scatter(x_medium,y_medium, label = 'Middle bar')
plt.scatter(x_medium_shade,y_medium_shade, label = 'Middle bar shade')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
plt.ylabel('N-S direction')
plt.xlabel('E-W direction')

plt.figure(2)
plt.scatter(x_polygon,y_polygon, label = 'Polygon')
plt.title('Wall shaded area') 
plt.ylabel('N-S direction')
plt.xlabel('E-W direction')
plt.show()
