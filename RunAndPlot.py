#!/usr/bin/python3
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *
from matplotlib.ticker import StrMethodFormatter
from mpl_toolkits import mplot3d
from os import system, path, remove
import sys
import re


# ---------------------------------------------------------------- 
#   Running this script
# ---------------------------------------------------------------- 

#   This script takes the following arguments:
#   0: run the computation only
#   1: plot the data only
#   2: run computation + plot data


# ---------------------------------------------------------------- 
#   Define parameters
# ---------------------------------------------------------------- 

# Select operation (SEL): 
# SEL = 1 : Mie scattering on a rectangular grid 
# SEL = 2 : Mie scattering on the surface of a sphere 
# SEL = 3 : Cross sections (fixed wavelength, varying radius)
# SEL = 4 : Cross sections (fixed radius, varying wavelength)
# SEL = 5 : Radiation pattern in the xy-, xz- or yz-plane

#-----------------------------------------
# General parameters

SEL     	= 1					# operation selection
fname	= "data.txt"		# Name of data file (for saving data and/or plotting)


#-----------------------------------------
# For a 0D, 1D, 2D, or 3D rectangular grid (SEL = 1)

if SEL == 1:
	r			= 250e-9			# Sphere radius [m]
	lam  		= 1000e-9		# Wavelength [m]
	n			= 60		    		# Number of Mie coefficients
	Lx			= 4000e-9	    # Grid length along x
	Ly			= 0		            # Grid length along y
	Lz			= 4000e-9      	# Grid length along z
	Cx		= 0                 	# Grid center along x
	Cy		= 0                	# Grid center along y
	Cz			= 0                 	# Grid center along z
	res		= 25e-9           # Grid resolution [m]

#-----------------------------------------
# For the surface of a sphere (SEL = 2)

if SEL == 2:
	r			= 250e-9			# Sphere radius [m]
	lam		= 1000e-9		# Wavelength [m]
	n 			= 60					# Number of Mie coefficients
	R			= 250.1e-9		# Points sphere radius [m]
	res 		= 2					# Points angular resolution [°]

#-----------------------------------------
# For the scattering cross section (fixed wavlength, varying radius) (SEL = 3)

if SEL == 3:
	lam		= 1000e-9		# Wavelength [m]
	n 			= 60					# Number of Mie coefficients
	rstart	= 10e-9			# Smallest sphere radius [m]
	rstop 	= 1800e-9 		# Largest sphere radius [m]
	npts 		= 1000				# Number of steps between rstart and rstop

#-----------------------------------------
# For the scattering cross section (fixed radius, varying wavlength) (SEL = 4)

if SEL == 4:
	r			= 250e-9			# Sphere radius [m]
	n 			= 60					# Number of Mie coefficients
	wstart	= 150e-9			# Smallest wavelength [m]
	wstop 	= 8000e-9		# Largest wavelength [m]
	npts	 	= 5000				# Number of steps between wstart and wstop


#-----------------------------------------
# For the radiation patterns in the xy-, xz- and yz-planes (SEL = 5)

if SEL == 5:
	r			= 800e-9			# Sphere radius [m]
	lam		= 1000e-9		# Wavelength [m]
	n 			= 60					# Number of Mie coefficients
	R			= r*1e6			# Points sphere radius [m]
	res 		= 1					# Points angular resolution [°]
	plane	= 2					# plane (1:xy, 2:xz, 3:yz) 









# =========================================================================================
# =========================================================================================
# =========================================================================================



if len(sys.argv) < 2:
	print("Please specify an argument:")
	print("---------------------------")
	print("0: run only")
	print("1: plot only")
	print("2: run + plot")
	print("---------------------------")
	sys.exit()




# ---------------------------------------------------------------- 
#   Run the simulation
# ---------------------------------------------------------------- 

if int(sys.argv[1]) != 1:

	print("Computation starting..")

	if 	path.isfile(fname):
		remove(fname)

	if SEL == 1:
		system("./MiePEC " + str(SEL) + " " + str(lam) + " " + str(n) + " " + str(r) + " " + str(Lx) + " " + str(Ly) + " " + str(Lz) + " " + str(Cx) + " " + str(Cy) + " " + str(Cz) + " " + str(res) + " " + fname)
	if SEL == 2:
		system("./MiePEC " + str(SEL) + " " + str(lam) + " " + str(n) + " " + str(r) + " " + str(R) + " " + str(res) + " " + fname)
	if SEL == 3:
		system("./MiePEC " + str(SEL) + " " + str(lam) + " " + str(n) + " " + str(rstart) + " " + str(rstop) + " " + str(npts) + " " + fname)
	if SEL == 4:
		system("./MiePEC " + str(SEL) + " " + str(r) + " " + str(n) + " " + str(wstart) + " " + str(wstop) + " " + str(npts) + " " + fname)
	if SEL == 5:
		system("./MiePEC " + str(SEL) + " " + str(lam) + " " + str(n) + " " + str(r) + " " + str(R) + " " + str(res) + " " + str(plane) + " " + fname)


if int(sys.argv[1]) == 0:
	sys.exit()


# ---------------------------------------------------------------- 
#   Load and plot the results
# ---------------------------------------------------------------- 



if 	path.isfile(fname) == 0:
	print("File does not exist!")
	sys.exit()

print("Plotting..")

if SEL == 1 or SEL == 2 or SEL == 5:
	
	data = loadtxt(fname, comments='#')

	if size(data) == 15:
		x = data[0]
		y = data[1]
		z = data[2]

		Ex0 = data[3] + 1j*data[4];
		Ey0 = data[5] + 1j*data[6];
		Ez0 = data[7] + 1j*data[8];

		Hx0 = data[9] + 1j*data[10];
		Hy0 = data[11]+ 1j*data[12];
		Hz0 = data[13]+ 1j*data[14];

	else:
		x = data[:,0]
		y = data[:,1]
		z = data[:,2]

		Ex0 = data[:,3] + 1j*data[:,4];
		Ey0 = data[:,5] + 1j*data[:,6];
		Ez0 = data[:,7] + 1j*data[:,8];

		Hx0 = data[:,9] + 1j*data[:,10];
		Hy0 = data[:,11]+ 1j*data[:,12];
		Hz0 = data[:,13]+ 1j*data[:,14];
		

	dim = 0
	if size(unique(x)) > 1:
		dim = dim + 1
	if size(unique(y)) > 1:
		dim = dim + 1
	if size(unique(z)) > 1:
		dim = dim + 1
		

elif SEL == 3:

	data 	= loadtxt(fname, comments='#')

	rv 		= data[:, 0]
	Csca		= data[:, 1]

elif SEL == 4:

	data 	= loadtxt(fname, comments='#')

	wv 		= data[:, 0]
	Csca		= data[:, 1]


# ---------------------------------------------------------------- 
#   Plot the results
# ---------------------------------------------------------------- 

if SEL == 1:
	
	# For a 0D grid
	if dim == 0:
		print("abs(Ex) = {}, abs(Ey) = {}, abs(Ez) = {}".format(abs(Ex0),abs(Ey0),abs(Ez0)))
		print("abs(Hx) = {}, abs(Hy) = {}, abs(Hz) = {}".format(abs(Hx0),abs(Hy0),abs(Hz0)))
		sys.exit()

	# For a 1D grid
	if dim == 1:
		if size(unique(x)) > 1:
			xlab = 'x [m]'
			xc = x
			xmin = min(x)
			xmax = max(x)
		elif size(unique(y)) > 1:
			xlab = 'y [m]'
			xc = y
			xmin = min(y)
			xmax = max(y)
		elif size(unique(z)) > 1:
			xlab = 'z [m]'
			xc = z
			xmin = min(z)
			xmax = max(z)


		rot = 30

		fig, axes = subplots(nrows=2, ncols=3)
		fig.tight_layout() 

		subplot(2,3,1)
		plot(xc,abs(Ex0))
		xlabel(xlab)
		ylabel("abs(Ex)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,2)
		plot(xc,abs(Ey0))
		xlabel(xlab)
		ylabel("abs(Ey)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,3)
		plot(xc,abs(Ez0))
		xlabel(xlab)
		ylabel("abs(Ez)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,4)
		plot(xc,abs(Hx0))
		xlabel(xlab)
		ylabel("abs(Hx)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,5)
		plot(xc,abs(Hy0))
		xlabel(xlab)
		ylabel("abs(Ey)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,6)
		plot(xc,abs(Hz0))
		xlabel(xlab)
		ylabel("abs(Ez)")
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		show()
		sys.exit()


	# For a 2D grid
	if dim == 2:
		if size(unique(x)) > 1 and size(unique(y)) > 1:
			Lx = size(unique(y))
			Ly = size(unique(x))
			xlab = 'x [m]'
			ylab = 'y [m]'
			xmin = min(x)
			xmax = max(x)
			ymin = min(y)
			ymax = max(y)
		elif size(unique(x)) > 1 and size(unique(z)) > 1:
			Lx = size(unique(z))
			Ly = size(unique(x))
			xlab = 'x [m]'
			ylab = 'z [m]'
			xmin = min(x)
			xmax = max(x)
			ymin = min(z)
			ymax = max(z)
		elif size(unique(y)) > 1 and size(unique(z)) > 1:
			Lx = size(unique(z))
			Ly = size(unique(y))
			xlab = 'y [m]'
			ylab = 'z [m]'
			xmin = min(y)
			xmax = max(y)
			ymin = min(z)
			ymax = max(z)


		Ex = reshape(Ex0, (Lx,Ly))
		Ey = reshape(Ey0, (Lx,Ly))
		Ez = reshape(Ez0, (Lx,Ly))

		Hx = reshape(Hx0, (Lx,Ly))
		Hy = reshape(Hy0, (Lx,Ly))
		Hz = reshape(Hz0, (Lx,Ly))

		rot = 30

		fig, axes = subplots(nrows=2, ncols=3)
		fig.tight_layout() 

		subplot(2,3,1)
		imshow(abs(Ex), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Ex)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,2)
		imshow(abs(Ey), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Ey)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,3)
		imshow(abs(Ez), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Ez)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,4)
		imshow(abs(Hx), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Hx)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,5)
		imshow(abs(Hy), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Hy)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		subplot(2,3,6)
		imshow(abs(Hz), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("abs(Hz)")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		figure()
		E = sqrt(abs(Ex)**2 + abs(Ey)**2 + abs(Ez)**2)
		imshow(abs(Ex), extent=[xmin, xmax, ymin, ymax], origin="lower", cmap=cm.jet)
		title("|Etot|")
		colorbar()
		xlabel(xlab)
		ylabel(ylab)
		gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2e}')) 
		xticks(rotation=rot)

		show()
		sys.exit()


	# For a 3D grid
	if dim == 3:
		print("3D grid plotting not implemented!")
		sys.exit()



if SEL == 2:

	E = sqrt(abs(Ex0)**2 + abs(Ey0)**2 + abs(Ez0)**2)

	ax = axes(projection='3d') #, proj_type = 'ortho', aspect='equal' 

	ax.scatter(x, y, z, c=abs(E), cmap=cm.jet, linewidth=0, s=30, alpha=1);

	ax.set_xlabel("x [m]")
	ax.set_ylabel("y [m]")
	ax.set_zlabel("z [m]")
	ax.elev = 0
	ax.azim = 0
	show()



if SEL == 3:

	with open(fname, 'r') as f:
		searchW = re.search('# wavelength = (.*), number', f.readline())
		lam = float(searchW.group(1))

	title("Normalized cross sections")
	plot(2*pi/lam*rv, Csca/(pi*rv**2), linewidth = 2, color='green', label='Csca')
	xlabel("$k r$")
	ylabel("$\sigma/(\pi r^2)$")
	grid()
	legend()

	show()


if SEL == 4:

	with open(fname, 'r') as f:
		searchR = re.search('# radius = (.*), number', f.readline())
		r = float(searchR.group(1))		


	title("Normalized cross sections")
	plot(2*pi/wv*r, Csca/(pi*r**2), linewidth = 2, color='green', label='Csca')
	xlabel("$k r$")
	ylabel("$\sigma/(\pi r^2)$")
	grid()
	legend()

	show()


if SEL == 5:

	Sx = 0.5*real(Ey0*conj(Hz0) - Ez0*conj(Hy0))
	Sy = 0.5*real(Ez0*conj(Hx0) - Ex0*conj(Hz0))
	Sz = 0.5*real(Ex0*conj(Hy0) - Ey0*conj(Hx0))

	if plane == 1:
		angle = arctan2(y,x)
		Sr = Sx*cos(angle) + Sy*sin(angle)		
		tt = "xy-plane (0 deg = +x)"
	elif plane == 2:
		angle = arctan2(x,z)
		Sr = Sz*cos(angle) + Sx*sin(angle)
		tt = "xz-plane (0 deg = +z)"
	elif plane == 3:
		angle = arctan2(y,z)
		Sr = Sz*cos(angle) + Sy*sin(angle)
		tt = "yz-plane (0 deg = +z)"


	polar(angle,10*log10(Sr/Sr.max()), color='blue', linewidth = 2, label=r"$10 log_{10}(I(\theta))$");
	title(tt)
	gca().set_rlim(-25,0)
	gca().set_rticks([0, -5, -10, -15,-20])
	legend()

	show()
