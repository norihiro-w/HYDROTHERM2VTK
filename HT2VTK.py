# Copyright (c) 2017 Norihiro Watanabe.
#
# HS2VTK is free software; you can redistribute and use in source and
# binary forms, with or without modification. HS2VTK is distributed in the
# hope that it will be useful but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import sys
import os

#------------------------------------------------------------------
# BEGIN: Common functions
#------------------------------------------------------------------
# Collect time steps in the given output file
def getTimeSteps(f):
	ts = []
	for line in f:
		if line.lstrip().find("Time Step No.") == 0:
			ts.append(int(line.split()[3].split(';')[0]))
	return ts

# Read floating values of the specified variable at the given time step no
def getValueArray(filename, t, varname):
	v = [0] * ncell
	f = open(filename)
	#TODO ny?
	while True:
		line = f.readline()
		if not f:
			break
		if line.find("Time Step No.")==-1:
			continue
		if int(line.split()[3].split(';')[0]) != t:
			continue
		while f:
			line = f.readline()
			if line.find("Time Step No.") != -1:
				break
			if line.find(varname)!=-1:
				break
		if not f:
			break

		line = f.readline()
		if line.find("Maximum value")!=-1:
			for i in range(0, 2):
				line = f.readline()
		else:
			for i in range(0,1):
				line = f.readline()
		n_cols = 0
		while n_cols < nx:
			for i in range(0,nz):
				cs = f.readline().split()
				for j in range(1,len(cs)):
					v[i*nx+ n_cols + (j-1)] = float(cs[j])
			n_cols = n_cols + len(cs)-1
			for i in range(0,4):
				line = f.readline()
		break

	f.close()

	# Correct ordering of 2d array
	v2 = list(v)
	for i in range(0, nz):
		for j in range(0,nx):
			v2[i*nx+j] = v[(nz-i-1)*nx + j]

	return v2


# Another version of reading floating values for some variables, e.g. permeability
def getValueArray2(filename, t, varname):
	v = [0] * ncell
	f = open(filename)
	#TODO ny?
	while True:
		line = f.readline()
		if not f:
			break
		if line.find("Time Step No.")==-1:
			continue
		if int(line.split()[3].split(';')[0]) != t:
			continue
		while f:
			line = f.readline()
			if line.find("Time Step No.") != -1:
				break
			if line.find(varname)!=-1:
				break
		if not f:
			break

		line = f.readline()
		if line.find("Maximum value")!=-1:
			line = f.readline()

		for i in range(0,nz):
			line = f.readline() # read "++ Row "
			j = 0
			while j < nx:
				cs = f.readline().split()
				for cv in cs:
					v[i * nx + j] = float(cv)
					j = j +1

		break

	f.close()

	v2 = list(v)
	for i in range(0, nz):
		for j in range(0,nx):
			v2[i*nx+j] = v[(nz-i-1)*nx + j]

	return v2

# Read an array of velocity values
def getVelocityArray(filename, t, vel_name):
	vx = getValueArray(filename, t, "X " + vel_name)
	vz = getValueArray(filename, t, "Z " + vel_name)

	v = [0.0] * len(vx) * 3
	for i in range(0, len(vx)):
		v[i*3 + 0] = vx[i]
		v[i*3 + 1] = 0.0
		v[i*3 + 2] = vz[i]

	return v

# Read and write scalar values to VTK
def outputValueArray2VTK(filename, t, varname, fs_vtk):
	vec = getValueArray(filename, t, varname)
	fs_vtk.write(varname.replace(" ", "_") + " 1 " + str(len(vec)) + " float\n")
	for v in vec:
		fs_vtk.write(str(v) + "\n")

# Read and write scalar values to VTK
def outputValueArray2VTK2(filename, t, varname, fs_vtk):
	vec = getValueArray2(filename, t, varname)
	fs_vtk.write(varname.replace(" ", "_") + " 1 " + str(len(vec)) + " float\n")
	for v in vec:
		fs_vtk.write(str(v) + "\n")

# Read and write vector values to VTK
def outputVelocityArray2VTK(filename, t, varname, fs_vtk):
	vec = getVelocityArray(filename, t, varname)
	fs_vtk.write(varname.replace(" ", "_") + " 3 " + str(len(vec)/3) + " float\n")
	for i in range(0, len(vec)/3):
		fs_vtk.write(str(vec[i*3]) + " " + str(vec[i*3+1]) + " " + str(vec[i*3+2]) + "\n")
#------------------------------------------------------------------
# END: Common functions
#------------------------------------------------------------------

#------------------------------------------------------------------
# Main program starts here
#------------------------------------------------------------------
argvs = sys.argv
argc = len(argvs)
if (argc != 4):
	print("Usage: + python %s <HYDROTHERM output dir> <HYDROTHERM output file suffix> <VTK file base name>" % argvs[0])
	quit()

# read and check the given arguments
hs_filedir=argvs[1]
hs_file_ex=argvs[2]
vtk_basename=argvs[3]
print ('-> HYDROTHERM output file dir: %s (Output_*.%s)' % (hs_filedir, hs_file_ex))
print ('-> VTK output basename: %s' % vtk_basename)

if not os.path.exists(hs_filedir):
	print ('***error: the given HYDROTHERM output file directory not found')
	quit()

#------------------------------------------------------------------
# read grid data
#------------------------------------------------------------------
hs_output_p = os.path.join(hs_filedir, "Out_pressure." + hs_file_ex)
if not os.path.exists(hs_output_p):
	print ('***error: pressure output not found')
	quit()

print ('-> read Out_pressure to get grid data and time step info')
fs_p = open(hs_output_p)
# skip headers
for i in range(0, 4):
	line = fs_p.readline()
# x
line = fs_p.readline()
cvals = []
for v in fs_p.readline().split():
	cvals.append(float(v))
while True:
	for i in range(0, 3):
		line = fs_p.readline()
	line = fs_p.readline()
	if len(line.strip())==0:
		break
	else:
		for v in fs_p.readline().split():
			cvals.append(float(v))
cxx = cvals
nx = len(cvals)

# y
for i in range(0, 3):
	line = fs_p.readline()
cvals = []
for v in fs_p.readline().split():
	cvals.append(float(v))
while True:
	for i in range(0, 3):
		line = fs_p.readline()
	line = fs_p.readline()
	if len(line.strip())==0:
		break
	else:
		for v in fs_p.readline().split():
			cvals.append(float(v))
cyy = cvals
ny = len(cvals)

# z
for i in range(0, 3):
	line = fs_p.readline()
cvals = []
for v in fs_p.readline().split():
	cvals.append(float(v))
while True:
	for i in range(0, 3):
		line = fs_p.readline()
	line = fs_p.readline()
	if len(line.strip())==0 or line.find("Time")!=-1:
		break
	else:
		for v in fs_p.readline().split():
			cvals.append(float(v))
czz = cvals
nz = len(cvals)
ncell = nx*ny*nz

fs_p.close()

print("nx, ny, nz = %d, %d, %d" % (nx, ny, nz))

fs_p = open(os.path.join(hs_filedir, "Out_pressure." + hs_file_ex))
ts = getTimeSteps(fs_p)
fs_p.close()
print("nr. time steps = %d" % len(ts))

# calculate nodal xx, yy, zz
npt = (nx+1)*(ny+1)*(nz+1)
#TODO assume
lx = 0
ly = 0
lz = 0

xx = []
xx.append(lx)
for cx in cxx:
	h = cx - lx
	xx.append(cx + h)
	lx = cx + h
yy = []
yy.append(ly)
for cx in cyy:
	h = cx - ly
	yy.append(cx + h)
	ly = cx + h
zz = []
zz.append(lz)
for cx in czz:
	h = cx - lz
	zz.append(cx + h)
	lz = cx + h

#------------------------------------------------------------------
# p, T, h, water vel, water mass flux, steam vel, steam mass flux
# water&steam: density, viscosity
# permeability, porosity
#------------------------------------------------------------------
n_celldata_var = 2 # p, T
output_h = os.path.exists(os.path.join(hs_filedir, "Out_enthalpy." + hs_file_ex))
if output_h:
	n_celldata_var = n_celldata_var + 1
	print("found enthalpy for output")
output_v = os.path.exists(os.path.join(hs_filedir, "Out_velocity." + hs_file_ex))
if output_v:
	n_celldata_var = n_celldata_var + 4 # vw, mw, vg, mg
	print("found velocities for output")
output_k = os.path.exists(os.path.join(hs_filedir, "Out_permeability." + hs_file_ex))
if output_k:
	n_celldata_var = n_celldata_var + 2 # kx, ky
	print("found permeability for output")

#------------------------------------------------------------------
# Output
#------------------------------------------------------------------
ts_no = 0
print("writing VTK files... ")
for t in ts:
	sys.stdout.write(str(ts_no) + " ")
	sys.stdout.flush()
	# output vtk header
	vtk_filepath = os.path.join(os.path.abspath(hs_filedir), vtk_basename + str(ts_no) + ".vtk")
	fs_vtk = open(vtk_filepath, 'w')
	fs_vtk.write("# vtk DataFile Version 2.0\n")
	fs_vtk.write("vtk output\n")
	fs_vtk.write("ASCII\n")
	fs_vtk.write("DATASET RECTILINEAR_GRID\n")
	fs_vtk.write("DIMENSIONS " + str(nx+1) + " " + str(ny+1) + " " + str(nz+1) + "\n")
	fs_vtk.write("X_COORDINATES " + str(nx+1) + " float\n")
	for v in xx:
		fs_vtk.write(str(v) + " ")
	fs_vtk.write("\n")
	fs_vtk.write("Y_COORDINATES " + str(ny+1) + " float\n")
	for v in yy:
		fs_vtk.write(str(v) + " ")
	fs_vtk.write("\n")
	fs_vtk.write("Z_COORDINATES " + str(nz+1) + " float\n")
	for v in zz:
		fs_vtk.write(str(v) + " ")
	fs_vtk.write("\n")

	fs_vtk.write("POINT_DATA " + str(npt) + "\n")
	fs_vtk.write("CELL_DATA " + str(ncell) + "\n")
	fs_vtk.write("FIELD FieldData " + str(n_celldata_var) + "\n")

	# p
	p = outputValueArray2VTK(os.path.join(hs_filedir, "Out_pressure." + hs_file_ex), t, "Pressure", fs_vtk)
	# T
	T = outputValueArray2VTK(os.path.join(hs_filedir, "Out_temperature." + hs_file_ex), t, "Temperature", fs_vtk)
	# h
	if output_h:
		h = outputValueArray2VTK(os.path.join(hs_filedir, "Out_enthalpy." + hs_file_ex), t, "Enthalpy", fs_vtk)
	if output_v:
		# v_w
		v_w = outputVelocityArray2VTK(os.path.join(hs_filedir, "Out_velocity." + hs_file_ex), t, "Water Interstitial Velocity", fs_vtk)
		# m_w
		m_w = outputVelocityArray2VTK(os.path.join(hs_filedir, "Out_velocity." + hs_file_ex), t, "Water Mass Flux", fs_vtk)
		# v_g
		v_w = outputVelocityArray2VTK(os.path.join(hs_filedir, "Out_velocity." + hs_file_ex), t, "Steam Interstitial Velocity", fs_vtk)
		# m_g
		m_w = outputVelocityArray2VTK(os.path.join(hs_filedir, "Out_velocity." + hs_file_ex), t, "Steam Mass Flux", fs_vtk)

	# permeability
	if output_k:
		kx = outputValueArray2VTK2(os.path.join(hs_filedir, "Out_permeability." + hs_file_ex), t, "X - Permeability", fs_vtk)
		kz = outputValueArray2VTK2(os.path.join(hs_filedir, "Out_permeability." + hs_file_ex), t, "Z - Permeability", fs_vtk)

	fs_vtk.close()

	ts_no = ts_no + 1

print("done.")
