import os
import Rappture
from Rappture.tools import executeCommand as RapptureExec
import shutil
import math
import numpy as np
import scipy.spatial.distance as scdistance
from tempfile import mkstemp
import sys
import re

def log(msg):
    try:
        driver.put('output.log(output_log)', msg, append=True)
        driver.put('output.log(output_log)', '\n----------\n\n', append=True)
    except NameError:
        pass

		
def memory_check(NAT, render_choice):

	# final value in MB
    user = os.getenv("USER","nobody")
    with open('ddaUser','w') as userQuota:
       userQuota.write('getquota user={0}\n'.format(user))

    with open('ddaUser','r') as userQuota:
       returncode,quotaStdout,quotaStderr = RapptureExec(['nc','-w','5','fshome.nanohub.org','301'],stdin=userQuota,streamOutput=False)
     
    try:
       os.remove('ddaUser')
    except:
       pass

    outline = quotaStdout + quotaStderr
    space_used    = float((outline.split(',')[-3]).split('=')[-1])/float(1024*1000)
    max_hardspace = float((outline.split(',')[-4]).split('=')[-1])/float(1024*1000)
    max_softspace = float((outline.split(',')[-5]).split('=')[-1])/float(1024*1000)
    free_mem = max_hardspace - space_used
    
    # 1 MB is reserved to prevent total hardspace usage. This reserve can be made larger if problems still arise.
    # multiply the NAT's mem-value by 2 since a temp file will be saved as a copy even though it gets deleted quickly.
    
    mem_to_use = (float(NAT) * 0.0000344244233274)*2 + 1
    if render_choice == "yes":
		mem_to_use = mem_to_use * 1.104

    check_space = free_mem - mem_to_use
    
    if (check_space < 0):
		return 0, mem_to_use, free_mem
    elif (check_space >= 0):
		return 1, mem_to_use, free_mem


def get_estimate(pip_path, input_filename):
	# Gives the true estimate of the number of dipoles in the shape
    estNAT = 0
    estInit = 1
    while (estNAT <= 000 and estInit <= 2500):
		estInit = estInit + 15
		run_estInit = estInit + int(estInit/float(10)) + 1
		returncode,pipStdout,pipStderr = RapptureExec([pip_path, '{0}'.format(run_estInit), input_filename, 'estimate_out'],streamOutput=False)
		if returncode != 0:		
			text = pipStdout + pipStderr
			log(text)
			sys.stderr.write(text)
			driver.result(returncode)
#			sys.exit(returncode)
			if ('{0}'.format(returncode) != '0'):
				log(returncode)
				log("PIP failed to generate the shape")

		with open('estimate_out','r') as estfile:
			reedline = estfile.readline()
			while not(re.search('ICOMP',reedline)):		
				prevline = reedline
				reedline = estfile.readline()
				if (re.search('NAT',reedline)):
					estNAT = int(reedline.split()[0])
					ex,why,zee = (int(prevline.split()[-3]),int(prevline.split()[-2]),int(prevline.split()[-1]))

	# ex, why, zee are the shape's initial max dimensions given by 'pip'
	# note that these are not the dipole counts
	# ex2, why2, zee2 are the dipoles found for the shape in each respective direction
	# note that these are not the max dimensions:
	# i.e. the z-value spans from -88 to 88, but the z-value is 117 due to spaces being not counted
    xyrat = ex / float(why)
    yzrat = why / float(zee)
    xzrat = ex / float(zee)
    					
    estPrism= ex * why * zee
    MAXDIST = max(ex,why,zee)
    
    if MAXDIST == ex:
		maxdirect = 'x'
    elif MAXDIST == why:
		maxdirect = 'y'
    elif MAXDIST == zee:
		maxdirect = 'z'
	
    estpercent = (estNAT/float(estPrism))

    return estpercent, estNAT, xyrat, yzrat, xzrat, maxdirect, MAXDIST


def getdist_lightshuttle(input_filename, input_data):
    # Quickly grab the Max Distance Value
    # between Blender shapes from the lightshuttle file.
    #
    # With uploaded files, take the time to recalculate what the explicit MaxVerticeDistance is from the obj file
    
    if input_data == "USE_BLENDER_FILE":
		sessionnum = os.getcwd().split('/')[-1]
		for root, _, files in os.walk(os.getcwd()):
			for fil in files:
				if fil.startswith("PolarLight") == True:
					your_shuttle = os.path.join(root,fil)
					with open(your_shuttle,'r') as blend_file:
						input_data1 = blend_file.readline()
						input_data2 = blend_file.readline()
						input_data3 = blend_file.readline()
						input_data4 = blend_file.readline()

    try:
		distance = float(input_data3.split()[-1])
		in_ex, in_why, in_zee = input_data4.split()[-3], input_data4.split()[-2], input_data4.split()[-1]
		in_ex, in_zee, in_why = float(in_ex), float(in_why), float(in_zee)
    except NameError:
		vert_data_array = []
		with open(input_filename,'r') as infile:
			for line in infile:
				if line.split()[0] == 'v':
					vert_data_array.append((float(line.split()[1]),float(line.split()[2]),float(line.split()[3])))
		points = np.array(vert_data_array)
		dist = scdistance.pdist(points)
		distance = dist.max()
		sqform = scdistance.squareform(dist)
		i = 0
		j = 0
		maxpoint = ""
		for vert in sqform:
			j = 0

			for entry in vert:
				if entry == distance:
					maxpoint = (i,j)
					break
				j = j + 1
			if maxpoint != "":
				break
			i = i + 1
					
		point1 = vert_data_array[maxpoint[0]]
		point2 = vert_data_array[maxpoint[1]]
		in_ex = float(abs(point1[0]-point2[0]))
		in_why = float(abs(point1[1]-point2[1]))
		in_zee = float(abs(point1[2]-point2[2]))
		

    return in_ex, in_why, in_zee




def send_remote(dipoles_in, input_filename, output_filename, multiple, pip_path):

    estPerc, estNAT, xyrat, yzrat, xzrat, maxdirect, maxdist = get_estimate(pip_path, input_filename)
    estPerc, estNAT, xyrat, yzrat, xzrat, maxdist = float(estPerc), int(estNAT), float(xyrat), float(yzrat), float(xzrat), float(maxdist)

    max_dipoles = float(dipoles_in)
    if maxdirect == 'x':
		max_dipoles2 = float(max_dipoles) * (1 / float(xyrat))
		max_dipoles3 = float(max_dipoles) * (1 / float(xzrat))
    elif maxdirect == 'y':
		max_dipoles2 = float(max_dipoles) * xyrat
		max_dipoles3 = float(max_dipoles) * (1 / float(yzrat))
    elif maxdirect == 'z':
		max_dipoles2 = float(max_dipoles) * xzrat
		max_dipoles3 = float(max_dipoles) * yzrat

    max_dipoles  = '{0}'.format(int(round(max_dipoles)))
    max_dipoles2 = '{0}'.format(int(round(max_dipoles2)))
    max_dipoles3 = '{0}'.format(int(round(max_dipoles3)))

    true_estNAT = float(estPerc * (float(max_dipoles) * float(max_dipoles2) * float(max_dipoles3)))
	# assume shape is squished into a cube for memory calculation purposes
    true_estdim = true_estNAT ** (1/float(3))

	# predict memory usage based on dipoles-to-use
	#y = e^(4.021238 + 0.007589x), where x = dipoles in, y = MB of RAM
    memory_usage = math.exp(4.021238 + 0.007589*true_estdim)

    venue_name = 'rcac_M'

    if (true_estNAT > 10000000):
		walltime = '1000'
    elif (true_estNAT > 100000000):
		walltime = '1440'
    else:
		walltime = '800'

 #   elif (memory_usage > 800):
#		walltime = '70'


    saved_home = os.getcwd()
    inpath1 = os.path.join(saved_home, input_filename)
    walltime = '12000'

    command = ["submit"]
    command.append("-M")
    command.append("-v")
    command.append(venue_name)
    command.append("-n")
    command.append("0")
    command.append("-w")
    command.append(walltime)
    if (multiple != 0):
      command.append("-p")
      command.append("@@ID=1-{0}".format(multiple))
    command.append("ddaconvert_r69-pip")
    command.append("{0}".format(dipoles_in))
    if (multiple == 0):
      command.append("{0}".format(inpath1))
      command.append("{0}".format(output_filename))
    else:
      output_filename ='temp_out_shape'
      command.append("{0}/temp_reorder@@ID.obj".format(saved_home))
      command.append("temp_out_shape@@ID")

    if not os.path.exists('submit_results'):
		os.makedirs('submit_results')
    sub_path = os.path.realpath('submit_results')

		       
    os.chdir(sub_path)

    exit_status, stdout, stderr = RapptureExec(command, streamOutput=True)
    os.chdir(saved_home)

    if ('{0}'.format(exit_status) != '0'):
		fail_message = "\nThe remote submission returned unsuccessfully.\n"		
		log(fail_message)

    if (multiple != 0):				
		for root, _, files in os.walk(sub_path):
			for fil in files:
				if fil.startswith("temp_out_shape") == True:
					outpath1 = os.path.join(root,fil)
					outpath2 = os.path.join(saved_home, fil)						
					os.rename(outpath1, outpath2)

	
    else:
		outpath1 = os.path.join(sub_path, output_filename)
		outpath2 = os.path.join(saved_home, output_filename)
		if os.path.exists(outpath1):
			os.rename(outpath1, outpath2)


    return exit_status



def dx_generate(pipoutfile, dx_scale):

    count_array = []
    read_bufferpip = 0
    minx = 99999999
    miny = 99999999
    minz = 99999999
    with open(pipoutfile, 'r') as pof:
		for line in pof:					
			if (re.search('NBX, NBY, NBZ=',line)):
				NBZ = line.split()[-1]
				NBY = line.split()[-2]
				NBX = line.split()[-3]				
			if (re.search('= NAT',line)):
				NAT = line.split()[0]				
			if (read_bufferpip == 1):
				try:					
					if (int(line.split()[1]) < minx):
						minx = line.split()[1]
					if (int(line.split()[2]) < miny):
						miny = line.split()[2]
					if (int(line.split()[3]) < minz):
						minz = line.split()[3]
				except (KeyError,ValueError):
					pass
			if (re.search('ICOMP',line)):
				read_bufferpip = 1
								
    count_array.append(NBX)
    count_array.append(NBY)
    count_array.append(NBZ)
    count_array.append(NAT)
    count_array.append(minx)
    count_array.append(miny)
    count_array.append(minz)
	# count_array contains: NBX, NBY, NBZ, NAT, minx, miny, minz

    # Create a "z-fast" count array of points for the number of points required for the grid.

    zfast = [[['-100' for nz in range(0, int(count_array[2]))] for ny in range(0, int(count_array[1]))]	for nx in range(0, int(count_array[0]))]

    valx = 0
    valy = 0
    valz = 0
    originx=0
    originy=0
    originz=0
    read_bufferpip = 0

    # Normalize the dipole values to be consistent with grid values (i.e. no coords less than 0)    
    with open(pipoutfile, 'r') as pof:
		for item in pof:
			if (read_bufferpip == 1):
				if (int(count_array[4]) < 0):
					valx = (int(item.split()[1])) - int(count_array[4])
					originx=int(count_array[4])
				else:
					valx = int(item.split()[1])
				if (int(count_array[5]) < 0):
					valy = (int(item.split()[2])) - int(count_array[5])
					originy=int(count_array[5])
				else:
					valy = int(item.split()[2])
				if (int(count_array[6]) < 0):
					valz = (int(item.split()[3])) - int(count_array[6])
					originz=int(count_array[6])
				else:
					valz = int(item.split()[3])

				try:
					zfast[valx][valy][valz] = '100'
				except IndexError:
					pass
					
			if (re.search('ICOMP',item)):
				read_bufferpip = 1

    dx_body = '\n'
    counter = 0
    for value in zfast:
		for v2 in value:
			for v3 in v2:				
				dx_body += '{0} '.format(v3)
				counter = counter + 1
				if counter == 3:					
					dx_body += '\n'
					counter = 0

	# write DX header - gridnum, gridorigin, griddelta, gridpoints
	# gridnum - x y z lengths for shape
	# gridorigin - origin to apply shape at, should probably always be 0,0,0
	# gridpoints - number of points in shape (dipoles, basically). Requires a box.
	# griddelta - scaling factor between grid points, setting to always be 4 Angstrom for now.

    dx_header= """
object 1 class gridpositions counts {grid_numx} {grid_numy} {grid_numz}
origin {gridorigin}
delta {griddelta} 0.0 0.0
delta 0.0 {griddelta} 0.0
delta 0.0 0.0 {griddelta}
object 2 class gridconnections counts {grid_numx} {grid_numy} {grid_numz}
object 3 class array type double rank 0 items {gridpoints} data follows
			   """
    gridnumx = int(count_array[0])
    gridnumy = int(count_array[1])
    gridnumz = int(count_array[2])
    gridnumpts = int(count_array[3])
    grid_fill = gridnumx * gridnumy * gridnumz

    dx_head = dx_header.format(
    grid_numx='{0}'.format(gridnumx),
    grid_numy='{0}'.format(gridnumy),
    grid_numz='{0}'.format(gridnumz),
    gridorigin='{0} {0} {0}'.format(originx,originy,originz) ,
    griddelta='{0}'.format(dx_scale),
    gridpoints='{0}'.format(grid_fill)
    )

    dx_foot = """
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
			  """
			
    driver.put('output.string(dxgen).about.label'.format(n), 'Data Explorer Filetype (.dx)')
    driver.put('output.string(dxgen).current', (dx_head + dx_body + dx_foot))

def convertShapeToVTK(filename,outputname):
	
	# Read the data from the Shapefile
	read_now = 0
	lookuptable = {}
	with open(filename,'r') as sf:
		for line in sf:

			if read_now == 1:
				rx = int(line.split()[1])
				ry = int(line.split()[2])
				rz = int(line.split()[3])
				lookuptable['{0},{1},{2}'.format(rx,ry,rz)] = "1 "
			
				
			if (re.search('NBX',line)):
				fdx = int(line.split()[-3])
				fdy = int(line.split()[-2])
				fdz = int(line.split()[-1])

			if (re.search('IX',line)):
				read_now = 1
			

	if (int(fdx/2)*2 < fdx):
		fdxcl = -int(fdx/2)
		fdxcu = int(fdx/2)
	else:
		fdxcl = -int(fdx/2)+1
		fdxcu = int(fdx/2)				
	if (int(fdy/2)*2 < fdy):
		fdycl = -int(fdy/2)
		fdycu = int(fdy/2)
	else:
		fdycl = -int(fdy/2)+1
		fdycu = int(fdy/2)	
	if (int(fdz/2)*2 < fdz):
		fdzcl = -int(fdz/2)
		fdzcu = int(fdz/2)
	else:
		fdzcl = -int(fdz/2)+1
		fdzcu = int(fdz/2)	
	
	# Then write the VTK
	with open(outputname,'w') as vtkf:
		vtkf.write("# vtk DataFile Version 3.0\n3D Rendering\nASCII\nDATASET RECTILINEAR_GRID\n")
		vtkf.write("DIMENSIONS {0} {1} {2}\n".format(fdx,fdy,fdz))
		vtkf.write("X_COORDINATES {0} float\n".format(fdx))
		for item in range(fdxcl,fdxcu+1):
			vtkf.write("{0} ".format(item))
		vtkf.write("\n")
		vtkf.write("Y_COORDINATES {0} float\n".format(fdy))
		for item in range(fdycl,fdycu+1):
			vtkf.write("{0} ".format(item))
		vtkf.write("\n")
		vtkf.write("Z_COORDINATES {0} float\n".format(fdz))
		for item in range(fdzcl,fdzcu+1):
			vtkf.write("{0} ".format(item))
		vtkf.write("\n")
		vtkf.write("POINT_DATA {0}\n".format(int(float(fdx)*float(fdy)*float(fdz))))
		vtkf.write("SCALARS component double 1\n")
		vtkf.write("LOOKUP_TABLE default\n")
		counter = 0
		for ez in range(fdzcl,fdzcu+1):
			for ey in range(fdycl,fdycu+1):
				for ex in range(fdxcl,fdxcu+1):
					counter += 1
					try:
						vtkf.write("{0}".format(lookuptable['{0},{1},{2}'.format(ex,ey,ez)]))
					except KeyError:
						vtkf.write("0 ")

					if counter == 9:
						counter = 0
						vtkf.write("\n")

def make_fullshape(pip_path,max_dipoles,input_filename,output_filename, remote_choice):

	if (remote_choice == '1'):
		returncode,pipStdout,pipStderr = RapptureExec([pip_path, max_dipoles, input_filename, output_filename+'_0'],streamOutput=False)
		if returncode != 0:		
			text = pipStdout + pipStderr
			log(text)
			sys.stderr.write(text)
			driver.result(returncode)
#			sys.exit(returncode)
			if ('{0}'.format(returncode) != '0'):
				log(returncode)
				log("PIP failed to generate the shape")
		log(pipStdout + pipStderr)
	elif (remote_choice == '3'):
		exit_code = send_remote(max_dipoles, input_filename, output_filename, 0, pip_path)
		if ('{0}'.format(exit_code) != '0'):
			driver.result(exit_code)
			log(exit_code)
			log("The remote submission failed")
	return 1
	
def devolve_glyphs(object_num, save_object_num):	
	# Rappture data-read
    # Prepare extra output plots of the combined image if only the combined image is needed.
    if save_object_num > object_num:
		for n in range (2, save_object_num+1):			
			shutil.copyfile(os.path.join(temp_folder,'vtkout1'), os.path.join(temp_folder,'vtkout{0}'.format(n)))
		object_num = save_object_num

    for n in range(1, object_num+2):
		mesh_data = [0,'.']
		field_data = [0,'.']
		read_buffer1 = 0
		read_buffer2 = 0
		with open('vtkout{0}'.format(n),'r') as read_file:
			read_buffx = 0
			read_buffy = 0
			read_buffz = 0
			coord_datax = []
			coord_datay = []
			coord_dataz = []
			for line in read_file:
				a = line

				if (re.search('Y_COORDINATES',a)):
					read_buffx = 0
				if (re.search('Z_COORDINATES',a)):
					read_buffy = 0
				if (re.search('POINT_DATA',a)):
					read_buffz = 0
				
				if (read_buffx == 1):
					for item in line.split():						
						coord_datax.append(item)
				if (read_buffy == 1):
					for item in line.split():						
						coord_datay.append(item)
				if (read_buffz == 1):
					for item in line.split():						
						coord_dataz.append(item)
				
				if (re.search('DATASET',a)):
					read_buffer1 = 1
				if (re.search('POINT_DATA',a)):
					read_buffer1 = 0			
					read_buffz = 0
				if (re.search('LOOKUP_TABLE',a)):
					read_buffer2 = 1
				if (re.search('X_COORDINATES',a)):
					read_buffx = 1
				if (re.search('Y_COORDINATES',a)):
					read_buffx = 0
					read_buffy = 1
				if (re.search('Z_COORDINATES',a)):
					read_buffy = 0
					read_buffz = 1

												
				if (read_buffer1 == 1):
					mesh_data.append(a)
				if (read_buffer2 == 1):
					if (not(re.search('LOOKUP_TABLE',a))):
						for item in a.split():
							field_data.append(item)

# X-fast
# Combine datas
		full_data = []
		indexfdata = 2
		for zval in coord_dataz:
			for yval in coord_datay:
				for xval in coord_datax:
					full_data.append((xval, yval, zval, field_data[indexfdata]))
#					full_data.append((xval, yval, zval, saved_fdata[indexfdata]))
					indexfdata = indexfdata + 1
					
		with open('indexed_data{0}'.format(n), 'w') as indata:
#			indata.write('\n xlen {0}, ylen{1}, zlength{2}\n'.format(len(coord_datax), len(coord_datay), len(coord_dataz)))
			for x, y, z, entry in full_data:
				if '{0}'.format(entry) == '1':				
					entry = '{0}'.format(n)
				indata.write('{0},{1},{2},{3}\n'.format(x,y,z,entry))
					
# Note that the last data set read is the full one, which is nice.
    total_dex = []
    for n in range(1, object_num+1):
		rewrite_dex = []
		with open('indexed_data{0}'.format(n), 'r') as dexed_data:
			current_dexline = dexed_data.readline()
			curx, cury, curz, curent = current_dexline.split(',')
			for zval in coord_dataz:
				for yval in coord_datay:
					for xval in coord_datax:
						if (curx == xval) and (cury == yval) and (curz == zval):
							rewrite_dex.append(curent)
							current_dexline = dexed_data.readline()
							if (current_dexline != ''):								
								curx, cury, curz, curent = current_dexline.split(',')
						else:
							rewrite_dex.append('0')
		total_dex.append(rewrite_dex)


    compiled_dex = []
    if len(total_dex) > 1:
	    for n in range(1, len(total_dex)):
			endix = 0
			for j in range(0, len(total_dex[n-1])):
				if (int(total_dex[n-1][j]) == 0) and (int(total_dex[n][j]) == 0):
					compiled_dex.append('0')
				elif (int(total_dex[n-1][j]) != 0) and (int(total_dex[n][j]) == 0):
					compiled_dex.append('{0}'.format(int(total_dex[n-1][j])))						
				elif (int(total_dex[n-1][j]) == 0) and (int(total_dex[n][j]) != 0):
					compiled_dex.append('{0}'.format(int(total_dex[n][j])))						
				# Shouldn't ever reach here...but may be useful
				else:
					compiled_dex.append('1')
					

			total_dex[n] = compiled_dex
			compiled_dex = []

    last_dex = total_dex[-1]
    with open('lastdex','w') as lastone:
		breakcount = 0
		for entry in last_dex:
			lastone.write('{0} '.format(entry))
			breakcount = breakcount+1
			if breakcount == 9:
				lastone.write('\n')
				breakcount = 0
				


def print_glyphs(object_num, save_object_num, runstat):
		
	# Rappture data-read

    # Prepare extra output plots of the combined image if only the combined image is needed.
    if (save_object_num > object_num) and (runstat != 1):
		for n in range (2, save_object_num+1):			
			shutil.copyfile(os.path.join(temp_folder,'vtkout1'), os.path.join(temp_folder,'vtkout{0}'.format(n)))
    if (save_object_num > object_num):		
		object_num = save_object_num

    if runstat == 0:
		object_num = 0

    for n in range(object_num+1, object_num+2):
		mesh_data = [0,'.']
		field_data = [0,'.']
		read_buffer1 = 0
		read_buffer2 = 0
#		colorval_list = ['2578cb','4025cb','a225cb','cb2582','cb2925','cb8d25','9dcb25','32cb25']
#		normalized_n = n - 1
#		if n > 7:
#			normalized_n = normalized_n % 8
#		colorval = colorval_list[normalized_n]
			
		with open('vtkout{0}'.format(n),'r') as read_file:
			for line in read_file:
				a = line
				if (re.search('DATASET',a)):
					read_buffer1 = 1
				if (re.search('POINT_DATA',a)):
					read_buffer1 = 0			
				if (re.search('LOOKUP_TABLE',a)):
					read_buffer2 = 1
												
				if (read_buffer1 == 1):
					mesh_data.append(a)
				if (read_buffer2 == 1):
					if (not(re.search('LOOKUP_TABLE',a))):
						field_data.append(a)
	
		with open('field_data','w') as data1_file:
			if runstat==1:
				with open('lastdex','r') as dexful:
					data1_file.write(dexful.read())
			else:	
				for item in field_data:
					if ((item != 0) and (item !='.')):
						data1_file.write(item)
		with open('mesh_data','w') as data2_file:
			for item in mesh_data:
				if ((item != 0) and (item !='.')):
					data2_file.write(item)
		with open('field_data','r') as data1_file:
			driver.put('output.field(myfield).about.label', '3D Rendering')
			driver.put('output.field(myfield).about.view', 'glyphs')			
			driver.put('output.field(myfield).component.mesh', 'output.mesh(mymesh)')
			driver.put('output.field(myfield).component.association', 'pointdata')
			driver.put('output.field(myfield).component.elemtype', 'scalars')
			driver.put('output.field(myfield).component.elemsize', '1')
			driver.put('output.field(myfield).component.style', '-shape sphere -scaleMode scalar -gscale 0.1')
#			driver.put('output.field(myfield).component.style', '-color #{0} -shape sphere -scaleMode scalar -gscale 0.1'.format(colorval))
			driver.put('output.field(myfield).component.values', data1_file.read())
			
		with open('mesh_data','r') as data2_file:
			driver.put('output.mesh(mymesh).about.label', 'Object Mesh')
			driver.put('output.mesh(mymesh).dim', '3')
			driver.put('output.mesh(mymesh).units', 'um')
			driver.put('output.mesh(mymesh).hide', 'yes')
			driver.put('output.mesh(mymesh).vtk', data2_file.read())
	
def check_overlap(invert_choice, count_dex, object_count):
	# Scan the output file for overlapping data points.
    # Apply the points to either the largest or smallest object containing the points.
    
	# Split the write_dex into a dynamic set of lists for each object.
	# XYZ coords from each Dipole from object number 'n' is put into an entry at list of lists entry 'n-1'
    decks = [[] for i in range(1, object_count+1)]
    with open('write_dex', 'r') as write_dex:
		for element in write_dex:
			decks[int(element.split()[6]) - 1].append((element.split()[1], element.split()[2], element.split()[3]))

	# Find the overlapping values that are present between any two objects.
    found_dex = {}
    for n in range(1, object_count+1):
		for j in range(1, object_count+1):			
			if j != n:								
				size_result = cmp( int(count_dex[n-1][1]), int(count_dex[j-1][1]) )		
				found_list = list(set(decks[n-1]).intersection(decks[j-1]))
				
				for item in found_list:					
					val1= '-1'
					val2= '{0}'.format(item[0])
					val3= '{0}'.format(item[1])
					val4= '{0}'.format(item[2])
					if ((size_result <= 0) and (invert_choice == 'no')) or ((size_result > 0) and (invert_choice == 'yes')):
						try:
							found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,n,n,n)]
						except KeyError:
							try:
								found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,j,j,j)]
							except KeyError:
								found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,n,n,n)] = '{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,n,n,n)
								
								
					if ((size_result <= 0) and (invert_choice == 'yes')) or ((size_result > 0) and (invert_choice == 'no')):

						try:
							found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,j,j,j)]
						except KeyError:
							try:
								found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,n,n,n)]
							except KeyError:
								found_dex['{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,j,j,j)] = '{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}'.format(val1,val2,val3,val4,j,j,j)


    return found_dex


			


def progress():
    try:
        stage = progress_stages.pop(0)
    except IndexError:
        return
    Rappture.Utils.progress(*stage)


if __name__ == '__main__':
    progress_stages = [(0, "Initializing Converter..."),
                       (50, "Converting..."),
                       (70, "Rendering..."),
                       (90, "Loading output files..."),
                   ]
    progress()

    # Open driver
    driver_name = sys.argv[1]
    driver = Rappture.library(driver_name)
    if driver is None:
        print "Error opening file " + driver_name
        exit(-1)

    driver_number = Rappture.tools.getDriverNumber(driver_name)

    # Get the tool root directory
    tool_path = sys.argv[2]
    if tool_path == "":
        tool_path = os.path.dirname(driver.get('tool.version.application.directory(tool)'))

#	Pre-Processing File Cleanup in case of Aborted sessions

    for filename in ('field_data', 'mesh_data', 'vtkout', 'vtrout_1.vtr', 'estimate_out',
					 'faceless_input', 'dupe_output_filename', 'write_dex'):
        if os.path.exists(filename):
			os.remove(filename)
    if os.path.exists('submit_results'):
		shutil.rmtree('submit_results', ignore_errors = True)

    for n in range(1,10):
		vtkfile = ('vtkout{0}'.format(n))
		temp_shape_n = 'temp_shape{0}.obj'.format(n)
		reorder_shape_n = 'temp_reorder{0}.obj'.format(n)
		temp_out_shape_n = 'temp_out_shape{0}'.format(n)
		indexed_data_n = 'indexed_data{0}'.format(n)
		for filename in(vtkfile,temp_shape_n,temp_out_shape_n, reorder_shape_n, indexed_data_n):
			if os.path.exists(filename):
				os.remove(filename)

			
    # Various variable and path declarations.
    read_buffer1 = 0
    read_buffer2 = 0
    bin_path = os.path.join(tool_path, 'bin')
    pip_path = os.path.join(tool_path, 'bin/pip')
    resolution_choice = driver.get('input.group(res_group).choice(res_choice).current')
    nmblender = float(driver.get('input.number(nmblender).current'))

    render_choice = driver.get('input.boolean(display_choice).current')
    pop_choice = driver.get('input.choice(pop_choice).current')
    pop_val1 = driver.get('input.group(max_dip).integer(max_dipoles).current')
    pop_val2 = driver.get('input.group(dipCar_group).integer(dipCar).current')
    pop_val3 = driver.get('input.group(dipvol_group).integer(dipvol).current')
    remote_choice = driver.get('input.choice(RUNSTATE).current')    
    dx_choice = driver.get('input.boolean(dx_choice).current')
    dx_scale = float(driver.get('input.group(dxscale).number(dxscaleval).current'))
    input_filename = 'shape_{0}.obj'.format(driver_number)
    output_filename = 'shape_{0}.dat'.format(driver_number)
    temp_folder = os.getcwd()
    mcheck = 1


    # Write .obj file
    input_data = driver.get(driver.get('input.loader.upload.to') + '.current')
    blender_choice = driver.get(driver.get('input.loader.upload.to') + '.current')
    sessionnum = os.getcwd().split('/')[-1]		
    if input_data == "USE_BLENDER_FILE":
		for root, _, files in os.walk(os.getcwd()):
			for fil in files:
				if fil.startswith("saved_mesh") == True:
					your_shuttle = os.path.join(root,fil)
					with open(your_shuttle,'r') as blend_file:
						input_data = blend_file.read()
					driver.put((driver.get('input.loader.upload.to') + '.current'), input_data)
					driver.put('input.loader.current', 'Uploaded data')
			
    with open(input_filename, 'w') as obj_file:
        obj_file.write(input_data)


    # Scan .obj file for multiple shapes and if found split them so that they can be processed by multi-pip
        
    object_count = 0    
    object_lock1 = 0

	# Objects are defined by find an object delimiter, followed by finding it actually has vertices
    with open(input_filename,'r') as obj_file:
		for line in obj_file:
			if (re.search('v ',line)) and (object_lock1 == 1):
				object_count += 1
				object_lock1 = 0
			if not(re.search('v ',line)) and (object_lock1 == 1):
				object_lock1 = 0

			if (re.search('o ',line)):
				object_lock1 = 1


#    material_regex = re.compile(r'usemtl ')
#    prepare_regex = re.compile(r's ')
#    prepare_lock1_regex = re.compile(r'f ')

    object_offset = 0
    save_catch = 0
    object_index = []

    check_data = driver.get('input.loader.current')
    tool_error = 0
    if (object_count == 0) and ((check_data == 'Uploaded data') or (input_data == 'USE_BLENDER_FILE')):
		tool_error = 1
	
	# Rewrite files that have more than 1 object
    if (object_count > 1) and (tool_error != 1):	
		n = 0
		face_index = []
		with open(input_filename,'r') as obj_file:

			while True:
				if not line:
					break									
				if (re.search('o ',line)):
					n += 1				
					Rappture.Utils.progress(*(int(n), "Reading Shapes from File..."))
					object_index.append((n,line))
					temp_shape_n = 'temp_shape{0}.obj'.format(n)
					with open('{0}'.format(temp_shape_n), 'w') as temp_shape:
						temp_shape.write("# Blender v2.70 (sub 0) OBJ File: ''\n")
						temp_shape.write("# www.blender.org\n")
						temp_shape.write(line)
						line = obj_file.readline()
				else:
					line = obj_file.readline()
					if not line:
						break		
					else:
						continue
										
				if n > 0:
					# Grab the face_index array for each shape and save it in a list.
					# Separate each object into its own file.
					# Note that the "line" variable is still read from input_filename, not the temp_shape
					with open('{0}'.format(temp_shape_n), 'a') as temp_shape:					
						while (not re.search('o ',line)):
							if not line:
								break
							if (not (re.search('usemtl ',line))) and (not( re.search('f ',line))):
								temp_shape.write(line)
							if (re.search('usemtl ',line)):
								temp_shape.write('usemtl None\n')

							# Make sure that none of the lines which can have text with an "f " in it are used
							# i.e. material names, objects names, etc. ... everything else should be numbers anyways
							if (not (re.search('usemtl ',line))) and (not (re.search('s ',line)))\
								and (not (re.search('o ',line))) and (not (re.search('v ',line))) and ((re.search('f ',line))):
								for item in line.split():
									if item != 'f':
										item = int(item)  
										face_index.append((n,item))
									temp_shape.write('{0} '.format(item))
								temp_shape.write('\n')
							line = obj_file.readline()

		#set face_index with
		offset_index = [0]
		for n in range(1, object_count + 1):
			temp_list = []
			for element in face_index:
				if element[0] == n:
					temp_list.append(element[1])
			offset_index.append(min(temp_list)-1)


	# Make a second pass through the temporary files in order to place empty vertices for all unused shapes.
	# This will maintain the 'global perspective' of the dipole draw for each used shape.
		
	# Create a vertice-only file of text bites for all the shapes, removes face data
		with open(input_filename,'r') as obj_file:    
			for line in obj_file:
				if (not(re.search('f ',line))):
					with open('faceless_input', 'a') as faceless_file:
						faceless_file.write(line)
				
		n = 0
		for n in range(1, object_count + 1):
			not_this_one = 1
			temp_shape_n = 'temp_shape{0}.obj'.format(n)
			with open('{0}'.format(temp_shape_n), 'a') as temp_shape:
				with open('faceless_input','r') as faceless_file:				
					for line in faceless_file:						
						if (re.search('o ',line)):
							not_this_one = 0
						if (('{0}'.format(object_index[n-1][1])) == ('{0}'.format(line))):
							not_this_one = 1
						if (not_this_one == 0):
							temp_shape.write(line)
	# Make a third pass to set the face values for each object to enumerate only with respect to itself
		n = 0
		for n in range(1, object_count + 1):
			temp_shape_n = 'temp_shape{0}.obj'.format(n)
			reorder_shape_n = 'temp_reorder{0}.obj'.format(n)
			face_offset = offset_index[n]
			with open('{0}'.format(temp_shape_n), 'r') as temp_shape:
				with open('{0}'.format(reorder_shape_n), 'a') as reorder_shape:
					for line in temp_shape:
						if (re.search('f ',line)):
							for item in line.split():
								if item != 'f':
									item = int(item) - face_offset  
								reorder_shape.write('{0} '.format(item))

							reorder_shape.write('\n')
						elif not(re.search('f ',line)):
							reorder_shape.write(line)


    progress()




    # Get the defining characteristics of the system for the 3 population methods:
    # The length ratios, the max axis (and its value in dipoles),
    # the percentage of the system to be filled, the corresponding number of dipoles
    if (tool_error != 1):
		estPerc, estNAT, xyrat, yzrat, xzrat, maxdirect, maxdist = get_estimate(pip_path, input_filename)
		estPerc, estNAT, xyrat, yzrat, xzrat, maxdist = float(estPerc), int(estNAT), float(xyrat), float(yzrat), float(xzrat), float(maxdist)


		#  If the user chose a preset, use that value instead of a custom value    
		#  Method 1:
		if pop_choice == "1":
			max_dipoles = float(pop_val1)
			try:
				inex, inwhy, inzee = getdist_lightshuttle(input_filename, blender_choice)
				max_in_dist = max(inex, inwhy, inzee)
				if resolution_choice == "Low(1)":
					max_dipoles = 1
				elif resolution_choice == "Medium(2)":
					max_dipoles = 2
				elif resolution_choice == "High(3)":
					max_dipoles = 3
			
				max_dipoles = max_dipoles * max_in_dist * nmblender
				# Set the value for nm between 2 dipoles on the backend for the user.
				dp2nm = (max_in_dist*nmblender)/max_dipoles
				
			except IndexError:
				tool_error = 1
				driver.put('output.string(output_file).current', "To scale by Blender units, the shape must be exported through the Toolworkflow's Blender\n", append=True)




		# Modify max dipoles to suit the program since entry form is user-friendly.
		# Used in Method 2:
		elif pop_choice == "2":
			length_ratio = float(pop_val2) / float(float(estNAT) ** (1/float(3)))
			max_dipoles = float(pop_val2)
			# Set the value for nm between 2 dipoles on the backend for the user.
			try:
				inex, inwhy, inzee = getdist_lightshuttle(input_filename, blender_choice)
				max_in_dist = max(inex, inwhy, inzee)
				dp2nm = (max_in_dist*nmblender)/max_dipoles
			except IndexError:
				dp2nm = (length_ratio*nmblender)/max_dipoles

		# Using Method 3:
		elif pop_choice == "3":
			length_ratio = (float(pop_val3) / float(estNAT)) ** (1/float(3))
			max_dipoles = maxdist * length_ratio * nmblender
			# Set the value for nm between 2 dipoles on the backend for the user.
			try:
				inex, inwhy, inzee = getdist_lightshuttle(input_filename, blender_choice)
				max_in_dist = max(inex, inwhy, inzee)
				dp2nm = (max_in_dist*nmblender)/max_dipoles
			except IndexError:
				dp2nm = (length_ratio*nmblender)/max_dipoles

		# Take an initial estimate for disk space for predicted large problems
		if maxdirect == 'x':
			max_dipoles2 = float(max_dipoles) * (1 / float(xyrat))
			max_dipoles3 = float(max_dipoles) * (1 / float(xzrat))
		elif maxdirect == 'y':
			max_dipoles2 = float(max_dipoles) * xyrat
			max_dipoles3 = float(max_dipoles) * (1 / float(yzrat))
		elif maxdirect == 'z':
			max_dipoles2 = float(max_dipoles) * xzrat
			max_dipoles3 = float(max_dipoles) * yzrat

		max_dipoles  = '{0}'.format(int(round(max_dipoles)))
		max_dipoles2 = '{0}'.format(int(round(max_dipoles2)))
		max_dipoles3 = '{0}'.format(int(round(max_dipoles3)))

		if (float(max_dipoles) >= 99):
			trNAT = float(estPerc * (float(max_dipoles) * float(max_dipoles2) * float(max_dipoles3)))
			mcheck, sizeMB, freemem = memory_check(trNAT, render_choice)
			if (sizeMB > 1000):
				remote_choice = "3"
				
			if mcheck == 0:
				fail_message = "\n\nThe simulation/conversion requested\n requires more disk space than the user has available.\n The disk space required is {0}MB.\n The disk available is {1}MB.\n".format(sizeMB, freemem)
				tool_error = 1
				log(fail_message)
				render_choice = 'no'
				dx_choice = 'no'
		else:
			mcheck = 1

   
    # Run pip just once, i.e. for one object or for files not using the current Blender object format.
    separate_choice = driver.get('input.boolean(separate_choice).current')

    if remote_choice == '3':
		render_choice == 'no'

    save_object_count = object_count
    if ((object_count <= 1) or (separate_choice == 'no')) and (mcheck == 1) and (tool_error != 1):
		if (remote_choice == '1'):
			returncode,pipStdout,pipStderr = RapptureExec([pip_path, max_dipoles, input_filename, output_filename],streamOutput=False)
			if returncode != 0:
				text = pipStdout + pipStderr
				log(text)
				sys.stderr.write(text)
				driver.result(returncode)
#				sys.exit(returncode)
				if ('{0}'.format(returncode) != '0'):
					log(returncode)
					log("fullscale-PIP failed to generate the shape")
			log(pipStdout + pipStderr)
		elif (remote_choice == '3'):
			exit_code = send_remote(max_dipoles, input_filename, output_filename, 0, pip_path)
			if ('{0}'.format(exit_code) != '0'):
				driver.result(exit_code)
				log(exit_code)
				log("The remote submission failed")
#				sys.exit(exit_code)

	# Run multi-pip with multiple objects from Blender
    runstat = 0	
    if ((object_count > 1) and (separate_choice == 'yes') and (mcheck ==1) and (tool_error != 1)):
		runstat = make_fullshape(pip_path, max_dipoles, input_filename, output_filename, remote_choice)
		n = 0
		if (remote_choice == '1'):
			while n < object_count:			
				n += 1
				reorder_shape_n = 'temp_reorder{0}.obj'.format(n)
				temp_out_shape_n = 'temp_out_shape{0}'.format(n)
	
				returncode,pipStdout,pipStderr = RapptureExec([pip_path, max_dipoles, reorder_shape_n, temp_out_shape_n],streamOutput=False)
				if returncode != 0:				
					text = pipStdout + pipStderr
					log(text)
					sys.stderr.write(text)
					driver.result(returncode)
					if ('{0}'.format(returncode) != '0'):
						log(returncode)
						log("multiple-shape-PIP failed to generate the shape")					
#					sys.exit(returncode)
				log(pipStdout + pipStderr)
		elif (remote_choice == '3'):
			exit_code = send_remote(max_dipoles, input_filename, 'temp_out_shape', object_count, pip_path)
			if ('{0}'.format(exit_code) != '0'):
				driver.result(exit_code)
				log(exit_code)
				log("The remote submission failed.")
#				sys.exit(exit_code)
				
	# Read and Combine the separated temporary output files into a single file for DDSCAT
		count_dex = []
		running_NAT = 0
		val1 = 0

		with open('write_dex','w') as write_dex:		
			for n in range(1, object_count + 1):
				temp_out_shape_n = 'temp_out_shape{0}'.format(n)
				line_count = 0

				with open(temp_out_shape_n, 'r') as temp_out:
					for line in temp_out:
						line_count = line_count + 1
						if line_count == 1:
							firstline = line
						if (re.search('NAT',line)):
							running_NAT = running_NAT + int(line.split()[0])
							count_dex.append((n,int(line.split()[0])))
							NATline = '   {0} = NAT\n'.format(running_NAT)
						if line_count == 3:
							threeline = line
						if line_count == 4:
							fourline = line
						if line_count == 5:
							fiveline = line
						if line_count == 6:
							sixline = line
						if line_count == 7:
							sevenline = line
						if line_count >= 8:
							val1 = val1 + 1
							val2 = line.split()[1]
							val3 = line.split()[2]
							val4 = line.split()[3]
							write_dex.write('{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}\n'.format(val1,val2,val3,val4,n,n,n))


		with open('dupe_output_filename', 'a') as outfile:
			with open('write_dex','r') as write_dex:
				outfile.write(firstline+NATline+threeline+fourline+fiveline+sixline+sevenline)
				for element in write_dex:
					outfile.write(element)

				
		# Check for overlapping values, subtracting the duplicated values from the larger shape by default. Can be inverted.
		Rappture.Utils.progress(*(int(60), "Removing Overlapped Values..."))
		invert_choice = driver.get('input.group(inverted).boolean(invert_choice).current')		
		write_dex_overlap = check_overlap(invert_choice, count_dex, object_count)
		if write_dex_overlap == {}:
			os.rename(os.path.join(temp_folder,'dupe_output_filename'), os.path.join(temp_folder,output_filename))
		if write_dex_overlap != {}:
			# Recalc NAT based on incoming removals
			try:
				recalcNAT = int(NATline.split()[0]) - len(write_dex_overlap)
				NATline = '   {0} = NAT\n'.format(recalcNAT)		
			except ValueError:
				pass

			# Prepare a truncated list of overlaps as a library for upcoming comparison
			trunc_overlap = {}
			for entry in write_dex_overlap:
				a = "{0},{1},{2}".format(int(entry.split()[1]),int(entry.split()[2]),int(entry.split()[3]))
				trunc_overlap[a] = (int(entry.split()[1]),int(entry.split()[2]),int(entry.split()[3]))
				
			# Rewrite without duplicate values, with NAT updated
			with open(output_filename, 'a') as outfilea:
				outfilea.write(firstline+NATline+threeline+fourline+fiveline+sixline+sevenline)
				with open('dupe_output_filename','r') as dupe_file:
					i_count = 0
					val1 = 0
					for line in dupe_file:
						i_count = i_count + 1
						if i_count >= 8:
							line_match = 0
							XYZ_set = "{0},{1},{2}".format(int(line.split()[1]),int(line.split()[2]),int(line.split()[3]))
							try:
								line_match = trunc_overlap[XYZ_set]
								line_match = 1
							except KeyError:
								line_match = 0

							if line_match == 0:
								val1 = val1+1
								outfilea.write('{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}\n'.format(val1,int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4]),int(line.split()[5]),int(line.split()[6])))
					for element in write_dex_overlap:
						val1 = val1+1
						outfilea.write('{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}\n'.format(val1,int(element.split()[1]),int(element.split()[2]),int(element.split()[3]),int(element.split()[4]),int(element.split()[5]),int(element.split()[6])))
	
		# Rewrite the separated object input files that were used to create the master file, from the updated master file without duplicates.
			for n in range(1, object_count + 1):			
				with open(output_filename, 'r') as outfileb:	

					temp_out_shape_n = 'temp_out_shape{0}'.format(n)
					os.remove(temp_out_shape_n)
					line_count = 0
					reduced_NAT = 0
					for line in outfileb:
						line_count += 1
						if line_count >= 8:
							DDD_set = (int(line.split()[4]),int(line.split()[5]),int(line.split()[6]))
							nnn_set = (n,n,n)
							if DDD_set == nnn_set:
								reduced_NAT = reduced_NAT+1
				with open(output_filename, 'r') as outfile2:	
					with open(temp_out_shape_n, 'a') as temp_out:				
						line_count=0
						val_n = 0
						for line in outfile2:
							line_count += 1
							if (line_count < 8) and (line_count != 2):
								temp_out.write(line)
							if line_count == 2:
								temp_out.write('  {0} = NAT\n'.format(reduced_NAT))
							if line_count >= 8:
								DDD_set = (int(line.split()[4]),int(line.split()[5]),int(line.split()[6]))
								nnn_set = (n,n,n)
								if DDD_set == nnn_set:
									val_n = val_n+1
									temp_out.write('{0:>7}{1:>6}{2:>6}{3:>6}{4:>5}{5:>2}{6:>2}\n'.format(val_n,int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4]),int(line.split()[5]),int(line.split()[6])))								

		
		
    progress()
	# If 3D Rendering is asked for, do conversions necessary to render in Rappture and/or VTK:
	
    old_ofile = output_filename
    if (render_choice == 'yes' and (tool_error != 1) and (remote_choice != "3")):
	# Fix for cases where a non-standard, acceptable input file is used so that they still run.
	# Fix for case where multiple objects are found, but you only want the main image.
		if (object_count <= 1) or (separate_choice == 'no'):
			object_count = 1
			os.rename(os.path.join(temp_folder,output_filename), os.path.join(temp_folder,'temp_out_shape1'))
			output_filename = 'temp_out_shape1'

		if (runstat == 1):
			os.rename(os.path.join(temp_folder,old_ofile+'_0'), os.path.join(temp_folder,'temp_out_shape{0}'.format(object_count+1)))
			object_count = object_count + 1

		for n in range(1, object_count+1):
			infilname = os.path.join(temp_folder,'temp_out_shape{0}'.format(n))	
			outfilname = os.path.join(temp_folder, 'vtkout{0}'.format(n))
			convertShapeToVTK(infilname,outfilname)


		object_count = save_object_count
		if runstat == 1:			
			devolve_glyphs(object_count, save_object_count)
		print_glyphs(object_count, save_object_count, runstat)



    progress()


    # Load and generate output files into interface
    # The rehandle scales ok to large size problems.
    shuttle_status = driver.get('input.boolean(load_choice).current')
    if dx_choice == 'yes' and (tool_error != 1) and (remote_choice != "3"):
		Rappture.Utils.progress(*(int(95), "Generating DX File..."))
		dx_generate(output_filename, dx_scale)
		
    if (tool_error == 1):
		driver.put('output.string(output_file).current', "Invalid Shape File Attempted, no objects found in file.", append=True)
		
    if os.path.exists(output_filename):
		with open(output_filename, 'r') as shape_file:
			try:
				driver.put('output.string(output_file).current', shape_file.read())
			except NameError:
				pass
    if os.path.exists(output_filename):
		if shuttle_status == 'yes':
			with open(output_filename, 'r') as shape_file:
				try:
					# delete any old shuttles to DDSCAT
					for root, _, files in os.walk('/tmp/'):
						for fil in files:
							if fil.endswith(".elttuhs"+sessionnum) == True:
								os.remove(os.path.join(root,fil))
							if fil.endswith(".nmelttuhs"+sessionnum) == True:
								os.remove(os.path.join(root,fil))
					
					fddac, temp_path_ddac = mkstemp()
					os.rename(temp_path_ddac, temp_path_ddac+'.elttuhs'+sessionnum)
					with open(temp_path_ddac+'.elttuhs'+sessionnum,'w') as shuttle:
						shuttle.write(shape_file.read())
					os.close(fddac)
					
					nmddac, temp_path_ddac2 = mkstemp()
					os.rename(temp_path_ddac2, temp_path_ddac2+'.nmelttuhs'+sessionnum)
					with open(temp_path_ddac2+'.nmelttuhs'+sessionnum,'w') as shuttle:
						shuttle.write('NMfor2Dipoles: {0}'.format(dp2nm))
						shuttle.write('\nobject_count: {0}'.format(object_count))
					os.close(nmddac)
					
				except NameError:
					pass


    if (remote_choice == "3") and (tool_error != 1):
		driver.put('output.log(output_log)', "\n The converted file may be too large display in Rappture, but completed successfully.\n Additionally, the created file will still be passed along to DDSCAT \n", append=True)
		driver.put('output.log(output_log)', "\n Rendering this file in 3D is disallowed (file is too large) \n", append=True)
		driver.put('output.log(output_log)', "\n Creating a DX file for this file is disallowed (file is too large) \n", append=True)
#	Post-Processing File Cleanup

    for filename in ('field_data', 'mesh_data', 'vtkout', 'vtrout_1.vtr', 'lastdex',output_filename+'_0', 'estimate_out',
					'vtrout.pvd', 'a1a2.pvd','a1a2_1.vtr',output_filename,input_filename,'faceless_input', 'dupe_output_filename',
					'write_dex'):
        if os.path.exists(filename):
			os.remove(filename)
#			1
    if os.path.exists('submit_results'):
		shutil.rmtree('submit_results', ignore_errors = True)

    for n in range(1,object_count + 2):
		vtkfile = ('vtkout{0}'.format(n))
		vtrfile = ('vtrout{0}'.format(n))
		vtroutfile = ('vtrout{0}_1.vtr'.format(n))    
		vtroutpvd = ('vtrout{0}.pvd'.format(n))
		temp_shape_n = 'temp_shape{0}.obj'.format(n)
		reorder_shape_n = 'temp_reorder{0}.obj'.format(n)
		temp_out_shape_n = 'temp_out_shape{0}'.format(n)
		indexed_data_n = 'indexed_data{0}'.format(n)
		for filename in(vtkfile,vtrfile,vtroutfile,vtroutpvd,temp_shape_n,temp_out_shape_n, reorder_shape_n, indexed_data_n):
			if os.path.exists(filename):
#				1
				os.remove(filename)
		


		

#	Exit with results displayed
    driver.result(0)

