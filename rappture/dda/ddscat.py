import os
import Rappture
from Rappture.tools import executeCommand as RapptureExec
import shutil
import sys
import re
import stat
import tempfile
import string
import threading
import time
import math
import zipfile
import numpy as np
from math import cos, sin

param_file_template = """' ========= Parameter file for v7.3 ==================='
'**** Preliminaries ****'
'{CMDTRQ}' = CMTORQ*6 (DOTORQ, NOTORQ) -- either do or skip torque calculations
'{CMDSOL}' = CMDSOL*6 (PBCGS2, PBCGST, GPBICG, QMRCCG, PETRKP) -- CCG method
'{CMDFFT}' = CMETHD*6 (GPFAFT, FFTMKL) -- FFT method
'{CALPHA}' = CALPHA*6 (GKDLDR, LATTDR, FLTRCD) -- DDA method
'{CBINFLAG}' = CBINFLAG (NOTBIN, ORIBIN, ALLBIN) -- specify binary output
'**** Initial Memory Allocation ****'
{dipole_dim1} {dipole_dim2} {dipole_dim3} = dimensioning allowance for target generation
'**** Target Geometry and Composition ****'
'{CSHAPE}' = CSHAPE*9 shape directive
{SHPAR1} {SHPAR2} {SHPAR3} {SHPAR4} {SHPAR5} {SHPAR6} {SHPAR7}  = shape parameters 1 - 3
{NCOMP} {NCOMP1} {NCOMP2} {NCOMP3}       = NCOMP = number of dielectric materials
{materials1}
{materials2}
{materials3}
{materials4}
{materials5}
{materials6}
{materials7}
{materials8}
{materials9}
'**** Additional Nearfield calculation? ****'
{NRFLD} = NRFLD (=0 to skip nearfield calc., =1 to calculate nearfield E)
{NRFLD_r1} {NRFLD_r2} {NRFLD_r3} {NRFLD_r4} {NRFLD_r5} {NRFLD_r6} (fract. extens. of calc. vol. in -x,+x,-y,+y,-z,+z)
'**** Error Tolerance ****'
{TOL} = TOL = MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)/(NORM OF AC|E>)
'**** maximum number of iterations allowed ****'
{MXITER} = MXITER
'**** Interaction cutoff parameter for PBC calculations ****'
{GAMMA} = GAMMA (1e-2 is normal, 3e-3 for greater accuracy)
'**** Angular resolution for calculation of <cos>, etc. ****'
{ETASCA} = ETASCA (number of angles is proportional to [(3+x)/ETASCA]^2 )
'**** Vacuum wavelengths (micron) ****'
{WAVINI} {WAVEND} {NWAV} '{WCDIVID}' = wavelengths (first,last,how many,how=LIN,INV,LOG)
'**** Refractive index of ambient medium'
{NAMBIENT} = NAMBIENT
'**** Effective Radii (micron) **** '
{AEFFINI} {AEFFEND} {NRAD} '{RCDIVID}' = aeff (first,last,how many,how=LIN,INV,LOG)
'**** Define Incident Polarizations ****'
({X1},{X2}) ({Y1},{Y2}) ({Z1},{Z2}) = Polarization state e01 (k along x axis)
{IORTH} = IORTH  (=1 to do only pol. state e01; =2 to also do orth. pol. state)
'**** Specify which output files to write ****'
{IWRKSC} = IWRKSC (=0 to suppress, =1 to write ".sca" file for each target orient.
'**** Prescribe Target Rotations ****'
{BETA}    {BETA}   1  = BETAMI, BETAMX, NBETA  (beta=rotation around a1)
{THET}    {THET}   1  = THETMI, THETMX, NTHETA (theta=angle between a1 and k)
{PHI}    {PHI}   1  = PHIMIN, PHIMAX, NPHI (phi=rotation angle of a1 around k)
'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'
0 0 0 = first IWAV, first IRAD, first IORI (0 0 0 to begin fresh)
'**** Select Elements of S_ij Matrix to Print ****'
6 = NSMELTS = number of elements of S_ij to print (not more than 9)
11 12 21 22 31 41 = indices ij of elements to print
'**** Specify Scattered Directions ****'
{FRAME_TYPE} = CMDFRM (LFRAME, TFRAME for Lab Frame or Target Frame)
{NPLANES}{NPLANE_TEXT}
{PLANE1}
{PLANE2}
{PERIODIC_SCATTERING_ORDERS}


"""

def log(msg):
    try:
        driver.put('output.log(output_log)', msg, append=True)
    except NameError:
        pass



def find_all(name, path):
    """Find all the files with the specified name in path.

    Returns a generator of file names, ordered by wavelength value.

    """
    for root, _, files in os.walk(path):
        if name in files:
            yield os.path.join(root, name)


def check_lightshuttle():
    """Find the incident light point that was saved in Blender.

    Returns the rotation values that need to be implemented prior to any manually DDSCAT-set rotations.

    """

    sessionnum = os.getcwd().split('/')[-1]
    for root, _, files in os.walk(os.getcwd()):
        for fil in files:
            if fil.startswith("PolarLight") == True:
                your_shuttle = os.path.join(root,fil)
                with open(your_shuttle,'r') as blend_file:
                    input_data1 = blend_file.readline()
                    input_data2 = blend_file.readline()

                xLF = [float(input_data1.split()[1]),float(input_data1.split()[2]),float(input_data1.split()[3])]
                yLF = [float(input_data2.split()[1]),float(input_data2.split()[2]),float(input_data2.split()[3])]

    if 'xLF' not in locals():
        xLF = [1, 0, 0]
    if 'yLF' not in locals():
        yLF = [0, 1, 0]

	# could start using numpy arrays...
    g = np.array(xLF)
    h = np.array(yLF)
    z = np.cross(g,h)
    zLF = z.tolist()    
	# zLF = [ xLF[1]*yLF[2] - xLF[2]*yLF[1], xLF[2]*yLF[0] - xLF[0]*yLF[2], xLF[0]*yLF[1] - xLF[1]*yLF[0] ]

    a1 = [1, 0, 0]
    a2 = [0, 1, 0]
    a3 = [0, 0, 1]

	# normalize the vectors being used
    magx = math.sqrt(sum(xLF[i]*xLF[i] for i in range(len(xLF))))
    xLF = [ xLF[i]/magx for i in range(len(xLF)) ]
    magy = math.sqrt(sum(yLF[i]*yLF[i] for i in range(len(yLF))))
    yLF = [ yLF[i]/magy for i in range(len(yLF)) ]
    magz = math.sqrt(sum(zLF[i]*zLF[i] for i in range(len(zLF))))
    zLF = [ zLF[i]/magz for i in range(len(zLF)) ]

    dotX = sum(xLF[i]*a1[i] for i in range(len(xLF)))
    dotY = sum(yLF[i]*a2[i] for i in range(len(yLF)))
    dotZ = sum(zLF[i]*a3[i] for i in range(len(zLF)))
    
    brotX   = math.degrees(math.acos(dotX))
    brotY   = math.degrees(math.acos(dotY))
    brotZ   = math.degrees(math.acos(dotZ))




    return brotX,brotY,brotZ



def memory_check(NAT):

    """
    Checks if the generation of a shape of dipole size = NAT will cause the user's disk quota
    to be overdrawn.
    Input: NAT in number of dipoles
    Output: Pass/Fail Status, how much space was requested (MB), how much space there is (MB)

    """

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
    mem_to_use = (float(NAT) * 0.00115) + 1
    check_space = free_mem - mem_to_use
    
    if (check_space < 0):
		return 0, mem_to_use, free_mem
    elif (check_space >= 0):
		return 1, mem_to_use, free_mem

def memory_check_filesize(fileSize):
	
    """

    Checks if the generation of a fileSize (in bytes) will cause the user's disk quota
    to be overdrawn.
    Input: fileSize in bytes
    Output: Pass/Fail Status, how much space was requested (MB), how much space there is (MB)

    """

	# incoming fileSize is in bytes
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
    mem_to_use = (((fileSize)/1024)/1024) + 1
    check_space = free_mem - mem_to_use
    
    if (check_space < 0):
		return 0, mem_to_use, free_mem
    elif (check_space >= 0):
		return 1, mem_to_use, free_mem


def get_mem(NAT, Field_status):
	
    """

    Returns the memory amount estimated for the current settings.
    Based on Number of Dipoles, Nearfield Calculation (on/off).

    """
    memory_usage = 0

    venue_name = ''
    if Field_status != '0':
		memory_usage = (float(NAT)*0.016 + 42)
    elif Field_status == '0':
		memory_usage = (float(NAT)*0.002 + 42)
		
    if (memory_usage <= 16000):
		venue_name = 'rcac_S'
    elif (memory_usage <= 32000):
		venue_name = 'rcac_M'
    elif (memory_usage <= 48000):
		venue_name = 'rcac_L'
    elif (memory_usage <= 64000):
		venue_name = 'rcac_XL'
    elif (memory_usage <= 128000):
		venue_name = 'rcac_XXL'
    elif (memory_usage <= 192000):
		venue_name = 'rcac_XXXL'
    elif (memory_usage > 192000):
		venue_name = 'invalid_problem_size'

    return memory_usage, venue_name 
		
def data_NF235_LIST(par1, par2, par3):
    """
    
    Selects the correct value for DDSCAT's volume extension memory requirements.
    Essentially a recalculation of the X,Y,Z dimensions of the shape for calculation purposes.

    Returns the MXNX, MXNY, MXNZ values needed for the memory allocation.
    
    """
    NF235_list = [1,2,3,4,5,6,8,9,10,12,15,16,18,20,24,25,27,     
         30,32,36,40,45,48,50,54,60,64,72,75,80,81,90,96,100,     
         108,120,125,128,135,144,150,160,162,180,192,200,216,225, 
         240,243,250,256,270,288,300,320,324,360,375,384,400,405, 
         432,450,480,486,500,512,540,576,600,625,640,648,675,720, 
         729,750,768,800,810,864,900,960,972,1000,1024,1080,1125, 
         1152,1200,1215,1250,1280,1296,1350,1440,1458,1500,1536,  
         1600,1620,1728,1800,1875,1920,1944,2000,2025,2048,2160,  
         2187,2250,2304,2400,2430,2500,2560,2592,2700,2880,2916,  
         3000,3072,3125,3200,3240,3375,3456,3600,3645,3750,3840,  
         3888,4000,4050,4096]

    current_par1 = 0
    current_par2 = 0
    current_par3 = 0
         
    while(par1 != NF235_list[current_par1]):
		if par1 > NF235_list[current_par1]:
			current_par1 = current_par1 + 1
		elif par1 <= NF235_list[current_par1]:
			par1 = NF235_list[current_par1]
    while(par2 != NF235_list[current_par2]):
		if par2 > NF235_list[current_par2]:
			current_par2 = current_par2 + 1
		elif par2 <= NF235_list[current_par2]:
			par2 = NF235_list[current_par2]
    while(par3 != NF235_list[current_par3]):
		if par3 > NF235_list[current_par3]:
			current_par3 = current_par3 + 1
		elif par3 <= NF235_list[current_par3]:
			par3 = NF235_list[current_par3]

    return par1, par2, par3

def remove_all_w(num_jobs):
    """ 
    
    Find all the w0** files and remove them. Each wavelength can generate a new w0** file of each type.

    """	
    for n in range(0,num_jobs):
		for filename in ('w{0}r000k000.sca'.format(str(n).zfill(3)), \
		'w{0}r000.avg'.format(str(n).zfill(3)), 'w{0}r000k000.fml'.format(str(n).zfill(3)),\
		'w{0}r000k000.E1'.format(str(n).zfill(3)), 'w{0}r000k000.E2'.format(str(n).zfill(3)), \
		'w{0}r000k000.pol1'.format(str(n).zfill(3)), 'w{0}r000k000.pol2'.format(str(n).zfill(3)),\
		'w{0}r000k000.EB1'.format(str(n).zfill(3)), 'w{0}r000k000.EB2'.format(str(n).zfill(3))):
			if os.path.exists(filename):
				os.remove(filename)

def find_stderr(sign, path):
    """
    
    Find all the files with the specified tag .stderr in the path.
    Returns (respectively): 
    A formatted log of errors found,
    a concatenated list of outputs, 
    a formatted list of outputs.
    
    """
    output = ''
    stdout_list_log = ''
    stdout_list_cat = ''
    wave_regex = re.compile(r'1 wavelengths from\s+(.+) to\s+\1')
    wavelength_vals = []
    avoid_dupes = []
    avoid_dupes2 = []
    
    for root, _, files in os.walk(path):
		for fil in files:
			if fil.endswith(".stderr") == True:
				with open(os.path.join(root, fil), 'r') as output_file:
					for line in output_file:
						wave_match = wave_regex.search(line)
						if wave_match is not None:
							wavelength = float(wave_match.group(1))
							wavelength_vals.append((wavelength,os.path.join(root, fil)))

    wavelength_vals.sort(key=lambda x: x[0])
    
    
    for wavelength, filepather in wavelength_vals:
		with open(filepather, 'r') as output_filer:
			if filepather not in avoid_dupes:				
				output += output_filer.read()
				output += '\n +++ Next Output File +++ \n'
				avoid_dupes.append(filepather)
			
    for wavelength, filepather in wavelength_vals:
		replacedfilepath = '{0}'.format(filepather.replace('.stderr','.stdout'))
		with open(replacedfilepath, 'r') as output_filer2:
			if replacedfilepath not in avoid_dupes2:
				output_red = output_filer2.read()
				stdout_list_log += '\n Wavelength (um): {0} \n'.format(wavelength)
				stdout_list_log += output_red
				stdout_list_cat += output_red
				stdout_list_log += '\n +++ Next Error File +++ \n'
				avoid_dupes2.append(replacedfilepath)



    return stdout_list_log, stdout_list_cat, output



def find_regex(regex, path):
    """Find all the files that match the specified regex in path.

    Returns a generator of file names.

    """
    regex_c = re.compile(regex)
    for root, _, files in os.walk(path):
        for name in files:
            if regex_c.match(name) is not None:
                yield os.path.join(root, name)

def parse_table(filename):
    """Parse a DDSCAT output table.

    filename  The file name of the input table.

    Returns a tuple (header, data) where:
    header  A string containing the table header.
    data    A list of table rows as strings.

    """
    with open(filename, 'r') as table:
        table_lines = table.readlines()

        # Find the last line of the header.
        i = len(table_lines) - 1
        for line in reversed(table_lines):
            if line != '' and line[0] not in string.digits:
                i += 1
                break
            i -= 1

        header = ''.join(table_lines[:i])
        data = table_lines[i:]
        return header, data

def parse_output_log():
    """Parse the Output being sent to the log

    Returns an error message if one is needed/found.
    Returns a '0' if no error is found.

    """
    parse_this = driver.get('output.log(output_log)')
    list_this = parse_this.split('\n')
    error_msg_val = ''
    append_now = 0

    for line in list_this:
		if (append_now == 1):
			error_msg_val = error_msg_val + '\n' + line
		if (re.search('sigterm',line) or re.search('forrtl',line)):
			error_msg_val = error_msg_val + '\n' + line
			append_now = 1

    if (parse_this == ''):
		error_msg_val = 'No Output Log was Generated!'

    if (append_now == 1):
		return error_msg_val
    elif (append_now == 0):
		return '0'


def collate_table(output_name, partial_names, sort_slice):
    """Create a DDSCAT output table from a list of partial tables.

    output_name    The file name of the output table.
    partial_names  A list of table file names to merge.
    sort_slice     In a fixed width table, the column range of the sort key.

    Returns Wavelength with Max Light Extinction.

	Also, secondarily handles selecting a maximum E-Field from a list of E-Fields when applicable.
	This is because the best time to catch such handling occurs when collating a qtable.

    """
    with open(output_name, 'w') as table:
        header = None
        header2 = None
        data = []
        data2 = []
        Efield_path = "0"
        max_field_tuple = [(0,0)]
        for partial_name in partial_names:
            header, partial_data = parse_table(partial_name)
            data += partial_data
            data2.append((partial_data,partial_name))
        data.sort(key=lambda row: float(row[sort_slice]))

        table.write(header)
        for row in data:
            table.write(row)

    with open(output_name, 'r') as table:
        if (output_name == 'qtable'):
			b = []
			read_buffer = 0
			max_Qext_name = '1'
			for line in table:
				a = line
				if (read_buffer==1):
					b.append((a.split()[1],a.split()[2]))
				if (re.search('wave       Q_ext',a)):
					read_buffer=1

			if b == []:
				b.append('0 1')
			try:				
				max_Qext = max(b, key=lambda x:float(x[1]))
			except ValueError:
				max_Qext = ('0 1')

			for item,key in data2:
				try:
					if (re.search('{0}'.format(max_Qext[0]),('{0}'.format(item[0])).split()[1])):						
						max_Qext_name = '{0}'.format(key)
				except IndexError:
					1
			Efield_path = max_Qext_name

			working_path = os.path.join(os.getcwd(), 'w000r000k000.E1')
			working_pathB = os.path.join(os.getcwd(), 'w000r000k000.EB1')
			
			Efield_path  = os.path.join(max_Qext_name.split('/qtable')[0],'w000r000k000.E1')
			EBfield_path = os.path.join(max_Qext_name.split('/qtable')[0],'w000r000k000.EB1')
			
			if os.path.exists(Efield_path):
				os.rename(Efield_path,working_path)
			if os.path.exists(EBfield_path):
				os.rename(EBfield_path,working_pathB)
				
			return max_Qext[0]
			
def local_maxqext_grab(output_name, partial_names, sort_slice, bfield):
    """
    Returns Wavelength with Max Light Extinction for locally run simulations.

	Also, secondarily handles selecting a maximum E-Field from a list of E-Fields when applicable.
	This is because the best time to catch such handling occurs when parsing the qtable.

    """

    b = []
    read_buffer = 0
    max_Qext_name = '1'
	
    with open ('qtable','r') as table:				
		for line in table:
			a = line
			if (read_buffer==1):
				b.append((a.split()[1],a.split()[2]))
			if (re.search('wave       Q_ext',a)):
				read_buffer=1

		if b == []:
			b.append('0 1')
		max_Qext = max(b, key=lambda x:float(x[1]))
		
		count = 0	
		save_count = 0
		for item, key in b:				
			if item == max_Qext[0]:
				save_count = count
			count = count + 1
		if len('{0}'.format(save_count)) == 1:
			save_count = '00{0}'.format(save_count)
		elif len('{0}'.format(save_count)) == 2:
			save_count = '0{0}'.format(save_count)
		Efield_to_use = 'w000r000k000.E1'
		EBfield_to_use = 'w000r000k000.EB1'
	
    return Efield_to_use, EBfield_to_use, max_Qext[0]		


def BuildVTKfiles(RawDataFile, squareval, getSecret, gsX, gsY, gsZ):
    """Using the Raw data from DDSCAT Build the 
		data files for the          : E-field, log E-field, E-field Vectors, 
		and if requested in addition: B-field,              B-field Vectors, Poynting Vectors

    """

    read_x = 0
    read_y = 0
    read_z = 0

    min_x = 999999999
    min_y = 999999999
    min_z = 999999999
    min_e = 999999999
    max_x = -999999999
    max_y = -999999999
    max_z = -999999999
    max_e = -999999999
    num_pts_x = 0
    num_pts_y = 0
    num_pts_z = 0
    
    exc = {}
    eyc = {}
    ezc = {}

    eec = []
    eev = []
    bbc = []
    bbv = []
    
    pv = []
    secretdata = []

    
    with open(RawDataFile,'r') as getdata:
		for line in getdata:
						
			if re.search('Xcoord',line) or (line == "") or (line == "\n"):
				pass		
			elif re.search('Dimensions',line):
				num_pts_x = int(line.split()[-3])
				num_pts_y = int(line.split()[-2])
				num_pts_z = int(line.split()[-1])
			else:
				xval_in   = round(float(line.split()[0]),4)
				yval_in   = round(float(line.split()[1]),4)
				zval_in   = round(float(line.split()[2]),4)
				
				eval_in   = round(float(line.split()[3]),6)

				exRval_in = float(line.split('(')[1].split(',')[0])
				exIval_in = float(line.split('(')[1].split(',')[0].split(')')[0])
				eyRval_in = float(line.split('(')[2].split(',')[0])
				eyIval_in = float(line.split('(')[2].split(',')[0].split(')')[0])
				ezRval_in = float(line.split('(')[3].split(',')[0])
				ezIval_in = float(line.split('(')[3].split(',')[0].split(')')[0])

				bon_in    = float(line.split('(')[3].split(',')[1].split(')')[1].split()[0])
				bval_in   = float(line.split('(')[3].split(',')[1].split(')')[1].split()[1])
				
				bxRval_in = float(line.split('(')[4].split(',')[0])
				bxIval_in = float(line.split('(')[4].split(',')[0].split(')')[0])
				byRval_in = float(line.split('(')[5].split(',')[0])
				byIval_in = float(line.split('(')[5].split(',')[0].split(')')[0])
				bzRval_in = float(line.split('(')[6].split(',')[0])
				bzIval_in = float(line.split('(')[6].split(',')[0].split(')')[0])
				
				px_in     = float(line.split('(')[6].split(',')[1].split(')')[1].split()[0])
				py_in     = float(line.split('(')[6].split(',')[1].split(')')[1].split()[1])
				pz_in     = float(line.split('(')[6].split(',')[1].split(')')[1].split()[2])


				exc[xval_in] = xval_in
				eyc[yval_in] = yval_in
				ezc[zval_in] = zval_in			

				eec.append(eval_in)
				
				if (squareval == "2"):
					exRval_in = exRval_in**2
					eyRval_in = eyRval_in**2
					ezRval_in = ezRval_in**2
					bxRval_in = bxRval_in**2
					byRval_in = byRval_in**2
					bzRval_in = bzRval_in**2
					px_in = px_in**2
					py_in = py_in**2
					pz_in = pz_in**2

					
				eev.append((exRval_in,eyRval_in,ezRval_in))
				bbc.append(bval_in)
				bbv.append((bxRval_in,byRval_in,bzRval_in))
    
				pv.append((px_in,py_in,pz_in))
				
				if (xval_in < min_x):
					min_x = xval_in
				if (xval_in > max_x):
					max_x = xval_in

				if (yval_in < min_y):
					min_y = yval_in
				if (yval_in > max_y):
					max_y = yval_in
					
				if (zval_in < min_z):
					min_z = zval_in
				if (zval_in > max_z):
					max_z = zval_in

				if (eval_in < min_e):
					min_e = eval_in
				if (eval_in > max_e):
					max_e = eval_in

				if getSecret == "On":
					if gsX == "-1000":
						setX = 0
					else:
						setX = 1
					if gsY == "-1000":
						setY = 0
					else:
						setY = 1
					if gsZ == "-1000":
						setZ = 0
					else:
						setZ = 1
					
					if ((round(float(gsX),4) == round(xval_in,4)) and setX == 1):
						wantX = 1
					else:
						wantX = 0
					if ((round(float(gsY),4) == round(yval_in,4)) and setY == 1):
						wantY = 1
					else:
						wantY = 0
					if ((round(float(gsZ),4) == round(zval_in,4)) and setZ == 1):
						wantZ = 1
					else:
						wantZ = 0
						
					if (wantX == setX) and (wantY == setY) and (wantZ == setZ):
						secretdata.append(line)
    # First, we'll write some of that vector data to some basic text for Rappture
    # Could be re-structured if there are too many lines coming from a large simulation
    exc = sorted(exc)
    eyc = sorted(eyc)
    ezc = sorted(ezc)
    if getSecret == "On":
		with open('secret_data','w') as sd:
			sd.write('         Xcoord           Ycoord           Zcoord       EField           EField-X(Re)     EField-X(Im)       EField-Y(Re)     EField-Y(Im)       EField-Z(Re)     EField-Z(Im)    B On/Off       Bfield           BField-X(Re)     BField-X(Im)       BField-Y(Re)     BField-Y(Im)       BField-Z(Re)     BField-Z(Im)   Poynting X           Poynting Y       Poynting Z      \n')
			for item in secretdata:
				sd.write(item)
    with open('EField_Vec','w') as evecfile:
		for item in eev:
			evecfile.write('{0} {1} {2}\n'.format(item[0],item[1],item[2]))

    if float(bon_in) == 1:	
		with open('BField_Vec','w') as bvecfile:
			for item in bbv:
				bvecfile.write('{0} {1} {2}\n'.format(item[0],item[1],item[2]))
				
		with open('Poynting_Vec','w') as pvecfile:
			for item in pv:
				pvecfile.write('{0} {1} {2}\n'.format(item[0],item[1],item[2]))		



    # Then write the VTKs
    # At least everything is sorted already!
    fdx = num_pts_x
    fdy = num_pts_y
    fdz = num_pts_z
    lc = 0
    with open('headerVTK','w') as vtkf:
		vtkf.write("# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET RECTILINEAR_GRID\n")
		vtkf.write("DIMENSIONS {0} {1} {2}\n".format(fdx,fdy,fdz))
		vtkf.write("X_COORDINATES {0} float\n".format(fdx))
		for item in exc:
			vtkf.write("{0} ".format(item))
			lc += 1
			if lc == 9:
				vtkf.write("\n")
				lc = 0
		vtkf.write("\n")
		lc=0
		vtkf.write("Y_COORDINATES {0} float\n".format(fdy))
		for item in eyc:
			vtkf.write("{0} ".format(item))
			lc += 1
			if lc == 9:
				vtkf.write("\n")
				lc = 0
		vtkf.write("\n")
		lc=0
		vtkf.write("Z_COORDINATES {0} float\n".format(fdz))
		for item in ezc:
			vtkf.write("{0} ".format(item))
			lc += 1
			if lc == 9:
				vtkf.write("\n")
				lc = 0
		vtkf.write("\n")
		lc=0
				
		fdxyz = (int(float(fdx)*float(fdy)*float(fdz)))
		vtkf.write("POINT_DATA {0}\n".format(fdxyz))

    # Since we throw compositions into the regular E-field view as one...
    # We need to have a different ending to the header there (than B-field, log-field).		

    shutil.copyfile('headerVTK','EField_VTK.vtk')
    shutil.copyfile('headerVTK','EField_VTK_logscale.vtk')
    if float(bon_in) == 1:		
		shutil.copyfile('headerVTK','BField_VTK.vtk')
		with open('BField_VTK.vtk','a') as vtkf:
			vtkf.write('FIELD FieldData 1\n')
			vtkf.write('Intensity 1 {0} float\n'.format(fdxyz))

    with open('EField_VTK.vtk','a') as vtkf:
		vtkf.write('FIELD FieldData 2\n')
		vtkf.write('Intensity 1 {0} float\n'.format(fdxyz))
    with open('EField_VTK_logscale.vtk','a') as vtkf:
		vtkf.write('FIELD FieldData 1\n')
		vtkf.write('Intensity 1 {0} float\n'.format(fdxyz))
		
    
    with open('EField_VTK.vtk','a') as efile, open('EField_VTK_logscale.vtk','a') as elogfile:
		for item in eec:
			efile.write('{0} '.format(item))
			elogfile.write('{0} '.format(math.log(item)))
			lc += 1
			if lc == 9:
				efile.write("\n")
				elogfile.write("\n")
				lc = 0
		efile.write('\n')
		elogfile.write('\n')		
		lc=0
    # Note, skipping compositions for now! 
    # It may be helpful to be loading E-fields without composition data...
    # ...so that we can load larger files and faster!
    #
    # Otherwise something like this will follow, but you need CompositionData still:
		#vtkf.write('FIELD FieldData 1\n')
		efile.write('Composition 1 {0} float\n'.format(fdxyz))

		# Read some composition data in				
		with open('composition.txt','r') as compo:
			lc = 0
			CompositionData = compo.read()
			for item in CompositionData.split():
				try:
					efile.write('{0} '.format(int(float(item))))
				except ValueError:
					pass
				lc += 1
				if lc == 9:
					efile.write("\n")

					lc = 0

			efile.write('\n')

			lc=0		

				
    if os.path.exists('BField_VTK.vtk'):
		with open('BField_VTK.vtk','a') as bfile:			
			for item in bbc:
				bfile.write('{0} '.format(item))
				lc += 1
				if lc == 9:
					bfile.write("\n")
					lc = 0
			bfile.write('\n')
	# Again, could also potentially place some Composition data in the B-Field plot at this point.

    os.remove('headerVTK')

				
    return min_x,min_y,min_z,min_e,max_x,max_y,max_z,max_e,num_pts_x,num_pts_y,num_pts_z
			

def get_timing_info(output_string):
    """Get the total time taken to run a DDSCAT job from a string containing
    its stderr output. Also return the wavelength processed if the file
    only contains one.

    Returns (wavelength, total time).

    """
    total_time = 0.0
    wavelength = None
    time_regex = re.compile(r'(\d+\.\d{3}) (= CPU time)')
    time_regex2 = re.compile(r'(\d+\.) (= CPU time)')
    wave_regex = re.compile(r' >DDSCAT \s+(.+) = WAVE =')
    for line in output_string.split('\n'):
        time_match = time_regex.search(line)
        if time_match is not None:
            total_time += float(time_match.group(1))
        time_match2 = time_regex2.search(line)
        if time_match2 is not None:
            total_time += float(time_match2.group(1))

        wave_match = wave_regex.search(line)
        if wave_match is not None:
            wavelength = float(wave_match.group(1))
    return wavelength, total_time
    
def get_timing_info_local(output_string):
    """Get the total time taken to run a DDSCAT job from a string containing
    its stderr output. Also return the wavelength processed if the file
    only contains one.

    Returns (wavelength, total time) array.

    """
    total_time = 0.0
    wavelength = None
    save_wavelength = -1.0
    wavelength_times = []
    buffer1 = 0
    buffer2 = 0
    time_regex = re.compile(r'(\d+\.\d{3}) (= CPU time)')
    time_regex2 = re.compile(r'(\d+\.) (= CPU time)')
    wave_regex = re.compile(r' >DDSCAT \s+(.+) = WAVE =')
    for line in output_string.split('\n'):

        wave_match = wave_regex.search(line)
        if wave_match is not None:
            wavelength = float(wave_match.group(1))
            buffer2 = buffer2 + 1            
        time_match = time_regex.search(line)
        if time_match is not None:
            save_wavelength = wavelength
            total_time += float(time_match.group(1))
            buffer1 = buffer1 + 1           
        time_match2 = time_regex2.search(line)
        if time_match2 is not None:
            save_wavelength = wavelength
            total_time += float(time_match2.group(1))
            buffer1 = buffer1 + 1           
        if (buffer1 != 0) and (buffer2 == 2):
            buffer1 = 0
            buffer2 = 1
            wavelength_times.append((save_wavelength, total_time))
            total_time = 0

    wavelength_times.append((save_wavelength,total_time))

    return wavelength_times
    

def gen_custom_diels(diel_num_list):
    """ Given that custom constant dielectrics are expected, write and name the files appropriately 
	    Inputs: 
	    Tuple list of dieletric file number expected (i.e. 1,3,5,7 for diel1, diel3, diel5, diel7)
	      and list of values used to build the corresponding files.
	    
	    Generates:
	    Dielectric constant value files for every item in diel_num_list.
	    Gives same dielectric file names as uploaded dielectric files.
    """

    for item, value in diel_num_list:
		if (os.path.getsize('custom_dielectric{0}'.format(item)) == 0):
			with open('custom_dielectric{0}'.format(item),'w') as cust_diel:
				cust_diel.write("Internally Generated Constant Refractive Index for Dielectric{0}\n".format(item))
				cust_diel.write("1 2 3 0 0 0 = columns for wave, Re(n), Im(n), eps1, eps2\n")
				wave_list = [x/float(1000) for x in range (1,1001)]
				for wave_entry in wave_list:				
					cust_diel.write("{0} {1} 1.000E-6\n".format(wave_entry, value))

 
# Nanobio node
# Courtesy of Nahil Sobh, University of Illinois

def rotate3D(theta_X, theta_Y, theta_Z):
    """
    ========================================
    Rotates the Coordinates Axes given the:
    1- Rotation Around X denoted by Theta_X
    2- Rotation Around Y denoted by Theta_Y
    3- Rotation Around Z denoted by Theta_Z
    ========================================
    """
    theta_X = np.radians(theta_X)
    theta_Y = np.radians(theta_Y)
    theta_Z = np.radians(theta_Z)
    
    Sx = np.sin(theta_X)
    Cx = np.cos(theta_X)
    
    Sy = np.sin(theta_Y)
    Cy = np.cos(theta_Y)
    
    Sz = np.sin(theta_Z)
    Cz = np.cos(theta_Z)


    Rx = np.array( [ [   1,   0,   0] , 
                     [   0,  Cx, -Sx] , 
                     [   0,  Sx,  Cx] ], dtype=np.float)
                     
    Ry = np.array( [ [  Cy,   0,  Sy] , 
                     [   0,   1,   0] , 
                     [ -Sy,   0,  Cy] ], dtype=np.float)   
                     
    Rz = np.array( [ [  Cz, -Sz,   0] , 
                     [  Sz,  Cz,   0] , 
                     [   0,   0,   1] ], dtype=np.float)   
                  
    Rxyz = np.dot(Rz,np.dot(Ry,Rx))
    a1 = Rxyz[:,0]
    a2 = Rxyz[:,1]
    a3 = Rxyz[:,2]
    phi   = np.degrees(np.arctan2( Rxyz[2,0],Rxyz[1,0] ))
    beta  = np.degrees(np.arctan2((-1 * Rxyz[0,2]),Rxyz[0,1]))
    theta = np.degrees(np.arctan2(np.sqrt(np.square(Rxyz[0,1])+np.square(Rxyz[0,2])),Rxyz[0,0]))
	
    return phi, beta, theta

# Deprecated function!
def percent_timer(num_jobs, wall_time, sub_path):    
    """ Populates a list with the simulation job IDs and increments percentage bar based on their completion/walltime.
    	Increments percentage bar based on jobs that have initiated a 'Running' state or are 'Done'.
    	
    	Current Status: Disabled due to threading not being a good solution.
    	
    """
    queue_id = [0]
    job_id = [0]
    first_running_id = [0]
    done_id=[0]

    while_buffer = 0
    n = 0
    current_sum = 20
	
    queue_regex = re.compile(r'\((\d+)\) Simulation Queued at')
    job_regex = re.compile(r'\((\d+)\) Job Submitted at')
    running_regex = re.compile(r'\((\d+)\) Simulation Running at')
    done_regex = re.compile(r'\((\d+)\) Simulation Done at')
	
    while (sorted(done_id) != sorted(queue_id)) or (while_buffer == 0):
		with open(os.path.join(sub_path,'submit_log'), 'rw') as submit_stream:
			for line in submit_stream:
				queue_match = queue_regex.search(line)
				job_match = job_regex.search(line)
				running_match = running_regex.search(line)
				done_match = done_regex.search(line)
				if (queue_match) and (queue_match.group(1) not in queue_id):
					queue_id.append(queue_match.group(1))
					while_buffer = 1
				if (job_match) and (job_match.group(1) not in job_id):
					job_id.append(job_match.group(1))
					while_buffer = 1					
				if (running_match) and (running_match.group(1) not in first_running_id):
					n = n+1					
					first_running_id.append(running_match.group(1))
					current_sum = 10 + n*(40/num_jobs)
					Rappture.Utils.progress(*(int(current_sum), "Running DDSCAT..."))
		
				queue_id = list(set(queue_id)|set(job_id))
	
				if (done_match) and (done_match.group(1) not in done_id):
					n = n+1					
					done_id.append(done_match.group(1))
					current_sum = 10 + n*(40/num_jobs)
					Rappture.Utils.progress(*(int(current_sum), "Running DDSCAT..."))
					
		


				


def progress():
    try:
        stage = progress_stages.pop(0)
    except IndexError:
        return
    Rappture.Utils.progress(*stage)


# ==================== The main Program starts here ======================

if __name__ == "__main__":
    progress_stages = [(0, "Initializing DDSCAT..."),
                       (10, "Running DDSCAT..."),
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


	#   Pre-processing cleanup in case of an aborted simulation

    for filename in ('custom_dielectric1','custom_dielectric2','custom_dielectric3','custom_dielectric4',
					'custom_dielectric5','custom_dielectric6','custom_dielectric7','custom_dielectric8',
					'custom_dielectric9', 'mtable','qtable', 'qtable2', 'qtable_data',
					'Styles', 'field_datafile.txt',
					'EField_VTK.vtk','BField_VTK.vtk','EField_VTK_logscale.vtk',
					'EField_Vec','BField_Vec','Poynting_Vec',
					'AuDiel.tab','AgDiel.tab','Pt_diel.tab','Pd_diel.tab','Cu_diel.txt','TiO2','SiO2',
					'shape.dat', 'target.out', 'ddscat.par', 'composition.txt',
					'w000r000k000.sca','w000r000.avg','w000r000k000.fml','zipfilename','zipfilename_VTK'):
        if os.path.exists(filename):
			os.remove(filename)
    if os.path.exists('submit_results'):
		shutil.rmtree('submit_results', ignore_errors = True)
    if os.path.exists('submit_results_efield'):
		shutil.rmtree('submit_results_efield', ignore_errors = True)


    # Set the Failure flag to initially not failed.
    ddscat_fail_flag = '0'
    ddscat_fail_message = ''
    

    parameter_groups = (
        ('input.phase(page1).group(options).group(hide_CSHAPE).choice({0})', ('CSHAPE',)),
        ('input.phase(page1).group(options).group(NCOMP_SET).integer({0})', ('NCOMP',)),
        ('input.phase(page1).group(options).group(NCOMP_SET1).integer({0})', ('NCOMP1',)),
        ('input.phase(page1).group(options).group(NCOMP_SET2).integer({0})', ('NCOMP2',)),
        ('input.phase(page1).group(options).group(NCOMP_SET3).integer({0})', ('NCOMP3',)),
        ('input.phase(page1).group(options).group(customSHPARs).number({0})', ('customDDIST',)),
        ('input.phase(page1).group(options).group(SHPARs).number({0})', ('SHPAR1','SHPAR2','SHPAR3','SHPAR4','SHPAR5','SHPAR6','DDIST')),
        ('input.phase(page3).group(advanced).choice({0})', ('CMDTRQ', 'CMDSOL', 'CALPHA', 'GAMMA')),
        ('input.phase(page3).group(advanced).number({0})', ('ETASCA', 'TOL')),
        ('input.phase(page3).group(advanced).integer({0})', ('NPLANES',)),
        ('input.phase(page3).group(NRFLD_HEAD).choice({0})', ('NRFLD',)),
        ('input.phase(page5).group(process).integer({0})', ('MXITER',)),
        ('input.phase(page3).group(NRFLD_increase).number({0})',
         ('NRFLD_r1', 'NRFLD_r2', 'NRFLD_r3', 'NRFLD_r4', 'NRFLD_r5', 'NRFLD_r6')),
        ('input.phase(page2).group(Wavelengths).number({0})', ('WAVINI', 'WAVEND')),
        ('input.phase(page2).group(Wavelengths).integer({0})', ('NWAV',)),
        ('input.phase(page2).group(Wavelengths).choice({0})', ('WCDIVID',)),
        ('input.phase(page2).group(Wavelengths).string({0})', ('WAV_table',)),
        ('input.phase(page1).group(options).group(Ambient).number({0})', ('NAMBIENT',)),
        ('input.phase(page1).boolean({0})', ('IORTH',)),
        ('input.phase(page1).group(options).group(Polarization).group(X).number({0})', ('X1', 'X2')),
        ('input.phase(page1).group(options).group(Polarization).group(Y).number({0})', ('Y1', 'Y2')),
        ('input.phase(page1).group(options).group(Polarization).group(Z).number({0})', ('Z1', 'Z2')),
        ('input.phase(page1).group(options).group(Rotations).group(Beta).number({0})', ('BETA',)),
        ('input.phase(page1).group(options).group(Rotations).group(Theta).number({0})', ('THET',)),
        ('input.phase(page1).group(options).group(Rotations).group(Phi).number({0})', ('PHI',)),
    )

    params = {}
    for group, param_names in parameter_groups:
        for param_name in param_names:
            params[param_name] = driver.get(group.format(param_name)+'.current')

    params['CMDFFT'] = 'GPFAFT'   # Do not use the Intel MKL library
    params['CBINFLAG'] = 'NOTBIN' # Do not write binary files
    params['IWRKSC'] = '0'          # Write a .sca file for each target orientation
    
    period_type = driver.get('input.phase(page3).choice(PERIOD).current')
    params['SHPAR7']=''
    params['FRAME_TYPE'] = 'LFRAME'
    params['NPLANE_TEXT']=' = NPLANES = number of scattering planes'
    params['PERIODIC_SCATTERING_ORDERS'] = ''

    brotX,brotY,brotZ = check_lightshuttle()
    LightType = driver.get('input.phase(page1).group(options).group(Rotations).group(ILight).boolean(ILIGHT).current')
    if LightType == 'no':
		brotX,brotY,brotZ = 0,0,0
    initialrotX = float(params['PHI'])
    initialrotY = float(params['BETA'])
    initialrotZ = float(params['THET'])
    rotX = brotX + initialrotX
    rotY = brotY + initialrotY
    rotZ = brotZ + initialrotZ
    rotPhi,rotBeta,rotTheta = rotate3D(rotX,rotY,rotZ)
    params['PHI'] = '{0}'.format(rotPhi)
    params['BETA'] = '{0}'.format(rotBeta)
    params['THET'] = '{0}'.format(rotTheta)
    cudiel = {}


	# Perform data rewrite for Polarization types
	# Technically this could be done in the tool.xml via an example uploader (future implementation).
    PolType = driver.get('input.phase(page1).group(options).group(Polarization_set).choice(Polar_choice).current')
    if PolType == '1':
		params['Y1']= '1'
		params['Y2']= '0'
		params['Z1']= '0'
		params['Z2']= '0'
    elif PolType == '2':
		params['Y1']= '0'
		params['Y2']= '0'
		params['Z1']= '1'
		params['Z2']= '0'
    elif PolType == '3':
		params['Y1']= '1'
		params['Y2']= '0'
		params['Z1']= '0'
		params['Z2']= '1'
    elif PolType == '4':
		params['Y1']= '1'
		params['Y2']= '0'
		params['Z1']= '0'
		params['Z2']= '-1'
    elif PolType == '5':
		params['Y1']= '1'
		params['Y2']= '0'
		params['Z1']= '0'
		params['Z2']= '0'
		params['IORTH']= 'yes'

	# Perform data rewrite for custom dielectrics being used.
    for i in range(1,10):		
		check_diel = ''
		check_diel = driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel{0}).loader(compload{0}).current'.format(i))
		if check_diel == 'Uploaded data':
			cudiel[i] = '1'
		else:
			cudiel[i] = driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel{0}).choice(CDIEL{0}).current'.format(i))


	# Note that the current version of DDSCAT (7.3) forces the custom shape file to be named 'shape.dat'
	# Thus, a copy of the shuttle file has to be made in order to preserve unique shuttles.
    sessionnum = os.getcwd().split('/')[-1]		
    check_variable = driver.get('input.phase(page1).group(options).group(uploader).loader(loaded).current')
    cshape_check = driver.get('input.phase(page1).group(options).group(hide_CSHAPE).choice(CSHAPE).current')
    input_shape_filename = 'shape.dat'

	# If the input file is an uploaded file:	
    # Place any custom shape files in the same directory as the parameter file
    if (check_variable == 'Uploaded data'):
		cshape_check = '8'
		driver.put('input.phase(page1).group(options).group(hide_CSHAPE).choice(CSHAPE).current','8')
		input_shape_filename = 'shape.dat'
		input_shape_data = driver.get(driver.get('input.phase(page1).group(options).group(uploader).loader(loaded).upload.to') + '.current')
		driver.put("input.phase(page1).group(options).group(uploader).string(UploadedFile).current","empty",append=0)
		with open(input_shape_filename, 'w') as in_shape_file:
			in_shape_file.write(input_shape_data)
		in_shape_file = ""
		input_shape_data = ""
		with open(input_shape_filename,'r') as in_shape_file:
			read_dip_counter = 0
			dipole_NAT = 0
			ddim_array_x = []
			ddim_array_y = []
			ddim_array_z = []
			for line in in_shape_file:
				if read_dip_counter == 1:
					dip_array = line.split()
					try:
						ddim_array_x.append(int(float(dip_array[1])))
						ddim_array_y.append(int(float(dip_array[2])))
						ddim_array_z.append(int(float(dip_array[3])))
						dipole_NAT = dipole_NAT + 1
					except IndexError:
						pass
					except ValueError:
						pass
				if 'ICOMP(x,y,z)' in line.split():
					read_dip_counter = 1

		try:
			ddip1 = (max(ddim_array_x) - min(ddim_array_x)) + 1
			ddip2 = (max(ddim_array_y) - min(ddim_array_y)) + 1
			ddip3 = (max(ddim_array_z) - min(ddim_array_z)) + 1
			NAT2x, NAT2y, NAT2z = data_NF235_LIST(int(float(ddip1)),int(float(ddip2)),int(float(ddip3)))
		except (ValueError,IndexError):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNo valid Input Object File was found.\n"
		fileSize = os.path.getsize('shape.dat')
		mcheck1, sizeMB1, freemem1 = memory_check_filesize(fileSize)
		mcheck2, sizeMB2, freemem2 = memory_check(dipole_NAT)
		sizeMB = sizeMB1 + sizeMB2
		mcheck = mcheck1 + mcheck2
		if (mcheck == 0):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\n\nThe simulation/conversion requested\n requires more disk space than the user has available.\n The disk space required is approximately {0}MB.\n The disk available is {1}MB.\n".format(sizeMB, freemem1)



	# If the input file is a shuttle file:
    if (cshape_check == '9') and (check_variable != 'Uploaded data'):
		check_variable = 'Uploaded data'
		for root, _, files in os.walk('/tmp/'):
			for fil in files:
				if fil.endswith(".elttuhs"+sessionnum) == True:
					
					this_shuttle = os.path.join(root,fil)
					
					# Check that the user has enough space to make the copy
					fileSize = os.path.getsize(this_shuttle)
					with open(this_shuttle,'r') as getThatNAT:
						read_dip_counter = 0
						dipole_NAT = 0
						ddim_array_x = []
						ddim_array_y = []
						ddim_array_z = []
						for line in getThatNAT:
							if read_dip_counter == 1:
								dip_array = line.split()
								try:
									ddim_array_x.append(int(float(dip_array[1])))
									ddim_array_y.append(int(float(dip_array[2])))
									ddim_array_z.append(int(float(dip_array[3])))
									dipole_NAT = dipole_NAT + 1
								except IndexError:
									pass
								except ValueError:
									pass
							if 'ICOMP(x,y,z)' in line.split():
								read_dip_counter = 1
					try:
						ddip1 = (max(ddim_array_x) - min(ddim_array_x)) + 1
						ddip2 = (max(ddim_array_y) - min(ddim_array_y)) + 1
						ddip3 = (max(ddim_array_z) - min(ddim_array_z)) + 1
						NAT2x, NAT2y, NAT2z = data_NF235_LIST(int(float(ddip1)),int(float(ddip2)),int(float(ddip3)))
					except ValueError:
						ddscat_fail_flag = '1'
						ddscat_fail_message += "\nNo valid Input Object File was found.\n"
						break
					mcheck1, sizeMB1, freemem1 = memory_check_filesize(fileSize)
					mcheck2, sizeMB2, freemem2 = memory_check(dipole_NAT)
					sizeMB = sizeMB1 + sizeMB2
					mcheck = mcheck1 + mcheck2
					if (mcheck == 0):
						ddscat_fail_flag = '1'
						ddscat_fail_message += "\n\nThe simulation/conversion requested\n requires more disk space than the user has available.\n The disk space required is approximately {0}MB.\n The disk available is {1}MB.\n".format(sizeMB, freemem1)
						break

					with open(this_shuttle,'r') as shuttle_file:
						with open('shape.dat', 'w') as in_shape_file:
							check_space = shuttle_file.readline()
							if check_space != '':
								input_data = check_space
							input_data += shuttle_file.read()
							in_shape_file.write(input_data)
#							driver.put('input.phase(page1).group(options).group(uploader).string(UploadedFile).current', input_data)
#							driver.put('input.phase(page1).group(options).group(uploader).loader(loaded).current', 'Uploaded data')
							input_data = ""
					shuttle_file = ""
					in_shape_file = ""
				if fil.endswith(".nmelttuhs"+sessionnum) == True:
					this_shuttel = os.path.join(root,fil)
					with open(this_shuttel,'r') as ts:
						a = ts.readline()
						b = ts.readline()
					params['customDDIST'] = a.split()[-1]
					driver.put('input.phase(page1).group(options).group(customSHPARs).number(customDDIST).current',params['customDDIST'])
#					params['NCOMP'] = b.split()[-1]

							
    if check_variable == 'Uploaded data':
		params['CSHAPE'] = '8'
		count_ncomp = int(params['NCOMP'])
		for n in range(count_ncomp+1, 10):
			check_ncomp = cudiel[n]
			if check_ncomp != 'None':
				cudiel[n] = '5'
		for n in range(2,int(params['NCOMP'])+1):				
			check_ncomp = cudiel[n]
			if check_ncomp == '5':
				count_ncomp = count_ncomp - 1
		params['NCOMP'] = ('{0}'.format(count_ncomp))


	# Grab the NAT and x,y,z lengths from the shape file
    if os.path.exists('shape.dat') and (os.path.getsize('shape.dat') != 0):
		with open('shape.dat','r') as shapein:
			shline = shapein.readline()
			if shline == '\n':
				shline = shapein.readline()




	# Convert SHPAR values to dipoles from (nm)
    DipolesPerNM = float(params['DDIST'])
    params['SHPAR1'] = (float(params['SHPAR1'])*(DipolesPerNM))
    params['SHPAR2'] = (float(params['SHPAR2'])*(DipolesPerNM))
    if (params['CSHAPE'] != '4'):		
		params['SHPAR3'] = (float(params['SHPAR3'])*(DipolesPerNM))
    else:
		params['SHPAR3'] = (float(params['SHPAR3']))
    params['SHPAR4'] = (float(params['SHPAR4'])*(DipolesPerNM))
    params['SHPAR5'] = (float(params['SHPAR5'])*(DipolesPerNM))
    params['SHPAR6'] = (float(params['SHPAR6'])*(DipolesPerNM))



    # Set the plane values accordingly. 
    # Currently deprecated, default is set to always use 1 plane.
    if params['NPLANES'] == '1':
		params['PLANE1'] = '0. 0. 180. 1 = phi, thetan_min, thetan_max (deg) for plane A'
		params['PLANE2'] = ''
    if params ['NPLANES'] == '2':
		params['PLANE1'] = '0. 0. 180. 5 = phi, thetan_min, thetan_max (deg) for plane A'
		params['PLANE2'] = '90. 0. 180. 5 = phi, thetan_min, thetan_max (deg) for plane B'


    # Custom dimensioning is applied based on valid sizings given in the NF235 list.
    par1, par2, par3 = data_NF235_LIST(int(round(params['SHPAR1'])), int(round(params['SHPAR2'])), int(round(params['SHPAR3'])))
    params['dipole_dim1']= int(par1)
    params['dipole_dim2']= int(par2)
    params['dipole_dim3']= int(par3)

    # Adjustment for cylinder type memory requirements
    if (params['CSHAPE'] == '4') or (params['CSHAPE'] == '5') or (params['CSHAPE'] == '6'):
		par3 = par2
		if (params['CSHAPE'] == '5'):
			par1 = par1 + par2
		par1,par2,par3 = data_NF235_LIST(int(par1),int(par2),int(par3))
		params['dipole_dim1']= int(par1)
		params['dipole_dim2']= int(par2)
		params['dipole_dim3']= int(par3)		

    if (params['CSHAPE'] != '8') and (params['CSHAPE'] != '9'):
		dipole_NAT = int(par1) * int(par2) * int(par3)




	# Begin Shape Check routine to confirm correct parameters for respective Shape options
    none_check1 = driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel2).choice(CDIEL2).current')
    none_check2 = driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel3).choice(CDIEL3).current')

    if params['CSHAPE'] == '1':
		params['CSHAPE'] = 'ELLIPSOID'
		params['NCOMP'] = ''
		params['NCOMP2'] = ''
		params['NCOMP3'] = ''
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0) or (params['SHPAR3'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR4']="'shape.dat'"
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
    elif params['CSHAPE'] == '2':
		params['CSHAPE'] = 'ANIELLIPS'
		params['NCOMP'] = ''
		params['NCOMP1'] = ''
		params['NCOMP2'] = ''
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR4']="'shape.dat'"
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
		if none_check1 == '5':
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNot all required dielectric materials have been allocated\n"		
		if none_check2 == '5':
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNot all required dielectric materials have been allocated\n"
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0) or (params['SHPAR3'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
			
    elif params['CSHAPE'] == '3':
		params['CSHAPE'] = 'CONELLIPS'
		params['NCOMP'] = ''
		params['NCOMP1'] = ''
		params['NCOMP3'] = ''
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR7']="'shape.dat'"
		if none_check1 == '5':
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNot all required dielectric materials have been allocated\n"
		if (int(params['SHPAR1']) < int(params['SHPAR4'])) \
		 or (int(params['SHPAR2']) < int(params['SHPAR5'])) \
		 or (int(params['SHPAR3']) < int(params['SHPAR6'])):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nThe first concentric ellipsoid specified\n must have larger or equal parameters (SHPAR 1-3)\n compared to the second ellipsoid (SHPAR 4-6).\n"
			
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0) or (params['SHPAR3'] == 0.0) or (params['SHPAR4'] == 0.0) or (params['SHPAR5'] == 0.0) or (params['SHPAR6'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
    elif params['CSHAPE'] == '4':
		params['CSHAPE'] = 'CYLINDER1'
		params['NCOMP'] = ''
		params['NCOMP2'] = ''
		params['NCOMP3'] = ''
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
		if ((params['SHPAR3'] != 1.0) and (params['SHPAR3'] != 2.0) and (params['SHPAR3'] != 3.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nThe third parameter for Cylinders must have a value of 1 or 2 or 3.\n"		
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR4']="'shape.dat'"
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
    elif params['CSHAPE'] == '5':
		params['CSHAPE'] = 'CYLNDRCAP'
		params['NCOMP'] = ''
		params['NCOMP2'] = ''
		params['NCOMP3'] = ''
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR3']="'shape.dat'"
	#		params['SHPAR4']=''
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
    elif params['CSHAPE'] == '6':
		params['CSHAPE'] = 'UNIAXICYL'
		params['NCOMP'] = ''
		params['NCOMP1'] = ''
		params['NCOMP3'] = ''
	#	if (period_type == '1') or (period_type =='2'):
	#		params['SHPAR3']="'shape.dat'"
	#		params['SHPAR4']=''
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
		if none_check1 == '5':
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNot all required dielectric materials have been allocated\n"
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
    elif params['CSHAPE'] == '7':
		params['CSHAPE'] = 'RCTGLPRSM'
		params['NCOMP'] = ''
		params['NCOMP2'] = ''
		params['NCOMP3'] = ''
	#	if (period_type == '1') or (period_type =='2'):
	#		params['CSHAPE'] = 'RCTGL_PBC'
	#		params['SHPAR3']="'shape.dat'"
	#		params['SHPAR4']=''
	#		params['SHPAR5']=''
	#		params['SHPAR6']=''
		if ((params['SHPAR1'] == 0.0) or (params['SHPAR2'] == 0.0) or (params['SHPAR3'] == 0.0)):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nA parameter was specified with a value of 0.\n"		
    elif params['CSHAPE'] == '8':
		params['CSHAPE'] = 'FROM_FILE'
		params['NCOMP1'] = ''
		params['NCOMP2'] = ''
		params['NCOMP3'] = ''
		params['SHPAR1'] = ''
		params['SHPAR2'] = ''
		params['SHPAR3'] = ''
		if (period_type == '1') or (period_type =='2'):
			DipolesPerNM = 1/float(params['customDDIST'])
			params['CSHAPE']='FRMFILPBC'
			params['SHPAR1']= driver.get('input.phase(page3).group(PERIOD_SHPARs).number(PERIOD_SHPAR1).current')
			params['SHPAR2']= driver.get('input.phase(page3).group(PERIOD_SHPARs).number(PERIOD_SHPAR2).current')
			if driver.get('input.phase(page3).group(PERIOD_SHPARs).boolean(snap_PS1).current') == "yes":
				try:
					params['SHPAR1']= ddip2
				except NameError:
					pass
			if driver.get('input.phase(page3).group(PERIOD_SHPARs).boolean(snap_PS2).current') == "yes":
				try:
					params['SHPAR2']= ddip3
				except NameError:
					pass				
			params['SHPAR1'] = (float(params['SHPAR1'])*(DipolesPerNM))
			params['SHPAR2'] = (float(params['SHPAR2'])*(DipolesPerNM))
			params['SHPAR3']="'shape.dat'"
			params['SHPAR4']=''
			params['SHPAR5']=''
			params['SHPAR6']=''
		if (period_type == '1'):
			params['SHPAR2']='0'

	# Reallocate the dipole memory dimensions using NF235 values	
		try:
			ddip1 = NAT2x
			ddip2 = NAT2y
			ddip3 = NAT2z
			par1,par2,par3 = data_NF235_LIST(int(float(ddip1)),int(float(ddip2)),int(float(ddip3)))
			params['dipole_dim1'] = int(par1)
			params['dipole_dim2'] = int(par2)
			params['dipole_dim3'] = int(par3)
		except (ValueError,IndexError,NameError):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nNo valid Input Object File was found.\n"


    # Set correct numeric value for IORTH
    # IORTH is no longer always 1.
    if params['IORTH'] == 'yes':
        params['IORTH'] = '2'

    else:
        params['IORTH'] = '1'

    # If the user inputs custom wavelengths, write them to a file
    if params['WCDIVID'] == 'TAB':
        try:
            with open('wave.tab', 'w') as tabfile:
                tabfile.write(params['WAV_table'])
        except IOError, e:
            log('ERROR: ' + e.strerror)

    

    # Place any custom dielectric files in the same directory as the other dielectrics
    material_dir = os.path.join(tool_path, 'data/diel')


    # Repeated for Custom Dielectrics 2-9, which reside together in a different tool section then Diel #1
    indielname={}
    indieldata={}
    infile={}

    for i in range(1,10):
	    indielname['input_diel_filename{0}'.format(i)] = 'custom_dielectric{0}'.format(i)
	    indieldata['input_diel_data{0}'.format(i)] = driver.get(driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel{0}).loader(compload{0}).upload.to'.format(i)) + '.current')
	    with open(indielname['input_diel_filename{0}'.format(i)],'w') as infile['in_diel_file{0}'.format(i)]:
	    	infile['in_diel_file{0}'.format(i)].write(indieldata['input_diel_data{0}'.format(i)])

    # Generate any constant custom dielectric files needed.
    chek_en1 = ''
    chek_enn = ''
    diel_num_list = []
    for en in range (1,10):
		chek_enn = driver.get('input.phase(page1).group(options).group(dielectrics1to9).group(truediel{0}).loader(compload{0}).current'.format(en))
		en_value = driver.get('input.phase(page1).group(options).group(mydielectrics1to9).group(minidiel{0}).number(customm_CDIEL{0}).current'.format(en))
		if '{0}'.format(chek_enn)=="Input Constant Custom Dielectric":
			diel_num_list.append((en, en_value))
    
    gen_custom_diels(diel_num_list)

   # Prepare the file paths relevant to current dielectric file choice(s), stored in prep_material/prepmat[]
    prepmat={}
    prepmatval={}
    for i in range (1,10):
		prepmat['prep_material{0}'.format(i)] = os.path.join(tool_path, 'data/diel', cudiel[i])

   #Prepare Logic Values for identifying custom input dielectric files
   # Logic values are: 1 or 10 = custom diel, 5 = no diel, name = library diel file

    for i in range(1,10):
		prepmatval['prep_material{0}_val'.format(i)] = cudiel[i]


   # Modify material name if Custom Dielectric Material is selected

    for i in range(1,10):
    	if ((prepmatval['prep_material{0}_val'.format(i)] == '1') or (prepmatval['prep_material{0}_val'.format(i)] == '10')):
    		prepmat['prep_material{0}'.format(i)] = os.path.realpath(indielname['input_diel_filename{0}'.format(i)])

    # Place dielectric file paths in the paramater file

    for i in range(1,10):
    	if ((prepmatval['prep_material{0}_val'.format(i)] == '1') or (prepmatval['prep_material{0}_val'.format(i)] == '10')):
    		params['materials{0}'.format(i)] = """'{0}' = dielectric file {1}""".format(indielname['input_diel_filename{0}'.format(i)], i)
    	else:
    		params['materials{0}'.format(i)] = """'{0}' = dielectric file {1}""".format(prepmatval['prep_material{0}_val'.format(i)], i)

    # diel_files = dielectric files to be passed for processing, initialized with the first diel file.
    # diel_files must be sent the full path, .par file must be sent just the name.
    # Default case handling:
    # - Zeroed out input for when dielectrics are not in use.
    # - Pass the 'dielectric_' renames if a custom dielectric is used.
    # - Pass the default dielectric names if defaults are used.
    diel_files = [prepmat['prep_material1']]
    for i in range(2,10):
    	if prepmatval['prep_material{0}_val'.format(i)] == '5':
		params['materials{0}'.format(i)] = ''
	else:
		diel_files.append(prepmat['prep_material{0}'.format(i)])


    # Add any custom shapefiles to the list of files to send to processing.
    check_variable = driver.get('input.phase(page1).group(options).group(uploader).loader(loaded).current')
    cshape_check = driver.get('input.phase(page1).group(options).group(hide_CSHAPE).choice(CSHAPE).current')
    if cshape_check == '9':
		check_variable = 'Uploaded data'
    if check_variable == 'Uploaded data':
		diel_files.append(os.path.realpath(input_shape_filename))

    progress()

    # Perform last-stage value grabbing to determine which type of
    # DDSCAT to run and how many cores, threads to use.
    # ddscat_check is 0 for local, 1 for remote.

    run_type = driver.get('input.phase(page5).group(process).choice(RUNSTATE).current')
    num_cores = 1
    wall_time = int(driver.get('input.phase(page5).group(process).integer(WALLTIME_M).current'))
    collect_timing = "yes"

    
    diel_files_submit = []
    for diel_file in diel_files:
        diel_files_submit.append("-i")
        diel_files_submit.append(diel_file)


    # If Nearfield is to be calculated, select the file type to send to our lite version of ddpostprocess
    # Note: this is not configured for IORTH=2 at all.
    #       Thus, plots of .E2 and .EB2 are currently ignored.
    
    nearfield_calculate = driver.get('input.phase(page3).group(NRFLD_HEAD).choice(NRFLD).current')
    if (nearfield_calculate == '1'):
		ddppfile_to_use = 'w000r000k000.E1'
    if (nearfield_calculate == '2'):
		ddppfile_to_use = 'w000r000k000.EB1'
    if (nearfield_calculate == '1') or (nearfield_calculate == '2'):

	# Zero out any extended E-field boundaries for TUCs that are touching.
	# Note that for Periodicity, SHPAR1 holds the y-value. SHPAR2 holds the z-value.		
		params['NRFLD_r1'] = '0.5'
		params['NRFLD_r2'] = '0.5'
		params['NRFLD_r3'] = '0.5'
		params['NRFLD_r4'] = '0.5'
		params['NRFLD_r5'] = '0.5'
		params['NRFLD_r6'] = '0.5'
		if ((period_type == '1') or (period_type == '2')) and (params['SHPAR1'] != '') and (params['SHPAR2'] != ''):
			if os.path.exists('shape.dat') and (os.path.getsize('shape.dat') != 0):
				if ('{0}'.format(float(params['SHPAR1'])) == NAT2y):
					params['NRFLD_r3'] = '0'
					params['NRFLD_r4'] = '0'
				if ('{0}'.format(float(params['SHPAR2'])) == NAT2z) and (period_type == '2'):
					params['NRFLD_r5'] = '0'
					params['NRFLD_r6'] = '0'

    # If periodic conditions are used, build the Periodic Parameter Conditions
    if period_type == '1':		
		params['FRAME_TYPE'] = 'TFRAME'
		params['NPLANES'] = ''
		params['PLANE1'] = ''
		params['PLANE2'] = ''
		params['NPLANE_TEXT']= ''
		params['PERIODIC_SCATTERING_ORDERS']= '1 = number of scattering cones\n0.  0.  180.  0.05 = OrderM zetamin zetamax dzeta for scattering cone 1'
    if period_type == '2':		
		params['FRAME_TYPE'] = 'TFRAME'
		params['NPLANES'] = ''
		params['PLANE1'] = ''
		params['PLANE2'] = ''
		params['NPLANE_TEXT']= ''
		params['PERIODIC_SCATTERING_ORDERS']= '1 = number of scattering orders\n0.  0. = OrderM OrderN for scattered radiation'



	# If a custom shapefile is used, space its dipoles properly.
    if params['CSHAPE'] == 'FROM_FILE':
		DipolesPerNM = 1/float(params['customDDIST'])


	#Prepare the aeff values from DDIST for DDSCAT use.
    distance = (1/(DipolesPerNM))
    # convert back to microns
    distance = (distance/1000)
	# covert error handling
    if ddscat_fail_flag == '1':
		dipole_NAT = 1
    volume = float(((distance**3)*(dipole_NAT)))
    pie = math.pi
    aeff = ((3*volume)/(4*pie))**(float(1)/3)

  

    # Prepare the AEFFINI, AEFFEND, NRAD, RCDIVID values for the par file based on user input.    
    params['AEFFINI']='{0}'.format(aeff)
    params['AEFFEND']='{0}'.format(aeff)
    params['NRAD']='1'
    params['RCDIVID']='LIN'


	# Perform conversions for user-specified wavelength in (nm) to microns
    start_wav_temp = float(params['WAVINI'])
    params['WAVINI'] = (start_wav_temp/float(1000))
    end_wav_temp = float(params['WAVEND'])
    params['WAVEND'] = (end_wav_temp/float(1000))



    # If requested, run the command for submitting DDSCAT to a cluster,
    # splitting multiple wavelengths into different jobs.
    # Job specification is implicit to number of wavelength splits. 
    # Currently defaulted to 1 job per wavelength.
    # The single job with the max light extinction is re-submitted for nearfield simulations.
    # In the case of only a single wavelength, only a single nearfield simulation is run.
    
    Rappture.Utils.progress(*(int(30), "Running DDSCAT..."))
    
	# First, check the amount CPU/Memory requested is available:
    if (run_type == 'remote_splitting') and (ddscat_fail_flag == '0'):
        Field_status = driver.get('input.phase(page3).group(NRFLD_HEAD).choice(NRFLD).current')
        memory_usage, venue_name = get_mem(dipole_NAT, Field_status)
        if (Field_status != '0'):
			mcheck, sizeMB, freemem = memory_check(dipole_NAT)
			if (not mcheck):
				ddscat_fail_flag = '1'
				ddscat_fail_message += "\n\nThe simulation/conversion requested\n requires more disk space than the user has available.\n The disk space required is {0}MB.\n The disk available is {1}MB.\n".format(sizeMB, freemem)
				
        cores_to_use = int(math.ceil(float(memory_usage)/float(4000)))
        if cores_to_use == 0:
			cores_to_use = 1
        if (cores_to_use <= 48):
			Field_status = '0'
			memory_usage_noe, venue_name_noe = get_mem(dipole_NAT, Field_status)
			normal_cores_to_use = int(math.ceil(float(memory_usage_noe)/float(4000)))
			if normal_cores_to_use == 0:
				normal_cores_to_use = 1
        if (cores_to_use > 48):
			cores_to_use = 0
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\n The job requested was predicted to use {0}MB of RAM while the maximum allowed is 192000MB of RAM\n".format(memory_usage)

		# Note that the above error capture method still holds even though we're using tiers instead of core counts.
		# This is because the limit is still 192 GB, which was 48x4GB cores previously.

	# Second, actually prepare and send the submission if no fail flag is set
    if (run_type == 'remote_splitting') and (ddscat_fail_flag == '0'):
        single_length = 0
        nearfield_set = driver.get('input.phase(page3).group(NRFLD_HEAD).choice(NRFLD).current')
        start_wav = float(params['WAVINI'])
        end_wav = float(params['WAVEND'])
        if (start_wav == end_wav) and (nearfield_set != '0'):
			maxWaveExt = float(start_wav)
			single_length = 1
			
        num_wavs = int(params['NWAV'])

		# reset cores for current submit command method which requires that cores be input as 0 for all submissions.
		# counting cores is deprecated, but may be useful in the future if DDSCAT gets better parallel functionality.
        normal_cores_to_use = 0
        cores_to_use = 0

        assert start_wav <= end_wav
        assert num_wavs >= 1

        # Calculate specific wavelengths to use based on the selected
        # interpolation method.
        wavs = []
        if params['WCDIVID'] == 'LIN':
            if num_wavs != 1:				
				step = (end_wav - start_wav) / (num_wavs - 1)
            elif num_wavs == 1:
				step = 0
            current_wav = start_wav
            for _ in range(num_wavs):
                wavs.append(current_wav)
                current_wav += step
            wavs[-1] = end_wav

        # There should be at least one wavelength per job.
        num_jobs = num_wavs

        # Group wavelengths by job number.
        wav_groups = [list() for _ in range(num_jobs)]
        while len(wavs) > 0:
            for job_num in range(num_jobs):
                wav_groups[job_num].append(wavs.pop())

        home_dir = os.getcwd()

        params['WCDIVID'] = 'TAB'
        params['NRFLD'] = '0'

        with open('ddscat.par', 'w') as param_file:
            param_file.write(param_file_template.format(**params))


        # Output the .par file parameters to screen
        with open('ddscat.par','r') as parfile:
            for line in parfile:
				sys.stdout.write(line)

        # Create a temporary directory for each job number and write a
        # wave.tab file.
        wav_files = []
        working_dirs = []
        for job_num in range(num_jobs):
            working_dir = tempfile.mkdtemp(dir=home_dir)
            working_dirs.append(working_dir)
            os.chdir(working_dir)

            with open('wave.tab', 'w') as wav_table:
                # DDSCAT expects the wave.tab file to start with a header line.
                wav_table.write('\n')
                for wav in wav_groups[job_num]:
                    wav_table.write(str(wav) + '\n')

            wav_files.append(os.path.join(working_dir, 'wave.tab'))
            os.chdir(home_dir)

        command = ["submit"]
        command.append("-M")
        command.append("-v")
        command.append(venue_name_noe)
        command.append("-n")
        command.append(str(normal_cores_to_use))
        command.append("-e")
        command.append("OMP_NUM_THREADS=1")
        command.append("-e")
        command.append("DDSCAT_DISABLE_TARGET_OUT=TRUE")
        command.append("-w")
        command.append(str(wall_time))
        command.append("-p")
        command.append("@@wav=%s" % (','.join(wav_files)))
        command.append("-i")
        command.append("@@wav")
        command.append("-i")
        command.append(os.path.join(home_dir,'ddscat.par'))
        command += diel_files_submit
        if (params['CSHAPE'] == '8') or (params['CSHAPE'] == '9') :
           command.append("-i")
           command.append(os.path.join(home_dir,'shape.dat'))
        command.append("ddscat-7.3.0-intel-14_openmp")

        # Run submit , percentage bar incrementer run in parallel with 'submit' command.

        # Make a temporary directory to return the submit results to, save its path name.
        # Initialize an empty log for saving the stdout from the submit command.
        saved_home = os.getcwd()
        if not os.path.exists('submit_results'):
			os.makedirs('submit_results')
			sub_path = os.path.realpath('submit_results')
        os.chdir(sub_path)                
        with open(os.path.realpath('submit_log'), 'w') as submit_stream:			
			1

	# Deprecated timer threading:
	
	#      submitThread = threading.Thread(target=percent_timer, args=(num_jobs, wall_time, sub_path))
	#      submitThread.daemon=True
	#      submitThread.start()

        if single_length == 0:
			exit_status, stdout, stderr = RapptureExec(command, streamOutput=True)
			submit_log = stdout + stderr
			with open('submit_log','w') as slog:
				slog.write(submit_log)
			os.chdir(saved_home)			
			stdout_log, stdout_list, capture_out = find_stderr(".*.stderr$",home_dir)
        else:
			exit_status, stdout, stderr = '0','skip',''
			os.chdir(saved_home)				
			stdout_log, stdout_list, capture_out = 'skip stdoutlog\n','','skip capout\n'
			
        
        if ('{0}'.format(exit_status) != '0')  or (stdout_list != '') or (capture_out==''):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nThe standard remote submission returned unsuccessfully.\n"	
			if venue_name == 'invalid_problem_size':
				ddscat_fail_message += "\nInvalid Problem Size - The simulation requested more than 192GB of memory.\n"	



        # Collect job output from not-single-wavelength-Efield simulations.
        if (ddscat_fail_flag != '1') and (single_length != 1):
			collate_table('mtable', find_all('mtable', home_dir), slice(0, 10))
			maxWaveExt = collate_table('qtable', find_all('qtable', home_dir), slice(10, 21))
			collate_table('qtable2', find_all('qtable2', home_dir), slice(10, 21))

       
        # Remove temporary directories.
        for working_dir in working_dirs:
            shutil.rmtree(working_dir, ignore_errors=True)
        

        # re-submit the single wavelength which the Nearfield should be calculated for

        if (nearfield_set != '0') and (ddscat_fail_flag == '0'):
			Rappture.Utils.progress(*(int(90), "Running DDSCAT for Nearfield Wavelength..."))
			params['WAVINI'] = '{0}'.format(float(maxWaveExt))
			params['WAVEND'] = '{0}'.format(float(maxWaveExt))
			params['NWAV'] = '1'
			params['WCDIVID'] = 'LIN'
			params['NRFLD'] = '{0}'.format(nearfield_set)
			with open('ddscat.par', 'w') as param_file:				
				param_file.write(param_file_template.format(**params))

        # Create the submit command.

			command = ["submit"]
			command.append("-M")
			command.append("-v")
			command.append(venue_name)
			command.append("-n")
			command.append(str(cores_to_use))
			command.append("-e")
			command.append("OMP_NUM_THREADS=1")
			command.append("-e")
			command.append("DDSCAT_DISABLE_TARGET_OUT=TRUE")
			command.append("-w")
			command.append(str(wall_time))
			command.append("-i")
			command.append(os.path.join(home_dir,'ddscat.par'))
			command += diel_files_submit
			if (params['CSHAPE'] == '8') or (params['CSHAPE'] == '9') :
			   command.append("-i")
			   command.append(os.path.join(home_dir,'shape.dat'))
			command.append("ddscat-7.3.0-intel-14_openmp")

        # Make a temporary directory to return the submit results to, save its path name.
        # Initialize an empty log for saving the stdout from the submit command.
			saved_home = os.getcwd()        
			if not os.path.exists('submit_results_efield'):				
				os.makedirs('submit_results_efield')
				sub_path = os.path.realpath('submit_results_efield')

			working_path = os.path.join(os.getcwd(), 'w000r000k000.E1')
			working_pathB = os.path.join(os.getcwd(), 'w000r000k000.EB1')
			os.chdir(sub_path)                
			Efield_path = os.path.join(os.getcwd(), 'w000r000k000.E1')
			EBfield_path = os.path.join(os.getcwd(), 'w000r000k000.EB1')

			exit_status, stdout, stderr = RapptureExec(command, streamOutput=True)
			submit_log = stdout + stderr
			with open('submit_log','w') as slog:
				slog.write(submit_log)
			os.chdir(saved_home)				

			stdout_log, stdout_list, capture_out = find_stderr(".*.stderr$",home_dir)

			
			if ('{0}'.format(exit_status) != '0')  or (stdout_list != '') or (capture_out == ''):				
				ddscat_fail_flag = '1'
				ddscat_fail_message += "\nThe Nearfield-generating remote submission returned unsuccessfully.\n"		
				if venue_name == 'invalid_problem_size':
					ddscat_fail_message += "\nInvalid Problem Size - The simulation requested more than 192GB of memory.\n"	

			if os.path.exists(Efield_path):
				os.rename(Efield_path,working_path)
			if os.path.exists(EBfield_path):
				os.rename(EBfield_path,working_pathB)


        # Collect all .stdout and .stderr files and put them in the output log.
        stdout_log, stdout_list, capture_out = find_stderr(".*.stderr$",home_dir)

        log(capture_out)

        # Collect job output for single-wavelength-Efield simulations.
        if (ddscat_fail_flag != '1') and (single_length == 1):
			collate_table('mtable', find_all('mtable', home_dir), slice(0, 10))
			maxWaveExt = collate_table('qtable', find_all('qtable', home_dir), slice(10, 21))
			collate_table('qtable2', find_all('qtable2', home_dir), slice(10, 21))

		# Collect timing output
        if collect_timing == 'yes':
            total_time = 0.0
            wavelength_times = []
            concat_stderr_home = []
            for stderr_path in find_regex(r'.*\.stderr', home_dir):
                with open(stderr_path, 'r') as stderr_file:
                    wavelength, time = get_timing_info(stderr_file.read())
                    stderr_read = stderr_file.read()

                concat_stderr_home.append((wavelength, stderr_read))
                wavelength_times.append((wavelength, time))
                total_time += time

            # Sort by wavelength.
            concat_stderr_home.sort(key=lambda x: x[0])	
            wavelength_times.sort(key=lambda x: x[0])	

            catch_neartime = []
            timing_info = 'Total time: {} sec.\n'.format(total_time)
            for wavelength, time in wavelength_times:
				if wavelength in catch_neartime:
					timing_info += 'Wavelength (Nearfield) {}: {} sec.\n'.format(wavelength, time)
				if wavelength not in catch_neartime:
					catch_neartime.append(wavelength)
					timing_info += 'Wavelength {}: {} sec.\n'.format(wavelength, time)




    # Run command for regular DDSCAT.
    elif (run_type == 'local') and (ddscat_fail_flag == '0'):
        Field_status = driver.get('input.phase(page3).group(NRFLD_HEAD).choice(NRFLD).current')
        memory_usage, venue_name = get_mem(dipole_NAT, Field_status)
        if (Field_status != '0'):
			mcheck, sizeMB, freemem = memory_check(dipole_NAT)
			if (not mcheck):
				ddscat_fail_flag = '1'
				ddscat_fail_message += "\n\nThe simulation/conversion requested\n requires more disk space than the user has available.\n The disk space required is {0}MB.\n The disk available is {1}MB.\n".format(sizeMB, freemem)
			if (memory_usage > 16000):
				ddscat_fail_flag = '1'
				ddscat_fail_message += "\n\nThe simulation/conversion requested\n requires more memory than the user has available.\n The memory required is {0}MB.\n The memory available is 16000MB.\n".format(memory_usage)


        with open('ddscat.par', 'w') as param_file:
            param_file.write(param_file_template.format(**params))
            
        # Output the .par file parameters to screen

        with open('ddscat.par','r') as parfile:
            for line in parfile:
				sys.stdout.write(line)

        sys.stdout.flush()

        # Copy material files to the working directory.
        for diel_file in diel_files:
			if diel_file.split('/')[-1] in ('AuDiel.tab','AgDiel.tab','Pt_diel.tab','Pd_diel.tab','Cu_diel.txt','TiO2.tab','SiO2.tab'):				
				shutil.copy(diel_file, '.')


        if (params['CSHAPE'] == '8') or (params['CSHAPE'] == '9'):
			runDDSCAT = ['ddscat','DDSCAT_DISABLE_TARGET_OUT=TRUE']
        else:
			runDDSCAT = ['ddscat']

        start_wav_local = float(params['WAVINI'])
        end_wav_local = float(params['WAVEND'])
        nearfield_local = float(params['NRFLD'])
        if ((start_wav_local == end_wav_local) or (nearfield_local == 0)):
			exit_code, stdout, stderr = RapptureExec(runDDSCAT, streamOutput=True)
#
        elif ((nearfield_local != 0) and (start_wav_local != end_wav_local)):
			params['NRFLD'] = '0'
			os.rename('ddscat.par','ddscat.par_save')
			with open('ddscat.par', 'w') as param_file:
				param_file.write(param_file_template.format(**params))
			exit_code, stdout, stderr1 = RapptureExec(runDDSCAT, streamOutput=True)
#			log(stderr)
			Efield_to_use, EBfield_to_use, maxWaveExt = local_maxqext_grab('qtable', find_all('qtable', os.getcwd()), slice(10,21), nearfield_calculate)
			os.rename('qtable','qtable_save')
			os.rename('qtable2','qtable2_save')
			os.rename('mtable','mtable_save')
			if nearfield_calculate == "1":
				params['NRFLD'] = '1'
			elif nearfield_calculate == "2":
				params['NRFLD'] = '2'
			params['WAVINI'] = '{0}'.format(maxWaveExt)
			params['WAVEND'] = '{0}'.format(maxWaveExt)
			params['NWAV'] = '1'
			with open('ddscat.par', 'w') as param_file:
				param_file.write(param_file_template.format(**params))
			exit_code, stdout, stderr = RapptureExec(runDDSCAT, streamOutput=True)
			stderr = stderr1 + stderr
			os.rename('ddscat.par_save','ddscat.par')
			os.rename('qtable_save','qtable')
			os.rename('qtable2_save','qtable2')
			os.rename('mtable_save','mtable')


        if ('{0}'.format(exit_code) != '0') or (stdout != ''):
			ddscat_fail_flag = '1'
			ddscat_fail_message += "\nDDSCAT failed to exit successfully.\n"		

        stdout_log = stdout

        if collect_timing == 'yes':
            total_time = 0.0
            wavelength_times = get_timing_info_local(stderr)

            for wavelength, time in wavelength_times:				
				total_time += time

            catch_neartime = []
            timing_info = 'Total time: {} sec.\n'.format(total_time)
            for wavelength, time in wavelength_times:
				if wavelength in catch_neartime:
					timing_info += 'Wavelength (Nearfield) {}: {} sec.\n'.format(wavelength, time)
				if wavelength not in catch_neartime:
					catch_neartime.append(wavelength)
					timing_info += 'Wavelength {}: {} sec.\n'.format(wavelength, time)


        log(stdout)
        log(stderr)
        
        home_dir = os.getcwd()
        if ddscat_fail_flag != '1':			
			Efield_to_use, EBfield_to_use, maxWaveExt = local_maxqext_grab('qtable', find_all('qtable', home_dir), slice(10,21), nearfield_calculate)
        if not (os.path.exists('w000r000k000.E1')):
			Efield_to_use = 'w000r000k000.E1'			
        if not (os.path.exists('w000r000k000.EB1')):
			EBfield_to_use = 'w000r000k000.EB1'		


    progress()




	# Output an error message to the output logs if DDSCAT seems to have failed.
	# Note that this output is placed here so that it is the first thing the user sees in the logs if it occurs.
    if run_type != 'local':		
		Efield_to_use = 'w000r000k000.E1'
		EBfield_to_use = 'w000r000k000.EB1'
		
    check_outlog = parse_output_log()
    if check_outlog != '0':
		exit_code = 1
		ddscat_fail_flag = '1'
		stdout_log = '{0}'.format(check_outlog)
		ddscat_fail_message += 'Not enough time or memory was allocated to complete the job.'
		
    
    if ddscat_fail_flag == '1':

	# Quickly test some name allocations that don't get set if DDSCAT didn't run.
		try:			
			stdout_log += ''
		except NameError:
			stdout_log = ''
		try:			
			checkEfieldname = Efield_to_use
		except NameError:
			Efield_to_use = 'filenotfound'
		try:			
			checkEBfieldname = EBfield_to_use
		except NameError:
			EBfield_to_use = 'filenotfound'
		try:			
			checktiming = timing_info
		except NameError:			
			timing_info = 'No timing computed, DDSCAT failed.'
			
				
		driver.put('output.string(failure).about.label','***SIMULATION STATUS***')
		driver.put('output.string(failure).current','<<< DDSCAT has encountered an error in processing. >>>\n\n Resultant data displayed for the settings attempted is void. \n If it was generated, please see the Output Log for more detailed information\n')
		driver.put('output.string(failure).current','\nReason: ', append=True)
		driver.put('output.string(failure).current',ddscat_fail_message, append=True)
		driver.put('output.string(failure).current','\n\nCrash Log:\n', append=True)
		driver.put('output.string(failure).current','{0}'.format(stdout_log), append=True)
		driver.put('output.string(failure).current', '\n# Nearfield Not Requested or Not Usable', append=True)

	# Check that the file to be used for the E-field is not an empty file
    Efield_fail_flag = '0'

    if (nearfield_calculate == "0" or nearfield_calculate == "1"):
		if os.path.exists(Efield_to_use):
			Efsize = os.stat(Efield_to_use).st_size
			if Efsize <= 0:
					Efield_fail_flag = '1'

    if (nearfield_calculate == "2"):
		if os.path.exists(EBfield_to_use):
			EBfsize = os.stat(EBfield_to_use).st_size
			if EBfsize <= 0:
					Efield_fail_flag = '1'
				
			
    if Efield_fail_flag == '1':
		efield_fail_message = "\nThe request for an Electric Field did not return from DDSCAT successfully."
		driver.put('output.string(failure).about.label','***SIMULATION STATUS***')
		driver.put('output.string(failure).current','<<< DDSCAT has encountered an error in processing. >>>\n\n Resultant data displayed for the settings attempted is void. \n If it was generated, please see the Output Log for more detailed information\n')
		driver.put('output.string(failure).current',efield_fail_message, append=True)
		driver.put('output.string(failure).current', '\n# Nearfield Not Requested or Not Usable', append=True)


	# Write a successful output status if successful.
    if (Efield_fail_flag != '1') and (ddscat_fail_flag != '1'):		
		driver.put('output.string(failure).about.label','***SIMULATION STATUS***')
		driver.put('output.string(failure).current','<<< DDSCAT succeeded in processing! >>>\n If it was generated, please see the Output Log for more detailed information\n')
		driver.put('output.string(failure).current','\nWavelength Value Considered for Nearfield Calculations: {0}\n'.format(float(maxWaveExt)), append=True)

		if run_type != 'local':
			memcore2 = int(normal_cores_to_use)*4000
			memdiff2 = float(memcore2 - memory_usage_noe)
			
			memcore1 = int(cores_to_use)*4000
			memdiff1 = float(memcore1 - memory_usage)

#			driver.put('output.string(failure).current','\nNumber of Cores Used: {0}\n'.format(normal_cores_to_use), append=True)
#			driver.put('output.string(failure).current','\nPredicted Non-Efield Memory Usage: {0} MB Needed, {1} MB Requested, {2} MB Surplus\n'.format(memory_usage_noe, memcore2, memdiff2), append=True)

#			driver.put('output.string(failure).current','\nNumber of Cores Used For Nearfield Calculation: {0}\n'.format(cores_to_use), append=True)
#			driver.put('output.string(failure).current','\nPredicted Nearfield Memory Usage: {0} MB Needed, {1} MB Requested, {2} MB Surplus\n'.format(memory_usage, memcore1, memdiff1), append=True)



	# Prepare plots for output, guarantees first outputs in menu. If there is no error message.
	#qtable
    driver.put('output.curve(qtable).about.label','Light Extinction EF vs. Wavelength')
    driver.put('output.curve(qtable).about.description','The plot corrosponding to the extinction behavior denoted numerically in qtable.')
    driver.put('output.curve(qtable).xaxis.label','Wavelength')
    driver.put('output.curve(qtable).xaxis.units','uM')
    driver.put('output.curve(qtable).yaxis.label','Light Extinction Efficiency Factor')
    #qtable2
    driver.put('output.curve(qtable2).about.label','Light Absorption EF vs. Wavelength')
    driver.put('output.curve(qtable2).about.description','The plot corrosponding to the absorption behavior denoted numerically in qtable.')
    driver.put('output.curve(qtable2).xaxis.label','Wavelength')
    driver.put('output.curve(qtable2).xaxis.units','uM')
    driver.put('output.curve(qtable2).yaxis.label','Light Absorption Efficiency Factor')
    #qtable3
    driver.put('output.curve(qtable3).about.label','Light Scattering EF vs. Wavelength')
    driver.put('output.curve(qtable3).about.description','The plot corrosponding to the scattering behavior denoted numerically in qtable.')
    driver.put('output.curve(qtable3).xaxis.label','Wavelength')
    driver.put('output.curve(qtable3).xaxis.units','uM')
    driver.put('output.curve(qtable3).yaxis.label','Light Scattering Efficiency Factor')
    #qtable4
    driver.put('output.curve(qtable4).about.label','Phase Lag EF vs. Wavelength')
    driver.put('output.curve(qtable4).about.description','The plot corrosponding to the behavior denoted numerically in qtable2.')
    driver.put('output.curve(qtable4).xaxis.label','Wavelength')
    driver.put('output.curve(qtable4).xaxis.units','uM')
    driver.put('output.curve(qtable4).yaxis.label','Phase Lag Efficiency Factor')
		
    # If the E-field is to be calculated, run our lite version of ddpostprocess and display the results.
    no_display_flag = 0
    if (os.path.exists('shape.dat')) and (not os.path.exists('target.out')):
		os.rename('shape.dat','target.out')
    if (os.path.exists('target.out')):
		if (os.path.getsize('target.out') > 200000000):
			no_display_flag = 1


    if ((nearfield_calculate == '1') or (nearfield_calculate == '2')) and (ddscat_fail_flag != '1') and (Efield_fail_flag != '1'):
		ddpp_path = os.path.join(tool_path, 'bin/myddpostprocess')
		squareval = driver.get('input.phase(page3).group(NRFLD_HEAD).group(NRFLD_LINE).choice(NRFLD_IVTR).current')
		ddStatus, ddStdout, ddStderr = RapptureExec([ddpp_path, ddppfile_to_use, squareval], streamOutput=False)
		
		if ddStatus != 0:
			text = ddStdout + ddStderr
			log(text)
			sys.stderr.write(text)
			driver.result(ddStatus)
			sys.exit(ddStatus)
				
		RawDataFile = os.path.join(os.getcwd(),'field_datafile.txt')
		getSecret = driver.get('input.phase(page2).group(secretmenu).choice(secretdata).current')
		getSecretX = driver.get('input.phase(page2).group(secretmenu).number(fixX).current')
		getSecretY = driver.get('input.phase(page2).group(secretmenu).number(fixY).current')
		getSecretZ = driver.get('input.phase(page2).group(secretmenu).number(fixZ).current')
		min_x,min_y,min_z,min_e,max_x,max_y,max_z,max_e,num_pts_x,num_pts_y,num_pts_z =	BuildVTKfiles(RawDataFile, squareval, getSecret, getSecretX,getSecretY,getSecretZ)

		#Cleanup and printout for E-field processing

		remove_all_w(int(params['NWAV']))

		### Visualization for E-Field

		style1path1 = os.path.join(tool_path, 'rappture/dda/Styles')
		style1path2 = os.path.join(os.getcwd(), 'Styles')
		shutil.copyfile(style1path1,style1path2)		
			
		if os.path.exists('EField_VTK.vtk'):				
			driver.put('output.field(2).about.label', 'Electric Field (3D Field)', '0')
			with open('EField_VTK.vtk', 'r') as VTKfile:
				driver.put('output.field(2).component.vtk', VTKfile.read(), '0', '0')
			with open(os.path.realpath('Styles'), 'r') as Stylefile:
				driver.put('output.field(2).component.style', Stylefile.read(),'0')
			# apply axis labeling (i.e. units)
			# Note that the units are actually (nm), but confusingly the axes always display an extra (10^-3) so 
			# to counteract this I am just writing (um) and it looks like (um) x (10^-3)
			driver.put('output.field(2).about.xaxis.label', 'X (um)', '0')
			driver.put('output.field(2).about.yaxis.label', 'Y (um)', '0')
			driver.put('output.field(2).about.zaxis.label', 'Z (um)', '0')

		if os.path.exists('EField_VTK_logscale.vtk'):
			driver.put('output.field(4).about.label', 'Electric Field in Log Scale (3D Field)', '0')
			with open('EField_VTK_logscale.vtk', 'r') as VTKfile:
				driver.put('output.field(4).component.vtk', VTKfile.read(), '0', '0')
			with open(os.path.realpath('Styles'), 'r') as Stylefile:
				driver.put('output.field(4).component.style', Stylefile.read(),'0')

			driver.put('output.field(4).about.xaxis.label', 'X (nm)', '0')
			driver.put('output.field(4).about.yaxis.label', 'Y (nm)', '0')
			driver.put('output.field(4).about.zaxis.label', 'Z (nm)', '0')
			
		if os.path.exists('BField_VTK.vtk') and (Field_status == "2"):				
			driver.put('output.field(8).about.label', 'Magnetic Field (3D Field)', '0')
			with open('BField_VTK.vtk', 'r') as VTKfile:
				driver.put('output.field(8).component.vtk', VTKfile.read(), '0', '0')
			with open(os.path.realpath('Styles'), 'r') as Stylefile:
				driver.put('output.field(8).component.style', Stylefile.read(),'0')
			# apply axis labeling (i.e. units)
			# Note that the units are actually (nm), but confusingly the axes always display an extra (10^-3) so 
			# to counteract this I am just writing (um) and it looks like (um) x (10^-3)
			driver.put('output.field(8).about.xaxis.label', 'X (um)', '0')
			driver.put('output.field(8).about.yaxis.label', 'Y (um)', '0')
			driver.put('output.field(8).about.zaxis.label', 'Z (um)', '0')


	# If Vector Field requested, draw it
		vecILINE = driver.get('input.phase(page3).group(NRFLD_HEAD).group(vectorinfo).choice(NRFLD_VECTOR).current')
		if (nearfield_calculate == '0'):
			vecILINE = "0"
#		vspacing = 0.001 * float(driver.get('input.phase(page3).group(NRFLD_HEAD).group(vectorinfo).number(vecspacing).current'))

		if (vecILINE == "1" or vecILINE == "2"):
#					if (vspacing != 0):
#						num_pts_x = int(round(round((max_x - min_x)/float(vspacing),4)))
#						num_pts_y = int(round(round((max_y - min_y)/float(vspacing),4)))
#						num_pts_z = int(round(round((max_z - min_z)/float(vspacing),4)))

			# Prepare the vector mesh
			# Mesh
				driver.put('output.mesh(mymesh).about.label', 'Object Mesh')
				driver.put('output.mesh(mymesh).dim', '3')
			#    driver.put('output.mesh(mymesh).units', 'um')
				driver.put('output.mesh(mymesh).hide', 'yes')
				driver.put('output.mesh(mymesh).grid.xaxis.min',min_x)
				driver.put('output.mesh(mymesh).grid.xaxis.max',max_x)
				driver.put('output.mesh(mymesh).grid.xaxis.numpoints',num_pts_x)
				driver.put('output.mesh(mymesh).grid.yaxis.min',min_y)
				driver.put('output.mesh(mymesh).grid.yaxis.max',max_y)
				driver.put('output.mesh(mymesh).grid.yaxis.numpoints',num_pts_y)
				driver.put('output.mesh(mymesh).grid.zaxis.min',min_z)
				driver.put('output.mesh(mymesh).grid.zaxis.max',max_z)		
				driver.put('output.mesh(mymesh).grid.zaxis.numpoints',num_pts_z)
				# Note that the units are actually (nm), but confusingly the axes always display an extra (10^-3) so 
				# to counteract this I am just writing (um) and it looks like (um) x (10^-3)
				driver.put('output.field(myfield4).about.label', 'E-Field Vector Rendering')
				driver.put('output.field(myfield4).about.xaxis.label', 'X (um)', '0')
				driver.put('output.field(myfield4).about.yaxis.label', 'Y (um)', '0')
				driver.put('output.field(myfield4).about.zaxis.label', 'Z (um)', '0')
				driver.put('output.field(myfield4).about.view', 'glyphs')			
				driver.put('output.field(myfield4).component.mesh', 'output.mesh(mymesh)')
				driver.put('output.field(myfield4).component.elemtype', 'vectors')
				driver.put('output.field(myfield4).component.elemsize', '3')	

				with open('EField_Vec', 'r') as VECfile:
					driver.put('output.field(myfield4).component.values', VECfile.read(),append=True)

				if os.path.exists('BField_Vec') and (Field_status == "2"):
					
				# Prepare the vector mesh
				# Mesh
					driver.put('output.mesh(mymeshB).about.label', 'Object Mesh')
					driver.put('output.mesh(mymeshB).dim', '3')
				#    driver.put('output.mesh(mymesh).units', 'um')
					driver.put('output.mesh(mymeshB).hide', 'yes')
					driver.put('output.mesh(mymeshB).grid.xaxis.min',min_x)
					driver.put('output.mesh(mymeshB).grid.xaxis.max',max_x)
					driver.put('output.mesh(mymeshB).grid.xaxis.numpoints',num_pts_x)
					driver.put('output.mesh(mymeshB).grid.yaxis.min',min_y)
					driver.put('output.mesh(mymeshB).grid.yaxis.max',max_y)
					driver.put('output.mesh(mymeshB).grid.yaxis.numpoints',num_pts_y)
					driver.put('output.mesh(mymeshB).grid.zaxis.min',min_z)
					driver.put('output.mesh(mymeshB).grid.zaxis.max',max_z)		
					driver.put('output.mesh(mymeshB).grid.zaxis.numpoints',num_pts_z)
					# Note that the units are actually (nm), but confusingly the axes always display an extra (10^-3) so 
					# to counteract this I am just writing (um) and it looks like (um) x (10^-3)
					driver.put('output.field(myfield4B).about.label', 'B-Field Vector Rendering')
					driver.put('output.field(myfield4B).about.xaxis.label', 'X (um)', '0')
					driver.put('output.field(myfield4B).about.yaxis.label', 'Y (um)', '0')
					driver.put('output.field(myfield4B).about.zaxis.label', 'Z (um)', '0')
					driver.put('output.field(myfield4B).about.view', 'glyphs')			
					driver.put('output.field(myfield4B).component.mesh', 'output.mesh(mymeshB)')
					driver.put('output.field(myfield4B).component.elemtype', 'vectors')
					driver.put('output.field(myfield4B).component.elemsize', '3')	
						
					
					with open('BField_Vec', 'r') as VECfile:
						driver.put('output.field(myfield4B).component.values', VECfile.read(),append=True)
						
						
				# Prepare the Poynting vector mesh
				# Mesh
					driver.put('output.mesh(mymeshBC).about.label', 'Object Mesh')
					driver.put('output.mesh(mymeshBC).dim', '3')
				#    driver.put('output.mesh(mymesh).units', 'um')
					driver.put('output.mesh(mymeshBC).hide', 'yes')
					driver.put('output.mesh(mymeshBC).grid.xaxis.min',min_x)
					driver.put('output.mesh(mymeshBC).grid.xaxis.max',max_x)
					driver.put('output.mesh(mymeshBC).grid.xaxis.numpoints',num_pts_x)
					driver.put('output.mesh(mymeshBC).grid.yaxis.min',min_y)
					driver.put('output.mesh(mymeshBC).grid.yaxis.max',max_y)
					driver.put('output.mesh(mymeshBC).grid.yaxis.numpoints',num_pts_y)
					driver.put('output.mesh(mymeshBC).grid.zaxis.min',min_z)
					driver.put('output.mesh(mymeshBC).grid.zaxis.max',max_z)		
					driver.put('output.mesh(mymeshBC).grid.zaxis.numpoints',num_pts_z)
					# Note that the units are actually (nm), but confusingly the axes always display an extra (10^-3) so 
					# to counteract this I am just writing (um) and it looks like (um) x (10^-3)
					driver.put('output.field(myfield4BC).about.label', 'Poynting Vector Rendering')
					driver.put('output.field(myfield4BC).about.xaxis.label', 'X (um)', '0')
					driver.put('output.field(myfield4BC).about.yaxis.label', 'Y (um)', '0')
					driver.put('output.field(myfield4BC).about.zaxis.label', 'Z (um)', '0')
					driver.put('output.field(myfield4BC).about.view', 'glyphs')			
					driver.put('output.field(myfield4BC).component.mesh', 'output.mesh(mymeshBC)')
					driver.put('output.field(myfield4BC).component.elemtype', 'vectors')
					driver.put('output.field(myfield4BC).component.elemsize', '3')	
						
						
					with open('Poynting_Vec', 'r') as VECfile:
						driver.put('output.field(myfield4BC).component.values', VECfile.read(),append=True)


# Remove other temp files here too		
#		os.remove('EField_VTK')				

			
		elif os.path.exists('EField_VTK.vtk') and (no_display_flag == 1):
			
			try:
				import zlib
				mode = zipfile.ZIP_DEFLATED
			except:
				mode = zipfile.ZIP_STORED

			try:
				zip = zipfile.ZipFile('zipfilename_VTK','w',mode)
				zip.write('EField_VTK.vtk')
				zip.close()
				
				driver.put('output.string(EField_VTK).current','zipfilename_VTK',type='file',compress=True)
			except NameError:
				pass

			driver.put('output.string(failure).current','\nThe file "EField_VTK" was too large to print properly,\n', append=True)
			driver.put('output.string(failure).current','however a zipped binary can still be downloaded via the download button.\n', append=True)
			driver.put('output.string(failure).current','After downloading, the extension must be changed from .dat to .zip and unzipped.', append=True)
			os.remove('EField_VTK.vtk')

		
    # Output timing information.
    if collect_timing == 'yes':
        prefix = 'output.string(timing)'
        driver.put(prefix + '.about.label', 'Timing Information')
        driver.put(prefix + '.current', timing_info)


    # Avoid writing files that crash Rappture:
    if (os.path.exists('target.out')):
		if (os.path.getsize('target.out') > 200000000):
			try:
				import zlib
				mode = zipfile.ZIP_DEFLATED
			except:
				mode = zipfile.ZIP_STORED

			try:
				zip = zipfile.ZipFile('zipfilename','w',mode)
				zip.write('target.out')
				zip.close()
				driver.put('output.string(target.out).about.label','target.out'+" (DDSCAT)")
				driver.put('output.string(target.out).current','zipfilename',type='file',compress=True)
			except NameError:
				pass

		
			os.remove('target.out')
			driver.put('output.string(failure).current','\nThe file "target.out" was too large to print properly (200MB+),\n', append=True)
			driver.put('output.string(failure).current','however a zipped binary can still be downloaded via the download button.\n', append=True)
			driver.put('output.string(failure).current','After downloading, the extension must be changed from .dat to .zip and unzipped.', append=True)


    for filename in ('mtable', 'qtable', 'qtable2', 'target.out','ddscat.par'):
        prefix = 'output.string(%s)' % (filename,)
        driver.put(prefix + '.about.label', filename + " (DDSCAT)")
        if os.path.exists(filename):
			with open(filename, 'r') as output_file:
				driver.put(prefix + '.current', output_file.read())
			output_file = ""

    # Read the strings from the qtable
    qtable_plot_data=[0,'.']
    qtable2_plot_data=[0,'.']
    read_buffer=0

	# If qtable doesn't actually exist, the plots below will crash.
	# The easiest way to avoid this is to make a fake qtable in this case.
    if not (os.path.exists('qtable')):
		with open('qtable','w') as fake_table:
			fake_table.write('')
    if not (os.path.exists('qtable2')):
		with open('qtable2','w') as fake_table:
			fake_table.write('')

    with open('qtable','r') as plot_input_file:
	for line in plot_input_file:
		a = line
		if (re.search('wave       Q_ext',a)):
			read_buffer=1
		if (read_buffer==1):
			qtable_plot_data.append(a)

    with open('qtable2','r') as plot_input_file:
		read_buffer = 0
		for line in plot_input_file:
			a = line
			if (read_buffer==1):
				qtable2_plot_data.append(a)
			if (re.search('wave      Q_pha',a)):
				read_buffer=1


    xy4 = ['0 0']
    for item in qtable2_plot_data:
		try:
			qpha = item.split()[2]
			wavepha = item.split()[1]
			if ((not math.isnan(float(qpha))) and (not math.isnan(float(wavepha)))):
				xy4.append(('\n{0} {1}').format(wavepha, qpha))						
		except (IndexError, AttributeError):
			pass

			

    # Write the data points from the Qtable to file
    temp_counter = 0
    item_number = 0
    with open('qtable_data', 'w') as data_file:
	        for item in qtable_plot_data:
			item_number = item_number + 1
			if temp_counter == 1:
	   	 		data_file.write(item)
			temp_counter = 1

    # Prepare the 3 sets of xy values to send to interface
    # Want line corrosponding to item_number = 3, continuing until index end
    xy1 = ['0 0']
    xy2 = ['0 0']
    xy3 = ['0 0']
    reader_fail_flag = '0'
    element_counter = 0
    for element in qtable_plot_data:
	if (element_counter >=3):
		temp_counter2 = 0
		for word in qtable_plot_data[element_counter].split():
			if (len(word) <= 11):				
				if (temp_counter2 == 1):
					xval = word
				if (temp_counter2 == 2):
					y1val = word
				if (temp_counter2 == 3):
					y2val = word
				if (temp_counter2 == 4):
					y3val = word
				temp_counter2 = temp_counter2 + 1
			if (len(word) > 11):
				reader_fail_flag = '1'
				xval = 0
				y1val = 0
				y2val = 0
				y3val = 0

		if ((not math.isnan(float(xval))) and (not math.isnan(float(y1val)))):
			xy1.append(('\n{0} {1}').format(xval, y1val))
		else:
			reader_fail_flag = '1'
		if ((not math.isnan(float(xval))) and (not math.isnan(float(y2val)))):
			xy2.append(('\n{0} {1}').format(xval, y2val))
		else:
			reader_fail_flag = '1'
		if ((not math.isnan(float(xval))) and (not math.isnan(float(y3val)))):
			xy3.append(('\n{0} {1}').format(xval, y3val))
		else:
			reader_fail_flag = '1'
	element_counter = element_counter + 1

    # Plot the three desired wavelength crossections.
    read_buffer_xy1 = 0
    read_buffer_xy2 = 0
    read_buffer_xy3 = 0
    read_buffer_xy4 = 0
    for pair in xy1:
	    if read_buffer_xy1 == 1:
	    	driver.put('output.curve(qtable).component.xy', pair, append=True)
	    read_buffer_xy1 = 1
    for pair in xy2:
	    if read_buffer_xy2 == 1:
       	        driver.put('output.curve(qtable2).component.xy', pair, append=True)
	    read_buffer_xy2 = 1
    for pair in xy3:
	    if read_buffer_xy3 == 1:
   	        driver.put('output.curve(qtable3).component.xy', pair, append=True)
	    read_buffer_xy3 = 1
    for pair in xy4:
	    if read_buffer_xy4 == 1:
   	        driver.put('output.curve(qtable4).component.xy', pair, append=True)
	    read_buffer_xy4 = 1


    if reader_fail_flag == '1' and ddscat_fail_message == '':				
		driver.put('output.string(failure).about.label','***SIMULATION STATUS***')
		driver.put('output.string(failure).current','<<< DDSCAT has encountered an error in processing. >>>\n\n Resultant data displayed for the settings attempted is void. \n')
		driver.put('output.string(failure).current','\nCrash Log:\n', append=True)
		driver.put('output.string(failure).current',' Inappropriate Values were returned for Light Absorption, Scattering, and/or Extinction', append=True)
		driver.put('output.string(failure).current','\n See mtable, qtable, qtable2 for more information on the Inappropriate Values', append=True)
		driver.put('output.string(failure).current', '\n# Nearfield Not Requested or Not Usable', append=True)
		reader_fail_flag = '2'

    if os.path.exists('submit_results/submit_log'):
		driver.put('output.string(sublog).about.label', 'Remote Submission Log')
		with open('submit_results/submit_log','r') as sublogfile:			
			driver.put('output.string(sublog).current', sublogfile.read())
		

		if os.path.exists('submit_results_efield/submit_log'):
			with open('submit_results_efield/submit_log','r') as sublogfile:			
				driver.put('output.string(sublog).current', '\n\n Nearfield: \n', append=True)
				driver.put('output.string(sublog).current', sublogfile.read(), append=True)

    email = driver.get('input.phase(page5).group(process).boolean(email).current')
    if email == "yes":
		command = ["submit"]
		command.append("--progress")
		command.append("silent")
		command.append("mail2self")
		command.append("-t")
		command.append("Please check Nanohub.org for your simulation results")
		command.append("-s")
		if ddscat_fail_flag == "1" or Efield_fail_flag == "1" or reader_fail_flag == "2":
			command.append("nanoDDSCAT+ Simulation #{0} Failed".format(driver_number))
		else:
			command.append("nanoDDSCAT+ Simulation #{0} Completed".format(driver_number))
		# Send out an email about the remote submission
		exit_status, stdout, stderr = RapptureExec(command, streamOutput=False)

    driver.result(0)

#   Remove created files from the working directory as post-processing cleanup.
    remove_all_w(int(params['NWAV']))
    for filename in ('custom_dielectric1','custom_dielectric2','custom_dielectric3','custom_dielectric4',
					'custom_dielectric5','custom_dielectric6','custom_dielectric7','custom_dielectric8',
					'custom_dielectric9', 'mtable','qtable', 'qtable2', 'qtable_data',
					'Styles', 'field_datafile.txt',
					'EField_VTK.vtk','BField_VTK.vtk','EField_VTK_logscale.vtk',
					'EField_Vec','BField_Vec','Poynting_Vec',
					'AuDiel.tab','AgDiel.tab','TiO2','SiO2','Pt_diel.tab','Pd_diel.tab','Cu_diel.txt',
					'shape.dat', 'ddscat.par','target.out','ddpostprocess.par','zipfilename',
					'zipfilename_VTK', 'composition.txt'):
        if os.path.exists(filename):
#			1
			os.remove(filename)

    if os.path.exists('submit_results'):
		shutil.rmtree('submit_results', ignore_errors = True)
    if os.path.exists('submit_results_efield'):
		shutil.rmtree('submit_results_efield', ignore_errors = True)
