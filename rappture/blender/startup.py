# blender --python startup.py
# AbderRahman N. Sobh - 2015
#
# Note: Script for generating farthest distances between vertices is
# credited to user "67853-Muffy" from blenderartists.org

import bpy
import sys
import os
import time
import glob
import tempfile
import shutil
import tarfile
from math import *

import bmesh
from bpy.props import FloatProperty, IntProperty, BoolProperty, StringProperty, CollectionProperty
from math import pi
from mathutils import Quaternion, Vector
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from bpy_extras.io_utils import ImportHelper

class Pyramid(bpy.types.Operator):
    """Add a 3-layer device mesh"""
    bl_idname = "object.device"
    bl_label = "Create a Device"
    bl_options = {'REGISTER', 'UNDO'}

    #steps = IntProperty(name="Steps",description="Number of items to create",min=3, max=3, default=3  )
    steps = 3
    Z_offset12 = FloatProperty(name="Z offset Bot-Mid",description="Z offset12",min=-100.0, max=100.0, default=0.0 )
    Z_offset23 = FloatProperty(name="Z offset Mid-Top",description="Z offset23",min=-100.0, max=100.0, default=0.0 )	
	
    dX1 = FloatProperty(name="Bottom X",description="Bottom X",min= 0.0, max=100.0, default=1.0 )
    dY1 = FloatProperty(name="Bottom Y",description="Bottom Y",min= 0.0, max=100.0, default=1.0 )
    dZ1 = FloatProperty(name="Bottom Z",description="Bottom Z",min= 0.0, max=100.0, default=1.0 )
	
    dX2 = FloatProperty(name="Mid X",description="Mid X",min= 0.0, max=100.0, default=1.0 )
    dY2 = FloatProperty(name="Mid Y",description="Mid Y",min= 0.0, max=100.0, default=1.0 )
    dZ2 = FloatProperty(name="Mid Z",description="Mid Z",min= 0.0, max=100.0, default=1.0 )

    dX3 = FloatProperty(name="Top X",description="Top X",min= 0.0, max=100.0, default=1.0 )
    dY3 = FloatProperty(name="Top Y",description="Top Y",min= 0.0, max=100.0, default=1.0 )
    dZ3 = FloatProperty(name="Top Z",description="Top Z",min= 0.0, max=100.0, default=1.0 )
	
    
    linked = BoolProperty(name="Linked",description="created objects are linked to the original",default =False)
    apply_t = BoolProperty(name="Apply Transform",description="apply object transform",default =False)

    def execute(self, context):
        steps = self.steps
        Z_offset12 = self.Z_offset12
        Z_offset23 = self.Z_offset23		
        dX1 = self.dX1
        dX2 = self.dX2
        dX3 = self.dX3		
        dY1 = self.dY1
        dY2 = self.dY2
        dY3 = self.dY3
        dZ1 = self.dZ1		
        dZ2 = self.dZ2		
        dZ3 = self.dZ3				
        linked = self.linked
        apply_t=self.apply_t

        dX = [dX1, dX2, dX3]
        dY = [dY1, dY2, dY3]
        dZ = [dZ1, dZ2, dZ3]
        dobj = [0,0,0]

        bpy.ops.mesh.primitive_cube_add(radius=0, view_align=False, enter_editmode=False, location=(0, 0, 0), rotation=(0, 0, 0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
        oritemp = bpy.context.object
		
        for x in range(0, steps):

            if linked :
                bpy.ops.object.duplicate_move_linked()
            else :
                bpy.ops.object.duplicate()

            bpy.context.object.dimensions = (dX[x], dY[x], dZ[x])

            dobj[x] = bpy.context.object                        
            if x == 0:
                dobj[x].location = dobj[0].location + Vector((0,0,(dZ[0]/float(2))))
            elif x == 1:
                dobj[x].location = dobj[1].location + Vector((0,0,(dZ[1]/float(2))))
            elif x == 2:
                dobj[x].location = dobj[2].location + Vector((0,0,(dZ[2]/float(2))))




        for x in range(0, steps):
            if x == 0:
                dobj[x].location = dobj[x].location
            if x == 1:
                dobj[x].location = dobj[x].location + Vector((0,0,(dZ[0]/float(2))+Z_offset12))
            elif x == 2:
#                dobj[x].location = dobj[x].location + Vector((0,0,(dZ[1]/float(2))+Z_offset23))
                dobj[x].location = dobj[x].location + dobj[0].location + Vector((0,0,(dZ[1]/float(2))+Z_offset23+Z_offset12))
				

            if apply_t and not linked:
                bpy.ops.object.transform_apply(scale =True)

        bpy.context.scene.objects.unlink(oritemp)
        oritemp.user_clear()
        bpy.data.objects.remove(oritemp)
            
        return {'FINISHED'}

   


#Find BBox Center
def findBBoxCenter(OBJ) :
    
    vX = 0.0
    vY = 0.0
    vZ = 0.0
    
    for v in OBJ.bound_box :
        vX += v[0]
        vY += v[1]
        vZ += v[2]
    
    #bbox center, object space    
    BBC = ( vX/8 , vY/8 , vZ/8 )
    
    return BBC



#distance from point
def FindDistanceFromPoint( SP , TP ) : #SamplePoint, TargetPoint
    
    #subtract target point from sample point
    #this makes sample point relative to target point
    V = [0,0,0]
    for c in range(3):
        V[c] = SP[c] - TP[c]

    #find magnitude of SamplePoint, now relative to Target Point
    mag = pow(V[0],2) + pow(V[1],2) + pow(V[2],2)
        
    return sqrt(mag)

class ToolPropsPanel(bpy.types.Panel):
	bl_label = "Export Mesh to DDSCAT Convert"
	bl_space_type = "VIEW_3D"
	bl_region_type = "TOOL_PROPS"
	
	def draw(self, context):
		self.layout.operator("center.export")
		self.layout.operator("resetrotate.export")
		self.layout.operator("rotate.export")
		self.layout.operator("lock.export")
		self.layout.operator("dda.export")


class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "center.export"
	bl_label = "1.) Center Scene at Origin"
	
	def execute(self, context):

		# Compute the Centroid of the scene.
		object_counter = 0
		centx = 0
		centy = 0
		centz = 0
		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				bpy.context.scene.objects.active = item
				bpy.ops.object.mode_set(mode = 'OBJECT')
				bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
#				bpy.ops.object.mode_set(mode = 'EDIT')
				object_counter = object_counter + 1
				centx = centx + item.location[0]
				centy = centy + item.location[1]
				centz = centz + item.location[2]
		centx = centx/float(object_counter)
		centy = centy/float(object_counter)
		centz = centz/float(object_counter)
		
		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				item.location = [item.location[0]-centx, item.location[1]-centy, item.location[2]-centz]
			
		return{'FINISHED'}

class ImportTAR(bpy.types.Operator, ImportHelper):
    """Load STL triangle mesh data from a .tar archive"""
    bl_idname = "import_mesh.tar"
    bl_label = "Import tar archive"
    bl_options = {'UNDO'}

    filename_ext = ".tar"

    filter_glob = StringProperty(
            default="*.tar",
            options={'HIDDEN'},
            )
    files = CollectionProperty(
            name="File Path",
            type=bpy.types.OperatorFileListElement,
            )
    directory = StringProperty(
            subtype='DIR_PATH',
            )

    def execute(self, context):
	# Having trouble calling the relative import for stl_utils, blender_utils
	# So I'm just going to call them directly from filepath

        SCRIPT_DIR = "/apps/share64/debian7/blender/blender-2.70a/2.70/scripts/addons/io_mesh_stl/"
        sys.path.append(os.path.normpath(SCRIPT_DIR))
        import stl_utils
        import blender_utils

        paths = [os.path.join(self.directory, name.name)
                 for name in self.files]
        subpaths = []
        truepaths = []

        if not paths:
            paths.append(self.filepath)

        if bpy.ops.object.mode_set.poll():
            bpy.ops.object.mode_set(mode='OBJECT')

        if bpy.ops.object.select_all.poll():
            bpy.ops.object.select_all(action='DESELECT')
            
        for path in paths:
            with tarfile.open(path,'r') as infile:
                subpaths = infile.getnames()
                delpath = subpaths.pop(0)
                infile.extractall()

        for item in subpaths:
            truepaths.append(os.path.realpath(item))

        for path in truepaths:
            objName = bpy.path.display_name(os.path.basename(path))
            tris, pts = stl_utils.read_stl(path)
            blender_utils.create_and_link_mesh(objName, tris, pts)
            
        shutil.rmtree(os.path.realpath(delpath))            

        return {'FINISHED'}



class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "center.export"
	bl_label = "1.) Center Scene at Origin"
	
	def execute(self, context):

		# Compute the Centroid of the scene.
		object_counter = 0
		centx = 0
		centy = 0
		centz = 0
		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				bpy.context.scene.objects.active = item
				bpy.ops.object.mode_set(mode = 'OBJECT')
				bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
#				bpy.ops.object.mode_set(mode = 'EDIT')
				object_counter = object_counter + 1
				centx = centx + item.location[0]
				centy = centy + item.location[1]
				centz = centz + item.location[2]
		centx = centx/float(object_counter)
		centy = centy/float(object_counter)
		centz = centz/float(object_counter)
		
		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				item.location = [item.location[0]-centx, item.location[1]-centy, item.location[2]-centz]
			
		return{'FINISHED'}


class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "resetrotate.export"
	bl_label = "2.) Reset Lights"
	

	def execute(self, context):
		for item in bpy.data.objects:
			item.select = False			
		bpy.data.objects['RotateIncidentLight'].select = True
		bpy.context.space_data.show_manipulator = True
		bpy.context.space_data.transform_manipulators = {'ROTATE'}
		bpy.data.objects['RotateIncidentLight'].rotation_euler[0] = 0
		bpy.data.objects['RotateIncidentLight'].rotation_euler[1] = 0
		bpy.data.objects['RotateIncidentLight'].rotation_euler[2] = 0
		return{'FINISHED'}


class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "rotate.export"
	bl_label = "3.) Rotate Lights/Axes"
	
	def execute(self, context):
		bpy.ops.object.mode_set(mode="OBJECT")
		for item in bpy.data.objects:
			item.select = False			
		bpy.data.objects['RotateIncidentLight'].select = True
		bpy.context.space_data.show_manipulator = True
		bpy.context.space_data.transform_manipulators = {'ROTATE'}
		return{'FINISHED'}


class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "lock.export"
	bl_label = "4.) Lock/Unlock Incident Light"
	
	def execute(self, context):
		bpy.ops.object.mode_set(mode="OBJECT")

		for item in bpy.data.objects:
			item.select = False			
		bpy.data.objects['RotateIncidentLight'].select = True


		bpy.context.space_data.show_manipulator = True
		for area in bpy.context.screen.areas:
			if area.type == "VIEW_3D":
				area.spaces[0].transform_orientation = 'LOCAL'
				if bpy.data.objects['RotateIncidentLight'].lock_rotation[1] == False:
					bpy.data.objects['RotateIncidentLight'].lock_rotation[1] = True
					bpy.data.objects['RotateIncidentLight'].lock_rotation[2] = True
					bpy.context.space_data.transform_orientation = 'LOCAL'
				elif bpy.data.objects['RotateIncidentLight'].lock_rotation[1] == True:
					bpy.data.objects['RotateIncidentLight'].lock_rotation[1] = False
					bpy.data.objects['RotateIncidentLight'].lock_rotation[2] = False
					bpy.context.space_data.transform_orientation = 'GLOBAL'
		bpy.context.space_data.transform_manipulators = {'ROTATE'}

		return{'FINISHED'}


class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "line.export"
	bl_label = "Backend.) Draw Longest Dimension Span"

	def execute(self, context):	

		# deselect everything
		for item in bpy.data.objects:
			item.select = False
		# select everything
		for item in bpy.data.objects:
			item.select = False

		# delete any previous draws of the MaxDist
		try:
			type(bpy.data.objects['MaxDist'])
			for item in bpy.data.objects:
				item.select = False

			bpy.data.objects['MaxDist'].select = True
			bpy.ops.object.mode_set(mode="OBJECT")
			bpy.ops.object.delete()
			return{'FINISHED'}

		except KeyError:
			pass

		# delete any previous draws of the SysDupe
		try:
			type(bpy.data.objects['SysDupe'])
			for item in bpy.data.objects:
				item.select = False

			bpy.data.objects['SysDupe'].select = True
			bpy.ops.object.mode_set(mode="OBJECT")
			bpy.ops.object.delete()

		except KeyError:
			pass

		try:
			bpy.ops.object.mode_set(mode="OBJECT")
		except RuntimeError:
			pass

		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				item.select = True
				bpy.context.scene.objects.active = item
			else:
				item.select = False

		# make a duplicate object of the system, combine all items, and pass the live selection to it
		bpy.ops.object.duplicate()
		try:
			bpy.ops.object.join()
		except RuntimeError:
			pass
		bpy.context.object.name = "SysDupe"
		for item in bpy.data.objects:
			item.select = False
		bpy.data.objects["SysDupe"].select = True
		bpy.context.scene.objects.active = bpy.data.objects['SysDupe']

		#Collect Object Data
		AO = bpy.data.objects["SysDupe"]
		AOV = AO.data.vertices

		#--- find farthest vertex from bounding box center
		# This is the outermost vertex in object space.
		# Find bounding box center.
		BBC = findBBoxCenter(AO)
		#Find farthest vertex from BBC
		MaxDist = 0.0
		MaxDistIndex = 0

		for v in AOV :
			d = FindDistanceFromPoint(v.co,BBC)
			if d > MaxDist :
				MaxDist = d
				MaxDistIndex = v.index


		#--- Vertex[MaxDistIndex] is the farthest from the center of the bounding box
		#--- Now, find farthest vertex from this vertex

		MaxDist2 = 0.0
		MaxDistIndex2 = 0

		for v in AOV :
			d = FindDistanceFromPoint( v.co , AOV[MaxDistIndex].co )
			if d > MaxDist2 :
				MaxDist2 = d
				MaxDistIndex2 = v.index
				

		AOV[MaxDistIndex].select = True
		AOV[MaxDistIndex2].select = True
		bpy.data.scenes[0].update()


		me = bpy.data.meshes.new('MaxDist')
		me.vertices.add(2)
		me.vertices[0].co = AOV[MaxDistIndex].co
		me.vertices[1].co = AOV[MaxDistIndex2].co

		# transform to worldspace coords
		mat = AO.matrix_world
		me.transform(mat)
		
		me.edges.add(1)
		me.edges[0].vertices = (0,1)

		ob = bpy.data.objects.new('MaxDist',me)
		scn = bpy.context.scene
		scn.objects.link(ob)

		# delete the SysDupe
		try:
			type(bpy.data.objects['SysDupe'])
			for item in bpy.data.objects:
				item.select = False

			bpy.data.objects['SysDupe'].select = True
			bpy.ops.object.mode_set(mode="OBJECT")
			bpy.ops.object.delete()

		except KeyError:
			pass

		for item in bpy.data.objects:
			item.select = False

		bpy.data.objects['MaxDist'].select = True
		bpy.context.scene.objects.active = bpy.data.objects['MaxDist']
		bpy.ops.object.mode_set(mode="EDIT")
		bpy.context.space_data.use_occlude_geometry = False


		return{'FINISHED'}
		
class OBJECT_OT_HelloButton(bpy.types.Operator):
	bl_idname = "dda.export"
	bl_label = "5.) Export .obj for DDAConvert"
	
	def execute(self, context):

		# If scene was not centered, center it
		bpy.ops.center.export()
			# Check if MaxDist was drawn, if so take its value.
			# If not, then draw it and take its value.
		try:
			entry2 = bpy.data.objects['MaxDist'].dimensions
			x2,y2,z2 = entry2
			x2,y2,z2 = float(x2),float(y2),float(z2)
			objectdistance = sqrt((x2)**2 + (y2)**2 + (z2)**2)

		except KeyError:
			bpy.ops.object.mode_set(mode="OBJECT")
			# call maxDist
			bpy.ops.line.export()
			# get Value
			entry2 = bpy.data.objects['MaxDist'].dimensions
			x2,y2,z2 = entry2
			x2,y2,z2 = float(x2),float(y2),float(z2)
			objectdistance = sqrt((x2)**2 + (y2)**2 + (z2)**2)		
		
		# delete any previous draws of the MaxDist
		try:
			type(bpy.data.objects['MaxDist'])
			for item in bpy.data.objects:
				item.select = False

			bpy.data.objects['MaxDist'].select = True
			bpy.ops.object.mode_set(mode="OBJECT")
			bpy.ops.object.delete()

		except KeyError:
			pass


		
		object_locations = []
		for item in bpy.data.objects:
			if item.name != "Incident Light" and item.name != "----Y Polarization Direction" and item.name != "RotateIncidentLight" and item.name != "MaxDist" and item.name != "SysDupe":
				bpy.context.scene.objects.active = item
				bpy.ops.object.mode_set(mode = 'EDIT')
				bpy.ops.object.mode_set(mode = 'EDIT')
				bpy.ops.mesh.select_all(action='SELECT')
				bpy.ops.mesh.quads_convert_to_tris()
				object_locations.append("{0} {1} {2}".format(item.location[0],item.location[1],item.location[2]))
			if item.name == "Incident Light":
				lightvectora = item.location * item.matrix_world
				lightvector = [round(lightvectora[0],2), round(lightvectora[1],2), -1 * round(lightvectora[2],2)]
			if item.name == "----Y Polarization Direction":
				polarvectora = item.location * item.matrix_world
				polarvector = [round(polarvectora[0],2), round(polarvectora[1],2), -1 * round(polarvectora[2],2)]

		
		curr_time = '{0}'.format(int(time.time()))


		# create a fresh shuttle file from the mesh file, clear out all other shuttle-associated files
		for fil in glob.glob(os.getcwd()+'/PolarLight*'):		
			os.remove(fil)		    		
		LIGHT_POLAR_PATH = os.path.join(os.getcwd(),'PolarLight'+curr_time+'.txt')
		with open(LIGHT_POLAR_PATH, 'w') as lightf:
			lightf.write('light: {0} {1} {2}'.format(lightvector[0], lightvector[1], lightvector[2]))
			lightf.write('\npolar: {0} {1} {2}'.format(polarvector[0], polarvector[1], polarvector[2]))
			lightf.write('\nmaxdist: {0}'.format(objectdistance))
			lightf.write('\nmdxyz: {0} {1} {2}'.format(x2,y2,z2))


		# create a fresh shuttle file from the mesh file, clear out all other shuttle-associated files
		for fil in glob.glob(os.getcwd()+'/saved_mesh*'):
			os.remove(fil)		    
		MESH_OBJ_PATH = os.path.join(os.getcwd(),'saved_mesh'+curr_time+'.obj')
		bpy.ops.export_scene.obj(filepath=MESH_OBJ_PATH)

		sessionnum = os.getcwd().split('/')[-1]

		MESH_MTL_PATH = MESH_OBJ_PATH.split('.obj')[0]+'.mtl'
		if os.path.exists(MESH_MTL_PATH):
			os.remove(MESH_MTL_PATH)


		return{'FINISHED'}


# make operator available in layout
def menu_func(self, context):
    self.layout.operator("object.device", 
        text="Device", 
        icon='SPACE3')

def menu_func_import(self, context):
    self.layout.operator("import_mesh.tar", 
        text="Import STL Archive (as .tar)")

# update layout on (un)registering addon
def register():  
    bpy.utils.register_module(__name__)  
    bpy.types.INFO_MT_mesh_add.append(menu_func) 
    bpy.types.INFO_MT_file_import.append(menu_func_import)


register()
bpy.ops.wm.read_homefile()
bpy.context.user_preferences.view.manipulator_size = 85

		
bpy.utils.register_module(__name__)
