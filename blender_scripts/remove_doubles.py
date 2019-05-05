
import bpy

scn = bpy.context.scene
mat = bpy.data.materials.get("Glass")
if mat is None:
    # create material
    mat = bpy.data.materials.new(name="Glass")
    
for f in range(scn.frame_start, scn.frame_end):
    fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj = bpy.context.selected_objects[0]
    
    original_type = bpy.context.area.type
    bpy.context.area.type = "VIEW_3D"

    bpy.context.scene.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.remove_doubles()
    bpy.ops.object.editmode_toggle()
    
    bpy.context.area.type = original_type
    bpy.ops.import_scene.obj(filepath=fpath, use_normals=False)
    bpy.ops.object.delete(use_global=False)
