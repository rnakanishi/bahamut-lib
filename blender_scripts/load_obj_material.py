
import bpy

scn = bpy.context.scene
mat = bpy.data.materials.get("Surface")
if mat is None:
    # create material
    mat = bpy.data.materials.new(name="Surface")
    
for f in range(scn.frame_start, scn.frame_end+1):
    fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/results/marching/dualcubes_{:04d}.obj'.format(f))
    # fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj = bpy.context.selected_objects[0]
    
    if obj.data.materials:
    # assig n to 1st material slot
        obj.data.materials[0] = mat
    else:
    # no slots
        obj.data.materials.append(mat)
        
    # key as visible on the current frame
    obj.keyframe_insert('hide',frame=f)
    obj.keyframe_insert('hide_render',frame=f)
     # hide it
    obj.hide = True
    obj.hide_render = True
    # key as hidden on the previous frame
    obj.keyframe_insert('hide',frame=f-1)
    obj.keyframe_insert('hide_render',frame=f-1)
    # key as hidden on the next frame
    obj.keyframe_insert('hide',frame=f+1)
    obj.keyframe_insert('hide_render',frame=f+1)
