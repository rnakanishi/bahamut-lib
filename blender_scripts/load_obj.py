import bpy

scn = bpy.context.scene
for f in range(scn.frame_start, scn.frame_end):
    fpath = bpy.path.abspath('/home/rnakanishi/Documentos/blender/ando/128-p1-a1/{:04d}.obj'.format(f))
#     fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/results/test/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj = bpy.context.selected_objects[0]
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

