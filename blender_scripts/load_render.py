import bpy

scn = bpy.context.scene
mat = bpy.data.materials.get("Glass")
if mat is None:
    # create material
    mat = bpy.data.materials.new(name="Glass")
    
for f in range(scn.frame_start, scn.frame_end+1):
    fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/results/64-p2-va2-la2/{:04d}.obj'.format(f))
    # fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj = bpy.context.selected_objects[0]    

    bpy.context.scene.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.remove_doubles()
    
    bpy.ops.object.modifier_add(type='SMOOTH')
    bpy.context.object.modifiers["Smooth"].iterations = 20

    bpy.ops.object.editmode_toggle()
    bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Smooth")
    bpy.ops.export_scene.obj(filepath=fpath, use_normals=False, use_materials=False, use_selection=True)

    if obj.data.materials:
    # assig n to 1st material slot
        obj.data.materials[0] = mat
    else:
    # no slots
        obj.data.materials.append(mat)
    
    # Render commands
    bpy.ops.object.shade_smooth()
    bpy.context.scene.render.image_settings.file_format='JPEG'
    bpy.context.scene.render.filepath = "/home/rnakanishi/Documents/blender/ando/frames/64-p2/{:04d}.jpg".format(f)
    bpy.context.scene.render.engine = "CYCLES"
    bpy.context.scene.render.layers[0].cycles.use_denoising = True
    bpy.context.scene.render.resolution_percentage = 20
    bpy.ops.render.render(use_viewport = True, write_still=True)

    bpy.ops.object.delete(use_global=False)
