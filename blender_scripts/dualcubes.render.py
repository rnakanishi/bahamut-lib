
import bpy

scn = bpy.context.scene
surfmat = bpy.data.materials.get("surface")
wiremat = bpy.data.materials.get("wireframe")

for f in range(scn.frame_start, scn.frame_end+1):
    fpath = bpy.path.abspath(
        '/home/rnakanishi/git/bahamut-lib/results/dual_marching/{:04d}.obj'.format(f))
    objname = 'dualcubes_{:04d}'.format(f)
    # fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj1 = bpy.context.selected_objects[0]
    obj1.location[0] = 1
    print(obj1.name)
    obj1.data.materials.append(surfmat)

    bpy.ops.object.duplicate(linked=0)
    obj2 = bpy.context.selected_objects[0]
    print(obj2.name)
    obj2.data.materials[0] = wiremat
    bpy.context.scene.objects.active = obj2
    mod = obj2.modifiers.new("Wireframe", type='WIREFRAME')
    mod.thickness = 0.01
    
    fpath = bpy.path.abspath(
        '/home/rnakanishi/git/bahamut-lib/results/weno_cube/{:04d}.obj'.format(f))
    objname = 'dualcubes_{:04d}'.format(f)
    # fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
    bpy.ops.import_scene.obj(filepath=fpath)
    obj3 = bpy.context.selected_objects[0]
    obj3.location[0] = -1
    print(obj3.name)
    obj3.data.materials.append(surfmat)

    bpy.ops.object.duplicate(linked=0)
    obj4 = bpy.context.selected_objects[0]
    print(obj4.name)
    obj4.data.materials[0] = wiremat
    bpy.context.scene.objects.active = obj4
    mod = obj4.modifiers.new("Wireframe", type='WIREFRAME')
    mod.thickness = 0.01

    bpy.context.scene.render.image_settings.file_format = 'JPEG'
    bpy.context.scene.render.filepath = "/home/rnakanishi/Documents/blender/marching/frames/cube/{:04d}.jpg".format(
       f)
    bpy.context.scene.render.engine = "CYCLES"
    bpy.context.scene.render.resolution_percentage = 75
    bpy.context.scene.render.layers[0].cycles.use_denoising = True
    bpy.ops.render.render(use_viewport=True, write_still=True)

    objs = bpy.data.objects
    objs.remove(objs[obj1.name])
    objs.remove(objs[obj2.name])
    objs.remove(objs[obj3.name])
    objs.remove(objs[obj4.name])

