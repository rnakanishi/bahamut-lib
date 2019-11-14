import bpy

scn = bpy.context.scene
surfmat = bpy.data.materials.get("surface")
wiremat = bpy.data.materials.get("wireframe")
# dull = bpy.data.materials.get("dull_surface")
dull = bpy.data.materials.get("shiny")

base = '/home/rnakanishi/git/bahamut-lib/results/particles/3d/enright/'
folders = [
    # 'pls128_1',
    'semiL128_2', 'weno128_2', 'pls128_2'
]
bpy.context.scene.render.engine = "BLENDER_EEVEE"

for folder in folders:
    for f in range(scn.frame_start, scn.frame_end+1):
        fpath = bpy.path.abspath(
            # '/home/rnakanishi/git/bahamut-lib/results/redistance_cube/{:04d}.obj'.format(f))
            base + folder + '/mesh/n{:04d}.obj'.format(f))
        # '/home/rnakanishi/git/bahamut-lib/results/mesh/3d/enright/pls80/{:04d}.obj'.format(f))
        objname = 'n{:04d}'.format(f)
        # objname = bpy.context.selected_objects[0].name

        # fpath = bpy.path.abspath('/home/rnakanishi/git/bahamut-lib/obj/{:04d}.obj'.format(f))
        bpy.ops.import_scene.obj(filepath=fpath)
        obj1 = bpy.context.selected_objects[0]
        # obj1 = bpy.context.scene.objects[objname]
        # obj1 = bpy.context.selected_objects[0]
        # print(bpy.context.selected_objects[0].name)

        # obj1.location[0] = 2
        print(obj1.name)
        obj1.data.materials[0] = bpy.data.materials.get("shiny")
        # obj1.data.materials.append(bpy.data.materials.get("shiny"))
        print('Material: ' + obj1.data.materials[0].name)

        bpy.context.scene.render.image_settings.file_format = 'JPEG'
        bpy.context.scene.render.filepath = "/home/rnakanishi/Documents/blender/enright/" + \
            folder+"/{:04d}.jpg".format(f)
        bpy.context.scene.render.resolution_percentage = 100
        # bpy.context.scene.render.layers[0].cycles.use_denoising = True
        bpy.context.scene.view_layers[0].cycles.use_denoising = True

        bpy.ops.render.render(use_viewport=True, write_still=True)

        objs = bpy.data.objects
        objs.remove(objs[obj1.name])
        # objs.remove(objs[obj2.name])
        # objs.remove(objs[obj3.name])
        # objs.remove(objs[obj4.name])
