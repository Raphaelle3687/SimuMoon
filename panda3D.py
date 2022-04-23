from math import cos, sin, pi
import imageio # install imageio: pip install imageio
from panda3d_viewer import Viewer, ViewerConfig
import numpy as np
EARTH_RADIUS = 6.371*1e6 #radius in meters

config = ViewerConfig()
config.set_window_size(1000, 1000)
config.enable_antialiasing(True, multisamples=4)
config.enable_shadow(True)
config.show_axes(False)
config.show_grid(False)
config.show_floor(False)

viewer = Viewer(window_type='offscreen', config=config)

viewer.set_background_color((0, 0, 0))
#viewer.enable_lights( False)
#viewer.enable_light(1, True)

nFiles = 150


positionsArrays = []
radiusesArrays = []

for i in range(nFiles):
    positions=np.load("data/positions"+str(i)+".npy")
    radiuses=np.load("data/radiuses"+str(i)+".npy")
    positions/=EARTH_RADIUS
    radiuses/=EARTH_RADIUS
    positionsArrays.append(positions)
    radiusesArrays.append(radiuses)

def redoGroup(positions, radiuses, redo=False):
    if redo:
        viewer.remove_group("root")
    viewer.append_group('root')

    for i in range(len(positions)):
        viewer.append_sphere('root', 'sphere_node'+str(i), radius=radiuses[i])
        if i==0:
            viewer.set_material('root', 'sphere_node'+str(i), color_rgba=(0, 0.1, 0.5, 1))
        else:
            viewer.set_material('root', 'sphere_node'+str(i), color_rgba=(100/255, 60/255, 40/255, 1))

        viewer.move_nodes('root', {'sphere_node'+str(i): (positions[i], (1, 0, 0, 0))})

x = 60
y = 0 
z = 60
viewer.reset_camera(pos=(x, y, z), look_at=(0, 0, 0))
redoGroup(positionsArrays[0], radiusesArrays[0])

with imageio.get_writer('sphere_anim.gif', mode='I') as writer:

    image_rgb = viewer.get_screenshot(requested_format='RGB')
    writer.append_data(image_rgb)

    for i in range(1, len(positionsArrays)):
        redoGroup(positionsArrays[i], radiusesArrays[i], redo=True)
        image_rgb = viewer.get_screenshot(requested_format='RGB')
        writer.append_data(image_rgb)
        