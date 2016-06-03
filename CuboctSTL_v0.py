from stl import mesh
import math
import numpy as np
from matplotlib import pyplot
from mpl_toolkits import mplot3d


def make_closed_lattice(strut_width, chamfer_factor, pitch, x, y, z,  ):
    '''
    This function creates a closed cuboct lattice.
    :param strut_width: float lattice strut width
    :param chamfer_factor: float node chamfer factor. Lower number corresponds to more node reinforcement
    :param pitch: float lattice pitch (distance between voxels)
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param z: integer number of items in the lattice in z direction
    :return: numpy stl mesh object of closed cuboct lattice
    '''

    # Make the voxel to be arrayed
    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    # Array the voxel into a lattice
    one_lattice = box_array(one_voxel, pitch, x, y, z)
    # Cap the open sides of the lattice
    closed_lattice = box_cap(one_lattice, strut_width, chamfer_factor, pitch, x, y, z)

    return closed_lattice


def voxel(strut_width, chamfer_factor, pitch):
    '''
    Creates the mesh of an open cuboct voxel.
    :param strut_width: float
    :param chamfer_factor: float
    :param pitch: float
    :return: numpy stl mesh object of voxel
    '''
    # Define list of voxel nodes
    vnodes = [
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor)
    ]
    # Place and orient all the nodes (can probably be done more efficiently)
    vnodes[1].rotate([0.5, 0, 0], math.radians(180))
    translate(vnodes[1], np.array([0, 0, 1]) * pitch)  # top node
    vnodes[2].rotate([0.5, 0, 0], math.radians(90))
    translate(vnodes[2], np.array([0, -0.5, 0.5]) * pitch)
    vnodes[3].rotate([0.5, 0, 0], math.radians(270))
    translate(vnodes[3], np.array([0, 0.5, 0.5]) * pitch)
    vnodes[4].rotate([0, 0.5, 0], math.radians(90))
    translate(vnodes[4], np.array([0.5, 0, 0.5]) * pitch)
    vnodes[5].rotate([0, 0.5, 0], math.radians(270))
    translate(vnodes[5], np.array([-0.5, 0, 0.5]) * pitch)

    # Define voxel struts using strut function

    strut1 = strut(strut_width, chamfer_factor, pitch)
    bottom_struts = arraypolar(strut1, [0, 0, 1], 4)
    strut1.rotate([1, 0, 0], math.radians(180))
    translate(strut1, np.array([0, 0, 1]) * pitch)
    top_struts = arraypolar(strut1, [0, 0, 1], 4)

    strut2 = side_strut(strut_width, chamfer_factor, pitch)
    side_struts = arraypolar(strut2, [0, 0, 0.5], 4)

    combined_geometry = bottom_struts + top_struts + side_struts +vnodes

    return combine_meshes(*combined_geometry)


def closed_voxel(strut_width, chamfer_factor, pitch):
    '''
    Creates the mesh of a closed cuboct voxel (i.e. with capped nodes).
    :param strut_width: float
    :param chamfer_factor: float
    :param pitch: float
    :return: numpy stl mesh object of voxel
    '''

    # Define list of voxel nodes
    vnodes = [
        capped_node(strut_width, chamfer_factor),
        capped_node(strut_width, chamfer_factor),
        capped_node(strut_width, chamfer_factor),
        capped_node(strut_width, chamfer_factor),
        capped_node(strut_width, chamfer_factor),
        capped_node(strut_width, chamfer_factor)
    ]
    # Place and orient all the nodes (can probably be done more efficiently)
    vnodes[1].rotate([0.5, 0, 0], math.radians(180))
    translate(vnodes[1], np.array([0, 0, 1]) * pitch)  # top node
    vnodes[2].rotate([0.5, 0, 0], math.radians(90))
    translate(vnodes[2], np.array([0, -0.5, 0.5]) * pitch)
    vnodes[3].rotate([0.5, 0, 0], math.radians(270))
    translate(vnodes[3], np.array([0, 0.5, 0.5]) * pitch)
    vnodes[4].rotate([0, 0.5, 0], math.radians(90))
    translate(vnodes[4], np.array([0.5, 0, 0.5]) * pitch)
    vnodes[5].rotate([0, 0.5, 0], math.radians(270))
    translate(vnodes[5], np.array([-0.5, 0, 0.5]) * pitch)

    # Define voxel struts using strut function
    strut0 = strut(strut_width, chamfer_factor, pitch)
    strut1 = strut(strut_width, chamfer_factor, pitch)
    bottom_struts = arraypolar(strut1, [0, 0, 0.5], 4)
    strut1.rotate([1, 0, 0], math.radians(180))
    translate(strut1, np.array([0, 0, 1]) * pitch)
    top_struts = arraypolar(strut1, [0, 0, 0.5], 4)

    strut2 = side_strut(strut_width, chamfer_factor, pitch)
    side_struts = arraypolar(strut2, [0, 0, 0.5], 4)

    combined_geometry = bottom_struts + top_struts + side_struts +vnodes

    return combine_meshes(*combined_geometry)


def cap(strutwidth, chamfactor):
    '''
    This function generates the octahedral cap geometry for cuboct. The outward normal is facing in negative z-direction
    :param strutwidth: float
    :param chamfactor: float
    :return: numpy stl mesh object of cuboct cap
    '''

    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2
    l_2 = strutwidth / 2 + chamheight
    l_3 = l_2 + strutwidth * np.cos(np.pi / 4)  # horizontal position of points

    # Make the bottom cap (to make closed node in voxel)
    bottomcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    point1b = [halfw, l_3, 0]
    point2b = [l_3, halfw, 0]
    point8b = [-halfw, l_3, 0]
    point3b = [l_3, -halfw, 0]
    point4b = [halfw, -l_3, 0]
    point5b = [-halfw, -l_3, 0]
    point6b = [-l_3, -halfw, 0]
    point7b = [-l_3, halfw, 0]
    bottomcap['vectors'][0] = np.array([point8b, point1b, point7b])
    bottomcap['vectors'][1] = np.array([point1b, point2b, point7b])
    bottomcap['vectors'][2] = np.array([point7b, point2b, point6b])
    bottomcap['vectors'][3] = np.array([point2b, point3b, point6b])
    bottomcap['vectors'][4] = np.array([point3b, point5b, point6b])
    bottomcap['vectors'][5] = np.array([point3b, point4b, point5b])
    bottom = mesh.Mesh(bottomcap)

    return bottom

def box_cap(open_lattice, strutwidth, chamfactor, pitch, x, y, z):
    '''
    This function applies caps to the outside of an open lattice mesh to create a closed mesh suitable for printing
    :param open_lattice: lattice mesh object
    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :param x: integer number of voxels in the lattice in x direction
    :param y: integer number of voxels in the lattice in y direction
    :param z: integer number of voxels in the lattice in z direction
    :return: numpy stl mesh object of closed lattice
    '''

    closed_lattice = [open_lattice]  # Assume want list structure

    # Generate the cap geometry
    cap_geo = cap(strutwidth, chamfactor)

    # ------cap bottom---------
    bottom_caps = rec_array(cap_geo, x, y, [1, 0, 0], [0, 1, 0], pitch, pitch)
    closed_lattice += bottom_caps

    # ------cap top---------
    # flip the cap geometry to array the top
    # NOTE for FUTURE WORK: it might be more efficient long term to just flip the normals and translate the whole plane
    cap_geo_top = mesh.Mesh(cap_geo.data.copy())
    cap_geo_top.rotate([1, 0, 0], math.radians(180))
    translate(cap_geo_top, np.array([0, 0, 1])*pitch*z)
    top_caps = rec_array(cap_geo_top, x, y, [1, 0, 0], [0, 1, 0], pitch, pitch)
    closed_lattice += top_caps # rec_array returns a list, so this works

    # ------cap negX (left) -----------
    # rotate and translate cap geometry
    cap_geo_left = mesh.Mesh(cap_geo.data.copy())
    cap_geo_left.rotate([0, 1, 0], math.radians(270))
    translate(cap_geo_left, np.array([-1, 0, 0]) * pitch / 2)
    translate(cap_geo_left, np.array([0, 0, 1]) * pitch / 2)
    left_side_caps = rec_array(cap_geo_left, y, z, [0, 1, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += left_side_caps

    # ------cap posX (right) -----------
    cap_geo_right = mesh.Mesh(cap_geo.data.copy())
    cap_geo_right.rotate([0, 1, 0], math.radians(90))
    translate(cap_geo_right, np.array([1, 0, 0]) * pitch * x)
    translate(cap_geo_right, np.array([0, 0, 1]) * pitch / 2)
    translate(cap_geo_right, np.array([-1, 0, 0]) * pitch / 2)
    right_side_caps = rec_array(cap_geo_right, y, z, [0, 1, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += right_side_caps

    # --------cap front (negY) ------------
    cap_geo_front = mesh.Mesh(cap_geo.data.copy())
    cap_geo_front.rotate([1, 0, 0], math.radians(90))
    translate(cap_geo_front, np.array([0, -1, 0]) * pitch / 2)
    translate(cap_geo_front, np.array([0, 0, 1]) * pitch / 2)
    front_caps = rec_array(cap_geo_front, x, z, [1, 0, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += front_caps

    # -------cap back (posY) --------------
    cap_geo_back = mesh.Mesh(cap_geo.data.copy())
    cap_geo_back.rotate([1, 0, 0], math.radians(270))
    translate(cap_geo_back, np.array([0, 1, 0]) * pitch * y)
    translate(cap_geo_back, np.array([0, -1, 0]) * pitch / 2)
    translate(cap_geo_back, np.array([0, 0, 1]) * pitch / 2)
    back_caps = rec_array(cap_geo_back, x, z, [1, 0, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += back_caps

    return combine_meshes(*closed_lattice)

def rec_array(mesh_object, x, y, x_vector, y_vector, x_pitch, y_pitch):
    '''
    This function arrays a given mesh object in two dimensions x and y.
    Note that x and y in this function need not be globally defined x and y.
    NOTE: returns a list of the arrayed mesh objects, not a unified mesh object
    :param mesh_object: numpy stl mesh object
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param x_vector: vector  ex. [1, 0 , 0]
    :param y_vector: vector  ex. [0, 0 , 1]
    :param x_pitch: pitch in x direction (distance between arrayed objects)
    :param y_pitch: pitch in y direction (distance between arrayed objects)
    :return: rectangular_array: list containing all arrayed objects
    '''

    rectangular_array = [mesh_object]

    if x == 1:  # Don't need to do anything if x dimension is 1
        x = 1
    else:
        for i in range(x):
            if i == 0:  # Don't need to make a copy of the original voxel
                i = 0
            else:
                new_obj = mesh.Mesh(mesh_object.data.copy())  # Make a copy of the voxel
                translate(new_obj, np.array(x_vector) * x_pitch * i)
                rectangular_array += [new_obj]
    # array in y direction
    if y == 1:  # Don't need to do anything if y dimension is 1
        y = 1
    else:
        xline = rectangular_array
        for j in range(y):
            if j == 0:  # Don't need to make a copy of the original voxel line
                j = 0
            else:
                for thing in xline:
                    new_obj = mesh.Mesh(thing.data.copy())  # Make a copy of the voxel
                    translate(new_obj, np.array(y_vector) * y_pitch * j)
                    rectangular_array = rectangular_array + [
                        new_obj]  # Can't use += because modifies copy too and creates an infinite loop

    # return the list of arrayed mesh objects
    return rectangular_array


def node(strutwidth, chamfactor):
    '''
    This function creates a mesh of an open cuboct node.
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of open node
    '''

    # Calculate commonly used values for geometry definition
    chamheight = strutwidth/ chamfactor
    halfw = strutwidth / 2

    # Define geometry of the top cap
    topcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    # topcap['vectors'][0] = np.array([[halfw, l_2, h],
                                       #  [l_2, halfw, h],
                                        #[0, 0, h]])
    # Calculate the height and halflength of top octogonal cap
    h = chamheight + (strutwidth * np.sin(np.pi / 4) + strutwidth / 2)
    l_2 = strutwidth / 2 + chamheight
    point1 = [halfw, l_2, h]
    point2 = [l_2, halfw, h]
    point3 = [l_2, -halfw, h]
    point4 = [halfw, -l_2, h]
    point5 = [-halfw, -l_2, h]
    point6 = [-l_2, -halfw, h]
    point7 = [-l_2, halfw, h]
    point8 = [-halfw, l_2, h]
    # Define with right hand rule to assure outward facing normal
    topcap['vectors'][0] = np.array([point8, point7, point1])
    topcap['vectors'][1] = np.array([point1, point7, point2])
    topcap['vectors'][2] = np.array([point2, point7, point6])
    topcap['vectors'][3] = np.array([point2, point6, point3])
    topcap['vectors'][4] = np.array([point3, point6, point5])
    topcap['vectors'][5] = np.array([point3, point5, point4])
    top = mesh.Mesh(topcap)

    # Define Geometry of the chamfered sides
    chamfersides = np.zeros(8, dtype=mesh.Mesh.dtype)

    hs = l_2 # height of side points of node
    l_3 = l_2 + strutwidth*np.cos(np.pi/4) # horizontal position of points
    point1s = [halfw, l_3, hs]
    point2s = [l_3, halfw, hs]
    point3s = [l_3, -halfw, hs]
    point4s = [halfw, -l_3, hs]
    point5s = [-halfw, -l_3, hs]
    point6s = [-l_3, -halfw, hs]
    point7s = [-l_3, halfw, hs]
    point8s = [-halfw, l_3, hs]
    s_point = [
        [halfw, l_3, hs],
        [l_3, halfw, hs],
        [l_3, -halfw, hs],
        [halfw, -l_3, hs],
        [-halfw, -l_3, hs],
        [-l_3, -halfw, hs],
        [-l_3, halfw, hs],
        [-halfw, l_3, hs]
    ]
    chamfersides['vectors'][0] = np. array([point1, point2s, point1s])
    chamfersides['vectors'][1] = np. array([point1, point2, point2s])
    chamfersides['vectors'][2] = np. array([point3, point4s, point3s])
    chamfersides['vectors'][3] = np. array([point3, point4, point4s])
    chamfersides['vectors'][4] = np. array([point5, point6s, point5s])
    chamfersides['vectors'][5] = np. array([point5, point6, point6s])
    chamfersides['vectors'][6] = np. array([point7, point8s, point7s])
    chamfersides['vectors'][7] = np. array([point7, point8, point8s])
    chamfersides_mesh = mesh.Mesh(chamfersides)


    # Define the rectangular sides
    sides = np.zeros(4, dtype=mesh.Mesh.dtype)

    point1b = [halfw, l_3, 0]
    point2b = [l_3, halfw,0]
    point8b = [-halfw, l_3, 0]
    sides['vectors'][0] = np. array([point1s, point2b, point1b])
    sides['vectors'][1] = np. array([point1s, point2s, point2b])
    sides['vectors'][2] = np. array([point8s, point1b, point8b])
    sides['vectors'][3] = np. array([point8s, point1s, point1b])
    sidesubmesh1 = mesh.Mesh(sides.copy())
    sidesubmesh2 = mesh.Mesh(sides.copy())
    sidesubmesh3 = mesh.Mesh(sides.copy())
    sidesubmesh4 = mesh.Mesh(sides.copy())
    sidesubmesh2.rotate([0.0, 0.0, 0.5], math.radians(90))
    sidesubmesh3.rotate([0.0, 0.0, 0.5], math.radians(180))
    sidesubmesh4.rotate([0.0, 0.0, 0.5], math.radians(270))

    # Make final mesh for the open node geometry
    finalnodemesh = mesh.Mesh(np.concatenate([
        chamfersides_mesh.data.copy(),
        top.data.copy(),
        sidesubmesh1.data.copy(),
        sidesubmesh2.data.copy(),
        sidesubmesh3.data.copy(),
        sidesubmesh4.data.copy(),
    ]))
    return finalnodemesh


def capped_node(strutwidth, chamfactor):
    '''
    This function creates a mesh of a closed cuboct node (i.e. has bottom cap).
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of closed node
    '''

    # Calculate commonly used values for geometry definition
    chamheight = strutwidth/ chamfactor
    halfw = strutwidth / 2

    # Define geometry of the top cap
    topcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    # topcap['vectors'][0] = np.array([[halfw, l_2, h],
                                       #  [l_2, halfw, h],
                                        #[0, 0, h]])
    # Calculate the height and halflength of top octogonal cap
    h = chamheight + (strutwidth * np.sin(np.pi / 4) + strutwidth / 2)
    l_2 = strutwidth / 2 + chamheight
    point1 = [halfw, l_2, h]
    point2 = [l_2, halfw, h]
    point3 = [l_2, -halfw, h]
    point4 = [halfw, -l_2, h]
    point5 = [-halfw, -l_2, h]
    point6 = [-l_2, -halfw, h]
    point7 = [-l_2, halfw, h]
    point8 = [-halfw, l_2, h]

    # Define with right hand rule to assure outward facing normal
    topcap['vectors'][0] = np.array([point8, point7, point1])
    topcap['vectors'][1] = np.array([point1, point7, point2])
    topcap['vectors'][2] = np.array([point2, point7, point6])
    topcap['vectors'][3] = np.array([point2, point6, point3])
    topcap['vectors'][4] = np.array([point3, point6, point5])
    topcap['vectors'][5] = np.array([point3, point5, point4])
    top = mesh.Mesh(topcap)

    # Define Geometry of the chamfered sides
    chamfersides = np.zeros(8, dtype=mesh.Mesh.dtype)

    hs = l_2 # height of side points of node
    l_3 = l_2 + strutwidth*np.cos(np.pi/4) # horizontal position of points
    point1s = [halfw, l_3, hs]
    point2s = [l_3, halfw, hs]
    point3s = [l_3, -halfw, hs]
    point4s = [halfw, -l_3, hs]
    point5s = [-halfw, -l_3, hs]
    point6s = [-l_3, -halfw, hs]
    point7s = [-l_3, halfw, hs]
    point8s = [-halfw, l_3, hs]
    s_point = [
        [halfw, l_3, hs],
        [l_3, halfw, hs],
        [l_3, -halfw, hs],
        [halfw, -l_3, hs],
        [-halfw, -l_3, hs],
        [-l_3, -halfw, hs],
        [-l_3, halfw, hs],
        [-halfw, l_3, hs]
    ]
    chamfersides['vectors'][0] = np.array([point1, point2s, point1s])
    chamfersides['vectors'][1] = np.array([point1, point2, point2s])
    chamfersides['vectors'][2] = np.array([point3, point4s, point3s])
    chamfersides['vectors'][3] = np.array([point3, point4, point4s])
    chamfersides['vectors'][4] = np.array([point5, point6s, point5s])
    chamfersides['vectors'][5] = np.array([point5, point6, point6s])
    chamfersides['vectors'][6] = np.array([point7, point8s, point7s])
    chamfersides['vectors'][7] = np.array([point7, point8, point8s])
    chamfersides_mesh = mesh.Mesh(chamfersides)


    # Define the rectangular sides
    sides = np.zeros(4, dtype=mesh.Mesh.dtype)

    point1b = [halfw, l_3, 0]
    point2b = [l_3, halfw,0]
    point8b = [-halfw, l_3, 0]
    sides['vectors'][0] = np.array([point1s, point2b, point1b])
    sides['vectors'][1] = np.array([point1s, point2s, point2b])
    sides['vectors'][2] = np.array([point8s, point1b, point8b])
    sides['vectors'][3] = np.array([point8s, point1s, point1b])

    sidesubmesh1 = mesh.Mesh(sides.copy())
    sidesubmesh2 = mesh.Mesh(sides.copy())
    sidesubmesh3 = mesh.Mesh(sides.copy())
    sidesubmesh4 = mesh.Mesh(sides.copy())
    sidesubmesh2.rotate([0.0, 0.0, 0.5], math.radians(90))
    sidesubmesh3.rotate([0.0, 0.0, 0.5], math.radians(180))
    sidesubmesh4.rotate([0.0, 0.0, 0.5], math.radians(270))

    # Make the bottom cap (to make closed node in voxel)
    bottomcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    point3b = [l_3, -halfw, 0]
    point4b = [halfw, -l_3, 0]
    point5b = [-halfw, -l_3, 0]
    point6b = [-l_3, -halfw, 0]
    point7b = [-l_3, halfw, 0]
    bottomcap['vectors'][0] = np.array([point8b, point1b, point7b])
    bottomcap['vectors'][1] = np.array([point1b, point2b, point7b])
    bottomcap['vectors'][2] = np.array([point7b, point2b, point6b])
    bottomcap['vectors'][3] = np.array([point2b, point3b, point6b])
    bottomcap['vectors'][4] = np.array([point3b, point5b, point6b])
    bottomcap['vectors'][5] = np.array([point3b, point4b, point5b])
    bottom = mesh.Mesh(bottomcap)

    # Make final mesh for the open node geometry
    finalnodemesh = mesh.Mesh(np.concatenate([
        chamfersides_mesh.data.copy(),
        top.data.copy(),
        sidesubmesh1.data.copy(),
        sidesubmesh2.data.copy(),
        sidesubmesh3.data.copy(),
        sidesubmesh4.data.copy(),
        bottom.data.copy()
    ]))
    return finalnodemesh


def strut(strutwidth, chamfactor,  pitch):
    '''
    This function creates the mesh of a cuboct strut. It corresponds to the strut on the bottom half of the cuboct voxel
    that would project onto the positive x axis.

    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :return: numpy mesh object of strut
    '''


    # Define connection points on bottom node
    # Geometry Parameters
    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2
    halfp = pitch / 2
    h = chamheight + (strutwidth * np.sin(np.pi / 4) + strutwidth / 2) # height of top cap
    l_2 = strutwidth / 2 + chamheight # horizontal position of points on topcap
    hs = l_2  # height of side points of node
    l_3 = l_2 + strutwidth * np.cos(np.pi / 4)  # horizontal position of points

    point2_copy = [l_2, halfw, h]
    point3_copy = [l_2, -halfw, h]
    point2s_copy = [l_3, halfw, hs]
    point3s_copy = [l_3, -halfw, hs]
    # new points to attach to on side node
    point2n = [halfp-h, halfw, halfp - l_2]
    point3n = [halfp-h, -halfw, halfp - l_2]
    point2sn = [halfp-hs, halfw, halfp - l_3]
    point3sn = [halfp-hs, -halfw, halfp - l_3]

    singlestrut_geo = np.zeros(8, dtype=mesh.Mesh.dtype)
    # This version of the strut definition attempts to fix the mesh normals

    singlestrut_geo['vectors'][0] = np.array([point2_copy, point2n, point2sn])
    singlestrut_geo['vectors'][1] = np.array([point2_copy, point2sn, point2s_copy])
    singlestrut_geo['vectors'][2] = np.array([point3_copy, point3sn,  point3n])
    singlestrut_geo['vectors'][3] = np.array([point3_copy, point3s_copy, point3sn])
    singlestrut_geo['vectors'][4] = np.array([point2_copy, point3n, point2n])
    singlestrut_geo['vectors'][5] = np.array([point2_copy, point3_copy, point3n])
    singlestrut_geo['vectors'][6] = np.array([point2s_copy, point2sn, point3sn])
    singlestrut_geo['vectors'][7] = np.array([point2s_copy, point3sn, point3s_copy])

    finalsinglestrut = mesh.Mesh(singlestrut_geo)
    return finalsinglestrut


def side_strut(strutwidth, chamfactor,  pitch):
    '''
        This function creates the mesh of a cuboct side strut. It corresponds to the strut on the horizontal center
        plane of the cuboct voxel. This function was created because meshing the rotated bottom strut into a unified
        mesh wasn't working for some reason (still unknown).
        :param strutwidth: float
        :param chamfactor: float
        :param pitch: float
        :return: numpy mesh object of side strut
        '''

    # Define connection points on bottom node
    # Geometry Parameters
    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2
    halfp = pitch / 2
    h = chamheight + (strutwidth * np.sin(np.pi / 4) + strutwidth / 2) # height of top cap
    l_2 = strutwidth / 2 + chamheight # horizontal position of points on topcap
    hs = l_2  # height of side points of node
    l_3 = l_2 + strutwidth * np.cos(np.pi / 4)  # horizontal position of points

    # translate the points we need down half a pitch
    point2_copy = [l_2, halfw, h - pitch/2]
    point3_copy = [l_2, -halfw, h - pitch/2]
    point2s_copy = [l_3, halfw, hs - pitch/2]
    point3s_copy = [l_3, -halfw, hs - pitch/2]
    # new points to attach to on side node
    point2n = [halfp-h, halfw, halfp - l_2 - pitch/2]
    point3n = [halfp-h, -halfw, halfp - l_2 - pitch/2]
    point2sn = [halfp-hs, halfw, halfp - l_3 - pitch/2]
    point3sn = [halfp-hs, -halfw, halfp - l_3 - pitch/2]

    # rotate all those points 90 degrees about the x axis. [ x, y, z] --> [ x, -z, y]

    point2_rotated = [point2_copy[0], -point2_copy[2], point2_copy[1]]
    point3_rotated = [point3_copy[0], -point3_copy[2], point3_copy[1]]
    point2s_rotated = [point2s_copy[0], -point2s_copy[2], point2s_copy[1]]
    point3s_rotated = [point3s_copy[0], -point3s_copy[2], point3s_copy[1]]
    point2n_rotated = [point2n[0], -point2n[2], point2n[1]]
    point3n_rotated = [point3n[0], -point3n[2], point3n[1]]
    point2sn_rotated = [point2sn[0], -point2sn[2], point2sn[1]]
    point3sn_rotated = [point3sn[0], -point3sn[2], point3sn[1]]

    singlestrut_geo = np.zeros(8, dtype=mesh.Mesh.dtype)
    # This version of the strut definition attempts to fix the mesh normals

    singlestrut_geo['vectors'][0] = np.array([point2_rotated, point2n_rotated, point2sn_rotated])
    singlestrut_geo['vectors'][1] = np.array([point2_rotated, point2sn_rotated, point2s_rotated])
    singlestrut_geo['vectors'][2] = np.array([point3_rotated, point3sn_rotated,  point3n_rotated])
    singlestrut_geo['vectors'][3] = np.array([point3_rotated, point3s_rotated, point3sn_rotated])
    singlestrut_geo['vectors'][4] = np.array([point2_rotated, point3n_rotated, point2n_rotated])
    singlestrut_geo['vectors'][5] = np.array([point2_rotated, point3_rotated, point3n_rotated])
    singlestrut_geo['vectors'][6] = np.array([point2s_rotated, point2sn_rotated, point3sn_rotated])
    singlestrut_geo['vectors'][7] = np.array([point2s_rotated, point3sn_rotated, point3s_rotated])

    # Move the strut back to side height
    singlestrut_geo['vectors'] += [0, 0, pitch/2]

    finalsinglestrut = mesh.Mesh(singlestrut_geo)

    return finalsinglestrut


def translate(meshobj, tvect):
    '''
    -------function from Daniel Cellucci's latticegen code--------
    translates mesh objects
    :param meshobj:
    :param tvect: numpy array ex. np.array([0, 0, 1])
    :return:
    '''
    vects = meshobj.vectors
    for vsi,vectset in enumerate(vects):
        for vi,vector in enumerate(vectset):
            meshobj.vectors[vsi][vi] = vector+tvect


def place_object(mesh_object, x_trans, y_trans, z_trans):
    '''
    This function translates a mesh object in x, y, and z direction.
    :param mesh_object: mesh object to translate
    :param x_trans: float. distance to translate in x
    :param y_trans: float. distance to translate in y
    :param z_trans: float. distance to translate in z
    :return:
    '''

    translate(mesh_object, np.array([1, 0, 0]) * x_trans)
    translate(mesh_object, np.array([0, 1, 0]) * y_trans)
    translate(mesh_object, np.array([0, 0, 1]) * z_trans)


def rotation_matrix( axis, theta):
    '''
    ----------THIS FUNCTION TAKEN FROM NUMPY STL ---------
    Generate a rotation matrix to Rotate the matrix over the given axis by
    the given theta (angle)

    Uses the Euler-Rodrigues formula for fast rotations:
    `https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula`_

    :param numpy.array axis: Axis to rotate over (x, y, z)
    :param float theta: Rotation angle in radians, use `math.radians` to
    convert degrees to radians if needed.
    '''
    axis = np.asarray(axis)
    # No need to rotate if there is no actual rotation
    if not axis.any():
        return np.zeros((3, 3))

    theta = np.asarray(theta)
    theta /= 2.

    axis = axis / math.sqrt( np.dot( axis, axis))

    a = math.cos(theta)
    b, c, d = - axis * math.sin(theta)
    angles = a, b, c, d
    powers = [x * y for x in angles for y in angles]
    aa, ab, ac, ad = powers[0:4]
    ba, bb, bc, bd = powers[4:8]
    ca, cb, cc, cd = powers[8:12]
    da, db, dc, dd = powers[12:16]

    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])


def rotate(object, axis, theta, point=None):
    '''
    ----------THIS FUNCTION TAKEN FROM NUMPY STL ---------
    Rotate the matrix over the given axis by the given theta (angle)

    Uses the `rotation_matrix`_ in the background.

    :param numpy.array axis: Axis to rotate over (x, y, z)
    :param float theta: Rotation angle in radians, use `math.radians` to
    convert degrees to radians if needed.
    :param numpy.array point: Rotation point so manual translation is not
    required
    '''
    # No need to rotate if there is no actual rotation
    if not theta:
        return

    point = np.asarray(point or [0] * 3)
    rot_matrix = rotation_matrix(axis, theta)

    # No need to rotate if there is no actual rotation
    if not rotation_matrix.any():
        return

    def _rotate(matrix):
        if point.any():
            # Translate while rotating
            return (matrix + point).dot(rot_matrix) - point
        else:
            # Simply apply the rotation
            return matrix.dot(rot_matrix)

    for i in range(3):
        object.vectors[:, i] = _rotate(object.vectors[:, i])


def arraypolar(m_obj, r_axis, num, rotation_point=None):
    '''
    This function arrays mesh objects in evenly spaced circular pattern.
    :param m_obj: numpy stl mesh object to array
    :param r_axis: rotation axis ex. [0, 0, 1]
    :param num: number of items to array
    :param rotation_point: point to place rotation axis if not center
    :return: list of arrayed objects
    '''
    array_objects = list()
    for i in range(num):
        obj = mesh.Mesh(m_obj.data.copy())
        obj.rotate(r_axis, math.radians((360 / num) * i), rotation_point)
        array_objects.append(obj)
    return array_objects


def combine_meshes(*args):
    '''
    This function combines a list or lists of mesh objects into a single mesh object
    :param args: list of mesh objects
    :return: numpy stl mesh object of combined geometries
    '''
    combined_data = np.concatenate([m_obj.data for m_obj in args])
    return mesh.Mesh(combined_data, remove_duplicate_polygons = True)


def box_array(voxel_mesh, pitch, x, y, z):
    '''
    This function cubically arrays a mesh object
    :param voxel_mesh:
    :param pitch:
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param z: integer number of items in the lattice in z direction
    :return: numpy stl mesh object of arrayed geometry
    '''
    lattice = [voxel_mesh]  # Assume want list structure

    # array in x direction
    if x == 1:  # Don't need to do anything if x dimension is 1
        x = 1
    else:
        for i in range(x):
            if i == 0:  # Don't need to make a copy of the original voxel
                i = 0
            else:
                new_obj = mesh.Mesh(voxel_mesh.data.copy())  # Make a copy of the voxel
                translate(new_obj, np.array([1, 0, 0]) * pitch* i)
                lattice += [new_obj]
    # array in y direction
    if y == 1:  # Don't need to do anything if y dimension is 1
        y = 1
    else:
        xline = lattice
        for j in range(y):
            if j == 0:  # Don't need to make a copy of the original voxel line
                j = 0
            else:
                for thing in xline:
                    new_obj = mesh.Mesh(thing.data.copy())  # Make a copy of the voxel
                    translate(new_obj, np.array([0, 1, 0]) * pitch * j)
                    lattice = lattice + [new_obj]  # Can't use += because modifies copy too and creates an infinite loop
    # array in z direction
    if z == 1:  # Don't need to do anything if the z dimension is 1
        z = 1
    else:
        xyplane = lattice
        for k in range(z):
            if k == 0:  # Don't need to make a copy of the original voxel plane
                k == 0
            else:
                for thing in xyplane:
                    new_obj = mesh.Mesh(thing.data.copy())  # Make a copy of the voxel
                    translate(new_obj, np.array([0, 0, 1]) * pitch * k)
                    lattice = lattice + [new_obj]  # Can't use += because modifies copy too and creates an infinite loop
    open_lattice = combine_meshes(*lattice)

    return open_lattice


def lattice_codedstructure(voxel_mesh, cap_mesh, pitch, template, closed=True):
    '''
    This function creates a lattice structure with individual voxel placement prescribed by a structure template.
    :param voxel_mesh: numpy stl mesh object of voxel geometry to be arrayed
    :param cap_mesh: numpy stl mesh object of the voxel cap
    :param pitch: float. lattice pitch
    :param template: three-dimensional numpy array containing a 1 for voxel, 0 for no voxel in that location
    first dimension is x, second dimension is y, third dimension is z
    :param closed: boolean. Set to false for an open lattice
    :return: numpy stl mesh object of lattice structure defined by template
    '''

    lattice = []

    # Determine the x, y, and z size of the template (bounding box size in voxels)
    [x_size, y_size, z_size] = template.shape

    for i in range(x_size):
        for j in range(y_size):
            for k in range (z_size):
                if template[i, j, k] == 1:  # If a voxel is supposed to be placed
                    new_obj = mesh.Mesh(voxel_mesh.data.copy())  # Make a copy of the voxel
                    # Move the new voxel to the correct place
                    place_object(new_obj, pitch*i, pitch*j, pitch*k)
                    #translate(new_obj, np.array([1, 0, 0]) * pitch * i)
                    #translate(new_obj, np.array([0, 1, 0]) * pitch * j)
                    #translate(new_obj, np.array([0, 0, 1]) * pitch * k)
                    lattice += [new_obj]
                    # check the connectivity and do optional capping
                    flag = 0  # Create a flag to detect minimum connectivity

                    # ----------check top--------
                    if k+1< z_size: # if you aren't on the edge of the template
                        if template[i, j, k + 1] == 0: # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the top
                                cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so that normal vector is correct
                                place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k+1))
                                lattice += [cap_geo_top]
                        elif template[i, j, k+1] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so that normal vector is correct
                            place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k+1))
                            lattice += [cap_geo_top]

                    # --------check bottom---------
                    if k - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j, k - 1] == 0:  # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the top
                                cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                                place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                                lattice += [cap_geo_bottom]
                        elif template[i, j, k + 1] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                            place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                            lattice += [cap_geo_bottom]

                    # ----------check right (positive X)----------
                    if i + 1 < x_size:  # if you aren't on the edge of the template
                        if template[i + 1, j, k ] == 0:  # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the right
                                cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_right.rotate([0, 1, 0], math.radians(90))
                                place_object(cap_geo_right, pitch * i + pitch/2, pitch * j, pitch * k + pitch/2 )
                                lattice += [cap_geo_right]
                        elif template[i + 1, j, k] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_right.rotate([0, 1, 0], math.radians(90))
                            place_object(cap_geo_right, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_right]

                    # ----------check left (negative X)---------
                    if i - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i - 1, j, k] == 0:  # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the right
                                cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_left.rotate([0, 1, 0], math.radians(270))
                                place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap_geo_left]
                        elif template[i - 1, j, k] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_left.rotate([0, 1, 0], math.radians(270))
                            place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_left]

                    # ----------check back (positive Y)---------
                    if j + 1 < y_size:  # if you aren't on the edge of the template
                        if template[i, j+1, k] == 0:  # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the right
                                cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_back.rotate([1, 0, 0], math.radians(270))
                                place_object(cap_geo_back, pitch * i, pitch * j + pitch/2, pitch * k + pitch / 2)
                                lattice += [cap_geo_back]
                        elif template[i, j+1, k] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_back.rotate([1, 0, 0], math.radians(270))
                            place_object(cap_geo_back, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_back]

                    # check front (negative Y)
                    if j - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j-1, k] == 0:  # if there isn't a voxel there
                            if closed == True:
                                # place a cap on the right
                                cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_front.rotate([1, 0, 0], math.radians(90))
                                place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap_geo_front]
                        elif template[i, j-1, k] == 1:
                            flag = 1
                        else:
                            print(" Something went horribly wrong with the template.")
                    else:
                        if closed == True:
                            # You are on the edge, so place a cap
                            cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_front.rotate([1, 0, 0], math.radians(90))
                            place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_front]

                    # If the flag wasn't thrown, this voxel doesn't have any connectivity. Show an error
                    if flag == 0:
                        print(" There is a voxel in your template with zero connectivity.")

    return combine_meshes(*lattice)


def create_test_template():
    '''
    Creates a template to test the lattice_codedstructure function.
    :return: template
    '''
    
    template = np.zeros((3, 3, 3))
    template[1, 1, 1] = 1
    template[1, 0, 1] = 1
    template[0, 1, 1] = 1
    template[1, 0, 0] = 1
    template[2, 1, 1] = 1
    template[2, 2, 1] = 1

    return template


def main():
    pitch = 30
    strut_width = 2
    chamfer_factor = 2.75
    x=3
    y=3
    z=3

    #lattice = make_closed_lattice(strut_width, chamfer_factor, pitch, x, y, z)

    #lattice.save('test_lattice.stl')

    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    new_voxel = mesh.Mesh(one_voxel.data.copy())
    place_object(new_voxel, pitch*2, pitch*2, pitch*2)

    template = create_test_template()
    capmesh = cap(strut_width, chamfer_factor)
    structure = lattice_codedstructure(one_voxel, capmesh, pitch, template)
    structure.save('coded_structure_test.stl')

    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(structure.vectors))


    # Auto scale to the mesh size
    scale = one_voxel.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()



main()
