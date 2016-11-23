from stl import mesh
import math
import numpy as np
from matplotlib import pyplot
from mpl_toolkits import mplot3d


def make_lattice(strut_width, chamfer_factor, pitch, x, y, z, closed=True):
    """
    This function creates a closed cuboct lattice.
    :param strut_width: float lattice strut width
    :param chamfer_factor: float node chamfer factor. Lower number corresponds to more node reinforcement
    :param pitch: float lattice pitch (distance between voxels)
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param z: integer number of items in the lattice in z direction
    :param closed: optional boolean parameter to determine whether to close the lattice. To make an open lattice with
    no caps, set to False
    :return: numpy stl mesh object of cuboct lattice
    """

    # Make the voxel to be arrayed
    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    # Array the voxel into a lattice
    one_lattice = box_array(one_voxel, pitch, x, y, z)

    if closed is True:
        # Cap the open sides of the lattice
        final_lattice = box_cap(one_lattice, strut_width, chamfer_factor, pitch, x, y, z)
    else:
        final_lattice = one_lattice

    return final_lattice


def voxel(strut_width, chamfer_factor, pitch):
    """
    Creates the mesh of an open cuboct voxel.
    :param strut_width: float
    :param chamfer_factor: float
    :param pitch: float
    :return: numpy stl mesh object of voxel
    """
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

def half_voxel(strut_width, chamfer_factor, pitch):
    """
    This code creates half-voxel geometry that is closed on the half-surface. Note that the half geometry created with
    function includes complete mid-plane voxels and struts. It removes the top node and top struts. Therefore, the top
    surface is NOT a flat surface.
    :param strut_width:
    :param chamfer_factor:
    :param pitch:
    :return:
    """


    # Define list of voxel nodes
    vnodes = [
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor),
        node(strut_width, chamfer_factor)
    ]
    # Place and orient all the nodes (can probably be done more efficiently)
    vnodes[1].rotate([0, 0.5, 0], math.radians(270))
    translate(vnodes[1], np.array([-0.5, 0, 0.5]) * pitch)
    vnodes[2].rotate([0.5, 0, 0], math.radians(90))
    translate(vnodes[2], np.array([0, -0.5, 0.5]) * pitch)
    vnodes[3].rotate([0.5, 0, 0], math.radians(270))
    translate(vnodes[3], np.array([0, 0.5, 0.5]) * pitch)
    vnodes[4].rotate([0, 0.5, 0], math.radians(90))
    translate(vnodes[4], np.array([0.5, 0, 0.5]) * pitch)

    # Define half voxel struts using strut function
    strut1 = strut(strut_width, chamfer_factor, pitch)
    bottom_struts = arraypolar(strut1, [0, 0, 1], 4)
    strut2 = side_strut(strut_width, chamfer_factor, pitch)
    side_struts = arraypolar(strut2, [0, 0, 0.5], 4)
    #  May want to alter code so that there is a flat surface on the half-voxel surface


    # Define connection points on bottom node
    # Geometry Parameters
    # Calculate commonly used values for geometry definition
    chamheight = strut_width / chamfer_factor
    half_w = strut_width / 2
    halfp = pitch / 2
    h = chamheight + (strut_width * np.sin(np.pi / 4) + strut_width / 2)  # height of top cap
    l_2 = strut_width / 2 + chamheight  # horizontal position of points on topcap
    hs = l_2  # height of side points of node
    l_3 = l_2 + strut_width * np.cos(np.pi / 4)  # horizontal position of points

    # new points to attach to on side node
    point2nc = [halfp - h, half_w, halfp + l_2]
    point3nc = [halfp - h, -half_w, halfp + l_2]
    point2snc = [halfp - hs, half_w, halfp +l_3]
    point3snc = [halfp - hs, -half_w, halfp + l_3]

    strutcap_geo = np.zeros(2, dtype=mesh.Mesh.dtype)

    strutcap_geo['vectors'][0] = np.array([point3nc, point3snc, point2snc])
    strutcap_geo['vectors'][1] = np.array([point2nc, point3nc, point2snc])

    strut_cap_geo = mesh.Mesh(strutcap_geo)

    #translate(strut_cap_geo, np.array([0.5, 0, 0.5])*pitch)
    #strut_cap_geo.rotate([0, 1, 0], math.radians(90), [ -pitch*0.5, 0, -pitch*0.5])  #unclear why making neg pitch works
    strut_caps = arraypolar(strut_cap_geo, [0, 0, 1], 4)

    combined_geometry = bottom_struts + side_struts + vnodes + strut_caps

    return combine_meshes(*combined_geometry)


def hybrid_voxel(strut_width, chamfer_factor, pitch, max_strut_width_interface):

    # Define list of voxel nodes
    vnodes = [
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface),
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface),
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface),
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface),
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface),
        hybrid_node(strut_width, chamfer_factor, max_strut_width_interface)
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

    combined_geometry = bottom_struts + top_struts + side_struts + vnodes

    return combine_meshes(*combined_geometry)

def closed_voxel(strut_width, chamfer_factor, pitch):
    """
    Creates the mesh of a closed cuboct voxel (i.e. with capped nodes).
    :param strut_width: float
    :param chamfer_factor: float
    :param pitch: float
    :return: numpy stl mesh object of voxel
    """

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



def cap_cuboct(strutwidth, chamfactor):
    """
    This function generates the octahedral cap geometry for cuboct. The outward normal is facing in negative z-direction
    :param strutwidth: float
    :param chamfactor: float
    :return: numpy stl mesh object of cuboct cap
    """

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
    """
    This function applies caps to the outside of an open lattice mesh to create a closed mesh suitable for printing
    :param open_lattice: lattice mesh object
    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :param x: integer number of voxels in the lattice in x direction
    :param y: integer number of voxels in the lattice in y direction
    :param z: integer number of voxels in the lattice in z direction
    :return: numpy stl mesh object of closed lattice
    """

    closed_lattice = [open_lattice]  # Assume want list structure

    # Generate the cap geometry
    cap_geo = cap_cuboct(strutwidth, chamfactor)

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

def box_cap_sides_only(open_lattice, strutwidth, chamfactor, pitch, x, y, z):
    """
    This function applies caps to the outside of an open lattice mesh to create a closed mesh suitable for printing
    :param open_lattice: lattice mesh object
    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :param x: integer number of voxels in the lattice in x direction
    :param y: integer number of voxels in the lattice in y direction
    :param z: integer number of voxels in the lattice in z direction
    :return: numpy stl mesh object of closed lattice
    """

    closed_lattice = [open_lattice]  # Assume want list structure

    # Generate the cap geometry
    cap_geo = cap_cuboct(strutwidth, chamfactor)

    # ------cap negX (left) -----------
    # rotate and translate cap geometry
    cap_geo_left = mesh.Mesh(cap_geo.data.copy())
    cap_geo_left.rotate([0, 1, 0], math.radians(270))
    translate(cap_geo_left, np.array([-1, 0, 0]) * pitch / 2)
    #translate(cap_geo_left, np.array([0, 0, 1]) * pitch / 2)
    left_side_caps = rec_array(cap_geo_left, y, z, [0, 1, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += left_side_caps

    # ------cap posX (right) -----------
    cap_geo_right = mesh.Mesh(cap_geo.data.copy())
    cap_geo_right.rotate([0, 1, 0], math.radians(90))
    translate(cap_geo_right, np.array([1, 0, 0]) * pitch * x)
    #translate(cap_geo_right, np.array([0, 0, 1]) * pitch / 2)
    translate(cap_geo_right, np.array([-1, 0, 0]) * pitch / 2)
    right_side_caps = rec_array(cap_geo_right, y, z, [0, 1, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += right_side_caps

    # --------cap front (negY) ------------
    cap_geo_front = mesh.Mesh(cap_geo.data.copy())
    cap_geo_front.rotate([1, 0, 0], math.radians(90))
    translate(cap_geo_front, np.array([0, -1, 0]) * pitch / 2)
    #translate(cap_geo_front, np.array([0, 0, 1]) * pitch / 2)
    front_caps = rec_array(cap_geo_front, x, z, [1, 0, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += front_caps

    # -------cap back (posY) --------------
    cap_geo_back = mesh.Mesh(cap_geo.data.copy())
    cap_geo_back.rotate([1, 0, 0], math.radians(270))
    translate(cap_geo_back, np.array([0, 1, 0]) * pitch * y)
    translate(cap_geo_back, np.array([0, -1, 0]) * pitch / 2)
    #translate(cap_geo_back, np.array([0, 0, 1]) * pitch / 2)
    back_caps = rec_array(cap_geo_back, x, z, [1, 0, 0], [0, 0, 1], pitch, pitch)
    closed_lattice += back_caps

    return combine_meshes(*closed_lattice)

def rec_array(mesh_object, x, y, x_vector, y_vector, x_pitch, y_pitch):
    """
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
    """

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
    """
    This function creates a mesh of an open cuboct node.
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of open node
    """

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
    """
    This function creates a mesh of a closed cuboct node (i.e. has bottom cap).
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of closed node
    """

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


def hybrid_node(strutwidth, chamfactor, max_strut_width_interface, closed=False):
    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2

    # Define geometry of the top cap
    topcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    # topcap['vectors'][0] = np.array([[halfw, l_2, h],
    #  [l_2, halfw, h],
    # [0, 0, h]])
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

    hs = l_2  # height of side points of node
    l_3 = l_2 + strutwidth * np.cos(np.pi / 4)  # horizontal position of points
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

    # Need to define geometry of base necessary to interface with other node
    chamheight_max = max_strut_width_interface / chamfactor
    halfw_max = max_strut_width_interface / 2
    l_2_max = max_strut_width_interface / 2 + chamheight_max
    l_3_max = l_2_max + max_strut_width_interface * np.cos(np.pi / 4)  # horizontal position of points



    # Define the rectangular sides
    sides = np.zeros(4, dtype=mesh.Mesh.dtype)

    point1b = [halfw_max, l_3_max, 0]
    point2b = [l_3_max, halfw_max, 0]
    point8b = [-halfw_max, l_3_max, 0]
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
    """
    This function creates a mesh of a closed cuboct node (i.e. has bottom cap).
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of closed node
    """

    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2

    # Define geometry of the top cap
    topcap = np.zeros(6, dtype=mesh.Mesh.dtype)
    # topcap['vectors'][0] = np.array([[halfw, l_2, h],
    #  [l_2, halfw, h],
    # [0, 0, h]])
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

    hs = l_2  # height of side points of node
    l_3 = l_2 + strutwidth * np.cos(np.pi / 4)  # horizontal position of points
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
    point2b = [l_3, halfw, 0]
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
    """
    This function creates the mesh of a cuboct strut. It corresponds to the strut on the bottom half of the cuboct voxel
    that would project onto the positive x axis.
    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :return: numpy mesh object of strut
    """

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
    """
    This function creates the mesh of a cuboct side strut. It corresponds to the strut on the horizontal center
    plane of the cuboct voxel. This function was created because meshing the rotated bottom strut into a unified
    mesh wasn't working for some reason (still unknown).
    :param strutwidth: float
    :param chamfactor: float
    :param pitch: float
    :return: numpy-stl mesh object of side strut
    """

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
    """
    -------function from Daniel Cellucci's latticegen code--------
    translates mesh objects
    :param meshobj:
    :param tvect: numpy array ex. np.array([0, 0, 1])
    :return:
    """
    vects = meshobj.vectors
    for vsi,vectset in enumerate(vects):
        for vi,vector in enumerate(vectset):
            meshobj.vectors[vsi][vi] = vector+tvect


def place_object(mesh_object, x_trans, y_trans, z_trans):
    """
    This function translates a mesh object in x, y, and z direction.
    :param mesh_object: mesh object to translate
    :param x_trans: float. distance to translate in x
    :param y_trans: float. distance to translate in y
    :param z_trans: float. distance to translate in z
    :return:
    """

    translate(mesh_object, np.array([1, 0, 0]) * x_trans)
    translate(mesh_object, np.array([0, 1, 0]) * y_trans)
    translate(mesh_object, np.array([0, 0, 1]) * z_trans)


def rotation_matrix( axis, theta):
    """
    ----------THIS FUNCTION TAKEN FROM NUMPY STL ---------
    Generate a rotation matrix to Rotate the matrix over the given axis by
    the given theta (angle)

    Uses the Euler-Rodrigues formula for fast rotations:
    `https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula`_

    :param numpy.array axis: Axis to rotate over (x, y, z)
    :param float theta: Rotation angle in radians, use `math.radians` to
    convert degrees to radians if needed.
    """
    axis = np.asarray(axis)
    # No need to rotate if there is no actual rotation
    if not axis.any():
        return np.zeros((3, 3))

    theta = np.asarray(theta)
    theta /= 2.

    axis = axis / math.sqrt( np.dot( axis, axis) )

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
    """
    ----------THIS FUNCTION TAKEN FROM NUMPY STL ---------
    Rotate the matrix over the given axis by the given theta (angle)

    Uses the `rotation_matrix`_ in the background.

    :param numpy.array axis: Axis to rotate over (x, y, z)
    :param float theta: Rotation angle in radians, use `math.radians` to
    convert degrees to radians if needed.
    :param numpy.array point: Rotation point so manual translation is not
    required
    """
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
    """
    This function arrays mesh objects in evenly spaced circular pattern.
    :param m_obj: numpy stl mesh object to array
    :param r_axis: rotation axis ex. [0, 0, 1]
    :param num: number of items to array
    :param rotation_point: point to place rotation axis if not center
    :return: list of arrayed objects
    """
    array_objects = list()
    for i in range(num):
        obj = mesh.Mesh(m_obj.data.copy())
        obj.rotate(r_axis, math.radians((360 / num) * i), rotation_point)
        array_objects.append(obj)
    return array_objects


def combine_meshes(*args):
    """
    This function combines a list or lists of mesh objects into a single mesh object
    :param args: list of mesh objects
    :return: numpy stl mesh object of combined geometries
    """
    combined_data = np.concatenate([m_obj.data for m_obj in args])
    return mesh.Mesh(combined_data, remove_duplicate_polygons = True)


def box_array(voxel_mesh, pitch, x, y, z):
    """
    This function cubically arrays a mesh object
    :param voxel_mesh:
    :param pitch:
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param z: integer number of items in the lattice in z direction
    :return: numpy stl mesh object of arrayed geometry
    """
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
                translate(new_obj, np.array([1, 0, 0]) * pitch * i)
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
    """
    This function creates a lattice structure with individual voxel placement prescribed by a structure template.
    :param voxel_mesh: numpy stl mesh object of voxel geometry to be arrayed
    :param cap_mesh: numpy stl mesh object of the voxel cap
    :param pitch: float. lattice pitch
    :param template: three-dimensional numpy array containing a 1 for voxel, 0 for no voxel in that location
    first dimension is x, second dimension is y, third dimension is z
    :param closed: boolean. Set to false for an open lattice
    :return: numpy stl mesh object of lattice structure defined by template
    """

    lattice = []

    # Determine the x, y, and z size of the template (bounding box size in voxels)
    [x_size, y_size, z_size] = template.shape

    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if template[i, j, k] == 1:  # If a voxel is supposed to be placed
                    new_obj = mesh.Mesh(voxel_mesh.data.copy())  # Make a copy of the voxel
                    # Move the new voxel to the correct place
                    place_object(new_obj, pitch*i, pitch*j, pitch*k)
                    lattice += [new_obj]
                    # check the connectivity and do optional capping
                    flag = 0  # Create a flag to detect minimum connectivity

                    # ----------check top--------
                    if k+1 < z_size:  # if you aren't on the edge of the template
                        if template[i, j, k + 1] == 0: # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so normal vectors correct
                                place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k+1))
                                lattice += [cap_geo_top]
                        elif template[i, j, k+1] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check above voxel  x = {0} y = {1} z = {2}'.format(i, j, k))
                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so that normal vector is correct
                            place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k+1))
                            lattice += [cap_geo_top]

                    # --------check bottom---------
                    if k - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j, k - 1] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                                place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                                lattice += [cap_geo_bottom]
                        elif template[i, j, k - 1] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check below voxel  x = {0} y = {1} z = {2}'.format(i, j, k))
                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                            place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                            lattice += [cap_geo_bottom]

                    # ----------check right (positive X)----------
                    if i + 1 < x_size:  # if you aren't on the edge of the template
                        if template[i + 1, j, k ] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_right.rotate([0, 1, 0], math.radians(90))
                                place_object(cap_geo_right, pitch * i + pitch/2, pitch * j, pitch * k + pitch/2 )
                                lattice += [cap_geo_right]
                        elif template[i + 1, j, k] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check right of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))
                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_right.rotate([0, 1, 0], math.radians(90))
                            place_object(cap_geo_right, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_right]

                    # ----------check left (negative X)---------
                    if i - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i - 1, j, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_left.rotate([0, 1, 0], math.radians(270))
                                place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap_geo_left]
                        elif template[i - 1, j, k] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check left of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))
                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_left.rotate([0, 1, 0], math.radians(270))
                            place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_left]

                    # ----------check back (positive Y)---------
                    if j + 1 < y_size:  # if you aren't on the edge of the template
                        if template[i, j+1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_back.rotate([1, 0, 0], math.radians(270))
                                place_object(cap_geo_back, pitch * i, pitch * j + pitch/2, pitch * k + pitch / 2)
                                lattice += [cap_geo_back]
                        elif template[i, j+1, k] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check behind (pos y) voxel  x = {0} y = {1} z = {2}'.format(i, j, k))
                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_back.rotate([1, 0, 0], math.radians(270))
                            place_object(cap_geo_back, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_back]

                    # check front (negative Y)
                    if j - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j-1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_front.rotate([1, 0, 0], math.radians(90))
                                place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap_geo_front]
                        elif template[i, j-1, k] == 1:
                            flag = 1
                        else:
                            print('Template Error. Check front(neg y)of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_front.rotate([1, 0, 0], math.radians(90))
                            place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_front]

                    # If the flag wasn't thrown, this voxel doesn't have any connectivity. Show an error
                    if flag == 0:
                        print(" There is a voxel in your template with zero connectivity.")

    return combine_meshes(*lattice)

def compression_specimen(strut_width, chamfer_factor, pitch, x, y, z):
    """
    This function creates a closed cuboct lattice, with a half plane of half-voxels on the top and bottom.
    i.e there will be z-1 complete voxels in the specimen.
    Also note that due to the fact that the entire mid-plane voxel node and strut are included in the half-voxel geo,
    the overall height will be slightly over z*pitch.
    :param strut_width: float lattice strut width
    :param chamfer_factor: float node chamfer factor. Lower number corresponds to more node reinforcement
    :param pitch: float lattice pitch (distance between voxels)
    :param x: integer number of items in the lattice in x direction
    :param y: integer number of items in the lattice in y direction
    :param z: integer number of items in the lattice in z direction
    :return: numpy stl mesh object of cuboct lattice
    """

    # Make the voxel to be arrayed
    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    # Array the voxel into a lattice and translate up one half-pitch
    one_lattice = box_array(one_voxel, pitch, x, y, z-1)
    translate(one_lattice, np.array([0, 0, 0.5])*pitch)

    # Add the half-voxels to the top and bottom
    half_vox1 = half_voxel(strut_width, chamfer_factor, pitch)
    translate(half_vox1, np.array([0, 0, 1]) * pitch * (z - 0.5))
    top_half_plane = rec_array(half_vox1, x, y, [1, 0, 0], [0, 1, 0], pitch, pitch)

    half_vox2 = half_voxel(strut_width, chamfer_factor, pitch)
    half_vox2.rotate([1, 0, 0], math.radians(180))
    translate(half_vox2, np.array([0, 0, 0.5]) * pitch)
    bottom_half_plane = rec_array(half_vox2, x, y, [1, 0, 0], [0, 1, 0], pitch, pitch)

    # Cap the open sides of the lattice
    final_lattice = box_cap_sides_only(one_lattice, strut_width, chamfer_factor, pitch, x, y, z+1)

    all_geometry = [final_lattice] + top_half_plane + bottom_half_plane

    return combine_meshes(*all_geometry)


def hybrid_codedstructure_legacy(template, pitch, cap_mesh, voxel_meshes,  closed=True):
    """
    This is legacy code. Use hybrid_codedstructure for updated code.



    This function creates a lattice structure by placing individual voxels at locations indicated by a template.
    It can place are arbitrary number of different types of voxels, though currently allows for only one type of cap.
    The template should contain
    :param template: three-dimensional numpy array with integer codes for locations of voxels. Enter a 0 for no voxel
    placement, and integers starting at 1 for voxels of different types. The template must contain integer values, and
    the voxels must be coded starting at 1 (for example, if you have two different types of voxels, they must be coded
    "1" and "2" respectively in the template. They cannot be coded "2" and "3" or other values because of how the
    current connectivity checking code operates).
    :param pitch: lattice pitch
    :param cap_mesh: numpy stl mesh object of the cap geometry
    :param voxel_meshes: list of voxels to be used. ex. if there are two voxel types, voxel_meshes = [voxel_1, voxel_2]
    The order of the voxel meshes must correspond to their code in the template (first mesh in list is code "1" in
    template, second is code"2", etc.)
    :param closed: boolean value. Default True. Set to false to leave the geometry open (or uncapped).
    :return: numpy stl mesh object of the coded structure
    """

    lattice = []

    # Determine the x, y, and z size of the template (bounding box size in voxels)
    [x_size, y_size, z_size] = template.shape

    # Determine the number of voxel types and assign their codes for reading the template
    number_voxel_types = len(voxel_meshes)
    codes = set(np.arange(1, number_voxel_types + 1, dtype=np.int))
    print codes


    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if template[i, j, k] in codes:  # If a voxel is supposed to be placed
                    # The appropriate voxel to be placed is the mesh in voxel_meshes at (voxel code - 1) index
                    new_obj = mesh.Mesh(voxel_meshes[template[i, j, k] - 1].data.copy())  # Make a copy of the voxel
                    # Move the new voxel to the correct place
                    place_object(new_obj, pitch * i, pitch * j, pitch * k)
                    lattice += [new_obj]

                    # Even if not closing the lattice, want to check connectivity to ensure no hanging voxels
                    # check the connectivity and do optional capping
                    flag = 0  # Create a flag to detect minimum connectivity

                    # ----------check top--------
                    if k + 1 < z_size:  # if you aren't on the edge of the template
                        if template[i, j, k + 1] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so normal vectors correct
                                place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k + 1))
                                lattice += [cap_geo_top]
                        elif template[i, j, k + 1] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check atop voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_top = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so that normal vector is correct
                            place_object(cap_geo_top, pitch * i, pitch * j, pitch * (k + 1))
                            lattice += [cap_geo_top]

                    # --------check bottom---------
                    if k - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j, k - 1] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                                place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                                lattice += [cap_geo_bottom]
                        elif template[i, j, k - 1] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check below voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_bottom = mesh.Mesh(cap_mesh.data.copy())
                            place_object(cap_geo_bottom, pitch * i, pitch * j, pitch * k)
                            lattice += [cap_geo_bottom]

                    # ----------check right (positive X)----------
                    if i + 1 < x_size:  # if you aren't on the edge of the template
                        if template[i + 1, j, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_right.rotate([0, 1, 0], math.radians(90))
                                place_object(cap_geo_right, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap_geo_right]
                        elif template[i+1, j, k] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check right of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_right = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_right.rotate([0, 1, 0], math.radians(90))
                            place_object(cap_geo_right, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_right]

                    # ----------check left (negative X)---------
                    if i - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i - 1, j, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_left.rotate([0, 1, 0], math.radians(270))
                                place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap_geo_left]
                        elif template[i - 1, j, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check left of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_left = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_left.rotate([0, 1, 0], math.radians(270))
                            place_object(cap_geo_left, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                            lattice += [cap_geo_left]

                    # ----------check back (positive Y)---------
                    if j + 1 < y_size:  # if you aren't on the edge of the template
                        if template[i, j + 1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_back.rotate([1, 0, 0], math.radians(270))
                                place_object(cap_geo_back, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap_geo_back]
                        elif template[i, j + 1, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check back of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_back = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_back.rotate([1, 0, 0], math.radians(270))
                            place_object(cap_geo_back, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_back]

                    # check front (negative Y)
                    if j - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j - 1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                                cap_geo_front.rotate([1, 0, 0], math.radians(90))
                                place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap_geo_front]
                        elif template[i, j - 1, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check front of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))


                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            cap_geo_front = mesh.Mesh(cap_mesh.data.copy())
                            cap_geo_front.rotate([1, 0, 0], math.radians(90))
                            place_object(cap_geo_front, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                            lattice += [cap_geo_front]

                    # If the flag wasn't thrown, this voxel doesn't have any connectivity. Show an error
                    if flag == 0:
                        print(" There is a voxel in your template with zero connectivity.")

    return combine_meshes(*lattice)


def hybrid_codedstructure(template, pitch, voxel_meshes, voxel_cap_geos,  closed=True):
    """
    This function creates a lattice structure by placing individual voxels at locations indicated by a template.
    It can place are arbitrary number of different types of voxels, and allows for definition of capping logic for each
    voxel type.
    Limitations:
    Current version cannot detect whether a voxel is next to a voxel with a closed face, therefore requiring capping
    on that side. In other words, if a voxel is next to ANY type of voxel, it will assume that no cap is necessary
    on that side.
    :param template: three-dimensional numpy array with integer codes for locations of voxels. Enter a 0 for no voxel
    placement, and integers starting at 1 for voxels of different types. The template must contain integer values, and
    the voxels must be coded starting at 1 (for example, if you have two different types of voxels, they must be coded
    "1" and "2" respectively in the template. They cannot be coded "2" and "3" or other values because of how the
    current connectivity checking code operates).
    :param pitch: lattice pitch
    :param voxel_meshes: list of voxels to be used. ex. if there are two voxel types, voxel_meshes = [voxel_1, voxel_2]
    The order of the voxel meshes must correspond to their code in the template (first mesh in list is code "1" in
    template, second is code"2", etc.)
    :param voxel_cap_geos: This parameter can either be a *single* mesh object of the bottom cap geometry to be used on
    all voxels. OR it can be a list of capping geometries for each respective voxel type, the order of capping geometry
    entry corresponding to the order of voxel meshes entered in voxel_meshes. If a single mesh object is entered into
    the list, it should be the bottom cap geometry and will be used for all sides. Else, enter a list of cap geometry
    mesh objects of the form [top_cap, bottom_cap, right_cap, left_cap, back_cap, front_cap]. If no cap exists for
    a given side, enter a 0. Example: if default= [top_cap, 0, right_cap, left_cap, back_cap, front_cap] and capmesh is
    a single mesh object, if voxel_cap_geos = [ capmesh, default], the first voxel type will have the capmesh geometry
    capping all sides, and the second voxel type will be capped with the top_cap mesh object on the top side, right_cap
    mesh object on the right side, etc. The second voxel type will not have capping on the bottom since a 0 is entered.
    :param closed: boolean value. Default True. Set to false to leave the geometry open (or uncapped).
    :return: numpy stl mesh object of the coded structure
    """

    lattice = []

    # Determine the x, y, and z size of the template (bounding box size in voxels)
    [x_size, y_size, z_size] = template.shape

    # Determine the number of voxel types and assign their codes for reading the template
    number_voxel_types = len(voxel_meshes)
    codes = set(np.arange(1, number_voxel_types + 1, dtype=np.int))
    print 'Detected voxel codes: '
    print codes



    if isinstance(voxel_cap_geos, list) is True: # if the voxel cap geometry is specified for each voxel separately

        if len(voxel_cap_geos) is not len(voxel_meshes):
            print 'You have not specified capping geometries for every voxel type.'

        for idx, val in enumerate(voxel_cap_geos):
            if isinstance(voxel_cap_geos[idx], list) is False:  # if it isn't a list
                # assume that the thing entered was the bottom cap geometry mesh for this type of voxel
                single_cap_geometry = mesh.Mesh(voxel_cap_geos[idx].data.copy())

                # make the default capping geometry
                cap_geo_top = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so normal vectors correct
                cap_geo_bottom = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_right = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_right.rotate([0, 1, 0], math.radians(90))
                cap_geo_left = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_left.rotate([0, 1, 0], math.radians(270))
                cap_geo_back = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_back.rotate([1, 0, 0], math.radians(270))
                cap_geo_front = mesh.Mesh(single_cap_geometry.data.copy())
                cap_geo_front.rotate([1, 0, 0], math.radians(90))

                default_caps = [cap_geo_top, cap_geo_bottom, cap_geo_right, cap_geo_left, cap_geo_back, cap_geo_front]

                voxel_cap_geos[idx] = default_caps

    else:  # user should have input a single mesh object of the bottom voxel cap (with correct normals)
        single_cap_geometry = mesh.Mesh(voxel_cap_geos.data.copy())
        voxel_cap_geos = []

        # make the default capping geometry
        cap_geo_top = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so normal vectors correct
        cap_geo_bottom = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_right = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_right.rotate([0, 1, 0], math.radians(90))
        cap_geo_left = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_left.rotate([0, 1, 0], math.radians(270))
        cap_geo_back = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_back.rotate([1, 0, 0], math.radians(270))
        cap_geo_front = mesh.Mesh(single_cap_geometry.data.copy())
        cap_geo_front.rotate([1, 0, 0], math.radians(90))

        default_caps = [cap_geo_top, cap_geo_bottom, cap_geo_right, cap_geo_left, cap_geo_back, cap_geo_front]

        for instance in voxel_meshes:
            voxel_cap_geos += [default_caps]

    # This is the reading of template, placing of voxels, and capping procedure
    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if template[i, j, k] in codes:  # If a voxel is supposed to be placed
                    # The appropriate voxel to be placed is the mesh in voxel_meshes at (voxel code - 1) index
                    new_obj = mesh.Mesh(voxel_meshes[template[i, j, k] - 1].data.copy())  # Make a copy of the voxel
                    # Move the new voxel to the correct place
                    place_object(new_obj, pitch * i, pitch * j, pitch * k)
                    lattice += [new_obj]

                    # Even if not closing the lattice, want to check connectivity to ensure no hanging voxels
                    # check the connectivity and do optional capping
                    flag = 0  # Create a flag to detect minimum connectivity

                    # ----------check top--------
                    if k + 1 < z_size:  # if you aren't on the edge of the template
                        if template[i, j, k + 1] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                # check if there is a cap to place there
                                if voxel_cap_geos[template[i, j, k] - 1][0] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][0].data.copy())
                                    place_object(cap, pitch * i, pitch * j, pitch * (k + 1))
                                    lattice += [cap]
                        elif template[i, j, k + 1] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check atop voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # check if there is a cap to place there
                            if voxel_cap_geos[template[i, j, k] - 1][0] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][0].data.copy())
                                place_object(cap, pitch * i, pitch * j, pitch * (k + 1))
                                lattice += [cap]

                    # --------check bottom---------
                    if k - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j, k - 1] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the top
                                # check if there is a specified cap geometry to place there
                                if voxel_cap_geos[template[i, j, k] - 1][1] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][1].data.copy())
                                    place_object(cap, pitch * i, pitch * j, pitch * k)
                                    lattice += [cap]
                        elif template[i, j, k - 1] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check below voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # check if there is a specified cap geometry to place there
                            if voxel_cap_geos[template[i, j, k] - 1][1] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][1].data.copy())
                                place_object(cap, pitch * i, pitch * j, pitch * k)
                                lattice += [cap]

                    # ----------check right (positive X)----------
                    if i + 1 < x_size:  # if you aren't on the edge of the template
                        if template[i + 1, j, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                # check if there is a specified cap geometry to place there
                                if voxel_cap_geos[template[i, j, k] - 1][2] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][2].data.copy())
                                    place_object(cap, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                                    lattice += [cap]
                        elif template[i+1, j, k] in codes:  # if a valid voxel code (there is a voxel there)
                            flag = 1
                        else:
                            print('Template Error. Check right of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            # check if there is a specified cap geometry to place there
                            if voxel_cap_geos[template[i, j, k] - 1][2] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][2].data.copy())
                                place_object(cap, pitch * i + pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap]

                    # ----------check left (negative X)---------
                    if i - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i - 1, j, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the left
                                # check if there is a specified cap geometry to place there
                                if voxel_cap_geos[template[i, j, k] - 1][3] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][3].data.copy())
                                    place_object(cap, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                                    lattice += [cap]
                        elif template[i - 1, j, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check left of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            # check if there is a specified cap geometry to place there
                            if voxel_cap_geos[template[i, j, k] - 1][3] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][3].data.copy())
                                place_object(cap, pitch * i - pitch / 2, pitch * j, pitch * k + pitch / 2)
                                lattice += [cap]

                    # ----------check back (positive Y)---------
                    if j + 1 < y_size:  # if you aren't on the edge of the template
                        if template[i, j + 1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                # check if there is a specified cap geometry to place there
                                if voxel_cap_geos[template[i, j, k] - 1][4] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][4].data.copy())
                                    place_object(cap, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                                    lattice += [cap]
                        elif template[i, j + 1, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check back of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))

                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            # check if there is a specified cap geometry to place there
                            if voxel_cap_geos[template[i, j, k] - 1][4] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][4].data.copy())
                                place_object(cap, pitch * i, pitch * j + pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap]

                    # -----------check front (negative Y)---------------
                    if j - 1 >= 0:  # if you aren't on the edge of the template
                        if template[i, j - 1, k] == 0:  # if there isn't a voxel there
                            if closed:
                                # place a cap on the right
                                # check if there is a specified cap geometry to place there
                                if voxel_cap_geos[template[i, j, k] - 1][5] is not 0:
                                    cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][5].data.copy())
                                    place_object(cap, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                                    lattice += [cap]
                        elif template[i, j - 1, k] in codes:
                            flag = 1
                        else:
                            print('Template Error. Check front of voxel  x = {0} y = {1} z = {2}'.format(i, j, k))


                    else:
                        if closed:
                            # You are on the edge, so place a cap
                            # check if there is a specified cap geometry to place there
                            if voxel_cap_geos[template[i, j, k] - 1][5] is not 0:
                                cap = mesh.Mesh(voxel_cap_geos[template[i, j, k] - 1][5].data.copy())
                                place_object(cap, pitch * i, pitch * j - pitch / 2, pitch * k + pitch / 2)
                                lattice += [cap]

                    # If the flag wasn't thrown, this voxel doesn't have any connectivity. Show an error
                    if flag == 0:
                        print(" There is a voxel in your template with zero connectivity.")

    return combine_meshes(*lattice)

def create_test_template():
    """
    Creates a template to test the lattice_codedstructure function.
    :return: template
    """

    template = np.zeros((3, 3, 3), dtype=np.int)
    template[1, 1, 1] = 2
    template[1, 0, 1] = 1
    template[0, 1, 1] = 3
    template[1, 0, 0] = 1
    template[2, 1, 1] = 1
    template[2, 2, 1] = 1

    return template


def ct_template():

    template = np.zeros((25, 10, 25), dtype=np.int)

    template[:, :, 0:12] = 1
    template[15:25, :, 12] = 1
    template[:, :, 13:25] = 1


    return template


def hybrid_template():
    x = 11
    y = 11
    z = 11

    template = np.zeros((x, y, z), dtype=np.int)
    template[:, :, :] = 1
    template[1, :, :] = 2
    template[3, :, :] = 2
    template[5, :, :] = 2
    template[7, :, :] = 2
    template[9, :, :] = 2
    template[:, 1, :] = 2
    template[:, 3, :] = 2
    template[:, 5, :] = 2
    template[:, 7, :] = 2
    template[:, 9, :] = 2


    return template


def compression_template():
    """
    Creates a template for a compression specimen with type 1 voxel, with top and bottom voxel planes being voxel type 2
    :return: numpy 3d array template
    """
    x = 10
    y = 10
    z = 10

    template = np.zeros((x,y,z), dtype=np.int)

    template[:, :, :] = 1
    template[:, :, 0] = 2
    template[:, :, z-1] = 3

    return template

def preview_mesh(*args):
    """
    This function plots numpy stl mesh objects entered into args. Note it will scale the preview plot based on the last mesh
    object entered.
    :param args: mesh objects to plot  ex.- preview_mesh(mesh1, mesh2, mesh3)
    :return:
    """
    print ("...preparing preview...")
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    for mesh_obj in args:
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh_obj.vectors))

    # Auto scale to the mesh size. Note it will choose the last mesh
    scale = mesh_obj.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def pitch_from_relden(relden, cf, sw):
    """
    This function calculates the pitch of cuboct of a given relative density, chamfer factor, and strut width.
    :param relden: float. Desired relative density
    :param cf: float. Chamfer factor of voxel
    :param sw: float. strut width of voxel
    :return: lattice pitch
    """
    chamheight = sw / cf
    l_2 = sw / 2 + chamheight
    l_3 = l_2 + sw * np.cos(math.radians(45))  # horizontal position of points
    l_4 = np.sqrt(2) * (l_3 - sw / 2)
    tan_theta = ((l_3 - l_2) / ((l_4 / 2) - (np.sqrt(2) * chamheight / 2)))
    v1 = l_2 * (sw * sw + 4 * sw * (l_3 - sw / 2) + 2 * (l_3 - sw / 2) * (l_3 - sw / 2))
    h = (l_4 / 2) * tan_theta
    hs = chamheight * tan_theta * np.sqrt(2) / 2
    v2 = ((l_4 * l_4 * h) - (2 * (chamheight * chamheight * hs))) / 3
    v3 = 4 * sw * (0.5 * (l_3 - l_2) * (l_3 - l_2) + (l_3 - l_2) * chamheight)
    v4 = sw * sw * (l_3 - l_2)
    node_volume = v1 + v2 + v3 + v4

    c1 = relden
    c2 = (-6) * np.sqrt(2)*sw *sw
    c3 = -6*node_volume + 12*sw*sw*np.sqrt(2)*(l_2 + l_3)
    return max(np.roots([c1, 0, c2, c3]))


def main():
    # The following code is an example of using the hybrid structure code
    pitch = 30
    strut_width = 2
    chamfer_factor = 2.75

    template = create_test_template()

    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    two_voxel = hybrid_voxel(strut_width * 0.5, chamfer_factor, pitch, strut_width)
    three_voxel = half_voxel(strut_width, chamfer_factor, pitch)

    capmesh = cap_cuboct(strut_width, chamfer_factor)

    # make the default capping geometry
    cap_geo_top = mesh.Mesh(capmesh.data.copy())
    cap_geo_top.rotate([1, 0, 0], math.radians(180))  # rotate so normal vectors correct
    cap_geo_bottom = mesh.Mesh(capmesh.data.copy())
    cap_geo_right = mesh.Mesh(capmesh.data.copy())
    cap_geo_right.rotate([0, 1, 0], math.radians(90))
    cap_geo_left = mesh.Mesh(capmesh.data.copy())
    cap_geo_left.rotate([0, 1, 0], math.radians(270))
    cap_geo_back = mesh.Mesh(capmesh.data.copy())
    cap_geo_back.rotate([1, 0, 0], math.radians(270))
    cap_geo_front = mesh.Mesh(capmesh.data.copy())
    cap_geo_front.rotate([1, 0, 0], math.radians(90))

    default_caps = [0, cap_geo_bottom, cap_geo_right, cap_geo_left, cap_geo_back, cap_geo_front]

    hybrid_structure = hybrid_codedstructure(template, pitch, [one_voxel, two_voxel, three_voxel],
                                             [capmesh, capmesh, default_caps])
    hybrid_structure.save('hybrid_structure_test.stl')
    preview_mesh(hybrid_structure)




if __name__ == "__main__":
    main()

