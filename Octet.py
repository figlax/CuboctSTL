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
    cap_geo = cap(strutwidth, chamfactor, pitch)

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


def node(strutwidth, chamfactor):
    """
    This function creates a mesh of an open full octet node.
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
    chamfersides = np.zeros(8, dtype=mesh.Mesh.dtype)
    chamfersides['vectors'][0] = np. array([point1, point2s, point1s])
    chamfersides['vectors'][1] = np. array([point1, point2, point2s])
    chamfersides['vectors'][2] = np. array([point3, point4s, point3s])
    chamfersides['vectors'][3] = np. array([point3, point4, point4s])
    chamfersides['vectors'][4] = np. array([point5, point6s, point5s])
    chamfersides['vectors'][5] = np. array([point5, point6, point6s])
    chamfersides['vectors'][6] = np. array([point7, point8s, point7s])
    chamfersides['vectors'][7] = np. array([point7, point8, point8s])
    chamfersides_mesh = mesh.Mesh(chamfersides)

    # Define second chamfer to bottom octet struts
    l_4 = np.sqrt(2)*(l_3 - halfw)  # Is this square root going to be a problem for error prop?
    point_A_prime = [l_3 - np.sqrt(2)*strutwidth/2, l_3, halfw]
    point_B_prime = [l_3,  l_3 - np.sqrt(2)*strutwidth/2, halfw]
    point_H_prime = [-(l_3 - np.sqrt(2)*strutwidth/2), l_3, halfw]
    point_A_prime_bottom = [l_3 - np.sqrt(2)*strutwidth/2, l_3, 0]
    point_B_prime_bottom = [l_3, l_3 - np.sqrt(2)*strutwidth/2, 0]
    point_H_prime_bottom = [-(l_3 - np.sqrt(2)*strutwidth/2), l_3, 0]

    sides = np.zeros(6, dtype=mesh.Mesh.dtype)
    sides['vectors'][0] = np.array([point1s, point2s, point_A_prime])
    sides['vectors'][1] = np.array([point2s, point_B_prime, point_A_prime])
    sides['vectors'][2] = np.array([point1s, point_A_prime, point8s])
    sides['vectors'][3] = np.array([point8s, point_A_prime, point_H_prime])
    sides['vectors'][4] = np.array([point_H_prime, point_A_prime, point_H_prime_bottom])
    sides['vectors'][5] = np.array([point_H_prime_bottom, point_A_prime, point_A_prime_bottom])
    sides_mesh = mesh.Mesh(sides.copy())
    all_sides = arraypolar(sides_mesh, [0, 0, 1], 4)

    # Make final mesh for the open node geometry
    combined_geometry = [chamfersides_mesh] + [top] + all_sides

    return combine_meshes(*combined_geometry)

def corner_node(strutwidth, chamfactor):
    """
    This function creates a mesh of an open quarter (first quadrant) octet node to be used on the corners
    of the octet voxel.
    :param strutwidth: float
    :param chamfactor: float
    :return: finalnodemesh: numpy stl mesh object of open node
    """

    # Calculate commonly used values for geometry definition
    chamheight = strutwidth / chamfactor
    halfw = strutwidth / 2

    # Define geometry of the top cap
    topcap = np.zeros(3, dtype=mesh.Mesh.dtype)
    # Calculate the height and half-length of top octogonal cap
    h = chamheight + (strutwidth * np.sin(np.pi / 4) + strutwidth / 2)
    l_2 = strutwidth / 2 + chamheight
    point1 = [halfw, l_2, h]
    point2 = [l_2, halfw, h]
    point3 = [l_2, 0, h]
    point4 = [0, 0, h]
    point5 = [0, l_2, h]

    # Define with right hand rule to assure outward facing normal
    topcap['vectors'][0] = np.array([point1, point5, point2])
    topcap['vectors'][1] = np.array([point2, point5, point3])
    topcap['vectors'][2] = np.array([point3, point5, point4])

    top = mesh.Mesh(topcap)

    # Define Geometry of the chamfered sides
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
    chamfersides = np.zeros(2, dtype=mesh.Mesh.dtype)
    chamfersides['vectors'][0] = np. array([point1, point2s, point1s])
    chamfersides['vectors'][1] = np. array([point1, point2, point2s])
    chamfersides_mesh = mesh.Mesh(chamfersides)

    # Define second chamfer to bottom octet struts
    l_4 = np.sqrt(2)*(l_3 - halfw)  # Is this square root going to be a problem for error prop?
    point_A_prime = [l_3 - np.sqrt(2)*strutwidth/2, l_3, halfw]
    point_B_prime = [l_3,  l_3 - np.sqrt(2)*strutwidth/2, halfw]
    point_H_prime_quarter = [0, l_3, halfw]
    point_c_prime_quarter = [l_3, 0, halfw]
    point_A_prime_bottom = [l_3 - np.sqrt(2)*strutwidth/2, l_3, 0]
    point_B_prime_bottom = [l_3, l_3 - np.sqrt(2)*strutwidth/2, 0]
    point_H_prime_quarter_bottom = [0, l_3, 0]
    point_c_prime_quarter_bottom = [l_3, 0, 0]

    sides = np.zeros(10, dtype=mesh.Mesh.dtype)
    sides['vectors'][0] = np.array([point1s, point2s, point_A_prime])
    sides['vectors'][1] = np.array([point2s, point_B_prime, point_A_prime])
    sides['vectors'][2] = np.array([point1s, point_A_prime, [0, l_3, hs]])
    sides['vectors'][3] = np.array([[0, l_3, hs], point_A_prime, point_H_prime_quarter])
    sides['vectors'][4] = np.array([point_H_prime_quarter, point_A_prime, point_H_prime_quarter_bottom])
    sides['vectors'][5] = np.array([point_H_prime_quarter_bottom, point_A_prime, point_A_prime_bottom])
    sides['vectors'][6] = np.array([[l_3, 0, hs], point_B_prime, point2s])
    sides['vectors'][7] = np.array([point_c_prime_quarter, point_B_prime, [l_3, 0, hs]])
    sides['vectors'][8] = np.array([point_c_prime_quarter, point_c_prime_quarter_bottom, point_B_prime])
    sides['vectors'][9] = np.array([point_c_prime_quarter_bottom, point_B_prime_bottom, point_B_prime])
    sides_mesh = mesh.Mesh(sides.copy())



    # Make final mesh for the open node geometry
    combined_geometry = [chamfersides_mesh] + [top] + [sides_mesh]

    return combine_meshes(*combined_geometry)


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

def corner(strut_width, chamfer_factor, pitch):
    """

    :param strut_width:
    :param chamfer_factor:
    :param pitch:
    :return:
    """

    # Make a base node for the corner
    cnode = corner_node(strut_width, chamfer_factor)

    # Geometry Parameters
    # Calculate commonly used values for geometry definition
    chamheight = strut_width / chamfer_factor
    halfw = strut_width / 2
    halfp = pitch / 2
    h = chamheight + (strut_width * np.sin(np.pi / 4) + strut_width / 2)  # height of top cap
    l_2 = strut_width / 2 + chamheight  # horizontal position of points on topcap
    hs = l_2  # height of side points of node
    l_3 = l_2 + strut_width * np.cos(np.pi / 4)  # horizontal position of points

    # re-use points from strut code, with alteration to create half-struts
    point2 = [l_2, halfw, h]
    point3 = [l_2, 0, h]
    point2s = [l_3, halfw, hs]
    point3s = [l_3, 0, hs]
    # new points to attach to on side node
    point2n = [halfp - h, halfw, halfp - l_2]
    point3n = [halfp - h, 0, halfp - l_2]
    point2sn = [halfp - hs, halfw, halfp - l_3]
    point3sn = [halfp - hs, 0, halfp - l_3]


    face_halfstrut_geo = np.zeros(6, dtype=mesh.Mesh.dtype)
    face_halfstrut_geo['vectors'][0] = np.array([point3, point3n, point2n])
    face_halfstrut_geo['vectors'][1] = np.array([point3, point2n, point2])
    face_halfstrut_geo['vectors'][2] = np.array([point2, point2n, point2sn])
    face_halfstrut_geo['vectors'][3] = np.array([point2, point2sn, point2s])
    face_halfstrut_geo['vectors'][4] = np.array([point2s, point2sn, point3s])
    face_halfstrut_geo['vectors'][5] = np.array([point3s, point2sn, point3sn])
    face_halfstrut = mesh.Mesh(face_halfstrut_geo)

    face_halfstrut2 = mesh.Mesh(face_halfstrut.data.copy())
    face_halfstrut2.rotate([1, 0, 1], math.radians(180))
    face_halfstrut2.rotate([0,0,1], math.radians(270))

    face_halfstrut3 = mesh.Mesh(face_halfstrut.data.copy())
    face_halfstrut3.rotate([1, 0, 1], math.radians(270))
    face_halfstrut3.rotate([0, 1, 0], math.radians(-45))
    face_halfstrut3.rotate([0, 0, 1], math.radians(-45))

    combined_geometry = [cnode] + [face_halfstrut] + [face_halfstrut2] + [face_halfstrut3]

    return combine_meshes(*combined_geometry)

def cap(strut_width, chamfer_factor, pitch):
    """
    Creates a mesh of the face capping geometry for the octet voxel.
    :param strut_width: float
    :param chamfer_factor: float
    :param pitch: float
    :return: numpy stl mesh object of cap geometry
    """
    # Calculate commonly used values for geometry definition
    chamheight = strut_width / chamfer_factor
    halfw = strut_width / 2
    l_2 = strut_width / 2 + chamheight
    l_3 = l_2 + strut_width * np.cos(np.pi / 4)  # horizontal position of points
    l_4 = np.sqrt(2) * (l_3 - halfw)  # Is this square root going to be a problem for error prop?

    # Define points for node cap
    face_strut_pos = l_3 - np.sqrt(2) * strut_width / 2
    point_A_prime_bottom = [face_strut_pos, l_3, 0]
    point_B_prime_bottom = [l_3, face_strut_pos, 0]
    point_C_prime_bottom = [l_3, -face_strut_pos, 0]
    point_D_prime_bottom = [face_strut_pos, -l_3, 0]
    point_E_prime_bottom = [-face_strut_pos, -l_3, 0]
    point_F_prime_bottom = [-l_3, -face_strut_pos, 0]
    point_G_prime_bottom = [-l_3, face_strut_pos, 0]
    point_H_prime_bottom = [-face_strut_pos, l_3, 0]
    # Use these points to cap area over the node
    node_cap_geo = np.zeros(6, dtype=mesh.Mesh.dtype)
    node_cap_geo['vectors'][0] = np.array([point_H_prime_bottom, point_A_prime_bottom, point_G_prime_bottom])
    node_cap_geo['vectors'][1] = np.array([point_A_prime_bottom, point_B_prime_bottom, point_G_prime_bottom])
    node_cap_geo['vectors'][2] = np.array([point_G_prime_bottom, point_B_prime_bottom, point_F_prime_bottom])
    node_cap_geo['vectors'][3] = np.array([point_B_prime_bottom, point_C_prime_bottom, point_F_prime_bottom])
    node_cap_geo['vectors'][4] = np.array([point_C_prime_bottom, point_E_prime_bottom, point_F_prime_bottom])
    node_cap_geo['vectors'][5] = np.array([point_C_prime_bottom, point_D_prime_bottom, point_E_prime_bottom])
    node_cap = mesh.Mesh(node_cap_geo)

    # Define corner node cap
    corner_cap_geo = np.zeros(3, dtype=mesh.Mesh.dtype)
    corner_cap_geo['vectors'][0] = np.array([[0,0,0], [0, l_3, 0], [l_3, 0, 0]])
    corner_cap_geo['vectors'][1] = np.array([[0, l_3, 0], point_B_prime_bottom, [l_3, 0, 0]])
    corner_cap_geo['vectors'][2] = np.array([[0, l_3, 0], point_A_prime_bottom, point_B_prime_bottom])
    corner_cap = mesh.Mesh(corner_cap_geo)

    # Define strut cap
    point_corner_A_prime_bottom = [-pitch/2 + point_A_prime_bottom[0], -pitch/2 + point_A_prime_bottom[1], 0]
    point_corner_B_prime_bottom = [-pitch/2 + point_B_prime_bottom[0], -pitch/2 + point_B_prime_bottom[1], 0]
    strut_cap_geo = np.zeros(2, dtype=mesh.Mesh.dtype)
    strut_cap_geo['vectors'][0] = np.array([point_corner_A_prime_bottom, point_F_prime_bottom, point_E_prime_bottom])
    strut_cap_geo['vectors'][1] = np.array([point_corner_A_prime_bottom, point_E_prime_bottom, point_corner_B_prime_bottom])
    strut_cap = mesh.Mesh(strut_cap_geo)

    translate(corner_cap, np.array([-1, -1, 0])* pitch/2)
    corners = arraypolar(corner_cap, [0, 0, 1], 4)
    struts = arraypolar(strut_cap, [0, 0, 1], 4)

    combined_geometry = corners + struts + [node_cap]

    return combine_meshes(*combined_geometry)


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

    # Define voxel corner nodes
    v_corner = corner(strut_width, chamfer_factor, pitch)
    translate(v_corner, np.array([-1, 0, 0]) * pitch / 2)
    translate(v_corner, np.array([0, -1, 0]) * pitch / 2)
    bottom_corners = arraypolar(v_corner, [0, 0, 1], 4)
    v_corner_top = corner(strut_width, chamfer_factor, pitch)
    v_corner_top.rotate([0, 1, 0], math.radians(180))
    translate(v_corner_top, np.array([1, 0, 0]) * pitch / 2)
    translate(v_corner_top, np.array([0, -1, 0]) * pitch / 2)
    translate(v_corner_top, np.array([0, 0, 1])* pitch)
    top_corners = arraypolar(v_corner_top, [0, 0, 1], 4)

    # Define voxel struts using strut function
    strut1 = strut(strut_width, chamfer_factor, pitch)
    bottom_struts = arraypolar(strut1, [0, 0, 1], 4)
    strut1.rotate([1, 0, 0], math.radians(180))
    translate(strut1, np.array([0, 0, 1]) * pitch)
    top_struts = arraypolar(strut1, [0, 0, 1], 4)

    strut2 = side_strut(strut_width, chamfer_factor, pitch)
    side_struts = arraypolar(strut2, [0, 0, 0.5], 4)

    combined_geometry = bottom_struts + top_struts + side_struts + vnodes + bottom_corners + top_corners

    return combine_meshes(*combined_geometry)

def arraypolar(m_obj, r_axis, num, rotation_point=None, angle=360):
    """
    This function arrays mesh objects in evenly spaced circular pattern.
    :param m_obj: numpy stl mesh object to array
    :param r_axis: rotation axis ex. [0, 0, 1]
    :param num: number of items to array
    :param rotation_point: point to place rotation axis if not center
    :param angle: angle to fill if not 360 degrees. Enter in degrees
    :return: list of arrayed objects
    """
    array_objects = list()
    for i in range(num):
        obj = mesh.Mesh(m_obj.data.copy())
        obj.rotate(r_axis, math.radians((angle / num) * i), rotation_point)
        array_objects.append(obj)
    return array_objects


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


def combine_meshes(*args):
    """
    This function combines a list or lists of mesh objects into a single mesh object
    :param args: list of mesh objects
    :return: numpy stl mesh object of combined geometries
    """
    combined_data = np.concatenate([m_obj.data for m_obj in args])
    return mesh.Mesh(combined_data, remove_duplicate_polygons = True)


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


def main():
    pitch = 10
    strut_width = 0.6
    chamfer_factor = 3
    x = 3
    y = 3
    z = 3

    test_corner = corner(strut_width, chamfer_factor, pitch)
    test_node = node(strut_width, chamfer_factor)
    test_voxel = voxel(strut_width, chamfer_factor, pitch)
    test_corner_node = corner_node(strut_width, chamfer_factor)
    test_lattice = make_lattice(strut_width, chamfer_factor, pitch, x, y, z)

    preview_mesh(test_lattice)
    test_voxel.save('octet_test_voxel.stl')


if __name__ == "__main__":
    main()