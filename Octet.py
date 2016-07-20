from stl import mesh
import math
import numpy as np
from matplotlib import pyplot
from mpl_toolkits import mplot3d

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
    point_A_prime = [strutwidth*(l_3-halfw)/4 + halfw, l_3, halfw]
    point_B_prime = [l_3,  strutwidth * (l_3 - halfw) / 4 + halfw, halfw]
    point_H_prime = [-(strutwidth*(l_3-halfw)/4+ halfw), l_3, halfw]
    point_A_prime_bottom = [strutwidth * (l_3 - halfw) / 4 + halfw, l_3, 0]
    point_B_prime_bottom = [l_3, strutwidth * (l_3 - halfw) / 4 + halfw, 0]
    point_H_prime_bottom = [-(strutwidth * (l_3 - halfw) / 4 + halfw), l_3, 0]

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
    point_A_prime = [strutwidth*(l_3-halfw)/4 + halfw, l_3, halfw]
    point_B_prime = [l_3,  strutwidth * (l_3 - halfw) / 4 + halfw, halfw]
    point_H_prime = [-(strutwidth*(l_3-halfw)/4+ halfw), l_3, halfw]
    point_A_prime_bottom = [strutwidth * (l_3 - halfw) / 4 + halfw, l_3, 0]
    point_B_prime_bottom = [l_3, strutwidth * (l_3 - halfw) / 4 + halfw, 0]
    point_H_prime_bottom = [-(strutwidth * (l_3 - halfw) / 4 + halfw), l_3, 0]

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
    x = 10
    y = 10
    z = 10

    test_node = node(strut_width, chamfer_factor)
    test_voxel = voxel(strut_width, chamfer_factor, pitch)
    preview_mesh(test_voxel)

if __name__ == "__main__":
    main()