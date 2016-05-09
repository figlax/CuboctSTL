from stl import mesh
import math
import numpy as np
from matplotlib import pyplot
from mpl_toolkits import mplot3d

def voxel(strut_width, chamfer_factor, pitch):
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
    bottom_struts = arraypolar2(strut1, [0, 0, 0.5], 4)
    strut1.rotate([1, 0, 0], math.radians(180))
    translate(strut1, np.array([0, 0, 1]) * pitch)
    top_struts = arraypolar2(strut1, [0, 0, 0.5], 4)

    strut2 = strut(strut_width, chamfer_factor, pitch)
    strut2.rotate([1, 0, 0], math.radians(90))
    translate(strut2, np.array([0, -0.5, 0.5]) * pitch)
    side_struts = arraypolar2(strut2, [0, 0, 0.5], 4)

    combined_geometry = bottom_struts + top_struts + side_struts +vnodes


    return combine_meshes(*combined_geometry)






def node(strutwidth, chamfactor):

    #strutwidth = 2;
    #chamfactor = 2 + 2/3;

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
    topcap['vectors'][0] = np.array([point1, point2, point3])
    topcap['vectors'][1] = np.array([point1, point3, point4])
    topcap['vectors'][2] = np.array([point1, point4, point5])
    topcap['vectors'][3] = np.array([point1, point5, point6])
    topcap['vectors'][4] = np.array([point1, point6, point7])
    topcap['vectors'][5] = np.array([point1, point7, point8])
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
    chamfersides['vectors'][0] = np. array([point1, point1s, point2s])
    chamfersides['vectors'][1] = np. array([point1, point2s, point2])
    chamfersides['vectors'][2] = np. array([point3, point3s, point4s])
    chamfersides['vectors'][3] = np. array([point3, point4s, point4])
    chamfersides['vectors'][4] = np. array([point5, point5s, point6s])
    chamfersides['vectors'][5] = np. array([point5, point6s, point6])
    chamfersides['vectors'][6] = np. array([point7, point7s, point8s])
    chamfersides['vectors'][7] = np. array([point7, point8s, point8])
    chamfersides_mesh = mesh.Mesh(chamfersides)


    # Define the rectangular sides
    sides = np.zeros(4, dtype=mesh.Mesh.dtype)

    point1b = [halfw, l_3, 0]
    point2b = [l_3, halfw,0]
    point8b = [-halfw, l_3, 0]
    sides['vectors'][0] = np. array([point1s, point1b, point2b])
    sides['vectors'][1] = np. array([point1s, point2b, point2s])
    sides['vectors'][2] = np. array([point8s, point8b, point1b])
    sides['vectors'][3] = np. array([point8s, point1b, point1s])
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

def strut(strutwidth, chamfactor,  pitch):
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
    singlestrut_geo['vectors'][0] = np.array([point2_copy, point2n, point2sn])
    singlestrut_geo['vectors'][1] = np.array([point2_copy, point2sn, point2s_copy])
    singlestrut_geo['vectors'][2] = np.array([point3_copy, point3n, point3sn])
    singlestrut_geo['vectors'][3] = np.array([point3_copy, point3sn, point3s_copy])
    singlestrut_geo['vectors'][4] = np.array([point2_copy, point2n, point3n])
    singlestrut_geo['vectors'][5] = np.array([point2_copy, point3n, point3_copy])
    singlestrut_geo['vectors'][6] = np.array([point2s_copy, point2sn, point3sn])
    singlestrut_geo['vectors'][7] = np.array([point2s_copy, point3sn, point3s_copy])

    finalsinglestrut = mesh.Mesh(singlestrut_geo)
    return finalsinglestrut

def translate(meshobj, tvect):
    vects = meshobj.vectors
    for vsi,vectset in enumerate(vects):
        for vi,vector in enumerate(vectset):
            meshobj.vectors[vsi][vi] = vector+tvect

def rotation_matrix( axis, theta):
    '''
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


def arraypolar(meshobjects_list, r_axis, num):
    # This function takes an array of mesh objects that will be arrayed meshobjects
    # r_axis is the axis of rotation
    # num is the number of items in the array mesh objects

    for i in range(0, num):
        meshobjects_list[i].rotate(r_axis, math.radians((360 / num) * i))
    return meshobjects_list

def arraypolar2(m_obj, r_axis, num, rotation_point=None):
    array_objects = list()
    array_objects.append(m_obj)

    for i in range(num):
        obj = mesh.Mesh(m_obj.data.copy())
        obj.rotate(r_axis, math.radians((360 / num) * i), rotation_point)
        array_objects.append(obj)
    return array_objects

def combine_meshes(*args):
    combined_data = np.concatenate([m_obj.data for m_obj in args])
    return mesh.Mesh(combined_data)

def lattice_array(voxel_mesh, pitch, x, y, z):
    lattice = list()
    lattice.append(voxel, mesh)

    if x == 1:
        x = 1
    else:
        for i in range(x):
            new_obj = mesh.Mesh(voxel_mesh.data.copy())  #Make a copy of the voxel
            translate(new_obj, [1, 0, 0]*pitch*i)
            lattice.append(new_obj)
    if y == 1:
        y = 1
    else:
        xline = lattice
        for j in range(y):
            for thing in xline:
                new_obj = mesh.Mesh(thing.data.copy())  # Make a copy of the voxel
                translate(new_obj, [1, 0, 0] * pitch * j)
                lattice.append(new_obj)

    






def main():
    pitch = 30
    strut_width = 2
    chamfer_factor = 2.666

    one_voxel = voxel(strut_width, chamfer_factor, pitch)

    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(one_voxel.vectors))


    # Auto scale to the mesh size
    scale = one_voxel.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)

    print mesh.Mesh.dtype
    # Show the plot to the screen
    pyplot.show()




main()