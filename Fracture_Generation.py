from CuboctSTL_v0 import *

def ct_template2():

    template = np.zeros((26, 10, 26), dtype=np.int)

    template[:, :, 0:12] = 1
    template[15:26, :, 12] = 1
    template[15:26, :, 13] = 1
    template[1:15, :, 12] = 2
    template[1:15, :, 13] = 3
    template[:, :, 14:26] = 1


    return template


def ct_template3():

    template = np.zeros((11, 5, 10), dtype=np.int)

    template[:, :, 0:4] = 1
    template[6:11, :, 4] = 1
    template[6:11, :, 5] = 1
    template[1:6, :, 4] = 2
    template[1:6, :, 5] = 3
    template[:, :, 6:10] = 1


    return template

def ct_template2_withholes():

    template = np.zeros((26, 10, 26), dtype=np.int)

    template[:, :, 0:4] = 1

    template[0:4, :, 4] = 1
    template[6:26, :, 4] = 1

    template[0:2, :, 5:9] = 1
    template[8:26, :, 5:9] = 1

    template[0:3, :, 9] = 1
    template[7:26, :, 9] = 1

    template[:, :, 10:12] = 1

    template[15:26, :, 12] = 1
    template[15:26, :, 13] = 1
    template[1:15, :, 12] = 2
    template[1:15, :, 13] = 3

    template[:, :, 14:16] = 1

    template[0:3, :, 16] = 1
    template[7:26, :, 16] = 1

    template[0:2, :, 17:21] = 1
    template[8:26, :, 17:21] = 1

    template[0:4, :, 21] = 1
    template[6:26, :, 21] = 1

    template[:, :, 22:26] = 1

    return template





def main():
    """
        This script generates the basic structure for a c(t) fracture specimen.
        It will generally need to be repaired in post-processing
        :return:
        """
    pitch = 15
    sw = 0.7
    cf = 5


    one_voxel = voxel(sw, cf, pitch)
    two_voxel = half_voxel(sw, cf, pitch)
    three_voxel = half_voxel(sw, cf, pitch)
    three_voxel.rotate([1,0,0], math.radians(180))
    translate(three_voxel, np.array([0, 0, pitch]))
    template_1 = ct_template3()
    capmesh = cap_cuboct(sw, cf)

    # make all side caps that we will use in defining capping for half-voxel
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

    # this definition of capping procedure will not cap the top of the voxel
    bottomhalf_caps = [0, cap_geo_bottom, cap_geo_right, cap_geo_left, cap_geo_back, cap_geo_front]
    tophalf_caps = [cap_geo_top, 0, cap_geo_right, cap_geo_left, cap_geo_back, cap_geo_front]

    structure = hybrid_codedstructure(template_1, pitch, [one_voxel, two_voxel, three_voxel],
                                      [capmesh, bottomhalf_caps, tophalf_caps])
    #auto_name=generate_ct_file_name(sw, cf, pitch, '_w_holes')
    auto_name = generate_ct_file_name(sw, cf, pitch)
    structure.save('generated_stl_files/ct_specimens/' + auto_name)

if __name__ == "__main__":
    main()


