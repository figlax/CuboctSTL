from CuboctSTL_v0 import *


def hybrid_template1():
    x = 10
    y = 10
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



def main():
    # The following code is an example of using the hybrid structure code
    strut_width = .600
    chamfer_factor = 5
    pitch = pitch_from_relden(0.15, chamfer_factor, strut_width)

    template = hybrid_template1()

    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    two_voxel = hybrid_voxel(strut_width * 0.75, chamfer_factor, pitch, strut_width)

    capmesh = cap_cuboct(strut_width, chamfer_factor)



    hybrid_structure = hybrid_codedstructure(template, pitch, [one_voxel, two_voxel],
                                             capmesh)
    hybrid_structure.save('generated_stl_files/hetcuboct1_RD0-15_75perReduct_SW0-6_cf5.stl')


if __name__ == "__main__":
    main()


