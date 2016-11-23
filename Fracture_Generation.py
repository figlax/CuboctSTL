from CuboctSTL_v0 import *

def ct_template2():

    template = np.zeros((25, 10, 25), dtype=np.int)

    template[:, :, 0:12] = 1
    template[15:25, :, 12] = 1
    template[:, :, 13:25] = 1


    return template

def main():
    pitch = 6.04
    sw = 0.6
    cf = 3

    template = ct_template2()
    # voxel_instance = voxel(sw, cf, pitch)
    # cap_1 = cap_cuboct(sw, cf)
    # structure = lattice_codedstructure(voxel_instance, cap_1, pitch, template )
    # structure.save('cuboct_ct_base.stl')

    one_voxel = voxel(sw, cf, pitch)
    two_voxel = hybrid_voxel(sw * 0.25, cf, pitch, sw)
    template_1 = hybrid_template()
    cap_1 = cap_cuboct(sw, cf)

    structure = hybrid_codedstructure(template_1, pitch, [two_voxel, one_voxel], cap_1)

    structure.save('generated_stl_files/cuboct_ctspecimen_2.stl')

if __name__ == "__main__":
    main()


