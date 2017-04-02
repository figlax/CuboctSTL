from CuboctSTL_v0 import *


def main():
    pitch = 15
    strut_width = 0.7
    chamfer_factor = 5
    x = 1
    y = 1
    z = 1

    lattice = make_lattice(strut_width, chamfer_factor, pitch, x, y, z)
    lattice.save('test_lattice.stl')

    # template = create_test_template()
    template = np.zeros((1, 1, 1), dtype=np.int)
    template[0, 0, 0] = 1
    capmesh = cap_cuboct(strut_width, chamfer_factor)
    one_voxel = voxel(strut_width, chamfer_factor, pitch)

    structure = lattice_codedstructure(one_voxel, capmesh, pitch, template)
    structure.save('coded_structure_test_onevox.stl')

    preview_mesh(structure)  # This renders a 3D plot of the geometry
    print(3/2)
    print('Words')

if __name__ == "__main__":
    main()