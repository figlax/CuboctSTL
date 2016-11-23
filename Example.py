from CuboctSTL_v0 import *

def main():
    pitch = 30
    strut_width = 2
    chamfer_factor = 2.75
    x=3
    y=3
    z=3

    lattice = make_lattice(strut_width, chamfer_factor, pitch, x, y, z)
    lattice.save('test_lattice.stl')

    template = create_test_template()
    capmesh = cap_cuboct(strut_width, chamfer_factor)
    one_voxel = voxel(strut_width, chamfer_factor, pitch)
    structure = lattice_codedstructure(one_voxel, capmesh, pitch, template)
    structure.save('coded_structure_test.stl')

    preview_mesh(lattice) #This renders a 3D plot of the geometry

if __name__ == "__main__":
    main()