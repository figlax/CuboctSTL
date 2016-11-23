from CuboctSTL_v0 import *

def main():
    pitch = 30
    strut_width = 2
    chamfer_factor = 5
    x=3
    y=3
    z=3

    one_voxel = half_voxel(strut_width, chamfer_factor, pitch)
    two_voxel = half_voxel(strut_width, chamfer_factor, pitch)
    two_voxel.rotate([1,0,0], math.radians(180))
    translate(two_voxel, np.array([0,0, pitch]))


    preview_mesh(one_voxel, two_voxel) #This renders a 3D plot of the geometry

if __name__ == "__main__":
    main()