from CuboctSTL_v0 import *


def main():
    # The following code is an example of using the hybrid structure code
    strut_width = 0.600
    chamfer_factor = 5
    rel=0.03
    pitch = pitch_from_relden(rel, chamfer_factor, strut_width)
    x_vox = 10
    y_vox = 10
    z_vox = 10


    capmesh = cap_cuboct(strut_width, chamfer_factor)

    compression = compression_specimen(strut_width, chamfer_factor, pitch, x_vox, y_vox, z_vox)
    auto_file_name = generate_file_name(strut_width, chamfer_factor, x_vox, y_vox, z_vox, pitch, rel, half='yes', extra_text='test')
    compression.save("generated_stl_files/" + auto_file_name)


if __name__ == "__main__":
    main()