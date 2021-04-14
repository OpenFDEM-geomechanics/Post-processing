import os
import sys
import math
import formatting_codes
from PIL import Image

'''
Stitching Photos together
Merge two images into one, displayed side by side
# https://stackoverflow.com/questions/10657383/stitching-photos-together
    Inputs:
        - folder: path for main folder
        - file1: path to first image file
        - file2: path to second image file
        - tstep: Time step No.
        - fin_Tstep: check if Final Time step
    Returns:
        - the merged Image object
        - Save Screen shot
'''


def merge_images(folder, file1, file2, tstep, fin_tstep):
    # print("Loading Sub Routine")
    # Load both images
    image1 = Image.open(file1)
    image2 = Image.open(file2)
    # Get size both images
    (width1, height1) = image1.size
    (width2, height2) = image2.size
    # Get ratio between the images
    ratio = float(height1) / float(height2)
    newsize =(int(width2 * ratio), int(height1))
    # Resize the images
    image2_resize = image2.resize(newsize)
    (width2, height2) = image2_resize.size
    # Create blank image with dimensions equal to combined image
    result_width = width1 + width2
    result_height = max(height1, height2)
    # Paste images in the newly created image
    result = Image.new('RGB', (result_width, result_height))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2_resize, box=(width1, 0))

    # Stitched folder
    Stitched = os.path.join(folder, "Stitched")
    # Create subfolder
    if not os.path.exists(Stitched):  # Check to see if the folder exists
        os.makedirs(Stitched)  # if not then makes the folder
    # Save Stitched shot to folder
    result.save(os.path.join(Stitched, ("Time_Step_" + str("%0.4d" % int(float(tstep))) + ".png")))

    '''
    # In case convert locks up your system
    # https://github.com/phw/peek/issues/112
    locate this file => /etc/ImageMagick-6/policy.xml
    < policy domain = "resource" name = "memory" value = "2GiB" / >
    < policy domain = "resource" name = "disk" value = "1GiB" / >
    '''

    # At final timestep
    # Create an animated gif from the stitched images
    if fin_tstep == "FinalFrame":
        # print(formatting_codes.red_text("\nCreating Animation"))
        os.system("convert -delay 100 -dispose Background " + str(Stitched) + "/*.png -loop 1 " + str(
            Stitched) + "/animation.gif")
    return result

if __name__ == '__main__':
    # print(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    merge_images(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])