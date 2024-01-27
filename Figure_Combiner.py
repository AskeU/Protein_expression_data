from PIL import Image
import os

# Folder path
folder_path = r'C:\Users\askung\OneDrive - Danmarks Tekniske Universitet\WP1\NGS\2023_11_09_p85a_copy'

# File paths
image1_path = os.path.join(folder_path, '2new_cytosolic.png')
image2_path = os.path.join(folder_path, '2new_full_cell.png')
image3_path = os.path.join(folder_path, '2new_AUC.png')

# Load the PNG images
image1 = Image.open(image1_path)
image2 = Image.open(image2_path)
image3 = Image.open(image3_path)

# Check the sizes of the images to ensure compatibility
width1, height1 = image1.size
width2, height2 = image2.size
width3, height3 = image3.size

# Ensure all images have the same height
if height1 != height2 or height1 != height3:
    raise ValueError("Image heights do not match")

# Set the truncation widths for each image (individual for both sides)
truncate_width1_left = 250
truncate_width1_right = 900

truncate_width2_left = 450
truncate_width2_right = 250

truncate_width3_left = 250
truncate_width3_right = 150

# Calculate the combined width and height
combined_width = (width1 - truncate_width1_left - truncate_width1_right) + \
                 (width2 - truncate_width2_left - truncate_width2_right) + \
                 (width3 - truncate_width3_left - truncate_width3_right)
combined_height = height1

# Create a new image with the combined width and same height
combined_image = Image.new("RGBA", (combined_width, combined_height))

# Paste the images side by side onto the new image with truncation
combined_image.paste(image1.crop((truncate_width1_left, 0, width1 - truncate_width1_right, height1)), (truncate_width1_left, 0))
combined_image.paste(image2.crop((truncate_width2_left, 0, width2 - truncate_width2_right, height2)),
                     (width1 - truncate_width1_right, 0))
combined_image.paste(image3.crop((truncate_width3_left, 0, width3 - truncate_width3_right, height3)),
                     (width1 + width2 - truncate_width1_right - truncate_width2_right-truncate_width2_left-truncate_width1_left+300, 0))

# Save the combined image
combined_image.save(os.path.join(folder_path, 'combined_image.png'))
