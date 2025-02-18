from PIL import Image
import os

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
run_name = 'multicomp-reactions/2023-06-19-run01/'

def combine_images(columns, space, images, output_id):
    rows = len(images) // columns
    if len(images) % columns:
        rows += 1
    width_max = max([Image.open(image).width for image in images])
    height_max = max([Image.open(image).height for image in images])
    background_width = width_max*columns + (space*columns)-space
    background_height = height_max*rows + (space*rows)-space
    background = Image.new('RGBA', (background_width, background_height), (255, 255, 255, 255))
    x = 0
    y = 0
    for i, image in enumerate(images):
        img = Image.open(image)
        x_offset = int((width_max-img.width)/2)
        y_offset = int((height_max-img.height)/2)
        background.paste(img, (x+x_offset, y+y_offset))
        x += width_max + space
        if (i+1) % columns == 0:
            y += height_max + space
            x = 0
    background.save(data_folder + run_name + f'results/interpolation_plots_combined/interpolation_across_catalyst_{output_id}.png')

all_ids = range(154)

# group images by 12
for i in range(0, len(all_ids), 12):
    image_list = [data_folder + run_name + f'results/interpolation_plots/interpolation_across_catalyst_{index}.png'
                  for index in all_ids[i:i+12]]
    combine_images(columns=3, space=1, images=image_list, output_id=i//12)