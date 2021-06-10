import sys
from shapely.geometry import Polygon
import json, os


def calculate_polyarea(post_processing, coordinates):
    # temp_coor = os.path.join(post_processing, 'temp_coor.csv')

    # with open(temp_coor) as f:
    #     coordinates = json.load(f)

    # print(coordinates)

    polygon = Polygon(coordinates)  # create polygon based on Point ID on boundary
    poly_area = float(polygon.area)  # calculate area of the polygon

    return poly_area

# Function chooser
func_arg = {"-polyarea": calculate_polyarea}

# print coordinates1
if __name__ == "__main__":
    func_arg[sys.argv[1]](sys.argv[2], sys.argv[3])
