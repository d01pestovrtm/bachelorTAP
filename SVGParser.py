import xml.etree.ElementTree as ET
import numpy as np

#TODO as class?
#TODO check stroke - different colors to different links
#TODO int or float????
def parse_svg(svg_file):
    """
    Parse the SVG file and extract polylines in the same order they were given.
    """
    tree = ET.parse(svg_file)
    root = tree.getroot()

    # SVG namespace handling
    namespace = {'svg': 'http://www.w3.org/2000/svg'}

    # Extract all polylines
    polylines = []
    for polyline in root.findall('.//svg:polyline', namespace):
        #points is a string
        points = polyline.attrib.get('points', '').strip()
        
        #Maybe imnplement ist using np.fromstring
        #points_array = np.fromstring(points.replace(' ', ','), sep=',', dtype=float).reshape(-1, 2)
        
	#int or float?????
        points = np.array(
            [tuple(map(float, point.split(','))) for point in points.split()]
        )
        polylines.append(points)

    return polylines


