import configparser
import SVGParser as parser
from detect_crossings import wirt_presentation
from find_generators import find_min_generating_set
from find_homomorphism import perform_combinatoric
import matrix as m

# first version without relations reduction and parallel executing
__version__ = "1.0.0"

if __name__ == "__main__":
	config = configparser.ConfigParser()
	config.read('config.ini')
	config_filename = config['Settings']['input_file']
	#config_filename = 'trefoil'
	strands = parser.parse_svg(config_filename)
	print("##### Computing Wirtinger Presentation #####")
	wirt_pres = wirt_presentation(strands)
	indices = find_min_generating_set(wirt_pres)
	
	Alex_polynomial = None	
	if int(config['Settings']['compute_AP']):
		Alex_polynomial = m.compute_Jacobi_matrix_for_AP(wirt_pres, 0, 0).det()
		print("##### Alexander polynomial: #####")
		print(Alex_polynomial)
		
	tAlex_polynomial = None
	if int(config['Settings']['compute_TAP']):
		config_group_size = int(config['Settings']['group_size'])
		config_group_type = config['Settings']['group_type']
		print("##### Finding epimorphism #####")
		imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
		#numerator
		if imgs != None:
			det_upper = m.compute_Jacobi_matrix_for_TAP(imgs[0], wirt_pres, 0, 0, config_group_size).det()
			det_lower = m.compute_denominator_matrix(imgs[0][0], config_group_size).det()
			tAlex_polynomial = det_upper / det_lower
			print("##### TAP #####")
			print(tAlex_polynomial)
	
	if config['Settings']['output_file'] != None:
		with open(config['Settings']['output_file'], "w") as f:
			f.write(f'Alexander polynomial of the given knot:\n{Alex_polynomial}\n')
			f.write(f'Twisted Alexander polynomial of the given knot:\n{tAlex_polynomial}\n')
