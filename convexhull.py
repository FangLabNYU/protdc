import sys
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt

# Draw a convex hull in the area where the noise is located

def read_clustering_result(result_file):
	'''Get the clustering result information and save it to the list
	'''
	x, y, categorys = [], [], []
	with open(result_file, 'r', encoding = 'utf-8') as fin:
		for line in fin:
			line = line.rstrip().split()
			x.append(float(line[3]))
			y.append(float(line[2]))
			categorys.append(int(line[4]))
	return x, y, categorys

def draw_convex_hull(file_name, categorys, x, y):
	'''Draw convex hull according to noise region
	'''
	# Get the seed, signal and noise information respectively
	seed_x, seed_y = [], []
	signal_x, signal_y = [], []
	noise_x, noise_y = [], []
	for i in range(len(x)):
		if int(categorys[i]) == 1:
			if x[i] == 1.0 and y[i] == 100.0:
				seed_x.append(x[i])
				seed_y.append(y[i])
			else:
				signal_x.append(x[i])
				signal_y.append(y[i])
		else:
			noise_x.append(x[i])
			noise_y.append(y[i])

	noise_data = np.column_stack((noise_x, noise_y))
	signal_data = np.column_stack((signal_x, signal_y))
	# Generate convex hull
	hull = ConvexHull(noise_data)				

	# draw it
	plt.scatter(seed_x, seed_y, c = 'red', label = 'Seed')
	plt.scatter(signal_x, signal_y, c= 'orange', label = 'Signal')
	plt.scatter(noise_x, noise_y, c = 'grey', label = 'Noise')

	for simplex in hull.simplices:
		plt.plot(noise_data[simplex, 0], noise_data[simplex, 1], color = 'blue')

	plt.title(file_name)
	plt.xlabel('Length ratio')
	plt.ylabel('Similarity (%)')
	plt.legend()
	plt.show()

if __name__ == '__main__':
	
	# Clustering results file
	result_file = sys.argv[1]
	file_name = result_file.split('/')[-1]
	x, y, categorys = read_clustering_result(result_file)
	draw_convex_hull(file_name, categorys, x, y)

