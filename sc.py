# -*- coding: utf-8 -*-
# @Author: weijie
# @Date:   2021-08-29 23:43:59
# @Last Modified by:   weijie
# @Last Modified time: 2021-11-24 18:30:09

import math
import numpy as np
from scipy.sparse.linalg import eigs
from sklearn.cluster import *
import sys

def read_file(file_name):
	"""Read protein opscan_output files
	There are two main indicators to indicate the distances of proteins, one is the similarity between proteins, and the other is the length ratio of amino acid sequences of proteins
	"""
	x = []
	y = []
	with open(file_name, "r", encoding = "utf-8") as fin:
		for line in fin.readlines():
			line = line.split()
			y.append(float(line[1]))
			x.append(float(line[2]))
	return x, y

def standadization(x):
	"""The method implemented in this example is zero-mean normalization
	"""
	n = len(x)
	avg = 0.0
	for i in range(n):
		avg = avg + x[i]
	avg = avg / n
	for i in range(n):
		x[i] = x[i] - avg
	sd = 0.0
	for i in range(n):
		sd = sd + x[i] * x[i]
	sd = sd / n
	sd = math.sqrt(sd)
	for i in range(n):
		x[i] = x[i] / sd

def euclidean_distance(x1, y1, x2, y2):
	'''Euclidean distance formula
	'''
	d = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
	return d

def build_full_connected_graph(x, y, n, sigma):
	"""Construct similarity matrix-adjacency matrix
	"""
	A = np.zeros((n,n))
	e = 2.718281828
	for i in range(n):
		for j in range(i, n):
			dist = euclidean_distance(x[i], y[i], x[j], y[j])
			A[i, j] = pow(e, -dist / (2.0 * sigma * sigma))
			A[j, i] = A[i, j]
	return A

def calculate_eigenvalue(A, max_eigenvalue_count):
	"""Compute the laplacian matrix eigenvalues and eigenvectors, which we usually do here for sparse matrices
	"""
	degreeMatrix = np.diag(1.0 / np.sum(A, axis=1))
	W_norm = np.dot(degreeMatrix, A)

	L_vals, L_vecs = eigs(W_norm, max_eigenvalue_count)
	L_vals = L_vals.real
	L_vecs = L_vecs.real
	L_vecs = L_vecs[:,np.argsort(L_vals)]
	L_vals = L_vals[np.argsort(L_vals)]
	L_vals = L_vals.tolist()
	L_vals.sort(reverse = True)
	return L_vals, L_vecs

def initial_step(file_name):
	'''Initialization steps, read data, normalization, and etc.
	'''
	x, y = read_file(file_name)
	original_x = x.copy()
	original_y = y.copy()
	n = len(x)
	standadization(x)
	standadization(y)
	x = np.array(x)
	y = np.array(y)
	data = np.c_[x, y]
	max_eigenvalue_count = 50
	return x, y, data, original_x, original_y, max_eigenvalue_count, n

def get_clusters_k(eigengap):
	'''According to the threshold to determine the number of clusters k, here after experiments to adjust the parameter is 1e-3, others can change this parameter according to their needs
	'''
	k = 2
	for i in range(len(eigengap)):
		if eigengap[i] >= 1e-3:
			k = i + 1
			break
	if k < 2:
		k = 2
	return k

def determin_signal_noise(k, n, labels, original_x, original_y, extra_limitation = False):
	'''Based on the clustering results, determine the signal and noise proteins
	'''
	cluster = [[] for i in range(k)]
	for i in range(n): 
		cluster[labels[i]].append(original_y[i])
	y_min = [min(cluster[i]) for i in range(k)]

	k = len(set(labels))
	cluster = [[] for i in range(k)]
	for i in range(n):
		cluster[labels[i]].append(original_y[i])
	size = [len(cluster[i]) for i in range(k)]
	y_centroid = [np.median(cluster[i]) for i in range(k)]
	maxSizeLabel = size.index(max(size))

	# 1 for signal, 0 for noise
	categorys = [0] * n

	# A cluster with a center of mass greater than the maximum cluster size is a signal, and the opposite is a nois
	for i in range(n):
		if y_centroid[labels[i]] > y_centroid[maxSizeLabel]:
			categorys[i] = 1

	# Find the x, y values of the highest point of noise
	noise_y = [original_y[i] for i in range(n) if categorys[i] == 0]
	noise_x = [original_x[i] for i in range(n) if categorys[i] == 0]
	noise_y_max = max(noise_y)
	cur_noise_x = noise_x[noise_y.index(noise_y_max)]

	# All points to the right of the highest point are noise
	for i in range(n):
		if original_y[i] <= noise_y_max and original_x[i] >= cur_noise_x:
			categorys[i] = 0
	
	return labels, categorys

def clustering_kernel_step(k, n, eigenvectors, max_eigenvalue_count, original_x, original_y):
	'''Clustering steps
	'''
	H = eigenvectors[:, (max_eigenvalue_count-k):(max_eigenvalue_count)]
	# KMeans
	labels = KMeans(n_clusters = k).fit(H).labels_
	labels, categorys = determin_signal_noise(k, n, labels, original_x, original_y)
	return labels, categorys

def partition(file_name, x, y, data, original_x, original_y, n, max_eigenvalue_count, sigma):
	"""The laplacian matrix is constructed and the eigenvalues and eigenvectors are calculated for multiple classification.
	"""
	A = build_full_connected_graph(x, y, n, sigma)
	eigenvalues, eigenvectors = calculate_eigenvalue(A, max_eigenvalue_count)
	# Find the eigengap
	eigengap = [abs(eigenvalues[i] - eigenvalues[i+1]) for i in range(len(eigenvalues) - 1)]

	# The number of clusters k is obtained according to a threshold value, and then signal and noise are determined
	k = get_clusters_k(eigengap)
	labels, categorys = clustering_kernel_step(k, n, eigenvectors, max_eigenvalue_count, original_x, original_y)
		
	return labels, categorys

def output_result(infile, categorys, categorys_outfile):
	'''Get the output result
	'''
	with open(infile, "r", encoding = 'utf-8') as fin:
		with open(categorys_outfile, "w", encoding = 'utf-8') as fout:
			i = 0
			for line in fin:
				fout.write(line.rstrip() + " " + str(categorys[i]) + '\n')
				i = i + 1

if __name__ == '__main__':
	
	# Input file for x, y data
	file_path = sys.argv[1]
	file_name = file_path.split("/")[-2] + '/' + file_path.split("/")[-1]
	# print(file_path.split("/")[-1])
	x, y, data, original_x, original_y, max_eigenvalue_count, n = initial_step(file_path)

	# Use of full connectivity diagrams
	# sigma is the parameter of the Gaussian function
	sigma = 0.2
	labels, categorys = partition(file_name, x, y, data, original_x, original_y, n, max_eigenvalue_count, sigma)

	output_result(file_path, categorys, sys.argv[2])
