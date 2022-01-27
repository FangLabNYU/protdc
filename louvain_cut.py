# -*- coding: utf-8 -*-
# @Author: weijie
# @Date:   2021-05-17 08:54:46
# @Last Modified by:   weijie
# @Last Modified time: 2021-06-03 16:21:01

import community
import networkx as nx
import sys
import math
import string

def build_graph_weighted(filename):
	'''Read network files to form Graph
	'''
	G = nx.Graph()
	with open(filename, "r", encoding = "utf-8") as fin:
		edges = [line.rstrip().split() for line in fin]
	for edge in edges:
		G.add_edge(edge[0], edge[1], weight = float(edge[2]))
	return edges, G

def write_result(cut_result_path, file_name, cut_information_path, partition, edges, N):
	'''Exporting classification results to different files
	'''
	hash_table = {}
	for edge in edges:
		start_label = partition[edge[0]]
		end_label = partition[edge[1]]
		# If two points belong to the same classification
		if start_label == end_label:
			# After the louvain algorithm cut, the weights of the edges are recalculated, using the Jaccard Index
			# Note: When calculating the signals of a point, the point "SEED" is included in the calculation here
			if edge[0] not in hash_table:
				hash_table[edge[0]] = set()
				hash_table[edge[0]].add(edge[0])
			hash_table[edge[0]].add(edge[1])
			if edge[1] not in hash_table:
				hash_table[edge[1]] = set()
				hash_table[edge[1]].add(edge[1])
		# Write data that does not belong to the same category to a different directory
		else:
			start_network_name = str(start_label)
			end_network_name = str(end_label)
			with open(cut_information_path + file_name, "a", encoding = "utf-8") as fout:
				fout.write("(" + edge[0] + "," + edge[1] + "," + edge[2] + ")" + " " + file_name + "_" + start_network_name + " " + file_name + "_" + end_network_name + "\n")
	
	# Write data belonging to the same category to a single file
	for edge in edges:
		start_label = partition[edge[0]]
		end_label = partition[edge[1]]
		if start_label == end_label:
			network_name = str(start_label)
			# compute Jaccard Index
			signals_s = hash_table[edge[0]]
			signals_e = hash_table[edge[1]]
			union_l = len(signals_s.union(signals_e))
			intersection_l = len(signals_s.intersection(signals_e))
			jaccard = round(intersection_l / union_l, 4)
			if jaccard > 0:
				with open(cut_result_path + file_name + "_" + network_name, "a", encoding = "utf-8") as fout:
					fout.write(edge[0] + " " + edge[1] + " " + str(jaccard) + "\n")

if __name__ == '__main__':

	# Graph generation by edge table
	edges, G = build_graph_weighted(sys.argv[1])

	# Use python_louvain algorithm to get community divisions
	partition = community.best_partition(G, resolution = 1)

	# print(partition)

	file_name = sys.argv[1].split("/")[-1]
	cut_result_path = sys.argv[2]
	cut_information_path = sys.argv[3]

	N = len(set(list(partition.values())))
	# print(N)
	# Divided into one category, it means the louvain algorithm can not be cut; otherwise output to the next directory
	if N > 1:
		write_result(cut_result_path, file_name, cut_information_path, partition, edges, N)
