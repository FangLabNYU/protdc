import sys

# Calculate the degree of the points in a network

# Input files, network edge file
network = sys.argv[1]
# Output file, degree of each node after normalized
degree_file = sys.argv[2]

OG_name = network.split('/')[-1]

with open(degree_file, 'a', encoding = 'utf-8') as fout:
	with open(network, 'r', encoding = 'utf-8') as fin:
		protein_weight_sum, edges = {}, set()
		for line in fin:
			line = line.rstrip().split()
			# In our calculations, the edges are undirected
			if (line[0], line[1]) not in edges and (line[1], line[0]) not in edges:
				edges.add((line[0], line[1]))
				# Get the sum of the weights of the edges connected to each point
				if line[0] not in protein_weight_sum:
					protein_weight_sum[line[0]] = 0
				protein_weight_sum[line[0]] = protein_weight_sum[line[0]] + float(line[2])
				if line[1] not in protein_weight_sum:
					protein_weight_sum[line[1]] = 0
				protein_weight_sum[line[1]] = protein_weight_sum[line[1]] + float(line[2])

		size = len(list(protein_weight_sum.keys()))
		if size >= 2:
			# Note that the use of "size-1" here, rather than "size", prevents the points from being particularly small 
			degrees = [round(v/(size-1), 4) for v in list(protein_weight_sum.values())]
			for p, v in protein_weight_sum.items():
				fout.write(p + '\t' + str(round(v/(size-1), 4)) + '\n')
