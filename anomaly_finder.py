#! /usr/bin/python

"""
This is just a python3 interface for the anomaly_finder.py script written
by CW Linkem (https://github.com/cwlinkem/anomaly_zone)

If you use this, please cite the original publication:
Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone 
in species trees and evidence for a misleading signal in higher-level skink phylogeny 
(Squamata: Scincidae). Systematic Biology, 65(3), 465-477.

"""

import sys
import os
import math
from math import exp
import dendropy
import getopt
from dendropy.utility import bitprocessing

def main():
	params = parseArgs()
	print("\n####################anomaly_finder.py######################\n")
	#read trees
	taxon_set = dendropy.TaxonNamespace()
	print("Reading primary tree...", params.tree)
	mrc = dendropy.Tree.get_from_path(src=params.tree, schema=params.ftype, taxon_namespace = taxon_set)
	if params.btrees:
		print("Reading bootstrap trees...", params.btrees)
		trees = dendropy.TreeList.get_from_path(params.btrees, params.ftype, taxon_namespace = taxon_set)
	else:
		print("No bootstrap trees provided.")

	#initialize some dictionaries
	master_dict=dict()
	mrc_dict=dict()
	split_dict=dict()
	split_occ_dict=dict()

	#if bootstrap trees provided, get results for them
	#master_dict:
	#key = pair of edges
	#value = list of anomaly zone results (1=True; 0=False)
	if params.btrees:
		print("Finding anomalous nodes in bootsteap trees...")
		for sp_tree in trees:
			#dendropy.treesplit.encode_splits(sp_tree) #deprecated
			sp_tree.encode_bipartitions()
			#split_bitmasks = sp_tree.split_bitmask_edge_map.keys()
			master_dict = get_nodes(sp_tree, master_dict)
			split_occ_dict=split_freq(sp_tree, taxon_set, split_occ_dict)	

	#sys.exit(0)
	#get results for main tree
	print("Finding anomalous nodes in primary tree...")
	mrc.encode_bipartitions()
	mrc_dict=get_nodes(mrc, mrc_dict) #get AZ results for primary tree
	split_dict=split_mapper(mrc, taxon_set)  #match a list of taxa to the split pattern from the bitmask for the whole tree

	print("\n######################ENCODING#############################\n")
	
	print("Writing results with splits encoded as:")
	labelTree = getBitStringLabels(mrc, len(taxon_set))
	labelTree.print_plot(show_internal_node_labels=True)
	
	print("(If tree is large, import the following into FigTree)\n")
	print(labelTree.as_string(schema="newick"))

	#if bootstrap trees provided:
	#for each internal edge pair:
	print("\n######################RESULTS##############################\n")
	print("Showing anomalous nodes only:\n")
	if params.btrees:
		print("ParentEdge\tDescendantEdge\tProp.Anomalous")
		seen=False
		for k,v in master_dict.items():
			taxlab1=list()
			taxlab2=list()
			#print(k[0])
			#if edge pair has results for main tree
			if k in mrc_dict.keys():
				#print("k:",k, "v:",v)
				#bitstring for parent node 
				parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
				#print(parentLabel)
				#list of descendent taxa
				parentList = taxon_set.bitmask_taxa_list(k[0])
				#print(parentList)
				#label for descendant edge
				descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
				#list of descendant taxa
				descendantList = taxon_set.bitmask_taxa_list(k[1])
				#print(descendantList)
				for i in parentList:
					taxlab1.append(i.label)
				for j in descendantList:
					taxlab2.append(j.label)
				if sum(v) > 0:
					seen=True
					print(parentLabel, descendantLabel,sum(v)/len(v), sep="\t")
		if seen==False:
			print("No anomalous nodes detected.")
	else:
		print("ParentEdge\tDescendantEdge")
		seen=False
		for k,v in mrc_dict.items():
			taxlab1=list()
			taxlab2=list()
			anom=False
			#print("k:",k, "v:",v)
			#bitstring for parent node 
			parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
			#print(parentLabel)
			#list of descendent taxa
			parentList = taxon_set.bitmask_taxa_list(k[0])
			#print(parentList)
			#label for descendant edge
			descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
			#list of descendant taxa
			descendantList = taxon_set.bitmask_taxa_list(k[1])
			#print(descendantList)
			for i in parentList:
				taxlab1.append(i.label)
			for j in descendantList:
				taxlab2.append(j.label)
			if sum(v) > 0:
				seen=True
				print(parentLabel, descendantLabel, sep="\t")
		if seen==False:
			print("No anomalous nodes detected.")
	print("\n########################OUTPUT############################\n")

	
	if params.btrees:
		print("Tree with anomalous pairs annotated as \"AZ Event : Prop. BS\"")
		outTreeBS = getAnomalyLabelsBS(mrc, master_dict)
		outTreeBS.print_plot(show_internal_node_labels=True)
		
		print("(If tree is large, import the following into FigTree)\n")
		print(outTreeBS.as_string(schema="newick"))
		with open(params.out, 'w') as fh:
			l=outTreeBS.as_string(schema="newick")
			l=l+"\n"
			fh.write(l)
		fh.close()
	else:
		print("Tree with anomalous edge pairs annotated (AZ pairs numbered)")
		outTree = getAnomalyLabels(mrc, mrc_dict)
		outTree.print_plot(show_internal_node_labels=True)
		
		print("(If tree is large, import the following into FigTree)\n")
		print(outTree.as_string(schema="newick"))
		with open(params.out, 'w') as fh:
			l=outTree.as_string(schema="newick")
			l=l+"\n"
			fh.write(l)
		fh.close()

	print("\n########################DONE!############################\n\n")
	
	


#function returns tree with AZ edge pairs labelled
def getAnomalyLabels(tree, d):
	tree2 = tree.clone(depth=1)
	labelDict = dict()
	az_event = 1
	for k,v in d.items():
		#if AZ test positive
		if sum(v) > 0:
			for edge in k:
				if str(edge) in labelDict:
					#print("exists")
					labelDict[str(edge)] = str(labelDict[str(edge)]) + "/" + str(az_event)
				else:
					labelDict[str(edge)] = str(az_event)
			#increment event number
			az_event += 1
	for node in tree2.levelorder_node_iter():
		if node.parent_node:
			lab=str(node.edge.split_bitmask)
			if lab in labelDict:
				l = "\"" + labelDict[str(lab)] + "\""
				node.label=l
			else:
				lab = "\"" + "0" + "\""
				node.label=lab
		else:
			lab = "\"" + "0" + "\""
			node.label=lab
	return(tree2)
		
#function returns tree with AZ edge pairs labelled
def getAnomalyLabelsBS(tree, d):
	tree2 = tree.clone(depth=1)
	labelDict = dict()
	bsDict = dict()
	az_event = 1
	for k,v in d.items():
		#if AZ test positive
		if sum(v) > 0:
			for edge in k:
				if str(edge) in labelDict:
					#print("exists")
					labelDict[str(edge)] = str(labelDict[str(edge)]) + "/" + str(az_event)
					bsDict[str(edge)] = str(bsDict[str(edge)]) + "/" + str(sum(v)/len(v))
				else:
					labelDict[str(edge)] = str(az_event)
					bsDict[str(edge)] = str(sum(v)/len(v))
			#increment event number
			az_event += 1
	for node in tree2.levelorder_node_iter():
		if node.parent_node:
			lab=str(node.edge.split_bitmask)
			if lab in labelDict:
				node.label=labelDict[str(lab)] + ":" + bsDict[str(lab)]
			else:
				node.label="0:0"
		else:
			node.label="0:0"
	return(tree2)
#function returns tree with bitstring node labels 
def getBitStringLabels(tree, l):
	tree2 = tree.clone(depth=1)
	for node in tree2.postorder_node_iter():
		lab = "\"" + str(bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=l)) + "\""
		node.label=lab
	return tree2

#For a pair of edges, return 1 if anomalous; return 0 if not anomalous
def anomaly_calc(nodes):
	lambda1=nodes[0]
	lambda2=nodes[1]
	if lambda1 is not None:
		Z1=math.log(2.0/3+((3*exp(2*lambda1)-2)/(18*(exp(3*lambda1)+0.000000000000001-exp(2*lambda1))))) #function a(x) from Degnan and Rosenberg (2006)
		if lambda2 <= Z1:
			return 1
		else:
			return 0
	else:
		return 0

def split_mapper(tree, taxon_set):
	split_dict={}
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			taxlist=[]
			tax = taxon_set.bitmask_taxa_list(node.edge.split_bitmask)
			for i in tax:
				taxlist.append(i.label)
			split_dict[bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))] = taxlist
	return split_dict

#function
def split_freq(tree,taxon_set, split_occ_dict):
	#for each node 
	for node in tree.postorder_node_iter():
		#if not root
		if node.parent_node is None:
			node.value = 1.0
		else:
			#split=dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set)) #deprecated
			split=bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))
			print(split)
			if split not in split_occ_dict.keys():
				split_occ_dict[split] = [1]
			else:
				split_occ_dict[split].append(1)
	return split_occ_dict
	
#Function tests for internal edge pairs under anomaly zone
#modifies a dict where key = pair of edges; calue = list of results (1=true; 0=false)		
def get_nodes(tree, master_dict):
	pair_list=[]
	#for each internal node
	#calculate
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			#root
			node.value = 1.0
		else:
			#if not root, and internal node
			if node.taxon is None:
				#get edge length of parent, and own edge length
				#x = parent; y = current
				node_pair = (node.parent_node.edge_length, node.edge_length)
				#calculate if pair of edges are under anomaly boundary
				anomalous=anomaly_calc(node_pair)
				#print(anomalous)
				edge_pair = (node.parent_node.edge.split_bitmask, node.edge.split_bitmask)
				#for edge pair, keep anomaly zone tests results
				if edge_pair not in master_dict.keys():
					master_dict[edge_pair] = [anomalous]
				else:
					master_dict[edge_pair].append(anomalous)	
	return master_dict		

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 't:f:b:hco:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.tree=None
		self.ftype="newick"
		self.btrees=None
		self.out="out.tre"


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")
			if o =="-c":
				self.display_citation()

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "t":
				self.tree = arg
			elif opt in ('h', 'help', 'c'):
				pass
			elif opt == "f":
				if arg == "newick" or arg == "nexus":
					self.ftype =arg
				else:
					self.display_help("Invalid option for -f",arg)
			elif opt == "b":
				self.btrees=arg
			elif opt == "o":
				self.out = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.tree:
			self.display_help("Error: Need treefile <-t>")
		if not self.ftype:
			print("No filetype (-f) provided: using default (\"newick\")")



	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nanomaly_finder.py\n")
		print ("\nUsage: ", sys.argv[0], "-t <treefile> -f <nexus or newick> \n")
		print ("""Description: An implementation of the code from Linkem et al (2016)
		
Citation: Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone in species trees and evidence for a misleading signal in higher-level skink phylogeny (Squamata: Scincidae). Systematic Biology, 65(3), 465-477.""")
		print("""
	Arguments:
		-t		: Tree file (branches scaled by coalescent units)
		-b		: Tree file containing bootstrap trees, if being used
		-f		: Format of tree file: "nexus" or "newick"
		-c		: Displays citation information
		-o		: Output file name for annotated tree (default: out.tre)
		-h		: Displays help menu

""")
		sys.exit()
	
	def display_citation(self):
		print("\nCitation: Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone in species trees and evidence for a misleading signal in higher-level skink phylogeny (Squamata: Scincidae). Systematic Biology, 65(3), 465-477.\n")
		
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

		
