'''
This version modifies version 5, by reassigning queries to subtree based on minimizing the topolocial tree distance from the nearest leaf by 
Hamming distance in the full tree to the seed sequence from the subtree. If the nearest leaf is not present in the subtree, it is considered 
to have Hamming distance of infinity.
'''

import sys
import os
import utils
import shutil
import json
import time
import argparse
import treeswift
import copy
import multiprocessing as mp
import threading
import itertools
from collections import Counter

def main(args):
    tree_path = args.tree
    output = args.outdir
    outFile = args.output
    aln = args.alignment
    n = args.subtreesize
    run = args.tmpfilenbr
    subtree_flag = args.subtreetype
    fragment_flag = args.fragmentflag
    q_aln = args.qalignment
    model = args.model
    info = args.info
    nbr_closest = args.votes

    # output path, ref, query, backbone tree, info
    t0 = time.perf_counter()
    tree = treeswift.read_tree_newick(tree_path)
    tree.resolve_polytomies()

    leaf_dict = tree.label_to_node(selection='leaves')
    print ('{} seconds loading tree'.format(time.perf_counter() - t0))

    if q_aln != "":
        ref_dict = utils.read_data(aln)
        q_dict = utils.read_data(q_aln)
    else:
        aln_dict = utils.read_data(aln)
        ref_dict, q_dict = utils.seperate(aln_dict, leaf_dict)
    print ('{} seconds loading alignment'.format(time.perf_counter() - t0))

    try:
        os.mkdir("tmp{}".format(run))
    except OSError as error:
    	pass
    try:
        os.mkdir(output)
    except OSError as error:
        pass

    query_votes_dict = dict()
    query_top_vote_dict = dict()
    
    tmp_output = "tmp{}/".format(run) + "/closest.txt"

    if fragment_flag == True:
        os.system("./fragment_hamming {} {} {} {} {} {}".format(aln, len(ref_dict), q_aln, len(q_dict), tmp_output, nbr_closest))
    else:    
        os.system("./hamming {} {} {} {} {} {}".format(aln, len(ref_dict), q_aln, len(q_dict), tmp_output, nbr_closest))
    print ('{} seconds finding closest leaves'.format(time.perf_counter() - t0))
    #for name, seq in q_dict.items():
    #    y = utils.find_closest_hamming(seq, ref_dict, 5, fragment_flag)
    #    print ('{} Closest sister taxa found: {}'.format(name, y[0]))
    #    print ('{} seconds new finding closest leaf'.format(time.perf_counter() - t0))
   
    f = open(tmp_output)
    for line in f:
        line = line.strip()
        y = line.split(',')
        name = y.pop(0)

        #print(name, y)
        for idx, taxon in enumerate(y):

            leaf, hamming = taxon.split(':')
            y[idx] = (leaf, int(hamming))

        y = sorted(y, key=lambda x: x[1])
        #print(y)
        for idx, taxon in enumerate(y):
            y[idx] = taxon[0]

        query_votes_dict[name] = y
        query_top_vote_dict[name] = y[0]
    f.close()

    print ('{} seconds processing closest leaves'.format(time.perf_counter() - t0))

    lf_votes = Counter()
    leaf_queries = dict()
    for name, y in query_votes_dict.items():

        lf_votes.update(y)
        for ind, leaf in enumerate(y):
            top_vote = False
            if ind == 0:
                top_vote = True
            if leaf not in leaf_queries:           
                leaf_queries[leaf] = {(name,top_vote)}
            else:
                leaf_queries[leaf].add((name,top_vote))

    print (len(leaf_queries))
    subtree_dict = dict()
    subtree_leaf_label_dict = dict()
    nbr_of_queries = len(q_dict)
    print (len(query_votes_dict), nbr_of_queries)
    most_common_index = 0
    strike_count = 0
    while len(query_votes_dict) > 0:
        (seed_label, node_votes) = lf_votes.most_common(most_common_index+1)[most_common_index]
        
        node_y = leaf_dict[seed_label]
        labels = utils.subtree_nodes_with_edge_length(tree, node_y, n)
        subtree = tree.extract_tree_with(labels)
        label_set = set(labels)

        queries_by_subtree = set()
        subtree_query_set = set()

        #gather any other queries that can be used with this subtree
        for label in labels:
            leaf_queries_remove_set = set()
            if label in leaf_queries:
                    
                for leaf_query, top_vote in leaf_queries[label]:
                
                    if leaf_query not in query_votes_dict:
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                        continue
                        
                    if top_vote:
                        subtree_query_set.add(leaf_query)
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                    
                leaf_queries[label].difference_update(leaf_queries_remove_set)
        queries_by_subtree.update(subtree_query_set)

        if len(queries_by_subtree)> 0:
            subtree_dict[subtree] = (seed_label, queries_by_subtree)
            subtree_leaf_label_dict[subtree] = subtree.label_to_node(selection='leaves')
            print ("{} queries in subtree".format(len(queries_by_subtree)))
        votes_b4 = len(list(lf_votes.elements()))
        
        for query in queries_by_subtree:
            if query in query_votes_dict:
                lf_votes.subtract(query_votes_dict[query])
                query_votes_dict.pop(query)
        if len(queries_by_subtree)> 0:        
            print ("votes before: {}, votes after: {}".format(votes_b4, len(list(lf_votes.elements()))))
            print ("queries left: {}".format(len(query_votes_dict)))
        if len(queries_by_subtree) == 0:
            most_common_index += 1
        else:
            most_common_index = 0;
            
    jplace = dict()
    placements = []

    utils.add_edge_nbrs(tree)
    jplace["tree"] = utils.newick_edge_tokens(tree)
    print ('{} seconds adding tokens'.format(time.perf_counter() - t0))

    placed_query_list = []
    
    new_subtree_dict = dict()
    for query, closest_label in query_top_vote_dict.items():
   
        best_subtree = None
        best_distance = 99999999999999999
        for subtree, value in subtree_dict.items():
            leaf_label_dict = subtree_leaf_label_dict[subtree]
            seed_label, _ = value
            if closest_label in leaf_label_dict:
                distance = subtree.distance_between(leaf_label_dict[closest_label], leaf_label_dict[seed_label])
                if distance < best_distance:
                   best_subtree = subtree
                   best_distance = distance
        if best_subtree in new_subtree_dict:
            new_subtree_dict[best_subtree].append(query)
        else:
            new_subtree_dict[best_subtree] = [query]
        
            

    print ('{} seconds assigning subtrees'.format(time.perf_counter() - t0))
    final_subtree_count = 0
    
    tmp_output = "tmp{}/query_subtree_dict".format(run)
    
    for subtree, query_list in new_subtree_dict.items():

        if len(query_list) == 0:
            continue

        final_subtree_count += 1
        
        tmp_tree = "tmp{}/tree".format(run)
        tmp_aln = "tmp{}/aln".format(run) + ".fa"
        tmp_qaln = "tmp{}/qaln".format(run) + "q.fa"
        tmp_output = "tmp{}/".format(run) + "/epa_result.jplace"
        tmp_dir = "tmp{}/".format(run)
        try:
            os.mkdir(tmp_dir)
        except OSError as error:
            pass

        f = open(tmp_aln, "w")
        fq = open(tmp_qaln, "w")
        for name in query_list:
            fq.write(">"+name+"\n")
            fq.write(q_dict[name]+"\n")
        tmp_leaf_dict = subtree.label_to_node(selection='leaves')
        for label, _ in tmp_leaf_dict.items():
            if label != '':
                f.write(">"+label+"\n")
                f.write(ref_dict[label]+"\n")

        f.close()
        fq.close()

        #print ('{} seconds building alignment'.format(time.perf_counter() - t0))

        subtree.resolve_polytomies()
        subtree.suppress_unifurcations()
        subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

        #print ('{} seconds writing subtree'.format(time.perf_counter() - t0))

        os.system("./epa-ng -m {} -t {} -w {} -s {} -q {} --redo -T 16".format(info, tmp_tree, tmp_dir, tmp_aln, tmp_qaln))

        #print ('{} seconds running epa-ng'.format(time.perf_counter() - t0))

        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)
        
        
        if len(place_json["placements"]) > 0:

            added_tree, edge_dict = utils.read_tree_newick_edge_tokens(place_json["tree"])

            for tmp_place in place_json["placements"]:
                placed_query_list.append(tmp_place["n"][0])
                for i in range(len(tmp_place["p"])):

                    edge_num = tmp_place["p"][i][0]
                    edge_distal = tmp_place["p"][i][3]

                    #print(edge_num)
                    #print(edge_distal)
                    #print ('{} loading jplace'.format(time.perf_counter() - t0))

                    right_n = edge_dict[str(edge_num)]
                    left_n = right_n.get_parent()

                    #left and right path_l and path_r are in added_tree
                    left, path_l = utils.find_closest(left_n, {left_n, right_n})
                    right, path_r = utils.find_closest(right_n, {left_n, right_n})

                    left = leaf_dict[left.get_label()]
                    right = leaf_dict[right.get_label()]
                    _, path = utils.find_closest(left, {left}, y=right)
                    # now left right and path are in tree

                    length = sum([x.get_edge_length() for x in path_l])+edge_distal
                    # left path length through subtree before placement node

                    target_edge = path[-1]

                    for j in range(len(path)):
                        length -= path[j].get_edge_length()
                        if length < 0:
                            target_edge = path[j]
                            break

                    tmp_place["p"][i][0] = 0

                    #print (tree)
                    label = target_edge.get_label()
                    
                    #print (label)
                    [taxon, target_edge_nbr] = label.split('%%',1)
                    tmp_place["p"][i][0] = target_edge.get_edge_length()+length
                    tmp_place["p"][i][1] = int(target_edge_nbr)
                    
                    #print(target_edge.get_edge_length()+length)

                    #print(int(target_edge_nbr))

                placements.append(tmp_place.copy())

        place_file.close()
    
    print(len(placed_query_list))

    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]


    output = open('{}/{}.jplace'.format(output,outFile), 'w')
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    print ('{} seconds building jplace'.format(time.perf_counter() - t0))
    print ('Final number of subtrees used:', final_subtree_count)
    #shutil.rmtree("tmp{}".format(run))
    


def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--info", type=str,
                        help="Path to model parameters", required=True, default=None)
    
    parser.add_argument("-t", "--tree", type=str,
                        help="Path to reference tree with estimated branch lengths", required=True, default=None)
    
    parser.add_argument("-d", "--outdir", type=str,
                        help="Directory path for output", required=True, default=None)
    
    parser.add_argument("-a", "--alignment", type=str,
                        help="Path for query and reference sequence alignment in fasta format", required=True, default=None)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file name", required=False, default="EPA-ng-BXR")
    
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")

    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=2000)
    
    parser.add_argument("-V", "--votes", type=int,
                        help="Integer number of votes per query sequence",
                        required=False, default=5)
    
    parser.add_argument("-s", "--subtreetype", type=str,
                        help="d (default) for edge weighted distances, n for node distances, h for hamming distances",
                        required=False, default="b")
    
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    
    parser.add_argument("-f", "--fragmentflag", type=str2bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=True)

    parser.add_argument("-v", "--version", action="version", version="1.0.0", help="show the version number and exit")
                       
    return parser.parse_args()


def str2bool(b):
    if isinstance(b, bool):
       return b
    if b.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif b.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')        

if __name__ == "__main__":
    main(parseArgs())
