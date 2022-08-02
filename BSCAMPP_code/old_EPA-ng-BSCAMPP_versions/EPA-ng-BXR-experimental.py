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

    files = []
    dist_mat_files = []
    try:
        os.mkdir("tmp{}".format(run))
    except OSError as error:
    	pass
    try:
        os.mkdir(output)
    except OSError as error:
    	pass

    try:
        tmp_dict = open("tmp{}/dict.json".format(run), 'r')
        sister_taxon_dict = json.load(tmp_dict)
        tmp_dict.close()
    except:
        sister_taxon_dict = dict()
        for name, seq in q_dict.items():
            y = utils.find_closest_hamming(seq, ref_dict, 1, fragment_flag)
            if y[0] in sister_taxon_dict:
                sister_taxon_dict[y[0]].append(name)
            else:
                sister_taxon_dict[y[0]] = [name]
            print ('{} Closest sister taxon found: {}'.format(name, y[0]))
            print ('{} seconds new finding closest leaf'.format(time.perf_counter() - t0))
        output = open("tmp{}/dict.json".format(run), 'w')
        json.dump(sister_taxon_dict,output)
        output.close()
    
    trees, query_decomp_dict = utils.build_subtrees(sister_taxon_dict, leaf_dict, tree, 200, n)

    print ('{} seconds decomposing tree'.format(time.perf_counter() - t0))

    jplace = dict()
    placements = []

    utils.add_edge_nbrs(tree)
    jplace["tree"] = utils.newick_edge_tokens(tree)
    print ('{} seconds adding tokens'.format(time.perf_counter() - t0))


    for tree_index in range(len(trees)):

        tmp_tree = "tmp{}/tree_{}".format(run,tree_index)
        tmp_aln = "tmp{}/aln{}".format(run,tree_index) + ".fa"
        tmp_qaln = "tmp{}/qaln{}".format(run,tree_index) + "q.fa"
        tmp_output = "tmp{}/{}".format(run, tree_index) + "/epa_result.jplace"
        tmp_dir = "tmp{}/{}".format(run, tree_index)
        try:
            os.mkdir(tmp_dir)
        except OSError as error:
            pass

        f = open(tmp_aln, "w")
        fq = open(tmp_qaln, "w")
        for name in query_decomp_dict[tree_index]:
            fq.write(">"+name+"\n")
            fq.write(q_dict[name]+"\n")

        temp_subtree_label_dict = trees[tree_index].label_to_node(selection='leaves')
        for label, node in temp_subtree_label_dict.items():
            if label != '':
                f.write(">"+label+"\n")
                f.write(ref_dict[label]+"\n")

        f.close()
        fq.close()

        print ('{} seconds building alignment'.format(time.perf_counter() - t0))

        subtree = trees[tree_index]

        print ('{} seconds extracting subtree'.format(time.perf_counter() - t0))
        subtree.resolve_polytomies()
        subtree.suppress_unifurcations()
        subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

        print ('{} seconds writing subtree'.format(time.perf_counter() - t0))

        os.system("./epa-ng -m {} -t {} -w {} -s {} -q {} --redo -T 16".format(info, tmp_tree, tmp_dir, tmp_aln, tmp_qaln))

        print ('{} seconds running epa-ng'.format(time.perf_counter() - t0))

        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)

        if len(place_json["placements"]) > 0:

            added_tree, edge_dict = utils.read_tree_newick_edge_tokens(place_json["tree"])

            for tmp_place in place_json["placements"]:
                for i in range(len(tmp_place["p"])):
                    edge_num = tmp_place["p"][i][0]
                    edge_distal = tmp_place["p"][i][3]

                    print(edge_num)
                    print(edge_distal)
                    print ('{} loading jplace'.format(time.perf_counter() - t0))

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
                    print (label)
                    [taxon, target_edge_nbr] = label.split('%%',1)
                    tmp_place["p"][i][0] = target_edge.get_edge_length()+length
                    tmp_place["p"][i][1] = int(target_edge_nbr)

                    print(target_edge.get_edge_length()+length)

                    print(int(target_edge_nbr))

                placements.append(tmp_place.copy())

        place_file.close()

    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]

    fragment = ""
    if fragment_flag == True:
        fragment = "-f"

    output = open('{}/{}.jplace'.format(output,outFile), 'w')
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    print ('{} seconds building jplace'.format(time.perf_counter() - t0))
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
                        required=False, default=10000)
    
    parser.add_argument("-s", "--subtreetype", type=str,
                        help="d (default) for edge weighted distances, n for node distances, h for hamming distances",
                        required=False, default="b")
    
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    
    parser.add_argument("-f", "--fragmentflag", type=bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=False)

    parser.add_argument("-v", "--version", action="version", version="1.0.0", help="show the version number and exit")
                       
    return parser.parse_args()

if __name__ == "__main__":
    main(parseArgs())
