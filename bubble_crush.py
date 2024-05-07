#parse a genomic assembly in GFA format and pop the bubbles
#author: Roland Faure

date = "2020-05-07"
version = "0.1"

import argparse
import sys
import shutil

def parse_args() :
    parser = argparse.ArgumentParser(description='Parse a genomic assembly in GFA format and pop the bubbles')

    parser.add_argument("-i", "--input-assembly", help="Original assembly in GFA format (required)", required=True)
    parser.add_argument("-o", "--output-assembly", help="Output assembly with bubbles crushed (required)", required=True)

    parser.add_argument("-a", "--absolute-threshold", help="Branch of a bubble less covered than this will be deleted. 0 to turn off. [5]", default=5, type=float)
    parser.add_argument("-r", "--relative-threshold", help="Branch of a bubble less covered than this compared to other will be deleted. 0 to turn off. 1 for crushing the less covered branch. [0]", default=0.5, type=float)

    parser.add_argument("-m", "--dont_merge", help="Do not merge the contigs after bubble popping", default=False, action="store_true")
    parser.add_argument("-t", "--transfer-coverage", help="Transfer the coverage of the deleted branch to the other branch. [False]", action="store_true")
    
    return parser.parse_args()


def reverse_complement(seq) :
    comp = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    return "".join([comp[i] for i in seq[::-1]])

#function that merges all adjacent contigs but that start from a GFA and not a list of segments and output a GFA
#does not use the complex contig/link data structure defined above, as it takes a lot of time
def merge_adjacent_contigs_GFA(gfa_in, gfa_out):

    #go through the GFA and index all contigs, and all links
    segments = [] #just a list of segment names
    depths = {} #associates a segment with its depth
    lengths = {} #associates a segment with its length
    links = {} #associates a segment with two lists of tuples ([(neighbor, end, CIGAR), ...], [(neighbor, end, CIGAR), ...])
    location_of_segments_in_GFA_file = {} #associates a segment with its location in the GFA file (to recover the sequence quickly)

    with open(gfa_in, 'r') as f:

        line = f.readline()
        while line :

            if line[0] == 'S' :
                length_line = len(line)
                line = line.strip().split('\t')
                location_of_segments_in_GFA_file[line[1]] = f.tell()-length_line
                segments.append(line[1])
                if line[1] not in links :
                    links[line[1]] = [[], []]
                lengths[line[1]] = len(line[2])
                #look for the depth
                for field in line[3:] :
                    if field[:2] == "DP" :
                        depths[line[1]] = float(field[5:])
                        break

            elif line[0] == 'L' :
                line = line.split('\t')
                if line[1] not in links :
                    links[line[1]] = [[], []]
                if line[3] not in links :
                    links[line[3]] = [[], []]

                side1 = 1
                side2 = 0
                if line[2] == '-' :
                    side1 = 0
                if line[4] == '-' :
                    side2 = 1

                links[line[1]][side1].append((line[3], side2, line[5].strip()))
                if line[3] != line[1] or side1 != side2 : #to add only one link if it is a self link
                    links[line[3]][side2].append((line[1], side1, line[5].strip()))

            line = f.readline()

    #merge the contigs
    old_segments_to_new_segments = {} #associates the old segment name with (the new segment name, endOfTheNewSegment) #if it is not at an end put -1 we dont care there will be no link
    new_segments = [] #each new segment is a list [(segment, orientation, CIGAR), ...]
    for segment in segments :

        for end in range(2):
            if segment not in old_segments_to_new_segments :
                # print("Merging ", segment, " at end ", end)
                #check if this is the end of a new contig
                if len(links[segment][end]) != 1 or len(links[links[segment][end][0][0]][links[segment][end][0][1]]) != 1 :

                    #then start a new contig and go on until end of new contig
                    new_segment = [(segment, 1-end, "0M")]
                    endNow = 1-end #end from which we are trying to extend
                    segmentNow = segment
                    already_merged_segments = set([segment])
                    while len(links[segmentNow][endNow]) == 1 and len(links[links[segmentNow][endNow][0][0]][links[segmentNow][endNow][0][1]]) == 1 and links[segmentNow][endNow][0][0] not in already_merged_segments:

                        new_segment.append((links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1], links[segmentNow][endNow][0][2]))

                        segmentNow, endNow = links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1]
                        already_merged_segments.add(segmentNow)

                    #now the new segment is created, fill the old_segments_to_new_segments
                    new_segment_name = "_".join([i[0] for i in new_segment])
                    for i, s in enumerate(new_segment):
                        if i== 0 and i == len(new_segment)-1 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 2)
                        elif i == 0 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 0)
                        elif i == len(new_segment)-1 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 1)
                        else :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, -1)

                    new_segments.append(new_segment)

    #some contigs were in loops with only contigs that have one neighbor left and right and still have not been merged, so do this now
    for segment in segments :

        if segment not in old_segments_to_new_segments : #then it is in a loop
            
            end = 0
            #then start a new contig and go on until end of new contig
            new_segment = [(segment, 1-end, "0M")]
            endNow = 1-end #end from which we are trying to extend
            segmentNow = segment
            already_merged_segments = set([segment])
            while len(links[segmentNow][endNow]) == 1 and len(links[links[segmentNow][endNow][0][0]][links[segmentNow][endNow][0][1]]) == 1 and links[segmentNow][endNow][0][0] not in already_merged_segments:

                new_segment.append((links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1], links[segmentNow][endNow][0][2]))

                segmentNow, endNow = links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1]
                already_merged_segments.add(segmentNow)

            #now the new segment is created, fill the old_segments_to_new_segments
            new_segment_name = "_".join([i[0] for i in new_segment])
            for i, s in enumerate(new_segment):
                if i== 0 and i == len(new_segment)-1 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 2)
                elif i == 0 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 0)
                elif i == len(new_segment)-1 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 1)
                else :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, -1)

            new_segments.append(new_segment)
        

    #write the new GFA
    # print("outputting gfa")
    gfa = open(gfa_in, 'r')
    with open(gfa_out, 'w') as f:

        L_lines = []
        number_out = 0

        for new_segment in new_segments :

            # print("Writing new segment ", number_out)
            number_out += 1

            #compute the sequence of the new segment
            seq = ""
            depth_total = 0
            length_total = 0
            # print("Looking for sequence of ", new_segment)
            for seg in new_segment :

                length_total += lengths[seg[0]]
                depth_total += depths[seg[0]]*lengths[seg[0]]

                line_subsegment = ""
                gfa.seek(location_of_segments_in_GFA_file[seg[0]])
                line_subsegment = gfa.readline()
                line_subsegment = line_subsegment.strip().split('\t')
                if len(line_subsegment) < 3 :
                    seq_subsegment = ""
                else :
                    seq_subsegment = line_subsegment[2]
                gfa.seek(0)

                #if the CIGAR contains something else than M, print a warning, sum the M, retract the D, add the I and convert that to M
                length_overlap = 0
                string_of_overlap = ""
                warning = False
                for i in seg[2] :
                    if i == 'M' :
                        length_overlap += int(string_of_overlap)
                        string_of_overlap = ""
                    elif i == 'D' :
                        length_overlap -= int(string_of_overlap)
                        string_of_overlap = ""
                        warning = True
                    elif i == 'I' :
                        length_overlap += int(string_of_overlap)
                        string_of_overlap = ""
                        warning = True
                    elif not i.isdigit() :
                        print("Warning: CIGAR contains something else than M, D or I. Ignoring ", seg[2])
                        warning = True
                        string_of_overlap = ""
                    else :
                        string_of_overlap += i

                if warning :
                    print("Warning: CIGAR contains something else than M, merging contigs is not well defined, turning ", seg[2], " into ", str(length_overlap) + "M" )
                
                if seg[1] == 0 :
                    seq_subsegment = reverse_complement(seq_subsegment)

                seq += seq_subsegment[length_overlap:]

                # print("After having added ", seg[0], " ", seq)

            #write the new segment
            new_segment_name = "_".join([i[0] for i in new_segment])
            f.write("S\t"+new_segment_name+"\t"+seq+"\tDP:f:"+ str(depth_total/length_total) +"\n")

            #now check if the new segment has links at its ends
            #first left end
            for l in links[new_segment[0][0]][1-new_segment[0][1]] :
                old_neighbor = l[0]
                new_neighbor_name = old_segments_to_new_segments[old_neighbor][0]
                new_neighbor_end = old_segments_to_new_segments[old_neighbor][1]
                if new_neighbor_end == 2: #means this is a non-merged segment, fall back to the original orientation
                    new_neighbor_end = l[1]

                if new_neighbor_end == -1 :
                    print("ERROR : new segment has a link at its left end with a segment that is not at an end")
                    sys.exit(1)

                if new_segment_name < new_neighbor_name : #to avoid writing the same link twice
                    L_line = "L\t"+new_segment_name+"\t-\t" + new_neighbor_name + "\t"+ str("+-"[new_neighbor_end]) +"\t" + l[2] + "\n"
                    L_lines.append(L_line)

            #then right end
            for l in links[new_segment[-1][0]][new_segment[-1][1]] :
                old_neighbor = l[0]
                new_neighbor_name = old_segments_to_new_segments[old_neighbor][0]
                new_neighbor_end = old_segments_to_new_segments[old_neighbor][1]
                if new_neighbor_end == 2: #means this is a non-merged segment, fall back to the original orientation
                    new_neighbor_end = l[1]

                if new_neighbor_end == -1 :
                    print("ERROR : new segment has a link at its right end with a segment that is not at an end")
                    sys.exit(1)

                if new_segment_name <= new_neighbor_name : #to avoid writing the same link twice
                    L_lines.append("L\t"+new_segment_name+"\t+\t" + new_neighbor_name + "\t"+ "+-"[new_neighbor_end] +"\t" + l[2] + "\n")

        for l in L_lines :
            f.write(l)

    gfa.close()
    f.close()


def main() :

    args = parse_args()

    #parse the GFA file
    set_of_contigs = {} #set of contigs: name -> (length, depth)
    links = {} #links: contigname -> [[(othercontigname, othercontigend), ...],[(othercontigname, othercontigend), ...]], first list being the links to the left, second list being the links to the right

    f = open(args.input_assembly, 'r')
    for line in f :
        ls = line.split()
        if ls[0] == "S" :
            #try to parse the depth of the contig
            depth = 0
            for i in ls :
                if i.startswith("DP") :
                    depth = float(i[5:])
            if depth == 0 :
                print("Warning: no depth information for contig ", ls[1])
            set_of_contigs[ls[1]] = (len(ls[2]), depth)

        elif ls[0] == "L" :
            end1 = 1
            if ls[2] == "-" :
                end1 = 0
            end2 = 0
            if ls[4] == "-" :
                end2 = 1

            if ls[1] not in links :
                links[ls[1]] = [[],[]]
            if ls[3] not in links :
                links[ls[3]] = [[],[]]

            links[ls[1]][end1].append((ls[3], end2))
            links[ls[3]][end2].append((ls[1], end1))

    f.close()

    #now go through the contigs and pop the bubbles
    contigs_to_delete = []
    for contig in set_of_contigs.keys() :
        if contig in contigs_to_delete :
            continue
        for end in range(2):
            if len(links[contig][end]) == 2 and links[contig][end][0][0] not in contigs_to_delete and links[contig][end][1][0] not in contigs_to_delete:

                #see if the two neighbors 1) have one link left and right and 2) are linked to the same contig on the other side
                neighbor1 = links[contig][end][0][0]
                end1 = links[contig][end][0][1]
                neighbor2 = links[contig][end][1][0]
                end2 = links[contig][end][1][1]

                if len(links[neighbor1][end1]) == 1 and len(links[neighbor2][end2]) == 1 and len(links[neighbor1][1-end1]) == 1 and len(links[neighbor2][1-end2]) == 1 and links[neighbor1][1-end1][0][0] == links[neighbor2][1-end2][0][0] :

                    #this is a bubble, see if we can delete one branch
                    depth1 = set_of_contigs[neighbor1][1]
                    depth2 = set_of_contigs[neighbor2][1]

                    if depth1 > depth2 :
                        if depth2 < args.absolute_threshold or depth2 <= args.relative_threshold * depth1 :
                            contigs_to_delete.append(neighbor2)
                            if args.transfer_coverage :
                                set_of_contigs[neighbor1] = (set_of_contigs[neighbor1][0], set_of_contigs[neighbor1][1] + depth2)
                    else :
                        if depth1 < args.absolute_threshold or depth1 <= args.relative_threshold * depth2 :
                            contigs_to_delete.append(neighbor1)
                            if args.transfer_coverage :
                                set_of_contigs[neighbor2] = (set_of_contigs[neighbor2][0], set_of_contigs[neighbor2][1] + depth1)

    #now write the new assembly
    tmp_out = args.input_assembly + ".tmp"
    f = open(tmp_out, 'w')
    fi = open(args.input_assembly, 'r')
    f.write("H\tVN:Z:1.0\n")
    for line in fi :
        ls = line.split()
        if ls[0] == "S" :
            if ls[1] not in contigs_to_delete :
                #output with the new depth
                for i in ls :
                    if i.startswith("DP") :
                        f.write("DP:f:" + str(set_of_contigs[ls[1]][1]) + "\t")
                    else :
                        f.write(i + "\t")
                f.write("\n")
        elif ls[0] == "L" :
            if ls[1] not in contigs_to_delete and ls[3] not in contigs_to_delete :
                f.write(line)

    f.close()
    fi.close()

    #merge the contigs
    if not args.dont_merge :
        merge_adjacent_contigs_GFA(tmp_out, args.output_assembly)
    else :
        shutil.move(tmp_out, args.output_assembly)
            

if __name__ == "__main__" :
    main()
