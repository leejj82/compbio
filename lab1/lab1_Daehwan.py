#!/usr/bin/env python

"""
CMSC828H, Fall 2010, Lab 1
Daehwan Kim - infphilo@umiacs.umd.edu
"""

import sys
import re
from datetime import datetime, date, time
from copy import deepcopy

def KnuthMorrisPratt(string, pat):
    failure = [-1]*len(pat)
    lens, lenp = len(string), len(pat)
   
    for j in range(1, lenp):
        i = failure[j-1]
               
        while i >= 0 and pat[j] != pat[i+1]:
            i = failure[i]

        if pat[j] == pat[i+1]:
            failure[j] = i+1

    i, j = 0, 0
    while i < lens and j < lenp:
        if string[i] == pat[j]:
            i += 1
            j += 1
        elif j == 0:
            i += 1
        else:
            j = failure[j-1] + 1

    return (i - j, j)

def two(num):
    if num < 10:
        return '0{0}'.format(num)
    else:
        return '{0}'.format(num)

def three(num):
    if num < 100:
        return '0{0}'.format(two(num))
    else:
        return '{0}'.format(num)

class Assembler:
    num_match = 40

    VACANT = 0
    INPLAY = 1
    ELIMINATED = 2

    fragment_length = 0
    fragment_data = ["" for i in range(300)]

    EDGE_UNUSED = 0
    EDGE_REDUNDANT = 1
    EDGE_UNITIG = 2
    EDGE_REDUCED = 3
    EDGE_USED = 4

    mate_pair_range = [1400, 2600]
    input_file_name = 'lab01.fasta'
    input_file_prefix = 'lab01'

    def print_step_string(self, string):
         print("\n----------------------------------------\n")
         print(string)
         print("\n----------------------------------------\n")

    def Reverse(self, string):
        reversed_string = ""

        for i in reversed(range(len(string))):
            if string[i] == 'a':
                reversed_string += 't'
            elif string[i] == 't':
                reversed_string += 'a'
            elif string[i] == 'g':
                reversed_string += 'c'
            else:
                reversed_string += 'g'

        return reversed_string

    
    def print_edge(self, edge):
        print('{0}  {1}  {2}  {3:6d}'.format(three(edge[0]+1), three(edge[1]+1), edge[2], edge[3]))

    def print_edge_list(self, edge_list):
        for edge in edge_list:
            if edge[0] < edge[1] and edge[4] == self.EDGE_UNUSED:
                self.print_edge(edge)

    def print_unitig(self, unitig, unitig_count = 0, output_stream_list = [], contained_list = [], contain_list = []):
        size = len(unitig)
        if size > 0 and (not contained_list or not contained_list[unitig[0][0]]):
            sum = self.fragment_length
            for fragment in unitig:
                sum += fragment[2]

                if contain_list:
                    size += len(contain_list[fragment[0]])

            output_string = "UNI  {0}  {1:3d}  {2:5d}".format(two(unitig_count), size, sum)
            # print(output_string)

            for output_stream in output_stream_list:
                output_stream.write(output_string + "\n")

            for fragment in unitig:
                output_string = "  {0}  {1}  {2:3d}".format(three(fragment[0] + 1), fragment[1], fragment[2])
                # print(output_string)

                for output_stream in output_stream_list:
                    output_stream.write(output_string + "\n")

                if contain_list and len(contain_list[fragment[0]]):
                    for contain_fragment in contain_list[fragment[0]]:
                        direction = contain_fragment[1]
                        if fragment[1] == 'R':
                            if direction == 'F':
                                direction = 'R'
                            else:
                                direction = 'F'

                        output_string = "  {0}  {1}  {2:3d}".format(three(contain_fragment[0] + 1), direction, 0)
                        # print(output_string)

                        for output_stream in output_stream_list:
                            output_stream.write(output_string + "\n")

    def print_unitig_list(self, unitig_list, output_stream_list = [], contained_list = [], contain_list = []):
        unitig_count = 0;
        for unitig in unitig_list:
            size = len(unitig)
            if size > 0:
                unitig_count = unitig_count + 1
                self.print_unitig(unitig, unitig_count, output_stream_list, contained_list, contain_list)

    def reverse_unitig(self, unitig):
        if len(unitig) == 0:
            return

        unitig = deepcopy(unitig)
        unitig.reverse()
        for i in reversed(range(len(unitig))):
            if unitig[i][1] == 'F':
                unitig[i][1] = 'R'
            else:
                unitig[i][1] = 'F'

            if i > 0:
                unitig[i][2] = unitig[i-1][2]

        unitig[0][2] = 0
        return unitig

    def Setting(self, input_file_name):
        self.input_file_name = input_file_name
        self.input_file_prefix = self.input_file_name.split('.')[0]

    def Read_fastaFile(self):
        input_file = open(self.input_file_name, 'r')

        num, prev_num = 0, 0
        string = ""
        for line in input_file:
            if line[0] == '>':
                num = int(line[1:])
                if string:
                    self.fragment_data[prev_num - 1] = string
                string = ""
                prev_num = num
                
            else:
                string += str(line[:-1])

        if string:
            self.fragment_data[num - 1] = string

        if len(self.fragment_data) > 0:
            self.fragment_length = len(self.fragment_data[0])

        input_file.close()


    def Overlap_Impl(self, a, b, first_found_return = True):
        result = []
        
        (offset, match) = KnuthMorrisPratt(a, b)
        if match >= self.num_match:
            result.append(('F', offset))
            if first_found_return:
                return result

        (offset, match) = KnuthMorrisPratt(b, a)
        if match >= self.num_match:
            if len(result) == 0 or result[-1][1] != -offset:
                result.append(('F', -offset))
                if first_found_return:
                    return result

        b = self.Reverse(b)

        (offset, match) = KnuthMorrisPratt(a, b)
        if match >= self.num_match:
            result.append(('R', offset))
            if first_found_return:
                return result

        (offset, match) = KnuthMorrisPratt(b, a)
        if match >= self.num_match:
            if len(result) == 0 or result[-1][0] != 'R' or result[-1][1] != -offset:
                result.append(('R', -offset))
                if first_found_return:
                    return result

        return result


    def Overlap(self):
        print >> sys.stderr, str(datetime.now()), "Overlap begin"
        overlap_file = open(self.input_file_prefix + '.olaps', 'w')

        for i in range(len(self.fragment_data)-1):
            for j in range(i+1, len(self.fragment_data)):
                result = self.Overlap_Impl(self.fragment_data[i], self.fragment_data[j])

                if result:
                    (dir, offset) = result[0]
                    overlap_file.write(' {0}  {1}  {2}  {3:4d}\n'.format(three(i+1), three(j+1), dir, offset))

        overlap_file.close()
        print >> sys.stderr, str(datetime.now()), "Overlap end"

    def Unitig(self):
        overlap_file = open(self.input_file_prefix + '.olaps', 'r')
        unitig_file = open(self.input_file_prefix + '.unis', 'w')

        edge_range = []
        edge_list = []
        max_index = 0
        contain_list = [[] for i in range(1000)]
        contained_list = [False for i in range(1000)]

        line_re = re.compile(r'(\d+)\s+(\d+)\s+(F|R)\s+(-?\d+)')

        for line in overlap_file:
            line_input = list(line_re.search(line[:-1]).groups())

            for num in [0, 1, 3]:
                line_input[num] = int(line_input[num])
                if num <= 1:
                    line_input[num] = line_input[num] - 1

            max_index = max(max_index, line_input[0], line_input[1])

            line_input.append(self.EDGE_UNUSED)
            edge_list.append(line_input)

            distance = line_input[3]

            if distance == 0:
                contained_list[line_input[1]] = True
                contain_list[line_input[0]].append([line_input[1], line_input[2]])

            if line_input[2] == 'F':
                distance = -distance

            edge_list.append([line_input[1], line_input[0], line_input[2], distance, self.EDGE_UNUSED])

        max_index = max_index + 1

        temp_edge_list = []
        for edge in edge_list:
            if contained_list[edge[0]] or contained_list[edge[1]]:
                  pass
            else:
                temp_edge_list.append(edge)

        edge_list = temp_edge_list
        edge_list.sort()

        count = 0
        begin, end = 0, 0
        for i in range(max_index):
            while count < len(edge_list) and edge_list[count][0] == i:
                end = end + 1
                count = count + 1
            edge_range.append([begin, end])
            edge_list[begin:end] = sorted(edge_list[begin:end], key=lambda x:abs(x[3]))
            begin = end

        mark = [self.VACANT for i in range(max_index)]

        for v in range(max_index):
            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                w = edge[1]
                mark[w] = self.INPLAY

            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                w = edge[1]

                begin = edge[3]
                end = begin + self.fragment_length - 1
                range_vw = [min(0, begin), max(self.fragment_length-1, end)]

                if mark[w] == self.INPLAY:
                    for edge2 in edge_list[edge_range[w][0]:edge_range[w][1]]:
                        x = edge2[1]

                        if mark[x] == self.INPLAY and v < x:
                            if edge[2] == 'F':
                                begin2 = edge2[3] + begin
                            else:
                                begin2 = begin - edge2[3]
                            end2 = begin2 + self.fragment_length - 1
                            range_vx = [min(0, begin2), max(self.fragment_length-1, end2)]

                            if range_vx[0] <= range_vw[0] and range_vx[1] >= range_vw[1]:
                                mark[x] = self.ELIMINATED

            for edge in edge_list[edge_range[v][0]:edge_range[v][1]]:
                w = edge[1]
                if mark[w] == self.ELIMINATED:
                    edge[4] = self.EDGE_REDUNDANT

                    for opposite_edge in edge_list[edge_range[w][0]:edge_range[w][1]]:
                        if opposite_edge[1] == v:
                            opposite_edge[4] = self.EDGE_REDUNDANT
                            break

                mark[w] = self.VACANT

        edge_list.sort()
        #self.print_edge_list(edge_list)

        # Unitig
        unitig_list = [[[i, 'F', 0]] for i in range(max_index)]
        vertex_included = [i for i in range(max_index)]

        def belong(v):
            v_temp = v
            while True:
                v_parent = vertex_included[v_temp]
                if v_parent == v_temp:
                    vertex_included[v] = v_parent
                    return v_parent

                v_temp = v_parent

        # both edge and unitig update = 0
        # unitig update = 1
        # no update = 2
        def merge_unitig(edge, option = 0):
            v = belong(edge[0])
            w = belong(edge[1])
            print_flag = False

            if print_flag:
                print_edge(edge)

            unitig_v = deepcopy(unitig_list[v])
            unitig_w = deepcopy(unitig_list[w])

            if unitig_v[0][0] != edge[0]:
                v_index = len(unitig_v) - 1
            else:
                v_index = 0

            if unitig_w[0][0] != edge[1]:
                w_index = len(unitig_w) - 1
            else:
                w_index = 0

            v_info, w_info = unitig_v[v_index], unitig_w[w_index]

            if v_info[0] != edge[0] or w_info[0] != edge[1]:
                print("Too unexpected!")
                return []

            if v_info[1] != 'F':
                unitig_v = self.reverse_unitig(unitig_v)

            if w_info[1] != edge[2]:
                unitig_w = self.reverse_unitig(unitig_w)

            if edge[3] > 0:
                if unitig_v[-1][0] != edge[0] or unitig_w[0][0] != edge[1] or unitig_v[-1][1] != 'F' or unitig_w[0][1] != edge[2]:
                    print("A bit unexpected!")
                    return []
                
                unitig_w[0][2] = edge[3]
                unitig_v.extend(unitig_w)
            else:
                if unitig_v[0][0] != edge[0] or unitig_w[-1][0] != edge[1] or unitig_v[0][1] != 'F' or unitig_w[-1][1] != edge[2]:
                    print("A bit unexpected!")
                    return []

                unitig_v[0][2] = -edge[3]
                unitig_w.extend(unitig_v)
                unitig_v = unitig_w

            if option == 0:
                edge[4] = self.EDGE_UNITIG

            if option == 0:
                for edge2 in edge_list[edge_range[edge[1]][0]:edge_range[edge[1]][1]]:
                    if edge2[1] == edge[0]:
                        edge2[4] = self.EDGE_UNITIG
                        break

            if option == 0 or option == 1:
                vertex_included[w] = v

                unitig_list[w] = []
                unitig_list[v] = unitig_v

            return unitig_v

        def find_edges(v):
            result_edge_list = []
            v = belong(v)

            v_list = [unitig_list[v][0][0]]
            if len(unitig_list[v]) > 1:
                v_list.append(unitig_list[v][-1][0])

            def get_edge_list(w):
                result = []
                for edge in edge_list[edge_range[w][0]:edge_range[w][1]]:
                    if edge[4] == self.EDGE_UNUSED:
                        result.append(edge)
                return result

            result_list = []
            for v_temp in v_list:
                result_list.append(get_edge_list(v_temp))
                
            return result_list

        def merge_unitigs():
            for v in range(max_index):
                v = belong(v)
                v_edge_list_list = find_edges(v)

                for v_edge_list in v_edge_list_list:
                    if len(v_edge_list) > 2:
                        continue

                    if len(v_edge_list) == 2 and belong(v_edge_list[0][1]) == belong(v_edge_list[1][1]):
                        continue
                
                    for e in v_edge_list:
                        temp_edge_list = []
                        for e2 in edge_list[edge_range[e[1]][0]:edge_range[e[1]][1]]:
                            if e2[4] == self.EDGE_UNUSED:
                                temp_edge_list.append(e2)

                        left_count, right_count, left = 0, 0, 1
                        for e2 in temp_edge_list:
                            if e2[3] < 0:
                                left_count += 1
                            elif e2[3] > 0:
                                right_count += 1
                                if e2[1] == e[0]:
                                    left = 0

                        if (left == 1 and left_count > 1) or (left == 0 and right_count > 1):
                            continue

                        merge_unitig(e)

        merge_unitigs()

        self.print_unitig_list(unitig_list, [unitig_file], contained_list, contain_list)
        # self.print_edge_list(edge_list)

        overlap_file.close()
        unitig_file.close()

    def ContigAndExtra(self):
        unitig_file = open(self.input_file_prefix + '.unis', 'r')
        contig_file = open(self.input_file_prefix + '.contig', 'w')
        extra_file = open(self.input_file_prefix + '.extra', 'w')

        unitig_list = []

        line_head_re = re.compile(r'^UNI\s+(\d+)\s+(\d+)\s+(\d+)$')
        line_content_re = re.compile(r'^\s+(\d+)\s+(F|R)\s+(\d+)$')

        input_unitig = []
        for line in unitig_file:
            if line[0] == 'U':
                line_input = list(line_head_re.search(line[:-1]).groups())
                if input_unitig:
                    unitig_list.append(input_unitig)
                    
                input_unitig = []
                
            else:
                line_input = list(line_content_re.search(line[:-1]).groups())
                
                for num in [0, 2]:
                    line_input[num] = int(line_input[num])
                    if num == 0:
                        line_input[num] = line_input[num] - 1

                input_unitig.append(line_input)

        if input_unitig:
            unitig_list.append(input_unitig)

        def make_raw_sequence(unitig):
            result = ""
            for i in range(len(unitig)):
                fragment = unitig[i]

                output_string = deepcopy(self.fragment_data[fragment[0]])
                if fragment[1] == 'R':
                    output_string = self.Reverse(output_string)

                if i == 0 or fragment[2] >= self.fragment_length:
                    result += output_string
                else:
                    result += output_string[self.fragment_length-fragment[2]:]

            return result

        edge_list = []
        edge_range = []

        # reduce the number of match among unitigs
        self.num_match = 17

        for i in range(len(unitig_list)):
            unitig1 = unitig_list[i]
            unitig1_raw_sequence = make_raw_sequence(unitig1)

            for j in range(i+1, len(unitig_list)):
                unitig2 = unitig_list[j]
                unitig2_raw_sequence = make_raw_sequence(unitig2)
                
                result = self.Overlap_Impl(unitig1_raw_sequence, unitig2_raw_sequence, False)

                for temp in result:
                    (dir, offset) = temp
   
                    edge = [i, j, dir, offset, self.EDGE_UNUSED]
                    edge_list.append(edge)

                    distance = edge[3]
                    if edge[2] == 'F':
                        distance = -distance

                    edge_list.append([edge[1], edge[0], edge[2], distance, self.EDGE_UNUSED])

        edge_list.sort()
           
        count = 0
        begin, end = 0, 0
        for i in range(len(unitig_list)):
            while count < len(edge_list) and edge_list[count][0] == i:
                end = end + 1
                count = count + 1
            edge_range.append([begin, end])
            begin = end

        def print_closeness():
            closeness = {}
            for ui in range(len(unitig_list)):
                for fragment in unitig_list[ui]:
                    num1 = fragment[0]
                    if num1 % 2 == 0:
                        num2 = num1 + 1
                    else:
                        num2 = num1 - 1

                    for uj in range(ui+1, len(unitig_list)):
                        for fragment2 in unitig_list[uj]:
                            if num2 == fragment2[0]:
                                index = min(ui, uj)*1000 + max(ui, uj)
                                if index not in closeness:
                                    closeness[index] = 1
                                else:
                                    closeness[index] = closeness[index] + 1

            for k, v in closeness.items():
                print("{0:3d}  {1:3d}\t{2:3d}".format(int(k/1000) + 1, k%1000 + 1, v))

        def num_mates(u, v):
            pos, neg, consumed = 0, 0, 0
            skip = [False for _ in range(len(self.fragment_data))]
            for f in range(len(u)):
                fragment = u[f]
                num1 = fragment[0]
                if num1 % 2 == 0:
                    num2 = num1 + 1
                else:
                    num2 = num1 - 1

                if skip[num1]:
                    assert skip[num2]
                    continue

                found = False
                for f2 in range(f + 1, len(u)):
                    fragment2 = u[f2]
                    if num2 == fragment2[0] and fragment[1] == 'F' and fragment2[1] == 'R':
                        found = True
                        consumed += 1
                        skip[num1] = True
                        skip[num2] = True
                        break
                if found:
                    continue

            for f in range(len(u)):
                fragment = u[f]
                num1 = fragment[0]
                if num1 % 2 == 0:
                    num2 = num1 + 1
                else:
                    num2 = num1 - 1

                if skip[num1]:
                    assert skip[num2]
                    continue
                
                for fragment2 in v:
                    if num2 == fragment2[0]:
                        if fragment[1] == 'F' and fragment2[1] == 'R':
                            pos += 1
                        else:
                            neg += 1

            return pos, neg, consumed

        # self.print_unitig_list(unitig_list)
        # print("")
        # self.print_edge_list(edge_list)
        # print("")
        # print_closeness()

        def merge_unitig(v, w, wd):
            if type(v).__name__ == 'int':
                v = unitig_list[v]

            if type(w).__name__ == 'int':
                w = deepcopy(unitig_list[w])

            if wd == 'R':
                w = self.reverse_unitig(w)

            reference = self.fragment_data[v[-1][0]]
            if v[-1][1] == 'R':
                reference = self.Reverse(reference)
                
            result = self.Overlap_Impl(reference, self.fragment_data[w[0][0]], False)
            assert(len(result) == 1)
            dir, offset = 0, 0
            (dir, offset) = result[0]

            w[0][1] = dir
            w[0][2] = offset
            v.extend(w)
            return v

        def merge_fragment(v, w):
            return merge_unitig(v, [[w, 'F', 0]], 'F')

        def build_contig(contig, unitig_list, depth = 0):
            _, _, consumed = num_mates(contig, [])
            if len(self.fragment_data) == consumed * 2:
                return contig
            
            reference = self.fragment_data[contig[-1][0]]
            if contig[-1][1] == 'R':
                reference = self.Reverse(reference)

            candidates = []
            for u in range(len(unitig_list)):
                for wd in ['F', 'R']:
                    w = deepcopy(unitig_list[u])
                    if wd == 'R':
                        w = self.reverse_unitig(w)
                    result = self.Overlap_Impl(reference, self.fragment_data[w[0][0]], False)

                    if len(result) <= 0:
                        continue
                    assert(len(result) == 1)
                    dir, offset = 0, 0
                    (dir, offset) = result[0]
                    if dir != w[0][1] or offset == 0:
                        continue

                    pos, neg, consumed = num_mates(contig, w)
                    if pos <= neg:
                        continue
                        
                    w[0][1] = dir
                    w[0][2] = offset
                    v = deepcopy(contig)
                    v.extend(w)

                    candidates.append([v, pos, neg, consumed, u, wd])

            def my_cmp(a, b):
                if a[1] != b[1]:
                    return b[1] - a[1]
                return a[2] - b[2]

            candidates = sorted(candidates, cmp=my_cmp)                
            for candidate, pos, neg, consumed, u, wd in candidates:
                print "Depth %d:" % depth, "contig vs.", u, wd, "closeness:", pos, neg, consumed, "after:", num_mates(candidate, [])[2]
                build_contig(candidate, unitig_list, depth + 1)

        build_contig(unitig_list[2], unitig_list)

        unitig = unitig_list[2]
        unitig = merge_unitig(unitig, 1, 'F')
        unitig = merge_unitig(unitig, 5, 'R')
        unitig = merge_unitig(unitig, 1, 'R')
        unitig = merge_unitig(unitig, 3, 'R')
        unitig = merge_unitig(unitig, 1, 'F')
        unitig = merge_unitig(unitig, 4, 'R')
        unitig = merge_unitig(unitig, 1, 'F')
        unitig = merge_unitig(unitig, 0, 'F')

        contig_file.write(">Contig\n")
        contig = make_raw_sequence(unitig)

        line_len = 60
        lines = [contig[i*line_len:(i+1)*line_len] for i in range(0, (len(contig) + line_len + 1) / line_len)]
        for line in lines:
            contig_file.write(line + "\n")

        repeat_list = [unitig_list[1][i][0] for i in range(len(unitig_list[1]))]
        fragment_list = [[] for i in range(len(self.fragment_data))]
        offset = 0
        for fragment in unitig:
            index = fragment[0]
            direction = fragment[1]
            offset += fragment[2]
            
            if index in repeat_list:
                continue
            
            if direction == 'F':
                fragment_list[index] = [index, offset, offset + self.fragment_length - 1]            
            else:
                fragment_list[index] = [index, offset + self.fragment_length - 1, offset]

        offset = 0
        for fragment in unitig:
            index = fragment[0]
            direction = fragment[1]
            offset += fragment[2]
            
            if index not in repeat_list:
                continue

            if index % 2 == 0:
                mate = index + 1
            else:
                mate = index - 1

            assert(fragment_list[mate])
            mate_fragment = fragment_list[mate]

            distance = 0
            if mate_fragment[1] < mate_fragment[2]:
                distance = offset - mate_fragment[2]
            else:
                distance = mate_fragment[2] - offset - self.fragment_length + 1

            if distance < self.mate_pair_range[0] or distance > self.mate_pair_range[1]:
                continue

            if direction == 'F':
                fragment_list[index] = [index, offset, offset + self.fragment_length - 1]            
            else:
                fragment_list[index] = [index, offset + self.fragment_length - 1, offset]


        for fragment in fragment_list:
            output_string = "{0}  {1:5d}  {2:5d}".format(three(fragment[0]+1), fragment[1]+1, fragment[2]+1)
            extra_file.write(output_string + "\n")

        unitig_file.close()
        contig_file.close()
        extra_file.close()


    def execute(self, input_file_name):
        self.Setting(input_file_name)
        self.Read_fastaFile()
        # self.Overlap()
        self.Unitig()
        self.ContigAndExtra()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        my_assembler = Assembler()
        my_assembler.execute(sys.argv[-1])
    else:
        print >> sys.stderr, "lab1.py input_fasta_file"
