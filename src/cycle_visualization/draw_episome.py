
import csv
import re
import sys, getopt

from graphic_elements import Rectangle, HorizontalLine, VerticalLine, Arrow, Text, CycleSection
from .. import hg19util as hg

class Segment():
    def __init__(self, start, end, viral, segment_name, chr_name):
        self.start = start
        self.end = end
        self.viral = viral
        self.segment_name = segment_name
        self.chr_name = chr_name


class EpisomeDrawer():


    def __init__(self, image_output=True, min_copy_number_to_show=0):
        self.BASE_URL = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=%s:%s-%s&hgsid=569969581_amo4eiY8UIG9UO5nhpNMBPEzqa3h'
        self.image_output = image_output
        self.file_names = []
        self.min_copy_number_to_show = min_copy_number_to_show
        self.reconstructed_cycles = []
        self.reconstructed_segments = []
        self.cycle_sections = []
        self.rectangles = []
        self.horizontal_lines = []
        self.vertical_lines = []
        self.arrows = []
        self.texts = []
        self.axes_labels_texts = []
        self.cycle_show = []

    #  readDataFile
    # Segment 1 chr2 195486244 195686719
    # Segment 11 gi|310698439|ref|NC_001526.2| 2 2594
    # Cycle=1;Copy_count=5.59055418897;Segments=8+,12+,11+,3+
    # Cycle=3;Copy_count=1.34109844616;Segments=2+,9-,10+,12+,11+
    def readDataFile(self, content):
        intervals = []
        segs = []
        seg_name_to_index_map = {}
        cycles = []
        directions = []
        cycles_names = []
        copy_counts = []
        tmpl = []
        chr_dict = {}  # Dictionary that maintains the chromosomes, etc.
        lines = content.split('\n')
        for linestr in lines:
            # print (linestr)
            if (linestr == [] or linestr == ''):
                pass
            elif re.match('Interval', linestr):
                _, interval_id, chr_name, start, end = linestr.strip().split()
                start = int(start)
                end = int(end)
                interval_id = int(interval_id)
                intervals.append([chr_name, start, end])
            elif (re.match('Segment', linestr)):
                line = re.split('[\t\s]+', linestr)
                # print(line[3], line[4])
                m = re.search('(gi\|[^\|]+)', line[2])
                if (m != None):
                    line[2] = m.group(1)
                if (chr_dict.get(line[2]) == None):
                    chr_dict.update({line[2]:[0, 0]})
                if (re.match('chr', line[2])):
                    viral = 0  # Human segment
                elif (re.match('gi|', line[2])):
                    viral = 1  # Viral segment
                segs.append((int(line[3]), int(line[4]), viral, line[1], line[2]))
                seg_name_to_index_map[line[1]] = len(segs) - 1

            elif (re.match('Cycle=', linestr)):
                cycle = []
                dirs = []
                c = linestr.split(';')
                c_dict = {cc.split('=')[0]:cc.split('=')[1] for cc in c}
                segment_counter = 0
                for x in c_dict['Segments'].split(','):
                    obj = re.match('(\w+)(.)', x)
                    if (obj):
                        mynum = obj.group(1)
                    if (obj and obj.group(2) == '+' and (mynum != '0' or segment_counter > 0)):  # discard initial 0
                        cycle.append(mynum)
                        dirs.append(1 if mynum != '0' else 0)
                    elif (obj and obj.group(2) == '-'):  # Note that we keep the 0 at the end,
                        cycle.append(mynum)
                        dirs.append(-1 if mynum != '0' else 0)
                    segment_counter += 1
                cycle_name = c_dict['Cycle']
                cycles_names.append(cycle_name)
                copy_count = float(c_dict['Copy_count'])
                if copy_count > self.min_copy_number_to_show:
                    copy_counts.append(copy_count)
                    cycles.append(cycle)
                    directions.append(dirs)
        # print(segs)
        # print(cycles)
        # print(chr_dict)
        return (intervals, segs, seg_name_to_index_map, cycles, directions, chr_dict, copy_counts, cycles_names)


    #==================================================================
    # compute_chr_offsets
    # Gives equal size to each of the chromosomes in the cycle, and to the viruses
    #==================================================================
    def compute_chr_offsets(self, chr_offs):
        numkeys = len(chr_offs)
        if (numkeys > 1):
            v_len = 0.2
            h_len = (1.0 - v_len) / (numkeys - 1)
            offs = 0.0
            # chr_offs['src']=[0,s_len]
            # offs=offs+s_len
            for ch in sorted(chr_offs.keys()):
                # print "Chr types:", ch
                if (re.match('chr', ch)):
                    chr_offs[ch] = [offs, h_len]
                    offs = offs + h_len
                elif (re.match('gi', ch)):  # Viral
                    chr_offs[ch] = [offs, v_len]
                    offs = offs + v_len
                else:  # Other viruses
                    chr_offs[ch] = [offs, v_len]
                    offs = offs + v_len
        elif (numkeys == 1):
            for ch in sorted(chr_offs.keys()):
                # print "Chr types:", ch
                if (re.match('gi', ch)):  # Viral
                    chr_offs[ch] = [0, 1]
                elif (re.match('chr', ch)):
                    chr_offs[ch] = [0, 1]
                else:  # Other viruses
                    chr_offs[ch] = [0, 1]
        return chr_offs

    def compute_intervals(self, segments):
        intervals = {}
        for seg in segments:
            intervals[seg[4]] = [1e20, 0]
        for seg in segments:
            chr_name = seg[4]
            intervals[chr_name][0] = min(intervals[chr_name][0], seg[0])
            intervals[chr_name][1] = max(intervals[chr_name][1], seg[1])
        intervals_list = []
        for chr_name, points in intervals.items():
            intervals_list.append([chr_name, points[0], points[1]])
        return intervals_list


    #==================================================================
    # findMaxIntervals
    # returns size of largest intervals divided by 3 .
    #==================================================================
    def findMaxIntervals(self, segs):
        maxIntvl = {}
        for s in segs:
            if (maxIntvl.get(s[4]) == None):
                maxIntvl.update({s[4]:0})
            maxIntvl[s[4]] = max(s[1] - s[0], maxIntvl[s[4]])
        for ch in maxIntvl.keys():
            maxIntvl[ch] = maxIntvl[ch] / 3.0
        return(maxIntvl)

    # insertPointInSortList(s,L).
    # Sorted list L is a 2-tuple, [coordinate, 0/1]. )
    # denotes begin, 1 denotes end
    def insertPointInSortList(self, s, L):
        i = 0
        done = False
        while(i < len(L) and not(done)):
            # print("i:",i,len(L))
            if(s[0] <= L[i][0]):
                L.insert(i, s)
                done = True
            i = i + 1
        if (not(done)):
            L.append(s)

    def makeListOfSegmentEndPoints(self, segs, typ):
        sortedL = []
        for s in segs:
            if (s[4] == typ):
                self.insertPointInSortList([s[0], 0], sortedL)
                self.insertPointInSortList([s[1], 1], sortedL)
        # print("sortedL:",typ,sortedL)
        return(sortedL)

    # #mergeIntervals
    # Input is a list of sorted end-points, and because we pop below, we go from the largest end-point
    # Each element in sortedL is a tuple [coordinate, 0/1]
    # It is like a balanced parentheses string. Keep going until begin intervals equales end intervals.
    # At the end, we should store the fist end-point (a begin interval), and the last end-point, as a merged interval
    def mergeIntervals(self, sortedL):
        mergeL = []
        while (not(sortedL == [])):
            x = sortedL.pop()  # Largest end-point, therefore x[0] will be 'interval-end'.
            n = 1
            while (n > 0):
                tmpobj = sortedL.pop()
                if(tmpobj[1] > 0):  # End interval
                    n = n + 1
                else:  # Begin interval
                    n = n - 1
            mergeL.insert(0, [tmpobj[0], x[0]])  # Now stored as begin-end.
        # print("mergeL", mergeL)
        return(mergeL)

    def compactIntervals(self, L, maxGap):
        span = 0
        compactL = []
        x = L.pop()
        while (not(L == [])):
            # print("x:", x)
            x1 = L.pop()
            if (x[0] - x1[1] < maxGap):
                x[0] = x1[0]
            else:
                compactL.insert(0, x)
                span = span + x[1] - x[0] + maxGap
                x = [x1[0], x1[1]]
        compactL.insert(0, x)
        span = span + x[1] - x[0] + maxGap
        # print("compactL:",compactL)
        return(compactL, span)


    def findMergedIntervalOffset(self, mergedIntvl, beg, gap):
        offset = 0;
        for (intbeg, intend) in mergedIntvl:
            # print("intbeg,intend:",intbeg,intend,beg)
            if (intend < beg):
                offset += intend - intbeg + gap
            elif (intbeg <= beg):
                offset += beg - intbeg
            # print("offset:", offset)
        return(offset)

    def convertSegmentCoordinates(self, segs, span, compact, maxGap, ilist):
        newlist = []
        for s in segs:
            len = s[1] - s[0]
            # offset = self.findMergedIntervalOffset(compact[s[4]], s[0], maxGap[s[4]])
            # print "chr_offs", s[4]
            # print "chr_offs", chr_offs[s[4]][0]
            # x1 = chr_offs[s[4]][0] + chr_offs[s[4]][1] * (offset + maxGap[s[4]] / 2.0) / (span[s[4]] + 0.0)
            # x2 = chr_offs[s[4]][0] + chr_offs[s[4]][1] * (offset + len + maxGap[s[4]] / 2.0) / (span[s[4]] + 0.0)
            chromosome_name = s[4]
            x1 = ilist.xpos(chromosome_name, s[0])
            x2 = ilist.xpos(chromosome_name, s[1])
            newlist.append([x1, x2, 0, s[0], s[1], s[3]])
        return(newlist)


    def drawAxes(self, span, compact, maxIntvl, chr_offs, cycles_section_bottom):
        for ch in compact.keys():
            for s in compact[ch]:
                # print "ch:", ch, s[0], s[1]
                len = s[1] - s[0]
                offset = self.findMergedIntervalOffset(compact[ch], s[0], maxIntvl[ch])
                x1 = chr_offs[ch][0] + chr_offs[ch][1] * (offset + maxIntvl[ch] / 2.0) / (span[ch] + 0.0)
                x2 = chr_offs[ch][0] + chr_offs[ch][1] * (offset + len + maxIntvl[ch] / 2.0) / (span[ch] + 0.0)
                self.horizontal_lines.append(HorizontalLine(x1, x2, cycles_section_bottom, height=2, css_style='axes_line'))

    def getLabelOffsets(self, x, prevX, prevY):
        ox = 0
        oy = 0
        if (prevX > 0 and x > prevX and x - prevX < 0.02):
            oy = prevY - 5
            ox = ox + 0.011
        else:
            oy = 0
        # print("labeloffset:", x, x - prevX, ox, oy)
        return(ox, oy)


    def drawAxesLabels(self, sorted, compact, maxIntvl, span, chr_offs, cycles_section_bottom):
        prevX = 0
        minGap = 0.02
        newList = []
        for ch in compact.keys():
            for s in compact[ch]:
                offset = self.findMergedIntervalOffset(compact[ch], s[0], maxIntvl[ch])
                offsetX = chr_offs[ch][0] + chr_offs[ch][1] * (offset + maxIntvl[ch] / 2) / (span[ch] + 0.0)
                newList.append([offsetX, s[0]])
                # print ("s[0]:", s[0])

                offset = self.findMergedIntervalOffset(compact[ch], s[1], maxIntvl[ch])
                offsetX = chr_offs[ch][0] + chr_offs[ch][1] * (offset + maxIntvl[ch] / 2) / (span[ch] + 0.0)
                newList.append([offsetX, s[1]])
        for el in newList:
            self.axes_labels_texts.append(Text(el[1], el[0], cycles_section_bottom, css_class='axes_label'))


    def drawAxesDottedLines(self, sorted, compact, maxIntvl, span, chr_offs, cycles_section_top, cycles_section_bottom):
        prevX = 0
        minGap = 0.02
        newList = []
        for ch in sorted.keys():
            for s in sorted[ch]:
                offset = self.findMergedIntervalOffset(compact[ch], s[0], maxIntvl[ch])
                offsetX = chr_offs[ch][0] + chr_offs[ch][1] * (offset + maxIntvl[ch] / 2) / (span[ch] + 0.0)
                newList.append([offsetX, s[0]])

        for el in newList:
            self.vertical_lines.append(VerticalLine(el[0], cycles_section_top, cycles_section_bottom, css_style='vertical_dotted_line'))


    def drawSections(self, intervals, ilist, cycles_section_top, cycles_section_bottom, auto_scale):
        dotted_lines = [(x, type) for x, type, chrom in ilist.offset_breaks() if type == ':']
        for x, line_type in dotted_lines:
            css_style = 'vertical_dotted_line'
            self.vertical_lines.append(VerticalLine(x, cycles_section_top - 5, cycles_section_bottom, width=1, css_style=css_style))

        i = 0
        chromosome_lines = [(0, '--', intervals[0][0])] + [(x, type, chrom) for x, type, chrom in ilist.offset_breaks() if type == '--' or type == '-']
        for x, line_type, name in chromosome_lines:
            if line_type == '-':
                css_style = 'vertical_species_line'
            else:
                css_style = 'vertical_dashed_line'
            end = 1.0 if i == len(chromosome_lines) - 1 else chromosome_lines[i + 1][0]
            self.texts.append(Text(name, (x + end) / 2 - 0.012, cycles_section_bottom + 7))
            bottom = cycles_section_bottom + 20 if line_type == '-' else cycles_section_bottom
            width = 2 if line_type == '-' else 1
            if x != 0:
                self.vertical_lines.append(VerticalLine(x, cycles_section_top - 5, bottom, width=width, css_style=css_style))
            i += 1

        self.horizontal_lines.append(HorizontalLine(0, 1, cycles_section_bottom, height=2, css_style='axes_line'))


    def get_chromosome(self, left, right, chr_offs):
        for i in range(len(chr_offs.keys())):
            keys_list = list(chr_offs.keys())
            ch = keys_list[i]
            if chr_offs[ch][0] <= left and ((i + 1 < len(keys_list) and right <= chr_offs[keys_list[i + 1]][0]) or (i + 1 == len(chr_offs.keys()))):
                return ch


    def drawCycles(self, segs, seg_name_to_index_map, cycles, directions, top_y, chr_offs, struct, copy_counts, cycles_names):
        height = 6
        yoff = top_y - 10
        delx = 0.02
        dely = 2
        # print(cycles)
        cycle_count = 0
        self.cycle_show.append([])
        for i in range(len(cycles)):
            cycle = cycles[i]
            dirs = directions[i]
            top_y_of_cycle = yoff
            bottom_y_of_cycle = yoff
            for elem in cycle:
                bottom_y_of_cycle += height + 2 * dely
            if len(cycle) > 0 and abs(dirs[len(cycle) - 1]) > 0:
                bottom_y_of_cycle += height + 2 * dely
            # yoff = bottom_y_of_cycle
            yoff += 16
            # print('start yoff', yoff)

            rectangles = []
            segments_labels = []
            css_class = '%s_cycle_%s' % (struct, cycle_count)
            cycle_title = 'Read %s' % (cycles_names[cycle_count])# + ': CN=%.2f' % float(copy_counts[cycle_count])
            self.cycle_show[struct].append((cycle_title, css_class))
            prev = -1000
            first = -1000
            segment_num = 0
            # print ("cycle %s :" % str(cycle_count + 1), cycle)
            for j in range(len(cycle)):
                elem = cycle[j]
                # print(cycle)
                # print(elem)
                # print(seg_name_to_index_map)
                if elem != '0':
                    seg = segs[seg_name_to_index_map[elem]]
                    # print ("cycle label:", seg)
                    # print ('cycle:', seg[0], " ", yoff - height / 2, " ", seg[1] - seg[0], ' ', height)
                    chromosome = self.get_chromosome(seg[0], seg[1], chr_offs)
                    tooltip_url = self.BASE_URL % (chromosome, seg[3], seg[4])
                    rectangles.append(Rectangle(seg[0], (yoff), seg[1] - seg[0], height, tooltip_url, 'rectangle ' + css_class))
                    segments_labels.append(Text(seg[-1], seg[0] - 0.04, yoff - 5, css_class=css_class))

                    if (segment_num + 1 < len(cycle) and (cycle[segment_num + 1] != '0' and dirs[j] < 0)) or cycle[len(cycle) - 1] != '0' or (segment_num > 0 and dirs[j] > 0):
                        self.arrows.append(Arrow(seg[0] - 0.01, yoff - height * 1.32, right_direction=dirs[j] > 0, css_class=css_class))
                    if (segment_num + 1 < len(cycle) and (cycle[segment_num + 1] != '0' and dirs[j] > 0)) or cycle[len(cycle) - 1] != '0' or (segment_num > 0 and dirs[j] < 0):
                        self.arrows.append(Arrow(seg[1] + 0.001, yoff - height * 1.32, right_direction=dirs[j] > 0, css_class=css_class))
                if (prev == -1000):
                    first = elem
                    firsty = yoff
                    prev = elem
                elif (elem != '0'):
                    s = segs[seg_name_to_index_map[prev]]
                    t = segs[seg_name_to_index_map[elem]]
                    prev_dir = dirs[j - 1] if j > 0 else dirs[j]
                    x1, y1, typ1, x2, y2, typ2 = self.convertEdgeOrientation(prev_dir, dirs[j], s, t, yoff - height - 2 * dely, yoff)
                    # print("edge coord", x1, y1, typ1, x2, y2, typ2)
                    self.drawEdge(x1, y1, typ1, x2, y2, typ2, height, 0.02 + segment_num * 0.001, 2, css_class)
                prev = elem
                yoff = yoff + height + 2 * dely
                segment_num += 1
            first_segment = segs[seg_name_to_index_map[first]]
            if (elem != '0'):
                last_segment = segs[seg_name_to_index_map[elem]]
                if dirs[len(cycle) - 1] > 0:
                    x1 = last_segment[1]
                    last_dir = 1
                else:
                    x1 = last_segment[0]
                    last_dir = -1
                if dirs[0] > 0:
                    x2 = first_segment[0]
                    first_dir = 1
                else:
                    x2 = first_segment[1]
                    first_dir = -1
                y1 = yoff - height - 2 * dely
                y2 = firsty
                self.drawLastEdge(x1, y1, last_dir, x2, y2, first_dir, height, 0.02, 2, css_class)
                yoff = yoff + height + 2 * dely
            cycle_count += 1
            # print ('bottom and top: ', yoff, top_y_of_cycle)
            # print ('yoff', yoff)
            background_opacity = '1' if cycle_count % 2 else '0';
            self.cycle_sections.append(CycleSection(top_y_of_cycle, bottom_y_of_cycle, cycle_title, rectangles, segments_labels, background_opacity=background_opacity, css_class=css_class))
            yoff = bottom_y_of_cycle + height


    def convertEdgeOrientation(self, el1, el2, s, t, y1, y2):
        # print("edge", el1, el2)
        if (el1 < 0):
            typ1 = 0
            x1 = s[0]
        else:
            typ1 = 1
            x1 = s[1]
        if (el2 < 0):
            typ2 = 1
            x2 = t[1]
        else:
            typ2 = 0
            x2 = t[0]
        if (x1 < x2):
            return(x1, y1, typ1, x2, y2, typ2)
        else:
            return(x2, y2, typ2, x1, y1, typ1)


    def drawLine(self, x1, y1, x2, y2, css_class):
        # print("draw line", x1, y1, x2, y2)
        llim = 0.001
        rlim = 0.999
        if (x1 < llim):
            x1 = llim
        elif (x1 > rlim):
            x1 = rlim
        if (x2 < llim):
            x2 = llim
        elif (x2 > rlim):
            x2 = rlim

        if x1 == x2:
            bot_y = y1 if y1 > y2 else y2
            top_y = y1 + y2 - bot_y
            self.vertical_lines.append(VerticalLine(x1, top_y, bot_y, 2, 'vertical_line ' + css_class))
        else:
            left_x = x1 if x1 < x2 else x2
            right_x = x1 + x2 - left_x
            self.horizontal_lines.append(HorizontalLine(left_x, right_x, y1, 2, 'horizontal_line ' + css_class))


    def drawEdge(self, x1, y1, t1, x2, y2, t2, h, delx, dely, css_class):
        # print("drawing edge:", x1, y1, t1, x2, y2, t2, h, delx, dely)
        if (y1 < y2):
            dely = -dely
            h = -h

        if (t1 == 1 and t2 == 0):
            if (x2 - x1 >= 2 * delx):  # right to left
                self.drawLine(x1, y1, (x1 + x2) / 2, y1, css_class)
                self.drawLine((x1 + x2) / 2, y2, x2, y2, css_class)
                self.drawLine((x1 + x2) / 2, y1, (x1 + x2) / 2, y2, css_class)
            else:
                self.drawLine(x1, y1, x1 + delx, y1, css_class)
                self.drawLine(x1 + delx, y1, x1 + delx, y1 - h / 2 - dely, css_class)
                self.drawLine(x2 - delx, y1 - h / 2 - dely, x1 + delx, y1 - h / 2 - dely, css_class)
                self.drawLine(x2 - delx, y1 - h / 2 - dely, x2 - delx, y2, css_class)
                self.drawLine(x2 - delx, y2, x2, y2, css_class)
        elif (t1 == 1 and t2 == 1):  # right to right
            self.drawLine(x1, y1, x2 + delx, y1, css_class)
            self.drawLine(x2, y2, x2 + delx, y2, css_class)
            self.drawLine(x2 + delx, y1, x2 + delx, y2, css_class)
        elif (t1 == 0 and t2 == 0):  # left to left
            self.drawLine(x1, y1, x1 - delx, y1, css_class)
            self.drawLine(x1 - delx, y2, x2, y2, css_class)
            self.drawLine(x1 - delx, y1, x1 - delx, y2, css_class)
        else:  # left to right
            self.drawLine(x1, y1, x1 - delx, y1, css_class)
            self.drawLine(max(x1 - delx, 0), y1 - h / 2 - dely, min(x2 + delx, 1.0), y1 - h / 2 - dely, css_class)
            self.drawLine(max(x1 - delx, 0), y1, max(x1 - delx, 0), y1 - h / 2 - dely, css_class)
            self.drawLine(min(x2 + delx, 1.0), y1 - h / 2 - dely, min(x2 + delx, 1.0), y2, css_class)
            self.drawLine(x2, y2, min(x2 + delx, 1.0), y2, css_class)


    def drawLastEdge(self, x1, y1, dir1, x2, y2, dir2, h, delx, dely, css_class):
        minx = min(x1, x2) - delx
        if dir1 < 0:
            self.drawLine(x1, y1, minx, y1, css_class)
        else:
            self.drawLine(x1, y1, x1 + delx, y1, css_class)
            self.drawLine(x1 + delx, y1, x1 + delx, y1 + h / 2 + dely, css_class)
            self.drawLine(x1 + delx, y1 + h / 2 + dely, minx, y1 + h / 2 + dely, css_class)
            self.drawLine(minx, y1 + h / 2 + dely, minx, y1, css_class)

        self.drawLine(minx, y1, minx, y2, css_class)

        if dir2 > 0:
            self.drawLine(minx, y2, x2, y2, css_class)
        else:
            self.drawLine(minx, y2, minx, y2 - h / 2 - dely, css_class)
            self.drawLine(minx, y2 - h / 2 - dely, x2 + delx, y2 - h / 2 - dely, css_class)
            self.drawLine(x2 + delx, y2 - h / 2 - dely, x2 + delx, y2, css_class)
            self.drawLine(x2, y2, x2 + delx, y2, css_class)


    def draw_episome(self, input_files, output_file=None, auto_scale=0):
        cycles_section_top = 30
        cycles_section_size = 0
        space_between_decompositions = 90 / (1 + auto_scale)
        bottoms = []
        for i in range(len(input_files)):
            if i != 0:
                cycles_section_size += space_between_decompositions
            input_content = input_files[i][1]
            intervals, segments, seg_name_to_index_map, cycles, directions, chr_offs, copy_counts, cycles_names = self.readDataFile(input_content)

            number_of_element = sum(len(x) for x in cycles)
            cycles_section_size += number_of_element * 10
            for cycle in cycles:
                if cycle[-1] != 0:
                    cycles_section_size += 10
            bottoms.append(cycles_section_top + cycles_section_size)
        tops = [cycles_section_top] + [bottom + space_between_decompositions for bottom in bottoms[:-1]]
        cycles_section_bottom = cycles_section_top + cycles_section_size


        # print ('cycle section_top:', cycles_section_top)
        # print ('cycle section_bottom:', cycles_section_bottom)
        for i in range(len(input_files)):
            input_content = input_files[i][1]
            self.file_names.append(Text('%s: %s' % (str(i + 1), input_files[i][0]), 0.5, tops[i] - 27))
            # print ('bottom:', bottoms[i])
            # print ('top:', tops[i])
            intervals, segments, seg_name_to_index_map, cycles, directions, chr_offs, copy_counts, cycles_names = self.readDataFile(input_content)
            if i == 0:
                self.reconstructed_cycles = [cname for cname in cycles_names]
                self.reconstructed_segments = [segment_count for segment_count in range(len(segments))]

            self.compute_chr_offsets(chr_offs)
            if len(intervals) == 0:
                intervals = self.compute_intervals(segments)
            ilist = hg.interval_list([hg.interval(chr_name, start_point, end_point) for chr_name, start_point, end_point in intervals])
            maxIntvl = self.findMaxIntervals(segments)
            sortedL = {}
            compact = {}
            span = {}
            for ch in maxIntvl.keys():
                sortedL[ch] = self.makeListOfSegmentEndPoints(segments, ch)
                sortedCopy = list(sortedL[ch])
                mergeL = self.mergeIntervals(sortedCopy)
                compact[ch], span[ch] = self.compactIntervals(mergeL, maxIntvl[ch])

            newsegs = self.convertSegmentCoordinates(segments, span, compact, maxIntvl, ilist)

            if i == 0:
                self.drawSections(intervals, ilist, cycles_section_top, cycles_section_bottom, auto_scale)

            # self.drawAxesLabels(sortedL, compact, maxIntvl, span, chr_offs, bottoms[i])
            # self.drawAxes(span, compact, maxIntvl, chr_offs, bottoms[i])
            # self.drawAxesDottedLines(sortedL, compact, maxIntvl, span, chr_offs, tops[i], bottoms[i])
            self.drawCycles(newsegs, seg_name_to_index_map, cycles, directions, tops[i], chr_offs, i, copy_counts, cycles_names)

        # fig.canvas.mpl_connect('motion_notify_event', on_plot_hover)

        # plt.show()

if __name__ == '__main__':
    inputs = []
    output_file = None
    auto_scale = 0
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:i:o:s", ["help", "ifile=", "ofile=", "autoscale"])
    except getopt.GetoptError:
        print ('draw_episome.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('draw_episome.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputs.append(arg)
        elif opt in ("-o", "--ofile"):
            output_file = arg
        elif opt in ("-s", "--autoscale"):
            auto_scale = 1
    if len(inputs) == 0:
        if len(sys.argv) < 2:
            print ('Error: Please specify an input file')
            print ('draw_episome.py -i <inputfile> -o <outputfile>')
            sys.exit()
        inputs = [sys.argv[-1]]

    if output_file == None:
        output_file = inputs[0] + '.png'
    ed = EpisomeDrawer()
    ed.draw_episome(input_files=inputs, output_file=output_file, auto_scale=auto_scale)
