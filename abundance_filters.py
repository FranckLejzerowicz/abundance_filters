#!/usr/bin/env python

import cProfile
import subprocess
from Bio import SeqIO
import time
import re
import glob
import sys,os
import argparse
from argparse import RawTextHelpFormatter
from numpy import array

__author__ = "Franck Lejzerowicz"
__copyright__ = "Copyright 2017, The Deep-Sea Microbiome Project"
__credits__ = ["Jan Pawlowski"]
__license__ = "GPL V3"
__version__ = "1.0"
__maintainer__ = "Franck Lejzerowicz"
__email__ = "franck.lejzerowicz@unige.ch"

def abundance_filters():
    """
    An abundance filter based on OTU/ISU tables and multiple possible flavours.
    Automatic (or option-based) recognition of the OTU table format (presence of header / number of left side OTU-data columns).
    Possibility to the apply the filter only to the selected samples (see methods 'choice' and 'replicates') and to return the filtered OTU table for these samples only.
    Filtering on a per-sample basis, or on a across-dataset basis (see option --sum).
    """
    parser=argparse.ArgumentParser(formatter_class = RawTextHelpFormatter)
    parser.add_argument('-i', nargs = '?', required = True, help = 'Input OTU/ISU table name (required)')
    parser.add_argument('-o', nargs = '?', default = argparse.SUPPRESS, help = "Output OTU/ISU table name (default = add arguments extensions to input name, and '.stats' for stats file)")
    parser.add_argument('-t', nargs = '?', type = int, default = 10, help = 'Abundance threshold, inclusive (default = 10)')
    parser.add_argument('-x', nargs = '*', default = [], help = 'Samples to apply filter on and to output, using:\n    * indices numbers (from 1 to number of samples)\n    * regex for sample(s) names)')
    parser.add_argument('-meth', nargs = '?', default = 'simple', choices = ['simple', 'minimum', 'presence', 'choice', 'replicates'], help = "Filtering method: [default = 'simple']\n    * 'simple': filter OTUs/ISUs not abundant enough\n    * 'minimum': filter OTU/ISU not abundant enough in a percent of samples. Use with '-mode' [default = 10]\n    * 'presence': filter OTUs/ISUs not present in a min percent of samples. Use with '-mode' [default = 10]\n    * 'choice': define groups of samples where OTUs/ISUs must be abundant enough. Use with -mode\n    * 'replicates': use groups of samples among which OTUs/ISUs must be abundant enough. Use with -mode)\n    [EACH METHOD CAN OPERATE PER SAMPLE OR ACROSS THE ENTIRE DATASET, USING '--sum']")
    parser.add_argument('-mode', nargs = '*', default = [], help = "Filtering mode, for '-meth':\n    * 'minimum': percent of samples (integer between 1 and 100) [default = 10]\n\t(remove OTU/ISU that do not reach threshold in the percent of samples)\n\te.g. '-meth minimum -mode 50' = OTU/ISU must be abundant enough in at least half of the samples\n    * 'presence': percent of samples (integer between 1 and 100) [default = 10]\n\t(remove OTU/ISU that is not present in the percent of samples)\n\te.g. '-meth presence -mode 50' = OTU/ISU must be present in at least half of the samples\n   (** Both '-meth choice' and '-meth replicates' could be combined with a minimum percent of these samples [default = 100])\n    * 'choice':\n\t- one (or more) space-separated regex to select by sample name\n\t- at least 2 hyphen-separated integers to select by sample indices (starts at 1)\n\t- both regex and indices\n\t[See '--select' for output]\n\te.g. '-meth choice -mode gut_.*_100ppm 30' = OTU/ISU must be abundant enough in 30 %s of the regex-matching samples\n    * 'replicates': file with tab-/comma-separated replicates samples names in rows (first column could be the new 'sample' name)\n\tOptional integers for:\n\t- number of replicate occurrences per group of sample replicates [default = 2]\n\t- number of occurrences in terms of groups of samples replicates (or percent if float between 0 and 1) [default = 1]\n\te.g. '-meth replicates -mode <repFile.txt> 3 10' = OTU/ISU must be abundant enough in >=3 replicates of >= 10 groups of replicates\n   ['minimum', 'choice' and 'replicates': see option '--only' to apply OTU filtering in selected samples only]" % '%'.replace(r"%", r"%%"))
    parser.add_argument('--sum', action = 'store_true', default = False, help = 'Compare to the threshold the SUM of the abundances in all the samples to be filtered.\n(default = not activated = filter per sample)')
    parser.add_argument('--only', action = 'store_true', default = False, help = "Remove ISUs/OTUs only in the samples where they do not reach the abundance threshold.\nAutomatic with '-meth sample' / Manual for methods 'minimum', 'choice' and 'replicates'.\n(default = remove OTU in all samples, or in all replicates of a group of replicates for method 'replicates').")
    parser.add_argument('--selection', action = 'store_true', default = False, help = "Return filtered data table containing only the samples selected using '-meth choice'.")
    parser.add_argument('-name-otus', nargs=1, default = [None], help = "Character string for pass-filter OTUs/ISUs prefix, with OTUs/ISUs indices as suffix\n(only if input table contains no name, else default = no name if option not used).")
    parser.add_argument('-meta', nargs = '?', type = int, default = argparse.SUPPRESS, help = 'Number of non-numeric column in the OTU/ISU table (default = Auto-detect - Warning if OTU/ISU names are numeric)')
    parser.add_argument('--head', action = 'store_true', default = False, help = 'First line is a header line (default = automatic detection).')
    parser.add_argument('-sep', nargs = '?', default = '\t', help = 'Fields separator for the input table (default = tabulation).')
    parser.add_argument('--merge', action = 'store_true', default = False, help = "For '-meth' replicates only: Merge the replicates into one sample (first column entries of the replicates file) and take the average among the replicates (default = don't).")
    parse=parser.parse_args()
    args=vars(parse)

    # parse arguments
    filin = args['i']
    thresh = args['t']
    meth = args['meth']
    mode = args['mode']
    across = args['sum']
    select = args['selection']
    only = args['only']
    sep = args['sep']
    head = args['head']
    merge = args['merge']
    if args.has_key('meta'):
        meta = args['meta']
    else:
        meta = None
    name_otus = args['name_otus'][0]

    # prepare output file names
    if args.has_key('o'):
        filout = args['o']
    else:
        filout = filin.rstrip('txt') + 't%s_%s' % (thresh, meth)
        if meth == 'replicates':
            if len(mode) == 3:
                filout += '%sreps_%sgrps' % (mode[1], mode[2])
            if len(mode) == 2:
                filout += '2reps_%sgrps' % (mode[1])
            else:
                filout += '2reps_1grps'
        if args['x']:
            filout += '_regex-%s' % '-'.join(args['x'])
        if across:
            filout += '_sum'
        if meth in ['choice', 'replicates']:
            if only:
                filout += '_only'
            if select:
                filout += '_select'
        filout += '.tsv'
    # - stats output
    statout = filout[:filout.rindex('.')] + '_stats.tsv'
    Rout = filout[:filout.rindex('.')] + '.R'
    # - table output

    # dict with indices and regex info for subsample selection
    regex = get_samples_indices(args['x'])

    # Warn if filtering method inconsistencies
    if meth == 'presence' and across:
        print "Warning: options '-meth presence' and '--sum' not compatible with (Impossible to sum presences to an abundance threshold value)"

    # main
    tableCode = check_table(filin, sep, head, meta)
    all_samples, samples = get_samples(filin, tableCode, sep, regex)
    meth_mode, samples, maxReps, max_nGrps, subsamples, sample_repname, replicates = get_method(meth, mode, across, only, tableCode, filin, samples, all_samples)
    oTable, lineIdx, filoutPath = write_header_table(filout, samples, all_samples, name_otus, select, merge, tableCode, subsamples, sample_repname, replicates)
    oStats, statoutPath = write_header_stats(args, statout, filout)
    stats, repData = filter_table(filin, tableCode, sep, regex, meth_mode, thresh, samples, all_samples, oTable, name_otus, select, merge, lineIdx, maxReps, max_nGrps)
    write_stats(stats, oStats, samples, all_samples, select)
    oTable.close()
    oStats.close()
    print 'Outputs:'
    print filoutPath
    print statoutPath
    if meth_mode[0] == filt_replicates:
        pdf = replicatesFigure(Rout, repData, samples)
        print os.path.abspath(Rout)
        print os.path.abspath(pdf)
    return 0


def get_repData(repData, reads_in_replicates, maxReps, nGrps):
    """
    Collect numbers of OTUs that occur in each possible number of groups of sample replicates, and for each class of number of replicates per group
    Useful for the make of the figure: see function replicatesFigure(). Hence, only active of method 'replicates' is used.
    """
    # from 2 replicates to the max number of replicates in a group
    for rNum in range(2,maxReps+1):
        n, n2 = get_number_presences(reads_in_replicates, rNum)
        # init exact and cumulative nested dicts with current number of replicates per group as keys
        if repData['Exact'].has_key(rNum) == False:
            repData['Exact'][rNum] = {}
            repData['Cumulative'][rNum] = {}
        # for each number of groups if replicates
        for gNum in range(1, int(nGrps+1)):
            # init another nested dict level to 0
            if repData['Exact'][rNum].has_key(gNum) == False:
                repData['Exact'][rNum][gNum] = 0
                repData['Cumulative'][rNum][gNum] = 0
            # increment the dicts for OTUs present in the required number of groups and replicates per group
            if n >= gNum:
                repData['Cumulative'][rNum][gNum]+=1
            if n2 == gNum:
                repData['Exact'][rNum][gNum]+=1
    return repData


def get_stats(stats, filtered, samples, all_samples, select):
    stats['otu_n'] += 1
    filt_0 = sum(filtered[0])
    if select:
        # kept
        filt_1 = sum([filtered[1][sample[0]] for sample in samples])
        samFilt = [filtered[1][sample[0]] for sample in all_samples if sample in samples]
        if sum(samFilt):
            filt_0 += [filtered[1][sample[0]] for sample in all_samples if sample not in samples]
    else:
        # kept
        filt_1 = sum(filtered[1])
    # fill stats dictoionnary
    if filt_1:
        stats['otu_1'].append(filt_1)
    if filt_0:
        stats['otu_0'].append(filt_0)

    for sample in sorted(all_samples, key=lambda x: int(x[0])):
        sample_idx = sample[0]
        sample_name = sample[1]
        if select:
            if sample in samples:
                curKept = filtered[1][sample_idx]
                curRemoved = filtered[0][sample_idx]
            else:
                curRemoved = filtered[0][sample_idx] + filtered[1][sample_idx]
        else:
            curKept = filtered[1][sample_idx]
            curRemoved = filtered[0][sample_idx]
        stats['s_reads_1'][sample_name] += curKept
        stats['s_reads_0'][sample_name] += curRemoved
        if curKept:
            stats['s_otu_1'][sample_name] += 1
        if curRemoved:
            stats['s_otu_0'][sample_name] += 1
    return stats


def write_stats(stats, o, samples, all_samples, select):
    otu_kept = stats['otu_1']
    otu_removed = stats['otu_0']
    otu_n = stats['otu_n']
    percent_otu = 100*(len(otu_kept)/float(otu_n))
    o.write('Entries kept by filter:\t%s\t%4.2f\t%%\n' % (len(otu_kept), percent_otu))
    o.write('Entries removed by filter:\t%s\t%4.2f\t%%\n' % (int(otu_n-len(otu_kept)), (100-percent_otu)))
    percent_reads = 100*(sum(otu_kept)/float(sum(otu_kept)+sum(otu_removed)))
    o.write('Reads kept by filter:\t%s\t%4.2f\t%%\n' % (sum(otu_kept), percent_reads))
    o.write('Reads removed by filter:\t%s\t%4.2f\t%%\n' % (sum(otu_removed), (100-percent_reads)))
    o.write('Sample index\tSample name\tOTUs kept\tOTUs removed\tReads kept\tReads removed\n')

    s_reads_0 = stats['s_reads_0']
    s_reads_1 = stats['s_reads_1']
    s_otu_0 = stats['s_otu_0']
    s_otu_1 = stats['s_otu_1']
    for sample in sorted(all_samples, key=lambda x: int(x[0])):
        o.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (sample[0], sample[1], s_otu_1[sample[1]], s_otu_0[sample[1]], s_reads_1[sample[1]], s_reads_0[sample[1]]))
    o.close()


def filter_table(filin, code, sep, regex, meth_mode, thresh, samples, all_samples, ot, name_otus, select, merge, lineIdx, maxReps, max_nGrps):
    """
    Main filtering function
    Parses the input table to filter and return a list with for each line (e.g. OTU),
    a tuple (line header, dict with 0/1 key for removed/kept per sample, OTU counter)
    Return a dict on replicates filtering: stays empty if replicates-filtering not used
    """
    # filtering method and method-specific arguments (see function get_method())
    curMethod = meth_mode[0]
    arguments = meth_mode[1:]
    # table format encoding (see function check_table())
    head_code = code[0]
    meta_code = code[1]
    OTU = 0
    repData = {'Exact': {}, 'Cumulative': {}}
    stats = {
        'otu_n': 0,
        'otu_1': [],
        'otu_0': [],
        's_otu_1': {},
        's_otu_0': {},
        's_reads_1': {},
        's_reads_0': {}}
    for sample in all_samples:
        stats['s_otu_1'][sample[1]] = 0
        stats['s_otu_0'][sample[1]] = 0
        stats['s_reads_1'][sample[1]] = 0
        stats['s_reads_0'][sample[1]] = 0
    # parse input OTU/ISU table
    for lindx, line in enumerate(open(filin, 'rU')):
        if lindx and lindx % 1000==0:
            print '%s lines processed' % lindx
        split_line = line.strip().split(sep)
        if len(''.join(split_line)): # skip empty line
            # line field that may correspond to count data or sample names
            reads = split_line[meta_code:]
            # line header (may contain OTU name / taxonomy / ...)
            meta = [x.strip() for x in split_line[:meta_code]]
            # if reads is count data
            if lindx or lindx == 0 and head_code == 0:
                # apply the selected filtering method with passed arguments as parameters
                filtered_line, reads_in_replicates = curMethod(reads, samples, thresh, arguments)
                OTU += 1
                write_line(filtered_line, meta, OTU, ot, samples, all_samples, name_otus, merge, lineIdx)
                stats = get_stats(stats, filtered_line, samples, all_samples, select)
                if curMethod == filt_replicates:
                    repData = get_repData(repData, reads_in_replicates, maxReps, max_nGrps)
    return stats, repData


def get_number_presences(presences_per_rep, rNum):
    n = 0
    n2 = 0
    for x in presences_per_rep:
        if x >= rNum:
            n+=1
            if x == rNum:
                n2+=1
    return n,n2


def replicatesFigure(Rout, repData, samples):
    """
    Make figure with R subprocess
    Number of OTUs per number of groups of sample replicates and per class of number of replcates in these groups
    """
    R = open(Rout, 'w')
    pdf = '%s.pdf' % '.'.join(Rout.split('.')[:-1])
    R.write('pdf("%s", height=5, width=14)\n' % pdf)
    R.write('par(mar = c(5,7,7,2), xpd = NA)\n')
    R.write('layout(matrix(c(1,2, 3), 1, 3), widths=c(1,1,0.33))\n')
    maxRep = len(repData['Exact'].keys())
    colors = ['rainbow(%s)[%s]' % (maxRep, x) for x in range(1, int(maxRep+1))]
    for typ in ['Cumulative', 'Exact']:
        curMax = max([max(repData[typ][x].values()) for x in repData[typ]])
        for rdx, rNum in enumerate(sorted(repData[typ])):
            if rdx:
                R.write('points(x=c(%s), y=c(%s), ylim = c(0,%s), type="l", lwd=2.5, col=%s)\n' % (','.join(map(str, sorted(repData[typ][rNum].keys()))), ','.join(map(str, [repData[typ][rNum][x] for x in sorted(repData[typ][rNum].keys())])), curMax, colors[rdx]))
            else:
                R.write('plot(x=c(%s), y=c(%s), type="l", lwd=2.5, ylim = c(0,%s), bty="n", axes = FALSE, ann = FALSE, col=%s)\n' % (','.join(map(str, sorted(repData[typ][rNum].keys()))), ','.join(map(str, [repData[typ][rNum][x] for x in sorted(repData[typ][rNum].keys())])), curMax, colors[rdx]))
                R.write('axis(1, cex.axis=1, las=1, cex.axis=1.5)\n')
                R.write('axis(2, cex.axis=1, las=1, cex.axis=1.5)\n')
                R.write('mtext(side = 1, "Number of grougs of sample replicates", line = 3.5, cex=1.5)\n')
                R.write('mtext(side = 2, "Number of OTUs", line = 3.5, cex=1.5)\n')
                R.write('mtext(side = 3, "%s", cex = 3)\n' % typ)
    R.write('plot.new()\n')
    R.write('legend("left", inset=c(-2,0), legend=c(%s), title="Number of\nreplicates\nper group", col=c(%s), lty=1, cex=2, lwd=2, bty="n")\n' % (','.join(map(str, sorted(repData['Exact'].keys()))), ','.join(colors)))
    #R.write('dev.off()\n')
    R.close()
    subprocess.call(['R', '-f', Rout, '--vanilla', '--slave'])
    return pdf

def get_filt(n1, n2, only, intSplt, reg_idx, thresh):
    filt = {0: [], 1: []}
    nLine = len(intSplt)
    if n1 > n2:
        # treat each sample separately
        if only:
            for vdx, val in enumerate(intSplt):
                if vdx in reg_idx:
                    if val > thresh:
                        filt[1].append(val)
                        filt[0].append(0)
                    else:
                        filt[0].append(val)
                        filt[1].append(0)
                else:
                    filt[1].append(val)
                    filt[0].append(0)
        # if '--only' option not activated: keep all samples counts
        else:
            filt[1] = intSplt
            filt[0] = [0]*nLine
    else:
        filt[0] = intSplt
        filt[1] = [0]*nLine
    return filt


def filt_replicates(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: if the OTU/ISU does not reach the threshold in:
    - a minimum number of replicates per group of replicates
    - a minimum number of groups of replcaites
    """
    intSplt = get_intSplt(splt)
    filt = {0: [], 1: []}
    # get arguments
    repsList = arguments[0]
    # required number replicates per group
    nReps = arguments[1]
    # refine required number of groups of replicates
    nGrps = arguments[2]
    if isinstance(nGrps, float):
        nGrps = int(len(repsList)*nGrps)
    across = arguments[3]
    only = arguments[4]
    sample_repname = arguments[5]
    idx_names_replicates = arguments[6]
    idx_replicates = arguments[7]
    flat_idx_replicates = arguments[8]

    # corresponding number of reads
    reads_replicates = [[intSplt[y] for y in x] for x in idx_replicates]
    reads_in_replicates = [len([y for y in x if y]) for x in reads_replicates if len([y for y in x if y]) >= 2]

    # if '--sum' option activated
    if across:
        # number of reads per group of replicates
        sums_grps = [sum(x) for x in reads_replicates]
        # groups with enough reads
        grps_above_thresh = [1 if x > thresh else 0 for x in sums_grps]
        # number of groups with enough reads
        nGrps_above_thresh = sum(grps_above_thresh)
        # initialize as keep all OTU/ISU reads
        filt[1] = intSplt
        filt[0] = [0]*len(splt)
        # check number of replicate groups that sums of reads across replicates are above the threshold
        # if enough:
        if nGrps_above_thresh >= nGrps:
            # if '--only' option activated: the filter should apply only the selected samples
            if only:
                # remove the reads of the OTU/ISU in the replicates of the group that do not sum above the threshold
                # NB: the samples not declared as replicate are not filtered at all
                for grpdx, grp in enumerate(grps_above_thresh):
                    for repdx in idx_replicates[grpdx]:
                        if grp == 0:
                            filt[1][repdx] = 0
                            filt[0][repdx] = intSplt[repdx]
        # if not enough:
        else:
            # if '--only' option activated: the filter should apply only the selected samples
            if only:
                # remove the reads of the OTU/ISU in the replicates of the group that do not sum above the threshold
                # NB: the samples not declared as replicate are not filtered out
                for idx in flat_idx_replicates:
                    filt[1][idx] = 0
                    filt[0][idx] = intSplt[idx]
            else:
                filt[0] = intSplt
                filt[1] = [0]*len(splt)
    # if '--sum' option not activated
    else:
        # list of binary lists to encode replicates with enough reads per group
        replicates_above_thresh = [[1 if reads_rep > thresh else 0 for reads_rep in reads_reps] for reads_reps in reads_replicates]
        # number of replicates with enough reads per group
        n_replicates_above_thresh_per_group = [sum(x) for x in replicates_above_thresh]
        # groups that have enough replicates with enough reads
        groups_above_nReps = [1 if x >= nReps else 0 for x in n_replicates_above_thresh_per_group]
        # number of groups that have enough replicates with enough reads
        nGrps_above_nReps = sum(groups_above_nReps)
        # init as everythin kept (as if '--only' was set to False)
        filt[1] = intSplt
        filt[0] = [0]*len(splt)
        # if enough groups that have enough replicates with enough reads
        if nGrps_above_nReps >= nGrps:
            # if '--only' option activated: the filter should apply only the selected samples
            if only:
                # for each group and for each replicate in the group
                for grpdx, grp in enumerate(groups_above_nReps):
                    for repdx in idx_replicates[grpdx]:
                        # update to zero the replicate that would not have enough replicates
                        if grp == 0:
                            filt[1][repdx] = 0
                            filt[0][repdx] = intSplt[repdx]
        # if not enough groups that have enough replicates with enough reads
        else:
            # if '--only' option activated: the filter should apply only the selected samples
            if only:
                # update to zero each replicate in the list if replicates
                for idx in flat_idx_replicates:
                    filt[1][idx] = 0
                    filt[0][idx] = intSplt[idx]

            # if '--only' option not activated: keep all samples counts
            else:
                # update to zero the entire OTU
                filt[0] = intSplt
                filt[1] = [0]*len(splt)
    return filt, reads_in_replicates

def get_replicates(filin, code, samples, table):
    """
    Parse the file with as lines the tab or comma separated replicate samples of a group of replicates
    Return a list of lists: replicates lists, per group of sample replicates
    """
    # init emply lists
    reps = []
    reps_in = []
    reps_out = []
    reps_heads = []
    # get the table samples names
    sNames = [x[1] for x in samples]
    # parse the replicates file
    for i in open(filin, 'rU'):
        # get the list of replicates of the current group of replicates
        i = '\t'.join(i.strip().split())
        if ',' in i:
            curReps = i.strip().split(',')
        elif '\t' in i:
            curReps = i.strip().split('\t')
        else:
            continue
        # collect the current replicates group in the main list
        reps.append(curReps)
        # collect the replicates of the file that are in/not in the table
        newReps = []
        for curRep in curReps:
            if curRep not in sNames:
                reps_out.append(curRep)
            else:
                newReps.append(curRep)
        if reps_out:
            curHead = curReps[0]
            curTail = curReps[1:]
            if sum([1 for x in curTail if curHead in x]) == len(curTail):
                reps_heads.append(curReps)
        # put the group of replicates appearing in the table in secondary list
        if newReps:
            reps_in.append(newReps)
    # if there are replicate samples not in the table
    sample_repname = 0
    if reps_out:
        if sorted(reps_out) == sorted([x[0] for x in reps_heads]):
            sample_repname = 1
            retReps = reps_heads
        else:
            print 'Sample name(s) provided in the replicates file not present in the data table (%s):' % table
            for rep in sorted(reps_out):
                print ' -', rep
            if len(reps_out) != len(list(set(reps_out))):
                print 'and some samples appear in several groups of replicates:'
                for rep in sorted(reps_out):
                    if reps_out.count(rep)>1:
                        print ' *', rep
            # ask the user to continue although some samples of the file are not in the table
            # and if yes, return the secondary list of groups of replicates
            userChoice = raw_input('Continue with the samples present in the table samples? <y/n>')
            userChoice = userChoice.upper()
            if userChoice.startswith('Y'):
                print 'using %s groups of replicates' % len(reps_in)
                for grpx, grp in enumerate(reps_in):
                    print grpx, ', '.join(grp)
                time.sleep(3)
                retReps = reps_in
            else:
                print 'Exiting'
                sys.exit()
    # return the main list of groups of replicates
    else:
        retReps = reps
    idx_names_replicates, idx_replicates, flat_idx_replicates = get_idx_repData(retReps, sample_repname, samples)
    return retReps, sample_repname, idx_names_replicates, idx_replicates, flat_idx_replicates


def get_idx_repData(repsData, argIdx, samples):
    # nested list of replicates lists in terms of sample indices instead of samples names
    idx_names_replicates = [[[r for r in samples if r[1] == rep][0] for rep in reps[argIdx:]] for reps in repsData]
    idx_replicates = [[[r[0] for r in samples if r[1] == rep][0] for rep in reps[argIdx:]] for reps in repsData]
    # flatten the list of lists into one list
    flat_idx_replicates = sum(idx_replicates, [])
    return  idx_names_replicates, idx_replicates, flat_idx_replicates


def filt_choice(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: if the OTU/ISU does not reach the threshold in:
    - the samples corresponding to the regex / index in the initial sample list
    - potentially, a minimum number of these samples
    - the sum of these selected samples
    """
    intSplt = get_intSplt(splt)
    # subsamples selection
    subsamples = arguments[0]
    # numbers on samples info (default to 100)
    percent = arguments[1]
    across = arguments[2]
    only = arguments[3]
    choice_idx, choice_intSplt, minSamples = get_sample_info(intSplt, subsamples, percent)
    # if '--sum' option activated
    if across:
        # sum of abundances in all selected samples
        sum_choice_intSplt = sum(choice_intSplt)
        # check if this sum if above the threshold
        filt = get_filt(choice_intSplt, thresh, only, intSplt, choice_idx, thresh)
    # if '--sum' option not activated
    else:
        # abundances of the subsamples that are above the threshold
        choice_vals_above_thresh = [x for x in choice_intSplt if x > thresh]
        # number of subsamples that have an abundance above the threshold
        n_vals_above_thresh = len(choice_vals_above_thresh)
        # if number of subsamples that have an abundance above the threshold is above the number of subsamples threshold
        filt = get_filt(n_vals_above_thresh, minSamples-1, only, intSplt, choice_idx, thresh)
    # for this filter also return the subsamples under a new key
    return filt, 0

def filt_presence(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: if the OTU/ISU does not reach the threshold in a:
    - a minimum number of samples
    - the sum of a minimum number of samples
    """
    percent = arguments[0]
    only = arguments[1]
    intSplt = get_intSplt(splt)
    reg_idx, reg_intSplt, minSamples = get_sample_info(intSplt, samples, percent)
    # transform entry abundance in binary presence/absence
    nSup = sum([1 for x in reg_intSplt if x])
    filt = get_filt(nSup, minSamples-1, only, intSplt, reg_idx, thresh)
    return filt, 0

def filt_minimum(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: if the OTU/ISU does not reach the threshold in:
    - a minimum number of samples
    - the sum of a minimum number of samples
    """
    intSplt = get_intSplt(splt)
    percent = arguments[0]
    across = arguments[1]
    only = arguments[2]
    reg_idx, reg_intSplt, minSamples = get_sample_info(intSplt, samples, percent)
    # if '--sum' option activated
    if across:
        # calculate the sum based on "nSamples" first highest numbers 
        sum_max_half = sum(sorted(reg_intSplt)[::-1][:len(samples)])
        filt = get_filt(sum_max_half, thresh, only, intSplt, reg_idx, thresh)
    # if '--sum' option not activated
    else:
        # numbers above the threshold
        n_samples_above_thresh = [x for x in reg_intSplt if x > thresh]
        # number of numbers above the thresholds
        sum_max_half_thresh = len(n_samples_above_thresh)
        filt = get_filt(sum_max_half_thresh, minSamples-1, only, intSplt, reg_idx, thresh)
    return filt, 0

def filt_dataset(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: across the entire dataset
    """
    intSplt = get_intSplt(splt)
    only = arguments[0]
    reg_idx, reg_intSplt, minSamples = get_sample_info(intSplt, samples)
    # calculate the sum on the subsamples only
    spltSum = sum(reg_intSplt)
    filt = get_filt(spltSum, thresh, only, intSplt, reg_idx, thresh)
    return filt, 0

def filt_sample(splt, samples, thresh, arguments):
    """
    Return the filtered version of the input OTU/ISU line
    Method: for each sample separately
    """
    intSplt = get_intSplt(splt)
    reg_idx, reg_intSplt, minSamples = get_sample_info(intSplt, samples)
    filt = get_filt(1, 0, True, intSplt, reg_idx, thresh)
    return filt, 0

def get_intSplt(splt):
    try:
        intSplt = map(int, splt)
    except ValueError:
        try:
            intSplt = map(float, splt)
        except ValueError:
            print 'Table contain non-numeric values\nExiting'
            sys.exit()
    return intSplt

def get_sample_info(intSplt, samples, percent=100):
    reg_idx = [s[0] for s in samples]
    reg_intSplt = [intSplt[x] for x in reg_idx]
    minSamples = len(samples) * (percent / 100.)
    return reg_idx, reg_intSplt, minSamples

def get_subsamples(samples, regex_indices):
    """
    Use the indices and/or regex info passed as -mode argument to the method 'choice' to make the selection of subsamples
    Return the list of subsamples selection in the same format as the list of passed samples
    """
    # init empy lists
    subsamples = []
    regex = []
    vals = []
    # for each indices/regex input
    for regex_index in regex_indices:
        # if it contains a hyphen (possible indices list)
        if '-' in regex_index:
            # check if all hyphen-separated fields are integers and collect them
            for ind in regex_index.split('-'):
                if ind:
                    try:
                        val = int(ind)
                        if val not in vals:
                            vals.append(val)
                    # if one is not an integer: empty the integer list and skip to next input
                    except ValueError:
                        regex.append(regex_index)
                        vals = []
                        break
        # if no hyphen: interpret as a regex
        else:
            regex.append(regex_index)
    # if there was a list of indices in the input
    if vals:
        # samples with this indices in the table samples
        vals = sorted(list(set(vals)))
        subsamples = [x for x in samples if (x[0]+1) in vals]
        # out-of-range input values (no possible sample)
        unused_indices = [x for x in vals if (x-1) not in [y[0] for y in samples]]
        if unused_indices:
            print '*** Warning: some sample indices are out of range: %s (ignored)' % ', '.join(map(str, unused_indices))
    # if there was at least one regex in the input
    if regex:
        # collect subsamples as samples that match any input regex
        detect = []
        for sample in samples:
            for rgx in regex:
                if re.search(rgx, sample[1]):
                    if sample not in subsamples:
                        subsamples.append(sample)
                        detect.append(rgx)
                    break
        # show if no subsample selected or which regex was unsuccessful at detecting a sample
        if len(subsamples)==0 or len(list(set(regex)^set(detect))):
            print "*** Warning: no sample detected using regex: '%s' (ignored):" % "', '".join(sorted(list(set(detect)^set(regex))))
    return subsamples

def get_method(meth, mode, across, only, tableCode, filin, samples, all_samples):
    """
    Define the filtering function to use based on user inputs by checking interactions of passed options
    Note that the abundance threshol argument always applies - see specific filtering functions
    Returns a list which first item corresponf to the function to call, followed by arguments for this function
    """
    maxReps = 0
    max_nGrps = 0
    subsamples = []
    replicates = []
    sample_repname = 0
    # abundance-based filters
    if meth == 'simple':
        if across:
            # across-dataset filtering / may apply on a set of selected samples
            filtering = [filt_dataset, only]
        else:
            # per-sample filtering / no argument
            filtering = [filt_sample]
    # occurence-based filters
    elif meth == 'minimum' or meth == 'presence':
        # these filtering methods can only take one "mode" argument: minimum percent of samples occurrences required
        if mode == []:
            # make default value if lacking
            mode = [10]
        if len(mode) == 1:
            # check if the passed value
            try:
                # is an integer
                val = int(mode[0])
                # is a percentage
                if val in range(1, 101):
                    if meth == 'presence':
                        # filtering function and its arguments: percent value of the samples
                        filtering = [filt_presence, val, only]
                    else:
                        # filtering function and its arguments: percent value of the samples; across dataset; only on samples not the entire OTU
                        filtering = [filt_minimum, val, across, only]
                else:
                    print "Option error: method 'minimum' incompatible with the input value of %s (outside range 1-100)\nExiting" % mode[0]
                    sys.exit()
            except ValueError:
                print "Option error: method 'minimum' incompatible with the value '%s'\nExiting" % mode[0]
                sys.exit()
        else:
            print "Too many arguments to -mode for -method 'minimum'\nExiting"
            sys.exit()
    # on-selected-samples filter
    elif meth == 'choice':
        # necessitates at least one mode argument: the indices/regex for samples selection
        if mode == []:
            print 'No indices/regex passed as argument...\nExiting'
            sys.exit()
        # init the percent of the selected samples to filter on
        val = 100
        # if there are 2 or more indices/regex
        if len(mode) >= 2:
            # check if the last argument
            try:
                # is an integer
                val = int(mode[-1])
                # is a percentage
                if val in range(1,101):
                    # then the other arguments are interpreted as indices/regex
                    mode = mode[:-1]
                else:
                    print "Warning: mode value %s for method 'choice' outside range 1-100: interpreted as indices/regex" % mode[-1]
            except ValueError:
                print "Warning: mode value %s for method 'choice' not a numeric: interpreted as indices/regex" % mode[-1]
        # get the chosen subsamples based on indices/regex
        samples = get_subsamples(samples, mode)
        subsamples = get_subsamples(samples, mode)
        # filtering function and its arguments: subsamples, percent value of the subsamples; across subsamples; only on samples not the entire OTU
        filtering = [filt_choice, subsamples, val, across, only]
    # replicates-based filtering
    elif meth == 'replicates':
        # necessitates at least one mode argument: the indices/regex for samples selection
        if mode == []:
            print "No regex passed as argument for the method 'replicates'.\nExiting"
            sys.exit()
        # cannot take more than 3 arguments
        elif len(mode)>3:
            print "Too many arguments to '-mode' for the method 'replicates'.\nExiting"
            sys.exit()
        # the passed first argument must be the file providing which samples are replicates
        if os.path.isfile(mode[0]):
            # init the minimum number of number of groups of sample replicates
            nGroups = 1
            # init the minimum number of replicates per group
            nReps = 2
            # replicates as a list of lists
            replicates, sample_repname, subsamples, idx_replicates, flat_idx_replicates = get_replicates(mode[0], tableCode, all_samples, filin)
            max_nGrps = len(replicates)
            # numbers of replicates in each group of replicates
            if sample_repname:
                nReplicates = [len(x)-1 for x in replicates]
            else:
                nReplicates = [len(x) for x in replicates]
            # max number of replicates in a group of replicates
            maxReps = max(nReplicates)
            # parse the other arguments
            for mdx, mod in enumerate(mode[1:]):
                if mdx:
                    # check if the nReps input argument
                    try:
                        # is an integer
                        nReps = int(mod)
                        # update max number of replicates in a group of replicates
                        if nReps > maxReps:
                            nReps = maxReps
                        # keep it to 2 minimum
                        if nReps == 1:
                            print 'Minimum number of replicates for filtering based on replicates is 2'
                            nReps = 2
                    except ValueError:
                        print "Option error: method 'replicates' incompatible with input %s (not integer)\nExiting" % mod
                        sys.exit()
                else:
                    # check if the nGroups input argument
                    try:
                        # is an integer
                        nGroups = int(mod)
                    except ValueError:
                        try:
                            # is a float
                            nGroups = float(mod)
                        except ValueError:
                            print "Option error: method 'replicates' incompatible with input %s (not numeric)\nExiting" % mod
                            sys.exit()
                    # upper max number of group of replicates
                    if nGroups > max_nGrps:
                        nGroups = max_nGrps
                # for the first argument
            # prints
            if isinstance(nGroups, float):
                print "Using %s as the number of replicates per sample, and %s as the percent of groups of sample replicates" % (nReps, nGroups)
            else:
                print "Using %s as the number of replicates per sample, and %s as the number of groups of sample replicates" % (nReps, nGroups)
            # filtering function and its arguments:
            #   replicates per group of replicate
            #   number of replicates per groups required
            #   number of groups of replicates required with the above number of replicates per group
            #   across subsamples
            #   only on samples not the entire OTU
            filtering = [filt_replicates, replicates, nReps, nGroups, across, only, sample_repname, subsamples, idx_replicates, flat_idx_replicates]
        else:
            print "Need file containing replicates information to run method 'replicates'"
            sys.exit()
    return filtering, samples, maxReps, max_nGrps, subsamples, sample_repname, replicates

def get_samples(filin, code, sep, regex):
    """
    Parses the first line if the input table to find samples
    Return all the samples found, as well as the samples that could be selected by user input (indices/regex)
    """
    # table format encoding (see function check_table())
    head_code = code[0]
    meta_code = code[1]
    # parse input OTU/ISU table
    for lindx, line in enumerate(open(filin, 'rU')):
        split_line = line.strip().split(sep)
        if len(''.join(split_line)): # skip empty line
            # line field that may correspond to count data or sample names
            reads = split_line[meta_code:]
            # line header (may contain OTU name / taxonomy / ...)
            meta = [x.strip() for x in split_line[:meta_code]]
            # only look at first, header line
            if lindx == 0:
                # if the first line is supposed to be a header
                if head_code:
                    # init samples as all samples
                    all_samples = [(xdx, x.strip()) for xdx, x in enumerate(reads)]
                    samples = [(xdx, x.strip()) for xdx, x in enumerate(reads)]
                    # refine samples according to user input
                    if regex[0]:
                        samples = [x for x in all_samples if x[0] in regex[0]]
                    for rgx in regex[1]:
                        samples = [x for x in all_samples if re.search(rgx, x[1])]
                # if the first line is data
                else:
                    # init samples as all samples, but based on a generic sample name with an index
                    samples = [(xdx, 'Sample#%s' % x) for xdx,x in enumerate(range(1, (len(split_line)-meta_code)))]
                    all_samples = [(xdx, 'Sample#%s' % x) for xdx,x in enumerate(range(1, (len(split_line)-meta_code)))]
                    # only possible to refine samples according to indices passed by user
                    if regex[0]:
                        samples = [x.strip() for x in all_samples if x[0] in regex[0]]
    return all_samples, samples



def write_line(filtered_line, meta, OTU, ot, samples, all_samples, name_otus, merge, lineIdx):
    """
    Write filtered line in output table
    """
    curFilt_1 = list(filtered_line[1])
    # only add data on the selected samples if option '--selection' is used
    curFilt_1 = [curFilt_1[x[0]] for x in all_samples]
    if sum(curFilt_1):
        line = []
        if meta:
            line.append('\t'.join(meta))
        elif name_otus:
            line.append('%s%s' % (name_otus, OTU))
        if merge:
            outStr = [sum([curFilt_1[y] for y in x])/float(len(x)) for x in lineIdx]
        else:
            outStr = curFilt_1
        line += map(str, outStr)
        ot.write('%s\n' % '\t'.join(line))


def check_table(filin, sep, head, meta):
    """
    Identify the type of OTU table: with/without header line, with/without OTU names (e.g. Taxonomy)
    by parsing two lines:
        the first line that could be a header
        the second line which always should be a data line representing the whole table structure for OTUs
    """
    tableCode = {}
    for idx,i in enumerate(open(filin, 'rU')):
        j = i.strip().split(sep)
        # not first line
        if idx:
            # get the lists of indices for (i) numeric fields (ii) non-numeric fields
            if meta == None:
                numbers2, strings2 = fill_fields(j)
            else:
                numbers2, strings2 = range(meta, len(j)), range(meta)
            if head:
                otu_taxonomy = 0
                if strings2:
                    otu_taxonomy = max(strings2)+1
                return 1, otu_taxonomy
            else:
                # first line contained numbers
                if numbers:
                    # first line was data line (i.e. no header)
                    if numbers == numbers2:
                        # only data (i.e. no text at all)
                        if strings2 == strings == []:
                            # -header -OTU name -taxonomy 
                            return 0, 0
                        # text and data
                        else:
                            # -header +OTU name -taxonomy 
                            if len(strings2) == 1:
                                return 0, 1
                            # -header +OTU name +taxonomy 
                            elif len(strings) == 2:
                                return 0, 2
                # first line was header
                otu_taxonomy = 0
                if strings2:
                    otu_taxonomy = max(strings2)+1
                return 1, otu_taxonomy
        # first line
        else:
            # count number of fields
            nFields = len(j)
            # get the lists of indices for (i) numeric fields (ii) non-numeric fields
            if meta == None:
                numbers, strings = fill_fields(j)
            else:
                numbers, strings = range(meta, len(j)), range(meta)
    if tableCode:
        return tableCode
    else:
        print 'No table structure recognized\nExiting'
        sys.exit()

def fill_fields(lineSplit):
    """Returns two lists of indices: for the numeric values / for the non-numeric values"""
    numbers = []
    strings = []
    for fdx, field in enumerate(lineSplit):
        try:
            val = int(field)
            numbers.append(fdx)
        except ValueError:
            try:
                val = float(field)
                numbers.append(fdx)
            except ValueError:
                strings.append(fdx)
    return numbers, strings

def get_samples_indices(reg):
    """
    Returns a dict with separate integer and non-integer inputs,
    for samples selection based on indices or regex
    """
    regex = {0: [], 1: []}
    for r in reg:
        try:
            val = int(r)
            regex[0].append(val)
        except ValueError:
            regex[1].append(r)
    return regex

def write_header_stats(args, statout, filout):
    """
    Stats file output (to be elaborated)
    """
    curTime = time.strftime("%D_%H:%M:%S", time.localtime())
    statoutPath = os.path.abspath(statout)
    o = open(statoutPath, 'w')
    o.write('# script: %s\n' % os.path.abspath(sys.argv[0]))
    o.write('# date: %s\n' % curTime)
    o.write('# input table: %s\n' % args['i'][0])
    o.write('# output table: %s\n' % filout)
    o.write('# abundance threshold: %s\n' % args['t'])
    o.write('# method: %s\n' % args['meth'])
    o.write('# mode: %s\n' % '\t'.join(args['mode']))
    o.write('# across dataset: %s\n' % args['sum'])
    if args['x']:
        o.write('# regex to subsample: %s\n' % ', '.join(args['x']))
    if args['meth'] in ['choice', 'replicates']:
        if args['meth'] == 'choice':
            o.write('# output table for chosen samples: %s\n' % args['selection'])
        o.write('# filter only in chosen samples: %s\n' % args['only'])
    return o, statoutPath


def write_header_table(filout, samples, all_samples, name_otus, select, merge, code, subsamples, sample_repname, replicates):
    """
    Write output table with header
    """
    # reconstitute the header line
    header = []
    lenMeta = code[1]
    cols = ['ID', 'Taxonomy']
    if lenMeta:
        for i in range(lenMeta):
            if len(cols) == i:
                header.append('\t')
            else:
                header.append(cols[i])
    elif name_otus:
        header.append(name_otus)
    if select:
        try:
            all_samples = sorted(sum(subsamples, []), key = lambda x: int(x[0]))
        except TypeError:
            all_samples = sorted(subsamples, key = lambda x: int(x[0]))
    ctrl = {}
    lineIdx = []
    for i in sorted(all_samples, key = lambda x: int(x[0])):
        if merge:
            if sample_repname:
                samrepsnames = list(replicates)
            else:
                print 'The replicates method did not find sample generic name in the replicates file\nExiting'
                sys.exit()
            for repsams in samrepsnames:
                newsam = repsams[0]
                oldsams = repsams[1:]
                if i[1] in oldsams and ctrl.has_key(newsam) == False:
                    header.append(newsam)
                    ctrl[newsam] = 1
                    lineIdx.append([x[0] for x in all_samples if x[1] in oldsams])
                    break
        else:
            header.append(i[1])
    # start writing

    filoutPath = os.path.abspath(filout)
    o = open(filoutPath, 'w')
    o.write('%s\n' % '\t'.join(header))
    return o, lineIdx, filoutPath


#cProfile.run('abundance_filters()')
abundance_filters()
