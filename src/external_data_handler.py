import logging
import re
import sys

from breakpoint_graph import *

complementary_nucleotide = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "t", "c": "g", "g": "c", "t": "a"}


def rev_complement(seq):
    if seq:
        return ''.join([complementary_nucleotide[a] for a in seq[::-1]])

    return seq


def vcf_var_to_bp_pair(chrom, pos, ref, alt):
    rlen = len(ref)
    if alt.startswith("[") or alt.startswith("]"):
        da = alt[0]
        db = alt.rstrip("acgtACGTN")[-1]
        if not da == db and (da == ']' or da == '['):
            logging.debug("Alt allele in VCF was not a breakend SV alt allele: " + str((chrom, pos, ref, alt)))
            return None, None

        hom_seq = re.split(']|\[', alt)[-1][rlen:]
        alt_end = re.split(']|\[', alt)[1]
        chrom2, pos2 = alt_end.rsplit(":")
        if da == ']':
            v1 = breakpoint_vertex(chrom2, int(pos2), 1)
        else:
            v1 = breakpoint_vertex(chrom2, int(pos2), -1)

        v2 = breakpoint_vertex(chrom, int(pos), -1)  # AA notation is -1 for v2 if it is forward

    elif alt.endswith("[") or alt.endswith("]"):
        db = alt[-1]
        da = alt.lstrip("acgtACGTN")[0]
        if not da == db and (da == ']' or da == '['):
            logging.debug("Alt allele in VCF was not a breakend SV alt allele: " + str((chrom, pos, ref, alt)))
            return None, None

        hom_seq = re.split(']|\[', alt)[0][rlen:]
        alt_end = re.split(']|\[', alt)[1]
        chrom2, pos2 = alt_end.rsplit(":")
        if db == ']':
            v2 = breakpoint_vertex(chrom2, int(pos2), 1)  # AA notation is 1 for v2 if it is reverse
        else:
            v2 = breakpoint_vertex(chrom2, int(pos2), -1)  # AA notation is -1 for v2 if it is forward

        v1 = breakpoint_vertex(chrom, int(pos), 1)

    else:
        logging.debug("Alt allele in VCF was not a breakend SV alt allele: " + str((chrom, pos, ref, alt)))
        return None, None

    return (v1, v2), hom_seq


# determine the number of read pairs supporting an SV call
# this does not cound split-reads (single-end), as that is not what AA uses.
def get_sv_read_pair_support_from_vcf(fd, header_fields):
    ilist = fd['INFO'].rsplit(";")
    info_dict = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in ilist if "=" in x}
    pr_support = None
    # check PE and SR (delly, lumpy) - SR is single-end version - not checking for now
    if 'PE' in info_dict:
        pr_support = int(info_dict['PE'])

    # check RP or REF (gridss) - REF is the single-end version - not checking for now
    elif 'RP' in info_dict:
        pr_support = int(info_dict['RP'])

    elif 'FORMAT' in header_fields:
        sname = header_fields[-1]
        fkeys = fd['FORMAT'].rsplit(":")
        slist = fd[sname].rsplit(":")
        if len(fkeys) != len(slist):
            logging.error("Format field defined a different number of fields than found in data! Skipping SV " + fd['ID'])
            logging.error(str(fkeys))
            logging.error(str(slist))
            return 0

        format_dict = dict(zip(fkeys, slist))
        # check format field for PR/SR (manta)
        if 'PR' in format_dict:
            pr_support = int(format_dict['PR'])

        # check format field for DR (svaba)
        elif 'DR' in format_dict:
            pr_support = int(format_dict['DR'])

    if pr_support is None:
        logging.warning("Could not find paired-end read support count for SV ID: " + fd['ID'])
        pr_support = 0

    return pr_support


# reads a .vcf (or .vcf.gz) and converts it to a list of dictionaries (keys are header values)
# returns the dictionary and a list of the header fields
def read_vcf(vcf_file, filter_by_pass=True):
    dlist = []
    if vcf_file.endswith('.gz'):
        import gzip
        opener = gzip.open
    else:
        opener = open

    field_warn = False
    with opener(vcf_file, 'rt') as infile:
        header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        for line in infile:
            if line.startswith("#CHROM"):
                header_fields = line[1:].rstrip().rsplit()

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fd = dict(zip(header_fields, fields))
                # check that it's a single-sample VCF which should have 9 header field names and 1 sample name
                # in the header.
                if "FORMAT" in fd and len(header_fields) > 10:
                    logging.error("VCF appears to contain multiple samples or non-standard fields. Please provide a"
                                  " single-sample VCF with appropriate formatting\n")
                    sys.exit(1)

                # this is because SVABA does not obey VCF4.2 format :(
                if len(fields) > len(header_fields):
                    field_warn = True
                    fd[header_fields[-1]] = fields[-1]

                if (not filter_by_pass) or (filter_by_pass and fd['FILTER'] == "PASS"):
                    dlist.append(fd)

    if field_warn:
        logging.warning("VCF file has more fields than specified by the header!")

    return dlist, header_fields


def sv_vcf_to_bplist(vcf_file, filter_by_pass=True):
    dlist, hf = read_vcf(vcf_file, filter_by_pass)
    vcf_dnlist = []
    seen_bp_set = set()
    for fd in dlist:
        bp_pair, hom_seq = vcf_var_to_bp_pair(fd['CHROM'], fd['POS'], fd['REF'], fd['ALT'])
        if bp_pair:
            if bp_pair not in seen_bp_set:
                seen_bp_set.add(bp_pair)
                seen_bp_set.add(bp_pair[::-1])
                support = get_sv_read_pair_support_from_vcf(fd, hf)
                if not hom_seq:
                    hom_seq, hom_len = None, None
                else:
                    hom_len = len(hom_seq)
                bref = breakpoint_edge(bp_pair[0], bp_pair[1], hom_seq=hom_seq, hom=hom_len, externally_called=True)
                bref.edge_type = bref.type()
                vcf_dnlist.append((bref, support))
                bref_rev = breakpoint_edge(bp_pair[1], bp_pair[0], hom_seq=rev_complement(hom_seq), hom=hom_len, externally_called=True)
                bref_rev.edge_type = bref.type()
                vcf_dnlist.append((bref_rev, support))


    logging.info("Read " + str(len(vcf_dnlist)) + " SV calls from " + vcf_file)
    return vcf_dnlist
