import logging
import re

from breakpoint_graph import *

def vcf_var_to_bp_pair(chrom, pos, ref, alt):
    rlen = len(ref)
    if alt.startswith("[") or alt.startswith("]"):
        if alt.startswith("["):
            v1s = 1
        else:
            v1s = -1

        d = alt.rstrip("acgtACGTN")[-1]
        if d == "[":
            v2s = -1
        elif d == "]":
            v2s = 1
        else:
            logging.warning("Alt allele in VCF was not a valid SV alt allele: " + str((chrom, pos, ref, alt)))
            return None, None

        hom_seq = re.split(']|\[', alt)[-1][rlen:]

    elif alt.endswith("[") or alt.endswith("]"):
        if alt.endswith("["):
            v2s = -1
        else:
            v2s = 1

        d = alt.lstrip("acgtACGTN")[-1]
        if d == "[":
            v1s = 1
        elif d == "]":
            v1s = -1
        else:
            logging.warning("Alt allele in VCF was not a valid SV alt allele: " + str((chrom, pos, ref, alt)))
            return None, None

        hom_seq = re.split(']|\[', alt)[0][rlen:]

    else:
        logging.warning("Alt allele in VCF was not a valid SV alt allele: " + str((chrom, pos, ref, alt)))
        return None, None

    dest = re.split(']|\[', alt)[1]
    chrom2, pos2 = dest.rsplit(":")
    v1 = breakpoint_vertex(chrom, int(pos), v1s)
    v2 = breakpoint_vertex(chrom2, int(pos2), v2s)
    return (v1, v2), hom_seq


# determine the number of read pairs supporting an SV call
# this does not cound split-reads (single-end), as that is not what AA uses.
def get_sv_read_pair_support_from_vcf(fd, header_fields):
    pass
    # check PE and SR (delly, lumpy)

    # check RP or REF (gridss)

    # check format field for PR/SR (manta)


# reads a .vcf (or .vcf.gz) and converts it to a list of dictionaries (keys are header values)
# returns the dictionary and a list of the header fields
def read_vcf(vcf_file):
    dlist = []
    if vcf_file.endswith('.gz'):
        import gzip
        opener = gzip.open
    else:
        opener = open

    with opener(vcf_file, 'rt') as infile:
        header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        for line in infile:
            if line.startswith("#CHROM"):
                header_fields = line[1:].rstrip().rsplit()

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                # TODO: check that it's a single-sample VCF
                fd = dict(zip(fields, header_fields))
                dlist.append(fd)

    return dlist, header_fields


def sv_vcf_to_bplist(vcf_file):
    dlist, hf = read_vcf(vcf_file)
    vcf_dnlist = []
    seen_bp_set = set()
    for fd in dlist:
        bp_pair, hom_seq = vcf_var_to_bp_pair(fd['CHROM'], fd['POS'], fd['REF'], fd['ALT'])
        if bp_pair:
            if bp_pair not in seen_bp_set:
                seen_bp_set.add(bp_pair)
                seen_bp_set.add(bp_pair[::-1])
                # get supports
                support = get_sv_read_pair_support_from_vcf(fd, hf)
                bref = breakpoint_edge(bp_pair[0], bp_pair[1], hom_seq=hom_seq, hom=len(hom_seq))



    # bp1 and bp2 are breakpoint vertices
    #                     brec_refine = self.refine_discordant_edge(breakpoint_edge(bp1c, bp2c))
    #                     brec = brec_refine[0]
    #
    #                         if self.edge_passes_filters(vl, brec):
    #                             dnlist0.append((brec, len([(v[1], v[0]) for v in vl])))

    # init_dnlist =

    return vcf_dnlist