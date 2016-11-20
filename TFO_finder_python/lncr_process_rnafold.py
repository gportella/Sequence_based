#! /opt/local/bin/python2.7
""" Read a multifasta file with lncrna (actually dna transcripts)
    and find putative TFOs
"""

import sys
import argparse
import re
import subprocess
from itertools import izip
from constants_nn import NN_DH, NN_DG, COEF_PH, DUPLEX_CONC, VAL_PH,\
    CONC_TFO, SPEC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import numpy as np

R = 0.0019858

ROLLO = '\n \
Open a fasta file, parse the sequences and predict \
possible T,C triplex \n\
'
WHOWHEN = 'Guillem Portella, v2.0, 08-2016'


def find_tfo_from_mfasta(fasta_records, b_rnafold=False, b_very_verbose=False):
    """ From an individual fasta entry record
        search for tfo.
    :rtype: a list of LncRna objects containing the TFO info
    """
    found_tfo = []
    for idx, record in enumerate(fasta_records):
        tfo_record = LncRna(record, idx, b_rnafold, b_very_verbose)
        if tfo_record.found:
            found_tfo.append(tfo_record)
    return found_tfo


class LncRna(object):
    """  Store and built all the information of TFOs

    """

    def __init__(self, record, index, b_rnafold=False, b_very_verbose=False):
        self.found = False
        self.struct_computed = False
        self.struct = []  # if b_rnafold is set true, stores the sec. structure
        self.lnctfos = []  # will be a list of Tfos objects
        self.find_tfo(record, index, b_rnafold, b_very_verbose)

    class Tfos(object):
        """ Empty, I know it sucks because it is just a container
            of data, pylint, but give me a break
        """

        def __init__(self):
            self.position = [0] * 2
            self.id_lncrna = ""
            self.seq = ""
            self.index_id = 0  # why do I need this one...?
            self.thermo = {"dh": 0, "dg": 0, "tm": 0, "c50": 0}

    def search_by_seq(self, string, min_len):
        """ Find chunks containing either U, T or C
        :rtype : a list of tuples, the len(list) is the number
                 of occurrences of the motif (e.g. number of putative tfos)
        """
        # this here should be moved somewhere and reuse the regular_exp
        # as they do not change
        filter_regex = '[U|T|C]{%s,}' % (min_len)
        regular_exp = re.compile(filter_regex)
        b_in = []
        b_end = []
        for hit in regular_exp.finditer(string):
            b_in.append(hit.start())
            b_end.append(hit.end())
        return izip(b_in, b_end)

    def find_tfo(self, local_rec, index, b_rnafold=False,
                 b_very_verbose=False):
        """ Find the tfo regions based on sequence and
            If b_rnafold is set true, it will only report regions
            that are supposed to be unpaired according to RNAfold

            Defaults are bRnafold = False, b_very_verbose = False

        """
        # at least chunks of 12 consecutive C/U
        lim_tfo = 12
        found_tfo = LncRna.Tfos()
        tfos_lim = self.search_by_seq(str(local_rec.seq), lim_tfo)
        # !!!! lncipedia gives you DNA, so you must transcribe it, see below
        for b_init, b_end in tfos_lim:
            found_tfo.seq = str(local_rec.seq[b_init:b_end].transcribe())
            found_tfo.thermo = compute_prediction(
                found_tfo.seq, b_very_verbose)
            # only keep it if it's worth it, add the rest of the data
            if found_tfo.thermo["tm"] > 318.15:
                if b_rnafold:
                    if not self.struct_computed \
                            and len(local_rec.seq) < 1200:
                        self.struct = rna_struct_extern(
                            str(local_rec.seq.transcribe()))
                        self.struct_computed = True
                    if all(i in ".," for i in self.struct[
                           b_init:b_end]) and len(self.struct) > 0:
                        # This returns true even if self.struct is empty
                        # which is actually a good thing
                        self.found = True
                        found_tfo.position = b_init, b_end
                        found_tfo.index_id = index
                        found_tfo.id_lncrna = local_rec.id
                        self.lnctfos.append(found_tfo)
                    elif len(self.struct) > 0:
                        new_seq_list = [i if x == "." or x == "," else "X"
                                        for x, i in
                                        izip(self.struct[b_init:b_end],
                                             found_tfo.seq)]
                        new_seq = "".join(new_seq_list)
                        tfos_lim_sub = self.search_by_seq(new_seq, lim_tfo)
                        for b_init_s, b_end_s in tfos_lim_sub:
                            found_tfo.seq = new_seq[b_init_s:b_end_s]
                            found_tfo.thermo = compute_prediction(
                                found_tfo.seq, b_very_verbose)
                            if found_tfo.thermo["tm"] > 318:
                                self.found = True
                                found_tfo.position = b_init_s, b_end_s
                                found_tfo.index_id = index
                                found_tfo.id_lncrna = local_rec.id
                                self.lnctfos.append(found_tfo)
                else:
                    # print self.struct[b_init:b_end]
                    self.found = True
                    found_tfo.position = b_init, b_end
                    found_tfo.index_id = index
                    found_tfo.id_lncrna = local_rec.id
                    self.lnctfos.append(found_tfo)

    def __str__(self):
        return "The object that returns a list of TFOs in a given SeqIO record"


def rna_struct_extern(rna_seq):
    """ Compute the secondary structure using a call to RNAFold
    :rype: a string with the secondary structure, . or , means likely unpaired
    """
    cmd = "/usr/local/bin/RNAfold -p <<< %s " % (rna_seq)
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell=True)
    out, err = proc.communicate()
    if err is not 0 and err is None:
        print "Error calling RNAfold:"
        print err
        sys.exit()
    else:
        return out.split("\n")[2].split(" ")[0]


def check_coefficients(t_dh, t_dg, c_ph):
    """ Be carefull and check the values of the coefficients
        which should have been imported as constants
    """

    # argsparse should have done that, but it does not hurt
    if len(t_dh) != 3:
        print "I need exactly three values for DH prediction"
        print "The order of coefficients is CC, CT/TC, TT"
        sys.exit()

    if len(t_dg) != 4:
        print "I need exactly three values for DG prediction"
        print "The order of coefficients is C, T, CC"
        sys.exit()

    if len(c_ph) != 2:
        print "I need exactly two values for pH correction"
        print "The order of coefficients is C, T, CC"
        sys.exit()


def rev_hoog(sequence):
    """ Return the reverse hoogsten complements
    """
    my_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([my_complement.get(nt.upper(), '')
                    for nt in sequence[::-1]])


def overcount(string, sub):
    """ Count occurrences of sub in string
    """
    count, start = 0, 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


def nearcount_str(tfo, b_very_verbose, b_mixed):
    """This function returns a
        matrix of CC, TC/CT, TT composition
    """

    if b_mixed is False:
        A = np.array([0, 0, 0])
        if SPEC == 'DNA':
            A[0] = overcount(tfo, 'CC')
            A[1] = overcount(tfo, 'TC')
            A[1] += overcount(tfo, 'CT')
            A[2] = overcount(tfo, 'TT')
        elif SPEC == 'RNA':
            A[0] = overcount(tfo, 'CC')
            A[1] = overcount(tfo, 'UC')
            A[1] += overcount(tfo, 'CU')
            A[2] = overcount(tfo, 'UU')
        if b_very_verbose:
            print "The matrix with the dinucleotide step occupancies is "
            if SPEC == 'DNA':
                print ' {0:>4s} {1:>7s} {2:>4s} ' .format('CC', 'TC/CT', 'TT')
            else:
                print ' {0:>4s} {1:>7s} {2:>4s} ' .format('CC', 'UC/CU', 'UU')
            print ' {0:4d} {1:7d} {2:4d} ' .format(int(A[0]),
                                                   int(A[1]), int(A[2]))
    else:
        A = np.array([0, 0, 0, 0])
        if SPEC == 'DNA':
            A[0] = overcount(tfo, 'C')
            A[1] = overcount(tfo, 'T')
            A[2] = overcount(tfo, 'CC')
            A[3] = 1
        elif SPEC == 'RNA':
            A[0] = overcount(tfo, 'C')
            A[1] = overcount(tfo, 'U')
            A[2] = overcount(tfo, 'CC')
            A[3] = 1
        if b_very_verbose:
            print "The matrix with the dinucleotide step occupancies is "
            if SPEC == 'DNA':
                print ' {0:>4s} {1:>7s} {2:>4s} ' .format('C', 'T', 'CC')
            else:
                print ' {0:>4s} {1:>7s} {2:>4s} ' .format('C', 'U', 'CC')
            print ' {0:4d} {1:7d} {2:4d} ' .format(int(A[0]),
                                                   int(A[1]), int(A[2]))
    return A


def apply_ph_correct(om_dg, p_dg):
    """ Apply corrections to DG based on the number of C and CC
        The corrections formula is C*(pH-5.6)*(a+b*CC), where C and CC are the
        number of C and CC in the sequence, entries om_dg[0] and om_dg[1]

        :rtype: numpy matrix with corrected dg
    """

    rdg = p_dg + om_dg[0] * (VAL_PH - 5.6) * \
        (COEF_PH[0] + COEF_PH[1] * om_dg[1])

    return rdg


def compute_tm(pred_dh, pred_dg):
    """ Compute the Tm based on DG, DH and the concentration
        tm = (298*DH / ( DH-DG - 298*R*log(4/Ct)))
        the Ct (concentration) is the second record in cond_mat

    :rtype: numpy matrix with Tm estimates
    """
    ref_t = 298
    tm_trip = (ref_t * pred_dh /
               (pred_dh - pred_dg - 298 * R * np.log(4 / CONC_TFO / 1e-6)))

    return tm_trip


def compute_c50(p_dg):
    """ Compute concentration of TFO needed to achieve 50% of
        triplex formation. If I am not mistaken, that should go down
        like this, where names of species imply concentration:
        triplex@c50 = duplex_initial / 2  == d/2
        Oo = initial TFO conc
        d/2 / (d/2*Oo-d/2) = K_eq
        Solving for Oo we have
        Oo = (1/K_eq) + d/2 = exp(DG/kT) + d/2
    """

    ref_t = 310
    # we use uM units, hence the 1e6 factor
    oligo_c50 = 1e6 * np.exp(p_dg / (ref_t * R)) + 0.5 * DUPLEX_CONC

    return oligo_c50


def compute_prediction(sequence, b_very_verbose=False):
    """ Given a sequence returna dict with dh, dg, tm, c50
    :rtype: A dict with the results
    """

    b_mixed = False
    om_dh = nearcount_str(sequence, b_very_verbose, b_mixed)
    b_mixed = True
    om_dg = nearcount_str(sequence, b_very_verbose, b_mixed)
    b_mixed = False
    t_dh = np.dot(om_dh, NN_DH)
    t_dg = np.dot(om_dg, NN_DG)
    t_dg = apply_ph_correct(om_dg, t_dg)
    t_tm = compute_tm(t_dh, t_dg)
    t_c50 = compute_c50(t_dg)
    results = {"dh": t_dh, "dg": t_dg, "tm": t_tm, "c50": t_c50}
    return results


def print_pred_short(pred_dh, pred_dg, pred_tm, c50):
    """ Print out the predictions of the model

    """
    print 'Prediction DH: {0:4.2f}'. format(float(pred_dh))
    print 'Prediction DG: {0:4.2f}'. format(float(pred_dg))
    print 'Prediction Tm: {0:4.2f}'. format(float(pred_tm))
    print 'Prediction Tm_deg: {0:4.2f}'. format(float(pred_tm) - 273.15)
    print 'C_50 at 37 deg {0:5.3f} duplex conc: {1:5.3f} uM '.\
        format(float(DUPLEX_CONC), float(c50))


def parse_arguments():
    """ Trying it pythonic
    :rtype: arguments dictionary
    """

    my_parser = argparse.ArgumentParser(description=ROLLO, epilog=WHOWHEN)
    # right now we only expect one fasta file as input
    # if you want you could change nargs to + and process more than one
    my_parser.add_argument('-f', metavar='tfo_file', type=str, nargs=1,
                           required=True,
                           help='A multisequence fasta input file')
    my_parser.add_argument('-o', metavar='output_file', type=str,
                           default="lncr_tfo.fna", nargs='?', required=False,
                           help='A multisequence fasta output file')
    my_parser.add_argument(
        '-fold',
        action='store_true',
        help='Filter based on secondary structure (e.g. unpaired region) ')
    my_parser.add_argument(
        '-v',
        action='store_true',
        help='Tell me what you read. ')
    my_parser.add_argument(
        '-vv',
        action='store_true',
        help='Be very verbose. ')

    args = vars(my_parser.parse_args())
    return args


def get_input_records(args):
    """ read in the arguments
    :rtype: the bioseq records
    """

    b_very_verbose = args['vv']
    # args['f'] contains a list of only one file at the moment
    seq_file = args['f'][0]
    # this is a str no matter what
    # assume a series of physiological CONDitions at the moment

    # Reading in the sequences
    handle = open(seq_file, "rU")
    if b_very_verbose:
        print 'Reading {0}' .format(seq_file)
    the_records = list(SeqIO.parse(handle, "fasta"))
    if b_very_verbose:
        print 'We have read', len(records)
    handle.close()
    return the_records


def write_output_records(args, p_lncrna_w_tfos):
    """ Write out stuff
    """
    out_file = args['o']

    total_found = 0
    out_records = []
    for lncrna in p_lncrna_w_tfos:
        for tfo in lncrna.lnctfos:
            args_id = [tfo.id_lncrna, str(round(tfo.thermo["dg"], 2)),
                       str(round(tfo.thermo["tm"] - 273.15, 2))]
            id_f = "|".join(val for val in args_id)
            seq_rec = Seq(tfo.seq, IUPAC.unambiguous_rna)
            back_transc = seq_rec.back_transcribe()
            complement = back_transc.complement()
            rec_sense = SeqRecord(complement, id=id_f, description="")
            # if it gets too large you should implement a different solution
            out_records.append(rec_sense)
            rec_antisense = SeqRecord(back_transc, id=id_f, description="")
            out_records.append(rec_antisense)
            total_found += 1

    SeqIO.write(out_records, out_file, "fasta")
    print "\nDone! Found a total of {0} possible TFOs\n".format(total_found)

if __name__ == "__main__":
    arguments = parse_arguments()
    records = get_input_records(arguments)
    lncrna_w_tfos = find_tfo_from_mfasta(records, b_rnafold=arguments['fold'],
                                         b_very_verbose=arguments['vv'])
    write_output_records(arguments, lncrna_w_tfos)

    sys.exit()
