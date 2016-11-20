
# TFO searching 

TFO finder
Searching putative triplex binding sites on DNA for the formation of RNA-DNA2 triplex.

Find triplex forming oligonucleotieds (TFOs) in a given (set of) sequence(s)
The sequences should be introduced using the fasta format
There are a couple of modes, mainly.

a) single stranded --> gives back the sequences of the binding partner of
                       the TFOs it finds. It gives back both the sense
                       and antisense. Works for DNA and RNA, but using
                       always the same set of rules. You can also request
                       using RNAfold to fold the single stranded, this only
                       makes sense for RNA. Futhermore, if you suply
                       constraints in the RNAfold format (modified fastq)
                       these can also be used. Only TFOs found in an unpaired
                       region will be reported.

b) double stranded --> gives back the TFO sequences that would bind to
                       the double stranded DNA (e.g. genome). It searches
                       both the fwd and rev strands. The output fasta file
                       containts the coordinates of the motif wrt the
                       + strand in the header. E.g., if it reports |1:10|-|
                       this means that the binding site for the TFO is
                       in the reverse strand, and its 5' starts at 10 and
                       the 3' finishes at 1.

TODO --> better document the output format

The code is organised such that you have the command line parsing
in one hpp file, the routines to compute the predicted stabilities of TFOs
in another one, and main class that does the job. The class is called
LncRna because that was the original intention, but the code grew to
treat double stranded DNA (tribially). Each object of the lncrna class
has one or severall Tfo objects, which store the TFO for each input
sequence. find_tfo.hpp, where the code is, should be well documented.
The class makes heavy usage of Tags to differentiate between function
calls that do the same but speciallized (eg. folding vs non-folding).

Quite important: do use -O3 for production, makes the codes way faster
although the compiler takes a tad longer

The code works fine, but there might be some hacks that require polishing
here and there.
