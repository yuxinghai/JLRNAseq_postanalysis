from optparse import OptionParser
import os,sys


def IR_result_parser(name):
    '''
for parser IRfinder result to 
Chr+"\t"+Start+"\t"+End+"\t"ensemble_id+"\t"+Symbol
'''
    with open(name) as f:
        list1 = []
        for lines in f:
            if lines.startswith("#"):
                continue
	    elif lines.startswith("Chr"):
		continue
            else:
                Chr,Start,End,Intron_GeneName_GeneID,X,Direction,ExcludedBases,p_diff,p_increased,p_decreased,A_IRratio,A_IRok,A_IntronCover,A_IntronDepth, A_SplicesMax,A_SplicesExact, B_IRratio,B_IRok,BIntronCover,B_IntronDepth,B_SplicesMax,B_SplicesExact,replicates,A1_IRratio,A2_IRratio,B1_IRratio,B2_IRratio=lines.strip().split("\t")
                Symbol=Intron_GeneName_GeneID.split("/")[0]
		ensemble_id=Intron_GeneName_GeneID.split("/")[1]
                line=Chr+"\t"+Start+"\t"+End+"\t"+ensemble_id+"\t"+Symbol
                list1.append(line)

    return list1


'''
only for get gene start and end cordinary
Chr, Start, End, Strand,ensemble_id, Symbol
'''
def genename_parser(name):
    list2 = []
    with open(name) as gene:
        gene.readline()
        for lines in gene:
            lines=lines.strip().split()
            #Chr, Start, End, Strand,ensemble_id, Symbol = lines.strip().split()
            list2.append(lines)
    return list2




def main():
    '''
    get Chr, Start, End, Symbol,0,,Strand,IRstart,IRend

    '''

    usage = "%prog [options]"+"\n"
    parser = OptionParser(usage, version="%prog ")
    parser.add_option("-i", "--input-file", action="store", type="string", dest="input_file",
                      help="input file gene name")
    parser.add_option("-r", "--refgene", action="store", type="string", dest="ref_gene_model",
                      help="Reference gene model in bed format.")
    parser.add_option("-o", "--output", action="store", type="string", dest="output_file",
                      help="output file gtf")
    (options, args) = parser.parse_args()

    if not (options.input_file and options.ref_gene_model):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(options.ref_gene_model):
        print >> sys.stderr, '\n\n' + options.ref_gene_model + " does NOT exists" + '\n'
        sys.exit(0)
    if not os.path.exists(options.input_file):
        print >> sys.stderr, '\n\n' + options.input_file + " does NOT exists" + '\n'
        sys.exit(0)


    list3 = []
    for lines in IR_result_parser(options.input_file):
        lines=lines.split("\t")
        for ano in genename_parser(options.ref_gene_model):
	    #print lines,ano
	    #exit(0)
            if (lines[3] == ano[4]) and (lines[4] == ano[5]):
		#gene symbol and ensemble id equal
                new_line = "\t".join([ano[0], ano[1], ano[2], ano[5], "0", ano[3], lines[1], lines[2]])
		#Chr, Start, End, Symbol,0,Strand,IRstart,IRend
                list3.append(new_line)


    with open(options.output_file, "w") as output:
        for line in list3:
            Chr, Start, End, Symbol, X, Strand, ThickStart, Thickend = line.split("\t")
            Start=int(Start)
            End=int(End)
            ThickStart=int(ThickStart)
            Thickend=int(Thickend)
            if Strand=="+":
                rela_start = ThickStart - Start
                rela_end = Thickend - Start
            else:
                rela_start = End - Thickend
                rela_end = End - ThickStart

	    #divide whole gene to 100 part
            pos_s = (100.0 / (End - Start)) * rela_start
            pos_s=int(round(pos_s))
            pos_e = (100.0 / (End - Start)) * rela_end
            pos_e=int(round(pos_e))
            if pos_e==pos_s:
	    #prevent some gene pos_s equal to pos_e
                pos_e+=1

            new = Chr + "\t" + str(pos_s) + "\t" + str(pos_e)+"\t" + Symbol + "\t" + Strand+"\n"
            output.write(new)
	    
if __name__ == '__main__':
    main()





















