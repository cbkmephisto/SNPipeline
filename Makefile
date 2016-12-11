# makefile for SNPipeline
# 2013-11-20 12:01 PM

GPP=g++
FLAGS=-c -pipe -O3 -g -funsafe-math-optimizations # -Wall -W
CC=$(GPP) $(FLAGS)
LK=$(GPP)
# TODO: change path to eigen3 header files folder, needed by SNPipeline.abg2bin
# for Ubuntu users on 16.10, before launching 'make', try
# sudo apt-get install libeigen3-dev
INCEIGEN=-I/opt/local/include/eigen3/ -I/usr/include/eigen3/

all:	bin/SNPipeline\
	bin/SNPipeline.mergingAB\
	bin/finalReportReformer\
	bin/SNPipeline.pooledAB2VCF41\
	bin/SNPipeline.VCF2AB\
	bin/SNPipeline.trimDown2Map\
	bin/SNPipeline.abgExcl\
	bin/gReplace\
	bin/ggrep\
	bin/abg2FImpute\
	bin/mapUniter\
	bin/abg2findhap\
	bin/fhout2fiout\
	bin/fhout2hapview\
	bin/fout2abg\
	bin/fout2haplotype\
	bin/abg2M\
	bin/bout2genotype\
	bin/bout2haplotype\
	bin/vcf2haplotype\
	bin/SNPipeline.abg2bin


############ link
bin/vcf2haplotype: obj/vcf2haplotype.o obj/logger.o obj/VCF_Record.o
	$(LK) $(INC) obj/vcf2haplotype.o obj/logger.o obj/VCF_Record.o -o bin/vcf2haplotype

bin/bout2haplotype: obj/bout2haplotype.o obj/logger.o obj/geneticMap.o obj/beagleOutput.o
	$(LK) $(INC) obj/bout2haplotype.o obj/logger.o obj/geneticMap.o obj/beagleOutput.o -o bin/bout2haplotype

bin/bout2genotype: obj/bout2genotype.o obj/logger.o obj/geneticMap.o obj/beagleOutput.o
	$(LK) $(INC) obj/bout2genotype.o obj/logger.o obj/geneticMap.o obj/beagleOutput.o -o bin/bout2genotype

bin/abg2M: obj/abg2M.o
	$(LK) $(INC) obj/abg2M.o -o bin/abg2M

bin/fout2haplotype: obj/fout2haplotype.o
	$(LK) $(INC) obj/fout2haplotype.o obj/logger.o -o bin/fout2haplotype

bin/fout2abg: obj/fout2abg.o
	$(LK) $(INC) obj/fout2abg.o obj/logger.o -o bin/fout2abg

bin/abg2findhap: obj/abg2findhap.o
	$(LK) $(INC) obj/abg2findhap.o obj/logger.o -o bin/abg2findhap

bin/fhout2fiout: obj/fhout2fiout.o
	$(LK) $(INC) obj/fhout2fiout.o obj/logger.o -o bin/fhout2fiout

bin/fhout2hapview: obj/fhout2hapview.o
	$(LK) $(INC) obj/fhout2hapview.o obj/logger.o -o bin/fhout2hapview

bin/mapUniter: obj/mapUniter.o
	$(LK) $(INC) obj/mapUniter.o obj/logger.o -o bin/mapUniter

bin/abg2FImpute: obj/abg2FImpute.o
	$(LK) $(INC) obj/abg2FImpute.o obj/logger.o -o bin/abg2FImpute

bin/gReplace: obj/gReplace.o
	$(LK) $(INC) obj/gReplace.o -o bin/gReplace

bin/ggrep: obj/ggrep.o
	$(LK) $(INC) obj/ggrep.o -o bin/ggrep

bin/SNPipeline.abgExcl: obj/SNPipeline.abgExcl.o
	$(LK) $(INC) obj/SNPipeline.abgExcl.o -o bin/SNPipeline.abgExcl

bin/SNPipeline.abg2bin: obj/SNPipeline.abg2bin.o obj/logger.o
	$(LK) $(INC) obj/SNPipeline.abg2bin.o obj/logger.o -o bin/SNPipeline.abg2bin

bin/SNPipeline.VCF2AB: obj/SNPipeline.VCF2AB.o obj/logger.o
	$(LK) $(INC) obj/SNPipeline.VCF2AB.o obj/logger.o -o bin/SNPipeline.VCF2AB

bin/SNPipeline.trimDown2Map: obj/SNPipeline.trimDown2Map.o obj/logger.o
	$(LK) $(INC) obj/logger.o obj/SNPipeline.trimDown2Map.o -o bin/SNPipeline.trimDown2Map

bin/SNPipeline.pooledAB2VCF41: obj/SNPipeline.pooledAB2VCF41.o obj/logger.o
	$(LK) $(INC) obj/logger.o obj/SNPipeline.pooledAB2VCF41.o -o bin/SNPipeline.pooledAB2VCF41

bin/finalReportReformer: obj/finalReportReformer.o obj/logger.o
	$(LK) $(INC) obj/logger.o obj/finalReportReformer.o -o bin/finalReportReformer

bin/SNPipeline: obj/SNPipeline.o obj/logger.o obj/data.o
	$(LK) $(INC) obj/logger.o obj/SNPipeline.o obj/data.o -o bin/SNPipeline

bin/SNPipeline.mergingAB: obj/SNPipeline.mergingAB.o obj/logger.o
	$(LK) $(INC) obj/logger.o obj/SNPipeline.mergingAB.o -o bin/SNPipeline.mergingAB

############# obj
obj/VCF_Record.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/VCF_Record.cpp -o obj/VCF_Record.o

obj/vcf2haplotype.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/vcf2haplotype.cpp -o obj/vcf2haplotype.o

obj/bout2haplotype.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/bout2haplotype.cpp -o obj/bout2haplotype.o

obj/geneticMap.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/geneticMap.cpp -o obj/geneticMap.o

obj/beagleOutput.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/beagleOutput.cpp -o obj/beagleOutput.o

obj/bout2genotype.o:
	$(CC) $(CFLAGS) $(INC) bout2genotype/bout2genotype.cpp -o obj/bout2genotype.o

obj/abg2M.o: abg2M/abg2M.cpp
	$(CC) $(CFLAGS) $(INC) abg2M/abg2M.cpp -o obj/abg2M.o

obj/fout2haplotype.o: fout2haplotype/fout2haplotype.cpp
	$(CC) $(CFLAGS) $(INC) fout2haplotype/fout2haplotype.cpp -o obj/fout2haplotype.o

obj/fout2abg.o: fout2abg/fout2abg.cpp
	$(CC) $(CFLAGS) $(INC) fout2abg/fout2abg.cpp -o obj/fout2abg.o

obj/abg2findhap.o: abg2findhap/abg2findhap.cpp
	$(CC) $(CFLAGS) $(INC) abg2findhap/abg2findhap.cpp -o obj/abg2findhap.o

obj/fhout2fiout.o: abg2findhap/fhout2fiout.cpp
	$(CC) $(CFLAGS) $(INC) abg2findhap/fhout2fiout.cpp -o obj/fhout2fiout.o

obj/fhout2hapview.o: abg2findhap/fhout2hapview.cpp
	$(CC) $(CFLAGS) $(INC) abg2findhap/fhout2hapview.cpp -o obj/fhout2hapview.o

obj/mapUniter.o: mapUniter/mapUniter.cpp
	$(CC) $(CFLAGS) $(INC) mapUniter/mapUniter.cpp -o obj/mapUniter.o

obj/abg2FImpute.o: abg2FImpute/abg2FImpute.cpp
	$(CC) $(CFLAGS) $(INC) abg2FImpute/abg2FImpute.cpp -o obj/abg2FImpute.o

obj/gReplace.o: gReplace/gReplace.cpp
	$(CC) $(CFLAGS) $(INC) gReplace/gReplace.cpp -o obj/gReplace.o

obj/ggrep.o: ggrep/ggrep.cpp
	$(CC) $(CFLAGS) $(INC) ggrep/ggrep.cpp -o obj/ggrep.o

obj/SNPipeline.abgExcl.o: SNPipeline.abgExcl/SNPipeline.abgExcl.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline.abgExcl/SNPipeline.abgExcl.cpp -o obj/SNPipeline.abgExcl.o

obj/SNPipeline.abg2bin.o: SNPipeline.abg2bin/SNPipeline.abg2bin.cpp
	$(CC) $(CFLAGS) $(INC) $(INCEIGEN) SNPipeline.abg2bin/SNPipeline.abg2bin.cpp -o obj/SNPipeline.abg2bin.o

obj/SNPipeline.VCF2AB.o: SNPipeline.VCF2AB/SNPipeline_VCF2AB.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline.VCF2AB/SNPipeline_VCF2AB.cpp -o obj/SNPipeline.VCF2AB.o

obj/SNPipeline.trimDown2Map.o: SNPipeline.trimDown2Map/main.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline.trimDown2Map/main.cpp -o obj/SNPipeline.trimDown2Map.o

obj/SNPipeline.pooledAB2VCF41.o: SNPipeline_pooledAB2VCF41.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline_pooledAB2VCF41.cpp -o obj/SNPipeline.pooledAB2VCF41.o

obj/finalReportReformer.o: finalReportReformer/main.cpp
	$(CC) $(CFLAGS) $(INC) finalReportReformer/main.cpp -o obj/finalReportReformer.o

obj/SNPipeline.mergingAB.o: SNPipeline.mergingAB.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline.mergingAB.cpp -o obj/SNPipeline.mergingAB.o

obj/SNPipeline.o: SNPipeline.cpp
	$(CC) $(CFLAGS) $(INC) SNPipeline.cpp -o obj/SNPipeline.o

#obj/finalReport.o: finalReport.cpp
#	$(CC) $(CFLAGS) $(INC) finalReport.cpp -o obj/finalReport.o

#obj/affyInput.o: affyInput.cpp
#	$(CC) $(CFLAGS) $(INC) affyInput.cpp -o obj/affyInput.o

obj/logger.o: logger.cpp
	$(CC) $(CFLAGS) $(INC) logger.cpp -o obj/logger.o

obj/data.o: data.cpp
	$(CC) $(CFLAGS) $(INC) data.cpp -o obj/data.o

###############

clean:
	rm -rf obj/*o bin/*
