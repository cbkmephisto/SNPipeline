# makefile for SNPipeline
# 2013-11-20 12:01 PM

GPP=g++
FLAGS=-c -pipe -O3 -g -funsafe-math-optimizations # -Wall -W
CC=$(GPP) $(FLAGS)
LK=$(GPP)
INC=-Ilibcbk
INCEIGEN=-I/opt/local/include/eigen3/

all:    bin/SNPipeline\
        bin/SNPipeline.mergingAB\
        bin/finalReportReformer\
        bin/SNPipeline.pooledAB2VCF41\
        bin/SNPipeline.VCF2AB\
        bin/SNPipeline.trimDown2Map\
        bin/SNPipeline.abgExcl\
		bin/gReplace\
		bin/ggrep\
        bin/SNPipeline.abg2bin


############ link
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
