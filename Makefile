# Run the below command once in every terminal
# export LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}

SEQCC=g++
MPICC=mpiCC
MPIRUN=mpirun -np 4
CFLAGS=-std=c++11 -Wall -O3
MPI_MAIN=main.cpp
SEQ_MAIN=ap_seq.cpp

ARGS-1 = \
	--taskid=1 \
	--verbose=1 \
	--startk=0 \
	--endk=5 \
	--inputpath=./test-1/test-input--1.gra \
	--headerpath=./test-1/test-header--1.dat \
	--outputpath=./test-1/test-output--1.txt \

ARGS-2 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=5 \
	--inputpath=./test-2/test-input--2.gra \
	--headerpath=./test-2/test-header--2.dat \
	--outputpath=./test-2/test-output--2.txt \

ARGS0 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=10 \
	--inputpath=./test0/test-input-0.gra \
	--headerpath=./test0/test-header-0.dat \
	--outputpath=./test0/test-output-0.txt \

ARGS1 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=10 \
	--inputpath=./test1/test-input-1.gra \
	--headerpath=./test1/test-header-1.dat \
	--outputpath=./test1/test-output-1.txt \

ARGS2 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=8 \
	--inputpath=./test2/test-input-2.gra \
	--headerpath=./test2/test-header-2.dat \
	--outputpath=./test2/test-output-2.txt \

ARGS3 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=5 \
	--inputpath=./test3/test-input-3.gra \
	--headerpath=./test3/test-header-3.dat \
	--outputpath=./test3/test-output-3.txt \

ARGS4 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=6 \
	--inputpath=./test4/test-input-4.gra \
	--headerpath=./test4/test-header-4.dat \
	--outputpath=./test4/test-output-4.txt \

ARGS5 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=8 \
	--inputpath=./test5/test-input-5.gra \
	--headerpath=./test5/test-header-5.dat \
	--outputpath=./test5/test-output-5.txt \

ARGS6 = \
	--taskid=1 \
	--verbose=1 \
	--startk=1 \
	--endk=29 \
	--inputpath=./test6/test-input-6.gra \
	--headerpath=./test6/test-header-6.dat \
	--outputpath=./test6/test-output-6.txt \

ARGS7 = \
	--taskid=1 \
	--verbose=1 \
	--startk=10 \
	--endk=25 \
	--inputpath=./test7/test-input-7.gra \
	--headerpath=./test7/test-header-7.dat \
	--outputpath=./test7/test-output-7.txt \

ARGS8 = \
	--taskid=1 \
	--verbose=1 \
	--startk=2 \
	--endk=6 \
	--inputpath=./test8/test-input-8.gra \
	--headerpath=./test8/test-header-8.dat \
	--outputpath=./test8/test-output-8.txt \

ARGS9 = \
	--taskid=1 \
	--verbose=1 \
	--startk=3 \
	--endk=3 \
	--inputpath=./test9/test-input-9.gra \
	--headerpath=./test9/test-header-9.dat \
	--outputpath=./test9/test-output-9.txt \

ARGS10 = \
	--taskid=1 \
	--verbose=1 \
	--startk=3 \
	--endk=3 \
	--inputpath=./test10/test-input-10.gra \
	--headerpath=./test10/test-header-10.dat \
	--outputpath=./test10/test-output-10.txt \

#####################################################

mpi_exec: mpi_main.o
	$(MPICC) $(CFLAGS) mpi_main.o -o a3

mpi_main.o: $(MPI_MAIN)
	$(MPICC) $(CFLAGS) -c $(MPI_MAIN) -o mpi_main.o

mpi_run-2: mpi_exec
	$(MPIRUN) a3 $(ARGS-2)
mpi_run-1: mpi_exec
	$(MPIRUN) a3 $(ARGS-1)
mpi_run0: mpi_exec
	$(MPIRUN) a3 $(ARGS0)
mpi_run1: mpi_exec
	$(MPIRUN) a3 $(ARGS1)
mpi_run2: mpi_exec
	$(MPIRUN) a3 $(ARGS2)
mpi_run3: mpi_exec
	$(MPIRUN) a3 $(ARGS3)
mpi_run4: mpi_exec
	$(MPIRUN) a3 $(ARGS4)
mpi_run5: mpi_exec
	$(MPIRUN) a3 $(ARGS5)
mpi_run6: mpi_exec
	$(MPIRUN) a3 $(ARGS6)
mpi_run7: mpi_exec
	$(MPIRUN) a3 $(ARGS7)
mpi_run8: mpi_exec
	$(MPIRUN) a3 $(ARGS8)
mpi_run9: mpi_exec
	$(MPIRUN) a3 $(ARGS9)
mpi_run10: mpi_exec
	$(MPIRUN) a3 $(ARGS10)

mpi_run: mpi_run-2 mpi_run-1 mpi_run0 mpi_run1 mpi_run2 mpi_run3 mpi_run4 mpi_run5 mpi_run6 mpi_run7 mpi_run8 mpi_run9 mpi_run10

#####################################################

seq_exec: seq_main.o
	$(SEQCC) $(CFLAGS) seq_main.o -o a3_seq

seq_main.o: $(SEQ_MAIN)
	$(SEQCC) $(CFLAGS) -c $(SEQ_MAIN) -o seq_main.o

seq_run-2: seq_exec
	./a3_seq $(ARGS-2)
seq_run-1: seq_exec
	./a3_seq $(ARGS-1)
seq_run0: seq_exec
	./a3_seq $(ARGS0)
seq_run1: seq_exec
	./a3_seq $(ARGS1)
seq_run2: seq_exec
	./a3_seq $(ARGS2)
seq_run3: seq_exec
	./a3_seq $(ARGS3)
seq_run4: seq_exec
	./a3_seq $(ARGS4)
seq_run5: seq_exec
	./a3_seq $(ARGS5)
seq_run6: seq_exec
	./a3_seq $(ARGS6)
seq_run7: seq_exec
	./a3_seq $(ARGS7)
seq_run8: seq_exec
	./a3_seq $(ARGS8)
seq_run9: seq_exec
	./a3_seq $(ARGS9)
seq_run10: seq_exec
	./a3_seq $(ARGS10)

seq_run: seq_run-2 seq_run-1 seq_run0 seq_run1 seq_run2 seq_run3 seq_run4 seq_run5 seq_run6 seq_run7 seq_run8 seq_run9 seq_run10

#####################################################

verify.o: verify.cpp
	$(SEQCC) $(CFLAGS) verify.cpp -o verify.o

check0: verify.o
	./verify.o test0/test-output-0.txt test0/output0_verbose.txt
check1: verify.o
	./verify.o test1/test-output-1.txt test1/output1_verbose.txt
check2: verify.o
	./verify.o test2/test-output-2.txt test2/output2_verbose.txt
check3: verify.o
	./verify.o test3/test-output-3.txt test3/output3_verbose.txt
check4: verify.o
	./verify.o test4/test-output-4.txt test4/output4_verbose.txt
check5: verify.o
	./verify.o test5/test-output-5.txt test5/output5_verbose.txt
check6: verify.o
	./verify.o test6/test-output-6.txt test6/output6_verbose.txt
check7: verify.o
	./verify.o test7/test-output-7.txt test7/output7_verbose.txt
check8: verify.o
	./verify.o test8/test-output-8.txt test8/output8_verbose.txt
check9: verify.o
	./verify.o test9/test-output-9.txt test9/output9_verbose.txt
check10: verify.o
	./verify.o test10/test-output-10.txt test10/output10_verbose.txt

check: check0 check1 check2 check3 check4 check5 check6 check7 check8 check9 check10

#####################################################

clean:
	rm a3* *.o
