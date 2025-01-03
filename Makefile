# Run the below command once in every terminal
# export LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}

SEQCC=g++
MPICC=mpiCC
MPIRUN=mpirun -np 3
CFLAGS=-std=c++17 -Wall -O3 -fopenmp
MPI_MAIN=main.cpp
SEQ_MAIN=ap_seq.cpp

ARGS-1 = \
	--taskid=1 \
	--verbose=0 \
	--startk=0 \
	--endk=5 \
	--inputpath=A2/test-1/test-input--1.gra \
	--headerpath=A2/test-1/test-header--1.dat \
	--outputpath=A2/test-1/test-output--1.txt \

ARGS-2 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=5 \
	--inputpath=A2/test-2/test-input--2.gra \
	--headerpath=A2/test-2/test-header--2.dat \
	--outputpath=A2/test-2/test-output--2.txt \

ARGS0 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=10 \
	--inputpath=A2/test0/test-input-0.gra \
	--headerpath=A2/test0/test-header-0.dat \
	--outputpath=A2/test0/test-output-0.txt \

ARGS1 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=10 \
	--inputpath=A2/test1/test-input-1.gra \
	--headerpath=A2/test1/test-header-1.dat \
	--outputpath=A2/test1/test-output-1.txt \

ARGS2 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=8 \
	--inputpath=A2/test2/test-input-2.gra \
	--headerpath=A2/test2/test-header-2.dat \
	--outputpath=A2/test2/test-output-2.txt \

ARGS3 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=5 \
	--inputpath=A2/test3/test-input-3.gra \
	--headerpath=A2/test3/test-header-3.dat \
	--outputpath=A2/test3/test-output-3.txt \

ARGS4 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=6 \
	--inputpath=A2/test4/test-input-4.gra \
	--headerpath=A2/test4/test-header-4.dat \
	--outputpath=A2/test4/test-output-4.txt \

ARGS5 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=8 \
	--inputpath=A2/test5/test-input-5.gra \
	--headerpath=A2/test5/test-header-5.dat \
	--outputpath=A2/test5/test-output-5.txt \

ARGS6 = \
	--taskid=1 \
	--verbose=0 \
	--startk=1 \
	--endk=29 \
	--inputpath=A2/test6/test-input-6.gra \
	--headerpath=A2/test6/test-header-6.dat \
	--outputpath=A2/test6/test-output-6.txt \

ARGS7 = \
	--taskid=1 \
	--verbose=0 \
	--startk=10 \
	--endk=25 \
	--inputpath=A2/test7/test-input-7.gra \
	--headerpath=A2/test7/test-header-7.dat \
	--outputpath=A2/test7/test-output-7.txt \

ARGS8 = \
	--taskid=1 \
	--verbose=0 \
	--startk=2 \
	--endk=6 \
	--inputpath=A2/test8/test-input-8.gra \
	--headerpath=A2/test8/test-header-8.dat \
	--outputpath=A2/test8/test-output-8.txt \

ARGS9 = \
	--taskid=1 \
	--verbose=0 \
	--startk=3 \
	--endk=3 \
	--inputpath=A2/test9/test-input-9.gra \
	--headerpath=A2/test9/test-header-9.dat \
	--outputpath=A2/test9/test-output-9.txt \

ARGS10 = \
	--taskid=1 \
	--verbose=0 \
	--startk=3 \
	--endk=3 \
	--inputpath=A2/test10/test-input-10.gra \
	--headerpath=A2/test10/test-header-10.dat \
	--outputpath=A2/test10/test-output-10.txt \

# A3: TASK 1 Test Cases

ARGS311 = \
	--taskid=1 \
	--verbose=0 \
	--startk=3 \
	--endk=10 \
	--inputpath=A3/test1/test-input-1.gra \
	--headerpath=A3/test1/test-header-1.dat \
	--outputpath=A3/test1/task1-output-1.txt \

ARGS312 = \
	--taskid=1 \
	--verbose=0 \
	--startk=3 \
	--endk=12 \
	--inputpath=A3/test2/test-input-2.gra \
	--headerpath=A3/test2/test-header-2.dat \
	--outputpath=A3/test2/task1-output-2.txt \

ARGS313 = \
	--taskid=1 \
	--verbose=0 \
	--startk=4 \
	--endk=8 \
	--inputpath=A3/test3/test-input-3.gra \
	--headerpath=A3/test3/test-header-3.dat \
	--outputpath=A3/test3/task1-output-3.txt \

ARGS314 = \
	--taskid=1 \
	--verbose=0 \
	--startk=2 \
	--endk=3 \
	--inputpath=A3/test4/test-input-4.gra \
	--headerpath=A3/test4/test-header-4.dat \
	--outputpath=A3/test4/task1-output-4.txt \

# A3: TASK 2 Test Cases

ARGS321 = \
	--taskid=2 \
	--verbose=0 \
	--endk=3 \
	--inputpath=A3/test1/test-input-1.gra \
	--headerpath=A3/test1/test-header-1.dat \
	--outputpath=A3/test1/task2-output-1.txt \
	--p=10 \

ARGS322 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=3 \
	--inputpath=A3/test2/test-input-2.gra \
	--headerpath=A3/test2/test-header-2.dat \
	--outputpath=A3/test2/task2-output-2.txt \
	--p=10 \

ARGS323 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=3 \
	--inputpath=A3/test3/test-input-3.gra \
	--headerpath=A3/test3/test-header-3.dat \
	--outputpath=A3/test3/task2-output-3.txt \
	--p=4 \

ARGS324 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=2 \
	--inputpath=A3/test4/test-input-4.gra \
	--headerpath=A3/test4/test-header-4.dat \
	--outputpath=A3/test4/task2-output-4.txt \
	--p=20 \

ARGS326 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=4 \
	--inputpath=A3/test6/test-input-6.gra \
	--headerpath=A3/test6/test-header-6.dat \
	--outputpath=A3/test6/task2-output-6.txt \
	--p=6 \

ARGS327 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=3 \
	--inputpath=A3/test7/test-input-7.gra \
	--headerpath=A3/test7/test-header-7.dat \
	--outputpath=A3/test7/task2-output-7.txt \
	--p=4 \

ARGS328 = \
	--taskid=2 \
	--verbose=0 \
	--startk=1 \
	--endk=5 \
	--inputpath=A3/test8/test-input-8.gra \
	--headerpath=A3/test8/test-header-8.dat \
	--outputpath=A3/test8/task2-output-8.txt \
	--p=6 \

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

mpi_run311: mpi_exec
	$(MPIRUN) a3 $(ARGS311)
mpi_run312: mpi_exec
	$(MPIRUN) a3 $(ARGS312)
mpi_run313: mpi_exec
	$(MPIRUN) a3 $(ARGS313)
mpi_run314: mpi_exec
	$(MPIRUN) a3 $(ARGS314)

mpi_run321: mpi_exec
	$(MPIRUN) a3 $(ARGS321)
mpi_run322: mpi_exec
	$(MPIRUN) a3 $(ARGS322)
mpi_run323: mpi_exec
	$(MPIRUN) a3 $(ARGS323)
mpi_run324: mpi_exec
	$(MPIRUN) a3 $(ARGS324)
mpi_run326: mpi_exec
	$(MPIRUN) a3 $(ARGS326)
mpi_run327: mpi_exec
	$(MPIRUN) a3 $(ARGS327)
mpi_run328: mpi_exec
	$(MPIRUN) a3 $(ARGS328)


mpi_run21: mpi_run-2 mpi_run-1 mpi_run0 mpi_run1 mpi_run2 mpi_run3 mpi_run4 mpi_run5
# mpi_run6 mpi_run7 mpi_run8 mpi_run9 mpi_run10
mpi_run31: mpi_run311 mpi_run312 mpi_run313 mpi_run314
mpi_run32: mpi_run321 mpi_run322 mpi_run323 mpi_run324
# mpi_run326 mpi_run327 mpi_run328

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
	./verify.o A2/test0/test-output-0.txt A2/test0/output0_verbose.txt
check1: verify.o
	./verify.o A2/test1/test-output-1.txt A2/test1/output1_verbose.txt
check2: verify.o
	./verify.o A2/test2/test-output-2.txt A2/test2/output2_verbose.txt
check3: verify.o
	./verify.o A2/test3/test-output-3.txt A2/test3/output3_verbose.txt
check4: verify.o
	./verify.o A2/test4/test-output-4.txt A2/test4/output4_verbose.txt
check5: verify.o
	./verify.o A2/test5/test-output-5.txt A2/test5/output5_verbose.txt
check6: verify.o
	./verify.o A2/test6/test-output-6.txt A2/test6/output6.txt
check7: verify.o
	./verify.o A2/test7/test-output-7.txt A2/test7/output7.txt
check8: verify.o
	./verify.o A2/test8/test-output-8.txt A2/test8/output8.txt
check9: verify.o
	./verify.o A2/test9/test-output-9.txt A2/test9/output9.txt
check10: verify.o
	./verify.o A2/test10/test-output-10.txt A2/test10/output10.txt

check311: verify.o
	./verify.o A3/test1/task1-output-1.txt A3/test1/task1_output1_verbose.txt
check312: verify.o
	./verify.o A3/test2/task1-output-2.txt A3/test2/task1_output2_verbose.txt
check313: verify.o
	./verify.o A3/test3/task1-output-3.txt A3/test3/task1_output3_verbose.txt
check314: verify.o
	./verify.o A3/test4/task1-output-4.txt A3/test4/task1_output4_verbose.txt

check321: verify.o
	./verify.o A3/test1/task2-output-1.txt A3/test1/task2_output1_verbose.txt
check322: verify.o
	./verify.o A3/test2/task2-output-2.txt A3/test2/task2_output2_verbose.txt
check323: verify.o
	./verify.o A3/test3/task2-output-3.txt A3/test3/task2_output3_verbose.txt
check324: verify.o
	./verify.o A3/test4/task2-output-4.txt A3/test4/task2_output4_verbose.txt
check326: verify.o
	./verify.o A3/test6/task2-output-6.txt A3/test6/task2_output6.txt
check327: verify.o
	./verify.o A3/test7/task2-output-7.txt A3/test7/task2_output7.txt
check328: verify.o
	./verify.o A3/test8/task2-output-8.txt A3/test8/task2_output8.txt


check21: check0 check1 check2 check3 check4 check5
# check6 check7 check8 check9 check10
check31: check311 check312 check313 check314
check32: check321 check322 check323 check324
# check326 check327 check328

#####################################################

clean:
	rm a3* *.o
