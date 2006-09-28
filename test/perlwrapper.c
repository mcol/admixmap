/*
simple C program to provide a parallel frontend for a perlscript

compile with something like:
pathcc -operlwrapper -I/usr/local/mpich/path/include perlwrapper.c -L/usr/local/mpich/path/lib64 -lmpich -lgmp
*/
#include "mpi.h"
#include <stdio.h>
#include <string.h>

int main(int argc, char**argv){

    if(argc < 2){
	printf("%s\n", "Please specify a perl script to run");
	return 1;;
    }
    const char* perlscript = argv[1];

    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

/*construct command*/
    char command[20];
    sprintf(command, "%s %s %d %d", "perl", perlscript, rank, size);

    printf("\n%s\n", command);
/*run main script with rank and size as args*/
    system(command);

/*wait till everyone is finished*/
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0){
	/* pull results into report file   */
	if(argc >2){
	    char command2[20];
	    sprintf(command, "%s %s", "perl", argv[2]);
	    system(command2);
	}
	else{
	    /*do it from here*/
	    ;
	}
    }


    MPI_Finalize();
    return 0;
}
