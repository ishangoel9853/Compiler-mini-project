//***********************************DISTRIBUTED K-MEANS CLUSTERING (USING MPI)**************************************

/* [numClusters]: no. objects assigned in each new cluster */
/*delta :  % of objects change their clusters */
/* **clusters: out: [numClusters][dimensions] */
/* **newClusters: [numClusters][dimensions] */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#define MAX_CHAR_PER_LINE 128

double wtime(void) {
    double          now_time;
    struct timeval  etstart;
    struct timezone tzp;

    if (gettimeofday(&etstart, &tzp) == -1)
        perror("Error: calling gettimeofday() not successful.\n");

    now_time = ((double)etstart.tv_sec) +
               ((double)etstart.tv_usec) / 1000000.0;
    return now_time;
}

float** file_read(char *filename, int  *numObjs, int  *dimensions) {
    float **objects;
    int     i, j, len;
    ssize_t numBytesRead;

    FILE *infile;
    char *line, *ret;
    int   lineLen;

    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: no such file (%s)\n", filename);
        return NULL;
    }

    /* finding the number of objects */
    lineLen = MAX_CHAR_PER_LINE;
    line = (char*) malloc(lineLen);
    assert(line != NULL);

    (*numObjs) = 0;
    while (fgets(line, lineLen, infile) != NULL) {
        while (strlen(line) == lineLen-1) {
            len = strlen(line);
            fseek(infile, -len, SEEK_CUR);

            lineLen += MAX_CHAR_PER_LINE;    /* increase lineLen */
            line = (char*) realloc(line, lineLen);
            assert(line != NULL);

            ret = fgets(line, lineLen, infile);
            assert(ret != NULL);
        }

        if (strtok(line, " \t\n") != 0)
            (*numObjs)++;
    }
    rewind(infile);


    /* find the no. objects of each object */
    (*dimensions) = 0;
    while (fgets(line, lineLen, infile) != NULL) {
        if (strtok(line, " \t\n") != 0) {
            while (strtok(NULL, " ,\t\n") != NULL) (*dimensions)++;
            break;
        }
    }
    rewind(infile);
    printf("File %s numObjs   = %d\n",filename,*numObjs);
    printf("File %s dimensions = %d\n",filename,*dimensions);


    len = (*numObjs) * (*dimensions);
    objects    = (float**)malloc((*numObjs) * sizeof(float*));
    assert(objects != NULL);
    objects[0] = (float*) malloc(len * sizeof(float));
    assert(objects[0] != NULL);
    for (i=1; i<(*numObjs); i++)
        objects[i] = objects[i-1] + (*dimensions);

    i = 0;
    while (fgets(line, lineLen, infile) != NULL) {
        if (strtok(line, " \t\n") == NULL) continue;
        for (j=0; j<(*dimensions); j++)
            objects[i][j] = atof(strtok(NULL, " ,\t\n"));
        i++;
    }

    fclose(infile);
    free(line);

    return objects;
}


int file_write(char      *filename,     /* input file name */
               int        numClusters,  /* no. clusters */
               int        numObjs,      /* no. data objects */
               int        dimensions,    /* no. coordinates (local) */
               float    **clusters,     /* [numClusters][dimensions] centers */
               int       *membership)   /* [numObjs] */
{
    FILE *fptr;
    int   i, j;
    char  outFileName[1024];


    sprintf(outFileName, "%s.cluster_centres", filename);
    printf("Writing coordinates of K=%d cluster centers to file \"%s\"\n",
           numClusters, outFileName);
    fptr = fopen(outFileName, "w");
    for (i=0; i<numClusters; i++) {
        fprintf(fptr, "%d ", i);
        for (j=0; j<dimensions; j++)
            fprintf(fptr, "%f ", clusters[i][j]);
        fprintf(fptr, "\n");
    }
    fclose(fptr);


    sprintf(outFileName, "%s.membership", filename);
    printf("Writing membership of N=%d data objects to file \"%s\"\n",
           numObjs, outFileName);
    fptr = fopen(outFileName, "w");
    for (i=0; i<numObjs; i++)
        fprintf(fptr, "%d %d\n", i, membership[i]);
    fclose(fptr);

    return 1;
}


float** mpi_read(char     *filename,      /* input file name */
                 int      *numObjs,       /* no. data objects (local) */
                 int      *dimensions,     /* no. coordinates */
                 MPI_Comm  comm)
{
    float    **objects;
    int        i, j, len, divd, rem;
    int        rank, nproc;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

        if (rank == 0) {
            objects = file_read(filename, numObjs, dimensions);
            if (objects == NULL) *numObjs = -1;
        }

        MPI_Bcast(numObjs,   1, MPI_INT, 0, comm);
        MPI_Bcast(dimensions, 1, MPI_INT, 0, comm);

        if (*numObjs == -1) {
            MPI_Finalize();
            exit(1);
        }

        divd = (*numObjs) / nproc;
        rem  = (*numObjs) % nproc;

        if (rank == 0) {
            int index = (rem > 0) ? divd+1 : divd;

            (*numObjs) = index;

            for (i=1; i<nproc; i++) {
                int msg_size = (i < rem) ? (divd+1) : divd;
                MPI_Send(objects[index], msg_size*(*dimensions), MPI_FLOAT,
                         i, i, comm);
                index += msg_size;
            }


            objects[0] = realloc(objects[0],
                                 (*numObjs)*(*dimensions)*sizeof(float));
            assert(objects[0] != NULL);
            objects    = realloc(objects, (*numObjs)*sizeof(float*));
            assert(objects != NULL);
        }
        else {
            (*numObjs) = (rank < rem) ? divd+1 : divd;

            objects    = (float**)malloc((*numObjs)            *sizeof(float*));
            assert(objects != NULL);
            objects[0] = (float*) malloc((*numObjs)*(*dimensions)*sizeof(float));
            assert(objects[0] != NULL);
            for (i=1; i<(*numObjs); i++)
                objects[i] = objects[i-1] + (*dimensions);

            MPI_Recv(objects[0], (*numObjs)*(*dimensions), MPI_FLOAT, 0,
                     rank, comm, &status);
        }


    return objects;
}



int mpi_write(char      *filename,     /* input file name */
              int        numClusters,  /* no. clusters */
              int        numObjs,      /* no. data objects */
              int        dimensions,    /* no. coordinates (local) */
              float    **clusters,     /* [numClusters][dimensions] centers */
              int       *membership,   /* [numObjs] */
              int        totalNumObjs, /* total no. data objects */
              MPI_Comm   comm)
{
    int        divd, rem, len, err;
    int        i, j, k, rank, nproc;
    char       outFileName[1024], fs_type[32], str[32], *delim;
    MPI_File   fh;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    delim = strchr(filename, ':');
    if (delim != NULL) {
        strncpy(fs_type, filename, delim-filename);
        fs_type[delim-filename] = '\0';
        delim++;
    }
    else
        delim = filename;


    if (rank == 0) {
        printf("Writing coordinates of K=%d cluster centers to file \"%s.cluster_centres\"\n",
               numClusters, delim);
        sprintf(outFileName, "%s.cluster_centres", filename);
        err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n", outFileName,errstr);
            MPI_Finalize();
            exit(1);
        }


            for (i=0; i<numClusters; i++) {
                sprintf(str, "%d ", i);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                for (j=0; j<dimensions; j++) {
                    sprintf(str, "%f ", clusters[i][j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
                MPI_File_write(fh, "\n", 1, MPI_CHAR, &status);
            }

        MPI_File_close(&fh);
    }

    if (rank == 0)
        printf("Writing membership of N=%d data objects to file \"%s.membership\"\n",
               totalNumObjs, delim);


        if (rank == 0) {
            int divd = totalNumObjs / nproc;
            int rem  = totalNumObjs % nproc;

            sprintf(outFileName, "%s.membership", filename);
            err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            if (err != MPI_SUCCESS) {
                char errstr[MPI_MAX_ERROR_STRING];
                int  errlen;
                MPI_Error_string(err, errstr, &errlen);
                printf("Error at opening file %s (%s)\n", outFileName,errstr);
                MPI_Finalize();
                exit(1);
            }


            for (j=0; j<numObjs; j++) {
                sprintf(str, "%d %d\n", j, membership[j]);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
            }

            k = numObjs;
            for (i=1; i<nproc; i++) {
                numObjs = (i < rem) ? divd+1 : divd;
                MPI_Recv(membership, numObjs, MPI_INT, i, i, comm, &status);

                for (j=0; j<numObjs; j++) {
                    sprintf(str, "%d %d\n", k++, membership[j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
            }
            MPI_File_close(&fh);
        }
        else {
            MPI_Send(membership, numObjs, MPI_INT, 0, rank, comm);
        }

    return 1;
}

__inline static
float euclid_dist_2(int    dimensions,  /* no. dimensions */
                    float *x1,   /* [dimensions] */
                    float *x2)   /* [dimensions] */
{
    int i;
    float ans=0.0;

    for (i=0; i<dimensions; i++)
        ans += (x1[i]-x2[i]) * (x1[i]-x2[i]);

    return(ans);
}


__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     dimensions,   /* no. coordinates */
                         float  *object,      /* [dimensions] */
                         float **clusters)    /* [numClusters][dimensions] */
{
    int   j, i;
    float dist, min_d;

    j    = 0;
    min_d = euclid_dist_2(dimensions, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(dimensions, object, clusters[i]);
        if (dist < min_d) {
            min_d = dist;
            j    = i;
        }
    }
    return(j);
}


int mpi_kmeans(float    **objects,     /* in: [numObjs][dimensions] */
               int        dimensions,   /* no. coordinates */
               int        numObjs,     /* no. objects */
               int        numClusters, /* no. clusters */
               float      threshold,   /* % objects change membership */
               int       *membership,  /* out: [numObjs] */
               float    **clusters,    /* out: [numClusters][dimensions] */
               MPI_Comm   comm)        /* MPI communicator */
{
    int      i, j, rank, index, loop=0, total_numObjs;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    int     *clusterSize;    /* [numClusters]: temp buffer for Allreduce */
    float    delta;          /* % of objects change their clusters */
    float    delta_tmp;
    float  **newClusters;    /* [numClusters][dimensions] */


     MPI_Comm_rank(comm, &rank);

    for (i=0; i<numObjs; i++) membership[i] = -1;

    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);
    clusterSize    = (int*) calloc(numClusters, sizeof(int));
    assert(clusterSize != NULL);

    newClusters    = (float**) malloc(numClusters *            sizeof(float*));
    assert(newClusters != NULL);
    newClusters[0] = (float*)  calloc(numClusters * dimensions, sizeof(float));
    assert(newClusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        newClusters[i] = newClusters[i-1] + dimensions;

    MPI_Allreduce(&numObjs, &total_numObjs, 1, MPI_INT, MPI_SUM, comm);
     printf("%2d: numObjs=%d total_numObjs=%d numClusters=%d dimensions=%d\n",rank,numObjs,total_numObjs,numClusters,dimensions);

    do {
        double curT = MPI_Wtime();
        delta = 0.0;
        for (i=0; i<numObjs; i++) {

            index = find_nearest_cluster(numClusters, dimensions, objects[i],
                                         clusters);


            if (membership[i] != index) delta += 1.0;


            membership[i] = index;


            newClusterSize[index]++;
            for (j=0; j<dimensions; j++)
                newClusters[index][j] += objects[i][j];
        }


        MPI_Allreduce(newClusters[0], clusters[0], numClusters*dimensions,
                      MPI_FLOAT, MPI_SUM, comm);
        MPI_Allreduce(newClusterSize, clusterSize, numClusters, MPI_INT,
                      MPI_SUM, comm);


        for (i=0; i<numClusters; i++) {
            for (j=0; j<dimensions; j++) {
                if (clusterSize[i] > 1)
                    clusters[i][j] /= clusterSize[i];
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }

        MPI_Allreduce(&delta, &delta_tmp, 1, MPI_FLOAT, MPI_SUM, comm);
        delta = delta_tmp / total_numObjs;


            double maxTime;
            curT = MPI_Wtime() - curT;
            MPI_Reduce(&curT, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
            if (rank == 0) printf("%2d: loop=%d time=%f sec\n",rank,loop,curT);

    } while (delta > threshold && loop++ < 500);

    if (rank == 0) printf("%2d: delta=%f threshold=%f loop=%d\n",rank,delta,threshold,loop);

    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);
    free(clusterSize);

    return 1;
}


static void usage(char *argv0, float threshold) {
    char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -i filename    : file containing data to be clustered\n"
        "       -n num_clusters: number of clusters (K must > 1)\n"
        "       -t threshold   : threshold value (default %.4f)\n"
        "       -o             : output timing results (default no)\n";

    fprintf(stderr, help, argv0, threshold);
}



int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     i, j;
           int     is_output_timing, is_print_usage;

           int     numClusters, dimensions, numObjs, totalNumObjs;
           int    *membership;    /* [numObjs] */
           char   *filename;
           float **objects;       /* [numObjs][dimensions] data objects */
           float **clusters;      /* [numClusters][dimensions] cluster center */
           float   threshold;
           double  timing, io_timing, clustering_timing;

           int        rank, nproc, mpi_namelen;
           char       mpi_name[MPI_MAX_PROCESSOR_NAME];
           MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Get_processor_name(mpi_name,&mpi_namelen);



    threshold        = 0.001;
    numClusters      = 0;
    is_output_timing = 0;
    is_print_usage   = 0;
    filename         = NULL;

    while ( (opt=getopt(argc,argv,"p:i:n:t:aoh"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 't': threshold=atof(optarg);
                      break;
            case 'n': numClusters = atoi(optarg);
                      break;
            case 'o': is_output_timing = 1;
                      break;

            case 'h': is_print_usage = 1;
                      break;
            default: is_print_usage = 1;
                      break;
        }
    }

    if (filename == 0 || numClusters <= 1 || is_print_usage == 1) {
        if (rank == 0) usage(argv[0], threshold);
        MPI_Finalize();
        exit(1);
    }

     printf("Proc %d of %d running on %s\n", rank, nproc, mpi_name);

    MPI_Barrier(MPI_COMM_WORLD);
    io_timing = MPI_Wtime();


    objects = mpi_read(filename, &numObjs, &dimensions,
                       MPI_COMM_WORLD);

        int num = (numObjs < 4) ? numObjs : 4;
        for (i=0; i<num; i++) {
            char strline[1024], strfloat[16];
            sprintf(strline,"%d: objects[%d]= ",rank,i);
            for (j=0; j<dimensions; j++) {
                sprintf(strfloat,"%10f",objects[i][j]);
                strcat(strline, strfloat);
            }
            strcat(strline, "\n");
            printf("%s",strline);
        }


    timing            = MPI_Wtime();
    io_timing         = timing - io_timing;
    clustering_timing = timing;

    clusters    = (float**) malloc(numClusters *             sizeof(float*));
    assert(clusters != NULL);
    clusters[0] = (float*)  malloc(numClusters * dimensions * sizeof(float));
    assert(clusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        clusters[i] = clusters[i-1] + dimensions;

    MPI_Allreduce(&numObjs, &totalNumObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        for (i=0; i<numClusters; i++)
            for (j=0; j<dimensions; j++)
                clusters[i][j] = objects[i][j];
    }
    MPI_Bcast(clusters[0], numClusters*dimensions, MPI_FLOAT, 0, MPI_COMM_WORLD);

    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);


    mpi_kmeans(objects, dimensions, numObjs, numClusters, threshold, membership,
               clusters, MPI_COMM_WORLD);

    free(objects[0]);
    free(objects);

    timing            = MPI_Wtime();
    clustering_timing = timing - clustering_timing;


    mpi_write( filename, numClusters, numObjs, dimensions,
              clusters, membership, totalNumObjs, MPI_COMM_WORLD);

    free(membership);
    free(clusters[0]);
    free(clusters);


    if (is_output_timing) {
        double max_io_timing, max_clustering_timing;

        io_timing += MPI_Wtime() - timing;


        MPI_Reduce(&io_timing, &max_io_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&clustering_timing, &max_clustering_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("\nPerforming **** Simple Kmeans  (MPI) ****\n");
            printf("Num of processes = %d\n", nproc);
            printf("Input file:        %s\n", filename);
            printf("numObjs          = %d\n", totalNumObjs);
            printf("dimensions        = %d\n", dimensions);
            printf("numClusters      = %d\n", numClusters);
            printf("threshold        = %.4f\n", threshold);

            printf("I/O time           = %10.4f sec\n", max_io_timing);
            printf("Computation timing = %10.4f sec\n", max_clustering_timing);
        }
    }

    MPI_Finalize();
    return(0);
}
