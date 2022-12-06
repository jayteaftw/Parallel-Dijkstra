/* assert */
#include <assert.h>
/* INFINITY */
#include <math.h>
/* FILE, fopen, fclose, fscanf, rewind */
#include <stdio.h>
/* EXIT_SUCCESS, malloc, calloc, free */
#include <stdlib.h>
/* time, CLOCKS_PER_SEC */
#include <time.h>

#include <omp.h>

#include <iostream>

#include <chrono>

using namespace std;

#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]

static void load(
    const char * const filename,
    int * const np,
    float ** const ap
)
{
    int i, j, n, ret;
    FILE * fp=NULL;
    float * a;

    /* open the file */
    fp = fopen(filename, "r");
    assert(fp);

    /* get the number of nodes in the graph */
    ret = fscanf(fp, "%d", &n);
    assert(1 == ret);

    /* allocate memory for local values */
    a = (float*) malloc(n*n*sizeof(*a));
    assert(a);

    /* read in roots local values */
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
        ret = fscanf(fp, "%f", &a(i,j));
        assert(1 == ret);
        }
    }

    /* close file */
    ret = fclose(fp);
    assert(!ret);

    /* record output values */
    *np = n;
    *ap = a;
}

void dijkstra_omp_saved(  const int s, 
                    const int n, 
                    const float * const a, 
                    float ** const lp,
                    int num_threads){

    int i, j;
    struct float_int {
        float l;
        int u;
    } min;
    char * m;
    float * l;
    int *threads_indx;
    float_int *sub_min;

    omp_set_num_threads(2);
    num_threads = omp_get_max_threads();
    sub_min = (float_int*)malloc(num_threads*sizeof(*sub_min));

    m = (char*) calloc(n, sizeof(*m)); //Cluster
    assert(m);

    l = (float*) malloc(n*sizeof(*l));
    assert(l);
    

    //Threads indicies 
    threads_indx = (int*)calloc(num_threads+1, sizeof(*threads_indx));
    assert(threads_indx);
    int block_size;
    block_size = n / num_threads;
    for(int i=1; i<num_threads+1;++i){
        if ((threads_indx[i-1]+block_size) <= n)
            threads_indx[i] = threads_indx[i-1]+block_size;
        else
            threads_indx[i] = n;
        
    }
    for(i=0;i<num_threads+1;i++){
        cout<<to_string(threads_indx[i])<<endl;
    }
    
    
    
    //Adds all the current path from current source
    for (i=0; i<n; ++i) {
        l[i] = a(i,s);
    }

    m[s] = 1;
    min.u = -1; /* avoid compiler warning */
    
      
    
    //Loops through each source  
    for (i=1; i<n; ++i) {
            /* find local minimum for threadID */
        #pragma omp parallel
        {    
            int tid = omp_get_thread_num();
            sub_min[tid].l = INFINITY;
            sub_min[tid].u = -1;
            for (j=threads_indx[tid]; j<threads_indx[tid+1]; ++j) {
                //cout<<to_string(tid)<<" "<<to_string(j)<<endl;
                if (!m[j] && l[j] < sub_min[tid].l) { //If not in cluster and val < min
                    sub_min[tid].l = l[j];
                    sub_min[tid].u = j;
                }
            }
        }
            

        min.l = INFINITY;
        //Find minimum amongst the threads
        for(int idx = 0; idx <num_threads; ++idx ){
            cout<<"idx "<<to_string(idx)<<" "<<to_string(sub_min[idx].l)<<endl;
            if (sub_min[idx].l < min.l){
                
                min.l = sub_min[idx].l;
                min.u = sub_min[idx].u;
            }
        }

        //Add node to cluster
        m[min.u] = 1;
        
        //Update nodes with min.l
        cout<<to_string(min.l)<<" "<<to_string(min.u)<<endl;    
        for (j=0; j<n; ++j) {
            //If node not in cluster, add current prefix to that node if the sum is less than current value
            if (!m[j] && min.l+a(j,min.u) < l[j]) 
                l[j] = min.l+a(j,min.u);
            if (j < 10)
                cout<<l[j]<<", ";
        }

    
    }
    cout<<"end"<<endl<<endl;;
    free(m);
    *lp = l;
}

void dijkstra_omp(  const int s, 
                    const int n, 
                    const float * const a, 
                    float ** const lp,
                    const int num_threads){

    int i, j;
    struct float_int {
        float l;
        int u;
    } min;
    char * m;
    float * l;
    int *threads_indx;
    float_int *sub_min; //sub cluster min

    
    sub_min = (float_int*)malloc(num_threads*sizeof(*sub_min));

    m = (char*) calloc(n, sizeof(*m)); //Cluster
    assert(m);

    l = (float*) malloc(n*sizeof(*l));
    assert(l);
    
    //Threads indicies 
    threads_indx = (int*)calloc(num_threads+1, sizeof(*threads_indx));
    assert(threads_indx);
    int block_size;
    block_size = (n / num_threads);

    if (n % num_threads != 0){
        block_size += 1;
    }

    for(i=1; i<num_threads;++i){
        if ((threads_indx[i-1]+block_size) <= n)
            threads_indx[i] = threads_indx[i-1]+block_size;
        else
            threads_indx[i] = n;
               
    }
    threads_indx[i] = n;

    for(i=0;i<num_threads+1;i++){
        cout<<to_string(threads_indx[i])<<endl;
    }
    
    
    //Adds all the current path from current source
    for (i=0; i<n; ++i) {
        l[i] = a(i,s);
    }

    m[s] = 1;
    min.u = -1; /* avoid compiler warning */
    
      
    
    //Loops through each source  
     
    i = 1;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        while(i < n){ 
        
            
            #pragma omp barrier
            sub_min[tid].l = INFINITY;
            sub_min[tid].u = -1;
            for (j=threads_indx[tid]; j<threads_indx[tid+1]; ++j) {
                //cout<<to_string(tid)<<" "<<to_string(j)<<endl;
                if (!m[j] && l[j] < sub_min[tid].l) { //If not in cluster and val < min
                    sub_min[tid].l = l[j];
                    sub_min[tid].u = j;
                }
            }
            
            #pragma omp barrier

            //Use Master thread (thread zero) to find minimum amongst all sub clusters
            if (tid == 0){
                //cout<<endl<<"here "<<to_string(tid)<<endl;
                min.l = INFINITY;
                //Find minimum amongst the threads
               // cout<<"n:"<<i<<endl;
                for(int idx = 0; idx <num_threads; ++idx ){
                    //cout<<"idx "<<to_string(idx)<<" "<<to_string(sub_min[idx].l)<<endl;
                    if (sub_min[idx].l < min.l){
                        min.l = sub_min[idx].l;
                        min.u = sub_min[idx].u;
                    }
                }
                //Add node to cluster
                m[min.u] = 1;
                
                //cout<<to_string(min.l)<<" "<<to_string(min.u)<<endl; 
                //Update node count
                i++;   
            }

            #pragma omp barrier

            //Update nodes with min.l
                
            for (j=threads_indx[tid]; j<threads_indx[tid+1]; ++j) {
                //If node not in cluster, add current prefix to that node if the sum is less than current value
                if (!m[j] && min.l+a(j,min.u) < l[j]) 
                    l[j] = min.l+a(j,min.u);
                //if (j < 10)
                    //cout<<l[j]<<", ";
            }
        }
    }
    free(m);
    *lp = l;
}




static void
print_time(const double seconds)
{
  printf("Search Time: %0.06fs\n", seconds);
}

static void
print_numbers(
  const char * const filename,
  const int n,
  const float * const numbers)
{
  int i;
  FILE * fout;

  /* open file */
  if(NULL == (fout = fopen(filename, "w"))) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write numbers to fout */
  for(i=0; i<n; ++i) {
    fprintf(fout, "%10.4f\n", numbers[i]);
  }

  fclose(fout);
}

int
main(int argc, char ** argv)
{
    int n;
    clock_t ts, te;
    float * a, * l;
    
    if(argc < 3){
        printf("Invalid number of arguments.\nUsage: dijkstra <graph> <num_sources> [<output_file>].\n");
        return EXIT_FAILURE;
    }
    /* initialize random seed: */
    srand (time(NULL));
    unsigned int seed = time(NULL);

    /* figure out number of random sources to search from */
    int nsources = atoi(argv[2]);
    assert(nsources > 0);

    /* load data */
    printf("Loading graph from %s.\n", argv[1]);
    load(argv[1], &n, &a);

    omp_set_num_threads(28);
    int num_threads = omp_get_max_threads();

    printf("Performing %d searches from random sources.\n", nsources);
    auto start = chrono::high_resolution_clock::now();
    for(int i=0; i < nsources; ++i){
        //dijkstra_omp(rand_r(&seed) % n, n, a, &l,num_threads );
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << duration.count() << " ms" << endl;

    if(argc >= 4){ 
        printf("Computing result for source 0.\n");
        dijkstra_omp(0, n, a, &l,num_threads);
        printf("Writing result to %s.\n", argv[3]);
        print_numbers(argv[3], n, l);
    }

    free(a);
    free(l);

  return EXIT_SUCCESS;
}
