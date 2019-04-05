# Compiler-mini-project
Analysing K-means clustering algorithm and providing several parallel implementations.

Sequential: 
  
  For Single_k implementations : 
    Compilation: g++ filename.cpp -o executable_name -lm
    Execution: ./executable_name input_filename.txt num_clusters(int) 
  
  For Multiple_k implementations:
    Compilation:  g++ filename.cpp -o executable_name -lm
    Execution: ./executable_name input_filename.txt max_numclusters(int)
    
Parallel: 
  
  For Single_k implementations : 
    Compilation: g++ filename.cpp -o executable_name -lm -fopenmp
    Execution: ./executable_name input_filename.txt num_clusters(int) 
  
  For Multiple_k implementations:
   Compilation:  g++ filename.cpp -o executable_name -fopenmp
   Execution: ./executable_name input_filename.txt max_numclusters(int)
