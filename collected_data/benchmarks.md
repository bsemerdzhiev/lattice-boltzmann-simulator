Sequential runs:
6'000 steps - 4'203ms
50'000 steps - 36'616ms
100'000 steps - 71'630ms

Simulation started
CPU time passed: 76927.375ms

 Performance counter stats for './LBMSim':

    76,490,618,357      task-clock                       #    0.994 CPUs utilized
             4,232      context-switches                 #   55.327 /sec
                 1      cpu-migrations                   #    0.013 /sec
             1,282      page-faults                      #   16.760 /sec
   554,540,041,494      instructions                     #    1.74  insn per cycle
   318,868,060,124      cycles                           #    4.169 GHz
    25,583,360,972      branches                         #  334.464 M/sec
        26,157,236      branch-misses                    #    0.10% of all branches
 #     32.7 %  tma_backend_bound
                                                  #     56.2 %  tma_bad_speculation
                                                  #      0.1 %  tma_frontend_bound
                                                  #     10.9 %  tma_retiring

      76.930698325 seconds time elapsed

      76.442323000 seconds user
       0.044983000 seconds sys


After AOS

borislav@Borko:~/Documents/lbm-sim/build$ perf stat ./LBMSim
Simulation started
CPU time passed: 33944.137ms

 Performance counter stats for './LBMSim':

    33,927,868,342      task-clock                       #    0.999 CPUs utilized
               574      context-switches                 #   16.918 /sec
                 1      cpu-migrations                   #    0.029 /sec
               711      page-faults                      #   20.956 /sec
   478,841,214,032      instructions                     #    3.55  insn per cycle
   134,719,093,844      cycles                           #    3.971 GHz
    24,386,409,490      branches                         #  718.772 M/sec
        24,974,004      branch-misses                    #    0.10% of all branches
 #     19.8 %  tma_backend_bound
                                                  #     13.3 %  tma_bad_speculation
                                                  #      3.0 %  tma_frontend_bound
                                                  #     63.9 %  tma_retiring

      33.945722225 seconds time elapsed

      33.925430000 seconds user
       0.001999000 seconds sys

