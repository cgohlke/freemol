EXTERN struct t_multi_query {
    int surf, volume, boxsize, gmmx, conn, logp;
    int nfile, *ext_files, files,files1;
    int stop_analyze, method;
    float dist[50], ang[50], dihed[50], pmr[50];
    float coffdist[50], coffang[50], coffdihed[50], coffpmr[50];
    float rms;
    char name1[60],name2[60];
    }  multi_query;

