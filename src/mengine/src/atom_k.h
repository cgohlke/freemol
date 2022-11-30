EXTERN struct t_atom_k {
        int natomtype;
        int type[MAXATOMTYPE], valency[MAXATOMTYPE];
        int tclass[MAXCLASS], tclass1[MAXCLASS], tclass2[MAXCLASS];
        char symbol[MAXATOMTYPE][5], description[MAXATOMTYPE][21];
        int  number[MAXATOMTYPE];
        float weight[MAXATOMTYPE];
        int   ligands[MAXATOMTYPE];
        }  atom_k;
EXTERN struct t_atom_def {
    int natomtype;
    int type[MAXATOMTYPE],valency[MAXATOMTYPE],number[MAXATOMTYPE],ligands[MAXATOMTYPE];
    float weight[MAXATOMTYPE];
    char  symbol[MAXATOMTYPE][5];
    }  atom_def;
        

