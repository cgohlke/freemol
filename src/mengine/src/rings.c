#define EXTERN extern

#include "pcwin.h"
#include "rings.h"
#include "utility.h"

int have_ring3(void);
int have_ring4(void);
int have_ring5(void);
int have_ring6(void);
int isbond(int, int);
int is_ring31(int);
int is_ring33(int *);
int is_ring41(int);
int is_ring42(int,int);
int is_ring43(int, int, int);
int is_ring44(int *);
int is_ring51(int);
int is_ring52(int, int);
int is_ring53(int, int, int);
int is_ring54(int, int, int, int);
int is_ring55(int *);
int is_ring61(int);
int is_ring62(int, int);
void get_ring62(int,int,int *);
int is_ring63(int,int,int);
int is_ring66(int *);
int is_cyclo5(int, int *);
void ksort(int, int *);
void message_alert(char *, char *);
int is_cyclo6(int, int *);
int find_rsize(int, int);
void get_rsize(int,int,int,int *);
int aromatic_5(int *array,long int *flags,int *atomnum,int *nconnect,int *input_charge,int **,int **);
int aromatic_6(int *,long int *flags,int *atomnum,int *nconnect,int *input_charge,int **,int **);
void get_rings(int natom,long int *flags,int *atomnum,int **iat,int **bo);
void allocate_rings(int niatom);
int check_ring1(int);
int has_double_bond(int,int *,int **);
int check_connectivity(int ia,int *atomnum,int *nconnect,int *input_charge);
int check_connectivity_yatom(int ia,int *atomnum,int *nconnect,int *input_charge);
int icompare(int,int *,int *);
int get_bondorder(int ia,int ib);
void free_rings(int);
// ===========================
void allocate_rings(int niatom)
{
  rings.nring3 = 0;
  rings.nring4 = 0;
  rings.nring5 = 0;
  rings.nring6 = 0;
  rings.r13 = imatrix(0,(niatom+10)/3, 0,3);
  rings.r14 = imatrix(0,(niatom+10)/4, 0,4);
  rings.r15 = imatrix(0,(niatom+10)/5, 0,5);
  rings.r16 = imatrix(0,(niatom+10)/2, 0,6);
}
// =========================
void free_rings(int niatom)
{
  free_imatrix(rings.r13 ,0,(niatom+10)/3, 0,3);
  free_imatrix(rings.r14 ,0,(niatom+10)/4, 0,4);
  free_imatrix(rings.r15 ,0,(niatom+10)/5, 0,5);
  free_imatrix(rings.r16 ,0,(niatom+10)/2, 0,6);
}
// ======================
int have_ring3()
{
  if (rings.nring3 > 0)
    return TRUE;
  else
    return FALSE;
}
int have_ring4()
{
  if (rings.nring4 > 0)
    return TRUE;
  else
    return FALSE;
}
int have_ring5()
{
  if (rings.nring5 > 0)
    return TRUE;
  else
    return FALSE;
}
int have_ring6()
{
  if (rings.nring6 > 0)
    return TRUE;
  else
    return FALSE;
}
/* --------------------------------------------------------- */
// iatom is atom number to search on
// isize = ring size
// num = number of ring if more than one ring containing atom iatom
// array = array of ring atoms
//
void get_rsize(int iatom, int isize, int num, int *array)
{
    int i, j, k, icount;
    icount = -1;
    
    if (isize == 3)
    {
       for (i=0; i < rings.nring3; i++)
       {
           for (j=0; j < 3; j++)
           {
               if (iatom == rings.r13[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r13[i][k];
                  return;
               }
           }
        }
    } else if (isize == 4)
    {
       for (i=0; i < rings.nring4; i++)
       {
           for (j=0; j < 4; j++)
           {
               if (iatom == rings.r14[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r14[i][k];
                  return;
               }
           }
        }
    } else if (isize == 5)
    {
       for (i=0; i < rings.nring5; i++)
       {
           for (j=0; j < 5; j++)
           {
               if (iatom == rings.r15[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r15[i][k];
                  return;
               }
           }
        }
    } else if (isize == 6)
    {
       for (i=0; i < rings.nring6; i++)
       {
           for (j=0; j < 6; j++)
           {
               if (iatom == rings.r16[i][j])
                  icount++;
               if (icount == num)
               {
                  for (k=0; k < isize; k++)
                      array[k] = rings.r16[i][k];
                  return;
               }
           }
        }
    } 
}
/* --------------------------------------------------------- */
int find_rsize(int isize,int iatom)
{
    int i, j, icount;
    
    icount = 0;
    if (isize == 3)
    {
        for (i=0; i < rings.nring3; i++)
        {
            for(j=0; j < 3; j++)
            {
                if (rings.r13[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 4)
    {
        for (i=0; i < rings.nring4; i++)
        {
            for(j=0; j < 4; j++)
            {
                if (rings.r14[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 5)
    {
        for (i=0; i < rings.nring5; i++)
        {
            for(j=0; j < 5; j++)
            {
                if (rings.r15[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    }else if (isize == 6)
    {
        for (i=0; i < rings.nring6; i++)
        {
            for(j=0; j < 6; j++)
            {
                if (rings.r16[i][j] == iatom)
                    icount++;
            }
        }
        return(icount);
    } 
    return(icount);
}
// =======================
int is_ring31(int ia)
{
    int i,j;

    for (i=0; i < rings.nring3; i++)
    {
        for (j=0; j < 3; j++)
        {
            if (ia == rings.r13[i][j])
            {
               return(TRUE);
            }
        }
    }
    return FALSE;
}
/* ============================ */
int is_ring41(int ia)
{
    int i,j, found;

    found = FALSE;
    for (i=0; i < rings.nring4; i++)
    {
        for (j=0; j < 4; j++)
        {
            if (ia == rings.r14[i][j])
            {
                found = TRUE;
                return (found);
            }
        }
    }
    return (found);
}
/* -------------------------------------------------------- */            
int is_ring42(int ia, int ib)
{
    int i,j, k, itmp, jtmp;
    if (ia < ib )
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
    for (i=0; i < rings.nring4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if ( itmp == rings.r14[i][j])
            {
                for (k=0; k < 4; k++)
                {
                    if (jtmp == rings.r14[i][k])
                       return(TRUE);
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */            
int is_ring43(int ia, int ib, int ic)
{
    int i, j,k,l;
    for (i=0; i < rings.nring4; i++)
    {
        for (j=0; j < 4; j++)
        {
            if (rings.r14[i][j] == ia)
            {
                for (k=0; k < 4; k++)
                {
                    if (rings.r14[i][k] == ib)
                    {
                        for (l=0; l < 4; l++)
                        {
                            if (rings.r14[i][l] == ic)
                                return TRUE;
                        }
                    }
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */    
int is_ring33(int *array)
{
    int i;
    for (i=0; i < rings.nring3; i++)
    {
        if (array[0] == rings.r13[i][0] && array[1] == rings.r13[i][1] && array[2] == rings.r13[i][2])
            return(TRUE);
    }
    return FALSE;
}
/* -------------------------------------------------------- */    
int is_ring44(int *array)
{
    int i;
    for (i=0; i < rings.nring4; i++)
    {
        if (array[0] == rings.r14[i][0] && array[1] == rings.r14[i][1] &&
            array[2] == rings.r14[i][2] && array[3] == rings.r14[i][3])
            return(TRUE);
    }
    return FALSE;
}
/* -------------------------------------------------------- */    
int is_ring51(int ia)
{
    int i,j;

    for(i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (rings.r15[i][j] == ia)
            {
                return(TRUE);
            }
        }
    }
    return(FALSE);
}
// ======================================
int is_cyclo5(int ia, int *array)
{
    int i,j;

    for(i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (rings.r15[i][j] == ia)
            {
	      array[0] = rings.r15[i][0];
	      array[1] = rings.r15[i][1];
	      array[2] = rings.r15[i][2];
	      array[3] = rings.r15[i][3];
	      array[4] = rings.r15[i][4];
              return(TRUE);
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */
int is_ring52(int ia, int ib)
{
    int i, j,k, itmp, jtmp;
    if (ia < ib)
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
    for (i=0; i < rings.nring5; i++)
    {
       for (j=0; j < 5; j++)
       {
           if (itmp == rings.r15[i][j])
           {
               for (k=0; k < 5; k++)
               {
                   if (jtmp == rings.r15[i][k])
                       return(TRUE);
               }
           }
       }
    }
    return(FALSE);
}     
/* -------------------------------------------------------- */    
int is_ring53(int ia, int ib, int ic)
{
    int i, j,k,l;
    for (i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (rings.r15[i][j] == ia)
            {
                for (k=0; k < 5; k++)
                {
                    if (rings.r15[i][k] == ib)
                    {
                        for (l=0; l < 5; l++)
                        {
                            if (rings.r15[i][l] == ic)
                                return(TRUE);
                        }
                    }
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */
int is_ring54(int ia, int ib, int ic, int id)
{
    int i,j,k,l,m;

    for (i=0; i < rings.nring5; i++)
    {
        for (j=0; j < 5; j++)
        {
            if (rings.r15[i][j] == ia)
            {
                for (k=0; k < 5; k++)
                {
                    if (rings.r15[i][k] == ib)
                    {
                        for (l=0; l < 5; l++)
                        {
                            if (rings.r15[i][l] == ic)
                            {
                                for (m=0; m < 5; m++)
                                {
                                    if (rings.r15[i][m] == id)
                                        return TRUE;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
int is_ring55(int *array)
{
    int i;
    for(i=0; i < rings.nring5; i++)
    {
        if ( array[0] == rings.r15[i][0] && array[1] == rings.r15[i][1] && array[2] == rings.r15[i][2] &&
             array[3] == rings.r15[i][3] && array[4] == rings.r15[i][4])
             return(TRUE);
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
int is_ring61(int ia)
{
    int i, j;
    for(i=0; i < rings.nring6; i++)
    {
        for(j=0; j < 6; j++)
        {
            if (ia == rings.r16[i][j])
               return (TRUE);
        }
    }
    return (FALSE);
}
/* -------------------------------------------------------- */
int is_ring62(int ia, int ib)
{
    int i, j,k, itmp, jtmp;
    if (ia < ib)
    {
        itmp = ia;
        jtmp = ib;
    }else
    {
        itmp = ib;
        jtmp = ia;
    }
    for (i=0; i < rings.nring6; i++)
    {
       for (j=0; j < 6; j++)
       {
           if (itmp == rings.r16[i][j])
           {
               for (k=0; k < 6; k++)
               {
                   if (jtmp == rings.r16[i][k])
                       return(TRUE);
               }
           }
       }
    }
    return(FALSE);
} 
// ============================================
void get_ring62(int ia,int ib,int *array)
{
  int i,j,k,l;
  for (i=0; i < rings.nring6; i++)
    {
      for (j=0; j < 6; j++)
	{
	  if (rings.r16[i][j] == ia)
	    {
	      for (k=0; k < 6; k++)
		{
		  if (rings.r16[i][k] == ib)
		    {
		      for (l=0; l < 6; l++)
			  array[l] = rings.r16[i][l]; 
		      return;
		    }
		}
	    }
	}
    }
}   
/* -------------------------------------------------------- */    
int is_ring63(int ia, int ib, int ic)
{
    int i, j, k, l;

    for (i=0; i < rings.nring6; i++)
    {
        for (j=0; j < 6; j++)
        {
            if (rings.r16[i][j] == ia)
            {
                for (k=0; k < 6; k++)
                {
                    if (rings.r16[i][k] == ib)
                    {
                        for (l=0; l < 6; l++)
                        {
                            if (rings.r16[i][l] == ic)
                                return(TRUE);
                        }
                    }
                }
            }
        }
    }
    return(FALSE);
}
/* -------------------------------------------------------- */    
int is_ring66(int *array)
{
    int i;
    for(i=0; i < rings.nring6; i++)
    {
        if ( array[0] == rings.r16[i][0] && array[1] == rings.r16[i][1] && array[2] == rings.r16[i][2] &&
             array[3] == rings.r16[i][3] && array[4] == rings.r16[i][4] && array[5] == rings.r16[i][5])
             return(TRUE);
    }
    return (FALSE);
}
/* -------------------------------------------------------- */    
void get_rings(int natom,long int *flags,int *atomnum, int **iat,int **bo)
{
   int i,j, k, l, m, n,unique;
   int jj,kk;
   int jatm,katm,latm,xatm;
   int array[7],array1[7],array2[7];

   rings.nring3 = 0;
   rings.nring4 = 0;
   rings.nring5 = 0;
   rings.nring6 = 0;

// get three membered rings
   for (i=1; i <= natom; i++)
     {
       for (j=0; j < MAXIAT; j++)
	 {
	   if (iat[i][j] != 0)
	     {
	       jatm = iat[i][j];
	       for (k=0; k < MAXIAT; k++)
		 {
		   if (iat[jatm][k] != 0 && iat[jatm][k] != i)
		     {
		       if (isbond(i,iat[jatm][k]))
			 {
			   array[0] = i;
			   array[1] = jatm;
			   array[2] = iat[jatm][k];
			   unique = TRUE;
                           for (jj=0; jj < 3; jj++)
                              array1[jj] = array[jj];
                           for (jj=0; jj < rings.nring3; jj++)
                           {
                               for(kk=0; kk < 3; kk++)
                                array2[kk] = rings.r13[jj][kk];
                               if (icompare(3,array1,array2))
                                unique = FALSE;
                           }
			   if (unique)
			     {  
			       if (rings.nring3 < (natom+10)/3)
				 {                            
				   rings.r13[rings.nring3][0] = array[0];
				   rings.r13[rings.nring3][1] = array[1];
				   rings.r13[rings.nring3][2] = array[2];
				   rings.nring3++;
				 }
			     }
			 }
		     }
		 }
	     }
	 }
     }
// get four membered rings
   for(i=1; i <= natom; i++)
   {
       for(j=0; j < MAXIAT; j++)
       {
           if (iat[i][j] != 0)
           {
               jatm = iat[i][j];
               for (k=0; k < MAXIAT; k++)
               {
                   if (iat[jatm][k] != 0 && iat[jatm][k] != i)
                   {
                       katm = iat[jatm][k];
                       if (!(isbond(katm,i)))
                       {
                           for(l=0; l < MAXIAT; l++)
                           {
                               if (iat[katm][l] != 0 && iat[katm][l] != jatm && iat[katm][l] != i)
                               {
                                   if (isbond(i,iat[katm][l]) )
                                   {
                                        array[0] = i;
                                        array[1] = jatm;
                                        array[2] = katm;
                                        array[3] = iat[katm][l];
                                        unique = TRUE;
                                        for (jj=0; jj < 4; jj++)
                                           array1[jj] = array[jj];
                                        for (jj=0; jj < rings.nring4; jj++)
                                        {
                                            for(kk=0; kk < 4; kk++)
                                              array2[kk] = rings.r14[jj][kk];
                                             if (icompare(4,array1,array2))
                                               unique = FALSE;
                                        }
                                        if ( unique )
                                        {          
					  if (rings.nring4 < (natom+10)/4)
					    {                            
					      rings.r14[rings.nring4][0] = array[0];
					      rings.r14[rings.nring4][1] = array[1];
					      rings.r14[rings.nring4][2] = array[2];
					      rings.r14[rings.nring4][3] = array[3];
					      rings.nring4++;
					    }
                                        }
                                   }
                               }
                           }
                       }
                   }
               }
           }
       }
   }
// get five membered rings i-jatm-katm-latm-
   for(i=1; i <= natom; i++)
   {
       for(j=0; j < MAXIAT; j++)
       {
           if (iat[i][j] != 0)
           {
               jatm = iat[i][j];
               for (k=0; k < MAXIAT; k++)
               {
                   if (iat[jatm][k] != 0 && iat[jatm][k] != i)
                   {
                       katm = iat[jatm][k];
                       if (!(isbond(katm,i)))
                       {
                           for(l=0; l < MAXIAT; l++)
                           {
                               if (iat[katm][l] != 0 && iat[katm][l] != jatm && !(isbond(iat[katm][l],i)))
                               {
				 latm = iat[katm][l];
				 for (m=0; m < MAXIAT; m++)
				   {
				     if (iat[latm][m] != 0 && iat[latm][m] != katm && iat[latm][m] != jatm)
				       {
					 if (isbond(i,iat[latm][m]))
					   {
					     array[0] = i;
					     array[1] = jatm;
					     array[2] = katm;
					     array[3] = latm;
					     array[4] = iat[latm][m];
                                             unique = TRUE;
                                             for (jj=0; jj < 5; jj++)
                                                array1[jj] = array[jj];
                                             for (jj=0; jj < rings.nring5; jj++)
                                             {
                                                 for(kk=0; kk < 5; kk++)
                                                   array2[kk] = rings.r15[jj][kk];
                                                 if (icompare(5,array1,array2))
                                                   unique = FALSE;
                                             }
					     if (unique)
					       {
						 if (rings.nring5 < (natom+10)/5)
						   {
						     rings.r15[rings.nring5][0] = array[0];
						     rings.r15[rings.nring5][1] = array[1];
						     rings.r15[rings.nring5][2] = array[2];
						     rings.r15[rings.nring5][3] = array[3];
						     rings.r15[rings.nring5][4] = array[4];
						     rings.nring5++;
						   }
					       }
					   }
				       }
				   }
			       }
			   }
		       }
		   }
	       }
	   }
       }
   }
// get six membered rings
   for(i=1; i <= natom; i++)
   {
       for(j=0; j < MAXIAT; j++)
       {
           if (iat[i][j] != 0)
           {
               jatm = iat[i][j];
               for (k=0; k < MAXIAT; k++)
               {
                   if (iat[jatm][k] != 0 && iat[jatm][k] != i)
                   {
                       katm = iat[jatm][k];
                       if (!(isbond(katm,i)))  // not three
                       {
                           for(l=0; l < MAXIAT; l++)
                           {
                               if (iat[katm][l] != 0 && iat[katm][l] != jatm && iat[katm][l] != i)
                               {
				 latm = iat[katm][l];
				 if (!(isbond(i,latm))) // not four
				   {
				     for (m=0; m < MAXIAT; m++)
				       {
					 if (iat[latm][m] != 0 && iat[latm][m] != katm && iat[latm][m] != jatm && iat[latm][m] != i)
					   {
					     xatm = iat[latm][m];
					     if (!(isbond(i,xatm))) // not five
					       {
						 for (n = 0; n < MAXIAT; n++)
						   {
						     if (iat[xatm][n] != 0 && iat[xatm][n] != latm && iat[xatm][n] != katm && iat[xatm][n] != jatm && iat[xatm][n] != i)
						       if (isbond(i,iat[xatm][n]))
							 {
							   array[0] = i;
							   array[1] = jatm;
							   array[2] = katm;
							   array[3] = latm;
							   array[4] = xatm;
							   array[5] = iat[xatm][n];
							   unique = TRUE;
							   for (jj=0; jj < 6; jj++)
							     array1[jj] = array[jj];
							   for (jj=0; jj < rings.nring6; jj++)
							     {
							       for(kk=0; kk < 6; kk++)
								 array2[kk] = rings.r16[jj][kk];
							       if (icompare(6,array1,array2))
								 unique = FALSE;
							     }
							   if (unique)
							     {
							       if (rings.nring6 < (natom+10)/2)
								 {
								   rings.r16[rings.nring6][0] = array[0];
								   rings.r16[rings.nring6][1] = array[1];
								   rings.r16[rings.nring6][2] = array[2];
								   rings.r16[rings.nring6][3] = array[3];
								   rings.r16[rings.nring6][4] = array[4];
								   rings.r16[rings.nring6][5] = array[5];
								   rings.nring6++;
								 }
							     }
							 }
						   }
					       }
					   }
				       }
				   }
			       }
			   }
		       }
		   }
	       }
	   }
       }
   }
}
/* =================================================== */
void ksort(int num, int array[])
{
    int i,temp;
    int found;

L_1:
    found = FALSE;
    for (i=0; i < num-1; i++)
    {
       if (array[i+1] < array[i])
       {
          temp = array[i];
          array[i] = array[i+1];
          array[i+1] = temp;
          found = TRUE;
       }
    }
    if (found == TRUE)
       goto L_1;
}
/* --------------------------------------------------------- */
int aromatic_5(int *array,long int *flags,int *atomnum,int *nconnect,int *input_charge,int **iat,int **bo)
{
  int i,j,npi,ihetero,found;
  int jatm,katm;
  int ia,ib,ic,idbl,jdbl;
  int db[5];
    
    npi = 0;
    ihetero = 0;
    for (i=0; i < 5; i++)
      db[i] = 0;
    // printf("array: %d %d %d %d %d\n",array[0],array[1],array[2],array[3],array[4]);
    for (i=0; i < 5; i++)
      {
	if (has_double_bond(array[i],nconnect,bo) == TRUE)
	  {
	    if (check_connectivity(array[i],atomnum,nconnect,input_charge) == FALSE)
	      {
		return FALSE;
	      }
	  } else
	  {
	    if (check_connectivity_yatom(array[i],atomnum,nconnect,input_charge) == FALSE)
	      {
		return FALSE;
	      }
	  }
      }
    // printf("all aromatic\n");
    //
    for (i=0; i < 5; i++)
      {
	jatm = array[i];
	if (atomnum[jatm] > 6)
	  {
	    found = FALSE;
	    for (j=0; j < nconnect[jatm]; j++)
	      {
		if (bo[jatm][j] == 2)
		  found = TRUE;
	      }
	    if (found == FALSE) ihetero = jatm;
	  }
      }
    for (i=0; i < 4; i++)
      {
	jatm = array[i];
	for (j=i+1; j < 5; j++)
	  {
	    katm = array[j];
	    if (isbond(jatm,katm) && get_bondorder(jatm,katm) == 2)
	      {
		npi++;
		if (i==0 && j == 1)
		  db[0] = 1;
		else if (i == 0 && j == 4)
		  db[4] = 1;
		else 
		  db[i] = 1;
	      }
	  }
      }
    //printf("aro 5: %d %d - %d %d %d %d %d\n",npi,ihetero,db[0],db[1],db[2],db[3],db[4]);
    ia = ib = ic = 0;
    if (npi == 2 && ihetero != 0) return TRUE;  // two double bonds and one hetero atom - aromatic
    if (npi == 0) return FALSE; // no double bonds in ring not aromatic
    if (npi == 1) // 
      {
        if (db[0] == 1)
        {
            ia = array[2];
            ib = array[3];
            ic = array[4];
        } else if (db[1] == 1)
        {
            ia = array[0];
            ib = array[3];
            ic = array[4];
        } else if (db[2] == 1)
        {
            ia = array[0];
            ib = array[1];
            ic = array[4];
        } else if (db[3] == 1)
        {
            ia = array[0];
            ib = array[1];
            ic = array[2];
        } else if (db[4] == 1)
        {
            ia = array[1];
            ib = array[2];
            ic = array[3];
        }
	if (ia == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < nconnect[ib]; i++)
	      {
		if (bo[ib][i] == 2 && is_ring61(iat[ib][i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < nconnect[ic]; i++)
	      {
		if (bo[ic][i] == 2 && is_ring61(iat[ic][i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  } else if (ib == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < nconnect[ia]; i++)
	      {
		if (bo[ia][i] == 2 && is_ring61(iat[ia][i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < nconnect[ic]; i++)
	      {
		if (bo[ic][i] == 2 && is_ring61(iat[ic][i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  } else if (ic == ihetero)
	  {
	    idbl = FALSE;
	    jdbl = FALSE;
	    for (i=0; i < nconnect[ib]; i++)
	      {
		if (bo[ib][i] == 2 && is_ring61(iat[ib][i]))
		  idbl = TRUE;
	      }
	    for (i=0; i < nconnect[ia]; i++)
	      {
		if (bo[ia][i] == 2 && is_ring61(iat[ia][i]))
		  jdbl = TRUE;
	      }
	    if (idbl && jdbl) return TRUE;
	  }
      }
    return FALSE;                  
}
/* --------------------------------------------------------- */
int aromatic_6(int *array,long int *flags,int *atomnum,int *nconnect,int *input_charge, int **iat,int **bo)
{
  int i,j,k,idbl,jdbl,num;
  int ia,ib,jatm,katm;
  int AROMATIC;
  int inarray[10];
  long int aromatic_mask;

    aromatic_mask = (1L << AROMATIC_MASK);
    AROMATIC = TRUE;   
    for (i=0; i < 6; i++)
    {
      if (check_connectivity(array[i],atomnum,nconnect,input_charge) == FALSE)
	AROMATIC = FALSE;
    }
    if (AROMATIC == FALSE)
      return AROMATIC;
    // hava all aromatic type atoms check number of double bonds
    AROMATIC = FALSE;
    num = 0;
    for (i=0; i < 10; i++)
      inarray[i] = FALSE;

    for (i=0; i < 5; i++)
      {
	if (get_bondorder(array[i],array[i+1]) == 2)
	  {
	    num++;
	    inarray[i] = TRUE;
	    inarray[i+1] = TRUE;
	  }
      }
    if (get_bondorder(array[0],array[5]) == 2)
      {
	 num++;
	 inarray[0] = TRUE;
	 inarray[5] = TRUE;
      }

    if (num == 3) return TRUE;
    //
    if (num == 0)
      return FALSE;
    else if (num >= 1)
      {
	for (i=0; i < 5; i++)
	  {
	    if (inarray[i] == FALSE)
	      {
		ia = array[i];
		for (j=i+1; j < 6; j++)
		  {
		    if (inarray[j] == FALSE)
		      {
			ib = array[j];
			if (isbond(ia,ib))
			  {
			    idbl = FALSE;
			    jdbl = FALSE;
			    jatm = katm = 0;
			    for (k = 0; k < MAXIAT; k++)
			      {
				if (iat[ia][k] != 0 && bo[ia][k] == 2)
				  {
				    if (atomnum[iat[ia][k]] == 6 || atomnum[iat[ia][k]]  == 7)
				      {
					if (is_ring61(iat[ia][k]) || is_ring51(iat[ia][k]))
					  {
					    idbl = TRUE;
					    jatm = iat[ia][k];
					  }
				      }
				  }
				if (iat[ib][k] != 0 && bo[ib][k] == 2)
				  {
				    if (atomnum[iat[ib][k]] == 6 || atomnum[iat[ib][k]]  == 7)
				      {
					if (is_ring61(iat[ib][k])|| is_ring51(iat[ib][k]))
					  {
					    jdbl = TRUE;
					    katm = iat[ib][k];
					  }
				      }
				  }
			      }
			     if (idbl == TRUE && (flags[jatm] & aromatic_mask) &&   // bad test - jatm and katm may not
				  jdbl == TRUE && (flags[katm] & aromatic_mask) ) // have been examined yet
			      {
				num++;
				if (num == 3) return TRUE;
			      }
			  }
			inarray[i] = TRUE;
			inarray[j] = TRUE;
		      }
		  }
	      }
	  }
      }
    return AROMATIC;
}
// ==================================
int check_ring1(int ia)
{
    if (is_ring61(ia)) return TRUE;
    if (is_ring51(ia)) return TRUE;
    if (is_ring31(ia)) return TRUE;
    if (is_ring41(ia)) return TRUE;
    return FALSE;
}
// ==============================
int check_connectivity(int ia,int *atomnum,int *nconnect,int *input_charge)
{
  // C (x3)
  if (atomnum[ia] == 6 && nconnect[ia] == 3)
    return TRUE;
  // N (x2)
  else if (atomnum[ia] == 7 && nconnect[ia] == 2)
    return TRUE;
  // P (x2)
  else if (atomnum[ia] == 15 && nconnect[ia] == 2)
    return TRUE;
  // N+ (x2)
  else if (atomnum[ia] == 7 && nconnect[ia] == 3 && input_charge[ia] == 1)
    return TRUE;
  // P+ (x2)
  else if (atomnum[ia] == 15 && nconnect[ia] == 3 && input_charge[ia] == 1)
    return TRUE;
  // O+ (x2)
  else if (atomnum[ia] == 8 && nconnect[ia] == 2 && input_charge[ia] == 1)
    return TRUE;
  // S+ (x2)
  else if (atomnum[ia] == 16 && nconnect[ia] == 2 && input_charge[ia] == 1)
    return TRUE;

  return FALSE;
}
// ==============================
int check_connectivity_yatom(int ia,int *atomnum,int *nconnect,int *input_charge)
{
  // C- (x3) cyclopentadiene anion
    if (atomnum[ia] == 6 && input_charge[ia] == -1 && nconnect[ia] == 3)
    return TRUE;
  // N- (x2)
  else if (atomnum[ia] == 7 && input_charge[ia] == -1 && nconnect[ia] == 2)
    return TRUE;
  // N  (x3)
  else if (atomnum[ia] == 7 && nconnect[ia] == 3)
    return TRUE;
  // O (x2)
  else if (atomnum[ia] == 8 && nconnect[ia] == 2)
    return TRUE;
  // S (x2)
  else if (atomnum[ia] == 16 && nconnect[ia] == 2)
    return TRUE;
  // P (x2)
  else if (atomnum[ia] == 15 && nconnect[ia] == 3)
    return TRUE;
  return FALSE;
}
//  ====================================
int has_double_bond(int ia,int *nconnect,int **bo)
 {
   int i;
   for (i=0; i < nconnect[ia]; i++)
     {
       if (bo[ia][i] == 2)
	 return TRUE;
     }
   return FALSE;
 }
/* =======================================  */
int icompare(int iatom, int *iarray1, int *iarray2)
{
    int i;
    int scratch1[50],scratch2[50];

    for (i=0; i < iatom; i++)
    {
        scratch1[i] = iarray1[i];
        scratch2[i] = iarray2[i];
    }
    ksort(iatom, scratch1);
    ksort(iatom, scratch2);
    for (i=0; i < iatom; i++)
    {
        if (scratch1[i] != scratch2[i])
            return FALSE;
    }
    return TRUE;
}

