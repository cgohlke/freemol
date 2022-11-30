#define EXTERN extern

#include "pcwin.h"

#include <errno.h>

void * malloc_filename ( char * , char * ) ;
int FetchRecord(FILE *, char *);
FILE * fopen_path ( char * , char * , char * ) ;
int is_bond(int, int);
void message_alert(char *, char *);
void clean_string(char *);
// ========================================
void clean_string(char *astring)
{
    int iz,i;
    iz = strlen(astring);
    for (i=0; i < iz; i++)
    {
        if (isprint(astring[i]) == 0)
          astring[i] = ' ';
    }
    for (i = iz -1; i >= 0; i--)
    {
        if (astring[i] != ' ')
        {
            astring[i+1] = '\0';
            break;
        }
    }
}
/*======================================================================*/
void * malloc_filename ( char *path , char *name )
{
  int ix, iz;
  char *filename;
  
  ix = 0;
  if (path != NULL)
     ix = strlen(path);
  iz = strlen(name);
  filename = malloc (ix + 1 + iz + 1) ;
  if ( ix  == 0 )
    {
      strcpy (filename,name);
    }
  else
    {
      sprintf (filename,"%s\\%s",path,name) ;
    }
  return ( (void *) filename ) ;
}
/*======================================================================*/
FILE * fopen_path ( char *path , char *name , char *mode )
{
  char *filename = malloc_filename(path,name) ;
  FILE *rval = fopen (filename,mode) ;
  if ( rval == NULL )
  {
     fprintf(pcmlogfile,"open of --%s-- failed\n",filename) ;
     fprintf(pcmlogfile,"Error %d %s\n",errno,strerror(errno));
     fclose(pcmlogfile);
     exit(0);
  }
  free (filename) ;
  return ( rval ) ;
}

/* =============================================  */
int FetchRecord(FILE *fp, char *buffer)
{
      int ch;
      char *ptr;
      
      ptr = buffer;
      do 
      {
         ch = getc(fp);
         if (ch >= ' ')
            *ptr++ = ch;
         else if (ch == '\n')
         {
            *ptr = 0;
            return TRUE;
         } else if (ch == '\r')
         {
            ch = getc(fp);
            if (ch != '\n')
                ungetc(ch,fp);
            *ptr = 0;
            return TRUE;
         } else if (ch == EOF)
         {
             *ptr = 0;
             return (ptr != buffer);
         } else
            *ptr++ = ch;
      } while (ptr < buffer+255);
      
      do 
      {
           ch = getc(fp);
       } while (ch !='\n' && ch != '\r' && ch != EOF);
       if (ch == '\r')
       {
           ch = getc(fp);
           if (ch != '\n')
              ungetc(ch,fp);
        }
        *ptr = 0;
        return TRUE;
} 



        
