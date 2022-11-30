#ifndef _H_TABULATOR
#define _H_TABULATOR

/*! \mainpage tabulator
 * \ref tabulator.h
 * 
 * File "tabulator.h" is Copyright (c) 2008 DeLano Scientific LLC,
 * Palo Alto, California, USA.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of DeLano Scientific LLC nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY DELANO SCIENTIFIC LLC ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL DELANO SCIENTIFIC LLC BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*! \file tabulator.h
 * \brief Tabulator: reader/writer utility for simple text tables.

 * Tabulator manages formatted whitespace-separated ASCII text tables
 * for use as tag records in MDL format SD files (and elsewhere).  The
 * source code is released under the BSD open-source license.

 - \ref behaviors
 - \ref example
 - \ref options
 - \ref conventions

 * \section behaviors Instance Behaviors
 
 * Tabulator instances simulate a two dimensional array of character
 * strings and thus have a type char***.  They can be accessing using
 * C array syntax or C pointer syntax.  When using pointer syntax, you
 * can trust that the pointer arrays are null terminated as diagrammed
 * below:

 \verbatim  
 "0 0" "0 1" "0 2" NULL
 "1 0" "1 1" "1 2" NULL
 NULL \endverbatim

 * such that you can dump the table contents from the tabulator "tab"
 * with the following nested loop:

\verbatim
{
  char **col, ***row = tab;
  while(col = *(row++)) {
    while(*col) {
      printf("%s ",*(col++))
    }
    printf("\n");
  }
} \endverbatim

 * Column header tags are stored in row -1, so you display the tags as follows:
 \verbatim
{
  char **col_tag = tab[-1];
  while(*col_tag) {
    printf("%s ",*(col_tag++));
  }
  printf("\n");
}\endverbatim

 * \section example Example Usage
 * - \ref writing
 * - \ref reading
 * - \ref table

 * \subsection writing Writing a Table
 * \verbatim 
{
   char ***tab = tabulator_new_from_header("AT_ID1 AT_ID2 DISTANCE COLOR DRAW_TYPE", 0);

   tabulator_add_row(&tab, "31 12 1 red 6");
   tabulator_add_row(&tab, "3 8 2 blue 7");
   tabulator_add_row(&tab, "1 2 2 red 3");

   {
     char *buffer = tabulator_as_table(tab);
     printf("%s", buffer);
   }

   tabulator_free(tab);
 } \endverbatim

 * \subsection reading Reading a Table

 * Assuming that "buffer" points at a char* text table in memory...

 \verbatim 
{
  char ***tab2 = tabulator_new_from_table_using_header(buffer, "COLOR AT_ID2 DRAW_TYPE", 0);

  if(tab2) {
    int n_row = tabulator_height(tab2);
    int i;
    for(i=0;i<n_row;i++) {
      printf("Row %d: COLOR=%s AT_ID2=%s DRAW_TYPE=%s\n",i, tab2[i][0], tab2[i][1], tab2[i][2]);
    }
    tabulator_free(tab2);
  }
} \endverbatim

 * \subsection table Sample Table

 * Note that column alignment is purely for benefit of human
 * perception.  The logical structure of the file is solely determined
 * by tokens, whitespace, double-quotes, and newlines.  Note that with the current implementation,
 * double-quotes are not supported characters within tokens (that can be fixed).

 \verbatim
 + AT_ID1 AT_ID2 DISTANCE COLOR DRAW_TYPE 
 |     31     12        1   red         6
 |      3      8        2    ""         7
 |      1      2        2  blue         3
 \endverbatim
 * \section options Optional Preprocessor Defines

 * Must define before #include "tabulator.h".
 
 * \verbatim #define TABULATOR_INCLUDE_IMPLEMENTATION \endverbatim
 * Activates inclusion of the implementation code.
 
 * To use tabulator, one and only .c file in your application must define
 * TABULATOR_INCLUDE_IMPLEMENTATION before including "tabulator.h".

 * \verbatim #define TABULATOR_INCLUDE_UNIT_TEST \endverbatim
 * Activates inclusion "main" with the tabulator unit testing code.  

 * For example, the tabulator unit test can be built by compiling the following main.c:
\verbatim
#define TABULATOR_INCLUDE_IMPLEMENTATION
#define TABULATOR_INCLUDE_UNIT_TEST

#include "tabulator.h"
\endverbatim

 * \section conventions Coding Conventions

 * Tabulator source code conventions as per DeLano Scientific
 * "lazy C" coding style:
 *
 * - ALL_CAPS_UNDERSCORE constants, macros, and defines.
 * 
 * - lowercase_underscore symbol names throughout.
 * 
 * - within member functions, the local pointer "I" used to refer to
 *   the instance (synonymous with C++ "this", or Python "self").
 *   \verbatim I->attribute = value; \endverbatim

*/

#include <stdio.h>

/*! \def TABULATOR_FLAG_STRICT
 * Constructor flag dictating failure if the input text table does not
 * already contain each of the requested header tags.
 */
#define TABULATOR_FLAG_STRICT          0x1

/*! \def TABULATOR_FLAG_DEBUG_STDERR
 * Constructor flag dictating that tabulator should export debug
 * information to stderr.
 */
#define TABULATOR_FLAG_DEBUG_STDERR    0x2

/*!
 * Constructs an new empty tabulator instance with structure matching 
 * the provided column header tags.

 * \param table should point at the first character of the text table (the '+' token).

 * \param header should point at a char* string containing a
 * whitespace separated list of column tags.

 * \return a char*** table whose members can be read directly using C
 * array syntax.
 */
char ***tabulator_new_from_header(char *header, int flags);

/*!
 * Constructs an new empty tabulator instance with the given size.

 * \return a char*** table whose members can be read directly using C
 * array syntax and written using \ref tabulator_set
 */
char ***tabulator_new_with_size(int height, int width, int flags);

/*!
 * Constructs a new tabulator instance from an existing table record.
 * \return a char*** table whose members can be read directly using
 * C array syntax.
 */
char ***tabulator_new_from_table(char *table, int flags);

/*!
 * Constructs a new tabulator instance from a stream.  
 * \return a char*** table whose members can be read directly using
 * C array syntax.
 */
char ***tabulator_new_from_file(FILE *input, int flags);

/*!
 * Constructs a new tabulator instance from an existing table with a
 * structure matching the headers provided.

 * \param table should point at the first character of the text table
 * (the '+' token).

 * \param flags OR'd bitmask of tabulator constructor flags.

 * \return a char*** table whose members can be read directly using C
 * array syntax.

 * When TABULATOR_FLAG_STRICT is specified, then this constructor will
 * return NULL if any of the requested column tags are missing from
 * the source table.
 */
char ***tabulator_new_from_table_using_header(char *table, char *header, int flags);

/*!
 * Constructs a new tabulator instance from a stream with a structure
 * matching the headers provided.

 * \param input should point at the first character of the table in
 * the file (the '+' token).

 * \param flags OR'd bitmask of tabulator constructor flags.

 * \return a char*** table whose members can be read directly using C
 * array syntax.

 * When TABULATOR_FLAG_STRICT is specified, then this constructor will
 * return NULL if any of the requested column tags are missing from
 * the source table.
 */
char ***tabulator_new_from_file_using_header(FILE *input, char *header, int flags);

/*! 
 * Destroys a tabulator instance.
 */
void tabulator_free(char ***tab);

/*! 
 * Returns the number of rows in the table, not including the header row.
 
 * \return table height or -1 on error.
 */
int tabulator_height(char ***tab);

/*! 
 * Returns the number of columns in the table.

 * \return table width or -1 on error.
 */
int tabulator_width(char ***tab);

/*! 
 * Copies an existing tabulator instance to a new instance with a
 * structure matching the provided column header tags.

 * When TABULATOR_FLAG_STRICT is specified, then this constructor will
 * return NULL if any of the requested column tags are missing from
 * the source table.

 */
char ***tabulator_copy_using_header(char ***tab_ptr, char *header, int flags);

/*! 
 * Extends an existing table with a new row consisting of the 
 * whitespace-delimited tokens provided.
 
 *\param tab_ptr is a char**** pointer which WILL be modified as the
 * expanded tabulator instance is relocated within the heap.

 *\param line is a whitespace separated list of tokens to be used in
 *the new row.
 
 * If the number of tokens is less than the table width, the remaining
 * columns are populated with empty tokens ("").

 */
int tabulator_add_row(char ****tab_ptr, char *line);

/*! 
 * Returns a pointer to a character buffer which will remain valid until 
 * the tabulator instance is destroyed or the call is repeated
 */
char   *tabulator_as_table(char ***tab);

/*! 
 * Replaces an existing entry in the tabulator instance. 

 *\returns true if successful, false if not

 * This routine does not release the memory associated with the previous entry.

 */
int tabulator_set(char ***tab, int row, int col, char *str);

#endif

#ifdef TABULATOR_INCLUDE_UNIT_TEST
#include <stdlib.h>
#include <unistd.h>

static int alloc_cnt = 0;
static int alloc_max = 0;
static void *w_malloc(size_t size) 
{ 
  alloc_cnt++;
  if(alloc_max<alloc_cnt) alloc_max = alloc_cnt;
  return malloc(size); 
}
static void *w_calloc(size_t count, size_t size)
{ 
  alloc_cnt++;
  if(alloc_max<alloc_cnt) alloc_max = alloc_cnt;
  return calloc(count,size); 
}
static void w_free(void *ptr) 
{ 
  alloc_cnt--;
  free(ptr);
}
#else
#define w_malloc malloc
#define w_calloc calloc
#define w_free free
#endif

#ifdef TABULATOR_INCLUDE_IMPLEMENTATION

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#ifndef NULL
#define NULL ((void*)0)
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

typedef struct {
  int flags;
  FILE *log, *debug;
  char *output, *text, **col;
  int text_size, text_max;
  int n_col, col_max;
  int n_row, row_max;
  char **row[2]; /* "row" must be last field in struct */
} tabulator;

#define TABULATOR_TABLE_TOKEN '+'
#define TABULATOR_LINE_TOKEN '|'


static void tabulator_copy_and_unquote(char **dst_ptr, char **src_ptr)
{
  char *src = *src_ptr;
  char *dst = *dst_ptr;
  int quote = false;
  int esc = false;
  char ch,m1=0,m2=0,m3=0;
  
  if(*src=='"') {
    quote = true;
    src++;
  }
  while(*src) {
    ch = *(src++);
    if(esc) {
      switch(ch) {
      case '\\':
        *(dst++)=ch;
        esc = false;
        break;
      case '"':
        *(dst++)=ch;
        esc = false;
        break;
      case 'n':
        *(dst++)='\n';
        esc = false;
        break;
      default:
        /* encoded characters */
        m3=m2;
        m2=m1;
        m1=ch;
        if(m3) { /* octal */
          *(dst++)=(m3-'0')*64 + (m2-'0')*8 + (m1-'0');
          esc = false;
        }
        break;
      }
    } else if(quote) {
      if(ch=='"') {
        quote = false;
      } else if(ch=='\\') {
        esc = true; m1=0;m2=0;m3=0;
      } else {
        *(dst++)=ch;
      }
    } else {
      if(ch=='"') {
        quote = true;
      } else if(ch=='\\') {
        esc = true; m1=0;m2=0;m3=0;
      } else {
        if(!((ch<33)||(ch>126)))
          *(dst++)=ch;
        else 
          break;
      }
    }
  }
  *(dst++) = 0;
  *src_ptr = src;
  *dst_ptr = dst;
}

static int tabulator_char_is_white(char ch)
{ 
  return ch&&((ch<33)||(ch>126));    
}

static int tabulator_char_is_newline(char ch) { return (ch==10)||(ch==13); }

char ***tabulator_new_from_header(char *header, int flags)
{
  int ok = true;
  tabulator *I = (tabulator*)w_calloc(sizeof(tabulator),1);
  
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "new_from_header(\"%s\",0x%x): I^%d",
            header, flags, I&&I);
  }
  
  if(I && header) {
    I->flags = flags;
    I->debug = (I->flags & TABULATOR_FLAG_DEBUG_STDERR) ? stderr : NULL;
    I->text_max = strlen(header) + 3;
    I->text_size = I->text_max;
    
    if( (ok = ok && (I->text = (char*)w_calloc(I->text_max , 1))) ) {
      
      I->row_max = 2;
      
      /* count header strings */
      { 
        int n_col = 0;
        char *src = header;
        char *dst = I->text + 1; 
        while(1) { /* pack header into a sequence of non-null strings */
          while(tabulator_char_is_white(*src)) src++;
          if(!*src) {
            break;
          }
          n_col++;
          tabulator_copy_and_unquote(&dst,&src);
        }
        I->n_col = n_col;
      }

      if(I->debug) {
        fprintf(I->debug, ", I->n_col=%d",I->n_col);
      }
      
      /* copy header strings */

      I->col_max = I->n_col + 1;
      if( (ok = ok && (I->col = (char**)w_calloc(sizeof(char*), I->col_max))) ) {
        
        char *src = header;
        char *dst = I->text + 1; 
        char **col = I->col;

        I->row[0] = col; /* becomes c2s[-1] */
        
        while(1) { /* pack header into a sequence of non-null strings */
          while(tabulator_char_is_white(*src)) src++;
          if(!*src) {
            break;
          }
          *(col++) = dst;
          tabulator_copy_and_unquote(&dst,&src);
        }
      }
    }
  }
  
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr,"\n");
  }
  if(!ok) {
    if(I) tabulator_free(I->row + 1);
    return NULL;
  } else {
    return (I->row + 1); /* point at the second entry in I->row */
  }
}


static char *tabulator_copy_line(char *dst, char *src, int buffer_size)
{
  if(buffer_size>0) {
    while(*src && (!tabulator_char_is_newline(*src))) {
      if(buffer_size--) {
        *(dst++) = *(src++);
      } else {
        break;
      }
    }
    *dst = 0;
    if(((*src)==10)&&((*src)==13)) src++; /* CRLF */
    src++; /* skip CR or LF */
  }
  return src;
}

static int tabulator_grow_text(tabulator *I, int chunk)
{
  int new_text_size = I->text_size + chunk;
  if(new_text_size > I->text_max) {
    /* grow text storage as necessary */
    int new_text_max = (new_text_size + (new_text_size>>1));
    char *new_text = (char*)realloc(I->text, new_text_max + 1);
    
    if(new_text) {
      char **col_ptr = I->col;
      char **col_stop = I->col + (I->n_col+1) * (I->n_row+1);
      ptrdiff_t diff = new_text - I->text;
      if(diff) {
        while(col_ptr != col_stop) {
          (*col_ptr) = *(col_ptr) ? (*col_ptr + diff) : (*col_ptr);
          col_ptr++;
        }
      }
      I->text = new_text;
      memset(I->text + I->text_max, 0, new_text_max - I->text_max);
      I->text_max = new_text_max;
    } else {
      return false;
    }
  }
  return true;
}

static int tabulator_grow_col(tabulator *I, int chunk)
{
  int new_col_size = (I->n_col+1) * (I->n_row + 1 + chunk);
  if(new_col_size > I->col_max) {
    int new_col_max = (new_col_size + (new_col_size>>1));
    char **new_col = (char**)realloc(I->col, sizeof(char*) * (new_col_max + 1));
    if(new_col) {
      char ***row = I->row;
      ptrdiff_t diff = new_col - I->col;
      if(diff) {
        while(*row) { *row = *row + diff; row++; }
      }
      I->col = new_col;
      memset(I->col + I->col_max, 0, sizeof(char**) * (new_col_max - I->col_max));
      I->col_max = new_col_max;
    } else {
      return false;
    }
  }
  return true;
}


char ***tabulator_new_with_size(int height, int width, int flags)
{
  int ok = true;
  tabulator *I = (tabulator*)w_calloc(sizeof(tabulator)+(sizeof(char**)*height),1);
  
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "new_with_size(%d,%d): I^%d\n",
            height,width, I&&I);
  }
  
  if(I) {
    I->flags = flags;
    I->debug = (I->flags & TABULATOR_FLAG_DEBUG_STDERR) ? stderr : NULL;
    I->text_max = 3;
    I->text_size = I->text_max;
    
    if( (ok = ok && (I->text = (char*)w_calloc(I->text_max, 1))) ) {
      
      I->n_row = height;
      I->n_col = width;
      I->row_max = I->n_row + 3;
      I->col_max = (width + 1) * (I->n_row + 2);

      if( (ok = ok && (I->col = (char**)w_calloc(sizeof(char*), I->col_max)))) {
        int i;
        for(i=0;i<=height;i++) {
          int j;
          char **col = I->col + i*(width+1);
          I->row[i] = col;
          for(j=0;j<width;j++) {
            *(col++) = I->text;
          }
        }
      }
    }
  }
  if(!ok) {
    if(I) tabulator_free(I->row + 1);
    return NULL;
  } else {
    return (I->row + 1); /* point at the second entry in I->row */
  }
}

char ***tabulator_copy_using_header(char ***tab, char *header, int flags)
{
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "copy_using_header(^%d,\"%s\",0x%x):\n",
            tab&&1, header, flags);
  }
  {
    char ***result = tabulator_new_from_header(header, flags);
    if(result) {
      int ok = true;
      tabulator *I = ((tabulator*)(result+1))-1;
      tabulator *S = ((tabulator*)(tab+1))-1;
      int *xref = (int*)w_malloc(sizeof(int)*I->n_col);
      int *used = (int*)w_calloc(sizeof(int),S->n_col);
      
      if(!(xref && used)) {
        ok = false;
      } else {

        memset(xref,-1,sizeof(int)*I->n_col);

        /* allocate row pointer storage */
        {
          tabulator *new_I = (tabulator*)realloc(I, sizeof(tabulator) +
                                                 sizeof(char***) * S->n_row+1);
          if(new_I) {
            I = new_I;
            result = (I->row + 1);

            I->row_max = S->n_row+1;
            I->n_row = S->n_row;
            memset(I->row+2, 0, sizeof(char**)*I->n_row);

          } else 
            ok = false;
        }
    
        /* allocate col pointer storage */

        if(ok) ok = tabulator_grow_col(I, 0);

        /* create the header cross-reference table */
        if(ok) { 
          int i,s;
          char **i_col = I->row[0];
          char **s_col = S->row[0];
          for(i=0; i<I->n_col; i++) {
            for(s=0; s<S->n_col; s++) {
              if(!used[s]) {
                if( !strcmp(i_col[i],s_col[s]) ) {
                  xref[i] = s;
                  used[s] = true;
                }
              }
            }
          }
        }
      }

      if(ok) {
        if(flags & TABULATOR_FLAG_STRICT) {
          int i;
          for(i=0; i<I->n_col; i++) {
            if(xref[i]<0) {
              ok = false;
              break;
            }
          }
        }
      }

      /* iterate through rows and cols */

      if(ok) {
        int row_idx = 1;
        char ***src_row = S->row + 1;
        char ***dst_row = I->row + 1;
      
        while(*src_row) {
          char **src_col = *(src_row++);
          char **dst_col =  I->col + (I->n_col+1) * (row_idx++);
          int i;
          *(dst_row++) = dst_col;
        
          for(i=0;i<I->n_col;i++) {
            if(xref[i]<0)
              dst_col[i] = I->text; /* missing / blank entry */
            else {
              char *st = src_col[xref[i]];
              int len = strlen(st) + 1;
              if(tabulator_grow_text(I,len)) {
                memcpy( (dst_col[i] = I->text + I->text_size), st, len);
                I->text_size += len;
              } else 
                ok = false;
            }
          }
        }
      }
      w_free(xref);
      w_free(used);
      if(result && !ok) {
        tabulator_free(result);
        result = NULL;
      }
    }
    return result;
  }
}

char ***tabulator_new_from_table(char *table, int flags)
{
  char ***result = NULL;

  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "new_from_table(^%d,0x%x):\n",
            table&&1, flags);
  }

  if(table && (*table == TABULATOR_TABLE_TOKEN)) { 
    int buf_size = 0;

    /* determine size of table (for allocating scratch space) */
    {
      char *src = table; 
      while(*src) {
        while(*src && (!tabulator_char_is_newline(*src))) src++; /* seek end of line */
        if(((*src)==10)&&((*src)==13)) src++; /* CRLF */
        src++;
        if(*src != TABULATOR_LINE_TOKEN) /* new line without start token ends table */
          break;
      }
      buf_size = 1 + (src - table);
    }

    /* load the table */
    if(buf_size) {
      char *scratch = (char*)w_malloc(buf_size);
      if(scratch) {
        char *src = table; 
        if(*src && (*src) == TABULATOR_TABLE_TOKEN) {
          src = tabulator_copy_line(scratch, src+1, buf_size);
          result = tabulator_new_from_header(scratch, flags);
        }
        if(result) {
          while(*src && (*src == TABULATOR_LINE_TOKEN)) {
            src = tabulator_copy_line(scratch, src+1, buf_size);
            tabulator_add_row(&result, scratch);
          }
        }
        w_free(scratch);
      }
    }
  }
  return result;
}

/*!
 * Constructs a new tabulator instance from a table on a stream.
 * \return a char*** table whose members can be read directly using
 * C array syntax.
 */
char ***tabulator_new_from_file(FILE *input, int flags)
{
  char ***result = NULL;
  int buffer_size = 4000;
  int bytes_free = buffer_size;
  char *buffer = (char*)malloc(buffer_size);
  char *newline = NULL;
  char *last = buffer;
  if(buffer) {
    buffer[0] = 0;
    while(1) {

      if(!fgets(last,bytes_free,input)) {
        /* EOF or no characters read -> incomplete/invalid table */
        buffer[0] = 0;
        break;
      }
      if(newline) {

        /* check for a blank line (last would be pointing at a newline) */

        while(*last) {   
          /* skip accidental/invisible whitespace */
          if((*last!=10)&&(*last!=13)&&(tabulator_char_is_white(*last)))
            last++;
          else
            break;
        }
        
        if( ((last[0]==10) && (!last[1]) && (newline[0]==10)) ||  /* LF LF (Unix) */
            ((last[0]==13) && (!last[1]) && (newline[0]==13)) ||  /* CR CR (Mac) */
            ((last[0]==13) && (last[1]==10) && (!last[2]) && 
             (newline[0]==13) && (newline[1]==10)) ) { /* CRLF CRLF (Win) */
          /* found blank line */
          break;
        }
      }

      /* otherwise, scan to end of string */

      while(*last) { 
        last++;
        bytes_free--;
      }

      /* and fine the newline (if present) */

      newline = NULL;
      switch(last-buffer) {
      case 0:
        break;
      case 1:
        if(last[-1]==10)
           newline = last-1;
        break;
      default:
        if((last[-2]==13) && last[-1]==10)
          newline = last-2;
        else if(last[-1]==10)
          newline = last-1;
        break;
      }
     
      /* get us more space if we need it, and update our variables */

      if(bytes_free<2) {
        char *new_buffer;
        int new_size =  buffer_size + (buffer_size>>1);

        new_buffer = (char*)realloc(buffer, new_size);
        if(!new_buffer) {
          buffer[0] = 0; /* error */
          break;
        }

        bytes_free += new_size - buffer_size;
        last = new_buffer + (last - buffer);
        if(newline) newline = new_buffer + (newline - buffer);

        buffer_size = new_size;
        buffer = new_buffer;
      }
    }

    if(buffer[0]) { /* read a table */
      result = tabulator_new_from_table(buffer,flags);
    }
    free(buffer);
  }
  return result;
}

char ***tabulator_new_from_file_using_header(FILE *input, char *header, int flags)
{
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "new_from_file_using_header(^%d,\"%s\",0x%x):\n",
            input&&1, header, flags);
  }
  {
    char ***result = NULL;
    char ***tmp = tabulator_new_from_file(input, flags);
    
    if(tmp) {
      result = tabulator_copy_using_header(tmp, header, flags);
      tabulator_free(tmp);
    }
    return result;
  }
}

char ***tabulator_new_from_table_using_header(char *table, char *header, int flags)
{
  if(flags & TABULATOR_FLAG_DEBUG_STDERR) {
    fprintf(stderr, "new_from_table_using_header(^%d,\"%s\",0x%x):\n",
            table&&1, header, flags);
  }
  {
    char ***result = NULL;
    char ***tmp = tabulator_new_from_table(table, flags);
    
    if(tmp) {
      result = tabulator_copy_using_header(tmp, header, flags);
      tabulator_free(tmp);
    }
    return result;
  }
}

int tabulator_height(char ***tab)
{
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;
    return I->n_row;
  } else {
    return -1;
  }
}

int tabulator_width(char ***tab)
{
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;
    return I->n_col;
  } else {
    return -1;
  }
}

void tabulator_free(char ***tab)
{
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;

    if(I->debug) {
      fprintf(I->debug, "free: I->text^%d, I->col^%d, I->output^%d\n",
              I->text&&1, I->col&&1, I->output&&1);
    }
    if(I->text) w_free(I->text);
    if(I->col) w_free(I->col);
    if(I->output) w_free(I->output);
    w_free(I);
  }
}

#if 0
static void tabulator_debug_dump(char ***tab) 
{
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;
    fprintf(stderr,"dump: I->n_col=%d, I->n_row=%d \n",I->n_col, I->n_row);
    if(I->n_col) {
      char ***row = I->row;
      while(*row) {
        char **col = *(row++);
        fprintf(stderr,"dump: ");
        while(*col) {
          fprintf(stderr,"[%s] ",*(col++));
        }
        fprintf(stderr,"\n");
      }
    }
  } else {
    fprintf(stderr,"dump: null tabulator\n");
  }
}
#endif

int tabulator_add_row(char ****tab_ptr, char *line)
{
  if(*tab_ptr) {
    tabulator *I = ((tabulator*)((*tab_ptr)+1))-1;
    int line_len = strlen(line);
    int ok = true;

    if(I->debug) {
      fprintf(I->debug, "add_row(^1,\"%s\"):\n",line);
    }

    ok = tabulator_grow_text(I, line_len + 1);

    if(ok) ok = tabulator_grow_col(I, 1);

    if(ok) {

      /* extend row storage as necessary */
      int new_row_size = (I->n_row + 3); 
      if(new_row_size > I->row_max) {      
        int new_row_max = (new_row_size + (new_row_size>>1));
        tabulator *new_I = (tabulator*)realloc(I, sizeof(tabulator) 
                                               + sizeof(char***) * new_row_max);
        if(new_I) {
          I = new_I;
          memset(I->row + I->row_max, 0, sizeof(char**) * (new_row_max - I->row_max));
          I->row_max = new_row_max;
        } else {
          ok = false;
        }
      }

      if(ok) {

        /* tokenize and link */
        char *new_text = I->text + I->text_size;
        char **col_start = I->col + (I->n_col+1) * (I->n_row+1);
        int col_cnt = I->n_col; 
        I->row[ 1 + I->n_row++ ] = col_start;

        /* copy the source text */
        memcpy(new_text, line, line_len);
        I->text_size += line_len + 1;

        /* link all n_columns to text token (handles blanks) */
        while(col_cnt--) {

          /* eliminate leading whitespace */
          while(*new_text && tabulator_char_is_white(*new_text)) new_text++; 

          if(!*new_text) {
            *(col_start++) = I->text; /* or substitute a blank string */
          } else {
            char *tmp = new_text;
            *(col_start++) = new_text; /* copy the token */
            tabulator_copy_and_unquote(&tmp, &new_text);
          }
        }
      }
    }
    (*tab_ptr) = ((char***)(I+1))-1;
    return ok;
  } else {
    return false;
  }
}

static int tabulator_escaped_strlen(char *src,int *quotes)
{
  if(!src[0]) {
    return 2;
  } else {
    int base_len = strlen(src);
    char ch;
    int quotes_required = false;
    int extra_len = 0;
      
    while(*src) {
      ch = *(src++);
      if(tabulator_char_is_white(ch)) {
        quotes_required = true;
        switch(ch) {
        case ' ': /* space */
          break;
        case '\n': /* newline */
          extra_len++;
          break;
        default:
          extra_len +=3; /* will use \xxx */
          break;
        }
      } else {
        switch(ch) {
        case '"':
        case '\\':
          extra_len++;
          break;
        }
      }
    }
    if(quotes_required) extra_len += 2;
    if(quotes) {
      *quotes = quotes_required;
    }
    return base_len + extra_len;
  }
}

static void tabulator_escaped_copy(char *dst, int width, char *src)
{
  int need_quotes = false;
  int len = tabulator_escaped_strlen(src,&need_quotes);
  dst += (width-len);
  if(!src[0]) {
    dst[0]='"';
    dst[1]='"'; 
  } else {
    if(need_quotes) {
      *(dst++) = '"';
    }
    while(*src) {
      unsigned char ch = *(src++);
      if(tabulator_char_is_white(ch)) {
        switch(ch) {
        case ' ': /* space */
          *(dst++)=ch;
          break;
        case '\n': /* newline */
          *(dst++)='\\';
          *(dst++)='n';
          break;
        default:
          *(dst++)='\\';
          *(dst++)=((ch>>6)&0x7)+'0';
          *(dst++)=((ch>>3)&0x7)+'0';
          *(dst++)=((ch>>0)&0x7)+'0';
          break;
        }
      } else {
        switch(ch) {
        case '"':
        case '\\':
          *(dst++)='\\';
          *(dst++)=ch;
          break;
        default:
          *(dst++)=ch;
          break;
        }
      }
    }
    if(need_quotes) {
      *(dst++) = '"';
    }
   }
}

char *tabulator_as_table(char ***tab)
{
  char *result = NULL;
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;
    int *max_width = (int*)w_calloc(sizeof(int),I->n_col);
    
    if(max_width) {
      if(I->debug) {
        fprintf(I->debug, "as_table(^1): I->n_row=%d, I->n_col=%d\n",
                I->n_row, I->n_col);
      }
    
      /* measure field widths */
      {
        char ***row = I->row;
        while(*row) {
          char **col = *(row++);
          int *mw = max_width;
          while(*col) {
            int fld_len = tabulator_escaped_strlen(*col,NULL);
            if(fld_len > *mw) *mw = fld_len;
            mw++;
            col++;
          }
        }
      }

      {
        int output_size = 0;

        /* compute total table size */
        {
          int i, line_width = 0;
          for(i=0;i<I->n_col;i++) {
            line_width += (++max_width[i]);
          }
          output_size = ( (line_width + 2) *  /* +2 for prefix & newline */
                          (I->n_row + 1) /* +1 = for headers */
                          + 1); /* and then +1 for terminal newline */
        }
      
        /* allocate buffer for the output */

        if(I->output) w_free(I->output);
        result = I->output = w_malloc(output_size+1);

        if(result) {
          char *dst = result;

          /* fill with spaces */
          memset(result, 32, output_size);

          /* terminate */
          result[output_size] = 0;

          /* copy fields into the output at proper locations */

          {
            char ***row = I->row; 
            char token = TABULATOR_TABLE_TOKEN;
            while(*row) {
              char **col = *(row++);
              int *mw = max_width;
              *(dst++) = token;
              while(*col) {
                tabulator_escaped_copy(dst, *mw, *(col++));
                dst += *(mw++);
              }
              *(dst++) = '\n'; /* row-ending newline */
              token = TABULATOR_LINE_TOKEN;
            }
          }
          *(dst++) = '\n'; /* table-ending newline */
        }
      }
      w_free(max_width);
    }
  }
  return result;
}

int tabulator_set(char ***tab, int row, int col, char *str)
{
  if(tab) {
    tabulator *I = ((tabulator*)(tab+1))-1;
    if( (row>-2) && (row<I->n_row) && (col>=0) && (col<I->n_col) ) {
      int len = strlen(str) + 1;
      if(tabulator_grow_text(I,len)) {
        char *dst = I->text + I->text_size;
        memcpy(dst, str, len);
        tab[row][col] = dst;
        I->text_size += len;
        return true;
      }
    }
  }
  return false;
}

#endif

#ifdef TABULATOR_INCLUDE_UNIT_TEST
int main(int argc, char **argv)
{

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header("",
                                            TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    tabulator_free(tab);
  }
  
  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header(" ONE ",
                                            TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header("TWO THREE",
                                            TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char buffer[255];
    char ***tab = tabulator_new_from_header(" ONE FIELD_TWO THREE FOUR ",
                                            TABULATOR_FLAG_DEBUG_STDERR);

    tabulator_debug_dump(tab);
    sprintf(buffer,"   asdf  245.6   ");
    tabulator_add_row(&tab, buffer);
    tabulator_debug_dump(tab);

    sprintf(buffer,"-45.6 CA 235 01");
    tabulator_add_row(&tab, buffer);
    tabulator_debug_dump(tab);

    sprintf(buffer,"blah1 blah2 blah3 blah4 blah5 blah6");
    tabulator_add_row(&tab, buffer);
    tabulator_debug_dump(tab);

    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));      
    tabulator_free(tab);    
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    fprintf(stderr, "tab^%d\n",tab&&1);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+\n|\n|\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ TEST",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ TEST TOO",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ A Ab Abc\n| 1\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ A Ab Abc\n| 1 2 3\n|4 5 6\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ Abc Ab A\n| 1 2 3 4 5 6\n|7 8\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ ONE TWO THREE\n| 1 2 3\n| 11 22  33 \n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));

    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_copy_using_header(tab, "",
                                                 TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_copy_using_header(tab, "TWO",
                                                 TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_copy_using_header(tab, "THREE TWO ONE",
                                                 TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }

    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "THREE ONE TWO",
                                                           TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }

    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "THREE FOUR TWO",
                                                           TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
 
    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "THREE FOUR TWO",
                                                           TABULATOR_FLAG_STRICT |
                                                           TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
   tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header(" \"ONE\" \"TWO DAY\" DONT\\\\DO\\\"IT ",
                                            TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_add_row(&tab, "1 2");
    tabulator_debug_dump(tab);
    if(tab) fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));      
    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "\"TWO DAY\" ONE",
                                                           TABULATOR_FLAG_STRICT |
                                                           TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
    tabulator_free(tab);
  }
 
  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header(" \"TWO\\037\" hi\\nyou ",
                                            TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_add_row(&tab, "1 2");
    tabulator_debug_dump(tab);
    if(tab) fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));      
    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           " hi\\nyou \"TWO\\037\"",
                                                           TABULATOR_FLAG_STRICT |
                                                           TABULATOR_FLAG_DEBUG_STDERR);
      tabulator_debug_dump(tab2);      
      if(tab2) fprintf(stderr, "tabulator_as_table(tab2):\n%s",tabulator_as_table(tab2));      
      tabulator_free(tab2);
    }
    tabulator_free(tab);
  }
 
  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ A \"Ab \" Abc\n| 1 2 3\n|4 5 6\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));

    fprintf(stderr, "%d\n", tabulator_set(tab,-1,1,"Abx"));
    fprintf(stderr, "%d\n", tabulator_set(tab,1,1,"xxx"));

    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));

    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_table("+ A \"\" Abc\n| 1\n",
                                           TABULATOR_FLAG_DEBUG_STDERR);
    tabulator_debug_dump(tab);
    fprintf(stderr, "tabulator_as_table(tab):\n%s",tabulator_as_table(tab));
    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header("AT_ID1 AT_ID2 DISTANCE COLOR DRAW_TYPE", 0);
    
    tabulator_add_row(&tab, "31 12 1 red 6");
    tabulator_add_row(&tab, "3 8 2 blue 7");
    tabulator_add_row(&tab, "1 2 2 red 3");
    
    fprintf(stderr,"%s", tabulator_as_table(tab));

    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "COLOR AT_ID2 DRAW_TYPE", 0);
      if(tab) {
        {
          int n_row = tabulator_height(tab2);
          int i;
          for(i=0;i<n_row;i++) {
            fprintf(stderr,"Row %d: COLOR=%s AT_ID2=%s DRAW_TYPE=%s\n",i, tab[i][0], tab[i][1], tab[i][2]);
          }
          tabulator_free(tab2);
        }
      }
    }

    tabulator_free(tab);
  }

  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  {
    char ***tab = tabulator_new_from_header("AT_ID1 AT_ID2 DISTANCE COLOR DRAW_TYPE", 0);
    
    tabulator_add_row(&tab, "31 12 1 red 6");
    tabulator_add_row(&tab, "3 8 2 blue 7");
    tabulator_add_row(&tab, "1 2 2 red 3");
    
    fprintf(stderr,"%s", tabulator_as_table(tab));

    fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
    {
      char ***tab2 = tabulator_new_from_table_using_header(tabulator_as_table(tab),
                                                           "COLOR AT_ID2 DRAW_TYPE", 0);
      if(tab) {
        {
          int n_row = tabulator_height(tab2);
          int i;
          for(i=0;i<n_row;i++) {
            fprintf(stderr,"Row %d: COLOR=%s AT_ID2=%s DRAW_TYPE=%s\n",i, tab[i][0], tab[i][1], tab[i][2]);
          }
          tabulator_free(tab2);
        }
      }
    }

    tabulator_free(tab);
  }


  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  /* confirm the ability to read and write tables of all different
     size, with all sorts of zero-terminated strings */

  {
    int a,cnt=0;
    for(a=0;a<5000;a++) {
      int width = random() & 0x0F;
      int height = random() & 0x3F;
      char *buf = malloc(0x10*width*(height+1));

      char ***tab = tabulator_new_with_size(height,width,0);
      if(buf && tab) {
        int i;
        for(i=-1;i<height;i++) {
          int j;
          for(j=0;j<width;j++) {
            int stlen = (random() & 0x3);
            char *str = buf + (((i+1)*width)+j)*0xF;
            int s;
            for(s=0;s<stlen;s++) {
              str[s] = (random()&0xFE)+1;
            }
            str[stlen] = 0;
            tabulator_set(tab,i,j,str);
          }
        }

        {
          char *tmp = tabulator_as_table(tab);
          char ***tab2 = tabulator_new_from_table(tmp,TABULATOR_FLAG_STRICT);
          char *tmp2 = tabulator_as_table(tab2);
          int valid = !strcmp(tmp,tmp2);
   
          if(valid) {
            for(i=-1;i<height;i++) {
              int j;
              for(j=0;j<width;j++) {
                char *str = buf + (((i+1)*width)+j)*0xF;
                valid = valid && (!strcmp(str,tab2[i][j]));
                tabulator_set(tab,i,j,str);
              }
            }
          }
          
          if(!valid) {
            fprintf(stderr,"INVALID: %dx%d\n",height,width);
            fprintf(stderr, "[\n%s]<\n%s>",tmp,tmp2);
            break;
          } 
          tabulator_free(tab2);
        }
        free(buf);
        tabulator_free(tab);
      }
      cnt++;
    }
    fprintf(stderr,"Confirmed read/write conversion for %d tables.\n",cnt);
  }


  fprintf(stderr,"---------------- TEST-tabulator.h: %d ----------------\n",__LINE__); 
  /* confirm the ability to read and write tables of all different
     size to and from a file with all sorts of zero-terminated strings */

  {
    int a,cnt=0;
    for(a=0;a<12;a++) {
      
      /* first write */
      int width = a;
      int height = a;
      FILE *f = fopen("tabulator.tmp","w");
      if(f) {
        for(width=0;width<a;width++) {
          for(height=0;height<a;height++) {
            char *buf = malloc(0x10*width*(height+1));
            
            char ***tab = tabulator_new_with_size(height,width,0);
            if(buf && tab) {
              int i;
              for(i=-1;i<height;i++) {
                int j;
                for(j=0;j<width;j++) {
                  int stlen = (random() & 0x3);
                  char *str = buf + (((i+1)*width)+j)*0xF;
                  int s;
                  for(s=0;s<stlen;s++) {
                    str[s] = (random()&0xFE)+1;
                  }
                  str[stlen] = 0;
                  tabulator_set(tab,i,j,str);
                }
              }
              {
                char *tmp = tabulator_as_table(tab);
                fwrite(tmp,strlen(tmp),1,f);
              }
              tabulator_free(tab);
            }
            free(buf);
          }
        }
        fclose(f);

        f=fopen("tabulator.tmp","r");
        if(f) {
          for(width=0;width<a;width++) {
            for(height=0;height<a;height++) {
              char ***tab2 = tabulator_new_from_file(f,TABULATOR_FLAG_STRICT);
              
              if(!tab2) {
                fprintf(stderr,"UNABLE TO READ TABLE FROM FILE %dx%d\n",height,width);
              } else if((tabulator_width(tab2)!=width)||
                        (tabulator_height(tab2)!=height)) {
                fprintf(stderr,"INVALID: %dx%d\n",height,width);
              }
              if(tab2) {
                cnt++;
                tabulator_free(tab2);
              }
            }
          }
          if(tabulator_new_from_file(f,TABULATOR_FLAG_STRICT)) { /* should always fail */
            fprintf(stderr,"SHOULD HAVE FAILED\n");
          }
          fclose(f);
        }

        f=fopen("tabulator.tmp","r");
        if(f) {
          for(width=0;width<a;width++) {
            for(height=0;height<a;height++) {
              char ***tab2 = tabulator_new_from_file_using_header(f,"1 2 3 4 5",0);
              
              if(!tab2) {
                fprintf(stderr,"UNABLE TO READ TABLE FROM FILE %dx%d\n",height,width);
              } else if((tabulator_width(tab2)!=5)||
                        (tabulator_height(tab2)!=height)) {
                fprintf(stderr,"INVALID: %dx%d\n",height,width);
              }
              if(tab2) {
                cnt++;
                tabulator_free(tab2);
              }
            }
          }
          if(tabulator_new_from_file(f,TABULATOR_FLAG_STRICT)) { /* should always fail */
            fprintf(stderr,"SHOULD HAVE FAILED\n");
          }
          fclose(f);
        }


      }
    }
    unlink("tabulator.tmp");
    fprintf(stderr,"Confirmed file read/write conversion %d times.\n",cnt);
  }

  fprintf(stderr,"HEAP BLOCKS ALLOCATED: max=%d, now=%d\n",alloc_max, alloc_cnt);
  return 0;
}

#endif




