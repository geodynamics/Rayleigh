#include <sys/stat.h>

// Simple directory creation routine

int cmkdir_(char *cc, int ll)
{
    int status;
   cc[ll--] = '\0';  // NULL terminate the string

   //printf("Creating Directory: %s\n",cc);
   status = mkdir(cc,0755); // 755 sets permissions on the directory

   return(status);
}

