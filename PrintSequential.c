#include "PrintSequential.h"

#include <stdio.h>
#include <string.h>

#include "petsc.h"


namespace PrintSequential
{
  std::string vformat(int & err, int & size, char const * format, va_list argList)
  {
    err = false;
    size = 0;

    // Try a allocation of the stack first, for speed.
    int buffSize = 128;
    char buff[buffSize];

    int numChars = vsnprintf(buff, buffSize-1, format, argList);

    if (numChars == -1)
    {
      printf("Error in printSequential, the format is invalid.\n");
      exit(1);
    }

    if (numChars < buffSize)
      return std::string(buff);

    err = true;
    size = numChars;
    return std::string("");
  };


  std::string vformat(int & numChars, char const * format, va_list argList)
  {
    // Need to allocate enough memory to hold resulting string.
    ++numChars;
    char * newBuff = new char[numChars];

    numChars = vsnprintf(newBuff, numChars, format, argList);

    if (numChars == -1)
    {
      printf("Error in printSequential, the format is invalid.\n");
      exit(1);
    }

    std::string result = std::string(newBuff);

    delete [] newBuff;

    return result;  

};
}

void printSequential(char const * mess)
{
  int leng = strlen(mess);
  // If string is not empty, need to add one to the length
  // so the ending 0x00 byte is transfered.
  if (leng)
    ++leng;
      
  int rank, worldSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, & worldSize);

  int procArrayLocal[worldSize];
  int procArrayGlobal[worldSize];

  for (int i=0; i<worldSize; ++i)
    procArrayLocal[i] = 0;
  procArrayLocal[rank] = leng;


  MPI_Reduce(procArrayLocal, procArrayGlobal, worldSize,
             MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    // rank zero node goes through the procArrayGlobal and requests
    // the data to be printed for each processor with a non-zero
    // procArrayGlobal.

    static char const * localMess = "Rank %d: %s";

    if (leng > 0)
      printf( localMess, rank, mess);

    // Now find the maximum length and make a char variable that big.
    int maxLeng = 0;
    for (int proc=1; proc<worldSize; ++proc)
      if (procArrayGlobal[proc] > maxLeng)
        maxLeng = procArrayGlobal[proc];
    if (maxLeng == 0)
      return;

    char maxMess[maxLeng];

    for (int proc=1; proc<worldSize; ++proc)
    {
      if (procArrayGlobal[proc] > 0)
      {
        MPI_Recv(maxMess, procArrayGlobal[proc], MPI_BYTE,
                 proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf( localMess, proc, maxMess);
      }
    }
  }
  else
  {
    // If this node wants a message printed, respond to the MPI_Reci
    // call otherwise your done.
    if (leng > 0)
      MPI_Send(mess, leng, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  }
}


void printSequentialv(char const * format, ...)
{
  if (strlen(format) == 0)
  {
    printSequential(format);
    return;
  }

  va_list argList;
  va_start(argList, format);
  int err, size;
  std::string mess = PrintSequential::vformat(err, size, format, argList);
  va_end(argList);

  if (!err)
  {
    printSequential(mess.c_str());
    return;
  }

  va_start(argList, format);
  mess = PrintSequential::vformat(size, format, argList);
  va_end(argList);
  printSequential(mess.c_str());
}
