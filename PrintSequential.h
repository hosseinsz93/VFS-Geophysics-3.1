#ifndef __PRINT_SEQUENTIAL__
#define __PRINT_SEQUENTIAL__

#include <memory>
#include <stdarg.h>

// To use:

// include this file to make PrintSequential available.

// To use, call printSequential(mess); where mess is a pointer to a
// string variable. this includes char mess[size];

// If there is no message to print, call printSequential("");

// Notes: no need to include the rank in your message,
//            printSequential prepends one for you.
//        printSequential doesn't add \n to messages.
//            If you want one, you need to include it in your message.

void printSequential(char const * mess);

// Specify the format and list of arguments just like using printf.
void printSequentialv(char const * format, ...);

#endif
