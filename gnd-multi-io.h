/*
 * gnd-multi-io.h
 *
 * brief  : io header for multi-platform
 * support: gcc, visual studio 10
 */

#ifndef GND_MULTI_IO_H_
#define GND_MULTI_IO_H_


#if defined(__linux__)

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define gnd_multi_open				open
#define gnd_multi_write				write
#define gnd_multi_read				read
#define gnd_multi_close				close



#elif defined(__MINGW32__)

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define gnd_multi_open				open
#define gnd_multi_write				write
#define gnd_multi_read				read
#define gnd_multi_close				close

#define S_IRWXG				S_IRWXU
#define S_IRGRP				S_IRUSR
#define S_IWGRP				S_IWUSR
#define S_IXGRP				S_IXUSR
#define S_IRWXO				0x00
#define S_IROTH				0x00
#define S_IWOTH				0x00
#define S_IXOTH				0x00



#elif defined(_MSC_VER) 		// visual studio cl

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <io.h>

#define gnd_multi_open		_open

// open flags
#define O_APPEND			_O_APPEND
#define O_BINARY			_O_BINARY
#define O_TRUNC				_O_TRUNC
#define O_CREAT				_O_CREAT
#define O_EXCL				_O_EXCL

#define O_WRONLY			_O_WRONLY
#define O_RDONLY			_O_RDONLY
#define O_RDWD				_O_RDWD

#define S_IRWXU				_S_IREAD | _S_IWRITE
#define S_IRUSR				_S_IREAD
#define S_IWUSR				_S_IWRITE
#define S_IXUSR				0x00
#define S_IRWXG				_S_IREAD | _S_IWRITE
#define S_IRGRP				_S_IREAD
#define S_IWGRP				_S_IWRITE
#define S_IXGRP				0x00
#define S_IRWXO				_S_IREAD | _S_IWRITE
#define S_IROTH				_S_IREAD
#define S_IWOTH				_S_IWRITE
#define S_IXOTH				0x00

#define gnd_multi_write		_write
#define gnd_multi_read		_read
#define gnd_multi_close		_close

#endif

#endif
