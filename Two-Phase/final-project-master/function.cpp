// getch_nonblock.c
#include <conio.h>

// Check if a key was pressed
int kbhit_c() {
    return _kbhit();
}

// Get the character (non-blocking)
int getch_c() {
    return _getch();
}
