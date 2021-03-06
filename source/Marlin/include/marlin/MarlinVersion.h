#ifndef MarlinVersion_h
#define MarlinVersion_h 1

// version macros
#define MARLIN_MAJOR_VERSION 0
#define MARLIN_MINOR_VERSION 11
#define MARLIN_PATCH_LEVEL 0
#define MARLIN_VERSION_GE( MAJV , MINV , PLEV)  ( (  MARLIN_MAJOR_VERSION > MAJV ) || ( (MARLIN_MAJOR_VERSION==MAJV) && ( MARLIN_MINOR_VERSION > MINV ) ) ||   ( (MARLIN_MAJOR_VERSION==MAJV) && ( MARLIN_MINOR_VERSION == MINV ) && ( MARLIN_PATCH_LEVEL >= PLEV ) ) )

#endif
