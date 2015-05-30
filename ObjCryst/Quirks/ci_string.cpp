#include "ci_string.h"
/// Case-insensitive string class
/// From: Guru of the Week #29
/// e.g.: http://gcc.gnu.org/onlinedocs/libstdc++/21_strings/gotw29a.txt
///
/// Public domain.
int strnicmp(const char *s1, const char *s2, int len)
{
   unsigned char c1, c2;
   while (len)
   {
      c1 = *s1; c2 = *s2;
      s1++; s2++;
      if (!c1) return c2 ? -1 : 0;
      if (!c2) return 1;
      if (c1 != c2)
      {
         c1 = tolower(c1);
         c2 = tolower(c2);
         if (c1 != c2) return c1 < c2 ? -1 : 1;
      }
      len--;
   }
   return 0;
}

bool ci_char_traits::eq( char c1, char c2 )
{return tolower(c1) == tolower(c2);}

bool ci_char_traits::ne( char c1, char c2 )
   {return tolower(c1) != tolower(c2);}

bool ci_char_traits::lt( char c1, char c2 )
   {return tolower(c1) < tolower(c2);}

int ci_char_traits::compare(const char* s1,const char* s2,size_t n )
{
   #ifdef _MSC_VER
   return _strnicmp(s1, s2, n);
   #else
	return strnicmp(s1, s2, n);
   #endif
}

const char* ci_char_traits::find( const char* s, int n, char a )
{
   while( n-- > 0 && tolower(*s) != tolower(a) ) ++s;
   return s;
}
