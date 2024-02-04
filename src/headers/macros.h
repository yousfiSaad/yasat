
#ifndef MACROS__
#define MACROS__



#define FOR(i, lim) \
  for (int i = 0; i < lim; i++)

#define FORC(i, lim) \
  for (char i = 0; i < lim; i++)

#define FORS(i, lim) \
  for (short i = 0; i < lim; i++)

#define FORI8(i, lim) \
  for (int8_t i = 0; i < lim; i++)

#define MAX_SOLUTIONS 2

#define FORRC(i, j, k)         \
  for (char i = 0; i < k; i++) \
    for (char j = 0; j < k; j++)

#define FORRS(i, j, k)         \
  for (short i = 0; i < k; i++) \
    for (short j = 0; j < k; j++)

#define FORRI8(i, j, k)          \
  for (int8_t i = 0; i < k; i++) \
    for (int8_t j = 0; j < k; j++)

#define FORR(i, j, k)         \
  for (int i = 0; i < k; i++) \
    for (int j = 0; j < k; j++)

#endif