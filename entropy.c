#include <math.h>
#include <stdint.h>
#include <htscodecs/utils.h>

#include "entropy.h"

double entropy0(unsigned char *buf, uint32_t size) {
    uint32_t F[256] = {0};
    hist8(buf, size, F);

    double e = 0;
    for (int i = 0; i < 256; i++) {
	if (F[i])
	    e += -fast_log((double)F[i]/size) * F[i];
    }

    return e / log(256);
}

double entropy1(unsigned char *buf, uint32_t size) {
    uint32_t F[256][256] = {{0}}, T[256] = {0};
    hist1_4(buf, size, F, T);

    double e = 0;
    for (int i = 0; i < 256; i++) {
	if (!T[i]) continue;
	for (int j = 0; j < 256; j++) {
	    if (F[i][j])
		e += -fast_log((double)F[i][j]/T[i]) * F[i][j];
	}
    }

    return e / log(256);
}
