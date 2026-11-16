/*
 * accel_gather.c: the harmonic-sum gather at the heart of accelsearch,
 * add_gathered_row(): out[k] += in[rinds[k]] for k < n.
 *
 * in is a row of a subharmonic's f-fdot plane and rinds maps each
 * fundamental r bin to its subharmonic bin.  The subharmonic's r step is
 * harmnum/numharm < 1 of the fundamental's, so rinds[k+1] - rinds[k] is 0
 * or 1, and any V consecutive outputs read from V consecutive inputs.  So
 * instead of a gather, the SIMD versions load V inputs from in + rinds[k]
 * and permute them into place: AVX2 and AVX-512 (behind the runtime
 * dispatch in simd_dispatch.h on baseline x86-64 builds) and NEON (always
 * there on aarch64).  They may read up to V-1 floats past the last input
 * used, so in must have ACCEL_ROWPAD readable floats after it.  Every
 * version adds the same floats in the same order, so the sums are the
 * same, bit for bit.  tests/accelsearch/gather_test.c checks that.
 */
#include "accel.h"
#include "simd_dispatch.h"

static void add_gathered_row_scalar(float *restrict out, const float *restrict in,
                                    const int *restrict rinds, int n)
{
    for (int kk = 0; kk < n; kk++)
        out[kk] += in[rinds[kk]];
}

#if defined(PRESTO_SIMD_DISPATCH) || defined(__AVX2__)
#include <immintrin.h>

#ifdef PRESTO_SIMD_DISPATCH
#define ACCEL_TARGET_V3 __attribute__((target("arch=x86-64-v3")))
#define ACCEL_TARGET_V4 __attribute__((target("arch=x86-64-v4")))
#else
#define ACCEL_TARGET_V3
#define ACCEL_TARGET_V4
#endif

#define ACCEL_HAVE_AVX2 1
ACCEL_TARGET_V3 static void add_gathered_row_avx2(float *restrict out,
                                                  const float *restrict in,
                                                  const int *restrict rinds, int n)
{
    int kk = 0;
    for (; kk + 8 <= n; kk += 8) {
        const int base = rinds[kk];
        const __m256i idx =
            _mm256_sub_epi32(_mm256_loadu_si256((const __m256i *) (rinds + kk)),
                             _mm256_set1_epi32(base));
        const __m256 vals = _mm256_permutevar8x32_ps(_mm256_loadu_ps(in + base), idx);
        _mm256_storeu_ps(out + kk, _mm256_add_ps(_mm256_loadu_ps(out + kk), vals));
    }
    for (; kk < n; kk++)
        out[kk] += in[rinds[kk]];
}

/* Like the coherentsearch kernels, use 512-bit vectors only with
 * AVX512_VBMI2 (Ice Lake and newer): see presto_simd_level(). */
#if defined(PRESTO_SIMD_DISPATCH) || \
    (defined(__AVX512F__) && defined(__AVX512VBMI2__))
#define ACCEL_HAVE_AVX512 1
ACCEL_TARGET_V4 static void add_gathered_row_avx512(float *restrict out,
                                                    const float *restrict in,
                                                    const int *restrict rinds, int n)
{
    int kk = 0;
    for (; kk + 16 <= n; kk += 16) {
        const int base = rinds[kk];
        const __m512i idx =
            _mm512_sub_epi32(_mm512_loadu_si512((const void *) (rinds + kk)),
                             _mm512_set1_epi32(base));
        const __m512 vals = _mm512_permutexvar_ps(idx, _mm512_loadu_ps(in + base));
        _mm512_storeu_ps(out + kk, _mm512_add_ps(_mm512_loadu_ps(out + kk), vals));
    }
    for (; kk < n; kk++)
        out[kk] += in[rinds[kk]];
}
#endif
#endif

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define ACCEL_HAVE_NEON 1
/* tbl selects bytes, so lane j of the index takes bytes 4 * rel + {0, 1,
 * 2, 3} of the loaded window (rel <= 3, so the multiply-add has no carries). */
static void add_gathered_row_neon(float *restrict out, const float *restrict in,
                                  const int *restrict rinds, int n)
{
    const uint32x4_t bytes = vdupq_n_u32(0x03020100u);
    int kk = 0;
    for (; kk + 4 <= n; kk += 4) {
        const int base = rinds[kk];
        const uint32x4_t rel =
            vreinterpretq_u32_s32(vsubq_s32(vld1q_s32(rinds + kk), vdupq_n_s32(base)));
        const uint8x16_t sel = vreinterpretq_u8_u32(vmlaq_n_u32(bytes, rel, 0x04040404u));
        const float32x4_t vals = vreinterpretq_f32_u8(
            vqtbl1q_u8(vld1q_u8((const uint8_t *) (in + base)), sel));
        vst1q_f32(out + kk, vaddq_f32(vld1q_f32(out + kk), vals));
    }
    for (; kk < n; kk++)
        out[kk] += in[rinds[kk]];
}
#endif

#ifdef PRESTO_SIMD_DISPATCH
PRESTO_NO_SANITIZE static void (*add_gathered_row_resolve(void))
 (float *restrict, const float *restrict, const int *restrict, int)
{
    const int lev = presto_simd_level();
    return lev == 4 ? add_gathered_row_avx512
        : lev == 3 ? add_gathered_row_avx2 : add_gathered_row_scalar;
}

void add_gathered_row(float *restrict out, const float *restrict in,
                      const int *restrict rinds, int n)
    __attribute__((ifunc("add_gathered_row_resolve")));
#else
void add_gathered_row(float *restrict out, const float *restrict in,
                      const int *restrict rinds, int n)
{
#if defined(ACCEL_HAVE_AVX512)
    add_gathered_row_avx512(out, in, rinds, n);
#elif defined(ACCEL_HAVE_AVX2)
    add_gathered_row_avx2(out, in, rinds, n);
#elif defined(ACCEL_HAVE_NEON)
    add_gathered_row_neon(out, in, rinds, n);
#else
    add_gathered_row_scalar(out, in, rinds, n);
#endif
}
#endif

const char *accel_simd_name(void)
/* Which add_gathered_row() this run uses, for -perfinfo */
{
#if defined(PRESTO_SIMD_DISPATCH)
    const int lev = presto_simd_level();
    return lev == 4 ? "x86-64-v4 (AVX-512), chosen at run time"
        : lev == 3 ? "x86-64-v3 (AVX2), chosen at run time"
        : "baseline (scalar), chosen at run time";
#elif defined(ACCEL_HAVE_AVX512)
    return "AVX-512, as compiled";
#elif defined(ACCEL_HAVE_AVX2)
    return "AVX2, as compiled";
#elif defined(ACCEL_HAVE_NEON)
    return "NEON, as compiled";
#else
    return "scalar, as compiled";
#endif
}

int accel_gather_kernels(const char **names, accel_gather_fn *fns, int maxnum)
/* For tests: every version of add_gathered_row() this CPU can run, the */
/* scalar one first.  Returns how many (at most maxnum).                */
{
    int num = 0;
#define ACCEL_ADD_KERNEL(NAME, FN) \
    if (num < maxnum) { names[num] = NAME; fns[num++] = FN; }
    ACCEL_ADD_KERNEL("scalar", add_gathered_row_scalar);
#ifdef PRESTO_SIMD_DISPATCH
    const int lev = presto_simd_level();
    if (lev >= 3)
        ACCEL_ADD_KERNEL("avx2", add_gathered_row_avx2);
    if (lev >= 4)
        ACCEL_ADD_KERNEL("avx512", add_gathered_row_avx512);
#else
#ifdef ACCEL_HAVE_AVX2
    ACCEL_ADD_KERNEL("avx2", add_gathered_row_avx2);
#endif
#ifdef ACCEL_HAVE_AVX512
    ACCEL_ADD_KERNEL("avx512", add_gathered_row_avx512);
#endif
#endif
#ifdef ACCEL_HAVE_NEON
    ACCEL_ADD_KERNEL("neon", add_gathered_row_neon);
#endif
    ACCEL_ADD_KERNEL("dispatched", add_gathered_row);
#undef ACCEL_ADD_KERNEL
    return num;
}
