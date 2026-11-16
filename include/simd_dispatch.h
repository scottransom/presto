/*
 * Runtime CPU dispatch for hot SIMD kernels (coherentsearch, accelsearch).
 *
 * Builds for a baseline x86-64 CPU (e.g. conda-forge's -march=nocona) also
 * compile their hot kernels for x86-64-v3 (AVX2+FMA) and x86-64-v4 (AVX-512),
 * and an ifunc resolver picks one when the program loads, using
 * presto_simd_level() below.  ifunc limits this to x86-64 ELF (Linux), and
 * builds that already target AVX2 (e.g. -march=native) don't need it.
 * Define PRESTO_NO_DISPATCH (or CS_NO_DISPATCH) to turn it off.
 */
#ifndef PRESTO_SIMD_DISPATCH_H
#define PRESTO_SIMD_DISPATCH_H

#if !defined(PRESTO_NO_DISPATCH) && !defined(CS_NO_DISPATCH) && \
    !defined(__AVX2__) && defined(__x86_64__) && defined(__ELF__) && \
    defined(__GNUC__) && defined(__has_attribute)
#if __has_attribute(ifunc) && __has_attribute(target)
#define PRESTO_SIMD_DISPATCH 1
#endif
#endif

/* The ifunc resolvers run while the program is being relocated, before
 * any sanitizer runtime is set up, so they must not be instrumented. */
#if defined(__has_attribute)
#if __has_attribute(no_sanitize)
#define PRESTO_NO_SANITIZE __attribute__((no_sanitize("address", "undefined")))
#endif
#endif
#ifndef PRESTO_NO_SANITIZE
#define PRESTO_NO_SANITIZE
#endif

/* The kernel version the CPU gets: 4 (x86-64-v4), 3 (x86-64-v3) or 0 (the
 * build's baseline).  v4 needs AVX512_VBMI2 too, i.e. Ice Lake, Sapphire
 * Rapids, Zen 4 or newer: on Skylake-SP and Cascade Lake, 512-bit vectors
 * lower the clock so much that the v4 kernels are slower than v3, and even
 * than the baseline (measured: coherentsearch 0.94x the baseline's speed on
 * a Xeon Silver 4114, against 1.35x on a Sapphire Rapids w5-3433). */
PRESTO_NO_SANITIZE static inline int presto_simd_level(void)
{
#ifdef PRESTO_SIMD_DISPATCH
    __builtin_cpu_init();
    if (__builtin_cpu_supports("x86-64-v4") && __builtin_cpu_supports("avx512vbmi2"))
        return 4;
    if (__builtin_cpu_supports("x86-64-v3"))
        return 3;
#endif
    return 0;
}

#endif
