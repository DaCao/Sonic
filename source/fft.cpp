//   fft.cpp - impelementation of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#include <fstream>
#include <iostream>
#include <math.h>
#include "../include/fft.h"

//#include <complex>
//#include "complex.h"


using namespace std;
//using ::complex;

//Convolution function which returns the frequency representation of the result.
//NFFT is the FFT size, nSIG is the size for the input, NFIL is the size of the filter.



// store all the precalculated factors and indexes in the array
int PermutedIndex [1024] = {0, 512, 256, 768, 128, 640, 384, 896, 64, 576, 320, 832, 192, 704, 448, 960, 32, 544, 288, 800, 160, 672, 416, 928, 96, 608, 352, 864, 224, 736, 480, 992, 16, 528, 272, 784, 144, 656, 400, 912, 80, 592, 336, 848, 208, 720, 464, 976, 48, 560, 304, 816, 176, 688, 432, 944, 112, 624, 368, 880, 240, 752, 496, 1008, 8, 520, 264, 776, 136, 648, 392, 904, 72, 584, 328, 840, 200, 712, 456, 968, 40, 552, 296, 808, 168, 680, 424, 936, 104, 616, 360, 872, 232, 744, 488, 1000, 24, 536, 280, 792, 152, 664, 408, 920, 88, 600, 344, 856, 216, 728, 472, 984, 56, 568, 312, 824, 184, 696, 440, 952, 120, 632, 376, 888, 248, 760, 504, 1016, 4, 516, 260, 772, 132, 644, 388, 900, 68, 580, 324, 836, 196, 708, 452, 964, 36, 548, 292, 804, 164, 676, 420, 932, 100, 612, 356, 868, 228, 740, 484, 996, 20, 532, 276, 788, 148, 660, 404, 916, 84, 596, 340, 852, 212, 724, 468, 980, 52, 564, 308, 820, 180, 692, 436, 948, 116, 628, 372, 884, 244, 756, 500, 1012, 12, 524, 268, 780, 140, 652, 396, 908, 76, 588, 332, 844, 204, 716, 460, 972, 44, 556, 300, 812, 172, 684, 428, 940, 108, 620, 364, 876, 236, 748, 492, 1004, 28, 540, 284, 796, 156, 668, 412, 924, 92, 604, 348, 860, 220, 732, 476, 988, 60, 572, 316, 828, 188, 700, 444, 956, 124, 636, 380, 892, 252, 764, 508, 1020, 2, 514, 258, 770, 130, 642, 386, 898, 66, 578, 322, 834, 194, 706, 450, 962, 34, 546, 290, 802, 162, 674, 418, 930, 98, 610, 354, 866, 226, 738, 482, 994, 18, 530, 274, 786, 146, 658, 402, 914, 82, 594, 338, 850, 210, 722, 466, 978, 50, 562, 306, 818, 178, 690, 434, 946, 114, 626, 370, 882, 242, 754, 498, 1010, 10, 522, 266, 778, 138, 650, 394, 906, 74, 586, 330, 842, 202, 714, 458, 970, 42, 554, 298, 810, 170, 682, 426, 938, 106, 618, 362, 874, 234, 746, 490, 1002, 26, 538, 282, 794, 154, 666, 410, 922, 90, 602, 346, 858, 218, 730, 474, 986, 58, 570, 314, 826, 186, 698, 442, 954, 122, 634, 378, 890, 250, 762, 506, 1018, 6, 518, 262, 774, 134, 646, 390, 902, 70, 582, 326, 838, 198, 710, 454, 966, 38, 550, 294, 806, 166, 678, 422, 934, 102, 614, 358, 870, 230, 742, 486, 998, 22, 534, 278, 790, 150, 662, 406, 918, 86, 598, 342, 854, 214, 726, 470, 982, 54, 566, 310, 822, 182, 694, 438, 950, 118, 630, 374, 886, 246, 758, 502, 1014, 14, 526, 270, 782, 142, 654, 398, 910, 78, 590, 334, 846, 206, 718, 462, 974, 46, 558, 302, 814, 174, 686, 430, 942, 110, 622, 366, 878, 238, 750, 494, 1006, 30, 542, 286, 798, 158, 670, 414, 926, 94, 606, 350, 862, 222, 734, 478, 990, 62, 574, 318, 830, 190, 702, 446, 958, 126, 638, 382, 894, 254, 766, 510, 1022, 1, 513, 257, 769, 129, 641, 385, 897, 65, 577, 321, 833, 193, 705, 449, 961, 33, 545, 289, 801, 161, 673, 417, 929, 97, 609, 353, 865, 225, 737, 481, 993, 17, 529, 273, 785, 145, 657, 401, 913, 81, 593, 337, 849, 209, 721, 465, 977, 49, 561, 305, 817, 177, 689, 433, 945, 113, 625, 369, 881, 241, 753, 497, 1009, 9, 521, 265, 777, 137, 649, 393, 905, 73, 585, 329, 841, 201, 713, 457, 969, 41, 553, 297, 809, 169, 681, 425, 937, 105, 617, 361, 873, 233, 745, 489, 1001, 25, 537, 281, 793, 153, 665, 409, 921, 89, 601, 345, 857, 217, 729, 473, 985, 57, 569, 313, 825, 185, 697, 441, 953, 121, 633, 377, 889, 249, 761, 505, 1017, 5, 517, 261, 773, 133, 645, 389, 901, 69, 581, 325, 837, 197, 709, 453, 965, 37, 549, 293, 805, 165, 677, 421, 933, 101, 613, 357, 869, 229, 741, 485, 997, 21, 533, 277, 789, 149, 661, 405, 917, 85, 597, 341, 853, 213, 725, 469, 981, 53, 565, 309, 821, 181, 693, 437, 949, 117, 629, 373, 885, 245, 757, 501, 1013, 13, 525, 269, 781, 141, 653, 397, 909, 77, 589, 333, 845, 205, 717, 461, 973, 45, 557, 301, 813, 173, 685, 429, 941, 109, 621, 365, 877, 237, 749, 493, 1005, 29, 541, 285, 797, 157, 669, 413, 925, 93, 605, 349, 861, 221, 733, 477, 989, 61, 573, 317, 829, 189, 701, 445, 957, 125, 637, 381, 893, 253, 765, 509, 1021, 3, 515, 259, 771, 131, 643, 387, 899, 67, 579, 323, 835, 195, 707, 451, 963, 35, 547, 291, 803, 163, 675, 419, 931, 99, 611, 355, 867, 227, 739, 483, 995, 19, 531, 275, 787, 147, 659, 403, 915, 83, 595, 339, 851, 211, 723, 467, 979, 51, 563, 307, 819, 179, 691, 435, 947, 115, 627, 371, 883, 243, 755, 499, 1011, 11, 523, 267, 779, 139, 651, 395, 907, 75, 587, 331, 843, 203, 715, 459, 971, 43, 555, 299, 811, 171, 683, 427, 939, 107, 619, 363, 875, 235, 747, 491, 1003, 27, 539, 283, 795, 155, 667, 411, 923, 91, 603, 347, 859, 219, 731, 475, 987, 59, 571, 315, 827, 187, 699, 443, 955, 123, 635, 379, 891, 251, 763, 507, 1019, 7, 519, 263, 775, 135, 647, 391, 903, 71, 583, 327, 839, 199, 711, 455, 967, 39, 551, 295, 807, 167, 679, 423, 935, 103, 615, 359, 871, 231, 743, 487, 999, 23, 535, 279, 791, 151, 663, 407, 919, 87, 599, 343, 855, 215, 727, 471, 983, 55, 567, 311, 823, 183, 695, 439, 951, 119, 631, 375, 887, 247, 759, 503, 1015, 15, 527, 271, 783, 143, 655, 399, 911, 79, 591, 335, 847, 207, 719, 463, 975, 47, 559, 303, 815, 175, 687, 431, 943, 111, 623, 367, 879, 239, 751, 495, 1007, 31, 543, 287, 799, 159, 671, 415, 927, 95, 607, 351, 863, 223, 735, 479, 991, 63, 575, 319, 831, 191, 703, 447, 959, 127, 639, 383, 895, 255, 767, 511, 1023};


// functions for radix 4 FFT
using namespace std;
typedef  unsigned long    ulong;
//using std::complex;
//typedef  complex<double>  Complex;

// sundiff()
template <typename Type>
static inline void sumdiff(Type &a, Type &b)
// {a, b}  <--| {a+b, a-b}
{ Type t=a-b; a+=b; b=t; }


template <typename Type>
static inline void sumdiff_r(Type &a, Type &b)
// {a, b}  <--| {b+a, b-a}
// Up to scaling (by a factor 2) the inverse of diffsum(a,b).
{ Type t=b-a; a+=b; b=t; }


template <typename Type>
static inline void sumdiff05(Type &a, Type &b)
// {a, b}  <--| {0.5*(a+b),  0.5*(a-b)}
{ Type t=(a-b)*0.5; a+=b; a*=0.5; b=t; }

template <typename Type>
static inline void sumdiff05_r(Type &a, Type &b)
// {a, b}  <--| {0.5*(a+b), 0.5*(b-a)}
{ Type t=(b-a)*0.5; a+=b; a*=0.5; b=t; }

template <typename Type>
static inline void diffsum(Type &a, Type &b)
// {a, b}  <--| {a-b, a+b}
// Up to scaling (by a factor 2) the inverse of sumdiff_r(a,b).
{ Type t=a-b; b+=a; a=t; }

//template <typename Type>
//static inline void diffsum_r(Type &d, Type &a)
//// {a, b}  <--| {b-a, a+b}
//{ Type t=b-a; b+=a; a=t; }


// ------ 4 arg versions:

template <typename Type>
static inline void sumdiff(Type a, Type b, Type &s, Type &d)
// {s, d}  <--| {a+b, a-b}
{ s=a+b; d=a-b; }

template <typename Type>
static inline void sumdiff05(Type a, Type b, Type &s, Type &d)
// {s, d}  <--| {0.5*(a+b), 0.5*(a-b)}
{ s=(a+b)*0.5; d=(a-b)*0.5; }


// ------ 3 arg versions used in split-radix FFTs:

template <typename Type>
static inline void sumdiff3(Type &a, Type b, Type &d)
// {a, b, d} <--| {a+b, b, a-b}  (used in split-radix FFTs)
// NEVER call like func(a,b,a) or func(a,b,b)
{ d=a-b; a+=b; }

template <typename Type>
static inline void sumdiff3_r(Type &a, Type b, Type &d)
// {a,b,d} <--| {a+b, b, b-a}  (used in split-radix FFTs)
// NEVER call like func(a,b,a) or func(a,b,b)
{ d=b-a; a+=b; }

template <typename Type>
static inline void diffsum3(Type a, Type &b, Type &s)
// {a, b, s} <--| {a, a-b, a+b}  (used in split-radix FFTs)
// NEVER call like func(a,b,a) or func(a,b,b)
{ s=a+b; b=a-b; }

template <typename Type>
static inline void diffsum3_r(Type a, Type &b, Type &s)
// {a, b, s} <--| {a, b-a, a+b}  (used in split-radix FFTs)
// NEVER call like func(a,b,a) or func(a,b,b)
{ s=a+b; b-=a; }



// cmult():
static inline void cmult(double c, double s,
                         double x, double y,
                         double &u, double &v)
// {u,v} <--| {x*c-y*s, x*s+y*c}
{
    u = x * c - y * s;
    v = y * c + x * s;
}
// -------------------------

static inline void cmult(double c, double s,
                         complex x, complex y,
                         complex &u, complex &v)
// {u,v} <--| {x*c-y*s, x*s+y*c}
// used in generated complex fhts
{
    u = x * c - y * s;
    v = y * c + x * s;
}
// -------------------------

static inline void cmult(complex c, complex s,
                         complex x, complex y,
                         complex &u, complex &v)
// {u,v} <--| {x*c-y*s, x*s+y*c}
// used in generated complex fhts
{
    u = x * c - y * s;
    v = y * c + x * s;
}
// -------------------------


static inline void csqr(double a, double b,
                        double &u, double &v)
// {u,v} <--| {a*a-b*b, 2*a*b}
{
    u = a * a - b * b;
    v = a * b;
    v += v;
}
// -------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////


void
fft8_dit_core_p1(double *fr, double *fi)
// 8-point decimation in time FFT
// isign = +1
// input data must be in revbin_permuted order
//.
// Cf. Nussbaumer p.148f
{
    // INPUT_RE:
    double t1r = fr[0] + fr[1];
    double t2r = fr[2] + fr[3];
    double t7r = t1r + t2r;
    double t3r = fr[4] - fr[5];
    double t4r = fr[4] + fr[5];
    double t5r = fr[6] + fr[7];
    double t8r = t4r + t5r;
    double t6r = fr[6] - fr[7];
    
    double m0r = t7r + t8r;
    double m1r = t7r - t8r;
    double m2r = t1r - t2r;
    double m3r = fr[0] - fr[1];
    double m4r = M_SQRT1_2 * (t3r - t6r);
    
#define m5i t6r
#define m6i t7r
#define m7i t8r
    m7i = M_SQRT1_2 * (t3r + t6r);
    m5i = t5r - t4r;
    m6i = fr[3] - fr[2];
    
    // INPUT_IM:
    double t1i = fi[0] + fi[1];
    double t2i = fi[2] + fi[3];
    double t7i = t1i + t2i;
    double t3i = fi[4] - fi[5];
    double t4i = fi[4] + fi[5];
    double t5i = fi[6] + fi[7];
    double t8i = t4i + t5i;
    double t6i = fi[6] - fi[7];
    
    double m0i = t7i + t8i;
    double m1i = t7i - t8i;
    double m2i = t1i - t2i;
    double m3i = fi[0] - fi[1];
    double m4i = M_SQRT1_2 * (t3i - t6i);
    
#define m5r t6i
#define m6r t7i
#define m7r t8i
    m7r = M_SQRT1_2 * (t3i + t6i);
    m5r = t4i - t5i;
    m6r = fi[2] - fi[3];
    
#define s1r t1r
#define s2r t2r
#define s3r t3r
#define s4r t4r
    s1r = m3r + m4r;
    s2r = m3r - m4r;
    s3r = m6r + m7r;
    s4r = m6r - m7r;
    
    // OUTPUT_RE:
    fr[0] = m0r;
    fr[7] = s1r + s3r;
    fr[6] = m2r + m5r;
    fr[5] = s2r - s4r;
    fr[4] = m1r;
    fr[3] = s2r + s4r;
    fr[2] = m2r - m5r;
    fr[1] = s1r - s3r;
    
#define s1i t1r
#define s2i t2r
#define s3i t3r
#define s4i t4r
    s1i = m3i + m4i;
    s2i = m3i - m4i;
    s3i = m6i - m7i;
    s4i = m6i + m7i;
    
    // OUTPUT_IM:
    fi[0] = m0i;
    fi[7] = s1i + s3i;
    fi[6] = m2i + m5i;
    fi[5] = s2i - s4i;
    fi[4] = m1i;
    fi[3] = s2i + s4i;
    fi[2] = m2i - m5i;
    fi[1] = s1i - s3i;
}
// -------------------------

#undef s1r
#undef s2r
#undef s3r
#undef s4r

#undef s1i
#undef s2i
#undef s3i
#undef s4i

#undef m5r
#undef m6r
#undef m7r

#undef m5i
#undef m6i
#undef m7i
// -------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////


static const ulong LX = 2;

void
fft_dit4_core_p1(double *fr, double *fi, ulong ldn)
// Auxiliary routine for fft_dit4()
// Decimation in time (DIT) radix-4 FFT
// Input data must be in revbin_permuted order
// ldn := base-2 logarithm of the array length
// Fixed isign = +1
{
    const ulong n = (1UL<<ldn);
    
    if ( n<=2 )
    {
        if ( n==2 )
        {
            sumdiff(fr[0], fr[1]);
            sumdiff(fi[0], fi[1]);
        }
        return;
    }
    
    
    ulong ldm = (ldn&1);
    if ( ldm!=0 )  // n is not a power of 4, need a radix-8 step
    {
        for (ulong i0=0; i0<n; i0+=8)  fft8_dit_core_p1(fr+i0, fi+i0);
    }
    else
    {
        for (ulong i0=0; i0<n; i0+=4)
        {
            double xr, yr, ur, vr, xi, yi, ui, vi;
            ulong i1 = i0 + 1;
            ulong i2 = i1 + 1;
            ulong i3 = i2 + 1;
            
            sumdiff(fr[i0], fr[i1], xr, ur);
            sumdiff(fr[i2], fr[i3], yr, vi);
            sumdiff(fi[i0], fi[i1], xi, ui);
            sumdiff(fi[i3], fi[i2], yi, vr);
            
            sumdiff(ui, vi, fi[i1], fi[i3]);
            sumdiff(xi, yi, fi[i0], fi[i2]);
            sumdiff(ur, vr, fr[i1], fr[i3]);
            sumdiff(xr, yr, fr[i0], fr[i2]);
        }
    }
    ldm += 2*LX;
    
    
    for ( ; ldm<=ldn; ldm+=LX)
    {
        ulong m = (1UL<<ldm);
        ulong m4 = (m>>LX);
        const double ph0 = 2.0 * M_PI / (double)m;
        double ph = 0.0;
        for (ulong j=0; j<m4; j++)
        {
            double c, s, c2, s2, c3, s3;
            //            SinCos(j*ph0, &s, &c);
            s = sin(ph); c = cos(ph); //SinCos(ph, &s, &c);
            ph += ph0;
            csqr(c, s, c2, s2);
            cmult(c, s, c2, s2, c3, s3);
            
            for (ulong r=0; r<n; r+=m)
            {
                ulong i0 = j + r;
                ulong i1 = i0 + m4;
                ulong i2 = i1 + m4;
                ulong i3 = i2 + m4;
                
                double xr, yr, ur, vr, xi, yi, ui, vi;
                cmult(c2, s2, fr[i1], fi[i1], xr, xi);
                
                sumdiff3_r(xr, fr[i0], ur);
                sumdiff3_r(xi, fi[i0], ui);
                
                cmult(c,  s,  fr[i2], fi[i2], yr, vr);
                cmult(c3, s3, fr[i3], fi[i3], vi, yi);
                
                sumdiff(yr, vi);
                sumdiff(yi, vr);
                
                sumdiff(ur, vr, fr[i1], fr[i3]);
                sumdiff(ui, vi, fi[i1], fi[i3]);
                sumdiff(xr, yr, fr[i0], fr[i2]);
                sumdiff(xi, yi, fi[i0], fi[i2]);
            }
        }
    }
}
// -------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////

// end


complex* CFFT::convolutionF(const complex *input,const complex *filter, long nSIG, long NFIL, long &NFFT)
{
    //Check for invalid inputs.
    if(input == NULL || filter == NULL)
    {
        cout << "Could not perform convolution on empty arrays!" << endl;
        return NULL;
    }
    
    bool NFFTChanged = false;
    //If NFFT not a power of 2, or it is smaller than signal or filter, prompt for new.
    while (log(NFFT) / log(2) != (int)(log(NFFT) / log(2)) || NFFT < nSIG || NFFT < NFIL)
    {
        cout << "Please input a valid NFFT, which is >= nSIG(" << nSIG << ") and >= NFIL(" << NFIL <<") : ";
        cin >> NFFT;
        NFFTChanged = true;
    }
    
    //Perform zero padding.
    complex *fInput, *fFilter;
    
    fInput = new complex[NFFT];
    for(int i = 0; i < nSIG; i++)
        fInput[i] = input[i];
    
    fFilter = new complex[NFFT];
    for(int i = 0; i < NFIL; i++)
        fFilter[i] = filter[i];
    
    //Store the output data.
    complex *output = new complex[NFFT];
    
    //Perform FFT on both input and filter.
    CFFT::Forward(fInput, (unsigned int)NFFT);
    CFFT::Forward(fFilter, (unsigned int)NFFT);
    
    for(int i = 0; i < NFFT; i++)
        output[i] = fInput[i] * fFilter[i];
    
    return output;
}

//Convolution function which returns the time representation of the result.
//NFFT is the FFT size, nSIG is the size for the input, NFIL is the size of the filter.
complex* CFFT::convolutionT(const complex *input,const complex *filter, long nSIG, long NFIL, long &NFFT)
{
    //Store the output data.
    complex *output = convolutionF(input, filter, nSIG, NFIL, NFFT);
    
    //Perform IFFT on the ouput.
    CFFT::Inverse(output, (unsigned int)NFFT);
    
    return output;
}

complex* CFFT::stereoConvMonoInputF(const complex *input,const complex *filterLeft,const complex *filterRight, long nSIG, long NFILL, long NFILR, long &NFFT)
{
    complex *result = stereoConvMonoInputT(input, filterLeft, filterRight, nSIG, NFILL, NFILR, NFFT);
    CFFT::Forward(result, (unsigned int)NFFT);
    return result;
}

//The size of the ouput will be 2 times of the size of FFT for the input signal and the NFFT value
//will be doubled after running the function
complex* CFFT::stereoConvMonoInputT(const complex *input,const complex *filterLeft,const complex *filterRight, long nSIG, long NFILL, long NFilR, long &NFFT)
{
    complex *tempLeft = new complex[NFFT];
    complex *tempRight = new complex[NFFT];
    complex *result = new complex[2*NFFT];
    
    tempLeft = CFFT::convolutionT(input, filterLeft, nSIG, NFILL, NFFT);
    tempRight = CFFT::convolutionT(input, filterRight, nSIG, NFILL, NFFT);
    NFFT = NFFT * 2;
    for (int i = 0; i < NFFT / 2; i++)
    {
        result[2 * i] = tempLeft[i];
        result[2 * i + 1] = tempRight[i];
    }
    
    delete tempLeft;
    delete tempRight;
    
    return result;
}

/*
 complex* CFFT::stereoConvStereoInputT(const complex *input, const complex *filterLeft, const complex *filterRight, long nSIG, long NFILL, long NFILR, long &NFFT)
 {
	NFFT = NFFT / 2;
	complex *leftTemp = new complex[NFFT];
	complex *rightTemp = new complex[NFFT];
 
	
 
	
 
	return 0;
 
 }
 
 */


//storing the an array into a text file
//filename is the file name you want to store the data into
//datatype represents the data you wanna store: real/real+imag/amplitude
void CFFT::storingData(complex *data, int NFFT,string temp ,char datatype)
{
    //string temp= filename;
    ofstream outputFile(temp.c_str());
    if (outputFile.is_open())
    {
        switch (datatype)
        {
            case 'r':
                //ofstream outputFile("real.txt");
                for (int i = 0; i < NFFT; i++)
                    outputFile << data[i].re() << endl;
                break;
                
            case 'c':
                for (int i = 0; i < NFFT; i++)
                    outputFile << data[i].re() << "    " << data[i].im() << endl;
                break;
            case 'a':
                for (int i = 0; i < NFFT; i++)
                    outputFile << sqrt(data[i].norm()) << endl;
                break;
        }
        outputFile.close();
    }
}



//   FORWARD FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
bool CFFT::Forward(const complex *const Input, complex *const Output, const unsigned int N)
{
    //   Check input parameters
    if (!Input || !Output || N < 1 || N & (N - 1))
        return false;
    //   Initialize data
    Rearrange(Input, Output, N);
    //   Call FFT implementation
    Perform(Output, N);
    //   Succeeded
    return true;
}

//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
//     Data - both input data and output
//     N    - length of input data
bool CFFT::Forward(complex *const Data, const unsigned int N)
{
    //   Check input parameters
    if (!Data || N < 1 || N & (N - 1))
        return false;
    //   Rearrange
    Rearrange(Data, N);
    //   Call FFT implementation
    Perform(Data, N);
    //   Succeeded
    return true;
}

//   INVERSE FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
//     Scale  - if to scale result
bool CFFT::Inverse(const complex *const Input, complex *const Output, const unsigned int N, const bool Scale /* = true */)
{
    //   Check input parameters
    if (!Input || !Output || N < 1 || N & (N - 1))
        return false;
    //   Initialize data
    Rearrange(Input, Output, N);
    //   Call FFT implementation
    Perform(Output, N, true);
    //   Scale if necessary
    if (Scale)
        CFFT::Scale(Output, N);
    //   Succeeded
    return true;
}

//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
//     Data  - both input data and output
//     N     - length of both input data and result
//     Scale - if to scale result
bool CFFT::Inverse(complex *const Data, const unsigned int N, const bool Scale /* = true */)
{
    //   Check input parameters
    if (!Data || N < 1 || N & (N - 1))
        return false;
    //   Rearrange
    Rearrange(Data, N);
    //   Call FFT implementation
    Perform(Data, N, true);
    //   Scale if necessary
    if (Scale)
        CFFT::Scale(Data, N);
    //   Succeeded
    return true;
}

//   Rearrange function
void CFFT::Rearrange(const complex *const Input, complex *const Output, const unsigned int N)
{
    //	//   Data entry position
    //	unsigned int Target = 0;
    //	//   Process all positions of input signal
    for (unsigned int Position = 0; Position < N; ++Position)
    {
        Output[Position] = Input[Position];
        //		//  Set data entry
        //		Output[Target] = Input[Position];
        //		//   Bit mask
        //		unsigned int Mask = N;
        //		//   While bit is set
        //		while (Target & (Mask >>= 1))
        //			//   Drop bit
        //			Target &= ~Mask;
        //		//   The current bit is 0 - set it
        //		Target |= Mask;
    }
}

//   Inplace version of rearrange function

void CFFT::Rearrange(complex *const Data, const unsigned int N)
{
    //	//   Swap position
    //	unsigned int Target = 0;
    //	//   Process all positions of input signal
    //	for (unsigned int Position = 0; Position < N; ++Position)
    //	{
    //		//   Only for not yet swapped entries
    //		if (Target > Position)
    //		{
    //			//   Swap entries
    //			const complex Temp(Data[Target]);
    //			Data[Target] = Data[Position];
    //			Data[Position] = Temp;
    //		}
    //		//   Bit mask
    //		unsigned int Mask = N;
    //		//   While bit is set
    //		while (Target & (Mask >>= 1))
    //			//   Drop bit
    //			Target &= ~Mask;
    //		//   The current bit is 0 - set it
    //		Target |= Mask;
    //	}
    
}

//   FFT implementation
void CFFT::Perform(complex *const Data, const unsigned int N, const bool Inverse /* = false */)
{
    
    //	const double pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;
    //	//   Iteration through dyads, quadruples, octads and so on...
    //	for (unsigned int Step = 1; Step < N; Step <<= 1)
    //	{
    //		//   Jump to the next entry of the same transform factor
    //		const unsigned int Jump = Step << 1;
    //		//   Angle increment
    //		const double delta = pi / double(Step);
    //		//   Auxiliary sin(delta / 2)
    //		const double Sine = sin(delta * .5);
    //		//   Multiplier for trigonometric recurrence
    //		const complex Multiplier(-2. * Sine * Sine, sin(delta));
    //		//   Start value for transform factor, fi = 0
    //		complex Factor(1.);
    //		//   Iteration through groups of different transform factor
    //		for (unsigned int Group = 0; Group < Step; ++Group)
    //		{
    //			//   Iteration within group
    //			for (unsigned int Pair = Group; Pair < N; Pair += Jump)
    //			{
    //				//   Match position
    //				const unsigned int Match = Pair + Step;
    //				//   Second term of two-point transform
    //				const complex Product(Factor * Data[Match]);
    //				//   Transform for fi + pi
    //				Data[Match] = Data[Pair] - Product;
    //				//   Transform for fi
    //				Data[Pair] += Product;
    //			}
    //			//   Successive transform factor via trigonometric recurrence
    //			Factor = Multiplier * Factor + Factor;
    //		}
    //	}
    
    double dataReal [N];
    double dataImag [N];
    for (int i = 0; i < N; i++) {
        dataReal[i] = Data[PermutedIndex[i]].re();
        dataImag[i] = Data[PermutedIndex[i]].im();
    }
    double *fr = dataReal;
    double *fi = dataImag;
    int ldn = log2(N);
    
    fft_dit4_core_p1(fr, fi, ldn);
    
    //    if(!Inverse){fft_dit4_core_p1(fr, fi, ldn);}
    //    else{ fft_dit4_core_p1(fi, fr, ldn);}
    
    
    for (int i = 0; i < N; i++){
        complex tmp(fr[i], -fi[i]);
        Data[i] = tmp;
    }
    
    
}

//   Scaling of inverse FFT result
void CFFT::Scale(complex *const Data, const unsigned int N)
{
    const double Factor = 1. / double(N);
    //   Scale all data entries
    for (unsigned int Position = 0; Position < N; ++Position)
        Data[Position] *= Factor;
}
