# Experiments with modified FFT and discrete cosine transform in Toit

The blog explores an MDCT implementation in toit and how it 
progressively is made to run faster. All timings is calculated on a M1MAX Mac Pro

When used in audio compression at sample rate 48Khz then the mdct is called for
each 10ms of audio, so that gives the timings some relational meaning.

### Naive approach
The first very naive implementation which is 
horribly inefficient implementation
```toit
import math
N ::= 480
mdct x/List -> List:
  X := List N
  N.repeat: | k |
    X[k] = 0
    (2*N).repeat: | n |
      X[k] += x[n]*(math.cos (math.PI/N*(n+0.5+N/2)*(k+0.5)));
  return X
```
Timing: 126 ms

### Modified naive approach
First, we want to eliminate all constant calculations
```toit
import math
N ::= 480
mdct x/List -> List:
  X := List N
  k_half := 0.5
  N.repeat: | k |
    X[k] = 0
    n_half := N/2+0.5
    (2*N).repeat: | n |
      X[k] += x[n]*(math.cos (math.PI/N*n_half*k_half));
      n_half += 1
    k_half += 1
  return X
```
Timing: 86 ms

So, this is much better with these simple rewrites. This is still an $O(n^2)$ implementation.
### Replacing repeat wit explicit for loops

```toit
import math
N ::= 480
mdct x/List -> List:
  X := List N
  k_half := 0.5
  for k := 0; k < N; k++:
    X[k] = 0
    n_half := N/2+0.5
    for n := 0; n < 2*N; n++:
      X[k] += x[n]*(math.cos (math.PI/N*n_half*k_half));
      n_half += 1
    k_half += 1
  return X
```
Timing: 67.6 ms

So, for now the conclusion is that `repeat` is somewhat slower than doing
explicit for loops and that keeping intermediate results in local variable helps.

Note, one would think that the `math.PI/N` and `2*N` could also benefit from
being stored in local variables, but these a globally constant and it is slower to
put them in local variables.

### Divide an conquer.
The MDCT algorithm is able to be implemented using FFT. So lets start by looking at implementing an FFT.

The FFT is simplest when the input is a power of 2. So let us assume that it is, just zero extend the 
original input up to 512 points and just use the first 480 results.

First we create a complex class to help us out
```toit
import math
class Complex:
  real/num
  imganiry/num

  constructor .real .imganiry:

  constructor.nth_root k/int N/int:
    real = math.cos 2*math.PI*k/N
    imganiry = math.sin 2*math.PI*k/N

  operator * other/Complex -> Complex:
    return Complex
        (real*other.real - imganiry*other.imganiry)
        (real*other.imganiry + imganiry*other.real)

  operator + other/Complex -> Complex:
    return Complex
      real + other.real
      imganiry + other.imganiry

  operator - other/Complex -> Complex:
    return Complex
      real - other.real
      imganiry - other.imganiry

  operator / scalar/num -> Complex:
    return Complex real/scalar imganiry/scalar
    
  static ONE/Complex ::= Complex 1 0  
```

Then we can implement a naive FFT based on Cooley-Tukey
```toit
fft x/List -> List:
  N ::= x.size
  if N == 1: return x

  odd := List N / 2
  even := List N / 2
  N.repeat:
    if it % 2 == 0: even[it / 2] = x[it]
    else: odd[(it - 1) / 2] = x[it]

  e := fft even
  o := fft odd

  fft_result := List N
  (N / 2).repeat: | k |
    p := e[k]
    q := (Complex.nth_root -k N) * o[k]
    fft_result[k] = p + q
    fft_result[k+N/2] = p - q

  return fft_result
```
Timing: 2.189 ms

The first thing to notice is that this has a lot of array copying. 
That can be optimized with sending more context to the fft

```toit
fft x/List start_index/int=0 stride/int=1 N/int=x.size -> List:
  if N == 1: return x[start_index..start_index+1]

  e := fft x start_index 2 * stride N / 2
  o := fft x start_index + stride 2 * stride N / 2

  fft_result := List N
  (N/2).repeat: | k |
    p := e[k]
    q := (Complex.nth_root -k N)*o[k]
    fft_result[k] = p + q
    fft_result[k+N/2] = p - q

  return fft_result
```
Timing: 1.872 ms

Instead of calculating the nth_root inside the loop with cos/sin, the primitive root
can be precalculated and the loop can use complex multiplication. The new implementation shifts to 
using logN as a context parameter and the exp(-2*PI/N) is precalculated for
all the N's we encounter (2,4,8,16,..,512)
```toit
ROOTS ::= (::
   result := List 10
   for n := 1; n < 10; n++: result[n] = Complex.nth_root -1 1 << n
   result
).call

fft x/List start_index/int=0 stride/int=1 logN/int=(math.log x.size 2).to_int -> List:
  if logN == 0: return x[start_index..start_index+1]

  e := fft x start_index 2 * stride logN - 1
  o := fft x start_index + stride 2 * stride logN - 1

  N ::= 1 << logN
  N_half ::= N / 2
  fft_result := List N
  root ::= ROOTS[logN]
  twiddle := Complex 1 0
  for k := 0; k < N_half; k++:
    p := e[k]
    q := twiddle * o[k]
    twiddle = twiddle * root
    fft_result[k] = p + q
    fft_result[k+N_half] = p - q

  return fft_result
```
Timing: 1.584 (1.516 with -O2)

The next experiment is to remove recursion. This requires some tricks reordering the input array
with BIT reverse indexing.
The algorithm now looks like
```toit
BIT_REVERSE ::= (::
  N ::= 512
  result := List N
  for i := 0; i < N; i++:
    bit_reversed := 0
    for bit := 0; bit < 9; bit++:
      if ( i & (1 << bit) ) != 0: bit_reversed += 1 << (8 - bit)
    result[i] = bit_reversed
  result
).call

fft x/List -> List:
  N ::= x.size
  A ::= List N
  logN ::= 9

  // bit reverse x
  for i := 0; i < N; i++:
    A[BIT_REVERSE[i]] = x[i]

  for stride := 1; stride <= logN; stride++:
    m ::= 1 << stride
    m_half ::= m >> 1
    root ::= ROOTS[stride]
    for k := 0; k < N; k += m:
      twiddle := Complex 1 0
      for j := 0; j < m_half; j++:
        K_PLUS_J ::= k+j
        p := A[K_PLUS_J]
        q := twiddle * A[K_PLUS_J + m_half]
        A[K_PLUS_J] = p + q
        A[K_PLUS_J + m_half] = p - q
        twiddle *= root
        
  return A
```
Timing: 1.725 ms (1.679 with -O2)

So, surprisingly it turns out that the iterative approach is slower than the 
recursive approach, even though there are much fewer array allocations.
