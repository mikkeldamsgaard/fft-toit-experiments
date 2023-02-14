import math

class Complex:
  real/num? := ?
  imaginary/num? := ?

  constructor .real .imaginary:

  constructor.nth_root k/int N/int:
    real = math.cos 2*math.PI*k/N
    imaginary = math.sin 2*math.PI*k/N

  operator * other/Complex -> Complex:
    return Complex
        (real*other.real - imaginary*other.imaginary)
        (real*other.imaginary + imaginary*other.real)

  operator + other/Complex -> Complex:
    return Complex
      real + other.real
      imaginary + other.imaginary

  operator - other/Complex -> Complex:
    return Complex
      real - other.real
      imaginary - other.imaginary

  operator / scalar/num -> Complex:
    return Complex real/scalar imaginary/scalar

  static ONE/Complex ::= Complex 1 0

  stringify -> string:
    return "$(%.2f real) $(%.2f imaginary)i"

  set other/Complex:
    real = other.real
    imaginary = other.imaginary

  mult_with other/Complex:
    r := (real*other.real - imaginary*other.imaginary)
    i := (real*other.imaginary + imaginary*other.real)
    real = r
    imaginary = i

  add_with other/Complex:
    real += other.real
    imaginary += other.imaginary

  minus_with other/Complex:
    real -= other.real
    imaginary -= other.imaginary

fft x/List -> List:
 return fft_iterative x

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


fft_iterative x/List -> List:
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

ROOTS ::= (::
   result := List 10
   for n := 1; n < 10; n++: result[n] = Complex.nth_root -1 1 << n
   result
).call

fft_recursive x/List start_index/int=0 stride/int=1 logN/int=(math.log x.size 2).to_int -> List:
  if logN == 0: return x[start_index..start_index+1]

  e := fft_recursive x start_index 2 * stride logN - 1
  o := fft_recursive x start_index + stride 2 * stride logN - 1

  N ::= 1 << logN
  N_half ::= N / 2
  fft_result := List N
  root ::= ROOTS[logN]
  twiddle := Complex.ONE
  for k := 0; k < N_half; k++:
    p := e[k]
    q := twiddle * o[k]
    twiddle = twiddle * root
    fft_result[k] = p + q
    fft_result[k+N_half] = p - q

  return fft_result

ifft x/List -> List:
  N ::= x.size
  pre := List N
  pre[0] = x[0]
  for i := 1; i < N; i++: pre[N - i] = x[i]
  post := fft_iterative pre
  for i := 0; i < N; i++: post[i] /= N
  return post

SAMPLE_SIZE ::= 480

main:
  block := List SAMPLE_SIZE
  SAMPLE_SIZE.repeat: |n|
    block[n] = math.cos 120*2*math.PI/(SAMPLE_SIZE)*n

  fft_block := List 512 (Complex 0 0)
  SAMPLE_SIZE.repeat: fft_block[it] = Complex block[it] 0

  ffted_block := fft fft_block
  start_forward := Time.monotonic_us
  10.repeat:
    ffted_block = fft fft_block
  end_forward := Time.monotonic_us
  iffted_block := ifft ffted_block

  distance_squared := 0
  SAMPLE_SIZE.repeat: | k |
    diff := iffted_block[k].real - block[k]
    distance_squared += diff*diff

  end := Time.monotonic_us
  if (math.sqrt distance_squared) > 1.0e-11:
    throw "Invalid algorithm $(math.sqrt distance_squared)"

  print "timing: forward: $(%.3f (end_forward*1.0 - start_forward)/10000)"
  y := List SAMPLE_SIZE*2
