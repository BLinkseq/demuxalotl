package main

import "math"

// ── PAMLD decoder (port of pheniqs pamld.cpp) ─────────────────────────────────

// score[q][e][o]: -10log10 P(observe o | expect e, Phred q); o==4 is N.
// Mirrors PhredScale::assemble_substitution_lookup, including q==0 scoring 0.
var score [128][4][5]float64

var baseCode [256]uint8 // A,C,G,T -> 0..3, everything else -> 4 (N)

func init() {
	for i := range baseCode {
		baseCode[i] = 4
	}
	for i, c := range "ACGT" {
		baseCode[c] = uint8(i)
		baseCode[c+32] = uint8(i)
	}
	for q := 1; q < 128; q++ {
		perr := pow10(-0.1 * float64(q))
		for e := range 4 {
			for o := range 5 {
				switch o {
				case 4:
					score[q][e][o] = 10 * math.Log10(4)
				case e:
					score[q][e][o] = -10 * math.Log10(1-perr)
				default:
					score[q][e][o] = float64(q)
				}
			}
		}
	}
}

var ln10 = math.Log(10)

// pow10 returns math.Pow(10, y) bit for bit for finite y <= 0, minus the
// per-call Log(10) and the special-case ladder; same Modf/Exp/Frexp/Ldexp steps
// as math.pow. demux_test.go checks it against math.Pow.
func pow10(y float64) float64 {
	switch y {
	case 0:
		return 1
	case -0.5:
		return 1 / math.Sqrt(10)
	}
	yi, yf := math.Modf(-y)
	if yf > 0.5 {
		yf--
		yi++
	}
	a1 := 1.0
	if yf != 0 {
		a1 = math.Exp(yf * ln10)
	}
	x1, xe := 0.625, 4 // Frexp(10)
	ae := 0
	for i := int64(yi); i != 0; i >>= 1 {
		if i&1 == 1 {
			a1 *= x1
			ae += xe
		}
		x1 *= x1
		xe <<= 1
		if x1 < .5 {
			x1 += x1
			xe--
		}
	}
	return math.Ldexp(1/a1, -ae)
}

// kahan is a compensated running sum, as in pheniqs: c holds the low-order
// bits lost by each add and feeds them back into the next one.
type kahan struct{ sum, c float64 }

func (k *kahan) add(x float64) {
	y := x - k.c
	t := k.sum + y
	k.c = (t - k.sum) - y
	k.sum = t
}
