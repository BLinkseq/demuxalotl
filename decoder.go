package main

import (
	_ "embed"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

type decoder struct {
	start, end int               // segment location in the padded R1
	codes      []uint8           // barcodes flattened, each end-start long
	ids        map[string]string // barcode sequence -> codec id without '@'
	seqs       []string
	conc       float64 // prior, (1-noise)/n
	noiseP     float64 // noise * random barcode probability
	randomP    float64 // 4^-length
	confidence float64
}

// decode returns the best barcode index and whether the call passed the PAMLD
// noise and confidence filters (false == QC fail).
//
// Every float operation matches the naive per-barcode loop (same operands, same
// order), so results are bit-identical; the speedup is only from hoisting the
// score lookups out of the barcode loop and scoring four barcodes at a time so
// their independent Kahan chains overlap in the CPU pipeline.
func (d *decoder) decode(obs, qual []byte) (best int, pass bool) {
	n := len(obs)
	var sc [64][4]float64 // sc[i][e]: score of expecting base e at position i
	for i := range n {
		o, q := baseCode[obs[i]], min(max(qual[i], 33)-33, 127)
		for e := range 4 {
			sc[i][e] = score[q][e][o]
		}
	}
	nb := len(d.seqs)
	var stack [128]float64
	cond := stack[:min(nb, len(stack))]
	if nb > len(stack) {
		cond = make([]float64, nb)
	}
	b := 0
	for ; b+4 <= nb; b += 4 {
		c0, c1, c2, c3 := d.codes[b*n:(b+1)*n], d.codes[(b+1)*n:(b+2)*n], d.codes[(b+2)*n:(b+3)*n], d.codes[(b+3)*n:(b+4)*n]
		var k0, k1, k2, k3 kahan
		for i := range n {
			s := &sc[i]
			k0.add(s[c0[i]])
			k1.add(s[c1[i]])
			k2.add(s[c2[i]])
			k3.add(s[c3[i]])
		}
		cond[b] = pow10(-0.1 * k0.sum)
		cond[b+1] = pow10(-0.1 * k1.sum)
		cond[b+2] = pow10(-0.1 * k2.sum)
		cond[b+3] = pow10(-0.1 * k3.sum)
	}
	for ; b < nb; b++ {
		var k kahan
		code := d.codes[b*n : (b+1)*n]
		for i := range n {
			k.add(sc[i][code[i]])
		}
		cond[b] = pow10(-0.1 * k.sum)
	}
	best = -1
	var bestP, bestCond float64
	var sigma kahan
	for b, c := range cond {
		p := float64(c * d.conc) // explicit conversion: no fused multiply-add on arm64, same bits everywhere
		sigma.add(p)
		if p > bestP {
			best, bestP, bestCond = b, p, c
		}
	}
	if best < 0 || bestCond <= d.randomP {
		return -1, false
	}
	sigma.add(d.noiseP)
	return best, bestP/sigma.sum > d.confidence
}

// ── config ────────────────────────────────────────────────────────────────────

type codec struct {
	Barcode []string `json:"barcode"`
}

type config struct {
	Template struct {
		Transform struct{ Token []string } `json:"transform"`
	} `json:"template"`
	Molecular []struct {
		Base       string                   `json:"base"`
		Noise      float64                  `json:"noise"`
		Confidence float64                  `json:"confidence_threshold"`
		Transform  struct{ Token []string } `json:"transform"`
	} `json:"molecular"`
	Decoder map[string]struct {
		Codec map[string]codec `json:"codec"`
	} `json:"decoder"`
}

// parseToken reads a pheniqs "read:start:end" token ("" end means to the end of the read).
func parseToken(t string) (start, end int, err error) {
	f := strings.Split(t, ":")
	if len(f) != 3 || f[0] != "0" {
		return 0, 0, fmt.Errorf("unsupported token %q, expected 0:start:end", t)
	}
	if start, err = strconv.Atoi(f[1]); err != nil {
		return
	}
	if f[2] != "" {
		end, err = strconv.Atoi(f[2])
	}
	return
}

// defaultSchema is used unless --schema names a file. It is a copy of
// harpy/templates/pheniqs.config.json (a test keeps the two in sync).
//
//go:embed schema.json
var defaultSchema []byte

// loadConfig builds the decoders from the schema file at path, or from the embedded
// default schema when path is empty.
func loadConfig(path string) (decs []*decoder, r1Start int, err error) {
	data := defaultSchema
	if path != "" {
		if data, err = os.ReadFile(path); err != nil {
			return nil, 0, err
		}
	}
	var cfg config
	if err = json.Unmarshal(data, &cfg); err != nil {
		if path != "" {
			err = fmt.Errorf("%s: %w", path, err)
		}
		return
	}
	if len(cfg.Template.Transform.Token) == 0 {
		return nil, 0, errors.New("config has no template transform token")
	}
	if r1Start, _, err = parseToken(cfg.Template.Transform.Token[0]); err != nil {
		return
	}
	if len(cfg.Molecular) != 4 {
		return nil, 0, fmt.Errorf("expected 4 molecular decoders (D,A,B,C), got %d", len(cfg.Molecular))
	}
	for _, m := range cfg.Molecular {
		d := &decoder{ids: map[string]string{}}
		if len(m.Transform.Token) != 1 {
			return nil, 0, errors.New("each molecular decoder needs exactly one token")
		}
		if d.start, d.end, err = parseToken(m.Transform.Token[0]); err != nil {
			return
		}
		if d.end-d.start > 64 {
			return nil, 0, errors.New("barcode segment longer than 64 bases")
		}
		cdc := cfg.Decoder[m.Base].Codec
		keys := make([]string, 0, len(cdc))
		for k := range cdc {
			keys = append(keys, k)
		}
		// ponytail: pheniqs keeps file order for tie-breaking, sorted ids match it for the shipped config
		sort.Strings(keys)
		for _, k := range keys {
			if len(cdc[k].Barcode) == 0 {
				continue
			}
			bc := strings.ToUpper(cdc[k].Barcode[0])
			if len(bc) != d.end-d.start {
				return nil, 0, fmt.Errorf("barcode %s (%s) is not %d bases", k, bc, d.end-d.start)
			}
			for i := range bc {
				d.codes = append(d.codes, baseCode[bc[i]])
			}
			d.seqs = append(d.seqs, bc)
			d.ids[bc] = strings.TrimPrefix(k, "@")
		}
		if len(d.seqs) == 0 {
			return nil, 0, fmt.Errorf("no barcodes for decoder base %q", m.Base)
		}
		d.conc = (1 - m.Noise) / float64(len(d.seqs))
		d.randomP = math.Pow(4, -float64(d.end-d.start))
		d.noiseP = m.Noise * d.randomP
		d.confidence = m.Confidence
		decs = append(decs, d)
	}
	return decs, r1Start, nil
}
