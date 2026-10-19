// demuxolotl - FASTQ in, FASTQ out replacement for the gih-stagger | pheniqs mux | gih-convert chain from harpy.
//
// Per read pair it:
// 1. finds the ME in R1 and applies the stagger padding,
// 2. decodes the D/A/B/C barcode segments with a port of the Pheniqs PAMLD decoder, driven by the same pheniqs.config.json
// 3. writes R1/R2 FASTQ with BX/VX tags
//
// R1 and R2 read names must match pair for pair (a /1 and /2 suffix is ignored), otherwise it exits with an error.
//
// Usage:
//
//	demuxolotl [options] <prefix> <pheniqs.json> <R1.fq[.gz]> <R2.fq[.gz]>
//
// Build:
//
//	go build -ldflags="-s -w" -o demuxolotl demux.go
package demux

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"sync"

	"github.com/shenwei356/bio/seq"
)

const (
	winStart  = 49
	batchSize = 2000
	empty     = "00"
	asciiQ40  = 'I' // Phred 40 + 33, quality of stagger pad bases
)

// ── stagger pads ──────────────────────────────────────────────────────────────

var padSeq = [8]string{"TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", ""}

// ── ME search (as gih-stagger) ────────────────────────────────────────────────

func findME(seq, me []byte, maxMismatch float64) int {
	meLen := len(me)
	winEnd := min(65, len(seq)-meLen)
	if len(seq) < meLen || winStart > winEnd {
		return -1
	}
seqscan:
	for pos := winStart; pos <= winEnd; pos++ {
		mm := 0.0
		for i := range meLen {
			if r := seq[pos+i]; r == 'N' {
				mm += 0.3
			} else if r != me[i] {
				mm++
			}
			if mm > maxMismatch {
				continue seqscan
			}
		}
		return pos
	}
	return -1
}

// ── per-pair processing ───────────────────────────────────────────────────────

type stats struct {
	total, discarded, qcfail, valid, corrected int64
}

func (s *stats) add(o stats) {
	s.total += o.total
	s.discarded += o.discarded
	s.qcfail += o.qcfail
	s.valid += o.valid
	s.corrected += o.corrected
}

// worker owns the scratch space for one goroutine.
type worker struct {
	decs        []*decoder
	r1Start     int
	me          []byte
	maxMismatch float64
	seq, qual   []byte
	bxSet       map[string]struct{}
	st          stats
	qxrx        bool
}

// process handles one pair. keep is false when the ME was not found.
// It appends the two FASTQ records to o1/o2.
func (w *worker) process(p *readPair, o1, o2 *bytes.Buffer) {
	w.st.total++
	fq1, fq2 := &p.fq1, &p.fq2
	mePos := findME(fq1.seq, w.me, w.maxMismatch)
	if mePos == -1 {
		w.st.discarded++
		return
	}
	plen := 7
	if mePos >= 51 && mePos <= 58 {
		plen = mePos - 51
	}
	pad := padSeq[plen]
	bcFrom := 0
	if plen == 6 {
		bcFrom = 1
	}
	// padded R1 = pad + barcode region + sequence after the ME
	w.seq, w.qual = w.seq[:0], w.qual[:0]
	for range pad {
		w.qual = append(w.qual, asciiQ40)
	}
	w.seq = append(w.seq, pad...)
	w.seq = append(w.seq, fq1.seq[bcFrom:mePos]...)
	w.seq = append(w.seq, fq1.seq[mePos+len(w.me):]...)
	w.qual = append(w.qual, fq1.qual[bcFrom:mePos]...)
	w.qual = append(w.qual, fq1.qual[mePos+len(w.me):]...)

	// decode D, A, B, C
	var orig [4][]byte
	var corr [4]string
	var bcqual [4][]byte
	qcFail, changed := false, false
	for i, d := range w.decs {
		orig[i] = w.seq[d.start:d.end]
		bcqual[i] = w.qual[d.start:d.end]
		idx, pass := d.decode(orig[i], w.qual[d.start:d.end])
		if !pass {
			qcFail = true
		}
		if idx >= 0 {
			corr[i] = d.seqs[idx]
			changed = changed || corr[i] != string(orig[i])
		}
	}
	// failed reads keep the observed barcode, passing reads the corrected one
	var id [4]string
	valid := 1
	for i, d := range w.decs {
		s := corr[i]
		if qcFail {
			s = string(orig[i])
		}
		var ok bool
		if id[i], ok = d.ids[s]; !ok {
			id[i], valid = empty, 0
		}
	}
	bx := "A" + id[1] + "C" + id[2] + "B" + id[3] + "D" + id[0]
	w.bxSet[bx] = struct{}{}
	// counts are per read, and both mates carry the same tags
	w.st.valid += int64(2 * valid)
	if qcFail {
		w.st.qcfail += 2
	} else if changed {
		w.st.corrected += 2
	}

	if !w.qxrx {
		writeFastq(o1, fq1.name, "/1", bx, valid, w.seq[w.r1Start:], w.qual[w.r1Start:])
		writeFastq(o2, fq2.name, "/2", bx, valid, fq2.seq, fq2.qual)
	} else {
		QX := make([]byte, 28)
		RX := make([]byte, 28)
		for i := range orig {
			QX = append(QX, orig[i]...)
			RX = append(RX, bcqual[i]...)
		}
		writeFastqQxRx(o1, fq1.name, "/1", bx, valid, QX, RX, w.seq[w.r1Start:], w.qual[w.r1Start:])
		writeFastqQxRx(o2, fq2.name, "/2", bx, valid, QX, RX, fq2.seq, fq2.qual)
	}
}

// ── main ──────────────────────────────────────────────────────────────────────

type result struct{ r1, r2 *bytes.Buffer }

func Demux(r1, r2, prefix, meSeq, schema string, maxMM int, qxrx bool, threads int) error {
	decoders, r1Start, err := loadConfig(schema)
	if err != nil {
		return fmt.Errorf("Error reading config: %v", err)
	}
	// pgzip.NewReaderN with a single block intermittently fails with "gzip: invalid checksum", so never go below 2
	gzBlocks := max(threads/2, 2)
	seq.ValidateSeq = false // reads may hold N and other IUPAC codes, and nothing downstream needs them checked

	//---- R1 fastq reader -----------------
	rdr1, err := newFastqReader(r1, gzBlocks)
	if err != nil {
		return fmt.Errorf("Error opening R1: %v", err)
	}
	defer rdr1.Close()

	//---- R2 fastq reader -----------------
	rdr2, err := newFastqReader(r2, gzBlocks)
	if err != nil {
		return fmt.Errorf("Error opening R2: %v", err)
	}
	defer rdr2.Close()

	//---- R1 fastq writer -----------------
	out1, err := newGzOut(prefix+".R1.fq.gz", max(threads/2, 1))
	if err != nil {
		return fmt.Errorf("Error opening R1 output: %v", err)
	}

	//---- R2 fastq writer -----------------
	out2, err := newGzOut(prefix+".R2.fq.gz", max(threads/2, 1))
	if err != nil {
		return fmt.Errorf("Error opening R2 output: %v", err)
	}

	// batches are decoded in parallel and written in input order
	type job struct {
		batch []readPair
		res   chan result
	}
	jobs := make(chan job, threads*2)
	ordered := make(chan chan result, threads*2)
	var wg, wwg sync.WaitGroup
	var mu sync.Mutex
	var tot stats
	bxSet := map[string]struct{}{}

	for range threads {
		wg.Go(func() {
			w := &worker{decs: decoders, r1Start: r1Start, me: []byte(strings.ToUpper(meSeq)), maxMismatch: float64(maxMM), bxSet: map[string]struct{}{}, qxrx: qxrx}
			for j := range jobs {
				r := result{new(bytes.Buffer), new(bytes.Buffer)}
				for i := range j.batch {
					w.process(&j.batch[i], r.r1, r.r2)
				}
				j.res <- r
			}
			mu.Lock()
			tot.add(w.st)
			for k := range w.bxSet {
				bxSet[k] = struct{}{}
			}
			mu.Unlock()
		})
	}
	var writeErr error
	wwg.Go(func() {
		for res := range ordered {
			r := <-res
			if writeErr == nil {
				_, e1 := out1.buf.Write(r.r1.Bytes())
				_, e2 := out2.buf.Write(r.r2.Bytes())
				writeErr = errors.Join(e1, e2)
			}
		}
	})

	submit := func(b []readPair) {
		res := make(chan result, 1)
		ordered <- res
		jobs <- job{b, res}
	}
	batch := make([]readPair, 0, batchSize)
	var npairs int64
	for {
		fq1, err1 := next(rdr1)
		fq2, err2 := next(rdr2)
		eof1, eof2 := err1 == io.EOF, err2 == io.EOF
		if eof1 && eof2 {
			break
		}
		if err1 != nil && !eof1 {
			return fmt.Errorf("Error: R1 malformed near pair %d: %v", npairs+1, err1)
		}
		if err2 != nil && !eof2 {
			return fmt.Errorf("Error: R2 malformed near pair %d: %v", npairs+1, err2)
		}
		if eof1 != eof2 {
			return fmt.Errorf("Error: R1/R2 read count mismatch near pair %d (R1 EOF=%v, R2 EOF=%v)", npairs+1, eof1, eof2)
		}
		if !sameName(fq1.name, fq2.name) {
			return fmt.Errorf("Error: R1/R2 read names differ at pair %d: %q vs %q (files out of sync?)", npairs+1, fq1.name, fq2.name)
		}
		npairs++
		batch = append(batch, readPair{fq1, fq2})
		if len(batch) == batchSize {
			submit(batch)
			batch = make([]readPair, 0, batchSize)
		}
	}
	if len(batch) > 0 {
		submit(batch)
	}
	close(jobs)
	wg.Wait()
	close(ordered)
	wwg.Wait()

	if err := errors.Join(writeErr, out1.close(), out2.close()); err != nil {
		return fmt.Errorf("Error writing output: %v", err)
	}

	fmt.Fprintf(os.Stderr, "Total read pairs processed:          %d\n", tot.total)
	fmt.Fprintf(os.Stderr, "Total unique barcodes:               %d\n", len(bxSet))
	fmt.Fprintf(os.Stderr, "Total reads with barcode failure:    %d\n", tot.qcfail)
	fmt.Fprintf(os.Stderr, "Total reads with valid barcodes:     %d\n", tot.valid)
	fmt.Fprintf(os.Stderr, "Total reads with corrected barcodes: %d\n", tot.corrected)
	fmt.Fprintf(os.Stderr, "Discarded (ME sequence not found):   %d\n", tot.discarded)

	return nil
}
