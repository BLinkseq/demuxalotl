package main

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

const cfgPath = "../../templates/pheniqs.config.json"

// read builds a stagger-01 (7 T pad) R1: padded[7:58] + ME + tail, with the given A/B/C segments.
func read(t *testing.T, a, b, c string, q byte) (readPair, string) {
	t.Helper()
	const me = "AGATGTGTATAAGAGACAG"
	padded := "TTTTTTT" + "GGG" + a + "ACA" + b + "ACA" + c + "ACAGTA" // 7+3+12+3+12+3+12+6 = 58
	if len(padded) != 58 {
		t.Fatalf("bad test read, padded is %d", len(padded))
	}
	r1 := padded[7:] + me + strings.Repeat("ACGT", 20)
	r2 := strings.Repeat("TGCA", 20)
	qual := func(n int) []byte { return bytes.Repeat([]byte{q}, n) }
	return readPair{
		fastqRecord{[]byte("r1"), []byte(r1), qual(len(r1))},
		fastqRecord{[]byte("r1"), []byte(r2), qual(len(r2))},
	}, padded
}

func run(t *testing.T, p readPair) (string, stats) {
	t.Helper()
	decs, r1Start, err := loadConfig(cfgPath)
	if err != nil {
		t.Fatal(err)
	}
	w := &worker{decs: decs, r1Start: r1Start, me: []byte("AGATGTGTATAAGAGACAG"), maxMismatch: 2, bxSet: map[string]struct{}{}}
	var o1, o2 bytes.Buffer
	w.process(&p, &o1, &o2)
	return o1.String(), w.st
}

func TestExactAndCorrected(t *testing.T) {
	d, _, _ := loadConfig(cfgPath)
	a, b, c := d[1].seqs[0], d[2].seqs[1], d[3].seqs[2] // segment @01, @02, @03

	p, padded := read(t, a, b, c, 'I')
	out, st := run(t, p)
	if !strings.Contains(out, "\tVX:i:1\tBX:Z:A01C02B03D01\n") || st.qcfail != 0 || st.corrected != 0 {
		t.Fatalf("exact barcode: %q %+v", out, st)
	}
	// R1 sequence is padded[58:], the 6 bases after the barcode region plus ME-trimmed tail
	if want := strings.Repeat("ACGT", 20) + "\n"; !strings.Contains(out, "\n"+want+"+\n") {
		t.Fatalf("R1 trim: %q (padded %s)", out, padded)
	}

	// one low-quality mismatch in A is corrected
	mut := []byte(a)
	mut[3] = map[byte]byte{'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[3]]
	p, _ = read(t, string(mut), b, c, 'I')
	p.fq1.qual[3+3] = '#' // Q2 at the mutated base (A starts at padded index 10 -> raw index 3)
	out, st = run(t, p)
	if !strings.Contains(out, "BX:Z:A01C02B03D01\n") || st.corrected != 2 || st.qcfail != 0 {
		t.Fatalf("corrected barcode: %q %+v", out, st)
	}
}

func TestFailKeepsObserved(t *testing.T) {
	d, _, _ := loadConfig(cfgPath)
	// six high-quality mismatches in A: no confident call, observed (invalid) barcode kept
	mut := []byte(d[1].seqs[0])
	for i := range 6 {
		mut[i] = map[byte]byte{'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[i]]
	}
	p, _ := read(t, string(mut), d[2].seqs[1], d[3].seqs[2], 'I')
	out, st := run(t, p)
	if !strings.Contains(out, "\tVX:i:0\tBX:Z:A00C02B03D01\n") || st.qcfail != 2 || st.valid != 0 {
		t.Fatalf("failed barcode: %q %+v", out, st)
	}
}

func TestNoME(t *testing.T) {
	p, _ := read(t, strings.Repeat("A", 12), strings.Repeat("C", 12), strings.Repeat("G", 12), 'I')
	p.fq1.seq = bytes.Repeat([]byte("A"), 150)
	out, st := run(t, p)
	if out != "" || st.discarded != 1 {
		t.Fatalf("no ME: %q %+v", out, st)
	}
}

// pow10 must equal math.Pow(10, y) bit for bit: decode() decisions depend on it.
func TestPow10MatchesMathPow(t *testing.T) {
	check := func(y float64) {
		if a, b := pow10(y), math.Pow(10, y); math.Float64bits(a) != math.Float64bits(b) {
			t.Fatalf("pow10(%v)=%v, math.Pow=%v", y, a, b)
		}
	}
	for _, y := range []float64{0, math.Copysign(0, -1), -0.5, -1, -2, -0.25, -0.75, -1.5, -10, -100, -300, -320, -400} {
		check(y)
	}
	r := rand.New(rand.NewSource(1))
	for range 20_000_000 {
		check(-r.Float64() * 400)
	}
	for s := 0.0; s < 4000; s += 0.01 { // sums of Phred scores, as decode() builds them
		check(-0.1 * s)
	}
	for range 2_000_000 {
		check(-0.1 * (float64(r.Intn(60)*12) + float64(r.Intn(5))*10*math.Log10(4)))
	}
}

// next must return every record of a good file and an error, not io.EOF, for a truncated one.
func TestNextTruncated(t *testing.T) {
	seq.ValidateSeq = false
	const rec = "@r1 desc\nACGTN\n+\n@+III\n" // quality starts with @ and +, both legal
	count := func(s string) (int, error) {
		r, err := fastx.NewReaderFromIO(seq.DNAredundant, strings.NewReader(s), "")
		if err != nil {
			return 0, err
		}
		n := 0
		for {
			x, err := next(r)
			if err == io.EOF {
				return n, nil
			}
			if err != nil {
				return n, err
			}
			if string(x.name) != "r1" || string(x.seq) != "ACGTN" || string(x.qual) != "@+III" {
				t.Fatalf("bad record %q %q %q", x.name, x.seq, x.qual)
			}
			n++
		}
	}
	if n, err := count(rec + rec + rec); n != 3 || err != nil {
		t.Fatalf("good file: n=%d err=%v", n, err)
	}
	for _, cut := range []int{len("@r1 de"), len("@r1 desc\n"), len("@r1 desc\nACG"), len("@r1 desc\nACGTN\n"), len("@r1 desc\nACGTN\n+\n"), len("@r1 desc\nACGTN\n+\n@+I")} {
		if _, err := count(rec + rec[:cut]); err == nil {
			t.Errorf("truncation after %d bytes of the last record was accepted", cut)
		}
	}
}

// TestMain lets the tests run the real main() in a child process, so exit codes can be checked.
func TestMain(m *testing.M) {
	if os.Getenv("GIH_DEMUX_MAIN") == "1" {
		main()
		os.Exit(0)
	}
	os.Exit(m.Run())
}

func runMain(t *testing.T, r1, r2 string) (int, string) {
	t.Helper()
	return runMainStats(t, r1, r2, "")
}

// runMainStats is runMain with the stats file path set ("" = a temp file).
func runMainStats(t *testing.T, r1, r2, stats string) (int, string) {
	t.Helper()
	dir := t.TempDir()
	p1, p2 := filepath.Join(dir, "R1.fq"), filepath.Join(dir, "R2.fq")
	for p, s := range map[string]string{p1: r1, p2: r2} {
		if err := os.WriteFile(p, []byte(s), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if stats == "" {
		stats = filepath.Join(dir, "stats.txt")
	}
	cmd := exec.Command(os.Args[0], "--threads", "2", "--stats", stats,
		cfgPath, p1, p2, filepath.Join(dir, "o1.fq.gz"), filepath.Join(dir, "o2.fq.gz"))
	cmd.Env = append(os.Environ(), "GIH_DEMUX_MAIN=1")
	var stderr strings.Builder
	cmd.Stderr = &stderr
	err := cmd.Run()
	if ee, ok := err.(*exec.ExitError); ok {
		return ee.ExitCode(), stderr.String()
	}
	if err != nil {
		t.Fatal(err)
	}
	return 0, stderr.String()
}

func fastq(n int, seq string) string {
	var sb strings.Builder
	for i := range n {
		fmt.Fprintf(&sb, "@r%d\n%s\n+\n%s\n", i, seq, strings.Repeat("I", len(seq)))
	}
	return sb.String()
}

func TestMainExitCodes(t *testing.T) {
	r1 := fastq(20, strings.Repeat("ACGT", 12)+"AGATGTGTATAAGAGACAG"+strings.Repeat("ACGT", 10))
	r2 := fastq(20, strings.Repeat("TGCA", 20))
	if code, stderr := runMain(t, r1, r2); code != 0 {
		t.Fatalf("good input: exit %d: %s", code, stderr)
	}
	// both files cut inside their last record: counts match, so only the truncation can raise the error
	for _, tail := range []string{"@r20", "@r20\nACGT", "@r20\nACGT\n+\n", "@r20\nACGT\n+\nII"} {
		if code, stderr := runMain(t, r1+tail, r2+tail); code == 0 {
			t.Errorf("both files cut after %q exited 0: %s", tail, stderr)
		}
	}
	if code, stderr := runMain(t, r1, fastq(19, strings.Repeat("TGCA", 20))); code == 0 || !strings.Contains(stderr, "mismatch") {
		t.Errorf("read count mismatch: exit %d: %s", code, stderr)
	}
}

// A stats file that cannot be written (disk full) must fail the run, not exit 0 with a missing stats file.
func TestMainStatsWriteFailure(t *testing.T) {
	if _, err := os.Stat("/dev/full"); err != nil {
		t.Skip("no /dev/full on this platform")
	}
	r1 := fastq(20, strings.Repeat("ACGT", 12)+"AGATGTGTATAAGAGACAG"+strings.Repeat("ACGT", 10))
	r2 := fastq(20, strings.Repeat("TGCA", 20))
	if code, _ := runMainStats(t, r1, r2, "/dev/full"); code == 0 {
		t.Error("exited 0 although the stats file could not be written")
	}
}

func TestSameName(t *testing.T) {
	for _, c := range []struct {
		r1, r2 string
		want   bool
	}{
		{"read1", "read1", true},
		{"read1/1", "read1/2", true},
		{"read1/1", "read1", true},
		{"read1", "read1/2", true},
		{"read1", "read2", false},
		{"read1/1", "read2/2", false},
		{"read1", "read10", false},
		{"read1", "Read1", false},
		{"read1/1", "read1/1", false}, // R2 carrying a /1 suffix means the files are swapped or mislabeled
		{"read1/2", "read1/2", false},
		{"", "read1", false},
		{"", "", true},
	} {
		if got := sameName([]byte(c.r1), []byte(c.r2)); got != c.want {
			t.Errorf("sameName(%q, %q) = %v, want %v", c.r1, c.r2, got, c.want)
		}
	}
}

// Names are checked on every pair, not only the first, and the error says where.
func TestMainRejectsNameMismatch(t *testing.T) {
	seq1 := strings.Repeat("ACGT", 12) + "AGATGTGTATAAGAGACAG" + strings.Repeat("ACGT", 10)
	seq2 := strings.Repeat("TGCA", 20)
	build := func(n int, seq, suffix string, rename func(i int) string) string {
		var sb strings.Builder
		for i := range n {
			fmt.Fprintf(&sb, "@%s%s desc\n%s\n+\n%s\n", rename(i), suffix, seq, strings.Repeat("I", len(seq)))
		}
		return sb.String()
	}
	same := func(i int) string { return fmt.Sprintf("read%d", i) }

	if code, stderr := runMain(t, build(20, seq1, "/1", same), build(20, seq2, "/2", same)); code != 0 {
		t.Fatalf("matching names with /1 and /2 suffixes: exit %d: %s", code, stderr)
	}
	for _, bad := range []int{0, 7, 19} { // first, middle, last pair
		renamed := func(i int) string {
			if i == bad {
				return "other"
			}
			return same(i)
		}
		code, stderr := runMain(t, build(20, seq1, "", same), build(20, seq2, "", renamed))
		want := fmt.Sprintf("differ at pair %d", bad+1)
		if code == 0 || !strings.Contains(stderr, want) {
			t.Errorf("mismatch at pair %d: exit %d, stderr %q, want it to contain %q", bad+1, code, stderr, want)
		}
	}
	// the same reads in a different order is the classic out-of-sync case
	swapped := func(i int) string { return same((i + 1) % 20) }
	if code, _ := runMain(t, build(20, seq1, "", same), build(20, seq2, "", swapped)); code == 0 {
		t.Error("R2 shifted by one read exited 0")
	}
}

// CASAVA 1.8+ puts "1:N:0:INDEX" after a space; it must never take part in the name comparison,
// even though it differs between R1 and R2 (read number, filter flag, index).
func TestCasavaHeaders(t *testing.T) {
	seq.ValidateSeq = false
	r, err := fastx.NewReaderFromIO(seq.DNAredundant, strings.NewReader("@A00:1:FC:1:1101:10:20/1 1:N:0:ATGGAGA\nACGT\n+\nIIII\n"), "")
	if err != nil {
		t.Fatal(err)
	}
	rec, err := next(r)
	if err != nil || string(rec.name) != "A00:1:FC:1:1101:10:20/1" {
		t.Fatalf("name = %q, err %v: the CASAVA part must be dropped and /1 kept for sameName", rec.name, err)
	}

	seq1 := strings.Repeat("ACGT", 12) + "AGATGTGTATAAGAGACAG" + strings.Repeat("ACGT", 10)
	seq2 := strings.Repeat("TGCA", 20)
	build := func(seq, header string) string {
		var sb strings.Builder
		for i := range 20 {
			fmt.Fprintf(&sb, "@A00:1:FC:1:1101:%d:20%s\n%s\n+\n%s\n", i, header, seq, strings.Repeat("I", len(seq)))
		}
		return sb.String()
	}
	for name, h := range map[string][2]string{
		"CASAVA only":      {" 1:N:0:ATGGAGA", " 2:N:0:ATGGAGA"},
		"CASAVA, new idx":  {" 1:N:0:ATGGAGA", " 2:Y:23:CCCCCCC"},
		"suffix + CASAVA":  {"/1 1:N:0:ATGGAGA", "/2 2:N:0:TTTTTTT"},
		"suffix only":      {"/1", "/2"},
		"no suffix at all": {"", ""},
	} {
		if code, stderr := runMain(t, build(seq1, h[0]), build(seq2, h[1])); code != 0 {
			t.Errorf("%s: exit %d: %s", name, code, stderr)
		}
	}
	// a different read is still caught when the CASAVA parts are identical
	if code, _ := runMain(t, build(seq1, " 1:N:0:ATGGAGA"), strings.Replace(build(seq2, " 2:N:0:ATGGAGA"), ":1101:7:20", ":1101:8:20", 1)); code == 0 {
		t.Error("different read ids with identical CASAVA parts exited 0")
	}
}
