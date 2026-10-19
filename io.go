package main

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/klauspost/pgzip"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// ── FASTQ reader ──────────────────────────────────────────────────────────────

type fastqRecord struct{ name, seq, qual []byte }

type readPair struct{ fq1, fq2 fastqRecord }

// next reads one record from r. fastx reuses its record on every Read, and the
// pair is handed to another goroutine, so the fields are copied (one allocation).
// Qualities stay ASCII. Returns io.EOF at the clean end of the file.
func next(r *fastx.Reader) (fastqRecord, error) {
	rec, err := r.Read()
	if err != nil {
		return fastqRecord{}, err
	}
	if len(rec.Seq.Seq) == 0 { // fastx returns a header-only tail as an empty record, so a truncated file would pass silently
		return fastqRecord{}, fmt.Errorf("record %q has no sequence, is the file truncated?", rec.ID)
	}
	n, sEnd := len(rec.ID), len(rec.ID)+len(rec.Seq.Seq)
	buf := make([]byte, 0, sEnd+len(rec.Seq.Qual))
	buf = append(buf, rec.ID...)
	buf = append(buf, rec.Seq.Seq...)
	buf = append(buf, rec.Seq.Qual...)
	return fastqRecord{buf[:n:n], buf[n:sEnd:sEnd], buf[sEnd:]}, nil
}

// sameName reports whether two read names belong to the same pair: identical once a
// trailing /1 is dropped from R1's name and a trailing /2 from R2's.
func sameName(r1, r2 []byte) bool {
	return bytes.Equal(bytes.TrimSuffix(r1, []byte("/1")), bytes.TrimSuffix(r2, []byte("/2")))
}

// newFastqReader wraps path in a fastx reader. The file is decompressed here, with
// exactly gzBlocks pgzip blocks, so fastx sees plain text and never decompresses itself.
func newFastqReader(path string, gzBlocks int) (*fastx.Reader, error) {
	fh, err := openFastq(path, gzBlocks)
	if err != nil {
		return nil, err
	}
	return fastx.NewReaderFromIO(seq.DNAredundant, fh, "")
}

func openFastq(path string, gzBlocks int) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if !strings.HasSuffix(path, ".gz") {
		return f, nil
	}
	gz, err := pgzip.NewReaderN(f, 3<<20, gzBlocks)
	if err != nil {
		f.Close()
		return nil, err
	}
	return struct {
		io.Reader
		io.Closer
	}{gz, closers{gz, f}}, nil
}

type closers []io.Closer

func (c closers) Close() (err error) {
	for _, x := range c {
		if e := x.Close(); e != nil {
			err = e
		}
	}
	return
}

// ----- FASTQ Writer ----------------------
func writeFastq(b *bytes.Buffer, name []byte, dir, bx string, vx int, seq, qual []byte) {
	b.WriteByte('@')
	b.Write(name)
	b.WriteString(dir)
	b.WriteString("\tVX:i:")
	b.WriteByte(byte('0' + vx))
	b.WriteString("\tBX:Z:")
	b.WriteString(bx)
	b.WriteByte('\n')
	b.Write(seq)
	b.WriteString("\n+\n")
	b.Write(qual)
	b.WriteByte('\n')
}

type gzOut struct {
	f   *os.File
	gz  *pgzip.Writer
	buf *bufio.Writer
}

func newGzOut(path string, threads int) (*gzOut, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	gz, err := pgzip.NewWriterLevel(f, 4)
	if err != nil {
		f.Close()
		return nil, err
	}
	gz.SetConcurrency(2<<20, threads)
	return &gzOut{f, gz, bufio.NewWriterSize(gz, 2<<20)}, nil
}

func (o *gzOut) close() error {
	return errors.Join(o.buf.Flush(), o.gz.Close(), o.f.Close())
}
