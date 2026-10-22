// Copyright © 2026 Pavel Dimens | Github: pdimens
package cmd

import (
	"demuxolotl/demux"
	"fmt"
	"os"
	"runtime"
	"strings"

	"github.com/spf13/cobra"
)

var VERSION string = "1.0"

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:                   "demuxolotl [options] prefix R1.fq[.gz] R2.fq[.gz]",
	Version:               VERSION,
	Short:                 "BLink-seq FASTQ barcode demultiplexer",
	Long:                  "Demultiplex linked-read barcodes from raw FASTQ BLink-seq data for one sample. Uses an embedded schema by default; if a custom decoding configuration is required, include a Pheniqs-style decoder JSON.",
	DisableFlagsInUseLine: true,
	SilenceUsage:          true,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			fmt.Printf("%s\n", cmd.UsageString())
			return fmt.Errorf("please provide positional arguments: an output file prefix and input FASTQ files")
		}
		if err := cobra.ExactArgs(3)(cmd, args); err != nil {
			return err
		}
		if err := filecheck(args[1]); err != nil {
			return err
		}
		if err := filecheck(args[2]); err != nil {
			return err
		}
		if _, err := cmd.Flags().GetInt("max-mismatch"); err != nil {
			return err
		}
		me_seq, err := cmd.Flags().GetString("me")
		if err != nil {
			return err
		}
		schema, err := cmd.Flags().GetString("schema")
		if schema != "use embedded" {
			if doesExist := fileExists(schema); doesExist != nil {
				return doesExist
			}
		}

		upperME := strings.ToUpper(me_seq)
		for _, char := range upperME {
			switch char {
			case 'A':
			case 'T':
			case 'C':
			case 'G':
			default:
				return fmt.Errorf("Invalid (non-nucleotide) character in ME sequence: %c\n%v", char, me_seq)
			}
		}
		return nil
	},
	RunE: func(cmd *cobra.Command, args []string) error {
		threads, _ := cmd.Flags().GetInt("threads")
		qxrx, _ := cmd.Flags().GetBool("qx-rx")
		maxmiss, _ := cmd.Flags().GetInt("max-mismatch")
		me_seq, _ := cmd.Flags().GetString("me")
		schema, _ := cmd.Flags().GetString("schema")
		threads = safethreads(threads)
		runtime.GOMAXPROCS(safethreads(threads))
		return demux.Demux(args[1], args[2], args[0], me_seq, schema, maxmiss, qxrx, threads)
	},
}

func init() {
	rootCmd.Flags().IntP("threads", "t", 4, "Threads to use")
	rootCmd.Flags().StringP("me", "m", "AGATGTGTATAAGAGACAG", "ME sequence to search for")
	rootCmd.Flags().IntP("max-mismatch", "x", 2, "Maximum allowable ME mismatch score")
	rootCmd.Flags().BoolP("qx-rx", "q", false, "Output tags QX:Z (barcode quality) and RX:Z (corrected barcode)")
	rootCmd.Flags().StringP("schema", "s", "use embedded", "Pheniqs-style JSON schema file with the barcodes and decoder settings")

	//rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}
