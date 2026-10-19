package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"runtime"
)

func filecheck(filename string) error {
	//var err error
	if err := fileExists(filename); err != nil {
		return err
	}
	if err := fileReadable(filename); err != nil {
		return err
	}
	return nil
}

// check if file exists. Returns a true if it does, otherwise false.
func fileExists(filename string) error {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return fmt.Errorf("\033[33;1m%s\033[0m does not exist", filename)
	}
	if err != nil {
		return fmt.Errorf("\033[33;1m%s\033[0m: %w", filename, err)
	}
	if info.IsDir() {
		return fmt.Errorf("\033[33;1m%s\033[0m is a directory", filename)
	}
	return nil
}

func fileReadable(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("\033[33;1m%s\033[0m does not have read persmissions.\n", filename)
	}
	file.Close()
	return nil
}

// check if executable 'e' is in path, returns true if it is, otherwise false
func checkIfExecInPath(e string) bool {
	// ignore output, only use error
	_, err := exec.LookPath(e)
	// return true if error nil, else false
	return err == nil
}

// Parse the thread count provided to `threads` at the CLI and clamp it between 1 and the system max
func safethreads(threads int) int {
	maxCores := runtime.NumCPU()
	// clamp between 1 and max system threads
	if threads > maxCores {
		fmt.Fprintf(
			os.Stderr,
			"Threads requested (%v) exceed total available (%v), scaling down to maximum available threads.\n",
			threads, maxCores,
		)
	}
	return min(maxCores, max(threads, 1))
}
