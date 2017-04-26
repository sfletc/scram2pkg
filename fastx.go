package scram2pkg

import (
	"bufio"
	"fmt"
	"github.com/montanaflynn/stats"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"
	"errors"
	"bytes"
	"compress/gzip"
)

func error_shutdown() {
	fmt.Println("\nExiting scram2")
	os.Exit(1)
}


//SeqLoad loads 1 or more small RNA seq. collapsed FASTA files (ie. from FASTX Toolkit).
//It returns a map with a read sequence as key and a mean_se struct (normalised read mean and standard error) as a value.
//Only reads present in all input files are returned.  Little format checking is  performed.  It is required
// that the input file is correctly formatted.
func SeqLoad(seq_files []string, file_type string, adapter string, min_len int, max_len int, min_count float64) map[string]*mean_se {
	t1 := time.Now()
	wg := &sync.WaitGroup{}
	no_of_files := len(seq_files)
	wg.Add(no_of_files)

	file_names := make(chan string, len(seq_files))
	for _, file_name := range seq_files {
		file_names <- file_name
	}
	close(file_names)

	srna_maps := make(chan map[string]float64, len(seq_files))

	for a := 0; a < len(seq_files); a++ {
		if file_type == "cfa" {
			go loadCfaFile(file_names, srna_maps, min_len, max_len, min_count, wg)
		} else if file_type == "fa" {
			go loadFastx(file_names, []byte(">"), adapter, srna_maps, min_len, max_len, min_count, wg)
		} else if file_type == "fq" {
			go loadFastx(file_names, []byte("@"), adapter, srna_maps, min_len, max_len, min_count, wg)
		}
	}
	go func(cs chan map[string]float64, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(srna_maps, wg)
	seq_map_all_counts := compileCounts(srna_maps)
	seq_map := calcMeanSe(seq_map_all_counts, no_of_files)
	t3 := time.Since(t1)
	fmt.Println("Read file set processed: ", t3)
	return seq_map
}

//Load a single collapsed read file and return map of read sequence as key and normalised RPMR count as value
func loadCfaFile(file_names chan string, srna_maps chan map[string]float64,
	min_len int, max_len int, min_count float64, wg *sync.WaitGroup) {
	t1 := time.Now()
	srna_map := make(map[string]float64)
	var count float64
	var total_count float64
	file_name := <-file_names
	f, err := os.Open(file_name)
	defer f.Close()
	if err !=nil {
		fmt.Println("\nCan't load collapsed read file " + file_name)
		os.Exit(1)
	}
	scanner := bufio.NewScanner(f)
	if file_name[len(file_name)-2:]=="gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + file_name)
			os.Exit(1)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	}
	seq_next:=false
	for scanner.Scan() {
		fasta_line := scanner.Text()
		if seq_next==true && len(fasta_line)==0 {
			fmt.Println("Read file format problem - blank line between header and sequence in " + file_name)
			error_shutdown()
		}
		if strings.HasPrefix(fasta_line,">") {
			header_line := strings.Split(fasta_line, "-")
			err := checkHeaderError(header_line, file_name)
			if err != nil{
				error_shutdown()
			}
			count, err = strconv.ParseFloat(strings.Split(fasta_line, "-")[1], 32)
			seq_next = true
		} else if count >= min_count && len(fasta_line) >= min_len && len(fasta_line) <= max_len {
			total_count += count
			srna_map[strings.ToUpper(fasta_line)] = count
			seq_next=false
		}
	}
	srna_map = rpmrNormalize(srna_map, total_count)
	srna_maps <- srna_map
	t2 := time.Since(t1)
	fmt.Println("Single read file "+ file_name+" loaded: ", t2)
	wg.Done()
}

//Load a single FASTA or FASTQ file and return map of read sequence as key and normalised RPMR count as value.
//Trim adapter from 3' end using up to 12 nt of 5' end of adapter as seed if required
func loadFastx(file_names chan string, first_char []byte, adapter string, srna_maps chan map[string]float64,
	min_len int, max_len int, min_count float64, wg *sync.WaitGroup) {
	t1 := time.Now()
	trim := false
	var seed string
	trim, seed = trimAdapter(adapter, trim, seed)

	srna_map := make(map[string]float64)
	var total_count float64
	file_name := <-file_names
	f, err := os.Open(file_name)
	if err != nil {
		fmt.Println("\nCan't load read file " + file_name)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if file_name[len(file_name)-2:] == "gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + file_name)
			os.Exit(1)
		}
		defer gz.Close()
		scanner = bufio.NewScanner(gz)
	}

	seq_next := false
	for scanner.Scan() {
		fasta_line := scanner.Bytes()
		switch {
		case bytes.Equal(fasta_line[:1], first_char):
			seq_next = true
		case seq_next == true && trim == true:
			srna_map, total_count, seq_next = addTrimmedRead(fasta_line, seed, min_len, max_len, srna_map, total_count, seq_next)
		case seq_next == true && len(fasta_line) >= min_len && len(fasta_line) <= max_len:
			srna_map, total_count, seq_next = addFullLengthRead(srna_map, fasta_line, total_count, seq_next)
		}
	}
	srna_map, total_count = removeReadsBelowMin(min_count, srna_map, total_count)
	srna_map = rpmrNormalize(srna_map, total_count)
	srna_maps <- srna_map
	t2 := time.Since(t1)
	fmt.Println("Single read file "+file_name+" loaded: ", t2)
	wg.Done()
}

//If adapter present, set seed for trimming (5' up to 12 nt of adapter sequence)
func trimAdapter(adapter string, trim bool, seed string) (bool, string) {
	if adapter != "nil" {
		trim = true
		switch {
		case len(adapter) < 12:
			seed = adapter
		default:
			seed = adapter[:11]
		}
	}
	return trim, seed
}

//Trim read and add it to the srna_map
func addTrimmedRead(fasta_line []byte, seed string, min_len int, max_len int, srna_map map[string]float64, total_count float64, seq_next bool) (map[string]float64, float64, bool) {
	read_slice := bytes.Split(fasta_line, []byte(seed))
	if len(read_slice) == 2 && len(read_slice[0]) >= min_len && len(read_slice[0]) <= max_len {
		if srna_count, ok := srna_map[string(read_slice[0])]; ok {
			srna_map[string(read_slice[0])] = srna_count + 1.0
			total_count += 1.0
		} else {
			srna_map[string(read_slice[0])] = 1.0
			total_count += 1.0
		}
	}
	seq_next = false
	return srna_map, total_count, seq_next
}

//Add full-length read to the srna map
func addFullLengthRead(srna_map map[string]float64, fasta_line []byte, total_count float64, seq_next bool) (map[string]float64, float64, bool) {
	if srna_count, ok := srna_map[string(fasta_line)]; ok {
		srna_map[string(fasta_line)] = srna_count + 1.0
		total_count += 1.0
	} else {
		srna_map[string(fasta_line)] = 1.0
		total_count += 1.0
	}
	seq_next = false
	return srna_map, total_count, seq_next
}

//Remove reads with count below the stated minimum for the srna_map
func removeReadsBelowMin(min_count float64, srna_map map[string]float64, total_count float64) (map[string]float64, float64) {
	if min_count > 1 {
		for srna, srna_count := range srna_map {
			if srna_count < min_count {
				delete(srna_map, srna)
				total_count -= srna_count
			}
		}
	}
	return srna_map, total_count
}

//Reads per million reads normalization of an input read library
func rpmrNormalize(srna_map map[string]float64, total_count float64) map[string]float64 {
	for srna, srna_count := range srna_map {
		srna_map[srna] = 1000000 * srna_count / total_count
	}
	return srna_map
}

//checks for error in collapsed fasta header
func checkHeaderError(header_line []string, file_name string) error {
	if len(header_line)<2 || len(header_line)>2 {
		return errors.New("\n"+ file_name+" is incorrectly formatted")
	}
	return nil
}

//compile_counts generates a map wit read seq as key and a slice of normalised counts for each read file
func compileCounts(srna_maps chan map[string]float64) map[string][]float64 {
	//map [srna:[count1,count2....], ...]
	seq_map_all_counts := make(map[string][]float64)
	for single_seq_map := range srna_maps {
		for srna, count := range single_seq_map {
			if _, ok := seq_map_all_counts[srna]; ok {
				seq_map_all_counts[srna] = append(seq_map_all_counts[srna], count)
			} else {
				seq_map_all_counts[srna] = []float64{count}
			}
		}
	}
	return seq_map_all_counts
}

//mean_se is a struct comprising a normalised mean and standard error for a read
type mean_se struct {
	Mean float64
	Se   float64
}

type read_counts struct {
	read string
	counts []float64
}

//calc_mean_se calculates the mean and standard error for each slice of counts
func calcMeanSe(seq_map_all_counts map[string][]float64, no_of_files int) map[string]*mean_se {
	seq_map := make(map[string]*mean_se)
	sqrt:=math.Sqrt(float64(no_of_files))
	for srna, counts := range seq_map_all_counts {
		if len(counts) < no_of_files{
			zeros := make([]float64,no_of_files-len(counts))
			counts=append(counts,zeros[:]...)
		}
		switch {
		case no_of_files > 1:
			counts_mean, _ := stats.Mean(counts)
			counts_stddev, _ := stats.StandardDeviationSample(counts)
			seq_map[srna] = &mean_se{counts_mean, counts_stddev / sqrt}
		default:
			seq_map[srna] = &mean_se{counts[0], 0.0}

		}
	}
	return seq_map
}

//mean_se is a struct comprising a normalised mean and standard error for a read
type header_ref struct {
	header string
	seq   string
}

//RefLoad loads a reference sequence DNA file (FASTA format).
//It returns a slice of header_ref structs (individual reference header and sequence).
func RefLoad(ref_file string) []*header_ref {
	t1 := time.Now()
	var ref_slice []*header_ref
	var single_header_ref *header_ref
	var header string
	var refSeq bytes.Buffer
	f, err := os.Open(ref_file)
	defer f.Close()
	if err != nil {
		fmt.Println("Problem opening fasta reference file "+ref_file)
		error_shutdown()
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fasta_line := scanner.Text()
		switch {
		case strings.HasPrefix(fasta_line, ">"):
			single_header_ref = &header_ref{header, refSeq.String()}
			ref_slice=append(ref_slice, single_header_ref)
			header=fasta_line[1:]
			refSeq.Reset()
		case len(fasta_line) != 0:
			refSeq.WriteString(strings.ToUpper(fasta_line))
		}
	}
	single_header_ref = &header_ref{header, refSeq.String()}
	ref_slice=append(ref_slice, single_header_ref)
	ref_slice=ref_slice[1:]

	t2 := time.Since(t1)
	fmt.Println("No. of reference sequences: ", len(ref_slice))
	fmt.Println("Reference file processed: ", t2)
	return ref_slice
}