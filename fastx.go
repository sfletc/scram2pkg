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
//Only reads present in all input files are returned.  No FASTA format checking is  performed.  It is required
// that the input file is correctly formatted.
func SeqLoad(seq_files []string, file_type string, min_len int, max_len int, min_count float64) map[string]*mean_se {
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
			go load_cfa_file_go(file_names, srna_maps, min_len, max_len, min_count, wg)
		} else if file_type == "fa" {
			go load_fx_file_go(file_names, []byte(">"), srna_maps, min_len, max_len, min_count, wg)
		} else if file_type == "fq" {
			go load_fx_file_go(file_names, []byte("@"), srna_maps, min_len, max_len, min_count, wg)
		}
	}
	go func(cs chan map[string]float64, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(srna_maps, wg)
	fmt.Println(srna_maps)
	seq_map_all_counts := compile_counts(srna_maps)
	seq_map := calc_mean_se(seq_map_all_counts, no_of_files)
	t3 := time.Since(t1)
	fmt.Println("Read file set processed: ", t3)
	return seq_map
}

//Load a single collapsed read file and return map of read sequence as key and normalised RPMR count as value
func load_cfa_file_go(file_names chan string, srna_maps chan map[string]float64,
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
	seq_next:=false
	for scanner.Scan() {
		fasta_line := scanner.Text()
		if seq_next==true && len(fasta_line)==0 {
			fmt.Println("Read file format problem - blank line between header and sequence in " + file_name)
			error_shutdown()
		}
		if strings.HasPrefix(fasta_line,">") {
			header_line := strings.Split(fasta_line, "-")
			err := check_header_error(header_line, file_name)
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
	for srna, srna_count := range srna_map {
		srna_map[srna] = 1000000 * srna_count / total_count
	}
	srna_maps <- srna_map
	t2 := time.Since(t1)
	fmt.Println("Single read file "+ file_name+" loaded: ", t2)
	wg.Done()
}


func load_fx_file_go(file_names chan string, first_char []byte, srna_maps chan map[string]float64,
	min_len int, max_len int, min_count float64, wg *sync.WaitGroup) {
	t1 := time.Now()
	srna_map := make(map[string]float64)

	//var count float64
	var total_count float64
	file_name := <-file_names
	f, err := os.Open(file_name)
	if err !=nil {
		fmt.Println("\nCan't load read file " + file_name)
		os.Exit(1)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if file_name[len(file_name)-2:]=="gz" {
		gz, err := gzip.NewReader(f)
		if err != nil {
			fmt.Println("\nCan't decompress read file " + file_name)
			os.Exit(1)
		}
		scanner = bufio.NewScanner(gz)
	}

	seq_next:=false
	for scanner.Scan() {
		fasta_line := scanner.Bytes()
		//fmt.Println(fasta_line,fasta_line[0:1],string(fasta_line[0]))
		switch {
		//case seq_next==true && len(fasta_line)==0:
		//	fmt.Println("Read file format problem - blank line between header and sequence in " + file_name)
		//	error_shutdown()

		case bytes.Equal(fasta_line[:1],first_char):
			seq_next=true
		case seq_next==true && len(fasta_line) >= min_len && len(fasta_line) <= max_len:
			if srna_count, ok := srna_map[string(fasta_line)]; ok{
				srna_map[string(fasta_line)] = srna_count+1.0
				total_count += 1.0
			} else {
				srna_map[string(fasta_line)]=1.0
				total_count += 1.0
			}
			seq_next=false
		}
	}
	//TODO: write test function for this
	if min_count > 1 {
		for srna, srna_count := range srna_map {
			if srna_count < min_count {
				delete(srna_map, srna)
				total_count-=srna_count
			}
		}
	}
	for srna, srna_count := range srna_map {
		srna_map[srna] = 1000000 * srna_count / total_count
	}
	//fmt.Println(srna_map)
	srna_maps <- srna_map
	t2 := time.Since(t1)
	fmt.Println("Single read file "+ file_name+" loaded: ", t2)
	wg.Done()
}





func check_header_error(header_line []string, file_name string) error {
	if len(header_line)<2 || len(header_line)>2 {
		return errors.New("\n"+ file_name+" is incorrectly formatted")
	}
	return nil
}

//compile_counts generates a map wit read seq as key and a slice of normalised counts for each read file
func compile_counts(srna_maps chan map[string]float64) map[string][]float64 {
	//map [srna:[count1,count2....], ...]
	seq_map_all_counts := make(map[string][]float64)
	//no_of_files := float64(len(seq_files))
	for single_seq_map := range srna_maps {
		for srna, count := range single_seq_map {
			if _, ok := seq_map_all_counts[srna]; ok {
				//counts += count / no_of_files
				seq_map_all_counts[srna] = append(seq_map_all_counts[srna], count)
			} else {
				counts := []float64{count}
				//counts = count / no_of_files
				seq_map_all_counts[srna] = counts
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

//calc_mean_se calculates the mean and standard error for each slice of counts
func calc_mean_se(seq_map_all_counts map[string][]float64, no_of_files int) map[string]*mean_se {
	seq_map := make(map[string]*mean_se)
	//TODO: THIS MEANS ONLY READS THAT PRESENT IN ALL FILES ARE RETAINED - MAY WANT AN OPTION ABOUT THIS!!
	for srna, counts := range seq_map_all_counts {
		switch {
		case len(counts) == no_of_files && no_of_files > 1:
			counts_mean, _ := stats.Mean(counts)
			counts_stddev, _ := stats.StandardDeviationSample(counts)
			counts_stderr := counts_stddev / math.Sqrt(float64(no_of_files))
			read_stats := &mean_se{counts_mean, counts_stderr}
			seq_map[srna] = read_stats
		case len(counts) == no_of_files && no_of_files == 1:
			read_stats := &mean_se{counts[0], 0.0}
			seq_map[srna] = read_stats
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

//func check_dna(line string){
//	dna := map[rune]bool {
//		'A':true,
//		'T':true,
//		'C':true,
//		'G':true,
//		'N':true,
//	}
//	runes := []rune(line)
//	for i:=0;i<=len(runes)-1;i++{
//		if _, ok := dna[runes[i]]; ok {
//			continue
//		} else {
//			fmt.Println("\nNon-DNA character in reference")
//			fmt.Println()
//			error_shutdown()
//		}
//	}
//}