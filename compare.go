package scram2pkg

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"path/filepath"
)

//CompareNoSplitCounts takes and alignment map and returns a map with the ref_header as key and the
//mean_se (mean and standard error) of aligned reads for that ref seq as value.  Read counts are NOT split by the number
//of times a read aligns to all reference sequences.
func CompareNoSplitCounts(alignment_map map[string]map[string][]int, seq_map map[string]*mean_se) map[string]mean_se {
	cdp_alignment_map := make(map[string]mean_se)
	for header, alignment := range alignment_map {
		var counts float64
		var counts_err []float64 //error_shutdown this
		for srna, pos := range alignment {
			counts += seq_map[srna].Mean * float64(len(pos))
			counts_err = append(counts_err, seq_map[srna].Se*float64(len(pos)))
		}
		var err_sq_sum float64
		for _, errs := range counts_err {
			err_sq_sum += errs * errs
		}
		cdp_alignment_map = calc_header_mean_se(counts_err, counts, cdp_alignment_map, header)
	}
	return cdp_alignment_map
}

//CompareSplitCounts takes and alignment map and returns a map with the ref_header as key and the
//mean_se (mean and standard error) of aligned reads for that ref seq as value.  Read counts are split by the number
//of times a read aligns to all reference sequences.
func CompareSplitCounts(alignment_map map[string]map[string][]int, seq_map map[string]*mean_se) map[string]mean_se {
	cdp_alignment_map := make(map[string]mean_se)
	//Calc. no. of times each read aligns
	srna_alignment_map := calc_times_read_aligns(alignment_map)
	for header, alignment := range alignment_map {
		var counts float64
		var counts_err []float64
		for srna, pos := range alignment {
			counts += seq_map[srna].Mean * float64(len(pos)) / float64(srna_alignment_map[srna])
			// should be err=sqrt(x^2/n) for each alignment ??
			// They are perfectly correlated (dependent), so maybe not?
			counts_err = append(counts_err,
				seq_map[srna].Se*float64(len(pos))/float64(srna_alignment_map[srna]))
			//counts_err = append(counts_err, (math.Sqrt((seq_map[srna].Se*seq_map[srna].Se)/
			// float64(srna_alignment_map[srna])))*float64(len(pos)))
		}
		cdp_alignment_map = calc_header_mean_se(counts_err, counts, cdp_alignment_map, header)
	}
	return cdp_alignment_map
}

//Calculates the number of times an aligned read aligns
func calc_times_read_aligns(alignment_map map[string]map[string][]int) map[string]int {
	srna_alignment_map := make(map[string]int)
	for _, alignment := range alignment_map {
		for srna, pos := range alignment {
			if _, ok := srna_alignment_map[srna]; ok {
				srna_alignment_map[srna] += len(pos)
			} else {
				srna_alignment_map[srna] = len(pos)
			}
		}
	}
	return srna_alignment_map
}

//Calculates the mean  and se of alignments to a header
func calc_header_mean_se(counts_err []float64, counts float64, cdp_alignment_map map[string]mean_se,
	header string) map[string]mean_se {
	var err_sq_sum float64
	for _, errs := range counts_err {
		err_sq_sum += errs * errs
	}
	ref_align_err := math.Sqrt(err_sq_sum)
	cdp_counts_and_err := mean_se{counts, ref_align_err}
	cdp_alignment_map[header] = cdp_counts_and_err
	return cdp_alignment_map
}

//Compare combines individual alignments for set sets of sequences (treatments).  It returns a map of ref header
//as key and a slice of set 1 mean/se and set2 mean/se as value.
func Compare(counts_map_1 map[string]mean_se, counts_map_2 map[string]mean_se) map[string][]float64 {
	cdp_final_map := make(map[string][]float64)
	//TODO: This only includes if BOTH maps have a non-zero count alignment - perhaps should mod?
	for header, count_stats := range counts_map_1 {
		if count_stats_2, ok := counts_map_2[header]; ok {
			cdp_final_map[header] = append(cdp_final_map[header], count_stats.Mean,
				count_stats.Se, count_stats_2.Mean, count_stats_2.Se)
		}
	}
	return cdp_final_map
}

func MirnaCompare(mirna_alignment_map_1 map[string]*mean_se_dup, mirna_alignment_map_2 map[string]*mean_se_dup, noSplit bool) map[string][]float64 {
	cdp_final_map := make(map[string][]float64)
	for header, count_stats := range mirna_alignment_map_1 {
		if count_stats_2, ok := mirna_alignment_map_2[header]; ok {
			switch {
			case noSplit == true:
				cdp_final_map[header] = append(cdp_final_map[header], count_stats.mean_se.Mean,
					count_stats.mean_se.Se, count_stats_2.mean_se.Mean, count_stats_2.mean_se.Se)
			default:
				cdp_final_map[header] = append(cdp_final_map[header], count_stats.mean_se.Mean/count_stats.dup,
					count_stats.mean_se.Se/count_stats.dup, count_stats_2.mean_se.Mean/count_stats_2.dup, count_stats_2.mean_se.Se/count_stats_2.dup)
			}
		}
	}
	return cdp_final_map
}



//CompareToCsv writes the output to a csv file.
func CompareToCsv(cdp_alignment_map map[string][]float64, nt int, out_prefix string) {
	alignments := [][]string{
		{"Header", "Mean count 1", "Std. err 1", "Mean count 2", "Std. err 2"},
	}
	for header, count_stats := range cdp_alignment_map {
		alignment := []string{header,
			strconv.FormatFloat(count_stats[0], 'f', 3, 64),
			strconv.FormatFloat(count_stats[1], 'f', 8, 64),
			strconv.FormatFloat(count_stats[2], 'f', 3, 64),
			strconv.FormatFloat(count_stats[3], 'f', 8, 64)}
		alignments = append(alignments, alignment)
	}
	var out_file string
	switch{
	case nt>0:
		out_file = out_prefix + "_" + strconv.Itoa(nt) + ".csv"
	default:
		out_file = out_prefix + "_miR.csv"
	}

	out_dir := filepath.Dir(out_file)
	os.MkdirAll(out_dir,0777)
	f, err := os.Create(out_file)
	if err != nil {
		fmt.Println("Can't save to file")
		error_shutdown()
	}
	w := csv.NewWriter(f)
	w.WriteAll(alignments)
	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
}
