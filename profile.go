package scramPkg

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"sync"
	"path/filepath"
)

//Details of an alignment for a discrete srna
type single_alignment struct {
	Seq   string  // Seq is an aligned read sequence
	timesAligned int // No of times the read has aligned
	Pos   int     // Pos is an aligned read position (from 5' fwd, starting at 1)
	Count float64 // Count is a normalised count
	Se    float64 // Se is the standard error
}

//collection for sorting
type single_alignments []*single_alignment

func (slice single_alignments) Len() int {
	return len(slice)
}
func (slice single_alignments) Less(i, j int) bool {
	return slice[i].Pos < slice[j].Pos
}
func (slice single_alignments) Swap(i, j int) {
	slice[i], slice[j] = slice[j], slice[i]
}

//ProfileNoSplit takes and alignment map and a sequence map as an input.  It returns a map of single alignments
//with a reference header as key and a single alignments struct as value.  Each single alignments struct is comprised of
//single_alignment structs (read seq, position, count, se).  The count for each read alignment is NOT split by the
//number of times a read aligns.
func ProfileNoSplit(alignment_map map[string]map[string][]int, seq_map map[string]*mean_se) map[string]*single_alignments {
	//switch to concurrent once functional
	srna_alignment_map := calc_times_read_aligns(alignment_map)
	profile_alignments_map := make(map[string]*single_alignments)
	for header, alignments := range alignment_map {
		var combined_alignments single_alignments
		for srna, positions := range alignments {
			for _, position := range positions {

				switch {
				case position > 0:
					alignment := single_alignment{srna, srna_alignment_map[srna],
						position,seq_map[srna].Mean, seq_map[srna].Se}
					combined_alignments = append(combined_alignments, &alignment)
				case position < 0:
					alignment := single_alignment{srna, srna_alignment_map[srna],
						0 - position, 0 - seq_map[srna].Mean, seq_map[srna].Se}
					combined_alignments = append(combined_alignments, &alignment)
				}
			}
		}
		sort.Sort(combined_alignments)
		profile_alignments_map[header] = &combined_alignments
	}
	return profile_alignments_map
}

type alignment_struct struct {
	header string
	alignments map[string][]int
}

type output_struct struct {
	header string
	combined_alignments *single_alignments
}


//ProfileSplit takes and alignment map and a sequence map as an input.  It returns a map of single alignments
//with a reference header as key and a single alignments struct as value.  Each single alignments struct is comprised of
//single_alignment structs (read seq, position, count, se).  The count for each read alignment is split by the
//number of times a read aligns.
func ProfileSplit(alignment_map map[string]map[string][]int, seq_map map[string]*mean_se) map[string]*single_alignments {
	alignmentsNo:=len(alignment_map)
	wg := &sync.WaitGroup{}
	wg.Add(alignmentsNo)

	srna_alignment_map := calc_times_read_aligns(alignment_map)



	alignments_for_goroutine := make(chan alignment_struct, alignmentsNo)
	output_from_goroutine := make(chan output_struct, alignmentsNo)

	for header, alignments := range alignment_map {
		//single_header_alignments :=alignment_struct{header,alignments}
		alignments_for_goroutine <- alignment_struct{header,alignments}

	}
	for a := 0; a < alignmentsNo; a++ {
		go profileSplitWorker(alignments_for_goroutine,output_from_goroutine,seq_map,srna_alignment_map,wg)
	}

	go func(cs chan output_struct, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(output_from_goroutine, wg)

	profile_alignments_map := make(map[string]*single_alignments)
	for single_output_struct := range output_from_goroutine {
		profile_alignments_map[single_output_struct.header]=single_output_struct.combined_alignments
	}

	return profile_alignments_map
}

func profileSplitWorker(alignments_for_goroutine chan alignment_struct, output_from_goroutine chan output_struct,
	seq_map map[string]*mean_se, srna_alignment_map map[string]int, wg *sync.WaitGroup){
	single_align := <- alignments_for_goroutine
	var combined_alignments single_alignments
	for srna, positions := range single_align.alignments {
		for _, position := range positions {
			split_count_mean := seq_map[srna].Mean / float64(srna_alignment_map[srna])
			split_se := seq_map[srna].Se / float64(srna_alignment_map[srna])
			switch {
			case position > 0:
				alignment := single_alignment{srna, srna_alignment_map[srna],
							      position, split_count_mean, split_se}
				combined_alignments = append(combined_alignments, &alignment)
			case position < 0:
				alignment := single_alignment{srna, srna_alignment_map[srna],
							      0 - position, 0 - split_count_mean, split_se}
				combined_alignments = append(combined_alignments, &alignment)
			}
		}
	}
	sort.Sort(combined_alignments)
	output_from_goroutine <- output_struct {single_align.header, &combined_alignments}
	wg.Done()
}

//ProfileToCsv writes the  den results to a csv file
func ProfileToCsv(profile_alignments_map map[string]*single_alignments, ref_slice []*header_ref, nt int, out_prefix string) {

	var strand string
	rows := [][]string{
		{"Header", "len", "sRNA","Position", "Strand", "Count", "Std. Err","Times aligned"},
	}
	for _, ref := range ref_slice {
		if alignments, ok := profile_alignments_map[ref.header]; ok {

			for _, alignment := range *alignments {

				switch {
				//TODO: fix this hack
				case alignment.Count<0:
					strand = "-"
					alignment.Count=0.0-alignment.Count
				default:
					strand = "+"
				}
				row := []string{ref.header, strconv.Itoa(len(ref.seq)),
					alignment.Seq, strconv.Itoa(alignment.Pos),
					strand,
					strconv.FormatFloat(alignment.Count, 'f', 3, 64),
					strconv.FormatFloat(alignment.Se, 'f', 8, 64),
					strconv.Itoa(alignment.timesAligned)}
				rows = append(rows, row)
			}
		}
	}
	out_file := out_prefix + "_" + strconv.Itoa(nt) + ".csv"
	out_dir := filepath.Dir(out_file)
	os.MkdirAll(out_dir,0777)
	f, err := os.Create(out_file)
	if err != nil {
		fmt.Println("Can't open csv file for writing")
		error_shutdown()
	}
	w := csv.NewWriter(f)
	w.WriteAll(rows)
	if err := w.Error(); err != nil {
		log.Fatalln("error writing csv:", err)
	}
}
